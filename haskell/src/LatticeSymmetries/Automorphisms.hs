{-# LANGUAGE DuplicateRecordFields #-}
{-# LANGUAGE OverloadedRecordDot #-}

module LatticeSymmetries.Automorphisms
  ( Hypergraph (..)
  , PermutationGroup (..)
  , mkPermutationGroup
  , emptyPermutationGroup
  , nullPermutationGroup
  , MultiplicationTable (..)
  , mkMultiplicationTable
  , shrinkMultiplicationTable
  , hypergraphAutomorphisms
  , isAutomorphism
  )
where

import Control.Exception (assert)
import Control.Monad.ST.Strict (runST)
import Data.Aeson (ToJSON)
import Data.List qualified
import Data.List.NonEmpty qualified as NonEmpty
import Data.Map.Strict qualified as Map
import Data.Set qualified as Set
import Data.Validity
import Data.Vector qualified as B
import Data.Vector.Algorithms.Search (binarySearch)
import Data.Vector.Generic ((!))
import Data.Vector.Generic qualified as G
import Data.Vector.Generic.Mutable qualified as GM
import Data.Vector.Unboxed qualified as U
import Data.Vector.Unboxed.Mutable qualified as UM
import Foreign (Ptr, castPtr)
import Foreign.C (CInt, CString)
import GHC.Records (HasField (..))
import LatticeSymmetries.Context
import LatticeSymmetries.Dense
import LatticeSymmetries.FFI (ls_hs_destroy_object, newCobject)
import LatticeSymmetries.Permutation
import LatticeSymmetries.Utils (MutablePtr, sortVectorBy)
import Prelude hiding (group, identity, permutations, second, toList)

data Hypergraph a = Hypergraph {vertices :: !(Set a), hyperedges :: !(Set (Set a))}
  deriving stock (Eq, Ord, Show)

newtype MultiplicationTable = MultiplicationTable {unMultiplicationTable :: DenseMatrix U.Vector Int}
  deriving stock (Show, Eq, Generic)
  deriving anyclass (NFData)
  deriving newtype (ToJSON)

-- | Permutation group
data PermutationGroup = PermutationGroup {permutations :: !(B.Vector Permutation), table :: MultiplicationTable}
  deriving stock (Show, Eq, Generic)
  deriving anyclass (NFData, ToJSON)

instance HasField "size" PermutationGroup Int where
  getField g = G.length g.permutations

mkPermutationGroup :: G.Vector v Permutation => v Permutation -> PermutationGroup
mkPermutationGroup p = PermutationGroup (G.convert p) (mkMultiplicationTable p)

emptyPermutationGroup :: PermutationGroup
emptyPermutationGroup = PermutationGroup G.empty (MultiplicationTable (DenseMatrix 0 0 G.empty))

-- | Check whether a permutation group is empty
nullPermutationGroup :: PermutationGroup -> Bool
nullPermutationGroup g = G.null g.permutations

mkMultiplicationTable :: G.Vector v Permutation => v Permutation -> MultiplicationTable
mkMultiplicationTable ps =
  MultiplicationTable . generateDenseMatrix n n $ \i j ->
    binarySearch' ps (ps ! i <> ps ! j)
  where
    !n = G.length ps

shrinkMultiplicationTable :: G.Vector v Int => MultiplicationTable -> v Int -> MultiplicationTable
shrinkMultiplicationTable (MultiplicationTable matrix) indices =
  MultiplicationTable . generateDenseMatrix n n $ \i j ->
    binarySearch' indices $ indexDenseMatrix matrix (indices ! i, indices ! j)
  where
    n = G.length indices

binarySearch' :: (HasCallStack, G.Vector v a, Ord a) => v a -> a -> Int
binarySearch' !v !x = runST $ do
  mv <- G.unsafeThaw v
  index <- binarySearch mv x
  unless (v ! index == x) $ error "element not found"
  pure index

data SearchTree a r
  = SearchBranch !a [SearchTree a r]
  | SearchLeaf !(Maybe r)
  deriving stock (Eq, Ord, Show, Functor)

newtype Partitioning a = Partitioning (NonEmpty (U.Vector a))
  deriving stock (Show, Eq)

allSingletons :: UM.Unbox a => Partitioning a -> Bool
allSingletons (Partitioning xs) = all ((== 1) . G.length) xs

distancePartition
  :: forall a
   . (UM.Unbox a, Ord a, Show a)
  => Hypergraph a
  -- ^ The hypergraph to partition
  -> a
  -- ^ A selected vertex
  -> Partitioning a
  -- ^ Groups of vertices by their distance to the selected vertex
distancePartition g v0 =
  Partitioning
    . NonEmpty.fromList
    . fmap (G.fromList . Set.toAscList)
    $ go (Set.singleton v0) (Set.toList g.hyperedges)
  where
    go :: Set a -> [Set a] -> [Set a]
    go !seen hyperedges
      | Set.null boundary = let r = g.vertices Set.\\ seen in [r | not (Set.null r)]
      | otherwise = boundary : go (Set.union seen combined) others
      where
        shouldIncludeHyperedge = any (`Set.member` seen) . Set.toList
        (toBeIncluded, others) = Data.List.partition shouldIncludeHyperedge hyperedges
        !combined = foldl' Set.union Set.empty toBeIncluded
        !boundary = combined Set.\\ seen

intersectSorted :: (G.Vector v a, Ord a) => v a -> v a -> v a
intersectSorted a b = runST $ do
  out <- GM.unsafeNew $ min (G.length a) (G.length b)
  let go !size !i !j
        | i < G.length a && j < G.length b =
            let x = G.unsafeIndex a i
                y = G.unsafeIndex b j
             in case compare x y of
                  EQ -> GM.unsafeWrite out size x >> go (size + 1) (i + 1) (j + 1)
                  LT -> go size (i + 1) j
                  GT -> go size i (j + 1)
        | otherwise = pure size
  size <- go 0 0 0
  G.unsafeFreeze (GM.take size out)
{-# SCC intersectSorted #-}

intersectAllToAll :: (UM.Unbox a, Ord a) => Partitioning a -> Partitioning a -> Partitioning a
intersectAllToAll (Partitioning p1) (Partitioning p2) =
  Partitioning
    . NonEmpty.fromList
    $ intersectSorted
      <$> NonEmpty.toList p1
      <*> NonEmpty.toList p2

createMappings :: U.Unbox a => Partitioning a -> Partitioning a -> [Mapping a]
createMappings (Partitioning (NonEmpty.toList -> src)) (Partitioning (NonEmpty.toList -> tgt)) =
  G.toList $ G.zipWith Mapping (G.concat src) (G.concat tgt)

pickOne :: (UM.Unbox a, Ord a) => Map a (Partitioning a) -> Partitioning a -> (a, Partitioning a)
pickOne dps (Partitioning (v :| rest))
  | G.null v = error "pickOne expects a normalized Partitioning"
  | otherwise = (x, intersectAllToAll p (dps Map.! x))
  where
    !x = G.head v
    p = Partitioning (G.tail v :| rest)

pickAll :: (UM.Unbox a, Ord a) => Map a (Partitioning a) -> Partitioning a -> [(a, Partitioning a)]
pickAll dps (Partitioning (xss :| rest)) = do
  (x, xs) <- picks (G.toList xss)
  let p = Partitioning (G.fromList xs :| rest)
  pure (x, intersectAllToAll p (dps Map.! x))
  where
    picks :: [a] -> [(a, [a])]
    picks [] = []
    picks (x : xs) = (x, xs) : [(y, x : ys) | (y, ys) <- picks xs]

normalizeWhenCompatible :: (UM.Unbox a, Ord a) => Partitioning a -> Partitioning a -> Maybe (Partitioning a, Partitioning a)
normalizeWhenCompatible (Partitioning (NonEmpty.toList -> src)) (Partitioning (NonEmpty.toList -> tgt))
  | fmap G.length src == fmap G.length tgt =
      let normalize = Partitioning . NonEmpty.fromList . filter (not . G.null)
       in -- (src', tgt') = unzip $ Data.List.sort $ zip (filter (not . null) src) (filter (not . null) tgt)
          Just (normalize src, normalize tgt)
  | otherwise = Nothing
{-# SCC normalizeWhenCompatible #-}

-- | Generate a SearchTree of graph automorphisms using distance partition
autsSearchTree :: Hypergraph Int -> SearchTree [Mapping Int] Permutation
autsSearchTree g = dfs [] (Partitioning (vs :| [])) (Partitioning (vs :| []))
  where
    dfs :: [Mapping Int] -> Partitioning Int -> Partitioning Int -> SearchTree [Mapping Int] Permutation
    dfs !mappings !srcPart !trgPart
      -- check compatibility at the final step
      | allSingletons srcPart =
          -- isCompatible xys = isAutomorphism g xys
          assert (allSingletons trgPart) $
            let p =
                  permutationFromMappings (Just (Set.size g.vertices)) $
                    mappings
                      <> createMappings srcPart trgPart
             in SearchLeaf $ if isAutomorphism g p then Just p else Nothing
      | otherwise =
          let (!x, !srcPart') = pickOne dps srcPart
              branches = do
                (!y, !trgPart') <- pickAll dps trgPart
                case normalizeWhenCompatible srcPart' trgPart' of
                  Just (srcPart'', trgPart'') -> pure $ dfs (Mapping x y : mappings) srcPart'' trgPart''
                  Nothing -> []
           in SearchBranch mappings branches
    vs = G.fromList $ Set.toAscList g.vertices
    !dps = Map.fromAscList $ (\v -> (v, distancePartition g v)) <$> G.toList vs
{-# SCC autsSearchTree #-}

transversalGeneratingSet :: SearchTree a Permutation -> B.Vector Permutation
transversalGeneratingSet = B.fromList . go
  where
    go (SearchLeaf (Just !p)) = [p | not (isIdentityPermutation p)]
    go (SearchBranch _ (t : ts)) = concatMap (force . take 1 . go) ts <> go t
    go _ = []

groupFromTransversalGeneratingSet :: B.Vector Permutation -> PermutationGroup
groupFromTransversalGeneratingSet tgs =
  mkPermutationGroup
    . sortVectorBy compare
    . B.fromList
    $ Data.List.foldr1 (flip (<>))
      <$> sequence transversals
  where
    k = (G.head tgs).length
    transversals =
      fmap ((identityPermutation k :) . B.toList . fmap fst) $
        B.groupBy (\g h -> snd g == snd h) $
          sortVectorBy (comparing snd) $
            fmap (\p -> (p, minimalSupport p)) tgs

hypergraphAutomorphisms :: (Permutation -> Bool) -> Hypergraph Int -> PermutationGroup
hypergraphAutomorphisms p =
  (mkPermutationGroup . G.filter p . (.permutations))
    . groupFromTransversalGeneratingSet
    . transversalGeneratingSet
    . autsSearchTree

isAutomorphism :: Hypergraph Int -> Permutation -> Bool
isAutomorphism g p = g.hyperedges == Set.map transformEdge g.hyperedges
  where
    i = unPermutation p
    transformEdge = Set.map (i G.!)
