{-# LANGUAGE DuplicateRecordFields #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE OverloadedRecordDot #-}

module LatticeSymmetries.Group
  ( PermutationGroup (..)
  , fromGenerators
  , pgLength
  , Symmetry
  , mkSymmetry
  , symmetrySector
  , symmetryPermutation
  , symmetryPhase
  , getPeriodicity
  , Symmetries (..)
  , mkSymmetries
  , mkSymmetriesFromRepresentation
  , areSymmetriesReal
  , nullSymmetries
  , emptySymmetries
  , symmetriesGetNumberBits

    -- ** Low-level interface for FFI
  , newCpermutation_group
  , destroyCpermutation_group

    -- * Automorphisms
  , Hypergraph (..)
  , hypergraphAutomorphisms
  , mkMultiplicationTable
  , abelianSubgroup
  , Representation (..)
  , groupRepresentations
  , cyclicGraph
  , cyclicGraph3
  , rectangularGraph
  -- , distancePartition
  -- , autsSearchTree
  -- , naiveExtractLeaves
  -- , isAutomorphism
  -- , transversalGeneratingSet
  -- , groupFromTransversalGeneratingSet
  , MultiplicationTable (..)
  , AbelianSubsetHistory (..)
  -- , abelianSubset
  , GroupElement (..)
  , shrinkMultiplicationTable
  , getGroupElements
  , selectPrimeElements
  , groupGenerators
  , groupElementCycle
  , groupElementOrder
  , Coset (..)
  , abelianization
  , commutatorSubgroup
  )
where

import Control.Exception (assert)
import Control.Monad.ST.Strict (ST, runST)
import Control.Parallel.Strategies (parListChunk, parMap, rdeepseq, rpar, rparWith, runEval)
import Data.Aeson (FromJSON (..), ToJSON (..), object, withObject, (.:), (.=))
import Data.Complex
import Data.IntSet qualified as IntSet
import Data.List qualified
import Data.List.NonEmpty qualified
import Data.List.NonEmpty qualified as NonEmpty
import Data.Map.Strict (Map)
import Data.Map.Strict qualified as Map
import Data.Ratio
import Data.Set qualified as Set
import Data.Stream.Monadic (Step (Done, Skip, Yield), Stream (Stream))
import Data.Vector qualified as B
import Data.Vector.Algorithms.Search (binarySearch)
import Data.Vector.Generic ((!))
import Data.Vector.Generic qualified as G
import Data.Vector.Generic.Mutable qualified as GM
import Data.Vector.Storable qualified as S
import Data.Vector.Unboxed qualified as U
import Data.Vector.Unboxed.Mutable qualified as UM
import Foreign.C.Types (CDouble)
import Foreign.Marshal (free, new)
import Foreign.Ptr
import GHC.Exts (IsList (..))
import LatticeSymmetries.Algebra (sortVectorBy)
import LatticeSymmetries.Benes
import LatticeSymmetries.Dense
import LatticeSymmetries.FFI
import LatticeSymmetries.Utils
import Prelude hiding (identity, permutations, second, toList)

newtype PermutationGroup = PermutationGroup (B.Vector Permutation)
  deriving stock (Show, Eq)

pgLength :: PermutationGroup -> Int
pgLength (PermutationGroup gs) = G.length gs

getPeriodicity :: Permutation -> Int
getPeriodicity p₀ = go 1 p₀
  where
    identity = identityPermutation (G.length (unPermutation p₀))
    go !n !p
      | p == identity = n
      | otherwise = go (n + 1) (p₀ <> p)

nullPermutationGroup :: PermutationGroup -> Bool
nullPermutationGroup (PermutationGroup gs) = G.null gs

data Symmetry = Symmetry
  { symmetryPermutation :: !Permutation
  , symmetryPhase :: !(Ratio Int)
  }
  deriving stock (Show, Eq, Ord)

symmetrySector :: Symmetry -> Int
symmetrySector s
  | denominator sector == 1 = numerator sector
  | otherwise = error $ "this should not have happened: invalid symmetry: " <> show s
  where
    sector = symmetryPhase s * fromIntegral (getPeriodicity (symmetryPermutation s))

instance FromJSON Symmetry where
  parseJSON = withObject "Symmetry" $ \v -> do
    mkSymmetry <$> v .: "permutation" <*> v .: "sector"

instance ToJSON Symmetry where
  toJSON s = object ["permutation" .= symmetryPermutation s, "sector" .= symmetrySector s]

mkSymmetry :: Permutation -> Int -> Symmetry
mkSymmetry p sector
  | sector < 0 || sector >= periodicity =
      error $
        "invalid sector: " <> show sector <> "; permutation has periodicity " <> show periodicity
  | otherwise = Symmetry p (sector % periodicity)
  where
    periodicity = getPeriodicity p

symmetryNumberSites :: Symmetry -> Int
symmetryNumberSites (Symmetry p _) = G.length (unPermutation p)

modOne :: (Integral a) => Ratio a -> Ratio a
modOne x
  | x >= 1 = x - fromIntegral (numerator x `div` denominator x)
  | otherwise = x

instance Semigroup Symmetry where
  (<>) (Symmetry pa λa) (Symmetry pb λb) = Symmetry (pa <> pb) (modOne (λa + λb))

data Symmetries = Symmetries
  { symmGroup :: !PermutationGroup
  , symmNetwork :: !BatchedBenesNetwork
  , symmCharactersReal :: !(S.Vector Double)
  , symmCharactersImag :: !(S.Vector Double)
  }
  deriving stock (Show, Eq)

getCharacter :: Symmetry -> Complex Double
getCharacter (symmetryPhase -> φ)
  | φ == 0 = 1 :+ 0
  | φ == 1 % 2 = (-1) :+ 0
  | φ == (-1) % 2 = (-1) :+ 0
  | φ == 1 % 4 = 0 :+ 1
  | φ == (-1) % 4 = 0 :+ (-1)
  | otherwise = cos (-2 * pi * realToFrac φ) :+ sin (-2 * pi * realToFrac φ)

mkSymmetries :: [Symmetry] -> Either Text Symmetries
mkSymmetries [] = Right emptySymmetries
mkSymmetries gs@(g : _)
  | all ((== n) . symmetryNumberSites) gs =
      if isConsistent
        then
          let permutations = G.fromList $ symmetryPermutation <$> symmetries
              permGroup = PermutationGroup permutations
              benesNetwork = mkBatchedBenesNetwork $ G.map toBenesNetwork permutations
              characters = getCharacter <$> symmetries
              charactersReal = G.fromList $ realPart <$> characters
              charactersImag = G.fromList $ imagPart <$> characters
           in Right $ Symmetries permGroup benesNetwork charactersReal charactersImag
        else Left "incompatible symmetries"
  | otherwise = Left "symmetries have different number of sites"
  where
    n = symmetryNumberSites g
    identity = mkSymmetry (identityPermutation n) 0
    symmetries = fromGenerators identity gs
    isConsistent = and $ do
      s₁ <- symmetries
      s₂ <- symmetries
      let s₃@(Symmetry p₃ λ₃) = s₁ <> s₂
      pure $
        Set.member s₃ set
          && denominator (λ₃ * fromIntegral (getPeriodicity p₃)) == 1
      where
        set = Set.fromList symmetries

mkSymmetriesFromRepresentation :: (HasCallStack) => PermutationGroup -> Representation -> Symmetries
mkSymmetriesFromRepresentation (PermutationGroup g) (Representation r) =
  either error id $ mkSymmetries . G.toList $ G.zipWith Symmetry g r

instance FromJSON Symmetries where
  parseJSON xs = do
    r <- mkSymmetries <$> parseJSON xs
    either (fail . toString) pure r

instance ToJSON Symmetries where
  toJSON symmetries = toJSON (toSymmetryList symmetries)

instance IsList Symmetries where
  type Item Symmetries = Symmetry
  toList = toSymmetryList
  fromList x = either error id (mkSymmetries x)

symmetriesGetNumberBits :: Symmetries -> Int
symmetriesGetNumberBits (Symmetries (PermutationGroup gs) _ _ _)
  | G.null gs = 0
  | otherwise = G.length . unPermutation . G.head $ gs

emptySymmetries :: Symmetries
emptySymmetries = Symmetries (PermutationGroup G.empty) emptyBatchedBenesNetwork G.empty G.empty

areSymmetriesReal :: Symmetries -> Bool
areSymmetriesReal = G.all (== 0) . symmCharactersImag

toSymmetryList :: (HasCallStack) => Symmetries -> [Symmetry]
toSymmetryList (Symmetries (PermutationGroup gs) _ λsRe λsIm) =
  Data.List.zipWith3 toSymmetry (G.toList gs) (G.toList λsRe) (G.toList λsIm)
  where
    toSymmetry :: (HasCallStack) => Permutation -> Double -> Double -> Symmetry
    toSymmetry g λRe λIm
      | r ≈ 1 && s' ≈ fromIntegral (round s' :: Int) = mkSymmetry g (round s')
      | otherwise = error $ "failed to reconstruct Symmetry from " <> show (toList g) <> " and λ = " <> show λRe <> " + " <> show λIm <> "j"
      where
        -- φ = sector / periodicity, but care needs to be taken because φ is approximate
        s' = φ * fromIntegral n
        n = getPeriodicity g
        -- polar returns the phase in (-π, π], we add 2π to get a positive value
        (r, _φ) = polar (λRe :+ λIm)
        φ = _φ / (2 * pi) + (if _φ < 0 then 1 else 0)

fromGenerators :: (Semigroup a, Ord a) => a -> [a] -> [a]
fromGenerators _ [] = []
fromGenerators identity gs = go Set.empty (Set.singleton identity)
  where
    go !interior !boundary
      | Set.null boundary = Set.toAscList interior
      | otherwise = go interior' boundary'
      where
        interior' = interior `Set.union` boundary
        boundary' = Set.fromList [h <> g | h <- Set.toList boundary, g <- gs] Set.\\ interior'

nullSymmetries :: Symmetries -> Bool
nullSymmetries = nullPermutationGroup . symmGroup

newCpermutation_group :: Symmetries -> IO (Ptr Cpermutation_group)
newCpermutation_group x =
  let (DenseMatrix numberShifts numberMasks masks) = bbnMasks (symmNetwork x)
      shifts = bbnShifts $ symmNetwork x
      numberBits = symmetriesGetNumberBits x
   in S.unsafeWith masks $ \masksPtr ->
        S.unsafeWith shifts $ \shiftsPtr ->
          S.unsafeWith (symmCharactersReal x) $ \eigvalsRealPtr ->
            S.unsafeWith (symmCharactersImag x) $ \eigvalsImagPtr ->
              new $
                Cpermutation_group
                  { cpermutation_group_refcount = 0
                  , cpermutation_group_number_bits = fromIntegral numberBits
                  , cpermutation_group_number_shifts = fromIntegral numberShifts
                  , cpermutation_group_number_masks = fromIntegral numberMasks
                  , cpermutation_group_masks = masksPtr
                  , cpermutation_group_shifts = shiftsPtr
                  , cpermutation_group_eigvals_re = (castPtr :: Ptr Double -> Ptr CDouble) eigvalsRealPtr
                  , cpermutation_group_eigvals_im = (castPtr :: Ptr Double -> Ptr CDouble) eigvalsImagPtr
                  , cpermutation_group_haskell_payload = nullPtr
                  }

destroyCpermutation_group :: Ptr Cpermutation_group -> IO ()
destroyCpermutation_group = free

data Hypergraph a = Hypergraph {vertices :: !(Set a), hyperedges :: !(Set (Set a))}
  deriving stock (Eq, Ord, Show)

data SearchTree a r
  = SearchBranch !a [SearchTree a r]
  | SearchLeaf !(Maybe r)
  deriving stock (Eq, Ord, Show, Functor)

newtype Partitioning a = Partitioning (NonEmpty (U.Vector a))
  deriving stock (Show, Eq)

allSingletons :: (UM.Unbox a) => Partitioning a -> Bool
allSingletons (Partitioning xs) = all ((== 1) . G.length) xs

distancePartition
  :: forall a
   . (UM.Unbox a, Ord a)
  => Hypergraph a
  -- ^ The hypergraph to partition
  -> a
  -- ^ A selected vertex
  -> Partitioning a
  -- ^ Groups of vertices by their distance to the selected vertex
distancePartition g v0 =
  Partitioning . NonEmpty.fromList . fmap (G.fromList . Set.toAscList) $
    go (Set.singleton v0) (Set.toList g.hyperedges)
  where
    go :: Set a -> [Set a] -> [Set a]
    go !seen hyperedges
      | Set.null boundary = [foldl' Set.union Set.empty others | not (null others)]
      | otherwise = boundary : go (Set.union seen combined) others
      where
        shouldIncludeHyperedge = any (`Set.member` seen) . Set.toList
        (toBeIncluded, others) = Data.List.partition shouldIncludeHyperedge hyperedges
        !combined = foldl' Set.union Set.empty toBeIncluded
        !boundary = combined Set.\\ seen

-- | c n is the cyclic graph on n vertices
cyclicGraph :: Int -> Hypergraph Int
cyclicGraph n
  | n >= 3 = Hypergraph (Set.fromList vs) (Set.fromList (Set.fromList <$> es))
  | otherwise = error "n must be at least 3"
  where
    vs = [0 .. n - 1]
    es = [[i, (i + 1) `mod` n] | i <- [0 .. n - 1]]

-- cyclic hypergraph with 3-vertex edges
cyclicGraph3 :: Int -> Hypergraph Int
cyclicGraph3 n
  | n >= 3 = Hypergraph (Set.fromList vs) (Set.fromList (Set.fromList <$> es))
  | otherwise = error "n must be at least 3"
  where
    vs = [0 .. (n - 1)]
    es = [[i, (i + 1) `mod` n, (i + 2) `mod` n] | i <- [0 .. (n - 1)]]

rectangularGraph :: Int -> Int -> Hypergraph Int
rectangularGraph n k = Hypergraph (Set.fromList vs) (Set.fromList (Set.fromList <$> es))
  where
    -- vs contains an extra element that does not enter any edge in es
    vs = [0 .. n * k]
    es = [[k * i + j, k * i + ((j + 1) `mod` k)] | i <- [0 .. n - 1], j <- [0 .. k - 1]] ++ [[k * i + j, k * ((i + 1) `mod` n) + j] | i <- [0 .. n - 1], j <- [0 .. k - 1]]

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
  Partitioning . NonEmpty.fromList $
    intersectSorted <$> NonEmpty.toList p1 <*> NonEmpty.toList p2

createMappings :: (U.Unbox a) => Partitioning a -> Partitioning a -> [Mapping a]
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
                    mappings <> createMappings srcPart trgPart
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

naiveExtractLeaves :: SearchTree a Permutation -> [Permutation]
naiveExtractLeaves = go
  where
    go (SearchLeaf p) = maybe [] pure p
    go (SearchBranch _ branches) = concatMap go branches

transversalGeneratingSet :: SearchTree a Permutation -> B.Vector Permutation
transversalGeneratingSet = B.fromList . go
  where
    go (SearchLeaf (Just !p)) = [p | not (isIdentityPermutation p)]
    go (SearchBranch _ (t : ts)) = concatMap (force . take 1 . go) ts <> go t
    -- \| d < (1 :: Int) = concat (parMap (rparWith rdeepseq) (take 1 . go (d + 1)) ts) <> go (d + 1) t
    go _ = []

groupFromTransversalGeneratingSet :: B.Vector Permutation -> PermutationGroup
groupFromTransversalGeneratingSet tgs =
  PermutationGroup . sortVectorBy compare . B.fromList $
    Data.List.foldr1 (flip (<>)) <$> sequence transversals
  where
    k = permutationLength (G.head tgs)
    transversals =
      fmap ((identityPermutation k :) . B.toList . fmap fst) $
        B.groupBy (\g h -> snd g == snd h) $
          sortVectorBy (comparing snd) $
            fmap (\p -> (p, minimalSupport p)) tgs

hypergraphAutomorphisms :: Hypergraph Int -> PermutationGroup
hypergraphAutomorphisms =
  groupFromTransversalGeneratingSet . transversalGeneratingSet . autsSearchTree

abelianSubgroup :: PermutationGroup -> MultiplicationTable -> (PermutationGroup, MultiplicationTable)
abelianSubgroup (PermutationGroup g) t = (g', t')
  where
    -- abelianSubset returns candidates; we consider the first 100 and select the largest
    indices = Data.List.maximumBy (comparing G.length) . take 100 $ abelianSubset t
    t' = shrinkMultiplicationTable t indices
    g' = PermutationGroup . G.map (g G.!) . G.convert $ indices

newtype MultiplicationTable = MultiplicationTable {unMultiplicationTable :: DenseMatrix U.Vector Int}
  deriving stock (Show)

mkMultiplicationTable :: PermutationGroup -> MultiplicationTable
mkMultiplicationTable (PermutationGroup ps) =
  MultiplicationTable . generateDenseMatrix n n $ \i j ->
    binarySearch' ps (ps ! i <> ps ! j)
  where
    !n = G.length ps

data AbelianSubsetHistory = AbelianSubsetHistory
  { included :: !(U.Vector Bool)
  , mask :: !(U.Vector Bool)
  }
  deriving stock (Show)

commMatrixRow :: MultiplicationTable -> Int -> U.Vector Bool
commMatrixRow (MultiplicationTable matrix) i = U.generate matrix.dmCols $ \k ->
  indexDenseMatrix matrix (i, k) == indexDenseMatrix matrix (k, i)

addToAbelianSubsetHistory :: MultiplicationTable -> AbelianSubsetHistory -> Int -> AbelianSubsetHistory
addToAbelianSubsetHistory t h i = AbelianSubsetHistory included' mask'
  where
    included' = G.modify (\v -> UM.write v i True) h.included
    mask' = G.zipWith (&&) h.mask (commMatrixRow t i)

mergeAbelianSubsetHistories :: [AbelianSubsetHistory] -> [AbelianSubsetHistory]
mergeAbelianSubsetHistories = fmap go . Data.List.groupBy (\a b -> a.mask == b.mask)
  where
    go :: [AbelianSubsetHistory] -> AbelianSubsetHistory
    go hs = runST $ case hs of
      (h : rest) -> do
        mh <- G.thaw h.included
        forM_ rest $ \h' -> G.iforM_ h'.included $ \i v ->
          UM.modify mh (|| v) i
        AbelianSubsetHistory <$> G.unsafeFreeze mh <*> pure h.mask
      [] -> error "should never happen by construction"

abelianSubset :: MultiplicationTable -> [U.Vector Int]
abelianSubset t = (\h -> sortVectorBy compare . G.map fst . G.filter snd . G.indexed $ h.included) <$> go h0
  where
    !n = dmRows . unMultiplicationTable $ t
    h0 = AbelianSubsetHistory (G.replicate n False) (G.replicate n True)

    getChoices h =
      G.toList . G.filter (not . (h.included G.!)) . G.map fst . G.filter snd . G.indexed $ h.mask

    go :: AbelianSubsetHistory -> [AbelianSubsetHistory]
    go h = case getChoices h of
      [] -> [h]
      choices ->
        let hs' =
              Data.List.sortOn (negate . G.length . G.filter id . (.included))
                . mergeAbelianSubsetHistories
                . Data.List.sortOn (.mask)
                $ addToAbelianSubsetHistory t h <$> choices
         in concatMap go hs'

data GroupElement = GroupElement {index :: !Int, table :: !MultiplicationTable}
  deriving stock (Show)

instance Semigroup GroupElement where
  -- NOTE: a.table and b.table should be the same
  a <> b = GroupElement index' a.table
    where
      index' = unsafeIndexDenseMatrix (unMultiplicationTable a.table) (a.index, b.index)

binarySearch' :: (HasCallStack, G.Vector v a, Ord a) => v a -> a -> Int
binarySearch' !v !x = runST $ do
  mv <- G.unsafeThaw v
  index <- binarySearch mv x
  unless (v ! index == x) $ error "element not found"
  pure index

shrinkMultiplicationTable :: MultiplicationTable -> U.Vector Int -> MultiplicationTable
shrinkMultiplicationTable (MultiplicationTable matrix) indices =
  MultiplicationTable . generateDenseMatrix n n $ \i j ->
    binarySearch' indices $ indexDenseMatrix matrix (indices ! i, indices ! j)
  where
    n = G.length indices

getGroupElements :: MultiplicationTable -> [GroupElement]
getGroupElements t@(MultiplicationTable (DenseMatrix n _ _)) =
  [GroupElement i t | i <- [0 .. n - 1]]

primeFactors :: Int -> [Int]
primeFactors = go 2
  where
    go :: Int -> Int -> [Int]
    go !d !n
      | n <= 1 = []
      | d * d > n = [n]
      | n `mod` d == 0 = d : go d (n `div` d)
      | otherwise = go (d + 1) n

combineFactors :: [Int] -> [(Int, Int)]
combineFactors xs = (\l -> (head l, length l)) <$> Data.List.NonEmpty.group xs

groupElementOrder :: GroupElement -> Int
groupElementOrder x0 = go (x0 <> x0)
  where
    go !x
      | x.index == x0.index = 1
      | otherwise = 1 + go (x <> x0)

groupElementCycle :: GroupElement -> [GroupElement]
groupElementCycle x0 = x0 : go (x0 <> x0)
  where
    go !x
      | x.index == x0.index = []
      | otherwise = x : go (x <> x0)

isIdentity :: GroupElement -> Bool
isIdentity x = (x <> x).index == x.index

selectPrimeElements :: [GroupElement] -> [(GroupElement, Int, Int)]
selectPrimeElements = mapMaybe getInfo
  where
    getInfo x =
      let !order = groupElementOrder x
       in case combineFactors $ primeFactors order of
            [(p, n)] -> Just (x, p, n)
            _ -> Nothing

groupGenerators :: [GroupElement] -> [GroupElement]
groupGenerators group0 = go factors0 [] (selectPrimeElements group0) []
  where
    factors0 = combineFactors . primeFactors . length $ group0

    go [] gens _ _ = gens
    go ((p, n) : otherFactors) gens candidates subgroup =
      go factors' gens' candidates' subgroup'
      where
        (g, _, n') =
          Data.List.maximumBy (comparing third) $
            filter ((== p) . second) candidates
        factors'
          | n == n' = otherFactors
          | otherwise = (p, n - n') : otherFactors
        gens' = g : gens
        subgroup'
          | null subgroup = groupElementCycle g
          | otherwise = [a <> b | a <- subgroup, b <- groupElementCycle g]
        candidates' =
          let seen = IntSet.fromList $ (.index) <$> subgroup'
           in flip filter candidates $ \(cg, cp, cn) ->
                not (IntSet.member cg.index seen)
                  && cn <= Data.List.maximum [fn | (fp, fn) <- factors', fp == cp]

        second (_, b, _) = b
        third (_, _, c) = c

newtype Representation = Representation (B.Vector (Ratio Int))

groupRepresentationsFromGenerators :: B.Vector GroupElement -> [Representation]
groupRepresentationsFromGenerators generators =
  fmap (Representation . G.fromList) . mapM phases $ G.toList orders
  where
    orders = G.map groupElementOrder generators
    phases n = [i % n | i <- [0 .. n - 1]]

groupRepresentations :: MultiplicationTable -> [Representation]
groupRepresentations =
  groupRepresentationsFromGenerators . G.fromList . groupGenerators . getGroupElements

newtype Coset = Coset {unCoset :: IntSet.IntSet}
  deriving stock (Eq, Ord, Show)

getTimes :: MultiplicationTable -> Int -> Int -> Int
getTimes (MultiplicationTable matrix) = curry (unsafeIndexDenseMatrix matrix)

getInvert :: MultiplicationTable -> Int -> Int
getInvert t@(MultiplicationTable matrix) =
  \i ->
    case G.find (\j -> i `times` j == identityIndex) (U.generate n id) of
      Just j -> j
      Nothing -> error "group element has no inverse"
  where
    !n = matrix.dmRows
    !times = getTimes t
    !identityIndex =
      case G.find (\j -> isIdentity (GroupElement j t)) (U.generate n id) of
        Just k -> k
        Nothing -> error "the group has no or multiple identities"

commutatorSubgroup :: MultiplicationTable -> IntSet
commutatorSubgroup t =
  IntSet.fromList
    [ g `times` h `times` (inverses ! g) `times` (inverses ! h)
    | g <- [0 .. n - 1]
    , h <- [0 .. n - 1]
    ]
  where
    !n = t.unMultiplicationTable.dmRows
    !times = getTimes t
    !invert = getInvert t
    !inverses = U.generate n invert

abelianization :: MultiplicationTable -> (B.Vector Coset, MultiplicationTable)
abelianization t = (cosets, cosetMultiplicationTable t cosets)
  where
    n = t.unMultiplicationTable.dmRows
    cgs = commutatorSubgroup t
    times = getTimes t
    coset g = Coset $ IntSet.fromList [g `times` h | h <- IntSet.toList cgs]
    !cosets = go [] IntSet.empty [0 .. n - 1]
    go !cs !seen (g : gs)
      | g `IntSet.member` seen = go cs seen gs
      | otherwise =
          let !c' = coset g
              seen' = unCoset c' `IntSet.union` seen
           in go (c' : cs) seen' gs
    go cs _ [] = G.fromList cs

-- !cosets = B.fromList . Set.toAscList . Set.fromList $ coset <$> [0 .. n - 1]

cosetMultiplicationTable :: MultiplicationTable -> B.Vector Coset -> MultiplicationTable
cosetMultiplicationTable t cosets =
  MultiplicationTable . generateDenseMatrix n n $ \i j ->
    let x = Data.List.head . IntSet.elems . unCoset $ cosets ! i
        y = Data.List.head . IntSet.elems . unCoset $ cosets ! j
        !z = getTimes t x y
     in case G.findIndex (\(Coset coset) -> z `IntSet.member` coset) cosets of
          Just k -> k
          Nothing -> error "G/[G, G] is not closed under <>"
  where
    n = G.length cosets

-- strongGeneratingSet :: SearchTree [Mapping Int] -> [[Mapping Int]]
-- strongGeneratingSet = go []
--   where
--     go gs (SearchTree True xys _) = xys : gs
--     go gs (SearchTree False xys ts) =
--       case listToMaybe $ reverse $ filter (\(Mapping x y) -> x /= y) xys of -- the first vertex that isn't fixed
--         Nothing -> Data.List.foldl' (\hs t -> go hs t) gs ts
--         Just (Mapping x y) ->
--           if y `elem` (x .^^ gs)
--             then gs
--             else find1New gs ts
--     find1New gs (t : ts) =
--       let hs = instrongTerminals gs t
--        in if take 1 gs /= take 1 hs -- we know a new element would be placed at the front
--             then hs
--             else find1New gs ts
--     find1New gs [] = gs

isAutomorphism :: Hypergraph Int -> Permutation -> Bool
isAutomorphism g p = g.hyperedges == Set.map transformEdge g.hyperedges
  where
    i = unPermutation p
    transformEdge = Set.map (i G.!)
