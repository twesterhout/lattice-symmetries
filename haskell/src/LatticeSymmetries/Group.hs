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
  , areSymmetriesReal
  , nullSymmetries
  , emptySymmetries
  , symmetriesGetNumberBits

    -- ** Low-level interface for FFI
  , newCpermutation_group
  , destroyCpermutation_group

    -- * Automorphisms
  , Hypergraph (..)
  , cyclicGraph
  , cyclicGraph3
  , rectangularGraph
  , distancePartition
  , autsSearchTree
  , isAutomorphism
  , transversalGeneratingSet
  , groupFromTransversalGeneratingSet
  )
where

import Control.Exception (assert)
import Data.Aeson (FromJSON (..), ToJSON (..), object, withObject, (.:), (.=))
import Data.Complex
import Data.List qualified
import Data.List.NonEmpty qualified as NonEmpty
import Data.Map.Strict (Map)
import Data.Map.Strict qualified as Map
import Data.Ratio
import Data.Set qualified as Set
import Data.Vector qualified as B
import Data.Vector.Generic qualified as G
import Data.Vector.Storable qualified as S
import Foreign.C.Types (CDouble)
import Foreign.Marshal (free, new)
import Foreign.Ptr
import GHC.Exts (IsList (..))
import LatticeSymmetries.Algebra (sortVectorBy)
import LatticeSymmetries.Benes
import LatticeSymmetries.Dense
import LatticeSymmetries.FFI
import LatticeSymmetries.Utils
import Prelude hiding (identity, permutations, toList)

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

data SearchTree a = SearchTree !Bool !a [SearchTree a]
  deriving stock (Eq, Ord, Show, Functor)

newtype Partitioning a = Partitioning (NonEmpty [a])
  deriving stock (Show, Eq)

allSingletons :: Partitioning a -> Bool
allSingletons (Partitioning xs) = all p xs
  where
    p [_] = True
    p _ = False

distancePartition
  :: forall a
   . (Ord a)
  => Hypergraph a
  -- ^ The hypergraph to partition
  -> a
  -- ^ A selected vertex
  -> Partitioning a
  -- ^ Groups of vertices by their distance to the selected vertex
distancePartition g v0 =
  Partitioning . NonEmpty.fromList . fmap Set.toList $
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
    vs = [0 .. n * k - 1]
    es = [[k * i + j, k * i + ((j + 1) `mod` k)] | i <- [0 .. n - 1], j <- [0 .. k - 1]] ++ [[k * i + j, k * ((i + 1) `mod` n) + j] | i <- [0 .. n - 1], j <- [0 .. k - 1]]

intersectAllToAll :: (Ord a) => Partitioning a -> Partitioning a -> Partitioning a
intersectAllToAll (Partitioning p1) (Partitioning p2) =
  Partitioning . NonEmpty.fromList $
    Data.List.intersect <$> NonEmpty.toList p1 <*> NonEmpty.toList p2

createMappings :: Partitioning a -> Partitioning a -> [Mapping a]
createMappings (Partitioning (NonEmpty.toList -> src)) (Partitioning (NonEmpty.toList -> tgt)) =
  zipWith Mapping (concat src) (concat tgt)

pickOne :: (Ord a) => Map a (Partitioning a) -> Partitioning a -> (a, Partitioning a)
pickOne _ (Partitioning ([] :| _)) = error "pickOne expects a normalized Partitioning"
pickOne dps (Partitioning ((x : xs) :| rest)) = (x, intersectAllToAll p (dps Map.! x))
  where
    p = Partitioning (xs :| rest)

pickAll :: (Ord a) => Map a (Partitioning a) -> Partitioning a -> [(a, Partitioning a)]
pickAll dps (Partitioning (xss :| rest)) = do
  (x, xs) <- picks xss
  let p = Partitioning (xs :| rest)
  pure (x, intersectAllToAll p (dps Map.! x))
  where
    picks :: [a] -> [(a, [a])]
    picks [] = []
    picks (x : xs) = (x, xs) : [(y, x : ys) | (y, ys) <- picks xs]

normalizeWhenCompatible :: (Ord a) => Partitioning a -> Partitioning a -> Maybe (Partitioning a, Partitioning a)
normalizeWhenCompatible (Partitioning (NonEmpty.toList -> src)) (Partitioning (NonEmpty.toList -> tgt))
  | fmap length src == fmap length tgt =
      let normalize = Partitioning . NonEmpty.fromList . filter (not . null)
       in -- (src', tgt') = unzip $ Data.List.sort $ zip (filter (not . null) src) (filter (not . null) tgt)
          Just (normalize src, normalize tgt)
  | otherwise = Nothing

-- | Generate a SearchTree of graph automorphisms using distance partition
autsSearchTree :: Hypergraph Int -> SearchTree [Mapping Int]
autsSearchTree g = dfs [] (Partitioning (vs :| [])) (Partitioning (vs :| []))
  where
    dfs :: [Mapping Int] -> Partitioning Int -> Partitioning Int -> SearchTree [Mapping Int]
    dfs mappings srcPart trgPart
      -- check compatibility at the final step
      | allSingletons srcPart =
          -- isCompatible xys = isAutomorphism g xys
          assert (allSingletons trgPart) $
            let p = mappings <> createMappings srcPart trgPart
             in SearchTree (isAutomorphism g p) p []
      | otherwise =
          let (x, srcPart') = pickOne dps srcPart
              branches = do
                (y, trgPart') <- pickAll dps trgPart
                case normalizeWhenCompatible srcPart' trgPart' of
                  Just (srcPart'', trgPart'') ->
                    pure $
                      dfs (Mapping x y : mappings) srcPart'' trgPart''
                  Nothing -> []
           in SearchTree False mappings branches
    vs = Set.toAscList g.vertices
    !dps = Map.fromAscList [(v, distancePartition g v) | v <- vs]

transversalGeneratingSet :: Int -> SearchTree [Mapping Int] -> B.Vector Permutation
transversalGeneratingSet k = B.fromList . go
  where
    go (SearchTree True xys _) =
      let !p = permutationFromMappings (Just k) xys
       in [p | not (isIdentityPermutation p)]
    go (SearchTree False _ []) = []
    go (SearchTree False _ (t : ts)) = concatMap (take 1 . go) ts <> go t

groupFromTransversalGeneratingSet :: Int -> B.Vector Permutation -> B.Vector Permutation
groupFromTransversalGeneratingSet k tgs =
  B.fromList $ Data.List.foldr1 (<>) <$> sequence transversals
  where
    transversals =
      fmap ((identityPermutation k :) . B.toList . fmap fst) $
        B.groupBy (\g h -> snd g == snd h) $
          sortVectorBy (comparing snd) $
            fmap (\p -> (p, minimalSupport p)) tgs

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

isAutomorphism :: Hypergraph Int -> [Mapping Int] -> Bool
isAutomorphism g mappings = g.hyperedges == Set.map transformEdge g.hyperedges
  where
    i =
      unPermutation $
        permutationFromMappings (Just (Set.size g.vertices)) mappings
    transformEdge = Set.map (i G.!)
