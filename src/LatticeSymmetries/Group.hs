module LatticeSymmetries.Group
  ( Permutation,
    unPermutation,
    mkPermutation,
    identityPermutation,
    fromGenerators,
    pgLength,
    Symmetry,
    mkSymmetry,
    SymmetriesHeader (..),
    mkSymmetries,
    symmetriesFromHeader,
    Symmetries (..),
    areSymmetriesReal,
    nullSymmetries,
    emptySymmetries,
    emptySymmetriesHeader,
    symmetriesGetNumberBits,
    -- symmetriesFromSpec,
    -- borrowCsymmetries,
    -- withReconstructedSymmetries,

    -- ** Low-level interface for FFI
    Csymmetry (..),
    newCsymmetry,
    destroyCsymmetry,
    withCsymmetry,
    Csymmetries (..),
    newCsymmetries,
    destroyCsymmetries,
    withCsymmetries,
    borrowCpermutation_group,
    releaseCpermutation_group,
  )
where

import Data.Aeson (FromJSON (..), ToJSON (..), object, withObject, (.:), (.=))
import Data.Complex
import qualified Data.List
import Data.Ratio
import qualified Data.Set as Set
import qualified Data.Vector as B
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Storable as S
import Foreign.C.Types (CDouble)
import Foreign.ForeignPtr
import Foreign.Marshal (free, new)
import Foreign.Ptr
import Foreign.StablePtr
import Foreign.Storable
import GHC.Exts (IsList (..))
import GHC.ForeignPtr (unsafeForeignPtrToPtr)
import LatticeSymmetries.Benes
import LatticeSymmetries.Dense
import LatticeSymmetries.FFI
-- import LatticeSymmetries.IO
import LatticeSymmetries.Utils
import System.IO.Unsafe (unsafePerformIO)
import Prelude hiding (identity, permutations, toList)

newtype PermutationGroup = PermutationGroup (B.Vector Permutation)
  deriving stock (Show, Eq)

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
  { symmetryPermutation :: !Permutation,
    symmetryPhase :: !(Ratio Int)
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

modOne :: Integral a => Ratio a -> Ratio a
modOne x
  | x >= 1 = x - fromIntegral (numerator x `div` denominator x)
  | otherwise = x

instance Semigroup Symmetry where
  (<>) (Symmetry pa λa) (Symmetry pb λb) = Symmetry (pa <> pb) (modOne (λa + λb))

data SymmetriesHeader = SymmetriesHeader
  { symmHeaderGroup :: !PermutationGroup,
    symmHeaderNetwork :: !BatchedBenesNetwork,
    symmHeaderCharactersReal :: !(S.Vector Double),
    symmHeaderCharactersImag :: !(S.Vector Double)
  }
  deriving stock (Show, Eq)

symmetriesGetNumberBits :: SymmetriesHeader -> Int
symmetriesGetNumberBits (SymmetriesHeader (PermutationGroup gs) _ _ _)
  | G.null gs = 0
  | otherwise = G.length . unPermutation . G.head $ gs

emptySymmetriesHeader :: SymmetriesHeader
emptySymmetriesHeader = SymmetriesHeader (PermutationGroup G.empty) emptyBatchedBenesNetwork G.empty G.empty

data Symmetries = Symmetries
  { symmetriesHeader :: !SymmetriesHeader,
    symmetriesContents :: !(ForeignPtr Cpermutation_group)
  }
  deriving stock (Show)

instance Eq Symmetries where
  (==) a b = symmetriesHeader a == symmetriesHeader b

instance FromJSON SymmetriesHeader where
  parseJSON xs = do
    gs <- parseJSON xs
    case mkSymmetriesHeader gs of
      Just s -> pure s
      Nothing -> fail "incompatible symmetries"

instance ToJSON SymmetriesHeader where
  toJSON symmetries = toJSON (toSymmetryList symmetries)

instance FromJSON Symmetries where
  parseJSON xs = mkSymmetries <$> parseJSON xs

instance ToJSON Symmetries where
  toJSON symmetries = toJSON (toSymmetryList (symmetriesHeader symmetries))

areSymmetriesReal :: SymmetriesHeader -> Bool
areSymmetriesReal = G.all (== 0) . symmHeaderCharactersImag

toSymmetryList :: SymmetriesHeader -> [Symmetry]
toSymmetryList (SymmetriesHeader (PermutationGroup gs) _ λsRe λsIm) =
  Data.List.zipWith3 _toSymmetry (G.toList gs) (G.toList λsRe) (G.toList λsIm)
  where
    _toSymmetry :: Permutation -> Double -> Double -> Symmetry
    _toSymmetry g λRe λIm
      | r ≈ 1 && s' ≈ fromIntegral (round s' :: Int) = mkSymmetry g (round s')
      | otherwise = error $ "failed to reconstruct Symmetry from " <> show (toList g) <> " and λ = " <> show λRe <> " + " <> show λIm <> "j"
      where
        -- φ = sector / periodicity, but care needs to be taken because φ is approximate
        s' = φ * fromIntegral n
        n = getPeriodicity g
        -- polar returns the phase in (-π, π], we add 2π to get a positive value
        (r, _φ) = polar (λRe :+ λIm)
        φ = _φ / (2 * pi) + (if _φ < 0 then 1 else 0)

pgLength :: PermutationGroup -> Int
pgLength (PermutationGroup gs) = G.length gs

fromGenerators :: (Semigroup a, Ord a) => a -> [a] -> [a]
fromGenerators _ [] = []
fromGenerators identity gs = go Set.empty (Set.singleton identity)
  where
    go !interior !boundary
      | Set.null boundary = Set.toAscList $ interior
      | otherwise = go interior' boundary'
      where
        interior' = interior `Set.union` boundary
        boundary' = Set.fromList [h <> g | h <- Set.toList boundary, g <- gs] Set.\\ interior'

nullSymmetries :: SymmetriesHeader -> Bool
nullSymmetries = nullPermutationGroup . symmHeaderGroup

emptySymmetries :: Symmetries
emptySymmetries = symmetriesFromHeader emptySymmetriesHeader

-- symmetriesGetNumberBits :: Symmetries -> Int
-- symmetriesGetNumberBits = getNumberBits . symmetriesHeader

mkSymmetries :: [Symmetry] -> Symmetries
mkSymmetries gs = case symmetriesFromHeader <$> mkSymmetriesHeader gs of
  Just s -> s
  Nothing -> error "incompatible symmetries"

mkSymmetriesHeader :: [Symmetry] -> Maybe SymmetriesHeader
mkSymmetriesHeader [] = Just emptySymmetriesHeader
mkSymmetriesHeader gs@(g : _)
  | all ((== n) . symmetryNumberSites) gs = case isConsistent of
      True ->
        let permutations = G.fromList $ symmetryPermutation <$> symmetries
            permGroup = PermutationGroup permutations
            benesNetwork = mkBatchedBenesNetwork $ G.map toBenesNetwork permutations
            charactersReal = G.fromList $ (\φ -> cos (-2 * pi * realToFrac φ)) <$> symmetryPhase <$> symmetries
            charactersImag = G.fromList $ (\φ -> sin (-2 * pi * realToFrac φ)) <$> symmetryPhase <$> symmetries
         in Just $ SymmetriesHeader permGroup benesNetwork charactersReal charactersImag
      False -> Nothing
  | otherwise = error "symmetries have different number of sites"
  where
    n = symmetryNumberSites g
    identity = mkSymmetry (identityPermutation n) 0
    symmetries = fromGenerators identity gs
    isConsistent = all id $ do
      s₁ <- symmetries
      s₂ <- symmetries
      let s₃@(Symmetry p₃ λ₃) = s₁ <> s₂
      pure $
        Set.member s₃ set
          && denominator (λ₃ * fromIntegral (getPeriodicity p₃)) == 1
      where
        set = Set.fromList symmetries

symmetriesFromHeader :: SymmetriesHeader -> Symmetries
symmetriesFromHeader x = unsafePerformIO $ do
  fp <- mallocForeignPtr
  let (DenseMatrix numberShifts numberMasks masks) = bbnMasks $ symmHeaderNetwork x
      shifts = bbnShifts $ symmHeaderNetwork x
      numberBits = symmetriesGetNumberBits x
  withForeignPtr fp $ \ptr ->
    S.unsafeWith masks $ \masksPtr ->
      S.unsafeWith shifts $ \shiftsPtr ->
        S.unsafeWith (symmHeaderCharactersReal x) $ \eigvalsRealPtr ->
          S.unsafeWith (symmHeaderCharactersImag x) $ \eigvalsImagPtr ->
            poke ptr $
              Cpermutation_group
                { cpermutation_group_refcount = 0,
                  cpermutation_group_number_bits = fromIntegral numberBits,
                  cpermutation_group_number_shifts = fromIntegral numberShifts,
                  cpermutation_group_number_masks = fromIntegral numberMasks,
                  cpermutation_group_masks = masksPtr,
                  cpermutation_group_shifts = shiftsPtr,
                  cpermutation_group_eigvals_re = (castPtr :: Ptr Double -> Ptr CDouble) eigvalsRealPtr,
                  cpermutation_group_eigvals_im = (castPtr :: Ptr Double -> Ptr CDouble) eigvalsImagPtr,
                  cpermutation_group_haskell_payload = nullPtr
                }
  pure $ Symmetries x fp
{-# NOINLINE symmetriesFromHeader #-}

-- symmetryFromSpec :: SymmetrySpec -> Symmetry
-- symmetryFromSpec (SymmetrySpec p k) =
--   mkSymmetry (mkPermutation (G.fromList (toList p))) k

-- symmetriesFromSpec :: [SymmetrySpec] -> Symmetries
-- symmetriesFromSpec = mkSymmetries . fmap symmetryFromSpec

-- borrowCsymmetries :: Symmetries -> IO (Ptr Cpermutation_group)
-- borrowCsymmetries symm = do
--   withForeignPtr (symmetriesContents symm) $ \ptr -> do
--     _ <- symmetriesIncRefCount ptr
--     symmetriesPokePayload ptr =<< newStablePtr symm
--   pure $ unsafeForeignPtrToPtr (symmetriesContents symm)

-- withReconstructedSymmetries :: Ptr Cpermutation_group -> (Symmetries -> IO a) -> IO a
-- withReconstructedSymmetries p action =
--   action =<< deRefStablePtr =<< symmetriesPeekPayload p

newtype {-# CTYPE "lattice_symmetries_haskell.h" "ls_hs_symmetry" #-} Csymmetry = Csymmetry
  { unCsymmetry :: StablePtr Symmetry
  }
  deriving stock (Eq)
  deriving newtype (Storable)

newCsymmetry :: Symmetry -> IO (Ptr Csymmetry)
newCsymmetry symm = do
  new =<< Csymmetry <$> newStablePtr symm

destroyCsymmetry :: Ptr Csymmetry -> IO ()
destroyCsymmetry p = do
  freeStablePtr =<< (unCsymmetry <$> peek p)
  free p

withCsymmetry :: Ptr Csymmetry -> (Symmetry -> IO a) -> IO a
withCsymmetry p f = f =<< deRefStablePtr . unCsymmetry =<< peek p

newtype {-# CTYPE "lattice_symmetries_haskell.h" "ls_hs_symmetries" #-} Csymmetries = Csymmetries
  { unCsymmetries :: StablePtr SymmetriesHeader
  }
  deriving stock (Eq)
  deriving newtype (Storable)

newCsymmetries :: SymmetriesHeader -> IO (Ptr Csymmetries)
newCsymmetries symm = new =<< Csymmetries <$> newStablePtr symm

destroyCsymmetries :: Ptr Csymmetries -> IO ()
destroyCsymmetries p = do
  freeStablePtr =<< (unCsymmetries <$> peek p)
  free p

withCsymmetries :: Ptr Csymmetries -> (SymmetriesHeader -> IO a) -> IO a
withCsymmetries p f = f =<< deRefStablePtr . unCsymmetries =<< peek p

borrowCpermutation_group :: SymmetriesHeader -> IO (Ptr Cpermutation_group)
borrowCpermutation_group x =
  let (DenseMatrix numberShifts numberMasks masks) = bbnMasks $ symmHeaderNetwork x
      shifts = bbnShifts $ symmHeaderNetwork x
      numberBits = symmetriesGetNumberBits x
   in (new =<<) $
        S.unsafeWith masks $ \masksPtr ->
          S.unsafeWith shifts $ \shiftsPtr ->
            S.unsafeWith (symmHeaderCharactersReal x) $ \eigvalsRealPtr ->
              S.unsafeWith (symmHeaderCharactersImag x) $ \eigvalsImagPtr ->
                pure $
                  Cpermutation_group
                    { cpermutation_group_refcount = 0,
                      cpermutation_group_number_bits = fromIntegral numberBits,
                      cpermutation_group_number_shifts = fromIntegral numberShifts,
                      cpermutation_group_number_masks = fromIntegral numberMasks,
                      cpermutation_group_masks = masksPtr,
                      cpermutation_group_shifts = shiftsPtr,
                      cpermutation_group_eigvals_re = (castPtr :: Ptr Double -> Ptr CDouble) eigvalsRealPtr,
                      cpermutation_group_eigvals_im = (castPtr :: Ptr Double -> Ptr CDouble) eigvalsImagPtr,
                      cpermutation_group_haskell_payload = nullPtr
                    }

releaseCpermutation_group :: Ptr Cpermutation_group -> IO ()
releaseCpermutation_group = free
