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
  )
where

import Data.Aeson (FromJSON (..), ToJSON (..), object, withObject, (.:), (.=))
import Data.Complex
import Data.List qualified
import Data.Ratio
import Data.Set qualified as Set
import Data.Vector qualified as B
import Data.Vector.Generic qualified as G
import Data.Vector.Storable qualified as S
import Foreign.C.Types (CDouble)
import Foreign.Marshal (free, new)
import Foreign.Ptr
import GHC.Exts (IsList (..))
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

modOne :: Integral a => Ratio a -> Ratio a
modOne x
  | x >= 1 = x - fromIntegral (numerator x `div` denominator x)
  | otherwise = x

instance Semigroup Symmetry where
  (<>) (Symmetry pa λa) (Symmetry pb λb) =
    let !p = pa <> pb
        !λ = modOne (λa + λb)
     in Symmetry p λ

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
  | φ == 1 % 4 = 0 :+ (-1)
  | φ == 1 % 2 = (-1) :+ 0
  | φ == 3 % 4 = 0 :+ 1
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

toSymmetryList :: HasCallStack => Symmetries -> [Symmetry]
toSymmetryList (Symmetries (PermutationGroup gs) _ λsRe λsIm) =
  Data.List.zipWith3 toSymmetry (G.toList gs) (G.toList λsRe) (G.toList λsIm)
  where
    toSymmetry :: HasCallStack => Permutation -> Double -> Double -> Symmetry
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
