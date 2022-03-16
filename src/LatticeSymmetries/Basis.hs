{-# LANGUAGE CApiFFI #-}
{-# LANGUAGE DefaultSignatures #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE QuantifiedConstraints #-}
{-# LANGUAGE RankNTypes #-}

module LatticeSymmetries.Basis
  ( BasisState (..),
    Basis (..),
    Factor,
    ParticleTy (..),
    IndexType (..),
    GeneratorType (..),
    SpinfulOccupation (..),
    Cbasis (..),
    IsBasis (..),
    createCbasis,
    destroyCbasis,
    Cbasis_kernels (..),
    createCbasis_kernels,
    destroyCbasis_kernels,
    stateIndex,
    binomialCoefficient,
    flattenIndex,
    getNumberBits,
    getNumberWords,
    isStateIndexIdentity,
    withParticleType,
    withReconstructedBasis,
  )
where

import Control.Exception.Safe (bracket)
-- import Data.Bits
import Data.Some
import qualified Data.Text as Text
import qualified Data.Vector.Storable as S
import Foreign.C.Types
import Foreign.ForeignPtr
import Foreign.Marshal.Alloc (alloca, free)
import Foreign.Marshal.Utils (fromBool, new, with)
import Foreign.Ptr
import Foreign.StablePtr
import Foreign.Storable
import LatticeSymmetries.Algebra (Algebra (..))
import LatticeSymmetries.BitString
import LatticeSymmetries.Generator
import LatticeSymmetries.NonbranchingTerm
import System.IO.Unsafe (unsafePerformIO)
import Text.PrettyPrint.ANSI.Leijen (Pretty (..))
import qualified Text.PrettyPrint.ANSI.Leijen as Pretty

data BasisState = BasisState {-# UNPACK #-} !Int !BitString
  deriving stock (Show, Eq, Ord)

data ParticleTy = SpinTy | SpinfulFermionTy | SpinlessFermionTy
  deriving stock (Show, Eq)

class
  ( Enum (GeneratorType t),
    Bounded (GeneratorType t),
    HasMatrixRepresentation (GeneratorType t),
    HasNonbranchingRepresentation (Generator Int (GeneratorType t)),
    Algebra (GeneratorType t),
    Ord (IndexType t),
    HasSiteIndex (IndexType t),
    Pretty (IndexType t),
    Pretty (GeneratorType t),
    Pretty (Generator (IndexType t) (GeneratorType t)),
    HasSiteIndex (IndexType t)
  ) =>
  IsBasis t

instance IsBasis 'SpinTy

instance IsBasis 'SpinfulFermionTy

instance IsBasis 'SpinlessFermionTy

data Basis (t :: ParticleTy) where
  SpinBasis :: {sbNumberSites :: !Int, sbHammingWeight :: !(Maybe Int)} -> Basis 'SpinTy
  SpinfulFermionicBasis :: {sfbNumberSites :: !Int, sfbOccupation :: !SpinfulOccupation} -> Basis 'SpinfulFermionTy
  SpinlessFermionicBasis :: {fbNumberSites :: !Int, fbOccupation :: !(Maybe Int)} -> Basis 'SpinlessFermionTy

deriving stock instance Show (Basis t)

deriving stock instance Eq (Basis t)

-- typedef struct {
--   void* elts;
--   uint64_t num_elts;
--
--   void* freer;
-- } chpl_external_array;
data Cexternal_array = Cexternal_array
  { external_array_elts :: !(Ptr ()),
    external_array_num_elts :: !Word64,
    external_array_freer :: !(Ptr ())
  }

instance Storable Cexternal_array where
  sizeOf _ = 24
  alignment _ = 8
  peek p =
    Cexternal_array
      <$> peekByteOff p 0
      <*> peekByteOff p 8
      <*> peekByteOff p 16
  poke p x = do
    pokeByteOff p 0 (external_array_elts x)
    pokeByteOff p 8 (external_array_num_elts x)
    pokeByteOff p 16 (external_array_freer x)

newtype ChapelArray = ChapelArray (ForeignPtr Cexternal_array)

data Representatives = Representatives {rStates :: !ChapelArray, rRanges :: !ChapelArray}

buildRepresentatives :: Basis t -> Representatives
buildRepresentatives = undefined

type family IndexType (t :: ParticleTy) where
  IndexType 'SpinTy = Int
  IndexType 'SpinfulFermionTy = (SpinIndex, Int)
  IndexType 'SpinlessFermionTy = Int

type family GeneratorType (t :: ParticleTy) where
  GeneratorType 'SpinTy = SpinGeneratorType
  GeneratorType 'SpinfulFermionTy = FermionGeneratorType
  GeneratorType 'SpinlessFermionTy = FermionGeneratorType

type Factor t = Generator (IndexType t) (GeneratorType t)

data SpinfulOccupation
  = SpinfulNoOccupation
  | SpinfulTotalParticles !Int
  | -- | @SpinfulPerSector numberDown numberUp@
    SpinfulPerSector !Int !Int
  deriving stock (Show, Eq)

newtype {-# CTYPE "lattice_symmetries_haskell.h" "ls_hs_particle_type" #-} Cparticle_type = Cparticle_type CInt
  deriving stock (Show, Eq)
  deriving newtype (Storable)

c_LS_HS_SPIN :: Cparticle_type
c_LS_HS_SPIN = Cparticle_type 0

c_LS_HS_SPINFUL_FERMION :: Cparticle_type
c_LS_HS_SPINFUL_FERMION = Cparticle_type 1

c_LS_HS_SPINLESS_FERMION :: Cparticle_type
c_LS_HS_SPINLESS_FERMION = Cparticle_type 2

type Cindex_kernel = CPtrdiff -> Ptr Word64 -> CPtrdiff -> Ptr CPtrdiff -> CPtrdiff -> Ptr () -> IO ()

foreign import ccall "dynamic"
  mkCindex_kernel :: FunPtr Cindex_kernel -> Cindex_kernel

data {-# CTYPE "lattice_symmetries_haskell.h" "ls_hs_basis_kernels" #-} Cbasis_kernels = Cbasis_kernels
  { cbasis_state_index_kernel :: {-# UNPACK #-} !(FunPtr Cindex_kernel),
    cbasis_state_index_data :: {-# UNPACK #-} !(Ptr ())
  }

instance Storable Cbasis_kernels where
  sizeOf _ = 16
  alignment _ = 8
  peek p =
    Cbasis_kernels
      <$> peekByteOff p 0
      <*> peekByteOff p 8
  poke p (Cbasis_kernels index_kernel index_data) = do
    pokeByteOff p 0 index_kernel
    pokeByteOff p 8 index_data

data {-# CTYPE "lattice_symmetries_haskell.h" "ls_hs_basis" #-} Cbasis = Cbasis
  { cbasis_number_sites :: {-# UNPACK #-} !CInt,
    cbasis_number_particles :: {-# UNPACK #-} !CInt,
    cbasis_number_up :: {-# UNPACK #-} !CInt,
    cbasis_particle_type :: {-# UNPACK #-} !Cparticle_type,
    cbasis_state_index_is_identity :: {-# UNPACK #-} !CBool,
    cbasis_requires_projection :: {-# UNPACK #-} !CBool,
    cbasis_kernels :: {-# UNPACK #-} !(Ptr Cbasis_kernels),
    cbasis_haskell_payload :: {-# UNPACK #-} !(Ptr ())
  }

instance Storable Cbasis where
  sizeOf _ = 40
  alignment _ = 8
  peek p =
    Cbasis
      <$> peekByteOff p 0
      <*> peekByteOff p 4
      <*> peekByteOff p 8
      <*> peekByteOff p 12
      <*> peekByteOff p 16
      <*> peekByteOff p 17
      <*> peekByteOff p 24
      <*> peekByteOff p 32
  poke p x = do
    pokeByteOff p 0 (cbasis_number_sites x)
    pokeByteOff p 4 (cbasis_number_particles x)
    pokeByteOff p 8 (cbasis_number_up x)
    pokeByteOff p 12 (cbasis_particle_type x)
    pokeByteOff p 16 (cbasis_state_index_is_identity x)
    pokeByteOff p 17 (cbasis_requires_projection x)
    pokeByteOff p 24 (cbasis_kernels x)
    pokeByteOff p 32 (cbasis_haskell_payload x)

stateIndex :: Basis t -> BasisState -> Maybe Int
stateIndex basis (BasisState _ (BitString α))
  | α <= fromIntegral (maxBound :: Word64) = unsafePerformIO $
    bracket (createCbasis_kernels basis) destroyCbasis_kernels $ \kernelsPtr ->
      with (fromIntegral α :: Word64) $ \spinsPtr ->
        alloca $ \indicesPtr -> do
          kernels <- peek kernelsPtr
          mkCindex_kernel
            (cbasis_state_index_kernel kernels)
            1
            spinsPtr
            1
            indicesPtr
            1
            (cbasis_state_index_data kernels)
          i <- peek indicesPtr
          if i < 0
            then pure Nothing
            else pure $ Just (fromIntegral i)
  | otherwise = Nothing

data Cbinomial

--   = Cbinomial
--   { cbinomial_number_bits :: {-# UNPACK #-} !CInt,
--     cbinomial_coefficients :: {-# UNPACK #-} !(Ptr Word64)
--   }

foreign import capi unsafe "lattice_symmetries_haskell.h ls_hs_internal_malloc_binomials"
  ls_hs_internal_malloc_binomials :: CInt -> IO (Ptr Cbinomial)

foreign import capi unsafe "lattice_symmetries_haskell.h ls_hs_internal_free_binomials"
  ls_hs_internal_free_binomials :: Ptr Cbinomial -> IO ()

foreign import capi unsafe "lattice_symmetries_haskell.h &ls_hs_internal_free_binomials"
  ls_hs_internal_free_binomials_ptr :: FunPtr (Ptr Cbinomial -> IO ())

foreign import capi unsafe "lattice_symmetries_haskell.h ls_hs_internal_compute_binomials"
  ls_hs_internal_compute_binomials :: Ptr Cbinomial -> IO ()

foreign import capi unsafe "lattice_symmetries_haskell.h ls_hs_internal_binomial"
  ls_hs_internal_binomial :: CInt -> CInt -> Ptr Cbinomial -> IO CInt

createBinomials :: HasCallStack => Int -> IO (Ptr Cbinomial)
createBinomials numberBits = do
  p <- ls_hs_internal_malloc_binomials (fromIntegral numberBits)
  when (p == nullPtr) $ error "could not allocate Cbinomial"
  ls_hs_internal_compute_binomials p
  pure p

globalBinomials :: ForeignPtr Cbinomial
globalBinomials =
  unsafePerformIO $
    newForeignPtr ls_hs_internal_free_binomials_ptr =<< createBinomials 64
{-# NOINLINE globalBinomials #-}

binomialCoefficient :: Int -> Int -> Int
binomialCoefficient n k
  | n <= 64 && k <= 64 = unsafePerformIO $
    withForeignPtr globalBinomials $ \p ->
      fromIntegral <$> ls_hs_internal_binomial (fromIntegral n) (fromIntegral k) p
  | otherwise = error "too big"

-- bracket (createBinomials n) ls_hs_internal_free_binomials $ \p ->
--   fromIntegral <$> ls_hs_internal_binomial (fromIntegral n) (fromIntegral k) p

foreign import capi "lattice_symmetries_haskell.h &ls_hs_state_index_combinadics_kernel"
  ls_hs_state_index_combinadics_kernel :: FunPtr Cindex_kernel

foreign import capi "lattice_symmetries_haskell.h &ls_hs_state_index_identity_kernel"
  ls_hs_state_index_identity_kernel :: FunPtr Cindex_kernel

-- autoDestroyStateIndexKernel :: HasCallStack => FunPtr Cindex_kernel -> Ptr () -> IO ()
-- autoDestroyStateIndexKernel kernel env
--   | kernel == ls_hs_state_index_combinadics_kernel = ls_hs_internal_free_binomials (castPtr env)
--   | kernel == ls_hs_state_index_identity_kernel && env == nullPtr = pure ()
--   | kernel == nullFunPtr && env == nullPtr = pure ()
--   | otherwise = error "failed to automatically deallocate state_index kernel"

-- class IsBasis a where
--   type IndexType a
--   type GeneratorType a
--
--   numberStates :: a -> Maybe Int
--   minStateEstimate :: a -> BasisState
--   maxStateEstimate :: a -> BasisState
--
--   -- getNumberBits :: a -> Int
--
--   flattenIndex :: a -> IndexType a -> Int
--   default flattenIndex :: IndexType a ~ Int => a -> IndexType a -> Int
--   flattenIndex _ = id

-- getNumberBits :: IsBasis a => a -> Int
-- getNumberBits = fromIntegral . cbasis_number_sites . toCbasis

flattenIndex :: Basis t -> IndexType t -> Int
flattenIndex basis i = case basis of
  SpinBasis _ _ -> i
  SpinfulFermionicBasis n _ ->
    case i of
      (SpinUp, k) -> k
      (SpinDown, k) -> n + k
  SpinlessFermionicBasis _ _ -> i

getNumberSites :: Basis t -> Int
getNumberSites basis = case basis of
  SpinBasis n _ -> n
  SpinfulFermionicBasis n _ -> n
  SpinlessFermionicBasis n _ -> n

getNumberBits :: Basis t -> Int
getNumberBits basis = case basis of
  SpinBasis n _ -> n
  SpinfulFermionicBasis n _ -> 2 * n
  SpinlessFermionicBasis n _ -> n

getNumberWords :: Basis t -> Int
getNumberWords b = (n + 63) `div` 64
  where
    n = getNumberBits b

getNumberParticles :: Basis t -> Maybe Int
getNumberParticles basis = case basis of
  SpinBasis n _ -> Just n
  SpinfulFermionicBasis _ occupation ->
    case occupation of
      SpinfulNoOccupation -> Nothing
      SpinfulTotalParticles p -> Just p
      SpinfulPerSector u d -> Just (u + d)
  SpinlessFermionicBasis _ p -> p

getNumberUp :: Basis t -> Maybe Int
getNumberUp basis = case basis of
  SpinBasis _ h -> h
  SpinfulFermionicBasis _ occupation ->
    case occupation of
      SpinfulNoOccupation -> Nothing
      SpinfulTotalParticles _ -> Nothing
      SpinfulPerSector u _ -> Just u
  SpinlessFermionicBasis _ _ -> Nothing

getParticleType :: Basis t -> Cparticle_type
getParticleType basis = case basis of
  SpinBasis _ _ -> c_LS_HS_SPIN
  SpinfulFermionicBasis _ _ -> c_LS_HS_SPINFUL_FERMION
  SpinlessFermionicBasis _ _ -> c_LS_HS_SPINLESS_FERMION

isStateIndexIdentity :: Basis t -> Bool
isStateIndexIdentity basis = case basis of
  SpinBasis _ Nothing -> True
  SpinfulFermionicBasis _ SpinfulNoOccupation -> True
  SpinlessFermionicBasis _ Nothing -> True
  _ -> False

optionalNatural :: Maybe Int -> CInt
optionalNatural x = case x of
  (Just n) -> fromIntegral n
  Nothing -> (-1)

createCbasis :: Basis t -> IO (Ptr Cbasis)
createCbasis basis = do
  p <- castStablePtrToPtr <$> newStablePtr basis
  kernels <- createCbasis_kernels basis
  new $
    Cbasis
      { cbasis_number_sites = fromIntegral (getNumberSites basis),
        cbasis_number_particles = optionalNatural (getNumberParticles basis),
        cbasis_number_up = optionalNatural (getNumberUp basis),
        cbasis_particle_type = getParticleType basis,
        cbasis_state_index_is_identity = fromBool (isStateIndexIdentity basis),
        cbasis_requires_projection = fromBool False,
        cbasis_kernels = kernels,
        cbasis_haskell_payload = p
      }

withParticleType :: HasCallStack => Cparticle_type -> (forall (t :: ParticleTy). IsBasis t => Proxy t -> a) -> a
withParticleType t action
  | t == c_LS_HS_SPIN = action (Proxy :: Proxy 'SpinTy)
  | t == c_LS_HS_SPINFUL_FERMION = action (Proxy :: Proxy 'SpinfulFermionTy)
  | t == c_LS_HS_SPINLESS_FERMION = action (Proxy :: Proxy 'SpinlessFermionTy)
  | otherwise = error "invalid particle type"

withReconstructedBasis :: forall a. Ptr Cbasis -> (forall (t :: ParticleTy). IsBasis t => Basis t -> IO a) -> IO a
withReconstructedBasis p action = do
  basis <- peek p
  let run :: forall (t :: ParticleTy). IsBasis t => Proxy t -> IO a
      run _ = do
        (x :: Basis t) <- deRefStablePtr $ castPtrToStablePtr (cbasis_haskell_payload basis)
        action x
  withParticleType (cbasis_particle_type basis) run

destroyCbasis :: Ptr Cbasis -> IO ()
destroyCbasis p
  | p == nullPtr = pure ()
  | otherwise = do
    b <- peek p
    destroyCbasis_kernels (cbasis_kernels b)
    let freeHaskellPayload :: forall (t :: ParticleTy). Proxy t -> IO ()
        freeHaskellPayload _ = freeStablePtr (castPtrToStablePtr (cbasis_haskell_payload b) :: StablePtr (Basis t))
    withParticleType (cbasis_particle_type b) freeHaskellPayload
    free p

createCbasis_kernels :: HasCallStack => Basis t -> IO (Ptr Cbasis_kernels)
createCbasis_kernels basis -- (SpinBasis n hammingWeight)
  | getNumberBits basis > 64 =
    new
      Cbasis_kernels
        { cbasis_state_index_kernel = nullFunPtr,
          cbasis_state_index_data = nullPtr
        }
  | isStateIndexIdentity basis =
    new
      Cbasis_kernels
        { cbasis_state_index_kernel = ls_hs_state_index_identity_kernel,
          cbasis_state_index_data = nullPtr
        }
  | otherwise = case basis of
    (SpinfulFermionicBasis _ (SpinfulPerSector _ _)) -> error "not yet implemented"
    _ -> do
      cache <- createBinomials 64
      new
        Cbasis_kernels
          { cbasis_state_index_kernel = ls_hs_state_index_combinadics_kernel,
            cbasis_state_index_data = castPtr cache
          }

destroyCbasis_kernels :: HasCallStack => Ptr Cbasis_kernels -> IO ()
destroyCbasis_kernels p
  | p == nullPtr = pure ()
  | otherwise = do
    let destroy (Cbasis_kernels kernel env)
          | kernel == ls_hs_state_index_combinadics_kernel = ls_hs_internal_free_binomials (castPtr env)
          | kernel == ls_hs_state_index_identity_kernel && env == nullPtr = pure ()
          | kernel == nullFunPtr && env == nullPtr = pure ()
          | otherwise = error "failed to automatically deallocate state_index kernel"
    destroy =<< peek p
    free p

-- instance IsBasis SpinBasis where
--   type IndexType SpinBasis = Int
--   type GeneratorType SpinBasis = SpinGeneratorType
--   numberStates (SpinBasis n hammingWeight)
--     | n >= 64 = Nothing
--     | otherwise = Just $ case hammingWeight of
--       Just h -> binomialCoefficient n h
--       Nothing -> bit n
--   minStateEstimate (SpinBasis n hammingWeight) = BasisState n . BitString $ case hammingWeight of
--     Just h -> bit (h + 1) - 1
--     Nothing -> zeroBits
--   maxStateEstimate (SpinBasis n hammingWeight) = BasisState n . BitString $ case hammingWeight of
--     Just h -> (bit (h + 1) - 1) `shiftL` (n - h)
--     Nothing -> bit (n + 1) - 1
--   getNumberBits (SpinBasis n _) = n

-- instance IsBasis SpinfulFermionicBasis where
--   type IndexType SpinfulFermionicBasis = (SpinIndex, Int)
--   type GeneratorType SpinfulFermionicBasis = FermionGeneratorType
--   numberStates (SpinfulFermionicBasis n occupation)
--     | n >= 64 = Nothing
--     | otherwise = Just $ case occupation of
--       SpinfulNoOccupation -> bit (2 * n)
--       SpinfulTotalParticles p -> binomialCoefficient (2 * n) p
--       SpinfulPerSector u d -> binomialCoefficient n u * binomialCoefficient n d
--   minStateEstimate (SpinfulFermionicBasis n occupation) =
--     case occupation of
--       SpinfulNoOccupation -> minStateEstimate (SpinlessFermionicBasis (2 * n) Nothing)
--       SpinfulTotalParticles p -> minStateEstimate (SpinlessFermionicBasis (2 * n) (Just p))
--       SpinfulPerSector down up ->
--         let (BasisState _ minDown) = minStateEstimate (SpinlessFermionicBasis n (Just down))
--             (BasisState _ minUp) = minStateEstimate (SpinlessFermionicBasis n (Just up))
--          in BasisState (2 * n) $ (minUp `shiftL` n) .|. minDown
--   maxStateEstimate (SpinfulFermionicBasis n occupation) =
--     case occupation of
--       SpinfulNoOccupation -> maxStateEstimate (SpinlessFermionicBasis (2 * n) Nothing)
--       SpinfulTotalParticles p -> maxStateEstimate (SpinlessFermionicBasis (2 * n) (Just p))
--       SpinfulPerSector down up ->
--         let (BasisState _ minDown) = maxStateEstimate (SpinlessFermionicBasis n (Just down))
--             (BasisState _ minUp) = maxStateEstimate (SpinlessFermionicBasis n (Just up))
--          in BasisState (2 * n) $ (minUp `shiftL` n) .|. minDown
--   getNumberBits (SpinfulFermionicBasis n _) = 2 * n
--   flattenIndex (SpinfulFermionicBasis n _) i = case i of
--     (SpinUp, k) -> k
--     (SpinDown, k) -> k + n

-- instance IsBasis SpinlessFermionicBasis where
--   type IndexType SpinlessFermionicBasis = Int
--   type GeneratorType SpinlessFermionicBasis = FermionGeneratorType
--   numberStates (SpinlessFermionicBasis n occupation) = numberStates (SpinBasis n occupation)
--   minStateEstimate (SpinlessFermionicBasis n p) = minStateEstimate (SpinBasis n p)
--   maxStateEstimate (SpinlessFermionicBasis n p) = maxStateEstimate (SpinBasis n p)
--   getNumberBits (SpinlessFermionicBasis n _) = n
