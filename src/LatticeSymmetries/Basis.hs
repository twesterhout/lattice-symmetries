{-# LANGUAGE CApiFFI #-}
{-# LANGUAGE DefaultSignatures #-}
{-# LANGUAGE QuantifiedConstraints #-}

module LatticeSymmetries.Basis
  ( BasisState (..),
    SpinBasis (..),
    SpinfulFermionicBasis (..),
    SpinlessFermionicBasis (..),
    SpinfulOccupation (..),
    IsBasis (..),
    Cbasis (..),
    stateIndex,
    binomialCoefficient,
    getNumberBits,
    getNumberWords,
    BitString (..),
    wordsFromBitString,
    writeBitString,
  )
where

import Control.Exception.Safe (bracket)
import Data.Bits
import Foreign.C.Types
import Foreign.Marshal.Alloc (alloca)
import Foreign.Marshal.Utils (fromBool, with)
import Foreign.Ptr
import Foreign.Storable
import LatticeSymmetries.Algebra (FermionGeneratorType, SpinGeneratorType, SpinIndex (..))
import LatticeSymmetries.BitString
import System.IO.Unsafe (unsafePerformIO)

-- writeWords :: (Bits i, Integral i) => Int -> Ptr Word64 -> i -> IO ()
-- writeWords count p x₀ = go 0 x₀
--   where
--     nextWord !x = case bitSizeMaybe x of
--       Just n -> if n > 64 then x `shiftR` 64 else 0
--       Nothing -> x `shiftR` 64
--     go !i !x
--       | i < count = do
--         pokeElemOff p i (fromIntegral x :: Word64)
--         go (i + 1) (nextWord x)
--       | otherwise = pure ()

data BasisState = BasisState {-# UNPACK #-} !Int {-# UNPACK #-} !BitString
  deriving stock (Show, Eq, Ord)

data SpinBasis = SpinBasis
  { sbNumberSites :: !Int,
    sbMagnetization :: !(Maybe Int)
  }
  deriving stock (Show, Eq)

data SpinfulOccupation
  = SpinfulNoOccupation
  | SpinfulTotalParticles !Int
  | -- | @SpinfulPerSector numberDown numberUp@
    SpinfulPerSector !Int !Int
  deriving stock (Show, Eq)

data SpinfulFermionicBasis = SpinfulFermionicBasis
  { sfbNumberSites :: !Int,
    sfbOccupation :: !SpinfulOccupation
  }
  deriving stock (Show, Eq)

data SpinlessFermionicBasis = SpinlessFermionicBasis
  { fbNumberSites :: !Int,
    fbOccupation :: !(Maybe Int)
  }
  deriving stock (Show, Eq)

newtype Cparticle_type = Cparticle_type CInt
  deriving stock (Show, Eq)
  deriving newtype (Storable)

c_LS_HS_SPIN :: Cparticle_type
c_LS_HS_SPIN = Cparticle_type 0

c_LS_HS_FERMION :: Cparticle_type
c_LS_HS_FERMION = Cparticle_type 1

-- |
-- Cstate_index_kernel
--   batch_size
--   spins
--   spins_stride
--   indices
--   indices_stride
--   private_kernel_data
type Cindex_kernel = CPtrdiff -> Ptr Word64 -> CPtrdiff -> Ptr CPtrdiff -> CPtrdiff -> Ptr () -> IO ()

foreign import ccall "dynamic"
  mkCindex_kernel :: FunPtr Cindex_kernel -> Cindex_kernel

data Cbasis = Cbasis
  { cbasis_number_sites :: {-# UNPACK #-} !CInt,
    cbasis_number_particles :: {-# UNPACK #-} !CInt,
    cbasis_number_up :: {-# UNPACK #-} !CInt,
    cbasis_particle_type :: {-# UNPACK #-} !Cparticle_type,
    cbasis_state_index_is_identity :: {-# UNPACK #-} !CBool
  }

instance Storable Cbasis where
  sizeOf _ = 20
  alignment _ = 4
  peek p =
    Cbasis
      <$> peekByteOff p 0
      <*> peekByteOff p 4
      <*> peekByteOff p 8
      <*> peekByteOff p 12
      <*> peekByteOff p 16
  poke p (Cbasis number_sites number_particles number_up particle_type state_index_is_identity) = do
    pokeByteOff p 0 number_sites
    pokeByteOff p 4 number_particles
    pokeByteOff p 8 number_up
    pokeByteOff p 12 particle_type
    pokeByteOff p 16 state_index_is_identity

data Cbasis_kernels = Cbasis_kernels
  { cbasis_state_index_kernel :: {-# UNPACK #-} !(FunPtr Cindex_kernel),
    cbasis_state_index_data :: {-# UNPACK #-} !(Ptr ())
  }

stateIndex :: IsBasis basis => basis -> BasisState -> Maybe Int
stateIndex basis (BasisState _ (BitString α))
  | α <= fromIntegral (maxBound :: Word64) = unsafePerformIO $
    bracket (createCbasis_kernels basis) (destroyCbasis_kernels basis) $ \kernels ->
      with (fromIntegral α :: Word64) $ \spinsPtr ->
        alloca $ \indicesPtr -> do
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

binomialCoefficient :: Int -> Int -> Int
binomialCoefficient n k = unsafePerformIO $
  bracket (createBinomials n) ls_hs_internal_free_binomials $ \p ->
    fromIntegral <$> ls_hs_internal_binomial (fromIntegral n) (fromIntegral k) p

foreign import capi "lattice_symmetries_haskell.h &ls_hs_state_index_combinadics_kernel"
  ls_hs_state_index_combinadics_kernel :: FunPtr Cindex_kernel

foreign import capi "lattice_symmetries_haskell.h &ls_hs_state_index_identity_kernel"
  ls_hs_state_index_identity_kernel :: FunPtr Cindex_kernel

autoDestroyStateIndexKernel :: HasCallStack => FunPtr Cindex_kernel -> Ptr () -> IO ()
autoDestroyStateIndexKernel kernel env
  | kernel == ls_hs_state_index_combinadics_kernel = ls_hs_internal_free_binomials (castPtr env)
  | kernel == ls_hs_state_index_identity_kernel && env == nullPtr = pure ()
  | kernel == nullFunPtr && env == nullPtr = pure ()
  | otherwise = error "failed to automatically deallocate state_index kernel"

class IsBasis a where
  type IndexType a
  type GeneratorType a
  toCbasis :: a -> Cbasis
  createCbasis_kernels :: HasCallStack => a -> IO Cbasis_kernels
  destroyCbasis_kernels :: a -> Cbasis_kernels -> IO ()

  minStateEstimate :: a -> BasisState
  maxStateEstimate :: a -> BasisState

  flattenIndex :: a -> IndexType a -> Int
  default flattenIndex :: IndexType a ~ Int => a -> IndexType a -> Int
  flattenIndex _ = id

getNumberBits :: IsBasis a => a -> Int
getNumberBits = fromIntegral . cbasis_number_sites . toCbasis

getNumberWords :: IsBasis a => a -> Int
getNumberWords b = (n + 63) `div` 64
  where
    n = getNumberBits b

instance IsBasis SpinBasis where
  type IndexType SpinBasis = Int
  type GeneratorType SpinBasis = SpinGeneratorType
  toCbasis (SpinBasis n (Just m)) =
    Cbasis
      { cbasis_number_sites = fromIntegral n,
        cbasis_number_particles = fromIntegral n,
        -- we have: cbasis_number_up + cbasis_number_down = n
        --          cbasis_number_up - cbasis_number_down = m
        cbasis_number_up = fromIntegral $ (n + m) `div` 2,
        cbasis_particle_type = c_LS_HS_SPIN,
        cbasis_state_index_is_identity = fromBool False
      }
  toCbasis (SpinBasis n Nothing) =
    Cbasis
      { cbasis_number_sites = fromIntegral n,
        cbasis_number_particles = fromIntegral n,
        cbasis_number_up = -1,
        cbasis_particle_type = c_LS_HS_SPIN,
        cbasis_state_index_is_identity = fromBool True
      }
  createCbasis_kernels (SpinBasis n magnetization)
    | n <= 64 =
      case magnetization of
        Just _ -> do
          cache <- createBinomials 64
          pure $
            Cbasis_kernels
              { cbasis_state_index_kernel = ls_hs_state_index_combinadics_kernel,
                cbasis_state_index_data = castPtr cache
              }
        Nothing ->
          pure $
            Cbasis_kernels
              { cbasis_state_index_kernel = ls_hs_state_index_identity_kernel,
                cbasis_state_index_data = nullPtr
              }
    | otherwise =
      pure $
        Cbasis_kernels
          { cbasis_state_index_kernel = nullFunPtr,
            cbasis_state_index_data = nullPtr
          }
  destroyCbasis_kernels _ (Cbasis_kernels stateIndexKernel stateIndexEnv) = do
    autoDestroyStateIndexKernel stateIndexKernel stateIndexEnv
  minStateEstimate (SpinBasis n Nothing) = BasisState n $ zeroBits
  minStateEstimate (SpinBasis n (Just m)) = BasisState n . BitString $ bit (m + 1) - 1
  maxStateEstimate (SpinBasis n Nothing) = BasisState n . BitString $ bit (n + 1) - 1
  maxStateEstimate (SpinBasis n (Just m)) = BasisState n . BitString $ (bit (m + 1) - 1) `shiftL` (n - m)

instance IsBasis SpinfulFermionicBasis where
  type IndexType SpinfulFermionicBasis = (SpinIndex, Int)
  type GeneratorType SpinfulFermionicBasis = FermionGeneratorType
  toCbasis (SpinfulFermionicBasis n SpinfulNoOccupation) =
    Cbasis
      { cbasis_number_sites = fromIntegral n,
        cbasis_number_particles = -1,
        cbasis_number_up = -1,
        cbasis_particle_type = c_LS_HS_FERMION,
        cbasis_state_index_is_identity = fromBool True
      }
  toCbasis (SpinfulFermionicBasis n (SpinfulTotalParticles p)) =
    Cbasis
      { cbasis_number_sites = fromIntegral n,
        cbasis_number_particles = fromIntegral p,
        cbasis_number_up = -1,
        cbasis_particle_type = c_LS_HS_FERMION,
        cbasis_state_index_is_identity = fromBool False
      }
  toCbasis (SpinfulFermionicBasis n (SpinfulPerSector u d)) =
    Cbasis
      { cbasis_number_sites = fromIntegral n,
        cbasis_number_particles = fromIntegral (u + d),
        cbasis_number_up = fromIntegral u,
        cbasis_particle_type = c_LS_HS_FERMION,
        cbasis_state_index_is_identity = fromBool False
      }
  createCbasis_kernels (SpinfulFermionicBasis n occupation) =
    case occupation of
      SpinfulNoOccupation -> createCbasis_kernels (SpinlessFermionicBasis (2 * n) Nothing)
      SpinfulTotalParticles p -> createCbasis_kernels (SpinlessFermionicBasis (2 * n) (Just p))
      SpinfulPerSector down up -> error "not yet implemented"
  destroyCbasis_kernels _ (Cbasis_kernels stateIndexKernel stateIndexEnv) = do
    autoDestroyStateIndexKernel stateIndexKernel stateIndexEnv
  minStateEstimate (SpinfulFermionicBasis n occupation) =
    case occupation of
      SpinfulNoOccupation -> minStateEstimate (SpinlessFermionicBasis (2 * n) Nothing)
      SpinfulTotalParticles p -> minStateEstimate (SpinlessFermionicBasis (2 * n) (Just p))
      SpinfulPerSector down up ->
        let (BasisState _ minDown) = minStateEstimate (SpinlessFermionicBasis n (Just down))
            (BasisState _ minUp) = minStateEstimate (SpinlessFermionicBasis n (Just up))
         in BasisState (2 * n) $ (minUp `shiftL` n) .|. minDown
  maxStateEstimate (SpinfulFermionicBasis n occupation) =
    case occupation of
      SpinfulNoOccupation -> maxStateEstimate (SpinlessFermionicBasis (2 * n) Nothing)
      SpinfulTotalParticles p -> maxStateEstimate (SpinlessFermionicBasis (2 * n) (Just p))
      SpinfulPerSector down up ->
        let (BasisState _ minDown) = maxStateEstimate (SpinlessFermionicBasis n (Just down))
            (BasisState _ minUp) = maxStateEstimate (SpinlessFermionicBasis n (Just up))
         in BasisState (2 * n) $ (minUp `shiftL` n) .|. minDown
  flattenIndex (SpinfulFermionicBasis n _) i = case i of
    (SpinUp, k) -> k
    (SpinDown, k) -> k `shiftL` n

instance IsBasis SpinlessFermionicBasis where
  type IndexType SpinlessFermionicBasis = Int
  type GeneratorType SpinlessFermionicBasis = FermionGeneratorType
  toCbasis (SpinlessFermionicBasis n (Just p)) =
    Cbasis
      { cbasis_number_sites = fromIntegral n,
        cbasis_number_particles = fromIntegral p,
        cbasis_number_up = 0,
        cbasis_particle_type = c_LS_HS_FERMION,
        cbasis_state_index_is_identity = fromBool False
      }
  toCbasis (SpinlessFermionicBasis n Nothing) =
    Cbasis
      { cbasis_number_sites = fromIntegral n,
        cbasis_number_particles = -1,
        cbasis_number_up = 0,
        cbasis_particle_type = c_LS_HS_FERMION,
        cbasis_state_index_is_identity = fromBool True
      }
  createCbasis_kernels (SpinlessFermionicBasis n occupation) = do
    kernels <- createCbasis_kernels (SpinBasis n occupation)
    pure $
      Cbasis_kernels
        { cbasis_state_index_kernel = cbasis_state_index_kernel kernels,
          cbasis_state_index_data = cbasis_state_index_data kernels
        }
  destroyCbasis_kernels _ (Cbasis_kernels stateIndexKernel stateIndexEnv) = do
    autoDestroyStateIndexKernel stateIndexKernel stateIndexEnv
  minStateEstimate (SpinlessFermionicBasis n p) = minStateEstimate (SpinBasis n p)
  maxStateEstimate (SpinlessFermionicBasis n p) = maxStateEstimate (SpinBasis n p)
