{-# LANGUAGE CApiFFI #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE OverloadedLists #-}
{-# LANGUAGE QuantifiedConstraints #-}
{-# LANGUAGE RankNTypes #-}

module LatticeSymmetries.Basis
  ( mkSpinBasis,
    mkSpinlessFermionicBasis,
    mkSpinfulFermionicBasis,
    BasisState (..),
    Basis (..),
    BasisHeader (..),
    Factor,
    ParticleTy (..),
    IndexType (..),
    GeneratorType (..),
    SpinfulOccupation (..),
    Cbasis (..),
    IsBasis (..),
    borrowCbasis,
    Cparticle_type (..),
    Cbasis_kernels (..),
    createCbasis_kernels,
    destroyCbasis_kernels,
    stateIndex,
    flattenIndex,
    getNumberBits,
    getNumberWords,
    isStateIndexIdentity,
    hasFixedHammingWeight,
    minStateEstimate,
    maxStateEstimate,
    withParticleType,
    withReconstructedBasis,
    ls_hs_create_basis,
    ls_hs_min_state_estimate,
    ls_hs_max_state_estimate,
    ls_hs_spin_chain_10_basis,
  )
where

import Control.Exception.Safe (bracket)
import Data.Bits
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
import GHC.ForeignPtr
import LatticeSymmetries.Algebra (Algebra (..))
import LatticeSymmetries.BitString
import LatticeSymmetries.FFI
import LatticeSymmetries.Generator
import LatticeSymmetries.Group
import LatticeSymmetries.IO
import LatticeSymmetries.NonbranchingTerm
import System.IO (hPutStrLn, stderr)
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

data BasisHeader (t :: ParticleTy) where
  SpinHeader :: !Int -> !(Maybe Int) -> !(Maybe Int) -> !Symmetries -> BasisHeader 'SpinTy
  SpinfulFermionHeader :: !Int -> !SpinfulOccupation -> BasisHeader 'SpinfulFermionTy
  SpinlessFermionHeader :: !Int -> !(Maybe Int) -> BasisHeader 'SpinlessFermionTy

deriving stock instance Show (BasisHeader t)

deriving stock instance Eq (BasisHeader t)

data Basis t = Basis {basisHeader :: !(BasisHeader t), basisContents :: !(ForeignPtr Cbasis)}
  deriving (Show)

instance Eq (Basis t) where
  (==) a b = basisHeader a == basisHeader b

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

mkSpinBasis ::
  HasCallStack =>
  -- | Number of sites
  Int ->
  -- | Hamming weight
  Maybe Int ->
  -- | Spin inversion
  Maybe Int ->
  -- | Lattice symmetries
  Symmetries ->
  Basis 'SpinTy
mkSpinBasis n _ _ _
  | n <= 0 = error $ "invalid number of spins: " <> show n
mkSpinBasis n (Just h) _ _
  | h < 0 || h > n = error $ "invalid Hamming weight: " <> show h
mkSpinBasis _ _ (Just i) _
  | i /= (-1) && i /= 1 = error $ "invalid spin inversion: " <> show i
mkSpinBasis n (Just h) (Just _) _
  | n /= 2 * h = error $ "invalid spin inversion: " <> show n <> " spins, but Hamming weight is " <> show h
mkSpinBasis n _ _ g
  | not (nullSymmetries g) && symmetriesGetNumberBits g /= n = error "invalid symmetries"
mkSpinBasis n h i g = basisFromHeader $ SpinHeader n h i g

mkSpinfulFermionicBasis ::
  HasCallStack =>
  -- | Number of sites
  Int ->
  -- | Number of particles
  SpinfulOccupation ->
  Basis 'SpinfulFermionTy
mkSpinfulFermionicBasis n _
  | n <= 0 = error $ "invalid number of sites: " <> show n
mkSpinfulFermionicBasis n (SpinfulTotalParticles p)
  | p < 0 || p > n = error $ "invalid number of particles: " <> show p
mkSpinfulFermionicBasis n (SpinfulPerSector u d)
  | u < 0 || u > n || d < 0 || d > n =
    error $
      "invalid number of particles: " <> show u <> " up, " <> show d <> " down"
mkSpinfulFermionicBasis n occupation = basisFromHeader $ SpinfulFermionHeader n occupation

mkSpinlessFermionicBasis ::
  HasCallStack =>
  -- | Number of sites
  Int ->
  -- | Number of particles
  Maybe Int ->
  Basis 'SpinlessFermionTy
mkSpinlessFermionicBasis n _
  | n <= 0 = error $ "invalid number of sites: " <> show n
mkSpinlessFermionicBasis n (Just h)
  | h < 0 || h > n = error $ "invalid number of particles: " <> show h
mkSpinlessFermionicBasis n h = basisFromHeader $ SpinlessFermionHeader n h

ls_hs_spin_chain_10_basis :: IO (Ptr Cbasis)
ls_hs_spin_chain_10_basis = do
  let symmetries =
        mkSymmetries
          [ mkSymmetry (mkPermutation [1, 2, 3, 4, 5, 6, 7, 8, 9, 0]) 0,
            mkSymmetry (mkPermutation [9, 8, 7, 6, 5, 4, 3, 2, 1, 0]) 0
          ]
      basis = mkSpinBasis 10 (Just 5) Nothing symmetries
  borrowCbasis basis

ls_hs_create_basis :: HasCallStack => Cparticle_type -> CInt -> CInt -> CInt -> IO (Ptr Cbasis)
ls_hs_create_basis particleType numberSites numberParticles numberUp
  | particleType == c_LS_HS_SPIN =
    let h = if numberUp == -1 then Nothing else Just (fromIntegral numberUp)
        basis = mkSpinBasis (fromIntegral numberSites) h Nothing emptySymmetries
     in borrowCbasis basis
  | particleType == c_LS_HS_SPINFUL_FERMION =
    let p
          | numberParticles == -1 = SpinfulNoOccupation
          | numberUp == -1 = SpinfulTotalParticles (fromIntegral numberParticles)
          | otherwise = SpinfulPerSector (fromIntegral numberUp) (fromIntegral (numberParticles - numberUp))
        basis = mkSpinfulFermionicBasis (fromIntegral numberSites) p
     in borrowCbasis basis
  | particleType == c_LS_HS_SPINLESS_FERMION =
    let p = if numberParticles == -1 then Nothing else Just (fromIntegral numberParticles)
        basis = mkSpinlessFermionicBasis (fromIntegral numberSites) p
     in borrowCbasis basis
  | otherwise = error $ "invalid particle type: " <> show particleType

ls_hs_max_state_estimate :: HasCallStack => Ptr Cbasis -> IO Word64
ls_hs_max_state_estimate p =
  withReconstructedBasis p $ \basis ->
    let (BasisState n (BitString x)) = maxStateEstimate (basisHeader basis)
     in if n > 64
          then error "maximal state is not representable as a 64-bit integer"
          else pure $ fromIntegral x

ls_hs_min_state_estimate :: HasCallStack => Ptr Cbasis -> IO Word64
ls_hs_min_state_estimate p =
  withReconstructedBasis p $ \basis ->
    let (BasisState n (BitString x)) = minStateEstimate (basisHeader basis)
     in if n > 64
          then error "minimal state is not representable as a 64-bit integer"
          else pure $ fromIntegral x

newtype ChapelArray = ChapelArray (ForeignPtr Cexternal_array)

data Representatives = Representatives {rStates :: !ChapelArray, rRanges :: !ChapelArray}

buildRepresentatives :: Basis t -> Representatives
buildRepresentatives = undefined

stateIndex :: Basis t -> BasisState -> Maybe Int
stateIndex basis (BasisState _ (BitString α))
  | α <= fromIntegral (maxBound :: Word64) = unsafePerformIO . withCbasis basis $ \basisPtr ->
    with (fromIntegral α :: Word64) $ \spinsPtr ->
      alloca $ \indicesPtr -> do
        (f, env) <- basisPeekStateIndexKernel basisPtr
        if f == nullFunPtr
          then pure Nothing
          else do
            mkCindex_kernel f 1 spinsPtr 1 indicesPtr 1 env
            i <- peek indicesPtr
            if i < 0
              then pure Nothing
              else pure $ Just (fromIntegral i)
  | otherwise = Nothing

flattenIndex :: Basis t -> IndexType t -> Int
flattenIndex basis i = case basisHeader basis of
  SpinHeader _ _ _ _ -> i
  SpinfulFermionHeader n _ ->
    case i of
      (SpinUp, k) -> k
      (SpinDown, k) -> n + k
  SpinlessFermionHeader _ _ -> i

getNumberSites :: BasisHeader t -> Int
getNumberSites x = case x of
  SpinHeader n _ _ _ -> n
  SpinfulFermionHeader n _ -> n
  SpinlessFermionHeader n _ -> n
{-# INLINE getNumberSites #-}

getNumberBits :: BasisHeader t -> Int
getNumberBits x = case x of
  SpinHeader n _ _ _ -> n
  SpinfulFermionHeader n _ -> 2 * n
  SpinlessFermionHeader n _ -> n
{-# INLINE getNumberBits #-}

getNumberWords :: BasisHeader t -> Int
getNumberWords x = (n + 63) `div` 64
  where
    n = getNumberBits x
{-# INLINE getNumberWords #-}

getNumberParticles :: BasisHeader t -> Maybe Int
getNumberParticles x = case x of
  SpinHeader n _ _ _ -> Just n
  SpinfulFermionHeader _ occupation ->
    case occupation of
      SpinfulNoOccupation -> Nothing
      SpinfulTotalParticles p -> Just p
      SpinfulPerSector u d -> Just (u + d)
  SpinlessFermionHeader _ p -> p
{-# INLINE getNumberParticles #-}

getNumberUp :: BasisHeader t -> Maybe Int
getNumberUp x = case x of
  SpinHeader _ h _ _ -> h
  SpinfulFermionHeader _ occupation ->
    case occupation of
      SpinfulNoOccupation -> Nothing
      SpinfulTotalParticles _ -> Nothing
      SpinfulPerSector u _ -> Just u
  SpinlessFermionHeader _ _ -> Nothing
{-# INLINE getNumberUp #-}

getParticleType :: BasisHeader t -> Cparticle_type
getParticleType x = case x of
  SpinHeader _ _ _ _ -> c_LS_HS_SPIN
  SpinfulFermionHeader _ _ -> c_LS_HS_SPINFUL_FERMION
  SpinlessFermionHeader _ _ -> c_LS_HS_SPINLESS_FERMION
{-# INLINE getParticleType #-}

isStateIndexIdentity :: BasisHeader t -> Bool
isStateIndexIdentity x = case x of
  SpinHeader _ Nothing Nothing g -> nullSymmetries g
  SpinfulFermionHeader _ SpinfulNoOccupation -> True
  SpinlessFermionHeader _ Nothing -> True
  _ -> False
{-# INLINE isStateIndexIdentity #-}

requiresProjection :: BasisHeader t -> Bool
requiresProjection x = case x of
  SpinHeader _ _ i g -> not (nullSymmetries g) || isJust i
  SpinfulFermionHeader _ _ -> False
  SpinlessFermionHeader _ _ -> False
{-# INLINE requiresProjection #-}

hasFixedHammingWeight :: BasisHeader t -> Bool
hasFixedHammingWeight x = case x of
  SpinHeader _ (Just _) _ _ -> True
  SpinfulFermionHeader _ (SpinfulTotalParticles _) -> True
  SpinlessFermionHeader _ (Just _) -> True
  _ -> False
{-# INLINE hasFixedHammingWeight #-}

maxStateEstimate :: BasisHeader t -> BasisState
maxStateEstimate x = case x of
  SpinHeader n (Just h) i _ ->
    BasisState n . BitString $
      if isNothing i
        then (bit h - 1) `shiftL` (n - h)
        else (bit h - 1) `shiftL` (n - h - 1)
  -- TODO: improve the bound for fixed spin inversion
  SpinHeader n Nothing _ _ -> BasisState n . BitString $ bit n - 1
  SpinfulFermionHeader n SpinfulNoOccupation -> maxStateEstimate (SpinlessFermionHeader (2 * n) Nothing)
  SpinfulFermionHeader n (SpinfulTotalParticles p) -> maxStateEstimate (SpinlessFermionHeader (2 * n) (Just p))
  SpinfulFermionHeader n (SpinfulPerSector down up) ->
    let (BasisState _ minDown) = maxStateEstimate (SpinlessFermionHeader n (Just down))
        (BasisState _ minUp) = maxStateEstimate (SpinlessFermionHeader n (Just up))
     in BasisState (2 * n) $ (minUp `shiftL` n) .|. minDown
  SpinlessFermionHeader n p -> maxStateEstimate (SpinHeader n p Nothing emptySymmetries)

minStateEstimate :: BasisHeader t -> BasisState
minStateEstimate x = case x of
  SpinHeader n (Just h) _ _ -> BasisState n . BitString $ bit h - 1
  -- TODO: improve the bound for fixed spin inversion
  SpinHeader n Nothing _ _ -> BasisState n . BitString $ zeroBits
  SpinfulFermionHeader n SpinfulNoOccupation -> minStateEstimate (SpinlessFermionHeader (2 * n) Nothing)
  SpinfulFermionHeader n (SpinfulTotalParticles p) -> minStateEstimate (SpinlessFermionHeader (2 * n) (Just p))
  SpinfulFermionHeader n (SpinfulPerSector down up) ->
    let (BasisState _ minDown) = minStateEstimate (SpinlessFermionHeader n (Just down))
        (BasisState _ minUp) = minStateEstimate (SpinlessFermionHeader n (Just up))
     in BasisState (2 * n) $ (minUp `shiftL` n) .|. minDown
  SpinlessFermionHeader n p -> minStateEstimate (SpinHeader n p Nothing emptySymmetries)

optionalNatural :: Maybe Int -> CInt
optionalNatural x = case x of
  (Just n) -> fromIntegral n
  Nothing -> (-1)
{-# INLINE optionalNatural #-}

basisFromHeader :: BasisHeader t -> Basis t
basisFromHeader x = unsafePerformIO $ do
  fp <- mallocForeignPtr
  kernels <- createCbasis_kernels x
  addForeignPtrConcFinalizer fp (destroyCbasis_kernels kernels)
  withForeignPtr fp $ \ptr ->
    poke ptr $
      Cbasis
        { cbasis_refcount = 0,
          cbasis_number_sites = fromIntegral (getNumberSites x),
          cbasis_number_particles = optionalNatural (getNumberParticles x),
          cbasis_number_up = optionalNatural (getNumberUp x),
          cbasis_particle_type = getParticleType x,
          cbasis_state_index_is_identity = fromBool (isStateIndexIdentity x),
          cbasis_requires_projection = fromBool (requiresProjection x),
          cbasis_kernels = kernels,
          cbasis_haskell_payload = nullPtr
        }
  pure $ Basis x fp
{-# NOINLINE basisFromHeader #-}

borrowCbasis :: Basis t -> IO (Ptr Cbasis)
borrowCbasis basis = do
  withForeignPtr (basisContents basis) $ \ptr -> do
    _ <- basisIncRefCount ptr
    basisPokePayload ptr =<< newStablePtr basis
  pure $ unsafeForeignPtrToPtr (basisContents basis)

withCbasis :: Basis t -> (Ptr Cbasis -> IO a) -> IO a
withCbasis x action = withForeignPtr (basisContents x) action

withParticleType :: HasCallStack => Cparticle_type -> (forall (t :: ParticleTy). IsBasis t => Proxy t -> a) -> a
withParticleType t action
  | t == c_LS_HS_SPIN = action (Proxy :: Proxy 'SpinTy)
  | t == c_LS_HS_SPINFUL_FERMION = action (Proxy :: Proxy 'SpinfulFermionTy)
  | t == c_LS_HS_SPINLESS_FERMION = action (Proxy :: Proxy 'SpinlessFermionTy)
  | otherwise = error "invalid particle type"

withReconstructedBasis :: forall a. Ptr Cbasis -> (forall (t :: ParticleTy). IsBasis t => Basis t -> IO a) -> IO a
withReconstructedBasis p action = do
  let run :: forall (t :: ParticleTy). IsBasis t => Proxy t -> IO a
      run _ = do
        (x :: Basis t) <- deRefStablePtr =<< basisPeekPayload p
        action x
  tp <- basisPeekParticleType p
  withParticleType tp run

data Ccombinadics_kernel_data

foreign import capi unsafe "lattice_symmetries_haskell.h ls_hs_internal_create_combinadics_kernel_data"
  ls_hs_internal_create_combinadics_kernel_data :: CInt -> CBool -> IO (Ptr Ccombinadics_kernel_data)

foreign import capi unsafe "lattice_symmetries_haskell.h ls_hs_internal_destroy_combinadics_kernel_data"
  ls_hs_internal_destroy_combinadics_kernel_data :: Ptr Ccombinadics_kernel_data -> IO ()

foreign import capi unsafe "lattice_symmetries_haskell.h &ls_hs_state_index_combinadics_kernel"
  ls_hs_state_index_combinadics_kernel :: FunPtr Cindex_kernel

foreign import capi unsafe "lattice_symmetries_haskell.h &ls_hs_state_index_identity_kernel"
  ls_hs_state_index_identity_kernel :: FunPtr Cindex_kernel

data Chalide_kernel_data

foreign import capi unsafe "lattice_symmetries_haskell.h ls_internal_create_halide_kernel_data"
  ls_internal_create_halide_kernel_data :: Ptr Cpermutation_group -> CInt -> IO (Ptr Chalide_kernel_data)

foreign import capi unsafe "lattice_symmetries_haskell.h ls_internal_destroy_halide_kernel_data"
  ls_internal_destroy_halide_kernel_data :: Ptr Chalide_kernel_data -> IO ()

foreign import capi unsafe "lattice_symmetries_haskell.h &ls_hs_is_representative_halide_kernel"
  ls_hs_is_representative_halide_kernel :: FunPtr Cis_representative_kernel

foreign import capi unsafe "lattice_symmetries_haskell.h &ls_hs_state_info_halide_kernel"
  ls_hs_state_info_halide_kernel :: FunPtr Cstate_info_kernel

setStateInfoKernel :: BasisHeader t -> Cbasis_kernels -> IO Cbasis_kernels
setStateInfoKernel (SpinHeader n h i g) k
  | n <= 64 = do
    kernelData <-
      withForeignPtr (symmetriesContents g) $ \gPtr ->
        ls_internal_create_halide_kernel_data gPtr (maybe 0 fromIntegral i)
    pure $
      k
        { cbasis_state_info_kernel = ls_hs_state_info_halide_kernel,
          cbasis_state_info_data = castPtr kernelData,
          cbasis_is_representative_kernel = ls_hs_is_representative_halide_kernel,
          cbasis_is_representative_data = castPtr kernelData
        }
setStateInfoKernel x k =
  pure $
    k
      { cbasis_state_info_kernel = nullFunPtr,
        cbasis_state_info_data = nullPtr,
        cbasis_is_representative_kernel = nullFunPtr,
        cbasis_is_representative_data = nullPtr
      }

setStateIndexKernel :: BasisHeader t -> Cbasis_kernels -> IO Cbasis_kernels
setStateIndexKernel x k
  | getNumberBits x > 64 =
    pure $
      k
        { cbasis_state_index_kernel = nullFunPtr,
          cbasis_state_index_data = nullPtr
        }
  | isStateIndexIdentity x =
    pure $
      k
        { cbasis_state_index_kernel = ls_hs_state_index_identity_kernel,
          cbasis_state_index_data = nullPtr
        }
  | hasFixedHammingWeight x = do
    kernelData <-
      ls_hs_internal_create_combinadics_kernel_data
        (fromIntegral (getNumberBits x))
        (fromBool False)
    pure $
      k
        { cbasis_state_index_kernel = ls_hs_state_index_combinadics_kernel,
          cbasis_state_index_data = castPtr kernelData
        }
  | otherwise = case x of
    (SpinfulFermionHeader n (SpinfulPerSector _ _)) -> do
      kernelData <-
        ls_hs_internal_create_combinadics_kernel_data
          (fromIntegral n) -- NOTE: not 2 * n !
          (fromBool True)
      pure $
        k
          { cbasis_state_index_kernel = ls_hs_state_index_combinadics_kernel,
            cbasis_state_index_data = castPtr kernelData
          }
    _ -> error "not implemented"

createCbasis_kernels :: HasCallStack => BasisHeader t -> IO (Ptr Cbasis_kernels)
createCbasis_kernels x =
  new <=< setStateInfoKernel x <=< setStateIndexKernel x $
    Cbasis_kernels
      { cbasis_state_info_kernel = nullFunPtr,
        cbasis_state_info_data = nullPtr,
        cbasis_is_representative_kernel = nullFunPtr,
        cbasis_is_representative_data = nullPtr,
        cbasis_state_index_kernel = nullFunPtr,
        cbasis_state_index_data = nullPtr
      }

destroyCindex_kernel :: HasCallStack => Cbasis_kernels -> IO ()
destroyCindex_kernel p
  | kernel == ls_hs_state_index_combinadics_kernel = ls_hs_internal_destroy_combinadics_kernel_data (castPtr env)
  | kernel == ls_hs_state_index_identity_kernel && env == nullPtr = pure ()
  | kernel == nullFunPtr && env == nullPtr = pure ()
  | otherwise = error "failed to automatically deallocate state_index kernel"
  where
    kernel = cbasis_state_index_kernel p
    env = cbasis_state_index_data p

destroyCstate_info_kernel :: HasCallStack => Cbasis_kernels -> IO ()
destroyCstate_info_kernel p
  | cbasis_state_info_kernel p == ls_hs_state_info_halide_kernel
      && cbasis_is_representative_kernel p == ls_hs_is_representative_halide_kernel =
    do
      ls_internal_destroy_halide_kernel_data (castPtr (cbasis_state_info_data p))
      when (cbasis_state_info_data p /= cbasis_is_representative_data p) $
        ls_internal_destroy_halide_kernel_data (castPtr (cbasis_is_representative_data p))
  | cbasis_state_info_kernel p == nullFunPtr
      && cbasis_state_info_data p == nullPtr
      && cbasis_is_representative_kernel p == nullFunPtr
      && cbasis_is_representative_data p == nullPtr =
    pure ()
  | otherwise = error "failed to automatically deallocate state_info and is_representative kernels"

destroyCbasis_kernels :: HasCallStack => Ptr Cbasis_kernels -> IO ()
destroyCbasis_kernels p = do
  kernels <- peek p
  destroyCindex_kernel kernels
  destroyCstate_info_kernel kernels

basisFromSpec :: BasisSpec -> Basis 'SpinTy
basisFromSpec (BasisSpec n h i gs) = mkSpinBasis n h i (symmetriesFromSpec gs)
