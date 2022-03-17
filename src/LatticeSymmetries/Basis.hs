{-# LANGUAGE CApiFFI #-}
{-# LANGUAGE DefaultSignatures #-}
{-# LANGUAGE GADTs #-}
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
    -- createCbasis,
    borrowCbasis,
    -- destroyCbasis,
    Cparticle_type (..),
    Cbasis_kernels (..),
    createCbasis_kernels,
    destroyCbasis_kernels,
    stateIndex,
    binomialCoefficient,
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
  SpinHeader :: !Int -> !(Maybe Int) -> BasisHeader 'SpinTy
  SpinfulFermionHeader :: !Int -> !SpinfulOccupation -> BasisHeader 'SpinfulFermionTy
  SpinlessFermionHeader :: !Int -> !(Maybe Int) -> BasisHeader 'SpinlessFermionTy

deriving stock instance Show (BasisHeader t)

deriving stock instance Eq (BasisHeader t)

data Basis t = Basis {basisHeader :: !(BasisHeader t), basisContents :: !(ForeignPtr Cbasis)}
  deriving (Show, Eq)

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
  Basis 'SpinTy
mkSpinBasis n _
  | n <= 0 = error $ "invalid number of spins: " <> show n
mkSpinBasis n (Just h)
  | h < 0 || h > n = error $ "invalid Hamming weight: " <> show h
mkSpinBasis n h = basisFromHeader $ SpinHeader n h

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

ls_hs_create_basis :: HasCallStack => Cparticle_type -> CInt -> CInt -> CInt -> IO (Ptr Cbasis)
ls_hs_create_basis particleType numberSites numberParticles numberUp
  | particleType == c_LS_HS_SPIN =
    let h = if numberUp == -1 then Nothing else Just (fromIntegral numberUp)
        basis = mkSpinBasis (fromIntegral numberSites) h
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

stateIndex :: Basis t -> BasisState -> Maybe Int
stateIndex basis (BasisState _ (BitString α))
  | α <= fromIntegral (maxBound :: Word64) = unsafePerformIO $
    bracket (createCbasis_kernels (basisHeader basis)) destroyCbasis_kernels $ \kernelsPtr ->
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

foreign import capi unsafe "lattice_symmetries_haskell.h &ls_hs_state_index_combinadics_kernel"
  ls_hs_state_index_combinadics_kernel :: FunPtr Cindex_kernel

foreign import capi unsafe "lattice_symmetries_haskell.h &ls_hs_state_index_identity_kernel"
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
flattenIndex basis i = case basisHeader basis of
  SpinHeader _ _ -> i
  SpinfulFermionHeader n _ ->
    case i of
      (SpinUp, k) -> k
      (SpinDown, k) -> n + k
  SpinlessFermionHeader _ _ -> i

getNumberSites :: BasisHeader t -> Int
getNumberSites x = case x of
  SpinHeader n _ -> n
  SpinfulFermionHeader n _ -> n
  SpinlessFermionHeader n _ -> n
{-# INLINE getNumberSites #-}

getNumberBits :: BasisHeader t -> Int
getNumberBits x = case x of
  SpinHeader n _ -> n
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
  SpinHeader n _ -> Just n
  SpinfulFermionHeader _ occupation ->
    case occupation of
      SpinfulNoOccupation -> Nothing
      SpinfulTotalParticles p -> Just p
      SpinfulPerSector u d -> Just (u + d)
  SpinlessFermionHeader _ p -> p
{-# INLINE getNumberParticles #-}

getNumberUp :: BasisHeader t -> Maybe Int
getNumberUp x = case x of
  SpinHeader _ h -> h
  SpinfulFermionHeader _ occupation ->
    case occupation of
      SpinfulNoOccupation -> Nothing
      SpinfulTotalParticles _ -> Nothing
      SpinfulPerSector u _ -> Just u
  SpinlessFermionHeader _ _ -> Nothing
{-# INLINE getNumberUp #-}

getParticleType :: BasisHeader t -> Cparticle_type
getParticleType x = case x of
  SpinHeader _ _ -> c_LS_HS_SPIN
  SpinfulFermionHeader _ _ -> c_LS_HS_SPINFUL_FERMION
  SpinlessFermionHeader _ _ -> c_LS_HS_SPINLESS_FERMION
{-# INLINE getParticleType #-}

isStateIndexIdentity :: BasisHeader t -> Bool
isStateIndexIdentity x = case x of
  SpinHeader _ Nothing -> True
  SpinfulFermionHeader _ SpinfulNoOccupation -> True
  SpinlessFermionHeader _ Nothing -> True
  _ -> False
{-# INLINE isStateIndexIdentity #-}

hasFixedHammingWeight :: BasisHeader t -> Bool
hasFixedHammingWeight x = case x of
  SpinHeader _ (Just _) -> True
  SpinfulFermionHeader _ (SpinfulTotalParticles _) -> True
  SpinlessFermionHeader _ (Just _) -> True
  _ -> False
{-# INLINE hasFixedHammingWeight #-}

maxStateEstimate :: BasisHeader t -> BasisState
maxStateEstimate x = case x of
  SpinHeader n (Just h) -> BasisState n . BitString $ (bit h - 1) `shiftL` (n - h)
  SpinHeader n Nothing -> BasisState n . BitString $ bit n - 1
  SpinfulFermionHeader n SpinfulNoOccupation -> maxStateEstimate (SpinlessFermionHeader (2 * n) Nothing)
  SpinfulFermionHeader n (SpinfulTotalParticles p) -> maxStateEstimate (SpinlessFermionHeader (2 * n) (Just p))
  SpinfulFermionHeader n (SpinfulPerSector down up) ->
    let (BasisState _ minDown) = maxStateEstimate (SpinlessFermionHeader n (Just down))
        (BasisState _ minUp) = maxStateEstimate (SpinlessFermionHeader n (Just up))
     in BasisState (2 * n) $ (minUp `shiftL` n) .|. minDown
  SpinlessFermionHeader n p -> maxStateEstimate (SpinHeader n p)

minStateEstimate :: BasisHeader t -> BasisState
minStateEstimate x = case x of
  SpinHeader n (Just h) -> BasisState n . BitString $ bit h - 1
  SpinHeader n Nothing -> BasisState n . BitString $ zeroBits
  SpinfulFermionHeader n SpinfulNoOccupation -> minStateEstimate (SpinlessFermionHeader (2 * n) Nothing)
  SpinfulFermionHeader n (SpinfulTotalParticles p) -> minStateEstimate (SpinlessFermionHeader (2 * n) (Just p))
  SpinfulFermionHeader n (SpinfulPerSector down up) ->
    let (BasisState _ minDown) = minStateEstimate (SpinlessFermionHeader n (Just down))
        (BasisState _ minUp) = minStateEstimate (SpinlessFermionHeader n (Just up))
     in BasisState (2 * n) $ (minUp `shiftL` n) .|. minDown
  SpinlessFermionHeader n p -> minStateEstimate (SpinHeader n p)

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
          cbasis_requires_projection = fromBool False,
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

-- createCbasis :: Basis t -> IO (Ptr Cbasis)
-- createCbasis basis = do
--   p <- castStablePtrToPtr <$> newStablePtr basis
--   kernels <- createCbasis_kernels basis
--   putStrLn $ "Created Cbasis: " <> show p
--   new $
--     Cbasis
--       { cbasis_number_sites = fromIntegral (getNumberSites basis),
--         cbasis_number_particles = optionalNatural (getNumberParticles basis),
--         cbasis_number_up = optionalNatural (getNumberUp basis),
--         cbasis_particle_type = getParticleType basis,
--         cbasis_state_index_is_identity = fromBool (isStateIndexIdentity basis),
--         cbasis_requires_projection = fromBool False,
--         cbasis_kernels = kernels,
--         cbasis_haskell_payload = p
--       }

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

-- destroyCbasis :: Ptr Cbasis -> IO ()
-- destroyCbasis p
--   | p == nullPtr = pure ()
--   | otherwise = do
--     b <- peek p
--     destroyCbasis_kernels (cbasis_kernels b)
--     let freeHaskellPayload :: forall (t :: ParticleTy). Proxy t -> IO ()
--         freeHaskellPayload _ = freeStablePtr (castPtrToStablePtr (cbasis_haskell_payload b) :: StablePtr (Basis t))
--     withParticleType (cbasis_particle_type b) freeHaskellPayload
--     free p

createCbasis_kernels :: HasCallStack => BasisHeader t -> IO (Ptr Cbasis_kernels)
createCbasis_kernels x
  | getNumberBits x > 64 =
    new
      Cbasis_kernels
        { cbasis_state_index_kernel = nullFunPtr,
          cbasis_state_index_data = nullPtr
        }
  | isStateIndexIdentity x =
    new
      Cbasis_kernels
        { cbasis_state_index_kernel = ls_hs_state_index_identity_kernel,
          cbasis_state_index_data = nullPtr
        }
  | otherwise = case x of
    (SpinfulFermionHeader _ (SpinfulPerSector _ _)) -> error "not yet implemented"
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
