{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE TemplateHaskell #-}

module LatticeSymmetries.FFI where

-- ( Cparticle_type (..)
-- , c_LS_HS_SPIN
-- , c_LS_HS_SPINFUL_FERMION
-- , c_LS_HS_SPINLESS_FERMION
-- , AtomicCInt (..)
-- , Cscalar
-- , Cbasis_kernels (..)
-- , Cindex_kernel
-- , Cstate_info_kernel
-- , Cis_representative_kernel
-- , Cpermutation_group (..)
-- , Cbasis (..)
-- , basisIncRefCount
-- , basisDecRefCount
-- , basisPeekParticleType
-- , basisPeekPayload
-- , basisPokePayload
-- , Creplace_index
-- , mkCreplace_index
-- , Cnonbranching_terms (..)
-- , Coperator (..)
-- , operatorIncRefCount
-- , operatorDecRefCount
-- , Cexternal_array (..)
-- , emptyExternalArray
-- , Cyaml_config (..)
-- ) where

import Foreign
import Foreign.C.Types
import LatticeSymmetries.ComplexRational (Cscalar)
import LatticeSymmetries.Context
import Language.C.Inline.Unsafe qualified as CU

-- newtype AtomicCInt = AtomicCInt CInt
--   deriving stock (Show, Eq)

-- newtype Cparticle_type = Cparticle_type CInt
--   deriving stock (Show, Eq)
--   deriving newtype (Storable)

-- type Cindex_kernel = CPtrdiff -> Ptr Word64 -> CPtrdiff -> Ptr CPtrdiff -> CPtrdiff -> Ptr () -> IO ()

-- type Cstate_info_kernel =
--   CPtrdiff
--   -> Ptr Word64
--   -> CPtrdiff
--   -> Ptr Word64
--   -> CPtrdiff
--   -> Ptr Cscalar
--   -> Ptr CDouble
--   -> Ptr ()
--   -> IO ()

-- type Cis_representative_kernel =
--   CPtrdiff
--   -> Ptr Word64
--   -> CPtrdiff
--   -> Ptr Word8
--   -> Ptr CDouble
--   -> Ptr ()
--   -> IO ()

-- data Cbasis_kernels = Cbasis_kernels
--   { cbasis_state_info_kernel :: {-# UNPACK #-} !(FunPtr Cstate_info_kernel)
--   , cbasis_state_info_data :: {-# UNPACK #-} !(Ptr ())
--   , cbasis_is_representative_kernel :: {-# UNPACK #-} !(FunPtr Cis_representative_kernel)
--   , cbasis_is_representative_data :: {-# UNPACK #-} !(Ptr ())
--   , cbasis_state_index_kernel :: {-# UNPACK #-} !(FunPtr Cindex_kernel)
--   , cbasis_state_index_data :: {-# UNPACK #-} !(Ptr ())
--   }
--   deriving stock (Show, Generic)
--
-- data Cpermutation_group = Cpermutation_group
--   { cpermutation_group_refcount :: {-# UNPACK #-} !CInt
--   , cpermutation_group_number_bits :: {-# UNPACK #-} !CInt
--   , cpermutation_group_number_shifts :: {-# UNPACK #-} !CInt
--   , cpermutation_group_number_masks :: {-# UNPACK #-} !CInt
--   , cpermutation_group_masks :: {-# UNPACK #-} !(Ptr Word64)
--   , cpermutation_group_shifts :: {-# UNPACK #-} !(Ptr Word64)
--   , cpermutation_group_eigvals_re :: {-# UNPACK #-} !(Ptr CDouble)
--   , cpermutation_group_eigvals_im :: {-# UNPACK #-} !(Ptr CDouble)
--   , cpermutation_group_haskell_payload :: {-# UNPACK #-} !(Ptr ())
--   }
--   deriving stock (Show, Generic)

-- data Cbasis = Cbasis
--   { cbasis_refcount :: {-# UNPACK #-} !AtomicCInt
--   , cbasis_number_sites :: {-# UNPACK #-} !CInt
--   , cbasis_number_particles :: {-# UNPACK #-} !CInt
--   , cbasis_number_up :: {-# UNPACK #-} !CInt
--   , cbasis_particle_type :: {-# UNPACK #-} !Cparticle_type
--   , cbasis_spin_inversion :: {-# UNPACK #-} !CInt
--   , cbasis_state_index_is_identity :: {-# UNPACK #-} !CBool
--   , cbasis_requires_projection :: {-# UNPACK #-} !CBool
--   , cbasis_kernels :: {-# UNPACK #-} !(Ptr Cbasis_kernels)
--   , cbasis_representatives :: {-# UNPACK #-} !Cexternal_array
--   , cbasis_haskell_payload :: {-# UNPACK #-} !(Ptr ())
--   }
--   deriving stock (Show, Generic)

data Cnonbranching_terms = Cnonbranching_terms
  { cnonbranching_terms_number_terms :: {-# UNPACK #-} !CInt
  , cnonbranching_terms_number_bits :: {-# UNPACK #-} !CInt
  , cnonbranching_terms_v :: {-# UNPACK #-} !(Ptr Cscalar)
  , cnonbranching_terms_m :: {-# UNPACK #-} !(Ptr Word64)
  , cnonbranching_terms_l :: {-# UNPACK #-} !(Ptr Word64)
  , cnonbranching_terms_r :: {-# UNPACK #-} !(Ptr Word64)
  , cnonbranching_terms_x :: {-# UNPACK #-} !(Ptr Word64)
  , cnonbranching_terms_s :: {-# UNPACK #-} !(Ptr Word64)
  }
  deriving stock (Show, Generic)

-- type Creplace_index = CInt -> CInt -> Ptr CInt -> Ptr CInt -> IO ()

-- data Coperator = Coperator
--   { coperator_refcount :: {-# UNPACK #-} !AtomicCInt
--   , coperator_basis :: {-# UNPACK #-} !(Ptr Cbasis)
--   , coperator_off_diag_terms :: {-# UNPACK #-} !(Ptr Cnonbranching_terms)
--   , coperator_diag_terms :: {-# UNPACK #-} !(Ptr Cnonbranching_terms)
--   , coperator_haskell_payload :: {-# UNPACK #-} !(Ptr ())
--   }
--   deriving stock (Show, Generic)

-- data Cexternal_array = Cexternal_array
--   { external_array_elts :: !(Ptr ())
--   , external_array_num_elts :: !Word64
--   , external_array_freer :: !(Ptr ())
--   }
--   deriving stock (Show, Generic)

-- data Cyaml_config = Cyaml_config
--   { cyaml_config_basis :: {-# UNPACK #-} !(Ptr Cbasis)
--   , cyaml_config_hamiltonian :: {-# UNPACK #-} !(Ptr Coperator)
--   , cyaml_config_number_observables :: {-# UNPACK #-} !CInt
--   , cyaml_config_observables :: {-# UNPACK #-} !(Ptr (Ptr Coperator))
--   }
--   deriving stock (Show, Generic)

importLS

ls_hs_destroy_object :: forall a. IsCobject a => (Ptr a -> IO ()) -> Ptr a -> IO ()
ls_hs_destroy_object custom p = do
  let p' = castPtr @a @Cobject p
  refcount <- [CU.exp| int { ls_hs_internal_object_dec_ref_count($(ls_hs_object* p')) } |]
  when (refcount == 1) $ do
    custom p
    freeStablePtr . castPtrToStablePtr =<< [CU.exp| void* { $(ls_hs_object* p')->haskell_payload } |]

newCobject :: IsCobject b => a -> IO (Ptr b)
newCobject = fmap castPtr . newCobjectSized 0

newCobjectSized :: Int -> a -> IO (Ptr Cobject)
newCobjectSized sizeBytes x = do
  payload <- castStablePtrToPtr <$> newStablePtr x
  p <- callocBytes $ max (fromIntegral [CU.pure| size_t { sizeof(ls_hs_object) } |]) sizeBytes
  [CU.block| void { ls_hs_internal_object_init($(ls_hs_object* p), 1, $(void* payload)); } |]
  pure p

withCobject :: forall b a r. IsCobject b => Ptr b -> (a -> IO r) -> IO r
withCobject (castPtr @b @Cobject -> p) f =
  f
    =<< (deRefStablePtr . castPtrToStablePtr)
    =<< [CU.exp| void* { $(ls_hs_object* p)->haskell_payload } |]

foldCobject :: forall b a r. IsCobject b => (a -> IO r) -> Ptr b -> IO r
foldCobject f p = withCobject p f

-- instance GStorable AtomicCInt where
--   gsizeOf _ = sizeOf (undefined :: CInt)
--   galignment _ = alignment (undefined :: CInt)
--   gpeekByteOff (castPtr @_ @CChar -> ptr) (fromIntegral -> offset) =
--     AtomicCInt
--       <$> [CU.exp| int { atomic_load($(char* ptr) + $(int offset)) } |]
--   gpokeByteOff (castPtr @_ @CChar -> ptr) (fromIntegral -> offset) (AtomicCInt value) =
--     [CU.exp| void { atomic_store($(char* ptr) + $(int offset), $(int value)) } |]
--
-- c_LS_HS_SPIN :: Cparticle_type
-- c_LS_HS_SPIN = Cparticle_type [CU.pure| int { LS_HS_SPIN } |]
--
-- c_LS_HS_SPINFUL_FERMION :: Cparticle_type
-- c_LS_HS_SPINFUL_FERMION = Cparticle_type [CU.pure| int { LS_HS_SPINFUL_FERMION } |]
--
-- c_LS_HS_SPINLESS_FERMION :: Cparticle_type
-- c_LS_HS_SPINLESS_FERMION = Cparticle_type [CU.pure| int { LS_HS_SPINLESS_FERMION } |]
--
-- basisIncRefCount :: Ptr Cbasis -> IO CInt
-- basisIncRefCount p =
--   [CU.exp| int { atomic_fetch_add(&($(ls_hs_basis* p)->refcount), 1) } |]
--
-- basisDecRefCount :: Ptr Cbasis -> IO CInt
-- basisDecRefCount p =
--   [CU.exp| int { atomic_fetch_sub(&($(ls_hs_basis* p)->refcount), 1) } |]
--
-- basisPeekParticleType :: Ptr Cbasis -> IO Cparticle_type
-- basisPeekParticleType p =
--   Cparticle_type
--     <$> [CU.exp| int { $(ls_hs_basis const* p)->particle_type } |]
--
-- basisPeekPayload :: Ptr Cbasis -> IO (StablePtr a)
-- basisPeekPayload p =
--   castPtrToStablePtr
--     <$> [CU.exp| void* { $(ls_hs_basis const* p)->haskell_payload } |]
--
-- basisPokePayload :: Ptr Cbasis -> StablePtr a -> IO ()
-- basisPokePayload p (castStablePtrToPtr -> x) =
--   [CU.block| void { $(ls_hs_basis* p)->haskell_payload = $(void* x); } |]
--
-- operatorIncRefCount :: Ptr Coperator -> IO CInt
-- operatorIncRefCount p =
--   [CU.exp| int { atomic_fetch_add(&($(ls_hs_operator* p)->refcount), 1) } |]
--
-- operatorDecRefCount :: Ptr Coperator -> IO CInt
-- operatorDecRefCount p =
--   [CU.exp| int { atomic_fetch_sub(&($(ls_hs_operator* p)->refcount), 1) } |]
--
-- emptyExternalArray :: Cexternal_array
-- emptyExternalArray = Cexternal_array nullPtr 0 nullPtr
--
-- foreign import ccall "dynamic"
--   mkCreplace_index :: FunPtr Creplace_index -> Creplace_index
--
-- deriving instance GStorable Cbasis_kernels
--
-- deriving instance GStorable Cbasis
--
-- deriving instance GStorable Cnonbranching_terms
--
-- deriving instance GStorable Coperator
--
-- deriving instance GStorable Cexternal_array
--
-- deriving instance GStorable Cpermutation_group
--
-- deriving instance GStorable Cyaml_config
