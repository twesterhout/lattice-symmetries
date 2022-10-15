{-# LANGUAGE CApiFFI #-}
module LatticeSymmetries.FFI
  ( Cparticle_type (..),
    c_LS_HS_SPIN,
    c_LS_HS_SPINFUL_FERMION,
    c_LS_HS_SPINLESS_FERMION,
    Cbasis_kernels (..),
    Cindex_kernel,
    mkCindex_kernel,
    Cstate_info_kernel,
    Cis_representative_kernel,
    Cbasis (..),
    basisIncRefCount,
    basisPeekParticleType,
    basisPeekStateIndexKernel,
    basisPeekPayload,
    basisPokePayload,
    Cexpr (..),
    Creplace_index,
    mkCreplace_index,
    exprIncRefCount,
    exprDecRefCount,
    exprPeekPayload,
    exprPokePayload,
    Coperator_kernel_data,
    Cnonbranching_terms (..),
    Coperator (..),
    Cscalar,
    operatorIncRefCount,
    operatorPeekParticleType,
    operatorPeekNumberOffDiagTerms,
    operatorPeekPayload,
    operatorPokePayload,
    Cexternal_array (..),
    emptyExternalArray,
    -- ls_hs_internal_destroy_external_array,
    Cpermutation_group (..),
    Cchpl_kernels (..),
    Cyaml_config (..),
    -- ls_hs_internal_get_chpl_kernels,
  ) where

import Foreign
import Foreign.C.Types
import LatticeSymmetries.ComplexRational (Cscalar)

#include "lattice_symmetries_haskell.h"

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

type Cstate_info_kernel = CPtrdiff -> Ptr Word64 -> CPtrdiff ->
                                      Ptr Word64 -> CPtrdiff ->
                                      Ptr Cscalar -> Ptr CDouble ->
                                      Ptr () -> IO ()

-- foreign import ccall "dynamic"
--   mkCstate_info_kernel :: FunPtr Cstate_info_kernel -> Cstate_info_kernel

type Cis_representative_kernel = CPtrdiff -> Ptr Word64 -> CPtrdiff ->
                                             Ptr Word8 ->
                                             Ptr CDouble ->
                                             Ptr () -> IO ()

-- foreign import ccall "dynamic"
--   mkCis_representative_kernel :: FunPtr Cis_representative_kernel -> Cis_representative_kernel

data {-# CTYPE "lattice_symmetries_haskell.h" "ls_hs_basis_kernels" #-} Cbasis_kernels = Cbasis_kernels
  { cbasis_state_info_kernel :: {-# UNPACK #-} !(FunPtr Cstate_info_kernel),
    cbasis_state_info_data :: {-# UNPACK #-} !(Ptr ()),
    cbasis_is_representative_kernel :: {-# UNPACK #-} !(FunPtr Cis_representative_kernel),
    cbasis_is_representative_data :: {-# UNPACK #-} !(Ptr ()),
    cbasis_state_index_kernel :: {-# UNPACK #-} !(FunPtr Cindex_kernel),
    cbasis_state_index_data :: {-# UNPACK #-} !(Ptr ())
  }

instance Storable Cbasis_kernels where
  sizeOf _ = #{size ls_hs_basis_kernels}
  {-# INLINE sizeOf #-}
  alignment _ = #{alignment ls_hs_basis_kernels}
  {-# INLINE alignment #-}
  peek p =
    Cbasis_kernels
      <$> #{peek ls_hs_basis_kernels, state_info_kernel} p
      <*> #{peek ls_hs_basis_kernels, state_info_data} p
      <*> #{peek ls_hs_basis_kernels, is_representative_kernel} p
      <*> #{peek ls_hs_basis_kernels, is_representative_data} p
      <*> #{peek ls_hs_basis_kernels, state_index_kernel} p
      <*> #{peek ls_hs_basis_kernels, state_index_data} p
  {-# INLINE peek #-}
  poke p x = do
    #{poke ls_hs_basis_kernels, state_info_kernel} p (cbasis_state_info_kernel x)
    #{poke ls_hs_basis_kernels, state_info_data} p (cbasis_state_info_data x)
    #{poke ls_hs_basis_kernels, is_representative_kernel} p (cbasis_is_representative_kernel x)
    #{poke ls_hs_basis_kernels, is_representative_data} p (cbasis_is_representative_data x)
    #{poke ls_hs_basis_kernels, state_index_kernel} p (cbasis_state_index_kernel x)
    #{poke ls_hs_basis_kernels, state_index_data} p (cbasis_state_index_data x)
  {-# INLINE poke #-}

data {-# CTYPE "lattice_symmetries_haskell.h" "ls_hs_basis" #-} Cbasis = Cbasis
  { cbasis_refcount :: {-# UNPACK #-} !CInt,
    cbasis_number_sites :: {-# UNPACK #-} !CInt,
    cbasis_number_particles :: {-# UNPACK #-} !CInt,
    cbasis_number_up :: {-# UNPACK #-} !CInt,
    cbasis_particle_type :: {-# UNPACK #-} !Cparticle_type,
    cbasis_state_index_is_identity :: {-# UNPACK #-} !CBool,
    cbasis_requires_projection :: {-# UNPACK #-} !CBool,
    cbasis_kernels :: {-# UNPACK #-} !(Ptr Cbasis_kernels),
    cbasis_representatives :: {-# UNPACK #-} !Cexternal_array,
    cbasis_haskell_payload :: {-# UNPACK #-} !(Ptr ())
  }

foreign import ccall safe "ls_hs_internal_read_refcount"
  peekRefCount :: Ptr CInt -> IO CInt

foreign import ccall safe "ls_hs_internal_write_refcount"
  pokeRefCount :: Ptr CInt -> CInt -> IO CInt

foreign import ccall safe "ls_hs_internal_inc_refcount"
  incRefCount :: Ptr CInt -> IO CInt

foreign import ccall unsafe "ls_hs_internal_dec_refcount"
  decRefCount :: Ptr CInt -> IO CInt

basisIncRefCount :: Ptr Cbasis -> IO CInt
basisIncRefCount p = incRefCount (p `plusPtr` #{offset ls_hs_basis, refcount})
{-# INLINE basisIncRefCount #-}

basisPeekParticleType :: Ptr Cbasis -> IO Cparticle_type
basisPeekParticleType p = #{peek ls_hs_basis, particle_type} p
{-# INLINE basisPeekParticleType #-}

basisPeekStateIndexKernel :: Ptr Cbasis -> IO (FunPtr Cindex_kernel, Ptr ())
basisPeekStateIndexKernel p = do
  kernels <- #{peek ls_hs_basis, kernels} p
  f <- #{peek ls_hs_basis_kernels, state_index_kernel} kernels
  env <- #{peek ls_hs_basis_kernels, state_index_data} kernels
  pure (f, env)
{-# INLINE basisPeekStateIndexKernel #-}

basisPeekPayload :: Ptr Cbasis -> IO (StablePtr a)
basisPeekPayload p = #{peek ls_hs_basis, haskell_payload} p
{-# INLINE basisPeekPayload #-}

basisPokePayload :: Ptr Cbasis -> StablePtr a -> IO ()
basisPokePayload p x = #{poke ls_hs_basis, haskell_payload} p x
{-# INLINE basisPokePayload #-}

instance Storable Cbasis where
  sizeOf _ = #{size ls_hs_basis}
  {-# INLINE sizeOf #-}
  alignment _ = #{alignment ls_hs_basis}
  {-# INLINE alignment #-}
  peek p =
    Cbasis
      <$> peekRefCount (p `plusPtr` #{offset ls_hs_basis, refcount})
      <*> #{peek ls_hs_basis, number_sites} p
      <*> #{peek ls_hs_basis, number_particles} p
      <*> #{peek ls_hs_basis, number_up} p
      <*> #{peek ls_hs_basis, particle_type} p
      <*> #{peek ls_hs_basis, state_index_is_identity} p
      <*> #{peek ls_hs_basis, requires_projection} p
      <*> #{peek ls_hs_basis, kernels} p
      <*> #{peek ls_hs_basis, representatives} p
      <*> #{peek ls_hs_basis, haskell_payload} p
  {-# INLINE peek #-}
  poke p x = do
    _ <- pokeRefCount (p `plusPtr` #{offset ls_hs_basis, refcount}) (cbasis_refcount x)
    #{poke ls_hs_basis, number_sites}            p (cbasis_number_sites x)
    #{poke ls_hs_basis, number_particles}        p (cbasis_number_particles x)
    #{poke ls_hs_basis, number_up}               p (cbasis_number_up x)
    #{poke ls_hs_basis, particle_type}           p (cbasis_particle_type x)
    #{poke ls_hs_basis, state_index_is_identity} p (cbasis_state_index_is_identity x)
    #{poke ls_hs_basis, requires_projection}     p (cbasis_requires_projection x)
    #{poke ls_hs_basis, kernels}                 p (cbasis_kernels x)
    #{poke ls_hs_basis, representatives}         p (cbasis_representatives x)
    #{poke ls_hs_basis, haskell_payload}         p (cbasis_haskell_payload x)
  {-# INLINE poke #-}

data {-# CTYPE "lattice_symmetries_haskell.h" "ls_hs_expr" #-} Cexpr = Cexpr
  { cexpr_refcount :: {-# UNPACK #-} !CInt,
    cexpr_haskell_payload :: {-# UNPACK #-} !(Ptr ())
  }

exprIncRefCount :: Ptr Cexpr -> IO CInt
exprIncRefCount p = incRefCount (p `plusPtr` #{offset ls_hs_expr, refcount})
{-# INLINE exprIncRefCount #-}

exprDecRefCount :: Ptr Cexpr -> IO CInt
exprDecRefCount p = decRefCount (p `plusPtr` #{offset ls_hs_expr, refcount})
{-# INLINE exprDecRefCount #-}

exprPeekPayload :: Ptr Cexpr -> IO (StablePtr a)
exprPeekPayload p = #{peek ls_hs_expr, haskell_payload} p
{-# INLINE exprPeekPayload #-}

exprPokePayload :: Ptr Cexpr -> StablePtr a -> IO ()
exprPokePayload p x = #{poke ls_hs_expr, haskell_payload} p x
{-# INLINE exprPokePayload #-}

instance Storable Cexpr where
  sizeOf _ = #{size ls_hs_expr}
  {-# INLINE sizeOf #-}
  alignment _ = #{alignment ls_hs_expr}
  {-# INLINE alignment #-}
  peek p =
    Cexpr
      <$> peekRefCount (p `plusPtr` #{offset ls_hs_expr, refcount})
      <*> #{peek ls_hs_expr, haskell_payload} p
  {-# INLINE peek #-}
  poke p x = do
    _ <- pokeRefCount (p `plusPtr` #{offset ls_hs_expr, refcount}) (cexpr_refcount x)
    #{poke ls_hs_expr, haskell_payload} p (cexpr_haskell_payload x)
  {-# INLINE poke #-}

type Creplace_index = CInt -> CInt -> Ptr CInt -> Ptr CInt -> IO ()

foreign import ccall "dynamic"
  mkCreplace_index :: FunPtr Creplace_index -> Creplace_index

data {-# CTYPE "lattice_symmetries_haskell.h" "ls_hs_operator" #-} Cnonbranching_terms =
  Cnonbranching_terms
    { cnonbranching_terms_number_terms :: {-# UNPACK #-} !CInt,
      cnonbranching_terms_number_bits :: {-# UNPACK #-} !CInt,
      cnonbranching_terms_v :: {-# UNPACK #-} !(Ptr Cscalar),
      cnonbranching_terms_m :: {-# UNPACK #-} !(Ptr Word64),
      cnonbranching_terms_l :: {-# UNPACK #-} !(Ptr Word64),
      cnonbranching_terms_r :: {-# UNPACK #-} !(Ptr Word64),
      cnonbranching_terms_x :: {-# UNPACK #-} !(Ptr Word64),
      cnonbranching_terms_s :: {-# UNPACK #-} !(Ptr Word64)
    }
  deriving stock (Show)

instance Storable Cnonbranching_terms where
  sizeOf _ = #{size ls_hs_nonbranching_terms}
  {-# INLINE sizeOf #-}
  alignment _ = #{alignment ls_hs_nonbranching_terms}
  {-# INLINE alignment #-}
  peek p =
    Cnonbranching_terms
      <$> #{peek ls_hs_nonbranching_terms, number_terms} p
      <*> #{peek ls_hs_nonbranching_terms, number_bits} p
      <*> #{peek ls_hs_nonbranching_terms, v} p
      <*> #{peek ls_hs_nonbranching_terms, m} p
      <*> #{peek ls_hs_nonbranching_terms, l} p
      <*> #{peek ls_hs_nonbranching_terms, r} p
      <*> #{peek ls_hs_nonbranching_terms, x} p
      <*> #{peek ls_hs_nonbranching_terms, s} p
  {-# INLINE peek #-}
  poke p x = do
    #{poke ls_hs_nonbranching_terms, number_terms} p (cnonbranching_terms_number_terms x)
    #{poke ls_hs_nonbranching_terms, number_bits} p (cnonbranching_terms_number_bits x)
    #{poke ls_hs_nonbranching_terms, v} p (cnonbranching_terms_v x)
    #{poke ls_hs_nonbranching_terms, m} p (cnonbranching_terms_m x)
    #{poke ls_hs_nonbranching_terms, l} p (cnonbranching_terms_l x)
    #{poke ls_hs_nonbranching_terms, r} p (cnonbranching_terms_r x)
    #{poke ls_hs_nonbranching_terms, x} p (cnonbranching_terms_x x)
    #{poke ls_hs_nonbranching_terms, s} p (cnonbranching_terms_s x)
  {-# INLINE poke #-}

data {-# CTYPE "lattice_symmetries_haskell.h" "ls_internal_operator_kernel_data" #-}
  Coperator_kernel_data

data {-# CTYPE "lattice_symmetries_haskell.h" "ls_hs_operator" #-} Coperator = Coperator
  { coperator_refcount :: {-# UNPACK #-} !CInt,
    coperator_basis :: {-# UNPACK #-} !(Ptr Cbasis),
    coperator_off_diag_terms :: {-# UNPACK #-} !(Ptr Cnonbranching_terms),
    coperator_diag_terms :: {-# UNPACK #-} !(Ptr Cnonbranching_terms),
    coperator_apply_off_diag_cxt :: {-# UNPACK #-} !(Ptr Coperator_kernel_data),
    coperator_apply_diag_cxt :: {-# UNPACK #-} !(Ptr Coperator_kernel_data),
    coperator_haskell_payload :: {-# UNPACK #-} !(Ptr ())
  }

operatorIncRefCount :: Ptr Coperator -> IO CInt
operatorIncRefCount p = incRefCount (p `plusPtr` #{offset ls_hs_operator, refcount})
{-# INLINE operatorIncRefCount #-}

operatorPeekParticleType :: Ptr Coperator -> IO Cparticle_type
operatorPeekParticleType p = #{peek ls_hs_operator, basis} p >>= basisPeekParticleType
{-# INLINE operatorPeekParticleType #-}

operatorPeekPayload :: Ptr Coperator -> IO (StablePtr a)
operatorPeekPayload p = #{peek ls_hs_operator, haskell_payload} p
{-# INLINE operatorPeekPayload #-}

operatorPokePayload :: Ptr Coperator -> StablePtr a -> IO ()
operatorPokePayload p x = #{poke ls_hs_operator, haskell_payload} p x
{-# INLINE operatorPokePayload #-}

operatorPeekNumberOffDiagTerms :: Ptr Coperator -> IO Int
operatorPeekNumberOffDiagTerms p = do
  termsPtr <- #{peek ls_hs_operator, off_diag_terms} p
  if termsPtr == nullPtr
    then pure 0
    else #{peek ls_hs_nonbranching_terms, number_terms} termsPtr
{-# INLINE operatorPeekNumberOffDiagTerms #-}

instance Storable Coperator where
  sizeOf _ = #{size ls_hs_operator}
  {-# INLINE sizeOf #-}
  alignment _ = #{alignment ls_hs_operator}
  {-# INLINE alignment #-}
  peek p =
    Coperator
      <$> peekRefCount (p `plusPtr` #{offset ls_hs_operator, refcount})
      <*> #{peek ls_hs_operator, basis} p
      <*> #{peek ls_hs_operator, off_diag_terms} p
      <*> #{peek ls_hs_operator, diag_terms} p
      <*> #{peek ls_hs_operator, apply_off_diag_cxt} p
      <*> #{peek ls_hs_operator, apply_diag_cxt} p
      <*> #{peek ls_hs_operator, haskell_payload} p
  {-# INLINE peek #-}
  poke p x = do
    _ <- pokeRefCount (p `plusPtr` #{offset ls_hs_operator, refcount}) (coperator_refcount x)
    #{poke ls_hs_operator, basis}               p (coperator_basis x)
    #{poke ls_hs_operator, off_diag_terms}      p (coperator_off_diag_terms x)
    #{poke ls_hs_operator, diag_terms}          p (coperator_diag_terms x)
    #{poke ls_hs_operator, apply_off_diag_cxt}  p (coperator_apply_off_diag_cxt x)
    #{poke ls_hs_operator, apply_diag_cxt}      p (coperator_apply_diag_cxt x)
    #{poke ls_hs_operator, haskell_payload}     p (coperator_haskell_payload x)
  {-# INLINE poke #-}

data {-# CTYPE "lattice_symmetries_haskell.h" "chpl_external_array" #-} Cexternal_array = Cexternal_array
  { external_array_elts :: !(Ptr ()),
    external_array_num_elts :: !Word64,
    external_array_freer :: !(Ptr ())
  }

emptyExternalArray :: Cexternal_array
emptyExternalArray = Cexternal_array nullPtr 0 nullPtr

instance Storable Cexternal_array where
  sizeOf _ = #{size chpl_external_array}
  {-# INLINE sizeOf #-}
  alignment _ = #{alignment chpl_external_array}
  {-# INLINE alignment #-}
  peek p =
    Cexternal_array
      <$> #{peek chpl_external_array, elts} p
      <*> #{peek chpl_external_array, num_elts} p
      <*> #{peek chpl_external_array, freer} p
  {-# INLINE peek #-}
  poke p x = do
    #{poke chpl_external_array, elts} p (external_array_elts x)
    #{poke chpl_external_array, num_elts} p (external_array_num_elts x)
    #{poke chpl_external_array, freer} p (external_array_freer x)
  {-# INLINE poke #-}

-- foreign import ccall unsafe "&ls_hs_internal_destroy_external_array"
--   ls_hs_internal_destroy_external_array :: FunPtr (Ptr Cexternal_array -> IO ())

data {-# CTYPE "lattice_symmetries_haskell.h" "ls_hs_permutation_group" #-} Cpermutation_group = Cpermutation_group
  { cpermutation_group_refcount :: {-# UNPACK #-} !CInt,
    cpermutation_group_number_bits :: {-# UNPACK #-} !CInt,
    cpermutation_group_number_shifts :: {-# UNPACK #-} !CInt,
    cpermutation_group_number_masks :: {-# UNPACK #-} !CInt,
    cpermutation_group_masks :: {-# UNPACK #-} !(Ptr Word64),
    cpermutation_group_shifts :: {-# UNPACK #-} !(Ptr Word64),
    cpermutation_group_eigvals_re :: {-# UNPACK #-} !(Ptr CDouble),
    cpermutation_group_eigvals_im :: {-# UNPACK #-} !(Ptr CDouble)
  }

instance Storable Cpermutation_group where
  sizeOf _ = #{size ls_hs_permutation_group}
  {-# INLINE sizeOf #-}
  alignment _ = #{alignment ls_hs_permutation_group}
  {-# INLINE alignment #-}
  peek p =
    Cpermutation_group
      <$> peekRefCount (p `plusPtr` #{offset ls_hs_permutation_group, refcount})
      <*> #{peek ls_hs_permutation_group, number_bits} p
      <*> #{peek ls_hs_permutation_group, number_shifts} p
      <*> #{peek ls_hs_permutation_group, number_masks} p
      <*> #{peek ls_hs_permutation_group, masks} p
      <*> #{peek ls_hs_permutation_group, shifts} p
      <*> #{peek ls_hs_permutation_group, eigvals_re} p
      <*> #{peek ls_hs_permutation_group, eigvals_im} p
  {-# INLINE peek #-}
  poke p x = do
    _ <- pokeRefCount (p `plusPtr` #{offset ls_hs_permutation_group, refcount}) (cpermutation_group_refcount x)
    #{poke ls_hs_permutation_group, number_bits}   p (cpermutation_group_number_bits x)
    #{poke ls_hs_permutation_group, number_shifts} p (cpermutation_group_number_shifts x)
    #{poke ls_hs_permutation_group, number_masks}  p (cpermutation_group_number_masks x)
    #{poke ls_hs_permutation_group, masks}         p (cpermutation_group_masks x)
    #{poke ls_hs_permutation_group, shifts}        p (cpermutation_group_shifts x)
    #{poke ls_hs_permutation_group, eigvals_re}    p (cpermutation_group_eigvals_re x)
    #{poke ls_hs_permutation_group, eigvals_im}    p (cpermutation_group_eigvals_im x)
  {-# INLINE poke #-}

data {-# CTYPE "lattice_symmetries_haskell.h" "ls_chpl_kernels" #-} Cchpl_kernels = Cchpl_kernels
  { cchpl_kernels_enumerate_states :: FunPtr (Cbasis -> Word64 -> Word64 -> Ptr Cexternal_array -> IO ())
  }

instance Storable Cchpl_kernels where
  sizeOf _ = #{size ls_chpl_kernels}
  {-# INLINE sizeOf #-}
  alignment _ = #{alignment ls_chpl_kernels}
  {-# INLINE alignment #-}
  peek p =
    Cchpl_kernels
      <$> #{peek ls_chpl_kernels, enumerate_states} p
  {-# INLINE peek #-}
  poke p x = do
    #{poke ls_chpl_kernels, enumerate_states} p (cchpl_kernels_enumerate_states x)
  {-# INLINE poke #-}

data {-# CTYPE "lattice_symmetries_haskell.h" "ls_hs_yaml_config" #-} Cyaml_config = Cyaml_config
  { cyaml_config_basis :: {-# UNPACK #-} !(Ptr Cbasis),
    cyaml_config_hamiltonian :: {-# UNPACK #-} !(Ptr Coperator),
    cyaml_config_number_observables :: {-# UNPACK #-} !CInt,
    cyaml_config_observables :: {-# UNPACK #-} !(Ptr (Ptr Coperator))
  }

instance Storable Cyaml_config where
  sizeOf _ = #{size ls_hs_yaml_config}
  {-# INLINE sizeOf #-}
  alignment _ = #{alignment ls_hs_yaml_config}
  {-# INLINE alignment #-}
  peek p =
    Cyaml_config
      <$> #{peek ls_hs_yaml_config, basis} p
      <*> #{peek ls_hs_yaml_config, hamiltonian} p
      <*> #{peek ls_hs_yaml_config, number_observables} p
      <*> #{peek ls_hs_yaml_config, observables} p
  {-# INLINE peek #-}
  poke p x = do
    #{poke ls_hs_yaml_config, basis}               p (cyaml_config_basis x)
    #{poke ls_hs_yaml_config, hamiltonian}         p (cyaml_config_hamiltonian x)
    #{poke ls_hs_yaml_config, number_observables}  p (cyaml_config_number_observables x)
    #{poke ls_hs_yaml_config, observables}         p (cyaml_config_observables x)
  {-# INLINE poke #-}

-- foreign import ccall unsafe "ls_hs_internal_get_chpl_kernels"
--   ls_hs_internal_get_chpl_kernels :: IO (Ptr Cchpl_kernels)
