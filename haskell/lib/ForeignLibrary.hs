{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE OverloadedLists #-}
{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE TemplateHaskell #-}
{-# OPTIONS_GHC -Wno-missed-extra-shared-lib #-}

module ForeignLibrary () where

import Foreign.C.String (CString)
import Foreign.C.Types (CBool (..), CInt (..), CPtrdiff (..))
import HeaderFileGeneration
import LatticeSymmetries.Basis
import LatticeSymmetries.Context
import LatticeSymmetries.Lowering
import LatticeSymmetries.Operator
import LatticeSymmetries.Some
import LatticeSymmetries.Utils
import Prelude

$(pure [])

typesTable
  [ ([t|()|], "void")
  , ([t|CBool|], "bool")
  , ([t|CInt|], "int")
  , ([t|Int64|], "int64_t")
  , ([t|Word64|], "uint64_t")
  , ([t|CPtrdiff|], "ptrdiff_t")
  , ([t|Double|], "double")
  , ([t|CString|], "char const *")
  , ([t|Cbasis|], "ls_hs_basis")
  , ([t|Cexpr|], "ls_hs_expr")
  , ([t|Coperator|], "ls_hs_operator")
  , ([t|RawIsRepresentativeKernel|], "ls_hs_is_representative_kernel_type_v2")
  , ([t|RawStateInfoKernel|], "ls_hs_state_info_kernel_type_v2")
  , ([t|RawStateToIndexKernel|], "ls_hs_state_to_index_kernel_type")
  ]

headerFile "lattice_symmetries_functions.h"

addVerbatimPrefix
  [ "#include \"lattice_symmetries_types.h\""
  , "#include <stdint.h>"
  , ""
  , "#if defined(__cplusplus)"
  , "extern \"C\" {"
  , "#endif"
  , ""
  , "/* python-cffi: START */"
  ]

addVerbatimSuffix
  [ "/* python-cffi: STOP */"
  , ""
  , "#if defined(__cplusplus)"
  , "} // extern \"C\""
  , "#endif"
  ]

addDeclarations
  [ ------------------------
    -- Utilities
    "ls_hs_destroy_string"
  , ------------------------
    -- Basis
    "ls_hs_basis_from_json"
  , "ls_hs_basis_to_json"
  , "ls_hs_destroy_basis"
  , "ls_hs_init_basis_info"
  , "ls_hs_init_is_representative_kernel"
  , "ls_hs_init_state_info_kernel"
  , "ls_hs_init_state_to_index_kernel"
  , "ls_hs_fixed_hamming_state_to_index"
  , "ls_hs_fixed_hamming_index_to_state"
  , "ls_hs_basis_permutation_group"
  , ------------------------
    -- Expr
    "ls_hs_expr_to_json"
  , "ls_hs_expr_from_json"
  , "ls_hs_destroy_expr"
  , "ls_hs_expr_to_string"
  , "ls_hs_expr_plus"
  , "ls_hs_expr_minus"
  , "ls_hs_expr_times"
  , "ls_hs_expr_negate"
  , "ls_hs_expr_scale"
  , "ls_hs_expr_equal"
  , "ls_hs_replace_indices"
  , "ls_hs_expr_adjoint"
  , "ls_hs_expr_is_hermitian"
  , "ls_hs_expr_is_real"
  , "ls_hs_expr_is_identity"
  , "ls_hs_expr_spin_inversion_invariant"
  , "ls_hs_expr_conserves_number_particles"
  , "ls_hs_expr_particle_type"
  , "ls_hs_expr_number_sites"
  , "ls_hs_expr_permutation_group"
  , "ls_hs_expr_abelian_permutation_group"
  , "ls_hs_expr_hilbert_space_sectors"
  , "ls_hs_expr_ground_state_sectors"
  , ------------------------
    -- Operator
    "ls_hs_create_operator"
  , "ls_hs_destroy_operator"
  , "ls_hs_operator_from_json"
  , "ls_hs_operator_to_json"
  , "ls_hs_operator_max_number_off_diag"
  ]
