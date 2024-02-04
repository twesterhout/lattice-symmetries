#pragma once

#include <lattice_symmetries_types.h>

ls_hs_basis_info const *ls_chpl_get_basis_info(ls_hs_basis const *basis);
ls_hs_is_representative_kernel_type_v2
ls_chpl_get_is_representative_kernel(ls_hs_basis const *basis);

void ls_chpl_invoke_is_representative_kernel(ls_hs_is_representative_kernel_type_v2 kernel,
                                             int64_t count, uint64_t const *basis_states,
                                             double *norms);
