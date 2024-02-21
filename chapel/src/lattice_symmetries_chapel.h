#pragma once

#include <stdint.h>
#include <stddef.h>
#include <lattice_symmetries_types.h>

/* python-cffi: START */
void ls_chpl_init(void);
void ls_chpl_finalize(void);
void ls_chpl_display_timings(void);

ls_hs_basis_info const *ls_chpl_get_basis_info(ls_hs_basis const *basis);

ls_hs_is_representative_kernel_type_v2 ls_chpl_get_is_representative_kernel(ls_hs_basis const *basis);
void ls_chpl_invoke_is_representative_kernel(ls_hs_is_representative_kernel_type_v2 kernel,
                                             int64_t count, uint64_t const *basis_states,
                                             double *norms);

ls_hs_state_to_index_kernel_type ls_chpl_get_state_to_index_kernel(ls_hs_basis const *basis);
void ls_chpl_invoke_state_to_index_kernel(ls_hs_state_to_index_kernel_type kernel,
                                          int64_t count, uint64_t const *basis_states,
                                          int64_t *indices);

ls_hs_state_info_kernel_type_v2 ls_chpl_get_state_info_kernel(ls_hs_basis const *basis);
void ls_chpl_invoke_state_info_kernel(ls_hs_state_info_kernel_type_v2 kernel,
                                      int64_t count, uint64_t const *basis_states,
                                      uint64_t *representatives, int32_t *indices);

void ls_chpl_local_enumerate_states(ls_hs_basis* p, uint64_t lower, uint64_t upper);

void ls_chpl_matrix_vector_product_f64(ls_hs_operator const* matrix, int numVectors,
                                       double const* xPtr, double* yPtr);

void ls_chpl_matrix_vector_product_c128(ls_hs_operator const* matrix, int numVectors,
                                        ls_hs_scalar const* xPtr, ls_hs_scalar* yPtr);

void ls_chpl_off_diag_to_csr_c128(ls_hs_operator const* matrix,
                                  chpl_external_array* matrixElements,
                                  chpl_external_array* rowOffsets,
                                  chpl_external_array* colIndices);

void ls_chpl_extract_diag_c128(ls_hs_operator const* matrix, chpl_external_array* diag);
/* python-cffi: STOP */
