#pragma once

#include <stdint.h>
#include <stddef.h>
#include <lattice_symmetries_types.h>

#ifdef CHPL_C_BACKEND
#define LS_CONST
#else
#define LS_CONST const
#endif

/* python-cffi: START */
void ls_chpl_init(void);
void ls_chpl_finalize(void);
void ls_chpl_display_timings(void);

ls_hs_basis_info const *ls_chpl_get_basis_info(ls_hs_basis const *basis);

ls_hs_is_representative_kernel_type_v2 ls_chpl_get_is_representative_kernel(ls_hs_basis const *basis);
void ls_chpl_invoke_is_representative_kernel(ls_hs_is_representative_kernel_type_v2 kernel,
                                             int64_t count, uint64_t const *basis_states,
                                             uint16_t *norms);

ls_hs_state_to_index_kernel_type ls_chpl_get_state_to_index_kernel(ls_hs_basis const *basis);
void ls_chpl_invoke_state_to_index_kernel(ls_hs_state_to_index_kernel_type kernel,
                                          int64_t count, uint64_t const *basis_states,
                                          int64_t *indices);

ls_hs_state_info_kernel_type_v2 ls_chpl_get_state_info_kernel(ls_hs_basis const *basis);
void ls_chpl_invoke_state_info_kernel(ls_hs_state_info_kernel_type_v2 kernel,
                                      int64_t count, uint64_t const *basis_states,
                                      uint64_t *representatives, int32_t *indices);

void ls_chpl_local_enumerate_states(ls_hs_basis* p, uint64_t lower, uint64_t upper);

void ls_chpl_matrix_vector_product_f64(ls_hs_operator LS_CONST* matrix, int32_t numVectors,
                                       double LS_CONST* xPtr, double* yPtr);

void ls_chpl_matrix_vector_product_c128(ls_hs_operator LS_CONST* matrix, int32_t numVectors,
                                        ls_hs_scalar LS_CONST* xPtr, ls_hs_scalar* yPtr);

void ls_chpl_off_diag_to_csr_c128(ls_hs_operator LS_CONST* matrix,
                                  chpl_external_array* matrixElements,
                                  chpl_external_array* rowOffsets,
                                  chpl_external_array* colIndices);

void ls_chpl_extract_diag_c128(ls_hs_operator LS_CONST* matrix, chpl_external_array* diag);

void ls_chpl_experimental_axpy_c128(int64_t size, double alpha_re, double alpha_im,
                                    ls_hs_scalar LS_CONST* xs, ls_hs_scalar LS_CONST* ys,
                                    ls_hs_scalar* out);

void ls_chpl_matrix_vector_product_csr_i32_f64(int64_t numberRows, int64_t numberCols, int64_t numberNonZero,
                                               double LS_CONST* matrixElements, int32_t LS_CONST* rowOffsets, int32_t LS_CONST* colIndices,
                                               double LS_CONST* x, double* y, int64_t numTasks);
void ls_chpl_matrix_vector_product_csr_i32_c128(int64_t numberRows, int64_t numberCols, int64_t numberNonZero,
                                                ls_hs_scalar LS_CONST* matrixElements, int32_t LS_CONST* rowOffsets, int32_t LS_CONST* colIndices,
                                                ls_hs_scalar LS_CONST* x, ls_hs_scalar* y, int64_t numTasks);
void ls_chpl_matrix_vector_product_csr_i64_f64(int64_t numberRows, int64_t numberCols, int64_t numberNonZero,
                                               double LS_CONST* matrixElements, int64_t LS_CONST* rowOffsets, int64_t LS_CONST* colIndices,
                                               double LS_CONST* x, double* y, int64_t numTasks);
void ls_chpl_matrix_vector_product_csr_i64_c128(int64_t numberRows, int64_t numberCols, int64_t numberNonZero,
                                                ls_hs_scalar LS_CONST* matrixElements, int64_t LS_CONST* rowOffsets, int64_t LS_CONST* colIndices,
                                                ls_hs_scalar LS_CONST* x, ls_hs_scalar* y, int64_t numTasks);
/* python-cffi: STOP */
