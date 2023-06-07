#pragma once

#if defined(__cplusplus)
#include <complex>
#include <cstddef>
#include <cstdint>
#else
#if defined(LS_C2CHAPEL) // c2chapel does not support C11 atomics
#define _Atomic
#else
#include <stdatomic.h>
#endif
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#endif

#if defined(LS_C2CHAPEL)
#include "chpl-external-array.h"
#else

#if defined __has_include
#if __has_include("chpl-external-array.h")
#include "chpl-external-array.h"
#else
/* python-cffi: START */
typedef struct {
  void *elts;
  uint64_t num_elts;

  void *freer;
} chpl_external_array;
/* python-cffi: STOP */
#endif
#endif

#endif

#if defined(__cplusplus)
extern "C" {
#endif

#if defined(__cplusplus)
using ls_hs_scalar = std::complex<double>;
#else
#if defined(LS_NO_STD_COMPLEX)
/* python-cffi: START */
typedef struct ls_hs_scalar {
  double _real;
  double _imag;
} ls_hs_scalar;
/* python-cffi: STOP */
#else
typedef _Complex double ls_hs_scalar;
#endif
#endif

typedef uint64_t ls_hs_bits;

/* python-cffi: START */
void ls_hs_init(void);
void ls_hs_exit(void);

void ls_hs_set_exception_handler(void (*handler)(char const *message));
void ls_hs_error(char const *message);
/* python-cffi: STOP */

__attribute__((noreturn)) void ls_hs_fatal_error(char const *func, int line,
                                                 char const *message);

#define LS_FATAL_ERROR(msg) ls_hs_fatal_error(__func__, __LINE__, msg)

#define LS_CHECK(cond, msg)                                                    \
  ((cond) ? ((void)0) : ls_hs_fatal_error(__func__, __LINE__, msg))

/* python-cffi: START */
typedef struct ls_hs_symmetry ls_hs_symmetry;

typedef struct ls_hs_symmetries ls_hs_symmetries;

typedef enum ls_hs_particle_type {
  LS_HS_SPIN,
  LS_HS_SPINFUL_FERMION,
  LS_HS_SPINLESS_FERMION
} ls_hs_particle_type;

typedef void (*ls_hs_internal_state_index_kernel_type)(
    ptrdiff_t batch_size, uint64_t const *alphas, ptrdiff_t alphas_stride,
    ptrdiff_t *indices, ptrdiff_t indices_stride, void const *private_data);

typedef void (*ls_hs_internal_is_representative_kernel_type)(
    ptrdiff_t batch_size, uint64_t const *alphas, ptrdiff_t alphas_stride,
    uint8_t *are_representatives, double *norms, void const *private_data);

typedef void (*ls_hs_internal_state_info_kernel_type)(
    ptrdiff_t batch_size, uint64_t const *alphas, ptrdiff_t alphas_stride,
    uint64_t *betas, ptrdiff_t betas_stride, ls_hs_scalar *characters,
    double *norms, void const *private_data);

typedef struct ls_hs_basis_kernels {
  ls_hs_internal_state_info_kernel_type state_info_kernel;
  void *state_info_data;
  ls_hs_internal_is_representative_kernel_type is_representative_kernel;
  void *is_representative_data;
  ls_hs_internal_state_index_kernel_type state_index_kernel;
  void *state_index_data;
} ls_hs_basis_kernels;

typedef struct ls_hs_permutation_group {
  _Atomic int refcount;
  int number_bits;
  int number_shifts;
  int number_masks;
  uint64_t *masks;
  uint64_t *shifts;
  double *eigvals_re;
  double *eigvals_im;
  void *haskell_payload;
} ls_hs_permutation_group;

typedef struct ls_hs_basis {
  _Atomic int refcount;
  int number_sites;
  int number_particles;
  int number_up;
  ls_hs_particle_type particle_type;
  int spin_inversion;
  bool state_index_is_identity;
  bool requires_projection;
  ls_hs_basis_kernels *kernels;
  chpl_external_array representatives;
  void *haskell_payload;
} ls_hs_basis;

typedef struct ls_hs_expr ls_hs_expr;

typedef void (*ls_hs_index_replacement_type)(int spin, int site, int *new_spin,
                                             int *new_site);

typedef struct ls_hs_nonbranching_terms {
  int number_terms;
  int number_bits;
  // number_words = ceil(number_bits / 64)
  ls_hs_scalar const *v; // array of shape [number_terms]
  uint64_t const *m;     // array of shape [number_terms, number_words]
  uint64_t const *l;     // array of shape [number_terms, number_words]
  uint64_t const *r;     // array of shape [number_terms, number_words]
  uint64_t const *x;     // array of shape [number_terms, number_words]
  uint64_t const *s;     // array of shape [number_terms, number_words]
  // all arrays are contiguous in row-major order
} ls_hs_nonbranching_terms;

typedef struct ls_hs_operator {
  _Atomic int refcount;
  ls_hs_basis const *basis;
  ls_hs_nonbranching_terms const *off_diag_terms;
  ls_hs_nonbranching_terms const *diag_terms;
  // ls_internal_operator_kernel_data const *apply_off_diag_cxt;
  // ls_internal_operator_kernel_data const *apply_diag_cxt;
  void *haskell_payload;
} ls_hs_operator;

typedef struct ls_hs_yaml_config {
  ls_hs_basis const *basis;
  ls_hs_operator const *hamiltonian;
  int number_observables;
  ls_hs_operator const *const *observables;
} ls_hs_yaml_config;

typedef struct ls_chpl_kernels {
  void (*enumerate_states)(ls_hs_basis const *, uint64_t, uint64_t,
                           chpl_external_array *);
  void (*operator_apply_off_diag)(ls_hs_operator *, int64_t, uint64_t *,
                                  chpl_external_array *, chpl_external_array *,
                                  chpl_external_array *, int64_t);
  void (*operator_apply_diag)(ls_hs_operator *, int64_t, uint64_t *,
                              chpl_external_array *, int64_t);
  void (*matrix_vector_product)(ls_hs_operator *, int, double const *,
                                double *);
} ls_chpl_kernels;
/* python-cffi: STOP */

// Hide these functions from c2chapel, because they are implemented on the
// Chapel side
#if !defined(LS_C2CHAPEL)
/* python-cffi: START */
void ls_chpl_init(void);
void ls_chpl_finalize(void);
/* python-cffi: STOP */
#endif

/* python-cffi: START */
ls_chpl_kernels const *ls_hs_internal_get_chpl_kernels(void);
/* python-cffi: STOP */
void ls_hs_internal_set_chpl_kernels(ls_chpl_kernels const *kernels);

void ls_hs_internal_destroy_external_array(chpl_external_array *arr);
/* python-cffi: STOP */

// TODO: remove all of this
typedef struct ls_hs_state_index_binary_search_data
    ls_hs_state_index_binary_search_data;

ls_hs_state_index_binary_search_data *
ls_hs_create_state_index_binary_search_kernel_data(
    chpl_external_array const *representatives, int number_bits);

void ls_hs_destroy_state_index_binary_search_kernel_data(
    ls_hs_state_index_binary_search_data *cache);

void ls_hs_state_index_binary_search_kernel(ptrdiff_t batch_size,
                                            uint64_t const *spins,
                                            ptrdiff_t spins_stride,
                                            ptrdiff_t *indices,
                                            ptrdiff_t indices_stride,
                                            void const *private_kernel_data);

/* python-cffi: START */
// The following should, in principle, not exist
void ls_hs_state_index(ls_hs_basis const *const basis,
                       ptrdiff_t const batch_size,
                       uint64_t const *const restrict spins,
                       ptrdiff_t const spins_stride,
                       ptrdiff_t *const restrict indices,
                       ptrdiff_t const indices_stride);
void ls_hs_is_representative(ls_hs_basis const *basis, ptrdiff_t batch_size,
                             uint64_t const *alphas, ptrdiff_t alphas_stride,
                             uint8_t *are_representatives, double *norms);
void ls_hs_state_info(ls_hs_basis const *basis, ptrdiff_t batch_size,
                      uint64_t const *alphas, ptrdiff_t alphas_stride,
                      uint64_t *betas, ptrdiff_t betas_stride,
                      ls_hs_scalar *characters, double *norms);
void ls_hs_build_representatives(ls_hs_basis *basis, uint64_t const lower,
                                 uint64_t const upper);
void ls_hs_unchecked_set_representatives(ls_hs_basis *basis,
                                         chpl_external_array const *states);
/* python-cffi: STOP */

#ifdef __cplusplus
} // extern "C"
#endif

