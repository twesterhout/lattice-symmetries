#pragma once

#if defined(__cplusplus)
#include <atomic>
#include <complex>
#include <cstddef>
#include <cstdint>
#else
#if defined(LS_C2CHAPEL) || defined(PYTHON_CFFI) // c2chapel and cffi do not support C11 atomics
#define _Atomic
#else
#include <stdatomic.h>
#endif
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#endif

#if defined(__cplusplus)
#define LS_HS_ATOMIC(type) std::atomic<type>
#elif defined(LS_C2CHAPEL) || defined(PYTHON_CFFI)
#define LS_HS_ATOMIC(type) type
#else
#define LS_HS_ATOMIC(type) type _Atomic
#endif

/* python-cffi: START */
typedef struct halide_buffer_t halide_buffer_t;
/* python-cffi: STOP */

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

/* python-cffi: START */
void ls_hs_init(void);
void ls_hs_exit(void);
/* python-cffi: STOP */

#if defined(PYTHON_CFFI)
// This is a fake definiion specifically for Python's cffi
/* python-cffi: START */
typedef struct ls_hs_object {
    int64_t refcount;
    void *haskell_payload;
} ls_hs_object;
/* python-cffi: STOP */
#else
typedef struct ls_hs_object {
    LS_HS_ATOMIC(int64_t) refcount;
    void *haskell_payload;
} ls_hs_object;
#endif

/* python-cffi: START */
typedef struct ls_hs_permutation ls_hs_permutation;
typedef struct ls_hs_permutation_group ls_hs_permutation_group;
typedef struct ls_hs_rep_element ls_hs_rep_element;
typedef struct ls_hs_representation ls_hs_representation;

typedef void (*ls_hs_is_representative_kernel_type_v2)(halide_buffer_t const *basis_states,
                                                       halide_buffer_t *norms);

typedef void (*ls_hs_state_info_kernel_type_v2)(halide_buffer_t const *basis_states,
                                                halide_buffer_t *representatives,
                                                halide_buffer_t *indices);

typedef void (*ls_hs_state_to_index_kernel_type)(halide_buffer_t const *basis_states,
                                                 halide_buffer_t *indices);

typedef enum ls_hs_particle_type {
    LS_HS_SPIN,
    LS_HS_SPINFUL_FERMION,
    LS_HS_SPINLESS_FERMION
} ls_hs_particle_type;

typedef struct ls_hs_basis_info {
    bool has_permutation_symmetries;
    bool requires_projection;
    bool is_state_index_identity;
    bool is_real;
    int number_bits;
    int number_words;
    int number_sites;
    int number_particles;
    int number_up;
    int hamming_weight;
    int spin_inversion;
    uint64_t min_state_estimate;
    uint64_t max_state_estimate;
    ls_hs_particle_type particle_type;
} ls_hs_basis_info;

typedef struct ls_hs_basis {
    // Reference count and Haskell payload
    ls_hs_object base;
    // Kernels for working with symmetry-adapted bases. Initialized lazily.
    LS_HS_ATOMIC(ls_hs_is_representative_kernel_type_v2) is_representative_kernel;
    LS_HS_ATOMIC(ls_hs_state_info_kernel_type_v2) state_info_kernel;
    LS_HS_ATOMIC(ls_hs_state_to_index_kernel_type) state_to_index_kernel;
    // Pre-computed representatives & norms
    chpl_external_array local_representatives; // uint64_t [N]
    chpl_external_array local_norms;           // uint16_t [N]
    // Auxiliary information. Initialized lazily upon request. By using atomic, we can ensure that
    // multiple threads don't write to it at the same time.
    LS_HS_ATOMIC(ls_hs_basis_info *) info;
} ls_hs_basis;

typedef struct ls_hs_expr {
    ls_hs_object base;
} ls_hs_expr;

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
    ls_hs_object base;
    ls_hs_basis *basis;
    ls_hs_expr *expr;
    ls_hs_nonbranching_terms *diag_terms;
    ls_hs_nonbranching_terms *off_diag_terms;
    int max_number_off_diag;
    int max_number_off_diag_estimate;
} ls_hs_operator;

typedef struct ls_chpl_batched_operator {
    ls_hs_operator const *matrix;
    int64_t batch_size;
    uint64_t *betas;
    ls_hs_scalar *coeffs;
    uint64_t *target_states;
    int64_t *offsets;
    uint64_t *temp_spins;
    ls_hs_scalar *temp_coeffs;
    double *temp_norms;
    uint8_t *locale_indices;
} ls_chpl_batched_operator;
/* python-cffi: STOP */

// Hide these functions from c2chapel, because they are implemented on the Chapel side
#if !defined(LS_C2CHAPEL)
void ls_chpl_init(void);
void ls_chpl_finalize(void);
#endif

// typedef ls_chpl_kernels ls_chpl_kernels;
// typedef struct ls_chpl_kernels {
//     void (*apply_diag_kernel_f64)(ls_hs_nonbranching_terms const *, int64_t, uint64_t const *,
//                                   chpl_external_array *, chpl_external_array *);
//     void (*apply_diag_kernel_c128)(ls_hs_nonbranching_terms const *, int64_t, uint64_t const *,
//                                    chpl_external_array *, chpl_external_array *);
//     void (*apply_off_diag_kernel)(ls_hs_nonbranching_terms const *, int64_t, uint64_t const *,
//                                   chpl_external_array *, chpl_external_array *);
//     void (*enumerate_states)(ls_hs_basis const *, uint64_t, uint64_t, chpl_external_array *);
//     // void (*operator_apply_off_diag)(ls_hs_operator *, int64_t, uint64_t *, chpl_external_array
//     *,
//     //                                 chpl_external_array *, chpl_external_array *, int64_t);
//     // void (*operator_apply_diag)(ls_hs_operator *, int64_t, uint64_t *, chpl_external_array *,
//     //                             int64_t);
//     // void (*matrix_vector_product_f64)(ls_hs_operator *, int, double const *, double *);
//     // void (*matrix_vector_product_c128)(ls_hs_operator *, int, ls_hs_scalar const *,
//     ls_hs_scalar
//     // *); void (*operator_to_csr)(ls_hs_operator *, chpl_external_array *, chpl_external_array
//     *,
//     //                         chpl_external_array *, int64_t);
//     void (*matrix_vector_product_csr_i32_c128)(int64_t, int64_t, int64_t, ls_hs_scalar const *,
//                                                int32_t const *, int32_t const *,
//                                                ls_hs_scalar const *, ls_hs_scalar *, int64_t);
// } ls_chpl_kernels;

// Implemented directly in C to avoid the overhead of going through Haskell's runtime.
ls_hs_basis_info const *ls_hs_get_basis_info(ls_hs_basis const *);

// ls_chpl_kernels const *ls_hs_internal_get_chpl_kernels(void);
// WARNING: ls_hs_internal_set_chpl_kernels is not thread-safe
// void ls_hs_internal_set_chpl_kernels(ls_chpl_kernels const *kernels);
void ls_hs_internal_destroy_external_array(chpl_external_array *arr);

void ls_hs_internal_object_init(ls_hs_object *, int refcount, void *haskell_payload);
int ls_hs_internal_object_inc_ref_count(ls_hs_object *);
int ls_hs_internal_object_dec_ref_count(ls_hs_object *);

ls_hs_is_representative_kernel_type_v2
ls_hs_internal_mk_is_representative_kernel(halide_buffer_t const *masks,
                                           halide_buffer_t const *eigvals_re,
                                           halide_buffer_t const *shifts, int spin_inversion);

ls_hs_state_info_kernel_type_v2 ls_hs_internal_mk_state_info_kernel(halide_buffer_t const *masks,
                                                                    halide_buffer_t const *shifts,
                                                                    int spin_inversion);

ls_hs_state_to_index_kernel_type
mk_fixed_hamming_state_to_index_kernel(int number_sites, int hamming_weight,
                                       halide_buffer_t const *binomials);

#ifdef __cplusplus
} // extern "C"
#endif
