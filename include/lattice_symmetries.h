// c2chapel generates a warning for #pragma once that becomes an error ...
// Instead of messing with CFLAGS, let's just use the old school approach
#ifndef LATTICE_SYMMETRIES_H
#define LATTICE_SYMMETRIES_H

// #include <stdbool.h>
// #include <stddef.h>
#include <stdint.h>
#include <assert.h>
// #include <stdio.h>
// #include <stdlib.h>

typedef struct ls_numpy_array_1d { uint8_t* data; void* handle; } ls_numpy_array_1d;
typedef void (*ls_alloc_numpy_array_1d_callback)(uint64_t size, ls_numpy_array_1d* out);

static inline void ls_invoke_alloc_numpy_array_1d_callback(ls_alloc_numpy_array_1d_callback const callback, uint64_t const size, ls_numpy_array_1d* const out)
{
    assert(callback != NULL && "ls_invoke_alloc_numpy_array_1d_callback: 'callback' is NULL");
    assert(out != NULL && "ls_invoke_alloc_numpy_array_1d_callback: 'out' is NULL");
    (*callback)(size, out);
}

extern void chpl_library_init(int argc, char *argv[]);
extern void chpl_library_finalize(void);

#endif // LATTICE_SYMMETRIES_H


#if 0
#if defined(LS_NO_STD_COMPLEX)
/* python-cffi: START */
typedef struct ls_scalar {
    double _real;
    double _imag;
} ls_scalar;
/* python-cffi: STOP */
#else
typedef _Complex double ls_scalar;
#endif

/* python-cffi: START */
void ls_init(void);
void ls_exit(void);
/* python-cffi: STOP */

#define LS_HS_MAX_BLOCK_SIZE 16

// typedef void (*is_representative_kernel_type)(halide_buffer_t const *basis_states, halide_buffer_t *norms);
int ls_internal_invoke_is_representative_kernel(void const *kernel, int64_t count, uint64_t const *basis_states, uint16_t *norms);

// typedef void (*state_info_kernel_type)(halide_buffer_t const *basis_states, halide_buffer_t *representatives, halide_buffer_t *indices);
int ls_internal_invoke_state_info_kernel(void const *kernel, int64_t count, uint64_t const *basis_states, uint64_t *representatives, int32_t *indices);

// typedef void (*state_to_index_kernel_type)(halide_buffer_t const *basis_states, halide_buffer_t *indices);
int ls_internal_invoke_state_to_index_kernel(void const *kernel, int64_t count, uint64_t const *basis_states, int64_t *indices);

// void ls_hs_internal_axpy(int64_t size, double alpha_re, double alpha_im, ls_hs_scalar const *x, ls_hs_scalar const *y, ls_hs_scalar *out);

/* python-cffi: START */
typedef enum ls_particle_type {
    LS_HS_SPIN,
    LS_HS_SPINFUL_FERMION,
    LS_HS_SPINLESS_FERMION
} ls_particle_type;

typedef struct ls_basis_info {
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
    int number_characters;
    ls_scalar const *characters;
    // ... pointers to kernels ...
} ls_basis_info;

typedef struct ls_nonbranching_terms {
    int number_terms;
    int number_bits;
    // number_words = ceil(number_bits / 64)
    ls_scalar const *v;    // array of shape [number_terms]
    uint64_t const *m;     // array of shape [number_terms, number_words]
    uint64_t const *l;     // array of shape [number_terms, number_words]
    uint64_t const *r;     // array of shape [number_terms, number_words]
    uint64_t const *x;     // array of shape [number_terms, number_words]
    uint64_t const *s;     // array of shape [number_terms, number_words]
                           // all arrays are contiguous in row-major order
} ls_nonbranching_terms;
/* python-cffi: STOP */

typedef struct ls_operator ls_operator;

void ls_chpl_display_timings(void);


void* ls_chpl_local_enumerate_states(ls_basis_info* basisInfo, int64_t *numStates);
void ls_chpl_local_enumerate_states_complete(void* result, int64_t count, uint64_t *representatives);

void ls_chpl_matrix_vector_product_f64(ls_operator LS_CONST* matrix, int32_t numVectors, double LS_CONST* xPtr, double* yPtr);
void ls_chpl_matrix_vector_product_c128(ls_operator LS_CONST* matrix, int32_t numVectors, ls_hs_scalar LS_CONST* xPtr, ls_hs_scalar* yPtr);
#endif
