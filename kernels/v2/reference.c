#include "lattice_symmetries_types.h"
#include <assert.h>
#include <dlfcn.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// void ls_hs_fatal_error(char const *func, int const line, char const *message) {
//     fprintf(stderr, "[Error]   [%s#%i] %s\n[Error]   Aborting ...", func, line, message);
//     abort();
// }

// typedef void (*error_handler_type)(char const *);
// static _Atomic error_handler_type ls_hs_internal_error_handler;
//
// static void default_error_handler(char const *message) {
//     fprintf(stderr, "[Error]   %s\n[Error]   Aborting ...", message);
//     abort();
// }

// void ls_hs_set_exception_handler(error_handler_type handler) {
//     if (handler == NULL) {
//         handler = default_error_handler;
//     }
//     atomic_exchange(&ls_hs_internal_error_handler, handler);
// }

// void ls_hs_error(char const *message) {
//     error_handler_type handler = atomic_load(&ls_hs_internal_error_handler);
//     if (handler == NULL) {
//         default_error_handler(message);
//     } else {
//         (*handler)(message);
//     }
// }

typedef void (*ls_hs_internal_chpl_free_func)(void *);

void ls_hs_internal_destroy_external_array(chpl_external_array *arr) {
    // LS_CHECK(arr != NULL, "trying to destroy a NULL chpl_external_array");
    if (arr->freer != NULL) {
        ls_hs_internal_chpl_free_func const free_func = (ls_hs_internal_chpl_free_func)arr->freer;
        (*free_func)(arr->elts);
    }
}

void ls_hs_internal_object_init(ls_hs_object *object, int refcount, void *haskell_payload) {
    // LS_CHECK(object != NULL, "trying to init a NULL object");
    atomic_store(&object->refcount, refcount);
    object->haskell_payload = haskell_payload;
}
int ls_hs_internal_object_inc_ref_count(ls_hs_object *object) {
    // LS_CHECK(object != NULL, "trying to increase refcount of a NULL object");
    return atomic_fetch_add(&object->refcount, 1);
    // fprintf(stderr, "refcount: %d -> %d\n", r, r + 1);
}
int ls_hs_internal_object_dec_ref_count(ls_hs_object *object) {
    // LS_CHECK(object != NULL, "trying to decrease refcount of a NULL object");
    return atomic_fetch_sub(&object->refcount, 1);
    // fprintf(stderr, "refcount: %d -> %d\n", r, r - 1);
}

// TODO: currently not thread-safe, fix it
// NOTE: Kernels are initialized from within Chapel's module initialization functions, so we leave
// it up to Chapel to ensure these are run by just one task on each locale.
// static ls_chpl_kernels global_chpl_kernels = {.enumerate_states = NULL,
//                                               .operator_apply_off_diag = NULL,
//                                               .operator_apply_diag = NULL,
//                                               .matrix_vector_product_f64 = NULL,
//                                               .matrix_vector_product_c128 = NULL};
//
// ls_chpl_kernels const *ls_hs_internal_get_chpl_kernels(void) { return &global_chpl_kernels; }
//
// void ls_hs_internal_set_chpl_kernels(ls_chpl_kernels const *kernels) {
//     LS_CHECK(kernels != NULL, "trying to set Chapel kernels from NULL");
//     global_chpl_kernels = *kernels;
// }

// void ls_hs_state_index(ls_hs_basis const *const basis, ptrdiff_t const batch_size,
//                        uint64_t const *const restrict spins, ptrdiff_t const spins_stride,
//                        ptrdiff_t *const restrict indices, ptrdiff_t const indices_stride) {
//     LS_CHECK(basis->kernels->state_index_kernel != NULL, "state_index_kernel is NULL");
//     (*basis->kernels->state_index_kernel)(batch_size, spins, spins_stride, indices,
//     indices_stride,
//                                           basis->kernels->state_index_data);
// }
//
// void ls_hs_is_representative(ls_hs_basis const *basis, ptrdiff_t batch_size, uint64_t const
// *alphas,
//                              ptrdiff_t alphas_stride, uint8_t *are_representatives, double
//                              *norms) {
//     LS_CHECK(basis->kernels->is_representative_kernel != NULL,
//              "is_representative_kernel is NULL, perhaps this basis requires no "
//              "projection?");
//     (*basis->kernels->is_representative_kernel)(batch_size, alphas, alphas_stride,
//                                                 are_representatives, norms,
//                                                 basis->kernels->is_representative_data);
// }
//
// void ls_hs_state_info(ls_hs_basis const *basis, ptrdiff_t batch_size, uint64_t const *alphas,
//                       ptrdiff_t alphas_stride, uint64_t *betas, ptrdiff_t betas_stride,
//                       ls_hs_scalar *characters, double *norms) {
//     LS_CHECK(basis->kernels->state_info_kernel != NULL,
//              "state_info_kernel is NULL, perhaps this basis requires no "
//              "projection?");
//     (*basis->kernels->state_info_kernel)(batch_size, alphas, alphas_stride, betas, betas_stride,
//                                          characters, norms, basis->kernels->state_info_data);
// }

// void ls_hs_build_representatives(ls_hs_basis *basis, uint64_t const lower, uint64_t const upper)
// {
//     ls_chpl_kernels const *kernels = ls_hs_internal_get_chpl_kernels();
//     LS_CHECK(kernels->enumerate_states != NULL,
//              "enumerate_states kernel is NULL, Chapel code was supposed to "
//              "initialize it");
//     if (basis->representatives.num_elts > 0) {
//         // Early exit: representatives have already been built
//         return;
//     }
//     (*kernels->enumerate_states)(basis, lower, upper, &basis->representatives);
//
//     // if (basis->kernels->state_index_kernel == NULL) {
//     const int number_bits =
//         (basis->particle_type == LS_HS_SPINFUL_FERMION ? 2 : 1) * basis->number_sites;
//     const int default_cache_bits = 22;
//     basis->kernels->state_index_data = ls_hs_create_state_index_binary_search_kernel_data(
//         &basis->representatives, number_bits, default_cache_bits);
//     basis->kernels->state_index_kernel = &ls_hs_state_index_binary_search_kernel;
//     // }
// }
//
// void ls_hs_unchecked_set_representatives(ls_hs_basis *basis, chpl_external_array const *states,
//                                          int const cache_bits) {
//     LS_CHECK(basis->representatives.num_elts == 0, "representatives have already been set");
//     LS_CHECK(basis->kernels != NULL, "basis->kernels not set");
//     LS_CHECK(basis->kernels->state_index_kernel == NULL, "state_index_kernel has already been
//     set"); basis->representatives = *states; const int number_bits =
//         (basis->particle_type == LS_HS_SPINFUL_FERMION ? 2 : 1) * basis->number_sites;
//     basis->kernels->state_index_data = ls_hs_create_state_index_binary_search_kernel_data(
//         &basis->representatives, number_bits, cache_bits);
//     basis->kernels->state_index_kernel = &ls_hs_state_index_binary_search_kernel;
// }
