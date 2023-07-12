#include "lattice_symmetries_types.h"
#include <assert.h>
#include <complex.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void ls_hs_fatal_error(char const *func, int const line, char const *message) {
  fprintf(stderr, "[Error]   [%s#%i] %s\n[Error]   Aborting ...", func, line,
          message);
  abort();
}

typedef void (*error_handler_type)(char const *);
static _Atomic error_handler_type ls_hs_internal_error_handler;

static void default_error_handler(char const *message) {
  fprintf(stderr, "[Error]   %s\n[Error]   Aborting ...", message);
  abort();
}

void ls_hs_set_exception_handler(error_handler_type handler) {
  if (handler == NULL) {
    handler = default_error_handler;
  }
  atomic_exchange(&ls_hs_internal_error_handler, handler);
}

void ls_hs_error(char const *message) {
  error_handler_type handler = atomic_load(&ls_hs_internal_error_handler);
  LS_CHECK(handler != NULL, "error handler is NULL");
  (*handler)(message);
}

typedef void (*ls_hs_internal_chpl_free_func)(void *);

void ls_hs_internal_destroy_external_array(chpl_external_array *arr) {
  // fprintf(stderr, "Finalizing chpl_external_array %p ...\n", arr);
  LS_CHECK(arr != NULL, "trying to destroy a NULL chpl_external_array");
  if (arr->freer != NULL) {
    ls_hs_internal_chpl_free_func const free_func =
        (ls_hs_internal_chpl_free_func)arr->freer;
    (*free_func)(arr->elts);
  }
}

int ls_hs_internal_read_refcount(_Atomic int const *refcount) {
  return atomic_load(refcount);
}

void ls_hs_internal_write_refcount(_Atomic int *refcount, int value) {
  atomic_store(refcount, value);
}

int ls_hs_internal_inc_refcount(_Atomic int *refcount) {
  return atomic_fetch_add(refcount, 1);
}

int ls_hs_internal_dec_refcount(_Atomic int *refcount) {
  return atomic_fetch_sub(refcount, 1);
}

#if 1
void ls_internal_operator_apply_diag_x1(ls_hs_operator const *const op,
                                        ptrdiff_t const batch_size,
                                        uint64_t const *restrict const alphas,
                                        double *restrict const ys,
                                        double const *restrict const xs) {
  // The diagonal is zero
  if (op->diag_terms == NULL || op->diag_terms->number_terms == 0) {
    memset(ys, 0, (size_t)batch_size * sizeof(double));
    return;
  }

  ls_hs_nonbranching_terms const *restrict terms = op->diag_terms;
  ptrdiff_t const number_terms = terms->number_terms;
  for (ptrdiff_t batch_idx = 0; batch_idx < batch_size; ++batch_idx) {
    double acc = 0;
    uint64_t const alpha = alphas[batch_idx];
    for (ptrdiff_t term_idx = 0; term_idx < number_terms; ++term_idx) {
      uint8_t const delta = (alpha & terms->m[term_idx]) == terms->r[term_idx];
      if (delta != 0) {
        int const sign =
            1 - 2 * (__builtin_popcountll(alpha & terms->s[term_idx]) % 2);
        double const factor =
            (xs != NULL) ? delta * sign * xs[batch_idx] : delta * sign;
        acc += creal(terms->v[term_idx]) * factor;
      }
    }
    ys[batch_idx] = acc;
  }
}

void ls_internal_operator_apply_off_diag_x1(ls_hs_operator const *const op,
                                            ptrdiff_t const batch_size,
                                            uint64_t const *const alphas,
                                            uint64_t *const betas,
                                            ls_hs_scalar *const coeffs,
                                            ptrdiff_t *const offsets,
                                            double const *const xs) {
  // Nothing to apply
  if (op->off_diag_terms == NULL || op->off_diag_terms->number_terms == 0) {
    memset(offsets, 0, (size_t)(batch_size + 1) * sizeof(ptrdiff_t));
    return;
  }

  ls_hs_nonbranching_terms const *restrict terms = op->off_diag_terms;
  ptrdiff_t const number_terms = terms->number_terms;

  offsets[0] = 0;
  ptrdiff_t offset = 0;
  for (ptrdiff_t batch_idx = 0; batch_idx < batch_size; ++batch_idx) {
    for (ptrdiff_t term_idx = 0; term_idx < number_terms; ++term_idx) {
      uint8_t const delta =
          (alphas[batch_idx] & terms->m[term_idx]) == terms->r[term_idx];
      if (delta != 0) {
        int const sign = 1 - 2 * (__builtin_popcountll(alphas[batch_idx] &
                                                       terms->s[term_idx]) %
                                  2);
        double const factor =
            (xs != NULL) ? delta * sign * xs[batch_idx] : delta * sign;
        coeffs[offset] = CMPLX(creal(terms->v[term_idx]) * factor,
                               cimag(terms->v[term_idx]) * factor);
        betas[offset] = alphas[batch_idx] ^ terms->x[term_idx];
        ++offset;
      }
    }
    offsets[batch_idx + 1] = offset;
  }
  return;
}
#endif

void ls_hs_state_index(ls_hs_basis const *const basis,
                       ptrdiff_t const batch_size,
                       uint64_t const *const restrict spins,
                       ptrdiff_t const spins_stride,
                       ptrdiff_t *const restrict indices,
                       ptrdiff_t const indices_stride) {
  LS_CHECK(basis->kernels->state_index_kernel != NULL,
           "state_index_kernel is NULL");
  (*basis->kernels->state_index_kernel)(batch_size, spins, spins_stride,
                                        indices, indices_stride,
                                        basis->kernels->state_index_data);
}

void ls_hs_is_representative(ls_hs_basis const *basis, ptrdiff_t batch_size,
                             uint64_t const *alphas, ptrdiff_t alphas_stride,
                             uint8_t *are_representatives, double *norms) {
  LS_CHECK(basis->kernels->is_representative_kernel != NULL,
           "is_representative_kernel is NULL, perhaps this basis requires no "
           "projection?");
  (*basis->kernels->is_representative_kernel)(
      batch_size, alphas, alphas_stride, are_representatives, norms,
      basis->kernels->is_representative_data);
}

void ls_hs_state_info(ls_hs_basis const *basis, ptrdiff_t batch_size,
                      uint64_t const *alphas, ptrdiff_t alphas_stride,
                      uint64_t *betas, ptrdiff_t betas_stride,
                      ls_hs_scalar *characters, double *norms) {
  LS_CHECK(basis->kernels->state_info_kernel != NULL,
           "state_info_kernel is NULL, perhaps this basis requires no "
           "projection?");
  (*basis->kernels->state_info_kernel)(batch_size, alphas, alphas_stride, betas,
                                       betas_stride, characters, norms,
                                       basis->kernels->state_info_data);
}

void ls_hs_build_representatives(ls_hs_basis *basis, uint64_t const lower,
                                 uint64_t const upper) {
  ls_chpl_kernels const *kernels = ls_hs_internal_get_chpl_kernels();
  LS_CHECK(kernels->enumerate_states != NULL,
           "enumerate_states kernel is NULL, Chapel code was supposed to "
           "initialize it");
  if (basis->representatives.num_elts > 0) {
    // Early exit: representatives have already been built
    return;
  }
  (*kernels->enumerate_states)(basis, lower, upper, &basis->representatives);

  // if (basis->kernels->state_index_kernel == NULL) {
  const int number_bits =
    (basis->particle_type == LS_HS_SPINFUL_FERMION ? 2 : 1) * basis->number_sites;
  const int default_cache_bits = 22;
  basis->kernels->state_index_data =
      ls_hs_create_state_index_binary_search_kernel_data(
          &basis->representatives, number_bits, default_cache_bits);
  basis->kernels->state_index_kernel = &ls_hs_state_index_binary_search_kernel;
  // }
}

void ls_hs_unchecked_set_representatives(ls_hs_basis *basis,
                                         chpl_external_array const *states,
					 int const cache_bits) {
  LS_CHECK(basis->representatives.num_elts == 0,
           "representatives have already been set");
  LS_CHECK(basis->kernels != NULL, "basis->kernels not set");
  LS_CHECK(basis->kernels->state_index_kernel == NULL,
           "state_index_kernel has already been set");
  basis->representatives = *states;
  const int number_bits =
    (basis->particle_type == LS_HS_SPINFUL_FERMION ? 2 : 1) * basis->number_sites;
  basis->kernels->state_index_data =
      ls_hs_create_state_index_binary_search_kernel_data(
          &basis->representatives, number_bits, cache_bits);
  basis->kernels->state_index_kernel = &ls_hs_state_index_binary_search_kernel;
}

// TODO: currently not thread-safe, fix it
static ls_chpl_kernels global_chpl_kernels = {.enumerate_states = NULL,
                                              .operator_apply_off_diag = NULL,
                                              .operator_apply_diag = NULL,
                                              .matrix_vector_product = NULL};
// static pthread_mutex_t global_chpl_kernels_mutex = PTHREAD_MUTEX_INITIALIZER;

ls_chpl_kernels const *ls_hs_internal_get_chpl_kernels() {
  return &global_chpl_kernels;
}
void ls_hs_internal_set_chpl_kernels(ls_chpl_kernels const *kernels) {
  global_chpl_kernels = *kernels;
}
