#include "lattice_symmetries_haskell.h"
#include <assert.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

typedef uint64_t ls_hs_bits;

typedef struct ls_hs_nonbranching_term {
  ls_hs_scalar v;
  ls_hs_bits m;
  ls_hs_bits l;
  ls_hs_bits r;
  ls_hs_bits x;
  ls_hs_bits s;
} ls_hs_nonbranching_term;

static inline int popcount(ls_hs_bits x) { return __builtin_popcountl(x); }

// Compute term|α⟩ = coeff|β⟩
void ls_hs_apply_nonbranching_term(ls_hs_nonbranching_term const *const term,
                                   ls_hs_bits const alpha,
                                   ls_hs_bits *const beta,
                                   ls_hs_scalar *const coeff) {
  int const delta = (alpha & term->m) == term->r;
  int const sign = 1 - 2 * (popcount(alpha & term->s) % 2);
  *beta = alpha ^ term->x;
  *coeff = term->v * (delta * sign);
}

ls_hs_binomials *ls_hs_internal_malloc_binomials(int const number_bits) {
  ls_hs_binomials *p = malloc(sizeof(ls_hs_binomials));
  if (p == NULL) {
    goto fail_1;
  }
  p->dimension = number_bits + 1; // NOTE: +1 because we could have n and k
                                  // running from 0 to number_bits inclusive
  p->coefficients = malloc(p->dimension * p->dimension * sizeof(uint64_t));
  if (p->coefficients == NULL) {
    goto fail_2;
  }
  return p;

fail_2:
  free(p);
fail_1:
  return NULL;
}

void ls_hs_internal_free_binomials(ls_hs_binomials *p) {
  if (p != NULL) {
    free(p->coefficients);
  }
  free(p);
}

void ls_hs_internal_compute_binomials(ls_hs_binomials *p) {
  int const dim = p->dimension;
  uint64_t *const coeff = p->coefficients;
  int n = 0;
  int k = 0;
  coeff[n * dim + k] = 1;
  for (int k = 1; k < dim; ++k) {
    coeff[n * dim + k] = 0;
  }
  for (n = 1; n < dim; ++n) {
    coeff[n * dim + 0] = 1;
    for (k = 1; k <= n; ++k) {
      coeff[n * dim + k] =
          coeff[(n - 1) * dim + (k - 1)] + coeff[(n - 1) * dim + k];
    }
    for (; k < n; ++k) {
      coeff[n * dim + k] = 0;
    }
  }
}

uint64_t ls_hs_internal_binomial(int const n, int const k,
                                 ls_hs_binomials const *cache) {
  if (k > n) {
    return 0;
  }
  assert(0 <= n && n < cache->dimension);
  assert(0 <= k && k < cache->dimension);
  return cache->coefficients[n * cache->dimension + k];
}

static uint64_t rank_via_combinadics(uint64_t alpha,
                                     ls_hs_binomials const *cache) {
  uint64_t i = 0;
  for (int k = 1; alpha != 0; ++k) {
    int c = __builtin_ctzl(alpha);
    alpha &= alpha - 1;
    i += ls_hs_internal_binomial(c, k, cache);
  }
  return i;
}

// -- |
// -- Cindex_kernel
// --   batch_size
// --   spins
// --   spins_stride
// --   indices
// --   indices_stride
// --   private_kernel_data
// type Cindex_kernel = CPtrdiff -> Ptr Word64 -> CPtrdiff -> Ptr CPtrdiff ->
// CPtrdiff -> Ptr () -> IO ()
void ls_hs_state_index_combinadics_kernel(ptrdiff_t const batch_size,
                                          uint64_t const *spins,
                                          ptrdiff_t const spins_stride,
                                          ptrdiff_t *const restrict indices,
                                          ptrdiff_t const indices_stride,
                                          void const *private_kernel_data) {
  ls_hs_binomials const *cache = private_kernel_data;
  for (ptrdiff_t batch_idx = 0; batch_idx < batch_size; ++batch_idx) {
    indices[batch_idx * indices_stride] =
        rank_via_combinadics(spins[batch_idx * spins_stride], cache);
  }
}

void ls_hs_state_index_identity_kernel(ptrdiff_t const batch_size,
                                       uint64_t const *spins,
                                       ptrdiff_t const spins_stride,
                                       ptrdiff_t *const restrict indices,
                                       ptrdiff_t const indices_stride,
                                       void const *private_kernel_data) {
  for (ptrdiff_t batch_idx = 0; batch_idx < batch_size; ++batch_idx) {
    indices[batch_idx * indices_stride] =
        (ptrdiff_t)spins[batch_idx * spins_stride];
  }
}
