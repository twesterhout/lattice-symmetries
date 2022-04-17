#include "lattice_symmetries_haskell.h"
#include <assert.h>
#include <complex.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef uint64_t ls_hs_bits;

typedef void (*free_stable_ptr_type)(void *);
static free_stable_ptr_type ls_hs_internal_free_stable_ptr;

void ls_hs_internal_set_free_stable_ptr(free_stable_ptr_type f) {
  ls_hs_internal_free_stable_ptr = f;
}

typedef void (*ls_hs_internal_chpl_free_func)(void *);

void ls_hs_fatal_error(char const *func, char const *message) {
  fprintf(stderr, "[%s] fatal error: %s\n", func, message);
  abort();
}
// void ls_hs_internal_destroy_external_array(chpl_external_array *p) {
//   ls_hs_internal_chpl_free_func free_func =
//       (ls_hs_internal_chpl_free_func)p->freer;
//   if (free_func != NULL) {
//     (*free_func)(p->elts);
//   }
// }

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

void ls_hs_destroy_basis_v2(ls_hs_basis *basis) {
  fprintf(stdout, "ls_hs_destroy_basis_v2 ...\n");
  if (ls_hs_internal_dec_refcount(&basis->refcount) == 1) {
    // Optionally free representatives
    if (basis->representatives.freer != NULL) {
      ls_hs_internal_chpl_free_func free_func =
          (ls_hs_internal_chpl_free_func)basis->representatives.freer;
      (*free_func)(basis->representatives.elts);
    }
    // Free Haskell payload
    LS_CHECK(ls_hs_internal_free_stable_ptr != NULL,
             "ls_hs_internal_free_stable_ptr not set");
    fprintf(stdout, "Calling hs_free_stable_ptr ...\n");
    (*ls_hs_internal_free_stable_ptr)(basis->haskell_payload);
  }
}

void ls_hs_destroy_operator_v2(ls_hs_operator *op) {
  fprintf(stdout, "ls_hs_destroy_basis_v2 ...\n");
  if (ls_hs_internal_dec_refcount(&op->refcount) == 1) {
    if (ls_hs_internal_free_stable_ptr == NULL) {
      fprintf(stderr, "ls_hs_internal_free_stable_ptr not set :(\n");
      abort();
    }
    fprintf(stdout, "Calling hs_free_stable_ptr ...\n");
    (*ls_hs_internal_free_stable_ptr)(op->haskell_payload);
  }
}

#if 0
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
#endif

#if 0
void strided_memset(void *const restrict dest, ptrdiff_t const count,
                    ptrdiff_t const stride, void const *const restrict element,
                    ptrdiff_t const size) {
  for (ptrdiff_t i = 0; i < count; ++i) {
    // Fill one element
    memcpy(dest + i * stride * size, element, size);
  }
}
#endif

static inline void bitstring_and(int const number_words,
                                 uint64_t const *const a,
                                 uint64_t const *const b,
                                 uint64_t *restrict const out) {
  for (int i = 0; i < number_words; ++i) {
    out[i] = a[i] & b[i];
  }
}

static inline void bitstring_xor(int const number_words,
                                 uint64_t const *const a,
                                 uint64_t const *const b,
                                 uint64_t *restrict const out) {
  for (int i = 0; i < number_words; ++i) {
    out[i] = a[i] ^ b[i];
  }
}

static inline bool bitstring_equal(int const number_words,
                                   uint64_t const *const a,
                                   uint64_t const *const b) {
  for (int i = 0; i < number_words; ++i) {
    if (a[i] != b[i]) {
      return false;
    }
  }
  return true;
}

static inline int bitstring_popcount(int const number_words,
                                     uint64_t const *const restrict a) {
  int acc = 0;
  for (int i = 0; i < number_words; ++i) {
    acc += __builtin_popcountl(a[i]);
  }
  return acc;
}

void ls_hs_operator_apply_diag_kernel(ls_hs_operator const *op,
                                      ptrdiff_t const batch_size,
                                      uint64_t const *restrict const alphas,
                                      ptrdiff_t const alphas_stride,
                                      ls_hs_scalar *restrict const coeffs) {
  // fprintf(stderr, "ls_hs_operator_apply_diag_kernel ...\n");
  if (op->diag_terms == NULL ||
      op->diag_terms->number_terms == 0) { // the diagonal is zero
    // fprintf(stderr, "The diagonal is zero ...\n");
    memset(coeffs, 0, (size_t)batch_size * sizeof(ls_hs_scalar));
    return;
  }

  int const number_words = (op->diag_terms->number_bits + 63) / 64;
  // fprintf(stderr, "number_words=%d\n", number_words);
  uint64_t *restrict const temp = malloc(number_words * sizeof(uint64_t));
  if (temp == NULL) {
    // fprintf(stderr, "%s\n", "failed to allocate memory");
    abort();
  }

  ls_hs_nonbranching_terms const *restrict terms = op->diag_terms;
  ptrdiff_t const number_terms = terms->number_terms;
  for (ptrdiff_t batch_idx = 0; batch_idx < batch_size; ++batch_idx) {
    ls_hs_scalar acc = 0;
    uint64_t const *restrict const alpha = alphas + batch_idx * alphas_stride;
    for (ptrdiff_t term_idx = 0; term_idx < number_terms; ++term_idx) {
      ls_hs_scalar const v = terms->v[term_idx];
      uint64_t const *restrict const m = terms->m + term_idx * number_words;
      uint64_t const *restrict const r = terms->r + term_idx * number_words;
      uint64_t const *restrict const s = terms->s + term_idx * number_words;

      bitstring_and(number_words, alpha, m, temp);
      int const delta = bitstring_equal(number_words, temp, r);
      bitstring_and(number_words, alpha, s, temp);
      int const sign = 1 - 2 * (bitstring_popcount(number_words, temp) % 2);
      // fprintf(stderr, "α=%zu, s=%zu, temp=%zu, popcount(temp)=%d, %d\n",
      //         alpha[0], s[0], temp[0], bitstring_popcount(number_words,
      //         temp),
      //         __builtin_popcountl(temp[0]));
      acc += v * (delta * sign);
      // fprintf(stderr, "acc += (%f + %fi) * (%d * %d)\n", crealf(v),
      // cimagf(v),
      //         delta, sign);
    }
    coeffs[batch_idx] = acc;
  }
  free(temp);
}

void ls_hs_operator_apply_off_diag_kernel(
    ls_hs_operator const *op, ptrdiff_t batch_size, uint64_t const *alphas,
    ptrdiff_t alphas_stride, uint64_t *betas, ptrdiff_t betas_stride,
    ls_hs_scalar *coeffs) {
  // fprintf(stderr, "ls_hs_operator_apply_off_diag_kernel ...\n");
  if (op->off_diag_terms == NULL ||
      op->off_diag_terms->number_terms == 0) { // nothing to apply
    // fprintf(stderr, "Nothing to do...\n");
    return;
  }

  int const number_words = (op->off_diag_terms->number_bits + 63) / 64;
  uint64_t *restrict const temp = malloc(number_words * sizeof(uint64_t));
  if (temp == NULL) {
    // fprintf(stderr, "%s\n", "failed to allocate memory");
    abort();
  }

  ls_hs_nonbranching_terms const *restrict terms = op->off_diag_terms;
  ptrdiff_t const number_terms = terms->number_terms;
  for (ptrdiff_t batch_idx = 0; batch_idx < batch_size; ++batch_idx) {
    // fprintf(stderr, "batch_idx=%zi\n", batch_idx);
    uint64_t const *restrict const alpha = alphas + batch_idx * alphas_stride;
    for (ptrdiff_t term_idx = 0; term_idx < number_terms; ++term_idx) {
      // fprintf(stderr, "term_idx=%zi\n", batch_idx);
      ls_hs_scalar const v = terms->v[term_idx];
      uint64_t const *restrict const m = terms->m + term_idx * number_words;
      uint64_t const *restrict const r = terms->r + term_idx * number_words;
      uint64_t const *restrict const x = terms->x + term_idx * number_words;
      uint64_t const *restrict const s = terms->s + term_idx * number_words;
      uint64_t *restrict const beta =
          betas + (batch_idx * number_terms + term_idx) * betas_stride;

      bitstring_and(number_words, alpha, m, temp);
      int const delta = bitstring_equal(number_words, temp, r);
      bitstring_and(number_words, alpha, s, temp);
      // fprintf(stderr, "α=%zu, s=%zu, temp=%zu, popcount(temp)=%d\n",
      // alpha[0],
      //         s[0], temp[0], bitstring_popcount(number_words, temp));
      int const sign = 1 - 2 * (bitstring_popcount(number_words, temp) % 2);
      coeffs[batch_idx * number_terms + term_idx] = v * (delta * sign);
      bitstring_xor(number_words, alpha, x, beta);
      // fprintf(stderr, "coeffs[%zi] = (%f + %fi) * (%d * %d)\n",
      //         batch_idx * number_terms + term_idx, crealf(v), cimagf(v),
      //         delta, sign);
    }
  }
  free(temp);
}

static void compute_binomials(int dim, uint64_t *coeff) {
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

ls_hs_combinadics_kernel_data *
ls_hs_internal_create_combinadics_kernel_data(int number_bits,
                                              bool is_per_sector) {
  ls_hs_combinadics_kernel_data *p =
      malloc(sizeof(ls_hs_combinadics_kernel_data));
  if (p == NULL) {
    goto fail_1;
  }
  p->dimension = number_bits + 1; // NOTE: +1 because we could have n and k
                                  // running from 0 to number_bits inclusive
  p->is_per_sector = is_per_sector;
  p->coefficients = malloc(p->dimension * p->dimension * sizeof(uint64_t));
  if (p->coefficients == NULL) {
    goto fail_2;
  }
  compute_binomials(p->dimension, p->coefficients);
  return p;

fail_2:
  free(p);
fail_1:
  return NULL;
}

void ls_hs_internal_destroy_combinadics_kernel_data(
    ls_hs_combinadics_kernel_data *p) {
  if (p != NULL) {
    free(p->coefficients);
  }
  free(p);
}

static inline uint64_t binomial(int const n, int const k,
                                ls_hs_combinadics_kernel_data const *cache) {
  if (k > n) {
    return 0;
  }
  assert(0 <= n && n < cache->dimension);
  assert(0 <= k && k < cache->dimension);
  return cache->coefficients[n * cache->dimension + k];
}

static inline ptrdiff_t
rank_via_combinadics(uint64_t alpha,
                     ls_hs_combinadics_kernel_data const *cache) {
  ptrdiff_t i = 0;
  for (int k = 1; alpha != 0; ++k) {
    int c = __builtin_ctzl(alpha);
    alpha &= alpha - 1;
    i += binomial(c, k, cache);
  }
  return i;
}

void ls_hs_state_index_combinadics_kernel(ptrdiff_t const batch_size,
                                          uint64_t const *spins,
                                          ptrdiff_t const spins_stride,
                                          ptrdiff_t *const restrict indices,
                                          ptrdiff_t const indices_stride,
                                          void const *private_kernel_data) {
  if (batch_size == 0) {
    return;
  }
  ls_hs_combinadics_kernel_data const *cache = private_kernel_data;

  int const number_bits = cache->dimension - 1;
  // The following are only needed for fermionic systems when number of
  // particles per spin sector is fixed
  uint64_t const mask =
      number_bits >= 64 ? ~(uint64_t)0 : ((uint64_t)1 << number_bits) - 1;
  int const hamming_weight_low = __builtin_popcountl(spins[0] & mask);
  ptrdiff_t const number_low = binomial(number_bits, hamming_weight_low, cache);

  for (ptrdiff_t batch_idx = 0; batch_idx < batch_size; ++batch_idx) {
    uint64_t alpha = spins[batch_idx * spins_stride];
    ptrdiff_t index;
    if (cache->is_per_sector) {
      uint64_t const low = alpha & mask;
      uint64_t const high = alpha >> number_bits;
      index = rank_via_combinadics(high, cache) * number_low +
              rank_via_combinadics(low, cache);
    } else {
      index = rank_via_combinadics(alpha, cache);
    }
    indices[batch_idx * indices_stride] = index;
  }
}

void ls_hs_state_index_identity_kernel(ptrdiff_t const batch_size,
                                       uint64_t const *spins,
                                       ptrdiff_t const spins_stride,
                                       ptrdiff_t *const restrict indices,
                                       ptrdiff_t const indices_stride,
                                       void const *private_kernel_data) {
  (void)private_kernel_data;
  for (ptrdiff_t batch_idx = 0; batch_idx < batch_size; ++batch_idx) {
    indices[batch_idx * indices_stride] =
        (ptrdiff_t)spins[batch_idx * spins_stride];
  }
}

void ls_hs_state_index(ls_hs_basis const *const basis,
                       ptrdiff_t const batch_size,
                       uint64_t const *const restrict spins,
                       ptrdiff_t const spins_stride,
                       ptrdiff_t *const restrict indices,
                       ptrdiff_t const indices_stride) {
  if (basis->kernels->state_index_kernel == NULL) {
    fprintf(stderr, "state_index_kernel is NULL\n");
    abort();
  }
  (*basis->kernels->state_index_kernel)(batch_size, spins, spins_stride,
                                        indices, indices_stride,
                                        basis->kernels->state_index_data);
}

void ls_hs_evaluate_wavefunction_via_statevector(
    ls_hs_basis_kernels const *const kernels, ptrdiff_t const batch_size,
    uint64_t const *const restrict alphas, ptrdiff_t const alphas_stride,
    void const *const restrict state_vector, size_t const element_size,
    void *const restrict coeffs) {
  ptrdiff_t *const indices = malloc(batch_size * sizeof(ptrdiff_t));
  if (indices == NULL) {
    fprintf(stderr, "failed to allocate space for indices\n");
    abort();
  }

  if (kernels->state_index_kernel == NULL) {
    free(indices);
    fprintf(stderr, "no kernel to compute indices\n");
    abort();
  }
  (*kernels->state_index_kernel)(batch_size, alphas, alphas_stride, indices, 1,
                                 kernels->state_index_data);
  for (ptrdiff_t batch_idx = 0; batch_idx < batch_size; ++batch_idx) {
    void const *src =
        (uint8_t const *)state_vector + indices[batch_idx] * element_size;
    void *dest = (uint8_t *)coeffs + batch_idx * element_size;
    memcpy(dest, src, element_size);
  }

  free(indices);
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
  if (basis->kernels->state_index_kernel == NULL) {
    basis->kernels->state_index_data =
        ls_hs_create_state_index_binary_search_kernel_data(
            &basis->representatives);
    basis->kernels->state_index_kernel =
        &ls_hs_state_index_binary_search_kernel;
  }
}

// TODO: currently not thread-safe, fix it
static ls_chpl_kernels global_chpl_kernels = {.enumerate_states = NULL};
// static pthread_mutex_t global_chpl_kernels_mutex = PTHREAD_MUTEX_INITIALIZER;

ls_chpl_kernels const *ls_hs_internal_get_chpl_kernels() {
  return &global_chpl_kernels;
}
void ls_hs_internal_set_chpl_kernels(ls_chpl_kernels const *kernels) {
  global_chpl_kernels = *kernels;
}
