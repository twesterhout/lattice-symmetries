#include "lattice_symmetries_haskell.h"
#include <assert.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct ls_hs_state_index_binary_search_data {
  ptrdiff_t number_states;
  uint64_t *representatives;
  int shift;
  int number_bits;
  ptrdiff_t number_offsets;
  uint64_t *offsets;
};

static int compare_uint64(const void *_a, const void *_b) {
  uint64_t const a = *(uint64_t const *)_a;
  uint64_t const b = *(uint64_t const *)_b;
  if (a == b) {
    return 0;
  }
  if (a < b) {
    return -1;
  }
  return 1;
}

static ptrdiff_t binary_search(ptrdiff_t const size, uint64_t const *haystack,
                               uint64_t const needle) {
  uint64_t const *element_ptr =
      bsearch(&needle, haystack, size, sizeof(uint64_t), &compare_uint64);
  if (element_ptr == NULL) {
    return -1;
  }
  return element_ptr - haystack;
}

ls_hs_state_index_binary_search_data *
ls_hs_create_state_index_binary_search_kernel_data(
    chpl_external_array const *representatives) {
  ls_hs_state_index_binary_search_data *cache =
      malloc(sizeof(ls_hs_state_index_binary_search_data));
  cache->number_states = representatives->num_elts;
  cache->representatives = representatives->elts;
  cache->shift = 0;
  cache->number_bits = 0;
  cache->number_offsets = 0;
  cache->offsets = NULL;
}

void ls_hs_destroy_state_index_binary_search_kernel_data(
    ls_hs_state_index_binary_search_data *cache) {
  free(cache);
}

void ls_hs_state_index_binary_search_kernel(ptrdiff_t const batch_size,
                                            uint64_t const *spins,
                                            ptrdiff_t const spins_stride,
                                            ptrdiff_t *const restrict indices,
                                            ptrdiff_t const indices_stride,
                                            void const *private_kernel_data) {
  ls_hs_state_index_binary_search_data const *cache = private_kernel_data;
  for (ptrdiff_t batch_idx = 0; batch_idx < batch_size; ++batch_idx) {
    uint64_t const spin = spins[batch_idx * spins_stride];
    ptrdiff_t const index =
        binary_search(cache->number_states, cache->representatives, spin);
    // LS_CHECK(index == -1 || cache->representatives[index] == spin, ":(");
    indices[batch_idx * indices_stride] = index;
  }
}
