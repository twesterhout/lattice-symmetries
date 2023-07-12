#include "lattice_symmetries_types.h"
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
  int number_bits;
  int shift;
  ptrdiff_t range_size;
  ptrdiff_t number_offsets;
  ptrdiff_t *offsets;
};

#if 0
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
  uint64_t const *element_ptr = bsearch(&needle, haystack, (size_t)size,
                                        sizeof(uint64_t), &compare_uint64);
  if (element_ptr == NULL) {
    return -1;
  }
  LS_CHECK(*element_ptr == needle, "search failed... :(");
  return element_ptr - haystack;
}
#endif

static void
generate_offset_ranges(ls_hs_state_index_binary_search_data *cache) {
  LS_CHECK(0 < cache->number_bits && cache->number_bits < 64, "invalid bits");
  LS_CHECK(cache->shift < 64, "invalid shift");
  // auto const        extract_relevant = [shift](auto const x) noexcept {
  // return x >> shift; };
  uint64_t const *first = cache->representatives;
  uint64_t const *const last = cache->representatives + cache->number_states;
  uint64_t const *const begin = first;

  ptrdiff_t const size = (ptrdiff_t)1 << cache->number_bits;
  cache->number_offsets = size + 1;
  cache->offsets = malloc((size_t)cache->number_offsets * sizeof(ptrdiff_t));
  LS_CHECK(cache->offsets, "malloc failed");
  for (ptrdiff_t i = 0; i < size; ++i) {
    cache->offsets[i] = first - begin;
    while (first != last && ((*first) >> cache->shift) == (uint64_t)i) {
      ++first;
    }
  }
  cache->offsets[size] = first - begin;
  LS_CHECK(first == last, "not all states checked");
}

static ptrdiff_t normalize_offset_ranges(ptrdiff_t const number_offsets,
                                         ptrdiff_t offsets[]) {
  ptrdiff_t max_range_size = 0;
  for (ptrdiff_t i = 0; i < number_offsets - 1; ++i) {
    ptrdiff_t const n = offsets[i + 1] - offsets[i];
    if (n > max_range_size) {
      max_range_size = n;
    }
  }

  ptrdiff_t number_states = offsets[number_offsets - 1];
  for (ptrdiff_t i = 0; i < number_offsets - 1; ++i) {
    if (offsets[i] > number_states - max_range_size) {
      offsets[i] = number_states - max_range_size;
    }
  }

  return max_range_size;
}

ls_hs_state_index_binary_search_data *
ls_hs_create_state_index_binary_search_kernel_data(
    chpl_external_array const *representatives, int const number_bits,
    int const prefix_bits) {
  ls_hs_state_index_binary_search_data *cache =
      malloc(sizeof(ls_hs_state_index_binary_search_data));
  LS_CHECK(cache != NULL, "malloc failed");
  cache->number_states = (ptrdiff_t)representatives->num_elts;
  cache->representatives = representatives->elts;
  cache->number_bits = prefix_bits;
  if (cache->number_bits > number_bits) {
    cache->number_bits = number_bits;
  }
  cache->shift = number_bits - cache->number_bits;
  cache->range_size = 0;
  cache->number_offsets = 0;
  cache->offsets = NULL;
  if (cache->number_bits > 0) {
    generate_offset_ranges(cache);
    cache->range_size =
        normalize_offset_ranges(cache->number_offsets, cache->offsets);
  }
  return cache;
}

void ls_hs_destroy_state_index_binary_search_kernel_data(
    ls_hs_state_index_binary_search_data *cache) {
  if (cache->offsets != NULL) {
    free(cache->offsets);
  }
  free(cache);
}

inline size_t bit_floor(size_t const i)
{
    int const num_bits = sizeof(i) * 8;
    return (size_t)1 << (num_bits - __builtin_clzl(i) - 1);
}
inline size_t bit_ceil(size_t const i)
{
    int const num_bits = sizeof(i) * 8;
    return (size_t)1 << (num_bits - __builtin_clzl(i - 1));
}

static inline
uint64_t const* branchless_lower_bound(uint64_t const* begin,
                                       uint64_t const* const end,
                                       uint64_t const value)
{
    size_t length = end - begin;
    if (length == 0) { return end; }

    size_t step = bit_floor(length);
    if (step != length && begin[step] < value)
    {
        length -= step + 1;
        if (length == 0) { return end; }
        step = bit_ceil(length);
        begin = end - step;
    }

    for (step /= 2; step != 0; step /= 2)
    {
        if (begin[step] < value) { begin += step; }
    }

    return begin + (*begin < value);
}

static inline void
branchless_binary_search_x1(uint64_t const haystack[],
                            ptrdiff_t const haystack_size,
                            uint64_t const needle, ptrdiff_t *index) {
  uint64_t const *begin = haystack + *index;
  uint64_t const *end = begin + haystack_size;
  begin = branchless_lower_bound(begin, end, needle);
  *index = (begin != end && *begin == needle) ? begin - haystack : -1;
}

static inline void
normal_binary_search_x1(uint64_t const haystack[],
                        ptrdiff_t const haystack_size,
                        uint64_t const needle, ptrdiff_t *index) {
  ptrdiff_t n = haystack_size;
  uint64_t const* first = haystack + *index;

  while (n > 0) {
    ptrdiff_t const half = n >> 1;
    uint64_t const* middle = first + half;
    if (*middle < needle) {
      first = middle + 1;
      n = n - half - 1;
    }
    else {
      n = half;
    }
  }

  ptrdiff_t k = first - haystack;
  return (k < haystack_size && *first == needle) ? k : -1;
}

static inline void
ls_hs_internal_binary_search_x1(uint64_t const haystack[],
                                ptrdiff_t const haystack_size,
                                uint64_t const needle, ptrdiff_t *index) {
  uint64_t const *base = haystack + *index;
  ptrdiff_t n = haystack_size;

  while (n > 1) {
    ptrdiff_t const half = n / 2;
    // __builtin_prefetch(base + half / 2, 0, 0);
    // __builtin_prefetch(base + half + half / 2, 0, 0);
    // fprintf(stderr, "i=%zi, n=%zi, half=%zi, needle=%zu, base[half]=%zu\n",
    // base - haystack, n, half, needle, base[half]);
    n -= half;
    base = (base[half] < needle) ? base + half : base;
  }

  base += *base < needle;
  *index = (*base == needle) ? base - haystack : -1;
}

static inline void
ls_hs_internal_binary_search_x8(uint64_t const haystack[],
                                ptrdiff_t const range_size, // haystack_size,
                                uint64_t const needles[], ptrdiff_t indices[]) {
  uint64_t const *base[8] = {haystack + indices[0], haystack + indices[1],
                             haystack + indices[2], haystack + indices[3],
                             haystack + indices[4], haystack + indices[5],
                             haystack + indices[6], haystack + indices[7]};
  ptrdiff_t n = range_size; // haystack_size;

  while (n > 1) {
    ptrdiff_t const half = n / 2;
    __builtin_prefetch(base[0] + half, 0, 0);
    __builtin_prefetch(base[1] + half, 0, 0);
    __builtin_prefetch(base[2] + half, 0, 0);
    __builtin_prefetch(base[3] + half, 0, 0);
    __builtin_prefetch(base[4] + half, 0, 0);
    __builtin_prefetch(base[5] + half, 0, 0);
    __builtin_prefetch(base[6] + half, 0, 0);
    __builtin_prefetch(base[7] + half, 0, 0);
    // __builtin_prefetch(base0 + half + half / 2, 0, 0);
    // __builtin_prefetch(base1 + half / 2, 0, 0);
    // __builtin_prefetch(base1 + half + half / 2, 0, 0);
    // fprintf(stderr, "i=%zi, n=%zi, half=%zi, needle=%zu, base[half]=%zu\n",
    // base - haystack, n, half, needle, base[half]);
    n -= half;
    base[0] = (base[0][half] < needles[0]) ? base[0] + half : base[0];
    base[1] = (base[1][half] < needles[1]) ? base[1] + half : base[1];
    base[2] = (base[2][half] < needles[2]) ? base[2] + half : base[2];
    base[3] = (base[3][half] < needles[3]) ? base[3] + half : base[3];
    base[4] = (base[4][half] < needles[4]) ? base[4] + half : base[4];
    base[5] = (base[5][half] < needles[5]) ? base[5] + half : base[5];
    base[6] = (base[6][half] < needles[6]) ? base[6] + half : base[6];
    base[7] = (base[7][half] < needles[7]) ? base[7] + half : base[7];
  }

  base[0] += *(base[0]) < needles[0];
  base[1] += *(base[1]) < needles[1];
  base[2] += *(base[2]) < needles[2];
  base[3] += *(base[3]) < needles[3];
  base[4] += *(base[4]) < needles[4];
  base[5] += *(base[5]) < needles[5];
  base[6] += *(base[6]) < needles[6];
  base[7] += *(base[7]) < needles[7];
  // fprintf(stderr, "i=%zi, n=%zi, half=%zi, needle=%zu, base[half]=%zu\n",
  // base - haystack, n, 0, needle, base[0]);
  indices[0] = (*(base[0]) == needles[0]) ? base[0] - haystack : -1;
  indices[1] = (*(base[1]) == needles[1]) ? base[1] - haystack : -1;
  indices[2] = (*(base[2]) == needles[2]) ? base[2] - haystack : -1;
  indices[3] = (*(base[3]) == needles[3]) ? base[3] - haystack : -1;
  indices[4] = (*(base[4]) == needles[4]) ? base[4] - haystack : -1;
  indices[5] = (*(base[5]) == needles[5]) ? base[5] - haystack : -1;
  indices[6] = (*(base[6]) == needles[6]) ? base[6] - haystack : -1;
  indices[7] = (*(base[7]) == needles[7]) ? base[7] - haystack : -1;
}

void ls_hs_state_index_binary_search_kernel(ptrdiff_t const batch_size,
                                            uint64_t const *spins,
                                            ptrdiff_t const spins_stride,
                                            ptrdiff_t *const restrict indices,
                                            ptrdiff_t const indices_stride,
                                            void const *private_kernel_data) {
  ls_hs_state_index_binary_search_data const *cache = private_kernel_data;
  LS_CHECK(indices_stride == 1, "expected indices_stride==1");
  LS_CHECK(spins_stride == 1, "expected spins_stride==1");

#if 1
  ptrdiff_t const block_size = 8;
  ptrdiff_t batch_idx = 0;
  // for (; batch_idx + block_size <= batch_size; batch_idx += block_size) {
  //   for (int i = 0; i < block_size; ++i) {
  //     uint64_t const spin = spins[batch_idx + i];
  //     indices[batch_idx + i] = cache->offsets[spin >> cache->shift];
  //   }
  //   ls_hs_internal_binary_search_x8(cache->representatives, cache->range_size,
  //                                   spins + batch_idx, indices + batch_idx);
  // }
  for (; batch_idx < batch_size; ++batch_idx) {
    uint64_t const spin = spins[batch_idx];
    indices[batch_idx] = cache->offsets[spin >> cache->shift];
    normal_binary_search_x1(cache->representatives, cache->range_size,
                            spins[batch_idx], indices + batch_idx);
    // branchless_binary_search_x1(cache->representatives, cache->range_size,
    //                             spins[batch_idx], indices + batch_idx);
    // ls_hs_internal_binary_search_x1(cache->representatives, cache->range_size,
    //                                 spins[batch_idx], indices + batch_idx);
  }

#else

  for (ptrdiff_t batch_idx = 0; batch_idx < batch_size; ++batch_idx) {
#if true
    uint64_t const spin = spins[batch_idx];
    uint64_t const i =
        (spin >> cache->shift); // & (((uint64_t)1 << cache->number_bits) - 1);
    ptrdiff_t const offset = cache->offsets[i];
    ptrdiff_t const count = cache->range_size; // offsets[i + 1] - offset;
    ptrdiff_t const index =
        binary_search(count, cache->representatives + offset, spin);
    // LS_CHECK(index == -1 || cache->representatives[index] == spin, ":(");
    indices[batch_idx] = offset + index;
#else
    ls_hs_internal_binary_search_x1(cache->representatives,
                                    cache->number_states, spins[batch_idx],
                                    indices + batch_idx);
#endif
  }
#endif
}
