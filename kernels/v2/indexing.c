#include "lattice_symmetries_types.h"
#include <HalideRuntime.h>
#include <assert.h>
#include <callback.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct state_to_index_binary_search_data {
    ptrdiff_t number_states;
    uint64_t const *representatives;
    unsigned number_bits;
    unsigned shift;
    ptrdiff_t range_size;
    ptrdiff_t number_offsets;
    ptrdiff_t *offsets;
};

static ptrdiff_t normalize_offset_ranges(ptrdiff_t const number_offsets, ptrdiff_t offsets[]) {
    ptrdiff_t max_range_size = 0;
    for (ptrdiff_t i = 0; i < number_offsets - 1; ++i) {
        ptrdiff_t const size = offsets[i + 1] - offsets[i];
        if (size > max_range_size) {
            max_range_size = size;
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

static void generate_offset_ranges(struct state_to_index_binary_search_data *cache) {
    LS_CHECK(0 < cache->number_bits && cache->number_bits < 64, "invalid number_bits");
    LS_CHECK(cache->shift < 64, "invalid shift");
    uint64_t const *first = cache->representatives;
    uint64_t const *const last = cache->representatives + cache->number_states;
    uint64_t const *const begin = first;

    ptrdiff_t const size = (ptrdiff_t)((uint64_t)1 << cache->number_bits);
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

    cache->range_size = normalize_offset_ranges(cache->number_offsets, cache->offsets);
}

static int64_t binary_search_x1(int64_t const haystack_size, uint64_t const haystack[],
                                uint64_t const needle, ptrdiff_t guess) {
    uint64_t const *base = haystack + guess;
    ptrdiff_t size = haystack_size;

    while (size > 1) {
        ptrdiff_t const half = size / 2;
        size -= half;
        base = (base[half] < needle) ? base + half : base;
    }

    base += *base < needle;
    return (*base == needle) ? base - haystack : -1;
}

static void
ls_hs_state_to_index_binary_search_kernel(halide_buffer_t const *basis_states,
                                          halide_buffer_t const *indices,
                                          struct state_to_index_binary_search_data const *ctx) {
    LS_CHECK(basis_states->dimensions == 1 && basis_states->dim[0].stride == 1 &&
                 indices->dimensions == 1 && indices->dim[0].stride == 1 &&
                 basis_states->dim[0].min == 0 && indices->dim[0].min == 0 &&
                 basis_states->dim[0].extent == indices->dim[0].extent,
             "basis_states and indices must be contiguous 1d tensors of matching shape");
    // fprintf(stderr, "ls_hs_state_to_index_binary_search_kernel...\n");
    // fprintf(stderr, "state_to_index_binary_search_data { number_states = %zi, representatives = %p, number_bits = %u, shift = %u, range_size = %zi, number_offsets = %zi, offsets = %p }\n", ctx->number_states, ctx->representatives, ctx->number_bits, ctx->shift, ctx->range_size, ctx->number_offsets, ctx->offsets);
    // for (ptrdiff_t k = 0; k < ctx->number_offsets; ++k) {
    //     fprintf(stderr, "offsets[%zi] = %zi\n", k, ctx->offsets[k]);
    // }

    ptrdiff_t const batch_size = basis_states->dim[0].extent;
    uint64_t const *spins = (uint64_t const *)basis_states->host;
    int64_t *out = (int64_t *)indices->host;

    for (ptrdiff_t batch_idx = 0; batch_idx < batch_size; ++batch_idx) {
        uint64_t const needle = spins[batch_idx];
        if (needle >= ((uint64_t)1 << ctx->number_bits)) {
            out[batch_idx] = -1;
        }
        else {
            int64_t const guess = ctx->offsets[needle >> ctx->shift];
            out[batch_idx] = binary_search_x1(ctx->range_size, ctx->representatives, needle, guess);
        }
    }
}

static void invoke_binary_search_state_to_index_kernel(void *data, va_alist alist) {
    struct state_to_index_binary_search_data const *ctx =
        (struct state_to_index_binary_search_data const *)data;
    va_start_void(alist);
    halide_buffer_t *basis_states = va_arg_ptr(alist, halide_buffer_t *);
    halide_buffer_t *indices = va_arg_ptr(alist, halide_buffer_t *);
    ls_hs_state_to_index_binary_search_kernel(basis_states, indices, ctx);
    va_return_void(alist);
}

ls_hs_state_to_index_kernel_type ls_hs_internal_mk_binary_search_state_to_index_kernel(
    int64_t const number_representatives, uint64_t const *representatives,
    unsigned const number_bits, unsigned const prefix_bits) {

    struct state_to_index_binary_search_data *cxt =
        malloc(sizeof(struct state_to_index_binary_search_data));
    LS_CHECK(cxt != NULL, "malloc failed");

    cxt->number_states = number_representatives;
    cxt->representatives = representatives;
    cxt->number_bits = prefix_bits;
    if (cxt->number_bits > number_bits) {
        cxt->number_bits = number_bits;
    }
    cxt->shift = number_bits - cxt->number_bits;
    cxt->range_size = 0;
    cxt->number_offsets = 0;
    cxt->offsets = NULL;
    if (cxt->number_bits > 0) {
        generate_offset_ranges(cxt);
    }

    callback_t closure = alloc_callback(&invoke_binary_search_state_to_index_kernel, cxt);
    return (ls_hs_state_to_index_kernel_type)closure;
}

void ls_hs_internal_destroy_binary_search_state_to_index_kernel(
    ls_hs_state_to_index_kernel_type closure) {
    if (closure == NULL) {
        return;
    }
    LS_CHECK(is_callback((void *)closure), "trying to destroy a normal function pointer. This "
                                           "should never happen. Please, submit a bug report.");

    struct state_to_index_binary_search_data *cxt = callback_data((callback_t)closure);
    LS_CHECK(cxt != NULL,
             "callback_data is NULL. This should never happen. Please, submit a bug report.");

    if (cxt->offsets != NULL) {
        free(cxt->offsets);
    }
    free(cxt);
    free_callback((callback_t)closure);
}
