#include <complex.h>
#include <lattice_symmetries/lattice_symmetries.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>

#define L1_CACHE_SIZE 64

static inline uint64_t max(uint64_t const a, uint64_t const b) { return (a > b) ? a : b; }
static inline uint64_t min(uint64_t const a, uint64_t const b) { return (a < b) ? a : b; }

static ls_error_code ls_batched_get_index_serial(ls_spin_basis const* const basis,
                                                 uint64_t const count, ls_bits64 const* const spins,
                                                 uint64_t const spins_stride, uint64_t* const out,
                                                 uint64_t const out_stride)
{
    for (uint64_t i = 0; i < count; ++i) {
        uint64_t      index;
        ls_error_code status = ls_get_index(basis, spins[i * spins_stride], &index);
        if (LATTICE_SYMMETRIES_UNLIKELY(status != LS_SUCCESS)) { return status; }
        out[i * out_stride] = index;
    }
    return LS_SUCCESS;
}

LATTICE_SYMMETRIES_EXPORT ls_error_code ls_batched_get_index(
    ls_spin_basis const* const basis, uint64_t const count, ls_bits64 const* const spins,
    uint64_t const spins_stride, uint64_t* const out, uint64_t const out_stride)
{
    _Alignas(L1_CACHE_SIZE) ls_error_code status = LS_SUCCESS;

    uint64_t const chunk_size    = max(count / 10, 100);
    uint64_t const number_chunks = (count + (chunk_size - 1)) / chunk_size;
#pragma omp parallel for default(none) num_threads(2) schedule(dynamic, 1)                         \
    firstprivate(basis, count, spins, spins_stride, out, out_stride, chunk_size, number_chunks)    \
        shared(status)
    for (uint64_t i = 0; i < number_chunks; ++i) {
        ls_error_code local_status; // NOLINT: initialized by atomic read
#pragma omp atomic read
        local_status = status;
        if (LATTICE_SYMMETRIES_UNLIKELY(local_status != LS_SUCCESS)) { continue; }

        uint64_t const  local_count = min(chunk_size, count - i * chunk_size);
        uint64_t const* local_spins = spins + i * chunk_size * spins_stride;
        uint64_t*       local_out   = out + i * chunk_size * out_stride;
        local_status = ls_batched_get_index_serial(basis, local_count, local_spins, spins_stride,
                                                   local_out, out_stride);
        if (LATTICE_SYMMETRIES_UNLIKELY(local_status != LS_SUCCESS)) {
#pragma omp atomic write
            status = local_status;
        }
    }
    return status;
}

LATTICE_SYMMETRIES_EXPORT void
ls_batched_get_state_info(ls_spin_basis const* const basis, uint64_t const count,
                          ls_bits512 const* const spins, uint64_t const spins_stride,
                          ls_bits512* const repr, uint64_t const repr_stride,
                          _Complex double* const eigenvalues, uint64_t const eigenvalues_stride,
                          double* const norm, uint64_t const norm_stride)
{
    uint64_t const chunk_size = max(count / (unsigned)omp_get_max_threads(), 100);
#pragma omp parallel for default(none) schedule(dynamic, chunk_size)                               \
    firstprivate(basis, chunk_size, count, spins, spins_stride, repr, repr_stride, eigenvalues,    \
                 eigenvalues_stride, norm, norm_stride)
    for (uint64_t i = 0; i < count; ++i) {
        ls_get_state_info(basis, spins + i * spins_stride, repr + i * repr_stride,
                          eigenvalues + i * eigenvalues_stride, norm + i * norm_stride);
    }
}

LATTICE_SYMMETRIES_EXPORT
void ls_batched_apply_symmetry(ls_symmetry const* symmetry, uint64_t const count,
                               uint64_t* const spins, uint64_t const stride)
{
    int num_threads = omp_get_max_threads();
    if (num_threads > 2) { num_threads = 2; }
#pragma omp parallel for default(none) num_threads(num_threads)                                    \
    firstprivate(count, spins, stride, symmetry)
    for (uint64_t i = 0; i < count; ++i) {
        ls_apply_symmetry(symmetry, (ls_bits512*)(spins + i * stride));
    }
}

typedef struct store_callback_ctx_t {
    uint64_t               size;
    ls_bits512* const      spins;
    _Complex double* const coeffs;
} store_callback_ctx_t;

static ls_error_code store_callback(ls_bits512 const* const bits, void const* const coeff,
                                    void* const _cxt)
{
    store_callback_ctx_t* cxt = _cxt;
    // LATTICE_SYMMETRIES_LOG_DEBUG(
    //     "Writing [%zu, %zu, %zu, %zu, %zu, %zu, %zu, %zu] (%f, %f) to %zu...\n", bits->words[0],
    //     bits->words[1], bits->words[2], bits->words[3], bits->words[4], bits->words[5],
    //     bits->words[6], bits->words[7], creal(*(_Complex double const*)coeff),
    //     cimag(*(_Complex double const*)coeff), cxt->size);
    cxt->spins[cxt->size]  = *bits;
    cxt->coeffs[cxt->size] = conj(*(_Complex double const*)coeff);
    ++cxt->size;
    return LS_SUCCESS;
}

LATTICE_SYMMETRIES_EXPORT
uint64_t ls_batched_operator_apply(ls_operator const* operator, uint64_t const count,
                                   ls_bits512 const* const spins, ls_bits512* const out_spins,
                                   _Complex double* out_coeffs, uint64_t* out_counts)
{
    uint64_t const max_buffer_size = ls_operator_max_buffer_size(operator);
    uint64_t const num_threads     = omp_in_parallel() ? 1 : omp_get_max_threads();
    uint64_t const chunk_size      = (count + (num_threads - 1)) / num_threads;
    uint64_t*      chunk_offsets   = malloc(num_threads * sizeof(uint64_t));
    uint64_t*      chunk_counts    = malloc(num_threads * sizeof(uint64_t));
    LATTICE_SYMMETRIES_CHECK(chunk_offsets != NULL, "failed to allocate storage for offsets");
    LATTICE_SYMMETRIES_CHECK(chunk_counts != NULL, "failed to allocate storage for counts");

#pragma omp parallel if (num_threads > 1) num_threads(num_threads) default(none)                   \
    firstprivate(operator, chunk_size, chunk_offsets, chunk_counts, count, spins, max_buffer_size, \
                 out_spins, out_coeffs, out_counts)
    {
        uint64_t const thread_id     = omp_get_thread_num();
        uint64_t const begin         = thread_id * chunk_size;
        uint64_t const end           = min(count, begin + chunk_size);
        uint64_t       total_written = 0U;
        for (uint64_t i = begin; i < end; ++i) {
            uint64_t const       offset = begin * max_buffer_size + total_written;
            store_callback_ctx_t cxt    = {0, out_spins + offset, out_coeffs + offset};
            ls_error_code const  status =
                ls_operator_apply(operator, spins + i, &store_callback, &cxt);
            LATTICE_SYMMETRIES_ASSERT(status == LS_SUCCESS, NULL);
            out_counts[i] = cxt.size;
            total_written += cxt.size;
        }
        chunk_offsets[thread_id] = begin * max_buffer_size;
        chunk_counts[thread_id]  = total_written;
    }

    LATTICE_SYMMETRIES_ASSERT(chunk_offsets[0] == 0, NULL);
    uint64_t offset = chunk_counts[0];
    for (uint64_t i = 1; i < num_threads; ++i) {
        memmove(out_spins + offset, out_spins + chunk_offsets[i],
                sizeof(ls_bits512) * chunk_counts[i]);
        memmove(out_coeffs + offset, out_coeffs + chunk_offsets[i],
                sizeof(_Complex double) * chunk_counts[i]);
        offset += chunk_counts[i];
    }
    free(chunk_offsets);
    free(chunk_counts);
    return offset;
}
