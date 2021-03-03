#include <lattice_symmetries/lattice_symmetries.h>
#include <omp.h>

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
    uint64_t const chunk_size    = max(count / (unsigned)omp_get_max_threads(), 100);
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
