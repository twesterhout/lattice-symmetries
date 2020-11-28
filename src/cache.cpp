// Copyright (c) 2019-2020, Tom Westerhout
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "cache.hpp"
#include "error_handling.hpp"
#include "kernels.hpp"
#include "state_info.hpp"

#if defined(__APPLE__)
#    include <libkern/OSByteOrder.h>
#    include <machine/endian.h>

#    define htobe16(x) OSSwapHostToBigInt16(x)
#    define htole16(x) OSSwapHostToLittleInt16(x)
#    define be16toh(x) OSSwapBigToHostInt16(x)
#    define le16toh(x) OSSwapLittleToHostInt16(x)

#    define htobe32(x) OSSwapHostToBigInt32(x)
#    define htole32(x) OSSwapHostToLittleInt32(x)
#    define be32toh(x) OSSwapBigToHostInt32(x)
#    define le32toh(x) OSSwapLittleToHostInt32(x)

#    define htobe64(x) OSSwapHostToBigInt64(x)
#    define htole64(x) OSSwapHostToLittleInt64(x)
#    define be64toh(x) OSSwapBigToHostInt64(x)
#    define le64toh(x) OSSwapLittleToHostInt64(x)
#else
#    include <endian.h>
#endif
#include <sys/stat.h>

#include <omp.h>
#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <numeric>

#include <unordered_map>

namespace lattice_symmetries {

namespace {
    auto save_histogram(char const* filename, std::unordered_map<uint64_t, uint64_t> const& data)
        -> void
    {
        std::vector<std::pair<uint64_t, uint64_t>> table;
        table.reserve(data.size());
        for (auto const& e : data) {
            table.push_back(e);
        }
        std::sort(std::begin(table), std::end(table),
                  [](auto const& x, auto const& y) { return x.first < y.first; });

        auto* file = std::fopen(filename, "w");
        LATTICE_SYMMETRIES_CHECK(file != nullptr, "failed to open file");
        for (auto const& e : table) {
            std::fprintf(file, "%zu\t%zu\n", e.first, e.second);
        }
        std::fclose(file);
    }

    auto generate_ranges_helper(tcb::span<uint64_t const> states, unsigned const bits,
                                unsigned const shift, std::vector<range_node_t>& ranges,
                                int64_t const offset) -> void
    {
        LATTICE_SYMMETRIES_ASSERT(0 < bits && bits <= 32, "invalid bits");
        constexpr auto empty            = range_node_t::make_empty();
        auto const     size             = uint64_t{1} << bits;
        auto const     mask             = size - 1;
        auto const     extract_relevant = [shift, mask](auto const x) noexcept {
            return (x >> shift) & mask;
        };

        ranges.reserve(ranges.size() + size);
        auto const*       first = states.data();
        auto const* const last  = first + states.size();
        auto const* const begin = first;
        for (auto i = uint64_t{0}; i < size; ++i) {
            auto element = range_node_t::make_empty();
            if (first != last && extract_relevant(*first) == i) {
                element.start = first - begin + offset;
                ++first;
                ++element.size;
                while (first != last && extract_relevant(*first) == i) {
                    ++first;
                    ++element.size;
                }
            }
            ranges.push_back(element);
        }
        LATTICE_SYMMETRIES_CHECK(first == last, "not all states checked");
    }

    auto generate_ranges_v2(tcb::span<uint64_t const> states, unsigned const bits,
                            unsigned const shift) -> std::vector<range_node_t>
    {
        constexpr auto            cutoff = 512U;
        std::vector<range_node_t> ranges;
        generate_ranges_helper(states, bits, shift, ranges, 0);

        auto next_shift = shift;
        for (;;) {
            auto       done = true;
            auto const size = ranges.size();

            next_shift -= bits;
            for (auto i = uint64_t{0}; i < size; ++i) {
                auto& r = ranges[i];
                if (r.is_range() && r.size > cutoff) {
                    done = false;
                    LATTICE_SYMMETRIES_CHECK(next_shift <= 64, nullptr);
                    auto const next_states = states.subspan(static_cast<uint64_t>(r.start),
                                                            static_cast<uint64_t>(r.size));
                    auto const offset      = r.start;

                    r.start = -static_cast<int64_t>(ranges.size());
                    r.size  = 0;
                    generate_ranges_helper(next_states, bits, next_shift, ranges, offset);
                }
            }
            if (done) { break; }
        }

        // std::unordered_map<uint64_t, uint64_t> histogram;
        // for (auto const& r : ranges) {
        //     if (r.is_range()) { ++histogram[static_cast<uint64_t>(r.size)]; }
        // }
        // save_histogram("ranges_histogram_v2.dat", histogram);
        return ranges;
    }

    auto generate_ranges(tcb::span<uint64_t const> states, unsigned const bits,
                         unsigned const shift)
    {
        LATTICE_SYMMETRIES_ASSERT(0 < bits && bits <= 32, "invalid bits");
        constexpr auto empty = std::make_pair(~uint64_t{0}, uint64_t{0});
        auto const     size  = uint64_t{1} << bits;

        std::vector<std::pair<uint64_t, uint64_t>> ranges;
        ranges.reserve(size);
        auto const*       first = states.data();
        auto const* const last  = first + states.size();
        auto const* const begin = first;
        for (auto i = uint64_t{0}; i < size; ++i) {
            auto element = empty;
            if (first != last && ((*first) >> shift) == i) {
                element.first = static_cast<uint64_t>(first - begin);
                ++first;
                ++element.second;
                while (first != last && ((*first) >> shift) == i) {
                    ++first;
                    ++element.second;
                }
            }
            ranges.push_back(element);
        }
        LATTICE_SYMMETRIES_CHECK(first == last, "not all states checked");

        // std::unordered_map<uint64_t, uint64_t> histogram;
        // for (auto const& r : ranges) {
        //     ++histogram[r.second];
        // }
        // save_histogram("ranges_histogram.dat", histogram);
        return ranges;
    }

    template <bool FixedHammingWeight> auto next_state(uint64_t const v) noexcept -> uint64_t
    {
        if constexpr (FixedHammingWeight) {
            auto const t = v | (v - 1U); // t gets v's least significant 0 bits set to 1
            // Next set to 1 the most significant bit to change,
            // set to 0 the least significant ones, and add the necessary 1 bits.
            return (t + 1U)
                   | (((~t & -~t) - 1U) >> (static_cast<unsigned>(__builtin_ctzl(v)) + 1U));
        }
        else {
            return v + 1;
        }
    }

    template <bool FixedHammingWeight>
    auto generate_states_task(uint64_t current, uint64_t const upper_bound,
                              basis_base_t const& header, small_basis_t const& payload,
                              std::vector<uint64_t>& states) -> void
    {
        if constexpr (FixedHammingWeight) {
            LATTICE_SYMMETRIES_ASSERT(popcount(current) == popcount(upper_bound),
                                      "current and upper_bound must have the same Hamming weight");
        }
        auto const handle = [&header, &payload, &states](uint64_t const x) noexcept -> void {
            // uint64_t             repr;
            // std::complex<double> character;
            // double               norm;
            // auto const is_repr = is_representative(payload.batched_symmetries, payload.other_symmetries, x);
            auto const is_repr = is_representative(header, payload, x);
            if (is_repr) { states.push_back(x); }
            // get_state_info(batched_symmetries, symmetries, x, repr, character, norm);
            // if (repr == x && norm > 0.0) { states.push_back(x); }
        };

        for (; current < upper_bound; current = next_state<FixedHammingWeight>(current)) {
            handle(current);
        }
        LATTICE_SYMMETRIES_ASSERT(current == upper_bound, "");
        handle(current);
    }

    auto generate_states_task(uint64_t current, uint64_t const upper_bound,
                              basis_base_t const& header, small_basis_t const& payload,
                              std::vector<uint64_t>& states) -> void
    {
        if (header.hamming_weight.has_value()) {
            return generate_states_task<true>(current, upper_bound, header, payload, states);
        }
        return generate_states_task<false>(current, upper_bound, header, payload, states);
    }

    template <bool FixedHammingWeight>
    auto split_into_tasks(uint64_t current, uint64_t const bound, uint64_t chunk_size)
        -> std::vector<std::pair<uint64_t, uint64_t>>
    {
        --chunk_size;
        auto const hamming_weight = popcount(current);
        auto       ranges         = std::vector<std::pair<uint64_t, uint64_t>>{};
        for (;;) {
            if (bound - current <= chunk_size) {
                ranges.emplace_back(current, bound);
                break;
            }
            auto next = FixedHammingWeight ? closest_hamming(current + chunk_size, hamming_weight)
                                           : current + chunk_size;
            LATTICE_SYMMETRIES_ASSERT(next >= current, "");
            if (next >= bound) {
                ranges.emplace_back(current, bound);
                break;
            }
            ranges.emplace_back(current, next);
            current = next_state<FixedHammingWeight>(next);
        }
        return ranges;
    }

    auto get_bounds(unsigned const number_spins, std::optional<unsigned> hamming_weight) noexcept
        -> std::pair<uint64_t, uint64_t>
    {
        if (hamming_weight.has_value()) {
            // special case 0 == all spins down
            if (*hamming_weight == 0U) { return {0U, 0U}; }
            // NOLINTNEXTLINE: special case 64 == all spins up
            if (*hamming_weight == 64U) { return {~uint64_t{0}, ~uint64_t{0}}; }
            auto const current = ~uint64_t{0} >> (64U - *hamming_weight);
            auto const bound   = number_spins > *hamming_weight
                                     ? (current << (number_spins - *hamming_weight))
                                     : current;
            return {current, bound};
        }
        auto const current = uint64_t{0};
        auto const bound   = ~uint64_t{0} >> (64U - number_spins);
        return {current, bound};
    }
} // namespace

LATTICE_SYMMETRIES_EXPORT
auto split_into_tasks(unsigned number_spins, std::optional<unsigned> hamming_weight,
                      uint64_t const chunk_size) -> std::vector<std::pair<uint64_t, uint64_t>>
{
    LATTICE_SYMMETRIES_ASSERT(0 < number_spins && number_spins <= 64, "invalid number of spins");
    LATTICE_SYMMETRIES_ASSERT(!hamming_weight.has_value() || *hamming_weight <= number_spins,
                              "invalid hamming weight");
    auto const [current, bound] = get_bounds(number_spins, hamming_weight);
    if (hamming_weight.has_value()) { return split_into_tasks<true>(current, bound, chunk_size); }
    return split_into_tasks<false>(current, bound, chunk_size);
}

auto closest_hamming(uint64_t x, unsigned const hamming_weight) noexcept -> uint64_t
{
    LATTICE_SYMMETRIES_ASSERT(hamming_weight <= 64, "invalid hamming weight");
    auto weight = popcount(x);
    if (weight > hamming_weight) {
        auto const max_value =
            hamming_weight == 0 ? uint64_t{0} : (~uint64_t{0} << (64U - hamming_weight));
        // Keep clearing lowest bits until we reach the desired Hamming weight.
        for (auto i = 0U; weight > hamming_weight; ++i) {
            LATTICE_SYMMETRIES_ASSERT(i < 64U, "index out of bounds");
            if (test_bit(x, i)) {
                clear_bit(x, i);
                --weight;
            }
        }
        if (x < max_value) { x = next_state<true>(x); }
    }
    else if (weight < hamming_weight) {
        // Keep setting lowest bits until we reach the desired Hamming weight.
        for (auto i = 0U; weight < hamming_weight; ++i) {
            LATTICE_SYMMETRIES_ASSERT(i < 64U, "index out of bounds");
            if (!test_bit(x, i)) {
                set_bit(x, i);
                ++weight;
            }
        }
    }
    return x;
}

namespace {
    auto generate_states(basis_base_t const& header, small_basis_t const& payload)
        -> std::vector<std::vector<uint64_t>>
    {
        LATTICE_SYMMETRIES_CHECK(0 < header.number_spins && header.number_spins <= 64,
                                 "invalid number of spins");
        LATTICE_SYMMETRIES_CHECK(!header.hamming_weight.has_value()
                                     || *header.hamming_weight <= header.number_spins,
                                 "invalid hamming weight");

        auto const chunk_size = [&header]() {
            auto const number_chunks    = 100U * static_cast<unsigned>(omp_get_max_threads());
            auto const [current, bound] = get_bounds(header.number_spins, header.hamming_weight);
            return std::max((bound - current) / number_chunks, uint64_t{1});
        }();
        auto ranges = lattice_symmetries::split_into_tasks(header.number_spins,
                                                           header.hamming_weight, chunk_size);
        auto states = std::vector<std::vector<uint64_t>>(ranges.size());
#pragma omp parallel for schedule(dynamic, 1) default(none) shared(header, payload, ranges, states)
        for (auto i = size_t{0}; i < ranges.size(); ++i) {
            auto const [current, bound] = ranges[i];
            // NOLINTNEXTLINE: 1024 * 1024 == 1048576 == 1MB
            states[i].reserve((1024 * 1024) / sizeof(uint64_t));
            generate_states_task(current, bound, header, payload, states[i]);
        }
        return states;
    }

    auto concatenate(std::vector<std::vector<uint64_t>> const& chunks)
    {
        auto r = std::vector<uint64_t>{};
        r.reserve(std::accumulate(std::begin(chunks), std::end(chunks), size_t{0},
                                  [](auto acc, auto const& x) { return acc + x.size(); }));
        std::for_each(std::begin(chunks), std::end(chunks),
                      [&r](auto const& x) { r.insert(std::end(r), std::begin(x), std::end(x)); });
        return r;
    }
} // namespace

// basis_cache_t::basis_cache_t(tcb::span<batched_small_symmetry_t const> batched,
//                              tcb::span<small_symmetry_t const> other, unsigned number_spins,
//                              std::optional<unsigned> hamming_weight,
//                              std::vector<uint64_t>   _unsafe_states)
//     : _shift{bits >= number_spins ? 0U : (number_spins - bits)}
//     , _states{_unsafe_states.empty()
//                   ? concatenate(generate_states(batched, other, number_spins, hamming_weight))
//                   : std::move(_unsafe_states)}
//     , _ranges{generate_ranges(_states, bits, _shift)}
// {}

basis_cache_t::basis_cache_t(basis_base_t const& header, small_basis_t const& payload,
                             std::vector<uint64_t> _unsafe_states)
    : _shift{bits >= header.number_spins ? 0U : (header.number_spins - bits)}
    , _shift_v2{bits_v2 >= header.number_spins ? 0U : (header.number_spins - bits_v2)}
    , _states{_unsafe_states.empty() ? concatenate(generate_states(header, payload))
                                     : std::move(_unsafe_states)}
    , _ranges{} // generate_ranges(_states, bits, _shift)}
    , _ranges_v2{generate_ranges_v2(_states, bits_v2, _shift_v2)}
{}

auto basis_cache_t::states() const noexcept -> tcb::span<uint64_t const> { return _states; }

auto basis_cache_t::number_states() const noexcept -> uint64_t { return _states.size(); }

auto basis_cache_t::index_v2(uint64_t const x) const noexcept -> outcome::result<uint64_t>
{
    using std::begin, std::end;
    auto const  mask  = (uint64_t{1} << bits_v2) - 1;
    auto        shift = _shift_v2;
    auto const* r     = _ranges_v2.data() + (x >> shift);
    while (!r->is_range()) {
        shift -= bits_v2;
        r = _ranges_v2.data() + (-r->start) + ((x >> shift) & mask);
    }
    auto const* first = _states.data() + r->start;
    auto const* last  = first + r->size;
    auto const* i     = search_sorted(first, last, x);
    if (i == last) { return outcome::failure(LS_NOT_A_REPRESENTATIVE); }
    return static_cast<uint64_t>(i - _states.data());
}

auto basis_cache_t::index(uint64_t const x) const noexcept -> outcome::result<uint64_t>
{
    if constexpr (false) {
        using std::begin, std::end;
        auto const& range = _ranges[x >> _shift];
        auto const* first = _states.data() + range.first;
        auto const* last  = first + range.second;
        auto const* i     = search_sorted(first, last, x);
        if (i == last) { return outcome::failure(LS_NOT_A_REPRESENTATIVE); }
        return static_cast<uint64_t>(i - _states.data());
    }
    else {
        return index_v2(x);
    }
}

namespace {
    struct close_file_fn_t {
        // NOLINTNEXTLINE: we're not using GSL, so no gsl::owner
        auto operator()(std::FILE* file) noexcept -> void { std::fclose(file); }
    };

    auto open_file(char const* filename, char const* mode) noexcept
        -> outcome::result<std::unique_ptr<std::FILE, close_file_fn_t>>
    {
        // NOLINTNEXTLINE: we're not using GSL, so no gsl::owner
        auto* p = std::fopen(filename, mode);
        if (p == nullptr) { return LS_COULD_NOT_OPEN_FILE; }
        return std::unique_ptr<std::FILE, close_file_fn_t>{p};
    }

    auto file_size(char const* filename) noexcept -> uint64_t
    {
        struct stat buf; // NOLINT: buf is initialized by stat
        stat(filename, &buf);
        return static_cast<uint64_t>(buf.st_size);
    }
} // namespace

auto save_states(tcb::span<uint64_t const> states, char const* filename) -> outcome::result<void>
{
    constexpr auto                 chunk_size = uint64_t{4096};
    constexpr std::array<char, 16> header     = {42, 42, 42, 42, 42, 42, 42, 42,
                                             42, 42, 42, 42, 42, 42, 42, 42};
    OUTCOME_TRY(stream, open_file(filename, "wb"));
    if (std::fwrite(header.data(), sizeof(char), std::size(header), stream.get())
        != std::size(header)) {
        return LS_FILE_IO_FAILED;
    }
    auto buffer = std::vector<uint64_t>(chunk_size);
    for (auto first = std::begin(states), last = std::end(states); first != last;) {
        auto const  count = std::min(chunk_size, static_cast<uint64_t>(std::distance(first, last)));
        auto const* next  = std::next(first, static_cast<int64_t>(count));
        std::transform(first, next, std::begin(buffer), [](auto const x) { return htole64(x); });
        if (std::fwrite(buffer.data(), sizeof(uint64_t), count, stream.get()) != count) {
            return LS_FILE_IO_FAILED;
        }
        // Move forward
        first = next;
    }
    return LS_SUCCESS;
}

auto load_states(char const* filename) -> outcome::result<std::vector<uint64_t>>
{
    OUTCOME_TRY(stream, open_file(filename, "rb"));

    auto           size        = file_size(filename);
    constexpr auto header_size = 16U;
    if (size < header_size) { return LS_CACHE_IS_CORRUPT; }
    size -= header_size;
    if (size % sizeof(uint64_t) != 0) { return LS_CACHE_IS_CORRUPT; }

    std::array<char, header_size> header; // NOLINT: header is initialized by fread
    if (std::fread(header.data(), sizeof(char), header_size, stream.get()) != std::size(header)) {
        return LS_FILE_IO_FAILED;
    }
    // NOLINTNEXTLINE: 42 is indeed a magic number, that's why it's used here
    if (!std::all_of(std::begin(header), std::end(header), [](auto const c) { return c == 42; })) {
        return LS_CACHE_IS_CORRUPT;
    }
    auto states = std::vector<uint64_t>(size / sizeof(uint64_t));
    if (std::fread(states.data(), sizeof(uint64_t), states.size(), stream.get()) != states.size()) {
        return LS_FILE_IO_FAILED;
    }
    std::transform(std::begin(states), std::end(states), std::begin(states),
                   [](auto const x) { return le64toh(x); });
    return outcome::success(std::move(states));
}

} // namespace lattice_symmetries
