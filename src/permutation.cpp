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

#include "permutation.hpp"
#include "bits.hpp"
#include "error_handling.hpp"
#include "macros.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <tuple>

namespace lattice_symmetries {

template <class Int> auto is_permutation_helper(tcb::span<Int const> xs) -> bool
{
    std::vector<Int> range(xs.size());
    std::iota(std::begin(range), std::end(range), Int{0});
    return std::is_permutation(std::begin(xs), std::end(xs), std::begin(range));
}

LATTICE_SYMMETRIES_EXPORT
auto is_permutation(tcb::span<unsigned const> xs) -> bool { return is_permutation_helper(xs); }
auto is_permutation(tcb::span<uint16_t const> xs) -> bool { return is_permutation_helper(xs); }

enum class status_t { success, invalid_argument, no_solution, swap_impossible };
enum class index_type_t { small, big };

struct solver_t {
  private:
    struct state_t {
        std::vector<uint16_t> source;
        std::vector<uint16_t> target;
        std::vector<uint16_t> inverse_source;
        std::vector<uint16_t> inverse_target;
    };

    struct scratch_t {
        ls_bits512                visited;
        ls_bits512                source_mask;
        ls_bits512                target_mask;
        std::vector<index_type_t> types;
        int                       delta;

        auto reset(tcb::span<std::pair<uint16_t, uint16_t> const> pairs) -> void
        {
            set_zero(visited);
            set_zero(source_mask);
            set_zero(target_mask);
            LATTICE_SYMMETRIES_ASSERT(!pairs.empty(), "");
            LATTICE_SYMMETRIES_ASSERT(pairs[0].first < pairs[0].second, "");
            delta = pairs[0].second - pairs[0].first;
            for (auto const& p : pairs) {
                LATTICE_SYMMETRIES_ASSERT(p.first + delta == p.second, "");
                LATTICE_SYMMETRIES_ASSERT(0 <= p.first && p.second < types.size(), "");
                LATTICE_SYMMETRIES_ASSERT(
                    !test_bit(visited, p.first) && !test_bit(visited, p.second), "");
                types[p.first]  = index_type_t::small;
                types[p.second] = index_type_t::big;
                set_bit(visited, p.first);
                set_bit(visited, p.second);
            }
            set_zero(visited);
        }
    };

    state_t   _info;
    scratch_t _cxt;

  public:
    solver_t(std::vector<uint16_t> _source, std::vector<uint16_t> _target)
        : _info{std::move(_source), std::move(_target), {}, {}}, _cxt{}
    {
        auto const n = _info.source.size();
        _info.inverse_source.resize(n);
        _info.inverse_target.resize(n);
        for (auto i = size_t{0}; i < n; ++i) {
            _info.inverse_source[_info.source[i]] = static_cast<uint16_t>(i);
            _info.inverse_target[_info.target[i]] = static_cast<uint16_t>(i);
        }
        _cxt.types.resize(n);
    }

  private:
    [[nodiscard]] auto size() const noexcept { return static_cast<unsigned>(_info.source.size()); }

    [[nodiscard]] auto get_index_type(uint16_t const index) const noexcept -> index_type_t
    {
        LATTICE_SYMMETRIES_ASSERT(index < size(), "invalid index");
        return _cxt.types[index];
    }
    static constexpr auto other_type(index_type_t const type) noexcept -> index_type_t
    {
        switch (type) {
        case index_type_t::small: return index_type_t::big;
        case index_type_t::big: return index_type_t::small;
        default: LATTICE_SYMMETRIES_UNREACHABLE;
        }
    }
    [[nodiscard]] auto other_index(uint16_t const index, index_type_t const type) const noexcept
        -> uint16_t
    {
        LATTICE_SYMMETRIES_ASSERT(index < size(), "invalid index");
        LATTICE_SYMMETRIES_ASSERT(type == get_index_type(index), "invalid type");
        switch (type) {
        case index_type_t::small: return static_cast<uint16_t>(index + _cxt.delta);
        case index_type_t::big: return static_cast<uint16_t>(index - _cxt.delta);
        default: LATTICE_SYMMETRIES_UNREACHABLE;
        }
    }
    [[nodiscard]] auto other_index(uint16_t const index) const noexcept
    {
        return other_index(index, get_index_type(index));
    }

    [[nodiscard]] constexpr auto already_visited(uint16_t const i) const noexcept
    {
        return test_bit(_cxt.visited, i);
    }
    auto mark_visited(uint16_t const i) noexcept -> void
    {
        set_bit(_cxt.visited, i);
        auto const other = other_index(i);
        if (other < size()) { set_bit(_cxt.visited, other); }
    }
    auto mark_not_visited(uint16_t const i) noexcept -> void
    {
        clear_bit(_cxt.visited, i);
        auto const other = other_index(i);
        if (other < size()) { clear_bit(_cxt.visited, other); }
    }

    [[nodiscard]] static auto _find_in(uint16_t const                  value,
                                       tcb::span<uint16_t const> const inverse) noexcept -> uint16_t
    {
        LATTICE_SYMMETRIES_ASSERT(value < inverse.size(), "");
        return inverse[value];
    }
    [[nodiscard]] auto index_in_source(uint16_t const value) const noexcept
    {
        return _find_in(value, _info.inverse_source);
    }
    [[nodiscard]] auto index_in_target(uint16_t const value) const noexcept
    {
        return _find_in(value, _info.inverse_target);
    }

    [[nodiscard]] auto _swap_in(uint16_t const index, index_type_t const type,
                                tcb::span<uint16_t> const perm,
                                tcb::span<uint16_t> const inverse_perm,
                                ls_bits512&               mask) const noexcept -> status_t
    {
        auto const other = other_index(index, type);
        if (other >= size()) { return status_t::swap_impossible; }
        std::swap(perm[index], perm[other]);
        inverse_perm[perm[index]] = index;
        inverse_perm[perm[other]] = other;
        toggle_bit(mask, std::min(index, other));
        return status_t::success;
    }
    [[nodiscard]] auto swap_in_source(uint16_t const index, index_type_t const type) noexcept
    {
        return _swap_in(index, type, _info.source, _info.inverse_source, _cxt.source_mask);
    }
    [[nodiscard]] auto swap_in_target(uint16_t const index, index_type_t const type) noexcept
    {
        return _swap_in(index, type, _info.target, _info.inverse_target, _cxt.target_mask);
    }

    auto solve_in_source(uint16_t i, bool const should_swap) noexcept -> status_t
    {
        if (already_visited(i)) {
            if (should_swap) { return status_t::no_solution; }
            return status_t::success;
        }

        auto type = get_index_type(i);
        if (should_swap) {
            auto const status = swap_in_source(i, type);
            if (status == status_t::swap_impossible) { return status; }
            LATTICE_SYMMETRIES_ASSERT(status == status_t::success, "");
        }

        mark_visited(i);
        if (!should_swap) {
            i    = other_index(i);
            type = other_type(type);
            if (i >= size()) { return status_t::success; }
        }

        auto const target_index          = index_in_target(_info.source[i]);
        auto const should_swap_in_target = get_index_type(target_index) != type;
        auto const status                = solve_in_target(target_index, should_swap_in_target);
        if (status == status_t::swap_impossible) {
            mark_not_visited(i);
            if (should_swap) {
                auto const c = swap_in_source(i, type);
                LATTICE_SYMMETRIES_ASSERT(c == status_t::success, "");
            }
        }
        return status;
    }

    auto solve_in_target(uint16_t i, bool const should_swap) noexcept -> status_t
    {
        auto type = get_index_type(i);
        if (should_swap) {
            auto const status = swap_in_target(i, type);
            if (status == status_t::swap_impossible) { return status; }
            LATTICE_SYMMETRIES_ASSERT(status == status_t::success, "");
        }
        else {
            i    = other_index(i);
            type = other_type(type);
            if (i >= size()) { return status_t::success; }
        }

        auto const source_index          = index_in_source(_info.target[i]);
        auto const should_swap_in_source = get_index_type(source_index) != type;
        return solve_in_source(source_index, should_swap_in_source);
    }

  public:
    auto solve_stage(tcb::span<std::pair<uint16_t, uint16_t> const> const pairs,
                     ls_bits512& source_mask, ls_bits512& target_mask) noexcept -> status_t
    {
        _cxt.reset(pairs);
        auto backup = _info;
        auto status = status_t::success;
        for (auto const& [i, j] : pairs) {
            if (!already_visited(i)) {
                status = solve_in_source(j, false);
                if (status == status_t::swap_impossible) { status = solve_in_source(j, true); }
            }
            if (status != status_t::success) { break; }
        }
        if (status == status_t::success) {
            source_mask = _cxt.source_mask;
            target_mask = _cxt.target_mask;
        }
        else {
            _info = backup;
        }
        return status;
    }

    auto solve() -> fat_benes_network_t
    {
        std::vector<unsigned> stages;
        stages.reserve(static_cast<size_t>(std::ceil(std::log2(_info.source.size()))));
        for (auto i = 1U; i <= _info.source.size() / 2; i *= 2) {
            stages.push_back(i);
        }
        if (stages.empty()) {
            LATTICE_SYMMETRIES_CHECK(_info.source == _info.target, "");
            return fat_benes_network_t{{}, {}, 0U};
        }

        std::vector<ls_bits512> left;
        std::vector<ls_bits512> right;
        left.reserve(stages.size());
        right.reserve(stages.size());
        for (auto const delta : stages) {
            left.emplace_back();
            right.emplace_back();

            std::vector<std::pair<uint16_t, uint16_t>> pairs;
            pairs.reserve(_info.source.size() / 2);
            for (auto i = 0U; i + delta <= _info.source.size(); i += 2 * delta) {
                for (auto j = 0U; j < delta; ++j) {
                    auto const index = i + j;
                    if (index >= _info.source.size() - delta) { break; }
                    pairs.emplace_back(index, index + delta);
                }
            }

            auto const status = solve_stage(pairs, left.back(), right.back());
            LATTICE_SYMMETRIES_CHECK(status == status_t::success, "");
        }
        LATTICE_SYMMETRIES_CHECK(_info.source == _info.target, "");

        auto stages_right = stages;
        LATTICE_SYMMETRIES_CHECK(is_zero(left.back()), "");
        left.pop_back();
        stages.pop_back();
        left.insert(left.end(), right.rbegin(), right.rend());
        stages.insert(stages.end(), stages_right.rbegin(), stages_right.rend());
        return fat_benes_network_t{std::move(left), std::move(stages),
                                   static_cast<unsigned>(_info.source.size())};
    }
};

inline auto next_pow_of_2(uint64_t const x) noexcept -> uint64_t
{
    return x <= 1U ? uint64_t{1}
                   // NOLINTNEXTLINE: 64 is the number of bits in uint64_t, not a magic constant
                   : uint64_t{1} << static_cast<unsigned>(64 - __builtin_clzl(x - 1U));
}

template <class Int>
auto compile_helper(tcb::span<Int const> const permutation) -> outcome::result<fat_benes_network_t>
{
    if (!is_permutation(permutation)) { return LS_INVALID_PERMUTATION; }
    // NOLINTNEXTLINE: 512 is the number of bits in ls_bits512, not a magic constant
    if (permutation.size() > 512U) { return LS_PERMUTATION_TOO_LONG; }
    if (permutation.empty()) { return fat_benes_network_t{{}, {}, 0U}; }
    auto const working_size = next_pow_of_2(permutation.size());

    std::vector<uint16_t> source(working_size);
    std::iota(std::begin(source), std::end(source), uint16_t{0});
    std::vector<uint16_t> target(working_size);
    auto i = std::copy(std::begin(permutation), std::end(permutation), std::begin(target));
    std::copy(std::next(std::begin(source), static_cast<ptrdiff_t>(permutation.size())),
              std::end(source), i);

    auto network = solver_t{std::move(source), std::move(target)}.solve();
    network.size = static_cast<unsigned>(permutation.size());
    return network;
}

LATTICE_SYMMETRIES_EXPORT
auto compile(tcb::span<unsigned const> const permutation) -> outcome::result<fat_benes_network_t>
{
    return compile_helper(permutation);
}

auto compile(tcb::span<uint16_t const> const permutation) -> outcome::result<fat_benes_network_t>
{
    return compile_helper(permutation);
}

#if 0
auto fat_benes_network_t::optimize() -> void
{
    auto       first = size_t{0};
    auto const last  = masks.size();
    for (; first != last && !is_zero(masks[first]); ++first) {}
    if (first == last) { return; }
    for (auto i = first + 1; i != last; ++i) {
        if (!is_zero(masks[i])) {
            masks[first]  = masks[i];
            deltas[first] = deltas[i];
            ++first;
        }
    }
    masks.resize(first);
    deltas.resize(first);
}

auto fat_benes_network_t::operator()(std::vector<int> x) const -> std::vector<int>
{
    for (auto _i = 0; _i < masks.size(); ++_i) {
        auto const& mask  = masks[_i];
        auto const  delta = deltas[_i];
        for (auto j = 0; j < x.size() - delta; ++j) {
            if (test_bit(mask, j)) { std::swap(x[j], x[j + delta]); }
        }
    }
    return x;
}
#endif

} // namespace lattice_symmetries
