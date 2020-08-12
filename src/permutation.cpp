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
#include "error_handling.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <tuple>

#include <iostream>

namespace lattice_symmetries {

template <class T> auto operator<<(std::ostream& out, std::vector<T> const& xs) -> std::ostream&
{
    auto const write = [&out](auto const& x) {
        if constexpr (std::is_same_v<T, uint8_t>) { out << static_cast<int>(x); }
        else {
            out << x;
        }
    };
    out << '[';
    if (!xs.empty()) {
        write(xs[0]);
        for (auto i = 1U; i < xs.size(); ++i) {
            out << ", ";
            write(xs[i]);
        }
    }
    out << ']';
    return out;
}

auto is_permutation(std::vector<int> const& xs) -> bool
{
    std::vector<int> range(xs.size());
    std::iota(std::begin(range), std::end(range), 0);
    return std::is_permutation(std::begin(xs), std::end(xs), std::begin(range));
}

struct solver_t {
  private:
    std::vector<int> source;
    std::vector<int> target;
    std::vector<int> inverse_source;
    std::vector<int> inverse_target;

    // scratch space
    bits512 _visited;
    bits512 _source_mask;
    bits512 _target_mask;

    enum class index_type_t { small, big };

  public:
    solver_t(std::vector<int> _source, std::vector<int> _target)
        : source{std::move(_source)}
        , target{std::move(_target)}
        , inverse_source{}
        , inverse_target{}
        , _visited{}
        , _source_mask{}
        , _target_mask{}
    {
        inverse_source.resize(source.size());
        inverse_target.resize(target.size());
        for (auto i = size_t{0}; i < source.size(); ++i) {
            inverse_source[source[i]] = i;
            inverse_target[target[i]] = i;
        }
    }

  private:
    /// Given an index of an element in `source` or `target` determines its type.
    static constexpr auto get_index_type(int index, int const delta) noexcept -> index_type_t
    {
        LATTICE_SYMMETRIES_ASSERT(index >= 0 && delta > 0, "invalid arguments");
        index %= 2 * delta;
        if (index < delta) { return index_type_t::small; }
        return index_type_t::big;
    }
    static constexpr auto other_index(int const index, int const delta) noexcept -> int
    {
        switch (get_index_type(index, delta)) {
        case index_type_t::small: return index + delta;
        case index_type_t::big: return index - delta;
        default: LATTICE_SYMMETRIES_UNREACHABLE;
        }
    }

    constexpr auto already_visited(int const i) noexcept { return test_bit(_visited, i); }
    constexpr auto mark_visited(int const i, int const delta) noexcept
    {
        set_bit(_visited, i);
        set_bit(_visited, other_index(i, delta));
    }

    static auto _find_in(int const value, int const delta, std::vector<int> const& inverse) noexcept
        -> std::tuple<int, index_type_t>
    {
        LATTICE_SYMMETRIES_ASSERT(0 <= value && value < inverse.size(), "");
        auto const index = inverse[static_cast<size_t>(value)];
        return std::make_tuple(index, get_index_type(index, delta));
    }
    auto find_in_source(int v, int d) const noexcept { return _find_in(v, d, inverse_source); }
    auto find_in_target(int v, int d) const noexcept { return _find_in(v, d, inverse_target); }

    static auto _swap_in(int const index, index_type_t const type, int const delta,
                         std::vector<int>& perm, std::vector<int>& inverse_perm, bits512& mask)
        -> void
    {
        auto const other = [&]() {
            switch (type) {
            case index_type_t::small: return index + delta;
            case index_type_t::big: return index - delta;
            default: LATTICE_SYMMETRIES_UNREACHABLE;
            }
        }();
        LATTICE_SYMMETRIES_ASSERT(0 <= index && index < perm.size(), "");
        LATTICE_SYMMETRIES_ASSERT(0 <= other && other < perm.size(), "");
        std::swap(perm[index], perm[other]);
        LATTICE_SYMMETRIES_ASSERT(0 <= perm[index] && perm[index] < inverse_perm.size(), "");
        LATTICE_SYMMETRIES_ASSERT(0 <= perm[other] && perm[other] < inverse_perm.size(), "");
        inverse_perm[perm[index]] = index;
        inverse_perm[perm[other]] = other;
        set_bit(mask, std::min(other, index));
    }
    auto swap_in_source(int const index, index_type_t const type, int const delta) noexcept
    {
        _swap_in(index, type, delta, source, inverse_source, _source_mask);
    }
    auto swap_in_target(int const index, index_type_t const type, int const delta) noexcept
    {
        _swap_in(index, type, delta, target, inverse_target, _target_mask);
    }

    auto solve_in_source(int const i, index_type_t const type, int const delta) noexcept -> status_t
    {
        if (already_visited(i)) {
            if (test_bit(_source_mask, std::min(i, other_index(i, delta)))) {
                return status_t::no_solution;
            }
            return status_t::success;
        }
        mark_visited(i, delta);

        auto const [target_index, target_index_type] = find_in_target(source[i], delta);
        if (target_index_type != type) {
            swap_in_target(target_index, target_index_type, delta);
            return solve_in_target(target_index, target_index_type, delta);
        }
        switch (target_index_type) {
        case index_type_t::small:
            return solve_in_target(target_index + delta, index_type_t::big, delta);
        case index_type_t::big:
            return solve_in_target(target_index - delta, index_type_t::small, delta);
        default: LATTICE_SYMMETRIES_UNREACHABLE;
        }
    }

    auto solve_in_target(int const i, index_type_t const type, int const delta) noexcept -> status_t
    {
        auto const [source_index, source_index_type] = find_in_source(target[i], delta);
        if (source_index_type != type) {
            swap_in_source(source_index, source_index_type, delta);
            return solve_in_source(source_index, source_index_type, delta);
        }
        switch (source_index_type) {
        case index_type_t::small:
            return solve_in_source(source_index + delta, index_type_t::big, delta);
        case index_type_t::big:
            return solve_in_source(source_index - delta, index_type_t::small, delta);
        default: LATTICE_SYMMETRIES_UNREACHABLE;
        }
    }

    auto solve_in_source_final(int const i, int const j) noexcept -> status_t
    {
        if (source[i] == target[i] && source[j] == target[j]) { return status_t::success; }
        if (source[i] == target[j] && source[j] == target[i]) {
            swap_in_source(i, index_type_t::small, j - i);
            return status_t::success;
        }
        return status_t::no_solution;
    }

    auto solve_stage(int const delta, bits512& source_mask, bits512& target_mask) noexcept
        -> status_t
    {
        // std::cout << "delta = " << delta << '\n';
        // std::cout << source << '\n';
        // std::cout << target << '\n';
        set_zero(_visited);
        set_zero(_source_mask);
        set_zero(_target_mask);
        auto const status = for_each_pair(
            static_cast<int>(source.size()), delta, [&](auto const i, auto const /*i + delta*/) {
                auto status = status_t::success;
                if (!already_visited(i)) {
                    status = solve_in_source(i, index_type_t::small, delta);
                }
                return status;
            });
        if (status == status_t::success) {
            source_mask = _source_mask;
            target_mask = _target_mask;
        }
        // std::cout << source << '\n';
        // std::cout << target << '\n';
        // std::cout << "delta = " << delta << " - done " << '\n';
        return status;
    }

  public:
    auto solve(int start = 1) -> std::optional<fat_benes_network_t>
    {
        std::vector<int> stages;
        stages.reserve(static_cast<size_t>(std::ceil(std::log2(source.size()))));
        for (auto i = start; i <= source.size() / 2; i *= 2) {
            stages.push_back(i);
        }

        std::vector<bits512> left;
        std::vector<bits512> right;
        left.reserve(stages.size());
        right.reserve(stages.size());
        for (auto const delta : stages) {
            left.emplace_back();
            right.emplace_back();
            auto const status = solve_stage(delta, left.back(), right.back());
            if (status != status_t::success) {
                LATTICE_SYMMETRIES_CHECK(status == status_t::no_solution, "weird error");
                return std::nullopt;
            }
        }
        if (source != target) { return std::nullopt; }
        auto stages_right = stages;
        return fat_benes_network_t{std::move(left), std::move(right), std::move(stages),
                                   std::move(stages_right)};
    }
};

auto compile(std::vector<int> const& permutation, int const initial_delta)
    -> std::optional<fat_benes_network_t>
{
    LATTICE_SYMMETRIES_CHECK(is_permutation(permutation), "not a permutation");
    std::vector<int> source(permutation.size());
    std::iota(std::begin(source), std::end(source), 0);
    return solver_t{std::move(source), permutation}.solve(initial_delta);
}

auto fat_benes_network_t::optimize() -> void
{
    auto const remove_zeros = [](auto& masks, auto& stages) {
        auto       first = size_t{0};
        auto const last  = masks.size();
        for (; first != last && !is_zero(masks[first]); ++first) {}
        if (first == last) { return; }
        for (auto i = first + 1; i != last; ++i) {
            if (!is_zero(masks[i])) {
                masks[first]  = masks[i];
                stages[first] = stages[i];
                ++first;
            }
        }
        masks.resize(first);
        stages.resize(first);
    };
    remove_zeros(fwd_masks, fwd_deltas);
    remove_zeros(bwd_masks, bwd_deltas);
}

auto fat_benes_network_t::operator()(std::vector<int> x) const -> std::vector<int>
{
    for (auto _i = 0; _i < fwd_masks.size(); ++_i) {
        auto const& mask  = this->fwd_masks[_i];
        auto const  delta = this->fwd_deltas[_i];
        for_each_pair(x.size(), delta, [&mask, &x](auto const i, auto const j) {
            if (test_bit(mask, i)) { std::swap(x[i], x[j]); }
            return status_t::success;
        });
    }
    for (auto _i = bwd_masks.size(); _i-- > 0;) {
        auto const& mask  = this->bwd_masks[_i];
        auto const  delta = this->bwd_deltas[_i];
        for_each_pair(x.size(), delta, [&mask, &x](auto const i, auto const j) {
            if (test_bit(mask, i)) { std::swap(x[i], x[j]); }
            return status_t::success;
        });
    }
    return x;
}

} // namespace lattice_symmetries
