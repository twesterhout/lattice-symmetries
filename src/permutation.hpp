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

#pragma once

#include "bits.hpp"
#include "error_handling.hpp"
#include <optional>
#include <vector>

#define LATTICE_SYMMETRIES_UNREACHABLE __builtin_unreachable()

namespace lattice_symmetries {

enum class status_t { success, invalid_argument, no_solution };

struct fat_benes_network_t {
    std::vector<bits512> fwd_masks;
    std::vector<bits512> bwd_masks;
    std::vector<int>     fwd_deltas;
    std::vector<int>     bwd_deltas;

    auto optimize() -> void;
    auto operator()(std::vector<int> x) const -> std::vector<int>;
};

auto compile(std::vector<int> const& permutation, int initial_delta = 1)
    -> std::optional<fat_benes_network_t>;

/// Returns true when `xs` is a permutation of `{0, ..., xs.size() - 1}`.
auto is_permutation(std::vector<int> const& xs) -> bool;

template <class Function>
auto for_each_pair(int const size, int const offset, Function fn) noexcept -> status_t
{
    LATTICE_SYMMETRIES_CHECK(size >= 0, "invalid size");
    LATTICE_SYMMETRIES_CHECK(offset > 0, "invalid offset");
    for (auto i = 0;; i += 2 * offset) {
        for (auto j = 0; j < offset; ++j) {
            auto const index = i + j;
            if (index >= size - offset) { return status_t::success; }
            if (status_t c = fn(index, index + offset); c != status_t::success) { return c; }
        }
    }
}

} // namespace lattice_symmetries
