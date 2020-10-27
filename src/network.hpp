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
#include "permutation.hpp"
#include <array>

namespace lattice_symmetries {

struct small_network_t {
    static constexpr auto max_depth = 11U;

    uint64_t masks[max_depth];
    uint16_t deltas[max_depth];
    bool     flip;
    uint64_t flip_mask;
    uint16_t depth;
    uint16_t width;

    small_network_t(fat_benes_network_t const& fat, bool invert) noexcept;
    auto operator()(uint64_t bits) const noexcept -> uint64_t;
    auto operator()(tcb::span<unsigned> bits) const noexcept -> void;
};

struct alignas(32) big_network_t {
    static constexpr auto max_depth = 17U;

    bits512  masks[max_depth];
    bits512  flip_mask;
    uint16_t deltas[max_depth];
    bool     flip;
    uint16_t depth;
    uint16_t width;

    big_network_t(fat_benes_network_t const& fat, bool invert) noexcept;
    auto operator()(bits512& bits) const noexcept -> void;
    auto operator()(tcb::span<unsigned> bits) const noexcept -> void;
};

auto reconstruct_permutation(small_network_t const& network) -> std::vector<uint16_t>;
auto reconstruct_permutation(big_network_t const& network) -> std::vector<uint16_t>;

auto compose(small_network_t const& x, small_network_t const& y)
    -> outcome::result<small_network_t>;
auto compose(big_network_t const& x, big_network_t const& y) -> outcome::result<small_network_t>;

// auto operator==(small_network_t const& x, small_network_t const& y) -> bool;
// auto operator!=(small_network_t const& x, small_network_t const& y) -> bool;
// auto operator==(big_network_t const& x, big_network_t const& y) -> bool;
// auto operator!=(big_network_t const& x, big_network_t const& y) -> bool;

struct alignas(32) batched_small_network_t {
    static constexpr auto max_depth  = 11U;
    static constexpr auto batch_size = 8U;

    uint64_t masks[max_depth][batch_size];
    uint16_t deltas[max_depth];
    uint64_t flip_masks[batch_size];
    // bool     flips[batch_size];
    // uint64_t flip_mask;
    uint16_t depth;
    uint16_t width;

    explicit batched_small_network_t(
        std::array<small_network_t const*, batch_size> const& networks) noexcept;
    auto operator()(uint64_t bits[batch_size]) const noexcept -> void;
};

} // namespace lattice_symmetries
