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

#include "error_handling.hpp"
#include "lattice_symmetries/lattice_symmetries.h"
#include <outcome.hpp>
#include <span.hpp>
#include <vector>

namespace lattice_symmetries {

/// Benes network IR (Intermediate Representation).
///
/// This is the product of compilation of permutations.
struct fat_benes_network_t {
    std::vector<ls_bits512> masks;
    std::vector<unsigned>   deltas;
    unsigned                size;
};

/// Returns true when `xs` is a permutation of `{0, ..., xs.size() - 1}`.
auto is_permutation(tcb::span<unsigned const> xs) -> bool;
auto is_permutation(tcb::span<uint16_t const> xs) -> bool;

/// Compiles the \p permutation to a standard Benes network.
auto compile(tcb::span<unsigned const> permutation) -> outcome::result<fat_benes_network_t>;
auto compile(tcb::span<uint16_t const> permutation) -> outcome::result<fat_benes_network_t>;

} // namespace lattice_symmetries
