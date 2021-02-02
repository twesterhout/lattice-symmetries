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

#include "../network.hpp"
#include <immintrin.h>
#include <cstdint>

#define LATTICE_SYMMETRIES_DECLARE()                                                               \
    auto benes_forward_64(uint64_t                       x[batch_size],                            \
                          batched_small_network_t const& network) noexcept->void
#define LATTICE_SYMMETRIES_DECLARE_128()                                                           \
    auto benes_forward_64_direct(__m128i& x0, __m128i& x1, __m128i& x2, __m128i& x3,               \
                                 batched_small_network_t const& network) noexcept->void
#define LATTICE_SYMMETRIES_DECLARE_256()                                                           \
    auto benes_forward_64_direct(__m256i& x0, __m256i& x1,                                         \
                                 batched_small_network_t const& network) noexcept->void

namespace lattice_symmetries {

inline constexpr int batch_size = 8;

auto benes_forward_64(uint64_t& x, small_network_t const& network) noexcept -> void;
LATTICE_SYMMETRIES_DECLARE();

namespace avx2 {
    LATTICE_SYMMETRIES_DECLARE();
    LATTICE_SYMMETRIES_DECLARE_256();
} // namespace avx2

namespace avx {
    LATTICE_SYMMETRIES_DECLARE();
    LATTICE_SYMMETRIES_DECLARE_128();
} // namespace avx

namespace sse4 {
    LATTICE_SYMMETRIES_DECLARE();
    LATTICE_SYMMETRIES_DECLARE_128();
} // namespace sse4

namespace sse2 {
    LATTICE_SYMMETRIES_DECLARE();
    LATTICE_SYMMETRIES_DECLARE_128();
} // namespace sse2

} // namespace lattice_symmetries

#undef LATTICE_SYMMETRIES_DECLARE
#undef LATTICE_SYMMETRIES_DECLARE_128
#undef LATTICE_SYMMETRIES_DECLARE_256
