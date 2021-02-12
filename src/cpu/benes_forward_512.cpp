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

#include "benes_forward_512.hpp"

#if LATTICE_SYMMETRIES_HAS_AVX2()
#    define ARCH avx2
#elif LATTICE_SYMMETRIES_HAS_AVX()
#    define ARCH avx
#elif LATTICE_SYMMETRIES_HAS_SSE4()
#    define ARCH sse4
#else
#    define ARCH sse2
#endif

namespace lattice_symmetries::ARCH {

namespace {
    auto bit_permute_step_512(__m128i& x0, __m128i& x1, __m128i& x2, __m128i& x3, __m128i m0,
                              __m128i m1, __m128i m2, __m128i m3, int const d) noexcept -> void
    {
        constexpr auto bits_in_word = 64;

        __m128i y0, y1, y2, y3; // NOLINT
        switch (d) {
        case 256: // NOLINT: 256 == 512 / 2 and is a multiple of 64, thus shifting is simplified
            // y <- (x ^ (x >> d)) & m
            y0 = _mm_xor_si128(x0, x2);
            y1 = _mm_xor_si128(x1, x3);
            y2 = _mm_and_si128(x2, m2);
            y3 = _mm_and_si128(x3, m3);
            y0 = _mm_and_si128(y0, m0);
            y1 = _mm_and_si128(y1, m1);
            // y <- y ^ (y << d)
            y2 = _mm_xor_si128(y0, y2);
            y3 = _mm_xor_si128(y1, y3);
            break;
        case 128: // NOLINT: 128 == 512 / 4 and is a multiple of 64, thus shifting is simplified
            // y <- (x ^ (x >> d)) & m
            y0 = _mm_xor_si128(x0, x1);
            y1 = _mm_xor_si128(x1, x2);
            y2 = _mm_xor_si128(x2, x3);
            y3 = _mm_and_si128(x3, m3);
            y0 = _mm_and_si128(y0, m0);
            y1 = _mm_and_si128(y1, m1);
            y2 = _mm_and_si128(y2, m2);
            // y <- y ^ (y << d)
            y3 = _mm_xor_si128(y3, y2);
            y2 = _mm_xor_si128(y2, y1);
            y1 = _mm_xor_si128(y1, y0);
            break;
        default:
            // y <- (x ^ (x >> d)) & m
            if (d == bits_in_word) { // NOLINT: number of bits in a word
                constexpr auto bytes = bits_in_word / 8;
                y0 = _mm_or_si128(_mm_slli_si128(x1, bytes), _mm_srli_si128(x0, bytes));
                y1 = _mm_or_si128(_mm_slli_si128(x2, bytes), _mm_srli_si128(x1, bytes));
                y2 = _mm_or_si128(_mm_slli_si128(x3, bytes), _mm_srli_si128(x2, bytes));
                y3 = _mm_srli_si128(x3, bytes);
            }
            else {
                LATTICE_SYMMETRIES_ASSERT(d < bits_in_word, "not implemented");
                y0 = _mm_or_si128(_mm_slli_epi64(x1, bits_in_word - d), _mm_srli_epi64(x0, d));
                y1 = _mm_or_si128(_mm_slli_epi64(x2, bits_in_word - d), _mm_srli_epi64(x1, d));
                y2 = _mm_or_si128(_mm_slli_epi64(x3, bits_in_word - d), _mm_srli_epi64(x2, d));
                y3 = _mm_srli_epi64(x3, d);
            }
            y0 = _mm_xor_si128(x0, y0);
            y1 = _mm_xor_si128(x1, y1);
            y2 = _mm_xor_si128(x2, y2);
            y3 = _mm_xor_si128(x3, y3);
            y0 = _mm_and_si128(y0, m0);
            y1 = _mm_and_si128(y1, m1);
            y2 = _mm_and_si128(y2, m2);
            y3 = _mm_and_si128(y3, m3);
            // y <- y ^ (y << d)
            if (d == bits_in_word) {
                constexpr auto bytes = bits_in_word / 8;

                m0 = _mm_slli_si128(y0, bytes);
                m1 = _mm_or_si128(_mm_slli_si128(y1, bytes), _mm_srli_si128(y0, bytes));
                m2 = _mm_or_si128(_mm_slli_si128(y2, bytes), _mm_srli_si128(y1, bytes));
                m3 = _mm_or_si128(_mm_slli_si128(y3, bytes), _mm_srli_si128(y2, bytes));
            }
            else {
                m0 = _mm_slli_epi64(y0, d);
                m1 = _mm_or_si128(_mm_slli_epi64(y1, d), _mm_srli_epi64(y0, bits_in_word - d));
                m2 = _mm_or_si128(_mm_slli_epi64(y2, d), _mm_srli_epi64(y1, bits_in_word - d));
                m3 = _mm_or_si128(_mm_slli_epi64(y3, d), _mm_srli_epi64(y2, bits_in_word - d));
            }
            y0 = _mm_xor_si128(y0, m0);
            y1 = _mm_xor_si128(y1, m1);
            y2 = _mm_xor_si128(y2, m2);
            y3 = _mm_xor_si128(y3, m3);
            break;
        } // end switch
        // x <- x ^ y
        x0 = _mm_xor_si128(x0, y0);
        x1 = _mm_xor_si128(x1, y1);
        x2 = _mm_xor_si128(x2, y2);
        x3 = _mm_xor_si128(x3, y3);
    }
} // namespace

auto benes_forward_512(ls_bits512& x, big_network_t const& network) noexcept -> void
{
    __m128i x0, x1, x2, x3;                                              // NOLINT
    __m128i m0, m1, m2, m3;                                              // NOLINT
    x0 = _mm_loadu_si128(reinterpret_cast<__m128i const*>(x.words));     // NOLINT
    x1 = _mm_loadu_si128(reinterpret_cast<__m128i const*>(x.words) + 1); // NOLINT
    x2 = _mm_loadu_si128(reinterpret_cast<__m128i const*>(x.words) + 2); // NOLINT
    x3 = _mm_loadu_si128(reinterpret_cast<__m128i const*>(x.words) + 3); // NOLINT
    for (auto i = 0U; i < network.depth; ++i) {
        m0 = _mm_loadu_si128(reinterpret_cast<__m128i const*>(network.masks[i].words));
        m1 = _mm_loadu_si128(reinterpret_cast<__m128i const*>(network.masks[i].words) + 1);
        m2 = _mm_loadu_si128(reinterpret_cast<__m128i const*>(network.masks[i].words) + 2);
        m3 = _mm_loadu_si128(reinterpret_cast<__m128i const*>(network.masks[i].words) + 3);
        bit_permute_step_512(x0, x1, x2, x3, m0, m1, m2, m3, network.deltas[i]);
    }
    _mm_storeu_si128(reinterpret_cast<__m128i*>(x.words), x0);     // NOLINT
    _mm_storeu_si128(reinterpret_cast<__m128i*>(x.words) + 1, x1); // NOLINT
    _mm_storeu_si128(reinterpret_cast<__m128i*>(x.words) + 2, x2); // NOLINT
    _mm_storeu_si128(reinterpret_cast<__m128i*>(x.words) + 3, x3); // NOLINT
}
} // namespace lattice_symmetries::ARCH

#if defined(LATTICE_SYMMETRIES_ADD_DISPATCH_CODE)
namespace lattice_symmetries {
auto benes_forward_512(ls_bits512& x, lattice_symmetries::big_network_t const& network) noexcept
    -> void
{
    LATTICE_SYMMETRIES_DISPATCH(benes_forward_512, x, network);
}
} // namespace lattice_symmetries
#endif
