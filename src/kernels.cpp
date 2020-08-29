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

#include "kernels.hpp"
#include <immintrin.h>
#include <iostream>

#if defined(__AVX2__)
#    define benes_forward_simd benes_forward_avx2
#    define benes_forward_512_simd benes_forward_512_avx2
#elif defined(__AVX__)
#    define benes_forward_simd benes_forward_avx
#    define benes_forward_512_simd benes_forward_512_avx
#elif defined(__SSE2__) || defined(__x86_64__)
#    define benes_forward_simd benes_forward_sse2
#    define benes_forward_512_simd benes_forward_512_sse2
#else
#    error "unsupported architecture"
#endif

namespace lattice_symmetries {

#if defined(__AVX2__)
namespace {
    /// Performs one step of the Butterfly network. It exchanges bits with distance
    /// \p d between them if the corresponding bits in the mask \p m are set.
    inline auto bit_permute_step(__m256i& x0, __m256i& x1, __m256i m0, __m256i m1,
                                 int const d) noexcept -> void
    {
        __m256i y0, y1;
        y0 = _mm256_srli_epi64(x0, d);
        y1 = _mm256_srli_epi64(x1, d);
        y0 = _mm256_xor_si256(x0, y0);
        y1 = _mm256_xor_si256(x1, y1);
        y0 = _mm256_and_si256(y0, m0);
        y1 = _mm256_and_si256(y1, m1);
        x0 = _mm256_xor_si256(x0, y0);
        x1 = _mm256_xor_si256(x1, y1);
        y0 = _mm256_slli_epi64(y0, d);
        y1 = _mm256_slli_epi64(y1, d);
        x0 = _mm256_xor_si256(x0, y0);
        x1 = _mm256_xor_si256(x1, y1);
    }
} // namespace

namespace detail {
    auto benes_forward_simd(uint64_t x[8], uint64_t const (*masks)[8], unsigned size,
                            uint16_t const deltas[]) noexcept -> void
    {
        __m256i x0, x1;
        __m256i m0, m1;
        x0 = _mm256_load_si256(reinterpret_cast<__m256i const*>(x));
        x1 = _mm256_load_si256(reinterpret_cast<__m256i const*>(x) + 1);
        for (auto i = 0U; i < size; ++i) {
            m0 = _mm256_load_si256(reinterpret_cast<__m256i const*>(masks[i]));
            m1 = _mm256_load_si256(reinterpret_cast<__m256i const*>(masks[i]) + 1);
            bit_permute_step(x0, x1, m0, m1, deltas[i]);
        }
        _mm256_store_si256(reinterpret_cast<__m256i*>(x), x0);
        _mm256_store_si256(reinterpret_cast<__m256i*>(x) + 1, x1);
    }

} // namespace detail

#else
namespace {
    /// Performs one step of the Butterfly network. It exchanges bits with distance
    /// \p d between them if the corresponding bits in the mask \p m are set.
    inline auto bit_permute_step(__m128i& x0, __m128i& x1, __m128i& x2, __m128i& x3, __m128i m0,
                                 __m128i m1, __m128i m2, __m128i m3, int const d) noexcept -> void
    {
        __m128i y0, y1, y2, y3;
        y0 = _mm_srli_epi64(x0, d);
        y1 = _mm_srli_epi64(x1, d);
        y2 = _mm_srli_epi64(x2, d);
        y3 = _mm_srli_epi64(x3, d);
        y0 = _mm_xor_si128(x0, y0);
        y1 = _mm_xor_si128(x1, y1);
        y2 = _mm_xor_si128(x2, y2);
        y3 = _mm_xor_si128(x3, y3);
        y0 = _mm_and_si128(y0, m0);
        y1 = _mm_and_si128(y1, m1);
        y2 = _mm_and_si128(y2, m2);
        y3 = _mm_and_si128(y3, m3);
        x0 = _mm_xor_si128(x0, y0);
        x1 = _mm_xor_si128(x1, y1);
        x2 = _mm_xor_si128(x2, y2);
        x3 = _mm_xor_si128(x3, y3);
        y0 = _mm_slli_epi64(y0, d);
        y1 = _mm_slli_epi64(y1, d);
        y2 = _mm_slli_epi64(y2, d);
        y3 = _mm_slli_epi64(y3, d);
        x0 = _mm_xor_si128(x0, y0);
        x1 = _mm_xor_si128(x1, y1);
        x2 = _mm_xor_si128(x2, y2);
        x3 = _mm_xor_si128(x3, y3);
    }

    auto bit_permute_step_512(__m128i& x0, __m128i& x1, __m128i& x2, __m128i& x3, __m128i m0,
                              __m128i m1, __m128i m2, __m128i m3, int const d) noexcept -> void
    {
        __m128i y0, y1, y2, y3;
        switch (d) {
        case 256:
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
        case 128:
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
            if (d == 64) {
                y0 = _mm_or_si128(_mm_slli_si128(x1, 8), _mm_srli_si128(x0, 8));
                y1 = _mm_or_si128(_mm_slli_si128(x2, 8), _mm_srli_si128(x1, 8));
                y2 = _mm_or_si128(_mm_slli_si128(x3, 8), _mm_srli_si128(x2, 8));
                y3 = _mm_srli_si128(x3, 8);
            }
            else {
                LATTICE_SYMMETRIES_ASSERT(d < 64, "not implemented");
                y0 = _mm_or_si128(_mm_slli_epi64(x1, 64 - d), _mm_srli_epi64(x0, d));
                y1 = _mm_or_si128(_mm_slli_epi64(x2, 64 - d), _mm_srli_epi64(x1, d));
                y2 = _mm_or_si128(_mm_slli_epi64(x3, 64 - d), _mm_srli_epi64(x2, d));
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
            if (d == 64) {
                m0 = _mm_slli_si128(y0, 8);
                m1 = _mm_or_si128(_mm_slli_si128(y1, 8), _mm_srli_si128(y0, 8));
                m2 = _mm_or_si128(_mm_slli_si128(y2, 8), _mm_srli_si128(y1, 8));
                m3 = _mm_or_si128(_mm_slli_si128(y3, 8), _mm_srli_si128(y2, 8));
            }
            else {
                m0 = _mm_slli_epi64(y0, d);
                m1 = _mm_or_si128(_mm_slli_epi64(y1, d), _mm_srli_epi64(y0, 64 - d));
                m2 = _mm_or_si128(_mm_slli_epi64(y2, d), _mm_srli_epi64(y1, 64 - d));
                m3 = _mm_or_si128(_mm_slli_epi64(y3, d), _mm_srli_epi64(y2, 64 - d));
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

namespace detail {
    auto benes_forward_simd(uint64_t x[8], uint64_t const (*masks)[8], unsigned size,
                            uint16_t const deltas[]) noexcept -> void
    {
        __m128i x0, x1, x2, x3;
        __m128i m0, m1, m2, m3;
        x0 = _mm_load_si128(reinterpret_cast<__m128i const*>(x));
        x1 = _mm_load_si128(reinterpret_cast<__m128i const*>(x) + 1);
        x2 = _mm_load_si128(reinterpret_cast<__m128i const*>(x) + 2);
        x3 = _mm_load_si128(reinterpret_cast<__m128i const*>(x) + 3);
        for (auto i = 0U; i < size; ++i) {
            m0 = _mm_load_si128(reinterpret_cast<__m128i const*>(masks[i]));
            m1 = _mm_load_si128(reinterpret_cast<__m128i const*>(masks[i]) + 1);
            m2 = _mm_load_si128(reinterpret_cast<__m128i const*>(masks[i]) + 2);
            m3 = _mm_load_si128(reinterpret_cast<__m128i const*>(masks[i]) + 3);
            bit_permute_step(x0, x1, x2, x3, m0, m1, m2, m3, deltas[i]);
        }
        _mm_store_si128(reinterpret_cast<__m128i*>(x), x0);
        _mm_store_si128(reinterpret_cast<__m128i*>(x) + 1, x1);
        _mm_store_si128(reinterpret_cast<__m128i*>(x) + 2, x2);
        _mm_store_si128(reinterpret_cast<__m128i*>(x) + 3, x3);
    }

    auto benes_forward_512_simd(bits512& x, bits512 const masks[], unsigned size,
                                uint16_t const deltas[], bool flip,
                                bits512 const& flip_mask) noexcept -> void
    {
        __m128i x0, x1, x2, x3;
        __m128i m0, m1, m2, m3;
        x0 = _mm_load_si128(reinterpret_cast<__m128i const*>(x.words));
        x1 = _mm_load_si128(reinterpret_cast<__m128i const*>(x.words) + 1);
        x2 = _mm_load_si128(reinterpret_cast<__m128i const*>(x.words) + 2);
        x3 = _mm_load_si128(reinterpret_cast<__m128i const*>(x.words) + 3);
        for (auto i = 0U; i < size; ++i) {
            m0 = _mm_load_si128(reinterpret_cast<__m128i const*>(masks[i].words));
            m1 = _mm_load_si128(reinterpret_cast<__m128i const*>(masks[i].words) + 1);
            m2 = _mm_load_si128(reinterpret_cast<__m128i const*>(masks[i].words) + 2);
            m3 = _mm_load_si128(reinterpret_cast<__m128i const*>(masks[i].words) + 3);
            bit_permute_step_512(x0, x1, x2, x3, m0, m1, m2, m3, deltas[i]);
        }
        if (flip) {
            m0 = _mm_load_si128(reinterpret_cast<__m128i const*>(flip_mask.words));
            x0 = _mm_xor_si128(x0, m0);
            m1 = _mm_load_si128(reinterpret_cast<__m128i const*>(flip_mask.words) + 1);
            x1 = _mm_xor_si128(x1, m1);
            m2 = _mm_load_si128(reinterpret_cast<__m128i const*>(flip_mask.words) + 2);
            x2 = _mm_xor_si128(x2, m2);
            m3 = _mm_load_si128(reinterpret_cast<__m128i const*>(flip_mask.words) + 3);
            x3 = _mm_xor_si128(x3, m3);
        }
        _mm_store_si128(reinterpret_cast<__m128i*>(x.words), x0);
        _mm_store_si128(reinterpret_cast<__m128i*>(x.words) + 1, x1);
        _mm_store_si128(reinterpret_cast<__m128i*>(x.words) + 2, x2);
        _mm_store_si128(reinterpret_cast<__m128i*>(x.words) + 3, x3);
    }
} // namespace detail
#endif

#if defined(LATTICE_SYMMETRIES_ADD_DISPATCH_CODE)
auto benes_forward(uint64_t x[8], uint64_t const (*masks)[8], unsigned size,
                   uint16_t const deltas[]) noexcept -> void
{
    if (__builtin_cpu_supports("avx2")) {
        // ...
        detail::benes_forward_avx2(x, masks, size, deltas);
    }
    else if (__builtin_cpu_supports("avx")) {
        detail::benes_forward_avx(x, masks, size, deltas);
    }
    else {
        detail::benes_forward_sse2(x, masks, size, deltas);
    }
}

auto benes_forward(bits512& x, bits512 const masks[], unsigned size, uint16_t const deltas[],
                   bool flip, bits512 const& flip_mask) noexcept -> void
{
    if (__builtin_cpu_supports("avx")) {
        detail::benes_forward_512_avx(x, masks, size, deltas, flip, flip_mask);
    }
    else {
        detail::benes_forward_512_sse2(x, masks, size, deltas, flip, flip_mask);
    }
}
#endif

} // namespace lattice_symmetries
