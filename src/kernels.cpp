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
#include <vectorclass.h>

#if defined(__AVX2__)
// This is one of the few cases when macros really simplify life
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#    define benes_forward_simd benes_forward_avx2
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#    define benes_forward_512_simd benes_forward_512_avx2
#elif defined(__AVX__)
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#    define benes_forward_simd benes_forward_avx
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#    define benes_forward_512_simd benes_forward_512_avx
#elif defined(__SSE2__) || defined(__x86_64__)
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#    define benes_forward_simd benes_forward_sse2
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#    define benes_forward_512_simd benes_forward_512_sse2
#else
#    error "unsupported architecture"
#endif

#if LATTICE_SYMMETRIES_HAS_AVX2()
#    define search_sorted_arch search_sorted_avx2
#elif LATTICE_SYMMETRIES_HAS_AVX()
#    define search_sorted_arch search_sorted_avx
#else
#    define search_sorted_arch search_sorted_sse2
#endif

namespace lattice_symmetries {
namespace vcl = VCL_NAMESPACE;

namespace {
    LATTICE_SYMMETRIES_FORCEINLINE
    auto bit_permute_step(vcl::Vec8uq& x, vcl::Vec8uq const& m, int const d) noexcept -> void
    {
        vcl::Vec8uq y;
        y = x >> d;
        y ^= x;
        y &= m;
        x ^= y;
        y <<= d;
        x ^= y;
    }
} // namespace

namespace detail {
    auto benes_forward_simd(void* _x, uint64_t const (*masks)[batch_size], unsigned size,
                            uint16_t const deltas[]) noexcept -> void
    {
        auto&       x = *static_cast<vcl::Vec8uq*>(_x);
        vcl::Vec8uq m;
        for (auto i = 0; i < static_cast<int>(size); ++i) {
            m.load(masks[i]);
            bit_permute_step(x, m, deltas[i]);
        }
    }

    auto benes_forward_simd(uint64_t x[batch_size], uint64_t const (*masks)[batch_size],
                            unsigned size, uint16_t const deltas[]) noexcept -> void
    {
        vcl::Vec8uq x_v;
        x_v.load(x);
        benes_forward_simd(&x_v, masks, size, deltas);
        x_v.store(x);
    }
} // namespace detail

namespace detail {
    auto search_sorted_arch(uint64_t const* first, uint64_t const* last,
                            uint64_t const key) noexcept -> uint64_t const*
    {
        constexpr auto cutoff_size = 32;
        constexpr auto vector_size = vcl::Vec8uq::size();
        static_assert(cutoff_size >= 2 * vector_size);
        auto const key_v         = vcl::Vec8uq{key};
        auto       size          = last - first;
        auto const original_size = size;
        if (LATTICE_SYMMETRIES_UNLIKELY(size <= 0)) { return last; }
        while (size > cutoff_size) {
            auto const half = size / 2;
            auto const x    = first[half];
            if (x == key) { return first + half; }
            if (x < key) {
                first += half + 1;
                size -= half + 1;
                // LATTICE_SYMMETRIES_CHECK(first + vector_size <= last, "");
                auto const mask = vcl::to_bits(vcl::Vec8uq{}.load(first) == key_v);
                if (mask != 0) { return first + __builtin_ctz(mask); }
            }
            else {
                // LATTICE_SYMMETRIES_CHECK(half >= vector_size, "");
                size = half - vector_size;
                // if (first + size + vector_size > last) {
                //     std::printf("first=%p, size=%zi, original_size=%zi, vector_size=%i, last=%p\n", first, size, original_size, vector_size, last);
                // }
                // LATTICE_SYMMETRIES_CHECK(first + size + vector_size <= last, "");
                auto const mask = vcl::to_bits(vcl::Vec8uq{}.load(first + size) == key_v);
                if (mask != 0) { return first + size + __builtin_ctz(mask); }
            }
        }
        for (; size >= vector_size; first += vector_size, size -= vector_size) {
            LATTICE_SYMMETRIES_CHECK(first + vector_size <= last, "");
            auto const mask = vcl::to_bits(vcl::Vec8uq{}.load(first) == key_v);
            if (mask != 0) { return first + __builtin_ctz(mask); }
        }
        if (size > 0) {
            LATTICE_SYMMETRIES_CHECK(first + size <= last, "");
            auto const mask = vcl::to_bits(vcl::Vec8uq{}.load_partial(size, first) == key_v);
            if (mask != 0) { return first + __builtin_ctz(mask); }
        }
        LATTICE_SYMMETRIES_ASSERT(first + size == last, "");
        return last;
    }
} // namespace detail

#if LATTICE_SYMMETRIES_HAS_AVX2()
#    if 0
namespace {
    /// Performs one step of the Butterfly network. It exchanges bits with distance
    /// \p d between them if the corresponding bits in the mask \p m are set.
    inline auto bit_permute_step(__m256i& x0, __m256i& x1, __m256i m0, __m256i m1,
                                 int const d) noexcept -> void
    {
        __m256i y0, y1; // NOLINT: it increases readability in this case :)
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
    auto benes_forward_avx2(__m256i& x0, __m256i& x1, uint64_t const (*masks)[batch_size],
                            unsigned size, uint16_t const deltas[]) noexcept -> void
    {
        __m256i m0, m1; // NOLINT
        for (auto i = 0; i < static_cast<int>(size); ++i) {
            m0 = _mm256_load_si256(reinterpret_cast<__m256i const*>(masks[i]));     // NOLINT
            m1 = _mm256_load_si256(reinterpret_cast<__m256i const*>(masks[i]) + 1); // NOLINT
            bit_permute_step(x0, x1, m0, m1, deltas[i]);
        }
    }

    auto benes_forward_avx2(uint64_t x[batch_size], uint64_t const (*masks)[batch_size],
                            unsigned size, uint16_t const deltas[]) noexcept -> void
    {
        auto x0 = _mm256_load_si256(reinterpret_cast<__m256i const*>(x));     // NOLINT
        auto x1 = _mm256_load_si256(reinterpret_cast<__m256i const*>(x) + 1); // NOLINT
        benes_forward_avx2(x0, x1, masks, size, deltas);
        _mm256_store_si256(reinterpret_cast<__m256i*>(x), x0);     // NOLINT
        _mm256_store_si256(reinterpret_cast<__m256i*>(x) + 1, x1); // NOLINT
    }
} // namespace detail
#    endif

#else // AVX or SSE2
namespace {
#    if 0
    /// Performs one step of the Butterfly network. It exchanges bits with distance
    /// \p d between them if the corresponding bits in the mask \p m are set.
    inline auto bit_permute_step(__m128i& x0, __m128i& x1, __m128i& x2, __m128i& x3, __m128i m0,
                                 __m128i m1, __m128i m2, __m128i m3, int const d) noexcept -> void
    {
        __m128i y0, y1, y2, y3; // NOLINT
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
#    endif

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

namespace detail {
#    if 0
    auto benes_forward_simd(__m128i& x0, __m128i& x1, __m128i& x2, __m128i& x3,
                            uint64_t const (*masks)[batch_size], unsigned size,
                            uint16_t const deltas[]) noexcept -> void
    {
        __m128i m0, m1, m2, m3; // NOLINT
        for (auto i = 0U; i < size; ++i) {
            m0 = _mm_load_si128(reinterpret_cast<__m128i const*>(masks[i]));     // NOLINT
            m1 = _mm_load_si128(reinterpret_cast<__m128i const*>(masks[i]) + 1); // NOLINT
            m2 = _mm_load_si128(reinterpret_cast<__m128i const*>(masks[i]) + 2); // NOLINT
            m3 = _mm_load_si128(reinterpret_cast<__m128i const*>(masks[i]) + 3); // NOLINT
            bit_permute_step(x0, x1, x2, x3, m0, m1, m2, m3, deltas[i]);
        }
    }

    auto benes_forward_simd(uint64_t x[batch_size], uint64_t const (*masks)[batch_size],
                            unsigned size, uint16_t const deltas[]) noexcept -> void
    {
        __m128i x0, x1, x2, x3; // NOLINT
        // Really don't have much choice but to use reinterpret_cast
        x0 = _mm_load_si128(reinterpret_cast<__m128i const*>(x));     // NOLINT
        x1 = _mm_load_si128(reinterpret_cast<__m128i const*>(x) + 1); // NOLINT
        x2 = _mm_load_si128(reinterpret_cast<__m128i const*>(x) + 2); // NOLINT
        x3 = _mm_load_si128(reinterpret_cast<__m128i const*>(x) + 3); // NOLINT
        benes_forward_simd(x0, x1, x2, x3, masks, size, deltas);
        _mm_store_si128(reinterpret_cast<__m128i*>(x), x0);     // NOLINT
        _mm_store_si128(reinterpret_cast<__m128i*>(x) + 1, x1); // NOLINT
        _mm_store_si128(reinterpret_cast<__m128i*>(x) + 2, x2); // NOLINT
        _mm_store_si128(reinterpret_cast<__m128i*>(x) + 3, x3); // NOLINT
    }
#    endif

    auto benes_forward_512_simd(ls_bits512& x, ls_bits512 const masks[], unsigned size,
                                uint16_t const deltas[]) noexcept -> void
    {
        __m128i x0, x1, x2, x3;                                              // NOLINT
        __m128i m0, m1, m2, m3;                                              // NOLINT
        x0 = _mm_loadu_si128(reinterpret_cast<__m128i const*>(x.words));     // NOLINT
        x1 = _mm_loadu_si128(reinterpret_cast<__m128i const*>(x.words) + 1); // NOLINT
        x2 = _mm_loadu_si128(reinterpret_cast<__m128i const*>(x.words) + 2); // NOLINT
        x3 = _mm_loadu_si128(reinterpret_cast<__m128i const*>(x.words) + 3); // NOLINT
        for (auto i = 0U; i < size; ++i) {
            m0 = _mm_loadu_si128(reinterpret_cast<__m128i const*>(masks[i].words));     // NOLINT
            m1 = _mm_loadu_si128(reinterpret_cast<__m128i const*>(masks[i].words) + 1); // NOLINT
            m2 = _mm_loadu_si128(reinterpret_cast<__m128i const*>(masks[i].words) + 2); // NOLINT
            m3 = _mm_loadu_si128(reinterpret_cast<__m128i const*>(masks[i].words) + 3); // NOLINT
            bit_permute_step_512(x0, x1, x2, x3, m0, m1, m2, m3, deltas[i]);
        }
        _mm_storeu_si128(reinterpret_cast<__m128i*>(x.words), x0);     // NOLINT
        _mm_storeu_si128(reinterpret_cast<__m128i*>(x.words) + 1, x1); // NOLINT
        _mm_storeu_si128(reinterpret_cast<__m128i*>(x.words) + 2, x2); // NOLINT
        _mm_storeu_si128(reinterpret_cast<__m128i*>(x.words) + 3, x3); // NOLINT
    }
} // namespace detail
#endif

#if defined(LATTICE_SYMMETRIES_ADD_DISPATCH_CODE)
constexpr auto bit_permute_step(uint64_t const x, uint64_t const m, unsigned const d) noexcept
    -> uint64_t
{
    auto const y = (x ^ (x >> d)) & m;
    return x ^ y ^ (y << d);
}

auto benes_forward(uint64_t& x, uint64_t const masks[], unsigned size,
                   uint16_t const deltas[]) noexcept -> void
{
    for (auto i = 0U; i < size; ++i) {
        x = bit_permute_step(x, masks[i], deltas[i]);
    }
}

auto benes_forward(uint64_t x[batch_size], uint64_t const (*masks)[batch_size], unsigned size,
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

auto benes_forward(ls_bits512& x, ls_bits512 const masks[], unsigned size,
                   uint16_t const deltas[]) noexcept -> void
{
    if (__builtin_cpu_supports("avx")) { detail::benes_forward_512_avx(x, masks, size, deltas); }
    else {
        detail::benes_forward_512_sse2(x, masks, size, deltas);
    }
}

auto search_sorted(uint64_t const* first, uint64_t const* last, uint64_t const key) noexcept
    -> uint64_t const*
{
    uint64_t const* r;
    if constexpr (true) {
        if (__builtin_cpu_supports("avx2")) { r = detail::search_sorted_avx2(first, last, key); }
        else if (__builtin_cpu_supports("avx")) {
            r = detail::search_sorted_avx(first, last, key);
        }
        else {
            r = detail::search_sorted_sse2(first, last, key);
        }
    }
    // else {
    // auto const* r2 = std::lower_bound(first, last, key);
    // if (r2 != last && *r2 != key) { r2 = last; }
    // LATTICE_SYMMETRIES_CHECK(r == r2, "");
    // }
    return r;
}
#endif

} // namespace lattice_symmetries
