// Copyright (c) 2019-2021, Tom Westerhout
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

#include "search_sorted.hpp"
#include "lattice_symmetries/lattice_symmetries.h"
#include <vectorclass.h>

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
namespace vcl = VCL_NAMESPACE;

LATTICE_SYMMETRIES_FORCEINLINE
auto better_to_bits(vcl::Vec8b const& v) noexcept -> unsigned
{
#if LATTICE_SYMMETRIES_HAS_AVX2()
    __m256i    v0 = v.get_low();
    __m256i    v1 = v.get_high();
    auto const r  = _mm256_movemask_pd(_mm256_castsi256_pd(v0))
                   | (_mm256_movemask_pd(_mm256_castsi256_pd(v1)) << 4);
#else
    __m128i v0 = v.get_low().get_low();
    __m128i v1 = v.get_low().get_high();
    __m128i v2 = v.get_high().get_low();
    __m128i v3 = v.get_high().get_high();
    // clang-format off
    auto const r = _mm_movemask_pd(_mm_castsi128_pd(v0))
           | (_mm_movemask_pd(_mm_castsi128_pd(v1)) << 2)
           | (_mm_movemask_pd(_mm_castsi128_pd(v2)) << 4)
           | (_mm_movemask_pd(_mm_castsi128_pd(v3)) << 6);
    // clang-format on
#endif
    return static_cast<unsigned>(r);
}

LATTICE_SYMMETRIES_FORCEINLINE
auto linear_search(uint64_t const* data, unsigned const size, vcl::Vec8uq const& key_v) noexcept
    -> unsigned
{
    if (LATTICE_SYMMETRIES_UNLIKELY(size == 0)) { return size; }
    constexpr auto vector_size = vcl::Vec8uq::size();
    auto const     count       = size / vector_size;
    auto const     rest        = size % vector_size;
    auto           offset      = 0U;
    for (auto i = 0U; i < count; ++i, offset += vector_size) {
        auto const value_v      = vcl::Vec8uq{}.load(data + offset);
        auto const equal_v      = value_v == key_v;
        auto const bits_equal_v = better_to_bits(equal_v);
        if (bits_equal_v != 0) {
            return offset + static_cast<unsigned>(__builtin_ctz(bits_equal_v));
        }
    }
    if (rest != 0) {
        auto const value_v      = vcl::Vec8uq{}.load_partial(static_cast<int>(rest), data + offset);
        auto const equal_v      = value_v == key_v;
        auto const bits_equal_v = better_to_bits(equal_v);
        if (bits_equal_v != 0) {
            auto const local_index = static_cast<unsigned>(__builtin_ctz(bits_equal_v));
            if (local_index < rest) { return offset + local_index; }
        }
    }
    return size;
}

auto search_sorted(uint64_t const* data, uint64_t size, uint64_t const key) noexcept -> uint64_t
{
    if (LATTICE_SYMMETRIES_UNLIKELY(size == 0)) { return size; }
    constexpr auto binary_search_threshold = 64U;
    constexpr auto bytes_in_cache_line     = 64U;
    constexpr auto words_in_cache_line     = bytes_in_cache_line / 8U;
    constexpr auto cache_line_mask         = (~uint64_t{0}) << 6U;
    static_assert(vcl::Vec8uq::size() == words_in_cache_line);

    // _mm_prefetch(static_cast<void const*>(data + size / 2), _MM_HINT_NTA);
    constexpr auto align_to_cache_line = [](uint64_t const* p) {
        return p;
        // auto const n = reinterpret_cast<uintptr_t>(p);
        // return reinterpret_cast<uint64_t const*>(n & cache_line_mask);
    };
    auto const  key_v         = vcl::Vec8uq{key};
    auto const* original_data = data;
    while (size > binary_search_threshold) {
        auto const* middle       = align_to_cache_line(data + (size / 2));
        auto const  value_v      = vcl::Vec8uq{}.load(middle);
        auto const  equal_v      = value_v == key_v;
        auto const  bits_equal_v = better_to_bits(equal_v);
        if (bits_equal_v != 0) {
            return static_cast<uint64_t>(middle - original_data + __builtin_ctz(bits_equal_v));
        }

        auto const less_than_v      = value_v < key_v;
        auto const bits_less_than_v = better_to_bits(less_than_v);
        if (bits_less_than_v != 0) {
            auto const* new_data = middle + words_in_cache_line;
            size                 = static_cast<uint64_t>((data + size) - new_data);
            data                 = new_data;
            // _mm_prefetch(static_cast<void const*>(data + size / 2), _MM_HINT_NTA);
        }
        else {
            size = static_cast<uint64_t>(middle - data);
            // _mm_prefetch(static_cast<void const*>(data + size / 2), _MM_HINT_NTA);
        }
    }
    data = align_to_cache_line(data);
    return static_cast<uint64_t>(data - original_data) + linear_search(data, size, key_v);
}

} // namespace lattice_symmetries::ARCH

#if defined(LATTICE_SYMMETRIES_ADD_DISPATCH_CODE)
extern "C" {
using func_type = auto (*)(uint64_t const*, uint64_t, uint64_t) noexcept -> uint64_t;

static auto resolve_search_sorted() noexcept -> func_type
{
    using namespace lattice_symmetries;
    if (ls_has_avx2()) { return &avx2::search_sorted; }
    if (ls_has_avx()) { return &avx::search_sorted; }
    if (ls_has_sse4()) { return &sse4::search_sorted; }
    return &sse2::search_sorted;
}

} // extern "C"

namespace lattice_symmetries {

#    if defined(__APPLE__) && __APPLE__
auto search_sorted(uint64_t const* data, uint64_t size, uint64_t key) noexcept -> uint64_t
{
    (*resolve_search_sorted())(data, size, key);
}
#    else
__attribute__((visibility("default"))) __attribute__((ifunc("resolve_search_sorted"))) auto
search_sorted(uint64_t const* data, uint64_t size, uint64_t key) noexcept -> uint64_t;
#    endif

} // namespace lattice_symmetries
#endif
