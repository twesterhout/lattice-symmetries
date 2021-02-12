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

#include "benes_forward_64.hpp"
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
auto bit_permute_step_64(vcl::Vec8uq& x, vcl::Vec8uq const& m, int const d) noexcept -> void
{
    vcl::Vec8uq y;
    y = x >> d;
    y ^= x;
    y &= m;
    x ^= y;
    y <<= d;
    x ^= y;
}

LATTICE_SYMMETRIES_FORCEINLINE
auto benes_forward_64_direct(vcl::Vec8uq& x, batched_small_network_t const& network) noexcept
    -> void
{
    vcl::Vec8uq m;
    for (auto i = 0; i < static_cast<int>(network.depth); ++i) {
        m.load(network.masks[i]);
        bit_permute_step_64(x, m, network.deltas[i]);
    }
}

#if LATTICE_SYMMETRIES_HAS_AVX2()
auto benes_forward_64_direct(__m256i& x0, __m256i& x1,
                             batched_small_network_t const& network) noexcept -> void
{
    auto x = vcl::Vec8uq{vcl::Vec4uq{x0}, vcl::Vec4uq{x1}};
    benes_forward_64_direct(x, network);
    x0 = static_cast<__m256i>(x.get_low());
    x1 = static_cast<__m256i>(x.get_high());
}
#else
auto benes_forward_64_direct(__m128i& x0, __m128i& x1, __m128i& x2, __m128i& x3,
                             batched_small_network_t const& network) noexcept -> void
{
    auto x = vcl::Vec8uq{vcl::Vec4uq{x0, x1}, vcl::Vec4uq{x2, x3}};
    benes_forward_64_direct(x, network);
    x0 = static_cast<__m128i>(x.get_low().get_low());
    x1 = static_cast<__m128i>(x.get_low().get_high());
    x2 = static_cast<__m128i>(x.get_high().get_low());
    x3 = static_cast<__m128i>(x.get_high().get_high());
}
#endif

auto benes_forward_64(uint64_t x[batch_size], batched_small_network_t const& network) noexcept
    -> void
{
    auto y = vcl::Vec8uq{}.load(x);
    benes_forward_64_direct(y, network);
    y.store(x);
}
} // namespace lattice_symmetries::ARCH

#if defined(LATTICE_SYMMETRIES_ADD_DISPATCH_CODE)
namespace lattice_symmetries {
constexpr auto bit_permute_step_64(uint64_t const x, uint64_t const m, unsigned const d) noexcept
    -> uint64_t
{
    auto const y = (x ^ (x >> d)) & m;
    return x ^ y ^ (y << d);
}

auto benes_forward_64(uint64_t& x, small_network_t const& network) noexcept -> void
{
    for (auto i = 0U; i < network.depth; ++i) {
        x = bit_permute_step_64(x, network.masks[i], network.deltas[i]);
    }
}

auto benes_forward_64(uint64_t x[batch_size], batched_small_network_t const& network) noexcept
    -> void
{
    LATTICE_SYMMETRIES_DISPATCH(benes_forward_64, x, network);
}
} // namespace lattice_symmetries
#endif
