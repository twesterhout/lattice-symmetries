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
#include <cstdint>

namespace lattice_symmetries {

constexpr auto operator==(ls_bits512 const& x, ls_bits512 const& y) noexcept -> bool
{
    return x.words[0] == y.words[0] && x.words[1] == y.words[1] && x.words[2] == y.words[2]
           && x.words[3] == y.words[3] && x.words[4] == y.words[4] && x.words[5] == y.words[5]
           && x.words[6] == y.words[6] && x.words[7] == y.words[7];
}
constexpr auto operator!=(ls_bits512 const& x, ls_bits512 const& y) noexcept -> bool
{
    return !(x == y);
}
constexpr auto operator<(ls_bits512 const& x, ls_bits512 const& y) noexcept -> bool
{
    for (auto i = 0; i < 8; ++i) {
        if (x.words[i] < y.words[i]) { return true; }
        else if (x.words[i] > y.words[i]) {
            return false;
        }
    }
    return false;
}
constexpr auto operator>(ls_bits512 const& x, ls_bits512 const& y) noexcept -> bool
{
    return y < x;
}
constexpr auto operator<=(ls_bits512 const& x, ls_bits512 const& y) noexcept -> bool
{
    return !(x > y);
}
constexpr auto operator>=(ls_bits512 const& x, ls_bits512 const& y) noexcept -> bool
{
    return !(x < y);
}

constexpr auto set_zero(ls_bits64& bits) noexcept -> void { bits = 0UL; }
constexpr auto set_zero(ls_bits512& bits) noexcept -> void
{
    for (auto& w : bits.words) {
        set_zero(w);
    }
}
constexpr auto is_zero(ls_bits64 const bits) noexcept -> bool { return bits == 0; }
constexpr auto is_zero(ls_bits512 const& bits) noexcept -> bool
{
    for (auto const word : bits.words) {
        if (word != 0) { return false; }
    }
    return true;
}

constexpr auto set_bit(ls_bits64& bits, unsigned const i) noexcept -> void
{
    LATTICE_SYMMETRIES_ASSERT(i < 64U, "index out of bounds");
    bits |= uint64_t{1} << uint64_t{i};
}
constexpr auto set_bit(ls_bits512& bits, unsigned const i) noexcept -> void
{
    LATTICE_SYMMETRIES_ASSERT(i < 512U, "index out of bounds");
    set_bit(bits.words[i / 64U], i % 64U);
}

constexpr auto clear_bit(ls_bits64& bits, unsigned const i) noexcept -> void
{
    LATTICE_SYMMETRIES_ASSERT(i < 64U, "index out of bounds");
    bits &= ~(uint64_t{1} << uint64_t{i});
}
constexpr auto clear_bit(ls_bits512& bits, unsigned const i) noexcept -> void
{
    LATTICE_SYMMETRIES_ASSERT(i < 512U, "index out of bounds");
    clear_bit(bits.words[i / 64U], i % 64U);
}

constexpr auto test_bit(ls_bits64 const bits, unsigned const i) noexcept -> bool
{
    LATTICE_SYMMETRIES_ASSERT(i < 64U, "index out of bounds");
    return static_cast<bool>((bits >> i) & 1U);
}
constexpr auto test_bit(ls_bits512 const& bits, unsigned const i) noexcept -> bool
{
    LATTICE_SYMMETRIES_ASSERT(i < 512U, "index out of bounds");
    return test_bit(bits.words[i / 64U], i % 64U);
}

constexpr auto toggle_bit(ls_bits64& bits, unsigned const i) noexcept -> void
{
    LATTICE_SYMMETRIES_ASSERT(i < 64U, "index out of bounds");
    bits ^= uint64_t{1} << uint64_t{i};
}
constexpr auto toggle_bit(ls_bits512& bits, unsigned const i) noexcept -> void
{
    LATTICE_SYMMETRIES_ASSERT(i < 512U, "index out of bounds");
    return toggle_bit(bits.words[i / 64U], i % 64U);
}

constexpr auto set_bit_to(uint64_t& bits, unsigned const i, bool const value) noexcept -> void
{
    clear_bit(bits, i);
    bits |= static_cast<uint64_t>(value) << i;
}

constexpr auto set_bit_to(ls_bits512& bits, unsigned const i, bool const value) noexcept -> void
{
    LATTICE_SYMMETRIES_ASSERT(i < 512U, "index out of bounds");
    set_bit_to(bits.words[i / 64U], i % 64U, value);
}

inline auto popcount(unsigned long const x) noexcept -> unsigned
{
    return static_cast<unsigned>(__builtin_popcountl(x));
}
inline auto popcount(unsigned long long const x) noexcept -> unsigned
{
    return static_cast<unsigned>(__builtin_popcountll(x));
}

} // namespace lattice_symmetries
