// Copyright (c) 2020, Tom Westerhout
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

#include <atomic>
#include <cstdint>

namespace lattice_symmetries {

//
//  boost/detail/atomic_count_std_atomic.hpp
//
//  atomic_count for std::atomic
//
//  Copyright 2013 Peter Dimov
//
//  Distributed under the Boost Software License, Version 1.0.
//  See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt
//
class atomic_count_t {
  public:
    atomic_count_t() noexcept : _value{0} {}
    explicit atomic_count_t(long const v) noexcept : _value{static_cast<std::int_least32_t>(v)} {}
    atomic_count_t(atomic_count_t const&) = delete;
    atomic_count_t(atomic_count_t&&)      = delete;
    auto operator=(atomic_count_t &&) -> atomic_count_t& = delete;
    auto operator=(atomic_count_t const&) -> atomic_count_t& = delete;

    auto operator++() noexcept -> long
    {
        return _value.fetch_add(1, std::memory_order_acq_rel) + 1;
    }
    auto operator--() noexcept -> long
    {
        return _value.fetch_sub(1, std::memory_order_acq_rel) - 1;
    }

    explicit operator long() const noexcept { return _value.load(std::memory_order_acquire); }

  private:
    std::atomic_int_least32_t _value;
};

inline auto load(atomic_count_t const& counter) noexcept -> unsigned
{
    return static_cast<unsigned int>(static_cast<long>(counter));
}

inline auto increment(atomic_count_t& counter) noexcept -> void { ++counter; }
inline auto decrement(atomic_count_t& counter) noexcept -> unsigned
{
    return static_cast<unsigned int>(--counter);
}

} // namespace lattice_symmetries
