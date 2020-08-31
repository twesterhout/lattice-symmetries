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

#include "lattice_symmetries/lattice_symmetries.h"
#include "macros.hpp"
#include <outcome.hpp>
#include <string_view>
#include <system_error>

namespace outcome = OUTCOME_V2_NAMESPACE;

namespace std {
template <> struct is_error_code_enum<ls_error_code> : true_type {};
} // namespace std

namespace lattice_symmetries {

[[noreturn]] auto check_fail(char const* expr, char const* file, unsigned line,
                             char const* function, char const* msg) noexcept -> void;

[[noreturn]] auto assert_fail(char const* expr, char const* file, unsigned line,
                              char const* function, char const* msg) noexcept -> void;

// Define a custom error code category derived from std::error_category
class ls_error_category : public std::error_category {
  public:
    [[nodiscard]] auto name() const noexcept -> const char* final;
    [[nodiscard]] auto message(int c) const -> std::string final;
    // [[nodiscard]] auto default_error_condition(int c) const noexcept -> std::error_condition final;
};

auto get_error_category() noexcept -> ls_error_category const&;

} // namespace lattice_symmetries

inline auto make_error_code(ls_error_code const e) noexcept -> std::error_code
{
    return {static_cast<int>(e), lattice_symmetries::get_error_category()};
}

#define LATTICE_SYMMETRIES_CHECK(cond, msg)                                                        \
    (LATTICE_SYMMETRIES_LIKELY(cond)                                                               \
         ? static_cast<void>(0)                                                                    \
         : ::lattice_symmetries::check_fail(#cond, __FILE__, __LINE__,                             \
                                            static_cast<char const*>(__FUNCTION__), msg))

#if !defined(NDEBUG)
#    define LATTICE_SYMMETRIES_ASSERT(cond, msg)                                                   \
        (LATTICE_SYMMETRIES_LIKELY(cond)                                                           \
             ? static_cast<void>(0)                                                                \
             : ::lattice_symmetries::assert_fail(#cond, __FILE__, __LINE__,                        \
                                                 static_cast<char const*>(__FUNCTION__), msg))
#else
#    define LATTICE_SYMMETRIES_ASSERT(cond, msg) static_cast<void>(cond)
#endif
