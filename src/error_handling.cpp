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

#include "error_handling.hpp"
#include <cstdio>
#include <cstdlib>

namespace lattice_symmetries {

[[noreturn]] auto assert_fail(char const* expr, char const* file, unsigned line,
                              char const* function, char const* msg) noexcept -> void
{
    // clang-format off
    auto const* bug_msg =
        "╔═════════════════════════════════════════════════════════════════╗\n"
        "║       Congratulations, you have found a bug in lbfgs-cpp!       ║\n"
        "║              Please, be so kind to submit it here               ║\n"
        "║     https://github.com/twesterhout/lattice_symmetries/issues    ║\n"
        "╚═════════════════════════════════════════════════════════════════╝";
    // clang-format on
    std::fprintf(stderr, "%s\n\n", bug_msg);
    std::fprintf(stderr,
                 "\x1b[1m\x1b[97mAssertion failed\x1b[0m at "
                 "%s:%u: %s: \"%s\" evaluated to false",
                 file, line, function, expr);
    if (msg != nullptr) { std::fprintf(stderr, ": \x1b[1m\x1b[97m%s\x1b[0m\n", msg); }
    else {
        std::fprintf(stderr, "\n");
    }
    std::fflush(stderr);
    std::abort();
}

auto ls_error_category::name() const noexcept -> const char* { return "ls_error_category"; }
auto ls_error_category::message(int c) const -> std::string
{
    switch (static_cast<ls_error_code>(c)) {
    default: return "unknown error";
    }
}

auto ls_error_category::default_error_condition(int c) const noexcept -> std::error_condition
{
    switch (static_cast<ls_error_code>(c)) {
    default: return std::error_condition(c, *this);
    }
}

auto get_error_category() noexcept -> ls_error_category const&
{
    static ls_error_category c;
    return c;
}

} // namespace lattice_symmetries
