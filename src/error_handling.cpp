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
#include <cstring>
#include <memory>

namespace lattice_symmetries {

[[noreturn]] auto check_fail(char const* expr, char const* file, unsigned line,
                             char const* function, char const* msg) noexcept -> void
{
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-vararg,hicpp-vararg)
    std::fprintf(stderr,
                 "\x1b[1m\x1b[97mAssertion failed\x1b[0m at "
                 "%s:%u: %s: \"%s\" evaluated to false",
                 file, line, function, expr);
    if (msg != nullptr && std::strlen(msg) > 0) {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-vararg,hicpp-vararg)
        std::fprintf(stderr, ": \x1b[1m\x1b[97m%s\x1b[0m\n", msg);
    }
    else {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-vararg,hicpp-vararg)
        std::fprintf(stderr, "\n");
    }
    std::fflush(stderr);
    std::abort();
}

[[noreturn]] auto assert_fail(char const* expr, char const* file, unsigned line,
                              char const* function, char const* msg) noexcept -> void
{
    // clang-format off
    auto const* bug_msg =
        "╔═════════════════════════════════════════════════════════════════╗\n"
        "║   Congratulations, you have found a bug in lattice_symmetries!  ║\n"
        "║              Please, be so kind to submit it here               ║\n"
        "║     https://github.com/twesterhout/lattice-symmetries/issues    ║\n"
        "╚═════════════════════════════════════════════════════════════════╝";
    // clang-format on
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-vararg,hicpp-vararg)
    std::fprintf(stderr, "%s\n\n", bug_msg);
    check_fail(expr, file, line, function, msg);
}

auto ls_error_category::name() const noexcept -> const char* { return "ls_error_category"; }

auto ls_error_category::message(int c) const -> std::string
{
    switch (static_cast<ls_error_code>(c)) {
    case LS_SUCCESS: return "no error";
    case LS_OUT_OF_MEMORY: return "failed to allocate memory";
    case LS_INVALID_ARGUMENT: return "argument is invalid";
    case LS_INVALID_HAMMING_WEIGHT:
        return "specified Hamming weight is invalid. This means that Hamming weight is negative "
               "or that it exceeds the number of spins. Note that -1 is the only allowed negative "
               "value for Hamming weight, which means do not restrict magnetization";
    case LS_INVALID_NUMBER_SPINS:
        return "specified number of spins is invalid. This means one of two things: specified "
               "number of spins is zero or exceeds 512; or that basis size is incompatible with "
               "indices used to specify interactions";
    case LS_INVALID_PERMUTATION:
        return "argument is not a valid permutation of {0, 1, ..., N-1}. Does your array "
               "contain duplicates? or did you start counting from 1 instead of 0?";
    case LS_INVALID_SECTOR: return "specified symmetry sector exceeds periodicity of the operator";
    case LS_INVALID_STATE:
        return "spin configuration does not belong to the basis. This typically happens when the "
               "norm of the basis element is zero.";
    case LS_INVALID_DATATYPE: return "unsupported datatype";
    case LS_PERMUTATION_TOO_LONG: return "permutation length exceeds 512";
    case LS_INCOMPATIBLE_SYMMETRIES:
        return "specified group generators are incompatible. This either means that symmetry "
               "sectors do not intersect or that there is no one-dimensional group representation";
    case LS_NOT_A_REPRESENTATIVE:
        return "specified basis configuration is not a basis representative. Did you forget to "
               "call ls_get_state_info to obtain the representative?";
    case LS_WRONG_BASIS_TYPE:
        return "basis has wrong type. This typically happens when you try to use an operation "
               "which only works for small systems. E.g. it is not feasible to compute basis "
               "element index for a system of 100 spins.";
    case LS_CACHE_NOT_BUILT:
        return "operation requires a list of basis representatives which has not been built yet. "
               "Did you forget to call ls_build first?";
    case LS_COULD_NOT_OPEN_FILE:
        return "could not open file. Specified path is most likely non-existent";
    case LS_FILE_IO_FAILED:
        return "file input/output operation failed. Did you run out of disk space?";
    case LS_CACHE_IS_CORRUPT:
        return "file does not contain a list of basis representatives. Did you specify wrong file? "
               "or is the file corrupt?";
    case LS_OPERATOR_IS_COMPLEX:
        return "operator is complex. Are you trying to apply a complex operator to a real vector?";
    case LS_DIMENSION_MISMATCH:
        return "dimension of the operator does not match dimension of the vector";
    case LS_SYSTEM_ERROR:
    default: return "unknown error";
    }
}

// auto ls_error_category::default_error_condition(int c) const noexcept -> std::error_condition
// {
//     switch (static_cast<ls_error_code>(c)) {
//     default: return std::error_condition(c, *this);
//     }
// }

auto get_error_category() noexcept -> ls_error_category const&
{
    static ls_error_category c;
    return c;
}

} // namespace lattice_symmetries

extern "C" LATTICE_SYMMETRIES_EXPORT auto ls_error_to_string(ls_error_code code) -> char const*
{
    using namespace lattice_symmetries;
    auto msg = get_error_category().message(code);
    auto p   = std::make_unique<char[]>(msg.size() + 1);
    std::strncpy(p.get(), msg.c_str(), msg.size() + 1);
    return p.release();
}

extern "C" LATTICE_SYMMETRIES_EXPORT auto ls_destroy_string(char const* message) -> void
{
    // We really do want const_cast here
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-const-cast)
    std::default_delete<char[]>{}(const_cast<char*>(message));
}
