// Copyright (c) 2020-2021, Tom Westerhout
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

#include "lattice_symmetries/lattice_symmetries.h"
#include <math.h>
#include <pthread.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h> // for gettimeofday
#include <time.h>

LATTICE_SYMMETRIES_NORETURN
static void default_check_fail_handler(char const* const expr, char const* const file,
                                       unsigned const line, char const* const function,
                                       char const* const msg)
{
    // clang-format off
    char const* const bug_msg =
        "╔═════════════════════════════════════════════════════════════════╗\n"
        "║   Congratulations, you have found a bug in lattice_symmetries!  ║\n"
        "║              Please, be so kind to submit it here               ║\n"
        "║     https://github.com/twesterhout/lattice-symmetries/issues    ║\n"
        "╚═════════════════════════════════════════════════════════════════╝";
    // clang-format on
    fprintf(stderr, "%s\n\n", bug_msg);
    fprintf(stderr,
            "\x1b[1m\x1b[97mAssertion failed\x1b[0m at "
            "%s:%u: %s: \"%s\" evaluated to false",
            file, line, function, expr);
    // msg is NULL-terminated
    if (msg != NULL && strlen(msg) > 0) { // Flawfinder: ignore
        fprintf(stderr, ": \x1b[1m\x1b[97m%s\x1b[0m\n", msg);
    }
    else {
        fprintf(stderr, "\n");
    }
    fflush(stderr);
    abort();
}

static ls_error_handler assert_fail_handler       = &default_check_fail_handler;
static pthread_mutex_t  assert_fail_handler_mutex = PTHREAD_MUTEX_INITIALIZER;
static ls_error_handler check_fail_handler        = &default_check_fail_handler;
static pthread_mutex_t  check_fail_handler_mutex  = PTHREAD_MUTEX_INITIALIZER;
static bool             is_logging_enabled        = false;

LATTICE_SYMMETRIES_EXPORT void ls_set_assert_fail_handler(ls_error_handler const func)
{
    pthread_mutex_lock(&assert_fail_handler_mutex);
    assert_fail_handler = func != NULL ? func : &default_check_fail_handler;
    pthread_mutex_unlock(&assert_fail_handler_mutex);
}

LATTICE_SYMMETRIES_EXPORT void ls_set_check_fail_handler(ls_error_handler const func)
{
    pthread_mutex_lock(&check_fail_handler_mutex);
    check_fail_handler = func != NULL ? func : &default_check_fail_handler;
    pthread_mutex_unlock(&check_fail_handler_mutex);
}

// cppcheck-suppress unusedFunction
LATTICE_SYMMETRIES_EXPORT void ls_assert_fail(char const* const expr, char const* const file,
                                              unsigned const line, char const* const function,
                                              char const* const msg)
{
    (*assert_fail_handler)(expr, file, line, function, msg);
    abort();
}

LATTICE_SYMMETRIES_EXPORT
void ls_check_fail(char const* const expr, char const* const file, unsigned const line,
                   char const* const function, char const* const msg)
{
    (*check_fail_handler)(expr, file, line, function, msg);
    abort();
}

LATTICE_SYMMETRIES_EXPORT char const* ls_error_to_string(ls_error_code code)
{
    switch (code) {
    case LS_SUCCESS: return "no error";
    case LS_OUT_OF_MEMORY: return "failed to allocate memory";
    case LS_INVALID_ARGUMENT: return "argument is invalid";
    case LS_INVALID_HAMMING_WEIGHT:
        return "specified Hamming weight is invalid. This means that Hamming weight is negative "
               "or that it exceeds the number of spins. Note that -1 is the only allowed negative "
               "value for Hamming weight, it means 'do not restrict magnetization'";
    case LS_INVALID_SPIN_INVERSION:
        return "specified spin_inversion is invalid. Expected either -1, 0, or 1";
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

LATTICE_SYMMETRIES_EXPORT void ls_destroy_string(char const* message) {}

LATTICE_SYMMETRIES_EXPORT bool ls_is_logging_enabled()
{
    return __atomic_load_n(&is_logging_enabled, __ATOMIC_ACQUIRE);
}

LATTICE_SYMMETRIES_EXPORT void ls_enable_logging()
{
    return __atomic_store_n(&is_logging_enabled, true, __ATOMIC_RELEASE);
}

LATTICE_SYMMETRIES_EXPORT void ls_disable_logging()
{
    return __atomic_store_n(&is_logging_enabled, false, __ATOMIC_RELEASE);
}

static void print_current_time(FILE* out)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    int milliseconds = rint(tv.tv_usec / 1000.0);
    if (milliseconds >= 1000) {
        milliseconds -= 1000;
        ++tv.tv_sec;
    }
    char       time_buffer[16];
    struct tm* time_info = localtime(&tv.tv_sec);
    strftime(time_buffer, 16, "%H:%M:%S", time_info);
    fprintf(out, "%s.%03d", time_buffer, milliseconds);
}

LATTICE_SYMMETRIES_EXPORT void ls_private_log_debug(char const* file, unsigned const line,
                                                    char const* const function,
                                                    char const* const fmt, ...)
{
    if (ls_is_logging_enabled()) {
        fprintf(stderr, "\x1b[1m\x1b[97m[Debug]\x1b[0m [");
        print_current_time(stderr);
        fprintf(stderr, "] [%s:%u:%s] ", file, line, function);
        va_list args;
        va_start(args, fmt);
        vfprintf(stderr, fmt, args);
        va_end(args);
    }
}
