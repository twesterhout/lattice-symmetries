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

#define LATTICE_SYMMETRIES_UNREACHABLE __builtin_unreachable()
#define LATTICE_SYMMETRIES_LIKELY(x) __builtin_expect(x, 1)
#define LATTICE_SYMMETRIES_UNLIKELY(x) __builtin_expect(x, 0)
#define LATTICE_SYMMETRIES_EXPORT __attribute__((visibility("default")))

#if !defined(LATTICE_SYMMETRIES_FORCEINLINE)
#    if defined(_MSC_VER)
#        define LATTICE_SYMMETRIES_FORCEINLINE __forceinline
#    elif defined(__GNUC__) && __GNUC__ > 3 // Clang also defines __GNUC__ (as 4)
#        define LATTICE_SYMMETRIES_FORCEINLINE inline __attribute__((__always_inline__))
#    else
#        define LATTICE_SYMMETRIES_FORCEINLINE inline
#    endif
#endif

#if !defined(LATTICE_SYMMETRIES_NOINLINE)
#    if defined(_MSC_VER)
#        define LATTICE_SYMMETRIES_NOINLINE __declspec(noinline)
#    elif defined(__GNUC__) && __GNUC__ > 3 // Clang also defines __GNUC__ (as 4)
#        if defined(__CUDACC__)             // nvcc doesn't always parse __noinline__,
#            define LATTICE_SYMMETRIES_NOINLINE __attribute__((noinline))
#        else
#            define LATTICE_SYMMETRIES_NOINLINE __attribute__((__noinline__))
#        endif
#    else
#        define LATTICE_SYMMETRIES_NOINLINE
#    endif
#endif

#if defined(__AVX2__)
#    define LATTICE_SYMMETRIES_HAS_AVX2() 1
#    define LATTICE_SYMMETRIES_HAS_AVX() 1
#elif defined(__AVX__)
#    define LATTICE_SYMMETRIES_HAS_AVX2() 0
#    define LATTICE_SYMMETRIES_HAS_AVX() 1
#elif defined(__SSE2__) || defined(__x86_64__)
#    define LATTICE_SYMMETRIES_HAS_AVX2() 0
#    define LATTICE_SYMMETRIES_HAS_AVX() 0
#else
#    error "unsupported architecture; lattice-symmetries currently only works on x86_64"
#endif

#if LATTICE_SYMMETRIES_HAS_AVX2()
#    define vcl vcl_avx2
#elif LATTICE_SYMMETRIES_HAS_AVX()
#    define vcl vcl_avx
#else
#    define vcl vcl_sse2
#endif
