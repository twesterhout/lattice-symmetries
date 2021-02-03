// Copyright (c) 2021, Tom Westerhout
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
#include <pthread.h>

typedef struct ls_arch_flags {
    bool has_avx2;
    bool has_avx;
    bool has_sse4;
} ls_arch_flags;

static ls_arch_flags  cpu_info              = {false, false, false};
static pthread_once_t cpu_info_once_control = PTHREAD_ONCE_INIT;

static void init_cpu_info(void)
{
    __builtin_cpu_init();
    cpu_info.has_avx2 = __builtin_cpu_supports("avx2") > 0;
    cpu_info.has_avx  = __builtin_cpu_supports("avx") > 0;
    cpu_info.has_sse4 =
        __builtin_cpu_supports("sse4.1") > 0 && __builtin_cpu_supports("sse4.2") > 0;
}

bool ls_has_avx2()
{
    pthread_once(&cpu_info_once_control, &init_cpu_info);
    return cpu_info.has_avx2;
}

bool ls_has_avx()
{
    pthread_once(&cpu_info_once_control, &init_cpu_info);
    return cpu_info.has_avx2;
}

bool ls_has_sse4()
{
    pthread_once(&cpu_info_once_control, &init_cpu_info);
    return cpu_info.has_avx2;
}
