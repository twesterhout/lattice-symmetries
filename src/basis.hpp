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

#include "intrusive_ptr.hpp"
#include "symmetry.hpp"
#include <memory>
#include <optional>

namespace lattice_symmetries {

struct basis_base_t {
    mutable atomic_count_t  refcount;
    unsigned                number_spins;
    std::optional<unsigned> hamming_weight;
    int                     spin_inversion;
    bool                    has_symmetries;
};

struct basis_cache_t;

struct small_basis_t {
    std::vector<batched_small_symmetry_t> batched_symmetries;
    std::vector<small_symmetry_t>         other_symmetries;
    std::unique_ptr<basis_cache_t>        cache;

    explicit small_basis_t(ls_group const& group);
};

struct big_basis_t {
    std::vector<big_symmetry_t> symmetries;

    explicit big_basis_t(ls_group const& group);
};

auto is_real(ls_spin_basis const& basis) noexcept -> bool;

} // namespace lattice_symmetries
