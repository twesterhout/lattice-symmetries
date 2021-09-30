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
#include <cstdlib>
#include <functional>
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
    std::vector<batched_small_symmetry_t>   batched_symmetries;
    std::optional<batched_small_symmetry_t> other_symmetries;
    unsigned                                number_other_symmetries;
    std::unique_ptr<basis_cache_t>          cache;

    explicit small_basis_t(ls_group const& group);
};

struct big_basis_t {
    std::vector<big_symmetry_t> symmetries;

    explicit big_basis_t(ls_group const& group);
};

auto is_real(ls_spin_basis const& basis) noexcept -> bool;

template <bool CallDestructor = true> struct free_deleter_fn_t {
    template <class T> auto operator()(T* ptr) const noexcept -> void
    {
        if (ptr != nullptr) {
            ptr->~T();
            std::free(ptr);
        }
    }
};

using state_info_kernel_type =
    std::function<void(uint64_t, void const*, void*, std::complex<double>*, double*)>;
using is_representative_kernel_type = std::function<void(uint64_t, void const*, uint8_t*, double*)>;

} // namespace lattice_symmetries

struct ls_flat_spin_basis;

struct ls_flat_group {
    template <class T>
    using buffer_t = std::unique_ptr<T, lattice_symmetries::free_deleter_fn_t<false>>;
    static constexpr uint64_t alignment = 64U;
    friend struct ls_flat_spin_basis;

    std::array<unsigned, 3> shape;
    buffer_t<uint64_t>      masks;
    buffer_t<uint64_t>      shifts;
    buffer_t<double>        eigenvalues_real;
    buffer_t<double>        eigenvalues_imag;
    buffer_t<unsigned>      sectors;
    buffer_t<unsigned>      periodicities;

  private:
    explicit ls_flat_group(std::array<unsigned, 3> _shape) noexcept;
    [[nodiscard]] auto has_allocation_succeeded() const noexcept -> bool;
};

struct ls_flat_spin_basis {
    using atomic_count_t                = lattice_symmetries::atomic_count_t;
    using state_info_kernel_type        = lattice_symmetries::state_info_kernel_type;
    using is_representative_kernel_type = lattice_symmetries::is_representative_kernel_type;
    using unique_ptr_type =
        std::unique_ptr<ls_flat_spin_basis, lattice_symmetries::free_deleter_fn_t<true>>;

    mutable atomic_count_t        refcount;
    unsigned                      number_spins;
    int                           hamming_weight;
    int                           spin_inversion;
    ls_flat_group                 group;
    state_info_kernel_type        state_info_kernel;
    is_representative_kernel_type is_representative_kernel;

  private:
    struct private_tag_type {};

  public:
    // The constructor is public, but unusable from outside of this class which is what we want.
    // (we can't make it private because alloc_aligned() in allocate() needs access to it)
    ls_flat_spin_basis(std::array<unsigned, 3> shape, private_tag_type /*tag*/) noexcept;

    static auto allocate(std::array<unsigned, 3> shape) noexcept -> unique_ptr_type;
};
