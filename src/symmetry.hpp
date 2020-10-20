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

#include "lattice_symmetries/lattice_symmetries.h"
#include "network.hpp"
#include <array>
#include <complex>
#include <variant>

namespace lattice_symmetries {

struct small_symmetry_t {
    small_network_t      network;
    unsigned             sector;
    unsigned             periodicity;
    std::complex<double> eigenvalue;

    small_symmetry_t(fat_benes_network_t const& _network, bool _flip, unsigned _sector,
                     unsigned _periodicity, std::complex<double> _eigenvalue) noexcept
        : network{_network, _flip}
        , sector{_sector}
        , periodicity{_periodicity}
        , eigenvalue{_eigenvalue}
    {}
};

struct big_symmetry_t {
    big_network_t        network;
    unsigned             sector;
    unsigned             periodicity;
    std::complex<double> eigenvalue;

    big_symmetry_t(fat_benes_network_t const& _network, bool _flip, unsigned _sector,
                   unsigned _periodicity, std::complex<double> _eigenvalue) noexcept
        : network{_network, _flip}
        , sector{_sector}
        , periodicity{_periodicity}
        , eigenvalue{_eigenvalue}
    {}
};

struct batched_small_symmetry_t {
    static constexpr auto batch_size = batched_small_network_t::batch_size;

    batched_small_network_t                      network;
    std::array<unsigned, batch_size>             sectors;
    std::array<unsigned, batch_size>             periodicities;
    std::array<std::complex<double>, batch_size> eigenvalues;

    explicit batched_small_symmetry_t(tcb::span<small_symmetry_t const> symmetries);
};

template <class Int> auto compute_periodicity(tcb::span<Int const> permutation) -> unsigned;
auto compute_eigenvalue(unsigned sector, unsigned periodicity) noexcept -> std::complex<double>;

auto get_state_info(tcb::span<batched_small_symmetry_t const> batched_symmetries,
                    tcb::span<small_symmetry_t const> symmetries, uint64_t bits,
                    uint64_t& representative, std::complex<double>& character,
                    double& norm) noexcept -> void;

auto is_representative(tcb::span<batched_small_symmetry_t const> batched_symmetries,
                       tcb::span<small_symmetry_t const> symmetries, uint64_t bits) noexcept
    -> bool;

auto get_state_info(std::vector<big_symmetry_t> const& symmetries, bits512 const& bits,
                    bits512& representative, std::complex<double>& character, double& norm) noexcept
    -> void;

auto is_real(small_symmetry_t const&) noexcept -> bool;
auto is_real(batched_small_symmetry_t const&) noexcept -> bool;
auto is_real(big_symmetry_t const&) noexcept -> bool;

} // namespace lattice_symmetries

struct ls_symmetry {
    std::variant<lattice_symmetries::small_symmetry_t, lattice_symmetries::big_symmetry_t> payload;

    template <class T, class... Args>
    ls_symmetry(std::in_place_type_t<T> tag, Args&&... args) noexcept
        : payload{tag, std::forward<Args>(args)...}
    {}
};
