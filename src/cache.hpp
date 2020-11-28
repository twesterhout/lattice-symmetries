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

#include "basis.hpp"
#include "symmetry.hpp"
#include <memory>
#include <optional>
#include <vector>

namespace lattice_symmetries {

auto closest_hamming(uint64_t x, unsigned hamming_weight) noexcept -> uint64_t;
auto split_into_tasks(unsigned number_spins, std::optional<unsigned> hamming_weight,
                      uint64_t chunk_size) -> std::vector<std::pair<uint64_t, uint64_t>>;
// auto generate_states(tcb::span<batched_small_symmetry_t const> batched,
//                      tcb::span<small_symmetry_t const> other, unsigned number_spins,
//                      std::optional<unsigned> hamming_weight) -> std::vector<std::vector<uint64_t>>;

struct range_node_t {
    int64_t start;
    int64_t size;

    [[nodiscard]] constexpr auto is_range() const noexcept -> bool { return start >= 0; }
    [[nodiscard]] constexpr auto is_pointer() const noexcept -> bool { return !is_range(); }

    static constexpr auto make_empty() noexcept -> range_node_t
    {
        return range_node_t{std::numeric_limits<int64_t>::max(), int64_t{0}};
    }
};

struct basis_cache_t {
  private:
    static constexpr auto bits    = 16U;
    static constexpr auto bits_v2 = 4U;

    unsigned                                   _shift;
    unsigned                                   _shift_v2;
    std::vector<uint64_t>                      _states;
    std::vector<std::pair<uint64_t, uint64_t>> _ranges;
    std::vector<range_node_t>                  _ranges_v2;

  public:
    // basis_cache_t(tcb::span<batched_small_symmetry_t const> batched,
    //               tcb::span<small_symmetry_t const> other, unsigned number_spins,
    //               std::optional<unsigned> hamming_weight,
    //               std::vector<uint64_t>   _unsafe_states = {});

    basis_cache_t(basis_base_t const& header, small_basis_t const& payload,
                  std::vector<uint64_t> _unsafe_states = {});

    [[nodiscard]] auto states() const noexcept -> tcb::span<uint64_t const>;
    [[nodiscard]] auto number_states() const noexcept -> uint64_t;
    [[nodiscard]] auto index_v2(uint64_t x) const noexcept -> outcome::result<uint64_t>;
    [[nodiscard]] auto index(uint64_t x) const noexcept -> outcome::result<uint64_t>;
};

auto save_states(tcb::span<uint64_t const> states, char const* filename) -> outcome::result<void>;
auto load_states(char const* filename) -> outcome::result<std::vector<uint64_t>>;

} // namespace lattice_symmetries
