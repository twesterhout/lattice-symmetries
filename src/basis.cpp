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

#include "basis.hpp"
#include "group.hpp"
#include "macros.hpp"
#include <algorithm>

namespace lattice_symmetries {

namespace {
    template <class T> auto extract(tcb::span<ls_symmetry const> symmetries) -> std::vector<T>
    {
        auto r = std::vector<T>{};
        r.reserve(symmetries.size());
        std::transform(std::begin(symmetries), std::end(symmetries), std::back_inserter(r),
                       [](auto const& x) { return std::get<T>(x.payload); });
        return r;
    }

    auto get_number_spins(ls_group const& group) noexcept -> std::optional<unsigned>
    {
        if (group.payload.empty()) { return std::nullopt; }
        return std::visit([](auto const& x) noexcept { return x.network.width; },
                          group.payload.front().payload);
    }
} // namespace

small_basis_t::small_basis_t(ls_group const& group)
    : batched_symmetries{}
    , other_symmetries{extract<small_symmetry_t>(group.payload)}
    , cache{nullptr}
{}

big_basis_t::big_basis_t(ls_group const& group) : symmetries{extract<big_symmetry_t>(group.payload)}
{}

} // namespace lattice_symmetries

using namespace lattice_symmetries;

struct ls_spin_basis {
    basis_base_t                             header;
    std::variant<small_basis_t, big_basis_t> payload;

    template <class T>
    explicit ls_spin_basis(std::in_place_type_t<T> tag, ls_group const& group,
                           unsigned const number_spins, std::optional<unsigned> hamming_weight)
        : header{{}, number_spins, hamming_weight, ls_get_group_size(&group) > 0}
        , payload{tag, group}
    {}

    ~ls_spin_basis()
    {
        LATTICE_SYMMETRIES_CHECK(load(header.refcount) == 0, "there remain references to object");
    }
};

struct ls_states {
    tcb::span<uint64_t const> payload;
    ls_spin_basis*            parent;

    ls_states(tcb::span<uint64_t const> states, ls_spin_basis const* owner)
        : payload{states}, parent{ls_copy_spin_basis(owner)}
    {}

    ~ls_states() { ls_destroy_spin_basis(parent); }
};

extern "C" ls_error_code ls_create_spin_basis(ls_spin_basis** ptr, ls_group const* group,
                                              unsigned const number_spins, int const hamming_weight)
{
    if (number_spins == 0 || number_spins > 512) { return LS_INVALID_NUMBER_SPINS; }
    if (auto n = get_number_spins(*group); n.has_value() && number_spins != *n) {
        return LS_INVALID_NUMBER_SPINS;
    }
    if (!(hamming_weight == -1
          || (0 <= hamming_weight && hamming_weight <= static_cast<int>(number_spins)))) {
        return LS_INVALID_HAMMING_WEIGHT;
    }

    auto need_big = [group, number_spins]() {
        if (ls_get_group_size(group) > 0) {
            return std::holds_alternative<big_symmetry_t>(group->payload.front().payload);
        }
        return number_spins > 64;
    }();
    auto const _hamming_weight =
        hamming_weight == -1 ? std::nullopt : std::optional<unsigned>{hamming_weight};
    auto p = need_big ? std::make_unique<ls_spin_basis>(std::in_place_type_t<big_basis_t>{}, *group,
                                                        number_spins, _hamming_weight)
                      : std::make_unique<ls_spin_basis>(std::in_place_type_t<small_basis_t>{},
                                                        *group, number_spins, _hamming_weight);
    increment(p->header.refcount);
    *ptr = p.release();
    return LS_SUCCESS;
}

extern "C" ls_spin_basis* ls_copy_spin_basis(ls_spin_basis const* basis)
{
    LATTICE_SYMMETRIES_ASSERT(load(basis->header.refcount) > 0,
                              "refcount cannot be increased from zero");
    increment(basis->header.refcount);
    return const_cast<ls_spin_basis*>(basis);
}

extern "C" void ls_destroy_spin_basis(ls_spin_basis* basis)
{
    if (decrement(basis->header.refcount) == 0) { std::default_delete<ls_spin_basis>{}(basis); }
}

extern "C" unsigned ls_get_number_spins(ls_spin_basis const* basis)
{
    return basis->header.number_spins;
}

extern "C" unsigned ls_get_number_bits(ls_spin_basis const* basis)
{
    if (std::holds_alternative<big_basis_t>(basis->payload)) { return 512U; }
    return 64U;
}

extern "C" int ls_get_hamming_weight(ls_spin_basis const* basis)
{
    auto const& m = basis->header.hamming_weight;
    return m.has_value() ? static_cast<int>(*m) : -1;
}

extern "C" bool ls_has_symmetries(ls_spin_basis const* basis)
{
    return basis->header.has_symmetries;
}

extern "C" ls_error_code ls_get_number_states(ls_spin_basis const* basis, uint64_t* out)
{
    auto p = std::get_if<small_basis_t>(&basis->payload);
    if (LATTICE_SYMMETRIES_UNLIKELY(p == nullptr)) { return LS_WRONG_BASIS_TYPE; }
    if (LATTICE_SYMMETRIES_UNLIKELY(p->cache == nullptr)) { return LS_CACHE_NOT_BUILT; }
    *out = p->cache->number_states();
    return LS_SUCCESS;
}

extern "C" ls_error_code ls_get_index(ls_spin_basis const* basis, uint64_t const bits[],
                                      uint64_t* index)
{
    auto p = std::get_if<small_basis_t>(&basis->payload);
    if (LATTICE_SYMMETRIES_UNLIKELY(p == nullptr)) { return LS_WRONG_BASIS_TYPE; }
    if (LATTICE_SYMMETRIES_UNLIKELY(p->cache == nullptr)) { return LS_CACHE_NOT_BUILT; }
    auto r = p->cache->index(bits[0]);
    if (LATTICE_SYMMETRIES_UNLIKELY(!r)) {
        if (r.error().category() == get_error_category()) {
            return static_cast<ls_error_code>(r.error().value());
        }
        return LS_SYSTEM_ERROR;
    }
    *index = r.value();
    return LS_SUCCESS;
}

extern "C" ls_error_code ls_build(ls_spin_basis* basis)
{
    auto p = std::get_if<small_basis_t>(&basis->payload);
    if (p == nullptr) { return LS_WRONG_BASIS_TYPE; }
    if (p->cache == nullptr) {
        p->cache = std::make_unique<basis_cache_t>(p->batched_symmetries, p->other_symmetries,
                                                   basis->header.number_spins,
                                                   basis->header.hamming_weight);
    }
    return LS_SUCCESS;
}

extern "C" void ls_get_state_info(ls_spin_basis* basis, uint64_t const bits[],
                                  uint64_t representative[], void* character, double* norm)
{
    struct visitor_t {
        uint64_t const* const bits;
        uint64_t* const       representative;
        std::complex<double>& character;
        double&               norm;

        auto operator()(small_basis_t const& payload) const noexcept
        {
            get_state_info(payload.batched_symmetries, payload.other_symmetries, *bits,
                           *representative, character, norm);
        }
        auto operator()(big_basis_t const& payload) const noexcept
        {
            get_state_info(payload.symmetries, *reinterpret_cast<bits512 const*>(bits),
                           *reinterpret_cast<bits512*>(representative), character, norm);
        }
    };
    std::visit(
        visitor_t{bits, representative, *reinterpret_cast<std::complex<double>*>(character), *norm},
        basis->payload);
}

extern "C" ls_error_code ls_get_states(ls_states** ptr, ls_spin_basis const* basis)
{
    auto small_basis = std::get_if<small_basis_t>(&basis->payload);
    if (LATTICE_SYMMETRIES_UNLIKELY(small_basis == nullptr)) { return LS_WRONG_BASIS_TYPE; }
    if (LATTICE_SYMMETRIES_UNLIKELY(small_basis->cache == nullptr)) { return LS_CACHE_NOT_BUILT; }
    auto const states = small_basis->cache->states();
    auto       p      = std::make_unique<ls_states>(states, basis);
    *ptr              = p.release();
    return LS_SUCCESS;
}

extern "C" void ls_destroy_states(ls_states* states) { std::default_delete<ls_states>{}(states); }

extern "C" uint64_t const* ls_states_get_data(ls_states const* states)
{
    return states->payload.data();
}

extern "C" uint64_t ls_states_get_size(ls_states const* states) { return states->payload.size(); }
