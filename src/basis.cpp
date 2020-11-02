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

    auto split_into_batches(tcb::span<small_symmetry_t const> symmetries)
        -> std::tuple<std::vector<batched_small_symmetry_t>, std::vector<small_symmetry_t>>
    {
        constexpr auto batch_size = batched_small_symmetry_t::batch_size;
        auto           offset     = 0UL;

        std::vector<batched_small_symmetry_t> batched;
        for (; offset + batch_size < symmetries.size(); offset += batch_size) {
            batched.emplace_back(symmetries.subspan(offset, batch_size));
        }

        std::vector<small_symmetry_t> other;
        std::copy(std::next(std::begin(symmetries), static_cast<ptrdiff_t>(offset)),
                  std::end(symmetries), std::back_inserter(other));
        return std::make_tuple(std::move(batched), std::move(other));
    }

    auto get_number_spins(ls_group const& group) noexcept -> std::optional<unsigned>
    {
        if (group.payload.empty()) { return std::nullopt; }
        return std::visit([](auto const& x) noexcept { return x.network.width; },
                          group.payload.front().payload);
    }
} // namespace

small_basis_t::small_basis_t(ls_group const& group) : cache{nullptr}
{
    auto symmetries = extract<small_symmetry_t>(group.payload);
    // other_symmetries = symmetries;
    std::tie(batched_symmetries, other_symmetries) = split_into_batches(symmetries);
}

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

    ls_spin_basis(ls_spin_basis const&) = delete;
    ls_spin_basis(ls_spin_basis&&)      = delete;
    auto operator=(ls_spin_basis const&) -> ls_spin_basis& = delete;
    auto operator=(ls_spin_basis&&) -> ls_spin_basis& = delete;

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

    ls_states(ls_states const&) = delete;
    ls_states(ls_states&&)      = delete;
    auto operator=(ls_states const&) -> ls_states& = delete;
    auto operator=(ls_states&&) -> ls_states& = delete;

    ~ls_states() { ls_destroy_spin_basis(parent); }
};

extern "C" LATTICE_SYMMETRIES_EXPORT ls_error_code ls_create_spin_basis(ls_spin_basis** ptr,
                                                                        ls_group const* group,
                                                                        unsigned const number_spins,
                                                                        int const hamming_weight)
{
    // NOLINTNEXTLINE: 512 is the max supported system size (i.e. number of bits in bits512)
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
        return number_spins > 64; // NOLINT: 64 is number of bits in uint64_t
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

extern "C" LATTICE_SYMMETRIES_EXPORT ls_spin_basis* ls_copy_spin_basis(ls_spin_basis const* basis)
{
    LATTICE_SYMMETRIES_ASSERT(load(basis->header.refcount) > 0,
                              "refcount cannot be increased from zero");
    increment(basis->header.refcount);
    // NOLINTNEXTLINE: We really do want const_cast here since the only non-const operation on
    // NOLINTNEXTLINE: ls_spin_basis is ls_build which may be called from on any instance
    return const_cast<ls_spin_basis*>(basis);
}

extern "C" LATTICE_SYMMETRIES_EXPORT void ls_destroy_spin_basis(ls_spin_basis* basis)
{
    if (decrement(basis->header.refcount) == 0) { std::default_delete<ls_spin_basis>{}(basis); }
}

extern "C" LATTICE_SYMMETRIES_EXPORT unsigned ls_get_number_spins(ls_spin_basis const* basis)
{
    return basis->header.number_spins;
}

extern "C" LATTICE_SYMMETRIES_EXPORT unsigned ls_get_number_bits(ls_spin_basis const* basis)
{
    if (std::holds_alternative<big_basis_t>(basis->payload)) {
        return 512U; // NOLINT: number of bits in bits512
    }
    return 64U; // NOLINT: number of bits in uint64_t
}

extern "C" LATTICE_SYMMETRIES_EXPORT int ls_get_hamming_weight(ls_spin_basis const* basis)
{
    auto const& m = basis->header.hamming_weight;
    return m.has_value() ? static_cast<int>(*m) : -1;
}

extern "C" LATTICE_SYMMETRIES_EXPORT bool ls_has_symmetries(ls_spin_basis const* basis)
{
    return basis->header.has_symmetries;
}

extern "C" LATTICE_SYMMETRIES_EXPORT ls_error_code ls_get_number_states(ls_spin_basis const* basis,
                                                                        uint64_t*            out)
{
    auto const* p = std::get_if<small_basis_t>(&basis->payload);
    if (LATTICE_SYMMETRIES_UNLIKELY(p == nullptr)) { return LS_WRONG_BASIS_TYPE; }
    if (LATTICE_SYMMETRIES_UNLIKELY(p->cache == nullptr)) { return LS_CACHE_NOT_BUILT; }
    *out = p->cache->number_states();
    return LS_SUCCESS;
}

extern "C" LATTICE_SYMMETRIES_EXPORT ls_error_code ls_get_index(ls_spin_basis const* basis,
                                                                uint64_t const       bits[],
                                                                uint64_t*            index)
{
    auto const* p = std::get_if<small_basis_t>(&basis->payload);
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

extern "C" LATTICE_SYMMETRIES_EXPORT ls_error_code ls_build(ls_spin_basis* basis)
{
    auto* p = std::get_if<small_basis_t>(&basis->payload);
    if (p == nullptr) { return LS_WRONG_BASIS_TYPE; }
    if (p->cache == nullptr) {
        p->cache = std::make_unique<basis_cache_t>(p->batched_symmetries, p->other_symmetries,
                                                   basis->header.number_spins,
                                                   basis->header.hamming_weight);
    }
    return LS_SUCCESS;
}

extern "C" LATTICE_SYMMETRIES_EXPORT void
ls_get_state_info(ls_spin_basis* basis, uint64_t const bits[], uint64_t representative[],
                  void* character, double* norm) // NOLINT: nope, norm can't be const
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
            // We do need reinterpret_casts here
            get_state_info(payload.symmetries, *reinterpret_cast<bits512 const*>(bits),   // NOLINT
                           *reinterpret_cast<bits512*>(representative), character, norm); // NOLINT
        }
    };
    std::visit(visitor_t{bits, representative,
                         *reinterpret_cast<std::complex<double>*>(character), // NOLINT
                         *norm},
               basis->payload);
}

extern "C" LATTICE_SYMMETRIES_EXPORT ls_error_code ls_get_states(ls_states**          ptr,
                                                                 ls_spin_basis const* basis)
{
    auto const* small_basis = std::get_if<small_basis_t>(&basis->payload);
    if (LATTICE_SYMMETRIES_UNLIKELY(small_basis == nullptr)) { return LS_WRONG_BASIS_TYPE; }
    if (LATTICE_SYMMETRIES_UNLIKELY(small_basis->cache == nullptr)) { return LS_CACHE_NOT_BUILT; }
    auto const states = small_basis->cache->states();
    auto       p      = std::make_unique<ls_states>(states, basis);
    *ptr              = p.release();
    return LS_SUCCESS;
}

extern "C" LATTICE_SYMMETRIES_EXPORT void ls_destroy_states(ls_states* states)
{
    std::default_delete<ls_states>{}(states);
}

extern "C" LATTICE_SYMMETRIES_EXPORT uint64_t const* ls_states_get_data(ls_states const* states)
{
    return states->payload.data();
}

extern "C" LATTICE_SYMMETRIES_EXPORT uint64_t ls_states_get_size(ls_states const* states)
{
    return states->payload.size();
}

extern "C" LATTICE_SYMMETRIES_EXPORT ls_error_code ls_save_cache(ls_spin_basis const* basis,
                                                                 char const*          filename)
{
    auto const* small_basis = std::get_if<small_basis_t>(&basis->payload);
    if (small_basis == nullptr) { return LS_WRONG_BASIS_TYPE; }
    if (small_basis->cache == nullptr) { return LS_CACHE_NOT_BUILT; }
    auto const states = small_basis->cache->states();
    auto const r      = save_states(states, filename);
    if (!r) {
        if (r.error().category() == get_error_category()) {
            return static_cast<ls_error_code>(r.error().value());
        }
        return LS_SYSTEM_ERROR;
    }
    return LS_SUCCESS;
}

extern "C" LATTICE_SYMMETRIES_EXPORT ls_error_code ls_load_cache(ls_spin_basis* basis,
                                                                 char const*    filename)
{
    auto* p = std::get_if<small_basis_t>(&basis->payload);
    if (p == nullptr) { return LS_WRONG_BASIS_TYPE; }
    // Cache already built
    if (p->cache != nullptr) { return LS_SUCCESS; }

    auto&& r = load_states(filename);
    if (!r) {
        if (r.error().category() == get_error_category()) {
            return static_cast<ls_error_code>(r.error().value());
        }
        return LS_SYSTEM_ERROR;
    }
    p->cache = std::make_unique<basis_cache_t>(p->batched_symmetries, p->other_symmetries,
                                               basis->header.number_spins,
                                               basis->header.hamming_weight, r.value());
    return LS_SUCCESS;
}

namespace lattice_symmetries {
auto is_real(ls_spin_basis const& basis) noexcept -> bool
{
    struct visitor_fn_t {
        auto operator()(small_basis_t const& x) const noexcept -> bool
        {
            return std::all_of(std::begin(x.batched_symmetries), std::end(x.batched_symmetries),
                               [](auto const& s) { return is_real(s); })
                   && std::all_of(std::begin(x.other_symmetries), std::end(x.other_symmetries),
                                  [](auto const& s) { return is_real(s); });
        }

        auto operator()(big_basis_t const& x) const noexcept -> bool
        {
            return std::all_of(std::begin(x.symmetries), std::end(x.symmetries),
                               [](auto const& s) { return is_real(s); });
        }
    };
    return std::visit(visitor_fn_t{}, basis.payload);
}
} // namespace lattice_symmetries
