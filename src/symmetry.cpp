#include "symmetry.hpp"
#include <algorithm>
#include <cstring>
#include <memory>
#include <numeric>

namespace lattice_symmetries {

template <class Proj>
auto get_projection(tcb::span<small_symmetry_t const> symmetries, Proj proj) noexcept(
    noexcept(std::declval<Proj const&>()(std::declval<small_symmetry_t const&>())))
{
    LATTICE_SYMMETRIES_CHECK(symmetries.size() == batched_small_symmetry_t::batch_size,
                             "symmetries has wrong length");
    using element_type = decltype(std::declval<Proj>()(std::declval<small_symmetry_t const&>()));
    // NOLINTNEXTLINE: batch is initialized by std::transform
    std::array<element_type, batched_small_symmetry_t::batch_size> batch;
    std::transform(std::begin(symmetries), std::end(symmetries), std::begin(batch),
                   std::cref(proj));
    return batch;
}

batched_small_symmetry_t::batched_small_symmetry_t(tcb::span<small_symmetry_t const> symmetries)
    : network{get_projection(symmetries, [](auto const& s) noexcept { return &s.network; })}
    , sectors{get_projection(symmetries, [](auto const& s) noexcept { return s.sector; })}
    , periodicities{get_projection(symmetries,
                                   [](auto const& s) noexcept { return s.periodicity; })}
    , eigenvalues_real{get_projection(symmetries,
                                      [](auto const& s) noexcept { return s.eigenvalue.real(); })}
    , eigenvalues_imag{
          get_projection(symmetries, [](auto const& s) noexcept { return s.eigenvalue.imag(); })}
{}

// \p permutation must be a valid permutation!
template <class Int> auto compute_periodicity(tcb::span<Int const> permutation) -> unsigned
{
    auto const is_id = [](auto const& xs) {
        auto y = Int{0};
        for (auto const x : xs) {
            if (x != y) { return false; }
            ++y;
        }
        return true;
    };
    std::vector<unsigned> a(permutation.size());
    std::vector<unsigned> b(permutation.size());
    std::iota(std::begin(a), std::end(a), Int{0});
    auto periodicity = 0U;
    for (;;) {
        ++periodicity;
        std::transform(std::begin(permutation), std::end(permutation), std::begin(b),
                       [&a](auto const i) { return a[i]; });
        if (is_id(b)) { break; }
        std::swap(a, b);
    }
    return periodicity;
}

template auto compute_periodicity(tcb::span<uint16_t const> permutation) -> unsigned;
template auto compute_periodicity(tcb::span<unsigned const> permutation) -> unsigned;

auto compute_eigenvalue(unsigned sector, unsigned periodicity) noexcept -> std::complex<double>
{
    constexpr auto pi = 3.141592653589793238462643383279502884L;
    if (sector == 0U) { return std::complex<double>{1.0, 0.0}; }
    if (2U * sector == periodicity) { return std::complex<double>{-1.0, 0.0}; }
    auto const arg = -static_cast<long double>(2U * sector) / static_cast<long double>(periodicity);
    auto       re  = static_cast<double>(std::cos(pi * arg));
    auto       im  = static_cast<double>(std::sin(pi * arg));
    constexpr auto cutoff = 1e-9;
    if (std::abs(re) < cutoff) { re = 0.0; }
    if (std::abs(im) < cutoff) { im = 0.0; }
    return {re, im};
}

auto is_real(small_symmetry_t const& symmetry) noexcept -> bool
{
    return symmetry.eigenvalue.imag() == 0.0;
}

auto is_real(batched_small_symmetry_t const& symmetry) noexcept -> bool
{
    return std::all_of(std::begin(symmetry.eigenvalues_imag), std::end(symmetry.eigenvalues_imag),
                       [](auto const& x) { return x == 0.0; });
}

auto is_real(big_symmetry_t const& symmetry) noexcept -> bool
{
    return symmetry.eigenvalue.imag() == 0.0;
}

} // namespace lattice_symmetries

extern "C" {

// cppcheck-suppress unusedFunction
LATTICE_SYMMETRIES_EXPORT ls_error_code ls_create_symmetry(ls_symmetry** ptr, unsigned const length,
                                                           unsigned const permutation[],
                                                           unsigned const sector)
{
    using namespace lattice_symmetries;
    if (ptr == nullptr) { return LS_INVALID_ARGUMENT; }
    auto r = compile(tcb::span{permutation, length});
    if (!r) {
        if (r.assume_error().category() == get_error_category()) {
            return static_cast<ls_error_code>(r.assume_error().value());
        }
        return LS_SYSTEM_ERROR;
    }
    auto periodicity = compute_periodicity(tcb::span{permutation, length});
    if (sector >= periodicity) { return LS_INVALID_SECTOR; }
    auto const eigenvalue = compute_eigenvalue(sector, periodicity);

    // NOLINTNEXTLINE: 64 spins is the max system size which can be represented by uint64_t
    if (length > 64U) {
        auto p = std::make_unique<ls_symmetry>(std::in_place_type_t<big_symmetry_t>{},
                                               std::move(r).assume_value(), sector, periodicity,
                                               eigenvalue);
        *ptr   = p.release();
    }
    else {
        auto p = std::make_unique<ls_symmetry>(std::in_place_type_t<small_symmetry_t>{},
                                               std::move(r).assume_value(), sector, periodicity,
                                               eigenvalue);
        *ptr   = p.release();
    }
    return LS_SUCCESS;
}

// cppcheck-suppress unusedFunction
LATTICE_SYMMETRIES_EXPORT void ls_destroy_symmetry(ls_symmetry* symmetry)
{
    LATTICE_SYMMETRIES_CHECK(symmetry != nullptr, "trying to destroy a nullptr");
    std::default_delete<ls_symmetry>{}(symmetry);
}

LATTICE_SYMMETRIES_EXPORT unsigned ls_get_sector(ls_symmetry const* symmetry)
{
    return std::visit([](auto const& x) noexcept { return x.sector; }, symmetry->payload);
}

LATTICE_SYMMETRIES_EXPORT unsigned ls_get_periodicity(ls_symmetry const* symmetry)
{
    return std::visit([](auto const& x) noexcept { return x.periodicity; }, symmetry->payload);
}

// cppcheck-suppress unusedFunction
LATTICE_SYMMETRIES_EXPORT void ls_get_eigenvalue(ls_symmetry const* symmetry, void* out)
{
    auto const value =
        std::visit([](auto const& x) noexcept -> std::complex<double> { return x.eigenvalue; },
                   symmetry->payload);
    std::memcpy(out, &value, sizeof(std::complex<double>));
}

// cppcheck-suppress unusedFunction
LATTICE_SYMMETRIES_EXPORT double ls_get_phase(ls_symmetry const* symmetry)
{
    return std::visit(
        [](auto const& x) noexcept {
            return static_cast<double>(x.sector) / static_cast<double>(x.periodicity);
        },
        symmetry->payload);
}

LATTICE_SYMMETRIES_EXPORT unsigned ls_symmetry_get_number_spins(ls_symmetry const* symmetry)
{
    return std::visit([](auto const& x) noexcept { return x.network.width; }, symmetry->payload);
}

LATTICE_SYMMETRIES_EXPORT unsigned ls_symmetry_get_network_depth(ls_symmetry const* symmetry)
{
    return std::visit([](auto const& x) noexcept { return x.network.depth; }, symmetry->payload);
}

} // extern "C"

namespace lattice_symmetries {
struct symmetry_get_network_masks_fn_t {
    void* const    out;
    uint64_t const stride;

    auto operator()(small_symmetry_t const& symmetry) const noexcept -> void
    {
        auto* const p = static_cast<uint64_t*>(out);
        for (auto j = 0U; j < symmetry.network.depth; ++j) {
            p[j * stride] = symmetry.network.masks[j];
        }
    }

    auto operator()(big_symmetry_t const& symmetry) const noexcept -> void
    {
        auto* const p = static_cast<ls_bits512*>(out);
        for (auto j = 0U; j < symmetry.network.depth; ++j) {
            p[j * stride] = symmetry.network.masks[j];
        }
    }
};
} // namespace lattice_symmetries

extern "C" {
LATTICE_SYMMETRIES_EXPORT void ls_symmetry_get_network_masks(ls_symmetry const* symmetry, void* out,
                                                             uint64_t const stride)
{
    lattice_symmetries::symmetry_get_network_masks_fn_t f{out, stride};
    std::visit(f, symmetry->payload);
}
} // extern "C"

namespace lattice_symmetries {
struct symmetry_apply_fn_t {
    ls_bits512* bits;

    auto operator()(small_symmetry_t const& symmetry) const noexcept -> void
    {
        bits->words[0] = symmetry.network(bits->words[0]);
    }

    auto operator()(big_symmetry_t const& symmetry) const noexcept -> void
    {
        symmetry.network(*bits);
    }
};
} // namespace lattice_symmetries

extern "C" {

// cppcheck-suppress unusedFunction
LATTICE_SYMMETRIES_EXPORT void ls_apply_symmetry(ls_symmetry const* symmetry, ls_bits512* bits)
{
    return std::visit(lattice_symmetries::symmetry_apply_fn_t{bits}, symmetry->payload);
}

} // extern "C"
