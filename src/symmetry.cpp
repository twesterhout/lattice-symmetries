#include "symmetry.hpp"
#include "macros.hpp"
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
    // NOLINTNEXTLINE: batch is initialized in std::transform
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
    , eigenvalues{get_projection(symmetries, [](auto const& s) noexcept { return s.eigenvalue; })}
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

auto get_state_info(tcb::span<batched_small_symmetry_t const> const batched_symmetries,
                    tcb::span<small_symmetry_t const> const symmetries, uint64_t bits,
                    uint64_t& representative, std::complex<double>& character,
                    double& norm) noexcept -> void
{
    if (batched_symmetries.empty() && symmetries.empty()) {
        representative = bits;
        character      = {1.0, 0.0};
        norm           = 1.0;
        return;
    }
    constexpr auto       batch_size = batched_small_symmetry_t::batch_size;
    alignas(32) uint64_t initial[batch_size]; // NOLINT: 32-byte alignment for AVX
    alignas(32) uint64_t buffer[batch_size];  // NOLINT: same
    std::fill(std::begin(initial), std::end(initial), bits);

    auto r = bits;
    auto n = 0.0;
    auto e = std::complex<double>{1.0};

    for (auto const& symmetry : batched_symmetries) {
        std::copy(std::begin(initial), std::end(initial), std::begin(buffer));
        symmetry.network(static_cast<uint64_t*>(buffer));
        for (auto i = 0U; i < batch_size; ++i) {
            if (buffer[i] < r) {
                r = buffer[i];
                e = symmetry.eigenvalues[i];
            }
            else if (buffer[i] == bits) {
                n += symmetry.eigenvalues[i].real();
            }
        }
    }
    for (auto const& symmetry : symmetries) {
        auto const y = symmetry.network(bits);
        if (y < r) {
            r = y;
            e = symmetry.eigenvalue;
        }
        else if (y == bits) {
            n += symmetry.eigenvalue.real();
        }
    }

    // We need to detect the case when norm is not zero, but only because of
    // inaccurate arithmetics
    constexpr auto norm_threshold = 1.0e-5;
    if (std::abs(n) <= norm_threshold) { n = 0.0; }
    LATTICE_SYMMETRIES_ASSERT(n >= 0.0, "");
    auto const group_size = batch_size * batched_symmetries.size() + symmetries.size();
    n                     = std::sqrt(n / static_cast<double>(group_size));

    // Save results
    representative = r;
    character      = e;
    norm           = n;
}

auto is_representative(tcb::span<batched_small_symmetry_t const> const batched_symmetries,
                       tcb::span<small_symmetry_t const> const symmetries, uint64_t bits) noexcept
    -> bool
{
    if (batched_symmetries.empty() && symmetries.empty()) { return true; }
    constexpr auto       batch_size = batched_small_symmetry_t::batch_size;
    alignas(32) uint64_t initial[batch_size]; // NOLINT: 32-byte alignment for AVX
    alignas(32) uint64_t buffer[batch_size];  // NOLINT: same
    std::fill(std::begin(initial), std::end(initial), bits);

    auto r = bits;
    auto n = 0.0;

    for (auto const& symmetry : batched_symmetries) {
        std::copy(std::begin(initial), std::end(initial), std::begin(buffer));
        symmetry.network(static_cast<uint64_t*>(buffer));
        for (auto i = 0U; i < batch_size; ++i) {
            if (buffer[i] < r) { return false; }
            if (buffer[i] == bits) { n += symmetry.eigenvalues[i].real(); }
        }
    }
    for (auto const& symmetry : symmetries) {
        auto const y = symmetry.network(bits);
        if (y < r) { return false; }
        if (y == bits) { n += symmetry.eigenvalue.real(); }
    }

    // We need to detect the case when norm is not zero, but only because of
    // inaccurate arithmetics
    constexpr auto norm_threshold = 1.0e-5;
    if (std::abs(n) <= norm_threshold) { n = 0.0; }
    LATTICE_SYMMETRIES_ASSERT(n >= 0.0, "");
    return n > 0.0;
}

auto get_state_info(std::vector<big_symmetry_t> const& symmetries, bits512 const& bits,
                    bits512& representative, std::complex<double>& character, double& norm) noexcept
    -> void
{
    if (symmetries.empty()) {
        representative = bits;
        character      = {1.0, 0.0};
        norm           = 1.0;
        return;
    }
    bits512 buffer; // NOLINT: buffer is initialized inside the loop before it is used
    auto    r = bits;
    auto    n = 0.0;
    auto    e = std::complex<double>{1.0};

    for (auto const& symmetry : symmetries) {
        buffer = bits;
        symmetry.network(buffer);
        if (buffer < r) {
            r = buffer;
            e = symmetry.eigenvalue;
        }
        else if (buffer == bits) {
            n += symmetry.eigenvalue.real();
        }
    }

    // We need to detect the case when norm is not zero, but only because of
    // inaccurate arithmetics
    constexpr auto norm_threshold = 1.0e-5;
    if (std::abs(n) <= norm_threshold) { n = 0.0; }
    LATTICE_SYMMETRIES_ASSERT(n >= 0.0, "");
    n = std::sqrt(n / static_cast<double>(symmetries.size()));

    // Save results
    representative = r;
    character      = e;
    norm           = n;
}

auto is_real(small_symmetry_t const& symmetry) noexcept -> bool
{
    return symmetry.eigenvalue.imag() == 0.0;
}

auto is_real(batched_small_symmetry_t const& symmetry) noexcept -> bool
{
    return std::all_of(std::begin(symmetry.eigenvalues), std::end(symmetry.eigenvalues),
                       [](auto const& x) { return x.imag() == 0.0; });
}

auto is_real(big_symmetry_t const& symmetry) noexcept -> bool
{
    return symmetry.eigenvalue.imag() == 0.0;
}

} // namespace lattice_symmetries

extern "C" LATTICE_SYMMETRIES_EXPORT ls_error_code ls_create_symmetry(ls_symmetry**  ptr,
                                                                      unsigned const length,
                                                                      unsigned const permutation[],
                                                                      unsigned const sector)
{
    using namespace lattice_symmetries;
    if (ptr == nullptr) { return LS_INVALID_ARGUMENT; }
    auto const r = compile(tcb::span{permutation, length});
    if (!r) {
        if (r.error().category() == get_error_category()) {
            return static_cast<ls_error_code>(r.error().value());
        }
        return LS_SYSTEM_ERROR;
    }
    auto periodicity = compute_periodicity(tcb::span{permutation, length});
    if (sector >= periodicity) { return LS_INVALID_SECTOR; }
    auto const eigenvalue = compute_eigenvalue(sector, periodicity);

    // NOLINTNEXTLINE: 64 spins is the max system size which can be represented by uint64_t
    if (length > 64U) {
        auto p = std::make_unique<ls_symmetry>(std::in_place_type_t<big_symmetry_t>{}, r.value(),
                                               sector, periodicity, eigenvalue);
        *ptr   = p.release();
    }
    else {
        auto p = std::make_unique<ls_symmetry>(std::in_place_type_t<small_symmetry_t>{}, r.value(),
                                               sector, periodicity, eigenvalue);
        *ptr   = p.release();
    }
    return LS_SUCCESS;
}

extern "C" LATTICE_SYMMETRIES_EXPORT void ls_destroy_symmetry(ls_symmetry* symmetry)
{
    LATTICE_SYMMETRIES_CHECK(symmetry != nullptr, "trying to destroy a nullptr");
    std::default_delete<ls_symmetry>{}(symmetry);
}

extern "C" LATTICE_SYMMETRIES_EXPORT unsigned ls_get_sector(ls_symmetry const* symmetry)
{
    LATTICE_SYMMETRIES_CHECK(symmetry != nullptr, "trying to dereference a nullptr");
    return std::visit([](auto const& x) noexcept { return x.sector; }, symmetry->payload);
}

extern "C" LATTICE_SYMMETRIES_EXPORT unsigned ls_get_periodicity(ls_symmetry const* symmetry)
{
    LATTICE_SYMMETRIES_CHECK(symmetry != nullptr, "trying to dereference a nullptr");
    return std::visit([](auto const& x) noexcept { return x.periodicity; }, symmetry->payload);
}

extern "C" LATTICE_SYMMETRIES_EXPORT void ls_get_eigenvalue(ls_symmetry const* symmetry, void* out)
{
    LATTICE_SYMMETRIES_CHECK(symmetry != nullptr, "trying to dereference a nullptr");
    auto const value =
        std::visit([](auto const& x) noexcept -> std::complex<double> { return x.eigenvalue; },
                   symmetry->payload);
    std::memcpy(out, &value, sizeof(std::complex<double>));
}

extern "C" LATTICE_SYMMETRIES_EXPORT double ls_get_phase(ls_symmetry const* symmetry)
{
    LATTICE_SYMMETRIES_CHECK(symmetry != nullptr, "trying to dereference a nullptr");
    return std::visit(
        [](auto const& x) noexcept {
            return static_cast<double>(x.sector) / static_cast<double>(x.periodicity);
        },
        symmetry->payload);
}

extern "C" LATTICE_SYMMETRIES_EXPORT unsigned
ls_symmetry_get_number_spins(ls_symmetry const* symmetry)
{
    return std::visit([](auto const& x) noexcept { return x.network.width; }, symmetry->payload);
}

namespace lattice_symmetries {
struct symmetry_apply_fn_t {
    uint64_t* bits;

    auto operator()(small_symmetry_t const& symmetry) const noexcept -> void
    {
        bits[0] = symmetry.network(bits[0]);
    }

    auto operator()(big_symmetry_t const& symmetry) const noexcept -> void
    {
        auto const n    = symmetry.network.width;
        auto const size = (n + 63U) / 64U;
        bits512    buffer; // NOLINT: buffer is initialized on the next line
        std::copy(bits, bits + size, static_cast<uint64_t*>(buffer.words));
        symmetry.network(buffer);
        std::copy(static_cast<uint64_t const*>(buffer.words),
                  static_cast<uint64_t const*>(buffer.words) + size, bits);
    }
};
} // namespace lattice_symmetries

extern "C" LATTICE_SYMMETRIES_EXPORT void ls_apply_symmetry(ls_symmetry const* symmetry,
                                                            uint64_t           bits[])
{
    return std::visit(lattice_symmetries::symmetry_apply_fn_t{bits}, symmetry->payload);
}
