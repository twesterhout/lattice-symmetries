#include "symmetry.hpp"
#include <algorithm>
#include <memory>
#include <numeric>

namespace lattice_symmetries {

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
    if (sector == 0U) { return std::complex<double>{1.0, 0.0}; }
    if (2U * sector == periodicity) { return std::complex<double>{-1.0, 0.0}; }
    auto const arg = -static_cast<long double>(2U * sector) / static_cast<long double>(periodicity);
    auto       re  = static_cast<double>(std::cos(M_PIl * arg));
    auto       im  = static_cast<double>(std::sin(M_PIl * arg));
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
    constexpr auto       batch_size = batched_small_symmetry_t::batch_size;
    alignas(32) uint64_t initial[batch_size];
    alignas(32) uint64_t buffer[batch_size];
    std::fill(std::begin(initial), std::end(initial), bits);

    auto r = bits;
    auto n = 0.0;
    auto e = std::complex<double>{1.0};

    for (auto const& symmetry : batched_symmetries) {
        std::copy(std::begin(initial), std::end(initial), std::begin(buffer));
        symmetry.network(buffer);
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

auto get_state_info(std::vector<big_symmetry_t> const& symmetries, bits512 const& bits,
                    bits512& representative, std::complex<double>& character, double& norm) noexcept
    -> void
{
    bits512 buffer;
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

extern "C" ls_error_code ls_create_symmetry(ls_symmetry** ptr, unsigned const length,
                                            unsigned const permutation[], bool const flip,
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
    if (flip && periodicity % 2U != 0U) { periodicity *= 2U; }
    if (sector >= periodicity) { return LS_INVALID_SECTOR; }
    auto const eigenvalue = compute_eigenvalue(sector, periodicity);

    if (length > 64U) {
        auto p = std::make_unique<ls_symmetry>(std::in_place_type_t<big_symmetry_t>{}, r.value(),
                                               flip, sector, periodicity, eigenvalue);
        *ptr   = p.release();
    }
    else {
        auto p = std::make_unique<ls_symmetry>(std::in_place_type_t<small_symmetry_t>{}, r.value(),
                                               flip, sector, periodicity, eigenvalue);
        *ptr   = p.release();
    }
    return LS_SUCCESS;
}

extern "C" void ls_destroy_symmetry(ls_symmetry* symmetry)
{
    LATTICE_SYMMETRIES_CHECK(symmetry != nullptr, "trying to destroy a nullptr");
    std::default_delete<ls_symmetry>{}(symmetry);
}

extern "C" unsigned ls_get_sector(ls_symmetry const* symmetry)
{
    LATTICE_SYMMETRIES_CHECK(symmetry != nullptr, "trying to dereference a nullptr");
    return std::visit([](auto const& x) noexcept { return x.sector; }, symmetry->payload);
}

extern "C" bool ls_get_flip(ls_symmetry const* symmetry)
{
    LATTICE_SYMMETRIES_CHECK(symmetry != nullptr, "trying to dereference a nullptr");
    return std::visit([](auto const& x) noexcept { return x.network.flip; }, symmetry->payload);
}

extern "C" unsigned ls_get_periodicity(ls_symmetry const* symmetry)
{
    LATTICE_SYMMETRIES_CHECK(symmetry != nullptr, "trying to dereference a nullptr");
    return std::visit([](auto const& x) noexcept { return x.periodicity; }, symmetry->payload);
}

extern "C" void ls_get_eigenvalue(ls_symmetry const* symmetry, void* out)
{
    LATTICE_SYMMETRIES_CHECK(symmetry != nullptr, "trying to dereference a nullptr");
    *reinterpret_cast<std::complex<double>*>(out) =
        std::visit([](auto const& x) noexcept { return x.eigenvalue; }, symmetry->payload);
}

extern "C" double ls_get_phase(ls_symmetry const* symmetry)
{
    LATTICE_SYMMETRIES_CHECK(symmetry != nullptr, "trying to dereference a nullptr");
    return std::visit(
        [](auto const& x) noexcept {
            return static_cast<double>(x.sector) / static_cast<double>(x.periodicity);
        },
        symmetry->payload);
}

extern "C" unsigned ls_symmetry_get_number_spins(ls_symmetry const* symmetry)
{
    return std::visit([](auto const& x) noexcept { return x.network.width; }, symmetry->payload);
}
