#include "state_info.hpp"
#include "../bits.hpp"
#include "benes_forward_512.hpp"
#include "benes_forward_64.hpp"
#include <vectorclass.h>

#if LATTICE_SYMMETRIES_HAS_AVX2()
#    define ARCH avx2
#elif LATTICE_SYMMETRIES_HAS_AVX()
#    define ARCH avx
#elif LATTICE_SYMMETRIES_HAS_SSE4()
#    define ARCH sse4
#else
#    define ARCH sse2
#endif

namespace lattice_symmetries::ARCH {
namespace vcl = VCL_NAMESPACE;

inline constexpr auto required_alignment = 32;

// We set 64-bit mask to
//
//     0000....0011....1111
//               ~~~~~~~~~~
//                   n
LATTICE_SYMMETRIES_FORCEINLINE auto get_flip_mask_64(unsigned const n) noexcept -> ls_bits64
{
    // Play nice and do not shift by 64 bits
    // NOLINTNEXTLINE: 64 is the number of bits in ls_bits64
    return n == 0U ? uint64_t{0} : ((~uint64_t{0}) >> (64U - n));
}

LATTICE_SYMMETRIES_FORCEINLINE auto get_flip_mask_512(unsigned n) noexcept -> ls_bits512
{
    ls_bits512 mask;
    auto       i = 0U;
    // NOLINTNEXTLINE: 64 is the number of bits in uint64_t
    for (; n >= 64U; ++i, n -= 64) {
        mask.words[i] = ~uint64_t{0};
    }
    if (n != 0U) {
        mask.words[i] = get_flip_mask_64(n);
        ++i;
    }
    for (; i < std::size(mask.words); ++i) {
        mask.words[i] = uint64_t{0};
    }
    return mask;
}

LATTICE_SYMMETRIES_FORCEINLINE auto operator^=(ls_bits512& x, ls_bits512 const& y) noexcept
    -> ls_bits512&
{
    for (auto i = 0; i < static_cast<int>(std::size(x.words)); ++i) {
        x.words[i] ^= y.words[i];
    }
    return x;
}

struct alignas(required_alignment) batch_acc_64_t {
    vcl::Vec8uq const original;
    vcl::Vec8uq       r;
    vcl::Vec8q        i;
    vcl::Vec8d        n;

    LATTICE_SYMMETRIES_FORCEINLINE
    explicit batch_acc_64_t(uint64_t const bits) noexcept
        : original{bits}, r{original}, i{std::numeric_limits<int64_t>::max()}, n{0.0}
    {}

    LATTICE_SYMMETRIES_FORCEINLINE
    auto update(vcl::Vec8uq const& x, vcl::Vec8q const& index, vcl::Vec8d const& real) noexcept
        -> void
    {
        n = vcl::if_add(x == original, n, real);

        auto const smaller = x < r;
        if (vcl::horizontal_or(smaller)) {
            r = vcl::select(smaller, x, r);
            i = vcl::select(smaller, index, i);
        }
    }

    LATTICE_SYMMETRIES_FORCEINLINE
    auto update_first_few(vcl::Vec8uq const& x, vcl::Vec8q const& index, vcl::Vec8d const& real,
                          unsigned const count) noexcept -> void
    {
        LATTICE_SYMMETRIES_ASSERT(count < 8, "");
        auto const include = vcl::Vec8uq{0, 1, 2, 3, 4, 5, 6, 7} < count;

        n = vcl::if_add((x == original) && include, n, real);

        auto const smaller = (x < r) && include;
        if (vcl::horizontal_or(smaller)) {
            r = vcl::select(smaller, x, r);
            i = vcl::select(smaller, index, i);
        }
    }

    LATTICE_SYMMETRIES_FORCEINLINE
    auto update_norm_only(vcl::Vec8uq const& x, vcl::Vec8d const& real) noexcept -> bool
    {
        auto const smaller = x < r;
        if (vcl::horizontal_or(smaller)) { return false; }
        n = vcl::if_add(x == original, n, real);
        return true;
    }

    LATTICE_SYMMETRIES_FORCEINLINE
    auto update_first_few_norm_only(vcl::Vec8uq const& x, vcl::Vec8d const& real,
                                    unsigned const count) noexcept -> bool
    {
        auto const include = vcl::Vec8uq{0, 1, 2, 3, 4, 5, 6, 7} < count;
        auto const smaller = (x < r) && include;
        if (vcl::horizontal_or(smaller)) { return false; }
        n = vcl::if_add((x == original) && include, n, real);
        return true;
    }

    [[nodiscard]] LATTICE_SYMMETRIES_FORCEINLINE auto reduce_norm_only() const noexcept -> double
    {
        return vcl::horizontal_add(n);
    }

    [[nodiscard]] LATTICE_SYMMETRIES_FORCEINLINE auto reduce() const noexcept
        -> std::tuple<uint64_t, int64_t, double>
    {
        // Compute norm
        auto const norm = reduce_norm_only();

        // Reduce from 512 to 256
        auto const r256_low   = r.get_low();
        auto const r256_high  = r.get_high();
        auto const smaller256 = r256_low < r256_high;
        auto const r256       = vcl::select(smaller256, r256_low, r256_high);
        auto const i256       = vcl::select(smaller256, i.get_low(), i.get_high());

        // Reduce from 526 to 128
        auto const r128_low   = r256.get_low();
        auto const r128_high  = r256.get_high();
        auto const smaller128 = r128_low < r128_high;
        auto const r128       = vcl::select(smaller128, r128_low, r128_high);
        auto const i128       = vcl::select(smaller128, i256.get_low(), i256.get_high());

        // Reduce from 128 to 64
        auto const r64_low   = r128.extract(0);
        auto const r64_high  = r128.extract(1);
        auto const smaller64 = r64_low < r64_high;
        if (smaller64) { return {r64_low, i128.extract(0), norm}; }
        return {r64_high, i128.extract(1), norm};
    }
};

LATTICE_SYMMETRIES_FORCEINLINE auto
apply_symmetry(vcl::Vec8uq& x, batched_small_symmetry_t const& symmetry) noexcept -> void
{
#if LATTICE_SYMMETRIES_HAS_AVX2()
    __m256i x0 = x.get_low();
    __m256i x1 = x.get_high();
    ARCH::benes_forward_64_direct(x0, x1, symmetry.network);
    x = vcl::Vec8uq{x0, x1};
#else
    __m128i x0 = x.get_low().get_low();
    __m128i x1 = x.get_low().get_high();
    __m128i x2 = x.get_high().get_low();
    __m128i x3 = x.get_high().get_high();
    ARCH::benes_forward_64_direct(x0, x1, x2, x3, symmetry.network);
    x = vcl::Vec8uq{{x0, x1}, {x2, x3}};
#endif
}

LATTICE_SYMMETRIES_FORCEINLINE auto apply_symmetry(ls_bits512&           x,
                                                   big_symmetry_t const& symmetry) noexcept -> void
{
    ARCH::benes_forward_512(x, symmetry.network);
}

auto get_state_info_64(basis_base_t const& basis_header, small_basis_t const& basis_body,
                       uint64_t bits, uint64_t& representative, std::complex<double>& character,
                       double& norm) noexcept -> void
{
    if (!basis_header.has_symmetries) {
        representative = bits;
        character      = {1.0, 0.0};
        norm           = 1.0;
        return;
    }
    vcl::Vec8uq flip_mask;
    vcl::Vec8d  flip_coeff;
    if (basis_header.spin_inversion != 0) {
        flip_mask  = vcl::Vec8uq{get_flip_mask_64(basis_header.number_spins)};
        flip_coeff = vcl::Vec8d{static_cast<double>(basis_header.spin_inversion)};
    }

    batch_acc_64_t acc{bits};

    vcl::Vec8q i_v{0, 1, 2, 3, 4, 5, 6, 7};
    vcl::Vec8q constant_8{8};
    for (auto const& symmetry : basis_body.batched_symmetries) {
        auto x = acc.original;
        apply_symmetry(x, symmetry);
        vcl::Vec8d real;
        real.load_a(symmetry.eigenvalues_real.data());
        acc.update(x, i_v, real);
        if (basis_header.spin_inversion != 0) {
            x ^= flip_mask;
            if (basis_header.spin_inversion != 1) { real = -real; }
            acc.update(x, -i_v, real);
        }
        i_v += constant_8;
    }
    if (basis_body.other_symmetries.has_value()) {
        auto const& symmetry = *basis_body.other_symmetries;

        auto x = acc.original;
        apply_symmetry(x, symmetry);
        vcl::Vec8d real;
        real.load_a(symmetry.eigenvalues_real.data());
        acc.update_first_few(x, i_v, real, basis_body.number_other_symmetries);
        if (basis_header.spin_inversion != 0) {
            x ^= flip_mask;
            if (basis_header.spin_inversion != 1) { real = -real; }
            acc.update_first_few(x, -i_v, real, basis_body.number_other_symmetries);
        }
    }
    auto [r, i, n] = acc.reduce();
    representative = r;

    if (i == std::numeric_limits<int64_t>::max()) { character = {1.0, 0.0}; }
    else {
        auto const i_abs       = static_cast<uint64_t>(std::abs(i));
        auto const batch_index = i_abs / 8;
        auto const rest_index  = i_abs % 8;

        auto const& s = batch_index == basis_body.batched_symmetries.size()
                            ? *basis_body.other_symmetries
                            : basis_body.batched_symmetries[batch_index];
        auto const e = std::complex{s.eigenvalues_real[rest_index], s.eigenvalues_imag[rest_index]};
        character    = i < 0 ? (static_cast<double>(basis_header.spin_inversion) * e) : e;
    }

    // We need to detect the case when norm is not zero, but only because of
    // inaccurate arithmetics
    constexpr auto norm_threshold = 1.0e-5;
    if (std::abs(n) <= norm_threshold) { n = 0.0; }
    LATTICE_SYMMETRIES_ASSERT(n >= 0.0, "");
    auto const group_size =
        (static_cast<unsigned>(basis_header.spin_inversion != 0) + 1)
        * (batch_size * basis_body.batched_symmetries.size() + basis_body.number_other_symmetries);
    n    = std::sqrt(n / static_cast<double>(group_size));
    norm = n;
}

auto is_representative_64(basis_base_t const& basis_header, small_basis_t const& basis_body,
                          uint64_t bits) noexcept -> bool
{
    if (!basis_header.has_symmetries) { return true; }

    auto const flip_mask  = vcl::Vec8uq{get_flip_mask_64(basis_header.number_spins)};
    auto const flip_coeff = vcl::Vec8d{static_cast<double>(basis_header.spin_inversion)};

    batch_acc_64_t acc{bits};
    for (auto const& symmetry : basis_body.batched_symmetries) {
        auto x = acc.original;
        apply_symmetry(x, symmetry);
        vcl::Vec8d real;
        real.load_a(symmetry.eigenvalues_real.data());

        if (!acc.update_norm_only(x, real)) { return false; }
        if (basis_header.spin_inversion != 0) {
            x ^= flip_mask;
            real *= flip_coeff;
            if (!acc.update_norm_only(x, real)) { return false; }
        }
    }
    if (basis_body.other_symmetries.has_value()) {
        auto const& symmetry = *basis_body.other_symmetries;
        auto const  count    = basis_body.number_other_symmetries;

        auto x = acc.original;
        apply_symmetry(x, symmetry);
        vcl::Vec8d real;
        real.load_a(symmetry.eigenvalues_real.data());
        if (!acc.update_first_few_norm_only(x, real, count)) { return false; }
        if (basis_header.spin_inversion != 0) {
            x ^= flip_mask;
            real *= flip_coeff;
            if (!acc.update_first_few_norm_only(x, real, count)) { return false; };
        }
    }

    auto n = acc.reduce_norm_only();
    // We need to detect the case when norm is not zero, but only because of
    // inaccurate arithmetics
    constexpr auto norm_threshold = 1.0e-5;
    if (std::abs(n) <= norm_threshold) { n = 0.0; }
    LATTICE_SYMMETRIES_ASSERT(n >= 0.0, "");
    return n > 0.0;
}

auto get_state_info_512(basis_base_t const& basis_header, big_basis_t const& basis_body,
                        ls_bits512 const& bits, ls_bits512& representative,
                        std::complex<double>& character, double& norm) noexcept -> void
{
    if (!basis_header.has_symmetries) {
        representative = bits;
        character      = {1.0, 0.0};
        norm           = 1.0;
        return;
    }
    auto const flip_mask  = get_flip_mask_512(basis_header.number_spins);
    auto const flip_coeff = static_cast<double>(basis_header.spin_inversion);

    ls_bits512 buffer; // NOLINT: buffer is initialized inside the loop before it is used
    auto       r = bits;
    auto       n = 0.0;
    auto       e = std::complex<double>{1.0};

    for (auto const& symmetry : basis_body.symmetries) {
        buffer = bits;
        symmetry.network(buffer);
        if (buffer < r) {
            r = buffer;
            e = symmetry.eigenvalue;
        }
        else if (buffer == bits) {
            n += symmetry.eigenvalue.real();
        }
        if (basis_header.spin_inversion != 0) {
            buffer ^= flip_mask;
            if (buffer < r) {
                r = buffer;
                e = flip_coeff * symmetry.eigenvalue;
            }
            else if (buffer == bits) {
                n += flip_coeff * symmetry.eigenvalue.real();
            }
        }
    }

    // We need to detect the case when norm is not zero, but only because of
    // inaccurate arithmetics
    constexpr auto norm_threshold = 1.0e-5;
    if (std::abs(n) <= norm_threshold) { n = 0.0; }
    LATTICE_SYMMETRIES_ASSERT(n >= 0.0, "");
    auto const group_size = (static_cast<unsigned>(basis_header.spin_inversion != 0) + 1)
                            * basis_body.symmetries.size();
    n = std::sqrt(n / static_cast<double>(group_size));

    // Save results
    representative = r;
    character      = e;
    norm           = n;
}
} // namespace lattice_symmetries::ARCH

#if defined(LATTICE_SYMMETRIES_ADD_DISPATCH_CODE)
extern "C" {
using get_state_info_64_type = auto (*)(lattice_symmetries::basis_base_t const&  basis_header,
                                        lattice_symmetries::small_basis_t const& basis_body,
                                        uint64_t bits, uint64_t& representative,
                                        std::complex<double>& character, double& norm) noexcept
                               -> void;

static auto resolve_get_state_info_64() noexcept -> get_state_info_64_type
{
    using namespace lattice_symmetries;
    return &sse4::get_state_info_64;
}

using is_representative_64_type = auto (*)(lattice_symmetries::basis_base_t const&  basis_header,
                                           lattice_symmetries::small_basis_t const& basis_body,
                                           uint64_t bits) noexcept -> bool;

static auto resolve_is_representative_64() noexcept -> is_representative_64_type
{
    using namespace lattice_symmetries;
    return &sse4::is_representative_64;
}

using get_state_info_512_type = auto (*)(lattice_symmetries::basis_base_t const& basis_header,
                                         lattice_symmetries::big_basis_t const&  basis_body,
                                         ls_bits512 const& bits, ls_bits512& representative,
                                         std::complex<double>& character, double& norm) noexcept
                                -> void;

static auto resolve_get_state_info_512() noexcept -> get_state_info_512_type
{
    using namespace lattice_symmetries;
    return &sse4::get_state_info_512;
}
} // extern "C"

namespace lattice_symmetries {
__attribute__((ifunc("resolve_get_state_info_64"))) auto
get_state_info_64(basis_base_t const& basis_header, small_basis_t const& basis_body, uint64_t bits,
                  uint64_t& representative, std::complex<double>& character, double& norm) noexcept
    -> void;

__attribute__((ifunc("resolve_is_representative_64"))) auto
is_representative_64(basis_base_t const& basis_header, small_basis_t const& basis_body,
                     uint64_t bits) noexcept -> bool;

__attribute__((ifunc("resolve_get_state_info_512"))) auto
get_state_info_512(basis_base_t const& basis_header, big_basis_t const& basis_body,
                   ls_bits512 const& bits, ls_bits512& representative,
                   std::complex<double>& character, double& norm) noexcept -> void;
} // namespace lattice_symmetries
#endif
