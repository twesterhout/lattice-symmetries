#include "state_info.hpp"
#include "kernels.hpp"
#include "macros.hpp"
#include <vectorclass.h>

#if LATTICE_SYMMETRIES_HAS_AVX2()
#    define batch_acc_64_arch_t batch_acc_64_avx2_t
#    define get_state_info_arch get_state_info_avx2
#    define is_representative_arch is_representative_avx2
#    define benes_forward_arch detail::benes_forward_avx2
#elif LATTICE_SYMMETRIES_HAS_AVX()
#    define batch_acc_64_arch_t batch_acc_64_avx_t
#    define get_state_info_arch get_state_info_avx
#    define is_representative_arch is_representative_avx
#    define benes_forward_arch detail::benes_forward_avx
#else
#    define batch_acc_64_arch_t batch_acc_64_sse2_t
#    define get_state_info_arch get_state_info_sse2
#    define is_representative_arch is_representative_sse2
#    define benes_forward_arch detail::benes_forward_sse2
#endif

namespace lattice_symmetries {
namespace {
    // We set 64-bit mask to
    //
    //     0000....0011....1111
    //               ~~~~~~~~~~
    //                   n
    auto get_flip_mask_64(unsigned const n) noexcept -> bits64
    {
        // Play nice and do not shift by 64 bits
        // NOLINTNEXTLINE: 64 is the number of bits in bits64
        return n == 0U ? uint64_t{0} : ((~uint64_t{0}) >> (64U - n));
    }

    auto get_flip_mask_512(unsigned n) noexcept -> bits512
    {
        bits512 mask;
        auto    i = 0U;
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

    LATTICE_SYMMETRIES_FORCEINLINE
    auto operator^=(bits512& x, bits512 const& y) noexcept -> bits512&
    {
        for (auto i = 0; i < static_cast<int>(std::size(x.words)); ++i) {
            x.words[i] ^= y.words[i];
        }
        return x;
    }

    struct alignas(32) batch_acc_64_arch_t {
        vcl::Vec8uq const original;
        vcl::Vec8uq       r;
        vcl::Vec8q        i;
        vcl::Vec8d        n;

        explicit batch_acc_64_arch_t(uint64_t const bits) noexcept
            : original{bits}, r{original}, i{std::numeric_limits<int64_t>::max()}, n{0.0}
        {}

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

        auto update_norm_only(vcl::Vec8uq const& x, vcl::Vec8d const& real) noexcept -> bool
        {
            auto const smaller = x < r;
            if (vcl::horizontal_or(smaller)) { return false; }
            n = vcl::if_add(x == original, n, real);
            return true;
        }

        auto update_first_few_norm_only(vcl::Vec8uq const& x, vcl::Vec8d const& real,
                                        unsigned const count) noexcept -> bool
        {
            auto const include = vcl::Vec8uq{0, 1, 2, 3, 4, 5, 6, 7} < count;
            auto const smaller = (x < r) && include;
            if (vcl::horizontal_or(smaller)) { return false; }
            n = vcl::if_add((x == original) && include, n, real);
            return true;
        }

        [[nodiscard]] auto reduce_norm_only() const noexcept -> double
        {
            return vcl::horizontal_add(n);
        }

        auto reduce(uint64_t& repr) const noexcept -> std::tuple<int64_t, double>;
    };

    LATTICE_SYMMETRIES_NOINLINE
    auto batch_acc_64_arch_t::reduce(uint64_t& repr) const noexcept -> std::tuple<int64_t, double>
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
        if (smaller64) {
            repr = r64_low;
            return {i128.extract(0), norm};
        }
        repr = r64_high;
        return {i128.extract(1), norm};
    }
} // namespace

namespace detail {
    auto get_state_info_arch(basis_base_t const& basis_header, small_basis_t const& basis_body,
                             uint64_t bits, uint64_t& representative,
                             std::complex<double>& character, double& norm) noexcept -> void
    {
        // std::printf("get_state_info_arch(%zu)\n", bits);
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

        batch_acc_64_arch_t acc{bits};

        vcl::Vec8q i_v{0, 1, 2, 3, 4, 5, 6, 7};
        vcl::Vec8q constant_8{8};
        for (auto const& symmetry : basis_body.batched_symmetries) {
            auto x = acc.original;
            benes_forward_arch(&x, symmetry.network.masks, symmetry.network.depth,
                               symmetry.network.deltas);
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
            benes_forward_arch(&x, symmetry.network.masks, symmetry.network.depth,
                               symmetry.network.deltas);
            vcl::Vec8d real;
            real.load_a(symmetry.eigenvalues_real.data());
            acc.update_first_few(x, i_v, real, basis_body.number_other_symmetries);
            if (basis_header.spin_inversion != 0) {
                x ^= flip_mask;
                if (basis_header.spin_inversion != 1) { real = -real; }
                acc.update_first_few(x, -i_v, real, basis_body.number_other_symmetries);
            }
        }
        auto [i, n] = acc.reduce(representative);

        if (i == std::numeric_limits<int64_t>::max()) { character = {1.0, 0.0}; }
        else {
            auto const i_abs       = static_cast<uint64_t>(std::abs(i));
            auto const batch_index = i_abs / 8;
            auto const rest_index  = i_abs % 8;

            auto const& s = batch_index == basis_body.batched_symmetries.size()
                                ? *basis_body.other_symmetries
                                : basis_body.batched_symmetries[batch_index];
            auto const  e =
                std::complex{s.eigenvalues_real[rest_index], s.eigenvalues_imag[rest_index]};
            character = i < 0 ? (-e) : e;
        }

        // We need to detect the case when norm is not zero, but only because of
        // inaccurate arithmetics
        constexpr auto norm_threshold = 1.0e-5;
        if (std::abs(n) <= norm_threshold) { n = 0.0; }
        LATTICE_SYMMETRIES_ASSERT(n >= 0.0, "");
        auto const group_size = (static_cast<unsigned>(basis_header.spin_inversion != 0) + 1)
                                * (batch_size * basis_body.batched_symmetries.size()
                                   + basis_body.number_other_symmetries);
        n    = std::sqrt(n / static_cast<double>(group_size));
        norm = n;
    }

    auto is_representative_arch(basis_base_t const& basis_header, small_basis_t const& basis_body,
                                uint64_t bits) noexcept -> bool
    {
        if (!basis_header.has_symmetries) { return true; }

        auto const flip_mask  = vcl::Vec8uq{get_flip_mask_64(basis_header.number_spins)};
        auto const flip_coeff = vcl::Vec8d{static_cast<double>(basis_header.spin_inversion)};

        batch_acc_64_arch_t acc{bits};
        for (auto const& symmetry : basis_body.batched_symmetries) {
            auto x = acc.original;
            benes_forward_arch(&x, symmetry.network.masks, symmetry.network.depth,
                               symmetry.network.deltas);
            vcl::Vec8d real;
            real.load_a(symmetry.eigenvalues_real.data());

            if (!acc.update_norm_only(x, real)) { return false; }
            if (basis_header.spin_inversion != 0) {
                x ^= flip_mask;
                real *= flip_coeff;
                if (!acc.update_norm_only(x, real)) { return false; };
            }
        }
        if (basis_body.other_symmetries.has_value()) {
            auto const& symmetry = *basis_body.other_symmetries;
            auto const  count    = basis_body.number_other_symmetries;

            auto x = acc.original;
            benes_forward_arch(&x, symmetry.network.masks, symmetry.network.depth,
                               symmetry.network.deltas);
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

    auto get_state_info_arch(basis_base_t const& basis_header, big_basis_t const& basis_body,
                             bits512 const& bits, bits512& representative,
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

        bits512 buffer; // NOLINT: buffer is initialized inside the loop before it is used
        auto    r = bits;
        auto    n = 0.0;
        auto    e = std::complex<double>{1.0};

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
} // namespace detail

#if defined(LATTICE_SYMMETRIES_ADD_DISPATCH_CODE)
auto get_state_info(basis_base_t const& basis_header, small_basis_t const& basis_body,
                    uint64_t bits, uint64_t& representative, std::complex<double>& character,
                    double& norm) noexcept -> void
{
    if (__builtin_cpu_supports("avx2")) {
        detail::get_state_info_avx2(basis_header, basis_body, bits, representative, character,
                                    norm);
    }
    else if (__builtin_cpu_supports("avx")) {
        detail::get_state_info_avx(basis_header, basis_body, bits, representative, character, norm);
    }
    else {
        detail::get_state_info_sse2(basis_header, basis_body, bits, representative, character,
                                    norm);
    }
}

auto is_representative(basis_base_t const& basis_header, small_basis_t const& basis_body,
                       uint64_t bits) noexcept -> bool
{
    if (__builtin_cpu_supports("avx2")) {
        return detail::is_representative_avx2(basis_header, basis_body, bits);
    }
    else if (__builtin_cpu_supports("avx")) {
        return detail::is_representative_avx(basis_header, basis_body, bits);
    }
    else {
        return detail::is_representative_sse2(basis_header, basis_body, bits);
    }
}

auto get_state_info(basis_base_t const& basis_header, big_basis_t const& basis_body,
                    bits512 const& bits, bits512& representative, std::complex<double>& character,
                    double& norm) noexcept -> void
{
    if (__builtin_cpu_supports("avx2")) {
        detail::get_state_info_avx2(basis_header, basis_body, bits, representative, character,
                                    norm);
    }
    else if (__builtin_cpu_supports("avx")) {
        detail::get_state_info_avx(basis_header, basis_body, bits, representative, character, norm);
    }
    else {
        detail::get_state_info_sse2(basis_header, basis_body, bits, representative, character,
                                    norm);
    }
}
#endif
} // namespace lattice_symmetries
