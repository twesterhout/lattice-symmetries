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
    auto get_flip_mask(unsigned const n) noexcept -> bits64
    {
        // Play nice and do not shift by 64 bits
        // NOLINTNEXTLINE: 64 is the number of bits in bits64
        return n == 0U ? uint64_t{0} : ((~uint64_t{0}) >> (64U - n));
    }

    struct alignas(32) batch_acc_64_arch_t {
        vcl::Vec8uq const original;
        vcl::Vec8uq       r;
        vcl::Vec8d        n;
        vcl::Vec8d        e_re;
        vcl::Vec8d        e_im;

        explicit batch_acc_64_arch_t(uint64_t const bits) noexcept
            : original{bits}, r{original}, n{0.0}, e_re{1.0}, e_im{0.0}
        {}

        auto update(vcl::Vec8uq const& x, vcl::Vec8d const& real, vcl::Vec8d const& imag) noexcept
            -> void
        {
            auto const smaller = x < r;

            r    = vcl::select(smaller, x, r);
            e_re = vcl::select(smaller, real, e_re);
            e_im = vcl::select(smaller, imag, e_im);
            n    = vcl::if_add(x == original, n, real);
        }

        auto update_norm_only(vcl::Vec8uq const& x, vcl::Vec8d const& real) noexcept -> bool
        {
            auto const smaller = x < r;
            if (vcl::horizontal_or(smaller)) { return false; }
            n = vcl::if_add(x == original, n, real);
            return true;
        }

        auto reduce_norm_only() const noexcept -> double { return vcl::horizontal_add(n); }

        auto reduce(uint64_t& repr, std::complex<double>& character, double& norm) noexcept -> void;
    };

    LATTICE_SYMMETRIES_NOINLINE
    auto batch_acc_64_arch_t::reduce(uint64_t& repr, std::complex<double>& character,
                                     double& norm) noexcept -> void
    {
        // Compute norm
        norm = reduce_norm_only();

        // Reduce from 512 to 256
        auto const r256_low   = r.get_low();
        auto const r256_high  = r.get_high();
        auto const smaller256 = r256_low < r256_high;
        auto const r256       = vcl::select(smaller256, r256_low, r256_high);
        auto const e_re256    = vcl::select(smaller256, e_re.get_low(), e_re.get_high());
        auto const e_im256    = vcl::select(smaller256, e_im.get_low(), e_im.get_high());

        // Reduce from 526 to 128
        auto const r128_low   = r256.get_low();
        auto const r128_high  = r256.get_high();
        auto const smaller128 = r128_low < r128_high;
        auto const r128       = vcl::select(smaller128, r128_low, r128_high);
        auto const e_re128    = vcl::select(smaller128, e_re256.get_low(), e_re256.get_high());
        auto const e_im128    = vcl::select(smaller128, e_im256.get_low(), e_im256.get_high());

        // Reduce from 128 to 64
        auto const r64_low   = r128.extract(0);
        auto const r64_high  = r128.extract(1);
        auto const smaller64 = r64_low < r64_high;
        if (smaller64) {
            repr      = r64_low;
            character = std::complex{e_re128.extract(0), e_im128.extract(0)};
        }
        else {
            repr      = r64_high;
            character = std::complex{e_re128.extract(1), e_im128.extract(1)};
        }
    }
} // namespace

namespace detail {
    auto get_state_info_arch(basis_base_t const& basis_header, small_basis_t const& basis_body,
                             uint64_t bits, uint64_t& representative,
                             std::complex<double>& character, double& norm) noexcept -> void
    {
        if (!basis_header.has_symmetries) {
            representative = bits;
            character      = {1.0, 0.0};
            norm           = 1.0;
            return;
        }
        auto const flip_mask         = get_flip_mask(basis_header.number_spins);
        auto const flip_mask_vector  = vcl::Vec8uq{flip_mask};
        auto const flip_coeff        = static_cast<double>(basis_header.spin_inversion);
        auto const flip_coeff_vector = vcl::Vec8d{flip_coeff};

        batch_acc_64_arch_t acc{bits};
        for (auto const& symmetry : basis_body.batched_symmetries) {
            auto x = acc.original;
            benes_forward_arch(&x, symmetry.network.masks, symmetry.network.depth,
                               symmetry.network.deltas);

            vcl::Vec8d real;
            vcl::Vec8d imag;
            real.load_a(symmetry.eigenvalues_real.data());
            imag.load_a(symmetry.eigenvalues_imag.data());
            acc.update(x, real, imag);
            if (basis_header.spin_inversion != 0) {
                x ^= flip_mask_vector;
                if (basis_header.spin_inversion != 1) {
                    real *= flip_coeff_vector;
                    imag *= flip_coeff_vector;
                }
                acc.update(x, real, imag);
            }
        }

        uint64_t             r;
        std::complex<double> e;
        double               n;
        acc.reduce(r, e, n);

        for (auto const& symmetry : basis_body.other_symmetries) {
            auto const y = symmetry.network(bits);
            if (y < r) {
                r = y;
                e = symmetry.eigenvalue;
            }
            else if (y == bits) {
                n += symmetry.eigenvalue.real();
            }
            if (basis_header.spin_inversion != 0) {
                auto const inverted = y ^ flip_mask;
                if (inverted < r) {
                    r = inverted;
                    e = flip_coeff * symmetry.eigenvalue;
                }
                else if (inverted == bits) {
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
                                * (batch_size * basis_body.batched_symmetries.size()
                                   + basis_body.other_symmetries.size());
        n = std::sqrt(n / static_cast<double>(group_size));

        // Save results
        representative = r;
        character      = e;
        norm           = n;
    }

    auto is_representative_arch(basis_base_t const& basis_header, small_basis_t const& basis_body,
                                uint64_t bits) noexcept -> bool
    {
        if (!basis_header.has_symmetries) { return true; }

        auto const flip_mask         = get_flip_mask(basis_header.number_spins);
        auto const flip_mask_vector  = vcl::Vec8uq{flip_mask};
        auto const flip_coeff        = static_cast<double>(basis_header.spin_inversion);
        auto const flip_coeff_vector = vcl::Vec8d{flip_coeff};

        batch_acc_64_arch_t acc{bits};
        for (auto const& symmetry : basis_body.batched_symmetries) {
            auto x = acc.original;
            benes_forward_arch(&x, symmetry.network.masks, symmetry.network.depth,
                               symmetry.network.deltas);
            vcl::Vec8d real;
            real.load_a(symmetry.eigenvalues_real.data());

            if (!acc.update_norm_only(x, real)) { return false; }
            if (basis_header.spin_inversion != 0) {
                x ^= flip_mask_vector;
                real *= flip_coeff_vector;
                if (!acc.update_norm_only(x, real)) { return false; };
            }
        }

        auto n = acc.reduce_norm_only();
        for (auto const& symmetry : basis_body.other_symmetries) {
            auto const y = symmetry.network(bits);
            if (y < bits) { return false; }
            if (y == bits) { n += symmetry.eigenvalue.real(); }
        }

        // We need to detect the case when norm is not zero, but only because of
        // inaccurate arithmetics
        constexpr auto norm_threshold = 1.0e-5;
        if (std::abs(n) <= norm_threshold) { n = 0.0; }
        LATTICE_SYMMETRIES_ASSERT(n >= 0.0, "");
        return n > 0.0;
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
#endif
} // namespace lattice_symmetries
