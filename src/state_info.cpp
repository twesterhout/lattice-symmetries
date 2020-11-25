#include "state_info.hpp"

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

    // auto init_flip_mask(bits512& mask, unsigned n) noexcept -> void
    // {
    //     auto i = 0U;
    //     // NOLINTNEXTLINE: 64 is the number of bits in uint64_t
    //     for (; n >= 64U; ++i, n -= 64) {
    //         mask.words[i] = ~uint64_t{0};
    //     }
    //     if (n != 0U) {
    //         init_flip_mask(mask.words[i], n);
    //         ++i;
    //     }
    //     for (; i < std::size(mask.words); ++i) {
    //         mask.words[i] = uint64_t{0};
    //     }
    // }
} // namespace

auto get_state_info_v2(basis_base_t const& basis_header, small_basis_t const& basis_body,
                       uint64_t bits, uint64_t& representative, std::complex<double>& character,
                       double& norm) noexcept -> void
{
    if (!basis_header.has_symmetries) {
        representative = bits;
        character      = {1.0, 0.0};
        norm           = 1.0;
        return;
    }
    constexpr auto       batch_size = batched_small_symmetry_t::batch_size;
    auto const           flip_mask  = get_flip_mask(basis_header.number_spins);
    auto const           flip_coeff = static_cast<double>(basis_header.spin_inversion);
    alignas(32) uint64_t initial[batch_size]; // NOLINT: 32-byte alignment for AVX
    alignas(32) uint64_t buffer[batch_size];  // NOLINT: same
    std::fill(std::begin(initial), std::end(initial), bits);

    auto r = bits;
    auto n = 0.0;
    auto e = std::complex<double>{1.0};

    for (auto const& symmetry : basis_body.batched_symmetries) {
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
            if (basis_header.spin_inversion != 0) {
                auto const inverted = buffer[i] ^ flip_mask;
                if (inverted < r) {
                    r = inverted;
                    e = flip_coeff * symmetry.eigenvalues[i];
                }
                else if (inverted == bits) {
                    n += flip_coeff * symmetry.eigenvalues[i].real();
                }
            }
        }
    }
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
    auto const group_size =
        (static_cast<unsigned>(basis_header.spin_inversion != 0) + 1)
        * (batch_size * basis_body.batched_symmetries.size() + basis_body.other_symmetries.size());
    n = std::sqrt(n / static_cast<double>(group_size));

    // Save results
    representative = r;
    character      = e;
    norm           = n;
}

} // namespace lattice_symmetries
