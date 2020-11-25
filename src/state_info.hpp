#include "basis.hpp"

namespace lattice_symmetries {

auto get_state_info_v2(basis_base_t const& basis_header, small_basis_t const& basis_body,
                       uint64_t bits, uint64_t& representative, std::complex<double>& character,
                       double& norm) noexcept -> void;

auto get_state_info(basis_base_t const& basis_header, small_basis_t const& basis_body,
                    uint64_t bits, uint64_t& representative, std::complex<double>& character,
                    double& norm) noexcept -> void;
auto is_representative(basis_base_t const& basis_header, small_basis_t const& basis_body,
                       uint64_t bits) noexcept -> bool;

auto get_state_info(basis_base_t const& basis_header, big_basis_t const& basis_body,
                    bits512 const& bits, bits512& representative, std::complex<double>& character,
                    double& norm) noexcept -> void;

namespace detail {
    auto get_state_info_avx2(basis_base_t const& basis_header, small_basis_t const& basis_body,
                             uint64_t bits, uint64_t& representative,
                             std::complex<double>& character, double& norm) noexcept -> void;
    auto get_state_info_avx(basis_base_t const& basis_header, small_basis_t const& basis_body,
                            uint64_t bits, uint64_t& representative,
                            std::complex<double>& character, double& norm) noexcept -> void;
    auto get_state_info_sse2(basis_base_t const& basis_header, small_basis_t const& basis_body,
                             uint64_t bits, uint64_t& representative,
                             std::complex<double>& character, double& norm) noexcept -> void;
    auto is_representative_avx2(basis_base_t const& basis_header, small_basis_t const& basis_body,
                                uint64_t bits) noexcept -> bool;
    auto is_representative_avx(basis_base_t const& basis_header, small_basis_t const& basis_body,
                               uint64_t bits) noexcept -> bool;
    auto is_representative_sse2(basis_base_t const& basis_header, small_basis_t const& basis_body,
                                uint64_t bits) noexcept -> bool;

    auto get_state_info_avx2(basis_base_t const& basis_header, big_basis_t const& basis_body,
                             bits512 const& bits, bits512& representative,
                             std::complex<double>& character, double& norm) noexcept -> void;
    auto get_state_info_avx(basis_base_t const& basis_header, big_basis_t const& basis_body,
                            bits512 const& bits, bits512& representative,
                            std::complex<double>& character, double& norm) noexcept -> void;
    auto get_state_info_sse2(basis_base_t const& basis_header, big_basis_t const& basis_body,
                             bits512 const& bits, bits512& representative,
                             std::complex<double>& character, double& norm) noexcept -> void;
} // namespace detail

} // namespace lattice_symmetries
