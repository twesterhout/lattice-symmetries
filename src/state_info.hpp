#include "basis.hpp"

namespace lattice_symmetries {

auto get_state_info_v2(basis_base_t const& basis_header, small_basis_t const& basis_body,
                       uint64_t bits, uint64_t& representative, std::complex<double>& character,
                       double& norm) noexcept -> void;

} // namespace lattice_symmetries
