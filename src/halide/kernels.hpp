#pragma once

#include "basis.hpp"

namespace lattice_symmetries {

auto make_state_info_kernel(ls_flat_spin_basis const& basis) noexcept -> state_info_kernel_type;
auto make_is_representative_kernel(ls_flat_spin_basis const& basis) noexcept
    -> is_representative_kernel_type;

} // namespace lattice_symmetries
