#include "../basis.hpp"

namespace lattice_symmetries {

#define LATTICE_SYMMETRIES_DECLARE()                                                               \
    auto get_state_info_64(basis_base_t const& basis_header, small_basis_t const& basis_body,      \
                           uint64_t bits, uint64_t& representative,                                \
                           std::complex<double>& character, double& norm) noexcept->void;          \
    auto is_representative_64(basis_base_t const& basis_header, small_basis_t const& basis_body,   \
                              uint64_t bits) noexcept->bool;                                       \
    auto get_state_info_512(basis_base_t const& basis_header, big_basis_t const& basis_body,       \
                            ls_bits512 const& bits, ls_bits512& representative,                    \
                            std::complex<double>& character, double& norm) noexcept->void;

#define LATTICE_SYMMETRIES_DECLARE_FOR_ARCH(arch)                                                  \
    namespace arch {                                                                               \
        LATTICE_SYMMETRIES_DECLARE()                                                               \
    } /* namespace arch */

LATTICE_SYMMETRIES_DECLARE()

LATTICE_SYMMETRIES_DECLARE_FOR_ARCH(avx2)
LATTICE_SYMMETRIES_DECLARE_FOR_ARCH(avx)
LATTICE_SYMMETRIES_DECLARE_FOR_ARCH(sse4)
LATTICE_SYMMETRIES_DECLARE_FOR_ARCH(sse2)

#undef LATTICE_SYMMETRIES_DECLARE
#undef LATTICE_SYMMETRIES_DECLARE_FOR_ARCH

} // namespace lattice_symmetries
