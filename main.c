#include "lattice_symmetries_haskell.h"
#include <stdio.h>

int main() {
  ls_hs_init();
  ls_enable_logging();
  ls_hs_spin_basis_v1 basis;
  ls_hs_operator_v1 operator;
  ls_hs_basis_and_hamiltonian_from_yaml("heisenberg_chain_4.yaml", &basis,
                                        &operator);
  printf("number_spins = %u\n", ls_get_number_spins(basis.payload));
  ls_hs_destroy_operator(&operator);
  ls_hs_destroy_spin_basis(&basis);
  ls_hs_exit();
}
