#include "lattice_symmetries_haskell.h"
#include <stdio.h>

int main() {
  ls_hs_init();
  ls_enable_logging();
  printf("foo() = %i\n", ls_hs_foo());
  ls_spin_basis *basis = ls_hs_load_basis_from_yaml("heisenberg_chain_4.yaml");
  ls_hs_exit();
}
