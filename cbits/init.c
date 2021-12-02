#include "lattice_symmetries_haskell.h"
#include <HsFFI.h>
#include <stdio.h>

void ls_hs_init(void) {
  LATTICE_SYMMETRIES_LOG_DEBUG("%s", "Calling hs_init...\n");
  hs_init(NULL, NULL);
  LATTICE_SYMMETRIES_LOG_DEBUG("%s", "Initialized RTS!\n");
}

void ls_hs_exit(void) {
  LATTICE_SYMMETRIES_LOG_DEBUG("%s", "Calling hs_exit...\n");
  hs_exit();
  LATTICE_SYMMETRIES_LOG_DEBUG("%s", "Deinitialized RTS!\n");
}

void *ls_hs_symbol_table[] = {&ls_hs_init,
                              &ls_hs_exit,
                              &ls_hs_basis_and_hamiltonian_from_yaml,
                              &ls_hs_destroy_spin_basis,
                              &ls_hs_destroy_operator,
                              &ls_hs_hdf5_create_dataset_u64,
                              &ls_hs_hdf5_write_1d_chunk_u64};
