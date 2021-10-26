#include "lattice_symmetries_haskell.h"
#include <HsFFI.h>
#include <stdio.h>

void ls_hs_init(void) { hs_init(NULL, NULL); }

void ls_hs_exit(void) { hs_exit(); }

void *ls_hs_symbol_table[] = {
    &ls_hs_init, &ls_hs_exit, &ls_hs_basis_and_hamiltonian_from_yaml,
    &ls_hs_destroy_spin_basis, &ls_hs_destroy_operator};
