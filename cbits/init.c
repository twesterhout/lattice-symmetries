#include "lattice_symmetries_haskell.h"
#include <HsFFI.h>
#include <stdio.h>

void ls_hs_init(void) { hs_init(NULL, NULL); }

void ls_hs_exit(void) { hs_exit(); }

void *ls_hs_symbol_table[] = {&ls_hs_init, &ls_hs_exit, &ls_hs_foo};
