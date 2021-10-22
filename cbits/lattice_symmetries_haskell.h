#pragma once

#if defined(__cplusplus)
#include <cstdint>
#else
#include <stdint.h>
#endif

#include <lattice_symmetries/lattice_symmetries.h>

#if defined(__cplusplus)
extern "C" {
#endif

void ls_hs_init(void);
void ls_hs_exit(void);
int ls_hs_foo(void);
ls_spin_basis *ls_hs_load_basis_from_yaml(char const *path);

#ifdef __cplusplus
} // extern "C"
#endif
