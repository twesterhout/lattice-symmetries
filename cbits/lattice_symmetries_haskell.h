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

typedef struct ls_hs_spin_basis_v1 {
  ls_spin_basis const *const payload;
  void *const context;
} ls_hs_spin_basis_v1;

typedef struct ls_hs_operator_v1 {
  ls_operator const *const payload;
  void *const context;
} ls_hs_operator_v1;

void ls_hs_basis_and_hamiltonian_from_yaml(char const *path,
                                           ls_hs_spin_basis_v1 *basis,
                                           ls_hs_operator_v1 *hamiltonian);
void ls_hs_destroy_spin_basis(ls_hs_spin_basis_v1 *basis);
void ls_hs_destroy_operator(ls_hs_operator_v1 *op);

#if 0
typedef struct ls_sparse_operator ls_sparse_operator;
typedef struct ls_term ls_term;

typedef struct ls_hs_basis {
  ls_flat_spin_basis const *payload;
  void *context;
} ls_hs_basis;

typedef struct ls_hs_operator {
  ls_sparse_operator const *payload;
  void *context;
} ls_hs_operator;

typedef void (*ls_spin_copy_fn)(uint64_t const * /*source*/,
                                uint64_t * /*destination*/);
typedef void (*ls_spin_fill_fn)(uint64_t const * /*2D source*/,
                                uint64_t /*count*/,
                                uint64_t * /*2D destination*/);

typedef struct ls_output_buffer {
  uint64_t *spins;
  _Complex double *coeffs;
  _Complex double *const diagonal;
  uint64_t const number_words;
  ls_spin_copy_fn const spin_copy;
  ls_spin_fill_fn const spin_fill;
} ls_output_buffer;

typedef struct ls_workspace {
  uint64_t *spins;
  _Complex double *characters;
  double *norms;
} ls_workspace;

void ls_hs_create_operator(ls_hs_operator *self,
                           ls_flat_spin_basis const *basis,
                           unsigned number_terms, ls_term const *term);
void ls_hs_destroy_operator(ls_hs_operator *self);
void ls_hs_load_hamiltonian_from_yaml(ls_hs_operator *self, char const *path);
void ls_hs_operator_get_basis(ls_hs_operator const *self, ls_hs_basis *basis);

uint64_t ls_hs_apply_operator(ls_hs_operator const *term, uint64_t const *spin,
                              ls_output_buffer *out, ls_workspace *workspace);
#endif

#ifdef __cplusplus
} // extern "C"
#endif
