#ifndef HELPERS_H
#define HELPERS_H

#include "LatticeSymmetries_stub.h"
#include "lattice_symmetries_haskell.h"
#include <stdint.h>

typedef struct ls_csr_matrix {
  unsigned dimension;
  unsigned number_nonzero;
  unsigned *offsets;
  unsigned *columns;
  _Complex double *off_diag_elements;
  _Complex double *diag_elements;
} ls_csr_matrix;

typedef struct ls_bit_index {
  uint8_t word;
  uint8_t bit;
} ls_bit_index;

typedef unsigned (*ls_term_gather_fn)(uint64_t const * /*source*/,
                                      ls_bit_index const * /*tuple*/);
typedef void (*ls_term_scatter_fn)(unsigned, ls_bit_index const * /*tuple*/,
                                   uint64_t * /*destination*/);
typedef void (*ls_spin_copy_fn)(uint64_t const * /*source*/,
                                uint64_t * /*destination*/);
typedef void (*ls_spin_fill_fn)(uint64_t const * /*2D source*/,
                                uint64_t /*count*/,
                                uint64_t * /*2D destination*/);

struct ls_term {
  ls_csr_matrix matrix;
  unsigned number_tuples;
  unsigned tuple_size;
  ls_bit_index *tuples;
  ls_term_gather_fn gather_fn;
  ls_term_scatter_fn scatter_fn;
};

struct ls_sparse_operator {
  ls_flat_spin_basis const *basis;
  unsigned number_terms;
  ls_term *terms;
};

unsigned ls_internal_term_gather_1(uint64_t const *, ls_bit_index const *);
unsigned ls_internal_term_gather_2(uint64_t const *, ls_bit_index const *);
unsigned ls_internal_term_gather_3(uint64_t const *, ls_bit_index const *);
unsigned ls_internal_term_gather_4(uint64_t const *, ls_bit_index const *);

void ls_internal_term_scatter_1(unsigned, ls_bit_index const *, uint64_t *);
void ls_internal_term_scatter_2(unsigned, ls_bit_index const *, uint64_t *);
void ls_internal_term_scatter_3(unsigned, ls_bit_index const *, uint64_t *);
void ls_internal_term_scatter_4(unsigned, ls_bit_index const *, uint64_t *);

void ls_internal_spin_copy_1(uint64_t const *, uint64_t *);
void ls_internal_spin_copy_8(uint64_t const *, uint64_t *);

void ls_internal_spin_fill_1(uint64_t const *, uint64_t, uint64_t *);
void ls_internal_spin_fill_8(uint64_t const *, uint64_t, uint64_t *);

void ls_csr_matrix_from_dense(unsigned const dimension,
                              _Complex double const *dense, unsigned *offsets,
                              unsigned *columns,
                              _Complex double *off_diag_elements,
                              _Complex double *diag_elements);
void ls_dense_from_csr_matrix(unsigned const dimension, unsigned const *offsets,
                              unsigned const *columns,
                              _Complex double const *off_diag_elements,
                              _Complex double const *diag_elements,
                              _Complex double *dense);

uint64_t ls_hs_apply_term(ls_term const *term, uint64_t const *spin,
                          ls_output_buffer *out);

#endif // HELPERS_H
