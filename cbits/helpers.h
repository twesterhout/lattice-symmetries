#ifndef HELPERS_H
#define HELPERS_H

#include <lattice_symmetries/lattice_symmetries.h>
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

typedef struct ls_term {
  ls_csr_matrix matrix;
  unsigned number_tuples;
  unsigned tuple_size;
  ls_bit_index *tuples;
  ls_term_gather_fn gather_fn;
  ls_term_scatter_fn scatter_fn;
} ls_term;

typedef struct ls_output_buffer {
  uint64_t *spins;
  _Complex double *coeffs;
  _Complex double *const diagonal;
  uint64_t number_words;
} ls_output_buffer;

typedef struct ls_workspace {
  uint64_t *spins;
  _Complex double *characters;
  double *norms;
} ls_workspace;

typedef struct ls_sparse_operator {
  ls_flat_spin_basis const *basis;
  unsigned number_terms;
  ls_term *terms;
} ls_sparse_operator;

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

#endif // HELPERS_H
