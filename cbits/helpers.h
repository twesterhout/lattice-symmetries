#ifndef HELPERS_H
#define HELPERS_H

// #include "LatticeSymmetries_stub.h"
// #include "lattice_symmetries_haskell.h"
#include <stdint.h>

typedef struct ls_csr_matrix {
  unsigned *offsets;
  unsigned *columns;
  _Complex double *off_diag_elements;
  _Complex double *diag_elements;
  unsigned dimension;
  unsigned number_nonzero;
} ls_csr_matrix;

typedef struct ls_bit_index {
  uint8_t word;
  uint8_t bit;
} ls_bit_index;

typedef struct ls_term {
  ls_csr_matrix const *const matrix;
  ls_bit_index const *const tuple;
  uint64_t const *const sign_mask;
  unsigned const tuple_size;
} ls_term;

void ls_internal_apply_operator(uint64_t number_terms, ls_term const terms[],
                                uint64_t number_spins, uint64_t const *spins,
                                uint64_t capacity, uint64_t number_words,
                                uint64_t *out_spins,
                                _Complex double *out_coeffs,
                                uint64_t *out_counts);

// typedef unsigned (*ls_term_gather_fn)(uint64_t const * /*source*/,
//                                       ls_bit_index const * /*tuple*/);
// typedef void (*ls_term_scatter_fn)(unsigned, ls_bit_index const * /*tuple*/,
//                                    uint64_t * /*destination*/);
// typedef void (*ls_spin_copy_fn)(uint64_t const * /*source*/,
//                                 uint64_t * /*destination*/);
// typedef void (*ls_spin_fill_fn)(uint64_t const * /*2D source*/,
//                                 uint64_t /*count*/,
//                                 uint64_t * /*2D destination*/);

// typedef struct ls_sparse_operator {
//   ls_flat_spin_basis const *basis;
//   unsigned number_terms;
//   ls_term *terms;
// } ls_sparse_operator;

// unsigned ls_internal_term_gather_1(uint64_t const *, ls_bit_index const *);
// unsigned ls_internal_term_gather_2(uint64_t const *, ls_bit_index const *);
// unsigned ls_internal_term_gather_3(uint64_t const *, ls_bit_index const *);
// unsigned ls_internal_term_gather_4(uint64_t const *, ls_bit_index const *);

// void ls_internal_term_scatter_1(unsigned, ls_bit_index const *, uint64_t *);
// void ls_internal_term_scatter_2(unsigned, ls_bit_index const *, uint64_t *);
// void ls_internal_term_scatter_3(unsigned, ls_bit_index const *, uint64_t *);
// void ls_internal_term_scatter_4(unsigned, ls_bit_index const *, uint64_t *);

// void ls_internal_spin_copy_1(uint64_t const *, uint64_t *);
// void ls_internal_spin_copy_8(uint64_t const *, uint64_t *);

// void ls_internal_spin_fill_1(uint64_t const *, uint64_t, uint64_t *);
// void ls_internal_spin_fill_8(uint64_t const *, uint64_t, uint64_t *);

// void ls_csr_matrix_from_dense(unsigned const dimension,
//                               _Complex double const *dense, unsigned
//                               *offsets, unsigned *columns, _Complex double
//                               *off_diag_elements, _Complex double
//                               *diag_elements);
// void ls_dense_from_csr_matrix(unsigned const dimension, unsigned const
// *offsets,
//                               unsigned const *columns,
//                               _Complex double const *off_diag_elements,
//                               _Complex double const *diag_elements,
//                               _Complex double *dense);

// uint64_t ls_hs_apply_term(ls_term const *term, uint64_t const *spin,
//                           ls_output_buffer *out);

#endif // HELPERS_H
