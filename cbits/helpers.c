#include "helpers.h"
#include <string.h>

void ls_csr_matrix_from_dense(unsigned const dimension,
                              _Complex double const *const dense,
                              unsigned *offsets, unsigned *columns,
                              _Complex double *off_diag_elements,
                              _Complex double *diag_elements) {
  offsets[0] = 0;
  for (unsigned i = 0; i < dimension; ++i) {
    unsigned nonzero_in_row = 0;
    for (unsigned j = 0; j < dimension; ++j) {
      _Complex double const element = dense[i * dimension + j];
      if (i == j) {
        diag_elements[i] = element;
      } else if (element != 0) {
        *columns = j;
        *off_diag_elements = element;
        ++nonzero_in_row;
        ++columns;
        ++off_diag_elements;
      }
    }
    offsets[i + 1] = offsets[i] + nonzero_in_row;
  }
}

void ls_dense_from_csr_matrix(unsigned const dimension,
                              unsigned const *const offsets,
                              unsigned const *const columns,
                              _Complex double const *const off_diag_elements,
                              _Complex double const *const diag_elements,
                              _Complex double *const dense) {
  memset(dense, 0,
         (size_t)dimension * (size_t)dimension * sizeof(_Complex double));
  for (unsigned i = 0; i < dimension; ++i) {
    for (unsigned k = offsets[i]; k < offsets[i + 1]; ++k) {
      unsigned const j = columns[k];
      dense[i * dimension + j] = off_diag_elements[k];
    }
    dense[i * dimension + i] = diag_elements[i];
  }
}
