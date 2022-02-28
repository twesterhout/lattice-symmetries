#include "csr.h"

typedef void (*ls_csr_binary_op_type)(element_type const *,
                                      element_type const *, element_type *);

static void ls_csr_binary_op(ls_csr_matrix const *A, ls_csr_matrix const *B,
                             ls_csr_matrix *C, ls_csr_binary_op_type op) {
  C->offsets[0] = 0;
  index_type nnz = 0;
  index_type const number_rows = A->dimension;
  for (unsigned i = 0; i < number_rows; ++i) {
    index_type A_pos = A->offsets[i];
    index_type B_pos = B->offsets[i];
    index_type const A_end = A->offsets[i + 1];
    index_type const B_end = B->offsets[i + 1];

    while (A_pos < A_end && B_pos < B_end) {
      index_type const A_j = A->columns[A_pos];
      index_type const B_j = B->columns[B_pos];
      element_type result;

      if (A_j == B_j) {
        op(A->off_diag_elements + A_pos, B->off_diag_elements + B_pos, &result);
        if (result != 0) {
          C->columns[nnz] = A_j;
          C->off_diag_elements[nnz] = result;
          ++nnz;
        }
        ++A_pos;
        ++B_pos;
      } else if (A_j < B_j) {
        element_type const zero = 0;
        op(A->off_diag_elements + A_pos, &zero, &result);
        if (result != 0) {
          C->columns[nnz] = A_j;
          C->off_diag_elements[nnz] = result;
          ++nnz;
        }
        ++A_pos;
      } else { // B_j < A_j
        element_type const zero = 0;
        op(&zero, B->off_diag_elements + B_pos, &result);
        if (result != 0) {
          C->columns[nnz] = B_j;
          C->off_diag_elements[nnz] = result;
          ++nnz;
        }
        ++B_pos;
      }
    }

    for (; A_pos < A_end; ++A_pos) {
      element_type const zero = 0;
      element_type result;
      op(A->off_diag_elements + A_pos, &zero, &result);
      if (result != 0) {
        C->columns[nnz] = A->columns[A_pos];
        C->off_diag_elements[nnz] = result;
        ++nnz;
      }
    }
    for (; B_pos < B_end; ++B_pos) {
      element_type const zero = 0;
      element_type result;
      op(&zero, B->off_diag_elements + B_pos, &result);
      if (result != 0) {
        C->columns[nnz] = B->columns[B_pos];
        C->off_diag_elements[nnz] = result;
        ++nnz;
      }
    }

    C->offsets[i + 1] = nnz;
  }
  for (index_type i = 0; i < A->dimension; ++i) {
    op(A->diag_elements + i, B->diag_elements + i, C->diag_elements + i);
  }
  C->dimension = A->dimension;
  C->number_nonzero = nnz;
}

static void plus_op(element_type const *a, element_type const *b,
                    element_type *c) {
  *c = *a + *b;
}

void ls_csr_plus(ls_csr_matrix const *A, ls_csr_matrix const *B,
                 ls_csr_matrix *C) {
  ls_csr_binary_op(A, B, C, &plus_op);
}

static void minus_op(element_type const *a, element_type const *b,
                     element_type *c) {
  *c = *a - *b;
}

void ls_csr_minus(ls_csr_matrix const *A, ls_csr_matrix const *B,
                  ls_csr_matrix *C) {
  ls_csr_binary_op(A, B, C, &minus_op);
}

static void times_op(element_type const *a, element_type const *b,
                     element_type *c) {
  *c = *a * *b;
}

void ls_csr_times(ls_csr_matrix const *A, ls_csr_matrix const *B,
                  ls_csr_matrix *C) {
  ls_csr_binary_op(A, B, C, &times_op);
}

void ls_csr_kron(ls_csr_matrix const *A, ls_csr_matrix const *B,
                 ls_csr_matrix *C) {
  C->offsets[0] = 0;
  index_type nnz = 0;
  for (index_type A_i = 0; A_i < A->dimension; ++A_i) {
    for (index_type B_i = 0; B_i < B->dimension; ++B_i) {
      index_type const i = A_i * A->dimension + B_i;
      for (index_type A_pos = A->offsets[A_i]; A_pos < A->offsets[A_i + 1];
           ++A_pos) {
        for (index_type B_pos = B->offsets[B_i]; B_pos < B->offsets[B_i + 1];
             ++B_pos) {
          C->columns[nnz] =
              A->columns[A_pos] * A->dimension + B->columns[B_pos];
          C->off_diag_elements[nnz] =
              A->off_diag_elements[A_pos] * B->off_diag_elements[B_pos];
          ++nnz;
        }
      }
      C->offsets[i + 1] = nnz;
    }
  }
  for (index_type A_i = 0; A_i < A->dimension; ++A_i) {
    for (index_type B_i = 0; B_i < B->dimension; ++B_i) {
      index_type const i = A_i * A->dimension + B_i;
      C->diag_elements[i] = A->diag_elements[A_i] * B->diag_elements[B_i];
    }
  }
  C->dimension = A->dimension * B->dimension;
  C->number_nonzero = nnz;
}
