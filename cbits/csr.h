#pragma once

#include "helpers.h"

typedef unsigned index_type;
typedef _Complex double element_type;

void ls_csr_plus(ls_csr_matrix const *A, ls_csr_matrix const *B,
                 ls_csr_matrix *C);

void ls_csr_minus(ls_csr_matrix const *A, ls_csr_matrix const *B,
                  ls_csr_matrix *C);

void ls_csr_times(ls_csr_matrix const *A, ls_csr_matrix const *B,
                  ls_csr_matrix *C);

void ls_csr_kron(ls_csr_matrix const *A, ls_csr_matrix const *B,
                 ls_csr_matrix *C);

void ls_csr_mult(ls_csr_matrix const *A, ls_csr_matrix const *B,
                 ls_csr_matrix *C);
