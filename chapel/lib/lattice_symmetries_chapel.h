#include "stdchpl.h"
#include "ctype.h"
#include "wctype.h"
#include "lattice_symmetries_haskell.h"
void chpl__init_BatchedOperator(int64_t _ln,
                                int32_t _fn);
void ls_chpl_operator_apply_diag(ls_hs_operator * matrixPtr,
                                 int64_t count,
                                 uint64_t * alphas,
                                 chpl_external_array * coeffs,
                                 int64_t numTasks);
void ls_chpl_operator_apply_off_diag(ls_hs_operator * matrixPtr,
                                     int64_t count,
                                     uint64_t * alphas,
                                     chpl_external_array * betas,
                                     chpl_external_array * coeffs,
                                     chpl_external_array * offsets,
                                     int64_t numTasks);
void chpl__init_CommonParameters(int64_t _ln,
                                 int32_t _fn);
void chpl__init_ConcurrentAccessor(int64_t _ln,
                                   int32_t _fn);
void chpl__init_DistributedMatrixVector(int64_t _ln,
                                        int32_t _fn);
void ls_chpl_matrix_vector_product(ls_hs_operator * matrixPtr,
                                   int32_t numVectors,
                                   _real64 * xPtr,
                                   _real64 * yPtr);
void chpl__init_FFI(int64_t _ln,
                    int32_t _fn);
void chpl__init_ForeignTypes(int64_t _ln,
                             int32_t _fn);
void chpl__init_LatticeSymmetries(int64_t _ln,
                                  int32_t _fn);
void ls_chpl_init_kernels(void);
void chpl__init_StatesEnumeration(int64_t _ln,
                                  int32_t _fn);
void ls_chpl_enumerate_representatives(ls_hs_basis * p,
                                       uint64_t lower,
                                       uint64_t upper,
                                       chpl_external_array * dest);
void chpl__init_Vector(int64_t _ln,
                       int32_t _fn);
