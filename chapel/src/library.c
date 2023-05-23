#include <stdint.h>
#include <stddef.h>

void ls_chpl_init(void);
void ls_chpl_finalize(void);

extern void chpl__init_BatchedOperator(int64_t _ln, int32_t _fn);
extern void chpl__init_CommonParameters(int64_t _ln, int32_t _fn);
extern void chpl__init_ConcurrentAccessor(int64_t _ln, int32_t _fn);
extern void chpl__init_DistributedMatrixVector(int64_t _ln, int32_t _fn);
extern void chpl__init_FFI(int64_t _ln, int32_t _fn);
extern void chpl__init_ForeignTypes(int64_t _ln, int32_t _fn);
extern void chpl__init_LatticeSymmetries(int64_t _ln, int32_t _fn);
extern void chpl__init_StatesEnumeration(int64_t _ln, int32_t _fn);
extern void chpl__init_Vector(int64_t _ln, int32_t _fn);
extern void chpl_library_init(int argc, char *argv[]);
extern void chpl_library_finalize(void);

void ls_chpl_init(void) {
  int const argc = 1;
  char const *argv[2] = {"lattice_symmetries", NULL};
  chpl_library_init(argc, (char**)argv);
  chpl__init_BatchedOperator(1, 2);
  chpl__init_CommonParameters(1, 2);
  chpl__init_ConcurrentAccessor(1, 2);
  chpl__init_DistributedMatrixVector(1, 2);
  chpl__init_FFI(1, 2);
  chpl__init_ForeignTypes(1, 2);
  chpl__init_StatesEnumeration(1, 2);
  chpl__init_Vector(1, 2);
  chpl__init_LatticeSymmetries(1, 2);
}

void ls_chpl_finalize(void) { chpl_library_finalize(); }
