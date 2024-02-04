#define _GNU_SOURCE
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <execinfo.h>
#include "library.h"

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


static int ls_internal_comp_int64(const void* _a, const void* _b, void* _ctx) {
  int64_t const a = *(int64_t const*)_a;
  int64_t const b = *(int64_t const*)_b;
  uint64_t const* keys = (uint64_t const*)_ctx;
  uint64_t const key_a = keys[a];
  uint64_t const key_b = keys[b];
  return (key_a > key_b) - (key_a < key_b);
}

void ls_internal_qsort_int64(int64_t* xs, size_t const count, uint64_t const* keys) {
  qsort_r(xs, count, sizeof(int64_t), ls_internal_comp_int64, (void*)keys);
}

void ls_internal_func_and_parent(void** func, void** parent) {
  void* buffer[100];
  int const size = backtrace(buffer, 100);
  char** names = backtrace_symbols(buffer, size);
  for (int i = 0; i < size; ++i) {
      fprintf(stderr, "buffer[%i] = %zu, %s\n", i, (uintptr_t)buffer[i], names[i]);
  }
  free(names);
  *func = buffer[0];
  *parent = buffer[1];
}


