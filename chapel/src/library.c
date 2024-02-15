#define _GNU_SOURCE

#include "lattice_symmetries_chapel.h"

#include <HalideRuntime.h>
#include <lattice_symmetries_functions.h>

#ifdef LS_CHPL_LIBRARY

extern void chpl__init_CommonParameters(int64_t _ln, int32_t _fn);
extern void chpl__init_FFI(int64_t _ln, int32_t _fn);
extern void chpl__init_Utils(int64_t _ln, int32_t _fn);
extern void chpl__init_ForeignTypes(int64_t _ln, int32_t _fn);
extern void chpl__init_Timing(int64_t _ln, int32_t _fn);
extern void chpl__init_Vector(int64_t _ln, int32_t _fn);
extern void chpl__init_BatchedOperator(int64_t _ln, int32_t _fn);
extern void chpl__init_ConcurrentQueue(int64_t _ln, int32_t _fn);
extern void chpl__init_StatesEnumeration(int64_t _ln, int32_t _fn);
extern void chpl__init_MatrixVectorProduct(int64_t _ln, int32_t _fn);
extern void chpl__init_LatticeSymmetries(int64_t _ln, int32_t _fn);

extern void chpl_library_init(int argc, char *argv[]);
extern void chpl_library_finalize(void);

void ls_chpl_init(void) {
  int const argc = 1;
  char const *argv[2] = {"lattice_symmetries", NULL};
  chpl_library_init(argc, (char**)argv);
  chpl__init_CommonParameters(1, 2);
  chpl__init_FFI(1, 2);
  chpl__init_Utils(1, 2);
  chpl__init_ForeignTypes(1, 2);
  chpl__init_Timing(1, 2);
  chpl__init_Vector(1, 2);
  chpl__init_BatchedOperator(1, 2);
  chpl__init_ConcurrentQueue(1, 2);
  chpl__init_StatesEnumeration(1, 2);
  chpl__init_MatrixVectorProduct(1, 2);
  chpl__init_LatticeSymmetries(1, 2);
}

void ls_chpl_finalize(void) { chpl_library_finalize(); }

#endif // LS_CHPL_LIBRARY

ls_hs_basis_info const *ls_chpl_get_basis_info(ls_hs_basis const *basis) {
    ls_hs_basis_info const *ptr = atomic_load(&basis->info);
    if (ptr == NULL) {
        ls_hs_init_basis_info(basis);
        ptr = atomic_load(&basis->info);
    }
    return ptr;
}

ls_hs_is_representative_kernel_type_v2
ls_chpl_get_is_representative_kernel(ls_hs_basis const *basis) {
    ls_hs_is_representative_kernel_type_v2 fun = atomic_load(&basis->is_representative_kernel);
    if (fun == NULL) {
        ls_hs_init_is_representative_kernel(basis);
        fun = atomic_load(&basis->is_representative_kernel);
    }
    return fun;
}

void ls_chpl_invoke_is_representative_kernel(ls_hs_is_representative_kernel_type_v2 kernel,
                                             int64_t const count, uint64_t const *basis_states,
                                             double *norms) {
    halide_dimension_t dims[1] = {
        (halide_dimension_t){.min = 0, .extent = (int32_t)count, .stride = 1}};
    halide_buffer_t basis_states_buf = (halide_buffer_t){
        .device = 0,
        .device_interface = 0,
        .host = (uint8_t *)basis_states,
        .flags = 0,
        .type = (struct halide_type_t){.code = halide_type_uint, .bits = 64, .lanes = 1},
        .dimensions = 1,
        .dim = dims,
        .padding = 0};
    halide_buffer_t norms_buf = (halide_buffer_t){
        .device = 0,
        .device_interface = 0,
        .host = (uint8_t *)norms,
        .flags = 0,
        .type = (struct halide_type_t){.code = halide_type_float, .bits = 64, .lanes = 1},
        .dimensions = 1,
        .dim = dims,
        .padding = 0};
    kernel(&basis_states_buf, &norms_buf);
}

ls_hs_state_to_index_kernel_type
ls_chpl_get_state_to_index_kernel(ls_hs_basis const *basis) {
    ls_hs_state_to_index_kernel_type fun = atomic_load(&basis->state_to_index_kernel);
    // fprintf(stderr, "ls_chpl_get_state_to_index_kernel: %p\n", (void*)fun);
    if (fun == NULL) {
        ls_hs_init_state_to_index_kernel(basis);
        fun = atomic_load(&basis->state_to_index_kernel);
    }
    return fun;
}

void ls_chpl_invoke_state_to_index_kernel(ls_hs_state_to_index_kernel_type kernel,
                                          int64_t const count, uint64_t const *basis_states,
                                          int64_t *indices) {
    halide_dimension_t dims[1] = {
        (halide_dimension_t){.min = 0, .extent = (int32_t)count, .stride = 1}};
    halide_buffer_t basis_states_buf = (halide_buffer_t){
        .device = 0,
        .device_interface = 0,
        .host = (uint8_t *)basis_states,
        .flags = 0,
        .type = (struct halide_type_t){.code = halide_type_uint, .bits = 64, .lanes = 1},
        .dimensions = 1,
        .dim = dims,
        .padding = 0};
    halide_buffer_t indices_buf = (halide_buffer_t){
        .device = 0,
        .device_interface = 0,
        .host = (uint8_t *)indices,
        .flags = 0,
        .type = (struct halide_type_t){.code = halide_type_int, .bits = 64, .lanes = 1},
        .dimensions = 1,
        .dim = dims,
        .padding = 0};
    kernel(&basis_states_buf, &indices_buf);
}
