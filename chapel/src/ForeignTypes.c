#include "ForeignTypes.h"
#include <HalideRuntime.h>
#include <lattice_symmetries_functions.h>

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
