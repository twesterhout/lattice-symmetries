#include <HalideRuntime.h>
#include <lattice_symmetries.h>


#define BUFFER(_rank, _dims, _type, _data) \
    (halide_buffer_t){ .device = 0, .device_interface = 0, .host = (uint8_t *)(_data), .flags = 0, .type = (_type), .dimensions = (_rank), .dim = (_dims), .padding = 0}
#define BUFFER_I64(_rank, _dims, _data) BUFFER(_rank, _dims, ((struct halide_type_t){.code = halide_type_int, .bits = 64, .lanes = 1}), _data)
#define BUFFER_U64(_rank, _dims, _data) BUFFER(_rank, _dims, ((struct halide_type_t){.code = halide_type_uint, .bits = 64, .lanes = 1}), _data)
#define BUFFER_F64(_rank, _dims, _data) BUFFER(_rank, _dims, ((struct halide_type_t){.code = halide_type_float, .bits = 64, .lanes = 1}), _data)

void ls_invoke_diag_matrix_kernel(ls_diag_terms const* terms, int32_t const count, uint64_t const* alphas, double *coeffs_re, double *coeffs_im)
{
    halide_dimension_t batch_dims[1] = {(halide_dimension_t){.min = 0, .extent = count, .stride = 1}};
    halide_buffer_t alphas_buf = BUFFER_I64(1, batch_dims, alphas);
    halide_buffer_t coeffs_re_buf = BUFFER_F64(1, batch_dims, coeffs_re);

    halide_dimension_t term_dims[1] = {(halide_dimension_t){.min = 0, .extent = terms->number_terms, .stride = 1}};
    halide_buffer_t v_re_buf = BUFFER_F64(1, term_dims, terms->v_re);

    LS_CHECK(terms->kernel != NULL, "kernel is NULL");
    if (terms->v_im == NULL && coeffs_im == NULL) {
        ls_diag_matrix_real_kernel const fn = terms->kernel;
        fn(&alphas_buf, &v_re_buf, &coeffs_re_buf);
    }
    else if (terms->v_im != NULL && coeffs_im != NULL) {
        halide_buffer_t coeffs_im_buf = BUFFER_F64(1, batch_dims, coeffs_im);
        halide_buffer_t v_im_buf = BUFFER_F64(1, term_dims, terms->v_im);
        ls_diag_matrix_complex_kernel const fn = terms->kernel;
        fn(&alphas_buf, &v_re_buf, &v_im_buf, &coeffs_re_buf, &coeffs_im_buf);
    }
    else {
        LS_FATAL_ERROR("incompatible v_im and coeffs_im");
    }
}

void ls_invoke_off_diag_matrix_kernel(ls_off_diag_terms const* terms, int32_t const count, uint64_t const* alphas, double *coeffs_re, double *coeffs_im)
{
    halide_dimension_t batch_dims[2] = {(halide_dimension_t){.min = 0, .extent = count, .stride = 1},
                                        (halide_dimension_t){.min = 0, .extent = terms->number_terms, .stride = count}};
    halide_buffer_t alphas_buf = BUFFER_I64(1, batch_dims, alphas);
    halide_buffer_t coeffs_re_buf = BUFFER_F64(2, batch_dims, coeffs_re);
    halide_dimension_t term_dims[2] = {(halide_dimension_t){.min = 0, .extent = terms->number_reduced, .stride = 1},
                                       (halide_dimension_t){.min = 0, .extent = terms->number_terms, .stride = terms->number_reduced}};
    halide_buffer_t v_re_buf = BUFFER_F64(2, term_dims, terms->v_re);

    LS_CHECK(terms->kernel != NULL, "kernel is NULL");
    if (terms->v_im == NULL && coeffs_im == NULL) {
        ls_off_diag_matrix_real_kernel const fn = terms->kernel;
        fn(&alphas_buf, &v_re_buf, &coeffs_re_buf);
    }
    else if (terms->v_im != NULL && coeffs_im != NULL) {
        halide_buffer_t coeffs_im_buf = BUFFER_F64(2, batch_dims, coeffs_im);
        halide_buffer_t v_im_buf = BUFFER_F64(2, term_dims, terms->v_im);
        ls_off_diag_matrix_complex_kernel const fn = terms->kernel;
        fn(&alphas_buf, &v_re_buf, &v_im_buf, &coeffs_re_buf, &coeffs_im_buf);
    }
    else {
        LS_FATAL_ERROR("incompatible v_im and coeffs_im");
    }
}

void ls_invoke_xored_state_to_index(ls_xored_state_to_index const* ctx, int32_t const count, uint64_t const* alphas, int32_t const number_terms, uint64_t const* masks, int64_t* indices)
{
    halide_dimension_t alphas_dims[1] = {(halide_dimension_t){.min = 0, .extent = count, .stride = 1}};
    halide_buffer_t alphas_buf = BUFFER_I64(1, alphas_dims, alphas);
    halide_dimension_t masks_dims[1] = {(halide_dimension_t){.min = 0, .extent = number_terms, .stride = 1}};
    halide_buffer_t masks_buf = BUFFER_I64(1, masks_dims, masks);
    halide_dimension_t states_dims[1] = {(halide_dimension_t){.min = 0, .extent = ctx->number_states, .stride = 1}};
    halide_buffer_t states_buf = BUFFER_I64(1, states_dims, ctx->basis_states);
    halide_dimension_t indices_dims[2] = {(halide_dimension_t){.min = 0, .extent = count, .stride = 1},
                                          (halide_dimension_t){.min = 0, .extent = number_terms, .stride = count}};
    halide_buffer_t indices_buf = BUFFER_I64(2, indices_dims, indices);

    LS_CHECK(ctx->kernel != NULL, "kernel is NULL");
    ls_xored_state_to_index_kernel const fn = ctx->kernel;
    fn(&alphas_buf, &masks_buf, &states_buf, &indices_buf);
}
