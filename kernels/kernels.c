#include "kernels.h"
#include "lattice_symmetries_haskell.h"
#include <HalideRuntime.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

// Macro to indicate x86_64
// #if (defined(_M_AMD64) || defined(_M_X64) || defined(__amd64))
// #define LS_X86_64() 1
// #else
// #define LS_X86_64() 0
// #endif

typedef int (*ls_internal_halide_state_info_kernel_type)(
    struct halide_buffer_t *_x_buffer, uint64_t _flip_mask,
    struct halide_buffer_t *_masks_buffer,
    struct halide_buffer_t *_eigvals_re_buffer,
    struct halide_buffer_t *_eigvals_im_buffer,
    struct halide_buffer_t *_shifts_buffer,
    struct halide_buffer_t *_representative_buffer,
    struct halide_buffer_t *_character_buffer,
    struct halide_buffer_t *_norm_buffer);

typedef int (*ls_internal_halide_is_representative_kernel_type)(
    struct halide_buffer_t *_x_buffer, uint64_t _flip_mask,
    struct halide_buffer_t *_masks_buffer,
    struct halide_buffer_t *_eigvals_re_buffer,
    struct halide_buffer_t *_shifts_buffer,
    struct halide_buffer_t *_is_representative_buffer,
    struct halide_buffer_t *_norm_buffer);

// typedef int (*ls_internal_halide_apply_off_diag_kernel_type)(
//     struct halide_buffer_t *_v_buffer, struct halide_buffer_t *_m_buffer,
//     struct halide_buffer_t *_r_buffer, struct halide_buffer_t *_x_buffer,
//     struct halide_buffer_t *_s_buffer, struct halide_buffer_t *_alpha_buffer,
//     struct halide_buffer_t *_beta_buffer,
//     struct halide_buffer_t *_coeff_buffer);

typedef struct ls_internal_halide_kernels_list {
  ls_internal_halide_state_info_kernel_type state_info_general;
  ls_internal_halide_state_info_kernel_type state_info_symmetric;
  ls_internal_halide_state_info_kernel_type state_info_antisymmetric;
  ls_internal_halide_is_representative_kernel_type is_representative_general;
  ls_internal_halide_is_representative_kernel_type is_representative_symmetric;
  ls_internal_halide_is_representative_kernel_type
      is_representative_antisymmetric;
  // ls_internal_halide_apply_off_diag_kernel_type apply_off_diag;
} ls_internal_halide_kernels_list;

// static void init_kernels_list(int const number_bits,
//                               ls_internal_halide_kernels_list *list) {
//   LS_CHECK(number_bits <= 64, "too many bits");
//   *list = (ls_internal_halide_kernels_list){
//       .state_info_general = &ls_internal_state_info_general_kernel,
//       .state_info_symmetric = &ls_internal_state_info_symmetric_kernel,
//       .state_info_antisymmetric = &ls_internal_state_info_antisymmetric_kernel,
//       .is_representative_general =
//           &ls_internal_is_representative_general_kernel,
//       .is_representative_symmetric =
//           &ls_internal_is_representative_symmetric_kernel,
//       .is_representative_antisymmetric =
//           &ls_internal_is_representative_antisymmetric_kernel,
//       // .apply_off_diag = &ls_internal_apply_off_diag_kernel,
//   };
//   // #if LS_X86_64() == 1
//   //   char const *arch = getenv("LATTICE_SYMMETRIES_ARCH");
//   //   if (arch == NULL) { // strcmp will invoke undefined behavior if we pass
//   //   it a
//   //                       // NULL pointer
//   //     arch = "";
//   //   }
//   //   __builtin_cpu_init();
//   //   if (strcmp(arch, "avx2") == 0) {
//   //     *list = LS_AVX2_KERNELS;
//   //   } else if (strcmp(arch, "avx") == 0) {
//   //     *list = LS_AVX_KERNELS;
//   //   } else if (strcmp(arch, "sse4_1") == 0) {
//   //     *list = LS_SSE4_KERNELS;
//   //   } else if (strcmp(arch, "generic") == 0) {
//   //     *list = LS_GENERIC_KERNELS;
//   //   } else if (__builtin_cpu_supports("avx2") > 0) {
//   //     *list = LS_AVX2_KERNELS;
//   //   } else if (__builtin_cpu_supports("avx") > 0) {
//   //     *list = LS_AVX_KERNELS;
//   //   } else if (__builtin_cpu_supports("sse4.1") > 0) {
//   //     *list = LS_SSE4_KERNELS;
//   //   } else {
//   //     *list = LS_GENERIC_KERNELS;
//   //   }
//   // #else
//   //   *list = LS_GENERIC_KERNELS;
//   // #endif
// }

static inline uint64_t get_flip_mask_64(unsigned const n) {
  return n == 0U ? (uint64_t)0 : ((~(uint64_t)0) >> (64U - n));
}

typedef struct ls_internal_halide_kernel_data {
  struct halide_buffer_t masks;
  struct halide_buffer_t eigvals_re;
  struct halide_buffer_t eigvals_im;
  struct halide_buffer_t shifts;
  struct halide_dimension_t masks_dims[2];
  struct halide_dimension_t shifts_dim;
  uint64_t flip_mask;
  ls_internal_halide_state_info_kernel_type state_info_kernel;
  ls_internal_halide_is_representative_kernel_type is_representative_kernel;
} ls_internal_halide_kernel_data;

ls_internal_halide_kernel_data *
ls_internal_create_halide_kernel_data(ls_hs_permutation_group const *g,
                                      int const spin_inversion) {
  ls_internal_halide_kernel_data *self =
      malloc(sizeof(ls_internal_halide_kernel_data));
  if (self == NULL) {
    return NULL;
  }

  self->masks = (struct halide_buffer_t){.device = 0,
                                         .device_interface = NULL,
                                         .host = (uint8_t *)g->masks,
                                         .flags = 0,
                                         .type = {halide_type_uint, 64, 1},
                                         .dimensions = 2,
                                         .dim = &(self->masks_dims[0]),
                                         .padding = NULL};
  self->eigvals_re =
      (struct halide_buffer_t){.device = 0,
                               .device_interface = NULL,
                               .host = (uint8_t *)g->eigvals_re,
                               .flags = 0,
                               .type = {halide_type_float, 64, 1},
                               .dimensions = 1,
                               .dim = &(self->masks_dims[1]),
                               .padding = NULL};
  self->eigvals_im =
      (struct halide_buffer_t){.device = 0,
                               .device_interface = NULL,
                               .host = (uint8_t *)g->eigvals_im,
                               .flags = 0,
                               .type = {halide_type_float, 64, 1},
                               .dimensions = 1,
                               .dim = &(self->masks_dims[1]),
                               .padding = NULL};
  self->shifts = (struct halide_buffer_t){
      .device = 0,
      .device_interface = NULL,
      .host = (uint8_t *)g->shifts,
      .flags = 0,
      .type = {halide_type_uint, 64, 1},
      .dimensions = 1,
      .dim = &(self->shifts_dim),
      .padding = NULL,
  };
  self->masks_dims[0] =
      (struct halide_dimension_t){.min = 0,
                                  .extent = (int32_t)g->number_shifts,
                                  .stride = (int32_t)g->number_masks,
                                  .flags = 0};
  self->masks_dims[1] = (struct halide_dimension_t){
      .min = 0, .extent = (int32_t)g->number_masks, .stride = 1, .flags = 0};
  self->shifts_dim =
      (struct halide_dimension_t){.min = 0,
                                  .extent = (int32_t)g->number_shifts,
                                  .stride = 1,
                                  .flags = 0},
  self->flip_mask = get_flip_mask_64((unsigned)g->number_bits);

  // ls_internal_halide_kernels_list list;
  // init_kernels_list(g->number_bits, &list);
  if (spin_inversion == 0) {
    self->state_info_kernel = &ls_internal_state_info_general_kernel;
    self->is_representative_kernel =
        &ls_internal_is_representative_general_kernel;
  } else if (spin_inversion == 1) {
    self->state_info_kernel = &ls_internal_state_info_symmetric_kernel;
    self->is_representative_kernel =
        &ls_internal_is_representative_symmetric_kernel;
  } else if (spin_inversion == -1) {
    self->state_info_kernel = &ls_internal_state_info_antisymmetric_kernel;
    self->is_representative_kernel =
        &ls_internal_is_representative_antisymmetric_kernel;
  } else {
    LS_FATAL_ERROR("invalid spin_inversion");
  }
  return self;
}

void ls_internal_destroy_halide_kernel_data(ls_internal_halide_kernel_data *p) {
  free(p);
}

void ls_hs_is_representative_halide_kernel(
    ptrdiff_t batch_size, uint64_t const *alphas, ptrdiff_t alphas_stride,
    uint8_t *are_representatives, double *norms, void const *private_data) {
  ls_internal_halide_kernel_data const *const state = private_data;

  halide_dimension_t alphas_dim =
      (halide_dimension_t){.min = 0,
                           .extent = (int32_t)batch_size,
                           .stride = (int32_t)alphas_stride,
                           .flags = 0};
  halide_buffer_t alphas_buf = (halide_buffer_t){
      .device = 0,
      .device_interface = NULL,
      .host = (uint8_t *)alphas,
      .flags = 0,
      .type = {halide_type_uint, 64, 1},
      .dimensions = 1,
      .dim = &alphas_dim,
      .padding = NULL,
  };

  halide_dimension_t is_repr_dim = (halide_dimension_t){
      .min = 0, .extent = (int32_t)batch_size, .stride = 1, .flags = 0};
  halide_buffer_t is_repr_buf = (halide_buffer_t){
      .device = 0,
      .device_interface = NULL,
      .host = are_representatives,
      .flags = 0,
      .type = {halide_type_uint, 8, 1},
      .dimensions = 1,
      .dim = &is_repr_dim,
      .padding = NULL,
  };

  halide_buffer_t norm_buf = (halide_buffer_t){
      .device = 0,
      .device_interface = NULL,
      .host = (uint8_t *)norms,
      .flags = 0,
      .type = {halide_type_float, 64, 1},
      .dimensions = 1,
      .dim = &is_repr_dim,
      .padding = NULL,
  };

  (*state->is_representative_kernel)(
      &alphas_buf, state->flip_mask, (halide_buffer_t *)&state->masks,
      (halide_buffer_t *)&state->eigvals_re, (halide_buffer_t *)&state->shifts,
      &is_repr_buf, &norm_buf);
}

void ls_hs_state_info_halide_kernel(ptrdiff_t batch_size,
                                    uint64_t const *alphas,
                                    ptrdiff_t alphas_stride, uint64_t *betas,
                                    ptrdiff_t betas_stride,
                                    ls_hs_scalar *characters, double *norms,
                                    void const *private_data) {
  ls_internal_halide_kernel_data const *const state = private_data;

  halide_dimension_t alphas_dim =
      (halide_dimension_t){.min = 0,
                           .extent = (int32_t)batch_size,
                           .stride = (int32_t)alphas_stride,
                           .flags = 0};
  halide_buffer_t alphas_buf = (halide_buffer_t){
      .device = 0,
      .device_interface = NULL,
      .host = (uint8_t *)alphas,
      .flags = 0,
      .type = {halide_type_uint, 64, 1},
      .dimensions = 1,
      .dim = &alphas_dim,
      .padding = NULL,
  };

  halide_dimension_t betas_dim =
      (halide_dimension_t){.min = 0,
                           .extent = (int32_t)batch_size,
                           .stride = (int32_t)betas_stride,
                           .flags = 0};
  halide_buffer_t betas_buf = (halide_buffer_t){
      .device = 0,
      .device_interface = NULL,
      .host = (uint8_t *)betas,
      .flags = 0,
      .type = {halide_type_uint, 64, 1},
      .dimensions = 1,
      .dim = &betas_dim,
      .padding = NULL,
  };

  halide_dimension_t characters_dims[2] = {
      (halide_dimension_t){
          .min = 0, .extent = (int32_t)batch_size, .stride = 2, .flags = 0},
      (halide_dimension_t){.min = 0, .extent = 2, .stride = 1, .flags = 0}};
  halide_buffer_t characters_buf = (halide_buffer_t){
      .device = 0,
      .device_interface = NULL,
      .host = (uint8_t *)characters,
      .flags = 0,
      .type = {halide_type_float, 64, 1},
      .dimensions = 2,
      .dim = characters_dims,
      .padding = NULL,
  };

  halide_dimension_t norms_dim = (halide_dimension_t){
      .min = 0, .extent = (int32_t)batch_size, .stride = 1, .flags = 0};
  halide_buffer_t norms_buf = (halide_buffer_t){
      .device = 0,
      .device_interface = NULL,
      .host = (uint8_t *)norms,
      .flags = 0,
      .type = {halide_type_float, 64, 1},
      .dimensions = 1,
      .dim = &norms_dim,
      .padding = NULL,
  };

  (*state->state_info_kernel)(
      &alphas_buf, state->flip_mask, (halide_buffer_t *)&state->masks,
      (halide_buffer_t *)&state->eigvals_re,
      (halide_buffer_t *)&state->eigvals_im, (halide_buffer_t *)&state->shifts,
      &betas_buf, &characters_buf, &norms_buf);
}

// typedef struct ls_hs_nonbranching_terms {
//   int number_terms;
//   int number_bits;
//   // number_words = ceil(number_bits / 64)
//   ls_hs_scalar const *v; // array of shape [number_terms]
//   uint64_t const *m;     // array of shape [number_terms, number_words]
//   uint64_t const *l;     // array of shape [number_terms, number_words]
//   uint64_t const *r;     // array of shape [number_terms, number_words]
//   uint64_t const *x;     // array of shape [number_terms, number_words]
//   uint64_t const *s;     // array of shape [number_terms, number_words]
//   // all arrays are contiguous in row-major order
// } ls_hs_nonbranching_terms;

#if 0
typedef struct ls_internal_operator_kernel_data {
  struct halide_buffer_t v;
  struct halide_buffer_t m;
  struct halide_buffer_t l;
  struct halide_buffer_t r;
  struct halide_buffer_t x;
  struct halide_buffer_t s;
  struct halide_dimension_t v_dims[2];
  struct halide_dimension_t m_dims[2];
  ls_internal_halide_apply_off_diag_kernel_type apply_off_diag_kernel;
  // ls_internal_halide_apply_diag_kernel_type apply_diag_kernel;
} ls_internal_operator_kernel_data;

static void
init_operator_kernel_data_buffers(ls_internal_operator_kernel_data *self,
                                  ls_hs_nonbranching_terms const *terms) {
  LS_CHECK(terms != NULL, "terms is NULL");
  self->v = (struct halide_buffer_t){.device = 0,
                                     .device_interface = NULL,
                                     .host = (uint8_t *)terms->v,
                                     .flags = 0,
                                     .type = {halide_type_float, 64, 1},
                                     .dimensions = 2,
                                     .dim = self->v_dims,
                                     .padding = NULL};
  self->m = (struct halide_buffer_t){.device = 0,
                                     .device_interface = NULL,
                                     .host = (uint8_t *)terms->m,
                                     .flags = 0,
                                     .type = {halide_type_uint, 64, 1},
                                     .dimensions = 2,
                                     .dim = self->m_dims,
                                     .padding = NULL};
  self->l = (struct halide_buffer_t){.device = 0,
                                     .device_interface = NULL,
                                     .host = (uint8_t *)terms->l,
                                     .flags = 0,
                                     .type = {halide_type_uint, 64, 1},
                                     .dimensions = 2,
                                     .dim = self->m_dims,
                                     .padding = NULL};
  self->r = (struct halide_buffer_t){.device = 0,
                                     .device_interface = NULL,
                                     .host = (uint8_t *)terms->r,
                                     .flags = 0,
                                     .type = {halide_type_uint, 64, 1},
                                     .dimensions = 2,
                                     .dim = self->m_dims,
                                     .padding = NULL};
  self->x = (struct halide_buffer_t){.device = 0,
                                     .device_interface = NULL,
                                     .host = (uint8_t *)terms->x,
                                     .flags = 0,
                                     .type = {halide_type_uint, 64, 1},
                                     .dimensions = 2,
                                     .dim = self->m_dims,
                                     .padding = NULL};
  self->s = (struct halide_buffer_t){.device = 0,
                                     .device_interface = NULL,
                                     .host = (uint8_t *)terms->s,
                                     .flags = 0,
                                     .type = {halide_type_uint, 64, 1},
                                     .dimensions = 2,
                                     .dim = self->m_dims,
                                     .padding = NULL};
  int32_t const number_terms = terms->number_terms;
  int32_t const number_words = (terms->number_bits + 63) / 64;

  self->v_dims[0] = (struct halide_dimension_t){
      .min = 0, .extent = 2, .stride = 1, .flags = 0};
  self->v_dims[1] = (struct halide_dimension_t){
      .min = 0, .extent = number_terms, .stride = 2, .flags = 0};

  self->m_dims[0] = (struct halide_dimension_t){
      .min = 0, .extent = number_words, .stride = 1, .flags = 0};
  self->m_dims[1] = (struct halide_dimension_t){
      .min = 0, .extent = number_terms, .stride = number_words, .flags = 0};
}

ls_internal_operator_kernel_data *ls_internal_create_apply_diag_kernel_data(
    ls_hs_nonbranching_terms const *diag_terms) {
  ls_internal_operator_kernel_data *self =
      malloc(sizeof(ls_internal_operator_kernel_data));
  LS_CHECK(self != NULL,
           "malloc failed when allocating ls_internal_operator_kernel_data");
  if (diag_terms == NULL) {
    // self->apply_diag_kernel = NULL;
    return self;
  }
  init_operator_kernel_data_buffers(self, diag_terms);
  ls_internal_halide_kernels_list list;
  init_kernels_list(diag_terms->number_bits, &list);
  // self->apply_diag_kernel = list.apply_diag;
  return self;
}

ls_internal_operator_kernel_data *ls_internal_create_apply_off_diag_kernel_data(
    ls_hs_nonbranching_terms const *off_diag_terms) {
  ls_internal_operator_kernel_data *self =
      malloc(sizeof(ls_internal_operator_kernel_data));
  LS_CHECK(self != NULL,
           "malloc failed when allocating ls_internal_operator_kernel_data");
  // If there are no off-diagonal terms, we don't have to initialize anything.
  // We just set the function pointer to NULL to indicate to the outside code
  // that nothing's been initialized.
  if (off_diag_terms == NULL) {
    self->apply_off_diag_kernel = NULL;
    return self;
  }
  init_operator_kernel_data_buffers(self, off_diag_terms);
  ls_internal_halide_kernels_list list;
  init_kernels_list(off_diag_terms->number_bits, &list);
  self->apply_off_diag_kernel = list.apply_off_diag;
  return self;
}

void ls_internal_destroy_operator_kernel_data(
    ls_internal_operator_kernel_data *p) {
  free(p);
}

void ls_internal_operator_apply_off_diag(
    ptrdiff_t batch_size, uint64_t const *alphas, ptrdiff_t const alphas_stride,
    uint64_t *betas, ptrdiff_t const betas_stride, ls_hs_scalar *coeffs,
    ls_internal_operator_kernel_data const *cxt) {
  LS_CHECK(cxt->apply_off_diag_kernel != NULL,
           "apply_off_diag has not been initialized");

  int32_t const number_terms = cxt->m_dims[1].extent;
  int32_t const number_words = cxt->m_dims[0].extent;
  // fprintf(stderr, "batch_size=%zi, number_terms=%i, number_words=%i\n",
  //         batch_size, number_terms, number_words);

  halide_dimension_t alphas_dims[2] = {
      (halide_dimension_t){
          .min = 0, .extent = number_words, .stride = 1, .flags = 0},
      (halide_dimension_t){.min = 0,
                           .extent = (int32_t)batch_size,
                           .stride = (int32_t)alphas_stride,
                           .flags = 0},
  };

  halide_buffer_t alphas_buf = (struct halide_buffer_t){
      .device = 0,
      .device_interface = NULL,
      .host = (uint8_t *)alphas,
      .flags = 0,
      .type = {halide_type_uint, 64, 1},
      .dimensions = 2,
      .dim = alphas_dims,
      .padding = NULL,
  };

  halide_dimension_t betas_dims[3] = {
      (halide_dimension_t){
          .min = 0, .extent = number_words, .stride = 1, .flags = 0},
      (halide_dimension_t){.min = 0,
                           .extent = number_terms,
                           .stride = (int32_t)betas_stride,
                           .flags = 0},
      (halide_dimension_t){.min = 0,
                           .extent = (int32_t)batch_size,
                           .stride = (int32_t)betas_stride * number_terms,
                           .flags = 0},
  };

  halide_buffer_t betas_buf = (struct halide_buffer_t){
      .device = 0,
      .device_interface = NULL,
      .host = (uint8_t *)betas,
      .flags = 0,
      .type = {halide_type_uint, 64, 1},
      .dimensions = 3,
      .dim = betas_dims,
      .padding = NULL,
  };

  halide_dimension_t coeffs_dims[3] = {
      (halide_dimension_t){.min = 0, .extent = 2, .stride = 1, .flags = 0},
      (halide_dimension_t){
          .min = 0, .extent = number_terms, .stride = 2, .flags = 0},
      (halide_dimension_t){.min = 0,
                           .extent = (int32_t)batch_size,
                           .stride = 2 * number_terms,
                           .flags = 0},
  };
  halide_buffer_t coeffs_buf = (struct halide_buffer_t){
      .device = 0,
      .device_interface = NULL,
      .host = (uint8_t *)coeffs,
      .flags = 0,
      .type = {halide_type_float, 64, 1},
      .dimensions = 3,
      .dim = coeffs_dims,
      .padding = NULL,
  };

  ls_internal_operator_kernel_data *const mut_cxt =
      (ls_internal_operator_kernel_data *)cxt;
  (*cxt->apply_off_diag_kernel)(&mut_cxt->v, &mut_cxt->m, &mut_cxt->r,
                                &mut_cxt->x, &mut_cxt->s, &alphas_buf,
                                &betas_buf, &coeffs_buf);
}
#endif
