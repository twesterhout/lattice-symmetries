#pragma once

#if defined(__cplusplus)
#include <complex>
#include <cstddef>
#include <cstdint>
#else
#include <stddef.h>
#include <stdint.h>
#endif

#include <lattice_symmetries/lattice_symmetries.h>

#if defined(__cplusplus)
extern "C" {
#endif

void ls_hs_init(void);
void ls_hs_exit(void);

typedef struct ls_hs_spin_basis_v1 {
  ls_spin_basis const *const payload;
  void *const context;
} ls_hs_spin_basis_v1;

typedef struct ls_hs_operator_v1 {
  ls_operator const *const payload;
  void *const context;
} ls_hs_operator_v1;

void ls_hs_basis_and_hamiltonian_from_yaml(char const *path,
                                           ls_hs_spin_basis_v1 *basis,
                                           ls_hs_operator_v1 *hamiltonian);
void ls_hs_destroy_spin_basis(ls_hs_spin_basis_v1 *basis);
void ls_hs_destroy_operator(ls_hs_operator_v1 *op);

ls_error_code ls_hs_operator_apply(ls_hs_operator_v1 const *op, uint64_t count,
                                   uint64_t const *spins, uint64_t *offsets,
                                   uint64_t *out_spins,
                                   _Complex double *out_coeffs);

unsigned ls_hs_hdf5_get_dataset_rank(char const *filename, char const *dataset);
void ls_hs_hdf5_get_dataset_shape(char const *filename, char const *dataset,
                                  uint64_t *shape);
void ls_hs_hdf5_create_dataset_u64(char const *filename, char const *dataset,
                                   unsigned dim, uint64_t const *shape);
void ls_hs_hdf5_create_dataset_f64(char const *filename, char const *dataset,
                                   unsigned dim, uint64_t const *shape);
void ls_hs_hdf5_create_dataset_f32(char const *filename, char const *dataset,
                                   unsigned dim, uint64_t const *shape);
void ls_hs_hdf5_create_dataset_c64(char const *filename, char const *dataset,
                                   unsigned dim, uint64_t const *shape);
void ls_hs_hdf5_create_dataset_c128(char const *filename, char const *dataset,
                                    unsigned dim, uint64_t const *shape);
void ls_hs_hdf5_write_chunk_u64(char const *filename, char const *dataset,
                                unsigned dim, uint64_t const *offset,
                                uint64_t const *shape, uint64_t const *data);
void ls_hs_hdf5_write_chunk_f64(char const *filename, char const *dataset,
                                unsigned dim, uint64_t const *offset,
                                uint64_t const *shape, double const *data);
void ls_hs_hdf5_read_chunk_u64(char const *filename, char const *dataset,
                               unsigned dim, uint64_t const *offset,
                               uint64_t const *shape, uint64_t *data);
void ls_hs_hdf5_read_chunk_f64(char const *filename, char const *dataset,
                               unsigned dim, uint64_t const *offset,
                               uint64_t const *shape, double *data);

typedef enum ls_hs_particle_type {
  LS_HS_SPIN,
  LS_HS_FERMION
} ls_hs_particle_type;

typedef struct ls_hs_basis {
  ls_hs_particle_type particle_type;
  int number_sites; // always >= 0
  int number_up;    // either >= 0 or == -1, -1 means "unspecified"
  int number_down;  // either >= 0 or == -1, -1 means "unspecified"
  bool state_index_is_identity;
  // if number_up >= 0 and number_down >= 0
  //   then
  //     number_up + number_down <= number_sites
  //        when particle_type == LS_HS_FERMION
  //     number_up + number_down == number_sites
  //        when particle_type == LS_HS_SPIN
} ls_hs_basis;

#if defined(__cplusplus)
using ls_hs_scalar = std::complex<double>;
#else
typedef _Complex double ls_hs_scalar;
#endif

typedef struct ls_hs_nonbranching_terms {
  int number_terms;
  int number_bits;
  // number_words = ceil(number_bits / 64)
  ls_hs_scalar const *v; // array of shape [number_terms]
  uint64_t const *m;     // array of shape [number_terms, number_words]
  uint64_t const *l;     // array of shape [number_terms, number_words]
  uint64_t const *r;     // array of shape [number_terms, number_words]
  uint64_t const *x;     // array of shape [number_terms, number_words]
  uint64_t const *s;     // array of shape [number_terms, number_words]
  // all arrays are contiguous in row-major order
} ls_hs_nonbranching_terms;

typedef struct ls_hs_operator {
  ls_hs_basis const *basis;
  ls_hs_nonbranching_terms const *off_diag_terms;
  ls_hs_nonbranching_terms const *diag_terms;
  bool needs_projection;
} ls_hs_operator;

// {{{ Binomials
typedef struct ls_hs_binomials {
  int dimension;
  uint64_t *coefficients; // array of shape [dimension, dimension]
                          // stored in row-major order
  // coefficients[n, k] corresponds to binomial(n, k)
} ls_hs_binomials;

ls_hs_binomials *ls_hs_internal_malloc_binomials(int number_bits);
void ls_hs_internal_free_binomials(ls_hs_binomials *p);
void ls_hs_internal_compute_binomials(ls_hs_binomials *p);
uint64_t ls_hs_internal_binomial(int n, int k, ls_hs_binomials const *cache);

void ls_hs_state_index_combinadics_kernel(ptrdiff_t batch_size,
                                          uint64_t const *spins,
                                          ptrdiff_t spins_stride,
                                          ptrdiff_t *indices,
                                          ptrdiff_t indices_stride,
                                          void const *private_kernel_data);
void ls_hs_state_index_identity_kernel(ptrdiff_t batch_size,
                                       uint64_t const *spins,
                                       ptrdiff_t spins_stride,
                                       ptrdiff_t *indices,
                                       ptrdiff_t indices_stride,
                                       void const *private_kernel_data);
// }}}

#if 0
typedef struct ls_sparse_operator ls_sparse_operator;
typedef struct ls_term ls_term;

typedef struct ls_hs_basis {
  ls_flat_spin_basis const *payload;
  void *context;
} ls_hs_basis;

typedef struct ls_hs_operator {
  ls_sparse_operator const *payload;
  void *context;
} ls_hs_operator;

typedef void (*ls_spin_copy_fn)(uint64_t const * /*source*/,
                                uint64_t * /*destination*/);
typedef void (*ls_spin_fill_fn)(uint64_t const * /*2D source*/,
                                uint64_t /*count*/,
                                uint64_t * /*2D destination*/);

typedef struct ls_output_buffer {
  uint64_t *spins;
  _Complex double *coeffs;
  _Complex double *const diagonal;
  uint64_t const number_words;
  ls_spin_copy_fn const spin_copy;
  ls_spin_fill_fn const spin_fill;
} ls_output_buffer;

typedef struct ls_workspace {
  uint64_t *spins;
  _Complex double *characters;
  double *norms;
} ls_workspace;

void ls_hs_create_operator(ls_hs_operator *self,
                           ls_flat_spin_basis const *basis,
                           unsigned number_terms, ls_term const *term);
void ls_hs_destroy_operator(ls_hs_operator *self);
void ls_hs_load_hamiltonian_from_yaml(ls_hs_operator *self, char const *path);
void ls_hs_operator_get_basis(ls_hs_operator const *self, ls_hs_basis *basis);

uint64_t ls_hs_apply_operator(ls_hs_operator const *term, uint64_t const *spin,
                              ls_output_buffer *out, ls_workspace *workspace);
#endif

#ifdef __cplusplus
} // extern "C"
#endif
