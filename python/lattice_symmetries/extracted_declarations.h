typedef struct {
  void *elts;
  uint64_t num_elts;

  void *freer;
} chpl_external_array;
typedef struct ls_hs_scalar {
  double _real;
  double _imag;
} ls_hs_scalar;
void ls_hs_init(void);
void ls_hs_exit(void);

void ls_hs_set_exception_handler(void (*handler)(char const *message));
void ls_hs_error(char const *message);
typedef struct ls_hs_symmetry ls_hs_symmetry;

typedef struct ls_hs_symmetries ls_hs_symmetries;

typedef enum ls_hs_particle_type {
  LS_HS_SPIN,
  LS_HS_SPINFUL_FERMION,
  LS_HS_SPINLESS_FERMION
} ls_hs_particle_type;

typedef void (*ls_hs_internal_state_index_kernel_type)(
    ptrdiff_t batch_size, uint64_t const *alphas, ptrdiff_t alphas_stride,
    ptrdiff_t *indices, ptrdiff_t indices_stride, void const *private_data);

typedef void (*ls_hs_internal_is_representative_kernel_type)(
    ptrdiff_t batch_size, uint64_t const *alphas, ptrdiff_t alphas_stride,
    uint8_t *are_representatives, double *norms, void const *private_data);

typedef void (*ls_hs_internal_state_info_kernel_type)(
    ptrdiff_t batch_size, uint64_t const *alphas, ptrdiff_t alphas_stride,
    uint64_t *betas, ptrdiff_t betas_stride, ls_hs_scalar *characters,
    double *norms, void const *private_data);

typedef struct ls_hs_basis_kernels {
  ls_hs_internal_state_info_kernel_type state_info_kernel;
  void *state_info_data;
  ls_hs_internal_is_representative_kernel_type is_representative_kernel;
  void *is_representative_data;
  ls_hs_internal_state_index_kernel_type state_index_kernel;
  void *state_index_data;
} ls_hs_basis_kernels;

typedef struct ls_hs_permutation_group {
  _Atomic int refcount;
  int number_bits;
  int number_shifts;
  int number_masks;
  uint64_t *masks;
  uint64_t *shifts;
  double *eigvals_re;
  double *eigvals_im;
  void *haskell_payload;
} ls_hs_permutation_group;

typedef struct ls_hs_basis {
  _Atomic int refcount;
  int number_sites;
  int number_particles;
  int number_up;
  ls_hs_particle_type particle_type;
  int spin_inversion;
  bool state_index_is_identity;
  bool requires_projection;
  ls_hs_basis_kernels *kernels;
  chpl_external_array representatives;
  void *haskell_payload;
} ls_hs_basis;

typedef struct ls_hs_expr ls_hs_expr;

typedef void (*ls_hs_index_replacement_type)(int spin, int site, int *new_spin,
                                             int *new_site);

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
  _Atomic int refcount;
  ls_hs_basis const *basis;
  ls_hs_nonbranching_terms const *off_diag_terms;
  ls_hs_nonbranching_terms const *diag_terms;
  // ls_internal_operator_kernel_data const *apply_off_diag_cxt;
  // ls_internal_operator_kernel_data const *apply_diag_cxt;
  void *haskell_payload;
} ls_hs_operator;

typedef struct ls_hs_yaml_config {
  ls_hs_basis const *basis;
  ls_hs_operator const *hamiltonian;
  int number_observables;
  ls_hs_operator const *const *observables;
} ls_hs_yaml_config;

typedef struct ls_chpl_kernels {
  void (*enumerate_states)(ls_hs_basis const *, uint64_t, uint64_t,
                           chpl_external_array *);
  void (*operator_apply_off_diag)(ls_hs_operator *, int64_t, uint64_t *,
                                  chpl_external_array *, chpl_external_array *,
                                  chpl_external_array *, int64_t);
  void (*operator_apply_diag)(ls_hs_operator *, int64_t, uint64_t *,
                              chpl_external_array *, int64_t);
  void (*matrix_vector_product)(ls_hs_operator *, int, double const *,
                                double *);
} ls_chpl_kernels;
void ls_chpl_init(void);
void ls_chpl_finalize(void);
ls_chpl_kernels const *ls_hs_internal_get_chpl_kernels(void);
void ls_hs_internal_set_chpl_kernels(ls_chpl_kernels const *kernels);
void ls_hs_internal_destroy_external_array(chpl_external_array *arr);
// The following should, in principle, not exist
void ls_hs_state_index(ls_hs_basis const *const basis,
                       ptrdiff_t const batch_size,
                       uint64_t const *const restrict spins,
                       ptrdiff_t const spins_stride,
                       ptrdiff_t *const restrict indices,
                       ptrdiff_t const indices_stride);
void ls_hs_is_representative(ls_hs_basis const *basis, ptrdiff_t batch_size,
                             uint64_t const *alphas, ptrdiff_t alphas_stride,
                             uint8_t *are_representatives, double *norms);
void ls_hs_state_info(ls_hs_basis const *basis, ptrdiff_t batch_size,
                      uint64_t const *alphas, ptrdiff_t alphas_stride,
                      uint64_t *betas, ptrdiff_t betas_stride,
                      ls_hs_scalar *characters, double *norms);
void ls_hs_build_representatives(ls_hs_basis *basis, uint64_t const lower,
                                 uint64_t const upper);
void ls_hs_unchecked_set_representatives(ls_hs_basis *basis,
                                         chpl_external_array const *states,
					 int cache_bits);
ls_hs_symmetry * ls_hs_symmetry_from_json(char const *);
void ls_hs_destroy_symmetry(ls_hs_symmetry *);
int ls_hs_symmetry_sector(ls_hs_symmetry const*);
double ls_hs_symmetry_phase(ls_hs_symmetry const*);
int ls_hs_symmetry_length(ls_hs_symmetry const*);
int * ls_hs_symmetry_permutation(ls_hs_symmetry const*);
void ls_hs_destroy_permutation(int *);
ls_hs_symmetries * ls_hs_symmetries_from_json(char const *);
void ls_hs_destroy_symmetries(ls_hs_symmetries *);
ls_hs_basis * ls_hs_clone_basis(ls_hs_basis const*);
void ls_hs_destroy_basis(ls_hs_basis *);
ls_hs_basis * ls_hs_basis_from_json(char const *);
char const * ls_hs_basis_to_json(ls_hs_basis const*);
void ls_hs_destroy_string(char const *);
uint64_t ls_hs_min_state_estimate(ls_hs_basis const*);
uint64_t ls_hs_max_state_estimate(ls_hs_basis const*);
bool ls_hs_basis_has_fixed_hamming_weight(ls_hs_basis const*);
bool ls_hs_basis_has_spin_inversion_symmetry(ls_hs_basis const*);
bool ls_hs_basis_has_permutation_symmetries(ls_hs_basis const*);
bool ls_hs_basis_requires_projection(ls_hs_basis const*);
void ls_hs_basis_build(ls_hs_basis const*);
bool ls_hs_basis_is_built(ls_hs_basis const*);
int ls_hs_basis_number_bits(ls_hs_basis const*);
int ls_hs_basis_number_words(ls_hs_basis const*);
char const * ls_hs_basis_state_to_string(ls_hs_basis const*, uint64_t const*);
ptrdiff_t ls_hs_fixed_hamming_state_to_index(uint64_t);
uint64_t ls_hs_fixed_hamming_index_to_state(ptrdiff_t, int);
char const * ls_hs_expr_to_json(ls_hs_expr const*);
ls_hs_expr * ls_hs_expr_from_json(char const *);
void ls_hs_destroy_expr(ls_hs_expr *);
char const * ls_hs_expr_to_string(ls_hs_expr const*);
ls_hs_expr * ls_hs_expr_plus(ls_hs_expr const*, ls_hs_expr const*);
ls_hs_expr * ls_hs_expr_minus(ls_hs_expr const*, ls_hs_expr const*);
ls_hs_expr * ls_hs_expr_times(ls_hs_expr const*, ls_hs_expr const*);
ls_hs_expr * ls_hs_expr_scale(ls_hs_scalar const*, ls_hs_expr const*);
ls_hs_expr * ls_hs_replace_indices(ls_hs_expr const*, ls_hs_index_replacement_type);
bool ls_hs_expr_equal(ls_hs_expr const*, ls_hs_expr const*);
ls_hs_expr * ls_hs_expr_adjoint(ls_hs_expr const*);
bool ls_hs_expr_is_hermitian(ls_hs_expr const*);
bool ls_hs_expr_is_real(ls_hs_expr const*);
bool ls_hs_expr_is_identity(ls_hs_expr const*);
ls_hs_operator * ls_hs_create_operator(ls_hs_basis const*, ls_hs_expr const*);
ls_hs_operator * ls_hs_clone_operator(ls_hs_operator const*);
void ls_hs_destroy_operator(ls_hs_operator *);
int ls_hs_operator_max_number_off_diag(ls_hs_operator *);
ls_hs_expr * ls_hs_operator_get_expr(ls_hs_operator const*);
ls_hs_basis * ls_hs_operator_get_basis(ls_hs_operator const*);
void ls_hs_prepare_hphi(ls_hs_operator const*, char const *);
void ls_hs_prepare_mvmc(ls_hs_operator const*, char const *);
ls_hs_yaml_config * ls_hs_load_yaml_config(char const *);
void ls_hs_destroy_yaml_config(ls_hs_yaml_config *);
