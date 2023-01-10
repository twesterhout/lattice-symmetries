from cffi import FFI
import glob
import os

def get_chapel_lib_path(prefix="lattice_symmetries"):
    print(os.getcwd())
    print(os.listdir(prefix))
    for f in glob.glob(os.path.join(prefix, "lattice-symmetries-chapel-*")):
        print("get_chapel_lib_path -> {}".format(f[len(prefix) + 1:]))
        return f[len(prefix) + 1:]
    return os.environ.get("LATTICE_SYMMETRIES_PATH")

def get_include_dirs():
    prefix = get_chapel_lib_path()
    if prefix is None:
        assert False
        return []
    else:
        return [os.path.join("lattice_symmetries", prefix, "include")]

def get_library_dirs():
    prefix = get_chapel_lib_path()
    if prefix is None:
        assert False
        return []
    else:
        return [os.path.join("lattice_symmetries", prefix, "lib")]

def get_runtime_dirs():
    prefix = get_chapel_lib_path()
    if prefix is None:
        assert False
        return []
    else:
        return [os.path.join("$ORIGIN", prefix, "lib")]

ffibuilder = FFI()

ffibuilder.cdef(
    """
    void set_python_exception_handler();

    void ls_hs_init(void);
    void ls_hs_exit(void);

    void ls_hs_set_exception_handler(void (*handler)(char const* message));
    void ls_hs_error(char const *message);
    void ls_hs_destroy_string(char *str);

    typedef struct ls_hs_scalar { double real; double imag; } ls_hs_scalar;
    typedef struct ls_hs_symmetries ls_hs_symmetries;
    typedef struct ls_hs_basis { int spin_inversion; ...; } ls_hs_basis;
    typedef struct ls_hs_expr ls_hs_expr;
    typedef struct ls_hs_operator ls_hs_operator;

    typedef struct {
      void *elts;
      uint64_t num_elts;

      void *freer;
    } chpl_external_array;

    void ls_hs_destroy_external_array(chpl_external_array *arr);

    typedef struct ls_chpl_kernels {
      void (*enumerate_states)(ls_hs_basis const *, uint64_t, uint64_t,
                               chpl_external_array *);
      void (*operator_apply_off_diag)(ls_hs_operator *, int64_t, uint64_t *,
                                      chpl_external_array *, chpl_external_array *,
                                      chpl_external_array *, int64_t);
      void (*operator_apply_diag)(ls_hs_operator *, int64_t, uint64_t *,
                                  chpl_external_array *, int64_t);
      void (*matrix_vector_product)(ls_hs_operator *, int, double const *, double *);
    } ls_chpl_kernels;
    ls_chpl_kernels const *ls_hs_internal_get_chpl_kernels();

    typedef struct ls_hs_symmetry ls_hs_symmetry;
    ls_hs_symmetry *ls_hs_symmetry_from_json(const char *json_string);
    void ls_hs_destroy_symmetry(ls_hs_symmetry *);

    int ls_hs_symmetry_sector(ls_hs_symmetry const *);
    int ls_hs_symmetry_length(ls_hs_symmetry const *);
    int *ls_hs_symmetry_permutation(ls_hs_symmetry const *);
    void ls_hs_destroy_permutation(int *);

    ls_hs_symmetries *ls_hs_symmetries_from_json(const char *json_string);
    void ls_hs_destroy_symmetries(ls_hs_symmetries *);

    ls_hs_basis *ls_hs_basis_from_json(char const *json_string);
    ls_hs_basis *ls_hs_clone_basis(ls_hs_basis const *);
    void ls_hs_destroy_basis(ls_hs_basis *);
    char const *ls_hs_basis_to_json(ls_hs_basis const *);

    uint64_t ls_hs_max_state_estimate(ls_hs_basis const *);
    uint64_t ls_hs_min_state_estimate(ls_hs_basis const *);
    int ls_hs_basis_number_bits(ls_hs_basis const* basis);
    int ls_hs_basis_number_words(ls_hs_basis const* basis);
    bool ls_hs_basis_has_fixed_hamming_weight(ls_hs_basis const *);
    bool ls_hs_basis_has_spin_inversion_symmetry(ls_hs_basis const *);
    bool ls_hs_basis_has_permutation_symmetries(ls_hs_basis const *);
    bool ls_hs_basis_requires_projection(ls_hs_basis const *);

    char const *ls_hs_basis_state_to_string(ls_hs_basis const *,
                                            uint64_t const *state);

    void ls_hs_basis_build(ls_hs_basis *basis);
    bool ls_hs_basis_is_built(ls_hs_basis const *basis);

    ptrdiff_t ls_hs_basis_number_states(ls_hs_basis const* basis);
    uint64_t const* ls_hs_basis_states(ls_hs_basis const* basis);

    void ls_hs_state_index(ls_hs_basis const *basis, ptrdiff_t batch_size,
                           uint64_t const *spins, ptrdiff_t spins_stride,
                           ptrdiff_t *indices, ptrdiff_t indices_stride);

    void ls_hs_state_info(ls_hs_basis const *basis, ptrdiff_t batch_size,
                          uint64_t const *alphas, ptrdiff_t alphas_stride,
                          uint64_t *betas, ptrdiff_t betas_stride,
                          ls_hs_scalar *characters, double *norms);

    typedef enum ls_hs_particle_type {
      LS_HS_SPIN,
      LS_HS_SPINFUL_FERMION,
      LS_HS_SPINLESS_FERMION
    } ls_hs_particle_type;

    typedef void (*ls_hs_index_replacement_type)(int spin, int site, int *new_spin,
                                                 int *new_site);

    extern "Python" void python_replace_indices(int, int, int*, int*);

    char const *ls_hs_expr_to_json(ls_hs_expr const *expr);
    ls_hs_expr *ls_hs_expr_from_json(char const *json_string);
    ls_hs_expr *ls_hs_replace_indices(ls_hs_expr const *expr,
                                      ls_hs_index_replacement_type callback);
    char const *ls_hs_expr_to_string(ls_hs_expr const *expr);
    void ls_hs_destroy_expr(ls_hs_expr *expr);

    ls_hs_expr *ls_hs_expr_plus(ls_hs_expr const *a, ls_hs_expr const *b);
    ls_hs_expr *ls_hs_expr_minus(ls_hs_expr const *a, ls_hs_expr const *b);
    ls_hs_expr *ls_hs_expr_times(ls_hs_expr const *a, ls_hs_expr const *b);
    ls_hs_expr *ls_hs_expr_scale(ls_hs_scalar const *z, ls_hs_expr const *a);
    ls_hs_expr *ls_hs_expr_adjoint(ls_hs_expr const *);

    bool ls_hs_expr_equal(ls_hs_expr const *a, ls_hs_expr const *b);
    bool ls_hs_expr_is_hermitian(ls_hs_expr const *);
    bool ls_hs_expr_is_real(ls_hs_expr const *);
    bool ls_hs_expr_is_identity(ls_hs_expr const *);

    ls_hs_operator *ls_hs_create_operator(ls_hs_basis const *, ls_hs_expr const *);
    ls_hs_operator *ls_hs_clone_operator(ls_hs_operator const *);
    void ls_hs_destroy_operator(ls_hs_operator *);
    ls_hs_expr const *ls_hs_operator_get_expr(ls_hs_operator const *);
    ls_hs_basis const *ls_hs_operator_get_basis(ls_hs_operator const *);

    typedef struct ls_hs_yaml_config {
      ls_hs_basis const *basis;
      ls_hs_operator const *hamiltonian;
      int number_observables;
      ls_hs_operator const *const *observables;
    } ls_hs_yaml_config;

    ls_hs_yaml_config *ls_hs_load_yaml_config(char const *);
    void ls_hs_destroy_yaml_config(ls_hs_yaml_config *);

    void ls_chpl_init(void);
    void ls_chpl_finalize(void);
"""
)

ffibuilder.set_source(
    "lattice_symmetries._ls_hs",
    """
    #define LS_NO_STD_COMPLEX
    #include "lattice_symmetries_haskell.h"
    #include <Python.h>

    static void python_exception_handler(char const* message) {
      PyGILState_STATE gstate = PyGILState_Ensure();
      // fprintf(stderr, "Calling PyErr_SetString '%s' ...\\n", message);
      PyErr_SetString(PyExc_RuntimeError, message);
      // fprintf(stderr, "Returning ...\\n");
      PyGILState_Release(gstate);
    }

    static void set_python_exception_handler() {
      fprintf(stderr, "set_python_exception_handler ...\\n");
      ls_hs_set_exception_handler(&python_exception_handler);
    }

    static ptrdiff_t ls_hs_basis_number_states(ls_hs_basis const* basis) {
      if (basis->representatives.elts == NULL) { return -1; }
      return (ptrdiff_t)basis->representatives.num_elts;
    }
    static uint64_t const* ls_hs_basis_states(ls_hs_basis const* basis) {
      return basis->representatives.elts;
    }
""",
    include_dirs=get_include_dirs(),
    # [# "../cbits", "/home/tom/src/distributed-matvec/lib/lattice_symmetries_chapel.h"],
    #               "/home/tom/src/distributed-matvec/third_party/include"],
    extra_compile_args=["-Wall", "-Wextra"],
    libraries=["lattice_symmetries_chapel", "lattice_symmetries_haskell"],
    library_dirs=get_library_dirs(),
    extra_link_args=["-Wl,-rpath,{}".format(d) for d in get_runtime_dirs()],
    # [
    #               "/home/tom/src/distributed-matvec/third_party/lib",
    #               # "../dist-newstyle/build/x86_64-linux/ghc-8.10.7/lattice-symmetries-haskell-0.1.0.0/f/lattice_symmetries_haskell/noopt/build/lattice_symmetries_haskell",
    #               # "../kernels/build",
    #               "/home/tom/src/distributed-matvec/lib"],
)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)

