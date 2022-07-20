from cffi import FFI

ffibuilder = FFI()

ffibuilder.cdef(
    """
    void set_python_exception_handler();

    void ls_hs_init(void);
    void ls_hs_exit(void);

    void ls_hs_set_exception_handler(void (*handler)(char const* message));
    void ls_hs_error(char const* message);

    typedef struct ls_hs_basis { ...; } ls_hs_basis;
    typedef struct ls_hs_operator { ...; } ls_hs_operator;
    typedef struct ls_hs_scalar { double real; double imag; } ls_hs_scalar;

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

    ls_hs_basis *ls_hs_clone_basis(ls_hs_basis const *);
    void ls_hs_destroy_basis_v2(ls_hs_basis *);
    uint64_t ls_hs_max_state_estimate(ls_hs_basis const *);
    uint64_t ls_hs_min_state_estimate(ls_hs_basis const *);
    ls_hs_basis *ls_hs_basis_from_json(char const *json_string);
    char const *ls_hs_basis_to_json(ls_hs_basis const *);
    void ls_hs_basis_build(ls_hs_basis *basis);
    bool ls_hs_basis_has_fixed_hamming_weight(ls_hs_basis const *);
    void ls_hs_state_index(ls_hs_basis const *basis, ptrdiff_t batch_size,
                           uint64_t const *spins, ptrdiff_t spins_stride,
                           ptrdiff_t *indices, ptrdiff_t indices_stride);

    int ls_hs_basis_number_bits(ls_hs_basis const* basis);
    int ls_hs_basis_number_words(ls_hs_basis const* basis);

    ptrdiff_t ls_hs_basis_number_states(ls_hs_basis const* basis);
    uint64_t const* ls_hs_basis_states(ls_hs_basis const* basis);

    ls_hs_operator *ls_hs_create_operator(ls_hs_basis const *, char const *, int, int, int const *);
    ls_hs_operator *ls_hs_operator_plus(ls_hs_operator const *, ls_hs_operator const *);
    ls_hs_operator *ls_hs_operator_minus(ls_hs_operator const *, ls_hs_operator const *);
    ls_hs_operator *ls_hs_operator_times(ls_hs_operator const *, ls_hs_operator const *);
    ls_hs_operator *ls_hs_operator_scale(ls_hs_scalar const *, ls_hs_operator const *);
    ls_hs_operator *ls_hs_operator_hermitian_conjugate(ls_hs_operator const *);

    bool ls_hs_operator_is_hermitian(ls_hs_operator const *);
    bool ls_hs_operator_is_identity(ls_hs_operator const *);

    char const *ls_hs_operator_pretty_terms(ls_hs_operator const *);
    void ls_hs_destroy_string(char const *);

    void ls_hs_destroy_operator_v2(ls_hs_operator *);

    char const *ls_hs_basis_state_to_string(ls_hs_basis const *,
                                            uint64_t const *state);

    void chpl_library_init_wrapper(void);
    void chpl_library_finalize(void);
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

    extern void chpl__init_StatesEnumeration(int64_t _ln, int32_t _fn);
    extern void chpl__init_BatchedOperator(int64_t _ln, int32_t _fn);
    extern void chpl__init_CommunicationQueue(int64_t _ln,
                                       int32_t _fn);
    extern void chpl__init_ConcurrentAccessor(int64_t _ln, int32_t _fn);
    extern void chpl__init_DistributedMatrixVector(int64_t _ln, int32_t _fn);
    extern void chpl__init_FFI(int64_t _ln, int32_t _fn);
    extern void chpl__init_ForeignTypes(int64_t _ln, int32_t _fn);
    extern void chpl__init_LatticeSymmetries(int64_t _ln, int32_t _fn);
    extern void chpl__init_MatrixVectorProduct(int64_t _ln, int32_t _fn);
    extern void chpl__init_StatesEnumeration(int64_t _ln, int32_t _fn);
    extern void chpl__init_Vector(int64_t _ln, int32_t _fn);

    extern void chpl_library_init(int argc, char* argv[]);
    extern void chpl_library_finalize(void);
    extern void ls_chpl_init_kernels(void);

    static void chpl_library_init_wrapper(void) {
        int const argc = 1;
        char* argv[1] = {"lattice_symmetries"};
        chpl_library_init(argc, argv);
        chpl__init_BatchedOperator(1, 2);
        chpl__init_CommunicationQueue(1, 2);
        chpl__init_ConcurrentAccessor(1, 2);
        chpl__init_DistributedMatrixVector(1, 2);
        chpl__init_FFI(1, 2);
        chpl__init_ForeignTypes(1, 2);
        chpl__init_LatticeSymmetries(1, 2);
        chpl__init_MatrixVectorProduct(1, 2);
        chpl__init_StatesEnumeration(1, 2);
        chpl__init_Vector(1, 2);
    }
""",
    include_dirs=["../cbits", "/home/tom/src/distributed-matvec/lib/lattice_symmetries_chapel.h"],
    extra_compile_args=["-Wall", "-Wextra"],
    libraries=["lattice_symmetries_chapel", "lattice_symmetries_haskell", "lattice_symmetries_core"],
    library_dirs=["../dist-newstyle/build/x86_64-linux/ghc-8.10.7/lattice-symmetries-haskell-0.1.0.0/f/lattice_symmetries_haskell/noopt/build/lattice_symmetries_haskell",
                  "../kernels/build",
                  "/home/tom/src/distributed-matvec/lib"],
)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)

