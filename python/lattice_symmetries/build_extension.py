from cffi import FFI

ffibuilder = FFI()

ffibuilder.cdef(
    """
    void set_python_exception_handler();

    void ls_hs_init(void);
    void ls_hs_exit(void);

    void ls_hs_set_exception_handler(void (*handler)(char const* message));
    void ls_hs_error(char const* message);

    typedef struct ls_hs_basis ls_hs_basis;

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
"""
)

ffibuilder.set_source(
    "lattice_symmetries._ls_hs",
    """
    #include "lattice_symmetries_haskell.h"
    #include <Python.h>

    static void python_exception_handler(char const* message) {
      PyGILState_STATE gstate = PyGILState_Ensure();
      fprintf(stderr, "Calling PyErr_SetString '%s' ...\\n", message);
      PyErr_SetString(PyExc_RuntimeError, message);
      fprintf(stderr, "Returning ...\\n");
      PyGILState_Release(gstate);
    }

    static void set_python_exception_handler() {
      fprintf(stderr, "set_python_exception_handler ...\\n");
      ls_hs_set_exception_handler(&python_exception_handler);
    }
""",
    include_dirs=["../cbits"],
    libraries=["lattice_symmetries_haskell", "lattice_symmetries_core"],
    library_dirs=["../dist-newstyle/build/x86_64-linux/ghc-8.10.7/lattice-symmetries-haskell-0.1.0.0/f/lattice_symmetries_haskell/noopt/build/lattice_symmetries_haskell", "../kernels/build"],
)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)

