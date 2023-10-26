from cffi import FFI
import glob
import os

# def get_chapel_lib_path(prefix="lattice_symmetries"):
#     print(os.getcwd())
#     print(os.listdir(prefix))
#     for f in glob.glob(os.path.join(prefix, "lattice-symmetries-chapel*")):
#         print("get_chapel_lib_path -> {}".format(f[len(prefix) + 1:]))
#         return f[len(prefix) + 1:]
#     return os.environ.get("LATTICE_SYMMETRIES_PATH")


def get_include_dirs():
    return []
    # kernels_path = os.environ.get("LS_KERNELS_PATH")
    # if kernels_path is not None:
    #     return [os.path.join(kernels_path, "include")]
    # else:
    #     return []

    # prefix = get_chapel_lib_path()
    # if prefix is None:
    #     assert False
    #     return []
    # else:
    #     return [os.path.join("lattice_symmetries", prefix, "include")]


def get_library_dirs():
    return []
    # prefix = get_chapel_lib_path()
    # if prefix is None:
    #     assert False
    #     return []
    # else:
    #     return [os.path.join("lattice_symmetries", prefix, "lib")]


def get_runtime_dirs():
    return []
    # prefix = get_chapel_lib_path()
    # if prefix is None:
    #     assert False
    #     return []
    # else:
    #     return [os.path.join("$ORIGIN", prefix, "lib")]


def get_generated_declarations():
    with open("lattice_symmetries/extracted_declarations.h", "r") as f:
        return f.read()


ffibuilder = FFI()

ffibuilder.cdef(
    get_generated_declarations()
    + "\n"
    + """
void set_python_exception_handler();
ptrdiff_t ls_hs_basis_number_states(ls_hs_basis const* basis);
uint64_t const* ls_hs_basis_states(ls_hs_basis const* basis);
extern "Python" void python_replace_indices(int, int, int*, int*);
extern "Python" void python_process_symmetries(ls_hs_symmetries*);
"""
)

ffibuilder.set_source(
    "lattice_symmetries._ls_hs",
    """
#define LS_NO_STD_COMPLEX
#include "lattice_symmetries_types.h"
#include "lattice_symmetries_functions.h"
#include <Python.h>

static void python_exception_handler(char const* message) {
    PyGILState_STATE gstate = PyGILState_Ensure();
    // fprintf(stderr, "Calling PyErr_SetString '%s' ...\\n", message);
    PyErr_SetString(PyExc_RuntimeError, message);
    // fprintf(stderr, "Returning ...\\n");
    PyGILState_Release(gstate);
}

static void set_python_exception_handler() {
    // fprintf(stderr, "set_python_exception_handler ...\\n");
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
