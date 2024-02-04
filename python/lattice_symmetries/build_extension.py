from cffi import FFI


def get_generated_declarations():
    with open("lattice_symmetries/extracted_declarations.h", "r") as f:
        return f.read()


ffibuilder = FFI()

ffibuilder.cdef(
    get_generated_declarations()
    + "\n"
    + """
void set_python_exception_handler();
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
    PyErr_SetString(PyExc_RuntimeError, message);
    PyGILState_Release(gstate);
}

static void set_python_exception_handler() {
    ls_hs_set_exception_handler(&python_exception_handler);
}
""",
    extra_compile_args=["-Wall", "-Wextra", "-DPYTHON_CFFI=1"],
    libraries=[
        # "lattice_symmetries_chapel",
        "lattice_symmetries_haskell"
    ],
)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
