from cffi import FFI


def get_generated_declarations():
    with open("lattice_symmetries/extracted_declarations.h", "r") as f:
        return f.read().replace("LS_CONST", "const")


ffibuilder = FFI()

ffibuilder.cdef(get_generated_declarations())

ffibuilder.set_source(
    "lattice_symmetries._ls_hs",
    """
#define LS_NO_STD_COMPLEX
#include "lattice_symmetries_types.h"
#include "lattice_symmetries_functions.h"
#include "lattice_symmetries_chapel.h"
""",
    extra_compile_args=["-Wall", "-Wextra", "-DPYTHON_CFFI=1"],  # , "-fsanitize=address"],
    # extra_link_args=["-fsanitize=address"],
    libraries=["lattice_symmetries_chapel", "lattice_symmetries_haskell"],
)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
