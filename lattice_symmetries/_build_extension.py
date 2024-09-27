from cffi import FFI


ffibuilder = FFI()

ffibuilder.cdef(
    """
typedef struct ls_numpy_array_1d { uint8_t* data; void* handle; } ls_numpy_array_1d;
typedef void (*ls_alloc_numpy_array_1d_callback)(uint64_t, ls_numpy_array_1d*);

extern "Python" void ls_alloc_numpy_array_1d(uint64_t size, ls_numpy_array_1d* out);

void hello(void);

void ls_chpl_init(void);
void chpl_library_finalize(void);
"""
)

ffibuilder.set_source(
    "lattice_symmetries._ls",
    """
#include <force_link_glibc_2.27.h>
#include <lattice_symmetries.h>
#include <lattice_symmetries_chapel.h>

void ls_chpl_init(void) {
  int const argc = 1;
  char const *argv[2] = {"lattice_symmetries", NULL};
  chpl_library_init(argc, (char**)argv);
  chpl__init_FFI(1, 2);
  chpl__init_Library(1, 2);
}
""",
    extra_compile_args=["-Wall", "-Wextra"],
    include_dirs=["include/"],
    libraries=["lattice_symmetries_chapel"],
    library_dirs=["."],
    extra_link_args=["-Wl,-rpath=$ORIGIN"],
)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
