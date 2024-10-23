from cffi import FFI


ffibuilder = FFI()

ffibuilder.cdef(
    """
typedef struct ls_numpy_array_1d { ...; } ls_numpy_array_1d;
typedef struct _complex64 { ...; } _complex64;
typedef struct _complex128 { ...; } _complex128;

typedef struct ls_diag_terms {
    void           *kernel;
    int             number_terms;
    double const   *v_re;
    double const   *v_im;
} ls_diag_terms;

typedef struct ls_off_diag_terms {
    void           *kernel;
    int             number_terms;
    int             number_reduced;
    double const   *v_re;
    double const   *v_im;
    uint64_t const *x;
} ls_off_diag_terms;

typedef void (*ls_alloc_numpy_array_1d_callback)(uint64_t, ls_numpy_array_1d*);

extern "Python" void ls_alloc_numpy_array_1d(uint64_t size, ls_numpy_array_1d* out);

void hello(void);
void the_ultimate_solution(void* alloc_numpy_array_1d,
                           ls_numpy_array_1d * result);

void ls_enumerate_states_fixed_hamming(int64_t numChunks,
                                       const int64_t * offsets,
                                       const uint64_t * values,
                                       uint64_t * dest);

void ls_matrix_apply_f32(ls_diag_terms * diag,
                         ls_off_diag_terms * off_diag,
                         int64_t numStates,
                         const uint64_t * representativesPtr,
                         int64_t numVectors,
                         const float * xPtr,
                         float * yPtr);
void ls_matrix_apply_f64(ls_diag_terms * diag,
                         ls_off_diag_terms * off_diag,
                         int64_t numStates,
                         const uint64_t * representativesPtr,
                         int64_t numVectors,
                         const double * xPtr,
                         double * yPtr);
void ls_matrix_apply_c64(ls_diag_terms * diag,
                         ls_off_diag_terms * off_diag,
                         int64_t numStates,
                         const uint64_t * representativesPtr,
                         int64_t numVectors,
                         const _complex64 * xPtr,
                         _complex64 * yPtr);
void ls_matrix_apply_c128(ls_diag_terms * diag,
                          ls_off_diag_terms * off_diag,
                          int64_t numStates,
                          const uint64_t * representativesPtr,
                          int64_t numVectors,
                          const _complex128 * xPtr,
                          _complex128 * yPtr);

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
