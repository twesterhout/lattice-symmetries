#ifndef LATTICE_SYMMETRIES_H
#define LATTICE_SYMMETRIES_H

/// \file
/// \brief Public API
///
/// Long description!

/// \mainpage Lattice Symmetries
///
/// \section intro_sec Introduction
///
/// This is the introduction.
///
/// \section install_sec Installation
///
/// \subsection from_source_sec Compiling from source
///
/// The suggested way of installing `lattice_symmetries` library is compiling it from source. There
/// are almost no external dependencies, so the process is quite simple. To compile the code, you
/// will need the following:
///
///   * C & C++ compiler (with C++17 support);
///   * CMake (3.15+);
///   * Git
///
/// Start by cloning the repository:
///
/// ```{.sh}
/// git clone https://github.com/twesterhout/lattice-symmetries.git
/// cd lattice-symmetries
/// ```
///
/// Next, create a directory where build artifacts will be stored (we do not support in-source
/// builds):
///
/// ```{.sh}
/// mkdir build
/// cd build
/// ```
///
/// Run the configure step which will determine the compilers to use, download dependencies etc.
///
/// ```{.sh}
/// cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=</where/to/install> ..
/// ```
///
/// Build the library:
///
/// ```{.sh}
/// cmake --build .
/// ```
///
/// And finally install it:
///
/// ```{.sh}
/// cmake --build . --target install
/// ```
///

#if defined(__cplusplus)
#    include <cstdint>
#else
#    include <stdbool.h>
#    include <stdint.h>
#endif

#if defined(__cplusplus)
extern "C" {
#endif

/// \defgroup errors Error handling
/// \brief Status codes and corresponding utility functions
///
/// LatticeSymmetries library uses status codes for reporting errors. #ls_error_code specifies all
/// possible status codes which can be returned by library functions. Names of the constants should
/// explain errors pretty well, however for higher-leven wrappers it is useful to convert these
/// status codes to human-readable messages. #ls_error_to_string function provides this
/// functionality.
///
/// \note Even though internally the library is written in C++17, exceptions are disabled during
///       compilation. I.e. we never throw exceptions, all errors are reported using status codes.
///
/// @{

/// \brief Status codes used by the library
///
typedef enum ls_error_code {
    LS_SUCCESS = 0,             ///< No error
    LS_OUT_OF_MEMORY,           ///< Memory allocation failed
    LS_INVALID_ARGUMENT,        ///< Argument to a function is invalid
    LS_INVALID_HAMMING_WEIGHT,  ///< Invalid Hamming weight
    LS_INVALID_NUMBER_SPINS,    ///< Invalid number of spins
    LS_INVALID_PERMUTATION,     ///< Argument is not a valid permutation
    LS_INVALID_SECTOR,          ///< Sector exceeds the periodicity of the operator
    LS_INVALID_STATE,           ///< Invalid basis state
    LS_PERMUTATION_TOO_LONG,    ///< Such long permutations are not supported
    LS_INCOMPATIBLE_SYMMETRIES, ///< Symmetries are incompatible
    LS_NOT_A_REPRESENTATIVE,    ///< Spin configuration is not a representative
    LS_WRONG_BASIS_TYPE,        ///< Expected a basis of different type
    LS_CACHE_NOT_BUILT,         ///< List of representatives is not yet built
    LS_COULD_NOT_OPEN_FILE,     ///< Failed to open file
    LS_FILE_IO_FAILED,          ///< File input/output failed
    LS_CACHE_IS_CORRUPT,        ///< File does not contain a list of representatives
    LS_OPERATOR_IS_COMPLEX,     ///< Trying to apply complex operator to real vector
    LS_DIMENSION_MISMATCH,      ///< Operator dimension does not match vector length
    LS_SYSTEM_ERROR,            ///< Unknown error
} ls_error_code;

/// \brief Given an error code, obtains the corresponding error message.
///
/// \warning C-string is dynamically allocated and must be freed using #ls_destroy_string.
///
/// \param code status code returned by one of the library functions.
/// \return human readable description of the error.
// NOLINTNEXTLINE(modernize-use-trailing-return-type)
char const* ls_error_to_string(ls_error_code code);
/// \brief Deallocates error message.
///
/// \param message string obtained from #ls_error_to_string.
void ls_destroy_string(char const* message);

/// @}
// end of errors group

/// \defgroup symmetries Symmetries
/// \brief Working with individual symmetries
///
/// @{

/// \brief Symmetry operator
///
/// All lattice symmetries are built from two primitive operations: permutation and spin inversion.
/// Furthermore, each symmetry has corresponding eigenvalue to which we restrict the Hilbert space.
// NOLINTNEXTLINE(modernize-use-using)
typedef struct ls_symmetry ls_symmetry;

/// \brief Allocates and constructs a symmetry.
///
/// \warning #ls_symmetry must be destructed and deallocated using #ls_destroy_symmetry.
///
/// \param ptr upon successful completion \p ptr set to point to the newly constructed object. Must
///            not be `nullptr`.
/// \param length length of the \p permutation array. Must not exceed 512.
/// \param permutation permutation of `{0, 1, ..., length - 1}`.
/// \param flip whether application of this symmetry inverts spins.
/// \param sector symmetry sector to which we restrict the problem. It must be between zero
//                (inclusive) and periodicity of the symmetry operator (exclusive).
/// \return #LS_SUCCESS on successful completion. #LS_NOT_A_PERMUTATION if \p permutation does not
///         form a valid permutation. #LS_PERMUTATION_TOO_LONG if \p length exceeds 512.
///         #LS_INVALID_SECTOR if \p sector exceeds the periodicity of the symmetry.
/// \see #ls_destroy_symmetry
// NOLINTNEXTLINE(modernize-use-trailing-return-type)
ls_error_code ls_create_symmetry(ls_symmetry** ptr, unsigned length, unsigned const permutation[],
                                 bool flip, unsigned sector);
/// \brief Destructs and deallocates the symmetry.
///
/// This function **must** be called on objects constructed using #ls_create_symmetry.
///
/// \param symmetry pointer to symmetry. Must not be `nullptr`.
/// \see #ls_create_symmetry
void ls_destroy_symmetry(ls_symmetry* symmetry);
/// \brief Get symmetry sector.
///
/// \param symmetry pointer to symmetry. Must not be `nullptr`.
/// \return symmetry sector.
/// \see #ls_get_periodicity, #ls_get_phase, #ls_get_eigenvalue
// NOLINTNEXTLINE(modernize-use-trailing-return-type)
unsigned ls_get_sector(ls_symmetry const* symmetry);
/// \brief Determine whether the symmetry applies spin inversion.
///
/// \param symmetry pointer to symmetry. Must not be `nullptr`.
/// \return `true` if \p symmetry applies spin inversion.
// NOLINTNEXTLINE(modernize-use-trailing-return-type)
bool ls_get_flip(ls_symmetry const* symmetry);
/// \brief Get phase of the eigenvalue.
///
/// Phase is simply `sector / periodicity`. This function is provided for convenience only.
///
/// \param symmetry pointer to symmetry. Must not be `nullptr`.
/// \return complex phase of the eigenvalue of \p symmetry.
/// \see #ls_get_sector, #ls_get_periodicity, #ls_get_eigenvalue
// NOLINTNEXTLINE(modernize-use-trailing-return-type)
double ls_get_phase(ls_symmetry const* symmetry);
/// \brief Get eigenvalue of the symmetry.
///
/// Eigenvalue is written to \p out.
///
/// \param symmetry pointer to symmetry. Must not be `nullptr`.
/// \param out pointer to `std::complex<double>` or any other structure binary compatible
///            with it (e.g. `double _Complex` or `double[2]`). Must not be `nullptr`.
void ls_get_eigenvalue(ls_symmetry const* symmetry, void* out);
/// \brief Get periodicity of the symmetry.
///
/// \param symmetry pointer to symmetry. Must not be `nullptr`.
/// \return periodicity of the symmetry.
// NOLINTNEXTLINE(modernize-use-trailing-return-type)
unsigned ls_get_periodicity(ls_symmetry const* symmetry);
/// Apply symmetry to a spin configuration
void ls_apply_symmetry(ls_symmetry const* symmetry, uint64_t bits[]);

///
// NOLINTNEXTLINE(modernize-use-trailing-return-type)
unsigned ls_symmetry_get_number_spins(ls_symmetry const* symmetry);

/// @}
// end of symmetries group

/// \defgroup group Symmetry group
/// @{

// NOLINTNEXTLINE(modernize-use-using)
typedef struct ls_group ls_group;

/// \brief Allocates and symmetry group.
///
/// After successful completion of this function \p ptr points to the newly constructed
/// #ls_group.
///
/// \note #ls_group must be destructed and deallocated using #ls_destroy_group.
///
/// \param ptr is set to point to the newly constructed object.
/// \param size length of the \p generators array.
/// \param generators array of pointers to individual symmetries which act as generators.
/// \return #LS_SUCCESS on successful completion. #LS_INCOMPATIBLE_SYMMETRIES if some symmetries are
///         incompatible with each other.
/// \see #ls_destroy_group
// NOLINTNEXTLINE(modernize-use-trailing-return-type)
ls_error_code ls_create_group(ls_group** ptr, unsigned size, ls_symmetry const* generators[]);
/// \brief Destructs and deallocates the symmetry group.
///
/// This function **must** be called on objects constructed using #ls_create_group.
///
/// \param group pointer to the group. Must not be `nullptr`.
void ls_destroy_group(ls_group* group);
/// \brief Get size of the symmetry sector.
///
/// \param group pointer to symmetry group. Must not be `nullptr`.
/// \return number of elements in the group.
// NOLINTNEXTLINE(modernize-use-trailing-return-type)
unsigned ls_get_group_size(ls_group const* group);

/// @}
// end of group group

/// \defgroup basis Spin basis
/// @{

// NOLINTNEXTLINE(modernize-use-using)
typedef struct ls_spin_basis ls_spin_basis;
// NOLINTNEXTLINE(modernize-use-using)
typedef struct ls_states ls_states;

/// \brief Allocates and constructs a basis.
///
/// After successful completion of this function \p ptr points to the newly constructed
/// #ls_spin_basis.
///
/// \note #ls_spin_basis must be destructed and deallocated using #ls_destroy_spin_basis.
///
/// \param ptr is set to point to the newly constructed object.
/// \param group symmetry group. Group may be empty, but must not be `nullptr`.
/// \param number_spins number of spins in the system. When \p group is not empty, number of spins
///                     may be deduced from it. In such a case it must match the value of \p
///                     number_spins.
/// \param hamming_weight allows to restrict to a sector with particular magnetisation. Hamming
///                       weight is the number of spins pointing upward. If one does not wish to
///                       restrict to a sector with given magnetisation, \p hamming_weight should
///                       be set to `-1`.
// NOLINTNEXTLINE(modernize-use-trailing-return-type)
ls_error_code ls_create_spin_basis(ls_spin_basis** ptr, ls_group const* group,
                                   unsigned number_spins, int hamming_weight);
// NOLINTNEXTLINE(modernize-use-trailing-return-type)
ls_spin_basis* ls_copy_spin_basis(ls_spin_basis const* basis);
void           ls_destroy_spin_basis(ls_spin_basis* basis);
// NOLINTNEXTLINE(modernize-use-trailing-return-type)
unsigned ls_get_number_spins(ls_spin_basis const* basis);
// NOLINTNEXTLINE(modernize-use-trailing-return-type)
unsigned ls_get_number_bits(ls_spin_basis const* basis);
// NOLINTNEXTLINE(modernize-use-trailing-return-type)
int ls_get_hamming_weight(ls_spin_basis const* basis);
// NOLINTNEXTLINE(modernize-use-trailing-return-type)
bool ls_has_symmetries(ls_spin_basis const* basis);
// NOLINTNEXTLINE(modernize-use-trailing-return-type)
ls_error_code ls_get_number_states(ls_spin_basis const* basis, uint64_t* out);
// NOLINTNEXTLINE(modernize-use-trailing-return-type)
ls_error_code ls_build(ls_spin_basis* basis);
void ls_get_state_info(ls_spin_basis* basis, uint64_t const bits[], uint64_t representative[],
                       void* character, double* norm);
// NOLINTNEXTLINE(modernize-use-trailing-return-type)
ls_error_code ls_get_index(ls_spin_basis const* basis, uint64_t const bits[], uint64_t* index);

// NOLINTNEXTLINE(modernize-use-trailing-return-type)
ls_error_code ls_get_states(ls_states** ptr, ls_spin_basis const* basis);
void          ls_destroy_states(ls_states* states);
// NOLINTNEXTLINE(modernize-use-trailing-return-type)
uint64_t const* ls_states_get_data(ls_states const* states);
// NOLINTNEXTLINE(modernize-use-trailing-return-type)
uint64_t ls_states_get_size(ls_states const* states);

// NOLINTNEXTLINE(modernize-use-trailing-return-type)
ls_error_code ls_save_cache(ls_spin_basis const* basis, char const* filename);
// NOLINTNEXTLINE(modernize-use-trailing-return-type)
ls_error_code ls_load_cache(ls_spin_basis* basis, char const* filename);

/// @}
// end of group basis

// NOLINTNEXTLINE(modernize-use-using)
typedef struct ls_interaction ls_interaction;
// NOLINTNEXTLINE(modernize-use-using)
typedef struct ls_operator ls_operator;

// NOLINTNEXTLINE(modernize-use-trailing-return-type)
ls_error_code ls_create_interaction1(ls_interaction** ptr, void const* matrix_2x2,
                                     unsigned number_nodes, uint16_t const* nodes);
// NOLINTNEXTLINE(modernize-use-trailing-return-type)
ls_error_code ls_create_interaction2(ls_interaction** ptr, void const* matrix_4x4,
                                     unsigned number_edges, uint16_t const (*edges)[2]);
// NOLINTNEXTLINE(modernize-use-trailing-return-type)
ls_error_code ls_create_interaction3(ls_interaction** ptr, void const* matrix_8x8,
                                     unsigned number_triangles, uint16_t const (*triangles)[3]);
// NOLINTNEXTLINE(modernize-use-trailing-return-type)
ls_error_code ls_create_interaction4(ls_interaction** ptr, void const* matrix_16x16,
                                     unsigned number_plaquettes, uint16_t const (*plaquettes)[4]);
void          ls_destroy_interaction(ls_interaction* interaction);
// NOLINTNEXTLINE(modernize-use-trailing-return-type)
bool ls_interaction_is_real(ls_interaction const* interaction);

// NOLINTNEXTLINE(modernize-use-trailing-return-type)
ls_error_code ls_create_operator(ls_operator** ptr, ls_spin_basis const* basis,
                                 unsigned number_terms, ls_interaction const* const terms[]);
void          ls_destroy_operator(ls_operator* op);
// NOLINTNEXTLINE(modernize-use-trailing-return-type)
ls_error_code ls_operator_matvec_f32(ls_operator const* op, uint64_t size, float const* x,
                                     float* y);
// NOLINTNEXTLINE(modernize-use-trailing-return-type)
ls_error_code ls_operator_matvec_f64(ls_operator const* op, uint64_t size, double const* x,
                                     double* y);

ls_error_code ls_operator_matmat_f64(ls_operator const* op, uint64_t size, uint64_t block_size,
                                     double const* x, uint64_t x_stride, double* y,
                                     uint64_t y_stride);

// NOLINTNEXTLINE(modernize-use-trailing-return-type)
ls_error_code ls_operator_matvec_c64(ls_operator const* op, uint64_t size, void const* x, void* y);
// NOLINTNEXTLINE(modernize-use-trailing-return-type)
ls_error_code ls_operator_matvec_c128(ls_operator const* op, uint64_t size, void const* x, void* y);

// NOLINTNEXTLINE(modernize-use-trailing-return-type)
ls_error_code ls_operator_apply_64(ls_operator const* op, uint64_t bits, uint64_t* out_size,
                                   void* out_coeffs, uint64_t* out_bits);
// NOLINTNEXTLINE(modernize-use-trailing-return-type)
ls_error_code ls_operator_apply_512(ls_operator const* op, uint64_t const bits[], uint64_t out_size,
                                    uint64_t (*out)[8]);

bool ls_operator_is_real(ls_operator const* op);

#if defined(__cplusplus)
} // extern "C"
#endif

#endif // LATTICE_SYMMETRIES_H
