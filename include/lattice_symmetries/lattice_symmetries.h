#ifndef LATTICE_SYMMETRIES_H
#define LATTICE_SYMMETRIES_H

/// \file
/// \brief Public API
///
/// Long description!

#if defined(__cplusplus)
#    include <cstdint>
#else
#    include <stdint.h>
#endif

#if defined(__cplusplus)
extern "C" {
#endif

/// \brief Status codes used by the library
///
enum ls_error_code {
    LS_SUCCESS = 0,             ///< No error
    LS_OUT_OF_MEMORY,           ///< Memory allocation failed
    LS_INVALID_ARGUMENT,        ///< Argument to a function is invalid
    LS_PERMUTATION_TOO_LONG,    ///< Such long permutations are not supported
    LS_NOT_A_PERMUTATION,       ///< Argument is not a valid permutation of `{0, 1, ..., N-1}`
    LS_INVALID_SECTOR,          ///< Sector exceeds the periodicity of the operator
    LS_INCOMPATIBLE_SYMMETRIES, ///< Symmetries are incompatible
    LS_NOT_A_REPRESENTATIVE,    ///< Spin configuration is not a representative
    LS_INVALID_NUMBER_SPINS,
    LS_INVALID_HAMMING_WEIGHT,
    LS_WRONG_BASIS_TYPE,
    LS_CACHE_NOT_BUILT,
    LS_COULD_NOT_OPEN_FILE,
    LS_FILE_IO_FAILED,
    LS_CACHE_IS_CORRUPT,
    LS_INVALID_STATE,
    LS_OPERATOR_IS_COMPLEX,
    LS_SYSTEM_ERROR, ///< Unknown error
};

/// Given an error code, obtains the corresponding error message.
char const* ls_error_to_string(ls_error_code code);
void        ls_destroy_string(char const* message);

/// \defgroup symmetries Individual symmetries
/// @{

/// \brief Symmetry operator
///
/// All lattice symmetries are built from two primitive operations: permutation and spin inversion.
/// Furthermore, each symmetry has corresponding eigenvalue to which we restrict the Hilbert space.
typedef struct ls_symmetry ls_symmetry;

/// \brief Allocates and constructs a symmetry.
///
/// After successful completion of this function \p ptr points to the newly constructed
/// #ls_symmetry.
///
/// \note #ls_symmetry must be destructed and deallocated using #ls_destroy_symmetry.
///
/// \param ptr is set to point to the newly constructed object.
/// \param length length of the \p permutation array.
/// \param permutation permutation of `{0, 1, ..., length - 1}`.
/// \param flip whether application of this symmetry inverts spins.
/// \param sector symmetry sector to which we restrict the problem.
/// \return #LS_SUCCESS on successful completion. #LS_NOT_A_PERMUTATION if \p permutation does not
///         form a valid permutation. #LS_PERMUTATION_TOO_LONG if \p length exceeds 512.
///         #LS_INVALID_SECTOR if \p sector exceeds the periodicity of the symmetry.
/// \see #ls_destroy_symmetry
ls_error_code ls_create_symmetry(ls_symmetry** ptr, unsigned length, unsigned const permutation[],
                                 bool flip, unsigned sector);
/// \brief Destructs and deallocates a symmetry.
///
/// This function **must** be called on objects constructed using #ls_create_symmetry.
///
/// \param symmetry pointer to symmetry. Must not be `nullptr`.
void ls_destroy_symmetry(ls_symmetry* symmetry);
/// \brief Get symmetry sector.
///
/// \param symmetry pointer to symmetry. Must not be `nullptr`.
/// \return symmetry sector.
unsigned ls_get_sector(ls_symmetry const* symmetry);
/// \brief Get whether the symmetry applies spin inversion.
///
/// \param symmetry pointer to symmetry. Must not be `nullptr`.
/// \return `true` when \p symmetry applies spin inversion.
bool ls_get_flip(ls_symmetry const* symmetry);
/// \brief Get phase of the eigenvalue.
///
/// \param symmetry pointer to symmetry. Must not be `nullptr`.
/// \return complex phase of the eigenvalue of \p symmetry.
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
unsigned ls_get_periodicity(ls_symmetry const* symmetry);
/// Apply symmetry to a spin configuration
void ls_apply_symmetry(ls_symmetry const* symmetry, uint64_t bits[]);

unsigned ls_symmetry_get_number_spins(ls_symmetry const* symmetry);

/// @}
// end of symmetries group

/// \defgroup group Symmetry group
/// @{

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
unsigned ls_get_group_size(ls_group const* group);

/// @}
// end of group group

/// \defgroup basis Spin basis
/// @{

typedef struct ls_spin_basis ls_spin_basis;
typedef struct ls_states     ls_states;

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
ls_error_code  ls_create_spin_basis(ls_spin_basis** ptr, ls_group const* group,
                                    unsigned number_spins, int hamming_weight);
ls_spin_basis* ls_copy_spin_basis(ls_spin_basis const* basis);
void           ls_destroy_spin_basis(ls_spin_basis* basis);
unsigned       ls_get_number_spins(ls_spin_basis const* basis);
unsigned       ls_get_number_bits(ls_spin_basis const* basis);
int            ls_get_hamming_weight(ls_spin_basis const* basis);
bool           ls_has_symmetries(ls_spin_basis const* basis);
ls_error_code  ls_get_number_states(ls_spin_basis const* basis, uint64_t* out);
ls_error_code  ls_build(ls_spin_basis* basis);
void ls_get_state_info(ls_spin_basis* basis, uint64_t const bits[], uint64_t representative[],
                       void* character, double* norm);
ls_error_code ls_get_index(ls_spin_basis const* basis, uint64_t const bits[], uint64_t* index);

ls_error_code   ls_get_states(ls_states** ptr, ls_spin_basis const* basis);
void            ls_destroy_states(ls_states* states);
uint64_t const* ls_states_get_data(ls_states const* states);
uint64_t        ls_states_get_size(ls_states const* states);

ls_error_code ls_save_cache(ls_spin_basis const* basis, char const* filename);
ls_error_code ls_load_cache(ls_spin_basis* basis, char const* filename);

/// @}
// end of group basis

typedef struct ls_interaction ls_interaction;
typedef struct ls_operator    ls_operator;

ls_error_code ls_create_interaction1(ls_interaction** ptr, void const* matrix_2x2,
                                     unsigned number_nodes, uint16_t const* nodes);
ls_error_code ls_create_interaction2(ls_interaction** ptr, void const* matrix_4x4,
                                     unsigned number_edges, uint16_t const (*edges)[2]);
ls_error_code ls_create_interaction3(ls_interaction** ptr, void const* matrix_8x8,
                                     unsigned number_triangles, uint16_t const (*triangles)[3]);
ls_error_code ls_create_interaction4(ls_interaction** ptr, void const* matrix_16x16,
                                     unsigned number_plaquettes, uint16_t const (*plaquettes)[4]);
void          ls_destroy_interaction(ls_interaction* interaction);
bool          ls_interaction_is_real(ls_interaction const* interaction);

ls_error_code ls_create_operator(ls_operator** ptr, ls_spin_basis const* basis,
                                 unsigned number_terms, ls_interaction const terms[]);
void          ls_destroy_operator(ls_operator* op);
bool          ls_operator_is_real(ls_operator const* op);
ls_error_code ls_operator_matvec_f32(ls_operator const* op, float const* x, float* y);
ls_error_code ls_operator_matvec_f64(ls_operator const* op, double const* x, double* y);
ls_error_code ls_operator_matvec_c64(ls_operator const* op, void const* x, void* y);
ls_error_code ls_operator_matvec_c128(ls_operator const* op, void const* x, void* y);

#if defined(__cplusplus)
} // extern "C"
#endif

#endif // LATTICE_SYMMETRIES_H
