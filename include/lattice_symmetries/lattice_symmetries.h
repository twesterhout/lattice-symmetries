// Copyright (c) 2019-2020, Tom Westerhout
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef LATTICE_SYMMETRIES_H
#define LATTICE_SYMMETRIES_H

#define LATTICE_SYMMETRIES_UNREACHABLE __builtin_unreachable()
#define LATTICE_SYMMETRIES_LIKELY(x) __builtin_expect(x, 1)
#define LATTICE_SYMMETRIES_UNLIKELY(x) __builtin_expect(x, 0)
#define LATTICE_SYMMETRIES_EXPORT __attribute__((visibility("default")))
#define LATTICE_SYMMETRIES_NORETURN __attribute__((noreturn))

#if defined(__clang__)
#    define LATTICE_SYMMETRIES_CLANG() 1
#    define LATTICE_SYMMETRIES_GCC() 0
#    define LATTICE_SYMMETRIES_MSVC() 0
#elif defined(__GNUC__) || defined(__GNUG__)
#    define LATTICE_SYMMETRIES_CLANG() 0
#    define LATTICE_SYMMETRIES_GCC() 1
#    define LATTICE_SYMMETRIES_MSVC() 0
#elif defined(_MSC_VER)
#    define LATTICE_SYMMETRIES_CLANG() 0
#    define LATTICE_SYMMETRIES_GCC() 0
#    define LATTICE_SYMMETRIES_MSVC() 1
#else
#    define LATTICE_SYMMETRIES_CLANG() 0
#    define LATTICE_SYMMETRIES_GCC() 0
#    define LATTICE_SYMMETRIES_MSVC() 0
#endif

#if !defined(LATTICE_SYMMETRIES_FORCEINLINE)
#    if defined(_MSC_VER)
#        define LATTICE_SYMMETRIES_FORCEINLINE __forceinline
#    elif defined(__GNUC__) && __GNUC__ > 3 // Clang also defines __GNUC__ (as 4)
#        define LATTICE_SYMMETRIES_FORCEINLINE inline __attribute__((__always_inline__))
#    else
#        define LATTICE_SYMMETRIES_FORCEINLINE inline
#    endif
#endif

#if !defined(LATTICE_SYMMETRIES_NOINLINE)
#    if defined(_MSC_VER)
#        define LATTICE_SYMMETRIES_NOINLINE __declspec(noinline)
#    elif defined(__GNUC__) && __GNUC__ > 3 // Clang also defines __GNUC__ (as 4)
#        if defined(__CUDACC__)             // nvcc doesn't always parse __noinline__,
#            define LATTICE_SYMMETRIES_NOINLINE __attribute__((noinline))
#        else
#            define LATTICE_SYMMETRIES_NOINLINE __attribute__((__noinline__))
#        endif
#    else
#        define LATTICE_SYMMETRIES_NOINLINE
#    endif
#endif

#if defined(__AVX2__)
#    define LATTICE_SYMMETRIES_HAS_AVX2() 1
#else
#    define LATTICE_SYMMETRIES_HAS_AVX2() 0
#endif

#if defined(__AVX__)
#    define LATTICE_SYMMETRIES_HAS_AVX() 1
#else
#    define LATTICE_SYMMETRIES_HAS_AVX() 0
#endif

#if defined(__SSE4_1__) && defined(__SSE4_2__)
#    define LATTICE_SYMMETRIES_HAS_SSE4() 1
#else
#    define LATTICE_SYMMETRIES_HAS_SSE4() 0
#endif

#if !defined(__SSE2__) && !defined(__x86_64__)
#    error "unsupported architecture; lattice-symmetries currently only works on x86_64"
#endif

#if defined(__cplusplus)
#    include <cstdint>
#else
#    include <stdbool.h>
#    include <stdint.h>
#endif

#if !defined(LATTICE_SYMMETRIES_COMPLEX128)
#    if defined(__cplusplus)
#        include <complex>
#        define LATTICE_SYMMETRIES_COMPLEX128 std::complex<double>
#    else
#        define LATTICE_SYMMETRIES_COMPLEX128 _Complex double
#    endif
#endif

#if defined(__cplusplus)
extern "C" {
#endif

typedef enum ls_error_code {
    LS_SUCCESS = 0,             ///< No error
    LS_OUT_OF_MEMORY,           ///< Memory allocation failed
    LS_INVALID_ARGUMENT,        ///< Argument to a function is invalid
    LS_INVALID_HAMMING_WEIGHT,  ///< Invalid Hamming weight
    LS_INVALID_SPIN_INVERSION,  ///< Invalid value for spin_inversion
    LS_INVALID_NUMBER_SPINS,    ///< Invalid number of spins
    LS_INVALID_PERMUTATION,     ///< Argument is not a valid permutation
    LS_INVALID_SECTOR,          ///< Sector exceeds the periodicity of the operator
    LS_INVALID_STATE,           ///< Invalid basis state
    LS_INVALID_DATATYPE,        ///< Invalid datatype
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

char const* ls_error_to_string(ls_error_code code);
void        ls_destroy_string(char const* message);

typedef void (*ls_error_handler)(char const* expr, char const* file, unsigned line,
                                 char const* function, char const* msg);

void ls_set_check_fail_handler(ls_error_handler func);
void ls_set_assert_fail_handler(ls_error_handler func);

LATTICE_SYMMETRIES_NORETURN
void ls_assert_fail(char const* expr, char const* file, unsigned line, char const* function,
                    char const* msg);

LATTICE_SYMMETRIES_NORETURN
void ls_check_fail(char const* expr, char const* file, unsigned line, char const* function,
                   char const* msg);

#define LATTICE_SYMMETRIES_CHECK(cond, msg)                                                        \
    (LATTICE_SYMMETRIES_LIKELY(cond)                                                               \
         ? ((void)0)                                                                               \
         : ls_check_fail(#cond, __FILE__, __LINE__, __FUNCTION__, msg))

#if !defined(NDEBUG)
#    define LATTICE_SYMMETRIES_ASSERT(cond, msg)                                                   \
        (LATTICE_SYMMETRIES_LIKELY(cond)                                                           \
             ? ((void)0)                                                                           \
             : ls_assert_fail(#cond, __FILE__, __LINE__, __FUNCTION__, msg))
#else
#    define LATTICE_SYMMETRIES_ASSERT(cond, msg) ((void)0)
#endif

bool ls_is_logging_enabled();
void ls_enable_logging();
void ls_disable_logging();

bool ls_has_avx2();
bool ls_has_avx();
bool ls_has_sse4();

#define LATTICE_SYMMETRIES_DISPATCH(func, ...)                                                     \
    if (ls_has_avx2()) { return ::lattice_symmetries::avx2::func(__VA_ARGS__); }                   \
    if (ls_has_avx()) { return ::lattice_symmetries::avx::func(__VA_ARGS__); }                     \
    if (ls_has_sse4()) { return ::lattice_symmetries::sse4::func(__VA_ARGS__); }                   \
    return ::lattice_symmetries::sse2::func(__VA_ARGS__)

// This is an internal function!
void ls_private_log_debug(char const* file, unsigned line, char const* function, char const* fmt,
                          ...);

#define LATTICE_SYMMETRIES_LOG_DEBUG(fmt, ...)                                                     \
    ls_private_log_debug(__FILE__, __LINE__, __FUNCTION__, fmt, __VA_ARGS__)

typedef uint64_t ls_bits64;

typedef struct ls_bits512 {
    ls_bits64 words[8];
} ls_bits512;

typedef struct ls_symmetry ls_symmetry;

ls_error_code ls_create_symmetry(ls_symmetry** ptr, unsigned length, unsigned const permutation[],
                                 unsigned sector);
void          ls_destroy_symmetry(ls_symmetry* symmetry);
unsigned      ls_get_sector(ls_symmetry const* symmetry);
double        ls_get_phase(ls_symmetry const* symmetry);
void          ls_get_eigenvalue(ls_symmetry const* symmetry, void* out);
unsigned      ls_get_periodicity(ls_symmetry const* symmetry);
unsigned      ls_symmetry_get_number_spins(ls_symmetry const* symmetry);
unsigned      ls_symmetry_get_network_depth(ls_symmetry const* symmetry);
void ls_symmetry_get_network_masks(ls_symmetry const* symmetry, void* out, uint64_t stride);
void ls_symmetry_get_network_shifts(ls_symmetry const* symmetry, unsigned* shifts);
void ls_apply_symmetry(ls_symmetry const* symmetry, ls_bits512* bits);
void ls_batched_apply_symmetry(ls_symmetry const* symmetry, uint64_t count, uint64_t* spins,
                               uint64_t stride);

typedef struct ls_group ls_group;

ls_error_code      ls_create_group(ls_group** ptr, unsigned size, ls_symmetry const* generators[]);
ls_error_code      ls_create_trivial_group(ls_group** ptr, unsigned number_spins);
void               ls_destroy_group(ls_group* group);
unsigned           ls_get_group_size(ls_group const* group);
ls_symmetry const* ls_group_get_symmetries(ls_group const* group);
int                ls_group_get_number_spins(ls_group const* group);
int                ls_group_get_network_depth(ls_group const* group);
int                ls_group_dump_symmetry_info(ls_group const* group, void* masks, unsigned* shifts,
                                               LATTICE_SYMMETRIES_COMPLEX128* eigenvalues);

typedef struct ls_spin_basis ls_spin_basis;

typedef struct ls_states ls_states;

ls_error_code  ls_create_spin_basis(ls_spin_basis** ptr, ls_group const* group,
                                    unsigned number_spins, int hamming_weight, int spin_inversion);
ls_spin_basis* ls_copy_spin_basis(ls_spin_basis const* basis);
void           ls_destroy_spin_basis(ls_spin_basis* basis);
unsigned       ls_get_number_spins(ls_spin_basis const* basis);
unsigned       ls_get_number_bits(ls_spin_basis const* basis);
int            ls_get_hamming_weight(ls_spin_basis const* basis);
int            ls_get_spin_inversion(ls_spin_basis const* basis);
bool           ls_has_symmetries(ls_spin_basis const* basis);

ls_error_code ls_get_number_states(ls_spin_basis const* basis, uint64_t* out);
ls_error_code ls_build(ls_spin_basis* basis);
ls_error_code ls_build_unsafe(ls_spin_basis* basis, uint64_t size,
                              uint64_t const representatives[]);
void          ls_get_state_info(ls_spin_basis const* basis, ls_bits512 const* bits,
                                ls_bits512* representative, void* character, double* norm);
void ls_batched_get_state_info(ls_spin_basis const* basis, uint64_t count, ls_bits512 const* spins,
                               uint64_t spins_stride, ls_bits512* repr, uint64_t repr_stride,
                               LATTICE_SYMMETRIES_COMPLEX128* eigenvalues,
                               uint64_t eigenvalues_stride, double* norm, uint64_t norm_stride);
ls_error_code ls_get_index(ls_spin_basis const* basis, uint64_t bits, uint64_t* index);
ls_error_code ls_batched_get_index(ls_spin_basis const* basis, uint64_t count,
                                   ls_bits64 const* spins, uint64_t spins_stride, uint64_t* out,
                                   uint64_t out_stride);

ls_error_code   ls_get_states(ls_states** ptr, ls_spin_basis const* basis);
void            ls_destroy_states(ls_states* states);
uint64_t const* ls_states_get_data(ls_states const* states);
uint64_t        ls_states_get_size(ls_states const* states);

ls_error_code ls_save_cache(ls_spin_basis const* basis, char const* filename);
ls_error_code ls_load_cache(ls_spin_basis* basis, char const* filename);

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
                                 unsigned number_terms, ls_interaction const* const terms[]);
void          ls_destroy_operator(ls_operator* op);

typedef enum {
    LS_FLOAT32,
    LS_FLOAT64,
    LS_COMPLEX64,
    LS_COMPLEX128,
} ls_datatype;

typedef ls_error_code (*ls_callback)(ls_bits512 const* bits, void const* coeff, void* cxt);

ls_error_code ls_operator_apply(ls_operator const* op, ls_bits512 const* bits, ls_callback func,
                                void* cxt);

uint64_t ls_batched_operator_apply(ls_operator const* op, uint64_t count, ls_bits512 const* spins,
                                   ls_bits512* out_spins, LATTICE_SYMMETRIES_COMPLEX128* out_coeffs,
                                   uint64_t* out_counts);

ls_error_code ls_operator_matmat(ls_operator const* op, ls_datatype dtype, uint64_t size,
                                 uint64_t block_size, void const* x, uint64_t x_stride, void* y,
                                 uint64_t y_stride);

ls_error_code ls_operator_expectation(ls_operator const* op, ls_datatype dtype, uint64_t size,
                                      uint64_t block_size, void const* x, uint64_t x_stride,
                                      void* out);

uint64_t ls_operator_max_buffer_size(ls_operator const* op);

bool ls_operator_is_real(ls_operator const* op);

#if defined(__cplusplus)
} // extern "C"
#endif

#if defined(__cplusplus)
#    include <system_error>

namespace std {
template <> struct is_error_code_enum<ls_error_code> : true_type {};
} // namespace std

namespace lattice_symmetries {
class ls_error_category : public std::error_category {
  public:
    [[nodiscard]] auto name() const noexcept -> const char* final;
    [[nodiscard]] auto message(int c) const -> std::string final;
};
auto get_error_category() noexcept -> ls_error_category const&;
} // namespace lattice_symmetries

inline auto make_error_code(ls_error_code const e) noexcept -> std::error_code
{
    return {static_cast<int>(e), lattice_symmetries::get_error_category()};
}
#endif

#endif // LATTICE_SYMMETRIES_H
