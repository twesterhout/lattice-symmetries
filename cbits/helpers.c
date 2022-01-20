#include "helpers.h"
#include <assert.h>
#include <string.h>

static inline void ls_internal_copy_spin(int const number_words,
                                         uint64_t const *const restrict source,
                                         uint64_t *const restrict destination) {
  for (int j = 0; j < number_words; ++j) {
    destination[j] = source[j];
  }
}

static inline void ls_internal_fill_with_spin(
    int const number_words, uint64_t const *const restrict source,
    uint64_t const count, uint64_t *const restrict destination) {
  for (uint64_t i = 0; i < count; ++i) {
    ls_internal_copy_spin(number_words, source, destination + i * number_words);
  }
}

static inline unsigned
ls_internal_gather_bits(int const n, ls_bit_index const *const restrict tuple,
                        uint64_t const *const restrict source) {
  unsigned r = 0;
  for (int i = 0; i < n; ++i) {
    unsigned const b = (source[tuple[i].word] >> tuple[i].bit) & 1U;
    // NOTE: Order is important here!
    // I.e. doing (r | b) << 1U will break stuff!
    r = (r << 1U) | b;
  }
  return r;
}

static inline void
ls_internal_scatter_bits(int const n, ls_bit_index const *const restrict tuple,
                         unsigned r, uint64_t *const restrict destination) {
  for (int i = n; i-- > 0;) {
    uint64_t const w = destination[tuple[i].word];
    uint8_t const b = tuple[i].bit;
    destination[tuple[i].word] =
        (w & ~((uint64_t)1 << b)) | ((uint64_t)(r & 1U) << b);
    r >>= 1U;
  }
}

#if 0
static void ls_csr_matrix_from_dense(unsigned const dimension,
                                     _Complex double const *const dense,
                                     unsigned *offsets, unsigned *columns,
                                     _Complex double *off_diag_elements,
                                     _Complex double *diag_elements) {
  offsets[0] = 0;
  for (unsigned i = 0; i < dimension; ++i) {
    unsigned nonzero_in_row = 0;
    for (unsigned j = 0; j < dimension; ++j) {
      _Complex double const element = dense[i * dimension + j];
      if (i == j) {
        diag_elements[i] = element;
      } else if (element != 0) {
        *columns = j;
        *off_diag_elements = element;
        ++nonzero_in_row;
        ++columns;
        ++off_diag_elements;
      }
    }
    offsets[i + 1] = offsets[i] + nonzero_in_row;
  }
}

static void ls_dense_from_csr_matrix(
    unsigned const dimension, unsigned const *const offsets,
    unsigned const *const columns,
    _Complex double const *const off_diag_elements,
    _Complex double const *const diag_elements, _Complex double *const dense) {
  memset(dense, 0,
         (size_t)dimension * (size_t)dimension * sizeof(_Complex double));
  for (unsigned i = 0; i < dimension; ++i) {
    for (unsigned k = offsets[i]; k < offsets[i + 1]; ++k) {
      unsigned const j = columns[k];
      dense[i * dimension + j] = off_diag_elements[k];
    }
    dense[i * dimension + i] = diag_elements[i];
  }
}
#endif

typedef struct ls_internal_output_buffer {
  uint64_t *const spins;         // 2D array of shape [capacity, number_words]
                                 // (stored in row-major order)
  _Complex double *const coeffs; // 1D array of shape [capacity]
  _Complex double diagonal;
  uint64_t size;
  uint64_t const capacity;
  uint64_t const number_words;
} ls_internal_output_buffer;

static inline void
ls_internal_apply_term(ls_term const *const term, uint64_t const *const spin,
                       ls_internal_output_buffer *const out) {
  unsigned const k =
      ls_internal_gather_bits(term->tuple_size, term->tuple, spin);
  // Handle the diagonal
  out->diagonal += term->matrix->diag_elements[k];
  // Handle off-diagonal elements
  unsigned const b = term->matrix->offsets[k];
  unsigned const e = term->matrix->offsets[k + 1];
  unsigned const count = e - b;
  assert(out->size + count <= out->capacity);
  uint64_t *const spins_dest = out->spins + out->size * out->number_words;
  _Complex double *const coeff_dest = out->coeffs + out->size;
  ls_internal_fill_with_spin(out->number_words, spin, count, spins_dest);
  for (unsigned offset = 0; offset < count; ++offset) {
    unsigned const j = b + offset;
    ls_internal_scatter_bits(term->tuple_size, term->tuple,
                             term->matrix->columns[j],
                             spins_dest + offset * out->number_words);
    coeff_dest[offset] = term->matrix->off_diag_elements[j];
  }
  out->size += count;
}

static inline void ls_internal_apply_many_terms(
    uint64_t const number_terms, ls_term const *const terms,
    uint64_t const *const spin, ls_internal_output_buffer *const out) {
  for (uint64_t i = 0; i < number_terms; ++i) {
    ls_internal_apply_term(terms + i, spin, out);
  }
  // NOTE: Let's not forget the diagonal contribution
  if (out->diagonal != 0) {
    ls_internal_copy_spin(out->number_words, spin,
                          out->spins + out->size * out->number_words);
    out->coeffs[out->size] = out->diagonal;
    ++(out->size);
  }
}

void ls_internal_apply_operator(
    uint64_t const number_terms, ls_term const *const terms,
    uint64_t const number_spins, uint64_t const *const spins,
    uint64_t const capacity, uint64_t const number_words, uint64_t *out_spins,
    _Complex double *out_coeffs, uint64_t *out_counts) {
  uint64_t offset = 0;
  for (uint64_t batch_index = 0; batch_index < number_spins; ++batch_index) {
    assert(offset <= capacity);
    ls_internal_output_buffer out = {.spins = out_spins + offset * number_words,
                                     .coeffs = out_coeffs + offset,
                                     .diagonal = 0,
                                     .size = 0,
                                     .capacity = capacity - offset,
                                     .number_words = number_words};
    uint64_t const *const current_spin = spins + batch_index * number_words;
    ls_internal_apply_many_terms(number_terms, terms, current_spin, &out);
    out_counts[batch_index] = out.size;
    offset += out.size;
  }
}

// uint64_t ls_hs_apply_term(ls_term const *const term, uint64_t const *const
// spin,
//                           ls_output_buffer *const out) {
//   _Complex double *const original_out_coeffs = out->coeffs;
//   // ls_term_fill_fn const fill_fn = get_fill_fn_for(out);
//   for (unsigned i = 0; i < term->number_tuples; ++i) {
//     ls_bit_index const *const tuple = term->tuples + i * term->tuple_size;
//     unsigned const k = (*term->gather_fn)(spin, tuple);
//     // Handle off-diagonal elements
//     unsigned const b = term->matrix.offsets[k];
//     unsigned const e = term->matrix.offsets[k + 1];
//     (*out->spin_fill)(spin, e - b, out->spins);
//     for (unsigned j = b; j < e;
//          ++j, out->spins += out->number_words, ++out->coeffs) {
//       (*term->scatter_fn)(term->matrix.columns[j], tuple, out->spins);
//       *out->coeffs = term->matrix.off_diag_elements[j];
//     }
//     // Handle the diagonal
//     *out->diagonal += term->matrix.diag_elements[k];
//   }
//   uint64_t const written = out->coeffs - original_out_coeffs;
//   return written;
// }

// uint64_t ls_hs_apply_operator(ls_hs_operator const *const operator,
//                               uint64_t const * const spin,
//                               ls_output_buffer *const out,
//                               ls_workspace *const workspace) {
//   ls_output_buffer temp = {workspace->spins, out->coeffs, out->diagonal,
//                            out->number_words};
//   for (unsigned i = 0; i < operator->payload->number_terms; ++i) {
//     ls_hs_apply_term(operator->payload->terms + i, spin, &temp);
//   }
//   uint64_t const count = temp.coeffs - out->coeffs;
//   // Append original spin to other_spins to include it in the calculation of
//   // norms and characters
//   (*out->spin_fill)(spin, 1, temp.spins);
//   *temp.coeffs = (_Complex double)1.0;
//
//   ls_flat_spin_basis_state_info(operator->payload->basis, count + 1,
//                                 workspace->spins, out->spins,
//                                 workspace->characters, workspace->norms);
//   double const current_norm = workspace->norms[count];
//   LATTICE_SYMMETRIES_CHECK(current_norm > 0, "");
//
//   uint64_t i = 0;
//   for (; i < count && (workspace->norms[i] != 0); ++i) {
//     out->coeffs[i] = out->coeffs[i] * (workspace->norms[i] / current_norm) *
//                      workspace->characters[i];
//   }
//   uint64_t offset = i;
//   // NOTE: This will be triggered only if our operator does not preserve
//   // symmetries of the basis
//   if (i != count) {
//     for (++i; i < count; ++i) {
//       if (workspace->norms[i] != 0) {
//         (*out->spin_copy)(out->spins + i * out->number_words,
//                           out->spins + offset * out->number_words);
//         out->coeffs[offset] = out->coeffs[i] *
//                               (workspace->norms[i] / current_norm) *
//                               workspace->characters[i];
//         ++offset;
//       }
//     }
//   }
//   // Increment pointers to final values
//   out->spins = out->spins + offset * out->number_words;
//   out->coeffs = out->coeffs + offset;
//   return offset;
// }

#if 0
#define DEFINE_TERM_FILL_FN(k)                                                 \
  void ls_internal_spin_copy_##k(uint64_t const *restrict source,              \
                                 uint64_t *restrict destination) {             \
    for (int j = 0; j < (k); ++j) {                                            \
      destination[j] = source[j];                                              \
    }                                                                          \
  }                                                                            \
  void ls_internal_spin_fill_##k(uint64_t const *restrict source,              \
                                 uint64_t count,                               \
                                 uint64_t *restrict destination) {             \
    for (uint64_t i = 0; i < count; ++i, destination += (k)) {                 \
      ls_internal_spin_copy_##k(source, destination);                          \
    }                                                                          \
  }

#define DEFINE_TERM_GATHER_FN(n)                                               \
  unsigned ls_internal_term_gather_##n(uint64_t const *source,                 \
                                       ls_bit_index const *tuple) {            \
    unsigned r = 0;                                                            \
    for (int i = 0; i < (n); ++i) {                                            \
      unsigned const b = (source[tuple[i].word] >> tuple[i].bit) & 1U;         \
      /* NOTE: Order is important here! I.e. doing (r | b) << 1U will break    \
       * stuff! */                                                             \
      r = (r << 1U) | b;                                                       \
    }                                                                          \
    return r;                                                                  \
  }

#define DEFINE_TERM_SCATTER_FN(n)                                              \
  void ls_internal_term_scatter_##n(unsigned r, ls_bit_index const *tuple,     \
                                    uint64_t *destination) {                   \
    for (int i = (n); i-- > 0;) {                                              \
      uint64_t const w = destination[tuple[i].word];                           \
      uint8_t const b = tuple[i].bit;                                          \
      destination[tuple[i].word] =                                             \
          (w & ~((uint64_t)1 << b)) | ((uint64_t)(r & 1U) << b);               \
      r >>= 1U;                                                                \
    }                                                                          \
  }

DEFINE_TERM_FILL_FN(1)
DEFINE_TERM_FILL_FN(8)

DEFINE_TERM_GATHER_FN(1)
DEFINE_TERM_GATHER_FN(2)
DEFINE_TERM_GATHER_FN(3)
DEFINE_TERM_GATHER_FN(4)

DEFINE_TERM_SCATTER_FN(1)
DEFINE_TERM_SCATTER_FN(2)
DEFINE_TERM_SCATTER_FN(3)
DEFINE_TERM_SCATTER_FN(4)
#endif
