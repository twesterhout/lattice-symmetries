#include "lattice_symmetries_haskell.h"
#include <HsFFI.h>
#include <Rts.h>
#include <stdio.h>

void ls_hs_init(void) {
  int argc = 1;
  char *argv[] = {"lattice_symmetries", NULL};
  // "+RTS", "-N1", "--install-signal-handlers=no", "-RTS", NULL};
  char **pargv = argv;

  // For some reason, options from argv are not processed properly, so we
  // manually set all RTS options using rts_opts field of RtsConfig
  RtsConfig conf = defaultRtsConfig;
  conf.rts_opts_enabled = RtsOptsAll;
  conf.rts_opts = "-N1 --install-signal-handlers=no";
  hs_init_ghc(&argc, &pargv, conf);

  ls_hs_internal_set_free_stable_ptr(&hs_free_stable_ptr);
}

void ls_hs_exit(void) {
  // LATTICE_SYMMETRIES_LOG_DEBUG("%s", "Calling hs_exit...\n");
  hs_exit();
  // LATTICE_SYMMETRIES_LOG_DEBUG("%s", "Deinitialized RTS!\n");
}

// typedef struct ls_hs_operator_apply_callback_context {
//   uint64_t count;
//   uint64_t *spins;
//   _Complex double *coeffs;
// } ls_hs_operator_apply_callback_context;

// static ls_error_code ls_hs_operator_apply_callback(ls_bits512 const *_bits,
//                                                    void const *_coeff,
//                                                    void *_cxt) {
//   ls_hs_operator_apply_callback_context *cxt = _cxt;
//   *(cxt->spins + cxt->count) = _bits->words[0];
//   *(cxt->coeffs + cxt->count) = *(_Complex double const *)_coeff;
//   ++cxt->count;
//   return LS_SUCCESS;
// }

// ls_error_code ls_hs_operator_apply(ls_hs_operator_v1 const *op,
//                                    uint64_t const count,
//                                    uint64_t const *restrict spins,
//                                    uint64_t *restrict offsets,
//                                    uint64_t *restrict out_spins,
//                                    _Complex double *restrict out_coeffs) {
//   ls_hs_operator_apply_callback_context cxt = {0, out_spins, out_coeffs};
//   for (uint64_t i = 0; i < count; ++i) {
//     offsets[i] = cxt.count;
//     ls_error_code const status =
//         ls_operator_apply(op->payload, (ls_bits512 const *)(spins + i),
//                           &ls_hs_operator_apply_callback, &cxt);
//     if (LATTICE_SYMMETRIES_UNLIKELY(status != LS_SUCCESS)) {
//       return status;
//     }
//   }
//   offsets[count] = cxt.count;
//   return LS_SUCCESS;
// }

// ls_error_code ls_hs_basis_index(ls_hs_spin_basis_v1 const *basis,
//                                 uint64_t const count,
//                                 uint64_t const *restrict spins,
//                                 uint64_t *restrict indices) {
//   for (uint64_t i = 0; i < count; ++i) {
//     ls_error_code const status =
//         ls_get_index(basis->payload, spins[i], indices + i);
//     if (LATTICE_SYMMETRIES_UNLIKELY(status != LS_SUCCESS)) {
//       return status;
//     }
//   }
//   return LS_SUCCESS;
// }
