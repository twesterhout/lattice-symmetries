module BatchedOperator {

use BitOps;
use CTypes;
use Time;
import OS.POSIX;

use FFI;
use ForeignTypes;
use StatesEnumeration;

// TODO: this is currently implemented inefficiently
/*
private proc localCompressMultiply(batchSize : int, numberTerms : int,
                                   sigmas : c_ptr(uint(64)), cs : c_ptr(complex(128)),
                                   xs : c_ptr(?eltType), offsets : c_ptr(int)) {
  var offset : int = 0;
  var i : int = 0;

  offsets[0] = 0;
  for batchIdx in 0 ..# batchSize {
    // foreach termIdx in 0 ..# numberTerms {
    //   cs[batchIdx * numberTerms + termIdx] *= xs[batchIdx];
    // }
    for termIdx in 0 ..# numberTerms {
      if (cs[i] != 0) {
        if (i > offset) {
          cs[offset] = cs[i];
          sigmas[offset] = sigmas[i];
        }
        cs[offset] *= xs[batchIdx];
        offset += 1;
      }
      i += 1;
    }
    offsets[batchIdx + 1] = offset;
  }
  // offsets[batchSize] = batchSize * numberTerms;
}
*/

proc ls_internal_operator_apply_diag_x1(
    const ref op : ls_hs_operator,
    batch_size : int,
    alphas : c_ptrConst(uint(64)),
    ys : c_ptr(real(64)),
    xs : c_ptr(real(64))) {
  // The diagonal is zero
  if (op.diag_terms == nil || op.diag_terms.deref().number_terms == 0) {
    POSIX.memset(ys, 0, batch_size:c_size_t * c_sizeof(real(64)));
    return;
  }

  const ref terms = op.diag_terms.deref();
  const number_terms = terms.number_terms;

  foreach batch_idx in 0 ..# batch_size {
    var acc : real(64) = 0;
    const alpha = alphas[batch_idx];
    for term_idx in 0 ..# number_terms {
      const delta = (alpha & terms.m[term_idx]) == terms.r[term_idx];
      if delta {
        const sign = 1 - 2 * (popCount(alpha & terms.s[term_idx]) % 2);
        const factor = if xs != nil then sign * xs[batch_idx] else sign;
        acc += terms.v[term_idx].re * factor;
      }
    }
    ys[batch_idx] = acc;
  }
}
// void ls_internal_operator_apply_diag_x1(ls_hs_operator const *const op,
//                                         ptrdiff_t const batch_size,
//                                         uint64_t const *restrict const alphas,
//                                         double *restrict const ys,
//                                         double const *restrict const xs) {
//   // The diagonal is zero
//   if (op->diag_terms == NULL || op->diag_terms->number_terms == 0) {
//     memset(ys, 0, (size_t)batch_size * sizeof(double));
//     return;
//   }
// 
//   ls_hs_nonbranching_terms const *restrict terms = op->diag_terms;
//   ptrdiff_t const number_terms = terms->number_terms;
//   for (ptrdiff_t batch_idx = 0; batch_idx < batch_size; ++batch_idx) {
//     double acc = 0;
//     uint64_t const alpha = alphas[batch_idx];
//     for (ptrdiff_t term_idx = 0; term_idx < number_terms; ++term_idx) {
//       uint8_t const delta = (alpha & terms->m[term_idx]) == terms->r[term_idx];
//       if (delta != 0) {
//         int const sign =
//             1 - 2 * (__builtin_popcountll(alpha & terms->s[term_idx]) % 2);
//         double const factor =
//             (xs != NULL) ? delta * sign * xs[batch_idx] : delta * sign;
//         acc += creal(terms->v[term_idx]) * factor;
//       }
//     }
//     ys[batch_idx] = acc;
//   }
// }
// 
proc ls_internal_operator_apply_off_diag_x1(
    const ref op : ls_hs_operator,
    batch_size : int,
    alphas : c_ptrConst(uint(64)),
    betas : c_ptr(uint(64)),
    coeffs : c_ptr(complex(128)),
    offsets : c_ptr(c_ptrdiff),
    xs : c_ptrConst(real(64))) {
  // Nothing to apply
  if (op.off_diag_terms == nil || op.off_diag_terms.deref().number_terms == 0) {
    POSIX.memset(offsets, 0, (batch_size + 1):c_size_t * c_sizeof(c_ptrdiff));
    return;
  }

  const ref terms = op.off_diag_terms.deref();
  const number_terms = terms.number_terms;

  offsets[0] = 0;
  var offset = 0;
  for batch_idx in 0 ..# batch_size {
    for term_idx in 0 ..# number_terms {
      const alpha = alphas[batch_idx];
      const delta = (alpha & terms.m[term_idx]) == terms.r[term_idx];
      if delta {
        const sign = 1 - 2 * (popCount(alpha & terms.s[term_idx]) % 2);
        const factor = if xs != nil then sign * xs[batch_idx] else sign;
        coeffs[offset] = terms.v[term_idx] * factor;
        betas[offset] = alpha ^ terms.x[term_idx];
        offset += 1;
      }
    }
    offsets[batch_idx + 1] = offset;
  }
}
// void ls_internal_operator_apply_off_diag_x1(ls_hs_operator const *const op,
//                                             ptrdiff_t const batch_size,
//                                             uint64_t const *const alphas,
//                                             uint64_t *const betas,
//                                             ls_hs_scalar *const coeffs,
//                                             ptrdiff_t *const offsets,
//                                             double const *const xs) {
//   // Nothing to apply
//   if (op->off_diag_terms == NULL || op->off_diag_terms->number_terms == 0) {
//     memset(offsets, 0, (size_t)(batch_size + 1) * sizeof(ptrdiff_t));
//     return;
//   }
// 
//   ls_hs_nonbranching_terms const *restrict terms = op->off_diag_terms;
//   ptrdiff_t const number_terms = terms->number_terms;
// 
//   offsets[0] = 0;
//   ptrdiff_t offset = 0;
//   for (ptrdiff_t batch_idx = 0; batch_idx < batch_size; ++batch_idx) {
//     for (ptrdiff_t term_idx = 0; term_idx < number_terms; ++term_idx) {
//       uint8_t const delta =
//           (alphas[batch_idx] & terms->m[term_idx]) == terms->r[term_idx];
//       if (delta != 0) {
//         int const sign = 1 - 2 * (__builtin_popcountll(alphas[batch_idx] &
//                                                        terms->s[term_idx]) %
//                                   2);
//         double const factor =
//             (xs != NULL) ? delta * sign * xs[batch_idx] : delta * sign;
//         coeffs[offset] = CMPLX(creal(terms->v[term_idx]) * factor,
//                                cimag(terms->v[term_idx]) * factor);
//         betas[offset] = alphas[batch_idx] ^ terms->x[term_idx];
//         ++offset;
//       }
//     }
//     offsets[batch_idx + 1] = offset;
//   }
//   return;
// }


config const batchedOperatorNewKernel = true;

record BatchedOperator {
  var _matrixPtr : c_ptr(Operator);
  var batchSize : int;
  var _numberOffDiagTerms : int;
  var _dom : domain(1);
  var _spins1 : [_dom] uint(64);
  var _spins2 : [_dom] uint(64);
  var _coeffs1 : [_dom] complex(128);
  var _coeffs2 : [_dom] complex(128);
  var _norms : [_dom] real(64);
  var _localeIdxs : [_dom] uint(8);
  var _offsets : [0 ..# batchSize + 1] int;

  var applyOffDiagTimer : stopwatch;
  var stateInfoTimer : stopwatch;
  var memcpyTimer : stopwatch;
  var coeffTimer : stopwatch;
  var keysTimer : stopwatch;

  proc init(const ref matrix : Operator, batchSize : int) {
    this._matrixPtr = c_const_ptrTo(matrix);
    this.batchSize = batchSize;
    this._numberOffDiagTerms = matrix.numberOffDiagTerms();
    const numberTerms = max(_numberOffDiagTerms, 1);
    this._dom = {0 ..# (batchSize * (numberTerms + 1))};
    // logDebug("BatchedOperator._dom = ", this._dom);
  }
  proc init=(const ref other : BatchedOperator) {
    assert(other.locale == here);
    this._matrixPtr = other._matrixPtr;
    this.batchSize = other.batchSize;
    this._numberOffDiagTerms = other._numberOffDiagTerms;
    this._dom = other._dom;
  }

  inline proc matrix ref { return _matrixPtr.deref(); }
  inline proc basis ref { return matrix.basis; }
  // proc deinit() {
    // logDebug("BatchedOperator spent ", compressTime, " in localCompressMultiply and ",
    //          offDiagTime, " in apply_off_diag");
  // }

  proc computeOffDiag(count : int,
                      alphas : c_ptr(uint(64)),
                      xs : c_ptr(?eltType))
      : (int, c_ptr(uint(64)), c_ptr(complex(128)), c_ptr(uint(8))) {

    assert(count <= batchSize);
    // Simple case when no symmetries are used
    if !matrix.basis.requiresProjection() {
      const betas = c_pointer_return(_spins1[0]);
      const cs = c_pointer_return(_coeffs1[0]);
      const offsets = c_pointer_return(_offsets[0]);
      const keys = c_pointer_return(_localeIdxs[0]);
      // We compute H|αᵢ⟩ = ∑ⱼ cᵢⱼ|βᵢⱼ⟩ for each |αᵢ⟩ where i ∈ {0 ..# count}
      // at this point j ∈ {0 ..# _numberOffDiagTerms}. Both cᵢⱼ and |βᵢⱼ⟩ are
      // conceptually 2-dimensional arrays, but we use flattened
      // representations of them.
      applyOffDiagTimer.start();
      ls_internal_operator_apply_off_diag_x1(
        matrix.payload.deref(),
        count,
        alphas,
        betas,
        cs,
        offsets,
        xs);
      applyOffDiagTimer.stop();

      keysTimer.start();
      const totalCount = offsets[count];
      foreach i in 0 ..# totalCount {
        keys[i] = localeIdxOf(betas[i]):uint(8);
      }
      keysTimer.stop();
      return (totalCount, betas, cs, keys);
    }
    // Another simple case when no permutations were given,
    // only the spin inversion
    if !matrix.basis.hasPermutationSymmetries() && matrix.basis.hasSpinInversionSymmetry() {
      const betas = c_pointer_return(_spins1[0]);
      const cs = c_pointer_return(_coeffs1[0]);
      const offsets = c_pointer_return(_offsets[0]);
      const keys = c_pointer_return(_localeIdxs[0]);
      // We compute H|αᵢ⟩ = ∑ⱼ cᵢⱼ|βᵢⱼ⟩ for each |αᵢ⟩ where i ∈ {0 ..# count}
      // at this point j ∈ {0 ..# _numberOffDiagTerms}. Both cᵢⱼ and |βᵢⱼ⟩ are
      // conceptually 2-dimensional arrays, but we use flattened
      // representations of them.
      applyOffDiagTimer.start();
      ls_internal_operator_apply_off_diag_x1(
        matrix.payload.deref(),
        count,
        alphas,
        betas,
        cs,
        offsets,
        xs);
      applyOffDiagTimer.stop();

      const totalCount = offsets[count];
      const mask = (1:uint(64) << matrix.basis.numberSites()) - 1;
      const character = matrix.basis.spinInversion;
      assert(character != 0);

      stateInfoTimer.start();
      foreach i in 0 ..# totalCount {
        const current = betas[i];
        const inverted = current ^ mask;
        if inverted < current {
          betas[i] = inverted;
          cs[i] *= character;
        }
      }
      stateInfoTimer.stop();

      keysTimer.start();
      foreach i in 0 ..# totalCount {
        keys[i] = localeIdxOf(betas[i]):uint(8);
      }
      keysTimer.stop();
      return (totalCount, betas, cs, keys);
    }

    // The tricky case when we have to project betas first
    const tempSpins = c_pointer_return(_spins2[0]);
    const tempCoeffs = c_pointer_return(_coeffs2[0]);
    const offsets = c_pointer_return(_offsets[0]);
    applyOffDiagTimer.start();
    ls_internal_operator_apply_off_diag_x1(
      matrix.payload.deref(),
      count,
      alphas,
      tempSpins,
      tempCoeffs,
      offsets,
      xs);
    applyOffDiagTimer.stop();
    const totalCount = offsets[count];
    // We are also interested in norms of alphas, so we append them to tempSpins
    memcpyTimer.start();
    assert(totalCount + count <= _dom.size);
    POSIX.memcpy(tempSpins + totalCount, alphas, count:c_size_t * c_sizeof(uint(64)));
    memcpyTimer.stop();

    const betas = c_pointer_return(_spins1[0]);
    const cs = c_pointer_return(_coeffs1[0]);
    const norms = c_pointer_return(_norms[0]);
    stateInfoTimer.start();
    ls_hs_state_info(
      basis.payload,
      totalCount + count,
      tempSpins, 1,
      betas, 1,
      cs,
      norms);
    stateInfoTimer.stop();

    coeffTimer.start();
    foreach i in 0 ..# count {
      foreach k in offsets[i] ..< offsets[i + 1] {
        cs[k] *= tempCoeffs[k] * norms[k] / norms[totalCount + i];
      }
    }
    coeffTimer.stop();

    keysTimer.start();
    const keys = c_pointer_return(_localeIdxs[0]);
    foreach i in 0 ..# totalCount {
      keys[i] = localeIdxOf(betas[i]):uint(8);
    }
    keysTimer.stop();

    return (totalCount, betas, cs, keys);
  }

}

export proc ls_chpl_operator_apply_diag(matrixPtr : c_ptr(ls_hs_operator),
                                        count : int,
                                        alphas : c_ptr(uint(64)),
                                        coeffs : c_ptr(chpl_external_array),
                                        numTasks : int) {
  logDebug("Calling ls_chpl_operator_apply_diag ...");
  var matrix = new Operator(matrixPtr, owning=false);
  if matrix.basis.numberWords != 1 then
    halt("bases with more than 64 bits are not yet implemented");
  if matrix.basis.requiresProjection() then
    halt("bases that require projection are not yet supported");

  var _cs : [0 ..# count] real(64) = noinit;
  ls_internal_operator_apply_diag_x1(
    matrix.payload.deref(), count, alphas, c_ptrTo(_cs[0]), nil);

  coeffs.deref() = convertToExternalArray(_cs);
  logDebug("Done! Returning ...");
}

export proc ls_chpl_operator_apply_off_diag(matrixPtr : c_ptr(ls_hs_operator),
                                            count : int,
                                            alphas : c_ptr(uint(64)),
                                            betas : c_ptr(chpl_external_array),
                                            coeffs : c_ptr(chpl_external_array),
                                            offsets : c_ptr(chpl_external_array),
                                            numTasks : int) {
  logDebug("Calling ls_chpl_operator_apply_off_diag ...");
  var matrix = new Operator(matrixPtr, owning=false);
  if matrix.basis.numberWords != 1 then
    halt("bases with more than 64 bits are not yet implemented");
  if matrix.basis.requiresProjection() then
    halt("bases that require projection are not yet supported");

  const numberOffDiagTerms = matrix.numberOffDiagTerms();
  var _betas   : [0 ..# count * numberOffDiagTerms] uint(64);
  var _cs      : [0 ..# count * numberOffDiagTerms] complex(128);
  var _offsets : [0 ..# count + 1] int;

  if numberOffDiagTerms != 0 {
    ls_internal_operator_apply_off_diag_x1(
      matrix.payload.deref(),
      count,
      alphas,
      c_ptrTo(_betas[0]),
      c_ptrTo(_cs[0]),
      c_ptrTo(_offsets[0]),
      nil);
    const betasSize = _offsets[count];
    betas.deref() = convertToExternalArray(_betas);
    coeffs.deref() = convertToExternalArray(_cs);
    offsets.deref() = convertToExternalArray(_offsets);
  }
  else {
    betas.deref() = new chpl_external_array(nil, 0, nil);
    coeffs.deref() = new chpl_external_array(nil, 0, nil);
    offsets.deref() = convertToExternalArray(_offsets);
  }
  logDebug("Done! Returning ...");
}

} // end module BatchedOperator
