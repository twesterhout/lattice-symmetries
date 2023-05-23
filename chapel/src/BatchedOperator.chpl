module BatchedOperator {

use CTypes;
use Time;

use FFI;
use ForeignTypes;
use StatesEnumeration;

// TODO: this is currently implemented inefficiently
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

  var applyOffDiagTimer : Timer;
  var stateInfoTimer : Timer;
  var memcpyTimer : Timer;
  var coeffTimer : Timer;
  var keysTimer : Timer;

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
        matrix.payload,
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
        matrix.payload,
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
      matrix.payload,
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
    c_memcpy(tempSpins + totalCount, alphas, count:c_size_t * c_sizeof(uint(64)));
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
  ls_internal_operator_apply_diag_x1(matrix.payload, count, alphas, c_ptrTo(_cs[0]), nil);

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
      matrix.payload,
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
