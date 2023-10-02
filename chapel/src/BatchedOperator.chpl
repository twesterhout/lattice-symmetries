module BatchedOperator {

use BitOps;
use CTypes;
use Time;
import OS.POSIX;
import Map;

use FFI;
use Utils;
use ForeignTypes;
use StatesEnumeration;
import CSR;

proc ls_internal_operator_apply_diag_x1(
    const ref op : ls_hs_operator,
    batch_size : int,
    alphas : c_ptrConst(uint(64)),
    ys : c_ptr(?eltType),
    xs : c_ptrConst(eltType)) {
  // ls_internal_operator_apply_diag_x1(c_ptrToConst(op), batch_size, alphas, ys, xs);
  // return;

  // The diagonal is zero
  if (op.diag_terms == nil || op.diag_terms.deref().number_terms == 0) {
    POSIX.memset(ys, 0, batch_size:c_size_t * c_sizeof(eltType));
    return;
  }

  const ref terms = op.diag_terms.deref();
  const number_terms = terms.number_terms;

  foreach batch_idx in 0 ..# batch_size {
    var acc : eltType = 0;
    const alpha = alphas[batch_idx];
    for term_idx in 0 ..# number_terms {
      const delta = (alpha & terms.m[term_idx]) == terms.r[term_idx];
      if delta {
        const sign = 1 - 2 * parity(alpha & terms.s[term_idx]):int;
        const factor = if xs != nil then sign * xs[batch_idx] else sign;
        acc += terms.v[term_idx].re * factor;
      }
    }
    ys[batch_idx] = acc;
  }
}

proc ls_internal_operator_apply_off_diag_x1(
    const ref op : ls_hs_operator,
    batch_size : int,
    alphas : c_ptrConst(uint(64)),
    betas : c_ptr(uint(64)),
    coeffs : c_ptr(complex(128)),
    offsets : c_ptr(c_ptrdiff),
    xs) {
  // ls_internal_operator_apply_off_diag_x1(c_ptrToConst(op), batch_size, alphas, betas, coeffs, offsets, xs);
  // return;

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
    const alpha = alphas[batch_idx];
    // we know that we can never get alpha again since all terms are off
    // diagonal, so using alpha as a special state is fine.
    var old_beta = alpha;
    var old_coeff : complex(128) = 0;
    for term_idx in 0 ..# number_terms {
      const delta = (alpha & terms.m[term_idx]) == terms.r[term_idx];
      if delta {
        const sign = 1 - 2 * parity(alpha & terms.s[term_idx]):int;
        const factor = if xs != nil then sign * xs[batch_idx] else sign;
        // assert(terms.v[term_idx] != 0);
        const coeff = terms.v[term_idx] * factor;
        const beta = alpha ^ terms.x[term_idx];
        if beta != old_beta {
          if old_beta != alpha {
            coeffs[offset] = old_coeff;
            betas[offset] = old_beta;
            offset += 1;
          }
          old_beta = beta;
          old_coeff = coeff;
        }
        else {
          old_coeff += coeff;
        }
        // if offset > minimal_offset && betas[offset - 1] == beta {
        //   coeffs[offset - 1] += coeff;
        // }
        // else {
        //   coeffs[offset] = coeff;
        //   betas[offset] = beta;
        //   offset += 1;
        // }
      }
    }
    if old_beta != alpha {
      coeffs[offset] = old_coeff;
      betas[offset] = old_beta;
      offset += 1;
    }
    offsets[batch_idx + 1] = offset;
  }
}

private proc computeKeys(totalCount : int,
                         betas : c_ptrConst(uint(64)),
                         keys : c_ptr(uint(8))) {
  if keys == nil then return;

  if numLocales == 1 then
    POSIX.memset(keys, 0, totalCount:c_size_t * c_sizeof(uint(8)));
  else
    foreach i in 0 ..# totalCount do
      keys[i] = localeIdxOf(betas[i]):uint(8);
}
private proc computeKeys(totalCount, betas, keys, ref keysTimer) {
  keysTimer.start();
  computeKeys(totalCount, betas, keys);
  keysTimer.stop();
}

private proc computeOffDiagNoProjection(const ref matrix : Operator,
                                        count : int,
                                        alphas : c_ptrConst(uint(64)),
                                        xs,
                                        const ref others : BatchedOperatorPointers,
                                        in timers : BatchedOperatorTimers?) {
  // We compute H|αᵢ⟩ = ∑ⱼ cᵢⱼ|βᵢⱼ⟩ for each |αᵢ⟩ where i ∈ {0 ..# count} and
  // j ∈ {0 ..# numberOffDiagTerms}. Both cᵢⱼ and |βᵢⱼ⟩ are conceptually
  // 2-dimensional arrays, but we use flattened representations of them.
  if timers != nil then timers!.applyOffDiagTimer.start();
  ls_internal_operator_apply_off_diag_x1(
    matrix.payload.deref(),
    count,
    alphas,
    others.betas,
    others.coeffs,
    others.offsets,
    xs);
  if timers != nil then timers!.applyOffDiagTimer.stop();
}
private proc computeOffDiagOnlyInversion(const ref matrix : Operator,
                                         count : int,
                                         alphas : c_ptrConst(uint(64)),
                                         xs,
                                         const ref others : BatchedOperatorPointers,
                                         in timers : BatchedOperatorTimers?) {
  computeOffDiagNoProjection(matrix, count, alphas, xs, others, timers);

  const totalCount = others.offsets[count];
  const mask = (1:uint(64) << matrix.basis.numberSites()) - 1;
  const character = matrix.basis.spinInversion;
  assert(character != 0);

  if timers != nil then timers!.stateInfoTimer.start();
  foreach i in 0 ..# totalCount {
    const current = others.betas[i];
    const inverted = current ^ mask;
    if inverted < current {
      others.betas[i] = inverted;
      others.coeffs[i] *= character;
    }
  }
  if timers != nil then timers!.stateInfoTimer.stop();
}
private proc computeOffDiagWithProjection(const ref matrix : Operator,
                                          count : int,
                                          alphas : c_ptrConst(uint(64)),
                                          xs,
                                          ref others : BatchedOperatorPointers,
                                          in timers : BatchedOperatorTimers?) {
  assert(others.tempNorms != nil && others.tempSpins != nil && others.tempCoeffs != nil);
  // computeOffDiagNoProjection stores its results in others.betas and others.coeffs,
  // but we need them in others.tempSpins and others.tempCoeffs
  others.betas <=> others.tempSpins;
  others.coeffs <=> others.tempCoeffs;
  computeOffDiagNoProjection(matrix, count, alphas, xs, others, timers);
  // revert the swaps
  others.betas <=> others.tempSpins;
  others.coeffs <=> others.tempCoeffs;

  const totalCount = others.offsets[count];
  // We are also interested in norms of alphas, so we append them to tempSpins
  if timers != nil then timers!.memcpyTimer.start();
  POSIX.memcpy(others.tempSpins + totalCount, alphas, count:c_size_t * c_sizeof(uint(64)));
  if timers != nil then timers!.memcpyTimer.stop();

  if timers != nil then timers!.stateInfoTimer.start();
  ls_hs_state_info(
    matrix.basis.payload,
    totalCount + count,
    others.tempSpins, 1,
    others.betas, 1,
    others.coeffs,
    others.tempNorms);
  if timers != nil then timers!.stateInfoTimer.stop();

  if timers != nil then timers!.coeffTimer.start();
  foreach i in 0 ..# count {
    foreach k in others.offsets[i] ..< others.offsets[i + 1] {
      others.coeffs[k] *= others.tempCoeffs[k] * others.tempNorms[k] / others.tempNorms[totalCount + i];
    }
  }
  if timers != nil then timers!.coeffTimer.stop();
}

// Given `count` basis vectors `alphas` with corresponding coefficients `xs`,
// we apply the operator to each basis vector.
//
// For each basis vector, we have:
//
//   O |α⟩ = ∑ⱼ cⱼ|βⱼ⟩
//
// where the number of cⱼ coefficients depends on the basis vector |α⟩, but is
// bounded by `_numberOffDiagTerms`.
private inline proc computeOffDiagGeneric(const ref matrix : Operator,
                                          count : int,
                                          alphas : c_ptrConst(uint(64)),
                                          xs,
                                          ref others : BatchedOperatorPointers,
                                          in timers : BatchedOperatorTimers?) {
  // Simple case when no symmetries are used
  if !matrix.basis.requiresProjection() {
    computeOffDiagNoProjection(matrix, count, alphas, xs, others, timers);
  }
  // Another simple case when no permutations were given, only the spin inversion
  else if !matrix.basis.hasPermutationSymmetries()
            && matrix.basis.hasSpinInversionSymmetry() {
    computeOffDiagOnlyInversion(matrix, count, alphas, xs, others, timers);
  }
  // The tricky case when we have to project betas first
  else {
    computeOffDiagWithProjection(matrix, count, alphas, xs, others, timers);
  }
}

extern {
  #include <stdlib.h>
  #include <stdint.h>

  // int qsort_s( void *ptr, size_t count, size_t size,
  //              int (*comp)(const void *, const void *, void *),
  //              void *context );

  int ls_internal_comp_int64(const void* _a, const void* _b, void* _ctx);
  void ls_internal_qsort_int64(int64_t* xs, size_t count, uint64_t const* keys);

  int ls_internal_comp_int64(const void* _a, const void* _b, void* _ctx) {
    int64_t const a = *(int64_t const*)_a;
    int64_t const b = *(int64_t const*)_b;
    uint64_t const* keys = (uint64_t const*)_ctx;
    uint64_t const key_a = keys[a];
    uint64_t const key_b = keys[b];
    return (key_a > key_b) - (key_a < key_b);
  }

  void ls_internal_qsort_int64(int64_t* xs, size_t const count, uint64_t const* keys) {
    qsort_r(xs, count, sizeof(int64_t), ls_internal_comp_int64, (void*)keys);
  }
}


private proc csrGeneratePart(count : int,
                             betas : c_ptrConst(uint(64)),
                             coeffs : c_ptrConst(complex(128)),
                             offsets : c_ptrConst(int),
                             diag : c_ptrConst(real(64)),
                             rowOffsets : c_ptr(?idxType),
                             colIndices : c_ptr(idxType),
                             matrixElements : c_ptr(?eltType),
                             const ref basis : Basis,
                             numberOffDiagTerms : int) {
  var order : [0 ..# numberOffDiagTerms] int;
  var indices : [0 ..# numberOffDiagTerms] int;
  var numberNonZero : idxType = 0;

  for rowIndex in (0:idxType) ..# (count:idxType) {
    rowOffsets[rowIndex] = numberNonZero;

    const b = offsets[rowIndex];
    const e = offsets[rowIndex + 1];
    const n = e - b;

    // Reset the order
    foreach i in 0 ..# n do order[i] = i;
    // Sort according to betas
    ls_internal_qsort_int64(c_ptrTo(order[0]), n, c_ptrToConst(betas[b]));
    // Convert betas to indices
    if numLocales == 1 && basis.isStateIndexIdentity()
      then POSIX.memcpy(c_ptrTo(indices[0]), c_ptrToConst(betas[b]), n:c_size_t * c_sizeof(int(64)));
      else ls_hs_state_index(basis.payload, n, c_ptrToConst(betas[b]), 1, c_ptrTo(indices[0]), 1);

    // Sum duplicates
    var diagonalWritten = false;
    var k : idxType = 0;
    while k < n {
      var colIndex = indices[order[k]]:idxType;
      if colIndex < 0 then
        halt("invalid index: " + colIndex:string + " for state " + betas[order[k]]:string);

      var acc = coeffs[b + order[k]];
      k += 1;
      while k < n && indices[order[k]] == colIndex {
        acc += coeffs[b + order[k]];
        k += 1;
      }
      if colIndex == rowIndex {
        acc += diag[rowIndex];
        diagonalWritten = true;
      }

      if acc != 0 {
        if !diagonalWritten && colIndex > rowIndex {
          const diagonal = diag[rowIndex];
          if diagonal != 0 {
            colIndices[numberNonZero] = rowIndex;
            matrixElements[numberNonZero] = diagonal;
            numberNonZero += 1;
          }
          diagonalWritten = true;
        }
        colIndices[numberNonZero] = colIndex;
        matrixElements[numberNonZero] = acc;
        numberNonZero += 1;
      }
    }
    if !diagonalWritten {
      const diagonal = diag[rowIndex];
      if diagonal != 0 {
        colIndices[numberNonZero] = rowIndex;
        matrixElements[numberNonZero] = diagonal;
        numberNonZero += 1;
      }
    }
  }
  rowOffsets[count] = numberNonZero;
}


class BatchedOperatorTimers {
  var totalTimer : stopwatch;
  var applyOffDiagTimer : stopwatch;
  var stateInfoTimer : stopwatch;
  var memcpyTimer : stopwatch;
  var coeffTimer : stopwatch;
  var keysTimer : stopwatch;
}

record BatchedOperatorPointers {
  var betas  : c_ptr(uint(64));
  var coeffs : c_ptr(complex(128));
  var offsets : c_ptr(int);
  var tempSpins : c_ptr(uint(64));
  var tempCoeffs : c_ptr(complex(128));
  var tempNorms : c_ptr(real(64));
  var localeIdxs : c_ptr(uint(8));
}

// A wrapper around the Operator class that allows applying the operator to a
// batch of basis vectors instead of one.
//
// This record pre-allocates all buffers such that repeated application cause no allocations.
//
// Operations on this record are not thread-safe, so each task should keep a separate copy of BatchedOperator.
record BatchedOperator {
  var _matrixPtr : c_ptrConst(Operator);
  var pointers : BatchedOperatorPointers;
  var batchSize : int;
  var numberOffDiagTerms : int;
  var _dom : domain(1);
  var _spins : [_dom] uint(64);
  var _coeffs : [_dom] complex(128);
  var _offsets : [0 ..# batchSize + 1] int;
  var _localeIdxs : [_dom] uint(8);
  var _domTemp : domain(1);
  var _spinsTemp : [_domTemp] uint(64);
  var _coeffsTemp : [_domTemp] complex(128);
  var _normsTemp : [_domTemp] real(64);
  var timers : owned BatchedOperatorTimers;

  proc initPointers() {
    pointers.betas = c_ptrTo(_spins);
    pointers.coeffs = c_ptrTo(_coeffs);
    pointers.offsets = c_ptrTo(_offsets);
    pointers.tempSpins = if _domTemp.size != 0 then c_ptrTo(_spinsTemp) else nil;
    pointers.tempCoeffs = if _domTemp.size != 0 then c_ptrTo(_coeffsTemp) else nil;
    pointers.tempNorms = if _domTemp.size != 0 then c_ptrTo(_normsTemp) else nil;
    pointers.localeIdxs = c_ptrTo(_localeIdxs);
  }

  proc init(const ref matrix : Operator, batchSize : int) {
    this._matrixPtr = c_addrOfConst(matrix);
    this.batchSize = batchSize;
    this.numberOffDiagTerms = matrix.numberOffDiagTerms();
    const numberTerms = max(numberOffDiagTerms, 1);
    this._dom = {0 ..# batchSize * numberTerms};
    this._domTemp =
      if _matrixPtr.deref().basis.hasPermutationSymmetries()
        // It is important to have numberTerms + 1 here. See the POSIX.memcpy
        // call in computeOffDiagWithProjection for why.
        then {0 ..# (batchSize * (numberTerms + 1))}
        else {0 .. -1};
    this.timers = new BatchedOperatorTimers();
    complete();
    initPointers();
  }
  proc init=(const ref other : BatchedOperator) {
    assert(other.locale == here);
    this._matrixPtr = other._matrixPtr;
    this.batchSize = other.batchSize;
    this.numberOffDiagTerms = other.numberOffDiagTerms;
    this._dom = other._dom;
    this._domTemp = other._domTemp;
    this.timers = new BatchedOperatorTimers();
    complete();
    initPointers();
  }

  inline proc matrix ref { return _matrixPtr.deref(); }
  inline proc basis ref { return matrix.basis; }

  proc computeOffDiag(count : int, alphas : c_ptrConst(uint(64)), xs)
      : (int, c_ptr(uint(64)), c_ptr(complex(128)), c_ptr(uint(8))) {

    assert(count <= batchSize, "buffer overflow in BatchedOperator");
    if count == 0 {
      pointers.offsets[0] = 0;
      return (0, nil, nil, nil);
    }
    if numberOffDiagTerms == 0 {
      POSIX.memset(pointers.offsets, 0, (count + 1):c_size_t * c_sizeof(int));
      return (0, nil, nil, nil);
    }

    timers.totalTimer.start();
    defer timers.totalTimer.stop();

    computeOffDiagGeneric(matrix, count, alphas, xs, pointers, timers.borrow());
    const totalCount = pointers.offsets[count];
    computeKeys(totalCount, pointers.betas, pointers.localeIdxs, timers.keysTimer);

    return (totalCount, pointers.betas, pointers.coeffs, pointers.localeIdxs);
  }

  proc getTimings() {
    var tree = timingTree(
      "BatchedOperator.computeOffDiag", timers.totalTimer.elapsed(),
      [ ("ls_internal_operator_apply_off_diag_x1", timers.applyOffDiagTimer.elapsed())
      , ("memcpy", timers.memcpyTimer.elapsed())
      , ("ls_hs_state_info", timers.stateInfoTimer.elapsed())
      , ("rescale coeffs", timers.coeffTimer.elapsed())
      , ("localeIdxOf", timers.keysTimer.elapsed())
      ]
    );
    return tree;
  }
}

export proc ls_chpl_operator_apply_diag(matrixPtr : c_ptr(ls_hs_operator),
                                        count : int,
                                        alphas : c_ptr(uint(64)),
                                        coeffs : c_ptr(chpl_external_array),
                                        numTasks : int) {
  var matrix = new Operator(matrixPtr, owning=false);
  if matrix.basis.numberWords != 1 then
    halt("bases with more than 64 bits are not yet implemented");

  var _cs : [0 ..# count] real(64) = noinit;
  ls_internal_operator_apply_diag_x1(
    matrix.payload.deref(), count, alphas, c_ptrTo(_cs[0]), nil);

  coeffs.deref() = convertToExternalArray(_cs);
}

export proc ls_chpl_operator_apply_off_diag(matrixPtr : c_ptr(ls_hs_operator),
                                            count : int,
                                            alphas : c_ptr(uint(64)),
                                            betas : c_ptr(chpl_external_array),
                                            coeffs : c_ptr(chpl_external_array),
                                            offsets : c_ptr(chpl_external_array),
                                            numTasks : int) {
  assert(matrixPtr != nil);
  var matrix = new Operator(matrixPtr, owning=false);
  var batchedOperator = new BatchedOperator(matrix, max(1, count));

  batchedOperator.computeOffDiag(count, alphas, nil);
  betas.deref() = convertToExternalArray(batchedOperator._spins);
  coeffs.deref() = convertToExternalArray(batchedOperator._coeffs);
  offsets.deref() = convertToExternalArray(batchedOperator._offsets);
}

export proc ls_chpl_operator_to_csr(matrixPtr : c_ptr(ls_hs_operator),
                                    rowOffsets : c_ptr(chpl_external_array),
                                    colIndices : c_ptr(chpl_external_array),
                                    matrixElements : c_ptr(chpl_external_array),
                                    numTasks : int) {
  assert(matrixPtr != nil);
  var matrix = new Operator(matrixPtr, owning=false);
  const ref basis = matrix.basis;
  const ref representatives = basis.representatives();
  const count = representatives.size;
  const representativesPtr = c_ptrToConst(representatives);
  var batchedOperator = new BatchedOperator(matrix, max(1, count));

  const (totalCount, betas, coeffs, _keys) =
    batchedOperator.computeOffDiag(count, representativesPtr, nil);

  var _diagonal : [0 ..# count] real(64) = noinit;
  ls_internal_operator_apply_diag_x1(
    matrixPtr.deref(), count, representativesPtr, c_ptrTo(_diagonal), nil);

  var _rowOffsets : [0 ..# count + 1] int(32);
  var _colIndices : [0 ..# count + totalCount] int(32);
  var _matrixElements : [0 ..# count + totalCount] complex(128);

  csrGeneratePart(count,
                  betas,
                  coeffs,
                  batchedOperator.pointers.offsets,
                  c_ptrTo(_diagonal),
                  c_ptrTo(_rowOffsets),
                  c_ptrTo(_colIndices),
                  c_ptrTo(_matrixElements),
                  basis,
                  batchedOperator.numberOffDiagTerms);

  rowOffsets.deref() = convertToExternalArray(_rowOffsets);
  colIndices.deref() = convertToExternalArray(_colIndices);
  matrixElements.deref() = convertToExternalArray(_matrixElements);
}

export proc ls_chpl_matrix_vector_product_csr_i32_c128(
    numberRows : int, numberCols : int,
    numberNonZero : int, matrixElements : c_ptrConst(complex(128)),
    rowOffsets : c_ptrConst(int(32)), colIndices : c_ptrConst(int(32)),
    x : c_ptrConst(complex(128)), y : c_ptr(complex(128)), numTasks : int) {
  const matrix = new CSR.CSR(eltType=complex(128), idxType=int(32),
                             matrixElements=matrixElements,
                             rowOffsets=rowOffsets,
                             colIndices=colIndices,
                             numberRows=numberRows,
                             numberCols=numberCols,
                             numberNonZero=numberNonZero);
  CSR.csrMatvec(matrix, x, y, if numTasks <= 0 then here.maxTaskPar else numTasks);
}

} // end module BatchedOperator
