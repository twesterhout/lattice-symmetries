module BatchedOperator {

use FFI;
use ForeignTypes;
use Timing;
use Utils;

import IO.FormattedIO.format;
import OS.POSIX;
import Reflection.getRoutineName;
use BitOps;
use CTypes;

// require "library.h";

proc applyDiagKernel(
    diag_terms : c_ptrConst(ls_hs_nonbranching_terms),
    batch_size : int,
    alphas : c_ptrConst(uint(64)),
    ys : c_ptr(?eltType),
    xs : c_ptrConst(eltType)) {
  // The diagonal is zero
  if (diag_terms == nil || diag_terms.deref().number_terms == 0) {
    POSIX.memset(ys, 0, batch_size:c_size_t * c_sizeof(eltType));
    return;
  }

  const ref terms = diag_terms.deref();
  const number_terms = terms.number_terms;

  foreach batch_idx in 0 ..# batch_size {
    var acc : eltType = 0;
    const alpha = alphas[batch_idx];
    for term_idx in 0 ..# number_terms {
      const delta = (alpha & terms.m[term_idx]) == terms.r[term_idx];
      if delta {
        const sign = 1 - 2 * parity(alpha & terms.s[term_idx]):eltType;
        const factor = if xs != nil then sign * xs[batch_idx] else sign;
        acc += terms.v[term_idx]:eltType * factor;
      }
    }
    ys[batch_idx] = acc;
  }
}

proc applyOffDiagKernel(off_diag_terms : c_ptrConst(ls_hs_nonbranching_terms),
                        chunk : range(int),
                        alphas : c_ptrConst(uint(64)),
                        betas : c_ptr(uint(64)),
                        coeffs : c_ptr(complex(128)),
                        offsets : c_ptr(c_ptrdiff),
                        targetStates : c_ptr(uint(64)),
                        estimatedNumberTerms : int,
                        param left : bool = false) {

  const batch_size = chunk.size;
  // Nothing to apply
  if (off_diag_terms == nil || off_diag_terms.deref().number_terms == 0) {
    POSIX.memset(offsets, 0, (batch_size + 1):c_size_t * c_sizeof(c_ptrdiff));
    return;
  }

  const ref terms = off_diag_terms.deref();
  const number_terms = terms.number_terms;

  offsets[0] = 0;
  var offset = 0;
  for batch_idx in 0 ..# batch_size {
    const oldOffset = offset;
    const alpha = alphas[batch_idx];

    var term_idx = 0;
    do {
      const tX = terms.x[term_idx];
      var acc : complex(128) = 0;

      do {
        const l_or_r = if left then terms.l[term_idx] else terms.r[term_idx];
        const delta = (alpha & terms.m[term_idx]) == l_or_r;
        if delta {
          const sign = 1 - 2 * parity(alpha & terms.s[term_idx]):real;
          const coeff = terms.v[term_idx] * sign; //  + (terms.v[term_idx].im * sign) * 1.0i;
          // const beta = alpha ^ tX;
          acc += coeff;

          // if beta != old_beta {
          //   if old_beta != alpha {
          //     coeffs[offset] = old_coeff;
          //     betas[offset] = old_beta;
          //     targetStates[offset] = chunk.low + batch_idx;
          //     offset += 1;
          //   }
          //   old_beta = beta;
          //   old_coeff = coeff;
          // }
          // else {
          //   old_coeff += coeff;
          // }
        }
        term_idx += 1;
      } while term_idx < number_terms && terms.x[term_idx] == tX;

      if acc != 0 {
        coeffs[offset] = acc;
        betas[offset] = alpha ^ tX;
        targetStates[offset] = chunk.low + batch_idx;
        offset += 1;
      }
    } while term_idx < number_terms;

    // TODO: this check is important, because estimatedNumberTerms is an estimate !!
    if offset - oldOffset > estimatedNumberTerms then
      halt(try! "buffer overflow in applyOffDiagKernel: estimatedNumberTerms=%i, but written %i".format(
                estimatedNumberTerms, offset - oldOffset));

    offsets[batch_idx + 1] = offset;

    /*
    const oldOffset = offset;
    const alpha = alphas[batch_idx];
    // we know that we can never get alpha again since all terms are off
    // diagonal, so using alpha as a special state is fine.
    var old_beta = alpha;
    var old_coeff : complex(128) = 0;
    for term_idx in 0 ..# number_terms {
      // Right: T|α> = v * δ_{α^m,r} * (−1)^{h(α^s)} |α⊕x>,
      // Left:  <α|T = v * δ_{α^m,l} * (−1)^{h(α^s)} <α⊕x|,
      const l_or_r = if left then terms.l[term_idx] else terms.r[term_idx];
      const delta = (alpha & terms.m[term_idx]) == l_or_r;
      if delta {
        const sign = 1 - 2 * parity(alpha & terms.s[term_idx]):real;
        const coeff = terms.v[term_idx].re * sign + (terms.v[term_idx].im * sign) * 1.0i;
        const beta = alpha ^ terms.x[term_idx];
        if beta != old_beta {
          if old_beta != alpha {
            coeffs[offset] = old_coeff;
            betas[offset] = old_beta;
            targetStates[offset] = chunk.low + batch_idx;
            offset += 1;
          }
          old_beta = beta;
          old_coeff = coeff;
        }
        else {
          old_coeff += coeff;
        }
      }
    }
    if old_beta != alpha {
      coeffs[offset] = old_coeff;
      betas[offset] = old_beta;
      targetStates[offset] = chunk.low + batch_idx;
      offset += 1;
    }
    if offset - oldOffset > estimatedNumberTerms {
      for k in oldOffset .. offset - 1 {
        writeln(betas[k], ", ", coeffs[k]);
      }
      writeln("alpha=", alpha);

      halt(try! "buffer overflow in applyOffDiagKernel: estimatedNumberTerms=%i, but written %i".format(
                estimatedNumberTerms, offset - oldOffset));
    }
    offsets[batch_idx + 1] = offset;
    */
  }
}

/*
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
*/

private proc computeOffDiagNoProjection(const ref matrix : ls_chpl_batched_operator,
                                        chunk,
                                        alphas : c_ptrConst(uint(64)),
                                        param left : bool) {
  const _timer = recordTime(getRoutineName());
  const count = chunk.size;
  if count > matrix.batch_size then
    halt(try! "buffer overflow: allocated space for %i elements, but count=%i".format(matrix.batch_size, count));
  if matrix.betas == nil || matrix.coeffs == nil || matrix.offsets == nil then
    halt("betas, coeffs, or offsets should be pre-allocated");

  const estimatedNumberTerms = matrix.matrix.deref().max_number_off_diag_estimate;
  applyOffDiagKernel(
    matrix.matrix.deref().off_diag_terms,
    chunk,
    alphas,
    matrix.betas,
    matrix.coeffs,
    matrix.offsets,
    matrix.target_states,
    estimatedNumberTerms,
    left);
}
private proc computeOffDiagOnlyInversion(const ref matrix : ls_chpl_batched_operator,
                                         chunk,
                                         alphas : c_ptrConst(uint(64)),
                                         param left : bool) {
  const _timer = recordTime("computeOffDiagOnlyInversion");
  // computeOffDiagNoProjection(matrix, chunk, alphas, xs, left);
  assert(false,
         "computeOffDiagOnlyInversion is not yet implemented");

  // const totalCount = matrix.offsets[count];
  // const mask = (1:uint(64) << matrix.basis.info.number_sites) - 1;
  // const character = matrix.basis.info.spin_inversion;
  // assert(character != 0);

  // if timers != nil then timers!.stateInfoTimer.start();
  // foreach i in 0 ..# totalCount {
  //   const current = others.betas[i];
  //   const inverted = current ^ mask;
  //   if inverted < current {
  //     others.betas[i] = inverted;
  //     others.coeffs[i] *= character;
  //   }
  // }
  // if timers != nil then timers!.stateInfoTimer.stop();
}
private proc computeOffDiagWithProjection(const ref matrix : ls_chpl_batched_operator,
                                          chunk,
                                          alphas : c_ptrConst(uint(64)),
                                          param left : bool) {
  const _timer = recordTime("computeOffDiagWithProjection");
  assert(false,
         "computeOffDiagWithProjection is not yet implemented");
  /*
  assert(others.tempNorms != nil && others.tempSpins != nil && others.tempCoeffs != nil);
  // computeOffDiagNoProjection stores its results in others.betas and others.coeffs,
  // but we need them in others.tempSpins and others.tempCoeffs
  others.betas <=> others.tempSpins;
  others.coeffs <=> others.tempCoeffs;
  computeOffDiagNoProjection(matrix, count, alphas, xs, others, left, timers);
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
      const c = if left then conj(others.coeffs[k]) else others.coeffs[k];
      others.coeffs[k] = c * others.tempCoeffs[k] * others.tempNorms[k] / others.tempNorms[totalCount + i];
    }
  }
  if timers != nil then timers!.coeffTimer.stop();
  */
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
private inline proc computeOffDiagGeneric(const ref basisInfo : ls_hs_basis_info,
                                          const ref matrix : ls_chpl_batched_operator,
                                          chunk,
                                          alphas : c_ptrConst(uint(64)),
                                          param left : bool) {
  // Simple case when no symmetries are used
  if !basisInfo.requires_projection {
    computeOffDiagNoProjection(matrix, chunk, alphas, left);
  }
  // Another simple case when no permutations were given, only the spin inversion
  else if !basisInfo.has_permutation_symmetries && basisInfo.spin_inversion != 0 {
    computeOffDiagOnlyInversion(matrix, chunk, alphas, left);
  }
  // The tricky case when we have to project betas first
  else {
    computeOffDiagWithProjection(matrix, chunk, alphas, left);
  }
}

/*
extern proc ls_internal_qsort_int64(xs : c_ptr(int(64)), count : c_size_t, keys : c_ptrConst(uint(64)));

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
  logDebug("numberOffDiagTerms=", numberOffDiagTerms);
  var order : [0 ..# numberOffDiagTerms] int;
  var indices : [0 ..# numberOffDiagTerms] int;
  var numberNonZero : idxType = 0;

  for rowIndex in (0:idxType) ..# (count:idxType) {
    rowOffsets[rowIndex] = numberNonZero;

    const b = offsets[rowIndex];
    const e = offsets[rowIndex + 1];
    const n = e - b;
    if numberOffDiagTerms == 0 then assert(n == 0);

    if n > 0 {
      // Reset the order
      foreach i in 0 ..# n do order[i] = i;
      // Sort according to betas
      ls_internal_qsort_int64(c_ptrTo(order[0]), n, c_ptrToConst(betas[b]));
      // Convert betas to indices
      if numLocales == 1 && basis.info.is_state_index_identity
        then POSIX.memcpy(c_ptrTo(indices[0]), c_ptrToConst(betas[b]), n:c_size_t * c_sizeof(int(64)));
        else ls_hs_state_index(basis.payload, n, c_ptrToConst(betas[b]), 1, c_ptrTo(indices[0]), 1);
    }

    // Sum duplicates
    var diagonalWritten = false;
    var k : idxType = 0;
    while k < n {
      var colIndex = indices[order[k]]:idxType;

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
        else if colIndex < 0 {
          halt("invalid index: " + colIndex:string + " for state " + betas[b + order[k - 1]]:string);
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
*/

/*
class BatchedOperatorTimers {
  var totalTimer : stopwatch;
  var applyOffDiagTimer : stopwatch;
  var stateInfoTimer : stopwatch;
  var memcpyTimer : stopwatch;
  var coeffTimer : stopwatch;
  var keysTimer : stopwatch;
}
*/

/*
record BatchedOperatorPointers {
  var betas  : c_ptr(uint(64));
  var coeffs : c_ptr(complex(128));
  var offsets : c_ptr(int);
  var tempSpins : c_ptr(uint(64));
  var tempCoeffs : c_ptr(complex(128));
  var tempNorms : c_ptr(real(64));
  var localeIdxs : c_ptr(uint(8));
}
*/

// A wrapper around the Operator class that allows applying the operator to a
// batch of basis vectors instead of one.
//
// This record pre-allocates all buffers such that repeated application cause no allocations.
//
// Operations on this record are not thread-safe, so each task should keep a separate copy of BatchedOperator.
record BatchedOperator {
  var raw : ls_chpl_batched_operator;
  var owning : bool;
  var maxTotalCount : int;

  proc matrix const ref : ls_hs_operator {
    assert(raw.matrix != nil);
    return raw.matrix.deref();
  }
  proc basis const ref : ls_hs_basis {
    assert(raw.basis != nil);
    return matrix.basis.deref();
  }
  proc basisInfo const ref : ls_hs_basis_info {
    return getBasisInfo(matrix.basis);
  }

  proc ref allocateBuffers() {
    const safeNumberTerms = max(matrix.max_number_off_diag, 1);
    const estimatedNumberTerms = max(matrix.max_number_off_diag_estimate, 1);
    // It is important to have `+ safeNumberTerms` here. See the POSIX.memcpy
    // call in computeOffDiagWithProjection for why.
    const capacity = (raw.batch_size - 1) * estimatedNumberTerms + safeNumberTerms // Always needed
                     + safeNumberTerms; // Specifically for symmetrized bases

    raw.betas = allocate(uint(64), capacity);
    raw.coeffs = allocate(complex(128), capacity);
    raw.target_states = allocate(uint(64), capacity);
    raw.offsets = allocate(int(64), raw.batch_size + 1);
    if basisInfo.has_permutation_symmetries {
      raw.temp_spins = allocate(uint(64), capacity);
      raw.temp_coeffs = allocate(complex(128), capacity);
      raw.temp_norms = allocate(real(64), capacity);
    }
    if numLocales > 1 {
      raw.locale_indices = allocate(uint(8), capacity);
    }
  }

  proc ref deallocateBuffers() {
    // logDebug("BatchedOperator.deallocateBuffers");
    if raw.betas != nil then deallocate(raw.betas);
    if raw.coeffs != nil then deallocate(raw.coeffs);
    if raw.target_states != nil then deallocate(raw.target_states);
    if raw.offsets != nil then deallocate(raw.offsets);
    if raw.locale_indices != nil then deallocate(raw.locale_indices);
    if raw.temp_spins != nil then deallocate(raw.temp_spins);
    if raw.temp_coeffs != nil then deallocate(raw.temp_coeffs);
    if raw.temp_norms != nil then deallocate(raw.temp_norms);
  }

  proc init(raw : ls_chpl_batched_operator) {
    this.raw = raw;
    this.owning = false;
  }
  proc init(matrix : c_ptrConst(ls_hs_operator), batch_size : int) {
    this.raw = new ls_chpl_batched_operator(matrix, batch_size, nil, nil, nil, nil, nil, nil, nil);
    this.owning = true;
    init this;
    allocateBuffers();
  }

  proc ref deinit() {
    // logDebug("BatchedOperator.deinit");
    if owning then deallocateBuffers();
  }

  proc ref computeOffDiag(chunk : range(int), alphas : c_ptrConst(uint(64)), param left : bool)
      : (int, c_ptr(uint(64)), c_ptr(complex(128)), c_ptr(uint(64))) {
    const count = chunk.size;
    if count == 0 || matrix.max_number_off_diag == 0 then
      return (0, nil, nil, nil);
    if count > raw.batch_size then
      halt(try! "buffer overflow in BatchedOperator: count=%i, batch_size=%i".format(count, raw.batch_size));

    computeOffDiagGeneric(basisInfo, raw, chunk, alphas, left);
    const totalCount = raw.offsets[count];
    // if totalCount > maxTotalCount then
    //   maxTotalCount = totalCount;
    // computeKeys(totalCount, pointers.betas, pointers.localeIdxs, timers.keysTimer);
    return (totalCount, raw.betas, raw.coeffs, raw.target_states);
  }
}

// class BatchedOperator {
//   var _matrixPtr : c_ptrConst(Operator);
//   var pointers : BatchedOperatorPointers;
//   var batchSize : int;
//   var numberOffDiagTerms : int;
//   var hasPermutationSymmetries : bool;
//   var timers : owned BatchedOperatorTimers;
// 
//   proc initPointers() {
//     const numberTerms = max(numberOffDiagTerms, 1);
//     // It is important to have numberTerms + 1 here. See the POSIX.memcpy
//     // call in computeOffDiagWithProjection for why.
//     const capacity = batchSize * (numberTerms + 1);
// 
//     pointers.offsets = allocate(int, batchSize + 1);
//     pointers.betas = allocate(uint(64), capacity);
//     pointers.coeffs = allocate(complex(128), capacity);
//     pointers.localeIdxs = allocate(uint(8), capacity);
//     if hasPermutationSymmetries {
//       pointers.tempSpins = allocate(uint(64), capacity);
//       pointers.tempCoeffs = allocate(complex(128), capacity);
//       pointers.tempNorms = allocate(real(64), capacity);
//     }
//   }
// 
//   proc deinitPointers() {
//     if pointers.offsets != nil then deallocate(pointers.offsets);
//     if pointers.betas != nil then deallocate(pointers.betas);
//     if pointers.coeffs != nil then deallocate(pointers.coeffs);
//     deallocate(pointers.localeIdxs);
//     if hasPermutationSymmetries {
//       deallocate(pointers.tempSpins);
//       deallocate(pointers.tempCoeffs);
//       deallocate(pointers.tempNorms);
//     }
//   }
// 
//   proc init(const ref matrix : Operator, batchSize : int) {
//     this._matrixPtr = c_addrOfConst(matrix);
//     this.batchSize = batchSize;
//     this.numberOffDiagTerms = matrix.numberOffDiagTerms();
//     this.hasPermutationSymmetries = _matrixPtr.deref().basis.info.has_permutation_symmetries;
//     this.timers = new BatchedOperatorTimers();
//     init this;
//     initPointers();
//   }
// 
//   proc deinit() { deinitPointers(); }
// 
//   inline proc matrix ref { return _matrixPtr.deref(); }
//   inline proc basis ref { return matrix.basis; }
// 
//   proc computeOffDiag(count : int, alphas : c_ptrConst(uint(64)), xs, param left : bool)
//       : (int, c_ptr(uint(64)), c_ptr(complex(128)), c_ptr(uint(8))) {
// 
//     assert(count <= batchSize, "buffer overflow in BatchedOperator");
//     if count == 0 {
//       pointers.offsets[0] = 0;
//       return (0, nil, nil, nil);
//     }
//     if numberOffDiagTerms == 0 {
//       POSIX.memset(pointers.offsets, 0, (count + 1):c_size_t * c_sizeof(int));
//       return (0, nil, nil, nil);
//     }
// 
//     timers.totalTimer.start();
//     defer timers.totalTimer.stop();
// 
//     computeOffDiagGeneric(matrix, count, alphas, xs, pointers, left, timers.borrow());
//     const totalCount = pointers.offsets[count];
//     computeKeys(totalCount, pointers.betas, pointers.localeIdxs, timers.keysTimer);
// 
//     return (totalCount, pointers.betas, pointers.coeffs, pointers.localeIdxs);
//   }
// 
//   proc getTimings() {
//     var tree = timingTree(
//       "BatchedOperator.computeOffDiag", timers.totalTimer.elapsed(),
//       [ ("applyOffDiagKernel", timers.applyOffDiagTimer.elapsed())
//       , ("memcpy", timers.memcpyTimer.elapsed())
//       , ("ls_hs_state_info", timers.stateInfoTimer.elapsed())
//       , ("rescale coeffs", timers.coeffTimer.elapsed())
//       , ("localeIdxOf", timers.keysTimer.elapsed())
//       ]
//     );
//     return tree;
//   }
// }

/*
export proc ls_chpl_operator_apply_diag(matrixPtr : c_ptrConst(ls_hs_operator),
                                        count : int,
                                        alphas : c_ptrConst(uint(64)),
                                        coeffs : c_ptr(chpl_external_array)) {
  assert(matrixPtr != nil);
  const matrix = matrixPtr.deref();
  const basisInfo = getBasisInfo(matrix.basis);
  if basisInfo.number_words != 1 then
    halt(try! ("ls_chpl_operator_apply_diag: number_words=%i, "
              + "but bases with more than 64 bits are not yet implemented").format(basisInfo.number_words));

  var cs = allocate(real(64), count);
  applyDiagKernel(matrix.diag_terms, count, alphas, cs, nil);
  coeffs.deref() = chpl_make_external_array_ptr_free(cs, count);
}

export proc ls_chpl_operator_apply_off_diag(matrixPtr : c_ptrConst(ls_hs_operator),
                                            count : int,
                                            alphas : c_ptrConst(uint(64)),
                                            betas : c_ptr(chpl_external_array),
                                            coeffs : c_ptr(chpl_external_array),
                                            offsets : c_ptr(chpl_external_array)) {
  assert(matrixPtr != nil);
  var batchedOperator = new BatchedOperator(matrixPtr, count);
  const ref basisInfo = batchedOperator.basisInfo;
  if basisInfo.number_words != 1 then
    halt(try! ("ls_chpl_operator_apply_off_diag: number_words=%i, "
              + "but bases with more than 64 bits are not yet implemented").format(basisInfo.number_words));

  const (size, _, _, _) = batchedOperator.computeOffDiag(0 ..# count, alphas, nil, left=false);
  betas.deref() = chpl_make_external_array_ptr_free(batchedOperator.raw.betas, size);
  batchedOperator.raw.betas = nil;
  coeffs.deref() = chpl_make_external_array_ptr_free(batchedOperator.raw.coeffs, size);
  batchedOperator.raw.coeffs = nil;
  offsets.deref() = chpl_make_external_array_ptr_free(batchedOperator.raw.offsets, count + 1);
  batchedOperator.raw.offsets = nil;
}
*/

/*
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
    batchedOperator.computeOffDiag(count, representativesPtr, nil, left=false);

  var _diagonal : [0 ..# count] real(64) = noinit;
  ls_internal_operator_apply_diag_x1(
    matrixPtr.deref(), count, representativesPtr, c_ptrTo(_diagonal), nil);

  var _rowOffsets : [0 ..# count + 1] int(32);
  var _colIndices : [0 ..# count + totalCount] int(32);
  var _matrixElements : [0 ..# count + totalCount] complex(128);

  logDebug(count);
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
*/

} // end module BatchedOperator
