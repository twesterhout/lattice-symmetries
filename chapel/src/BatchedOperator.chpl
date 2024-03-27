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

proc applyOffDiagKernel(const ref matrix : ls_chpl_batched_operator,
                        chunk : range(int),
                        alphas : c_ptrConst(uint(64)),
                        param left : bool = false) {

  const batch_size = chunk.size;
  if batch_size > matrix.batch_size then
    halt(try! "buffer overflow: allocated space for %i elements, but chunk.size=%i".format(matrix.batch_size, batch_size));
  if matrix.betas == nil || matrix.coeffs == nil || matrix.offsets == nil then
    halt("betas, coeffs, and offsets should be pre-allocated");

  const off_diag_terms = matrix.matrix.deref().off_diag_terms;
  // Nothing to apply
  if (off_diag_terms == nil || off_diag_terms.deref().number_terms == 0) {
    POSIX.memset(matrix.offsets, 0, batch_size:c_size_t * c_sizeof(c_ptrdiff));
    return;
  }

  const ref terms = off_diag_terms.deref();
  const number_terms = terms.number_terms;
  const estimatedNumberTerms = matrix.matrix.deref().max_number_off_diag_estimate;
  const ref basisInfo = getBasisInfo(matrix.matrix.deref().basis);
  const spinInversionMask = (1:uint(64) << basisInfo.number_sites) - 1;
  const spinInversionCharacter = basisInfo.spin_inversion;

  // matrix.offsets[0] = 0;
  if batch_size > 0 && alphas == nil then halt("applyOffDiagKernel received null alphas");
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
          const coeff = terms.v[term_idx] * sign;
          acc += coeff;
          // writeln(try! "acc += %r".format(coeff.re));
        }
        term_idx += 1;
      } while term_idx < number_terms && terms.x[term_idx] == tX;

      if acc != 0 {
        var beta = alpha ^ tX;
        if spinInversionCharacter != 0 {
          const inverted = beta ^ spinInversionMask;
          if inverted < beta {
            beta = inverted;
            acc *= spinInversionCharacter;
            // writeln(try! "acc *= %i".format(spinInversionCharacter));
          }
        }
        matrix.coeffs[offset] = acc;
        matrix.betas[offset] = beta;
        if matrix.target_indices != nil then
          matrix.target_indices[offset] = chunk.low + batch_idx;
        // writeln(try! "y[%i] += %r <%u|x>".format(chunk.low + batch_idx, acc.re, beta));
        offset += 1;
      }
    } while term_idx < number_terms;

    // TODO: this check is important, because estimatedNumberTerms is an estimate !!
    if offset - oldOffset > estimatedNumberTerms then
      halt(try! "buffer overflow in applyOffDiagKernel: estimatedNumberTerms=%i, but written %i".format(
                estimatedNumberTerms, offset - oldOffset));

    matrix.offsets[batch_idx] = offset;
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

// Given `count` basis vectors `alphas` with corresponding coefficients `xs`,
// we apply the operator to each basis vector.
//
// For each basis vector, we have:
//
//   O |α⟩ = ∑ⱼ cⱼ|βⱼ⟩
//
// where the number of cⱼ coefficients depends on the basis vector |α⟩, but is
// bounded by `_numberOffDiagTerms`.
proc _computeOffDiag(const ref basisInfo : ls_hs_basis_info,
                     ref matrix : ls_chpl_batched_operator,
                     chunk,
                     alphas : c_ptrConst(uint(64)),
                     norms : c_ptrConst(uint(16)),
                     param left : bool) {
  const _timer = recordTime(getRoutineName());
  // Simple case when no symmetries are used
  if !basisInfo.has_permutation_symmetries {
    applyOffDiagKernel(matrix, chunk, alphas, left);
  }
  // The tricky case when we have to project betas first
  else {
    // logDebug("_computeOffDiag");
    if matrix.temp_spins == nil || matrix.temp_group_indices == nil then
      halt("temp_spins, temp_norms, and temp_coeffs should be pre-allocated");
    if basisInfo.characters == nil then
      halt("info->characters should not be NULL...");

    // applyOffDiagKernel stores its results in matrix.betas and matrix.coeffs,
    // but we need them in matrix.temp_spins.
    matrix.betas <=> matrix.temp_spins;
    applyOffDiagKernel(matrix, chunk, alphas, left);
    // undo the swaps
    matrix.betas <=> matrix.temp_spins;
    // logDebug("swapping worked");

    const totalCount = matrix.offsets[chunk.size - 1];
    // for k in 0 ..# chunk.size {
    //   logDebug("offsets[", k, "]=", matrix.offsets[k]);
    // }
    // logDebug("totalCount=", totalCount);
    const kernel = ls_chpl_get_state_info_kernel(matrix.matrix.deref().basis);
    if kernel == nil then
      halt("ls_chpl_get_state_info_kernel returned NULL");
    const countWithPadding = roundUpToMaxBlockSize(totalCount);
    if countWithPadding > 0 then
      ls_chpl_invoke_state_info_kernel(kernel, countWithPadding, matrix.temp_spins,
                                      matrix.betas, matrix.temp_group_indices);
    // logDebug("state info worked");

    const numCharacters = basisInfo.number_characters;
    const characters : c_ptrConst(complex(128)) = basisInfo.characters;
    if chunk.size > 0 {
      if norms == nil then halt("_computeOffDiag expects pre-computed norms");
      // logDebug("start norms...");
      foreach k in 0 ..< matrix.offsets[0] {
        const c = matrix.coeffs[k];
        const g = matrix.temp_group_indices[k];
        const character =
          if g > numCharacters then characters[g - numCharacters] * basisInfo.spin_inversion
                               else characters[g];
        // logDebug(try! "k=%i, c=%r + %r im, characters[%i] = %r + %r im, norms[0] = %r".format(k, c.re, c.im, g, characters[g].re, characters[g].im, norms[0]));
        matrix.coeffs[k] = c * character / sqrt(norms[0]:real);
      }
      foreach i in 1 ..< chunk.size {
        foreach k in matrix.offsets[i - 1] ..< matrix.offsets[i] {
          // assert(matrix.coeffs != nil && matrix.temp_group_indices != nil);
          const c = matrix.coeffs[k];
          const g = matrix.temp_group_indices[k];
          const character =
            if g > numCharacters then characters[g - numCharacters] * basisInfo.spin_inversion
                                else characters[g];
          // logDebug(try! "i=%i, k=%i, c=%r + %r im, characters[%i] = %r + %r im, norms[i] = %r".format(i, k, c.re, c.im, g, characters[g].re, characters[g].im, norms[i]));
          // assert(characters != nil && norms != nil);
          matrix.coeffs[k] = c * character / sqrt(norms[i]:real);
        }
      }
      // logDebug("done with norms...");
    }
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
    // NOTE: We round capacity up to a multiple of LS_HS_MAX_BLOCK_SIZE...
    const capacity = roundUpToMaxBlockSize((raw.batch_size - 1) * estimatedNumberTerms + safeNumberTerms);

    raw.betas = allocate(uint(64), capacity);
    raw.coeffs = allocate(complex(128), capacity);
    raw.target_indices = allocate(uint(64), capacity);
    raw.offsets = allocate(int(64), raw.batch_size);
    if basisInfo.has_permutation_symmetries {
      raw.temp_spins = allocate(uint(64), capacity);
      raw.temp_group_indices = allocate(int(32), capacity);
    }
    if numLocales > 1 {
      raw.locale_indices = allocate(uint(8), capacity);
    }
  }

  proc ref deallocateBuffers() {
    // logDebug("BatchedOperator.deallocateBuffers");
    if raw.betas != nil then deallocate(raw.betas);
    if raw.coeffs != nil then deallocate(raw.coeffs);
    if raw.target_indices != nil then deallocate(raw.target_indices);
    if raw.offsets != nil then deallocate(raw.offsets);
    if raw.locale_indices != nil then deallocate(raw.locale_indices);
    if raw.temp_spins != nil then deallocate(raw.temp_spins);
    if raw.temp_group_indices != nil then deallocate(raw.temp_group_indices);
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

  proc ref computeOffDiag(chunk : range(int), alphas : c_ptrConst(uint(64)), norms : c_ptrConst(uint(16)), param left : bool)
      : (int, c_ptr(uint(64)), c_ptr(complex(128)), c_ptr(uint(64))) {

    const count = chunk.size;
    if count == 0 || matrix.max_number_off_diag == 0 then
      return (0, nil, nil, nil);
    // if count > raw.batch_size then
    //   halt(try! "buffer overflow in BatchedOperator: count=%i, batch_size=%i".format(count, raw.batch_size));

    _computeOffDiag(basisInfo, raw, chunk, alphas, norms, left);
    const totalCount = raw.offsets[count - 1];
    // if totalCount > maxTotalCount then
    //   maxTotalCount = totalCount;
    // computeKeys(totalCount, pointers.betas, pointers.localeIdxs, timers.keysTimer);
    return (totalCount, raw.betas, raw.coeffs, raw.target_indices);
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
