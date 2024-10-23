use CTypes;
use FFI;
private use OS.POSIX;
import IO.FormattedIO.format;
import BitOps.ctz, BitOps.parity;
import RangeChunk.chunks;
import Math;

config const kMatrixVectorNumChunks : int = 1; // here.maxTaskPar;
config const kMatrixVectorInnerChunkSize : int = 1024;

export proc hello() {
  writeln("Hello world!");
}

export proc the_ultimate_solution(alloc_numpy_array_1d : c_fn_ptr, ref result : ls_numpy_array_1d) {
  const size = 12;
  ls_invoke_alloc_numpy_array_1d_callback(alloc_numpy_array_1d, size:c_size_t, result);
  foreach i in 0 ..# size do
    result.data[i] = 42;
}

/* Get the next integer with the same Hamming weight.

   Semantically equivalent to
   ```
   const m = popcount(v);
   v += 1;
   while (popcount(v) != m) { v += 1; }
   return v;
   ```
 */
private inline proc nextStateFixedHamming(v: uint(64)): uint(64) {
  const t = v | (v - 1);
  return (t + 1) | (((~t & (t + 1)) - 1) >> (ctz(v) + 1));
}

export proc ls_enumerate_states_fixed_hamming(numChunks : int(64),
                                              offsets : c_ptrConst(int(64)),
                                              values : c_ptrConst(uint(64)),
                                              dest : c_ptr(uint(64))) {
  forall chunkIdx in 0 ..# numChunks {
    var state = values[chunkIdx];
    for i in offsets[chunkIdx] .. offsets[chunkIdx + 1] - 1 {
      dest[i] = state;
      state = nextStateFixedHamming(state);
    }
  }
}


// proc _applyDiagKernel(diag_terms : c_ptrConst(ls_nonbranching_terms),
//                       const ref alphas : [] uint(64),
//                       const ref xs : [] ?eltType,
//                       ref ys : [] eltType,
//                       param numVectors : int,
//                       param left : bool = false) {
//   const ref terms = diag_terms.deref();
//   const numTerms = terms.number_terms;
//   const vs = terms.v:c_ptrConst(eltType); // TODO: are we sure about this?
//   var acc = new c_array(complex(128), numVectors);
// 
//   foreach (alpha, batchIdx) in zip(alphas, 0 ..) {
//     for param vectorIdx in 0 ..# numVectors do acc[vectorIdx] = 0; // unroll!
// 
//     for termIdx in 0 ..# numTerms {
//       const l_or_r = if left then terms.l[termIdx] else terms.r[termIdx];
//       const delta = (alpha & terms.m[termIdx]) == l_or_r;
//       if delta {
//         const sign = 1 - 2 * parity(alpha & terms.s[termIdx]):eltType;
//         const factor = sign * vs[termIdx]; // TODO: check types!
// 
//         for param vectorIdx in 0 ..# numVectors do acc[vectorIdx] += factor; // unroll!
//       }
//     }
// 
//     for param vectorIdx in 0 ..# numVectors do // unroll!
//       ys[vectorIdx, batchIdx] = acc[vectorIdx]:eltType * xs[vectorIdx, batchIdx];
//   }
// }

// proc applyDiagKernel(diag_terms : c_ptrConst(ls_nonbranching_terms),
//                      const ref alphas : [] uint(64),
//                      const ref xs : [] ?eltType,
//                      ref ys : [] eltType,
//                      param left : bool = false)
//     where alphas.domain.rank == 1 && ys.domain.rank == 2 && xs.domain.rank == 2 {
// 
//   if alphas.dim(0).size != ys.dim(1).size || alphas.dim(0).size != xs.dim(1).size then
//     halt(try! "dimension mismatch: %i vs. %i vs. %i".format(alphas.dim(0).size, ys.dim(1).size, xs.dim(1).size));
// 
//   const numVectors = xs.dim(0).size;
//   if numVectors != ys.dim(0).size then
//     halt(try! "dimension mismatch: %i vs. %i".format(numVectors, ys.dim(0).size));
// 
//   // The diagonal is zero
//   if (diag_terms == nil || diag_terms.deref().number_terms == 0) {
//     ys = 0;
//     return;
//   }
// 
//   select numVectors {
//     when 1 do _applyDiagKernel(diag_terms, alphas, xs, ys, numVectors=1, left=left);
//     when 2 do _applyDiagKernel(diag_terms, alphas, xs, ys, numVectors=2, left=left);
//     when 3 do _applyDiagKernel(diag_terms, alphas, xs, ys, numVectors=3, left=left);
//     when 4 do _applyDiagKernel(diag_terms, alphas, xs, ys, numVectors=4, left=left);
//     when 5 do _applyDiagKernel(diag_terms, alphas, xs, ys, numVectors=5, left=left);
//     when 6 do _applyDiagKernel(diag_terms, alphas, xs, ys, numVectors=6, left=left);
//     when 7 do _applyDiagKernel(diag_terms, alphas, xs, ys, numVectors=7, left=left);
//     when 8 do _applyDiagKernel(diag_terms, alphas, xs, ys, numVectors=8, left=left);
//     otherwise do halt(try! "numVectors=%i not supported yet".format(numVectors));
//   }
// }

// proc applyDiag(const ref matrix : ls_chpl_batched_operator,
//                const ref representatives : [] uint(64),
//                const ref x : [] ?eltType,
//                ref y : [] eltType,
//                numChunks : int = min(kMatrixVectorDiagonalNumChunks, representatives.size)) {
//   // const _timer = recordTime(getRoutineName());
//   const totalSize = representatives.size;
//   const ranges : [0 ..# numChunks] range(int, boundKind.both, strideKind.one) = chunks(0 ..# totalSize, numChunks);
//   const diag_terms = matrix.diag_terms;
// 
//   forall r in ranges do
//     applyDiagKernel(diag_terms, representatives[r], x[.., r], y[.., r]);
// }

// proc ls_matrix_apply_diag(const ref terms : ls_diag_terms,
//                           numStates : int(64),
//                           representativesPtr : c_ptrConst(uint(64)),
//                           numVectors : int(64),
//                           xPtr : c_ptrConst(?eltType),
//                           yPtr : c_ptr(eltType)) {
//   if terms.number_terms == 0 {
//     POSIX.memset(yPtr, 0, numVectors * numStates * c_sizeof(eltType));
//     return;
//   }
// 
//   const numChunks = min(kMatrixVectorDiagonalNumChunks, numStates);
//   const ranges : [0 ..# numChunks] range(int, boundKind.both, strideKind.one) = chunks(0 ..# numStates, numChunks);
//   forall r in ranges do
//     if r.size > 0 {
//       const innerChunkSize = min(r.size, kMatrixVectorDiagonalInnerChunkSize);
//       var coeffsReal = allocate(real(64), innerChunkSize);
//       defer deallocate(coeffsReal);
//       var coeffsImag = allocate(real(64), innerChunkSize);
//       defer deallocate(coeffsImag);
// 
//       var offset = 0;
//       while offset < r.size {
//         const count = min(innerChunkSize, r.size - offset);
//         const totalOffset = r.low + offset;
//         if isComplex(eltType)
//           then ls_invoke_diag_matrix_kernel(c_ptrToConst(terms), count.safeCast(int(32)), representativesPtr + totalOffset, coeffsReal, coeffsImag);
//           else ls_invoke_diag_matrix_kernel(c_ptrToConst(terms), count.safeCast(int(32)), representativesPtr + totalOffset, coeffsReal, nil);
// 
//         foreach k in 0 ..# numVectors do
//           foreach i in 0 ..# count do
//             if isComplex(eltType)
//               then yPtr[k * numStates + totalOffset + i] = (coeffsReal[i] + coeffsImag[i]:imag(64)):eltType * xPtr[k * numStates + totalOffset + i];
//               else yPtr[k * numStates + totalOffset + i] = coeffsReal[i]:eltType * xPtr[k * numStates + totalOffset + i];
//         offset += count;
//       }
//     }
// }
// 
// export proc ls_matrix_apply_diag_f32(const ref terms : ls_diag_terms, numStates : int(64), representativesPtr : c_ptrConst(uint(64)), numVectors : int(64), xPtr : c_ptrConst(real(32)), yPtr : c_ptr(real(32))) { ls_matrix_apply_diag(terms, numStates, representativesPtr, numVectors, xPtr, yPtr); }
// export proc ls_matrix_apply_diag_f64(const ref terms : ls_diag_terms, numStates : int(64), representativesPtr : c_ptrConst(uint(64)), numVectors : int(64), xPtr : c_ptrConst(real(64)), yPtr : c_ptr(real(64))) { ls_matrix_apply_diag(terms, numStates, representativesPtr, numVectors, xPtr, yPtr); }
// export proc ls_matrix_apply_diag_c64(const ref terms : ls_diag_terms, numStates : int(64), representativesPtr : c_ptrConst(uint(64)), numVectors : int(64), xPtr : c_ptrConst(complex(64)), yPtr : c_ptr(complex(64))) { ls_matrix_apply_diag(terms, numStates, representativesPtr, numVectors, xPtr, yPtr); }
// export proc ls_matrix_apply_diag_c128(const ref terms : ls_diag_terms, numStates : int(64), representativesPtr : c_ptrConst(uint(64)), numVectors : int(64), xPtr : c_ptrConst(complex(128)), yPtr : c_ptr(complex(128))) { ls_matrix_apply_diag(terms, numStates, representativesPtr, numVectors, xPtr, yPtr); }
// 

iter innerChunkIter(r : range, innerChunkSize : int)
{
  var offset = 0;
  while offset < r.size {
    const count = min(innerChunkSize, r.size - offset).safeCast(int(32));
    const totalOffset = r.low + offset;
    yield (totalOffset, count);
    offset += count;
  }
}

inline proc safeDeallocate(ptr) { if ptr != nil then deallocate(ptr); }

proc taskLoop(r : range,
              const ref diag : ls_diag_terms,
              const ref off_diag : ls_off_diag_terms,
              const ref state_to_index : ls_xored_state_to_index,
              numStates : int(64),
              representativesPtr : c_ptrConst(uint(64)),
              numVectors : int(64),
              xPtr : c_ptrConst(?eltType),
              yPtr : c_ptr(eltType)) {
  const innerChunkSize = min(r.size, kMatrixVectorInnerChunkSize);

  // A buffer for pre-computed matrix elements
  const bufferSize = 2 * innerChunkSize * (off_diag.number_terms + 1);
  var buffer = allocate(real(64), bufferSize);
  defer safeDeallocate(buffer);

  POSIX.memset(buffer, 0, bufferSize:c_size_t * c_sizeof(real(64)));
  const coeffsReal = buffer;
  const diagCoeffsReal = coeffsReal + innerChunkSize * off_diag.number_terms;
  const coeffsImag = if isComplex(eltType) then diagCoeffsReal + innerChunkSize else nil:c_ptr(real(64));
  const diagCoeffsImag = if isComplex(eltType) then coeffsImag + innerChunkSize * off_diag.number_terms else nil:c_ptr(real(64));

  // A buffer for pre-computed indices
  param isIndexToStateIdentity = true;
  var indexBuffer =
    if isIndexToStateIdentity
      then nil:c_ptr(int(64))
      else allocate(int(64), innerChunkSize * off_diag.number_terms);
  defer safeDeallocate(indexBuffer);
  if indexBuffer != nil then
    POSIX.memset(indexBuffer, 0, (innerChunkSize * off_diag.number_terms):c_size_t * c_sizeof(int(64)));

  // A buffer for 

  type accType = if isComplex(eltType) then complex(128) else real(64);
  for (offset, count) in innerChunkIter(r, innerChunkSize) {
    if diag.number_terms > 0 then
      ls_invoke_diag_matrix_kernel(c_ptrToConst(diag), count, representativesPtr + offset, diagCoeffsReal, diagCoeffsImag);
    if off_diag.number_terms > 0 {
      ls_invoke_off_diag_matrix_kernel(c_ptrToConst(off_diag), count, representativesPtr + offset, coeffsReal, coeffsImag);
      // TODO: Optionally project representatives ^ off_diag.x to new representatives
      // updating coeffsReal and coeffsImag inplace
      // if !requiresProjection then
      //   ls_invoke_xored_state_info_kernel(c_ptrToConst(state_info), count, representativesPtr + offset, off_diag.number_terms, off_diag.x,
      //                                     /*outputs ...*/nil);

      if !isIndexToStateIdentity then
        ls_invoke_xored_state_to_index_kernel(c_ptrToConst(state_to_index), count, representativesPtr + offset,
                                                                            off_diag.number_terms, off_diag.x, indexBuffer);
    }

    foreach k in 0 ..# numVectors do
      foreach i in 0 ..# count {
        var acc : accType =
          xPtr[k * numStates + offset + i]:accType
            * (if isComplex(eltType)
                then diagCoeffsReal[i] + diagCoeffsImag[i]:imag(64)
                else diagCoeffsReal[i]);

        for t in 0 ..# off_diag.number_terms {
          const beta =
            if isIndexToStateIdentity
              then (representativesPtr[offset + i] ^ off_diag.x[t]):int(64)
              else indexBuffer[t * count + i];

          const x = xPtr[k * numStates + beta]:accType;
          if isComplex(eltType) {
            const coeff = coeffsReal[t * count + i] + coeffsImag[t * count + i]:imag(64);
            acc += x * coeff;
          }
          else {
            const coeff = coeffsReal[t * count + i];
            acc = Math.fma(x, coeff, acc);
          }
        }

        yPtr[k * numStates + offset + i] = acc:eltType;
      }
  }
}

proc ls_matrix_apply(const ref diag : ls_diag_terms,
                     const ref off_diag : ls_off_diag_terms,
                     numStates : int(64),
                     representativesPtr : c_ptrConst(uint(64)),
                     numVectors : int(64),
                     xPtr : c_ptrConst(?eltType),
                     yPtr : c_ptr(eltType)) {
  const numChunks = min(kMatrixVectorNumChunks, numStates);
  const ranges : [0 ..# numChunks] range(int, boundKind.both, strideKind.one) = chunks(0 ..# numStates, numChunks);
  const state_to_index = new ls_xored_state_to_index();
  forall r in ranges do
    taskLoop(r, diag, off_diag, state_to_index, numStates, representativesPtr, numVectors, xPtr, yPtr);
}

export proc ls_matrix_apply_f32(const ref diag : ls_diag_terms, const ref off_diag : ls_off_diag_terms, numStates : int(64), representativesPtr : c_ptrConst(uint(64)), numVectors : int(64), xPtr : c_ptrConst(real(32)), yPtr : c_ptr(real(32))) { ls_matrix_apply(diag, off_diag, numStates, representativesPtr, numVectors, xPtr, yPtr); }
export proc ls_matrix_apply_f64(const ref diag : ls_diag_terms, const ref off_diag : ls_off_diag_terms, numStates : int(64), representativesPtr : c_ptrConst(uint(64)), numVectors : int(64), xPtr : c_ptrConst(real(64)), yPtr : c_ptr(real(64))) { ls_matrix_apply(diag, off_diag, numStates, representativesPtr, numVectors, xPtr, yPtr); }
export proc ls_matrix_apply_c64(const ref diag : ls_diag_terms, const ref off_diag : ls_off_diag_terms, numStates : int(64), representativesPtr : c_ptrConst(uint(64)), numVectors : int(64), xPtr : c_ptrConst(complex(64)), yPtr : c_ptr(complex(64))) { ls_matrix_apply(diag, off_diag, numStates, representativesPtr, numVectors, xPtr, yPtr); }
export proc ls_matrix_apply_c128(const ref diag : ls_diag_terms, const ref off_diag : ls_off_diag_terms, numStates : int(64), representativesPtr : c_ptrConst(uint(64)), numVectors : int(64), xPtr : c_ptrConst(complex(128)), yPtr : c_ptr(complex(128))) { ls_matrix_apply(diag, off_diag, numStates, representativesPtr, numVectors, xPtr, yPtr); }

/*
proc ls_matrix_apply_off_diag(const ref terms : ls_off_diag_terms,
                              numStates : int(64),
                              representativesPtr : c_ptrConst(uint(64)),
                              numVectors : int(64),
                              xPtr : c_ptrConst(?eltType),
                              yPtr : c_ptr(eltType)) {
  const D1 = {0 ..# numStates};
  const D2 = {0 ..# numVectors, 0 ..# numStates};
  // NOTE: the casts from c_ptrConst to c_ptr are fine since we save the arrays in
  // const variables afterwards thus regaining const correctness
  const representatives = makeArrayFromPtr(representativesPtr:c_ptr(uint(64)), D1);
  const x = makeArrayFromPtr(xPtr:c_ptr(eltType), D2);
  var y = makeArrayFromPtr(yPtr, D2);

  const numChunks = min(kMatrixVectorDiagonalNumChunks, numStates);
  const ranges : [0 ..# numChunks] range(int, boundKind.both, strideKind.one) = chunks(0 ..# numStates, numChunks);

  forall r in ranges do
    if r.size > 0 {
      const innerChunkSize = min(r.size, kMatrixVectorDiagonalInnerChunkSize);
      var coeffsReal = allocate(real(64), innerChunkSize);
      defer deallocate(coeffsReal);
      var coeffsImag = allocate(real(64), innerChunkSize);
      defer deallocate(coeffsImag);

      var offset = 0;
      while offset < r.size {
        const count = min(innerChunkSize, r.size - offset);
        const totalOffset = r.low + offset;
        if isComplex(eltType)
          then ls_invoke_diag_matrix_kernel(c_ptrToConst(terms), count, representativesPtr + totalOffset, coeffsReal, coeffsImag);
          else ls_invoke_diag_matrix_kernel(c_ptrToConst(terms), count, representativesPtr + totalOffset, coeffsReal, nil);

        foreach k in 0 ..# numVectors do
          foreach i in 0 ..# count do
            if isComplex(eltType)
              then yPtr[k * numStates + totalOffset + i] = (coeffsReal[i] + coeffsImag[i]:imag(64)):eltType * xPtr[k * numStates + totalOffset + i];
              else yPtr[k * numStates + totalOffset + i] = coeffsReal[i]:eltType * xPtr[k * numStates + totalOffset + i];
        offset += count;
      }
    }
}
*/

/*
proc applyOffDiagKernel(const ref matrix : ls_chpl_batched_operator,
                        chunk : range(int),
                        alphas : c_ptrConst(uint(64)),
                        param left : bool = false) {

  const batch_size = chunk.size;
  if batch_size > matrix.batch_size then
    halt(try! "buffer overflow: allocated space for %i elements, but chunk.size=%i".format(matrix.batch_size, batch_size));
  if matrix.betas == nil || matrix.coeffs == nil || matrix.offsets == nil then
    halt("betas, coeffs, and offsets should be pre-allocated");

  const off_diag_terms = matrix.off_diag_terms;
  // Nothing to apply
  if (off_diag_terms == nil || off_diag_terms.deref().number_terms == 0) {
    POSIX.memset(matrix.offsets, 0, batch_size:c_size_t * c_sizeof(int(64)));
    return;
  }

  const ref terms = off_diag_terms.deref();
  const numberTerms = terms.number_terms;
  const vs = terms.v:c_ptrConst(eltType); // TODO: are we sure about this?
  const spinInversionMask = matrix.spin_inversion_mask;
  const spinInversionCharacter = matrix.spin_inversion;

  if batch_size > 0 && alphas == nil then
    halt("applyOffDiagKernel received null alphas");

  var offset = 0;
  for batch_idx in 0 ..# batch_size {
    const oldOffset = offset;
    const alpha = alphas[batch_idx];

    var termIdx = 0;
    do {
      const tX = terms.x[termIdx];
      var acc : complex(128) = 0;

      do {
        const l_or_r = if left then terms.l[termIdx] else terms.r[termIdx];
        const delta = (alpha & terms.m[termIdx]) == l_or_r;
        if delta {
          const sign = 1 - 2 * parity(alpha & terms.s[termIdx]):real;
          const coeff = vs[termIdx] * sign;
          acc += coeff:complex(128);
        }
        termIdx += 1;
      } while termIdx < numberTerms && terms.x[termIdx] == tX;

      if acc != 0 {
        var beta = alpha ^ tX;
        if spinInversionCharacter != 0 {
          const inverted = beta ^ spinInversionMask;
          if inverted < beta {
            beta = inverted;
            acc *= spinInversionCharacter:real(64);
          }
        }

        matrix.coeffs[offset] = acc;
        matrix.betas[offset] = beta;
        if matrix.target_indices != nil then
          matrix.target_indices[offset] = chunk.low + batch_idx;
        offset += 1;
      }
    } while termIdx < numberTerms;

    matrix.offsets[batch_idx] = offset;
  }
}

private proc localProcessExperimental(const ref basis : Basis,
                                      xs : c_ptrConst(?coeffType),
                                      norms : c_ptrConst(uint(16)),
                                      size : int,
                                      basisStates : c_ptrConst(uint(64)),
                                      coeffs : c_ptrConst(complex(128)),
                                      targetIndices : c_ptrConst(uint(64)),
                                      minTargetIndex : uint(64),
                                      numDistinctTargetIndices : int,
                                      targetCoeffs : c_ptr(coeffType),
                                      indicesBuffer : c_ptr(int)) {
  const _timer = recordTime(getRoutineName());

  local {
    // Reset accumulators
    POSIX.memset(targetCoeffs, 0, numDistinctTargetIndices:c_size_t * c_sizeof(coeffType));

    var indices : c_ptr(int(64));
    if numLocales == 1 && basis.info.is_state_index_identity {
      // Special case when we don't have to call ls_hs_state_index
      indices = basisStates:c_ptr(int(64));
    }
    else {
      indices = indicesBuffer;
      const sizeWithPadding = roundUpToMaxBlockSize(size);
      basisStatesToIndices(basis, sizeWithPadding, basisStates, indices);
    }

    var k = 0;
    while k < size {
      var acc : coeffType =
        if indices[k] >=0
          then coeffs[k]:coeffType * xs[indices[k]] * sqrt(norms[indices[k]]:real)
          else 0;
      var targetIndex = targetIndices[k];
      k += 1;
      while k < size && targetIndices[k] == targetIndex {
        if indices[k] >= 0 then
          acc += coeffs[k]:coeffType * xs[indices[k]] * sqrt(norms[indices[k]]:real);
        k += 1;
      }
      targetCoeffs[targetIndex - minTargetIndex] += acc;
    }
  }
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
*/
