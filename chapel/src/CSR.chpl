module CSR {

private use Timing;

private use CTypes;

private proc mergePathSearch(diagonal : int, a : c_ptrConst(?),
                             sizeA : int, sizeB : int) : (int, int)
{
    var xMin = max(diagonal - sizeB, 0);
    var xMax = min(diagonal, sizeA);

    while xMin < xMax {
        const pivot = (xMin + xMax) / 2;
        if a[pivot] <= diagonal - pivot - 1
          then xMin = pivot + 1;
          else xMax = pivot;
    }

    return (min(xMin, sizeA), diagonal - xMin);
}

record CSR {
  type eltType;
  type idxType;
  var matrixElements : c_ptrConst(eltType);
  var rowOffsets : c_ptrConst(idxType);
  var colIndices : c_ptrConst(idxType);
  var numberRows : int;
  var numberCols : int;
  var numberNonZero : int;
}

proc csrSortIndices(ref matrix : CSR(?eltType, ?idxType)) {

}

proc csrMatvec(const ref matrix : CSR(?matrixEltType, ?idxType),
               xPtr : c_ptrConst(?vectorEltType),
               yPtr : c_ptr(?outEltType),
               numTasks : int = here.maxTaskPar) {
  const _timer = recordTime("csrMatvec");
  // Always use the widest available type for accumulation
  type accType = if isComplex(matrixEltType) || isComplex(vectorEltType) then complex(128) else real(64);

  // Temporary storage for inter-thread fix-up after load-balanced work
  // The last row-id each worked on by each thread when it finished its path segment
  var rowCarry : [0 ..# numTasks] int;
  // The running total within each thread when it finished its path segment
  var valueCarry : [0 ..# numTasks] accType;

  coforall threadIdx in 0 ..# numTasks {
    const numberMergeItems = matrix.numberRows + matrix.numberNonZero;
    const numberItemsPerTask = (numberMergeItems + numTasks - 1) / numTasks;
    const startDiagonal = min(numberItemsPerTask * threadIdx, numberMergeItems);
    const endDiagonal = min(startDiagonal + numberItemsPerTask, numberMergeItems);
    var (x, y) = mergePathSearch(startDiagonal, matrix.rowOffsets + 1,
                                 matrix.numberRows, matrix.numberNonZero);
    const (endX, endY) = mergePathSearch(endDiagonal, matrix.rowOffsets + 1,
                                         matrix.numberRows, matrix.numberNonZero);
    // Consume whole rows
    while x < endX {
      var acc : accType = 0;
      while y < matrix.rowOffsets[x + 1] {
        acc += matrix.matrixElements[y] * xPtr[matrix.colIndices[y]];
        y += 1;
      }
      yPtr[x] = acc:outEltType;
      x += 1;
    }
    // Consume partial portion of thread's last row
    var acc : accType = 0;
    while y < endY {
      acc += matrix.matrixElements[y] * xPtr[matrix.colIndices[y]];
      y += 1;
    }
    // Save carry-outs
    rowCarry[threadIdx] = x;
    valueCarry[threadIdx] = acc;
  }

  // Carry-out fix-up (rows spanning multiple threads)
  for threadIdx in 0 ..# numTasks - 1 do
    if rowCarry[threadIdx] < matrix.numberRows then
      yPtr[rowCarry[threadIdx]] += valueCarry[threadIdx]:outEltType;
}

proc ls_chpl_matrix_vector_product_csr(numberRows : int,
                                       numberCols : int,
                                       numberNonZero : int,
                                       matrixElements : c_ptrConst(?eltType),
                                       rowOffsets : c_ptrConst(?idxType),
                                       colIndices : c_ptrConst(idxType),
                                       x : c_ptrConst(eltType),
                                       y : c_ptr(eltType),
                                       numTasks : int) {
  const matrix = new CSR(eltType=eltType, idxType=idxType,
                         matrixElements=matrixElements,
                         rowOffsets=rowOffsets,
                         colIndices=colIndices,
                         numberRows=numberRows,
                         numberCols=numberCols,
                         numberNonZero=numberNonZero);
  csrMatvec(matrix, x, y, if numTasks <= 0 then here.maxTaskPar else numTasks);
}

export proc ls_chpl_matrix_vector_product_csr_i32_f64(numberRows : int, numberCols : int, numberNonZero : int, matrixElements : c_ptrConst(real(64)), rowOffsets : c_ptrConst(int(32)), colIndices : c_ptrConst(int(32)), x : c_ptrConst(real(64)), y : c_ptr(real(64)), numTasks : int) {
  ls_chpl_matrix_vector_product_csr(numberRows, numberCols, numberNonZero, matrixElements, rowOffsets, colIndices, x, y, numTasks);
}

export proc ls_chpl_matrix_vector_product_csr_i32_c128(numberRows : int, numberCols : int, numberNonZero : int, matrixElements : c_ptrConst(complex(128)), rowOffsets : c_ptrConst(int(32)), colIndices : c_ptrConst(int(32)), x : c_ptrConst(complex(128)), y : c_ptr(complex(128)), numTasks : int) {
  ls_chpl_matrix_vector_product_csr(numberRows, numberCols, numberNonZero, matrixElements, rowOffsets, colIndices, x, y, numTasks);
}

export proc ls_chpl_matrix_vector_product_csr_i64_f64(numberRows : int, numberCols : int, numberNonZero : int, matrixElements : c_ptrConst(real(64)), rowOffsets : c_ptrConst(int(64)), colIndices : c_ptrConst(int(64)), x : c_ptrConst(real(64)), y : c_ptr(real(64)), numTasks : int) {
  ls_chpl_matrix_vector_product_csr(numberRows, numberCols, numberNonZero, matrixElements, rowOffsets, colIndices, x, y, numTasks);
}

export proc ls_chpl_matrix_vector_product_csr_i64_c128(numberRows : int, numberCols : int, numberNonZero : int, matrixElements : c_ptrConst(complex(128)), rowOffsets : c_ptrConst(int(64)), colIndices : c_ptrConst(int(64)), x : c_ptrConst(complex(128)), y : c_ptr(complex(128)), numTasks : int) {
  ls_chpl_matrix_vector_product_csr(numberRows, numberCols, numberNonZero, matrixElements, rowOffsets, colIndices, x, y, numTasks);
}

} // end module CSR
