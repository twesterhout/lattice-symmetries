module CSR {

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

proc csrMatvec(const ref matrix : CSR(?matrixEltType, ?idxType),
               xPtr : c_ptrConst(?vectorEltType),
               yPtr : c_ptr(?outEltType),
               numTasks : int = here.maxTaskPar) {

  // Temporary storage for inter-thread fix-up after load-balanced work
  // The last row-id each worked on by each thread when it finished its path segment
  var rowCarry : [0 ..# numTasks] int;
  // The running total within each thread when it finished its path segment
  var valueCarry : [0 ..# numTasks] complex(128);

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
      var acc : complex(128) = 0;
      while y < matrix.rowOffsets[x + 1] {
        acc += matrix.matrixElements[y] * xPtr[matrix.colIndices[y]];
        y += 1;
      }
      yPtr[x] = acc:outEltType;
      x += 1;
    }
    // Consume partial portion of thread's last row
    var acc : complex(128) = 0.0;
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
      yPtr[rowCarry[threadIdx]] += valueCarry[threadIdx];
}

} // ena module CSR
