module MatrixVectorProduct {

proc perLocaleDiagonal(matrix : Operator,
                       const ref x : [] ?eltType,
                       ref y : [] eltType,
                       const ref representatives : [] uint(64),
                       in timings : TimingTree? = nil,
                       numChunks : int = min(matrixVectorDiagonalNumChunks,
                                             representatives.size)) {
  var timer = new stopwatch();
  timer.start();

  const totalSize = representatives.size;
  // const batchSize = (totalSize + numChunks - 1) / numChunks;
  const ranges : [0 ..# numChunks] range(int, boundKind.both, strideKind.one) =
    chunks(0 ..# totalSize, numChunks);
  // var workspace : [0 ..# batchSize] complex(128) = noinit;
  forall r in ranges {
    applyDiagKernel(
      matrix.payload.deref(),
      indices.size,
      c_ptrToConst(representatives[indices.low]),
      c_ptrTo(y[indices.low]),
      c_ptrToConst(x[indices.low]));
    localDiagonalBatch(r, matrix, x, y, representatives);
  }

  timer.stop();
  if timings != nil then
    timings!.addChild("localDiagonal", timer.elapsed());
}

}
