module CommonParameters {
  config const kDisplayTimings : bool = false;
  config const kHashedToBlockNumChunks = 2 * here.maxTaskPar;
  config const kBlockToHashedNumChunks = 2 * numLocales * here.maxTaskPar;
  config const kIsRepresentativeBatchSize : int = 10240;
  config const kEnumerateStatesNumChunks : int = 8 * numLocales * here.maxTaskPar;
  config const kCacheNumberBits : int = 26;

  // MatrixVectorProduct
  config const kMatrixVectorDiagonalNumChunks : int = numLocales * here.maxTaskPar;
  config const kRemoteBufferSize = 10000;
  config const kNumTasks = here.maxTaskPar;

  // OperatorToCsr
  config const kToCsrNumChunks = 100 * here.maxTaskPar;

  // For testing array equality
  config const kAbsTol = 1e-13;
  config const kRelTol = 1e-11;
}
