module CommonParameters {
  config const kDisplayTimings : bool = false;
  config const kHashedToBlockNumChunks = 2 * here.maxTaskPar;
  config const kBlockToHashedNumChunks = 2 * numLocales * here.maxTaskPar;
  config const kIsRepresentativeBatchSize : int = 10240;
  config const kEnumerateStatesNumChunks : int = 8 * numLocales * here.maxTaskPar;
  config const kCacheNumberBits : int = 26;
}
