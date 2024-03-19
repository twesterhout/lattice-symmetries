module MatrixVectorProduct {

private use BatchedOperator;
private use CSR;
private use CommonParameters;
private use ConcurrentQueue;
private use FFI;
private use ForeignTypes;
private use Timing;
private use Utils;
private use Vector;

private import Reflection.getRoutineName;
private import Sort;
private use AllLocalesBarriers;
private use CTypes;
private use ChapelLocks;
private use DynamicIters;
private use IO.FormattedIO;
private use OS.POSIX;
private use RangeChunk;
private use Time;

proc perLocaleDiagonal(const ref matrix : Operator,
                       const ref x : [] ?eltType,
                       ref y : [] eltType,
                       const ref representatives : [] uint(64),
                       numChunks : int = min(kMatrixVectorDiagonalNumChunks, representatives.size)) {
  assert(matrix.locale == here);
  const _timer = recordTime(getRoutineName());

  const totalSize = representatives.size;
  const ranges : [0 ..# numChunks] range(int, boundKind.both, strideKind.one) = chunks(0 ..# totalSize, numChunks);
  const diag_terms = matrix.payload.deref().diag_terms;
  forall r in ranges do
    applyDiagKernel(
      diag_terms,
      r.size,
      c_ptrToConst(representatives[r.low]),
      c_ptrTo(y[r.low]),
      c_ptrToConst(x[r.low]));
}

record WorkerBuffersView {
  type coeffType;
  var localeIdx : int;
  var capacity : int;
  var basisStates : c_ptr(uint(64));
  var coeffs : c_ptr(coeffType);

  proc init(type coeffType) {
    this.coeffType = coeffType;
    this.localeIdx = -1;
    this.capacity = 0;
    this.basisStates = nil;
    this.coeffs = nil;
  }
  proc init(type coeffType, capacity : int) {
    this.coeffType = coeffType;
    this.localeIdx = here.id;
    this.capacity = capacity;
    this.basisStates = nil;
    this.coeffs = nil;
  }
  proc init(const ref other : WorkerBuffersView(?coeffType)) {
    this.coeffType = other.coeffType;
    this.localeIdx = other.localeIdx;
    this.capacity = other.capacity;
    this.basisStates = other.basisStates;
    this.coeffs = other.coeffs;
  }
}

record WorkerBuffers {
  type coeffType;
  var raw : WorkerBuffersView(coeffType);

  proc init(type coeffType) {
    this.coeffType = coeffType;
    this.raw = new WorkerBuffersView(coeffType);
  }
  proc init(raw : WorkerBuffersView(?coeffType)) {
    this.coeffType = coeffType;
    this.raw = raw;
    init this;

    if raw.capacity > 0 {
      assert(this.raw.localeIdx == here.id);
      this.raw.basisStates = allocate(uint(64), raw.capacity);
      this.raw.coeffs = allocate(coeffType, raw.capacity);
    }
  }

  proc deinit() {
    if raw.capacity > 0 {
      deallocate(raw.basisStates);
      deallocate(raw.coeffs);
    }
  }

  forwarding raw;
}

inline proc encodeAsWord(x : int) { return x:uint(64); }
inline proc decodeFromWord(w : uint(64), ref x : int) { x = w:int; }


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

record BufferFilledMessage {
  var srcLocaleIdx : int(16);
  var srcTaskIdx : int(16);
  var size : int(32);
}

proc encodeAsWord(x : BufferFilledMessage) : uint(64) {
  return (x.srcLocaleIdx:uint(64) << (32 + 16))
       | (x.srcTaskIdx:uint(64) << 32)
       | x.size:uint(64);
}
proc decodeFromWord(w : uint(64), ref x : BufferFilledMessage) {
  x = new BufferFilledMessage(
    srcLocaleIdx=(w >> (32 + 16)):int(16),
    srcTaskIdx=(w >> 32):int(16),
    size=w:int(32));
}

record RemoteLocaleStateView {
  var localBuffersPtr : c_ptr(void);
  var localQueuePtr : c_ptr(void);
}

var globalLocaleStateStore : [0 ..# numLocales] RemoteLocaleStateView;


record LocaleState {
  type coeffType;

  // Contains range(int) ranges of 'representatives' that we should process.
  var _chunksDom : domain(1);
  var _chunks : [_chunksDom] range(int);
  // We push indices of chunks onto a queue such that workers can race to pop from it.
  // var chunkIndicesQueue : ConcurrentQueueWrapper(int);
  var chunkIndex : chpl__processorAtomicType(int);

  // Allocated buffers to which tasks from other locales write their data
  var _localBuffersDom : domain(2);
  var localBuffers : [_localBuffersDom] WorkerBuffers(coeffType);

  // Remote buffers to which we will send data
  var _remoteBuffersDom : domain(2);
  var remoteBuffers : [_remoteBuffersDom] WorkerBuffersView(coeffType);

  // Contains messages from remote locales that tell us that we should run
  // localProcess on data in localBuffers
  var localQueue : ConcurrentQueueWrapper(BufferFilledMessage);
  var numberCompleted : atomic(int);

  var _remoteQueuesDom : domain(1);
  var remoteQueues : [_remoteQueuesDom] ConcurrentQueueWrapper(BufferFilledMessage);

  var matrix : Operator;
  var basis : Basis;
  var localeIdxFn;
  var numberStates : int;
  var xPtr : c_ptrConst(coeffType);
  var yPtr : c_ptr(coeffType);
  var representativesPtr : c_ptrConst(uint(64));
  var normsPtr : c_ptrConst(uint(16));

  proc init(const ref matrix : Operator,
            const ref representatives : [] uint(64),
            const ref norms : [] uint(16),
            xs : c_ptrConst(?coeffType),
            ys : c_ptr(coeffType),
            numChunks : int,
            numTasks : int,
            remoteBufferSize : int) {
    const _timer = recordTime("LocaleState_init");
    this.coeffType = coeffType;

    this._chunksDom = {0 ..# numChunks};
    this._chunks = chunks(0 ..# representatives.size, numChunks);

    const D1 = if numLocales > 1 then {0 ..# numLocales, 0 ..# numTasks} else {0 .. -1, 0 .. -1};
    this._localBuffersDom = D1;
    this.localBuffers = [_i in D1] new WorkerBuffers(new WorkerBuffersView(coeffType, remoteBufferSize));

    const D2 = if numLocales > 1 then {0 ..# numTasks, 0 ..# numLocales} else {0 .. -1, 0 .. -1};
    this._remoteBuffersDom = D2;
    this.remoteBuffers = [_i in D2] new WorkerBuffersView(coeffType);

    this.localQueue = new ConcurrentQueueWrapper(capacity=if numLocales > 1 then numLocales * numTasks else 1,
                                                 nullElement=new BufferFilledMessage(-1, -1, 0));

    const D3 = if numLocales > 1 then {0 ..# numLocales} else {0 .. -1};
    this._remoteQueuesDom = D3;
    this.remoteQueues = [i in D3] new ConcurrentQueueWrapper(BufferFilledMessage, nil:c_ptr(ConcurrentQueue), i, false);

    this.matrix = matrix;
    this.basis = new Basis(this.matrix.payload.deref().basis, owning=false);
    this.localeIdxFn = this.basis.getLocaleIdxFn();
    this.numberStates = representatives.size;
    this.xPtr = xs;
    this.yPtr = ys;
    this.representativesPtr = safe_c_ptrToConst(representatives);
    this.normsPtr = safe_c_ptrToConst(norms);
    init this;

    chunkIndex.write(0);
    numberCompleted.write(0);

    // pre-compile Basis kernels
    if numberStates > 0 && (!basis.info.is_state_index_identity || numLocales != 1) {
      const _stateToIndexKernel = ls_chpl_get_state_to_index_kernel(basis.payload);
    }
  }
  proc init(const ref matrix : Operator,
            const ref representatives : [] uint(64),
            const ref norms : [] uint(16),
            const ref xs : [] ?coeffType,
            ref ys : [] coeffType,
            numChunks : int,
            numTasks : int,
            remoteBufferSize : int) {
    this.init(matrix, representatives, norms, safe_c_ptrToConst(xs), safe_c_ptrTo(ys), numChunks, numTasks, remoteBufferSize);
  }
}

proc ref LocaleState.getRemoteView() {
  return new RemoteLocaleStateView(
    localBuffersPtr=c_ptrTo(localBuffers[0, 0]):c_ptr(void),
    localQueuePtr=c_ptrTo(localQueue.queuePtr)
  );
}

proc ref LocaleState.completeInitialization(const ref globalStore) {
  const _timer = recordTime(getRoutineName());
  if numLocales > 1 {
    const numTasks = remoteBuffers.dim(0).size;
    for localeIdx in 0 ..# numLocales {
      const view : RemoteLocaleStateView = globalStore[localeIdx];

      for taskIdx in 0 ..# numTasks {
        var remoteBufferView : WorkerBuffersView(coeffType);
        const rvfPtr = view.localBuffersPtr:c_ptrConst(WorkerBuffers(coeffType)) + localeIdx * numTasks + taskIdx;
        on Locales[localeIdx] do
          remoteBufferView = rvfPtr.deref().raw;

        remoteBuffers[taskIdx, localeIdx] = remoteBufferView;
      }

      ref queue = remoteQueues[localeIdx];
      assert(queue.queuePtr == nil);
      assert(queue.owning == false);
      assert(queue.localeId == localeIdx);
      queue.queuePtr = view.localQueuePtr:c_ptr(ConcurrentQueue);
    }
  }
}

proc ref LocaleState.getNextChunkIndex() {
  const i = chunkIndex.fetchAdd(1);
  return if i < _chunksDom.size then i else -1;
}

proc LocaleState.allCompleted() : bool {
  const numTasks = localBuffers.dim(1).size;
  return numberCompleted.read() == (numLocales - 1) * numTasks;
}

proc ref LocaleState.processReceivedData(msg : BufferFilledMessage) {
  assert(false);
  if msg.size == -1 {
    // srcTaskIdx on srcLocaleIdx indicates to us that we should not expect any more data from it
    // isCompleted[msg.srcLocaleIdx, msg.srcTaskIdx].write(true);
    numberCompleted.add(1);
  }
  else {
    assert(false);
    ref buffer = localBuffers[msg.srcLocaleIdx, msg.srcTaskIdx];
    assert(msg.size <= buffer.capacity);
    // localProcess(basis, accessor,
    //              buffer.basisStates,
    //              buffer.coeffs, // :c_ptrConst(buffer.coeffs.eltType),
    //              msg.size);

    // Send a message back to srcTaskIdx on srcLocaleIdx saying that the buffer is free again
    // ref queue = remoteQueues[msg.srcLocaleIdx, msg.srcTaskIdx];
    // queue.push(new BufferEmptiedMessage(localeIdx=msg.srcLocaleIdx, taskIdx=msg.srcTaskIdx));
  }
}
proc ref LocaleState.tryProcessReceivedData(maxCount : int = numLocales) {
  var count : int = 0;
  do {
    var msg : BufferFilledMessage;
    const success = localQueue.tryPop(msg);
    if success {
      processReceivedData(msg);
      count += 1;
    }
  } while success && count < maxCount;
  // logDebug("tryProcessReceivedData: count=", count);
}

record Worker {
  var batchedOperator : BatchedOperator;
  var localeStatePtr;
  var taskIdx : int;
  var indices;
  var accumulators;

  proc init(taskIdx : int, ref localeState, maxChunkSize : int) {
    const _timer = recordTime("Worker_init");
    this.batchedOperator = new BatchedOperator(localeState.matrix.payload, maxChunkSize);
    this.localeStatePtr = c_addrOf(localeState);
    this.taskIdx = taskIdx;
    this.indices = allocate(int,
      roundUpToMaxBlockSize(maxChunkSize * (1 + max(1, localeState.matrix.max_number_off_diag))));
    this.accumulators = allocate(localeState.coeffType, maxChunkSize);
  }

  proc deinit() {
    deallocate(indices);
    deallocate(accumulators);
    // logDebug("maxTotalCount=", batchedOperator.maxTotalCount);
  }

  proc ref run() {
    const _timer = recordTime(getRoutineName());
    ref localeState = localeStatePtr.deref();
    const ref basis = localeState.basis;
    const ref localeIdxFn = localeState.localeIdxFn;

    while true {
      const chunkIdx = localeState.getNextChunkIndex();
      if chunkIdx == -1 { break; }

      const chunk : range(int) = localeState._chunks[chunkIdx];
      const (n, basisStatesPtr, coeffsPtr, targetIndicesPtr) =
        batchedOperator.computeOffDiag(
          chunk,
          localeState.representativesPtr + chunk.low,
          localeState.normsPtr + chunk.low,
          left=false);

      if numLocales > 1 {
        assert(false, "not implemented");
        // localeIdxFn(n, basisStatesPtr, batchedOperator.raw.locale_indices);

        // var radixOffsets : c_array(int, 257);
        // var submitted : c_array(bool, 256);
        // stableRadixOneStep(n, batchedOperator.raw.locale_indices, radixOffsets, basisStatesPtr, coeffsPtr);

        // for
      }
      else if n > 0 {
        // Process all the generated data locally
        localProcessExperimental(localeState.basis,
                                 localeState.xPtr,
                                 localeState.normsPtr,
                                 n,
                                 basisStatesPtr,
                                 coeffsPtr:c_ptrConst(complex(128)),
                                 targetIndicesPtr,
                                 chunk.low,
                                 chunk.size,
                                 accumulators,
                                 indices);

        foreach k in 0 ..# chunk.size {
          localeState.yPtr[chunk.low + k] += accumulators[k];
        }
      }
    } // end while true

    if numLocales > 1 {
      // no more chunks to produce, but we could still have data to analyze
      // logDebug("1 tryProcessReceivedData");
      localeState.tryProcessReceivedData();

      // tell other locales that we're done producing data
      for i in 0 ..# numLocales {
        // logDebug("queue.push");
        ref queue = localeState.remoteQueues[i];
        queue.push(new BufferFilledMessage(here.id.safeCast(int(16)), taskIdx.safeCast(int(16)), -1:int(32)));
      }

      while !localeState.allCompleted() {
        // logDebug("2 tryProcessReceivedData");
        localeState.tryProcessReceivedData();
      }
    }
  }
}

proc perLocaleOffDiagonal(matrix : Operator,
                          const ref xs : [] ?eltType,
                          ref ys : [] eltType,
                          const ref representatives : [] uint(64),
                          const ref norms : [] uint(16)) {
  const _timer = recordTime(getRoutineName());
  const numTasks = kNumTasks;
  const remoteBufferSize = max(kRemoteBufferSize, matrix.max_number_off_diag_estimate);
  const numChunks =
    min((representatives.size * matrix.max_number_off_diag_estimate + remoteBufferSize - 1) / remoteBufferSize,
        representatives.size);

  var localeState = new LocaleState(matrix, representatives, norms, xs, ys, numChunks, numTasks, remoteBufferSize);
  if numLocales > 1 {
    globalLocaleStateStore[here.id] = localeState.getRemoteView();
    allLocalesBarrier.barrier();
    localeState.completeInitialization(globalLocaleStateStore);
    // TODO: do we need this?
    allLocalesBarrier.barrier();
  }

  const maxChunkSize = max reduce [c in localeState._chunks] c.size;
  coforall taskIdx in 0 ..# numTasks with (ref localeState) {
    var worker = new Worker(taskIdx, localeState, maxChunkSize);
    worker.run();
  }

  // TODO: do we need this?
  allLocalesBarrier.barrier();
}

proc extractDiag(matrix : Operator,
                 type eltType,
                 const ref representatives : [] uint(64),
                 const ref norms : [] uint(16)) {
  const _timer = recordTime(getRoutineName());
  assert(matrix.locale == here);
  assert(representatives.locale == here);
  assert(norms.locale == here);

  var diag : [0 ..# representatives.size] eltType;

  const numChunks = min(kNumTasks, representatives.size);
  const ranges = chunks(0 ..# representatives.size, numChunks);
  forall chunk in ranges {
    applyDiagKernel(
      matrix.payload.deref().diag_terms,
      chunk.size,
      c_ptrToConst(representatives[chunk.low]),
      c_ptrTo(diag[chunk.low]),
      nil);
  }

  return diag;
}

proc convertOffDiagToCsr(matrix : Operator,
                         type eltType,
                         const ref representatives : [] uint(64),
                         const ref norms : [] uint(16)) {
  const _timer = recordTime(getRoutineName());
  assert(matrix.locale == here);
  assert(representatives.locale == here);
  assert(norms.locale == here);

  const ref basisInfo = matrix.basis.info;

  const numTasks = kNumTasks;
  const numChunks = min(representatives.size, kToCsrNumChunks);
  var localeState = new LocaleState(matrix, representatives, norms,
                                    nil:c_ptrConst(complex(128)), nil:c_ptr(complex(128)),
                                    numChunks, numTasks, remoteBufferSize=0);
  const maxChunkSize = max reduce [c in localeState._chunks] c.size;
  const estimatedNumberTerms = matrix.max_number_off_diag_estimate;

  const chunkDests : [localeState._chunks.domain] range(int);
  var capacity = 0;
  for (d, c) in zip(chunkDests, localeState._chunks) {
    const n = roundUpToMaxBlockSize(c.size * estimatedNumberTerms);
    d = capacity ..# n;
    capacity += n;
  }
  assert(capacity >= representatives.size * estimatedNumberTerms);

  var matrixElementsPtr = allocate(eltType, capacity);
  POSIX.memset(matrixElementsPtr, 0, capacity:c_size_t * c_sizeof(eltType));
  var rowOffsetsPtr = allocate(int(64), representatives.size + 1);
  POSIX.memset(rowOffsetsPtr, 0, (representatives.size + 1):c_size_t * c_sizeof(int(64)));
  var colIndicesPtr = allocate(int(64), capacity);
  POSIX.memset(colIndicesPtr, 0, capacity:c_size_t * c_sizeof(int(64)));
  var numberNonZero : atomic int;

  coforall taskIdx in 0 ..# numTasks with (ref localeState) {

    var batchedOperator = new BatchedOperator(matrix.payload, maxChunkSize);
    var buffer : [0 ..# estimatedNumberTerms] (complex(128), int);

    while true {
      const chunkIdx = localeState.getNextChunkIndex();
      if chunkIdx == -1 { break; }
      const chunk : range(int) = localeState._chunks[chunkIdx];
      const dest : range(int) = chunkDests[chunkIdx];

      const copy = batchedOperator.raw;
      if basisInfo.is_state_index_identity then
        batchedOperator.raw.betas = (colIndicesPtr + dest.low):c_ptr(uint(64));
      batchedOperator.raw.coeffs = matrixElementsPtr + dest.low;
      batchedOperator.raw.offsets = (rowOffsetsPtr + 1) + chunk.low;

      const (totalCount, _, _, _) =
        batchedOperator.computeOffDiag(
          chunk,
          localeState.representativesPtr + chunk.low,
          localeState.normsPtr + chunk.low,
          left=false);

      if !basisInfo.is_state_index_identity && totalCount > 0 {
        const countWithPadding = roundUpToMaxBlockSize(totalCount);
        basisStatesToIndices(matrix.basis, countWithPadding, batchedOperator.raw.betas, colIndicesPtr + dest.low);
      }

      for k in 0 ..# totalCount {
        const i = colIndicesPtr[dest.low + k];
        batchedOperator.raw.coeffs[k] *= sqrt(localeState.normsPtr[i]:real);
      }

      // Sort elements
      for k in 0 ..# chunk.size {
        const lo = if k == 0 then 0 else batchedOperator.raw.offsets[k - 1];
        const hi = batchedOperator.raw.offsets[k];
        const coeffsPtr = batchedOperator.raw.coeffs + lo;
        const indicesPtr = colIndicesPtr + dest.low + lo;

        for i in 0 ..# (hi - lo) do
          buffer[i] = (coeffsPtr[i], indicesPtr[i]);
        record Comparator { inline proc compare(a, b) { return (new Sort.DefaultComparator()).compare(a[1], b[1]); } }
        Sort.QuickSort.quickSort(buffer[0 ..# (hi - lo)], comparator=new Comparator());
        foreach i in 0 ..# (hi - lo) {
          coeffsPtr[i] = buffer[i][0];
          indicesPtr[i] = buffer[i][1];
        }
      }

      numberNonZero.add(totalCount);

      // NOTE: important!
      batchedOperator.raw = copy;
    } // end while true

  }

  var offset = 0;
  for chunkIdx in 0 ..# numChunks {
    const chunk = localeState._chunks[chunkIdx];
    const dest = chunkDests[chunkIdx];
    const chunkSize = rowOffsetsPtr[1 + chunk.high];
    // writeln(chunk, ", ", chunkSize, ", ", offset);
    foreach k in chunk {
      rowOffsetsPtr[1 + k] += offset;
    }
    assert(offset <= dest.low);
    if offset < dest.low {
      POSIX.memmove(matrixElementsPtr + offset, matrixElementsPtr + dest.low,
                    chunkSize:c_size_t * c_sizeof(eltType));
      POSIX.memmove(colIndicesPtr + offset, colIndicesPtr + dest.low,
                    chunkSize:c_size_t * c_sizeof(int(64)));
    }
    offset += chunkSize;
  }

  return new CSR(complex(128), int(64), matrixElementsPtr, rowOffsetsPtr, colIndicesPtr,
                 representatives.size, representatives.size, numberNonZero.read());
}

proc perLocaleMatrixVector(matrix : Operator,
                           const ref x : [] ?eltType,
                           ref y : [] eltType,
                           const ref representatives : [] uint(64),
                           const ref norms : [] uint(16)) {
  const _timer = recordTime(getRoutineName());
  assert(matrix.locale == here);
  assert(x.locale == here);
  assert(y.locale == here);
  assert(representatives.locale == here);
  perLocaleDiagonal(matrix, x, y, representatives);
  if matrix.payload.deref().max_number_off_diag > 0 {
    perLocaleOffDiagonal(matrix, x, y, representatives, norms);
  }
}

proc ls_chpl_matrix_vector_product(matrixPtr : c_ptrConst(ls_hs_operator),
                                   numVectors : c_int,
                                   xPtr : c_ptrConst(?eltType),
                                   yPtr : c_ptr(eltType)) {
  if matrixPtr == nil || xPtr == nil || yPtr == nil then
    halt("matrixPtr, xPtr, and yPtr must not be nil");

  var matrix = new Operator(matrixPtr, owning=false);
  const ref basis = matrix.basis;
  const numberStates = basis.payload.deref().local_representatives.num_elts;
  const basisStatesPtr = basis.payload.deref().local_representatives.elts:c_ptrConst(uint(64));
  const normsPtr = basis.payload.deref().local_norms.elts:c_ptrConst(uint(16));
  if !basis.payload.deref().is_built then halt("basis states have not been built");

  // NOTE: the cast from c_ptrConst to c_ptr is fine since we save the array in
  // a const variable afterwards thus regaining const correctness
  const representatives = makeArrayFromPtr(basisStatesPtr:c_ptr(uint(64)), {0 ..# numberStates});
  const norms = makeArrayFromPtr(normsPtr:c_ptr(uint(16)), {0 ..# numberStates});
  for k in 0 ..# (numVectors:int) {
    // NOTE: same thing applies to the cast here
    const x = makeArrayFromPtr(xPtr:c_ptr(eltType) + k * numberStates, {0 ..# numberStates});
    var y = makeArrayFromPtr(yPtr + k * numberStates, {0 ..# numberStates});
    perLocaleMatrixVector(matrix, x, y, representatives, norms);
  }
}

export proc ls_chpl_matrix_vector_product_f64(matrixPtr : c_ptrConst(ls_hs_operator),
                                              numVectors : c_int,
                                              xPtr : c_ptrConst(real(64)),
                                              yPtr : c_ptr(real(64))) {
  ls_chpl_matrix_vector_product(matrixPtr, numVectors, xPtr, yPtr);
}

export proc ls_chpl_matrix_vector_product_c128(matrixPtr : c_ptrConst(ls_hs_operator),
                                               numVectors : c_int,
                                               xPtr : c_ptrConst(complex(128)),
                                               yPtr : c_ptr(complex(128))) {
  ls_chpl_matrix_vector_product(matrixPtr, numVectors, xPtr, yPtr);
}

export proc ls_chpl_off_diag_to_csr_c128(matrixPtr : c_ptrConst(ls_hs_operator),
                                         matrixElements : c_ptr(chpl_external_array),
                                         rowOffsets : c_ptr(chpl_external_array),
                                         colIndices : c_ptr(chpl_external_array)) {

  const matrix = new Operator(matrixPtr, owning=false);
  const ref basis = matrix.basis;
  const numberStates = basis.payload.deref().local_representatives.num_elts;
  const basisStatesPtr = basis.payload.deref().local_representatives.elts:c_ptrConst(uint(64));
  const normsPtr = basis.payload.deref().local_norms.elts:c_ptrConst(uint(16));
  if !basis.payload.deref().is_built then halt("basis states have not been built");

  const basisStates = makeArrayFromPtr(basisStatesPtr:c_ptr(uint(64)), {0 ..# numberStates});
  const norms = makeArrayFromPtr(normsPtr:c_ptr(uint(16)), {0 ..# numberStates});
  var csr = convertOffDiagToCsr(matrix, complex(128), basisStates, norms);

  matrixElements.deref() = chpl_make_external_array_ptr_free(csr.matrixElements:c_ptr(void), csr.numberNonZero);
  rowOffsets.deref() = chpl_make_external_array_ptr_free(csr.rowOffsets:c_ptr(void), csr.numberRows + 1);
  colIndices.deref() = chpl_make_external_array_ptr_free(csr.colIndices:c_ptr(void), csr.numberNonZero);
}

export proc ls_chpl_extract_diag_c128(matrixPtr : c_ptrConst(ls_hs_operator),
                                      diag : c_ptr(chpl_external_array)) {
  const matrix = new Operator(matrixPtr, owning=false);
  const ref basis = matrix.basis;
  const numberStates = basis.payload.deref().local_representatives.num_elts;
  const basisStatesPtr = basis.payload.deref().local_representatives.elts:c_ptrConst(uint(64));
  const normsPtr = basis.payload.deref().local_norms.elts:c_ptrConst(uint(16));
  if !basis.payload.deref().is_built then halt("basis states have not been built");

  const basisStates = makeArrayFromPtr(basisStatesPtr:c_ptr(uint(64)), {0 ..# numberStates});
  const norms = makeArrayFromPtr(normsPtr:c_ptr(uint(16)), {0 ..# numberStates});

  var arr = extractDiag(matrix, complex(128), basisStates, norms);
  diag.deref() =
    if arr.size > 0 then convertToExternalArray(arr)
                    else chpl_make_external_array_ptr(nil:c_ptr(void), 0);
}

}
