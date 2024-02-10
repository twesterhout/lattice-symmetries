module MatrixVectorProduct {

private use BatchedOperator;
private use CommonParameters;
private use ConcurrentAccessor;
private use ConcurrentQueue;
private use FFI;
private use ForeignTypes;
private use Timing;
private use Vector;

private import Reflection.getRoutineName;
private use AllLocalesBarriers;
private use CTypes;
private use ChapelLocks;
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
  var _timer = recordTime("perLocaleDiagonal");

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


// var localProcessCounts : Vector(int);
// localProcessCounts.reserve(7000);
// var localProcessCountsLock : chpl_LocalSpinlock;

private proc localProcessExperimental(const ref basis : Basis,
                                      xs : c_ptrConst(?coeffType),
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
      const _timer = recordTime("allocate and basisStatesToIndices");
      indices = indicesBuffer;
      basisStatesToIndices(basis, size, basisStates, indices);
    }

    foreach k in 0 ..# size {
      const x = xs[indices[k]:int];
      if x != 0 then
        targetCoeffs[targetIndices[k] - minTargetIndex] += coeffs[k]:coeffType * x;
    }
  }
}

private proc localProcess(const ref basis : Basis,
                          ref accessor : ConcurrentAccessor(?coeffType),
                          basisStates : c_ptrConst(uint(64)),
                          coeffs : c_ptr(?t),
                          size : int) {
  const _timer = recordTime("localProcess");

  // localProcessCountsLock.lock();
  // localProcessCounts.pushBack(size);
  // localProcessCountsLock.unlock();

  local {
    // count == 0 has to be handled separately because c_ptrTo(indices) fails
    // when the size of indices is 0.
    if size == 0 then return;

    // Special case when we don't have to call ls_hs_state_index
    if numLocales == 1 && basis.info.is_state_index_identity {
      foreach k in 0 ..# size {
        const i = basisStates[k]:int;
        const c = coeffs[k]:coeffType;
        accessor.localAdd(i, c);
      }
    }
    else {
      var indices = allocate(int, size);
      defer deallocate(indices);
      basisStatesToIndices(basis, size, basisStates, indices);

      /*
      const _t = recordTime("localProcess.unstableRadixOneStep");
      var keys = allocate(uint(8), size);
      defer deallocate(keys);
      var offsets : c_array(int, 257);

      foreach k in 0 ..# size do
        keys[k] = (indices[k] >> 12):uint(8);

      unstableRadixOneStep(size, keys, offsets, indices, coeffs);

      const dataPtr = accessor._data:c_ptr(void):c_ptr(coeffType);
      for j in 0 ..# 256 {
        if offsets[j + 1] - offsets[j] > 0 {
          ref lock = accessor._locks[j];
          while lock.l.read() || lock.l.testAndSet(memoryOrder.acquire) {
            // nothing
          }

          foreach k in offsets[j] .. offsets[j + 1] - 1 {
            const i = indices[k];
            const c = coeffs[k]:coeffType;
            dataPtr[i] += c;
          }

          lock.l.clear(memoryOrder.release);
        }
      }
      */

      {
        const _t = recordTime("localProcess.localAdd");
        foreach k in 0 ..# size {
          const i = indices[k];
          const c = coeffs[k]:coeffType;
          if c != 0 {
            // Importantly, the user could have made a mistake and given us an
            // operator which does not respect the basis symmetries. Then we could
            // have that a |σ⟩ was generated that doesn't belong to our basis. In
            // this case, we should throw an error.
            if i >= 0 then accessor.localAdd(i, c);
                      else halt("invalid index: " + i:string +
                                " for state " + basisStates[k]:string +
                                " with coeff " + c:string);
          }
        }
      }
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
  var accessor : ConcurrentAccessor(coeffType);
  var numberStates : int;
  var xPtr : c_ptrConst(coeffType);
  var representativesPtr : c_ptrConst(uint(64));

  proc init(const ref matrix : Operator,
            const ref representatives : [] uint(64),
            const ref xs : [] ?coeffType,
            ref ys : [] coeffType,
            numChunks : int,
            numTasks : int,
            remoteBufferSize : int) {
    this.coeffType = coeffType;

    this._chunksDom = {0 ..# numChunks};
    this._chunks = chunks(0 ..# representatives.size, numChunks);
    // this.chunkIndicesQueue = new ConcurrentQueueWrapper(numChunks, nullElement=-1);

    const D1 = {0 ..# numLocales, 0 ..# numTasks};
    this._localBuffersDom = D1;
    this.localBuffers = [_i in D1] new WorkerBuffers(new WorkerBuffersView(coeffType, remoteBufferSize));

    const D2 = {0 ..# numTasks, 0 ..# numLocales};
    this._remoteBuffersDom = D2;
    this.remoteBuffers = [_i in D2] new WorkerBuffersView(coeffType);

    this.localQueue = new ConcurrentQueueWrapper(capacity=numLocales * numTasks,
                                                 nullElement=new BufferFilledMessage(-1, -1, 0));

    const D3 = {0 ..# numLocales};
    this._remoteQueuesDom = D3;
    this.remoteQueues = [i in D3] new ConcurrentQueueWrapper(BufferFilledMessage, nil:c_ptr(ConcurrentQueue), i, false);

    this.matrix = matrix;
    this.basis = new Basis(this.matrix.payload.deref().basis, owning=false);
    this.localeIdxFn = this.basis.getLocaleIdxFn();
    this.accessor = new ConcurrentAccessor(ys);
    this.numberStates = representatives.size;
    this.xPtr = c_ptrToConst(xs);
    this.representativesPtr = c_ptrToConst(representatives);
    init this;

    // for i in _chunksDom {
    //   const success = chunkIndicesQueue.tryPush(i);
    //   if !success then
    //     halt("failed to push to chunkIndicesQueue");
    // }
    chunkIndex.write(0);
    numberCompleted.write(0);

    // pre-compile Basis kernels
    if !basis.info.is_state_index_identity || numLocales != 1 {
      const _stateToIndexKernel = ls_chpl_get_state_to_index_kernel(basis.payload);
    }
  }
}

proc ref LocaleState.getRemoteView() {
  return new RemoteLocaleStateView(
    localBuffersPtr=c_ptrTo(localBuffers[0, 0]):c_ptr(void),
    localQueuePtr=c_ptrTo(localQueue.queuePtr)
  );
}

proc ref LocaleState.completeInitialization(const ref globalStore) {
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

proc LocaleState.getNextChunkIndex() {
  // var i : int;
  // const success = chunkIndicesQueue.tryPop(i);
  // return if success then i else -1;
  const i = chunkIndex.fetchAdd(1);
  return if i < _chunksDom.size then i else -1;
}

proc LocaleState.allCompleted() : bool {
  // logDebug("allCompleted: numberCompleted.read()=", numberCompleted.read());
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
    ref buffer = localBuffers[msg.srcLocaleIdx, msg.srcTaskIdx];
    assert(msg.size <= buffer.capacity);
    localProcess(basis, accessor,
                 buffer.basisStates,
                 buffer.coeffs, // :c_ptrConst(buffer.coeffs.eltType),
                 msg.size);

    // Send a message back to srcTaskIdx on srcLocaleIdx saying that the buffer is free again
    assert(false);
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
    this.batchedOperator = new BatchedOperator(localeState.matrix.payload, maxChunkSize);
    this.localeStatePtr = c_addrOf(localeState);
    this.taskIdx = taskIdx;
    this.indices = allocate(int, maxChunkSize * (1 + max(1, localeState.matrix.max_number_off_diag)));
    this.accumulators = allocate(localeState.coeffType, maxChunkSize);
  }

  proc deinit() {
    deallocate(indices);
    deallocate(accumulators);
  }

  proc ref run() {
    ref localeState = localeStatePtr.deref();
    const ref basis = localeState.basis;
    const ref localeIdxFn = localeState.localeIdxFn;

    while true {
      // logDebug("getNextChunkIndex");
      const chunkIdx = localeState.getNextChunkIndex();
      if chunkIdx == -1 { break; }

      const chunk : range(int) = localeState._chunks[chunkIdx];
      // logDebug("chunk=", chunk);
      const (n, basisStatesPtr, coeffsPtr, targetIndicesPtr) =
        batchedOperator.computeOffDiag(
          chunk,
          localeState.representativesPtr + chunk.low,
          nil, // localeState.xPtr + chunk.low,
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
        // const accumulators = allocate(localeState.coeffType, chunk.size);
        // defer deallocate(accumulators);

        localProcessExperimental(localeState.basis,
                                 localeState.xPtr,
                                 n,
                                 basisStatesPtr,
                                 coeffsPtr:c_ptrConst(complex(128)),
                                 targetIndicesPtr,
                                 chunk.low,
                                 chunk.size,
                                 accumulators,
				 indices);
        foreach k in 0 ..# chunk.size do
          if accumulators[k] != 0 then
            localeState.accessor.localAdd(chunk.low + k, accumulators[k]);
      }
    } // end while true

    // no more chunks to produce, but we could still have data to analyze
    // logDebug("1 tryProcessReceivedData");
    localeState.tryProcessReceivedData();

    // tell other locales that we're done producing data
    if numLocales > 1 {
      for i in 0 ..# numLocales {
        // logDebug("queue.push");
        ref queue = localeState.remoteQueues[i];
        queue.push(new BufferFilledMessage(here.id.safeCast(int(16)), taskIdx.safeCast(int(16)), -1:int(32)));
      }
    }

    while !localeState.allCompleted() {
      // logDebug("2 tryProcessReceivedData");
      localeState.tryProcessReceivedData();
    }
  }
}

proc perLocaleOffDiagonal(matrix : Operator,
                          const ref xs : [] ?eltType,
                          ref ys : [] eltType,
                          const ref representatives : [] uint(64)) {
  const _timer = recordTime("perLocaleOffDiagonal");
  const numTasks = kNumTasks;
  const remoteBufferSize = max(kRemoteBufferSize, matrix.max_number_off_diag);
  const numChunks =
    min((representatives.size * matrix.max_number_off_diag + remoteBufferSize - 1) / remoteBufferSize,
        representatives.size);

  var localeState = new LocaleState(matrix, representatives, xs, ys, numChunks, numTasks, remoteBufferSize);
  globalLocaleStateStore[here.id] = localeState.getRemoteView();
  allLocalesBarrier.barrier();
  localeState.completeInitialization(globalLocaleStateStore);

  // TODO: do we need this?
  allLocalesBarrier.barrier();

  const maxChunkSize = max reduce [c in localeState._chunks] c.size;
  coforall taskIdx in 0 ..# numTasks with (ref localeState) {
    var worker = new Worker(taskIdx, localeState, maxChunkSize);
    worker.run();
  }

  // TODO: do we need this?
  allLocalesBarrier.barrier();
}

proc perLocaleMatrixVector(matrix : Operator,
                           const ref x : [] ?eltType,
                           ref y : [] eltType,
                           const ref representatives : [] uint(64)) {
  const _timer = recordTime("perLocaleMatrixVector");
  assert(matrix.locale == here);
  assert(x.locale == here);
  assert(y.locale == here);
  assert(representatives.locale == here);
  // logDebug("== START perLocaleMatrixVector");
  perLocaleDiagonal(matrix, x, y, representatives);
  if matrix.payload.deref().max_number_off_diag > 0 then
    perLocaleOffDiagonal(matrix, x, y, representatives);
  // logDebug("== FINISH perLocaleMatrixVector");
}

}
