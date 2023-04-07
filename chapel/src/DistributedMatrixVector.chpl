module DistributedMatrixVector {

use CTypes;
use RangeChunk;
use AllLocalesBarriers;
use Time;
use DynamicIters;
use CommDiagnostics;
use ChapelLocks;
import Random;

use CommonParameters;
use FFI;
use ForeignTypes;
use ConcurrentAccessor;
use BatchedOperator;
// use CommunicationQueue;

config const kVerboseComm : bool = false; 
config const kVerboseGetTiming : bool = false; 
config const kUseQueue : bool = false;


private proc meanAndErrString(timings : [] real) {
  const mean = (+ reduce timings) / timings.size:real;
  const variance =
    (1.0 / timings.size:real)
      * (+ reduce ([i in timings.domain]
                     (timings[i] - mean) * (timings[i] - mean)));
  const err = round(100 * sqrt(variance)) / 100;
  return mean:string + " ± " + err:string;
}

/* 
 */
private proc localDiagonalBatch(indices : range(int, BoundedRangeType.bounded, false),
                                matrix : Operator, const ref x : [] ?eltType, ref y : [] eltType,
                                const ref representatives : [] uint(64)) {
  const batchSize = indices.size;
  // assert(workspace.size >= batchSize);
  // if workspace.size < batchSize then
  //   workspace.domain = {0 ..# batchSize};
  ls_internal_operator_apply_diag_x1(
    matrix.payload, batchSize, c_const_ptrTo(representatives[indices.low]),
    c_ptrTo(y[indices.low]), c_const_ptrTo(x[indices.low]));
  // ls_hs_operator_apply_diag_kernel(
  //   matrix.payload, batchSize,
  //   c_const_ptrTo(representatives[indices.low]), 1,
  //   c_ptrTo(workspace));
  // foreach i in indices {
  //   y[i] = x[i] * workspace[i - indices.low]:eltType;
  // }
}

config const matrixVectorDiagonalNumChunks : int = 10 * here.maxTaskPar;
config const matrixVectorOffDiagonalNumChunks : int = 32 * here.maxTaskPar;
config const matrixVectorMainLoopNumTasks : int = here.maxTaskPar;

private proc localDiagonal(matrix : Operator, const ref x : [] ?eltType, ref y : [] eltType,
                           const ref representatives : [] uint(64),
                           numChunks : int = min(matrixVectorDiagonalNumChunks,
                                                 representatives.size)) {
  const totalSize = representatives.size;
  // const batchSize = (totalSize + numChunks - 1) / numChunks;
  var ranges : [0 ..# numChunks] range(int, BoundedRangeType.bounded, false) =
    chunks(0 ..# totalSize, numChunks);
  // var workspace : [0 ..# batchSize] complex(128) = noinit;
  forall r in ranges {
    localDiagonalBatch(r, matrix, x, y, representatives);
  }
}

private proc localProcess(basisPtr : c_ptr(Basis), accessorPtr : c_ptr(ConcurrentAccessor(?coeffType)),
                          basisStates : c_ptr(uint(64)), coeffs : c_ptr(?t), size : int) {
  var timer = new stopwatch();
  timer.start();
  var allocateTimer = new stopwatch();
  var indexingTimer = new stopwatch();
  var accessTimer = new stopwatch();
  local {
    // count == 0 has to be handled separately because c_ptrTo(indices) fails
    // when the size of indices is 0.
    if size == 0 then return (0, 0, 0, 0);

    // Special case when we don't have to call ls_hs_state_index
    if basisPtr.deref().isStateIndexIdentity() {
      accessTimer.start();
      ref accessor = accessorPtr.deref();
      foreach k in 0 ..# size {
        const i = basisStates[k]:int;
        const c = coeffs[k]:coeffType;
        accessor.localAdd(i, c);
      }
      accessTimer.stop();
    }
    else {
      allocateTimer.start();
      var indices : [0 ..# size] int = noinit;
      allocateTimer.stop();

      indexingTimer.start();
      ls_hs_state_index(basisPtr.deref().payload, size, basisStates, 1, c_ptrTo(indices[0]), 1);
      indexingTimer.stop();

      accessTimer.start();
      ref accessor = accessorPtr.deref();
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
      accessTimer.stop();
    }
  }
  timer.stop();
  return (timer.elapsed(), allocateTimer.elapsed(),
          indexingTimer.elapsed(), accessTimer.elapsed());
}


/*
var globalPtrStore : [LocaleSpace] (c_ptr(Basis), c_ptr(ConcurrentAccessor(real(64))));

private proc localOffDiagonal(matrix : Operator, const ref x : [] ?eltType, ref y : [] eltType,
                              const ref representatives : [] uint(64),
                              numChunks : int = min(matrixVectorOffDiagonalNumChunks,
                                                    representatives.size)) {
  var timer = new stopwatch();
  timer.start();
  var initTime : real;
  var computeOffDiagTime : atomic real;
  var stagingAddTime : atomic real;
  var stagingFlushTime : atomic real;
  var queueDrainTime : real;
  var initTimer = new stopwatch();
  initTimer.start();

  const chunkSize = (representatives.size + numChunks - 1) / numChunks;
  logDebug("Local dimension=", representatives.size, ", chunkSize=", chunkSize);
  var accessor = new ConcurrentAccessor(y);
  globalPtrStore[here.id] = (c_const_ptrTo(matrix.basis), c_ptrTo(accessor));
  allLocalesBarrier.barrier();

  var queue = new CommunicationQueue(eltType, globalPtrStore);

  const ranges : [0 ..# numChunks] range(int) =
    chunks(0 ..# representatives.size, numChunks);
  initTimer.stop();
  initTime += initTimer.elapsed();
  forall rangeIdx in dynamic(0 ..# numChunks, chunkSize=1,
                             numTasks=matrixVectorMainLoopNumTasks)
                     with (ref queue,
                           var batchedOperator = new BatchedOperator(matrix, chunkSize),
                           var staging = new StagingBuffers(queue)) {
    const r : range(int) = ranges[rangeIdx];
    // logDebug("Processing ", r, " ...");

    var timer = new stopwatch();
    timer.start();
    const (n, basisStatesPtr, coeffsPtr, keysPtr) = batchedOperator.computeOffDiag(
        r.size, c_const_ptrTo(representatives[r.low]), c_const_ptrTo(x[r.low]));
    chpl_task_yield();
    var radixOffsets : c_array(int, 257);
    radixOneStep(n, keysPtr, radixOffsets, basisStatesPtr, coeffsPtr);
    chpl_task_yield();
    // const n = offsetsPtr[r.size];
    timer.stop();
    computeOffDiagTime.add(timer.elapsed(), memoryOrder.relaxed);

    timer.clear();
    timer.start();
    var numberToDo = numLocales;
    var isComplete : [0 ..# numLocales] bool = false;
    while numberToDo > 0 {
      const prev = numberToDo;
      for localeIdx in 0 ..# numLocales {
        if isComplete[localeIdx] then continue;

        const k = radixOffsets[localeIdx];
        const n = radixOffsets[localeIdx + 1] - k;
        const isDone = if n != 0 then queue.tryEnqueue(localeIdx, n, basisStatesPtr + k, coeffsPtr + k)
                                 else true;
        if isDone {
          isComplete[localeIdx] = true;
          numberToDo -= 1;
        }
        // else
        //   chpl_task_yield();
      }
      if numberToDo == prev then chpl_task_yield();
    }
    // staging.add(n, basisStatesPtr, coeffsPtr);
    timer.stop();
    stagingAddTime.add(timer.elapsed(), memoryOrder.relaxed);
  }

  var queueTimer = new stopwatch();
  queueTimer.start();
  queue.drain();
  queueTimer.stop();
  queueDrainTime += queueTimer.elapsed();

  allLocalesBarrier.barrier();
  timer.stop();
  logDebug("localOffDiagonal took ", timer.elapsed(), "\n",
           "  ├─ ", initTime, " in initialization\n",
           "  ├─ ", computeOffDiagTime, " in computeOffDiag (total)\n",
           "  ├─ ", stagingAddTime, " in staging.add (total)\n",
           "  │   └─ ", queue.enqueueTime, " in queue.enqueue\n",
           "  │       ├─ ", queue.localProcessTimeHere, " in localProcess on here\n",
           "  │       ├─ ", queue.lockWaitingTimeHere, " waiting for queue.lock\n",
           "  │       └─ ", queue.enqueueUnsafeTime, " in queue._enqueueUnsafe\n",
           "  │           └─ ", queue.flushBufferTime, " in queue._flushBuffer\n",
           "  │               ├─ ", queue.lockWaitingTimeRemote, " waiting for remoteBuffer.lock\n",
           "  │               ├─ ", queue.flushBufferPutTime, " in remote PUTs\n",
           "  │               └─ ", queue.flushBufferRemoteTime, " in remote tasks\n",
           "  │                   └─ ", queue.localProcessTimeRemote, " in localProcess on remote\n",
           "  └─ ", queueDrainTime, " in queue.drain\n");
}
*/

record PartitionInfo {
    var _countOrOffset : int;
    var nextOffset : int;

    inline proc count ref { return _countOrOffset; }
    inline proc offset ref { return _countOrOffset; }
};

proc partitionBy(in first : c_ptr(?eltType), last : c_ptr(eltType), predicate) {
    while true {
      if first == last then return last;
      if !predicate(first.deref()) then break;
      first += 1;
    }

    var it = first + 1;
    while it != last {
      if predicate(it.deref()) {
        first.deref() <=> it.deref();
        first += 1;
      }
      it += 1;
    }
    return first;
}

inline proc swapElements(a : int, b : int, arr : c_ptr(?t1)) {
  arr[a] <=> arr[b];
}
inline proc swapElements(a : int, b : int, arr1 : c_ptr(?t1), arr2 : c_ptr(?t2)) {
  swapElements(a, b, arr1);
  swapElements(a, b, arr2);
}

proc radixOneStep(numKeys : int, keys : c_ptr(uint(8)), offsets : c_array(int, 257), arrs...?numArrs)
{
  var partitions : c_array(PartitionInfo, 256);
  foreach i in 0 ..# numKeys {
    partitions[keys[i]:int].count += 1;
  }

  var remainingPartitions : c_array(uint(8), 256);
  var numPartitions : int;
  var total : int;
  for i in 0 ..# 256 {
    const count = partitions[i].count;
    if count > 0 {
      partitions[i].offset = total;
      total += count;
      remainingPartitions[numPartitions] = i:uint(8);
      numPartitions += 1;
    }
    partitions[i].nextOffset = total;
  }

  var lastRemaining = remainingPartitions:c_ptr(uint(8)) + numPartitions;
  var endPartition = remainingPartitions:c_ptr(uint(8)) + 1;
  while lastRemaining - endPartition > 0 {
    record Func {
      inline proc this(partitionIdx : uint(8)) {
        ref beginOffset = partitions[partitionIdx:int].offset;
        ref endOffset = partitions[partitionIdx:int].nextOffset;
        if beginOffset == endOffset then return false;

        for i in beginOffset .. endOffset - 1 {
          ref offset = partitions[keys[i]:int].offset;
          keys[i] <=> keys[offset];
          swapElements(i, offset, (...arrs));
          offset += 1;
        }
        return beginOffset != endOffset;
      }
    }
    lastRemaining = partitionBy(remainingPartitions:c_ptr(uint(8)), lastRemaining, new Func());
  }

  offsets[0] = 0;
  foreach i in 1 ..# 256 {
    offsets[i] = partitions[i - 1].nextOffset;
  }
}

record _LocalBuffer {
  type coeffType;
  var destLocaleIdx : int;
  var srcLocaleIdx : int;
  var capacity : int;
  var size : int;
  var basisStates : c_ptr(uint(64));
  var coeffs : c_ptr(coeffType);
  var isFull : c_ptr(chpl__processorAtomicType(bool));
  var isEmpty : chpl__processorAtomicType(bool);
  var isEOF : chpl__processorAtomicType(bool);

  proc postinit() {
    assert(destLocaleIdx == here.id);
    isEmpty.write(true);
    isEOF.write(false);
    basisStates = c_malloc(uint(64), capacity);
    coeffs = c_malloc(coeffType, capacity);
  }

  proc deinit() {
    c_free(basisStates);
    c_free(coeffs);
  }
}

config const kDeadlock = false;

record _RemoteBuffer {
  type coeffType;
  var destLocaleIdx : int;
  var srcLocaleIdx : int;
  var capacity : int;
  var size : c_ptr(int);
  var basisStates : c_ptr(uint(64));
  var coeffs : c_ptr(coeffType);
  var isFull : chpl__processorAtomicType(bool);
  var isEmpty : c_ptr(chpl__processorAtomicType(bool));
  var isEOF : c_ptr(chpl__processorAtomicType(bool));

  var _putTime : real;
  var _waitTime : real;
  var _submitTime : real;

  proc postinit() {
    isFull.write(false);
  }

  inline proc put(localBasisStates : c_ptr(uint(64)),
                  localCoeffs : c_ptr(coeffType),
                  count : int) {
    assert(here.id == srcLocaleIdx);
    assert(count <= capacity);
    // assert(!isFull.read());
    // startVerboseCommHere();
    PUT(localBasisStates, destLocaleIdx, basisStates, count:c_size_t * c_sizeof(uint(64)));
    PUT(localCoeffs, destLocaleIdx, coeffs, count:c_size_t * c_sizeof(coeffType));
    // stopVerboseCommHere();
  }

  inline proc put(localBasisStates : [] uint(64),
                  localCoeffs : [] coeffType,
                  size : int) {
    put(c_const_ptrTo(localBasisStates[0]), c_const_ptrTo(localCoeffs[0]), size);
  }

  proc submit(basisStatesPtr : c_ptr(uint(64)),
              coeffsPtr : c_ptr(coeffType),
              count : int,
              inTheMeantime) {
    var timer = new stopwatch();
    timer.start();

    var waitTimer = new stopwatch();
    waitTimer.start();
    // Wait for the buffer to become empty
    while isFull.read() {
      inTheMeantime();
      if kDeadlock {
        chpl_task_yield();
      }
      // tryProcessLocal(taskIdx, newLocalBuffers, matrix, accessor);
    }
    waitTimer.stop();
    _waitTime += waitTimer.elapsed();

    put(basisStatesPtr, coeffsPtr, count);
    isFull.write(true);
    const atomicPtr = isEmpty;
    const sizePtr = size;
    on Locales[destLocaleIdx] {
      sizePtr.deref() = count;
      atomicStoreBool(atomicPtr, false);
    }

    timer.stop();
    _submitTime += timer.elapsed();
  }

  proc submit(basisStatesPtr : c_ptr(uint(64)),
              coeffsPtr : c_ptr(coeffType),
              count : int) {
    record InTheMeantime {
      inline proc this() {
        return;
      }
    }
    submit(basisStatesPtr, coeffsPtr, count, new InTheMeantime());
  }

  proc finish() {
    // startVerboseCommHere();
    const atomicPtr = isEOF;
    on Locales[destLocaleIdx] {
      atomicStoreBool(atomicPtr, true);
      // pragma "local fn" pragma "fast-on safe extern function"
      // extern proc atomic_store_bool(ref obj : chpl__processorAtomicType(bool), value : bool) : void;

      // atomic_store_bool(atomicPtr.deref(), true);
      // isEOF.deref().write(true);
    }
    // stopVerboseCommHere();
  }
}




class GlobalPtrStore {
  var arr : [LocaleSpace] (c_ptr(Basis),
                           c_ptr(ConcurrentAccessor(real(64))),
                           c_ptr(_RemoteBuffer(complex(128))),
                           c_ptr(_LocalBuffer(complex(128))),
                           int);
}

var globalPtrStoreNoQueue = new GlobalPtrStore();
// var globalPtrStoreNoQueue : [LocaleSpace] (c_ptr(Basis),
//                                            c_ptr(ConcurrentAccessor(real(64))),
//                                            c_ptr(_RemoteBuffer(complex(128))),
//                                            c_ptr(_LocalBuffer(complex(128))),
//                                            int);

config const kRemoteBufferSize = 150000;
config const kNumTasks = here.maxTaskPar;
config const kNumConsumerTasks = 1;
config const kVerbose = false;
config const kUseConsumer : bool = false;

extern proc chpl_task_getId(): chpl_taskID_t;

inline proc atomicStoreBool(p : c_ptr(chpl__processorAtomicType(bool)), value : bool) {
  pragma "local fn" pragma "fast-on safe extern function"
  extern proc atomic_store_bool(ref _obj : chpl__processorAtomicType(bool),
                                _value : bool) : void;

  atomic_store_bool(p.deref(), value);
}

/*
proc tryProcessLocal(taskIdx : int, srcLocaleIdx : int, newLocalBuffers, matrix, accessor) {
  ref localBuffer = newLocalBuffers[srcLocaleIdx, taskIdx];
  if !localBuffer.isEmpty.read() {
    // if kVerbose then
    //   logDebug("335: Calling localProcess for newLocalBuffers[", srcLocaleIdx, ", ", taskIdx, "]...");
    localProcess(c_const_ptrTo(matrix.basis), c_const_ptrTo(accessor),
                 localBuffer.basisStates,
                 localBuffer.coeffs,
                 localBuffer.size);
    localBuffer.isEmpty.write(true);
    const atomicPtr = localBuffer.isFull;
    on Locales[localBuffer.srcLocaleIdx] {
      atomicStoreBool(atomicPtr, false);
    }
    // if kVerbose then
    //   logDebug("335: done with localProcess for newLocalBuffers[", srcLocaleIdx, ", ", taskIdx, "]");
    return true;
  }
  return false;
}
proc tryProcessLocal(taskIdx : int, newLocalBuffers, matrix, accessor) {
  var hasDoneWork : bool = false;
  for srcLocaleIdx in 0 ..# numLocales {
    hasDoneWork = hasDoneWork ||
      tryProcessLocal(taskIdx, srcLocaleIdx, newLocalBuffers, matrix, accessor);
  }
  return hasDoneWork;
}
*/

/*
record LocalOffDiagonalTimers {
  var total : stopwatch;
  var initialization : stopwatch;
  var computeOffDiag : atomic real;
  var submit : real;
  var wait : real;
  var put : real;
}
*/

proc _offDiagMakeLocalBuffers(numTasks : int, remoteBufferSize : int) {
  var localBuffers : [0 ..# numLocales, 0 ..# numTasks] _LocalBuffer(complex(128))
    = [(srcLocaleIdx, taskIdx) in {0 ..# numLocales, 0 ..# numTasks}]
        new _LocalBuffer(complex(128), here.id, srcLocaleIdx, remoteBufferSize);
  return localBuffers;
}

proc _offDiagMakeRemoteBuffers(numTasks : int, remoteBufferSize : int) {
  var remoteBuffers : [0 ..# numLocales, 0 ..# numTasks] _RemoteBuffer(complex(128))
    = [(destLocaleIdx, taskIdx) in {0 ..# numLocales, 0 ..# numTasks}]
        new _RemoteBuffer(complex(128), destLocaleIdx, here.id, remoteBufferSize);
  return remoteBuffers;
}

proc _offDiagInitLocalBuffers(numTasks : int, ref localBuffers, const ref ptrStore) {
  const destLocaleIdx = localBuffers.locale.id;
  coforall loc in Locales do on loc {
    const srcLocaleIdx = loc.id;
    const (_basis, _accessor, _remoteBufferPtr,
           _localBufferPtr, _numChunks) = ptrStore[srcLocaleIdx];
    for taskIdx in 0 ..# numTasks {
      ref remoteBuffer = (_remoteBufferPtr + destLocaleIdx * numTasks + taskIdx).deref();
      localBuffers[srcLocaleIdx, taskIdx].isFull = c_ptrTo(remoteBuffer.isFull);
    }
  }
}

proc _offDiagInitRemoteBuffers(numTasks : int, ref remoteBuffers, const ref ptrStore) {
  const srcLocaleIdx = remoteBuffers.locale.id;
  coforall loc in Locales do on loc {
    const destLocaleIdx = loc.id;
    const (_basis, _accessor, _remoteBufferPtr,
           _localBufferPtr, _numChunks) = ptrStore[destLocaleIdx];
    for taskIdx in 0 ..# numTasks {
      ref myLocalBuffer = (_localBufferPtr + srcLocaleIdx * numTasks + taskIdx).deref();
      ref remoteBuffer = remoteBuffers[destLocaleIdx, taskIdx];
      remoteBuffer.basisStates = myLocalBuffer.basisStates;
      remoteBuffer.coeffs = myLocalBuffer.coeffs;
      remoteBuffer.size = c_ptrTo(myLocalBuffer.size);
      remoteBuffer.isEmpty = c_ptrTo(myLocalBuffer.isEmpty);
      remoteBuffer.isEOF = c_ptrTo(myLocalBuffer.isEOF);
    }
  }
}

record Producer {
  type eltType;
  var _taskIdx : int;
  var numChunks : int;
  var numProducerTasks : int;

  var batchedOperator : BatchedOperator;
  var basisPtr : c_ptr(Basis);
  var accessorPtr : c_ptr(ConcurrentAccessor(eltType));
  var representativesPtr : c_ptr(uint(64));
  var xPtr : c_ptr(eltType);

  var rangesPtr : c_ptr(range(int));
  var moreWorkPtr : c_ptr(atomic bool);
  var currentChunkIdxPtr : c_ptr(atomic int);

  var remoteBuffersPtr : c_ptr(_RemoteBuffer(complex(128)));

  var runTimer : stopwatch;
  var computeOffDiagTimer : stopwatch;
  var radixOneStepTimer : stopwatch;
  var localProcessTimer : stopwatch;
  var localProcessAllocating : real;
  var localProcessIndexing : real;
  var localProcessAccessing : real;
  var putSize : int;
  var putTimer : stopwatch;
  var fastOnTimer : stopwatch;
  var submitTimer : stopwatch;

  proc init(taskIdx : int, numChunks : int, in batchedOperator : BatchedOperator,
            ref accessor : ConcurrentAccessor(?eltType),
            const ref representatives : [] uint(64),
            const ref x : [] eltType,
            const ref ranges : [] range(int),
            ref moreWork : atomic bool,
            ref currentChunkIdx : atomic int,
            ref remoteBuffers : [] _RemoteBuffer(complex(128))) {
    // logDebug("Creating Producer(", taskIdx, ")...");
    this.eltType = eltType;
    this._taskIdx = taskIdx;
    this.numChunks = numChunks;
    this.numProducerTasks = remoteBuffers.domain.dim(1).size;
    this.batchedOperator = batchedOperator;
    this.basisPtr = c_const_ptrTo(this.batchedOperator._matrixPtr.deref().basis);
    this.accessorPtr = c_ptrTo(accessor);
    this.representativesPtr = c_const_ptrTo(representatives);
    this.xPtr = c_const_ptrTo(x);
    this.rangesPtr = c_const_ptrTo(ranges);
    this.moreWorkPtr = c_ptrTo(moreWork);
    this.currentChunkIdxPtr = c_ptrTo(currentChunkIdx);
    this.remoteBuffersPtr = c_ptrTo(remoteBuffers);

    this.runTimer = new stopwatch();
    this.computeOffDiagTimer = new stopwatch();
    this.radixOneStepTimer = new stopwatch();
    this.localProcessTimer = new stopwatch();
    this.localProcessAllocating = 0;
    this.localProcessIndexing = 0;
    this.localProcessAccessing = 0;
    this.putSize = 0;
    this.putTimer = new stopwatch();
    this.fastOnTimer = new stopwatch();
    this.submitTimer = new stopwatch();
    // logDebug("Done creating Producer(", taskIdx, ")...");
  }

  inline proc remoteBuffers(localeIdx : int, taskIdx : int) ref {
    assert(0 <= localeIdx && localeIdx < numLocales);
    assert(0 <= taskIdx && taskIdx < numProducerTasks);
    return remoteBuffersPtr[localeIdx * numProducerTasks + taskIdx];
  }

  proc bandwidth() : real {
    const sentBytes = putSize * (c_sizeof(uint(64)) + c_sizeof(complex(128))):int;
    const sentTime = putTimer.elapsed();
    return sentBytes / (1024.0 * 1024.0 * 1024.0) / sentTime;
  }

  proc trySubmit(ref remoteBuffer,
                 basisStatesPtr : c_ptr(uint(64)),
                 coeffsPtr : c_ptr(complex(128)),
                 count : int) {
    if remoteBuffer.isFull.read() then
      return false;

    putTimer.start();
    remoteBuffer.put(basisStatesPtr, coeffsPtr, count);
    putTimer.stop();
    putSize += count;

    remoteBuffer.isFull.write(true);

    const atomicPtr = remoteBuffer.isEmpty;
    const sizePtr = remoteBuffer.size;
    fastOnTimer.start();
    on Locales[remoteBuffer.destLocaleIdx] {
      sizePtr.deref() = count;
      atomicStoreBool(atomicPtr, false);
    }
    fastOnTimer.stop();
    return true;
  }

  proc run() {
    assert(numLocales <= 256);
    runTimer.start();

    var radixOffsets : c_array(int, 257);
    var submitted : [0 ..# numLocales] bool;

    while moreWorkPtr.deref().read() {
      const rangeIdx = currentChunkIdxPtr.deref().fetchAdd(1);
      // Multiple threads passed moreWork.read() at once.
      // All whose fetchAdd() was after the one
      // that grabbed the final chunk just break.
      if rangeIdx >= numChunks then
        break;
      // Final rangeIdx -- signal that to everybody
      if rangeIdx == numChunks - 1 then
        moreWorkPtr.deref().write(false);

      // Compute a r.size rows of the matrix
      computeOffDiagTimer.start();
      const r : range(int) = rangesPtr[rangeIdx];
      const (n, basisStatesPtr, coeffsPtr, keysPtr) = batchedOperator.computeOffDiag(
          r.size, representativesPtr + r.low, xPtr + r.low);
      computeOffDiagTimer.stop();

      radixOneStepTimer.start();
      radixOneStep(n, keysPtr, radixOffsets, basisStatesPtr, coeffsPtr);
      radixOneStepTimer.stop();

      c_memset(c_ptrTo(submitted[0]), 0, numLocales:c_size_t * c_sizeof(bool));
      var remaining = numLocales;

      {
        const destLocaleIdx = here.id;
        const k = radixOffsets[destLocaleIdx];
        const n = radixOffsets[destLocaleIdx + 1] - k;
        ref remoteBuffer = remoteBuffers[destLocaleIdx, _taskIdx];
        localProcessTimer.start();
        assert(!remoteBuffer.isFull.read());
        const (_, allocTime, indexTime, accessTime) =
          localProcess(basisPtr,
                       accessorPtr,
                       basisStatesPtr + k,
                       coeffsPtr + k,
                       n);
        localProcessTimer.stop();
        localProcessAllocating += allocTime;
        localProcessIndexing += indexTime;
        localProcessAccessing += accessTime;

        submitted[destLocaleIdx] = true;
        remaining -= 1;
      }
      submitTimer.start();
      while remaining > 0 {
        for destLocaleIdx in 0 ..# numLocales {
          if submitted[destLocaleIdx] then continue;

          const k = radixOffsets[destLocaleIdx];
          const n = radixOffsets[destLocaleIdx + 1] - k;
          ref remoteBuffer = remoteBuffers[destLocaleIdx, _taskIdx];
          if trySubmit(remoteBuffer, basisStatesPtr + k, coeffsPtr + k, n) {
            submitted[destLocaleIdx] = true;
            remaining -= 1;
          }
        }
      }
      submitTimer.stop();
    }

    runTimer.stop();
  }
}

config const kShuffle : bool = false;

record Consumer {
  type eltType;

  var _taskIdx : int;
  var numConsumerTasks : int;
  var numProducerTasks : int;

  var slots;

  var basisPtr : c_ptr(Basis);
  var accessorPtr : c_ptr(ConcurrentAccessor(eltType));

  var totalNumberChunks : int;
  var numProcessedPtr : c_ptr(atomic int);

  var localBuffersPtr : c_ptr(_LocalBuffer(complex(128)));

  var runTimer : stopwatch;
  var localProcessTimer : stopwatch;
  var localProcessAllocating : real;
  var localProcessIndexing : real;
  var localProcessAccessing : real;
  var fastOnTimer : stopwatch;

  proc init(taskIdx : int, numConsumerTasks : int, numProducerTasks : int,
            const ref basis : Basis, ref accessor : ConcurrentAccessor(?eltType),
            totalNumberChunks : int, ref numProcessedChunks : atomic int,
            ref localBuffers : [] _LocalBuffer(complex(128))) {
    // logDebug("Creating Consumer(", taskIdx, ")...");
    this.eltType = eltType;
    this._taskIdx = taskIdx;
    this.numConsumerTasks = numConsumerTasks;
    this.numProducerTasks = numProducerTasks;
   
    // logDebug(0 ..# (numLocales - 1) * numProducerTasks);
    const everything = 0 ..# (numLocales - 1) * numProducerTasks;
    if taskIdx < everything.size {
      const r = chunk(everything, numConsumerTasks, taskIdx);
      var pairs : [0 ..# r.size] (int, int);
      var size = 0;
      var offset = 0;
      for localeIdx in 0 ..# numLocales {
        if localeIdx == here.id then continue;
        for otherTaskIdx in 0 ..# numProducerTasks {
          if r.first <= offset && offset <= r.last {
            pairs[size] = (localeIdx, otherTaskIdx);
            size += 1;
          }
          offset += 1;
        }
      }
      if kShuffle then
        Random.shuffle(pairs);

      this.slots = pairs;
      // logDebug("Slots: ", slots);
    }
    else {
      var pairs : [0 ..# 0] (int, int);
      this.slots = pairs;
    }

    this.basisPtr = c_const_ptrTo(basis);
    this.accessorPtr = c_ptrTo(accessor);
    this.totalNumberChunks = totalNumberChunks;
    this.numProcessedPtr = c_ptrTo(numProcessedChunks);
    this.localBuffersPtr = c_ptrTo(localBuffers[0, 0]);

    this.runTimer = new stopwatch();
    this.localProcessTimer = new stopwatch();
    this.fastOnTimer = new stopwatch();
  }

  proc localBuffers(localeIdx : int, taskIdx : int) ref {
    assert(0 <= localeIdx && localeIdx < numLocales);
    assert(0 <= taskIdx && taskIdx < numProducerTasks);
    return localBuffersPtr[localeIdx * numProducerTasks + taskIdx];
  }

  proc run() {
    runTimer.start();
    while numProcessedPtr.deref().read() < totalNumberChunks {
      var hasDoneWork = false;
      for (localeIdx, otherTaskIdx) in slots {
        ref localBuffer = localBuffers[localeIdx, otherTaskIdx];
        if !localBuffer.isEmpty.read() {
          local {
            localProcessTimer.start();
            numProcessedPtr.deref().add(1);
            const (_, allocTime, indexTime, accessTime) =
              localProcess(basisPtr, accessorPtr,
                           localBuffer.basisStates,
                           localBuffer.coeffs,
                           localBuffer.size);
            localBuffer.isEmpty.write(true);
            localProcessTimer.stop();
            localProcessAllocating += allocTime;
            localProcessIndexing += indexTime;
            localProcessAccessing += accessTime;
          }
          const atomicPtr = localBuffer.isFull;

          fastOnTimer.start();
          on Locales[localBuffer.srcLocaleIdx] {
            atomicStoreBool(atomicPtr, false);
          }
          fastOnTimer.stop();

          hasDoneWork = true;
        }
      }
    }
    runTimer.stop();
  }
}


private proc localOffDiagonalNoQueue(matrix : Operator, const ref x : [] ?eltType, ref y : [] eltType,
                                     const ref representatives : [] uint(64)) {
  // logDebug("Calling localOffDiagonalNoQueue...");
  var totalTimer = new stopwatch();
  var initializationTimer = new stopwatch();
  totalTimer.start();
  initializationTimer.start();

  const numTasks = kNumTasks;
  // If there is only one locale, there's no need to process tasks from other locales :)
  const numConsumerTasks =
    if numLocales == 1 then 0
                       else max(min(kNumConsumerTasks, numTasks - 1), 1);
  const numProducerTasks = numTasks - numConsumerTasks;

  const remoteBufferSize = max(kRemoteBufferSize, matrix.numberOffDiagTerms());
  // Other locales need access to our newLocalBuffers, newRemoteBuffers, etc.
  // so we create and store them in a common location (on locale 0) before a
  // barrier.
  var newLocalBuffers = _offDiagMakeLocalBuffers(numProducerTasks, remoteBufferSize);
  var newRemoteBuffers = _offDiagMakeRemoteBuffers(numProducerTasks, remoteBufferSize);
  var accessor = new ConcurrentAccessor(y);

  const numChunks =
    min(max((representatives.size * matrix.numberOffDiagTerms()
               + remoteBufferSize - 1) / remoteBufferSize,
            10 * numProducerTasks),
        representatives.size);
  // logDebug(matrix.basis, ", ", (newRemoteBuffers.dim(0).size, newRemoteBuffers.dim(1).size),
  //                        ", ", (newLocalBuffers.dim(0).size, newLocalBuffers.dim(1).size));
  // logDebug((c_const_ptrTo(matrix.basis), c_ptrTo(accessor),
  //          c_ptrTo(newRemoteBuffers[0, 0]), c_ptrTo(newLocalBuffers[0, 0]), numChunks));
  // logDebug("before: ", globalPtrStoreNoQueue.arr[here.id]);
  globalPtrStoreNoQueue.arr[here.id] = (c_const_ptrTo(matrix.basis),
                                    c_ptrTo(accessor),
                                    c_ptrTo(newRemoteBuffers[0, 0]),
                                    c_ptrTo(newLocalBuffers[0, 0]),
                                    numChunks);
  // logDebug("after: ", globalPtrStoreNoQueue.arr[here.id]);
  allLocalesBarrier.barrier();

  // logDebug("Check #1");
  const ptrStore : [0 ..# numLocales] globalPtrStoreNoQueue.arr.eltType = globalPtrStoreNoQueue.arr;
  // logDebug("Check #2");
  _offDiagInitLocalBuffers(numProducerTasks, newLocalBuffers, ptrStore);
  // logDebug("Check #3");
  _offDiagInitRemoteBuffers(numProducerTasks, newRemoteBuffers, ptrStore);
  // logDebug("Check #4");

  const ranges : [0 ..# numChunks] range(int) = chunks(0 ..# representatives.size, numChunks);
  // logDebug("Check #5");
  const batchedOperatorChunkSize = (representatives.size + numChunks - 1) / numChunks;
  var moreWork : atomic bool = true;
  var curChunkIdx : atomic int = 0;
  // logDebug("Check #6");

  var totalNumberChunks = 0;
  for localeIdx in 0 ..# numLocales {
    if localeIdx != here.id {
      const (_basis, _accessor, _remoteBufferPtr,
             _localBufferPtr, _numChunks) = ptrStore[localeIdx];
      totalNumberChunks += _numChunks;
    }
  }
  // logDebug("Check #7");
  
  var numProcessedChunks : atomic int = 0;

  // logDebug("numChunks = ", numChunks, ", chunkSize = ",
  //          batchedOperatorChunkSize, ", numTasks = ", numTasks,
  //          ", representatives.size = ", representatives.size,
  //          ", numberOffDiagTerms = ", matrix.numberOffDiagTerms());

  var producerRunTime : [0 ..# numProducerTasks] real;
  var producerComputeOffDiagTime : [0 ..# numProducerTasks] real;
  var producerApplyOffDiagTime : [0 ..# numProducerTasks] real;
  var producerMemcpyTime : [0 ..# numProducerTasks] real;
  var producerStateInfoTime : [0 ..# numProducerTasks] real;
  var producerCoeffTime : [0 ..# numProducerTasks] real;
  var producerKeysTime : [0 ..# numProducerTasks] real;
  var producerRadixOneStepTime : [0 ..# numProducerTasks] real;
  var producerLocalProcessTime : [0 ..# numProducerTasks] real;
  var producerAllocTime : [0 ..# numProducerTasks] real;
  var producerIndexTime : [0 ..# numProducerTasks] real;
  var producerAccessTime : [0 ..# numProducerTasks] real;
  var producerSubmitTime : [0 ..# numProducerTasks] real;
  var producerFastOnTime : [0 ..# numProducerTasks] real;
  var producerPutTime : [0 ..# numProducerTasks] real;
  var producerBandwidth : [0 ..# numProducerTasks] real;

  var consumerRunTime : [0 ..# numConsumerTasks] real;
  var consumerLocalProcessTime : [0 ..# numConsumerTasks] real;
  var consumerAllocTime : [0 ..# numConsumerTasks] real;
  var consumerIndexTime : [0 ..# numConsumerTasks] real;
  var consumerAccessTime : [0 ..# numConsumerTasks] real;
  var consumerFastOnTime : [0 ..# numConsumerTasks] real;

  // logDebug("Check #8");
  allLocalesBarrier.barrier();
  initializationTimer.stop();

  coforall taskIdx in 0 ..# numTasks with (ref accessor) {
    if taskIdx < numProducerTasks {
      // I'm a producer
      var producer = new Producer(
        taskIdx,
        numChunks,
        new BatchedOperator(matrix, batchedOperatorChunkSize),
        accessor,
        representatives,
        x,
        ranges,
        moreWork,
        curChunkIdx,
        newRemoteBuffers);
      producer.run();

      producerRunTime[taskIdx] = producer.runTimer.elapsed();
      producerComputeOffDiagTime[taskIdx] = producer.computeOffDiagTimer.elapsed();
      producerApplyOffDiagTime[taskIdx] = producer.batchedOperator.applyOffDiagTimer.elapsed();
      producerMemcpyTime[taskIdx] = producer.batchedOperator.memcpyTimer.elapsed();
      producerStateInfoTime[taskIdx] = producer.batchedOperator.stateInfoTimer.elapsed();
      producerCoeffTime[taskIdx] = producer.batchedOperator.coeffTimer.elapsed();
      producerKeysTime[taskIdx] = producer.batchedOperator.keysTimer.elapsed();
      producerRadixOneStepTime[taskIdx] = producer.radixOneStepTimer.elapsed();
      producerLocalProcessTime[taskIdx] = producer.localProcessTimer.elapsed();
      producerAllocTime[taskIdx] = producer.localProcessAllocating;
      producerIndexTime[taskIdx] = producer.localProcessIndexing;
      producerAccessTime[taskIdx] = producer.localProcessAccessing;
      producerSubmitTime[taskIdx] = producer.submitTimer.elapsed();
      producerFastOnTime[taskIdx] = producer.fastOnTimer.elapsed();
      producerPutTime[taskIdx] = producer.putTimer.elapsed();
      producerBandwidth[taskIdx] = producer.bandwidth();
    }
    else {
      // I'm a consumer
      var consumer = new Consumer(
        taskIdx - numProducerTasks,
        numConsumerTasks,
        numProducerTasks,
        matrix.basis,
        accessor,
        totalNumberChunks,
        numProcessedChunks,
        newLocalBuffers);
      consumer.run();

      const i = taskIdx - numProducerTasks;
      consumerRunTime[i] = consumer.runTimer.elapsed();
      consumerLocalProcessTime[i] = consumer.localProcessTimer.elapsed();
      consumerAllocTime[i] = consumer.localProcessAllocating;
      consumerIndexTime[i] = consumer.localProcessIndexing;
      consumerAccessTime[i] = consumer.localProcessAccessing;
      consumerFastOnTime[i] = consumer.fastOnTimer.elapsed();
    }
  }

  allLocalesBarrier.barrier();

  for srcLocaleIdx in 0 ..# numLocales {
    for taskIdx in 0 ..# numProducerTasks {
      // assert(newLocalBuffers[srcLocaleIdx, taskIdx].isEOF.read());
      assert(newLocalBuffers[srcLocaleIdx, taskIdx].isEmpty.read());
    }
  }
  for srcLocaleIdx in 0 ..# numLocales {
    for taskIdx in 0 ..# numProducerTasks {
      assert(!newRemoteBuffers[srcLocaleIdx, taskIdx].isFull.read());
    }
  }

  totalTimer.stop();
  if kDisplayTimings then
    logDebug("localOffDiagonalNoQueue: ", totalTimer.elapsed(), "\n",
             " ├─ ", initializationTimer.elapsed(), " in initialization\n",
             " ├─ producers: ", meanAndErrString(producerRunTime), "\n",
             " │   ├─ computeOffDiag: ", meanAndErrString(producerComputeOffDiagTime), "\n",
             " │   │   ├─ applyOffDiag: ", meanAndErrString(producerApplyOffDiagTime), "\n",
             " │   │   ├─ memcpy:       ", meanAndErrString(producerMemcpyTime), "\n",
             " │   │   ├─ stateInfo:    ", meanAndErrString(producerStateInfoTime), "\n",
             " │   │   ├─ rescale:      ", meanAndErrString(producerCoeffTime), "\n",
             " │   │   └─ localeIdxOf:  ", meanAndErrString(producerKeysTime), "\n",
             " │   ├─ radixOneStep:   ", meanAndErrString(producerRadixOneStepTime), "\n",
             " │   ├─ localProcess:   ", meanAndErrString(producerLocalProcessTime), "\n",
             " │   │   ├─ allocating: ", meanAndErrString(producerAllocTime), "\n",
             " │   │   ├─ indexing:   ", meanAndErrString(producerIndexTime), "\n",
             " │   │   └─ accessing:  ", meanAndErrString(producerAccessTime), "\n",
             " │   └─ submit:         ", meanAndErrString(producerSubmitTime), "\n",
             " │       ├─ PUT:    ", meanAndErrString(producerPutTime), "\n",
             " │       └─ fastOn: ", meanAndErrString(producerFastOnTime), "\n",
             " └─ consumers: ", meanAndErrString(consumerRunTime), "\n",
             "     ├─ localProcess: ", meanAndErrString(consumerLocalProcessTime), "\n",
             "     │   ├─ allocating: ", meanAndErrString(consumerAllocTime), "\n",
             "     │   ├─ indexing:   ", meanAndErrString(consumerIndexTime), "\n",
             "     │   └─ accessing:  ", meanAndErrString(consumerAccessTime), "\n",
             "     └─ fastOn:       ", meanAndErrString(consumerFastOnTime), "\n",
             "    (bandwidth in GB/s: ", meanAndErrString(producerBandwidth), ")");
}

proc localMatrixVector(matrix : Operator, const ref x : [] ?eltType, ref y : [] eltType,
                               const ref representatives : [] uint(64)) {
  // logDebug("Calling localMatrixVector...");
  assert(matrix.locale == here);
  assert(x.locale == here);
  assert(y.locale == here);
  assert(representatives.locale == here);
  if matrix.numberDiagTerms() > 0 then
    localDiagonal(matrix, x, y, representatives);
  // logDebug("Done with diagonal");
  // if kUseQueue then
  //   localOffDiagonal(matrix, x, y, representatives);
  // else
  if matrix.numberOffDiagTerms() > 0 then
    localOffDiagonalNoQueue(matrix, x, y, representatives);
}

proc matrixVectorProduct(const ref matrix : Operator,
                         const ref x,
                         ref y,
                         const ref representatives) {
  // logDebug("Calling matrixVectorProduct...");
  coforall loc in Locales with (ref y) do on loc {
    // var (_, myMatrix) = loadConfigFromYaml(matrixFilename, hamiltonian=true);
    const myMatrix = matrix;
    const ref myX = x.getBlock(loc.id)[0, ..];
    ref myY = y.getBlock(loc.id)[0, ..];
    const ref myBasisStates = representatives.getBlock(loc.id);
    // logDebug("Setting representatives...");
    myMatrix.basis.uncheckedSetRepresentatives(myBasisStates);
    // logDebug("Done setting representatives");
    var timer = new stopwatch();
    timer.start();
    localMatrixVector(myMatrix, myX, myY, myBasisStates);
    timer.stop();
    if kDisplayTimings then
      logDebug("Spent ", timer.elapsed(), " in localMatrixVector");
  }
}

export proc ls_chpl_matrix_vector_product(matrixPtr : c_ptr(ls_hs_operator), numVectors : c_int,
                                          xPtr : c_ptr(real(64)), yPtr : c_ptr(real(64))) {
  // logDebug("Calling ls_chpl_matrix_vector_product ...");
  var matrix = new Operator(matrixPtr, owning=false);
  if matrix.basis.numberWords != 1 then
    halt("bases with more than 64 bits are not yet implemented");
  if numVectors != 1 then
    halt("applying the Operator to more than 1 vector is not yet implemented");

  const ref representatives = matrix.basis.representatives();
  const numStates = representatives.size;
  var x = makeArrayFromPtr(xPtr, (numStates,));
  var y = makeArrayFromPtr(yPtr, (numStates,));
  localMatrixVector(matrix, x, y, representatives);
  // logDebug("Done!");
}

/*
proc localMatrixVectorPart(
      H : Operator,
      const ref x : [] ?eltType,
      ref y : [] eltType,
      const ref representatives : [] uint(64)) {
  y = 0;

  // Diagonal contribution
  forall (sigma, i) in zip(representatives, 0..) {

  }


  // Off-diagonal contribution

  coforall loc in Locales do on loc {
    var accessor = new unmanaged ConcurrentAccessor();

    for i,  in zip
    
  }


  for i in 0 ..# x.size {
    // representatives is a distributed vector of basis states. It is pre-computed
    const |σᵢ⟩ = representatives[i];
    for (cⱼ, |σⱼ⟩) in H|σᵢ⟩ {
      const t = x[i] * cⱼ;
      const j = indexOf(|σⱼ⟩);
      y[j] += t;
    }
  }

  y = 0; // initialize y with zeros
  for i in 0 ..# x.size {
    // representatives is a distributed vector of basis states. It is pre-computed
    const |σᵢ⟩ = representatives[i];
    for (cⱼ, |σⱼ⟩) in H|σᵢ⟩ {
      const t = x[i] * cⱼ;
      const j = indexOf(|σⱼ⟩);
      y[j] += t;
    }
  }
}
*/

}
