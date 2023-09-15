module DistributedMatrixVector {

use CTypes;
use RangeChunk;
use AllLocalesBarriers;
use List;
use Time;
use DynamicIters;
use CommDiagnostics;
use ChapelLocks;
import MemDiagnostics;
import Random;
import OS.POSIX;

use CommonParameters;
use FFI;
use Utils;
use ForeignTypes;
use ConcurrentAccessor;
use BatchedOperator;

config const matrixVectorDiagonalNumChunks : int = here.maxTaskPar;

private proc localDiagonalBatch(indices : range(int, boundKind.both, strideKind.one),
                                matrix : Operator, const ref x : [] ?eltType, ref y : [] eltType,
                                const ref representatives : [] uint(64)) {
  ls_internal_operator_apply_diag_x1(
    matrix.payload.deref(),
    indices.size,
    c_ptrToConst(representatives[indices.low]),
    c_ptrTo(y[indices.low]),
    c_ptrToConst(x[indices.low]));
}

private proc localDiagonal(matrix : Operator,
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
    localDiagonalBatch(r, matrix, x, y, representatives);
  }

  timer.stop();
  if timings != nil then
    timings!.addChild("localDiagonal", timer.elapsed());
}


record LocalProcessTimer {
  var total : stopwatch;
  var allocation : stopwatch;
  var indexing : stopwatch;
  var accessing : stopwatch;

  proc getTimings() {
    return timingTree(
      "localProcess", total.elapsed(),
      [ ("allocation", allocation.elapsed())
      , ("ls_hs_state_index", indexing.elapsed())
      , ("localAdd", accessing.elapsed())
      ]
    );
  }
}

private proc localProcess(basisPtr : c_ptrConst(Basis),
                          accessorPtr : c_ptr(ConcurrentAccessor(?coeffType)),
                          basisStates : c_ptrConst(uint(64)),
                          coeffs : c_ptrConst(?t),
                          size : int,
                          ref timer : LocalProcessTimer) {
  local {
    // count == 0 has to be handled separately because c_ptrTo(indices) fails
    // when the size of indices is 0.
    if size == 0 then return;
    timer.total.start();

    ref accessor = accessorPtr.deref();
    // Special case when we don't have to call ls_hs_state_index
    if numLocales == 1 && basisPtr.deref().isStateIndexIdentity() {
      timer.accessing.start();
      foreach k in 0 ..# size {
        const i = basisStates[k]:int;
        const c = coeffs[k]:coeffType;
        accessor.localAdd(i, c);
      }
      timer.accessing.stop();
    }
    else {
      timer.allocation.start();
      var indices : [0 ..# size] int = noinit;
      timer.allocation.stop();

      timer.indexing.start();
      ls_hs_state_index(basisPtr.deref().payload, size, basisStates, 1, c_ptrTo(indices[0]), 1);
      timer.indexing.stop();

      timer.accessing.start();
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
      timer.accessing.stop();
    }

    timer.total.stop();
  }
}

record PartitionInfo {
    var _countOrOffset : int;
    var nextOffset : int;

    inline proc count ref { return _countOrOffset; }
    inline proc offset ref { return _countOrOffset; }
};

inline proc partitionBy(in first : c_ptr(?eltType), last : c_ptr(eltType), predicate) {
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

proc radixOneStep(numKeys : int, keys : c_ptr(uint(8)), offsets : c_array(int, 257), arrs...?numArrs) {
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
  var incrementNumProcessed : chpl__processorAtomicType(bool);

  proc postinit() {
    assert(destLocaleIdx == here.id);
    incrementNumProcessed.write(true);
    isEmpty.write(true);
    basisStates = allocate(uint(64), capacity);
    coeffs = allocate(coeffType, capacity);
  }

  proc deinit() {
    deallocate(basisStates);
    deallocate(coeffs);
  }
}

// record RemoteBufferTimer {
//   var put : stopwatch;
//   var putSize : int;
//   var wait : stopwatch;
//   var waitCalls : int;
//   var waitCount : int;
//   var submit : stopwatch;
// }

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
  var incrementNumProcessed : c_ptr(chpl__processorAtomicType(bool));

  // var timer : RemoteBufferTimer;

  proc postinit() {
    isFull.write(false);
  }

  inline proc put(localBasisStates : c_ptrConst(uint(64)),
                  localCoeffs : c_ptrConst(coeffType),
                  count : int) {
    // timer.put.start();
    assert(here.id == srcLocaleIdx);
    assert(count <= capacity);
    PUT(localBasisStates, destLocaleIdx, basisStates, count:c_size_t * c_sizeof(uint(64)));
    PUT(localCoeffs, destLocaleIdx, coeffs, count:c_size_t * c_sizeof(coeffType));
    // timer.putSize += count * (c_sizeof(coeffTimer) + c_sizeof(uint(64))):int;
    // timer.put.stop();
  }

  inline proc put(localBasisStates : [] uint(64),
                  localCoeffs : [] coeffType,
                  size : int) {
    put(c_ptrToConst(localBasisStates[0]), c_ptrToConst(localCoeffs[0]), size);
  }

  // proc submit(basisStatesPtr : c_ptrConst(uint(64)),
  //             coeffsPtr : c_ptrConst(coeffType),
  //             count : int,
  //             inTheMeantime) {
  //   timer.submit.start();

  //   timer.wait.start();
  //   timer.waitCalls += 1; // Count number of times we had to wait
  //   // Wait for the buffer to become empty
  //   while isFull.read() {
  //     inTheMeantime();
  //     timer.waitCount += 1; // Count busy-waiting iterations
  //   }
  //   timer.wait.stop();

  //   put(basisStatesPtr, coeffsPtr, count);
  //   isFull.write(true);
  //   const atomicPtr = isEmpty;
  //   const sizePtr = size;
  //   on Locales[destLocaleIdx] {
  //     sizePtr.deref() = count;
  //     atomicStoreBool(atomicPtr, false);
  //   }

  //   timer.submit.stop();
  // }

  // proc submit(basisStatesPtr : c_ptr(uint(64)),
  //             coeffsPtr : c_ptr(coeffType),
  //             count : int) {
  //   record InTheMeantime {
  //     inline proc this() { return; }
  //   }
  //   submit(basisStatesPtr, coeffsPtr, count, new InTheMeantime());
  // }

  // proc finish() {
  //   const atomicPtr = isEOF;
  //   on Locales[destLocaleIdx] {
  //     atomicStoreBool(atomicPtr, true);
  //   }
  // }
}

class GlobalPtrStore {
  var arr : [LocaleSpace] (c_ptr(_RemoteBuffer(complex(128))),
                           c_ptr(_LocalBuffer(complex(128))),
                           int);

}

var globalPtrStoreNoQueue = new GlobalPtrStore();
// var globalPtrStoreNoQueue : [LocaleSpace] (c_ptr(Basis),
//                                            c_ptr(ConcurrentAccessor(real(64))),
//                                            c_ptr(_RemoteBuffer(complex(128))),
//                                            c_ptr(_LocalBuffer(complex(128))),
//                                            int);

config const kRemoteBufferSize = 10000;
config const kNumTasks = here.maxTaskPar;
config const kNumConsumerTasks = 24;

// extern proc chpl_task_getId(): chpl_taskID_t;

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
    const (remoteBufferPtr, _localBufferPtr, _numChunks) = ptrStore[srcLocaleIdx];
    for taskIdx in 0 ..# numTasks {
      ref remoteBuffer = (remoteBufferPtr + destLocaleIdx * numTasks + taskIdx).deref();
      localBuffers[srcLocaleIdx, taskIdx].isFull = c_ptrTo(remoteBuffer.isFull);
    }
  }
}

proc _offDiagInitRemoteBuffers(numTasks : int, ref remoteBuffers, const ref ptrStore) {
  const srcLocaleIdx = remoteBuffers.locale.id;
  coforall loc in Locales do on loc {
    const destLocaleIdx = loc.id;
    const (_remoteBufferPtr, localBufferPtr, _numChunks) = ptrStore[destLocaleIdx];
    for taskIdx in 0 ..# numTasks {
      ref myLocalBuffer = (localBufferPtr + srcLocaleIdx * numTasks + taskIdx).deref();
      ref remoteBuffer = remoteBuffers[destLocaleIdx, taskIdx];
      remoteBuffer.basisStates = myLocalBuffer.basisStates;
      remoteBuffer.coeffs = myLocalBuffer.coeffs;
      remoteBuffer.size = c_ptrTo(myLocalBuffer.size);
      remoteBuffer.isEmpty = c_ptrTo(myLocalBuffer.isEmpty);
      remoteBuffer.incrementNumProcessed = c_ptrTo(myLocalBuffer.incrementNumProcessed);
    }
  }
}

record ProducerTimer {
  var run : stopwatch;
  var computeOffDiag : stopwatch;
  var radixOneStep : stopwatch;
  var localProcess : LocalProcessTimer;
  var putSize : int;
  var put : stopwatch;
  var fastOn : stopwatch;
  var submit : stopwatch;
  var submitCalls : int;
  var submitIterations : int;
  var overflowCount : int;
}

record Producer {
  type eltType;
  var _taskIdx : int;
  var numChunks : int;
  var numProducerTasks : int;

  var batchedOperator : BatchedOperator;
  var basisPtr : c_ptrConst(Basis);
  var accessorPtr : c_ptr(ConcurrentAccessor(eltType));
  var representativesPtr : c_ptrConst(uint(64));
  var xPtr : c_ptrConst(eltType);

  var rangesPtr : c_ptrConst(range(int));
  var moreWorkPtr : c_ptr(atomic bool);
  var currentChunkIdxPtr : c_ptr(atomic int);

  var remoteBuffersPtr : c_ptr(_RemoteBuffer(complex(128)));

  var timer : ProducerTimer;
  // var runTimer : stopwatch;
  // var computeOffDiagTimer : stopwatch;
  // var radixOneStepTimer : stopwatch;
  // var localProcessTimer : LocalProcessTimer;
  // var putSize : int;
  // var putTimer : stopwatch;
  // var fastOnTimer : stopwatch;
  // var submitTimer : stopwatch;

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
    this.basisPtr = c_ptrToConst(this.batchedOperator._matrixPtr.deref().basis);
    this.accessorPtr = c_ptrTo(accessor);
    this.representativesPtr = c_ptrToConst(representatives);
    this.xPtr = c_ptrToConst(x);
    this.rangesPtr = c_ptrToConst(ranges);
    this.moreWorkPtr = c_ptrTo(moreWork);
    this.currentChunkIdxPtr = c_ptrTo(currentChunkIdx);
    this.remoteBuffersPtr = c_ptrTo(remoteBuffers);

    // this.runTimer = new stopwatch();
    // this.computeOffDiagTimer = new stopwatch();
    // this.radixOneStepTimer = new stopwatch();
    // this.localProcessTimer = new LocalProcessTimer();
    // this.putSize = 0;
    // this.putTimer = new stopwatch();
    // this.fastOnTimer = new stopwatch();
    // this.submitTimer = new stopwatch();
    // logDebug("Done creating Producer(", taskIdx, ")...");
  }

  inline proc remoteBuffers(localeIdx : int, taskIdx : int) ref {
    assert(0 <= localeIdx && localeIdx < numLocales);
    assert(0 <= taskIdx && taskIdx < numProducerTasks);
    return remoteBuffersPtr[localeIdx * numProducerTasks + taskIdx];
  }

  proc bandwidth() : real {
    const sentBytes = timer.putSize * (c_sizeof(uint(64)) + c_sizeof(complex(128))):int;
    const sentTime = timer.put.elapsed();
    return 8 * sentBytes / (1024.0 * 1024.0 * 1024.0) / sentTime;
  }

  proc trySubmit(ref remoteBuffer,
                 basisStatesPtr : c_ptr(uint(64)),
                 coeffsPtr : c_ptr(complex(128)),
                 count : int,
                 increment : bool = true) {
    if remoteBuffer.isFull.read() then
      return false;

    timer.put.start();
    remoteBuffer.put(basisStatesPtr, coeffsPtr, count);
    timer.putSize += count;
    timer.put.stop();

    remoteBuffer.isFull.write(true);

    const atomicPtr = remoteBuffer.isEmpty;
    const sizePtr = remoteBuffer.size;
    timer.fastOn.start();
    if increment {
      on Locales[remoteBuffer.destLocaleIdx] {
        sizePtr.deref() = count;
        atomicStoreBool(atomicPtr, false);
      }
    }
    else {
      const incrementAtomicPtr = remoteBuffer.incrementNumProcessed;
      on Locales[remoteBuffer.destLocaleIdx] {
        sizePtr.deref() = count;
        atomicStoreBool(incrementAtomicPtr, false);
        atomicStoreBool(atomicPtr, false);
      }
    }
    timer.fastOn.stop();
    return true;
  }

  proc run() {
    assert(numLocales <= 256);
    timer.run.start();

    var radixOffsets : c_array(int, 257);
    var radixExtraOffsets : [0 ..# numLocales] int;
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
      timer.computeOffDiag.start();
      const r : range(int) = rangesPtr[rangeIdx];
      const (n, basisStatesPtr, coeffsPtr, keysPtr) =
        batchedOperator.computeOffDiag(
          r.size,
          representativesPtr + r.low,
          (xPtr + r.low):c_ptrConst(eltType));
      timer.computeOffDiag.stop();

      timer.radixOneStep.start();
      if numLocales == 1 { // We don't have to reshuffle stuff
        radixOffsets[0] = 0;
        radixOffsets[1] = n;
      }
      else {
        radixOneStep(n, keysPtr, radixOffsets, basisStatesPtr, coeffsPtr);
      }
      timer.radixOneStep.stop();

      timer.submit.start();
      POSIX.memset(c_ptrTo(submitted[0]), 0, numLocales:c_size_t * c_sizeof(bool));
      POSIX.memset(c_ptrTo(radixExtraOffsets[0]), 0, numLocales:c_size_t * c_sizeof(int));
      var remaining = numLocales;
      timer.submit.stop();

      {
        const destLocaleIdx = here.id;
        const k = radixOffsets[destLocaleIdx];
        const n = radixOffsets[destLocaleIdx + 1] - k;
        ref remoteBuffer = remoteBuffers[destLocaleIdx, _taskIdx];

        assert(!remoteBuffer.isFull.read());
        localProcess(basisPtr,
                     accessorPtr,
                     basisStatesPtr + k,
                     (coeffsPtr + k):c_ptrConst(complex(128)),
                     n,
                     timer.localProcess);

        submitted[destLocaleIdx] = true;
        remaining -= 1;
      }

      timer.submit.start();
      timer.submitCalls += 1;
      while remaining > 0 {
        timer.submitIterations += 1;
        for destLocaleIdx in 0 ..# numLocales {
          if submitted[destLocaleIdx] then continue;

          const k = radixOffsets[destLocaleIdx] + radixExtraOffsets[destLocaleIdx];
          const n = radixOffsets[destLocaleIdx + 1] - k;
          ref remoteBuffer = remoteBuffers[destLocaleIdx, _taskIdx];
          const capacity = remoteBuffer.capacity;
          if n > capacity {
            if trySubmit(remoteBuffer, basisStatesPtr + k, coeffsPtr + k,
                         capacity, increment=false) {
              radixExtraOffsets[destLocaleIdx] += capacity;
              timer.overflowCount += 1;
            }
          }
          else {
            if trySubmit(remoteBuffer, basisStatesPtr + k, coeffsPtr + k, n) {
              submitted[destLocaleIdx] = true;
              remaining -= 1;
            }
          }
        }
      }
      timer.submit.stop();
    }

    timer.run.stop();
  }

  proc getTimings() {
    var tree = timingTree("Producer.run", timer.run.elapsed());
    tree.addChild(batchedOperator.getTimings());
    tree.addChild(timer.localProcess.getTimings());
    tree.addChild(timingTree(
      "radixOneStep", timer.radixOneStep.elapsed()
    ));
    tree.addChild(timingTree(
      "submit", timer.submit.elapsed(),
      [ ("extra iterations per submit", timer.submitIterations:real / timer.submitCalls:real)
      , ("number overflows", timer.overflowCount:real)

      , ("put", timer.put.elapsed())
      , ("fastOn", timer.fastOn.elapsed())
      ]
    ));
    tree.getChild("submit").getChild("put").addChild(timingTree(
      "bandwidth (Gb/s)", bandwidth()
    ));
    return tree;
  }
}

config const kShuffle : bool = false;
config const kFactor : int = numLocales; // round(1.75 * numLocales):int;

record Consumer {
  type eltType;

  var _taskIdx : int;
  var numConsumerTasks : int;
  var numProducerTasks : int;

  var slots;

  var basisPtr : c_ptrConst(Basis);
  var accessorPtr : c_ptr(ConcurrentAccessor(eltType));

  var totalNumberChunks : int;
  var numProcessedPtr : c_ptr(atomic int);

  var localBuffersPtr : c_ptr(_LocalBuffer(complex(128)));

  var runTimer : stopwatch;
  var localProcessTimer : LocalProcessTimer;
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

    this.basisPtr = c_ptrToConst(basis);
    this.accessorPtr = c_ptrTo(accessor);
    this.totalNumberChunks = totalNumberChunks;
    this.numProcessedPtr = c_ptrTo(numProcessedChunks);
    this.localBuffersPtr = c_ptrTo(localBuffers[0, 0]);

    this.runTimer = new stopwatch();
    this.localProcessTimer = new LocalProcessTimer();
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
            if localBuffer.incrementNumProcessed.read() then
              numProcessedPtr.deref().add(1);
            else
              localBuffer.incrementNumProcessed.write(true);
            localProcess(basisPtr, accessorPtr,
                         localBuffer.basisStates,
                         localBuffer.coeffs:c_ptrConst(complex(128)),
                         localBuffer.size,
                         localProcessTimer);
            localBuffer.isEmpty.write(true);
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

  proc getTimings() {
    var tree = timingTree("Consumer.run", runTimer.elapsed());
    tree.addChild(localProcessTimer.getTimings());
    tree.addChild(timingTree("fastOn", fastOnTimer.elapsed()));
    return tree;
  }
}


private proc localOffDiagonalNoQueue(matrix : Operator,
                                     const ref x : [] ?eltType,
                                     ref y : [] eltType,
                                     const ref representatives : [] uint(64),
                                     in timings : TimingTree? = nil) {
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
               + remoteBufferSize - 1) / (remoteBufferSize * kFactor),
            10 * numProducerTasks),
        representatives.size);
  // logDebug(matrix.basis, ", ", (newRemoteBuffers.dim(0).size, newRemoteBuffers.dim(1).size),
  //                        ", ", (newLocalBuffers.dim(0).size, newLocalBuffers.dim(1).size));
  // logDebug((c_const_ptrTo(matrix.basis), c_ptrTo(accessor),
  //          c_ptrTo(newRemoteBuffers[0, 0]), c_ptrTo(newLocalBuffers[0, 0]), numChunks));
  // logDebug("before: ", globalPtrStoreNoQueue.arr[here.id]);
  globalPtrStoreNoQueue.arr[here.id] = (c_ptrTo(newRemoteBuffers[0, 0]),
                                        c_ptrTo(newLocalBuffers[0, 0]),
                                        numChunks);
  // logDebug("after: ", globalPtrStoreNoQueue.arr[here.id]);
  allLocalesBarrier.barrier();

  // logDebug("Check #1");
  const ptrStore = globalPtrStoreNoQueue.arr;
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
      const (_remoteBufferPtr, _localBufferPtr, numChunks) = ptrStore[localeIdx];
      totalNumberChunks += numChunks;
    }
  }
  // logDebug("Check #7");

  var numProcessedChunks : atomic int = 0;

  // logDebug("remoteBufferSize: ", remoteBufferSize,
  //          "; numChunks: ", numChunks,
  //          "; batchedOperatorChunkSize: ", batchedOperatorChunkSize);

  // logDebug("numChunks = ", numChunks, ", chunkSize = ",
  //          batchedOperatorChunkSize, ", numTasks = ", numTasks,
  //          ", representatives.size = ", representatives.size,
  //          ", numberOffDiagTerms = ", matrix.numberOffDiagTerms());

  // var producerTimers : [0 ..# numProducerTasks] ProducerTimer;
  // var producerRunTime : [0 ..# numProducerTasks] real;
  // var producerComputeOffDiagTime : [0 ..# numProducerTasks] real;
  // var producerApplyOffDiagTime : [0 ..# numProducerTasks] real;
  // var producerMemcpyTime : [0 ..# numProducerTasks] real;
  // var producerStateInfoTime : [0 ..# numProducerTasks] real;
  // var producerCoeffTime : [0 ..# numProducerTasks] real;
  // var producerKeysTime : [0 ..# numProducerTasks] real;
  // var producerRadixOneStepTime : [0 ..# numProducerTasks] real;
  // var producerLocalProcessTime : [0 ..# numProducerTasks] real;
  // var producerAllocTime : [0 ..# numProducerTasks] real;
  // var producerIndexTime : [0 ..# numProducerTasks] real;
  // var producerAccessTime : [0 ..# numProducerTasks] real;
  // var producerSubmitTime : [0 ..# numProducerTasks] real;
  // var producerFastOnTime : [0 ..# numProducerTasks] real;
  // var producerPutTime : [0 ..# numProducerTasks] real;
  // var producerBandwidth : [0 ..# numProducerTasks] real;

  // var consumerRunTime : [0 ..# numConsumerTasks] real;
  // var consumerLocalProcessTime : [0 ..# numConsumerTasks] real;
  // var consumerAllocTime : [0 ..# numConsumerTasks] real;
  // var consumerIndexTime : [0 ..# numConsumerTasks] real;
  // var consumerAccessTime : [0 ..# numConsumerTasks] real;
  // var consumerFastOnTime : [0 ..# numConsumerTasks] real;
  var producerTimings : list(TimingTree, parSafe=true);
  var consumerTimings : list(TimingTree, parSafe=true);

  // logDebug("Check #8");
  allLocalesBarrier.barrier();
  initializationTimer.stop();

  coforall taskIdx in 0 ..# numTasks with (ref accessor, ref producerTimings, ref consumerTimings) {
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

      if kDisplayTimings then
        producerTimings.pushBack(producer.getTimings());

      //  logDebug(producer.getTimings());
      // producerRunTime[taskIdx] = producer.runTimer.elapsed();
      // producerComputeOffDiagTime[taskIdx] = producer.computeOffDiagTimer.elapsed();
      // producerApplyOffDiagTime[taskIdx] = producer.batchedOperator.applyOffDiagTimer.elapsed();
      // producerMemcpyTime[taskIdx] = producer.batchedOperator.memcpyTimer.elapsed();
      // producerStateInfoTime[taskIdx] = producer.batchedOperator.stateInfoTimer.elapsed();
      // producerCoeffTime[taskIdx] = producer.batchedOperator.coeffTimer.elapsed();
      // producerKeysTime[taskIdx] = producer.batchedOperator.keysTimer.elapsed();
      // producerRadixOneStepTime[taskIdx] = producer.radixOneStepTimer.elapsed();
      // producerLocalProcessTime[taskIdx] = producer.localProcessTimer.total.elapsed();
      // producerAllocTime[taskIdx] = producer.localProcessTimer.allocation.elapsed();
      // producerIndexTime[taskIdx] = producer.localProcessTimer.indexing.elapsed();
      // producerAccessTime[taskIdx] = producer.localProcessTimer.accessing.elapsed();
      // producerSubmitTime[taskIdx] = producer.submitTimer.elapsed();
      // producerFastOnTime[taskIdx] = producer.fastOnTimer.elapsed();
      // producerPutTime[taskIdx] = producer.putTimer.elapsed();
      // producerBandwidth[taskIdx] = producer.bandwidth();
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

      if kDisplayTimings then
        consumerTimings.pushBack(consumer.getTimings());
      // if kDisplayTimings then
      //   logDebug(consumer.getTimings());

      // const i = taskIdx - numProducerTasks;
      // consumerRunTime[i] = consumer.runTimer.elapsed();
      // consumerLocalProcessTime[i] = consumer.localProcessTimer.total.elapsed();
      // consumerAllocTime[i] = consumer.localProcessTimer.allocation.elapsed();
      // consumerIndexTime[i] = consumer.localProcessTimer.indexing.elapsed();
      // consumerAccessTime[i] = consumer.localProcessTimer.accessing.elapsed();
      // consumerFastOnTime[i] = consumer.fastOnTimer.elapsed();
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
  if timings != nil && kDisplayTimings {
    var tree = timingTree("localOffDiagonalNoQueue", totalTimer.elapsed());
    tree.addChild(timingTree("initialization", initializationTimer.elapsed()));
    if !producerTimings.isEmpty() then
      tree.addChild(combineTimingTrees(producerTimings.toArray()));
    if !consumerTimings.isEmpty() then
      tree.addChild(combineTimingTrees(consumerTimings.toArray()));
    timings!.addChild(tree);
  }
    // logDebug("localOffDiagonalNoQueue: ", totalTimer.elapsed(), "\n",
    //          " ├─ ", initializationTimer.elapsed(), " in initialization\n"
             // " ├─ producers: ", meanAndErrString(producerRunTime), "\n",
             // " │   ├─ computeOffDiag: ", meanAndErrString(producerComputeOffDiagTime), "\n",
             // " │   │   ├─ applyOffDiag: ", meanAndErrString(producerApplyOffDiagTime), "\n",
             // " │   │   ├─ memcpy:       ", meanAndErrString(producerMemcpyTime), "\n",
             // " │   │   ├─ stateInfo:    ", meanAndErrString(producerStateInfoTime), "\n",
             // " │   │   ├─ rescale:      ", meanAndErrString(producerCoeffTime), "\n",
             // " │   │   └─ localeIdxOf:  ", meanAndErrString(producerKeysTime), "\n",
             // " │   ├─ radixOneStep:   ", meanAndErrString(producerRadixOneStepTime), "\n",
             // " │   ├─ localProcess:   ", meanAndErrString(producerLocalProcessTime), "\n",
             // " │   │   ├─ allocating: ", meanAndErrString(producerAllocTime), "\n",
             // " │   │   ├─ indexing:   ", meanAndErrString(producerIndexTime), "\n",
             // " │   │   └─ accessing:  ", meanAndErrString(producerAccessTime), "\n",
             // " │   └─ submit:         ", meanAndErrString(producerSubmitTime), "\n",
             // " │       ├─ PUT:    ", meanAndErrString(producerPutTime), "\n",
             // " │       └─ fastOn: ", meanAndErrString(producerFastOnTime), "\n",
             // " └─ consumers: ", meanAndErrString(consumerRunTime), "\n",
             // "     ├─ localProcess: ", meanAndErrString(consumerLocalProcessTime), "\n",
             // "     │   ├─ allocating: ", meanAndErrString(consumerAllocTime), "\n",
             // "     │   ├─ indexing:   ", meanAndErrString(consumerIndexTime), "\n",
             // "     │   └─ accessing:  ", meanAndErrString(consumerAccessTime), "\n",
             // "     └─ fastOn:       ", meanAndErrString(consumerFastOnTime), "\n"
             // "    (bandwidth in GB/s: ", meanAndErrString(producerBandwidth), ")");
             // );
}

proc localMatrixVector(matrix : Operator,
                       const ref x : [] ?eltType,
                       ref y : [] eltType,
                       const ref representatives : [] uint(64),
                       in timings : TimingTree? = nil) {
  assert(matrix.locale == here);
  assert(x.locale == here);
  assert(y.locale == here);
  assert(representatives.locale == here);
  if matrix.numberDiagTerms() > 0 then
    localDiagonal(matrix, x, y, representatives, timings);
  else
    y = 0;
  if matrix.numberOffDiagTerms() > 0 then
    localOffDiagonalNoQueue(matrix, x, y, representatives, timings);
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
    var timings : TimingTree? = nil;
    if kDisplayTimings then
      timings = timingTree("matrixVectorProduct", 0.0);

    var timer = new stopwatch();
    timer.start();
    localMatrixVector(myMatrix, myX, myY, myBasisStates, timings);
    timer.stop();

    if kDisplayTimings {
      timings!.stat = timer.elapsed();
      logDebug(timings!);
    }
  }
}

export proc ls_chpl_matrix_vector_product_f64(matrixPtr : c_ptr(ls_hs_operator), numVectors : c_int,
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

export proc ls_chpl_matrix_vector_product_c128(matrixPtr : c_ptr(ls_hs_operator), numVectors : c_int,
                                               xPtr : c_ptr(complex(128)), yPtr : c_ptr(complex(128))) {
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
