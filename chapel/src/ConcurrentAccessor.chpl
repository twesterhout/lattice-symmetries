module ConcurrentAccessor {

use FFI;
use CTypes;
use ChapelLocks;

// config const concurrentAccessorNumLocks = 5 * here.maxTaskPar;

// A simple wrapper around an array `y` which synchronizes accesses. The only
// operation which is allowed is `y[i] += x` <-> `accessor.localAdd(i, x)`. We use
// multiple locks such that multiple threads have lower chance of trying to
// access the same block at the same time. Doing so is still thread-safe, but
// could negatively impact perforamnce especially on machines with high core
// count.
record ConcurrentAccessor {
  type eltType;
  var _data : c_ptr(chpl__processorAtomicType(eltType));
  var _numElts : int;

  proc init(ref arr : [] ?t)
      where arr.domain.rank == 1 && arr.domain.strides == strideKind.one {
    if arr.locale != here then
      halt("ConcurrentAccessor can only be constructed around a local array");
    this.eltType = t;
    // We go via c_void_ptr, because otherwise chpl issues a warning that we're
    // casting to a non-equivalent pointer type
    this._data = c_ptrTo(arr[arr.domain.low]):c_void_ptr:c_ptr(chpl__processorAtomicType(eltType));
    this._numElts = arr.size;
  }

  inline proc localAdd(i : int, x : eltType) {
    _data[i].add(x, order=memoryOrder.relaxed);
  }
}

}
