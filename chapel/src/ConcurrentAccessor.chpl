module ConcurrentAccessor {

private use FFI;
private use CTypes;
private use ChapelLocks;
private use IO.FormattedIO;

// config const concurrentAccessorNumLocks = 5 * here.maxTaskPar;

private proc realType(type t) type {
  if t == real(32) then return t;
  if t == real(64) then return t;
  if t == complex(64) then return real(32);
  if t == complex(128) then return real(64);
  compilerError("ConcurrentAccessor does not support " + t:string);
}

private proc isComplex(type t) param {
  return t == complex(64) || t == complex(128);
}

// A simple wrapper around an array `y` which synchronizes accesses. The only
// operation which is allowed is `y[i] += x` <-> `accessor.localAdd(i, x)`. We use
// multiple locks such that multiple threads have lower chance of trying to
// access the same block at the same time. Doing so is still thread-safe, but
// could negatively impact perforamnce especially on machines with high core
// count.
record ConcurrentAccessor {
  type eltType;
  var _data : c_ptr(chpl__processorAtomicType(realType(eltType)));
  var _numElts : int;

  var _locks : c_array(chpl_LocalSpinlock, 256);

  proc init(ref arr : [] ?t)
      where arr.domain.rank == 1 && arr.domain.strides == strideKind.one {
    if arr.locale != here then
      halt("ConcurrentAccessor can only be constructed around a local array");
    this.eltType = t;
    // We go via c_void_ptr, because otherwise chpl issues a warning that we're
    // casting to a non-equivalent pointer type
    this._data = c_ptrTo(arr[arr.domain.low]):c_ptr(void)
                                             :c_ptr(chpl__processorAtomicType(realType(eltType)));
    this._numElts = arr.size;
  }

  inline proc ref lock(i : int) { _locks[i:int].lock(); }
  inline proc ref unlock(i : int) { _locks[i:int].unlock(); }

  inline proc localAdd(i : int, x : eltType) {
    // assert(0 <= i && i < _numElts,
    //        try! "localAdd: index out of bounds: %i, but _numElts=%i".format(i, _numElts));
    if !isComplex(eltType) then
      _data[i].add(x, order=memoryOrder.relaxed);
    else {
      _data[2 * i].add(x.re, order=memoryOrder.relaxed);
      _data[2 * i + 1].add(x.im, order=memoryOrder.relaxed);
    }
  }
}

}
