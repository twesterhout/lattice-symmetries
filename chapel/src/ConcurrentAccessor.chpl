module ConcurrentAccessor {

use FFI;
use CTypes;
use ChapelLocks;

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

  inline proc localAdd(i : int, x : eltType) {
    if !isComplex(eltType) then
      _data[i].add(x, order=memoryOrder.relaxed);
    else {
      _data[2 * i].add(x.re, order=memoryOrder.relaxed);
      _data[2 * i + 1].add(x.im, order=memoryOrder.relaxed);
    }
  }
}

}
