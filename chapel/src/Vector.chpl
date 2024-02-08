module Vector {

// TODO: probably shouldn't use it...
use ArrayViewSlice;
use BlockDist;
use CTypes;

use Utils;

record Vector {
  type eltType;
  // _dom specifies the capacity, and _size is the current size
  var _dom : domain(1, idxType=int, strides=strideKind.one);
  var _arr : [_dom] eltType;
  var _size : int;

  proc init(type eltType) {
    this.eltType = eltType;
    this._size = 0;
  }
  proc init(arr : [] ?t) {
    this.eltType = t;
    this._dom = {0 ..# arr.size};
    this._arr = arr;
    this._size = arr.size;
  }

  forwarding _arr only this;

  proc ref reserve(capacity : int) {
    if capacity > _dom.size then
      _dom = {0 ..# capacity};
  }

  proc ref resize(newSize : int) {
    if newSize > _size then
      reserve(newSize);
    _size = newSize;
  }

  inline proc size { return _size; }

  proc ref defaultGrow(factor : real = 1.5) {
    const currentCapacity = _dom.size;
    const newCapacity =
      max(currentCapacity + 1, round(factor * currentCapacity):int);
    reserve(newCapacity);
  }

  inline proc ref pushBack(x : eltType) {
    if _size == _dom.size then
      defaultGrow();
    _arr[_size] = x;
    _size += 1;
  }

  proc ref append(xs : [] eltType) {
    if _size + xs.size > _dom.size then
      reserve(_size + xs.size);
    _arr[_size ..# xs.size] = xs;
    _size += xs.size;
  }
  proc ref append(const ref xs : Vector(eltType)) {
    append(xs._arr[0 ..# xs._size]);
  }

  proc ref shrink() {
    if _size < _dom.size then
      _dom = {0 ..# _size};
  }

  proc ref clear() {
    resize(0);
    shrink();
  }

  inline proc ref data : c_ptr(eltType) { return c_ptrTo(_arr[_dom.low]); }

  pragma "reference to const when const this"
  pragma "fn returns aliasing array"
  proc toArray() {
    pragma "no auto destroy" var d = {0 ..# _size};
    d._value._free_when_no_arrs = true;
    d._value.definedConst = true;
    var a = new unmanaged ArrayViewSliceArr(
        eltType=_arr.eltType,
        _DomPid=d._pid, dom=d._instance,
        _ArrPid=_arr._pid, _ArrInstance=_arr._value);
    d._value.add_arr(a, locking=false, addToList=false);
    return _newArray(a);
  }
}

proc isVector(type x : Vector(?)) param { return true; }
proc isVector(type x) param { return false; }

// proc _getDataPtrs(type eltType, const ref counts, ref data) {
//   const size = counts.size;
//   // logDebug("_getDataPtrs(", eltType:string, ", counts=", counts, ", ", data.domain, ")");
//   // assert(data.size == size);
//   // assert(data.size > 0);
//   // assert(data[0].size > 0);
//   var dataPtrs : [0 ..# size] c_ptr(eltType);
//   if size == 0 then return dataPtrs;
// 
//   const dataPtrsPtr = c_ptrTo(dataPtrs);
//   const localeIdx = here.id;
//   forall (c, i) in zip(counts, 0 ..) {
//     assert(here == c.locale);
//     ref myData = data[i];
//     // if myData.domain.rank == 1 then
//     //   assert(myData.domain.low == 0);
//     const myDataPtr = if myData.size > 0
//                         then c_ptrTo(data[i][myData.domain.low])
//                         else nil;
//     dataPtrs[i] = myDataPtr;
//     // if here.id != localeIdx then
//     //   PUT(myDataPtr, localeIdx, dataPtrsPtr + i, c_sizeof(c_ptr(eltType)));
//     // else
//     //   dataPtrsPtr[i] = myDataPtr;
//   }
//   return dataPtrs;
// }


class LocBlockArr {
  type eltType;
  param rank : int;
  var dom : domain(rank);
  var data : [dom] eltType;

  forwarding data only this;

  proc dataPtr() {
    if dom.size == 0 then return nil;
    return c_ptrTo(data[dom.low]);
  }
}

// proc getBlockPtrs(type eltType, const ref arr, ref ptrs) {
//   // var ptrs : [0 ..# arr.size] c_ptr(eltType);
//   for i in arr.dim(0) {
//     ptrs[i] = arr[i].dataPtr();
//     // ref locBlock = arr[i];
//     // if locBlock.dom.size > 0 {
//     //   ref x = locBlock.data[locBlock.dom.low];
//     //   ptrs[i] = __primitive("_wide_get_addr", x):c_ptr(eltType);
//     // }
//   }
//   // return ptrs;
// }

// proc getOuterDom(numBlocks : int, param distribute : bool) {
//   const box = {0 ..# numBlocks};
//   return if distribute then box dmapped blockDist(box, Locales)
//                        else box;
// }
// 
// iter makeInnerDoms(shape : int ...?n, lastDimSizes : [] int) {
//   for n in lastDimSizes do
//     yield sizeToDomain((...shape), n);
// }

class BlockVector {
  type eltType;
  param rank : int;
  var _outerDom;
  var _locBlocks : [_outerDom] unmanaged LocBlockArr(eltType, rank);
  var _dataPtrs;

  proc finalizeInitialization(counts) {
    forall (locBlock, dataPtr, shape) in zip(_locBlocks, _dataPtrs, counts) {
      locBlock.dom = if isTuple(shape) then sizeToDomain((...shape)) else sizeToDomain(shape);
      dataPtr = locBlock.dataPtr();
    }
  }

  proc deinit() {
    forall locBlock in _locBlocks {
      delete locBlock;
    }
  }

  proc init(type t, param rank : int) {
    this.eltType = t;
    this.rank = rank;
    const dom = LocaleSpace dmapped blockDist(LocaleSpace, Locales);
    this._outerDom = dom;
    this._locBlocks = [i in dom] new unmanaged LocBlockArr(t, rank);
    this._dataPtrs = [i in dom.dim(0)] nil:c_ptr(t);
  }

  proc init(type t, counts : [] int) where counts.domain.rank == 1 {
    assert(counts.size == numLocales);
    init(t, 1);
    finalizeInitialization(counts);
  }
  proc init(type t, counts : [] ?shape) where counts.domain.rank == 1 && isTuple(shape) {
    assert(counts.size == numLocales);
    init(t, shape.size);
    finalizeInitialization(counts);
  }
  proc init(chunks : [] ?t) where chunks.domain.rank == 1 && isVector(t) {
    assert(chunks.size == numLocales);
    init(t.eltType, 1);
    finalizeInitialization([c in chunks] c.size);
    forall (locBlock, chunk) in zip(_locBlocks, chunks) {
      locBlock.data = chunk.toArray();
    }
  }

  inline proc count(i) { return _locBlocks[i].count; }

  inline proc counts { return _locBlocks.count; }

  inline proc _data ref { return _locBlocks; }

  // inline proc getBlockDomain(blockIdx : int)
  //     where innerRank == 1 {
  //   return {0 ..# count(blockIdx)};
  // }
  // inline proc getBlockDomain(blockIdx : int)
  //     where innerRank == 2 {
  //   return {0 ..# _locBlocks[blockIdx].data.shape[0],
  //           0 ..# _locBlocks[blockIdx].count};
  // }

  // pragma "fn returns aliasing array"
  pragma "reference to const when const this"
  proc getBlock(loc : locale = here) ref {
    return _locBlocks[loc.id].data;
    // pragma "no auto destroy" var d = getBlockDomain(blockIdx);
    // d._value._free_when_no_arrs = true;
    // d._value.definedConst = true;
    // var a = new unmanaged ArrayViewSliceArr(
    //     eltType=eltType,
    //     _DomPid=d._pid, dom=d._instance,
    //     _ArrPid=_locBlocks[blockIdx].data._pid,
    //     _ArrInstance=_locBlocks[blockIdx].data._value);
    // d._value.add_arr(a, locking=false, addToList=false);
    // return _newArray(a);
  }

  // pragma "fn returns aliasing array"
  pragma "reference to const when const this"
  inline proc this(loc : locale) ref { return getBlock(loc); }

  pragma "reference to const when const this"
  inline proc this(i : int, j) ref {
    return _locBlocks[i].data[j];
  }

  // pragma "fn returns aliasing array"
  // inline proc this(i : int, j : range(int)) ref where innerRank == 1 {
    // pragma "no auto destroy" var d = {j};
    // d._value._free_when_no_arrs = true;
    // d._value.definedConst = true;
    // var a = new unmanaged ArrayViewSliceArr(
    //     eltType=eltType,
    //     _DomPid=d._pid, dom=d._instance,
    //     _ArrPid=_locBlocks[i].data._pid,
    //     _ArrInstance=_locBlocks[i].data._value);
    // d._value.add_arr(a, locking=false, addToList=false);
    // return _newArray(a);
  // }

  inline proc this(i : int, j...) ref where rank != 1 {
    return _locBlocks[i].data[(...j)];
  }
  // inline proc this(i : int, j, k) ref where innerRank == 2 {
  //   return _locBlocks[i].data[j, k];
  // }
}

proc similar(const ref other : BlockVector(?eltType, ?rank)) {
  return similar(other.eltType, other);
}
proc similar(type t, const ref other : BlockVector(?eltType, ?rank)) {
  if rank == 1 then
    return new BlockVector(t, other.counts);
  else if rank == 2 then
    return new BlockVector(t, other.innerDom.shape[0], other.counts);
  else
    compilerError("unsupported rank: " + rank:string);
}

} // end module Vector
