module Vector {

use FFI;
// TODO: probably shouldn't use it...
use ArrayViewSlice;
use BlockDist;
use CTypes;

record Vector {
  type eltType;
  // _dom specifies the capacity, and _size is the current size
  var _dom : domain(1, idxType=int, stridable=false);
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
  // proc this(i : int) {
  //   if i >= _size then
  //     halt("index " + i:string + " is out of bounds for domain {0 ..# " + _size:string + "}");
  //   return _arr[i];
  // }

  proc reserve(capacity : int) {
    if capacity > _dom.size then
      _dom = {0 ..# capacity};
  }

  proc resize(newSize : int) {
    if newSize > _size then
      reserve(newSize);
    _size = newSize;
  }

  inline proc size { return _size; }

  proc defaultGrow(factor : real = 1.5) {
    const currentCapacity = _dom.size;
    const newCapacity =
      max(currentCapacity + 1, round(factor * currentCapacity):int);
    reserve(newCapacity);
  }

  inline proc pushBack(x : eltType) {
    if _size == _dom.size then
      defaultGrow();
    _arr[_size] = x;
    _size += 1;
  }

  proc append(xs : [] eltType) {
    if _size + xs.size > _dom.size then
      reserve(_size + xs.size);
    _arr[_size ..# xs.size] = xs;
    _size += xs.size;
  }
  proc append(const ref xs : Vector(eltType)) {
    append(xs._arr[0 ..# xs._size]);
  }

  proc shrink() {
    if _size < _dom.size then
      _dom = {0 ..# _size};
  }

  inline proc data : c_ptr(eltType) { return c_ptrTo(_arr[_dom.low]); }

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

proc isVector(type x : Vector) param { return true; }
proc isVector(type x) param { return false; }

proc _getDataPtrs(type eltType, const ref counts, ref data) {
  const size = counts.size;
  // logDebug("_getDataPtrs(", eltType:string, ", counts=", counts, ", ", data.domain, ")");
  // assert(data.size == size);
  // assert(data.size > 0);
  // assert(data[0].size > 0);
  var dataPtrs : [0 ..# size] c_ptr(eltType);
  if size == 0 then return dataPtrs;

  const dataPtrsPtr = c_ptrTo(dataPtrs);
  const localeIdx = here.id;
  forall (c, i) in zip(counts, 0 ..) {
    assert(here == c.locale);
    ref myData = data[i];
    // if myData.domain.rank == 1 then
    //   assert(myData.domain.low == 0);
    const myDataPtr = if myData.size > 0
                        then c_ptrTo(data[i][myData.domain.low])
                        else nil;
    dataPtrs[i] = myDataPtr;
    // if here.id != localeIdx then
    //   PUT(myDataPtr, localeIdx, dataPtrsPtr + i, c_sizeof(c_ptr(eltType)));
    // else
    //   dataPtrsPtr[i] = myDataPtr;
  }
  return dataPtrs;
}


class LocBlockVector {
  type eltType;
  param rank : int;
  var dom : domain(rank);
  var data : [dom] eltType;
  var count : int;

  forwarding data only this;
}

proc getBlockPtrs(type eltType, const ref arr, ref ptrs) {
  // var ptrs : [0 ..# arr.size] c_ptr(eltType);
  for i in arr.dim(0) {
    ref locBlock = arr[i];
    if locBlock.dom.size > 0 {
      ref x = locBlock.data[locBlock.dom.low];
      ptrs[i] = __primitive("_wide_get_addr", x):c_ptr(eltType);
    }
  }
  // return ptrs;
}

proc getOuterDom(numBlocks : int, param distribute : bool) {
  const box = {0 ..# numBlocks};
  return if distribute then box dmapped Block(box, Locales)
                       else box;
}

proc getInnerDom(counts : [] int) {
  const maxNumElts = max reduce counts;
  return {0 ..# maxNumElts};
}
proc getInnerDom(batchSize : int, counts : [] int) {
  const maxNumElts = max reduce counts;
  return {0 ..# batchSize, 0 ..# maxNumElts};
}

class BlockVector {
  type eltType;
  param innerRank : int;
  var outerDom;
  var _locBlocks : [outerDom] unmanaged LocBlockVector(eltType, innerRank);
  var _dataPtrs;

  proc finalizeInitialization(innerDom : domain(innerRank), counts) {
    forall (i, n) in zip(outerDom, counts) with (in innerDom) {
      _locBlocks[i].dom = innerDom;
      _locBlocks[i].count = n;
      _dataPtrs[i] = c_ptrTo(_locBlocks[i].data[innerDom.low]);
    }
    // logDebug("this is going to fail...");
    // _dataPtrs = getBlockPtrs(eltType, _locBlocks);
    // getBlockPtrs(eltType, _locBlocks, _dataPtrs);
    // logDebug("hm... nope, it didn't");
  }
  inline proc finalizeInitialization(counts : [] int) {
    finalizeInitialization(getInnerDom(counts), counts);
  }
  inline proc finalizeInitialization(batchSize : int, counts : [] int) {
    finalizeInitialization(getInnerDom(batchSize, counts), counts);
  }

  proc deinit() {
    forall locBlock in _locBlocks {
      delete locBlock;
    }
  }

  proc init(type t, param rank : int, numBlocks : int, param distribute : bool) {
    this.eltType = t;
    this.innerRank = rank;
    const dom = getOuterDom(numBlocks, distribute);
    this.outerDom = dom;
    this._locBlocks = [i in dom] new unmanaged LocBlockVector(t, rank);
    this._dataPtrs = [i in 0 ..# numBlocks] nil:c_ptr(t);
  }

  proc init(type t, counts : [] int, param distribute : bool = true) {
    init(t, 1, counts.size, distribute);
    finalizeInitialization(counts);
  }

  proc init(type t, batchSize : int, counts : [] int, param distribute : bool = true) {
    init(t, 2, counts.size, distribute);
    finalizeInitialization(batchSize, counts);
  }

  proc init(chunks : [] ?t, param distribute : bool = true)
      where chunks.domain.rank == 1 && isVector(t) {
    init(t.eltType, 1, chunks.size, distribute);
    const counts : [0 ..# chunks.size] int = [c in chunks] c.size;
    finalizeInitialization(counts);
    forall (locBlock, chunk) in zip(this._locBlocks, chunks) {
      locBlock.data[0 ..# chunk.size] = chunk.toArray();
    }
  }

  // proc init(chunks : [] ?eltType, counts : [] int)
  //     where !isArray(eltType) {
  //   this.eltType = eltType;
  //   this.innerRank = chunks.rank - 1;
  //   this._outerDom = counts.domain;
  //   this._counts = counts;
  //   const maxNumElts = max reduce _counts;
  //   this._innerDom = {0 ..# maxNumElts};
  //   this._dataPtrsDom = {0 ..# counts.size};
  //   this.complete();
  //   this._dataPtrs = _getDataPtrs(eltType, _counts, _data);
  //   forall (arr, count, i) in zip(this._data, this._counts, 0 ..) {
  //     arr[0 ..# count] = chunks[i, 0 ..# count];
  //   }
  // }

  inline proc numBlocks { return outerDom.size; }
  inline proc count(i) { return _locBlocks[i].count; }

  inline proc counts { return _locBlocks.count; }

  inline proc innerDom { return _locBlocks[0].dom; }

  inline proc _data ref { return _locBlocks; }

  inline proc getBlockDomain(blockIdx : int)
      where innerRank == 1 {
    return {0 ..# count(blockIdx)};
  }
  inline proc getBlockDomain(blockIdx : int)
      where innerRank == 2 {
    return {0 ..# _locBlocks[blockIdx].data.shape[0],
            0 ..# _locBlocks[blockIdx].count};
  }

  pragma "reference to const when const this"
  pragma "fn returns aliasing array"
  proc getBlock(blockIdx : int) {
    pragma "no auto destroy" var d = getBlockDomain(blockIdx);
    d._value._free_when_no_arrs = true;
    d._value.definedConst = true;
    var a = new unmanaged ArrayViewSliceArr(
        eltType=eltType,
        _DomPid=d._pid, dom=d._instance,
        _ArrPid=_locBlocks[blockIdx].data._pid,
        _ArrInstance=_locBlocks[blockIdx].data._value);
    d._value.add_arr(a, locking=false, addToList=false);
    return _newArray(a);
  }

  pragma "reference to const when const this"
  pragma "fn returns aliasing array"
  inline proc this(loc : locale) { return getBlock(loc.id); }

  inline proc this(i : int, j : int) ref where innerRank == 1 {
    return _locBlocks[i].data[j];
  }

  pragma "reference to const when const this"
  pragma "fn returns aliasing array"
  inline proc this(i : int, j : range(int)) ref where innerRank == 1 {
    pragma "no auto destroy" var d = {j};
    d._value._free_when_no_arrs = true;
    d._value.definedConst = true;
    var a = new unmanaged ArrayViewSliceArr(
        eltType=eltType,
        _DomPid=d._pid, dom=d._instance,
        _ArrPid=_locBlocks[i].data._pid,
        _ArrInstance=_locBlocks[i].data._value);
    d._value.add_arr(a, locking=false, addToList=false);
    return _newArray(a);
  }

  inline proc this(i : int, j...) ref where innerRank != 1 {
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
    return new BlockVector(t, other.numBlocks, other.counts);
  else
    compilerError("unsupported rank: " + rank:string);
}

} // end module Vector
