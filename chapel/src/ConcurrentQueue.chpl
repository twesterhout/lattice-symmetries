module ConcurrentQueue {

  // private use BitOps;
  // private use ChapelLocks;
  private use CTypes;

  require "ConcurrentQueue.h";

  extern record ConcurrentQueue {};

  extern proc ls_chpl_create_ConcurrentQueue(capacity : uint(32), nullElement : uint(64)) : c_ptr(ConcurrentQueue);

  extern proc ls_chpl_destroy_ConcurrentQueue(queue : c_ptr(ConcurrentQueue));

  extern proc ls_chpl_ConcurrentQueue_size(queue : c_ptr(ConcurrentQueue)) : int(64);

  pragma "fast-on safe extern function"
  extern proc ls_chpl_ConcurrentQueue_try_push(queue : c_ptr(ConcurrentQueue), value : uint(64)) : bool;

  pragma "fast-on safe extern function"
  extern proc ls_chpl_ConcurrentQueue_push(queue : c_ptr(ConcurrentQueue), value : uint(64));

  pragma "fast-on safe extern function"
  extern proc ls_chpl_ConcurrentQueue_try_pop(queue : c_ptr(ConcurrentQueue), value : c_ptr(uint(64))) : bool;

  pragma "fast-on safe extern function"
  extern proc ls_chpl_ConcurrentQueue_pop(queue : c_ptr(ConcurrentQueue), value : c_ptr(uint(64)));

  record ConcurrentQueueWrapper {
    type eltType;
    var queuePtr : c_ptr(ConcurrentQueue);
    var localeId : int;
    var owning : bool;

    proc init(type eltType, queuePtr : c_ptr(ConcurrentQueue) = nil, localeId : int = -1, owning : bool = false) {
      this.eltType = eltType;
      this.queuePtr = queuePtr;
      this.localeId = localeId;
      this.owning = owning;
    }
    proc init(capacity : int, nullElement : ?eltType) {
      this.eltType = eltType;
      this.queuePtr = ls_chpl_create_ConcurrentQueue(capacity.safeCast(uint(32)), encodeAsWord(nullElement));
      this.localeId = here.id;
      this.owning = true;
      // init this;
      // assert(queuePtr != nil);
    }
    proc deinit() {
      if owning then
        on Locales[localeId] do
          ls_chpl_destroy_ConcurrentQueue(queuePtr);
    }

    proc tryPush(value : eltType) : bool {
      const rvfValue : uint(64) = encodeAsWord(value);
      const rvfQueue = queuePtr;
      var success : bool;
      if here.id == localeId then
        success = ls_chpl_ConcurrentQueue_try_push(rvfQueue, rvfValue);
      else
        on Locales[localeId] do
          success = ls_chpl_ConcurrentQueue_try_push(rvfQueue, rvfValue);
      return success;
    }

    proc push(value : eltType) {
      const rvfValue : uint(64) = encodeAsWord(value);
      const rvfQueue = queuePtr;
      if here.id == localeId then
        ls_chpl_ConcurrentQueue_push(rvfQueue, rvfValue);
      else
        on Locales[localeId] do
          ls_chpl_ConcurrentQueue_push(rvfQueue, rvfValue);
    }

    proc tryPop(ref value) : bool {
      if here.id != localeId then halt("tryPop should be executed locally");

      var encodedValue : uint(64);
      const success = ls_chpl_ConcurrentQueue_try_pop(queuePtr, c_ptrTo(encodedValue));
      if success then
        decodeFromWord(encodedValue, value);
      return success;
    }

    proc pop() : eltType {
      if here.id != localeId then halt("pop should be executed locally");

      var encodedValue : uint(64);
      ls_chpl_ConcurrentQueue_pop(queuePtr, c_ptrTo(encodedValue));
      var value : eltType;
      decodeFromWord(encodedValue, value);
      return value;
    }

    proc size() {
      if here.id != localeId then halt("size should be executed locally");

      return ls_chpl_ConcurrentQueue_size(queuePtr);
    }
  }

  /*
  proc nextPowerOfTwo(x : int(?n)) {
    assert(x > 0, "expected a positive number");
    return
      if x == 1
        then 1
        else 1 << (n - clz(x - 1));
  }

  class Queue {
    type eltType;

    var head : chpl__processorAtomicType(int);
    var tail : chpl__processorAtomicType(int);
    var dom;
    var buffer : [dom] chpl__processorAtomicType(eltType);
    var capacity : int;
    var mask : int;
    var nullElement : eltType;

    proc init(in capacity : int, nullElement : ?eltType) {
      assert(capacity > 0, "expected a positive capacity");
      this.eltType = eltType;
      capacity = nextPowerOfTwo(capacity);
      this.dom = {0 ..# capacity};
      this.capacity = capacity;
      this.mask = capacity - 1;
      this.nullElement = nullElement;
      init this;

      head.write(0);
      tail.write(0);
      for i in dom do
        buffer[i].write(nullElement);
    }

proc push(value : eltType) : bool {
  var localHead : int;
  do {
    // Writers write at 'head' and readers read at 'tail'
    //      -+---+---+---+---+---+-
    //       |   |   |   |   |   |
    //      -+---+---+---+---+---+-
    //         ^           ^
    //        tail        head
    localHead = head.read();
    const localTail = tail.read();

    // The queue is full
    // assert(localHead - localTail <= capacity, "buffer overflow?");
    if localHead - localTail == capacity then
      return false;

  } while buffer[localHead & mask].read() != nullElement || !head.compareAndSwap(localHead, localHead + 1); // buffer[localHead & mask].read() != nullElement || !head.compareAndSwap(localHead, localHead + 1);
  // - buffer[localHead & mask] != nullElement means that a reader has
  //   incremented tail, but hasn't yet had the chance to reset the value.
  // - head.compareAndSwap(localHead, localHead + 1) returns true if the store of localHead + 1 was successful.
  //   If it returns false, it means that some other thread beat us to it, so we restart the loop.

  const old = buffer[localHead & mask].exchange(value);
  // assert(old == nullElement);
  return true;
}

    proc pop() : (bool, eltType) {
      var localTail : int;
      do {
        localTail = tail.read();
        const localHead = head.read();

        // The queue is empty
        if localTail >= localHead then
          return (false, nullElement);

      } while buffer[localTail & mask].read() == nullElement || !tail.compareAndSwap(localTail, localTail + 1);
      // - buffer[localTail & mask] == nullElement means that a writer has
      //   incremented head, but hasn't yet had the chance to reset the value.
      // - tail.compareAndSwap(localTail, localTail + 1) returns true if the store of localTail + 1 was successful.
      //   If it returns false, it means that some other thread beat us to it, so we restart the loop.

      localTail = localTail & mask;
      const old = buffer[localTail].exchange(nullElement);
      assert(old != nullElement);
      return (true, old);
    }
  }
  */

}
