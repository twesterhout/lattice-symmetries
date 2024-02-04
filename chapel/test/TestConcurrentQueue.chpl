use ConcurrentQueue;

use CTypes;
use Random;
use BlockDist;
use CommDiagnostics;


config const kArraySize = 1000;
config const kCapacity = 100;
config const kVerboseComm = false;

proc test_queue() {
  var dom = blockDist.createDomain({0 ..# kArraySize});
  var arr : [dom] uint(64);
  coforall loc in Locales with (ref arr) do on loc {
    ref block = arr[arr.localSubdomain()];

    // var localChunk : [0 ..# block.size] block.eltType;
    // assert (localChunk.isRectangular());
    fillRandom(block, min=1:block.eltType, max=1000:block.eltType, seed=42 + loc.id);

    // block = localChunk;
  }
  const totalCount = arr.size;
  const totalSum = + reduce arr;

  // var queue = new Queue(capacity=100, nullElement=0xFFFFFFFFFFFFFFFF:uint(64));
  const c_queue = ls_chpl_create_ConcurrentQueue(100:uint(32), 0xFFFFFFFFFFFFFFFF:uint(64));
  defer ls_chpl_destroy_ConcurrentQueue(c_queue);

  var count : atomic int;
  count.write(0);
  var sum : atomic arr.eltType;
  sum.write(0);

  if kVerboseComm then startVerboseComm();

  coforall loc in Locales with (const ref arr) do on loc {
    const myQueue = c_queue;

    if loc.id == 0 {

      assert(here.maxTaskPar > 1);
      coforall taskIdx in 0 ..# here.maxTaskPar {
        // First task is still pushing elements to the queue
        if taskIdx == 0 {
          for i in arr.localSubdomain() do
            ls_chpl_ConcurrentQueue_push(myQueue, arr[i]);
        }
        // Other tasks are popping from the queue
        else {
          while count.read() < totalCount {
            var x : uint(64);
            const success = ls_chpl_ConcurrentQueue_try_pop(myQueue, c_addrOf(x));
            if success {
              sum.add(x);
              count.add(1);
            }
          }
        }
      }

    }
    else {
      // Pushing all elements to the queue
      forall i in arr.localSubdomain() {
        const x = arr[i];
        on Locales[0] do
          ls_chpl_ConcurrentQueue_push(myQueue, x);
      }
    }
  }

  if kVerboseComm then stopVerboseComm();

  writeln(totalSum, ", ", sum.read());
  assert(totalSum == sum.read());
}

proc main() {
  test_queue();
}
