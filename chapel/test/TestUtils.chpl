use Utils; // initRuntime / deinitRuntime
use ForeignTypes;

use CTypes;
import OS.POSIX;

proc test_prefixSum() {
  assert(prefixSum([1, 2, 3, 1, 5]).equals([0, 1, 3, 6, 7]));
  assert(prefixSum([1, 2, 3, 1, 5], inclusive=true).equals([0, 1, 3, 6, 7, 12]));
  var emptyArr : [0 .. -1] int;
  assert(prefixSum(emptyArr, inclusive=false).equals(emptyArr));
  assert(prefixSum(emptyArr, inclusive=true).equals([0]));

  var arr1 : [0 ..# 2, 0 ..# 3] int;
  arr1[0, ..] = [1, 2, 3];
  arr1[1, ..] = [4, 5, 6];
  var expected1 : [0 ..# 2, 0 ..# 3] int;
  expected1[0, ..] = [0, 0, 0];
  expected1[1, ..] = [1, 2, 3];
  assert(prefixSum(arr1, dim=0).equals(expected1));
  var expected2 : [0 ..# 2, 0 ..# 3] int;
  expected2[0, ..] = [0, 1, 3];
  expected2[1, ..] = [0, 4, 9];
  assert(prefixSum(arr1, dim=1).equals(expected2));
}

proc test_sum() {
  var arr1 : [0 ..# 2, 0 ..# 3] int;
  arr1[0, ..] = [1, 2, 3];
  arr1[1, ..] = [4, 5, 6];
  assert(sum(arr1, dim=0).equals([5, 7, 9]));
  assert(sum(arr1, dim=1).equals([6, 15]));
}

proc test_radixOneStep() {
  var keys = [1, 2, 3, 1, 1, 3, 2, 5]:uint(8);
  var offsets1 : c_array(int, 257);
  POSIX.memset(offsets1, 0, 257:c_size_t * c_sizeof(int));
  var arr1 = [1, 2, 3, 4, 5, 6, 7, 8];
  var arr2 = [10, 20, 30, 40, 50, 60, 70, 80];
  unstableRadixOneStep(keys.size, c_ptrTo(keys), offsets1, c_ptrTo(arr1));
  assert(keys.equals([1, 1, 1, 2, 2, 3, 3, 5]));
  assert(offsets1[0] == 0);
  assert(offsets1[1] == 0);
  assert(offsets1[2] == 3);
  assert(offsets1[3] == 5);
  assert(offsets1[4] == 7);
  assert(offsets1[5] == 7);
  assert(offsets1[6] == 8);
  for i in 7 .. 256 do
    assert(offsets1[i] == 8);
  // this may break because the partitioning is unstable
  assert(arr1.equals([1, 5, 4, 2, 7, 3, 6, 8]));

  keys = [1, 2, 3, 1, 1, 3, 2, 5]:uint(8);
  POSIX.memset(offsets1, 0, 257:c_size_t * c_sizeof(int));
  arr1 = [1, 2, 3, 4, 5, 6, 7, 8];
  arr2 = [10, 20, 30, 40, 50, 60, 70, 80];
  stableRadixOneStep(keys.size, c_ptrToConst(keys), offsets1, c_ptrTo(arr1), c_ptrTo(arr2));
  assert(offsets1[0] == 0);
  assert(offsets1[1] == 0);
  assert(offsets1[2] == 3);
  assert(offsets1[3] == 5);
  assert(offsets1[4] == 7);
  assert(offsets1[5] == 7);
  assert(offsets1[6] == 8);
  for i in 7 .. 256 do
    assert(offsets1[i] == 8);
  assert(arr1.equals([1, 4, 5, 2, 7, 3, 6, 8]));
  assert(arr2.equals([10, 40, 50, 20, 70, 30, 60, 80]));

  // Test when all keys are the same
  const keys2 = [0, 0, 0, 0, 0]:uint(8);
  var offsets2 : c_array(int, 257);
  POSIX.memset(offsets2, 0, 257:c_size_t * c_sizeof(int));
  var arr3 = [1, 2, 3, 4, 5];
  var arr4 = [10, 20, 30, 40, 50];
  stableRadixOneStep(keys2.size, c_ptrToConst(keys2), offsets2, c_ptrTo(arr3), c_ptrTo(arr4));
  assert(offsets2[0] == 0);
  assert(offsets2[1] == 5);
  for i in 2 .. 256 do
    assert(offsets2[i] == 5);
  assert(arr3.equals([1, 2, 3, 4, 5]));
  assert(arr4.equals([10, 20, 30, 40, 50]));
}

proc test_refcounts() {
  const op = new Operator(
      "{"
    + "\"basis\": { \"number_spins\": 2 },"
    + "\"expression\": { \"expression\": \"2 (σ⁺₀ σ⁻₁ + σ⁺₁ σ⁻₀) + σᶻ₀ σᶻ₁\" }"
    + "}");
}

proc main() {
  initRuntime();
  defer deinitRuntime();

  test_prefixSum();
  test_sum();
  test_radixOneStep();
  test_refcounts();

  return 0;
}
