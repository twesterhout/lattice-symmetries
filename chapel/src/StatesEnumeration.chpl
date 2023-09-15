module StatesEnumeration {

use CommonParameters;
use FFI;
use ForeignTypes;
use Vector;

use BitOps;
use BlockDist;
use ChplConfig;
use CTypes;
use CyclicDist;
use DynamicIters;
use RangeChunk;
use Time;
import Communication;
import OS.POSIX;

/* Get the next integer with the same Hamming weight.

   Semantically equivalent to
   ```
   const m = popcount(v);
   v += 1;
   while (popcount(v) != m) { v += 1; }
   return v;
   ```
 */
private inline proc nextStateFixedHamming(v: uint(64)): uint(64) {
  const t = v | (v - 1);
  return (t + 1) | (((~t & (t + 1)) - 1) >> (ctz(v) + 1));
}

private inline proc nextStateGeneral(v: uint(64)): uint(64) {
  return v + 1;
}

private inline proc nextState(v: uint(64), param isHammingWeightFixed : bool): uint(64) {
  return if isHammingWeightFixed then nextStateFixedHamming(v) else nextStateGeneral(v);
}

private inline iter _manyNextStateIter(in v: uint(64), bound: uint(64),
                                       param isHammingWeightFixed : bool) {
  while true {
    yield v;
    // If v == bound, incrementing v further is unsafe because it can overflow.
    if v == bound then break;
    v = nextState(v, isHammingWeightFixed);
  }
}
private inline iter manyNextState(v: uint(64), bound: uint(64),
                                  isHammingWeightFixed : bool) {
  assert(v <= bound, "v is greater than bound");
  if isHammingWeightFixed then
    for x in _manyNextStateIter(v, bound, true) do
      yield x;
  else
    for x in _manyNextStateIter(v, bound, false) do
      yield x;
}
private proc manyNextState(v: uint(64), bound: uint(64), ref buffer: [] uint(64),
                           isHammingWeightFixed : bool) : int {
  if buffer.size == 0 then return 0;

  var offset = 0;
  for (x, i) in zip(manyNextState(v, bound, isHammingWeightFixed), 0 ..) {
    buffer[offset] = x;
    offset += 1;
    if offset == buffer.size then break;
  }

  return offset;
}

/* TODO: finish later because the version in Haskell is slow.
private proc computeBinomials(dim : int = 64) {
  assert(dim <= 64);
  var coeff : [0 ..# dim, 0 ..# dim] uint(64);
  coeff[0, 0] = 1;
  for n in 1 ..< dim {
    coeff[n, 0] = 1;
    for k in 1 .. n do
      coeff[n, k] = coeff[n - 1, k - 1] + coeff[n - 1, k];
  }
  return coeff;
}

record BinomialsCache {
  var binomials;

  proc init() {
    this.binomials = computeBinomials();
  }

  proc this(n : int, k : int) {
    if (n <= 0 || k > n) then
      return 0;
    return binomials[n, k];
  }

  proc fixedHammingStateToIndex(in alpha : uint(64)) {
    var i = 0;
    var k = 0;
    while alpha != 0 {
      const c = ctz(alpha);
      alpha = alpha & (alpha - 1);

    }
  }
}
*/

private inline proc unprojectedStateToIndex(basisState : uint(64),
                                            isHammingWeightFixed : bool) : int {
  return if isHammingWeightFixed
           then ls_hs_fixed_hamming_state_to_index(basisState)
           else basisState:int;
}
private inline proc unprojectedIndexToState(stateIndex : int,
                                            hammingWeight : int) : uint(64) {
  return if hammingWeight != -1
           then ls_hs_fixed_hamming_index_to_state(stateIndex, hammingWeight:c_int)
           else stateIndex:uint(64);
}

/* We want to split `r` into `numChunks` non-overlapping ranges but such that for each range
   `chunk` we have that `popcount(chunk.low) == popcount(chunk.high) == popcount(r.low)` if
   `isHammingWeightFixed` was set to `true`.
 */
private proc determineEnumerationRanges(r : range(uint(64)), in numChunks : int,
                                        isHammingWeightFixed : bool) {
  // const ranges = determineEnumerationRanges(globalRange, min(numChunks, globalRange.size),
  //                                           isHammingWeightFixed);

  var timer = new stopwatch();
  timer.start();
  const hammingWeight = if isHammingWeightFixed then popCount(r.low):int else -1;
  const lowIdx = unprojectedStateToIndex(r.low, isHammingWeightFixed);
  const highIdx = unprojectedStateToIndex(r.high, isHammingWeightFixed);
  const totalSize = highIdx - lowIdx + 1;
  if numChunks > totalSize then
    numChunks = totalSize;
  var ranges : [0 ..# numChunks] range(uint(64));
  for (r, i) in zip(chunks(lowIdx .. highIdx, numChunks), 0..) {
    ranges[i] = unprojectedIndexToState(r.low, hammingWeight)
                  .. unprojectedIndexToState(r.high, hammingWeight);
  }
  timer.stop();
  if kDisplayTimings then
    logDebug("determineEnumerationRanges(", r, ") took ", timer.elapsed());
  return ranges;
}
private proc determineEnumerationRanges(r : range(uint(64)), numChunks : int, basis : Basis) {
  if basis.isSpinBasis() ||
     basis.isSpinlessFermionicBasis() ||
     (basis.isSpinfulFermionicBasis() && basis.numberUp() == -1) {
    return determineEnumerationRanges(r, numChunks, basis.isHammingWeightFixed());
  }
  if basis.isSpinfulFermionicBasis() {
    // Both the number of ↑ and ↓ spins are set.
    const numberSites = basis.numberSites();
    const mask = (1 << numberSites) - 1; // isolate the lower numberSites bits
    const numberUp = basis.numberUp();
    const numberDown = basis.numberParticles() - numberUp;

    const minA = r.low & mask;
    const maxA = r.high & mask;
    const minB = (r.low >> numberSites) & mask;
    const maxB = (r.high >> numberSites) & mask;
    assert(popCount(minA) == numberUp && popCount(maxA) == numberUp &&
             popCount(minB) == numberDown && popCount(maxB) == numberDown,
           "invalid r.low=" + r.low:string + " and r.high=" + r.high:string +
           " for a spinless fermion system with " + numberUp:string +
           " spins up and " + numberDown:string + " spins down");

    // TODO: This algorithm will be inefficient when we have only 1 spin down and many spins up.
    // In such cases, we should split into subranges using minA .. maxA range.
    const rangesB = determineEnumerationRanges(minB .. maxB, numChunks, true);
    var ranges : [0 ..# rangesB.size] range(uint(64));

    for (dest, rB) in zip(ranges, rangesB) {
      const low = (rB.low << numberSites) | minA;
      const high = (rB.high << numberSites) | maxA;
      dest = low .. high;
    }
    return ranges;
  }
  halt("unsupported basis type");
}

/* Hash function which we use to map spin configurations to locale indices.

   Typical usage:
   ```chapel
   var localeIndex = (hash64_01(x) % numLocales:uint):int;
   ```
*/
private inline proc hash64_01(in x: uint(64)): uint(64) {
  x = (x ^ (x >> 30)) * (0xbf58476d1ce4e5b9:uint(64));
  x = (x ^ (x >> 27)) * (0x94d049bb133111eb:uint(64));
  x = x ^ (x >> 31);
  return x;
}

inline proc localeIdxOf(basisState : uint(64)) : int 
    where CHPL_COMM == "" {
  return 0;
}
inline proc localeIdxOf(basisState : uint(64)) : int 
    where CHPL_COMM != "" {
  return (hash64_01(basisState) % numLocales:uint):int;
}

private proc _enumStatesComputeMasksAndCounts(const ref states, ref outMasks) {
  const totalCount = states.size;
  var counts : [0 ..# numLocales] int;

  outMasks.resize(totalCount);
  if numLocales == 1 { // all data belongs to locale 0
    if totalCount > 0 then
      POSIX.memset(outMasks.data, 0, totalCount:c_size_t * c_sizeof(outMasks.eltType));
    counts[0] = totalCount;
  }
  else {
    foreach i in 0 ..# totalCount {
      const key = localeIdxOf(states[i]);
      outMasks[i] = key:uint(8);
      counts[key] += 1;
    }
  }
  return counts;
}

private proc assertBoundsHaveSameHammingWeight(r : range(?t), isHammingWeightFixed : bool) {
  assert(!isHammingWeightFixed || popCount(r.low) == popCount(r.high),
    "r.low=" + r.low:string + " and r.high=" + r.high:string
      + " have different Hamming weight: " + popCount(r.low):string
      + " vs. " + popCount(r.high):string);
}


private proc _enumerateStatesProjected(r : range(uint(64)), const ref basis : Basis,
                                       ref outStates : Vector(uint(64))) {
  if r.size == 0 then return;
  const isHammingWeightFixed = basis.isHammingWeightFixed();
  assertBoundsHaveSameHammingWeight(r, isHammingWeightFixed);

  var buffer: [0 ..# kIsRepresentativeBatchSize] uint(64);
  var flags: [0 ..# kIsRepresentativeBatchSize] uint(8);
  var norms: [0 ..# kIsRepresentativeBatchSize] real(64);
  var lower = r.low;
  const upper = r.high;

  while true {
    const written = manyNextState(lower, upper, buffer, isHammingWeightFixed);
    isRepresentative(basis, buffer[0 ..# written], flags[0 ..# written], norms[0 ..# written]);

    for i in 0 ..# written do
      if flags[i]:bool && norms[i] > 0 then
        outStates.pushBack(buffer[i]);

    const last = buffer[written - 1];
    if last == upper then break;

    if isHammingWeightFixed { lower = nextState(last, true); }
    else { lower = nextState(last, false); }
  }
}
private proc _enumerateStatesUnprojected(r : range(uint(64)), const ref basis : Basis,
                                         ref outStates : Vector(uint(64))) {
  const isHammingWeightFixed = basis.isHammingWeightFixed();
  const hasSpinInversionSymmetry = basis.hasSpinInversionSymmetry();
  assertBoundsHaveSameHammingWeight(r, isHammingWeightFixed);
  var low = r.low;
  var high = r.high;
  if hasSpinInversionSymmetry {
    const numberSites = basis.numberSites();
    const mask = (1:uint(64) << numberSites) - 1; // isolate the lower numberSites bits
    high = min(high, high ^ mask);
    assert(popCount(high) == popCount(low));
  }
  const lowIdx = unprojectedStateToIndex(low, isHammingWeightFixed);
  const highIdx = unprojectedStateToIndex(high, isHammingWeightFixed);
  const totalCount = highIdx - lowIdx + 1;
  if totalCount > 0 {
    outStates.resize(totalCount);
    manyNextState(low, high, outStates._arr, isHammingWeightFixed);
  }
}
private proc _enumerateStatesUnprojectedSpinfulFermion(r : range(uint(64)),
                                                       const ref basis : Basis,
                                                       ref outStates : Vector(uint(64))) {
  // logDebug("_enumerateStatesUnprojectedSpinfulFermion ...");
  assert(basis.isSpinfulFermionicBasis());
  const numberSites = basis.numberSites();
  const mask = (1 << numberSites) - 1; // isolate the lower numberSites bits
  const numberUp = basis.numberUp();
  const numberDown = basis.numberParticles() - numberUp;
  assert(numberUp >= 0 && numberDown >= 0);
  // logDebug("numberSites = ", numberSites, ", numberUp = ", numberUp, ", numberDown = ", numberDown);

  const minA = r.low & mask;
  const maxA = r.high & mask;
  // logDebug("minA = ", minA, ", maxA = ", maxA);
  assert(popCount(minA) == numberUp && popCount(maxA) == numberUp);
  const countA =
    unprojectedStateToIndex(maxA, isHammingWeightFixed=true) -
      unprojectedStateToIndex(minA, isHammingWeightFixed=true) + 1;

  const minB = (r.low >> numberSites) & mask;
  const maxB = (r.high >> numberSites) & mask;
  assert(popCount(minB) == numberDown && popCount(maxB) == numberDown);
  const countB =
    unprojectedStateToIndex(maxB, isHammingWeightFixed=true) -
      unprojectedStateToIndex(minB, isHammingWeightFixed=true) + 1;

  outStates.resize(countA * countB);
  var offset = 0;
  for xB in manyNextState(minB, maxB, isHammingWeightFixed=true) {
    for xA in manyNextState(minA, maxA, isHammingWeightFixed=true) {
      const x = (xB << numberSites) | xA;
      outStates[offset] = x;
      offset += 1;
    }
  }
}

private proc _enumerateStates(r : range(uint(64)), const ref basis : Basis,
                              ref outStates : Vector(uint(64))) {
  if basis.isSpinBasis() {
    if basis.hasPermutationSymmetries() then
      _enumerateStatesProjected(r, basis, outStates);
    else
      _enumerateStatesUnprojected(r, basis, outStates);
  }
  else if basis.isSpinfulFermionicBasis() {
    if basis.numberUp() == -1 then
      _enumerateStatesUnprojected(r, basis, outStates);
    else
      _enumerateStatesUnprojectedSpinfulFermion(r, basis, outStates);
  }
  else {
    assert(basis.isSpinlessFermionicBasis());
    _enumerateStatesUnprojected(r, basis, outStates);
  }
}

proc prefixSum(arr : [] ?eltType) {
  var sums : [0 ..# arr.size] eltType;
  var total : eltType = 0;
  for i in sums.domain {
    sums[i] = total;
    total += arr[i];
  }
  return sums;
}

proc permuteBasedOnMasks(arrSize : int, masks : c_ptrConst(?maskType), arr : c_ptrConst(?eltType),
                         counts : [] int, destOffsets : [] int,
                         destPtrs : [] c_ptr(eltType)) {
  var offsets : [0 ..# numLocales] int = prefixSum(counts);
  var src : [0 ..# arrSize] eltType = noinit;
  for i in 0 ..# arrSize {
    const key = masks[i]:int;
    src[offsets[key]] = arr[i];
    offsets[key] += 1;
  }

  var copyTimer = new stopwatch();
  copyTimer.start();
  var i = 0;
  for localeIdx in 0 ..# numLocales {
    if counts[localeIdx] > 0 {
      const srcPtr = c_ptrTo(src[i]);
      const destPtr = destPtrs[localeIdx] + destOffsets[localeIdx];
      const destSize = counts[localeIdx]:c_size_t * c_sizeof(eltType);
      PUT(srcPtr, localeIdx, destPtr, destSize);
      i += counts[localeIdx];
    }
  }
  copyTimer.stop();

  return copyTimer.elapsed();
}

inline proc _enumStatesMakeBuckets(numChunks : int) {
  const dom = {0 ..# numChunks} dmapped Cyclic(startIdx=0);
  var buckets : [dom] (Vector(uint(64)), Vector(uint(8)));
  return buckets;
}

proc _enumStatesComputeCounts(ref buckets,
                              const ref ranges : [] range(uint(64)),
                              const ref basis : Basis) {
  assert(here.id == 0);

  const numChunks = ranges.size;
  var counts : [0 ..# numChunks, 0 ..# numLocales] int;
  const countsPtr = c_ptrTo(counts[0, 0]);

  // const ref serializedBasis = basis.json_repr;
  coforall loc in Locales do on loc {
    var outerTimer = new stopwatch();
    outerTimer.start();

    const myBasis = basis;
    // new Basis(serializedBasis);
    const myRanges : [0 ..# numChunks] range(uint(64)) = ranges;
    const mySubdomain = buckets.localSubdomain();

    var enumerateStatesTime : atomic real = 0;
    var computeMasksAndCountsTime : atomic real = 0;
    var putsTime : atomic real = 0;

    forall chunkIdx in dynamic(mySubdomain, chunkSize=1) {
      ref (outStates, outMasks) = buckets.localAccess(chunkIdx);
      var timer = new stopwatch();
      timer.start();
      // This is the actual computation!
      _enumerateStates(myRanges[chunkIdx], myBasis, outStates);
      timer.stop();
      enumerateStatesTime.add(timer.elapsed(), memoryOrder.relaxed);

      timer.clear();
      timer.start();
      const myCounts = _enumStatesComputeMasksAndCounts(outStates, outMasks);
      timer.stop();
      computeMasksAndCountsTime.add(timer.elapsed(), memoryOrder.relaxed);

      timer.clear();
      timer.start();
      // const myCounts = _enumerateStatesUnprojectedNew(myRanges[chunkIdx], myBasis,
      //                                                 outStates, outMasks);
      // Copy computed chunks back to locale 0. We use low-level PUT, because
      // Chapel generates bad code otherwise
      const myCountsPtr = c_ptrToConst(myCounts[0]);
      const destPtr = countsPtr + chunkIdx * numLocales;
      const destSize = numLocales:c_size_t * c_sizeof(int);
      PUT(myCountsPtr, 0, destPtr, destSize);

      timer.stop();
      putsTime.add(timer.elapsed(), memoryOrder.relaxed);
    }

    outerTimer.stop();
    if kDisplayTimings then
      logDebug("_enumStatesComputeCounts: ", outerTimer.elapsed(), "\n",
               "  ├─ ", enumerateStatesTime, " in _enumerateStates\n",
               "  ├─ ", computeMasksAndCountsTime, " in _enumStatesComputeMasksAndCounts\n",
               "  └─ ", putsTime, " in PUT\n",
               "     (parallel speedup: ",
                        (enumerateStatesTime.read() + computeMasksAndCountsTime.read()
                          + putsTime.read()) / outerTimer.elapsed(),
                     ")");
  }

  return counts;
}

proc _enumStatesCountsToOffsets(counts, out _totalCounts) {
  const numChunks = counts.shape[0];

  var offsets : [0 ..# numChunks, 0 ..# numLocales] int;
  var totalCounts : [0 ..# numLocales] int;
  foreach localeIdx in 0 ..# numLocales {
    var total = 0;
    for chunkIdx in 0 ..# numChunks {
      offsets[chunkIdx, localeIdx] = total;
      total += counts[chunkIdx, localeIdx];
    }
    totalCounts[localeIdx] = total;
  }

  _totalCounts = totalCounts;
  return offsets;
}

proc _enumStatesMakeMasksOffsets(counts) {
  const numChunks = counts.shape[0];
  var masksOffsets : [0 ..# numChunks + 1] int;
  var total = 0;
  for chunkIdx in 0 ..# numChunks {
    masksOffsets[chunkIdx] = total;
    foreach k in 0 ..# counts.shape[1] do
      total += counts[chunkIdx, k];
    // total += (+ reduce counts[chunkIdx, ..]);
  }
  masksOffsets[numChunks] = total;
  return masksOffsets;
}

proc _enumStatesMakeMasks(counts, totalCounts) {
  const numMasks = + reduce totalCounts;
  const masksBox = {0 ..# numMasks};
  const masksDom = masksBox dmapped Block(boundingBox=masksBox);
  var masks : [masksDom] uint(8);
  return masks;
}

proc _enumStatesPrepareDescriptors(const ref masks, counts, totalCounts) {
  const numChunks = counts.shape[0];
  const masksOffsets = _enumStatesMakeMasksOffsets(counts);
  const numMasks = masksOffsets[numChunks];
  const ref masksDom = masks.domain;

  var masksPtrs : [0 ..# numLocales] c_ptr(uint(8));
  const mainLocaleIdx = here.id;
  const masksPtrsPtr = c_ptrTo(masksPtrs[0]);
  forall loc in Locales do on loc {
    const ref x = masks.localAccess(masks.localSubdomain(loc).low);
    const ptr = __primitive("_wide_get_addr", x):c_ptr(uint(8));
    Communication.put(dest = masksPtrsPtr + loc.id,
                      src = c_addrOfConst(ptr),
                      destLocID = mainLocaleIdx,
                      numBytes = c_sizeof(c_ptr(uint(8))));
  }

  var masksDescriptors : [0 ..# numChunks] (int, c_ptr(uint(8)));
  for i in 0 ..# numChunks {
    const firstIndex = masksOffsets[i];
    const lastIndex = if i < numChunks - 1
                        then masksOffsets[i + 1] - 1
                        else numMasks - 1;
    const loc = masksDom.distribution.dsiIndexToLocale(firstIndex);
    // The range firstIndex .. lastIndex crosses the locale boundary in masks
    // array.
    if loc != masksDom.distribution.dsiIndexToLocale(lastIndex) then
      masksDescriptors[i] = (firstIndex, nil);
    else {
      const ptrOffset = firstIndex - masksDom.localSubdomain(loc).low;
      masksDescriptors[i] = (loc.id, masksPtrs[loc.id] + ptrOffset);
    }
  }
  return masksDescriptors;
}

proc _enumStatesDistribute(const ref buckets, ref masks,
                           const ref counts,
                           const ref offsets,
                           const ref basisStatesPtrs,
                           const ref masksDescriptors) {
  var copyTimes : [0 ..# numLocales] real;
  var maskCopyTimes : [0 ..# numLocales] real;

  var distributeTimer = new stopwatch();
  distributeTimer.start();

  const basisStatesPtrsPtr = c_ptrToConst(basisStatesPtrs);
  const mainLocaleIdx = here.id;
  const numChunks = counts.shape[0];
  coforall loc in Locales do on loc {
    const mySubdomain = buckets.localSubdomain();
    const myCounts : [0 ..# numChunks, 0 ..# numLocales] int = counts;
    const myOffsets : [0 ..# numChunks, 0 ..# numLocales] int = offsets;
    const myMasksDescriptors : [0 ..# numChunks] (int, c_ptr(uint(8))) = masksDescriptors;
    // Simple assignment fails with a segmentation fault when compiling with
    // CHPL_COMM=none... no idea why, but the following is a workaround
    var myDestPtrs : [0 ..# numLocales] c_ptr(uint(64)) = noinit;
    GET(c_ptrTo(myDestPtrs), mainLocaleIdx, basisStatesPtrsPtr,
        numLocales:c_size_t * c_sizeof(myDestPtrs.eltType));

    var myCopyTime : atomic real;
    var myMaskCopyTime : atomic real;
    forall bucketIdx in mySubdomain {
      const ref (myStates, myMasks) = buckets.localAccess(bucketIdx);

      if myStates.size > 0 {
        // Distributing states
        const copyTime =
          permuteBasedOnMasks(myStates.size,
                              c_ptrToConst(myMasks[0]),
                              c_ptrToConst(myStates[0]),
                              myCounts[bucketIdx, ..],
                              myOffsets[bucketIdx, ..],
                              myDestPtrs);
        myCopyTime.add(copyTime, memoryOrder.relaxed);

        // Distributing masks
        var timer = new stopwatch();
        timer.start();
        const (localeOrOffset, targetPtr) = myMasksDescriptors[bucketIdx];
        if targetPtr != nil {
          const targetLocaleIdx = localeOrOffset;
          PUT(c_ptrToConst(myMasks._arr[0]), targetLocaleIdx, targetPtr,
              myMasks.size:c_size_t * c_sizeof(uint(8)));
        }
        else {
          const targetOffset = localeOrOffset;
          masks[targetOffset ..# myMasks.size] = myMasks.toArray();
        }
        timer.stop();
        myMaskCopyTime.add(timer.elapsed(), memoryOrder.relaxed);
      }
    }

    copyTimes[loc.id] = myCopyTime.read();
    maskCopyTimes[loc.id] = myMaskCopyTime.read();
  }
  distributeTimer.stop();

  return (distributeTimer.elapsed(), copyTimes, maskCopyTimes);
}

proc enumerateStates(ranges : [] range(uint(64)), const ref basis : Basis, out _masks) {
  var timer = new stopwatch();
  timer.start();

  // We distribute ranges among locales using Cyclic distribution to ensure
  // an even workload. For each range, a vector of basis states and a vector of
  // masks is computed. Masks indicate on which locale a basis state should live.
  var bucketsTimer = new stopwatch();
  bucketsTimer.start();
  const numChunks = ranges.size;
  var buckets = _enumStatesMakeBuckets(numChunks);
  bucketsTimer.stop();

  // How many states coming from a certain chunk live on a certain locale.
  // Each chunk computes its own row in parallel and then does a remote PUT here.
  var countsTimer = new stopwatch();
  countsTimer.start();
  const counts = _enumStatesComputeCounts(buckets, ranges, basis);
  countsTimer.stop();

  // Transform counts into offsets such that each task known where to copy data.
  // Total counts tell us how much memory we have to allocate on each locale.
  const totalCounts;
  const offsets = _enumStatesCountsToOffsets(counts, totalCounts);

  // Allocate space for states
  var allocateTimer = new stopwatch();
  allocateTimer.start();
  var basisStates = new BlockVector(uint(64), totalCounts, distribute=true);
  allocateTimer.stop();
  // logDebug("634: the following is going to fail");
  // Simple assignment fails with a segmentation fault when compiling with
  // CHPL_COMM=none. This is a workaround
  var basisStatesPtrs : [0 ..# numLocales] c_ptr(uint(64)) = noinit;
  // if enableSegFault then
  //   basisStatesPtrs = basisStates._dataPtrs;
  // else
  POSIX.memcpy(c_ptrTo(basisStatesPtrs), c_ptrToConst(basisStates._dataPtrs),
           numLocales:c_size_t * c_sizeof(c_ptr(uint(64))));
  // logDebug("634: nope it didn't");

  // Allocate space for masks
  var masksTimer = new stopwatch();
  masksTimer.start();
  var masks = _enumStatesMakeMasks(counts, totalCounts);
  masksTimer.stop();
  var descriptorsTimer = new stopwatch();
  descriptorsTimer.start();
  var masksDescriptors =
    _enumStatesPrepareDescriptors(masks, counts, totalCounts);
  descriptorsTimer.stop();

  // Distribute buckets to basisStates and maks
  const (distributeTime, copyTimes, maskCopyTimes) =
    _enumStatesDistribute(buckets, masks, counts, offsets, basisStatesPtrs, masksDescriptors);

  timer.stop();
  if kDisplayTimings then
    logDebug("enumerateStates: ", timer.elapsed(), "\n",
             "  ├─ ", countsTimer.elapsed(), " in _enumStatesComputeCounts\n",
             "  ├─ ", bucketsTimer.elapsed(), " in _enumStatesMakeBuckets\n",
             "  ├─ ", allocateTimer.elapsed(), " in allocation\n",
             "  ├─ ", masksTimer.elapsed(), " in _enumStatesMakeMasks\n",
             "  ├─ ", descriptorsTimer.elapsed(), " in _enumStatesPrepareDescriptors\n",
             "  └─ ", distributeTime, " shuffling stuff around\n",
             "      ├─ ", copyTimes, " copying states\n",
             "      └─ ", maskCopyTimes, " copying masks");
  _masks = masks;
  return basisStates;
}
proc enumerateStates(globalRange : range(uint(64)), basis : Basis, out masks,
                     in numChunks : int = kEnumerateStatesNumChunks) {
  // If we don't have to project, the enumeration ranges can be split equally,
  // so there's no need to use more chunks than there are cores
  if !basis.requiresProjection() then
    numChunks = numLocales * here.maxTaskPar;
  const ranges = determineEnumerationRanges(globalRange, min(numChunks, globalRange.size), basis);
  return enumerateStates(ranges, basis, masks);
}
proc enumerateStates(const ref basis : Basis, out masks,
                     numChunks : int = kEnumerateStatesNumChunks) {
  const lower = basis.minStateEstimate();
  const upper = basis.maxStateEstimate();
  // logDebug("computed lower=", lower, ", upper=", upper);
  return enumerateStates(lower .. upper, basis, masks, numChunks);
}


export proc ls_chpl_enumerate_representatives(p : c_ptr(ls_hs_basis),
                                              lower : uint(64),
                                              upper : uint(64),
                                              dest : c_ptr(chpl_external_array)) {
  // logDebug("ls_chpl_enumerate_representatives ...");
  const basis = new Basis(p, owning=false);
  // var rs = localEnumerateRepresentatives(basis, lower, upper);
  const masks;
  // logDebug("calling enumerateStates(", lower, ", ", upper, ") ...");
  var rs = enumerateStates(basis, masks);
  // logDebug("enumerateStates done!");
  // ref v = rs[here];
  var v = rs[here];
  // writeln(v.type:string);
  // writeln(getExternalArrayType(v):string);
  // writeln(v._value.isDefaultRectangular():string);
  dest.deref() = convertToExternalArray(v);
}

}
