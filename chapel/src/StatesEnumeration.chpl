module StatesEnumeration {

private use CommonParameters;
private use FFI;
private use ForeignTypes;
private use Vector;
private use Timing;
private use Utils;

private import Communication;
private import IO.FormattedIO.format;
private import OS.POSIX;
private use BitOps;
private use BlockDist;
private use CTypes;
private use ChplConfig;
private use CyclicDist;
private use DynamicIters;
private use RangeChunk;

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

/* Get the next integer
 */
private inline proc nextStateGeneral(v: uint(64)): uint(64) {
  return v + 1;
}

/* Call nextStateFixedHamming or nextStateGeneral depending on the isHammingWeightFixed
 */
private inline proc nextState(v: uint(64), param isHammingWeightFixed : bool): uint(64) {
  return if isHammingWeightFixed then nextStateFixedHamming(v) else nextStateGeneral(v);
}

/* Iterate on `nextState` until we reach the upper bound `bound`.
 */
private inline iter _manyNextStateIter(in v: uint(64), bound: uint(64),
                                       param isHammingWeightFixed : bool) {
  while true {
    yield v;
    // If v == bound, incrementing v further is unsafe because it can overflow.
    if v >= bound then break;
    v = nextState(v, isHammingWeightFixed);
  }
}
/* Same as `_manyNextStateIter`, but `isHammingWeightFixed` is now a run-time
   parameter instead of a compile-time one.
 */
private inline iter manyNextState(v: uint(64), bound: uint(64),
                                  isHammingWeightFixed : bool) {
  if v > bound then
    halt(try! "v (%i) is greater than bound (%i)".format(v, bound));

  if isHammingWeightFixed then
    for x in _manyNextStateIter(v, bound, true) do
      yield x;
  else
    for x in _manyNextStateIter(v, bound, false) do
      yield x;
}
/* Same as `manyNextStateIter`, but instead of iterating, fill the buffer `buffer`.
   Fills `buffer` until it becomes full or until `bound` is reached.
   Returns the number of elements written to `buffer`.
 */
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

private proc unprojectedStateToIndex(basisStates : [] uint(64),
                                     isHammingWeightFixed : bool) {
  var indices : [0 ..# basisStates.size] int(64);
  if isHammingWeightFixed {
    ls_hs_fixed_hamming_state_to_index(basisStates.size, c_ptrToConst(basisStates), c_ptrTo(indices));
  }
  else {
    indices = basisStates:int(64);
  }
  return indices;
}
private proc unprojectedStateToIndex(_t : (uint(64), uint(64)), isHammingWeightFixed : bool) {
  const (low, high) = _t;
  const indices = unprojectedStateToIndex([low, high], isHammingWeightFixed);
  return (indices[0], indices[1]);
}
private proc unprojectedStateToIndex((low, high) : (uint(64), uint(64)), isHammingWeightFixed : bool) {
  const indices = unprojectedStateToIndex([low, high], isHammingWeightFixed);
  return (indices[0], indices[1]);
}
private proc unprojectedIndexToState(indices : [] int(64), numberSites : int, hammingWeight : int) {
  var basisStates : [0 ..# indices.size] uint(64);
  if hammingWeight != -1 {
    ls_hs_fixed_hamming_index_to_state(numberSites:c_int, hammingWeight:c_int, indices.size,
                                       c_ptrToConst(indices), c_ptrTo(basisStates));
  }
  else {
    basisStates = indices:uint(64);
  }
  return basisStates;
}

private proc assertBoundsHaveSameHammingWeight(r : range(?t), isHammingWeightFixed : bool) {
  if isHammingWeightFixed && popCount(r.low) != popCount(r.high) then
    halt(try! "r.low=%u and r.high=%u have different Hamming weight: %i != %i".format(
         r.low, r.high, popCount(r.low), popCount(r.high)));
}

/* We want to split `r` into `numChunks` non-overlapping ranges but such that for each range
   `chunk` we have that `popcount(chunk.low) == popcount(chunk.high) == popcount(r.low)` if
   `isHammingWeightFixed` was set to `true`.
 */
proc determineEnumerationRanges(r : range(uint(64)), in numChunks : int,
                                numberSites : int, hammingWeight : int) {
  var _timer = recordTime("determineEnumerationRanges");
  assertBoundsHaveSameHammingWeight(r, hammingWeight != -1);

  // Convert basisState bounds to index bounds
  const (lowIdx, highIdx) = unprojectedStateToIndex((r.low, r.high), hammingWeight != -1);
  const totalSize = highIdx - lowIdx + 1;
  if numChunks > totalSize then
    numChunks = totalSize;

  // Create ranges of indices
  var lowIndices : [0 ..# numChunks] int(64);
  var highIndices : [0 ..# numChunks] int(64);
  for (r, i) in zip(chunks(lowIdx .. highIdx, numChunks), 0..) {
    lowIndices[i] = r.low;
    highIndices[i] = r.high;
  }

  // Convert index ranges to ranges of basisStates
  var ranges : [0 ..# numChunks] range(uint(64));
  for (l, h, i) in zip(unprojectedIndexToState(lowIndices, numberSites, hammingWeight),
                       unprojectedIndexToState(highIndices, numberSites, hammingWeight),
                       0..) {
    ranges[i] = l .. h;
  }
  return ranges;
}
proc determineEnumerationRanges(r : range(uint(64)), numChunks : int,
                                const ref basisInfo : ls_hs_basis_info) {
  if basisInfo.particle_type == LS_HS_SPIN ||
     basisInfo.particle_type == LS_HS_SPINLESS_FERMION ||
     (basisInfo.particle_type == LS_HS_SPINFUL_FERMION && basisInfo.number_up == -1) {
    return determineEnumerationRanges(r, numChunks, basisInfo.number_bits, basisInfo.hamming_weight);
  }
  if basisInfo.particle_type == LS_HS_SPINFUL_FERMION {
    // Both the number of ↑ and ↓ spins are set.
    const numberSites = basisInfo.number_sites;
    const mask = (1:uint(64) << numberSites) - 1; // isolate the lower numberSites bits
    const numberUp = basisInfo.number_up;
    const numberDown = basisInfo.number_particles - numberUp;
    if numberUp < 0 || numberDown < 0 || numberUp + numberDown > 2 * numberSites then
      halt(try! "invalid combination of (number_up, number_particles, number_sites): (%i, %i, %i)".format(
           basisInfo.number_up, basisInfo.number_particles, basisInfo.number_sites));

    const minA = r.low & mask;
    const maxA = r.high & mask;
    const minB = (r.low >> numberSites) & mask;
    const maxB = (r.high >> numberSites) & mask;
    if popCount(minA) != numberUp || popCount(maxA) != numberUp ||
       popCount(minB) != numberDown || popCount(maxB) != numberDown then
      halt(try! "invalid r.low=%u and r.high=%u for a spinless fermion system with %i spins up and %i spins down".format(
           r.low, r.high, numberUp, numberDown));

    // TODO: This algorithm will be inefficient when we have only 1 spin down and many spins up.
    // In such cases, we should split into subranges using minA .. maxA range.
    const rangesB = determineEnumerationRanges(minB .. maxB, numChunks, basisInfo.number_sites, numberDown);
    var ranges : [0 ..# rangesB.size] range(uint(64));

    for (dest, rB) in zip(ranges, rangesB) {
      const low = (rB.low << numberSites) | minA;
      const high = (rB.high << numberSites) | maxA;
      dest = low .. high;
    }
    return ranges;
  }
  // This should never happen
  halt("unsupported basis type: particle_type=" + basisInfo.particle_type:string);
}

/* Hash function which we use to map spin configurations to locale indices.

   Typical usage:
   ```chapel
   var localeIndex = (hash64_01(x) % numLocales:uint):int;
   ```
*/
/*
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
*/

/* Given a vector of states `states`, write to `outMasks` the target locale for
   each state. The target locale is computed using a hash function.
   Furthermore, a histogram of target locales is computed, i.e., how many states
   belong to which locale. This information is returned as a 1D array.
 */
/*
private proc _enumStatesComputeMasksAndCounts(const ref states, ref outMasks) {
  const totalCount = states.size;
  var counts : [0 ..# numLocales] int;

  outMasks.resize(totalCount);
  // TODO: get rid of this clause altogether
  if numLocales == 1 { // all data belongs to locale 0
    if totalCount > 0 then
      POSIX.memset(outMasks.data, 0, totalCount:c_size_t * c_sizeof(outMasks.eltType));
    counts[0] = totalCount;
  }
  else {
    foreach i in 0 ..# totalCount {
      // TODO: generalize to an arbitrary hash function
      const key = localeIdxOf(states[i]);
      outMasks[i] = key:uint(8);
      counts[key] += 1;
    }
  }
  return counts;
}
*/

private proc _enumerateStatesProjected(r : range(uint(64)),
                                       const ref basis : Basis,
                                       ref outStates : Vector(uint(64)),
                                       ref outNorms : Vector(uint(16))) {
  var _timer = recordTime("_enumerateStatesProjected");
  if r.size == 0 then return;
  const ref basisInfo = basis.info;
  const isHammingWeightFixed = basisInfo.hamming_weight != -1;
  assertBoundsHaveSameHammingWeight(r, isHammingWeightFixed);

  var buffer: [0 ..# kIsRepresentativeBatchSize] uint(64);
  var norms: [0 ..# kIsRepresentativeBatchSize] real(64);
  var lower = r.low;
  const upper = r.high;

  while true {
    const written = manyNextState(lower, upper, buffer, isHammingWeightFixed);
    const last = buffer[written - 1];
    assert(last <= upper,
           "manyNextState went beyond the upper bound: "
           + last:string + " > " + upper:string);

    // TODO: get rid of array slices for better performance
    isRepresentative(basis, buffer[0 ..# written], norms[0 ..# written]);

    for i in 0 ..# written do
      if norms[i] > 0 {
        outStates.pushBack(buffer[i]);
        outNorms.pushBack(norms[i]:uint(16));
      }

    if last == upper then break;

    if isHammingWeightFixed { lower = nextState(last, true); }
    else { lower = nextState(last, false); }
  }
}
private proc _enumerateStatesUnprojected(r : range(uint(64)),
                                         const ref basis : Basis,
                                         ref outStates : Vector(uint(64)),
                                         ref outNorms : Vector(uint(16))) {
  var _timer = recordTime("_enumerateStatesUnprojected");
  const ref basisInfo = basis.info;
  const isHammingWeightFixed = basisInfo.hamming_weight != -1;
  const hasSpinInversionSymmetry = basisInfo.spin_inversion != 0;
  assertBoundsHaveSameHammingWeight(r, isHammingWeightFixed);
  const low = r.low;
  var high = r.high;
  if hasSpinInversionSymmetry {
    const numberSites = basisInfo.number_sites;
    const mask = (1:uint(64) << numberSites) - 1; // isolate the lower numberSites bits
    high = min(high, high ^ mask);
    if isHammingWeightFixed && popCount(high) != popCount(low) then
      halt(try! "Hamming weight changed upon a spin flip of %u".format(high));
  }

  const (lowIdx, highIdx) = unprojectedStateToIndex((low, high), isHammingWeightFixed);
  const totalCount = highIdx - lowIdx + 1;
  if totalCount > 0 {
    outStates.resize(totalCount);
    outNorms.resize(totalCount);
    manyNextState(low, high, outStates._arr, isHammingWeightFixed);
    // All norms are 1, because
    // - either no projection takes place (i.e., hasSpinInversionSymmetry==false)
    // - or we apply spin inversion, but x is always not equal to invert(x)
    outNorms._arr[0 ..# totalCount] = 1;
  }
}
private proc _enumerateStatesUnprojectedSpinfulFermion(r : range(uint(64)),
                                                       const ref basis : Basis,
                                                       ref outStates : Vector(uint(64)),
                                                       ref outNorms : Vector(uint(16))) {
  const ref basisInfo = basis.info;
  assert(basisInfo.particle_type == LS_HS_SPINFUL_FERMION,
         "_enumerateStatesUnprojectedSpinfulFermion only works for LS_HS_SPINFUL_FERMION, but particle_type="
         + basisInfo.particle_type:string);
  const numberSites = basisInfo.number_sites;
  const mask = (1 << numberSites) - 1; // isolate the lower numberSites bits
  const numberUp = basisInfo.number_up;
  const numberDown = basisInfo.number_particles - numberUp;
  assert(0 <= numberUp && 0 <= numberDown && numberUp + numberDown <= 2 * numberSites,
          "invalid combination of (number_up, number_particles, number_sites): "
          + (basisInfo.number_up, basisInfo.number_particles, basisInfo.number_sites):string);

  const minA = r.low & mask;
  const maxA = r.high & mask;
  assertBoundsHaveSameHammingWeight(minA .. maxA, isHammingWeightFixed=true);
  const (minAIdx, maxAIdx) = unprojectedStateToIndex((minA, maxA), isHammingWeightFixed=true);
  const countA = maxAIdx - minAIdx + 1;

  const minB = (r.low >> numberSites) & mask;
  const maxB = (r.high >> numberSites) & mask;
  assertBoundsHaveSameHammingWeight(minB .. maxB, isHammingWeightFixed=true);
  const (minBIdx, maxBIdx) = unprojectedStateToIndex((minB, maxB), isHammingWeightFixed=true);
  const countB = maxBIdx - minBIdx + 1;

  outStates.resize(countA * countB);
  outNorms.resize(countA * countB);
  var offset = 0;
  for xB in manyNextState(minB, maxB, isHammingWeightFixed=true) {
    for xA in manyNextState(minA, maxA, isHammingWeightFixed=true) {
      const x = (xB << numberSites) | xA;
      outStates[offset] = x;
      outNorms[offset] = 1;
      offset += 1;
    }
  }
}
proc _enumerateStatesGeneric(r : range(uint(64)),
                             const ref basis : Basis,
                             ref outStates : Vector(uint(64)),
                             ref outNorms : Vector(uint(16))) {
  const ref basisInfo = basis.info;
  select basisInfo.particle_type {
    when LS_HS_SPIN {
      if basisInfo.has_permutation_symmetries then
        _enumerateStatesProjected(r, basis, outStates, outNorms);
      else
        _enumerateStatesUnprojected(r, basis, outStates, outNorms);
    }
    when LS_HS_SPINLESS_FERMION {
      _enumerateStatesUnprojected(r, basis, outStates, outNorms);
    }
    when LS_HS_SPINFUL_FERMION {
      if basisInfo.number_up == -1 then
        _enumerateStatesUnprojected(r, basis, outStates, outNorms);
      else
        _enumerateStatesUnprojectedSpinfulFermion(r, basis, outStates, outNorms);
    }
    otherwise do halt("invalid particle_type: " + basisInfo.particle_type:string);
  }
}

/*
proc distributeByKeys(n : int, keysPtr : c_ptrConst(uint(8)),
                      arrPtr : c_ptr(?eltType), rDestPtrs : [] c_ptr(eltType)) {
  const _timer = recordTime("distributeByKeys");
  assert(numLocales <= 256,
         "distributeByKeys currently supports up to 256 locales");

  var offsets : c_array(int, 257);
  if numLocales == 1 { // We don't have to reshuffle stuff
    offsets[0] = 0;
    offsets[1] = n;
  }
  else {
    radixOneStep(n, keysPtr, offsets, arrPtr);
  }

  for localeIdx in 0 ..# numLocales {
    const count = offsets[localeIdx + 1] - offsets[localeIdx];
    if count > 0 {
      Communication.put(dest = rDestPtrs[localeIdx],
                        src = c_ptrTo(arrPtr[offsets[localeIdx]]),
                        destLocID = localeIdx,
                        numBytes = count:c_size_t * c_sizeof(eltType));
    }
  }
}
*/

class Bucket {
  var basisStates : Vector(uint(64));
  var norms : Vector(uint(16));
  var keys : Vector(uint(8));
  var offsets : c_array(int, 257);
}

/*
private proc enumStatesMakeBuckets(numChunks : int) {
  const dom = {0 ..# numChunks} dmapped cyclicDist(startIdx=0);
  var buckets : [dom] (Vector(uint(64)), Vector(real(64)), Vector(uint(8)));
  return buckets;
}
*/

proc enumStatesFillBuckets(ref buckets : [] owned Bucket,
                           const ref ranges : [] range(uint(64)),
                           const ref basis : Basis) {
  assert(buckets.size == ranges.size,
         "the number of buckets (" + buckets.size:string
         + ") does not match the number of ranges (" + ranges.size:string + ")");
  assert(here.id == 0, "enumStatesFillBuckets expects to be called from Locale 0");

  const numChunks = ranges.size;
  coforall loc in Locales do on loc {
    const myBasis = basis;
    const myRanges : [0 ..# numChunks] range(uint(64)) = ranges;
    const mySubdomain = buckets.localSubdomain();

    forall chunkIdx in dynamic(mySubdomain, chunkSize=1) {
      ref bucket = buckets.localAccess(chunkIdx);
      // Fill in basisStates and norms
      _enumerateStatesGeneric(myRanges[chunkIdx], myBasis, bucket.basisStates, bucket.norms);
      const count = bucket.basisStates.size;
      if numLocales > 1 {
        // Fill in keys
        const localeIdxFn = basis.getLocaleIdxFn();
        bucket.keys.resize(count);
        localeIdxFn(count, bucket.basisStates.data, bucket.keys.data);
        // Shuffle based on keys
        stableRadixOneStep(count, bucket.keys.data:c_ptrConst(uint(8)),
                           // TODO: explicit cast is **very** important here, see
                           // https://github.com/chapel-lang/chapel/issues/24318
                           bucket.offsets:c_ptr(int),
                           bucket.basisStates.data, bucket.norms.data);
      }
      else {
        bucket.offsets[0] = 0;
        bucket.offsets[1] = count;
        // TODO: should we fill the remaining offsets as well?
      }
    }
  }
}

proc enumStatesPerBucketCounts(const ref buckets : [] owned Bucket) {
  var perBucketLocaleCounts : [0 ..# buckets.size, 0 ..# numLocales] int;
  const perBucketLocaleCountsPtr = c_ptrTo(perBucketLocaleCounts[0, 0]);
  const localeIdx = here.id;
  forall (bucket, i) in zip(buckets, 0..) {
    var myCounts : [0 ..# numLocales] int;
    foreach i in 0 ..# numLocales do
      myCounts[i] = bucket.offsets[i + 1] - bucket.offsets[i];
    // TODO: replace with a simple assignment operator when rank-changing
    // copies become fast.
    Communication.put(dest = perBucketLocaleCountsPtr + i * numLocales,
                      src = c_ptrToConst(myCounts),
                      destLocID = localeIdx,
                      numBytes = numLocales:c_size_t * c_sizeof(int));
  }
  return perBucketLocaleCounts;
}

proc enumerateStates(ranges : [] range(uint(64)), const ref basis : Basis, out _basisStates, out _norms, out _keys) {
  // We distribute ranges among locales using Cyclic distribution to ensure
  // an even workload. For each range, a vector of basis states and a vector of
  // masks is computed. Masks indicate on which locale a basis state should live.
  const dom = {0 ..# ranges.size} dmapped cyclicDist(startIdx=0);
  var buckets : [dom] owned Bucket = [i in dom] new Bucket();
  enumStatesFillBuckets(buckets, ranges, basis);
  // How many states coming from a certain chunk should be sent to a particular locale.
  const perBucketLocaleCounts = enumStatesPerBucketCounts(buckets); // [0 ..# buckets.size, 0 ..# numLocales]

  // Allocate space for states
  const perLocaleCounts = sum(perBucketLocaleCounts, dim=0); // [0 ..# numLocales]
  const totalCount = sum(perLocaleCounts);
  var basisStates = new BlockVector(uint(64), perLocaleCounts);
  var norms = new BlockVector(uint(16), perLocaleCounts);

  // Perform transfers to target locales
  const basisStatesPtrs : [0 ..# numLocales] c_ptr(uint(64)) = basisStates._dataPtrs;
  const normsPtrs : [0 ..# numLocales] c_ptr(uint(16)) = norms._dataPtrs;
  const perLocaleOffsets = prefixSum(perBucketLocaleCounts, dim=0); // [0 ..# buckets.size, 0 ..# numLocales]
  forall (bucket, _i) in zip(buckets, 0..) {
    assert(here == bucket.locale);
    // TODO: check that these are bulk transferred
    const myBasisStatesPtrs = basisStatesPtrs;
    const myNormsPtrs = normsPtrs;
    // TODO: this one is probably not bulk transferred :(
    const myPerLocaleOffsets = perLocaleOffsets[_i, ..];

    for targetLocaleIdx in 0 ..# numLocales {
      const count = bucket.offsets[targetLocaleIdx + 1] - bucket.offsets[targetLocaleIdx];
      const srcOffset = bucket.offsets[targetLocaleIdx];
      const destOffset = myPerLocaleOffsets[targetLocaleIdx];
      Communication.put(dest = myBasisStatesPtrs[targetLocaleIdx] + destOffset,
                        src = bucket.basisStates.data + srcOffset,
                        destLocID = targetLocaleIdx,
                        numBytes = count:c_size_t * c_sizeof(uint(64)));
      Communication.put(dest = myNormsPtrs[targetLocaleIdx] + destOffset,
                        src = bucket.norms.data + srcOffset,
                        destLocID = targetLocaleIdx,
                        numBytes = count:c_size_t * c_sizeof(uint(16)));
    }

    // we're done with bucket.basisStates and bucket.norms, so we can free up some memory
    bucket.basisStates.clear();
    bucket.norms.clear();
  }

  // Allocate space for keys
  assert(numLocales <= max(uint(8)), "locale indices do not fit into an uint(8)");
  var keys = blockDist.createArray(0 ..# totalCount, uint(8));

  // Perform transfers to target locales
  if numLocales > 1 {
    const perBucketCounts = sum(perBucketLocaleCounts, dim=1); // [0 ..# buckets.size]
    const perBucketOffsets = prefixSum(perBucketCounts); // [0 ..# buckets.size]
    forall (bucket, i) in zip(buckets, 0..) {
      assert(here == bucket.locale);
      const myOffset = perBucketOffsets[i];
      keys[myOffset ..# bucket.keys.size] = bucket.keys.toArray();

      // we're done with bucket.keys, so we can free the memory
      bucket.keys.clear();
    }
  }

  _basisStates = basisStates;
  _norms = norms;
  _keys = keys;
}
proc enumerateStates(const ref basis : Basis, out basisStates, out norms, out keys,
                     in numChunks = kEnumerateStatesNumChunks) {
  const ref basisInfo = basis.info;
  const r = basisInfo.min_state_estimate .. basisInfo.max_state_estimate;
  // If we don't have to project, the enumeration ranges can be split equally,
  // so there's no need to use more chunks than there are cores
  if !basisInfo.requires_projection then
    numChunks = numLocales * here.maxTaskPar;
  numChunks = min(r.size, numChunks);
  const chunks = determineEnumerationRanges(r, numChunks, basisInfo);
  enumerateStates(chunks, basis, basisStates, norms, keys);
}

export proc ls_chpl_local_enumerate_states(p : c_ptr(ls_hs_basis),
                                           lower : uint(64),
                                           upper : uint(64)) {
  const basis = new Basis(p, owning=false);
  // var rs = localEnumerateRepresentatives(basis, lower, upper);
  const basisStates;
  const norms;
  const keys;

  const ref basisInfo = basis.info;
  const r = lower .. upper;
  // If we don't have to project, the enumeration ranges can be split equally,
  // so there's no need to use more chunks than there are cores
  const numChunks = min(r.size,
                        if !basisInfo.requires_projection
                          then numLocales * here.maxTaskPar
                          else kEnumerateStatesNumChunks);
  const chunks = determineEnumerationRanges(r, numChunks, basisInfo);
  enumerateStates(chunks, basis, basisStates, norms, keys);


  p.deref().local_representatives = convertToExternalArray(basisStates[here]);
  p.deref().local_norms = convertToExternalArray(norms[here]);
}

/*
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
*/

/*
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
*/

/*
proc _enumStatesMakeMasks(counts, totalCounts) {
  const numMasks = + reduce totalCounts;
  const masksBox = {0 ..# numMasks};
  const masksDom = masksBox dmapped blockDist(boundingBox=masksBox);
  var masks : [masksDom] uint(8);
  return masks;
}
*/

/*
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
*/

/*
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
  coforall loc in Locales with (ref copyTimes, ref maskCopyTimes) do on loc {
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
*/

/*
proc enumerateStates(globalRange : range(uint(64)), basis : Basis, out masks,
                     in numChunks : int = kEnumerateStatesNumChunks) {
  // If we don't have to project, the enumeration ranges can be split equally,
  // so there's no need to use more chunks than there are cores
  if !basis.info.requires_projection then
    numChunks = numLocales * here.maxTaskPar;
  const ranges = determineEnumerationRanges(globalRange, min(numChunks, globalRange.size), basis.info);
  return enumerateStates(ranges, basis, masks);
}
proc enumerateStates(const ref basis : Basis, out masks,
                     numChunks : int = kEnumerateStatesNumChunks) {
  
  const lower = basis.info.min_state_estimate;
  const upper = basis.info.max_state_estimate;
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
*/

}
