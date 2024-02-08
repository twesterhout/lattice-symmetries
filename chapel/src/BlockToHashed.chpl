module BlockToHashed {
  use CommonParameters;
  use FFI;
  import StatesEnumeration.permuteBasedOnMasks;
  import Vector.BlockVector;

  import Communication;
  use CTypes;
  use RangeChunk;
  use Time;

  proc blockArrHistograms(const ref arr : [] ?i, maxElement : i)
      where isBlockDist(masks.domain.distribution) &&
            isIntegral(i) {
    const ref arrTargetLocales = arr.targetLocales();
    const size = maxElement:int + 1;
    var counts : [0 ..# arrTargetLocales.size, 0 ..# size] int;
    const countsPtr = c_ptrTo(counts[0, 0]);
    const mainLocaleIdx = here.id;

    coforall loc in arrTargetLocales do on loc {
      const mySubdomain = arr.localSubdomain();
      var myCounts : [0 ..# size] int;
      forall x in arr.localAccess(mySubdomain) with (+ reduce myCounts) {
        assert(0 <= x && x <= maxElement);
        myCounts[x:int] += 1;
      }

      Communication.put(dest = countsPtr + loc.id * size,
                        src = c_ptrToConst(myCounts[0]),
                        destLocID = mainLocaleIdx,
                        numBytes = size:c_size_t * c_sizeof(int));
    }
    return counts;
  }

  private proc blockToHashedTaskCounts(masksSize : int, masksPtr : c_ptrConst(?i),
                                       numChunksPerLocale : int) {
    var counts : [0 ..# numChunksPerLocale, 0 ..# numLocales] int;
    const ranges : [0 ..# numChunksPerLocale] range(int) =
      chunks(0 ..# masksSize, numChunksPerLocale);

    forall (r, chunkIdx) in zip(ranges, 0 ..# numChunksPerLocale) {
      foreach i in r do
        counts[chunkIdx, masksPtr[i]:int] += 1;
    }
    return counts;
  }

  proc blockToHashedTaskOffsets(counts, perLocaleOffsets) {
    const numChunksPerLocale = counts.shape[0];
    var offsets : [0 ..# numChunksPerLocale, 0 ..# numLocales] int;
    for destLocaleIdx in 0 ..# numLocales {
      var total = perLocaleOffsets[destLocaleIdx];
      for chunkIdx in 0 ..# numChunksPerLocale {
          const count = counts[chunkIdx, destLocaleIdx];
          offsets[chunkIdx, destLocaleIdx] = total;
          total += count;
      }
    }
    return offsets;
  }

  proc blockToHashedNumChunksPerLocale(masks, numChunks) {
    const minBlockSize = min reduce [loc in Locales] masks.localSubdomain(loc).size;
    const suggested = max(1, numChunks / numLocales);
    return min(minBlockSize, suggested);
  }

  proc arrFromBlockToHashed(const ref arr : [] ?eltType, const ref keys : [] ?i)
      where isBlockDist(arr.domain.distribution) &&
            isBlockDist(keys.domain.distribution) &&
            (arr.rank == 1 || arr.rank == 2) {

    param rank = arr.domain.rank;
    const maxElement = max reduce keys;
    assert(maxElement < numLocales,
           "maxElement=%i, but numLocales=%i".format(maxElement, numLocales));
    const perLocaleLocaleCounts = blockArrHistograms(keys, (numLocales - 1):i);
    const perLocaleOffsets = prefixSum(perLocaleCounts, dim=0);

    const batchSize = if rank == 1 then 1 else arr.dim(0).size;
    const destShapes : [0 ..# numLocales] = [n in sum(perLocaleLocaleCounts, dim=0)] if rank == 1 then n else (batchSize, n);
    var destArr = new BlockVector(eltType, destShapes);

    const mainLocaleIdx = here.id;
    const perLocaleOffsetsPtr = c_ptrToConst(perLocaleOffsets[0]);
    coforall loc in Locales do on loc {

      var myPerLocaleOffsets : [0 ..# numLocales] int;
      Communication.get(dest = c_ptrTo(myPerLocaleOffsets[0]),
                        src = perLocaleOffsetsPtr + loc.id * numLocales,
                        srcLocID = mainLocaleIdx,
                        numBytes = numLocales:c_size_t * c_sizeof(int));

      var myDestPtrs : [0 ..# numLocales] c_ptr(eltType);
      Communication.get(dest = c_ptrTo(myDestPtrs[0]),
                        src = destPtrsPtr,
                        srcLocID = mainLocaleIdx,
                        numBytes = numLocales:c_size_t * c_sizeof(c_ptr(eltType)));

      const mySubdomain = masks.localSubdomain();
      const myMasksPtr = c_ptrToConst(masks.localAccess(mySubdomain.low));
      const myMasksSize = mySubdomain.size;
      const myPerTaskCounts = blockToHashedTaskCounts(myMasksSize, myMasksPtr, numChunksPerLocale);
      const myPerTaskOffsets = blockToHashedTaskOffsets(myPerTaskCounts, myPerLocaleOffsets);
      const ranges : [0 ..# numChunksPerLocale] range(int) = chunks(0 ..# myMasksSize, numChunksPerLocale);

      var myPermuteTime : atomic real;
      for batchIdx in 0 ..# batchSize {
        const myArrPtr = if rank == 1
                           then c_ptrToConst(arr.localAccess(mySubdomain.low))
                           else c_ptrToConst(arr.localAccess(batchIdx, mySubdomain.low));

        forall (r, chunkIdx) in zip(ranges, 0 ..# numChunksPerLocale) {
          var myTimer = new stopwatch();
          myTimer.start();
          permuteBasedOnMasks(r.size,
                              myMasksPtr + r.low,
                              myArrPtr + r.low,
                              myPerTaskCounts[chunkIdx, ..],
                              myPerTaskOffsets[chunkIdx, ..],
                              myDestPtrs);
          myTimer.stop();
          myPermuteTime.add(myTimer.elapsed(), memoryOrder.relaxed);
        }
        myDestPtrs += batchStride;
      }

      const _myPermuteTime = myPermuteTime.read();
      Communication.put(dest = permuteTimePtr + loc.id,
                        src = c_ptrToConst(_myPermuteTime),
                        destLocID = mainLocaleIdx,
                        numBytes = c_sizeof(real));
    }
    timer.distribute.stop();
    timer.total.stop();

    if kDisplayTimings then logDebug(timer);
    return destArr;
  }

}
