module BlockToHashed {
  use FFI;
  use Vector;
  use CommonParameters;

  use CTypes;
  use RangeChunk;
  use Time;

  private proc _blockToHashedLocaleCounts(const ref masks : [] ?i) where isIntegral(i) {
    var counts : [0 ..# numLocales, 0 ..# numLocales] int;
    const countsPtr = c_ptrTo(counts[counts.domain.low]);
    const mainLocaleIdx = here.id;
    coforall loc in Locales do on loc {
      const mySubdomain = masks.localSubdomain();
      // const ref myMasks = masks[masks.localSubdomain()];
      var myCounts : [0 ..# numLocales] int;
      forall key in masks.localAccess(mySubdomain) with (+ reduce myCounts) {
        myCounts[key:int] += 1;
      }

      const putOffset = loc.id * numLocales;
      const putSize = numLocales:c_size_t * c_sizeof(int);
      PUT(c_ptrTo(myCounts[0]), mainLocaleIdx, countsPtr + putOffset, putSize);
    }
    return counts;
  }

  private proc _blockToHashedLocaleOffsets(counts) {
    var offsets : [0 ..# numLocales, 0 ..# numLocales] int;
    foreach destLocaleIdx in 0 ..# numLocales {
      var total = 0;
      for srcLocaleIdx in 0 ..# numLocales {
        const count = counts[srcLocaleIdx, destLocaleIdx];
        offsets[srcLocaleIdx, destLocaleIdx] = total;
        total += count;
      }
    }
    return offsets;
  }

  private proc _blockToHashedTaskCounts(masksSize : int, masksPtr : c_ptr(?i),
                                        numChunksPerLocale : int) {
    // if kVerboseComm then startVerboseCommHere();
    var counts : [0 ..# numChunksPerLocale, 0 ..# numLocales] int;
    const ranges : [0 ..# numChunksPerLocale] range(int) =
      chunks(0 ..# masksSize, numChunksPerLocale);

    forall (r, chunkIdx) in zip(ranges, 0 ..# numChunksPerLocale) {
      foreach i in r do
        counts[chunkIdx, masksPtr[i]:int] += 1;
    }
    return counts;
  }

  proc _blockToHashedTaskOffsets(counts, perLocaleOffsets) {
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

  proc _blockToHashedMakeDestArr(arr : [] ?eltType, counts)
      where arr.domain.rank == 1 {
    const destCounts = [i in 0 ..# numLocales] (+ reduce counts[.., i]);
    return new BlockVector(eltType, destCounts);
  }
  proc _blockToHashedMakeDestArr(arr : [] ?eltType, counts)
      where arr.domain.rank == 2 {
    const destCounts = [i in 0 ..# numLocales] (+ reduce counts[.., i]);
    return new BlockVector(eltType, arr.shape[0], destCounts);
  }

  proc _blockToHashedNumChunksPerLocale(masks, numChunks) {
    const minBlockSize = min reduce [loc in Locales] masks.localSubdomain(loc).size;
    const suggested = max(1, numChunks / numLocales);
    return min(minBlockSize, suggested);
  }

  proc arrFromBlockToHashed(const ref arr : [] ?eltType, const ref masks : [] ?i,
                            numChunksPerLocale : int =
                              _blockToHashedNumChunksPerLocale(masks,
                                kBlockToHashedNumChunks)) {
    var timer = new Timer();
    var countsTimer = new Timer();
    var makeDestArrTimer = new Timer();
    var distributeTimer = new Timer();
    var permuteTime : [0 ..# numLocales] real;
    const permuteTimePtr = c_ptrTo(permuteTime[0]);
    timer.start();

    countsTimer.start();
    const perLocaleCounts = _blockToHashedLocaleCounts(masks);
    const perLocaleOffsets = _blockToHashedLocaleOffsets(perLocaleCounts);
    const perLocaleOffsetsPtr = c_const_ptrTo(perLocaleOffsets[perLocaleOffsets.domain.low]);
    countsTimer.stop();

    makeDestArrTimer.start();
    var destArr = _blockToHashedMakeDestArr(arr, perLocaleCounts);
    const destPtrsPtr = c_ptrTo(destArr._dataPtrs);
    makeDestArrTimer.stop();

    distributeTimer.start();
    param rank = arr.domain.rank;
    const batchSize = if rank == 1
                        then 1
                        else arr.shape[0];
    const batchStride = if rank == 1
                          then 0
                          else destArr.innerDom.shape[1];
    const mainLocaleIdx = here.id;
    coforall loc in Locales do on loc {
      // const _myPtr : c_ptr(eltType) = c_ptrTo(destArr[loc.id, destArr.innerDom.low]);
      // if kUseLowLevelComm then
      //   PUT(c_const_ptrTo(_myPtr), 0, destPtrsPtr + loc.id, c_sizeof(c_ptr(eltType)));
      // else
      //   destPtrs[loc.id] = _myPtr;

      // allLocalesBarrier.barrier();

      var myPerLocaleOffsets : [0 ..# numLocales] int;
      GET(c_ptrTo(myPerLocaleOffsets[0]), mainLocaleIdx,
          perLocaleOffsetsPtr + loc.id * numLocales,
          numLocales:c_size_t * c_sizeof(int));

      var myDestPtrs : [0 ..# numLocales] c_ptr(eltType);
      GET(c_ptrTo(myDestPtrs[0]), mainLocaleIdx, destPtrsPtr,
          numLocales:c_size_t * c_sizeof(c_ptr(eltType)));

      const mySubdomain = masks.localSubdomain();
      const myMasksPtr = c_const_ptrTo(masks.localAccess(mySubdomain.low));
      const myMasksSize = mySubdomain.size;

      const myPerTaskCounts =
        _blockToHashedTaskCounts(myMasksSize, myMasksPtr, numChunksPerLocale);
      const myPerTaskOffsets = _blockToHashedTaskOffsets(myPerTaskCounts, myPerLocaleOffsets);
      // const (perTaskCount, perTaskOffset) =
      //   getPerTaskCountAndOffset(myMasksSize, myMasksPtr, myPerLocaleOffset);
      const ranges : [0 ..# numChunksPerLocale] range(int) =
        chunks(0 ..# myMasksSize, numChunksPerLocale);

      // var computeTimer = new Timer();
      // computeTimer.start();
      // var permuteSmallPrepareTime : atomic real;
      // var permuteSmallComputeTime : atomic real;
      var myPermuteTime : atomic real;

      for batchIdx in 0 ..# batchSize {
        const myArrPtr = if rank == 1
                           then c_const_ptrTo(arr.localAccess(mySubdomain.low))
                           else c_const_ptrTo(arr.localAccess(batchIdx, mySubdomain.low));

        forall (r, chunkIdx) in zip(ranges, 0 ..# numChunksPerLocale) {
          var myTimer = new Timer();
          myTimer.start();
          permuteBasedOnMasks(r.size,
                              myMasksPtr + r.low,
                              myArrPtr + r.low,
                              myPerTaskCounts[chunkIdx, ..],
                              myPerTaskOffsets[chunkIdx, ..],
                              myDestPtrs);
          myTimer.stop();
          myPermuteTime.add(myTimer.elapsed(), memoryOrder.relaxed);
          // proc permuteBasedOnMasks(arrSize : int, masks : c_ptr(?maskType), arr : c_ptr(?eltType),
          //                          counts : [] int, destOffsets : [] int,
          //                          destPtrs : [] c_ptr(eltType)) {
          // const _copy =
          //   permuteSmall(r.size,
          //                myMasksPtr + r.low,
          //                myArrPtr + r.low,
          //                perTaskCount[chunkIdx, ..],
          //                perTaskOffset[chunkIdx, ..],
          //                myDestPtrs);
          // permuteSmallPrepareTime.add(_prep);
          // permuteSmallComputeTime.add(_compute);
          // permuteSmallCopyTime.add(_copy, memoryOrder.relaxed);
        }
        myDestPtrs += batchStride;
      }

      const _myPermuteTime = myPermuteTime.read();
      PUT(c_const_ptrTo(_myPermuteTime), mainLocaleIdx,
          permuteTimePtr + loc.id,
          c_sizeof(real));
      // computeTimer.stop();
      // logDebug("arrFromBlockToHashed main loop took ", computeTimer.elapsed(), "\n",
      //          "                      from which ", permuteSmallCopyTime.read(),
      //          " were spent in remote PUTs (divide by the number of tasks because ",
      //          "they ran in parallel)");
    }
    distributeTimer.stop();

    timer.stop();
    logDebug("arrFromBlockToHashed took ", timer.elapsed(), "\n",
             "  ├─ ", countsTimer.elapsed(), " were spent computing counts\n",
             "  ├─ ", makeDestArrTimer.elapsed(), " allocating destArr\n",
             "  └─ ", distributeTimer.elapsed(), " in the main loop\n",
             "      └─ ", permuteTime, " in permuteBasedOnMasks");
    // if kVerboseComm then stopVerboseComm();
    return destArr;
  }

}
