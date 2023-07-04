module BlockToHashed {
  use CommonParameters;
  use FFI;
  import StatesEnumeration.permuteBasedOnMasks;
  import Vector.BlockVector;

  import Communication;
  use CTypes;
  use RangeChunk;
  use Time;

  private proc blockToHashedLocaleCounts(const ref masks : [] ?i) where isIntegral(i) {
    var counts : [0 ..# numLocales, 0 ..# numLocales] int;
    const countsPtr = c_ptrTo(counts[counts.domain.low]);
    const mainLocaleIdx = here.id;
    coforall loc in Locales do on loc {
      const mySubdomain = masks.localSubdomain();
      var myCounts : [0 ..# numLocales] int;
      forall key in masks.localAccess(mySubdomain) with (+ reduce myCounts) do
        myCounts[key:int] += 1;

      Communication.put(dest = countsPtr + loc.id * numLocales,
                        src = c_ptrTo(myCounts[0]),
                        destLocID = mainLocaleIdx,
                        numBytes = numLocales:c_size_t * c_sizeof(int));
    }
    return counts;
  }

  private proc blockToHashedLocaleOffsets(counts) {
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

  private proc blockToHashedTaskCounts(masksSize : int, masksPtr : c_ptr(?i),
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

  proc blockToHashedMakeDestArr(arr : [] ?eltType, counts) {
    const destCounts = [i in 0 ..# numLocales] (+ reduce counts[.., i]);
    return
      if arr.domain.rank == 1
        then new BlockVector(eltType, destCounts)
        else new BlockVector(eltType, arr.shape[0], destCounts);
  }

  proc blockToHashedNumChunksPerLocale(masks, numChunks) {
    const minBlockSize = min reduce [loc in Locales] masks.localSubdomain(loc).size;
    const suggested = max(1, numChunks / numLocales);
    return min(minBlockSize, suggested);
  }

  record BlockToHashedTimer {
    var total : stopwatch;
    var counts : stopwatch;
    var makeDestArr : stopwatch;
    var distribute : stopwatch;
    var permute : [0 ..# numLocales] real;

    proc writeThis(f) throws {
      const tree = timingTree("arrFromBlockToHashed", total.elapsed(),
        [ ("counts", counts.elapsed())
        , ("makeDestArr", makeDestArr.elapsed())
        , ("distribute", distribute.elapsed())
        ]);
      tree.children[2].addChild(timingTree(
        "permute (mean over locales, sum over tasks)", meanAndErr(permute)[0]
      ));

      f.write(tree);
    }
  }


  proc arrFromBlockToHashed(const ref arr : [] ?eltType, const ref masks : [] ?i,
                            numChunksPerLocale : int =
                              blockToHashedNumChunksPerLocale(masks, kBlockToHashedNumChunks)) {
    var timer = new BlockToHashedTimer();
    const permuteTimePtr = c_ptrTo(timer.permute[0]);
    timer.total.start();

    timer.counts.start();
    const perLocaleCounts = blockToHashedLocaleCounts(masks);
    const perLocaleOffsets = blockToHashedLocaleOffsets(perLocaleCounts);
    const perLocaleOffsetsPtr = c_const_ptrTo(perLocaleOffsets[perLocaleOffsets.domain.low]);
    timer.counts.stop();

    timer.makeDestArr.start();
    var destArr = blockToHashedMakeDestArr(arr, perLocaleCounts);
    const destPtrsPtr = c_ptrTo(destArr._dataPtrs);
    timer.makeDestArr.stop();

    timer.distribute.start();
    param rank = arr.domain.rank;
    const batchSize = if rank == 1
                        then 1
                        else arr.shape[0];
    const batchStride = if rank == 1
                          then 0
                          else destArr.innerDom.shape[1];
    const mainLocaleIdx = here.id;
    coforall loc in Locales do on loc {
      var myPerLocaleOffsets : [0 ..# numLocales] int = noinit;
      // NOTE: Actually, we want myPerLocaleOffsets = perLocaleOffsets[loc.id, ..],
      // but Chapel fails to generate a bulk copy operation and copies element-by-element.
      Communication.get(dest = c_ptrTo(myPerLocaleOffsets[0]),
                        src = perLocaleOffsetsPtr + loc.id * numLocales,
                        srcLocID = mainLocaleIdx,
                        numBytes = numLocales:c_size_t * c_sizeof(int));

      var myDestPtrs : [0 ..# numLocales] c_ptr(eltType) = noinit;
      // TODO: we're not doing rank-changing copies here, so myDestPtrs = destPtrs
      // should do the trick... Check it!
      Communication.get(dest = c_ptrTo(myDestPtrs[0]),
                        src = destPtrsPtr,
                        srcLocID = mainLocaleIdx,
                        numBytes = numLocales:c_size_t * c_sizeof(c_ptr(eltType)));

      const mySubdomain = masks.localSubdomain();
      const myMasksPtr = c_const_ptrTo(masks.localAccess(mySubdomain.low));
      const myMasksSize = mySubdomain.size;
      const myPerTaskCounts = blockToHashedTaskCounts(myMasksSize, myMasksPtr, numChunksPerLocale);
      const myPerTaskOffsets = blockToHashedTaskOffsets(myPerTaskCounts, myPerLocaleOffsets);
      const ranges : [0 ..# numChunksPerLocale] range(int) = chunks(0 ..# myMasksSize, numChunksPerLocale);

      var myPermuteTime : atomic real;
      for batchIdx in 0 ..# batchSize {
        const myArrPtr = if rank == 1
                           then c_const_ptrTo(arr.localAccess(mySubdomain.low))
                           else c_const_ptrTo(arr.localAccess(batchIdx, mySubdomain.low));

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
                        src = c_const_ptrTo(_myPermuteTime),
                        destLocID = mainLocaleIdx,
                        numBytes = c_sizeof(real));
    }
    timer.distribute.stop();
    timer.total.stop();

    if kDisplayTimings then logDebug(timer);
    return destArr;
  }

}
