module HashedToBlock {
  use CommonParameters;
  use FFI;

  use BlockDist;
  use CTypes;
  use RangeChunk;
  use Time;

  private proc hashedToBlockComputeCounts(const ref masks : [] ?i, numChunks : int) {
    var counts : [0 ..# numLocales, 0 ..# numChunks, 0 ..# numLocales] int;
    const countsPtr = c_ptrTo(counts[counts.domain.low]);
    coforall loc in Locales do on loc {
      const mySubdomain = masks.localSubdomain(loc);
      const myRanges : [0 ..# numChunks] range(int) =
        chunks(mySubdomain.dim(0), numChunks);
      var myCounts : [0 ..# numChunks, 0 ..# numLocales] int;
      forall (r, chunkIdx) in zip(myRanges, 0 ..# numChunks) {
        foreach key in masks.localAccess(r) {
          myCounts[chunkIdx, key:int] += 1;
        }
      }
      const myCountsPtr = c_ptrToConst(myCounts[myCounts.domain.low]);
      const destPtr = countsPtr + loc.id * numChunks * numLocales;
      const destSize = (numChunks * numLocales):c_size_t * c_sizeof(int);
      PUT(myCountsPtr, 0, destPtr, destSize);
    }
    return counts;
  }

  private proc hashedToBlockComputeSrcOffsets(counts) {
    const numChunks = counts.shape[1];
    var offsets : [0 ..# numLocales, 0 ..# numChunks, 0 ..# numLocales] int;
    forall localeIdx in 0 ..# numLocales {
      var total = 0;
      for i in 0 ..# numLocales {
        for j in 0 ..# numChunks {
          offsets[i, j, localeIdx] = total;
          total += counts[i, j, localeIdx];
        }
      }
    }
    return offsets;
  }

  private proc hashedToBlockMakeDestArray(arr, masks) {
    param rank = arr.innerRank;
    const batchSize = if rank == 1 then 1
                                   else arr.innerDom.shape[0];
    const boundingBox = if rank == 1 then {0 ..# masks.size}
                                     else {0 ..# batchSize, 0 ..# masks.size};
    const targetLocales = if rank == 1 then Locales
                                       else reshape(Locales, {0 ..# 1, 0 ..# numLocales});
    const destArrDom = boundingBox dmapped Block(boundingBox, targetLocales);
    var destArr : [destArrDom] arr.eltType;
    return destArr;
  }

  private proc hashedToBlockNumChunks(masks, numChunks) {
    const minChunkSize = min reduce [loc in Locales] masks.localSubdomain(loc).size;
    return min(numChunks, minChunkSize);
  }

  record HashedToBlockTimer {
    var total : stopwatch;
    var counts : stopwatch;
    var distribute : stopwatch;

    proc writeThis(f) throws {
      f.write(timingTree(
        "arrFromHashedToBlock", total.elapsed(),
        [ ("counts", counts.elapsed())
        , ("distribute", distribute.elapsed())
        ]
      ));
    }
  }

  proc arrFromHashedToBlock(const ref arr, const ref masks,
                            numChunks = hashedToBlockNumChunks(masks, kHashedToBlockNumChunks)) {
    var timer = new HashedToBlockTimer();
    timer.total.start();

    timer.counts.start();
    const counts = hashedToBlockComputeCounts(masks, numChunks);
    timer.counts.stop();

    const srcOffsets = hashedToBlockComputeSrcOffsets(counts);
    var destArr = hashedToBlockMakeDestArray(arr, masks);
    const countsPtr = c_ptrToConst(counts[counts.domain.low]);
    const srcOffsetsPtr = c_ptrToConst(srcOffsets[srcOffsets.domain.low]);
    const mainLocaleIdx = here.id;
    const arrPtrsPtr = c_ptrToConst(arr._dataPtrs[0]);
    param rank = arr.innerRank;
    const batchSize = if rank == 1 then 1 else arr.innerDom.shape[0];
    const batchStride = if rank == 1 then 0 else arr.innerDom.shape[1];

    timer.distribute.start();
    coforall loc in Locales do on loc {
      const myDestSubdomain = destArr.localSubdomain();
      const myMasksSubdomain = masks.localSubdomain();
      // const mySize = myMasksSubdomain.size;
      const myRanges : [0 ..# numChunks] range(int) =
        chunks(myMasksSubdomain.dim(0), numChunks);

      var myCounts : [0 ..# numChunks, 0 ..# numLocales] int;
      var mySrcOffsets : [myCounts.domain] int;
      var myArrPtrs : [0 ..# numLocales] c_ptr(arr.eltType);
      cobegin {
        GET(c_ptrTo(myCounts[myCounts.domain.low]), mainLocaleIdx,
            countsPtr + loc.id * (numChunks * numLocales),
            (numChunks * numLocales):c_size_t * c_sizeof(int));
        GET(c_ptrTo(mySrcOffsets[mySrcOffsets.domain.low]), mainLocaleIdx,
            srcOffsetsPtr + loc.id * (numChunks * numLocales),
            (numChunks * numLocales):c_size_t * c_sizeof(int));
        GET(c_ptrTo(myArrPtrs[myArrPtrs.domain.low]), mainLocaleIdx,
            arrPtrsPtr, numLocales:c_size_t * c_sizeof(myArrPtrs.eltType));
      }
     
      for batchIdx in 0 ..# batchSize {
        forall (r, chunkIdx) in zip(myRanges, 0 ..# numChunks) {
          var maxCount = max reduce myCounts[chunkIdx, ..];
          var srcArr : [0 ..# numLocales, 0 ..# maxCount] arr.eltType;
          for srcLocaleIdx in 0 ..# numLocales {
            const n = myCounts[chunkIdx, srcLocaleIdx];
            const k = mySrcOffsets[chunkIdx, srcLocaleIdx];
            GET(c_ptrTo(srcArr[srcLocaleIdx, 0]), srcLocaleIdx,
                myArrPtrs[srcLocaleIdx] + batchIdx * batchStride + k,
                n:c_size_t * c_sizeof(srcArr.eltType));
          }

          var written : [0 ..# numLocales] int;
          for destIdx in r {
            const key = masks.localAccess(destIdx):int;
            ref offset = written[key];
            if rank == 1 then
              destArr.localAccess(destIdx) = srcArr[key, offset];
            else
              destArr.localAccess(batchIdx, destIdx) = srcArr[key, offset];
            offset += 1;
          }
        }
      }
    }
    timer.distribute.stop();
    timer.total.stop();
    if kDisplayTimings then logDebug(timer);

    return destArr;
  }

} // module HashedToBlock
