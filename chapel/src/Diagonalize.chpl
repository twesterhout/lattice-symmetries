module Diagonalize {
  use AllLocalesBarriers;
  use CTypes;
  use CommDiagnostics;
  import Random;

  use BlockToHashed;
  use HashedToBlock;
  use LatticeSymmetries;
  use MyHDF5;
  use PRIMME;

  /*
  // A buffer located on locale 0 to help with the broadcast                      
  var tmpBuffDom = {0..0:c_int};
  var tmpBuff: [tmpBuffDom] real;

  proc broadcastReal(buffer: c_ptr(real), count: c_ptr(c_int)) {
    const n = count.deref(),
          inds = 0..<n;

    if here.id == 0 {
      // grow the temp buff if it's not big enough                                
      if n > tmpBuffDom.size then
        tmpBuffDom = {inds};

      // copy locale 0's data into the buffer                                     
      forall i in inds do
        tmpBuff[i] = buffer[i];
    }

    // wait until locale 0's got tmpBuff set up before proceeding                 
    allLocalesBarrier.barrier();

    // Locale 0 already has the data so doesn't need to do anything               
    if (here.id != 0) then
      forall i in inds do
        buffer[i] = tmpBuff[i];
  }

  export proc ls_chpl_broadcast_real(sendBuf : c_void_ptr, count : c_ptr(c_int),
                                     primme : c_ptr(primme_params), ierr : c_ptr(c_int)) {
    if primme.deref().broadcastReal_type != primme_op_double then
      halt("broadcastReal is implemented for double precision only");
    broadcastReal(sendBuf : c_ptr(real), count);
    ierr.deref() = 0;
  }


  record AtomicBuf {
    var dom : domain(1);
    var arr : [dom] atomic real;
  }

  // A buffer of atomics on locale 0 for computing the reduction                  
  // var atomicBuffDom: domain(1); // = {0..0:c_int};
  var atomicBuff = new AtomicBuf();
  // : [atomicBuffDom] atomic real;

  // var D : domain(1);
  // var A : [D] real;

  proc globalSumReal(sendBuf: c_ptr(real), recvBuf: c_ptr(real),
                     count: c_ptr(c_int)) {
    // startVerboseCommHere();
    const n = count.deref();
    const inds = 0 ..# n:int;
    // writeln("calling globalSumReal(", sendBuf, ", ", recvBuf, ", ", n, ") from ", here);

    const sendArr = [i in inds] sendBuf[i];
    writeln(here, ": sendBuf=", ([i in inds] sendBuf[i]):string);
    writeln(here, ": sendArr=", sendArr);
    writeln(here, ": recvBuf=", ([i in inds] recvBuf[i]):string);

    // ref atomicBuffRef = atomicBuff;

    // grow the temp buff if it's not big enough                                  
    if here.id == 0 then
      if n > atomicBuff.dom.size {
        // atomicBuffDom = {inds};
        // writeln(atomicBuffDom);
        // writeln(atomicBuff);
        atomicBuff.dom = {inds};

        // D = {0 ..# n:int};
      }

    // Make sure locale 0 has had the chance to resize before proceeding          
    allLocalesBarrier.barrier();
    assert(atomicBuff.dom.size >= n);
    assert(atomicBuff.arr.size >= n);
    // writeln(here, ": D = ", D, ", A = ", A);

    // allLocalesBarrier.barrier();
    // assert(D.size >= n);
    // assert(A.domain.size >= n);
    // assert(atomicBuffRef.size >= n);

    // have all locales atomically add their results to the atomicBuff            
    for i in inds {
      writeln("adding " + sendBuf[i]:string + " to atomicBuff.arr[" + i:string + "]");
      const elt = sendBuf[i];
      assert(elt == sendArr[i]);
      atomicBuff.arr[i].add(sendBuf[i]);
    }

    // Make sure all locales have accumulated their contributions                 
    allLocalesBarrier.barrier();

    // Have each locale copy the results out into its buffer                      
    for i in inds do
      recvBuf[i] = atomicBuff.arr[i].read();

    writeln(here, ": recvBuf=", ([i in inds] recvBuf[i]):string);
    // stopVerboseCommHere();
  }

  export proc ls_chpl_global_sum_real(sendBuf : c_void_ptr, recvBuf : c_void_ptr,
                                      count : c_ptr(c_int), primme : c_ptr(primme_params),
                                      ierr : c_ptr(c_int)) {
    if primme.deref().globalSumReal_type != primme_op_double then
      halt("globalSum is implemented for double precision only");
    globalSumReal(sendBuf : c_ptr(real), recvBuf : c_ptr(real), count);
    ierr.deref() = 0;
  }

  */

  inline proc borrowOperator(primme : c_ptr(primme_params)) {
    const p = primme.deref().matrix : c_ptr(ls_hs_operator);
    return new Operator(p, owning=false);
  }

  export proc ls_chpl_primme_matvec(_x : c_void_ptr, _ldx : c_ptr(int(64)),
                                    _y : c_void_ptr, _ldy : c_ptr(int(64)),
                                    _blockSize : c_ptr(c_int), primme : c_ptr(primme_params),
                                    _ierr : c_ptr(c_int)) {
    // logDebug("Calling ls_chpl_primme_matvec ...");
    const blockSize = _blockSize.deref():int;
    // if blockSize != 1 then
    //   halt("currently only blockSize=1 is supported");
    const ldx = _ldx.deref();
    const ldy = _ldy.deref();
    const n = primme.deref().nLocal;
    assert(ldx >= n);
    assert(ldy >= n);

    // const precision = primme.deref().matrixMatvec_type;
    // const dtype = 
    type eltType = real;
    const matrix = borrowOperator(primme);
    const ref representatives = matrix.basis.representatives();

    for k in 0 ..# blockSize {
      ref X = makeArrayFromPtr(_x : c_ptr(eltType) + ldx * k, (n,));
      ref Y = makeArrayFromPtr(_y : c_ptr(eltType) + ldy * k, (n,));
      localMatrixVector(matrix, X, Y, representatives);
    }

    _ierr.deref() = 0;
    // logDebug("Done with ls_chpl_primme_matvec ...");
  }

  config const input : string = "data/heisenberg_chain_10.yaml";
  config const kOutput : string = "exact_diagonalization_output.h5";
  config const numEvals : int = 1;
  config const kEps : real = 1e-6;

  // Advanced:
  config const kMaxBasisSize : int = 0;
  config const kMinRestartSize : int = 0;
  config const kMaxBlockSize : int = 1;

  proc configurePrimme(ref params, dimension, matrix, basisStates) {
      params.n = dimension;
      params.matrixMatvec = c_ptrTo(ls_chpl_primme_matvec);
      params.numEvals = numEvals:c_int;
      params.target = primme_smallest;
      params.eps = kEps;

      params.numProcs = numLocales:c_int;
      params.procID = here.id:c_int;
      params.nLocal = basisStates.size;
      params.globalSumReal = c_ptrTo(primmeGlobalSumReal);
      params.globalSumReal_type = primme_op_double;
      params.broadcastReal = c_ptrTo(primmeBroadcastReal);
      params.broadcastReal_type = primme_op_double;

      // params.initSize = ;
      params.maxBasisSize = kMaxBasisSize:c_int;
      params.minRestartSize = kMinRestartSize:c_int;
      params.maxBlockSize = kMaxBlockSize:c_int;

      // params.commInfo = ;
      params.matrix = matrix.payload;
      // params.massMatrix = ;
      // params.preconditioner = ;
      // params.convtest = ;
      // params.monitor = ;

      // params.ldevecs = ;
      // params.numOrthoConst = ;
      // params.dynamicMethodSwitch = ;
      // params.locking = ;
      // params.maxMatvecs = ;
      // params.maxOuterIterations = ;
      // params.iseed = ;
      // params.aNorm = ;
      // params.BNorm = ;
      // params.invBNorm = ;
      params.printLevel = 5;
      // params.outputFile = ;
      // params.ShiftsForPreconditioner = ;
      // params.initBasisMode = ;
      // params.projectionParams = ;
      // params.restartingParams = ;
      // params.correctionParams = ;
      // params.stats = ;
      // params.convTestFun = ;
      // params.ldOPs = ;
      // params.monitorFun = ;
      // params.orth = ;

      primme_set_method(PRIMME_DEFAULT_MIN_MATVECS, params);
  }

  proc makeBasisStates(basis : Basis, out basisStates, out masks) {
    if doesObjectExist(kOutput, "basis/representatives") {
      logDebug("Reading representatives from '" + kOutput + "' ...");
      const blockBasisStates =
        readDatasetAsBlocks(kOutput, "basis/representatives", rank = 1, eltType = uint(64));
      const tempMasks = [x in blockBasisStates] localeIdxOf(x):uint(8);
      basisStates = arrFromBlockToHashed(blockBasisStates, tempMasks);
      masks = tempMasks;
    }
    else {
      logDebug("Generating basis states ...");
      const tempMasks;
      basisStates = enumerateStates(basis, tempMasks);
      logDebug("Writing representatives to '" + kOutput + "' ...");
      writeDatasetAsBlocks(kOutput, "basis/representatives",
                           arrFromHashedToBlock(basisStates, tempMasks));
      // basisStates = tempBasisStates;
      masks = tempMasks;
    }
  }

  proc saveEigenvectors(evals : [] real(64),
                        evecs : BlockVector,
                        resNorms : [] real(64),
                        masks : [] uint(8)) {
    writeDatasetAsBlocks(kOutput, "hamiltonian/eigenvectors",
                         arrFromHashedToBlock(evecs, masks));
    writeDataset(kOutput, "hamiltonian/eigenvalues", evals);
    writeDataset(kOutput, "hamiltonian/residuals", resNorms);
  }

  proc main() {
    initRuntime();
    defer deinitRuntime();

    // Parse the input configuration file
    const (basis, matrix, observables) =
      loadConfigFromYaml(input, hamiltonian=true, observables=true);

    if numEvals < 1 then
      halt("invalid numEvals: " + numEvals:string);

    // if !matrix.isHermitian then
    //   halt("Hamiltonian is not Hermitian");

    // if !matrix.isReal then
    //   halt("Hamiltonian is not real");

    // Create the output file
    makeFile(kOutput);
    makeGroup(kOutput, "basis");
    makeGroup(kOutput, "hamiltonian");
    makeGroup(kOutput, "observables");

    // Build the basis
    const basisStates;
    const masks;
    makeBasisStates(basis, basisStates, masks);

    writeln(basisStates);
    writeln(masks);

    const dimension = + reduce basisStates.counts;
    logDebug("Hilbert space dimension: ", dimension);

    // Allocate space for eigenvectors
    var evecs = new BlockVector(real(64), numEvals, basisStates.counts);
    var evals : [0 ..# numEvals] real(64);
    var resNorms : [0 ..# numEvals] real(64);

    logDebug("Diagonalizing the Hamiltonian ...");
    coforall loc in Locales with (ref evecs)
                            do on loc {
      const myMatrix = matrix;
      const ref myBasisStates = basisStates.getBlock(loc.id);
      myMatrix.basis.uncheckedSetRepresentatives(myBasisStates);

      ref myEvecs = evecs.getBlock(loc.id);
      var myEvals : [0 ..# numEvals] real(64);
      var myResNorms : [0 ..# numEvals] real(64);

      // Configuring PRIMME
      var params : primme_params;
      primme_initialize(params);
      defer primme_free(params);

      configurePrimme(params, dimension, myMatrix, myBasisStates);
      const ierr = dprimme(c_ptrTo(myEvals),
                           c_ptrTo(myEvecs[myEvecs.domain.low]),
                           c_ptrTo(myResNorms),
                           params);
      if ierr != 0 then
        halt("PRIMME failed with error code: " + ierr:string);

      if loc.id == 0 {
        evals = myEvals;
        resNorms = myResNorms;
      }
    }

    logDebug("Obtained eigenvalues: ", evals);
    logDebug("Residual norms:       ", resNorms);

    saveEigenvectors(evals, evecs, resNorms, masks);

  }
}
