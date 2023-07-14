use LatticeSymmetries;
use MyHDF5;
use HashedToBlock;
use BlockToHashed;
use Time;

proc localLoadVectors(filename : string, x : string = "/x", y : string = "/y") {
  var input = readDataset(filename, x, real(64), rank = 2)[0, ..];
  var output = readDataset(filename, y, real(64), rank = 2)[0, ..];
  return (input, output);
}

config const kHamiltonian = "data/heisenberg_chain_10.yaml";
config const kVectors = "data/matvec/heisenberg_chain_10.h5";
config const kAbsTol = 1e-13;
config const kRelTol = 1e-11;

proc approxEqual(a : real, b : real, atol = kAbsTol, rtol = kRelTol) {
  return abs(a - b) <= max(atol, rtol * max(abs(a), abs(b)));
}
proc approxEqual(a : [], b : [], atol = kAbsTol, rtol = kRelTol) {
  return [i in a.domain] approxEqual(a[i], b[i], atol, rtol);
}

proc main() {
  initRuntime();
  defer deinitRuntime();

  const (_, matrix) = loadConfigFromYaml(kHamiltonian, hamiltonian=true);

  const masks;
  const basisStates = enumerateStates(matrix.basis, masks);

  const x = arrFromBlockToHashed(readDatasetAsBlocks(kVectors, "/x"), masks);

  var z = similar(x);

  var timer = new stopwatch();
  timer.start();
  matrixVectorProduct(matrix, x, z, basisStates);
  timer.stop();

  const yBlock = readDatasetAsBlocks(kVectors, "/y");
  const zBlock = arrFromHashedToBlock(z, masks);

  const closeEnough = && reduce approxEqual(yBlock, zBlock);
  writeln(closeEnough);
  if (!closeEnough) {
    var maxErrorCount = 10;
    for i in yBlock.domain {
      if !approxEqual(zBlock[i], yBlock[i]) && maxErrorCount > 0 {
        writeln("at ", i, ": ",
                zBlock[i], " (computed) != ",
                yBlock[i], " (expected); Δ = ",
                abs(zBlock[i] - yBlock[i]));
        maxErrorCount -= 1;
      }
    }
  }
  assert(closeEnough);
  writeln(timer.elapsed());
}
