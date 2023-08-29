use LatticeSymmetries;
import Random;
import Time;


config const kHamiltonian = "data/heisenberg_chain_10.yaml";
config const kRandomSeed = 42;
config const kRunMatrixVectorProduct = true;

proc main() {
  initRuntime();
  defer deinitRuntime();

  const (_, matrix) = loadConfigFromYaml(kHamiltonian, hamiltonian=true);

  const masks;
  const basisStates = enumerateStates(matrix.basis, masks);
  const dimension = + reduce basisStates.counts;
  logDebug("Hilbert space dimension: ", dimension);
  logDebug("masks.size(): ", masks.localSubdomain().size);

  var timer = new Time.stopwatch();
  timer.start();
  var x = new BlockVector(real(64), 1, basisStates.counts);
  coforall loc in Locales with (ref x) do on loc {
    ref myX = x.getBlock(loc.id);
    Random.fillRandom(myX, seed = kRandomSeed + loc.id);
  }
  var z = similar(x);
  timer.stop();
  logDebug("Allocation and filling: ", timer.elapsed());

  if kRunMatrixVectorProduct then
    matrixVectorProduct(matrix, x, z, basisStates);
}
