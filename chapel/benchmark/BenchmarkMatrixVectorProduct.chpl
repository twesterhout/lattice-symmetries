use LatticeSymmetries;
import Random;

config const kHamiltonian = "data/heisenberg_chain_10.yaml";
config const kRandomSeed = 42;

proc main() {
  initRuntime();
  defer deinitRuntime();

  const (_, matrix) = loadConfigFromYaml(kHamiltonian, hamiltonian=true);
  const masks;
  const basisStates = enumerateStates(matrix.basis, masks);

  var x = new BlockVector(real(64), 1, basisStates.counts);
  coforall loc in Locales with (ref x) do on loc {
    ref myX = x.getBlock(loc.id);
    Random.fillRandom(myX, seed = kRandomSeed + loc.id);
  }

  var z = similar(x);

  matrixVectorProduct(matrix, x, z, basisStates);
}
