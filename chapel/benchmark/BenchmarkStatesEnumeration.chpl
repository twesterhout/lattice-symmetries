use LatticeSymmetries;

config const kHamiltonian = "data/heisenberg_chain_10.yaml";

proc main() {
  initRuntime();
  defer deinitRuntime();

  const (_, matrix) = loadConfigFromYaml(kHamiltonian, hamiltonian=true);
  const masks;
  const basisStates = enumerateStates(matrix.basis, masks);
  const dimension = + reduce basisStates.counts;
  logDebug("Hilbert space dimension: ", dimension);
}
