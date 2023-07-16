use LatticeSymmetries;
use HashedToBlock;
use BlockToHashed;
use Time;

config const kHamiltonian = "data/heisenberg_chain_10.yaml";

proc main() {
  initRuntime();
  defer deinitRuntime();

  const (_, matrix) = loadConfigFromYaml(kHamiltonian, hamiltonian=true);

  const masks;
  const basisStates = enumerateStates(matrix.basis, masks);

  const statesBlock = arrFromHashedToBlock(basisStates, masks);

  const statesHashed = arrFromBlockToHashed(statesBlock, masks);

  forall loc in Locales do on loc {
    const ref expected = basisStates.getBlock(loc.id);
    const ref predicted = statesHashed.getBlock(loc.id);
    const closeEnough = && reduce (expected == predicted);
    if (!closeEnough) {
      var maxErrorCount = 10;
      for i in expected.domain {
        if expected[i] != predicted[i] && maxErrorCount > 0 {
          writeln("at ", i, ": ",
                  predicted[i], " (computed) != ",
                  expected[i], " (expected)");
          maxErrorCount -= 1;
        }
      }
    }
    assert(closeEnough);
  }
}
