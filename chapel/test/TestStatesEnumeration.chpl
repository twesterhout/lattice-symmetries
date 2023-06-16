use LatticeSymmetries;
use BlockToHashed;
use HashedToBlock;
use MyHDF5;
use Time;

config const kHamiltonian = "data/heisenberg_kagome_12.yaml";
config const kRepresentatives = "data/heisenberg_kagome_12.h5";

proc main() {
  initRuntime();
  defer deinitRuntime();

  const basis = loadConfigFromYaml(kHamiltonian);

  var timer = new stopwatch();
  timer.start();
  const masks;
  const basisStates = enumerateStates(basis, masks);
  timer.stop();

  const predicted = arrFromHashedToBlock(basisStates, masks);
  const reference =
    readDatasetAsBlocks(kRepresentatives, "/representatives",
                        rank = 1, eltType = uint(64));

  const theSame = && reduce [i in reference.domain] reference[i] == predicted[i];
  if !theSame {
    var maxErrorCount = 10;
    for i in reference.domain {
      if reference[i] != predicted[i] && maxErrorCount > 0 {
        writeln("at index ", i, ": ", reference[i], " != ", predicted[i]);
        maxErrorCount -= 1;
      }
    }
  }
  writeln(timer.elapsed());
  return 0;
}
