import HDF5Extensions;
use ForeignTypes;
use MatrixVectorProduct;
use StatesEnumeration;
use Utils;
import Timing;

use IO;
use JSON;

config const kHamiltonian = "data/heisenberg_chain_8.json";
config const kVectors = "data/matvec/heisenberg_chain_8.h5";

proc localLoadVectors(filename : string, x : string = "/x", y : string = "/y") {
  var input = HDF5Extensions.readDataset(filename, x, real(64), rank = 2)[0, ..];
  var output = HDF5Extensions.readDataset(filename, y, real(64), rank = 2)[0, ..];
  return (input, output);
}

proc loadOperator(filename : string) {
  var f = open(filename, ioMode.r);
  var r = f.reader();
  return new Operator(r.readAll(string));
}

proc main() {
  initRuntime();
  defer deinitRuntime();

  const matrix = loadOperator(kHamiltonian);
  const _k = ls_chpl_get_state_to_index_kernel(matrix.basis.payload);
  const (x, y) = localLoadVectors(kVectors);

  const basisStates;
  const norms;
  const keys;
  enumerateStates(matrix.basis, basisStates, norms, keys);

  var z : [x.domain] x.eltType;
  perLocaleMatrixVector(matrix, x, z, basisStates[here]);

  checkArraysEqual(z, y);

  stdout.withSerializer(new jsonSerializer()).writeln(Timing.summarize());

  // const counts = localProcessCounts.toArray();
  // writeln((+ reduce counts):real / counts.size:real);
  // writeln(min reduce counts);
  // writeln(max reduce counts);
}
