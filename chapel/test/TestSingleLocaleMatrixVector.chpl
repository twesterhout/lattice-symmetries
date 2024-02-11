import HDF5Extensions;
use ForeignTypes;
use MatrixVectorProduct;
use StatesEnumeration;
use Utils;
import Timing;

use IO;
use JSON;
use List;

config const kHamiltonian = "../test/test_4_2_10_10_expr.json";
config const kVectors = "../test/test_4_2_10_10_arrays.json";
config const kVerbose = false;

record TestInfo {
  var x_real : list(real);
  var x_imag : list(real);
  var y_real : list(real);
  var y_imag : list(real);
}

proc main() {
  initRuntime();
  defer deinitRuntime();

  const matrix = loadOperatorFromFile(kHamiltonian);
  // writeln((new OperatorForeignInterface()).to_json(matrix.payload));

  const basisStates, norms, keys;
  enumerateStates(matrix.basis, basisStates, norms, keys);

  const testData = open(kVectors, ioMode.r).reader(deserializer=new jsonDeserializer()).read(TestInfo);
  const x = [(re, im) in zip(testData.x_real, testData.x_imag)] re + im * 1.0i;
  const y = [(re, im) in zip(testData.y_real, testData.y_imag)] re + im * 1.0i;

  var z : [x.domain] x.eltType;
  perLocaleMatrixVector(matrix, x, z, basisStates[here]);

  if kVerbose {
    writeln(y);
    writeln(z);
  }

  checkArraysEqual(z, y);
}
