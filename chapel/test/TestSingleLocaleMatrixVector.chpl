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
  var states: list(uint);
  var x_real : list(real);
  var x_imag : list(real);
  var y_real : list(real);
  var y_imag : list(real);
}

/*
proc unpackToFullBasis(const ref basis : Basis,
                       const ref basisStates : [] uint(64),
                       const ref norms : [] uint(16),
                       const ref xs : [] ?eltType) {

  const numberBits = basis.info.number_bits;
  const spinInversion = basis.info.spin_inversion;
  assert(basis.info.is_state_index_identity);

  const numberTargetStates = 1:uint(64) << numberBits;
  var ys : [0 ..# numberTargetStates] eltType;

  const groupSize = if spinInversion != 0 then 2 else 1;
  const spinInversionMask = (1:uint(64) << numberBits) - 1;

  forall (basisState, norm, x, i) in zip(basisStates, norms, xs, 0..) {
    ys[basisState:int] = x * math.sqrt(norm / groupSize);
    if spinInversion != 0 {
      const invertedBasisState = basisState ^ spinInversionMask;
      ys[invertedBasisState:int] = x * spinInversion * math.sqrt(norm / groupSize);
    }
  }

  return ys;
}
*/


proc main() {
  initRuntime();
  defer deinitRuntime();

  const testData = open(kVectors, ioMode.r).reader(deserializer=new jsonDeserializer()).read(TestInfo);
  const x = [(re, im) in zip(testData.x_real, testData.x_imag)] re + im * 1.0i;
  const y = [(re, im) in zip(testData.y_real, testData.y_imag)] re + im * 1.0i;

  const matrix = loadOperatorFromFile(kHamiltonian);
  const ref basis = matrix.basis;

  const basisStates, norms, keys;
  enumerateStates(basis, basisStates, norms, keys);
  if kVerbose {
    writeln("basisStates = ", basisStates[here]);
    writeln("expectedBasisStates = ", testData.states.toArray());
  }
  checkArraysEqual(testData.states.toArray(), basisStates[here]);
  matrix.basis.payload.deref().local_representatives = convertToExternalArray(basisStates[here]);
  matrix.basis.payload.deref().local_norms = convertToExternalArray(norms[here]);

  var z : [x.domain] x.eltType;
  perLocaleMatrixVector(matrix, x, z, basisStates[here]);
  if kVerbose {
    writeln("x = ", x);
    writeln("y = ", y);
    writeln("z = ", z);
  }
  checkArraysEqual(z, y);
}
