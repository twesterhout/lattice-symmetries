import Timing;
use CSR;
use ForeignTypes;
use MatrixVectorProduct;
use StatesEnumeration;
use Utils;

use CTypes;
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

proc main() {
  initRuntime();
  defer deinitRuntime();

  const matrix = loadOperatorFromFile(kHamiltonian);
  const ref basis = matrix.basis;
  // writeln((new OperatorForeignInterface()).to_json(matrix.payload));

  const basisStates, norms, keys;
  enumerateStates(basis, basisStates, norms, keys);

  const testData = open(kVectors, ioMode.r).reader(deserializer=new jsonDeserializer()).read(TestInfo);
  const x = [(re, im) in zip(testData.x_real, testData.x_imag)] re + im * 1.0i;
  const y = [(re, im) in zip(testData.y_real, testData.y_imag)] re + im * 1.0i;
  var z : [x.domain] x.eltType;

  const csrMatrix = convertOffDiagToCsr(matrix, complex(128), basisStates[here], norms[here]);
  const diag = extractDiag(matrix, complex(128), basisStates[here], norms[here]);

  // const capacity = basisStates[here].size * (matrix.max_number_off_diag_estimate + 1);
  // writeln(makeArrayFromPtr(csrMatrix.matrixElements:c_ptr(complex(128)), {0 ..# capacity}));
  // writeln(makeArrayFromPtr(csrMatrix.colIndices:c_ptr(int(64)), {0 ..# capacity}));
  // writeln(makeArrayFromPtr(csrMatrix.rowOffsets:c_ptr(int(64)), {0 ..# (basisStates[here].size + 1)}));

  csrMatvec(csrMatrix, c_ptrToConst(x), c_ptrTo(z));
  z += diag * x;

  if kVerbose {
    writeln("basisStates = ", basisStates[here]);
    writeln("expectedBasisStates = ", testData.states.toArray());
    writeln("x = ", x);
    writeln("y = ", y);
    writeln("z = ", z);
  }

  checkArraysEqual(testData.states.toArray(), basisStates[here]);
  if matrix.basis.info.spin_inversion == -1 { // we do not trust QuSpin when spin_inversion==-1
    const absZ = abs(z);
    const absY = abs(y);
    checkArraysEqual(absZ, absY);
  }
  else {
    checkArraysEqual(z, y);
  }

}
