use Utils;
use BatchedOperator;
use ForeignTypes;
use MatrixVectorProduct;

use CTypes;
import OS.POSIX;

config const kRun = true;

proc test_BatchedOperator() {
  const op = new Operator(
      "{"
    + "\"basis\": { \"number_spins\": 2 },"
    + "\"expression\": { \"expression\": \"2 (σ⁺₀ σ⁻₁ + σ⁺₁ σ⁻₀) + σᶻ₀ σᶻ₁\" }"
    + "}");
  // 1
  //  -1   2
  //   2  -1
  //         1

  if kRun {
    var alphas : [0 ..# 2] uint(64) = [0, 1];
    var betas : chpl_external_array;
    var coeffs : chpl_external_array;
    var offsets : chpl_external_array;
    ls_chpl_operator_apply_off_diag(op.payload, alphas.size, c_ptrToConst(alphas),
                                    c_ptrTo(betas), c_ptrTo(coeffs), c_ptrTo(offsets));
    defer chpl_free_external_array(betas);
    defer chpl_free_external_array(coeffs);
    defer chpl_free_external_array(offsets);

    checkArraysEqual(makeArrayFromExternArray(betas, uint(64)), [2:uint]);
    checkArraysEqual(makeArrayFromExternArray(coeffs, complex(128)), [2:complex(128)]);
    checkArraysEqual(makeArrayFromExternArray(offsets, int(64)), [0, 0, 1]);
  }
}

proc test_smallMatrixVectorProduct() {
  const op = new Operator(
      "{"
    + "\"basis\": { \"number_spins\": 2 },"
    + "\"expression\": { \"expression\": \"2 (σ⁺₀ σ⁻₁ + σ⁺₁ σ⁻₀) + σᶻ₀ σᶻ₁\" }"
    + "}");
  // 1
  //  -1   2
  //   2  -1
  //         1
  var representatives : [0 ..# 4] uint(64) = [0, 1, 2, 3]:uint;
  var x : [0 ..# 4] real(64) = [0, 1, 0, 0];
  var y : [0 ..# 4] real(64);

  perLocaleMatrixVector(op, x, y, representatives);
  writeln(y);
}

proc test_smallMatrixVectorProduct2() {
  const op = new Operator(
      "{"
    + "\"basis\": { \"number_spins\": 4, \"hamming_weight\": 2 },"
    + "\"expression\": { \"expression\": \"2 (σ⁺₀ σ⁻₁ + σ⁺₁ σ⁻₀) + σᶻ₀ σᶻ₁ + 2 (σ⁺₁ σ⁻₂ + σ⁺₂ σ⁻₁) + σᶻ₁ σᶻ₂ + 2 (σ⁺₂ σ⁻₃ + σ⁺₃ σ⁻₂) + σᶻ₂ σᶻ₃ + 2 (σ⁺₃ σ⁻₀ + σ⁺₀ σ⁻₃) + σᶻ₃ σᶻ₀\" }"
    + "}");
  var representatives : [0 ..# 6] uint(64) = [3, 5, 6, 9, 10, 12]:uint;
  var x : [0 ..# 6] real(64) = [0, 1, 0, 0, 0, 0];
  var y : [0 ..# 6] real(64);

  perLocaleMatrixVector(op, x, y, representatives);
  writeln(y);
}

proc test_smallMatrixVectorProduct3() {
  const op = new Operator(
      "{"
    + "\"basis\": { \"number_spins\": 4, \"hamming_weight\": 2 },"
    + "\"expression\": { \"expression\": \"2 (σ⁺₀ σ⁻₁ + σ⁺₁ σ⁻₀) + σᶻ₀ σᶻ₁ + 2 (σ⁺₁ σ⁻₂ + σ⁺₂ σ⁻₁) + σᶻ₁ σᶻ₂ + 2 (σ⁺₂ σ⁻₃ + σ⁺₃ σ⁻₂) + σᶻ₂ σᶻ₃ + 2 (σ⁺₃ σ⁻₀ + σ⁺₀ σ⁻₃) + σᶻ₃ σᶻ₀\" }"
    + "}");
  var representatives : [0 ..# 6] uint(64) = [3, 5, 6, 9, 10, 12]:uint;
  var x : [0 ..# 6] real(64) = [0, 1, 0, 0, 0, 0];
  var y : [0 ..# 6] real(64);

  perLocaleMatrixVector(op, x, y, representatives);
  writeln(y);
}

proc main() {
  initRuntime();
  defer deinitRuntime();

  test_BatchedOperator();
  test_smallMatrixVectorProduct();
  // test_smallMatrixVectorProduct2();

  return 0;
}
