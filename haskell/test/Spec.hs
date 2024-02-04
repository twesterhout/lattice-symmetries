module Main (main) where

import LatticeSymmetries.AutomorphismsSpec qualified
import LatticeSymmetries.BenesSpec qualified
import LatticeSymmetries.BitStringSpec qualified
import LatticeSymmetries.ComplexRationalSpec qualified
import LatticeSymmetries.DenseSpec qualified
import LatticeSymmetries.ExprSpec qualified
import LatticeSymmetries.GeneratorSpec qualified
import LatticeSymmetries.GroupSpec qualified
import LatticeSymmetries.LoweringSpec qualified
import LatticeSymmetries.NonbranchingTermSpec qualified
import LatticeSymmetries.OperatorSpec qualified
import LatticeSymmetries.ParserSpec qualified
import LatticeSymmetries.PermutationSpec qualified
import Test.Hspec

main :: IO ()
main = hspec $ do
  describe "ComplexRational" LatticeSymmetries.ComplexRationalSpec.spec
  describe "BitString" LatticeSymmetries.BitStringSpec.spec
  describe "DenseMatrix" LatticeSymmetries.DenseSpec.spec
  describe "Permutation" LatticeSymmetries.PermutationSpec.spec
  describe "Automorphisms" LatticeSymmetries.AutomorphismsSpec.spec
  describe "Benes" LatticeSymmetries.BenesSpec.spec
  describe "Generator" LatticeSymmetries.GeneratorSpec.spec
  describe "Parser" LatticeSymmetries.ParserSpec.spec
  describe "Expr" LatticeSymmetries.ExprSpec.spec
  describe "Lowering" LatticeSymmetries.LoweringSpec.spec
  describe "NonbranchingTerm" LatticeSymmetries.NonbranchingTermSpec.spec
  describe "Operator" LatticeSymmetries.OperatorSpec.spec
  describe "Group" LatticeSymmetries.GroupSpec.spec
