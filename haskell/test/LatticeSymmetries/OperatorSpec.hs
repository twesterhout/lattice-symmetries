{-# LANGUAGE OverloadedLists #-}
{-# LANGUAGE OverloadedRecordDot #-}
{-# LANGUAGE OverloadedStrings #-}

module LatticeSymmetries.OperatorSpec (spec) where

import Data.Text.IO qualified
import Data.Vector.Generic qualified as G
import LatticeSymmetries.Basis
import LatticeSymmetries.Benes (Permutation (unPermutation))
import LatticeSymmetries.Expr (Expr, mapIndices, mkExpr, simplifyExpr)
import LatticeSymmetries.Generator (ParticleTag (SpinTag, SpinfulFermionTag))
import LatticeSymmetries.Group (PermutationGroup (PermutationGroup), emptySymmetries, hypergraphAutomorphisms)
import LatticeSymmetries.Operator
import LatticeSymmetries.Utils (toPrettyText)
import Test.Hspec

spec :: Spec
spec = do
  describe "operatorToHypergraph" $ do
    it "handles simple hopping terms" $ do
      let basis = SpinfulFermionBasis 2 SpinfulNoOccupation
          expr = mkExpr SpinfulFermionTag "- (c†₁↑ c₀↑ + c†₀↑ c₁↑ + c†₁↓ c₀↓ + c†₀↓ c₁↓)"
          operator = Operator basis expr
          (PermutationGroup automorphisms) = operatorSymmetryGroup operator
      G.forM_ automorphisms $ \p ->
        (applyPermutation p operator).opTerms `shouldBe` expr
    it "not all hypergraph symmetries are Hamiltonian symmetries" $ do
      let basis = SpinBasis 2 Nothing Nothing emptySymmetries
          expr = mkExpr SpinTag "σ⁺₀ σ⁻₁"
          operator = Operator basis expr
          (PermutationGroup automorphisms) = operatorSymmetryGroup operator
      G.forM_ automorphisms $ \p ->
        (applyPermutation p operator).opTerms `shouldBe` expr
