{-# LANGUAGE OverloadedRecordDot #-}
{-# LANGUAGE OverloadedStrings #-}

module LatticeSymmetries.OperatorSpec (spec) where

import Data.Vector.Generic qualified as G
import LatticeSymmetries.Automorphisms
import LatticeSymmetries.Basis
import LatticeSymmetries.Expr (mkExpr)
import LatticeSymmetries.Generator (ParticleTag (SpinTag, SpinfulFermionTag))
import LatticeSymmetries.Group
import LatticeSymmetries.Operator
import Test.Hspec
import Utils

spec :: Spec
spec = do
  describe "operatorToHypergraph" $ pure ()

-- it "handles simple hopping terms" $ do
--   expr <- extractRight $ mkExpr SpinfulFermionTag "- (c†₁↑ c₀↑ + c†₀↑ c₁↑ + c†₁↓ c₀↓ + c†₀↓ c₁↓)"
--   let basis = SpinfulFermionBasis 2 SpinfulNoOccupation
--   let operator = Operator basis expr
--       (PermutationGroup automorphisms) = operatorSymmetryGroup operator
--   G.forM_ automorphisms $ \p -> (applyPermutation p operator).opTerms `shouldBe` expr
-- it "not all hypergraph symmetries are Hamiltonian symmetries" $ do
--   expr <- extractRight $ mkExpr SpinTag "σ⁺₀ σ⁻₁"
--   let basis = SpinBasis 2 Nothing Nothing emptySymmetries
--       operator = Operator basis expr
--       (PermutationGroup automorphisms) = operatorSymmetryGroup operator
--   G.forM_ automorphisms $ \p -> (applyPermutation p operator).opTerms `shouldBe` expr
