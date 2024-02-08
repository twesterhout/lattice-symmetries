{-# LANGUAGE OverloadedRecordDot #-}
{-# LANGUAGE OverloadedStrings #-}

module LatticeSymmetries.OperatorSpec (spec) where

import Control.Exception (bracket)
import Data.Aeson qualified as Aeson
import Data.ByteString (packCString)
import LatticeSymmetries.Basis
import LatticeSymmetries.Expr
import LatticeSymmetries.Generator
import LatticeSymmetries.Group
import LatticeSymmetries.Operator
import LatticeSymmetries.Some
import LatticeSymmetries.Utils
import Test.Hspec
import Utils

spec :: Spec
spec = do
  describe "ls_hs_create_operator / ls_hs_destroy_operator" $ do
    it "constructs ls_hs_operator" $ do
      let basis = SpinBasis 5 (Just 2) Nothing emptyRepresentation
          expr = either error id $ mkExpr SpinTag "σ⁺₀ σ⁻₁ σ⁺₂ σ⁻₃ σ⁺₄"
          operator = Operator basis expr
      bracket (newCbasis basis) ls_hs_destroy_basis $ \cBasis ->
        bracket (newCexpr (SomeExpr SpinTag expr)) ls_hs_destroy_expr $ \cExpr ->
          decodeResult (ls_hs_create_operator cBasis cExpr) ls_hs_destroy_operator $ \cOperator -> do
            withCoperator cOperator (`shouldBe` SomeOperator SpinTag operator)
  -- ls_hs_operator_max_number_off_diag cOperator `shouldReturn` 1
  describe "ls_hs_operator_from_json / ls_hs_operator_to_json" $ do
    it "roundtrips" $ do
      let basis = SpinBasis 5 (Just 2) Nothing emptyRepresentation
          expr = either error id $ mkExpr SpinTag "σ⁺₀ σ⁻₁ σ⁺₂ σ⁻₃ σ⁺₄"
          operator = Operator basis expr

      bracket (newCencoded (SomeOperator SpinTag operator)) ls_hs_destroy_string $ \jsonString ->
        decodeResult (ls_hs_operator_from_json jsonString) ls_hs_destroy_operator $ \cOperator ->
          withCoperator cOperator (`shouldBe` SomeOperator SpinTag operator)

      bracket (newCoperator Nothing Nothing (SomeOperator SpinTag operator)) ls_hs_destroy_operator $ \cOperator ->
        bracket (ls_hs_operator_to_json cOperator) ls_hs_destroy_string $ \jsonString ->
          decodeUtf8 @Text <$> packCString jsonString
            `shouldReturn` (decodeUtf8 . toStrict) (Aeson.encode (SomeOperator SpinTag operator))

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
