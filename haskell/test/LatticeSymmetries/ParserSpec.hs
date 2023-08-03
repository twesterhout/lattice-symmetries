{-# LANGUAGE OverloadedStrings #-}

module LatticeSymmetries.ParserSpec (spec) where

import LatticeSymmetries.Expr
import LatticeSymmetries.Generator
import LatticeSymmetries.Utils
import Test.Hspec

spec :: Spec
spec = do
  describe "mkExpr" $ do
    it "handles simple operators" $ do
      toPrettyText (mkExpr SpinTag "σˣ₅") `shouldBe` "σ⁺₅ + σ⁻₅"
      toPrettyText (mkExpr SpinTag "σʸ₈₉₄₃") `shouldBe` "-ⅈ σ⁺₈₉₄₃ + ⅈ σ⁻₈₉₄₃"
      toPrettyText (mkExpr SpinTag "\\sigma^y_90") `shouldBe` "-ⅈ σ⁺₉₀ + ⅈ σ⁻₉₀"
      toPrettyText (mkExpr SpinfulFermionTag "4.0 n₀↑ n₀↓") `shouldBe` "4.0 n₀↑ n₀↓"
    it "handles complex expressions" $ do
      toPrettyText (mkExpr SpinfulFermionTag "5 (c^\\dagger_10\\up - c\\dagger2\\down) n3\\up")
        `shouldBe` "5.0 n₃↑ c†₁₀↑ - 5.0 n₃↑ c†₂↓"
      toPrettyText (mkExpr SpinlessFermionTag "3.25 (c†23-5(c12-n2))")
        `shouldBe` "16.25 n₂ - 16.25 c₁₂ + 3.25 c†₂₃"
