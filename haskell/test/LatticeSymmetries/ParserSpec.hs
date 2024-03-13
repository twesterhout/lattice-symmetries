{-# LANGUAGE OverloadedStrings #-}

module LatticeSymmetries.ParserSpec (spec) where

import Control.Arrow qualified
import Data.Text.IO qualified
import LatticeSymmetries.ComplexRational (prettyRational)
import LatticeSymmetries.Expr
import LatticeSymmetries.Generator
import LatticeSymmetries.Parser
import LatticeSymmetries.Utils
import Prettyprinter hiding (space)
import Test.Hspec
import Test.Hspec.QuickCheck
import Text.Megaparsec (errorBundlePretty, parse)
import Text.Megaparsec.Char (space)
import Text.Megaparsec.Char.Lexer (signed)
import Utils

spec :: Spec
spec = do
  let testParse :: (Eq a, Show a) => Parser a -> Text -> a -> Expectation
      testParse p s y = case Control.Arrow.left (toText . errorBundlePretty) (parse p "" s) of
        Left e -> error e
        Right x -> x `shouldBe` y
  describe "pReal" $ do
    prop "parses real numbers" $ \x -> testParse (signed space pReal) (renderDoc (prettyRational x)) x
  describe "pComplexRational" $ do
    prop "parses complex numbers" $ \x -> testParse (signed space pCoeff) (toPrettyText x) x
  describe "mkExpr" $ do
    let p :: (Pretty a, Show a) => Either Text a -> IO Text
        p x@(Right _) = fmap toPrettyText . extractRight $ x
        p x@(Left e) = do
          Data.Text.IO.putStrLn e
          fmap toPrettyText . extractRight $ x
    it "handles simple operators" $ do
      p (mkExpr SpinTag "σˣ₅") `shouldReturn` "σ⁺₅ + σ⁻₅"
      p (mkExpr SpinTag "σʸ₈₉₄₃") `shouldReturn` "-ⅈ σ⁺₈₉₄₃ + ⅈ σ⁻₈₉₄₃"
      p (mkExpr SpinTag "\\sigma^y_90") `shouldReturn` "-ⅈ σ⁺₉₀ + ⅈ σ⁻₉₀"
      p (mkExpr SpinfulFermionTag "4.0 n₀↑ n₀↓") `shouldReturn` "4.0 n₀↑ n₀↓"
      p (mkExpr SpinfulFermionTag "4.0 n₀↑ I") `shouldReturn` "4.0 n₀↑"
      p (mkExpr SpinTag "4.0 I") `shouldReturn` "4.0 I"
    it "handles complex expressions" $ do
      p (mkExpr SpinfulFermionTag "5 (c^\\dagger_10\\up - c\\dagger2\\down) n3\\up") `shouldReturn` "5.0 n₃↑ c†₁₀↑ - 5.0 n₃↑ c†₂↓"
      p (mkExpr SpinlessFermionTag "3.25 (c†23-5(c12-n2))") `shouldReturn` "16.25 n₂ - 16.25 c₁₂ + 3.25 c†₂₃"
      p (mkExpr SpinfulFermionTag "- (c†₁↑ c₀↑ + c†₀↑ c₁↑ + c†₁↓ c₀↓ + c†₀↓ c₁↓)") `shouldReturn` "-c†₀↑ c₁↑ + c₀↑ c†₁↑ - c†₀↓ c₁↓ + c₀↓ c†₁↓"
      p (mkExpr SpinfulFermionTag "(3.0 + 3.0ⅈ) n₁↑ + (3.0 + 1.0ⅈ) n₁↑ n₂↓ c†₀↑ + (-2/3 - 7/3ⅈ) I + (1/3 + 0.5ⅈ) I n₁↑ I") `shouldReturn` "(-2/3 - 7/3ⅈ) I + (3.0 + 1.0ⅈ) c†₀↑ n₁↑ n₂↓ + (10/3 + 3.5ⅈ) n₁↑"
    it "handles complex numbers" $ do
      p (mkExpr SpinTag "j σ⁺₅") `shouldReturn` "ⅈ σ⁺₅"
      p (mkExpr SpinTag "-2im σ⁺₅") `shouldReturn` "-2.0ⅈ σ⁺₅"
      p (mkExpr SpinTag "(-3.1 + 2ⅈ) σ⁺₅") `shouldReturn` "(-3.1 + 2.0ⅈ) σ⁺₅"
      p (mkExpr SpinTag "(3.1 - 2ⅈ) σ⁺₅") `shouldReturn` "(3.1 - 2.0ⅈ) σ⁺₅"
