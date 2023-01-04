{-# LANGUAGE OverloadedLists #-}

module LatticeSymmetries.ParserSpec (spec) where

import Data.Ratio ((%))
import LatticeSymmetries.Algebra
import LatticeSymmetries.ComplexRational
import LatticeSymmetries.Generator
import LatticeSymmetries.Parser
import Test.Hspec
import Text.Parsec (parse)

spec :: Spec
spec = do
  describe "pSpinOperator" $ do
    it "parses simple operators" $ do
      parse pSpinOperator "" ("σᶻ₁₀" :: Text)
        `shouldBe` Right (Expr [Scaled 1 [Generator 10 SpinZ]])
      parse pSpinOperator "" ("S⁺₀" :: Text)
        `shouldBe` Right (Expr [Scaled (fromRational (1 % 2)) [Generator 0 SpinPlus]])
    it "parses composite operators" $ do
      parse pSpinOperator "" ("σˣ₅" :: Text)
        `shouldBe` Right (Expr [Scaled 1 [Generator 5 SpinPlus], Scaled 1 [Generator 5 SpinMinus]])
      parse pSpinOperator "" ("σʸ₈₉₄₃" :: Text)
        `shouldBe` Right
          ( Expr
              [ Scaled (ComplexRational 0 (-1)) [Generator 8943 SpinPlus],
                Scaled (ComplexRational 0 1) [Generator 8943 SpinMinus]
              ]
          )
    it "parses scaled operators" $ do
      parse pNumber "" ("(4.5 + 3ⅈ)" :: Text)
        `shouldBe` Right (ComplexRational (9 % 2) 3)
      parse pNumber "" ("4.5" :: Text)
        `shouldBe` Right (ComplexRational (9 % 2) 0)
      parse pNumber "" ("5" :: Text)
        `shouldBe` Right (ComplexRational 5 0)
      parse (pOperatorString pSpinOperator) "" ("4.5 σᶻ₁₀" :: Text)
        `shouldBe` Right (Expr [Scaled (fromRational (9 % 2)) [Generator 10 SpinZ]])
      parse (pOperatorString pSpinOperator) "" ("4.5  × σᶻ₁₀" :: Text)
        `shouldBe` Right (Expr [Scaled (fromRational (9 % 2)) [Generator 10 SpinZ]])
      parse (pOperatorString pSpinOperator) "" ("4.5×σᶻ₁₀" :: Text)
        `shouldBe` Right (Expr [Scaled (fromRational (9 % 2)) [Generator 10 SpinZ]])
      parse (pOperatorString pSpinOperator) "" ("(4.5 + 3ⅈ) σᶻ₁₀ σ⁺₂₀" :: Text)
        `shouldBe` Right
          ( Expr
              [ Scaled
                  (ComplexRational (9 % 2) 3)
                  [Generator 10 SpinZ, Generator 20 SpinPlus]
              ]
          )
      parse (pExpr pSpinOperator) "" ("4.5 σᶻ₁₀ + σ⁺₂₀" :: Text)
        `shouldBe` Right
          ( Expr
              [ Scaled (ComplexRational (9 % 2) 0) [Generator 10 SpinZ],
                Scaled 1 [Generator 20 SpinPlus]
              ]
          )
      parse (pExpr pSpinOperator) "" ("(σ⁺₂₀ + σ⁻₈) σᶻ₅" :: Text)
        `shouldBe` Right
          ( Expr
              [ Scaled 1 [Generator 5 SpinZ, Generator 8 SpinMinus],
                Scaled 1 [Generator 5 SpinZ, Generator 20 SpinPlus]
              ]
          )
      parse (pExpr pSpinOperator) "" ("4.5 σᶻ₁₀ - (σ⁺₂₀ + (8 - 2ⅈ) σ⁻₈) σᶻ₅" :: Text)
        `shouldBe` Right
          ( Expr
              [ Scaled (ComplexRational (-8) 2) [Generator 5 SpinZ, Generator 8 SpinMinus],
                Scaled (-1) [Generator 5 SpinZ, Generator 20 SpinPlus],
                Scaled (ComplexRational (9 % 2) 0) [Generator 10 SpinZ]
              ]
          )
      mkExpr' SpinTag "σᶻ₉ σᶻ₁₀ + σ⁺₅ σ⁺₆"
        `shouldBe` Right
          ( Expr
              [ Scaled 1 [Generator 5 SpinPlus, Generator 6 SpinPlus],
                Scaled 1 [Generator 9 SpinZ, Generator 10 SpinZ]
              ]
          )
  it "" $ do
    replicateSiteIndices [[0, 1], [1, 2], [2, 3]] ("σ⁺₁ σ⁻₂" :: Expr 'SpinTy)
      `shouldBe` "σ⁺₀ σ⁻₁ + σ⁺₁ σ⁻₂ + σ⁺₂ σ⁻₃"
