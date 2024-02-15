{-# LANGUAGE OverloadedLists #-}

module LatticeSymmetries.ExprSpec (spec) where

import Control.Exception (bracket)
import Data.Aeson qualified as Aeson
import Data.ByteString (packCString, useAsCString)
import Foreign (fromBool)
import Foreign.C (CString)
import LatticeSymmetries.Algebra
import LatticeSymmetries.Automorphisms
import LatticeSymmetries.ComplexRational (ComplexRational (ComplexRational))
import LatticeSymmetries.Expr
import LatticeSymmetries.Generator
import LatticeSymmetries.NonbranchingTerm (HasNonbranchingRepresentation (nonbranchingRepresentation))
import LatticeSymmetries.Some
import LatticeSymmetries.Utils (decodeCString, ls_hs_destroy_string, toPrettyText)
import Test.Hspec
import Test.Hspec.QuickCheck
import Utils
import Prelude hiding (Product, Sum)

decodeCexprResult :: IO CString -> IO SomeExpr
decodeCexprResult action =
  bracket action ls_hs_destroy_string $ \cStr ->
    bracket (extractRight =<< extractRight =<< decodeCString @(Either Text _) cStr) ls_hs_destroy_expr $ foldCexpr pure

spec :: Spec
spec = do
  describe "ls_hs_expr_from_json / ls_hs_expr_to_json" $ do
    modifyMaxSize (const 100) $ prop "round trips" $ \(someExpr :: SomeExpr) ->
      useAsCString (toStrict (Aeson.encode someExpr)) $ \exprEncoded ->
        bracket (ls_hs_expr_from_json exprEncoded) ls_hs_destroy_string $ \cExprStr -> do
          cExpr <- extractRight =<< extractRight =<< decodeCString @(Either Text _) cExprStr
          bracket (ls_hs_expr_to_json cExpr) ls_hs_destroy_string $ \cStr -> do
            Just e <- Aeson.decodeStrict <$> packCString cStr
            -- when (mapSomeExpr simplifyExpr someExpr /= e) $ do
            --   Text.IO.putStrLn $ toPrettyText (mapSomeExpr simplifyExpr someExpr) <> " | " <> toPrettyText e
            e `shouldBe` mapSomeExpr simplifyExpr someExpr
  describe "ls_hs_expr_to_string" $ do
    prop "pretty prints" $ \(someExpr :: SomeExpr) ->
      bracket (newCexpr someExpr) ls_hs_destroy_expr $ \cExpr ->
        bracket (ls_hs_expr_to_string cExpr) ls_hs_destroy_string $ \cStr ->
          decodeUtf8 <$> packCString cStr `shouldReturn` toPrettyText someExpr
  describe "Num" $ do
    let check :: (IsBasis t, Show (Expr t)) => ParticleTag t -> Expr t -> Expr t -> Expr t -> Expectation
        check _tag a b c = do
          simplifyExpr (negate (simplifyExpr (negate a))) `shouldBe` simplifyExpr a
          simplifyExpr (a + b) `shouldBe` simplifyExpr (b + a)
          simplifyExpr (a - b) `shouldBe` simplifyExpr (negate (b - a))
          simplifyExpr (simplifyExpr (a + b) + c) `shouldBe` simplifyExpr (a + simplifyExpr (b + c))
          simplifyExpr (simplifyExpr (a + b) * c) `shouldBe` simplifyExpr (a * c + b * c)
    prop "SpinTy" $ check SpinTag
    prop "SpinlessFermionTy" $ check SpinlessFermionTag
    prop "SpinfulFermionTy" $ check SpinfulFermionTag
  describe "ls_hs_{plus,minus,times,scale}" $ do
    let check :: (IsBasis t, Show (Expr t)) => ParticleTag t -> Expr t -> Expr t -> Expectation
        check tag a b = do
          bracket (newCexpr (SomeExpr tag a)) ls_hs_destroy_expr $ \cA ->
            bracket (newCexpr (SomeExpr tag b)) ls_hs_destroy_expr $ \cB -> do
              ls_hs_expr_equal cA cA `shouldReturn` fromBool True
              ls_hs_expr_equal cB cB `shouldReturn` fromBool True
              decodeCexprResult (ls_hs_expr_plus cA cB) `shouldReturn` SomeExpr tag (simplifyExpr (a + b))
              decodeCexprResult (ls_hs_expr_minus cA cB) `shouldReturn` SomeExpr tag (simplifyExpr (a - b))
              decodeCexprResult (ls_hs_expr_times cA cB) `shouldReturn` SomeExpr tag (simplifyExpr (a * b))
              bracket (ls_hs_expr_negate cA) ls_hs_destroy_expr $
                foldCexpr (`shouldBe` SomeExpr tag (-a))
              bracket (ls_hs_expr_scale 2.0 (-5.0) cA) ls_hs_destroy_expr $
                foldCexpr (`shouldBe` SomeExpr tag (ComplexRational 2 (-5) `scale` a))
    prop "SpinTy" $ check SpinTag
    prop "SpinlessFermionTy" $ check SpinlessFermionTag
    prop "SpinfulFermionTy" $ check SpinfulFermionTag
  describe "exprPermutationGroup" $ do
    it "computes permutation groups" $ do
      let e = mkExpr SpinTag "2 (σ⁺₀ σ⁻₁ + σ⁺₁ σ⁻₀) + σᶻ₀ σᶻ₁"
      (.permutations) . exprPermutationGroup Nothing <$> e `shouldBe` Right [[0, 1], [1, 0]]
  describe "HasNonbranchingRepresentation" $ do
    prop "identities match" $ do
      nonbranchingRepresentation (Generator (1 :: Int) SpinIdentity)
        `shouldBe` nonbranchingRepresentation (Generator (2 :: Int) FermionIdentity)
  describe "isInvariantUponSpinInversion" $ do
    it "works for SpinTy" $ do
      isInvariantUponSpinInversion <$> mkExpr SpinTag "2 (σ⁺₀ σ⁻₁ + σ⁺₁ σ⁻₀) + σᶻ₀ σᶻ₁"
        `shouldBe` Right True
      isInvariantUponSpinInversion <$> mkExpr SpinTag "σ⁺₀ σ⁻₁"
        `shouldBe` Right False
      isInvariantUponSpinInversion <$> mkExpr SpinTag "σᶻ₁"
        `shouldBe` Right False
      isInvariantUponSpinInversion <$> mkExpr SpinTag "S^x_10"
        `shouldBe` Right True
  describe "conservesNumberParticles" $ do
    it "works for SpinTy" $ do
      conservesNumberParticles <$> mkExpr SpinTag "2 (σ⁺₀ σ⁻₁ + σ⁺₁ σ⁻₀) + σᶻ₀ σᶻ₁"
        `shouldBe` Right True
      conservesNumberParticles <$> mkExpr SpinTag "σᶻ₀ σᶻ₁"
        `shouldBe` Right True
      conservesNumberParticles <$> mkExpr SpinTag "2 (σ⁺₀ σ⁻₁ + σ⁺₁ σ⁻₀)"
        `shouldBe` Right True
      conservesNumberParticles <$> mkExpr SpinTag "σ⁺₀ σ⁺₁ + σ⁻₁ σ⁻₀"
        `shouldBe` Right False
