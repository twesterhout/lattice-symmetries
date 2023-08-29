{-# LANGUAGE OverloadedLists #-}
{-# LANGUAGE OverloadedRecordDot #-}

module LatticeSymmetries.ConversionSpec (spec) where

import Data.Text.IO qualified as Text
import LatticeSymmetries.Basis
import LatticeSymmetries.Conversion
import LatticeSymmetries.Expr
import LatticeSymmetries.Generator
import LatticeSymmetries.Group
import LatticeSymmetries.Utils
import Test.Hspec

spec :: Spec
spec = do
  describe "Conversion" $ do
    it "prepareHPhi" $ do
      let
        basis = SpinBasis 4 (Just 2) Nothing emptySymmetries
        expr =
          spinsToFermions $
            mkExpr SpinTag "σˣ₀ σˣ₁ + σʸ₀ σʸ₁ + σᶻ₀ σᶻ₁ + σˣ₁ σˣ₂ + σʸ₁ σʸ₂ + σᶻ₁ σᶻ₂ + σˣ₂ σˣ₃ + σʸ₂ σʸ₃ + σᶻ₂ σᶻ₃ + σˣ₃ σˣ₀ + σʸ₃ σʸ₀ + σᶻ₃ σᶻ₀"
        (hphi, others) = extractInteractions expr
      isHermitianExpr expr `shouldBe` True
      Text.putStrLn $ toPrettyText others
      Text.putStrLn $ show hphi.exchange
      others `shouldBe` Expr []
      prepareHPhi "/tmp/lattice-symmetries-haskell/hphi" basis hphi
    it "preparemVMC" $ do
      let
        basis = SpinfulFermionBasis 2 (SpinfulTotalParticles 2)
        expr = mkExpr SpinfulFermionTag "- 2 (c†₁↑ c₀↑ + c†₀↑ c₁↑ + c†₁↓ c₀↓ + c†₀↓ c₁↓) + 4.0 n₀↑ n₀↓ + 4.0 n₁↑ n₁↓"
        (int, others) = extractInteractions expr
      isHermitianExpr expr `shouldBe` True
      others `shouldBe` Expr []
      prepareVMC "/tmp/lattice-symmetries-haskell/vmc" basis int
