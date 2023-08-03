module LatticeSymmetries.ConversionSpec (spec) where

import Data.Text qualified as Text
import Data.Text.IO qualified as Text
import LatticeSymmetries.Basis
import LatticeSymmetries.Conversion
import LatticeSymmetries.Expr
import LatticeSymmetries.Generator
import LatticeSymmetries.Group
import LatticeSymmetries.Operator
import LatticeSymmetries.Parser
import LatticeSymmetries.Utils (toPrettyText)
import Test.Hspec

spec :: Spec
spec = do
  describe "Conversion" $ do
    -- it "writeTransFile" $ do
    --   let expr = mkExpr SpinfulFermionTag "c†↑₁ c↑₀ + c†↑₀ c↑₁ + 0.2 n↑₃ n↑₄ + n↑₃ n↓₃ + c†↑₁ c↑₀ c↑₂ c†↑₃"
    --   -- isHermitianExpr expr `shouldBe` True
    --   Text.putStrLn $ toPrettyText expr
    --   let (transfer, others) = extractTransFile expr
    --       (coulomb, others') = extractCoulombIntraFile others
    --       (twoBody, others'') = extractInterAllFile others'
    --   Text.putStrLn transfer
    --   Text.putStrLn coulomb
    --   Text.putStrLn twoBody
    --   Text.putStrLn $ toPrettyText others'
    --   True `shouldBe` True
    -- it "spinsToFermions" $ do
    --   let expr = mkExpr SpinTag "σˣ₀ σˣ₁ + σʸ₀ σʸ₁ + σᶻ₀ σᶻ₁"
    --   Text.putStrLn $ toPrettyText (spinsToFermions expr)
    it "prepareHPhi" $ do
      let
        -- basis = SpinBasis 3 (Just 1) Nothing emptySymmetries
        basis = SpinfulFermionBasis 2 (SpinfulTotalParticles 2)
        -- expr = spinsToFermions $ mkExpr SpinTag "σˣ₀ σˣ₁ + σʸ₀ σʸ₁ + σᶻ₀ σᶻ₁"
        -- expr = spinsToFermions $ mkExpr SpinTag "σˣ₀ σˣ₁ + σʸ₀ σʸ₁ + σᶻ₀ σᶻ₁ + σˣ₁ σˣ₂ + σʸ₁ σʸ₂ + σᶻ₁ σᶻ₂ + σˣ₀ σˣ₂ + σʸ₀ σʸ₂ + σᶻ₀ σᶻ₂"
        -- expr = mkExpr SpinfulFermionTag " c†↑₁ c↑₀ - c†↑₀ c↑₁ + 4.0 n↑₀ n↓₀ + 4.0 n↑₁ n↓₁"
        expr = mkExpr SpinfulFermionTag "- 2 (c†₁↑ c₀↑ + c†₀↑ c₁↑ + c†₁↓ c₀↓ + c†₀↓ c₁↓) + 4.0 n₀↑ n₀↓ + 4.0 n₁↑ n₁↓"
        (hphi, others) = extractInteractions expr
      isHermitianExpr expr `shouldBe` True
      Text.putStrLn (toPrettyText expr)
      Text.putStrLn (toPrettyText others)
      prepareHPhi basis hphi "/tmp/hubbard"
    it "preparemVMC" $ do
      let
        basis = SpinfulFermionBasis 2 (SpinfulTotalParticles 2)
        expr = mkExpr SpinfulFermionTag "- 2 (c†₁↑ c₀↑ + c†₀↑ c₁↑ + c†₁↓ c₀↓ + c†₀↓ c₁↓) + 4.0 n₀↑ n₀↓ + 4.0 n₁↑ n₁↓"
        (int, _) = extractInteractions expr
      isHermitianExpr expr `shouldBe` True
      -- Text.putStrLn (toPrettyText others)
      prepareVMC basis int "/tmp/vmc"

-- let
--   -- basis = SpinBasis 3 (Just 1) Nothing emptySymmetries
--   basis = SpinfulFermionBasis 2 (SpinfulTotalParticles 2)
--   -- expr = spinsToFermions $ mkExpr SpinTag "σˣ₀ σˣ₁ + σʸ₀ σʸ₁ + σᶻ₀ σᶻ₁"
--   -- expr = spinsToFermions $ mkExpr SpinTag "σˣ₀ σˣ₁ + σʸ₀ σʸ₁ + σᶻ₀ σᶻ₁ + σˣ₁ σˣ₂ + σʸ₁ σʸ₂ + σᶻ₁ σᶻ₂ + σˣ₀ σˣ₂ + σʸ₀ σʸ₂ + σᶻ₀ σᶻ₂"
--   -- expr = mkExpr SpinfulFermionTag " c†↑₁ c↑₀ - c†↑₀ c↑₁ + 4.0 n↑₀ n↓₀ + 4.0 n↑₁ n↓₁"
--   expr = mkExpr SpinfulFermionTag "- 2 (c†₁↑ c₀↑ + c†₀↑ c₁↑ + c†₁↓ c₀↓ + c†₀↓ c₁↓) + 4.0 n₀↑ n₀↓ + 4.0 n₁↑ n₁↓"
--   (hphi, others) = extractHPhi basis expr
-- isHermitianExpr expr `shouldBe` True
-- Text.putStrLn (toPrettyText expr)
-- Text.putStrLn (toPrettyText others)
-- prepareHPhi hphi "/tmp/hubbard"
