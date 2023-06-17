module LatticeSymmetries.GeneratorSpec (spec) where

import qualified Data.Aeson as Aeson
import LatticeSymmetries.Generator
import Test.Hspec

spec :: Spec
spec = do
  describe "FromJSON ParticleTy" $
    it "parses ParticleTy" $ do
      forM_ ([SpinTy, SpinlessFermionTy, SpinfulFermionTy] :: [ParticleTy]) $ \tp ->
        Aeson.decode (Aeson.encode tp) `shouldBe` Just tp
      Aeson.decode "\"spin\"" `shouldBe` Just SpinTy
      Aeson.decode "\"spinless-fermion\"" `shouldBe` Just SpinlessFermionTy
      Aeson.decode "\"spinful-fermion\"" `shouldBe` Just SpinfulFermionTy
