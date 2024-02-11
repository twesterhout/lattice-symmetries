module LatticeSymmetries.GeneratorSpec (spec) where

import Data.Aeson qualified as Aeson
import Data.Bits
import LatticeSymmetries.BitString
import LatticeSymmetries.Generator
import Test.Hspec
import Test.Hspec.QuickCheck
import Test.QuickCheck
import Utils

instance Arbitrary (BasisState t) where
  arbitrary = do
    n <- chooseInt (1, 500)
    bits <- fromInteger <$> chooseInteger (0, bit n - 1)
    pure $ BasisState n (BitString bits)

spec :: Spec
spec = do
  describe "ParticleTy" $ do
    it "FromJSON works" $ do
      Aeson.decode "\"spin\"" `shouldBe` Just SpinTy
      Aeson.decode "\"spinless-fermion\"" `shouldBe` Just SpinlessFermionTy
      Aeson.decode "\"spinful-fermion\"" `shouldBe` Just SpinfulFermionTy
  describe "ToJSON / FromJSON (ParticleTy, SpinGeneratorType, FermionGeneratorType)" $ do
    prop "ParticleTy" $ shouldRoundTrip @ParticleTy
  describe "BasisState" $ do
    prop "invertBasisState" $ \x@(BasisState n bits) -> do
      let y@(BasisState n' bits') = invertBasisState x
      n' `shouldBe` n
      (bits' Data.Bits..&. bits) `shouldBe` zeroBits
      invertBasisState y `shouldBe` x
