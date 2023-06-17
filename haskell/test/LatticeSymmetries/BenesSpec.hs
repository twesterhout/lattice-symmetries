{-# LANGUAGE OverloadedLists #-}

module LatticeSymmetries.BenesSpec (spec) where

import qualified Data.Aeson as Aeson
import Data.Bits
import Data.Vector (Vector)
import qualified Data.Vector.Generic as G
import LatticeSymmetries.Benes
import Test.Hspec
import Test.Hspec.QuickCheck
import Test.QuickCheck

data PermutationTestData = PermutationTestData !Permutation !(Vector Integer)
  deriving stock (Show)

instance Arbitrary PermutationTestData where
  arbitrary = sized $ \size -> do
    n <- chooseInt (1, 500)
    p <- mkPermutation . G.fromList <$> shuffle [0 .. n - 1]
    v <- G.replicateM size (chooseInteger (0, bit n - 1))
    pure $ PermutationTestData p v

checkPermuteBits :: PermutationTestData -> Expectation
checkPermuteBits (PermutationTestData p v) = do
  let !network = toBenesNetwork p
  G.forM_ v $ \x ->
    permuteBits network x `shouldBe` permuteBits' p x

spec :: Spec
spec = do
  describe "permuteBits" $ do
    it "should permute bits" $ do
      let p₁ = toBenesNetwork (mkPermutation [0, 1, 2])
      permuteBits p₁ 0b100 `shouldBe` 0b100
      let p₂ = toBenesNetwork (mkPermutation [1, 2, 0])
      permuteBits p₂ 0b100 `shouldBe` 0b010
      permuteBits p₂ 0b101 `shouldBe` 0b110
      let v₃ = mkPermutation [1, 2, 3, 0, 5, 6, 7, 4, 9, 10, 11, 8]
          p₃ = toBenesNetwork v₃
      permuteBits p₃ 0b100111010010 `shouldBe` 0b110011100001
      permuteBits p₃ 0b100111010010 `shouldBe` permuteBits' v₃ 0b100111010010
    modifyMaxSize (const 100) $ prop "correctly permutes bits" $ checkPermuteBits
  describe "ToJSON Permutation" $ do
    it "encodes Permutation" $ do
      Aeson.encode (mkPermutation [0, 1, 2, 3]) `shouldBe` "[0,1,2,3]"
      Aeson.encode (mkPermutation [3, 5, 0, 1, 2, 4]) `shouldBe` "[3,5,0,1,2,4]"
  describe "FromJSON Permutation" $ do
    it "parses Permutation" $ do
      Aeson.decode "[3, 0, 2, 1]" `shouldBe` Just (mkPermutation [3, 0, 2, 1])
      Aeson.decode "[3, 0, 2, 1, 4, 5, 6]" `shouldBe` Just (mkPermutation [3, 0, 2, 1, 4, 5, 6])
      Aeson.decode "{ \"oops\": [2, 1, 0] }" `shouldBe` (Nothing :: Maybe Permutation)
