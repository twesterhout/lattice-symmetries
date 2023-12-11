{-# LANGUAGE OverloadedLists #-}

module LatticeSymmetries.BenesSpec (spec) where

import Data.Aeson qualified as Aeson
import Data.Bits
import Data.Vector (Vector)
import Data.Vector.Generic qualified as G
import LatticeSymmetries.Benes
import Test.Hspec
import Test.Hspec.QuickCheck
import Test.QuickCheck

data PermutationTestData = PermutationTestData !Permutation !(Vector Integer)
  deriving stock (Show)

instance Arbitrary PermutationTestData where
  arbitrary = sized $ \size -> do
    n <- chooseInt (1, 500)
    p <- either error id . mkPermutation . G.fromList <$> shuffle [0 .. n - 1]
    v <- G.replicateM size (chooseInteger (0, bit n - 1))
    pure $ PermutationTestData p v

spec :: Spec
spec = do
  describe "permuteBits" $ do
    it "should permute bits" $ do
      let p₁ = permutationToBenesNetwork <$> mkPermutation [0, 1, 2]
      benesPermuteBits <$> p₁ <*> pure 0b100 `shouldBe` Right 0b100
      let p₂ = permutationToBenesNetwork <$> mkPermutation [1, 2, 0]
      benesPermuteBits <$> p₂ <*> pure 0b100 `shouldBe` Right 0b010
      benesPermuteBits <$> p₂ <*> pure 0b101 `shouldBe` Right 0b110
      let v₃ = mkPermutation [1, 2, 3, 0, 5, 6, 7, 4, 9, 10, 11, 8]
          p₃ = permutationToBenesNetwork <$> v₃
      benesPermuteBits <$> p₃ <*> pure 0b100111010010 `shouldBe` Right 0b110011100001
      benesPermuteBits <$> p₃ <*> pure 0b100111010010 `shouldBe` permuteBits <$> v₃ <*> pure 0b100111010010
    prop "correctly permutes bits" $ \(PermutationTestData p v) ->
      let !network = permutationToBenesNetwork p
       in G.forM_ v $ \x ->
            benesPermuteBits network x `shouldBe` permuteBits p x
  describe "ToJSON/FromJSON Permutation" $ do
    it "encodes Permutation" $ do
      Aeson.encode (either error id $ mkPermutation [0, 1, 2, 3]) `shouldBe` "[0,1,2,3]"
      Aeson.encode (either error id $ mkPermutation [3, 5, 0, 1, 2, 4]) `shouldBe` "[3,5,0,1,2,4]"
    it "parses Permutation" $ do
      Aeson.decode "[3, 0, 2, 1]" `shouldBe` either error Just (mkPermutation [3, 0, 2, 1])
      Aeson.decode "[3, 0, 2, 1, 4, 5, 6]" `shouldBe` either error Just (mkPermutation [3, 0, 2, 1, 4, 5, 6])
      Aeson.decode "{ \"oops\": [2, 1, 0] }" `shouldBe` (Nothing :: Maybe Permutation)
    prop "encoding roundtrips" $ \(PermutationTestData p _) ->
      Aeson.decode (Aeson.encode p) `shouldBe` Just p
