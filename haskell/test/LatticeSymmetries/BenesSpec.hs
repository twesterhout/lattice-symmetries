{-# LANGUAGE OverloadedLists #-}

module LatticeSymmetries.BenesSpec (spec) where

import Data.Aeson qualified as Aeson
import Data.Bits
import Data.Vector (Vector)
import Data.Vector.Generic qualified as G
import LatticeSymmetries.Benes
import LatticeSymmetries.BitString
import LatticeSymmetries.Permutation
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
      benesPermuteBits <$> p₁ <*> pure (BitString 0b100) `shouldBe` Right (BitString 0b100)
      let p₂ = permutationToBenesNetwork <$> mkPermutation [1, 2, 0]
      benesPermuteBits <$> p₂ <*> pure (BitString 0b100) `shouldBe` Right (BitString 0b010)
      benesPermuteBits <$> p₂ <*> pure (BitString 0b101) `shouldBe` Right (BitString 0b110)
      let v₃ = mkPermutation [1, 2, 3, 0, 5, 6, 7, 4, 9, 10, 11, 8]
          p₃ = permutationToBenesNetwork <$> v₃
      benesPermuteBits <$> p₃ <*> pure (BitString 0b100111010010) `shouldBe` Right (BitString 0b110011100001)
      benesPermuteBits <$> p₃ <*> pure (BitString 0b100111010010) `shouldBe` permuteBits <$> v₃ <*> pure (BitString 0b100111010010)
      let v₄ = mkPermutation [0]
          p₄ = permutationToBenesNetwork <$> v₄
      print p₄
      benesPermuteBits <$> p₄ <*> pure (BitString 0b1) `shouldBe` Right (BitString 0b1)
      benesPermuteBits <$> p₄ <*> pure (BitString 0b1) `shouldBe` permuteBits <$> v₄ <*> pure (BitString 0b1)
    prop "correctly permutes bits" $ \(PermutationTestData p v) ->
      let !network = permutationToBenesNetwork p
       in G.forM_ v $ \x ->
            benesPermuteBits network (BitString x) `shouldBe` permuteBits p (BitString x)
  describe "mkBatchedBenesNetwork" $ do
    it "should build a BatchedBenesNetwork" $ do
      let !p = permutationToBenesNetwork <$> mkPermutation [0]
      mkBatchedBenesNetwork . G.singleton <$> p `shouldBe` Right emptyBatchedBenesNetwork
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
