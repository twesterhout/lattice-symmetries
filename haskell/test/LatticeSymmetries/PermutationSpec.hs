{-# LANGUAGE OverloadedLists #-}

module LatticeSymmetries.PermutationSpec (spec) where

import Data.Aeson qualified as Aeson
import Data.Bits (Bits (bit, setBit, testBit, zeroBits, (.&.)))
import Data.Vector qualified as B
import Data.Vector.Generic qualified as G
import LatticeSymmetries.Permutation
import Test.Hspec
import Test.Hspec.QuickCheck (modifyMaxSize, prop)
import Utils (decodeResult)
import Data.ByteString (useAsCString)
import Control.Exception (bracket)
import LatticeSymmetries.Utils (ls_hs_destroy_string, decodeCString)
import LatticeSymmetries.FFI (withCobject)
import LatticeSymmetries.Context

spec :: Spec
spec = do
  describe "mkPermutation" $ do
    it "constructs permutations" $ do
      mkPermutation [0] `shouldSatisfy` isRight
      mkPermutation [0, 1, 2, 3] `shouldSatisfy` isRight
    it "recognizes invalid permutations" $ do
      mkPermutation [1, 2, 3, 2] `shouldSatisfy` isLeft
      mkPermutation [1, 2, 4, 5] `shouldSatisfy` isLeft
    it "does not support zero-length permutations" $ mkPermutation [] `shouldSatisfy` isLeft
    it "does not support permutations starting with 1" $ mkPermutation [1, 2, 3] `shouldSatisfy` isLeft
  describe "getPeriodicity" $ do
    modifyMaxSize (const 50) $ prop "computes periodicities" $ \p -> do
      let n = getPeriodicity p
      B.findIndex isIdentityPermutation (B.iterateN n (<> p) p) `shouldBe` Just (n - 1)
  describe "permuteVector" $ prop "permutes vector elements" $ \p ->
    permuteVector p (B.generate p.length id) `shouldBe` G.convert (unPermutation p)
  describe "permuteBits" $ prop "permutes bits" $ \p x0 -> do
    let x = x0 .&. bit (p.length - 1)
        v = B.generate p.length (testBit x)
        v' = permuteVector p v
        x' = G.ifoldl' (\acc i -> bool acc (setBit acc i)) zeroBits v'
    permuteBits p x `shouldBe` x'
  describe "identityPermutation" $ prop "isIdentityPermutation" $ \n ->
    when (n > 0) $ identityPermutation n `shouldSatisfy` isIdentityPermutation
  describe "permutationFromMappings" $ prop "round trips" $ \p -> do
    let mappingsAll = G.toList . G.map (uncurry Mapping) $ G.indexed (unPermutation p)
        mappings = filter (\(Mapping a b) -> a /= b) mappingsAll
    permutationFromMappings (Just p.length) mappings `shouldBe` p
    permutationFromMappings Nothing mappingsAll `shouldBe` p
  describe "minimalSupport" $ do
    it "computes something"
      $ minimalSupport
      <$> mkPermutation [0, 2, 1]
      `shouldBe` Right (Just 1)
    prop "computes minimal support" $ \p -> do
      let mappings = G.toList . G.map (uncurry Mapping) . G.filter (uncurry (/=)) $ G.indexed (unPermutation p)
      case mappings of
        [] -> minimalSupport p `shouldBe` Nothing
        (Mapping i _ : _) -> minimalSupport p `shouldBe` Just i
  describe "ToJSON/FromJSON Permutation" $ do
    prop "round trips" $ \(p :: Permutation) ->
      Aeson.decode (Aeson.encode p) `shouldBe` Just p
  describe "FFI" $ do
    prop "ls_hs_create_permutation / ls_hs_destroy_permutation / ls_hs_permutation_info" $ \(p :: Permutation) ->
      useAsCString (toStrict (Aeson.encode p)) $ \pEncoded ->
        decodeResult (ls_hs_create_permutation pEncoded) ls_hs_destroy_permutation $ \cPermutation -> do
          -- testing round trip
          withCobject @Cpermutation cPermutation (`shouldBe` p)
          -- testing PermutationInfo
          bracket (ls_hs_permutation_info cPermutation) ls_hs_destroy_string $ \cStr ->
            decodeCString cStr `shouldReturn` Right (PermutationInfo (unPermutation p) (getPeriodicity p))
