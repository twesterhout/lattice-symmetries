{-# LANGUAGE OverloadedLists #-}
{-# LANGUAGE OverloadedRecordDot #-}

module LatticeSymmetries.GroupSpec (spec) where

import qualified Data.Aeson as Aeson
import Data.Bits
import Data.Vector (Vector)
import qualified Data.Vector.Generic as G
import LatticeSymmetries.Benes
import LatticeSymmetries.Group
import Test.Hspec
import Test.Hspec.QuickCheck
import Test.QuickCheck

spec :: Spec
spec = do
  -- describe "mkSymmetries" $ do
  --   it "should build symmetry groups" $ do
  --     print $
  --       mkSymmetries
  --         [ mkSymmetry (mkPermutation [1, 2, 3, 0]) 0,
  --           mkSymmetry (mkPermutation [3, 2, 1, 0]) 1
  --         ]
  --     True `shouldBe` True
  describe "FromJSON Symmetry" $ do
    it "parses Symmetry" $ do
      Aeson.decode "{\"permutation\": [1, 2, 3, 0], \"sector\": 2}" `shouldBe` Just (mkSymmetry [1, 2, 3, 0] 2)
      Aeson.decode "{\"permutation\": [0, 1], \"sector\": 0}" `shouldBe` Just (mkSymmetry [0, 1] 0)
      -- Empty permutations are not supported
      -- Aeson.decode "{\"permutation\": [], \"sector\": 0}" `shouldBe` Just (mkSymmetry [] 0)
  describe "FromJSON Symmetries" $ do
    it "parses Symmetries" $ do
      Aeson.decode "[{\"permutation\": [1, 2, 0], \"sector\": 1}]"
        `shouldBe` Just ([mkSymmetry [1, 2, 0] 1] :: Symmetries)
      Aeson.decode
        "[{\"permutation\": [1, 2, 3, 0], \"sector\": 0}, \
        \ {\"permutation\": [3, 2, 1, 0], \"sector\": 0}]"
        `shouldBe` Just ([mkSymmetry [1, 2, 3, 0] 0, mkSymmetry [3, 2, 1, 0] 0] :: Symmetries)
  describe "mkSymmetries" $ do
    it "builds cycles" $ do
      let g1 = mkSymmetry [1,2,3,0,5,6,7,4,9,10,11,8,13,14,15,12] 1
      fmap (.symmCharactersReal) (mkSymmetries [g1]) `shouldBe` Right [1, 0, -1, 0]
      fmap (.symmCharactersImag) (mkSymmetries [g1]) `shouldBe` Right [0, -1, 0, 1]
