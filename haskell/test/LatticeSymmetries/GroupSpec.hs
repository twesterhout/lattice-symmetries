{-# LANGUAGE OverloadedLists #-}
{-# LANGUAGE OverloadedRecordDot #-}

module LatticeSymmetries.GroupSpec (spec) where

import Data.Aeson qualified as Aeson
import Data.Bits
import Data.Complex
import Data.Ratio
import Data.Vector (Vector)
import Data.Vector.Generic qualified as G
import LatticeSymmetries.Benes
import LatticeSymmetries.Group
import Test.Hspec
import Test.Hspec.QuickCheck
import Test.QuickCheck
import Utils

syms :: [Either Text Symmetry] -> IO Symmetries
syms xs = fmap compileGroupRepresentation . (extractRight . groupRepresentationFromGenerators) =<< sequence (extractRight <$> xs)


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
  prop "getCharacter" $ \(p :: Int) (q :: Int) ->
    when (q /= 0) $ do
      let x = abs p % abs q - fromIntegral (abs p `div` abs q)
          z = (Symmetry (identityPermutation 1) x).character
      realPart z `shouldBeApprox` cos (-2 * pi * realToFrac @_ @Double x)
      imagPart z `shouldBeApprox` sin (-2 * pi * realToFrac @_ @Double x)

  describe "FromJSON Symmetry" $ do
    it "parses Symmetry" $ do
      Aeson.decode "{\"permutation\": [1, 2, 3, 0], \"sector\": 2}"
        `shouldBe` Just (either error id $ mkSymmetry [1, 2, 3, 0] 2)
      Aeson.decode "{\"permutation\": [0, 1], \"sector\": 0}"
        `shouldBe` Just (either error id $ mkSymmetry [0, 1] 0)
  -- Empty permutations are not supported
  -- Aeson.decode "{\"permutation\": [], \"sector\": 0}" `shouldBe` Just (mkSymmetry [] 0)
  -- describe "FromJSON Symmetries" $ do
  --   it "parses Symmetries" $ do
  --     expected1 <- syms [mkSymmetry [1, 2, 0] 1]
  --     Aeson.decode "[{\"permutation\": [1, 2, 0], \"sector\": 1}]" `shouldBe` Just expected1
  --     expected2 <- syms [mkSymmetry [1, 2, 3, 0] 0, mkSymmetry [3, 2, 1, 0] 0]
  --     Aeson.decode
  --       "[{\"permutation\": [1, 2, 3, 0], \"sector\": 0}, {\"permutation\": [3, 2, 1, 0], \"sector\": 0}]"
  --       `shouldBe` Just expected2
  describe "mkSymmetries" $ do
    it "builds cycles" $ do
      g <- syms [mkSymmetry [1, 2, 3, 0, 5, 6, 7, 4, 9, 10, 11, 8, 13, 14, 15, 12] 1]
      g.symmCharactersReal `shouldBe` [1, 0, -1, 0]
      g.symmCharactersImag `shouldBe` [0, -1, 0, 1]
