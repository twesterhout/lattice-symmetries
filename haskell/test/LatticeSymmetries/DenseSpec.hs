{-# LANGUAGE RankNTypes #-}

module LatticeSymmetries.DenseSpec (spec) where

import Data.Vector (Vector)
import qualified Data.Vector as B
import LatticeSymmetries.Dense
import Test.Hspec

spec :: Spec
spec = do
  describe "IsList DenseMatrix" $ do
    it "converts lists to matrices" $ do
      denseMatrixFromList [[1, 2], [3, 4], [5, 6]] `shouldBe` DenseMatrix 3 2 (B.generate 6 (+ 1))
      denseMatrixFromList [] `shouldBe` DenseMatrix 0 0 (B.empty :: Vector Int)
      denseMatrixFromList [[], []] `shouldBe` DenseMatrix 2 0 (B.empty :: Vector Int)
      denseMatrixFromList [[1, 2, 3]] `shouldBe` DenseMatrix 1 3 (B.generate 3 (+ 1))
    it "converts matrices to lists" $ do
      denseMatrixToList (DenseMatrix 3 2 (B.generate 6 (+ 1))) `shouldBe` [[1, 2], [3, 4], [5, 6]]
      denseMatrixToList (DenseMatrix 0 0 B.empty) `shouldBe` ([] :: [[Int]])
      denseMatrixToList (DenseMatrix 2 0 B.empty) `shouldBe` ([[], []] :: [[Int]])
      denseMatrixToList (DenseMatrix 1 3 (B.generate 3 (+ 1))) `shouldBe` [[1, 2, 3]]
