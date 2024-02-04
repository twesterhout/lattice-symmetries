{-# LANGUAGE OverloadedLists #-}

module LatticeSymmetries.AutomorphismsSpec (spec) where

import Data.Set qualified as Set
import LatticeSymmetries.Automorphisms
import Test.Hspec
import Utils ()

-- | c n is the cyclic graph on n vertices
cyclicGraph :: Int -> Hypergraph Int
cyclicGraph n
  | n >= 3 = Hypergraph (Set.fromList vs) (Set.fromList (Set.fromList <$> es))
  | otherwise = error "n must be at least 3"
  where
    vs = [0 .. n - 1]
    es = [[i, (i + 1) `mod` n] | i <- [0 .. n - 1]]

-- cyclic hypergraph with 3-vertex edges
cyclicGraph3 :: Int -> Hypergraph Int
cyclicGraph3 n
  | n >= 3 = Hypergraph (Set.fromList vs) (Set.fromList (Set.fromList <$> es))
  | otherwise = error "n must be at least 3"
  where
    vs = [0 .. (n - 1)]
    es = [[i, (i + 1) `mod` n, (i + 2) `mod` n] | i <- [0 .. (n - 1)]]

rectangularGraph :: Int -> Int -> Hypergraph Int
rectangularGraph n k = Hypergraph (Set.fromList vs) (Set.fromList (Set.fromList <$> es))
  where
    -- vs contains an extra element that does not enter any edge in es
    vs = [0 .. n * k]
    es = [[k * i + j, k * i + ((j + 1) `mod` k)] | i <- [0 .. n - 1], j <- [0 .. k - 1]] ++ [[k * i + j, k * ((i + 1) `mod` n) + j] | i <- [0 .. n - 1], j <- [0 .. k - 1]]

spec :: Spec
spec = do
  describe "hypergraphAutomorphisms" $ do
    it "works for cyclic graphs" $ do
      forM_ ([3 .. 10] :: [Int]) $ \n ->
        (hypergraphAutomorphisms (const True) (cyclicGraph n)).size `shouldBe` 2 * n
      forM_ ([50, 100] :: [Int]) $ \n ->
        (hypergraphAutomorphisms (const True) (cyclicGraph n)).size `shouldBe` 2 * n
    it "works for cyclic graphs with 3-vertex edges" $ do
      forM_ ([5 .. 10] :: [Int]) $ \n ->
        (hypergraphAutomorphisms (const True) (cyclicGraph3 n)).size `shouldBe` 2 * n
      forM_ ([50, 100] :: [Int]) $ \n ->
        (hypergraphAutomorphisms (const True) (cyclicGraph3 n)).size `shouldBe` 2 * n
    it "works for rectangularGraphs" $ do
      (hypergraphAutomorphisms (const True) (rectangularGraph 4 4)).size `shouldBe` 384
      (hypergraphAutomorphisms (const True) (rectangularGraph 3 4)).size `shouldBe` 48
