{-# LANGUAGE OverloadedLists #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE OverloadedRecordDot #-}

module LatticeSymmetries.GroupSpec (spec) where

import Data.Aeson qualified as Aeson
import Data.Bits
import Data.Complex
import Data.Ratio
import Data.Text.IO qualified
import Data.Vector (Vector)
import Data.Vector.Generic qualified as G
import LatticeSymmetries.Benes
import LatticeSymmetries.Expr (exprPermutationGroup, mkExpr, Expr, mapIndices)
import LatticeSymmetries.Generator
import LatticeSymmetries.Group
import LatticeSymmetries.Lowering
import LatticeSymmetries.Permutation
import LatticeSymmetries.Utils (toPrettyText)
import Test.Hspec
import Test.Hspec.QuickCheck
import Test.QuickCheck
import Utils

-- syms :: [Either Text Symmetry] -> IO Symmetries
-- syms xs = fmap (compileGroupRepresentation . unRepresentation) . (extractRight . groupRepresentationFromGenerators) =<< sequence (extractRight <$> xs)

reversePermutation :: Permutation -> Permutation
reversePermutation = Permutation . G.reverse . unPermutation

linearChainRepresentations :: Bool -> Int -> Vector (Representation Permutation)
linearChainRepresentations periodic size
  | not periodic =
      [ Representation [RepElement [0 .. size - 1] 0, RepElement (reversePermutation [0 .. size - 1]) 0]
      , Representation [RepElement [0 .. size - 1] 0, RepElement (reversePermutation [0 .. size - 1]) (1 % 2)]
      ]
  | odd size =
      let t = either error id . mkPermutation $ [1 .. size - 1] <> [0]
          p = either error id . mkPermutation $ G.reverse [0 .. size - 1]
       in [ either error id (fromGenerators [RepElement t 0, RepElement p 0])
          , either error id (fromGenerators [RepElement t 0, RepElement p (1 % 2)])
          ]
  | otherwise =
      let t = either error id . mkPermutation $ [1 .. size - 1] <> [0]
          p = either error id . mkPermutation $ G.reverse [0 .. size - 1]
       in [ either error id (fromGenerators [RepElement t 0, RepElement p 0])
          , either error id (fromGenerators [RepElement t (1 % 2), RepElement p 0])
          , either error id (fromGenerators [RepElement t (1 % 2), RepElement p (1 % 2)])
          , either error id (fromGenerators [RepElement t 0, RepElement p (1 % 2)])
          ]

heisenbergOnGraph :: Vector (Int, Int) -> Expr SpinTy
heisenbergOnGraph edges = G.foldl1' (+) (G.map replace edges)
  where
    expr = either error id $ mkExpr SpinTag "S+0 S-1 + S-0 S+1"
    replace (i, j) = flip mapIndices expr $ \case
      0 -> i
      1 -> j
      _ -> error "should not happen"

heisenbergOnChain :: Bool -> Int -> Expr SpinTy
heisenbergOnChain periodic size
  | size < 2 = error "not supported"
  | periodic && size > 2 = heisenbergOnGraph . G.fromList $ [(i, (i + 1) `mod` size) | i <- [0 .. size - 1]]
  | otherwise = heisenbergOnGraph . G.fromList $ [(i, i + 1) | i <- [0 .. size - 2]]

-- heisenbergOnSquare :: Bool -> Int -> Int -> Expr SpinTy
-- heisenbergOnSquare periodic width height
--   | width < 2 || height < 2 = error "not supported"
--   | periodic && width > 2 && height > 2 =
--       heisenbergOnGraph . G.fromList $ [(i, (i + 1) `mod` size) | i <- [0 .. size - 1]]
--   | otherwise = heisenbergOnGraph . G.fromList $ [(i, i + 1) | i <- [0 .. size - 2]]

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
          z = (RepElement (identityPermutation 1) x).character
      realPart z `shouldBeApprox` cos (-2 * pi * realToFrac @_ @Double x)
      imagPart z `shouldBeApprox` sin (-2 * pi * realToFrac @_ @Double x)

  describe "FromJSON Symmetry" $ do
    prop "round trips" $ shouldRoundTrip @Symmetry
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
  -- describe "mkSymmetries" $ do
  --   it "builds cycles" $ do
  --     g <- syms [mkSymmetry [1, 2, 3, 0, 5, 6, 7, 4, 9, 10, 11, 8, 13, 14, 15, 12] 1]
  --     g.symmCharactersReal `shouldBe` [1, 0, -1, 0]
  --     g.symmCharactersImag `shouldBe` [0, -1, 0, 1]

  describe "groupRepresentations" $ do
    prop "non-periodic chains" $ \size ->
      when (size >= 2) $ do
        let e = heisenbergOnChain False size
        groupRepresentations (exprPermutationGroup Nothing e) `shouldBe` linearChainRepresentations False size
    prop "periodic chains" $ \size ->
      when (size > 2) $ do
        let e = heisenbergOnChain True size
        groupRepresentations (exprPermutationGroup Nothing e) `shouldBe` linearChainRepresentations True size

    it "3-site non-periodic chain" $ do
      e <- extractRight $ mkExpr SpinTag "S+0 S-1 + S+1 S-0 + S+1 S-2 + S+2 S-1"
      groupRepresentations (exprPermutationGroup Nothing e)
        `shouldBe` [ Representation [RepElement [0, 1, 2] 0, RepElement [2, 1, 0] 0]
                   , Representation [RepElement [0, 1, 2] 0, RepElement [2, 1, 0] (1 % 2)]
                   ]
    it "3-site periodic chain" $ do
      e <- extractRight $ mkExpr SpinTag "S+0 S-1 + S+1 S-0 + S+1 S-2 + S+2 S-1 + S+2 S-0 + S-2 S+0"
      groupRepresentations (exprPermutationGroup Nothing e)
        `shouldBe` linearChainRepresentations True 3
