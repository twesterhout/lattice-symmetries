module LatticeSymmetries.NonbranchingTermSpec (spec) where

import LatticeSymmetries.Generator
import LatticeSymmetries.NonbranchingTerm
import Test.Hspec

spec :: Spec
spec = do
  describe "Semigroup NonbranchingTerm" $ do
    it ".." $ do
      let a = nonbranchingRepresentation (Generator (3 :: Int) FermionCreate)
          b = nonbranchingRepresentation (Generator (4 :: Int) FermionAnnihilate)
          c = nonbranchingRepresentation (Generator (3 :: Int) FermionAnnihilate)
          d = nonbranchingRepresentation (Generator (4 :: Int) FermionIdentity)
      -- print a
      -- print b
      case a <> a of
        (NonbranchingTerm v _ _ _ _ _) -> v `shouldBe` 0
      case b <> b of
        (NonbranchingTerm v _ _ _ _ _) -> v `shouldBe` 0
      (a <> c) `shouldBe` (nonbranchingRepresentation (Generator (3 :: Int) FermionCount))
      (d <> d) `shouldBe` d
      (a <> d) `shouldBe` a
      (d <> a) `shouldBe` a
      (b <> d) `shouldBe` b
      (d <> b) `shouldBe` b

-- print (a <> b)
-- print (b <> a)
