{-# LANGUAGE RankNTypes #-}

module LatticeSymmetries.ComplexRationalSpec (spec) where

import LatticeSymmetries.ComplexRational
import LatticeSymmetries.Context (Cscalar)
import LatticeSymmetries.Utils
import Test.Hspec
import Test.Hspec.QuickCheck

checkOp :: (forall a. Num a => a -> a -> a) -> (Cscalar, Cscalar) -> Expectation
checkOp op (z₁, z₂) = toComplexDouble r' `shouldSatisfy` (≈ r)
  where
    z₁' = fromComplexDouble z₁ :: ComplexRational
    z₂' = fromComplexDouble z₂ :: ComplexRational
    r = z₁ `op` z₂
    r' = z₁' `op` z₂'

checkDivide :: (Cscalar, Cscalar) -> Expectation
checkDivide (z₁, z₂)
  | z₂ == 0 = pure ()
  | otherwise = toComplexDouble r' `shouldSatisfy` (≈ r)
  where
    z₁' = fromComplexDouble z₁ :: ComplexRational
    z₂' = fromComplexDouble z₂ :: ComplexRational
    r = z₁ / z₂
    r' = z₁' / z₂'

spec :: Spec
spec = do
  describe "Num ComplexRational" $ do
    prop "adds numbers" $ checkOp (+)
    prop "subtracts numbers" $ checkOp (-)
    prop "multiplies numbers" $ checkOp (*)
    prop "divides numbers" checkDivide
