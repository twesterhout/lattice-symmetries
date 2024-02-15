{-# LANGUAGE UndecidableInstances #-}

module Utils
  ( compareWith
  , shouldBeApprox
  , extractRight
  , decodeResult
  , shouldRoundTrip
  )
where

import Control.Exception (bracket, throwIO)
import Data.Aeson (FromJSON, ToJSON)
import Data.Aeson qualified as Aeson
import Data.Bits (Bits (bit))
import Data.Ratio ((%))
import Data.Vector.Generic qualified as G
import Foreign.C (CString)
import Foreign.Ptr (Ptr)
import GHC.Exts (IsList (..))
import GHC.Stack
import LatticeSymmetries.Algebra
import LatticeSymmetries.BitString
import LatticeSymmetries.ComplexRational
import LatticeSymmetries.Context
import LatticeSymmetries.Expr
import LatticeSymmetries.Generator
import LatticeSymmetries.Group
import LatticeSymmetries.Permutation
import LatticeSymmetries.Some
import LatticeSymmetries.Utils
import Test.HUnit.Lang (FailureReason (..), HUnitFailure (..))
import Test.Hspec
import Test.QuickCheck
import Prelude hiding (Product, Sum, toList)

compareWith :: (HasCallStack, Show a) => (a -> a -> Bool) -> a -> a -> Expectation
compareWith comparator result expected =
  unless (comparator result expected) $ do
    throwIO (HUnitFailure location $ ExpectedButGot Nothing expectedMsg actualMsg)
  where
    expectedMsg = show expected
    actualMsg = show result
    location = case reverse (toList callStack) of
      (_, loc) : _ -> Just loc
      [] -> Nothing

shouldBeApprox :: (Ord a, ApproxEq a, Show a) => a -> a -> Expectation
shouldBeApprox = compareWith (≈)

extractRight :: (HasCallStack, Show a, Show b) => Either a b -> IO b
extractRight (Right x) = pure x
extractRight x@(Left _) = throwIO (HUnitFailure location $ ExpectedButGot Nothing "Right" (show x))
  where
    location = case reverse (toList callStack) of
      (_, loc) : _ -> Just loc
      [] -> Nothing

instance Arbitrary Permutation where
  arbitrary = sized $ \k -> do
    n <- chooseInt (1, max 1 k)
    either error id . mkPermutation . G.fromList <$> shuffle [0 .. n - 1]

instance Arbitrary (RepElement Permutation) where
  arbitrary = do
    p <- arbitrary
    let n = getPeriodicity p
    s <- chooseInt (0, n - 1)
    pure $ RepElement p (s % n)

instance Arbitrary BitString where
  arbitrary = sized $ \k -> BitString . fromIntegral <$> chooseInteger (0, bit k - 1)

instance Arbitrary ComplexRational where
  arbitrary = ComplexRational <$> arbitrary <*> arbitrary

instance Arbitrary ParticleTy where
  arbitrary = oneof [pure SpinTy, pure SpinfulFermionTy, pure SpinlessFermionTy]

instance Arbitrary SpinIndex where
  arbitrary = chooseEnum (SpinUp, SpinDown)

instance Arbitrary SpinGeneratorType where
  arbitrary = chooseEnum (SpinIdentity, SpinMinus)

instance Arbitrary FermionGeneratorType where
  arbitrary = chooseEnum (FermionIdentity, FermionAnnihilate)

instance Arbitrary g => Arbitrary (Generator Int g) where
  arbitrary = sized $ \n ->
    Generator <$> chooseInt (0, max 0 (n - 1)) <*> arbitrary

instance Arbitrary g => Arbitrary (Generator (SpinIndex, Int) g) where
  arbitrary = sized $ \n ->
    Generator <$> ((,) <$> arbitrary <*> chooseInt (0, max 0 (n - 1))) <*> arbitrary

instance (Arbitrary c, Arbitrary g) => Arbitrary (Scaled c g) where
  arbitrary = Scaled <$> arbitrary <*> arbitrary

instance Arbitrary g => Arbitrary (Product g) where
  arbitrary = do
    k <- chooseInt (0, 5)
    Product . G.fromListN k <$> vector k

instance Arbitrary a => Arbitrary (NonEmpty a) where
  arbitrary = (:|) <$> arbitrary <*> arbitrary

instance Arbitrary g => Arbitrary (Sum g) where
  arbitrary = do
    k <- chooseInt (1, 10)
    Sum . G.fromListN k <$> vector k

instance (Arbitrary (IndexType t), Arbitrary (GeneratorType t), Arbitrary (Generator (IndexType t) (GeneratorType t))) => Arbitrary (Expr t) where
  arbitrary = Expr <$> arbitrary

instance Arbitrary SomeExpr where
  arbitrary =
    oneof
      [ SomeExpr SpinTag <$> arbitrary
      , SomeExpr SpinlessFermionTag <$> arbitrary
      , SomeExpr SpinfulFermionTag <$> arbitrary
      ]

shouldRoundTrip :: (ToJSON a, FromJSON a, Eq a, Show a) => a -> Expectation
shouldRoundTrip x = Aeson.decode (Aeson.encode x) `shouldBe` Just x

decodeResult :: forall a r. IsCobject a => IO CString -> (Ptr a -> IO ()) -> (Ptr a -> IO r) -> IO r
decodeResult action deleter f =
  bracket action ls_hs_destroy_string $ \cStr ->
    bracket (extractRight =<< extractRight =<< decodeCString @(Either Text (Ptr a)) cStr) deleter f
