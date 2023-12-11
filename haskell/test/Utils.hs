{-# LANGUAGE DefaultSignatures #-}
{-# LANGUAGE OverloadedRecordDot #-}
{-# LANGUAGE UndecidableInstances #-}

module Utils
  ( compareWith
  , shouldBeApprox
  , extractRight
  )
where

import Control.Exception (throwIO)
import Control.Monad (unless, void)
import GHC.Exts (IsList (..))
import GHC.Stack
import LatticeSymmetries.Utils
import System.IO (stderr)
import Test.HUnit.Lang (FailureReason (..), HUnitFailure (..))
import Test.Hspec
import Prelude hiding (Eq (..), toList)

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
