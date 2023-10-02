{-# LANGUAGE DefaultSignatures #-}
{-# LANGUAGE OverloadedRecordDot #-}
{-# LANGUAGE UndecidableInstances #-}

module Utils
  ( compareWith
  , shouldBeApprox
  )
where

import Control.Exception (throwIO)
import Control.Monad (unless, void)
import GHC.Exts (IsList (..))
import GHC.Stack
import System.IO (stderr)
import Test.HUnit.Lang (FailureReason (..), HUnitFailure (..))
import Test.Hspec
import Prelude hiding (Eq (..), toList)
import LatticeSymmetries.Utils

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
