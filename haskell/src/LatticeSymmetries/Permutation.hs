{-# LANGUAGE DefaultSignatures #-}
{-# LANGUAGE OverloadedRecordDot #-}
{-# LANGUAGE TemplateHaskell #-}

module LatticeSymmetries.Permutation
  ( Permutation (..)
  , mkPermutation
  , HasPeriodicity (..)
  , permuteVector
  , permuteBits
  , identityPermutation
  , isIdentityPermutation
  , Mapping (..)
  , permutationFromMappings
  , minimalSupport
  )
where

import Control.Monad.ST
import Data.Aeson (FromJSON (..), ToJSON (..))
import Data.Bits
import Data.IntSet qualified as IntSet
import Data.List qualified as L
import Data.Set qualified as Set
import Data.Validity (check)
import Data.Vector.Generic ((!))
import Data.Vector.Generic qualified as G
import Data.Vector.Unboxed qualified as U
import Data.Vector.Unboxed.Deriving (derivingUnbox)
import Data.Vector.Unboxed.Mutable qualified as UM
import GHC.Exts (IsList (..))
import GHC.Records (HasField (getField))
import LatticeSymmetries.BitString
import LatticeSymmetries.Utils
import Prelude hiding (cycle, identity)

-- | A permutation of numbers @[0 .. N - 1]@.
newtype Permutation = Permutation {unPermutation :: U.Vector Int}
  deriving stock (Show, Eq, Ord, Generic)
  deriving anyclass (NFData)

instance ToJSON Permutation where
  toJSON (Permutation p) = toJSON p

instance FromJSON Permutation where
  parseJSON = (eitherToParser . mkPermutation) <=< parseJSON

instance IsList Permutation where
  type Item Permutation = Int
  toList (Permutation p) = G.toList p
  fromList p = either error id $ mkPermutation (G.fromList p)

class HasPeriodicity a where
  getPeriodicity :: a -> Int
  default getPeriodicity :: (Semigroup a, Eq a) => a -> Int
  getPeriodicity g0 = go 1 (g0 <> g0)
    where
      go !i !curr
        | curr == g0 = i
        | otherwise = go (i + 1) (curr <> g0)

fastGetPeriodicity :: Permutation -> Int
fastGetPeriodicity (Permutation p)
  | G.null p = 0
  | otherwise = L.foldl1' lcm $ go (IntSet.fromDistinctAscList [1 .. G.length p - 1]) 1 0 0
  where
    go !todo !size !x0 !i
      | p G.! i == x0 =
          (size :) $
            if IntSet.null todo
              then []
              else
                let (x0', todo') = IntSet.deleteFindMin todo
                 in go todo' (1 :: Int) x0' x0'
      | otherwise =
          let i' = p G.! i
              todo' = IntSet.delete i' todo
           in go todo' (size + 1) x0 i'

instance HasPeriodicity Permutation where
  getPeriodicity = fastGetPeriodicity

-- | Rearrange elements of the input vector according to the given permutation.
--
-- Throws an error in case the length of the input vector does not match the length of the permutation.
permuteVector
  :: (HasCallStack, G.Vector v a)
  => Permutation
  -- ^ Specifies how to order elements
  -> v a
  -- ^ Input vector
  -> v a
  -- ^ Rearranged vector
permuteVector (Permutation p) xs
  | G.length p == G.length xs = G.generate (G.length p) (G.unsafeIndex xs . G.unsafeIndex p)
  | otherwise = error $ "length mismatch: " <> show (G.length p) <> " != " <> show (G.length xs)

-- | Rearrange bits of the input 'BitString' according to the given permutation.
permuteBits :: Permutation -> BitString -> BitString
permuteBits (Permutation p) x = go 0 zeroBits
  where
    go !i !y
      | i < G.length p =
          let y' = if testBit x (p ! i) then setBit y i else y
           in go (i + 1) y'
      | otherwise = y

-- | Get the length of the permutation. If we are given a permutation of numbers @[0 .. N-1]@, then
-- this function will return @N@.
instance HasField "length" Permutation Int where
  getField (Permutation p) = G.length p

-- | Generate the identity permutation of given length.
identityPermutation
  :: HasCallStack
  => Int
  -- ^ Length of the permutation
  -> Permutation
identityPermutation size
  | size > 0 = Permutation $ G.generate size id
  | otherwise = error $ "invalid size: " <> show size

-- | Create a permutation from vector.
mkPermutation :: U.Vector Int -> Either Text Permutation
mkPermutation p = prettyValidate (Permutation p)

instance Validity Permutation where
  validate (Permutation p) =
    mconcat
      [ check (not (G.null p)) "p is not empty"
      , check
          (Set.toAscList (Set.fromList (G.toList p)) == [0 .. G.length p - 1])
          (show p <> "is a permutation of [0 .. " <> show (G.length p - 1) <> "]")
      ]

instance Semigroup Permutation where
  (<>) x (Permutation ys) = Permutation $ permuteVector x ys

data Mapping a = Mapping !a !a
  deriving stock (Show, Eq, Ord)

derivingUnbox
  "Mapping"
  [t|forall a. UM.Unbox a => Mapping a -> (a, a)|]
  [|\(Mapping a b) -> (a, b)|]
  [|\(!a, !b) -> Mapping a b|]

-- | Create a permutation from a list of index mappings, e.g. [1->2, 3->5, 2->1].
-- Assumes that the list of mappings is well-formed.
permutationFromMappings :: HasCallStack => Maybe Int -> [Mapping Int] -> Permutation
permutationFromMappings maybeN mappings = runST $ do
  p <- UM.generate n id
  forM_ mappings $ \(Mapping x y) ->
    UM.write p x y
  either error id . mkPermutation <$> U.unsafeFreeze p
  where
    !n = fromMaybe ((+ 1) . L.maximum $ (\(Mapping x y) -> max x y) <$> mappings) maybeN

minimalSupport :: Permutation -> Maybe Int
minimalSupport = G.findIndex (uncurry (/=)) . G.indexed . unPermutation

isIdentityPermutation :: Permutation -> Bool
isIdentityPermutation = G.all (uncurry (==)) . G.indexed . unPermutation
