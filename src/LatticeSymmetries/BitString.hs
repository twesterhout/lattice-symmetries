module LatticeSymmetries.BitString
  ( BitString (..),
    readBitString,
    writeBitString,
    wordsFromBitString,
  )
where

import Control.Monad.Primitive
import Control.Monad.ST
import Data.Bits
import qualified Data.Primitive.Ptr as P
import Data.Vector.Generic (Vector)
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Generic.Mutable as GM
import Foreign.Ptr

newtype BitString = BitString Integer
  deriving stock (Show, Eq, Ord)
  deriving newtype (Bits)

writeBitString :: PrimMonad m => Int -> Ptr Word64 -> BitString -> m ()
writeBitString numberWords p (BitString x₀) = do
  let go !i !x
        | i < numberWords = do
          P.writeOffPtr p i (fromIntegral x :: Word64)
          go (i + 1) (x `shiftR` 64)
        | otherwise = pure ()
  go 0 x₀

readBitString :: PrimMonad m => Int -> Ptr Word64 -> m BitString
readBitString numberWords p = do
  let go !i !x
        | i < numberWords = do
          y <- BitString <$> (fromIntegral :: Word64 -> Integer) <$> P.readOffPtr p i
          go (i + 1) ((y `shiftL` (64 * i)) .|. x)
        | otherwise = pure x
  go 0 zeroBits

wordsFromBitString :: Vector v Word64 => Int -> BitString -> v Word64
wordsFromBitString numberWords (BitString x₀) = runST $ do
  mvector <- GM.unsafeNew numberWords
  let go !i !x
        | i < numberWords = do
          GM.unsafeWrite mvector i (fromIntegral x :: Word64)
          go (i + 1) (x `shiftR` 64)
        | otherwise = pure ()
  go 0 x₀
  G.unsafeFreeze mvector
