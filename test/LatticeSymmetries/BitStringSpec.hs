module LatticeSymmetries.BitStringSpec (spec) where

import qualified Data.Vector.Storable as S
import qualified Data.Vector.Storable.Mutable as SM
import LatticeSymmetries.BitString
import System.IO.Unsafe (unsafePerformIO)
import Test.Hspec
import Test.Hspec.QuickCheck
import Test.QuickCheck (Arbitrary (..))
import Prelude hiding (words)

log2 :: Integer -> Int
log2 !n
  | n == 1 = 0
  | n > 1 = 1 + log2 (div n 2)
  | otherwise = error "invalid n"

numberWords :: BitString -> Int
numberWords (BitString n)
  | n == 0 = 0
  | otherwise = 1 + log2 n

instance Arbitrary BitString where
  arbitrary = (BitString . abs) <$> arbitrary

spec :: Spec
spec = do
  describe "readBitString <-> writeBitString" $ do
    prop "writes to vector & reads from memory" $ \x ->
      let n = numberWords x
          words = wordsFromBitString n x
          y = unsafePerformIO $
            S.unsafeWith words $ \wordsPtr ->
              readBitString n wordsPtr
       in y `shouldBe` x
    prop "writes to memory & reads from memory" $ \x ->
      let n = numberWords x
          y = unsafePerformIO $ do
            buf <- SM.new n
            SM.unsafeWith buf $ \bufPtr ->
              writeBitString n bufPtr x
            SM.unsafeWith buf $ \bufPtr ->
              readBitString n bufPtr
       in y `shouldBe` x
