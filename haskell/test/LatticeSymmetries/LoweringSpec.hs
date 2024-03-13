module LatticeSymmetries.LoweringSpec (spec) where

import Control.Exception (bracket)
import Data.Bits
import Data.Maybe (fromJust)
import Data.Ratio ((%))
import Data.Vector qualified as B
import Data.Vector.Storable qualified as S
import Data.Vector.Storable.Mutable qualified as SM
import LatticeSymmetries.Basis
import LatticeSymmetries.BitString
import LatticeSymmetries.Group
import LatticeSymmetries.Lowering
import LatticeSymmetries.Permutation
import LatticeSymmetries.Utils (prettyValidate)
import Test.Hspec
import Test.Hspec.QuickCheck
import Test.QuickCheck
import Test.QuickCheck.Gen (chooseInt64)
import Utils (extractRight)
import Prelude hiding (group)

data LinearChainRepresentation = LinearChainRepresentation !Int !(Representation Permutation)
  deriving stock (Show, Eq)

data LinearChainTestData = LinearChainTestData !Int !(Representation Permutation) !(S.Vector Word64)
  deriving stock (Show, Eq)

mkLinearChainRepresentation :: Int -> Int -> Either Text LinearChainRepresentation
mkLinearChainRepresentation n sector = do
  p <- mkPermutation . fromList $ ([1 .. n - 1] <> [0])
  t <- prettyValidate $ RepElement p (sector % n)
  s <- fromGenerators [t]
  pure $ LinearChainRepresentation n s

instance Arbitrary LinearChainRepresentation where
  arbitrary = do
    n <- chooseInt (2, 64)
    sector <- chooseInt (0, n - 1)
    pure . fromRight (error "should be Right") $ mkLinearChainRepresentation n sector

instance Arbitrary LinearChainTestData where
  arbitrary = do
    (LinearChainRepresentation n group) <- arbitrary
    basisStates <-
      if n < 6
        then pure [0 .. bit n - 1]
        else do
          k <- chooseInt (1, 64)
          replicateM k (chooseInt (0, bit n - 1))
    pure $ LinearChainTestData n group (fromList $ fromIntegral <$> basisStates)

data FixedHammingBitStrings = FixedHammingBitStrings !Int !Int !(S.Vector Int64) !(S.Vector Word64)
  deriving stock (Show, Eq)

instance Arbitrary FixedHammingBitStrings where
  arbitrary = do
    n <- chooseInt (1, 63)
    k <- chooseInt (0, n)
    batchSize <- chooseInt (1, 100)
    indices <- S.replicateM batchSize $ chooseInt64 (0, fromIntegral (fromJust (binomial n k)) - 1)
    let states = S.map (fixedHammingIndexToState k . fromIntegral) indices
    pure $ FixedHammingBitStrings n k indices states

referenceIsRepresentative :: Int -> Maybe Int -> Representation Permutation -> S.Vector Word64 -> S.Vector Word16
referenceIsRepresentative n spinInversion group = S.map (maybe 0 fromIntegral . isRepresentativeSlow group spinInversion . BasisState n . BitString . fromIntegral)

referenceStateInfo :: Int -> Maybe Int -> Representation Permutation -> S.Vector Word64 -> (S.Vector Word64, S.Vector Int32)
referenceStateInfo n spinInversion group = bimap S.convert S.convert . B.unzip . B.map f . B.convert
  where
    f :: Word64 -> (Word64, Int32)
    f w =
      let (BasisState _ (BitString r), i) = stateInfoSlow group spinInversion $ BasisState n (BitString (fromIntegral w))
       in (fromIntegral r, fromIntegral i)

-- referenceStateInfo :: Int -> Representation Permutation -> S.Vector Word64 -> (S.Vector Word64, S.Vector Int32)
-- referenceStateInfo n group states = let (rs, is) = unzip $ fmap f (S.toList states) in (S.fromList rs, S.fromList is)
--   where
--     f x =
--       let (BasisState _ (BitString rep), index) = stateInfoSlow group Nothing (BasisState n . BitString $ fromIntegral x)
--        in (fromIntegral rep, fromIntegral index)

spec :: Spec
spec = do
  describe "isRepresentative" $ do
    it "example 1" $ do
      forM_ [3, 4] $ \n ->
        forM_ [0, 1, 2] $ \k -> do
          (LinearChainRepresentation _ group) <- extractRight $ mkLinearChainRepresentation n k
          let basisStates :: S.Vector Word64
              basisStates = fromList [0 .. bit n - 1]
          forM_ [Nothing, Just 1, Just (-1)] $ \i ->
            bracket (createIsRepresentativeKernel_v2 group i) destroyIsRepresentativeKernel_v2 $ \kernelPtr ->
              invokeIsRepresentativeKernel kernelPtr basisStates
                `shouldReturn` referenceIsRepresentative n i group basisStates
    modifyMaxSuccess (const 20) $
      prop "small linear chains" $
        \(LinearChainTestData n group basisStates) -> do
          forM_ [Nothing, Just 1, Just (-1)] $ \i ->
            bracket (createIsRepresentativeKernel_v2 group i) destroyIsRepresentativeKernel_v2 $ \kernelPtr ->
              invokeIsRepresentativeKernel kernelPtr basisStates
                `shouldReturn` referenceIsRepresentative n i group basisStates

  describe "stateInfo_v3" $ do
    it "example 1" $ do
      (LinearChainRepresentation n group) <- extractRight (mkLinearChainRepresentation 3 0)
      let basisStates :: S.Vector Word64
          basisStates = fromList [0 .. bit n - 1]
      forM_ [Just 1 {-[Nothing, Just 1, Just (-1)]-}] $ \i ->
        bracket (createStateInfoKernel_v3 group i) destroyStateInfoKernel_v3 $ \kernelPtr ->
          invokeStateInfoKernel kernelPtr basisStates
            `shouldReturn` referenceStateInfo n i group basisStates

  -- modifyMaxSuccess (const 20) $
  --   prop "small linear chains" $
  --     \(LinearChainTestData n group basisStates) -> do
  --       forM_ [Nothing, Just 1, Just (-1)] $ \i ->
  --         bracket (createStateInfoKernel_v3 group i) destroyStateInfoKernel_v3 $ \kernelPtr ->
  --           invokeStateInfoKernel kernelPtr basisStates
  --             `shouldReturn` referenceStateInfo n i group basisStates

  -- describe "fixedHammingStateToIndex" $ do
  --   modifyMaxSuccess (const 10) $
  --     prop "computes indices correctly" $
  --       \(FixedHammingBitStrings numberSites hammingWeight indices states) -> do
  --         bracket
  --           (createFixedHammingStateToIndexKernel numberSites hammingWeight)
  --           destroyFixedHammingStateToIndexKernel
  --           $ \kernelPtr ->
  --             invokeFixedHammingStateToIndexKernel kernelPtr states `shouldReturn` indices

  describe "ls_hs_fixed_hamming_state_to_index" $ do
    prop "computes indices correctly" $
      \(FixedHammingBitStrings _numberSites _hammingWeight indices states) -> do
        out <- SM.new (S.length states)
        S.unsafeWith states $ \statesPtr ->
          SM.unsafeWith out $ \indicesPtr -> do
            ls_hs_fixed_hamming_state_to_index
              (fromIntegral (S.length states))
              statesPtr
              indicesPtr
            S.freeze out `shouldReturn` indices

  describe "ls_hs_fixed_hamming_index_to_state" $ do
    prop "computes indices correctly" $
      \(FixedHammingBitStrings numberSites hammingWeight indices states) -> do
        out <- SM.new (S.length states)
        S.unsafeWith indices $ \indicesPtr ->
          SM.unsafeWith out $ \outPtr -> do
            ls_hs_fixed_hamming_index_to_state
              (fromIntegral numberSites)
              (fromIntegral hammingWeight)
              (fromIntegral (S.length states))
              indicesPtr
              outPtr
            S.freeze out `shouldReturn` states
