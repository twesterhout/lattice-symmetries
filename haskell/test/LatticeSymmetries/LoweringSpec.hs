module LatticeSymmetries.LoweringSpec (spec) where

import Control.Exception (bracket)
import Data.Bits
import Data.Maybe (fromJust)
import Data.Vector qualified as B
import Data.Vector.Generic qualified as G
import Data.Vector.Storable qualified as S
import Data.Vector.Storable.Mutable qualified as SM
import LatticeSymmetries.Basis
import LatticeSymmetries.BitString
import LatticeSymmetries.Dense (DenseMatrix (DenseMatrix))
import LatticeSymmetries.Generator
import LatticeSymmetries.Group
import LatticeSymmetries.Lowering
import LatticeSymmetries.Permutation
import Test.Hspec
import Test.Hspec.QuickCheck
import Test.QuickCheck
import Test.QuickCheck.Gen (chooseInt64)
import Prelude hiding (group)

data LinearChainRepresentation = LinearChainRepresentation !Int !(Representation Permutation)
  deriving stock (Show, Eq)

data LinearChainTestData = LinearChainTestData !Int !(Representation Permutation) !(S.Vector Word64)
  deriving stock (Show, Eq)

mkLinearChainRepresentation :: Int -> Int -> Either Text LinearChainRepresentation
mkLinearChainRepresentation n sector = do
  p <- mkPermutation . fromList $ ([1 .. n - 1] <> [0])
  t <- mkSymmetry p sector
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

referenceIsRepresentative :: Int -> Representation Permutation -> S.Vector Word64 -> S.Vector Double
referenceIsRepresentative n group = S.map (maybe 0 realToFrac . isRepresentativeSlow group Nothing . BasisState n . BitString . fromIntegral)

referenceStateInfo :: Int -> Representation Permutation -> S.Vector Word64 -> (S.Vector Word64, S.Vector Int32)
referenceStateInfo n group states = let (rs, is) = unzip $ fmap f (S.toList states) in (S.fromList rs, S.fromList is)
  where
    f x =
      let (BasisState _ (BitString rep), index) = stateInfoSlow group Nothing (BasisState n . BitString $ fromIntegral x)
       in (fromIntegral rep, fromIntegral index)

spec :: Spec
spec = do
  describe "isRepresentative" $ do
    it "example 1" $ do
      let (LinearChainRepresentation _ group) = fromRight (error "oops") $ mkLinearChainRepresentation 4 0
      let basisStates :: S.Vector Word64
          basisStates = fromList [0 .. bit 4 - 1]
          isRep :: S.Vector Double
          isRep = fromList [4, 1, 0, 1, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 4]
      referenceIsRepresentative 4 group basisStates `shouldBe` isRep
      kernelPtr <- createIsRepresentativeKernel_v2 group Nothing
      invokeIsRepresentativeKernel kernelPtr basisStates `shouldReturn` isRep
      destroyIsRepresentativeKernel_v2 kernelPtr
    -- TODO: compilation doesn't work for 1 spin
    -- it "example 2" $ do
    --   let (LinearChainRepresentation n group) = fromRight (error "oops") $ mkLinearChainRepresentation 1 0
    --   let basisStates :: S.Vector Word64
    --       basisStates = fromList [0 .. bit n - 1]
    --       isRep :: S.Vector Double
    --       isRep = fromList [1, 1]
    --   S.map (maybe 0 realToFrac . isRepresentativeSlow group Nothing . BasisState n . BitString . fromIntegral) basisStates `shouldBe` isRep
    --   kernelPtr <- createIsRepresentativeKernel_v2 group Nothing
    --   invokeIsRepresentativeKernel kernelPtr basisStates `shouldReturn` isRep
    it "example 3" $ do
      let (LinearChainRepresentation n group) = fromRight (error "oops") $ mkLinearChainRepresentation 3 0
      let basisStates :: S.Vector Word64
          basisStates = fromList [0 .. bit n - 1]
          isRep :: S.Vector Double
          isRep = fromList [3, 1, 0, 1, 0, 0, 0, 3]
      referenceIsRepresentative n group basisStates `shouldBe` isRep
      kernelPtr <- createIsRepresentativeKernel_v2 group Nothing
      invokeIsRepresentativeKernel kernelPtr basisStates `shouldReturn` isRep
      destroyIsRepresentativeKernel_v2 kernelPtr

  --   modifyMaxSuccess (const 10)
  --     $ prop "small linear chains"
  --     $ \(LinearChainTestData n group basisStates) -> do
  --       kernelPtr <- createIsRepresentativeKernel_v2 (unRepresentation group) Nothing
  --       invokeIsRepresentativeKernel kernelPtr basisStates
  --         `shouldReturn` referenceIsRepresentative n (unRepresentation group) basisStates

  -- describe "stateInfo" $ do
  --   it "example 1" $ do
  --     let (LinearChainRepresentation n group) = fromRight (error "oops") $ mkLinearChainRepresentation 3 0
  --     let basisStates :: S.Vector Word64
  --         basisStates = fromList [0 .. bit n - 1]
  --         (reps, indices) = referenceStateInfo n (unRepresentation group) basisStates
  --     kernelPtr <- createStateInfoKernel_v2 (unRepresentation group) Nothing
  --     invokeStateInfoKernel kernelPtr basisStates `shouldReturn` (reps, indices)
  --     destroyStateInfoKernel_v2 kernelPtr

  --   modifyMaxSuccess (const 10)
  --     $ prop "small linear chains"
  --     $ \(LinearChainTestData n group basisStates) -> do
  --       kernelPtr <- createStateInfoKernel_v2 (unRepresentation group) Nothing
  --       invokeStateInfoKernel kernelPtr basisStates `shouldReturn` referenceStateInfo n (unRepresentation group) basisStates
  --       destroyStateInfoKernel_v2 kernelPtr
  describe "fixedHammingStateToIndex" $ do
    -- it "compiles" $ do
    --   kernelPtr <- createFixedHammingStateToIndexKernel 4 2
    --   let states = S.fromList [0b11, 0b110, 0b101, 0b1100]
    --   invokeFixedHammingStateToIndexKernel kernelPtr states `shouldReturn` S.fromList [0, 2, 1, 5]
    --   destroyFixedHammingStateToIndexKernel kernelPtr
    --   True `shouldBe` True
    modifyMaxSuccess (const 10) $
      prop "computes indices correctly" $
        \(FixedHammingBitStrings numberSites hammingWeight indices states) -> do
          bracket
            (createFixedHammingStateToIndexKernel numberSites hammingWeight)
            destroyFixedHammingStateToIndexKernel
            $ \kernelPtr ->
              invokeFixedHammingStateToIndexKernel kernelPtr states `shouldReturn` indices
  describe "ls_hs_fixed_hamming_state_to_index" $ do
    prop "computes indices correctly" $
      \(FixedHammingBitStrings numberSites hammingWeight indices states) -> do
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
