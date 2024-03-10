{-# LANGUAGE OverloadedRecordDot #-}

import Control.Exception
import Data.Maybe (fromJust)
import Data.Vector.Storable qualified as S
import Data.Vector.Storable.Mutable qualified as SM
import Foreign (Bits (bit), castPtr)
import Gauge.Benchmark
import Gauge.Main
import Language.Halide hiding ((==))
import LatticeSymmetries.Group (fromGenerators, mkSymmetry)
import LatticeSymmetries.Lowering
import LatticeSymmetries.Permutation (mkPermutation)
import System.Random
import Prelude hiding (group)

benchmarkFixedHammingStateToIndex :: IO ()
benchmarkFixedHammingStateToIndex = do
  let pureGen = mkStdGen 140
      numberSites = 36
      hammingWeight = 18

      indices :: S.Vector Int64
      !indices =
        S.fromList
          . take 1000
          . unfoldr (Just . uniformR (0, fromIntegral (fromJust (binomial numberSites hammingWeight)) - 1))
          $ pureGen
      !states = S.map (fixedHammingIndexToState hammingWeight . fromIntegral) indices

  out <- SM.new @_ @Int64 (S.length states)
  bracket (createFixedHammingStateToIndexKernel numberSites hammingWeight) destroyFixedHammingStateToIndexKernel $
    \kernelPtr -> do
      withHalideBuffer @1 @Word64 states $ \basisStatesBuf ->
        withHalideBuffer @1 @Int64 out $ \indicesBuf -> do
          let !b = nfIO $ toFun_fixed_hamming_state_to_index_kernel kernelPtr (castPtr basisStatesBuf) (castPtr indicesBuf)
          benchmark b
      print . (indices ==) =<< S.freeze out

  S.unsafeWith states $ \basisStatesBuf ->
    SM.unsafeWith out $ \indicesBuf -> do
      let !b = nfIO $ ls_hs_fixed_hamming_state_to_index (fromIntegral (S.length states)) basisStatesBuf indicesBuf
      benchmark b
  print . (indices ==) =<< S.freeze out

benchmarkStateInfo :: IO ()
benchmarkStateInfo = do
  let pureGen = mkStdGen 140
      numberSites = 40

      basisStates :: S.Vector Word64
      !basisStates =
        S.fromList
          . take 1024
          . unfoldr (Just . uniformR (0, bit numberSites - 1))
          $ pureGen

      Right p = mkPermutation . fromList $ ([1 .. numberSites - 1] <> [0])
      Right t = mkSymmetry p 0
      Right group = fromGenerators [t]

  representatives <- SM.new @_ @Word64 (S.length basisStates)
  indices <- SM.new @_ @Int32 (S.length basisStates)

  -- bracket (createStateInfoKernel_v2 group Nothing) destroyStateInfoKernel_v2 $ \kernelPtr ->
  --   withHalideBuffer @1 @Word64 basisStates $ \basisStatesBuf ->
  --     withHalideBuffer @1 @Word64 representatives $ \representativesBuf ->
  --       withHalideBuffer @1 @Int32 indices $ \indicesBuf -> do
  --         let !b = nfIO $ toFun_state_info_kernel kernelPtr (castPtr basisStatesBuf) (castPtr representativesBuf) (castPtr indicesBuf)
  --         benchmark b

  -- bracket (createStateInfoKernel_v3 group Nothing) destroyStateInfoKernel_v3 $ \kernelPtr ->
  --   withHalideBuffer @1 @Word64 basisStates $ \basisStatesBuf ->
  --     withHalideBuffer @1 @Word64 representatives $ \representativesBuf ->
  --       withHalideBuffer @1 @Int32 indices $ \indicesBuf -> do
  --         let !b = nfIO $ toFun_state_info_kernel kernelPtr (castPtr basisStatesBuf) (castPtr representativesBuf) (castPtr indicesBuf)
  --         benchmark b

  norms <- SM.new @_ @Word16 (S.length basisStates)
  bracket (createIsRepresentativeKernel_v2 group Nothing) destroyIsRepresentativeKernel_v2 $ \kernelPtr ->
    withHalideBuffer @1 @Word64 basisStates $ \basisStatesBuf ->
      withHalideBuffer @1 @Word16 norms $ \normsBuf -> do
        let !b = nfIO $ toFun_is_representative_kernel kernelPtr (castPtr basisStatesBuf) (castPtr normsBuf)
        benchmark b

main :: IO ()
main = do
  benchmarkStateInfo
