{-# LANGUAGE OverloadedRecordDot #-}

import Control.Exception
import Data.Maybe (fromJust)
import Data.Vector.Storable qualified as S
import Data.Vector.Storable.Mutable qualified as SM
import Foreign (castPtr)
import Gauge.Benchmark
import Gauge.Main
import Language.Halide hiding ((==))
import LatticeSymmetries.Basis
import LatticeSymmetries.Dense (DenseMatrix (DenseMatrix))
import LatticeSymmetries.Lowering
import System.Random

main :: IO ()
main = do
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


