module LatticeSymmetries.Dense
  ( DenseMatrix (..),
    dmShape,
    denseMatrixFromList,
    denseMatrixToList,
    indexDenseMatrix,
    extractDiagonal,
    countOffDiagNonZero,
  )
where

import Data.Vector.Generic ((!))
import qualified Data.Vector.Generic as G
import qualified GHC.Exts as GHC (IsList (..))

-- | Dense matrix in row-major order (C layout)
data DenseMatrix v a = DenseMatrix {dmRows :: !Int, dmCols :: !Int, dmData :: !(v a)}
  deriving stock (Show, Eq, Generic)
  deriving anyclass (NFData)

-- | Get matrix shape
dmShape :: DenseMatrix v a -> (Int, Int)
dmShape m = (dmRows m, dmCols m)

instance (G.Vector v a, Num a) => Num (DenseMatrix v a) where
  (+) a b = DenseMatrix (dmRows a) (dmCols a) $ G.zipWith (+) (dmData a) (dmData b)
  (-) a b = DenseMatrix (dmRows a) (dmCols a) $ G.zipWith (-) (dmData a) (dmData b)
  (*) a b = DenseMatrix (dmRows a) (dmCols a) $ G.zipWith (*) (dmData a) (dmData b)
  abs a = DenseMatrix (dmRows a) (dmCols a) $ G.map abs (dmData a)
  signum a = DenseMatrix (dmRows a) (dmCols a) $ G.map signum (dmData a)
  fromInteger _ = error "Num instance of DenseMatrix does not implement fromInteger"

denseMatrixFromList :: (HasCallStack, G.Vector v a) => [[a]] -> DenseMatrix v a
denseMatrixFromList rs
  | G.length elements == nRows * nCols = DenseMatrix nRows nCols elements
  | otherwise = error "nested list has irregular shape"
  where
    !nRows = length rs
    !nCols = case rs of
      [] -> 0
      (r : _) -> length r
    elements = G.fromListN (nRows * nCols) $ mconcat rs

denseMatrixToList :: G.Vector v a => DenseMatrix v a -> [[a]]
denseMatrixToList (DenseMatrix _ nCols v) = go (G.toList v)
  where
    go elements = case splitAt nCols elements of
      (row, []) -> [row]
      (row, rest) -> row : go rest

instance (G.Vector v a) => GHC.IsList (DenseMatrix v a) where
  type Item (DenseMatrix v a) = [a]
  fromList = denseMatrixFromList
  toList = denseMatrixToList

countOffDiagNonZero :: (G.Vector v a, Eq a, Num a) => DenseMatrix v a -> Int
countOffDiagNonZero m@(DenseMatrix nRows nCols _) = go1 0 0
  where
    go2 !i !j !acc
      | i == j = go2 i (j + 1) acc
      | j < nCols =
        if indexDenseMatrix m (i, j) /= 0
          then go2 i (j + 1) (acc + 1)
          else go2 i (j + 1) acc
      | otherwise = acc
    go1 !i !acc
      | i < nRows = let acc' = go2 i 0 acc in go1 (i + 1) acc'
      | otherwise = acc

indexDenseMatrix :: (HasCallStack, G.Vector v a) => DenseMatrix v a -> (Int, Int) -> a
indexDenseMatrix (DenseMatrix nRows nCols v) (i, j)
  | i < nRows && j < nCols = v ! (nCols * i + j)
  | otherwise = error $ "invalid index " <> show (i, j) <> " for a matrix of shape " <> show (nRows, nCols)

extractDiagonal :: (HasCallStack, G.Vector v a, G.Vector v Int) => DenseMatrix v a -> v a
extractDiagonal m@(DenseMatrix nRows nCols _)
  | nRows == nCols = G.map (\i -> indexDenseMatrix m (i, i)) (G.enumFromN 0 nRows)
  | otherwise = error "cannot extract the diagonal of a non-square matrix"
