{-# LANGUAGE OverloadedRecordDot #-}

import Data.Aeson qualified as Aeson
import Data.Vector.Generic qualified as G
import LatticeSymmetries.Basis
import LatticeSymmetries.Dense (DenseMatrix (DenseMatrix))
import LatticeSymmetries.Expr
import LatticeSymmetries.Lowering
import LatticeSymmetries.NonbranchingTerm
import LatticeSymmetries.Operator
import LatticeSymmetries.Some

main :: IO ()
main = do
  r <-
    Aeson.eitherDecodeFileStrict @SomeOperator "/home/tom/Projects/lattice-symmetries/chapel/data/nikita/expression_25.json"
      -- "/home/tom/Projects/lattice-symmetries/test/random_quspin/test_10_6_10_40_expr.json"

  case r of
    Right op -> withSomeOperator op $ \operator@(Operator basis _) -> do
      let offDiag = G.filter (not . nbtIsDiagonal) (getNonbranchingTerms operator)
          numberBits = getNumberBits basis
      print $ unsafeEstimateMaxNumberOffDiag numberBits (getHammingWeight basis) offDiag
      print $ getMaxNumberOffDiag numberBits (getHammingWeight basis) offDiag
    Left e -> putStrLn e
