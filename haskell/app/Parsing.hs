{-# LANGUAGE OverloadedRecordDot #-}

import LatticeSymmetries.Basis
import LatticeSymmetries.Dense (DenseMatrix (DenseMatrix))
import LatticeSymmetries.Lowering
import LatticeSymmetries.Some
import qualified Data.Aeson as Aeson
import LatticeSymmetries.Expr (estimateNumberSites)
import LatticeSymmetries.Operator (SomeOperator(SomeOperator), Operator (Operator), foldSomeOperator)

main :: IO ()
main = do
  r <-
    Aeson.eitherDecodeFileStrict @SomeOperator "/home/tom/Projects/lattice-symmetries/chapel/data/nikita/expression_25.json"
  case r of
    Right expr -> print (foldSomeOperator (\(Operator _ e) -> estimateNumberSites e) expr)
    Left e -> putStrLn e
