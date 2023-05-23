module Main where

import LatticeSymmetries.Basis
import LatticeSymmetries.Operator
import LatticeSymmetries.Parser
import Text.PrettyPrint.ANSI.Leijen (hardline, pretty, putDoc)

main :: IO ()
main = do
  let basis = SpinBasis 4 (Just 2)
      indices = [[i, (i + 1) `div` 4] | i <- [0 .. 3]]
      a = operatorFromString basis "σᶻ₀ σᶻ₁" indices
      b = operatorFromString basis "σˣ₀ σˣ₁" indices
      c = operatorFromString basis "σʸ₀ σʸ₁" indices
  putStrLn "Starting ..."
  putDoc $ pretty (opTerms (a + b + c)) <> hardline
  putStrLn "Done!"
