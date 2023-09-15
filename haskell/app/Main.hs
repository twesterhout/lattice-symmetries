import Data.Vector.Generic qualified as G
import LatticeSymmetries.Group

main :: IO ()
main = do
  let !g =
        groupFromTransversalGeneratingSet . transversalGeneratingSet $
          autsSearchTree (rectangularGraph 16 16)
  let !t = mkMultiplicationTable g
  -- print (G.length g)
  print $ G.length <$> take 16 (abelianSubset t)
