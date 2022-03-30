module LatticeSymmetries.Group
  ( Permutation,
    unPermutation,
    mkPermutation,
    identityPermutation,
    fromGenerators,
    pgLength,
  )
where

import Control.Monad.ST
import qualified Data.Set as S
import qualified Data.Vector as B
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Generic.Mutable as GM
import qualified Data.Vector.Unboxed as U
import LatticeSymmetries.Benes

newtype PermutationGroup = PermutationGroup (B.Vector Permutation)
  deriving stock (Show, Eq)

pgLength :: PermutationGroup -> Int
pgLength (PermutationGroup gs) = G.length gs

fromGenerators :: [Permutation] -> PermutationGroup
fromGenerators [] = PermutationGroup G.empty
fromGenerators gs@(p : _) = go S.empty (S.singleton (identityPermutation (G.length v)))
  where
    v = unPermutation p
    go !interior !boundary
      | S.null boundary = PermutationGroup . G.fromList . S.toAscList $ interior
      | otherwise = go interior' boundary'
      where
        interior' = interior `S.union` boundary
        boundary' = S.fromList [h <> g | h <- S.toList boundary, g <- gs] S.\\ interior'
