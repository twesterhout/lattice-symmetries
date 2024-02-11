-- |
-- Module      : LatticeSymmetries.NonbranchingTerm
-- Description : Non-branching terms in the Hamiltonian
-- Copyright   : (c) Tom Westerhout, 2022
-- Stability   : experimental
module LatticeSymmetries.NonbranchingTerm
  ( NonbranchingTerm (..)
  , HasNonbranchingRepresentation (..)
  , nbtIsDiagonal
  , getMaxNumberOffDiag
  , unsafeEstimateMaxNumberOffDiag
  )
where

import Control.Monad.ST (runST)
import Data.Bits
import Data.Set qualified as Set
import Data.Vector (Vector)
import Data.Vector.Algorithms.Intro qualified
import Data.Vector.Generic qualified as G
import Data.Vector.Generic.Mutable qualified as GM
import LatticeSymmetries.BitString
import LatticeSymmetries.ComplexRational

-- | Non-branching term of the Hamiltonian.
--
-- Suppose that \(T\) is such a term. The "non-branching" part means that when we apply \(T\) to a
-- basis state \(|\alpha\rangle\), we get just one other basis state back, i.e.
-- \( T|\alpha\rangle = c|\beta\rangle \).
--
-- We use the representation from <http://arxiv.org/abs/2203.04158>. See Appendix A in the paper.
data NonbranchingTerm = NonbranchingTerm
  -- (v, m, l, r, x, s)
  { nbtV :: !ComplexRational
  , nbtM :: !BitString
  , nbtL :: !BitString
  , nbtR :: !BitString
  , nbtX :: !BitString
  , nbtS :: !BitString
  }
  deriving stock (Show, Eq)

-- | Implements a .&. complement b without using complement function
andWithComplement :: BitString -> BitString -> BitString
andWithComplement a b = a .&. (b `xor` fakeOnes)
  where
    fakeOnes = go (BitString 0xFFFFFFFFFFFFFFFF) 64
    go !n !i = if n >= a then n else go ((n `shiftL` i) .|. n) (2 * i)

-- | Composition of operators.
instance Semigroup NonbranchingTerm where
  (<>) (NonbranchingTerm vₐ mₐ lₐ rₐ _ sₐ) (NonbranchingTerm vᵦ mᵦ lᵦ rᵦ xᵦ sᵦ) = NonbranchingTerm v m l r x s
    where
      v = if (rₐ `xor` lᵦ) .&. mₐ .&. mᵦ /= zeroBits then 0 else ((-1) ^ p) * vₐ * vᵦ
      m = mₐ .|. mᵦ
      -- What we want is
      --   r = rᵦ .|. (rₐ .&. complement mᵦ)
      --   l = lₐ .|. (lᵦ .&. complement mₐ)
      -- however, complement does not work for Integer, because that would create an infinitely long
      -- bitstring
      r = rᵦ .|. andWithComplement rₐ mᵦ
      l = lₐ .|. andWithComplement lᵦ mₐ
      x = l `xor` r
      s = andWithComplement (sₐ `xor` sᵦ) m
      z = (r .&. sᵦ) `xor` xᵦ
      p = popCount $ (r .&. sᵦ) `xor` (z .&. sₐ)

-- | Specifies that an operator is non-branching and can be cast into 'NonbranchingTerm'.
--
-- Our low-level C kernels deal with non-branching terms for efficiency, but in the Haskell land we
-- want to deal with mathematical expressions (i.e. polynomials of creation/annihilation operators).
-- This typeclass connects the two worlds.
--
-- @g@ will typically be 'LatticeSymmetries.Generator.SpinGeneratorType', 'LatticeSymmetries.Generator.FermionGeneratorType' or 'LatticeSymmetries.Generator.Generator'.
class HasNonbranchingRepresentation g where
  nonbranchingRepresentation :: g -> NonbranchingTerm

-- | Checks whether an operator is diagonal.
nbtIsDiagonal :: NonbranchingTerm -> Bool
nbtIsDiagonal t = nbtX t == zeroBits

step :: HasCallStack => (NonbranchingTerm -> BitString) -> Int -> Vector NonbranchingTerm -> Vector (Vector NonbranchingTerm)
step !p !i !v = G.filter (not . G.null) $ fromListN 2 [ups <> nones, downs <> nones]
  where
    hasConstraint t = testBit (p t) i || testBit t.nbtM i
    unsafeGetConstraint t = testBit (p t) i
    (ups, downs, nones) = runST $ do
      mv <- G.unsafeThaw v
      k1 <- GM.unstablePartition hasConstraint mv
      if k1 == 0
        then pure (G.empty, G.empty, v)
        else do
          k2 <- GM.unstablePartition unsafeGetConstraint (GM.take k1 mv)
          (,,) <$> G.unsafeFreeze (GM.take k2 mv)
              <*> G.unsafeFreeze (GM.slice k2 (k1 - k2) mv)
              <*> G.unsafeFreeze (GM.drop k1 mv)

getNumberUnique :: (G.Vector v a, Ord a) => v a -> Int
getNumberUnique = Set.size . Set.fromList . G.toList

unsafeEstimateMaxNumberOffDiag' :: (NonbranchingTerm -> BitString) -> Int -> Vector NonbranchingTerm -> Int
unsafeEstimateMaxNumberOffDiag' p numberBits v0 = if G.null v0 then 0 else go 0 (G.singleton v0)
  where
    go :: Int -> Vector (Vector NonbranchingTerm) -> Int
    go !i !vs
      | i < numberBits = go (i + 1) $ runST $ do
          let expanded = G.map (\x -> (getNumberUnique $ (.nbtX) <$> x, x)) $ G.concatMap (step p i) vs
          mv <- G.unsafeThaw expanded -- G.concatMap (step p i) vs
          let k = min (GM.length mv) 16
          Data.Vector.Algorithms.Intro.selectBy (comparing (negate . fst)) mv k
          G.map snd . G.take k <$> G.unsafeFreeze mv
      | otherwise = G.maximum $ getNumberUnique . G.map (.nbtX) <$> vs

unsafeEstimateMaxNumberOffDiag :: Int -> Maybe Int -> Vector NonbranchingTerm -> Int
unsafeEstimateMaxNumberOffDiag numberBits _hammingWeigh terms = max l r
  where
    !l = unsafeEstimateMaxNumberOffDiag' (.nbtL) numberBits terms
    !r = unsafeEstimateMaxNumberOffDiag' (.nbtR) numberBits terms

getMaxNumberOffDiag :: Int -> Maybe Int -> Vector NonbranchingTerm -> Int
getMaxNumberOffDiag _numberBits _hammingWeight terms = getNumberUnique $ (.nbtX) <$> terms
