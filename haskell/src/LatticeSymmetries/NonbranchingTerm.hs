-- |
-- Module      : LatticeSymmetries.NonbranchingTerm
-- Description : Non-branching terms in the Hamiltonian
-- Copyright   : (c) Tom Westerhout, 2022
-- Stability   : experimental
module LatticeSymmetries.NonbranchingTerm
  ( NonbranchingTerm (..)
  , HasNonbranchingRepresentation (..)
  , nbtIsDiagonal
  )
where

import Data.Bits (Bits (complement, popCount, shiftL, zeroBits, (.&.), (.|.)))
import LatticeSymmetries.BitString (BitString (BitString))
import LatticeSymmetries.ComplexRational (ComplexRational)

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
      s = (sₐ `xor` sᵦ) .&. complement m
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
