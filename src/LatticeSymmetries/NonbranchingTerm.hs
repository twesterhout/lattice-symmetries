module LatticeSymmetries.NonbranchingTerm
  ( NonbranchingTerm (..),
    HasNonbranchingRepresentation (..),
    isNonbranchingTermDiagonal,
    applyNonbranchingTerm,
  )
where

import Data.Bits
import LatticeSymmetries.BitString
import LatticeSymmetries.ComplexRational

data NonbranchingTerm = NonbranchingTerm
  -- (v, m, l, r, x, s)
  { nbtV :: !ComplexRational,
    nbtM :: !BitString,
    nbtL :: !BitString,
    nbtR :: !BitString,
    nbtX :: !BitString,
    nbtS :: !BitString
  }
  deriving stock (Show, Eq)

-- | implements a .&. complement b without using complement function
andWithComplement :: BitString -> BitString -> BitString
andWithComplement a b = a .&. (b `xor` fakeOnes)
  where
    fakeOnes = go (BitString 0xFFFFFFFFFFFFFFFF) 64
    go !n !i = if n >= a then n else go ((n `shiftL` i) .|. n) (2 * i)

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

class HasNonbranchingRepresentation g where
  nonbranchingRepresentation :: g -> NonbranchingTerm

isNonbranchingTermDiagonal :: NonbranchingTerm -> Bool
isNonbranchingTermDiagonal t = nbtX t == zeroBits

applyNonbranchingTerm :: NonbranchingTerm -> BitString -> (ComplexRational, BitString)
applyNonbranchingTerm (NonbranchingTerm v m _ r x s) α = (c, β)
  where
    δ = fromEnum $ (α .&. m) == r
    sign = 1 - 2 * (popCount (α .&. s) `mod` 2)
    β = α `xor` x
    c = v * fromIntegral (δ * sign)
