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
import Data.Vector.Unboxed qualified as U
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

-- countingMerge :: G.Vector v a => (a -> a -> Ordering) -> v a -> v a -> (Int, v a)
-- countingMerge cmp va vb = runST $ do
--   let !na = G.length va
--       !nb = G.length vb
--   dest <- GM.unsafeNew (na + nb)
--   let write !count !i !x
--         | i > 0 = do
--             y <- GM.unsafeRead dest (i - 1)
--             GM.unsafeWrite dest i x
--             case cmp y x of
--               EQ -> pure count
--               _ -> pure (count + 1)
--         | otherwise = do
--             GM.unsafeWrite dest i x
--             pure (count + 1)
--   let go !ia !ib !count
--         | ia < na && ib < nb = do
--             let !a = G.unsafeIndex va ia
--                 !b = G.unsafeIndex vb ib
--             case cmp a b of
--               LT -> go (ia + 1) ib =<< write count (ia + ib) a
--               GT -> go ia (ib + 1) =<< write count (ia + ib) b
--               EQ -> do
--                 count' <- write count (ia + ib) a
--                 GM.unsafeWrite dest (ia + ib + 1) b
--                 go (ia + 1) (ib + 1) count'
--         | ia < na = go (ia + 1) ib =<< write count (ia + ib) (G.unsafeIndex va ia)
--         | ib < nb = go ia (ib + 1) =<< write count (ia + ib) (G.unsafeIndex vb ib)
--         | otherwise = pure count
--   count <- go 0 0 0
--   (count,) <$> G.unsafeFreeze dest
-- {-# SCC countingMerge #-}

mergeX :: U.Vector (Word64, Word64, Word64) -> U.Vector (Word64, Word64, Word64) -> U.Vector (Word64, Word64, Word64)
mergeX !va !vb = runST $ do
  let !na = G.length va
      !nb = G.length vb
  dest <- GM.unsafeNew (na + nb)
  let remainderA !ia !ib
        | ia < na = G.unsafeIndexM va ia >>= GM.unsafeWrite dest (ia + ib) >> remainderA (ia + 1) ib
        | otherwise = pure ()
  let remainderB !ia !ib
        | ib < nb = G.unsafeIndexM vb ib >>= GM.unsafeWrite dest (ia + ib) >> remainderB ia (ib + 1)
        | otherwise = pure ()
  let go !ia !ib
        | ia < na && ib < nb = do
            let a@(!_, !_, !xa) = G.unsafeIndex va ia
                b@(!_, !_, !xb) = G.unsafeIndex vb ib
            case compare xa xb of
              LT -> GM.unsafeWrite dest (ia + ib) a >> go (ia + 1) ib
              GT -> GM.unsafeWrite dest (ia + ib) b >> go ia (ib + 1)
              EQ -> GM.unsafeWrite dest (ia + ib) a >> GM.unsafeWrite dest (ia + ib + 1) b >> go (ia + 1) (ib + 1)
        | ia < na = remainderA ia ib
        | ib < nb = remainderB ia ib
        | otherwise = pure ()
  go 0 0
  G.unsafeFreeze dest
{-# SCC mergeX #-}

getCountX :: U.Vector (Word64, Word64, Word64) -> Int
getCountX !v
  | G.null v = 0
  | otherwise = let !r = go 1 (G.unsafeIndex v 0) 1 in r
  where
    go !n (!_, !_, !xa) !i
      | i < G.length v =
          let b@(!_, !_, !xb) = G.unsafeIndex v i
           in if xa == xb then go n b (i + 1) else go (n + 1) b (i + 1)
      | otherwise = n
{-# SCC getCountX #-}

stepU :: Int -> (Int, U.Vector (Word64, Word64, Word64)) -> Vector (Int, U.Vector (Word64, Word64, Word64))
stepU !i (!n, !v)
  | n <= 1 = G.singleton (n, v)
  | otherwise = runST $ do
      let (!upsAndDowns, !nones) = U.partition hasConstraint v
      mv <- G.unsafeThaw upsAndDowns
      k2 <- GM.unstablePartition unsafeGetConstraint mv
      let (ups, downs) = GM.splitAt k2 mv
      Data.Vector.Algorithms.Intro.sortBy comparingX ups
      Data.Vector.Algorithms.Intro.sortBy comparingX downs
      v1 <- mergeX <$> G.unsafeFreeze ups <*> pure nones
      v2 <- mergeX <$> G.unsafeFreeze downs <*> pure nones
      pure $ fromListN 2 [(getCountX v1, v1), (getCountX v2, v2)]
  where
    hasConstraint (!r, !m, !_) = testBit r i || testBit m i
    unsafeGetConstraint (!r, !_, !_) = testBit r i
    comparingX (!_, !_, !x1) (!_, !_, !x2) = compare x1 x2
{-# SCC stepU #-}

getNumberUnique :: (G.Vector v a, Ord a) => v a -> Int
getNumberUnique = Set.size . Set.fromList . G.toList

unsafeEstimateMaxNumberOffDiag' :: Int -> U.Vector (Word64, Word64, Word64) -> Int
unsafeEstimateMaxNumberOffDiag' numberBits v0
  | G.null v0 = 0
  | otherwise = go 0 $ G.singleton (getNumberUnique (G.map getX v0), v0)
  where
    getX (!_, !_, x) = x
    go :: Int -> Vector (Int, U.Vector (Word64, Word64, Word64)) -> Int
    go !i !vs
      | i < numberBits = go (i + 1) $ runST $ do
          let expanded = G.concatMap (stepU i) vs
          mv <- G.unsafeThaw expanded
          let k = min (GM.length mv) 64
          Data.Vector.Algorithms.Intro.selectBy (\(n1, _) (n2, _) -> compare n2 n1) mv k
          G.take k <$> G.unsafeFreeze mv
      | otherwise = G.maximum $ G.map fst vs

unsafeEstimateMaxNumberOffDiag :: Int -> Maybe Int -> Vector NonbranchingTerm -> Int
unsafeEstimateMaxNumberOffDiag numberBits hammingWeight terms
  | numberBits > 64 = getMaxNumberOffDiag numberBits hammingWeight terms
  | otherwise = max l r
  where
    toW :: BitString -> Word64
    toW = fromIntegral . unBitString
    !l =
      unsafeEstimateMaxNumberOffDiag' numberBits . G.convert $
        G.map (\t -> (toW t.nbtL, toW t.nbtM, toW t.nbtX)) terms
    !r =
      unsafeEstimateMaxNumberOffDiag' numberBits . G.convert $
        G.map (\t -> (toW t.nbtR, toW t.nbtM, toW t.nbtX)) terms

getMaxNumberOffDiag :: Int -> Maybe Int -> Vector NonbranchingTerm -> Int
getMaxNumberOffDiag _numberBits _hammingWeight terms = getNumberUnique $ (.nbtX) <$> terms
