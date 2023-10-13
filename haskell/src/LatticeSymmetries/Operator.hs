{-# LANGUAGE CApiFFI #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE OverloadedRecordDot #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE UndecidableInstances #-}

module LatticeSymmetries.Operator
  ( Operator (..)
  , mkOperator
  , SomeOperator (..)
  , withSomeOperator
  , foldSomeOperator
  , maxNumberOffDiag
  , operatorSymmetryGroup
  , operatorAbelianRepresentations
  , isInvariant
  , applyPermutation

    -- ** Helpers for FFI
  , newCoperator
  , cloneCoperator
  , destroyCoperator
  , withCoperator
  )
where

import Data.Primitive.Ptr qualified as P
import Data.Set qualified as Set
import Data.Vector (Vector)
import Data.Vector qualified as B
import Data.Vector.Generic qualified as G
import Foreign.Marshal
import Foreign.Ptr
import Foreign.StablePtr
import Foreign.Storable
import GHC.Show qualified
import LatticeSymmetries.Algebra
import LatticeSymmetries.Basis
import LatticeSymmetries.Benes (Permutation (unPermutation))
import LatticeSymmetries.BitString
import LatticeSymmetries.ComplexRational
import LatticeSymmetries.Expr
import LatticeSymmetries.FFI
import LatticeSymmetries.Generator
import LatticeSymmetries.Group (Hypergraph (Hypergraph), PermutationGroup (PermutationGroup), Symmetries, abelianSubgroup, groupRepresentations, hypergraphAutomorphisms, mkMultiplicationTable, mkSymmetriesFromRepresentation)
import LatticeSymmetries.NonbranchingTerm
import LatticeSymmetries.Utils
import Prelude hiding (Product, Sum)

data Operator (t :: ParticleTy) = Operator
  { opBasis :: !(Basis t)
  , opTerms :: !(Expr t)
  }

deriving stock instance
  (Show (Basis t), Show (IndexType t), Show (GeneratorType t))
  => Show (Operator t)

data SomeOperator where
  SomeOperator :: (IsBasis t) => ParticleTag t -> Operator t -> SomeOperator

instance Show SomeOperator where
  show (SomeOperator SpinTag x) = show x
  show (SomeOperator SpinlessFermionTag x) = show x
  show (SomeOperator SpinfulFermionTag x) = show x

withSomeOperator :: SomeOperator -> (forall t. (IsBasis t) => Operator t -> a) -> a
withSomeOperator x f = case x of
  SomeOperator _ operator -> f operator
{-# INLINE withSomeOperator #-}

foldSomeOperator :: (forall t. (IsBasis t) => Operator t -> a) -> SomeOperator -> a
foldSomeOperator f x = case x of
  SomeOperator _ operator -> f operator
{-# INLINE foldSomeOperator #-}

mkOperator :: (IsBasis t) => Basis t -> Expr t -> Operator t
mkOperator = Operator

getNonbranchingTerms
  :: (HasNonbranchingRepresentation (Generator Int (GeneratorType t)))
  => Operator t
  -> Vector NonbranchingTerm
getNonbranchingTerms operator =
  case nonbranchingRepresentation <$> opTermsFlat operator of
    -- NOTE: sorting based on nbtX is important!
    -- Terms with the same nbtX generate the same spin configuration. If they
    -- appear one after another, we can eliminate the duplicates easily. This
    -- is done in the off_diag kernel in Chapel.
    -- If duplicates are not eliminated, we might run into bufer overflows
    -- since we allocate buffers assuming that there are not duplicates.
    (Sum v) -> sortVectorBy (comparing (.nbtX)) $ G.filter ((/= 0) . (.nbtV)) v
  where
    opTermsFlat :: Operator t -> Polynomial ComplexRational (Generator Int (GeneratorType t))
    opTermsFlat (Operator basis terms) =
      fmap (fmap (fmap (\(Generator i g) -> Generator (flattenIndex basis i) g))) (unExpr terms)

-- typedef struct ls_hs_nonbranching_terms {
--   int number_terms;
--   int number_bits;
--   // number_words = ceil(number_bits / 64)
--   ls_hs_scalar const *v; // array of shape [number_terms]
--   uint64_t const *m;     // array of shape [number_terms, number_words]
--   uint64_t const *l;     // array of shape [number_terms, number_words]
--   uint64_t const *r;     // array of shape [number_terms, number_words]
--   uint64_t const *x;     // array of shape [number_terms, number_words]
--   uint64_t const *s;     // array of shape [number_terms, number_words]
--   // all arrays are contiguous in row-major order
-- } ls_hs_nonbranching_terms;

-- writeWords :: (Bits i, Integral i) => Int -> Ptr Word64 -> i -> IO ()
-- writeWords count p x₀ = go 0 x₀
--   where
--     nextWord !x = case bitSizeMaybe x of
--       Just n -> if n > 64 then x `shiftR` 64 else 0
--       Nothing -> x `shiftR` 64
--     go !i !x
--       | i < count = do
--         pokeElemOff p i (fromIntegral x :: Word64)
--         go (i + 1) (nextWord x)
--       | otherwise = pure ()

createCnonbranching_terms :: Int -> Vector NonbranchingTerm -> IO (Ptr Cnonbranching_terms)
createCnonbranching_terms numberBits terms
  | G.null terms = pure nullPtr
  | otherwise = do
      let numberTerms = G.length terms
          numberWords = (numberBits + 63) `div` 64
      -- NOTE: the following is not exception-safe :/
      -- TODO: fix it
      vPtr <- mallocBytes (numberTerms * sizeOf (undefined :: Cscalar))
      let mkArray = mallocBytes (numberTerms * numberWords * sizeOf (undefined :: Word64))
      mPtr <- mkArray
      lPtr <- mkArray
      rPtr <- mkArray
      xPtr <- mkArray
      sPtr <- mkArray
      loopM 0 (< numberTerms) (+ 1) $ \i -> do
        let (NonbranchingTerm v m l r x s) = (G.!) terms i
        pokeElemOff vPtr i (toComplexDouble v)
        writeBitString numberWords (P.advancePtr mPtr (i * numberWords)) m
        writeBitString numberWords (P.advancePtr lPtr (i * numberWords)) l
        writeBitString numberWords (P.advancePtr rPtr (i * numberWords)) r
        writeBitString numberWords (P.advancePtr xPtr (i * numberWords)) x
        writeBitString numberWords (P.advancePtr sPtr (i * numberWords)) s
      let c_terms = Cnonbranching_terms (fromIntegral numberTerms) (fromIntegral numberBits) vPtr mPtr lPtr rPtr xPtr sPtr
      p <- malloc
      poke p c_terms
      pure p

destroyCnonbranching_terms :: Ptr Cnonbranching_terms -> IO ()
destroyCnonbranching_terms p
  | p /= nullPtr = do
      Cnonbranching_terms _ _ vPtr mPtr lPtr rPtr xPtr sPtr <- peek p
      free vPtr
      free mPtr
      free lPtr
      free rPtr
      free xPtr
      free sPtr
      free p
  | otherwise = pure ()

newCoperator :: (IsBasis t, HasCallStack) => Maybe (Ptr Cbasis) -> Operator t -> IO (Ptr Coperator)
newCoperator maybeBasisPtr x = do
  let (diag, offDiag) = G.partition nbtIsDiagonal (getNonbranchingTerms x)
      numberBits = getNumberBits . opBasis $ x
  diag_terms <- createCnonbranching_terms numberBits diag
  off_diag_terms <- createCnonbranching_terms numberBits offDiag
  basis <- case maybeBasisPtr of
    Just p -> cloneCbasis p
    Nothing -> newCbasis (opBasis x)
  payload <- newStablePtr (SomeOperator (getParticleTag . opBasis $ x) x)
  new $
    Coperator
      { coperator_refcount = AtomicCInt 1
      , coperator_basis = basis
      , coperator_off_diag_terms = off_diag_terms
      , coperator_diag_terms = diag_terms
      , -- coperator_apply_off_diag_cxt = nullPtr,
        -- coperator_apply_diag_cxt = nullPtr,
        coperator_haskell_payload = castStablePtrToPtr payload
      }

withCoperator :: Ptr Coperator -> (SomeOperator -> IO a) -> IO a
withCoperator p f = f =<< (deRefStablePtr . castPtrToStablePtr . coperator_haskell_payload) =<< peek p

cloneCoperator :: (HasCallStack) => Ptr Coperator -> IO (Ptr Coperator)
cloneCoperator p = do
  _ <- operatorIncRefCount p
  pure p

destroyCoperator :: (HasCallStack) => Ptr Coperator -> IO ()
destroyCoperator p = do
  refcount <- operatorDecRefCount p
  when (refcount == 1) $ do
    x <- peek p
    destroyCbasis (coperator_basis x)
    destroyCnonbranching_terms (coperator_off_diag_terms x)
    destroyCnonbranching_terms (coperator_diag_terms x)
    freeStablePtr . castPtrToStablePtr $ coperator_haskell_payload x
    free p

maxNumberOffDiag :: (IsBasis t) => Operator t -> Int
maxNumberOffDiag op = Set.size $ Set.fromList . G.toList . G.map nbtX $ offDiag
  where
    (_, offDiag) = G.partition nbtIsDiagonal (getNonbranchingTerms op)

operatorToHypergraph :: (IsBasis t) => Operator t -> Hypergraph Int
operatorToHypergraph (Operator basis expr) = flattenHypergraph $ exprToHypergraph expr
  where
    flattenHypergraph (Hypergraph vertices hyperedges) =
      Hypergraph
        (Set.map (flattenIndex basis) vertices)
        (Set.map (Set.map (flattenIndex basis)) hyperedges)

applyPermutation :: (IsBasis t) => Permutation -> Operator t -> Operator t
applyPermutation (unPermutation -> p) (Operator basis expr) =
  Operator basis . simplifyExpr $ mapIndices remap expr
  where
    remap i = unFlattenIndex basis $ p G.! flattenIndex basis i

isInvariant :: (IsBasis t) => Operator t -> Permutation -> Bool
isInvariant operator permutation =
  (applyPermutation permutation operator).opTerms == operator.opTerms

operatorSymmetryGroup :: (IsBasis t) => Operator t -> PermutationGroup
operatorSymmetryGroup operator = PermutationGroup $ B.filter (isInvariant operator) automorphisms
  where
    (PermutationGroup automorphisms) = hypergraphAutomorphisms . operatorToHypergraph $ operator

operatorAbelianRepresentations :: (IsBasis t) => Operator t -> [Symmetries]
operatorAbelianRepresentations operator = mkSymmetriesFromRepresentation g' <$> reps
  where
    g = operatorSymmetryGroup operator
    (g', t') = abelianSubgroup g (mkMultiplicationTable g)
    reps = groupRepresentations t'
