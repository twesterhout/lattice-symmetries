{-# LANGUAGE CApiFFI #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE OverloadedRecordDot #-}
{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE UndecidableInstances #-}

module LatticeSymmetries.Operator
  ( Operator (..)
  , mkOperator
  , SomeOperator (..)
  , withSomeOperator
  , foldSomeOperator
  , getNonbranchingTerms
  , getHilbertSpaceSectors
  -- , maxNumberOffDiag
  -- , operatorSymmetryGroup
  -- , operatorAbelianRepresentations
  -- , isInvariant
  -- , applyPermutation

    -- ** Helpers for FFI

  -- , newCoperator
  -- , cloneCoperator
  -- , destroyCoperator
  , withCoperator
  , foldCoperator
  , newCoperator
  , ls_hs_create_operator
  , ls_hs_destroy_operator
  , ls_hs_operator_max_number_off_diag
  , ls_hs_operator_to_json
  , ls_hs_operator_from_json
  , ls_hs_expr_hilbert_space_sectors
  )
where

import Data.Aeson (FromJSON (parseJSON), ToJSON (toJSON), object, withObject, (.:), (.=))
import Data.Primitive.Ptr qualified as P
import Data.Vector (Vector)
import Data.Vector.Generic qualified as G
import Foreign.C (CBool (..), CInt (..), CSize (..), CString)
import Foreign.Marshal
import Foreign.Ptr
import Foreign.StablePtr
import Foreign.Storable
import GHC.Show qualified
import Language.C.Inline.Unsafe qualified as CU
import LatticeSymmetries.Algebra
import LatticeSymmetries.Basis
import LatticeSymmetries.BitString
import LatticeSymmetries.ComplexRational
import LatticeSymmetries.Context
import LatticeSymmetries.Expr
import LatticeSymmetries.FFI
import LatticeSymmetries.Generator
import LatticeSymmetries.Group (Representation (Representation), abelianRepresentations, abelianSubgroup, emptyRepresentation)
import LatticeSymmetries.NonbranchingTerm
import LatticeSymmetries.Permutation (Permutation (Permutation))
import LatticeSymmetries.Some
import LatticeSymmetries.Utils
import Prelude hiding (Product, Sum, group)

importLS

data Operator (t :: ParticleTy) = Operator
  { opBasis :: !(Basis t)
  , opTerms :: !(Expr t)
  }

deriving stock instance
  (Show (Basis t), Show (IndexType t), Show (GeneratorType t))
  => Show (Operator t)

deriving stock instance
  (Eq (Basis t), Eq (IndexType t), Eq (GeneratorType t))
  => Eq (Operator t)

data SomeOperator where
  SomeOperator :: IsBasis t => ParticleTag t -> Operator t -> SomeOperator

instance Show SomeOperator where
  show (SomeOperator SpinTag x) = show x
  show (SomeOperator SpinlessFermionTag x) = show x
  show (SomeOperator SpinfulFermionTag x) = show x

instance Eq SomeOperator where
  (SomeOperator SpinTag x) == (SomeOperator SpinTag y) = x == y
  (SomeOperator SpinfulFermionTag x) == (SomeOperator SpinfulFermionTag y) = x == y
  (SomeOperator SpinlessFermionTag x) == (SomeOperator SpinlessFermionTag y) = x == y
  _ == _ = False

withSomeOperator :: SomeOperator -> (forall t. IsBasis t => Operator t -> a) -> a
withSomeOperator x f = case x of
  SomeOperator _ operator -> f operator
{-# INLINE withSomeOperator #-}

foldSomeOperator :: (forall t. IsBasis t => Operator t -> a) -> SomeOperator -> a
foldSomeOperator f x = case x of
  SomeOperator _ operator -> f operator
{-# INLINE foldSomeOperator #-}

mkOperator :: (MonadFail m, IsBasis t) => Basis t -> Expr t -> m (Operator t)
mkOperator b e = pure $ Operator b e

mkSomeOperator :: MonadFail m => SomeBasis -> SomeExpr -> m SomeOperator
mkSomeOperator (SomeBasis b@(SpinBasis {})) (SomeExpr SpinTag e) = SomeOperator SpinTag <$> mkOperator b e
mkSomeOperator (SomeBasis b@(SpinBasis {})) (SomeExpr t _) =
  fail $
    "could not construct operator from " <> show (getParticleTag b) <> " basis and " <> show t <> " expression"
mkSomeOperator (SomeBasis b) (SomeExpr t _) = fail $ "not implemented: " <> show (getParticleTag b) <> ", " <> show t

instance FromJSON SomeOperator where
  parseJSON = withObject "Operator" $ \v -> do
    basis <- v .: "basis"
    expr <- v .: "expression"
    mkSomeOperator basis expr

instance ToJSON SomeOperator where
  toJSON (SomeOperator t (Operator basis expr)) =
    object ["basis" .= SomeBasis basis, "expression" .= SomeExpr t expr]

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

-- newCbasis :: IsBasis t => Basis t -> IO (Ptr Cbasis)
-- newCbasis x = do
--   payload <- castStablePtrToPtr <$> newStablePtr (SomeBasis x)
--   -- NOTE: important to initialize memory to 0 such that we don't have to manually initialize fields
--   p <- callocBytes $ fromIntegral [CU.pure| size_t { sizeof(ls_hs_basis) } |]
--   [CU.block| void { ls_hs_internal_object_init(&$(ls_hs_basis* p)->base, 1, $(void* payload)); } |]
--   pure p
newCnonbranching_terms :: Int -> Vector NonbranchingTerm -> IO (Ptr Cnonbranching_terms)
newCnonbranching_terms numberBits terms = do
  let numberWords = (numberBits + 63) `div` 64
      numberTerms = G.length terms
      c_numberTerms = fromIntegral numberTerms
      c_numberBits = fromIntegral numberBits
      mkArray
        | numberTerms > 0 = callocBytes (numberTerms * numberWords * sizeOf (undefined :: Word64))
        | otherwise = pure nullPtr
  p <- callocBytes $ fromIntegral [CU.pure| size_t { sizeof(ls_hs_nonbranching_terms) } |]
  vPtr <-
    if numberTerms > 0
      then callocBytes (numberTerms * sizeOf (undefined :: Cscalar))
      else pure nullPtr
  mPtr <- mkArray
  lPtr <- mkArray
  rPtr <- mkArray
  xPtr <- mkArray
  sPtr <- mkArray
  [CU.block| void {
    ls_hs_nonbranching_terms* p = $(ls_hs_nonbranching_terms* p);
    p->number_terms = $(int c_numberTerms);
    p->number_bits = $(int c_numberBits);
    p->v = $(ls_hs_scalar const* vPtr);
    p->m = $(uint64_t const* mPtr);
    p->l = $(uint64_t const* lPtr);
    p->r = $(uint64_t const* rPtr);
    p->x = $(uint64_t const* xPtr);
    p->s = $(uint64_t const* sPtr);
  } |]
  loopM 0 (< numberTerms) (+ 1) $ \i -> do
    let (NonbranchingTerm v m l r x s) = terms G.! i
    pokeElemOff vPtr i (toComplexDouble v)
    writeBitString numberWords (P.advancePtr mPtr (i * numberWords)) m
    writeBitString numberWords (P.advancePtr lPtr (i * numberWords)) l
    writeBitString numberWords (P.advancePtr rPtr (i * numberWords)) r
    writeBitString numberWords (P.advancePtr xPtr (i * numberWords)) x
    writeBitString numberWords (P.advancePtr sPtr (i * numberWords)) s
  pure p

destroyCnonbranching_terms :: Ptr Cnonbranching_terms -> IO ()
destroyCnonbranching_terms p
  | p /= nullPtr = do
      let freeIfNotNull x
            | x /= nullPtr = free x
            | otherwise = pure ()
      freeIfNotNull =<< [CU.exp| ls_hs_scalar const* { $(ls_hs_nonbranching_terms* p)->v } |]
      freeIfNotNull =<< [CU.exp| uint64_t const* { $(ls_hs_nonbranching_terms* p)->m } |]
      freeIfNotNull =<< [CU.exp| uint64_t const* { $(ls_hs_nonbranching_terms* p)->l } |]
      freeIfNotNull =<< [CU.exp| uint64_t const* { $(ls_hs_nonbranching_terms* p)->r } |]
      freeIfNotNull =<< [CU.exp| uint64_t const* { $(ls_hs_nonbranching_terms* p)->x } |]
      freeIfNotNull =<< [CU.exp| uint64_t const* { $(ls_hs_nonbranching_terms* p)->s } |]
      free p
  | otherwise = pure ()

getNonbranchingTerms
  :: forall t
   . (Typeable t, HasNonbranchingRepresentation (Generator Int (GeneratorType t)))
  => Operator t
  -> Vector NonbranchingTerm
getNonbranchingTerms operator =
  case nonbranchingRepresentation <$> opTermsFlat of
    -- NOTE: sorting based on nbtX is important!
    -- Terms with the same nbtX generate the same spin configuration. If they
    -- appear one after another, we can eliminate the duplicates easily. This
    -- is done in the off_diag kernel in Chapel.
    -- If duplicates are not eliminated, we might run into bufer overflows
    -- since we allocate buffers assuming that there are not duplicates.
    (Sum v) -> sortVectorBy (comparing (.nbtX)) $ G.filter ((/= 0) . (.nbtV)) v
  where
    numberSites = getNumberSites operator.opBasis
    flattenGenerator (Generator i g) = Generator (flattenIndex @t numberSites i) g
    opTermsFlat = fmap (fmap flattenGenerator) <$> unExpr operator.opTerms

getHilbertSpaceSectors :: forall t. (HasCallStack, IsBasis t) => Expr t -> Vector (Basis t)
getHilbertSpaceSectors expr = sortOnDim $ case particleDispatch @t of
  SpinTag -> G.fromList $ do
    h <- if hasU1 then Just <$> [0 .. n] else [Nothing]
    i <-
      -- Only apply spin inversion when exactly half the spins are up
      if hasZ2 && maybe True ((== n) . (2 *)) h
        then [Just 1, Just (-1)]
        else [Nothing]
    g <-
      -- Only apply symmetries when n > 1 and h /= 0 and h /= n
      if n > 1 && maybe True (\h' -> h' /= 0 || h' /= n) h
        then G.toList . abelianRepresentations . abelianSubgroup . exprPermutationGroup (Just n) $ expr
        else [emptyRepresentation]
    pure $ SpinBasis n h i g
  SpinlessFermionTag -> G.fromList $ do
    SpinlessFermionBasis n <$> (if hasU1 then Just <$> [0 .. n] else [Nothing])
  SpinfulFermionTag -> undefined
  where
    n = estimateNumberSites expr
    hasU1 = conservesNumberParticles expr
    hasZ2 = isInvariantUponSpinInversion expr
    sortOnDim = fmap snd . G.reverse . sortVectorBy (comparing fst) . fmap (\x -> (estimateHilbertSpaceDimension x, x))

-- newCoperator :: (IsBasis t, HasCallStack) => Maybe (Ptr Cbasis) -> Operator t -> IO (Ptr Coperator)
-- newCoperator maybeBasisPtr x = do
--   let (diag, offDiag) = G.partition nbtIsDiagonal (getNonbranchingTerms x)
--       numberBits = getNumberBits . opBasis $ x
--   diag_terms <- createCnonbranching_terms numberBits diag
--   off_diag_terms <- createCnonbranching_terms numberBits offDiag
--   basis <- case maybeBasisPtr of
--     Just p -> cloneCbasis p
--     Nothing -> newCbasis (opBasis x)
--   payload <- newStablePtr (SomeOperator (getParticleTag . opBasis $ x) x)
--   new $
--     Coperator
--       { coperator_refcount = AtomicCInt 1
--       , coperator_basis = basis
--       , coperator_off_diag_terms = off_diag_terms
--       , coperator_diag_terms = diag_terms
--       , -- coperator_apply_off_diag_cxt = nullPtr,
--         -- coperator_apply_diag_cxt = nullPtr,
--         coperator_haskell_payload = castStablePtrToPtr payload
--       }

withCoperator :: Ptr Coperator -> (SomeOperator -> IO a) -> IO a
withCoperator p f =
  f
    =<< (deRefStablePtr . castPtrToStablePtr)
    =<< [CU.exp| void* { $(ls_hs_operator* p)->base.haskell_payload } |]

foldCoperator :: (SomeOperator -> IO a) -> Ptr Coperator -> IO a
foldCoperator f p = withCoperator p f

newCoperator :: Maybe (Ptr Cbasis) -> Maybe (Ptr Cexpr) -> SomeOperator -> IO (Ptr Coperator)
newCoperator cBasis cExpr operator@(SomeOperator t op) = do
  let (diag, offDiag) = G.partition nbtIsDiagonal (foldSomeOperator getNonbranchingTerms operator)
      numberBits = foldSomeOperator (getNumberBits . opBasis) operator
  diag_terms <- newCnonbranching_terms numberBits diag
  off_diag_terms <- newCnonbranching_terms numberBits offDiag

  -- NOTE: important to initialize memory to 0 such that we don't have to manually initialize fields
  p <- callocBytes $ fromIntegral [CU.pure| size_t { sizeof(ls_hs_operator) } |]

  payload <- castStablePtrToPtr <$> newStablePtr operator
  cBasis' <- maybe (newCbasis op.opBasis) pure cBasis
  cExpr' <- maybe (newCexpr (SomeExpr t op.opTerms)) pure cExpr
  let incRefCountBasis = fromBool (isJust cBasis)
      incRefCountExpr = fromBool (isJust cExpr)
      maxNumberOffDiag =
        fromIntegral $
          getMaxNumberOffDiag (getNumberBits op.opBasis) (getHammingWeight op.opBasis) offDiag
      -- NOTE: disable computing a better upper bound on the number of terms because it's taking too long...
      -- TODO: fix me 😈
      maxNumberOffDiagEstimate = maxNumberOffDiag
  -- fromIntegral $ unsafeEstimateMaxNumberOffDiag (getNumberBits op.opBasis) (getHammingWeight op.opBasis) offDiag
  [CU.block| void {
    ls_hs_operator* p = $(ls_hs_operator* p);
    ls_hs_internal_object_init(&p->base, 1, $(void* payload));
    p->basis = $(ls_hs_basis* cBasis');
    if ($(bool incRefCountBasis)) {
      ls_hs_internal_object_inc_ref_count(&p->basis->base);
    }
    p->expr = $(ls_hs_expr* cExpr');
    if ($(bool incRefCountExpr)) {
      ls_hs_internal_object_inc_ref_count(&p->expr->base);
    }
    p->diag_terms = $(ls_hs_nonbranching_terms* diag_terms);
    p->off_diag_terms = $(ls_hs_nonbranching_terms* off_diag_terms);
    p->max_number_off_diag = $(int maxNumberOffDiag);
    p->max_number_off_diag_estimate = $(int maxNumberOffDiagEstimate);
  } |]
  pure p

ls_hs_create_operator :: Ptr Cbasis -> Ptr Cexpr -> IO CString
ls_hs_create_operator cBasis cExpr =
  withCbasis cBasis $ \basis ->
    withCexpr cExpr $ \expr ->
      newCencoded
        =<< rightM (newCoperator (Just cBasis) (Just cExpr)) (mkSomeOperator @(Either Text) basis expr)

ls_hs_destroy_operator :: Ptr Coperator -> IO ()
ls_hs_destroy_operator = ls_hs_destroy_object $ \p -> do
  destroyCnonbranching_terms =<< [CU.exp| ls_hs_nonbranching_terms* { $(ls_hs_operator* p)->diag_terms } |]
  destroyCnonbranching_terms =<< [CU.exp| ls_hs_nonbranching_terms* { $(ls_hs_operator* p)->off_diag_terms } |]
  ls_hs_destroy_expr =<< [CU.exp| ls_hs_expr* { $(ls_hs_operator* p)->expr } |]
  ls_hs_destroy_basis =<< [CU.exp| ls_hs_basis* { $(ls_hs_operator* p)->basis } |]

ls_hs_operator_max_number_off_diag :: Ptr Coperator -> IO CInt
ls_hs_operator_max_number_off_diag = foldCoperator $ foldSomeOperator $ \op -> do
  pure . fromIntegral $
    unsafeEstimateMaxNumberOffDiag (getNumberBits op.opBasis) (getHammingWeight op.opBasis) $
      G.filter (not . nbtIsDiagonal) (getNonbranchingTerms op)

ls_hs_operator_from_json :: CString -> IO CString
ls_hs_operator_from_json = newCencoded <=< rightM (newCoperator Nothing Nothing) <=< decodeCString @SomeOperator

ls_hs_operator_to_json :: Ptr Coperator -> IO CString
ls_hs_operator_to_json = foldCoperator newCencoded

ls_hs_expr_hilbert_space_sectors :: Ptr Cexpr -> IO CString
ls_hs_expr_hilbert_space_sectors = foldCexpr $ foldSomeExpr $ \expr ->
  newCencoded =<< G.mapM newCbasis (getHilbertSpaceSectors expr)

-- cloneCoperator :: (HasCallStack) => Ptr Coperator -> IO (Ptr Coperator)
-- cloneCoperator p = do
--   _ <- operatorIncRefCount p
--   pure p
--
-- destroyCoperator :: (HasCallStack) => Ptr Coperator -> IO ()
-- destroyCoperator p = do
--   refcount <- operatorDecRefCount p
--   when (refcount == 1) $ do
--     x <- peek p
--     destroyCbasis (coperator_basis x)
--     destroyCnonbranching_terms (coperator_off_diag_terms x)
--     destroyCnonbranching_terms (coperator_diag_terms x)
--     freeStablePtr . castPtrToStablePtr $ coperator_haskell_payload x
--     free p

-- maxNumberOffDiag :: IsBasis t => Operator t -> Int
-- maxNumberOffDiag op = Set.size $ Set.fromList . G.toList . G.map nbtX $ offDiag
--   where
--     (_, offDiag) = G.partition nbtIsDiagonal (getNonbranchingTerms op)

-- operatorToHypergraph :: IsBasis t => Operator t -> Hypergraph Int
-- operatorToHypergraph (Operator basis expr) = flattenHypergraph $ exprToHypergraph expr
--   where
--     flattenHypergraph (Hypergraph vertices hyperedges) =
--       Hypergraph
--         (Set.map (flattenIndex basis) vertices)
--         (Set.map (Set.map (flattenIndex basis)) hyperedges)
--
-- applyPermutation :: IsBasis t => Permutation -> Operator t -> Operator t
-- applyPermutation (unPermutation -> p) (Operator basis expr) =
--   Operator basis . simplifyExpr $ mapIndices remap expr
--   where
--     remap i = unFlattenIndex basis $ p G.! flattenIndex basis i
--
-- isInvariant :: IsBasis t => Operator t -> Permutation -> Bool
-- isInvariant operator permutation =
--   (applyPermutation permutation operator).opTerms == operator.opTerms
--
-- operatorSymmetryGroup :: IsBasis t => Operator t -> PermutationGroup
-- operatorSymmetryGroup operator = hypergraphAutomorphisms (isInvariant operator) . operatorToHypergraph $ operator
--
-- operatorAbelianRepresentations :: IsBasis t => Operator t -> [Vector Symmetry]
-- operatorAbelianRepresentations operator = abelianRepresentations g (mkMultiplicationTable g)
--   where
--     g = operatorSymmetryGroup operator
