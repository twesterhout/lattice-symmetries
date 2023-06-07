{-# LANGUAGE CApiFFI #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE UndecidableInstances #-}

module LatticeSymmetries.Operator
  ( Operator (..),
    mkOperator,
    SomeOperator (..),
    withSomeOperator,
    foldSomeOperator,
    maxNumberOffDiag,

    -- ** Helpers for FFI
    newCoperator,
    cloneCoperator,
    destroyCoperator,
    withCoperator,
  )
where

import qualified Data.Primitive.Ptr as P
import qualified Data.Set as Set
import Data.Vector (Vector)
import qualified Data.Vector.Generic as G
import Foreign.ForeignPtr
import Foreign.Marshal
import Foreign.Ptr
import Foreign.StablePtr
import Foreign.Storable
import GHC.ForeignPtr
import qualified GHC.Show
import LatticeSymmetries.Algebra
import LatticeSymmetries.Basis
import LatticeSymmetries.BitString
import LatticeSymmetries.ComplexRational
import LatticeSymmetries.Expr
import LatticeSymmetries.FFI
import LatticeSymmetries.Generator
import LatticeSymmetries.NonbranchingTerm
import LatticeSymmetries.Utils
import System.IO.Unsafe (unsafePerformIO)
import Prelude hiding (Product, Sum)

data Operator (t :: ParticleTy) = Operator
  { opBasis :: !(Basis t),
    opTerms :: !(Expr t)
    -- Polynomial ComplexRational (Generator (IndexType t) (GeneratorType t)))
  }

deriving stock instance
  (Show (Basis t), Show (IndexType t), Show (GeneratorType t)) =>
  Show (Operator t)

-- data Operator (t :: ParticleTy) = Operator
--   { opHeader :: !(OperatorHeader t),
--     opContents :: !(ForeignPtr Coperator)
--   }

-- deriving stock instance (Show (OperatorHeader t)) => Show (Operator t)

data SomeOperator where
  SomeOperator :: IsBasis t => ParticleTag t -> Operator t -> SomeOperator

instance Show SomeOperator where
  show (SomeOperator SpinTag x) = show x
  show (SomeOperator SpinlessFermionTag x) = show x
  show (SomeOperator SpinfulFermionTag x) = show x

withSomeOperator :: SomeOperator -> (forall t. IsBasis t => Operator t -> a) -> a
withSomeOperator x f = case x of
  SomeOperator _ operator -> f operator
{-# INLINE withSomeOperator #-}

foldSomeOperator :: (forall t. IsBasis t => Operator t -> a) -> SomeOperator -> a
foldSomeOperator f x = case x of
  SomeOperator _ operator -> f operator
{-# INLINE foldSomeOperator #-}

-- withSameBasis :: IsBasis t => OperatorHeader t -> OperatorHeader t -> (Basis t -> a) -> a
-- withSameBasis a b f
--   | opBasis a == opBasis b = f (opBasis a)
--   | otherwise = error "expected operators defined on the same basis"

-- instance IsBasis t => Num (OperatorHeader t) where
--   (+) a b = withSameBasis a b $ \basis -> OperatorHeader basis (opTerms a + opTerms b)
--   (-) a b = withSameBasis a b $ \basis -> OperatorHeader basis (opTerms a - opTerms b)
--   (*) a b = withSameBasis a b $ \basis -> OperatorHeader basis (opTerms a * opTerms b)
--   negate a = OperatorHeader (opBasis a) (negate (opTerms a))
--   abs a = OperatorHeader (opBasis a) (abs (opTerms a))
--   signum a = OperatorHeader (opBasis a) (signum (opTerms a))
--   fromInteger _ = error "Num instance of OperatorHeader does not implement fromInteger"

-- instance IsBasis t => CanScale ComplexRational (OperatorHeader t) where
--   scale c a = OperatorHeader (opBasis a) (c `scale` opTerms a)

-- instance IsBasis t => Num (Operator t) where
--   (+) a b = operatorFromHeader $ opHeader a + opHeader b
--   (-) a b = operatorFromHeader $ opHeader a - opHeader b
--   (*) a b = operatorFromHeader $ opHeader a * opHeader b
--   negate a = operatorFromHeader $ negate (opHeader a)
--   abs a = operatorFromHeader $ abs (opHeader a)
--   signum a = operatorFromHeader $ signum (opHeader a)
--   fromInteger n = operatorFromHeader $ fromInteger n

-- instance IsBasis t => CanScale ComplexRational (Operator t) where
--   scale c a = operatorFromHeader $ scale c (opHeader a)

mkOperator :: IsBasis t => Basis t -> Expr t -> Operator t
mkOperator basis expr = Operator basis expr
{-# NOINLINE mkOperator #-}

-- getNumberTerms :: Operator c basis -> Int
-- getNumberTerms operator = let (Sum v) = opTerms operator in G.length v

-- opIsIdentity :: IsBasis t => Operator t -> Bool
-- opIsIdentity operator = G.length terms == 1 && isId (G.head terms)
--   where
--     isId (Scaled _ (Product p)) = G.length p == 1 && isIdentity (G.head p)
--     terms = unSum . opTerms . opHeader $ operator

-- opIsHermitian :: IsBasis t => Operator t -> Bool
-- opIsHermitian operator = terms == conjugateExpr terms
--   where
--     terms = opTerms $ opHeader operator

-- conjugateOperator :: IsBasis t => Operator t -> Operator t
-- conjugateOperator operator =
--   operatorFromHeader $
--     OperatorHeader basis (conjugateExpr terms)
--   where
--     basis = opBasis $ opHeader operator
--     terms = opTerms $ opHeader operator

-- opHermitianConjugate :: IsBasis t => Operator t -> Operator t
-- opHermitianConjugate operator = operatorFromHeader $
--   where

-- opIsHermitian :: IsBasis t => Symmetric

getNonbranchingTerms ::
  HasNonbranchingRepresentation (Generator Int (GeneratorType t)) =>
  Operator t ->
  Vector NonbranchingTerm
getNonbranchingTerms operator =
  case nonbranchingRepresentation <$> opTermsFlat operator of
    (Sum v) -> v
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
      let c_terms =
            Cnonbranching_terms
              (fromIntegral numberTerms)
              (fromIntegral numberBits)
              vPtr
              mPtr
              lPtr
              rPtr
              xPtr
              sPtr
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

--
-- typedef struct ls_hs_operator {
--   ls_hs_basis const *basis;
--   ls_hs_nonbranching_terms const *off_diag_terms;
--   ls_hs_nonbranching_terms const *diag_terms;
--   bool needs_projection;
-- } ls_hs_operator;

-- foreign import ccall safe "ls_hs_destroy_basis_v2"
--   ls_hs_destroy_basis :: Ptr Cbasis -> IO ()

-- foreign import ccall unsafe "ls_internal_create_apply_off_diag_kernel_data"
--   ls_internal_create_apply_off_diag_kernel_data :: Ptr Cnonbranching_terms -> IO (Ptr Coperator_kernel_data)

-- foreign import ccall unsafe "ls_internal_create_apply_diag_kernel_data"
--   ls_internal_create_apply_diag_kernel_data :: Ptr Cnonbranching_terms -> IO (Ptr Coperator_kernel_data)

-- foreign import ccall unsafe "ls_internal_destroy_operator_kernel_data"
--   ls_internal_destroy_operator_kernel_data :: Ptr Coperator_kernel_data -> IO ()

newCoperator :: (IsBasis t, HasCallStack) => Maybe (Ptr Cbasis) -> Operator t -> IO (Ptr Coperator)
newCoperator maybeBasisPtr x = do
  -- logDebug' $ "newCoperator"
  let (diag, offDiag) = G.partition nbtIsDiagonal (getNonbranchingTerms x)
      numberBits = getNumberBits . opBasis $ x
  diag_terms <- createCnonbranching_terms numberBits diag
  off_diag_terms <- createCnonbranching_terms numberBits offDiag
  basis <- case maybeBasisPtr of
    Just p -> cloneCbasis p
    Nothing -> newCbasis (opBasis x)
  -- logDebug' $ "basis is " <> show basis
  -- borrowCbasis (opBasis x)
  payload <- newStablePtr (SomeOperator (getParticleTag . opBasis $ x) x)
  r <-
    new $
      Coperator
        { coperator_refcount = 1,
          coperator_basis = basis,
          coperator_off_diag_terms = off_diag_terms,
          coperator_diag_terms = diag_terms,
          -- coperator_apply_off_diag_cxt = nullPtr,
          -- coperator_apply_diag_cxt = nullPtr,
          coperator_haskell_payload = castStablePtrToPtr payload
        }
  pure r

withCoperator :: Ptr Coperator -> (SomeOperator -> IO a) -> IO a
withCoperator p f = f =<< (deRefStablePtr . castPtrToStablePtr . coperator_haskell_payload) =<< peek p

cloneCoperator :: HasCallStack => Ptr Coperator -> IO (Ptr Coperator)
cloneCoperator p = do
  basis <- coperator_basis <$> peek p
  _ <- operatorIncRefCount p
  -- payload <- newStablePtr =<< deRefStablePtr (castPtrToStablePtr . coperator_haskell_payload $ x)
  -- new $ x {coperator_haskell_payload = castStablePtrToPtr payload}
  pure p

destroyCoperator :: HasCallStack => Ptr Coperator -> IO ()
destroyCoperator p = do
  refcount <- operatorDecRefCount p
  when (refcount == 1) $ do
    x <- peek p
    destroyCbasis (coperator_basis x)
    -- logDebug' $ "destroying off-diag terms ..."
    destroyCnonbranching_terms (coperator_off_diag_terms x)
    -- logDebug' $ "destroying diag terms ..."
    destroyCnonbranching_terms (coperator_diag_terms x)
    -- logDebug' $ "freeing stable ptr ..."
    freeStablePtr . castPtrToStablePtr $ coperator_haskell_payload x
    -- logDebug' $ "freeing ptr ..."
    free p

-- addForeignPtrConcFinalizer fp (ls_internal_destroy_operator_kernel_data apply_off_diag_cxt)
-- addForeignPtrConcFinalizer fp (ls_internal_destroy_operator_kernel_data apply_diag_cxt)
-- pure $ Operator x fp

-- operatorFromHeader ::
--   HasNonbranchingRepresentation (Generator Int (GeneratorType t)) =>
--   OperatorHeader t ->
--   Operator t
-- operatorFromHeader x = unsafePerformIO $ do
--   let (diag, offDiag) = G.partition nbtIsDiagonal (getNonbranchingTerms x)
--       numberBits = getNumberBits . basisHeader . opBasis $ x
--   diag_terms <- createCnonbranching_terms numberBits diag
--   off_diag_terms <- createCnonbranching_terms numberBits offDiag
--   -- apply_diag_cxt <- ls_internal_create_apply_diag_kernel_data diag_terms
--   -- apply_off_diag_cxt <- ls_internal_create_apply_off_diag_kernel_data off_diag_terms
--   basis <- borrowCbasis (opBasis x)
--
--   fp <- mallocForeignPtr
--   addForeignPtrConcFinalizer fp (destroyCbasis basis)
--   -- addForeignPtrConcFinalizer fp (ls_internal_destroy_operator_kernel_data apply_off_diag_cxt)
--   -- addForeignPtrConcFinalizer fp (ls_internal_destroy_operator_kernel_data apply_diag_cxt)
--   addForeignPtrConcFinalizer fp (destroyCnonbranching_terms off_diag_terms)
--   addForeignPtrConcFinalizer fp (destroyCnonbranching_terms diag_terms)
--   withForeignPtr fp $ \ptr ->
--     poke ptr $
--       Coperator
--         { coperator_refcount = 0,
--           coperator_basis = basis,
--           coperator_off_diag_terms = off_diag_terms,
--           coperator_diag_terms = diag_terms,
--           coperator_apply_off_diag_cxt = nullPtr,
--           coperator_apply_diag_cxt = nullPtr,
--           coperator_haskell_payload = nullPtr
--         }
--   pure $ Operator x fp
-- {-# NOINLINE operatorFromHeader #-}

-- borrowCoperator :: Operator t -> IO (Ptr Coperator)
-- borrowCoperator operator = do
--   payload <- newStablePtr operator
--   withForeignPtr (opContents operator) $ \ptr -> do
--     _ <- operatorIncRefCount ptr
--     operatorPokePayload ptr payload
--   pure $ unsafeForeignPtrToPtr (opContents operator)

-- destroyCoperator :: Ptr Coperator -> IO ()
-- destroyCoperator ptr = do
--   refcount <- operatorDecRefCount ptr
--   when (refcount == 1) $
--     freeStablePtr =<< operatorPeekPayload ptr

--
-- createCoperator ::
--   forall t.
--   ( HasNonbranchingRepresentation (Generator Int (GeneratorType t))
--   ) =>
--   Operator t ->
--   IO (Ptr Coperator)
-- createCoperator operator = do
--   let (diag, offDiag) = G.partition isNonbranchingTermDiagonal (getNonbranchingTerms operator)
--       numberBits = getNumberBits (opBasis operator)
--   basis <- createCbasis (opBasis operator)
--   diag_terms <- createCnonbranching_terms numberBits diag
--   off_diag_terms <- createCnonbranching_terms numberBits offDiag
--   p <- castStablePtrToPtr <$> newStablePtr operator
--   putStrLn $ "Created Coperator: " <> show p
--   new $
--     Coperator
--       basis
--       off_diag_terms
--       diag_terms
--       (fromBool False)
--       p

-- destroyCoperator :: HasCallStack => Ptr Coperator -> IO ()
-- destroyCoperator p
--   | p /= nullPtr = do
--     (Coperator basis off_diag_terms diag_terms _ payload) <- peek p
--     putStrLn $ "Destroying Coperator: " <> show payload
--     !particleType <- cbasis_particle_type <$> peek basis
--     destroyCbasis basis
--     destroyCnonbranching_terms off_diag_terms
--     destroyCnonbranching_terms diag_terms
--     let freeHaskellPayload :: forall (t :: ParticleTy). Proxy t -> IO ()
--         freeHaskellPayload _ =
--           freeStablePtr (castPtrToStablePtr payload :: StablePtr (Operator t))
--     withParticleType particleType freeHaskellPayload
--     free p
--   | otherwise = error "should not happen"

-- withReconstructedOperator ::
--   forall a.
--   Ptr Coperator ->
--   (forall (t :: ParticleTy). IsBasis t => Operator t -> IO a) ->
--   IO a
-- withReconstructedOperator p action = do
--   let run :: forall (t :: ParticleTy). IsBasis t => Proxy t -> IO a
--       run _ = (action :: Operator t -> IO a) =<< deRefStablePtr =<< operatorPeekPayload p
--   tp <- operatorPeekParticleType p
--   withParticleType tp run

-- withSameTypeAs ::
--   forall t a.
--   HasCallStack =>
--   Operator t ->
--   ((forall t'. Operator t' -> a) -> a) ->
--   (Operator t -> a) ->
--   a
-- withSameTypeAs a _with _action = _with f
--   where
--     f :: forall t'. Operator t' -> a
--     f b = case (basisHeader . opBasis . opHeader $ a, basisHeader . opBasis . opHeader $ b) of
--       ((SpinHeader _ _ _ _), (SpinHeader _ _ _ _)) -> _action b
--       ((SpinfulFermionHeader _ _), (SpinfulFermionHeader _ _)) -> _action b
--       ((SpinlessFermionHeader _ _), (SpinlessFermionHeader _ _)) -> _action b
--       _ -> error "operators have different types"

maxNumberOffDiag :: IsBasis t => Operator t -> Int
maxNumberOffDiag op = Set.size $ Set.fromList . G.toList . G.map nbtX $ offDiag
  where
    (_, offDiag) = G.partition nbtIsDiagonal (getNonbranchingTerms op)

-- applyOperator ::
--   forall t.
--   HasNonbranchingRepresentation (Generator Int (GeneratorType t)) =>
--   Operator t ->
--   BasisState ->
--   Vector (ComplexRational, BasisState)
-- applyOperator operator (BasisState _ x) = unsafePerformIO $
--   bracket (createCoperator operator) destroyCoperator $ \c_operator -> do
--     numberTerms <- peekNumberOffDiagTerms c_operator
--     hPutStrLn stderr $ printf "numberTerms = %d" numberTerms
--     let numberWords = getNumberWords (opBasis operator)
--         -- numberTerms = getNumberTerms operator
--         αs = wordsFromBitString numberWords x
--         αsStride = fromIntegral numberWords
--         βsStride = fromIntegral numberWords
--         batchSize = 1
--     βs <- SM.unsafeNew ((numberTerms + 1) * numberWords)
--     cs <- SM.unsafeNew (numberTerms + 1)
--     S.unsafeWith αs $ \αsPtr ->
--       SM.unsafeWith βs $ \βsPtr ->
--         SM.unsafeWith cs $ \csPtr ->
--           do
--             writeBitString
--               numberWords
--               (P.advancePtr βsPtr (numberTerms * numberWords))
--               x
--             ls_hs_operator_apply_diag_kernel
--               c_operator
--               batchSize
--               αsPtr
--               αsStride
--               (plusPtr csPtr (numberTerms * sizeOf (undefined :: Cscalar)))
--             ls_hs_operator_apply_off_diag_kernel
--               c_operator
--               batchSize
--               αsPtr
--               αsStride
--               βsPtr
--               βsStride
--               csPtr
--     (hPutStrLn stderr . show) =<< G.unsafeFreeze cs
--     cs' <- G.map fromComplexDouble <$> G.convert <$> G.unsafeFreeze cs
--     (hPutStrLn stderr . show) cs'
--     (hPutStrLn stderr . show) =<< G.unsafeFreeze βs
--     let numberBits = getNumberBits (opBasis operator)
--     βs' <-
--       SM.unsafeWith βs $ \p ->
--         G.generateM (numberTerms + 1) $ \i ->
--           BasisState numberBits
--             <$> readBitString numberWords (P.advancePtr p (i * numberWords))
--     hPutStrLn stderr $ "Done"
--     (hPutStrLn stderr . show) βs'
--     let !r = G.filter (\(c, _) -> c /= 0) $ G.zip cs' βs'
--     (hPutStrLn stderr . show) r
--     pure $ r

-- void ls_hs_operator_apply_diag_kernel(ls_hs_operator const *op,
--                                       ptrdiff_t batch_size,
--                                       uint64_t const *alphas,
--                                       ptrdiff_t alphas_stride,
--                                       ls_hs_scalar *coeffs);
-- foreign import capi unsafe "lattice_symmetries_haskell.h ls_hs_operator_apply_diag_kernel"
--   ls_hs_operator_apply_diag_kernel :: Ptr Coperator -> CPtrdiff -> Ptr Word64 -> CPtrdiff -> Ptr Cscalar -> IO ()

-- void ls_hs_operator_apply_off_diag_kernel(
--     ls_hs_operator const *op, ptrdiff_t batch_size, uint64_t const *alphas,
--     ptrdiff_t alphas_stride, uint64_t *betas, ptrdiff_t betas_stride,
--     ls_hs_scalar *coeffs);
-- foreign import capi unsafe "lattice_symmetries_haskell.h ls_hs_operator_apply_off_diag_kernel"
--   ls_hs_operator_apply_off_diag_kernel :: Ptr Coperator -> CPtrdiff -> Ptr Word64 -> CPtrdiff -> Ptr Word64 -> CPtrdiff -> Ptr Cscalar -> IO ()

-- foreign import capi unsafe "lattice_symmetries_haskell.h ls_hs_evaluate_wavefunction_via_statevector"
--   ls_hs_evaluate_wavefunction_via_statevector :: Ptr Cbasis_kernels -> CPtrdiff -> Ptr Word64 -> CPtrdiff -> Ptr () -> CSize -> Ptr () -> IO ()
-- ls_hs_create_spin_basis_from_yaml :: HasCallStack => CString -> IO (Ptr Cbasis)
-- ls_hs_create_spin_basis_from_yaml cFilename = do
--   filename <- peekCString cFilename
--   logDebug' $ "Loading Basis from " <> show filename <> " ..."
--   r <- decodeFileWithWarnings filename
--   case r of
--     Left e -> do
--       logError' $ toText $ prettyPrintParseException e
--       pure nullPtr
--     Right (warnings, (BasisOnlyConfig (basis :: Basis 'SpinTy))) -> do
--       mapM_ (logWarning' . show) warnings
--       borrowCbasis basis
