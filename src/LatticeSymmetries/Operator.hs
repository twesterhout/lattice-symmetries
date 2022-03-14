{-# LANGUAGE CApiFFI #-}
{-# LANGUAGE UndecidableInstances #-}

module LatticeSymmetries.Operator where

import Control.Exception.Safe (bracket)
import Data.Bits
import Data.Complex
import qualified Data.Primitive.Ptr as P
import Data.Vector (Vector)
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Storable as S
import qualified Data.Vector.Storable.Mutable as SM
import Foreign.C.Types
import Foreign.Marshal.Alloc
import Foreign.Marshal.Utils (fromBool, new)
import Foreign.Ptr
import Foreign.StablePtr
import Foreign.Storable
import LatticeSymmetries.Algebra
import LatticeSymmetries.Basis
import LatticeSymmetries.BitString
import LatticeSymmetries.ComplexRational
import LatticeSymmetries.Generator
import LatticeSymmetries.NonbranchingTerm
import LatticeSymmetries.Utils (loopM)
import System.IO (hPutStrLn)
import System.IO.Unsafe (unsafePerformIO)
import Text.Printf
import Prelude hiding (Sum)

data Operator (t :: ParticleTy) = Operator
  { opBasis :: !(Basis t),
    opTerms :: !(Polynomial ComplexRational (Generator (IndexType t) (GeneratorType t)))
  }

deriving instance (Show (IndexType t), Show (GeneratorType t)) => Show (Operator t)

opTermsFlat :: Operator t -> Polynomial ComplexRational (Generator Int (GeneratorType t))
opTermsFlat operator =
  fmap (fmap (fmap (\(Generator i g) -> Generator (flattenIndex basis i) g))) (opTerms operator)
  where
    basis = opBasis operator

instance
  ( Enum (GeneratorType t),
    Bounded (GeneratorType t),
    HasMatrixRepresentation (GeneratorType t),
    Algebra (GeneratorType t),
    Ord (IndexType t)
  ) =>
  Num (Operator t)
  where
  (+) a b
    | opBasis a == opBasis b = Operator (opBasis a) (simplify $ opTerms a + opTerms b)

-- getNumberTerms :: Operator c basis -> Int
-- getNumberTerms operator = let (Sum v) = opTerms operator in G.length v

getNonbranchingTerms ::
  HasNonbranchingRepresentation (Generator Int (GeneratorType t)) =>
  Operator t ->
  Vector NonbranchingTerm
getNonbranchingTerms operator =
  case nonbranchingRepresentation <$> opTermsFlat operator of
    (Sum v) -> v

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

data Cnonbranching_terms = Cnonbranching_terms
  { cnonbranching_terms_number_terms :: {-# UNPACK #-} !CInt,
    cnonbranching_terms_number_bits :: {-# UNPACK #-} !CInt,
    cnonbranching_terms_v :: {-# UNPACK #-} !(Ptr Cscalar),
    cnonbranching_terms_m :: {-# UNPACK #-} !(Ptr Word64),
    cnonbranching_terms_l :: {-# UNPACK #-} !(Ptr Word64),
    cnonbranching_terms_r :: {-# UNPACK #-} !(Ptr Word64),
    cnonbranching_terms_x :: {-# UNPACK #-} !(Ptr Word64),
    cnonbranching_terms_s :: {-# UNPACK #-} !(Ptr Word64)
  }
  deriving stock (Show)

instance Storable Cnonbranching_terms where
  sizeOf _ = 56 -- 2 * 4 + 6 * 8
  alignment _ = 8
  peek p =
    Cnonbranching_terms
      <$> peekByteOff p 0
      <*> peekByteOff p 4
      <*> peekByteOff p 8
      <*> peekByteOff p 16
      <*> peekByteOff p 24
      <*> peekByteOff p 32
      <*> peekByteOff p 40
      <*> peekByteOff p 48
  poke p (Cnonbranching_terms nTerms nBits v m l r x s) = do
    pokeByteOff p 0 nTerms
    pokeByteOff p 4 nBits
    pokeByteOff p 8 v
    pokeByteOff p 16 m
    pokeByteOff p 24 l
    pokeByteOff p 32 r
    pokeByteOff p 40 x
    pokeByteOff p 48 s

--
-- typedef struct ls_hs_operator {
--   ls_hs_basis const *basis;
--   ls_hs_nonbranching_terms const *off_diag_terms;
--   ls_hs_nonbranching_terms const *diag_terms;
--   bool needs_projection;
-- } ls_hs_operator;

data Coperator = Coperator
  { coperator_basis :: {-# UNPACK #-} !(Ptr Cbasis),
    coperator_off_diag_terms :: {-# UNPACK #-} !(Ptr Cnonbranching_terms),
    coperator_diag_terms :: {-# UNPACK #-} !(Ptr Cnonbranching_terms),
    coperator_needs_projection :: {-# UNPACK #-} !CBool,
    coperator_haskell_payload :: {-# UNPACK #-} !(Ptr ())
  }

instance Storable Coperator where
  sizeOf _ = 40
  alignment _ = 8
  peek p =
    Coperator
      <$> peekByteOff p 0
      <*> peekByteOff p 8
      <*> peekByteOff p 16
      <*> peekByteOff p 24
      <*> peekByteOff p 32
  poke p (Coperator basis off_diag_terms diag_terms needs_projection payload) = do
    pokeByteOff p 0 basis
    pokeByteOff p 8 off_diag_terms
    pokeByteOff p 16 diag_terms
    pokeByteOff p 24 needs_projection
    pokeByteOff p 32 payload

createCoperator ::
  forall t.
  ( HasNonbranchingRepresentation (Generator Int (GeneratorType t))
  ) =>
  Operator t ->
  IO (Ptr Coperator)
createCoperator operator = do
  let (diag, offDiag) = G.partition isNonbranchingTermDiagonal (getNonbranchingTerms operator)
      numberBits = getNumberBits (opBasis operator)
  basis <- createCbasis (opBasis operator)
  diag_terms <- createCnonbranching_terms numberBits diag
  off_diag_terms <- createCnonbranching_terms numberBits offDiag
  p <- castStablePtrToPtr <$> newStablePtr operator
  new $
    Coperator
      basis
      off_diag_terms
      diag_terms
      (fromBool False)
      p

destroyCoperator :: HasCallStack => Ptr Coperator -> IO ()
destroyCoperator p
  | p /= nullPtr = do
    (Coperator basis off_diag_terms diag_terms _ payload) <- peek p
    !particleType <- cbasis_particle_type <$> peek basis
    destroyCbasis basis
    destroyCnonbranching_terms off_diag_terms
    destroyCnonbranching_terms diag_terms
    let freeHaskellPayload :: forall (t :: ParticleTy). Proxy t -> IO ()
        freeHaskellPayload _ =
          freeStablePtr (castPtrToStablePtr payload :: StablePtr (Operator t))
    withParticleType particleType freeHaskellPayload
    free p
  | otherwise = error "should not happen"

peekNumberOffDiagTerms :: Ptr Coperator -> IO Int
peekNumberOffDiagTerms p = do
  p' <- coperator_off_diag_terms <$> peek p
  if p' /= nullPtr
    then (fromIntegral . cnonbranching_terms_number_terms) <$> peek p'
    else pure 0

applyOperator ::
  forall t.
  HasNonbranchingRepresentation (Generator Int (GeneratorType t)) =>
  Operator t ->
  BasisState ->
  Vector (ComplexRational, BasisState)
applyOperator operator (BasisState _ x) = unsafePerformIO $
  bracket (createCoperator operator) destroyCoperator $ \c_operator -> do
    numberTerms <- peekNumberOffDiagTerms c_operator
    hPutStrLn stderr $ printf "numberTerms = %d" numberTerms
    let numberWords = getNumberWords (opBasis operator)
        -- numberTerms = getNumberTerms operator
        αs = wordsFromBitString numberWords x
        αsStride = fromIntegral numberWords
        βsStride = fromIntegral numberWords
        batchSize = 1
    βs <- SM.unsafeNew ((numberTerms + 1) * numberWords)
    cs <- SM.unsafeNew (numberTerms + 1)
    S.unsafeWith αs $ \αsPtr ->
      SM.unsafeWith βs $ \βsPtr ->
        SM.unsafeWith cs $ \csPtr ->
          do
            writeBitString
              numberWords
              (P.advancePtr βsPtr (numberTerms * numberWords))
              x
            ls_hs_operator_apply_diag_kernel
              c_operator
              batchSize
              αsPtr
              αsStride
              (plusPtr csPtr (numberTerms * sizeOf (undefined :: Cscalar)))
            ls_hs_operator_apply_off_diag_kernel
              c_operator
              batchSize
              αsPtr
              αsStride
              βsPtr
              βsStride
              csPtr
    (hPutStrLn stderr . show) =<< G.unsafeFreeze cs
    cs' <- G.map fromComplexDouble <$> G.convert <$> G.unsafeFreeze cs
    (hPutStrLn stderr . show) cs'
    (hPutStrLn stderr . show) =<< G.unsafeFreeze βs
    let numberBits = getNumberBits (opBasis operator)
    βs' <-
      SM.unsafeWith βs $ \p ->
        G.generateM (numberTerms + 1) $ \i ->
          BasisState numberBits
            <$> readBitString numberWords (P.advancePtr p (i * numberWords))
    hPutStrLn stderr $ "Done"
    (hPutStrLn stderr . show) βs'
    let !r = G.filter (\(c, _) -> c /= 0) $ G.zip cs' βs'
    (hPutStrLn stderr . show) r
    pure $ r

-- void ls_hs_operator_apply_diag_kernel(ls_hs_operator const *op,
--                                       ptrdiff_t batch_size,
--                                       uint64_t const *alphas,
--                                       ptrdiff_t alphas_stride,
--                                       ls_hs_scalar *coeffs);
foreign import capi unsafe "lattice_symmetries_haskell.h ls_hs_operator_apply_diag_kernel"
  ls_hs_operator_apply_diag_kernel :: Ptr Coperator -> CPtrdiff -> Ptr Word64 -> CPtrdiff -> Ptr Cscalar -> IO ()

-- void ls_hs_operator_apply_off_diag_kernel(
--     ls_hs_operator const *op, ptrdiff_t batch_size, uint64_t const *alphas,
--     ptrdiff_t alphas_stride, uint64_t *betas, ptrdiff_t betas_stride,
--     ls_hs_scalar *coeffs);
foreign import capi unsafe "lattice_symmetries_haskell.h ls_hs_operator_apply_off_diag_kernel"
  ls_hs_operator_apply_off_diag_kernel :: Ptr Coperator -> CPtrdiff -> Ptr Word64 -> CPtrdiff -> Ptr Word64 -> CPtrdiff -> Ptr Cscalar -> IO ()
