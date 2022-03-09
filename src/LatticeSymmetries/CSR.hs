{-# LANGUAGE CApiFFI #-}
{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE UndecidableInstances #-}

module LatticeSymmetries.CSR where

-- import Data.Binary (Binary (..))

-- import Data.Vector.Binary

-- import Control.Exception.Safe (bracket, impureThrow, throwIO)
import Control.Monad.Primitive (PrimMonad, PrimState)
import Control.Monad.ST
-- import qualified Control.Monad.ST.Unsafe (unsafeIOToST)
-- import Data.Aeson
-- import Data.Aeson.Types (typeMismatch)

-- import Data.Scientific (toRealFloat)

import Data.Bits (Bits, toIntegralSized)
import Data.Complex
import qualified Data.List
import Data.Type.Equality
-- import qualified Data.List.NonEmpty as NonEmpty

-- import qualified Data.Vector.Fusion.Stream.Monadic as Stream

-- import qualified Data.Vector.Storable as S
import qualified Data.Vector
import qualified Data.Vector as B
import qualified Data.Vector.Algorithms.Intro as Intro
import Data.Vector.Fusion.Bundle (Bundle)
import qualified Data.Vector.Fusion.Bundle as Bundle (inplace)
import qualified Data.Vector.Fusion.Bundle.Monadic as Bundle
import Data.Vector.Fusion.Bundle.Size (Size (..), toMax)
import Data.Vector.Fusion.Stream.Monadic (Step (..), Stream (..))
import qualified Data.Vector.Fusion.Util (unId)
import Data.Vector.Generic ((!))
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Generic.Mutable as GM
import qualified Data.Vector.Storable as S
import qualified Data.Vector.Unboxed
-- import qualified Data.Vector.Storable.Mutable as SM
-- import Data.Vector.Unboxed (Unbox)
-- import qualified Data.Vector.Unboxed as U
-- import Data.Yaml (decodeFileWithWarnings)
-- import Foreign.C.String (CString, peekCString)
import Foreign.C.Types (CInt (..), CUInt (..))
-- import Foreign.ForeignPtr
-- import Foreign.ForeignPtr.Unsafe (unsafeForeignPtrToPtr)
-- import Foreign.Marshal.Alloc (alloca, free, malloc, mallocBytes)
-- import Foreign.Marshal.Array (newArray, withArrayLen)
import Foreign.Marshal.Utils (with)
import Foreign.Ptr (FunPtr, Ptr, castPtr)
-- import Foreign.StablePtr
import Foreign.Storable (Storable (..))
import qualified GHC.Exts as GHC (IsList (..))
import GHC.TypeLits
import qualified GHC.TypeLits as GHC
-- import qualified GHC.ForeignPtr as GHC (Finalizers (..), ForeignPtr (..), ForeignPtrContents (..))
-- import GHC.Generics
-- import qualified GHC.IORef as GHC (atomicSwapIORef)
-- import GHC.Prim
-- import qualified GHC.Ptr as GHC (Ptr (..))
-- import qualified Language.C.Inline as C
-- import qualified Language.C.Inline.Unsafe as CU
import LatticeSymmetries.Dense
-- import LatticeSymmetries.IO
-- import LatticeSymmetries.Types
import qualified System.IO.Unsafe
import qualified Unsafe.Coerce
-- import qualified System.Mem.Weak
import Prelude hiding (group, product, sort)

-- typedef struct ls_csr_matrix {
--   unsigned *offsets;
--   unsigned *columns;
--   _Complex double *off_diag_elements;
--   _Complex double *diag_elements;
--   unsigned dimension;
--   unsigned number_nonzero;
-- } ls_csr_matrix;
data Ccsr_matrix = Ccsr_matrix
  { c_csr_matrix_offsets :: !(Ptr CUInt),
    c_csr_matrix_columns :: !(Ptr CUInt),
    c_csr_matrix_off_diag_elements :: !(Ptr (Complex Double)),
    c_csr_matrix_diag_elements :: !(Ptr (Complex Double)),
    c_csr_matrix_dimension :: !CUInt,
    c_csr_matrix_number_nonzero :: !CUInt
  }

instance Storable Ccsr_matrix where
  alignment _ = 8
  sizeOf _ = 40
  peek p =
    Ccsr_matrix
      <$> peekByteOff p 0
      <*> peekByteOff p 8
      <*> peekByteOff p 16
      <*> peekByteOff p 24
      <*> peekByteOff p 32
      <*> peekByteOff p 36
  poke p x = do
    pokeByteOff p 0 (c_csr_matrix_offsets x)
    pokeByteOff p 8 (c_csr_matrix_columns x)
    pokeByteOff p 16 (c_csr_matrix_off_diag_elements x)
    pokeByteOff p 24 (c_csr_matrix_diag_elements x)
    pokeByteOff p 32 (c_csr_matrix_dimension x)
    pokeByteOff p 36 (c_csr_matrix_number_nonzero x)

data CsrMatrix
  = CsrMatrix
      !(S.Vector CUInt)
      !(S.Vector CUInt)
      !(S.Vector (Complex Double))
      !(S.Vector (Complex Double))
  deriving stock (Show, Eq)

withCsrMatrix :: CsrMatrix -> (Ptr Ccsr_matrix -> IO a) -> IO a
withCsrMatrix m@(CsrMatrix offsets columns offDiagElems diagElems) action =
  S.unsafeWith offsets $ \offsetsPtr ->
    S.unsafeWith columns $ \columnsPtr ->
      S.unsafeWith offDiagElems $ \offDiagElemsPtr ->
        S.unsafeWith diagElems $ \diagElemsPtr ->
          with (Ccsr_matrix offsetsPtr columnsPtr offDiagElemsPtr diagElemsPtr c_dimension c_number_non_zero) action
  where
    c_dimension = fromIntegral (csrDim m)
    c_number_non_zero = fromIntegral (csrNumberNonZero m)

unsafeNewCsrMatrix :: Int -> Int -> IO CsrMatrix
unsafeNewCsrMatrix dimension nnz =
  pure $
    CsrMatrix
      (G.replicate (dimension + 1) 0)
      (G.replicate nnz 0)
      (G.replicate nnz 0)
      (G.replicate dimension 0)

csrMatrixShrink :: HasCallStack => Int -> CsrMatrix -> CsrMatrix
csrMatrixShrink nnz m@(CsrMatrix offsets columns offDiagElems diagElems)
  | nnz == csrNumberNonZero m = m
  | nnz < csrNumberNonZero m = CsrMatrix offsets (G.take nnz columns) (G.take nnz offDiagElems) diagElems
  | otherwise = error "general resizing not implemented"

csrDim :: CsrMatrix -> Int
csrDim (CsrMatrix _ _ _ diagElems) = G.length diagElems

csrNumberNonZero :: CsrMatrix -> Int
csrNumberNonZero (CsrMatrix _ _ offDiagElems _) = G.length offDiagElems

csrMatrixFromDense ::
  (HasCallStack, G.Vector v (Complex Double), G.Vector v Int) =>
  DenseMatrix v (Complex Double) ->
  CsrMatrix
csrMatrixFromDense m@(DenseMatrix nRows nCols v)
  | nRows == nCols = runST $ do
    let nnz = countOffDiagNonZero m
        diagElems = G.convert $ extractDiagonal m
    offsetsBuffer <- GM.new (nRows + 1)
    columnsBuffer <- GM.new nnz
    offDiagElemsBuffer <- GM.new nnz
    let go2 !i !j !offset
          | i == j = go2 i (j + 1) offset
          | j < nCols =
            let x = indexDenseMatrix m (i, j)
             in if x /= 0
                  then do
                    GM.write columnsBuffer offset (fromIntegral j)
                    GM.write offDiagElemsBuffer offset x
                    go2 i (j + 1) (offset + 1)
                  else go2 i (j + 1) offset
          | otherwise = pure offset
        go1 !i !offset
          | i < nRows = do
            offset' <- go2 i 0 offset
            GM.write offsetsBuffer (i + 1) (fromIntegral offset')
            go1 (i + 1) offset'
          | otherwise = pure offset
    GM.write offsetsBuffer 0 0
    offset <- go1 0 0
    unless (nnz == offset) $ error "this is probably a bug"
    CsrMatrix
      <$> G.unsafeFreeze offsetsBuffer
      <*> G.unsafeFreeze columnsBuffer
      <*> G.unsafeFreeze offDiagElemsBuffer
      <*> pure diagElems
  | otherwise = error "expected a square matrix"

foreign import capi unsafe "csr.h ls_csr_plus"
  ls_csr_plus :: Ptr Ccsr_matrix -> Ptr Ccsr_matrix -> Ptr Ccsr_matrix -> IO ()

foreign import capi unsafe "csr.h ls_csr_minus"
  ls_csr_minus :: Ptr Ccsr_matrix -> Ptr Ccsr_matrix -> Ptr Ccsr_matrix -> IO ()

foreign import capi unsafe "csr.h ls_csr_times"
  ls_csr_times :: Ptr Ccsr_matrix -> Ptr Ccsr_matrix -> Ptr Ccsr_matrix -> IO ()

foreign import capi unsafe "csr.h ls_csr_kron"
  ls_csr_kron :: Ptr Ccsr_matrix -> Ptr Ccsr_matrix -> Ptr Ccsr_matrix -> IO ()

binaryOp :: (Ptr Ccsr_matrix -> Ptr Ccsr_matrix -> Ptr Ccsr_matrix -> IO ()) -> CsrMatrix -> CsrMatrix -> CsrMatrix
binaryOp opKernel a b = System.IO.Unsafe.unsafePerformIO $ do
  let pessimisticNNZ = csrNumberNonZero a + csrNumberNonZero b
  c <- unsafeNewCsrMatrix (csrDim a) pessimisticNNZ
  nnz <-
    withCsrMatrix a $ \aPtr ->
      withCsrMatrix b $ \bPtr ->
        withCsrMatrix c $ \cPtr -> do
          opKernel aPtr bPtr cPtr
          fromIntegral . c_csr_matrix_number_nonzero <$> peek cPtr
  pure $ csrMatrixShrink nnz c

instance Num CsrMatrix where
  (+) = binaryOp ls_csr_plus
  (-) = binaryOp ls_csr_minus
  (*) = binaryOp ls_csr_times
  abs (CsrMatrix offsets columns offDiagElems diagElems) = CsrMatrix offsets columns (G.map abs offDiagElems) (G.map abs diagElems)
  signum (CsrMatrix offsets columns offDiagElems diagElems) = CsrMatrix offsets columns (G.map signum offDiagElems) (G.map signum diagElems)
  fromInteger _ = error "Num instance of CsrMatrix does not implement fromInteger"

csrScale :: Complex Double -> CsrMatrix -> CsrMatrix
csrScale c (CsrMatrix offsets columns offDiagElems diagElems) =
  CsrMatrix offsets columns (G.map (c *) offDiagElems) (G.map (c *) diagElems)

csrKron :: CsrMatrix -> CsrMatrix -> CsrMatrix
csrKron a b = System.IO.Unsafe.unsafePerformIO $ do
  let dimension = csrDim a * csrDim b
      nnz = csrNumberNonZero a * csrNumberNonZero b
  c <- unsafeNewCsrMatrix dimension nnz
  nnz' <-
    withCsrMatrix a $ \aPtr ->
      withCsrMatrix b $ \bPtr ->
        withCsrMatrix c $ \cPtr -> do
          ls_csr_kron aPtr bPtr cPtr
          fromIntegral . c_csr_matrix_number_nonzero <$> peek cPtr
  unless (nnz == nnz') $
    error $ "this is probably a bug: " <> show nnz <> " /= " <> show nnz'
  pure c

csrKronMany :: HasCallStack => [CsrMatrix] -> CsrMatrix
csrKronMany [] = error "expected a non-empty list of matrices"
csrKronMany xs = Data.List.foldl1' csrKron xs
