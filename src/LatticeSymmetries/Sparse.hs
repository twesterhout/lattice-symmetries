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

module LatticeSymmetries.Sparse where

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
import Data.Type.Equality
-- import qualified Data.List
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
import Data.Vector.Generic ((!))
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Generic.Mutable as GM
import qualified Data.Vector.Storable
import qualified Data.Vector.Unboxed
-- import qualified Data.Vector.Storable.Mutable as SM
-- import Data.Vector.Unboxed (Unbox)
-- import qualified Data.Vector.Unboxed as U
-- import Data.Yaml (decodeFileWithWarnings)
-- import Foreign.C.String (CString, peekCString)
-- import Foreign.C.Types (CInt (..), CUInt (..))
-- import Foreign.ForeignPtr
-- import Foreign.ForeignPtr.Unsafe (unsafeForeignPtrToPtr)
-- import Foreign.Marshal.Alloc (alloca, free, malloc, mallocBytes)
-- import Foreign.Marshal.Array (newArray, withArrayLen)
-- import Foreign.Marshal.Utils (new, with)
-- import Foreign.Ptr (FunPtr, Ptr, castPtr)
-- import Foreign.StablePtr
-- import Foreign.Storable (Storable (..))
import qualified GHC.Exts as GHC (IsList (..))
import GHC.TypeLits
import qualified GHC.TypeLits as GHC
import qualified Unsafe.Coerce
-- import qualified GHC.ForeignPtr as GHC (Finalizers (..), ForeignPtr (..), ForeignPtrContents (..))
-- import GHC.Generics
-- import qualified GHC.IORef as GHC (atomicSwapIORef)
-- import GHC.Prim
-- import qualified GHC.Ptr as GHC (Ptr (..))
-- import qualified Language.C.Inline as C
-- import qualified Language.C.Inline.Unsafe as CU
-- import LatticeSymmetries.Context
-- import LatticeSymmetries.IO
-- import LatticeSymmetries.Types
-- import qualified System.IO.Unsafe
-- import qualified System.Mem.Weak
import Prelude hiding (group, product, sort)

-- C.context (C.baseCtx <> C.bsCtx <> C.funCtx <> lsCtx)
-- C.include "<lattice_symmetries/lattice_symmetries.h>"
-- C.include "helpers.h"

type KnownDenseMatrix v r c a = (G.Vector v a, KnownNat r, KnownNat c)

type KnownCOO v i r c a = (G.Vector v (i, i, a), KnownNat r, KnownNat c)

type KnownCSR v i r c a = (G.Vector v i, G.Vector v a, KnownNat r, KnownNat c)

natToInt :: forall n. KnownNat n => Int
natToInt = fromIntegral $ GHC.TypeLits.natVal (Proxy @n)

-- | Dense matrix in row-major order (C layout)
data DenseMatrix v (r :: Nat) (c :: Nat) a = DenseMatrix {dmData :: !(v a)}
  deriving stock (Show, Eq, Generic)

type StorableDenseMatrix r c a = DenseMatrix Data.Vector.Storable.Vector r c a

type UnboxedDenseMatrix r c a = DenseMatrix Data.Vector.Unboxed.Vector r c a

type BoxedDenseMatrix r c a = DenseMatrix Data.Vector.Vector r c a

-- | Get number of rows in the matrix
dmRows :: forall r c a v. KnownDenseMatrix v r c a => DenseMatrix v r c a -> Int
dmRows _ = natToInt @r
{-# INLINE dmRows #-}

-- | Get number of columns in the matrix
dmCols :: forall r c a v. KnownDenseMatrix v r c a => DenseMatrix v r c a -> Int
dmCols _ = natToInt @c
{-# INLINE dmCols #-}

-- | Get matrix shape
dmShape :: forall r c a v. KnownDenseMatrix v r c a => DenseMatrix v r c a -> (Int, Int)
dmShape m = (dmRows m, dmCols m)
{-# INLINE dmShape #-}

instance (KnownDenseMatrix v r c a, Num a) => Num (DenseMatrix v r c a) where
  (+) a b = DenseMatrix $ G.zipWith (+) (dmData a) (dmData b)
  (-) a b = DenseMatrix $ G.zipWith (-) (dmData a) (dmData b)
  (*) a b = DenseMatrix $ G.zipWith (*) (dmData a) (dmData b)
  abs a = DenseMatrix $ G.map abs (dmData a)
  signum _ = error "Num instance for DenseMatrix does not implement signum"
  fromInteger z = DenseMatrix $ G.replicate (natToInt @r * natToInt @c) (fromInteger z)

-- | Sparse matrix in Coordinate format
data COO v i (r :: Nat) (c :: Nat) a = COO {cooData :: !(v (i, i, a))}
  deriving stock (Generic)

deriving instance Show (v (i, i, a)) => Show (COO v i r c a)

deriving instance Eq (v (i, i, a)) => Eq (COO v i r c a)

-- | Get number of rows in the matrix
cooRows :: forall r c a i v. KnownCOO v i r c a => COO v i r c a -> Int
cooRows _ = natToInt @r
{-# INLINE cooRows #-}

-- | Get number of columns in the matrix
cooCols :: forall r c a i v. KnownCOO v i r c a => COO v i r c a -> Int
cooCols _ = natToInt @c
{-# INLINE cooCols #-}

-- | Get matrix shape
cooShape :: forall r c a i v. KnownCOO v i r c a => COO v i r c a -> (Int, Int)
cooShape m = (cooRows m, cooCols m)
{-# INLINE cooShape #-}

data CSR v i (r :: Nat) (c :: Nat) a = CSR
  { csrOffsets :: !(v i),
    csrIndices :: !(v i),
    csrData :: !(v a)
  }
  deriving stock (Generic)

deriving instance (Show (v i), Show (v a)) => Show (CSR v i r c a)

deriving instance (Eq (v i), Eq (v a)) => Eq (CSR v i r c a)

data SomeCSR v i a where
  SomeCSR :: (KnownNat r, KnownNat c) => CSR v i r c a -> SomeCSR v i a

sameShape ::
  forall r₁ c₁ r₂ c₂ a₁ a₂ v₁ v₂ i₁ i₂.
  (KnownCSR v₁ i₁ r₁ c₁ a₁, KnownCSR v₂ i₂ r₂ c₂ a₂) =>
  CSR v₁ i₁ r₁ c₁ a₁ ->
  CSR v₂ i₂ r₂ c₂ a₂ ->
  Maybe ('(r₁, c₁) :~: '(r₂, c₂))
sameShape a b =
  case sameNat (Proxy @r₁) (Proxy @r₂) of
    Just Refl -> case sameNat (Proxy @c₁) (Proxy @c₂) of
      Just Refl -> Just Refl
      Nothing -> Nothing
    Nothing -> Nothing

withSomeCsr ::
  (G.Vector v i, G.Vector v a) =>
  SomeCSR v i a ->
  (forall r c. (KnownNat r, KnownNat c) => CSR v i r c a -> b) ->
  b
withSomeCsr (SomeCSR m) f = f m
{-# INLINE withSomeCsr #-}

-- | Binary search in a vector. Return index, if found.
binarySearch :: (G.Vector v a, Ord a) => v a -> a -> Maybe Int
binarySearch v z = go 0 (G.length v)
  where
    {-# INLINE go #-}
    go !l !u
      | l < u =
        -- NOTE: we assume that the vector is short enought such that @u + l@ does not overflow
        let !m = (u + l) `div` 2
         in case compare (v ! m) z of
              LT -> go (m + 1) u
              EQ -> Just m
              GT -> go l m
      | otherwise = Nothing
{-# INLINE binarySearch #-}

csrIndex :: (KnownCSR v i r c a, Integral i, Num a) => CSR v i r c a -> (Int, Int) -> a
csrIndex matrix (i, j) =
  case binarySearch (G.slice l (u - l) (csrIndices matrix)) (fromIntegral j) of
    Just k -> csrData matrix ! (l + k)
    Nothing -> 0
  where
    !l = fromIntegral $ csrOffsets matrix ! i
    !u = fromIntegral $ csrOffsets matrix ! (i + 1)
{-# INLINE csrIndex #-}

csrNumberNonZero :: KnownCSR v i r c a => CSR v i r c a -> Int
csrNumberNonZero = G.length . csrIndices
{-# INLINE csrNumberNonZero #-}

csrRows :: forall r c a i v. KnownNat r => CSR v i r c a -> Int
csrRows _ = natToInt @r
{-# INLINE csrRows #-}

csrCols :: forall r c a i v. KnownNat c => CSR v i r c a -> Int
csrCols _ = natToInt @c
{-# INLINE csrCols #-}

csrShape :: (KnownNat r, KnownNat c) => CSR v i r c a -> (Int, Int)
csrShape csr = (csrRows csr, csrCols csr)
{-# INLINE csrShape #-}

csrTraverseIndex ::
  forall i₁ i₂ f v r c a.
  (Applicative f, KnownCSR v i₁ r c a, KnownCSR v i₂ r c a) =>
  (i₁ -> f i₂) ->
  CSR v i₁ r c a ->
  f (CSR v i₂ r c a)
csrTraverseIndex f (CSR offsets indices elements) = construct <$> offsets' <*> indices'
  where
    construct o i = CSR o i elements
    offsets' = G.fromListN n <$> traverse f (G.toList offsets)
    indices' = G.fromListN n <$> traverse f (G.toList indices)
    !n = G.length indices
{-# INLINE csrTraverseIndex #-}

-- | Change the underlying integral type. Returns a @Maybe@ because casts may overflow.
csrReIndex ::
  (KnownCSR v i₁ r c a, KnownCSR v i₂ r c a, Bits i₁, Integral i₁, Bits i₂, Integral i₂) =>
  CSR v i₁ r c a ->
  Maybe (CSR v i₂ r c a)
csrReIndex = csrTraverseIndex toIntegralSized
{-# INLINE csrReIndex #-}

-- | Change the underlying vector type.
csrReVector ::
  (KnownCSR v₁ i r c a, KnownCSR v₂ i r c a) =>
  CSR v₁ i r c a ->
  CSR v₂ i r c a
csrReVector (CSR offsets indices elements) =
  CSR (G.convert offsets) (G.convert indices) (G.convert elements)
{-# INLINE csrReVector #-}

cooToBundle :: KnownCOO v i r c a => COO v i r c a -> Bundle v (i, i, a)
cooToBundle = G.stream . cooData

cooFromBundle :: KnownCOO v i r c a => Bundle v (i, i, a) -> COO v i r c a
cooFromBundle = COO . G.unstream . Bundle.reVector

computeShape :: (Integral i, G.Vector v (i, i, a)) => v (i, i, a) -> (Int, Int)
computeShape v
  | G.null v = (0, 0)
  | otherwise = (fromIntegral maxRow + 1, fromIntegral maxColumn + 1)
  where
    (maxRow, maxColumn) = G.foldl' combine (0, 0) v
    combine (!r, !c) (!i, !j, _) = let !r' = max i r; !c' = max j c in (r', c')

data CombineNeighborsHelper a
  = CombineNeighborsFirst
  | CombineNeighborsPrevious !a
  | CombineNeighborsDone

combineNeighborsImpl :: Monad m => (a -> a -> Bool) -> (a -> a -> a) -> Stream m a -> Stream m a
{-# INLINE combineNeighborsImpl #-}
combineNeighborsImpl equal combine (Stream step s₀) = Stream step' (CombineNeighborsFirst, s₀)
  where
    {-# INLINE step' #-}
    step' (CombineNeighborsFirst, s) = do
      r <- step s
      case r of
        Yield a s' -> pure $ Skip (CombineNeighborsPrevious a, s')
        Skip s' -> pure $ Skip (CombineNeighborsFirst, s')
        Done -> pure $ Done
    step' (CombineNeighborsPrevious a, s) = do
      r <- step s
      case r of
        Yield b s' ->
          if equal a b
            then pure $ Skip (CombineNeighborsPrevious (combine a b), s')
            else pure $ Yield a (CombineNeighborsPrevious b, s')
        Skip s' -> pure $ Skip (CombineNeighborsPrevious a, s')
        Done -> pure $ Yield a (CombineNeighborsDone, s)
    step' (CombineNeighborsDone, _) = pure $ Done

combineNeighbors :: G.Vector v a => (a -> a -> Bool) -> (a -> a -> a) -> v a -> v a
combineNeighbors equal combine =
  G.unstream . Bundle.inplace (combineNeighborsImpl equal combine) toMax . G.stream

cooNormalize :: (KnownCOO v i r c a, Integral i, Num a) => COO v i r c a -> COO v i r c a
cooNormalize (COO v) = COO v'
  where
    comparison (i₁, j₁, _) (i₂, j₂, _) = compare (i₁, j₁) (i₂, j₂)
    group =
      combineNeighbors
        (\(i₁, j₁, _) (i₂, j₂, _) -> i₁ == i₂ && j₁ == j₂)
        (\(i, j, a) (_, _, b) -> (i, j, a + b))
    sort x = runST $ do
      buffer <- G.thaw x
      Intro.sortBy comparison buffer
      G.unsafeFreeze buffer
    v' = group $ sort v

-- preprocessCoo :: Num a => [(Int, Int, a)] -> [(Int, Int, a)]
-- preprocessCoo =
--   fmap (Data.List.foldl1' (\(!i, !j, !x) (_, _, y) -> (i, j, x + y)))
--     . Data.List.groupBy (\a b -> key a == key b)
--     . sortOn key
--   where
--     key (i, j, _) = (i, j)

unsafeCooToCsr ::
  (KnownCOO v1 i r c a, KnownCSR v2 i r c a, G.Vector v1 a, G.Vector v1 i, Integral i) =>
  COO v1 i r c a ->
  CSR v2 i r c a
unsafeCooToCsr coo@(COO coordinates) = CSR offsets indices elements
  where
    indices = G.convert $ G.map (\(_, j, _) -> j) coordinates
    elements = G.convert $ G.map (\(_, _, x) -> x) coordinates
    n = cooRows coo
    offsets = runST $ do
      rs <- GM.replicate (n + 1) 0
      G.forM_ coordinates $ \(!i, _, _) ->
        GM.modify rs (+ 1) (fromIntegral i + 1)
      _ <- loopM 0 (< n) (+ 1) $ \ !i -> do
        r <- GM.read rs i
        GM.modify rs (+ r) (i + 1)
      G.unsafeFreeze rs

cooToCsr ::
  forall v1 v2 i r c a.
  (KnownCOO v1 i r c a, KnownCSR v2 i r c a, G.Vector v1 a, G.Vector v1 i, Integral i, Num a) =>
  COO v1 i r c a ->
  CSR v2 i r c a
cooToCsr = unsafeCooToCsr . cooNormalize

csrOffsetsToIndices :: (G.Vector v i, Integral i) => v i -> v i
csrOffsetsToIndices offsets = runST $ do
  let nRows = G.length offsets - 1
      nnz = fromIntegral $ G.last offsets
  indices <- GM.new nnz
  _ <- loopM 0 (< nRows) (+ 1) $ \i ->
    let !b = fromIntegral $ offsets ! i
        !e = fromIntegral $ offsets ! (i + 1)
     in GM.set (GM.slice b (e - b) indices) (fromIntegral i)
  G.unsafeFreeze indices

csrToCoo ::
  forall v1 v2 i r c a.
  (KnownCSR v1 i r c a, KnownCOO v2 i r c a, G.Vector v2 a, G.Vector v2 i, Integral i) =>
  CSR v1 i r c a ->
  COO v2 i r c a
csrToCoo csr = COO v
  where
    v =
      G.zip3
        (G.convert $ csrOffsetsToIndices (csrOffsets csr))
        (G.convert $ csrIndices csr)
        (G.convert $ csrData csr)

denseToCoo ::
  forall v1 v2 i r c a.
  (KnownDenseMatrix v1 r c a, KnownCOO v2 i r c a, Eq a, Num a, Integral i) =>
  DenseMatrix v1 r c a ->
  COO v2 i r c a
denseToCoo dense = COO v
  where
    numberRows = dmRows dense
    numberCols = dmCols dense
    v =
      G.unstream
        . Bundle.filter (\(_, _, x) -> x /= 0)
        . Bundle.map (\(i, j) -> (fromIntegral i, fromIntegral j, indexDenseMatrix dense i j))
        $ cartesian
          (Bundle.enumFromStepN 0 1 numberRows)
          (Bundle.enumFromStepN 0 1 numberCols)

denseToCsr ::
  forall v i r c a.
  (KnownDenseMatrix v r c a, KnownCSR v i r c a, Eq a, Num a, Integral i) =>
  DenseMatrix v r c a ->
  CSR v i r c a
denseToCsr = cooToCsr . denseToCoo @v @Data.Vector.Vector

cooToDense ::
  forall v1 v2 i r c a.
  (HasCallStack, KnownCOO v1 i r c a, KnownDenseMatrix v2 r c a, Integral i) =>
  COO v1 i r c a ->
  DenseMatrix v2 r c a
cooToDense coo
  | (n, m) == (natToInt @r, natToInt @c) = runST $ do
    elements <- GM.new (n * m)
    G.forM_ (cooData coo) $ \(i, j, x) ->
      GM.write elements (fromIntegral i * m + fromIntegral j) x
    DenseMatrix <$> G.unsafeFreeze elements
  | otherwise = error "incompatible shape"
  where
    (n, m) = cooShape coo

csrToDense ::
  forall v i r c a.
  (KnownDenseMatrix v r c a, KnownCSR v i r c a, G.Vector v (i, i, a), Integral i) =>
  CSR v i r c a ->
  DenseMatrix v r c a
csrToDense = cooToDense . csrToCoo @v @Data.Vector.Vector

-- csrFromList :: KnownCSR v i r c a => [(i, i, a)] -> CSR v i r c a
-- csrFromList = undefined

someCsrShape :: SomeCSR v i a -> (Int, Int)
someCsrShape (SomeCSR m) = csrShape m

cartesian :: Bundle v a -> Bundle v b -> Bundle v (a, b)
cartesian a b =
  Bundle.fromStream
    (cartesianImpl (Bundle.elements a) (Bundle.elements b))
    (Bundle.sSize a `product` Bundle.sSize b)
  where
    product Unknown _ = Unknown
    product _ Unknown = Unknown
    product (Exact x) (Exact y) = Exact (x * y)
    product (Exact x) (Max y) = Max (x * y)
    product (Max x) (Exact y) = Max (x * y)
    product (Max x) (Max y) = Max (x * y)

cartesianImpl :: Monad m => Stream m a -> Stream m b -> Stream m (a, b)
cartesianImpl (Stream step1 s1₀) (Stream step2 s2₀) = Stream step' (Nothing, s1₀, s2₀)
  where
    step' (Just a, s1, s2) = do
      r₂ <- step2 s2
      case r₂ of
        Yield b s2' -> pure $ Yield (a, b) (Just a, s1, s2')
        Skip s2' -> pure $ Skip (Just a, s1, s2')
        -- NOTE: It's important to pass s2₀ here!
        Done -> pure $ Skip (Nothing, s1, s2₀)
    step' (Nothing, s1, s2) = do
      r₁ <- step1 s1
      case r₁ of
        Yield a s1' -> pure $ Skip (Just a, s1', s2)
        Skip s1' -> pure $ Skip (Nothing, s1', s2)
        Done -> pure $ Done

-- mergeImpl :: Stream m a -> Stream m a -> Stream m a
-- mergeImpl (Stream step1 s1₀) (Stream step2 s2₀) = Stream step' (s1₀, s2₀)
--   where
--     step' (s1, s2) = do

cooKron ::
  (KnownCOO v i r1 c1 a, KnownCOO v i r2 c2 a, Num i, Num a) =>
  COO v i r1 c1 a ->
  COO v i r2 c2 a ->
  COO v i (r1 GHC.* r2) (c1 GHC.* c2) a
cooKron a b = COO v
  where
    combine ((i₁, j₁, x₁), (i₂, j₂, x₂)) =
      (i₁ * fromIntegral n₁ + i₂, j₁ * fromIntegral m₁ + j₂, x₁ * x₂)
    v = G.unstream . Bundle.map combine $ cartesian (cooToBundle a) (cooToBundle b)
    (n₁, m₁) = cooShape a

-- cooPlus :: (Integral i, Unbox i, Num a, Unbox a) => COO v i r c a -> COO v i r c a -> COO v i r c a
-- cooPlus a b = undefined

withSomeNat :: forall r. HasCallStack => Int -> (forall n. KnownNat n => Proxy n -> r) -> r
withSomeNat n f = case GHC.TypeLits.someNatVal (fromIntegral n) of
  Just (SomeNat (Proxy :: Proxy n)) -> f (Proxy @n)
  Nothing -> error "negative natural number"

csrKron ::
  forall v i r1 c1 r2 c2 a.
  ( KnownCSR v i r1 c1 a,
    KnownCSR v i r2 c2 a,
    KnownNat (r1 GHC.* r2),
    KnownNat (c1 GHC.* c2),
    Integral i,
    Num a
  ) =>
  CSR v i r1 c1 a ->
  CSR v i r2 c2 a ->
  CSR v i (r1 GHC.* r2) (c1 GHC.* c2) a
csrKron a b = (cooToCsr @B.Vector) $ cooKron (csrToCoo a) (csrToCoo b)

loopM :: Monad m => i -> (i -> Bool) -> (i -> i) -> (i -> m ()) -> m i
loopM i₀ cond inc action = go i₀
  where
    go !i
      | cond i = do () <- action i; go (inc i)
      | otherwise = return i
{-# INLINE loopM #-}

iFoldM :: Monad m => i -> (i -> Bool) -> (i -> i) -> a -> (a -> i -> m a) -> m a
iFoldM i₀ cond inc x₀ action = go x₀ i₀
  where
    go !x !i
      | cond i = do !x' <- action x i; go x' (inc i)
      | otherwise = pure x
{-# INLINE iFoldM #-}

instance (KnownCSR v i r c a, Integral i, Eq a, Num a) => Num (CSR v i r c a) where
  (+) = csrBinaryOp (+)
  (-) = csrBinaryOp (-)
  (*) = csrBinaryOp (*)
  abs m = m {csrData = G.map abs (csrData m)}
  signum = error "Num instance for CSR does not implement signum"
  fromInteger = error "Num instance for CSR does not implement fromInteger"

csrBinaryOp ::
  (KnownCSR v i r c a, KnownCSR v i r c b, Integral i, Eq b, Num a, Num b) =>
  (a -> a -> b) ->
  CSR v i r c a ->
  CSR v i r c a ->
  CSR v i r c b
csrBinaryOp op a b = runST $ do
  outOffsets <- GM.new (csrRows a + 1)
  let estimatedNumberNonZero = csrNumberNonZero a + csrNumberNonZero b
  outIndices <- GM.new estimatedNumberNonZero
  outData <- GM.new estimatedNumberNonZero
  let consumeBoth !kA !kB !nnz =
        let !r = op (csrData a ! kA) (csrData b ! kB)
         in if r /= 0
              then do
                GM.write outIndices nnz (csrIndices a ! kA)
                GM.write outData nnz r
                pure (nnz + 1)
              else pure nnz
      consumeLeft !kA !nnz =
        let !r = op (csrData a ! kA) 0
         in if r /= 0
              then do
                GM.write outIndices nnz (csrIndices a ! kA)
                GM.write outData nnz r
                pure (nnz + 1)
              else pure nnz
      consumeRight !kB !nnz =
        let !r = op 0 (csrData b ! kB)
         in if r /= 0
              then do
                GM.write outIndices nnz (csrIndices b ! kB)
                GM.write outData nnz r
                pure (nnz + 1)
              else pure nnz
  GM.write outOffsets 0 0
  nnz <- iFoldM 0 (< csrRows a) (+ 1) 0 $ \_nnz i -> do
    let !beginA = fromIntegral $ csrOffsets a ! i
        !beginB = fromIntegral $ csrOffsets b ! i
        !endA = fromIntegral $ csrOffsets a ! (i + 1)
        !endB = fromIntegral $ csrOffsets b ! (i + 1)
        loop !kA !kB !nnz
          | kA < endA && kB < endB =
            let !jA = csrIndices a ! kA
                !jB = csrIndices b ! kB
             in case compare jA jB of
                  EQ -> consumeBoth kA kB nnz >>= loop (kA + 1) (kB + 1)
                  LT -> consumeLeft kA nnz >>= loop (kA + 1) kB
                  GT -> consumeRight kB nnz >>= loop kA (kB + 1)
          | kA < endA = tailLeft kA nnz
          | otherwise = tailRight kB nnz
        tailLeft !kA !nnz
          | kA < endA = consumeLeft kA nnz >>= tailLeft (kA + 1)
          | otherwise = pure nnz
        tailRight !kB !nnz
          | kB < endB = consumeRight kB nnz >>= tailRight (kB + 1)
          | otherwise = pure nnz
    nnz <- loop beginA beginB _nnz
    GM.write outOffsets (i + 1) (fromIntegral nnz)
    pure nnz
  CSR
    <$> G.unsafeFreeze outOffsets
    <*> G.unsafeFreeze (GM.slice 0 nnz outIndices)
    <*> G.unsafeFreeze (GM.slice 0 nnz outData)

-- template <class I, class T, class T2, class binary_op>
-- void csr_binop_csr_canonical(const I n_row, const I n_col,
--                              const I Ap[], const I Aj[], const T Ax[],
--                              const I Bp[], const I Bj[], const T Bx[],
--                                    I Cp[],       I Cj[],       T2 Cx[],
--                              const binary_op& op)
-- {
--     //Method that works for canonical CSR matrices
--
--     Cp[0] = 0;
--     I nnz = 0;
--
--     for(I i = 0; i < n_row; i++){
--         I A_pos = Ap[i];
--         I B_pos = Bp[i];
--         I A_end = Ap[i+1];
--         I B_end = Bp[i+1];
--
--         //while not finished with either row
--         while(A_pos < A_end && B_pos < B_end){
--             I A_j = Aj[A_pos];
--             I B_j = Bj[B_pos];
--
--             if(A_j == B_j){
--                 T result = op(Ax[A_pos],Bx[B_pos]);
--                 if(result != 0){
--                     Cj[nnz] = A_j;
--                     Cx[nnz] = result;
--                     nnz++;
--                 }
--                 A_pos++;
--                 B_pos++;
--             } else if (A_j < B_j) {
--                 T result = op(Ax[A_pos],0);
--                 if (result != 0){
--                     Cj[nnz] = A_j;
--                     Cx[nnz] = result;
--                     nnz++;
--                 }
--                 A_pos++;
--             } else {
--                 //B_j < A_j
--                 T result = op(0,Bx[B_pos]);
--                 if (result != 0){
--                     Cj[nnz] = B_j;
--                     Cx[nnz] = result;
--                     nnz++;
--                 }
--                 B_pos++;
--             }
--         }
--
--         //tail
--         while(A_pos < A_end){
--             T result = op(Ax[A_pos],0);
--             if (result != 0){
--                 Cj[nnz] = Aj[A_pos];
--                 Cx[nnz] = result;
--                 nnz++;
--             }
--             A_pos++;
--         }
--         while(B_pos < B_end){
--             T result = op(0,Bx[B_pos]);
--             if (result != 0){
--                 Cj[nnz] = Bj[B_pos];
--                 Cx[nnz] = result;
--                 nnz++;
--             }
--             B_pos++;
--         }
--
--         Cp[i+1] = nnz;
--     }
-- }

maxNumNonZeroAfterMatMul ::
  forall r k c a i v.
  (KnownCSR v i r k a, KnownCSR v i k c a, G.Vector v Int, Integral i) =>
  CSR v i r k a ->
  CSR v i k c a ->
  Int
maxNumNonZeroAfterMatMul a b = runST $ do
  let nRows = csrRows a
      nCols = csrCols b
  mask <- G.unsafeThaw (G.replicate nCols (-1) :: v Int)
  iFoldM 0 (< nRows) (+ 1) 0 $ \nnz i ->
    let jjBegin = fromIntegral $ csrOffsets a ! i
        jjEnd = fromIntegral $ csrOffsets a ! (i + 1)
     in iFoldM jjBegin (< jjEnd) (+ 1) nnz $ \nnzRow jj ->
          let j = fromIntegral $ csrIndices a ! jj
              kkBegin = fromIntegral $ csrOffsets b ! j
              kkEnd = fromIntegral $ csrOffsets b ! (j + 1)
           in iFoldM kkBegin (< kkEnd) (+ 1) nnzRow $ \acc kk -> do
                let k = fromIntegral $ csrIndices b ! kk
                m <- GM.read mask k
                if m /= i
                  then do GM.write mask k i; pure (acc + 1)
                  else pure acc

sortManyByKey ::
  (Ord c, PrimMonad m, G.Vector v a, G.Vector v b) =>
  (a -> b -> c) ->
  G.Mutable v (PrimState m) a ->
  G.Mutable v (PrimState m) b ->
  m ()
sortManyByKey key a b = do
  mbuffer <-
    (G.unsafeThaw =<<) $
      Data.Vector.zip
        <$> (G.convert <$> G.unsafeFreeze a)
        <*> (G.convert <$> G.unsafeFreeze b)
  Intro.sortBy (comparing (uncurry key)) mbuffer
  buffer <- G.unsafeFreeze mbuffer
  G.copy a $ G.convert (G.map fst buffer)
  G.copy b $ G.convert (G.map snd buffer)

csrMatMul ::
  forall r k c a i v.
  (KnownCSR v i r k a, KnownCSR v i k c a, G.Vector v Int, Integral i, Num a) =>
  CSR v i r k a ->
  CSR v i k c a ->
  CSR v i r c a
csrMatMul a b =
  -- cooToCsr . csrToCoo $
  runST $ do
    let !bufferSize = maxNumNonZeroAfterMatMul a b
        !nRows = csrRows a
        !nCols = csrCols b
    outOffsets <- GM.new (nRows + 1)
    outIndices <- GM.new bufferSize
    outData <- GM.new bufferSize
    next <- G.unsafeThaw (G.replicate nCols (-1) :: v Int)
    sums <- G.unsafeThaw (G.replicate nCols 0 :: v a)
    nnz <- iFoldM 0 (< nRows) (+ 1) 0 $ \ !nnz !i -> do
      let !jjBegin = fromIntegral $ csrOffsets a ! i
          !jjEnd = fromIntegral $ csrOffsets a ! (i + 1)
      (head, length) <- iFoldM jjBegin (< jjEnd) (+ 1) (-2 :: Int, 0 :: Int) $
        \(!_head, !_length) !jj ->
          let !j = fromIntegral $ csrIndices a ! jj
              !v = csrData a ! jj
              !kkBegin = fromIntegral $ csrOffsets b ! j
              !kkEnd = fromIntegral $ csrOffsets b ! (j + 1)
           in iFoldM kkBegin (< kkEnd) (+ 1) (_head, _length) $ \(!head, !length) !kk -> do
                let !k = fromIntegral $ csrIndices b ! kk
                GM.modify sums (+ v * (csrData b ! kk)) k
                _next <- GM.read next k
                if _next == -1
                  then do GM.write next k head; pure (k, length + 1)
                  else pure (head, length)
      (_, nnz') <- iFoldM 0 (< length) (+ 1) (head, nnz) $ \(!head, !nnz) _ -> do
        GM.write outIndices nnz (fromIntegral head)
        GM.write outData nnz =<< GM.read sums head
        head' <- GM.read next head
        GM.write next head (-1)
        GM.write sums head 0
        pure (head', nnz + 1)
      GM.write outOffsets (i + 1) (fromIntegral nnz')
      pure nnz'
    outOffsets' <- G.unsafeFreeze outOffsets
    loopM 0 (< nRows) (+ 1) $ \ !i -> do
      let !jjBegin = fromIntegral $ outOffsets' ! i
          !jjEnd = fromIntegral $ outOffsets' ! (i + 1)
          key !j _ = j
      sortManyByKey
        (\j _ -> j)
        (GM.slice jjBegin (jjEnd - jjBegin) outIndices)
        (GM.slice jjBegin (jjEnd - jjBegin) outData)
    CSR
      <$> pure outOffsets'
      <*> G.unsafeFreeze (GM.slice 0 nnz outIndices)
      <*> G.unsafeFreeze (GM.slice 0 nnz outData)

-- void ls_csr_matrix_from_dense(unsigned const dimension,
--                               _Complex double const *const dense,
--                               unsigned *offsets, unsigned *columns,
--                               _Complex double *off_diag_elements,
--                               _Complex double *diag_elements) {
--   offsets[0] = 0;
--   for (unsigned i = 0; i < dimension; ++i) {
--     unsigned nonzero_in_row = 0;
--     for (unsigned j = 0; j < dimension; ++j) {
--       _Complex double const element = dense[i * dimension + j];
--       if (i == j) {
--         diag_elements[i] = element;
--       } else if (element != 0) {
--         *columns = j;
--         *off_diag_elements = element;
--         ++nonzero_in_row;
--         ++columns;
--         ++off_diag_elements;
--       }
--     }
--     offsets[i + 1] = offsets[i] + nonzero_in_row;
--   }
-- }

isDenseMatrixSquare :: KnownDenseMatrix v r c a => DenseMatrix v r c a -> Bool
isDenseMatrixSquare m = dmRows m == dmCols m

isDenseMatrixEmpty :: KnownDenseMatrix v r c a => DenseMatrix v r c a -> Bool
isDenseMatrixEmpty m = dmRows m * dmCols m == 0

indexDenseMatrix :: KnownDenseMatrix v r c a => DenseMatrix v r c a -> Int -> Int -> a
indexDenseMatrix m i j = dmData m ! (c * i + j)
  where
    c = dmCols m

isDenseMatrixHermitian :: (KnownDenseMatrix v r c (Complex a), Num a, Eq a) => DenseMatrix v r c (Complex a) -> Bool
isDenseMatrixHermitian matrix = isDenseMatrixSquare matrix && go 0 0
  where
    (r, c) = dmShape matrix
    -- Iterate over upper triangle (including the diagonal) of the matrix
    go :: Int -> Int -> Bool
    go !i !j
      | j == c = let !i' = i + 1 in (i' >= r) || go i' i'
      | otherwise =
        let !p = indexDenseMatrix matrix i j == conjugate (indexDenseMatrix matrix j i)
         in p && go i (j + 1)

denseMatrixFromList ::
  forall r c a v.
  (G.Vector v a, KnownNat r, KnownNat c) =>
  [[a]] ->
  Either Text (DenseMatrix v r c a)
denseMatrixFromList rs
  | length rs /= natToInt @r = Left $ "expected " <> show (natToInt @r) <> " rows"
  | any ((/= natToInt @c) . length) rs = Left $ "expected " <> show (natToInt @c) <> " columns"
  | otherwise = Right . DenseMatrix . G.fromList . mconcat $ rs

instance (G.Vector v a, KnownNat r, KnownNat c) => GHC.IsList (DenseMatrix v r c a) where
  type Item (DenseMatrix v r c a) = [a]
  fromList rows = case denseMatrixFromList rows of
    Right m -> m
    Left msg -> error msg

cooFromList ::
  forall r c a i v.
  (KnownCOO v i r c a, Integral i) =>
  [(i, i, a)] ->
  Either Text (COO v i r c a)
cooFromList coordinates
  | all valid coordinates = Right . COO . G.fromList $ coordinates
  | otherwise = Left "index out of bounds"
  where
    valid (i, j, _) = 0 <= i && i < r && 0 <= j && j < c
    r = fromIntegral $ natToInt @r
    c = fromIntegral $ natToInt @c

instance (KnownCOO v i r c a, Integral i) => GHC.IsList (COO v i r c a) where
  type Item (COO v i r c a) = (i, i, a)
  fromList coordinates = case cooFromList coordinates of
    Right m -> m
    Left msg -> error msg

instance (KnownCSR v i r c a, Integral i, Num a) => GHC.IsList (CSR v i r c a) where
  type Item (CSR v i r c a) = (i, i, a)
  fromList coordinates = cooToCsr @Data.Vector.Vector @v $ fromList coordinates
  toList = G.toList . cooData . csrToCoo @v @Data.Vector.Vector

instance (G.Vector v i, G.Vector v a, Integral i, Num a) => GHC.IsList (SomeCSR v i a) where
  type Item (SomeCSR v i a) = (i, i, a)
  fromList coordinates =
    withSomeNat n $ \(Proxy :: Proxy r) ->
      withSomeNat m $ \(Proxy :: Proxy c) ->
        SomeCSR $ fromList @(CSR v i r c a) coordinates
    where
      (n, m) = computeShape (Data.Vector.fromList coordinates)

-- data {-# CTYPE "ls_bit_index" #-} BitIndex = BitIndex !Word8 !Word8
--   deriving stock (Show, Eq, Generic)

--  deriving anyclass (Binary)

-- instance Binary CUInt where
--   put (CUInt x) = Data.Binary.put x
--   get = CUInt <$> Data.Binary.get

-- data SparseSquareMatrix = SparseSquareMatrix
--   { ssmDimension :: {-# UNPACK #-} !Int,
--     ssmOffsets :: {-# UNPACK #-} !(S.Vector CUInt),
--     ssmColumns :: {-# UNPACK #-} !(S.Vector CUInt),
--     ssmOffDiagElements :: {-# UNPACK #-} !(S.Vector (Complex Double)),
--     ssmDiagElements :: {-# UNPACK #-} !(S.Vector (Complex Double))
--   }
--   deriving stock (Show, Eq, Generic)

--   deriving anyclass (Binary)

-- isSparseMatrixHermitian :: SparseSquareMatrix -> Bool
-- isSparseMatrixHermitian = isDenseMatrixHermitian . sparseToDense

-- withCsparse_matrix :: SparseSquareMatrix -> (Csparse_matrix -> IO a) -> IO a
-- withCsparse_matrix matrix action =
--   S.unsafeWith (scsrOffsets matrix) $ \offsetsPtr ->
--     S.unsafeWith (ssmColumns matrix) $ \columnsPtr ->
--       S.unsafeWith (ssmOffDiagElements matrix) $ \offDiagElementsPtr ->
--         S.unsafeWith (ssmDiagElements matrix) $ \diagElementsPtr ->
--           action $
--             Csparse_matrix
--               (fromIntegral . ssmDimension $ matrix)
--               (fromIntegral . S.length . ssmColumns $ matrix)
--               offsetsPtr
--               columnsPtr
--               offDiagElementsPtr
--               diagElementsPtr

-- typedef struct ls_csr_matrix {
--     unsigned         dimension;
--     unsigned         number_nonzero;
--     unsigned*        offsets;
--     unsigned*        columns;
--     _Complex double* off_diag_elements;
--     _Complex double* diag_elements;
-- } ls_csr_matrix;
-- data {-# CTYPE "helpers.h" "ls_csr_matrix" #-} Csparse_matrix
--   = Csparse_matrix
--       {-# UNPACK #-} !CUInt
--       {-# UNPACK #-} !CUInt
--       {-# UNPACK #-} !(Ptr CUInt)
--       {-# UNPACK #-} !(Ptr CUInt)
--       {-# UNPACK #-} !(Ptr (Complex Double))
--       {-# UNPACK #-} !(Ptr (Complex Double))

-- instance Storable Csparse_matrix where
--   sizeOf _ = 40
--   alignment _ = 8
--   peek p =
--     Csparse_matrix
--       <$> peekByteOff p 0
--       <*> peekByteOff p 4
--       <*> peekByteOff p 8
--       <*> peekByteOff p 16
--       <*> peekByteOff p 24
--       <*> peekByteOff p 32
--   poke p (Csparse_matrix dimension number_nonzero offsets columns off_diag_elements diag_elements) = do
--     pokeByteOff p 0 dimension
--     pokeByteOff p 4 number_nonzero
--     pokeByteOff p 8 offsets
--     pokeByteOff p 16 columns
--     pokeByteOff p 24 off_diag_elements
--     pokeByteOff p 32 diag_elements

-- trueCsparse_matrixSizeOf :: Int
-- trueCsparse_matrixSizeOf = fromIntegral [CU.pure| unsigned int { sizeof(ls_csr_matrix) } |]

-- trueCsparse_matrixAlignment :: Int
-- trueCsparse_matrixAlignment = fromIntegral [CU.pure| unsigned int { __alignof__(ls_csr_matrix) } |]

-- typedef struct ls_bit_index {
--     uint8_t word;
--     uint8_t bit;
-- } ls_bit_index;
-- data {-# CTYPE "helpers.h" "ls_bit_index" #-} Cbit_index
--   = Cbit_index {-# UNPACK #-} !Word8 {-# UNPACK #-} !Word8
--   deriving (Show, Eq, Generic)

--  deriving anyclass (Binary)

-- instance Storable Cbit_index where
--   sizeOf _ = 2
--   alignment _ = 1
--   peek p = Cbit_index <$> peekByteOff p 0 <*> peekByteOff p 1
--   poke p (Cbit_index word bit) = pokeByteOff p 0 word >> pokeByteOff p 1 bit

-- trueCbit_indexSizeOf :: Int
-- trueCbit_indexSizeOf = fromIntegral [CU.pure| unsigned int { sizeof(ls_bit_index) } |]

-- trueCbit_indexAlignment :: Int
-- trueCbit_indexAlignment = fromIntegral [CU.pure| unsigned int { __alignof__(ls_bit_index) } |]

-- typedef unsigned (*ls_term_gather_fn)(uint64_t const* /*source*/, ls_bit_index const* /*tuple*/);
-- typedef void (*ls_term_scatter_fn)(unsigned, ls_bit_index const* /*tuple*/,
--                                    uint64_t* /*destination*/);
-- typedef struct ls_term {
--     ls_csr_matrix      matrix;
--     unsigned           number_tuples;
--     unsigned           tuple_size;
--     ls_bit_index*      tuples;
--     ls_term_gather_fn  gather_fn;
--     ls_term_scatter_fn scatter_fn;
-- } ls_term;
-- data {-# CTYPE "helpers.h" "ls_term" #-} Cterm
--   = Cterm
--       {-# UNPACK #-} !Csparse_matrix
--       {-# UNPACK #-} !CUInt
--       {-# UNPACK #-} !CUInt
--       {-# UNPACK #-} !(Ptr Cbit_index)
--       {-# UNPACK #-} !(FunPtr (Ptr Word64 -> Ptr Cbit_index -> IO CUInt))
--       {-# UNPACK #-} !(FunPtr (CUInt -> Ptr Cbit_index -> Ptr Word64 -> IO ()))

-- instance Storable Cterm where
--   sizeOf _ = 72
--   alignment _ = 8
--   peek p =
--     Cterm
--       <$> peekByteOff p 0
--       <*> peekByteOff p 40
--       <*> peekByteOff p 44
--       <*> peekByteOff p 48
--       <*> peekByteOff p 56
--       <*> peekByteOff p 64
--   poke p (Cterm matrix number_tuples tuple_size tuples gather_fn scatter_fn) = do
--     pokeByteOff p 0 matrix
--     pokeByteOff p 40 number_tuples
--     pokeByteOff p 44 tuple_size
--     pokeByteOff p 48 tuples
--     pokeByteOff p 56 gather_fn
--     pokeByteOff p 64 scatter_fn

-- trueCtermSizeOf :: Int
-- trueCtermSizeOf = fromIntegral [CU.pure| unsigned int { sizeof(ls_term) } |]

-- trueCtermAlignment :: Int
-- trueCtermAlignment = fromIntegral [CU.pure| unsigned int { __alignof__(ls_term) } |]

-- data SitesList = SitesList
--   { slNumberTuples :: !Int,
--     slTupleSize :: !Int,
--     slData :: !(S.Vector Cbit_index)
--   }
--   deriving stock (Show, Eq, Generic)

--  deriving anyclass (Binary)

-- sitesListFromList :: [[Int]] -> Either Text SitesList
-- sitesListFromList rows = do
--   (DenseMatrix (numberTuples, tupleSize) v) <- denseMatrixFromList rows
--   when (S.any (< 0) v) $
--     Left "sites list cannot have negative elements"
--   when (S.any (> fromIntegral (maxBound :: Word16)) v) $
--     Left "sites list cannot such large elements"
--   Right $ SitesList numberTuples tupleSize (S.map toBitIndex v)
--   where
--     toBitIndex x = Cbit_index (fromIntegral (x `div` 64)) (fromIntegral (x `mod` 64))

-- newtype SitesList = SitesList [[Int]]

-- data OperatorTerm = OperatorTerm {otMatrix :: !SparseSquareMatrix, otSites :: !SitesList}
--   deriving stock (Show, Eq, Generic)

--  deriving anyclass (Binary)

-- data {-# CTYPE "ls_hs_operator_term" #-} OperatorTermWrapper
--   = OperatorTermWrapper
--       {-# UNPACK #-} !(Ptr Cterm)
--       {-# UNPACK #-} !(StablePtr OperatorTerm)

-- instance Storable OperatorTermWrapper where
--   {-# INLINE sizeOf #-}
--   sizeOf _ = 16
--   {-# INLINE alignment #-}
--   alignment _ = 8
--   {-# INLINE peek #-}
--   peek p = OperatorTermWrapper <$> peekByteOff p 0 <*> peekByteOff p 8
--   {-# INLINE poke #-}
--   poke p (OperatorTermWrapper term stable) =
--     pokeByteOff p 0 term >> pokeByteOff p 8 stable

-- vectorFromPtr :: Storable a => Int -> Ptr a -> IO (S.Vector a)
-- vectorFromPtr n p = S.freeze =<< SM.unsafeFromForeignPtr0 <$> newForeignPtr_ p <*> pure n

-- denseMatrixFromPtr :: Int -> Int -> Ptr a -> IO (DenseMatrix v r c a)
-- denseMatrixFromPtr = undefined

-- ls_hs_create_operator_term_from_dense ::
--   CUInt ->
--   Ptr (Complex Double) ->
--   CUInt ->
--   CUInt ->
--   Ptr Word16 ->
--   IO (Ptr OperatorTermWrapper)
-- ls_hs_create_operator_term_from_dense dimension matrixData numberTuples tupleSize tuplesData = do
--   let dimension' = fromIntegral dimension
--       numberTuples' = fromIntegral numberTuples
--       tupleSize' = fromIntegral tupleSize
--   matrixContents <- vectorFromPtr (dimension' * dimension') matrixData
--   let matrix = case denseToSparse $ DenseMatrix (dimension', dimension') matrixContents of
--         Right m -> m
--         Left e -> error e
--       toBitIndex x = Cbit_index (fromIntegral (x `div` 64)) (fromIntegral (x `mod` 64))
--   sitesContents <- vectorFromPtr (numberTuples' * tupleSize') tuplesData
--   let sites = SitesList numberTuples' tupleSize' (S.map toBitIndex sitesContents)
--   when (2 ^ (slTupleSize sites) /= ssmDimension matrix) $
--     error $ "wrong matrix dimension"
--   let term = OperatorTerm matrix sites
--   wrapper <- OperatorTermWrapper <$> allocateCterm term <*> newStablePtr term
--   new wrapper

-- allocateCterm :: OperatorTerm -> IO (Ptr Cterm)
-- allocateCterm term = do
--   p <- malloc
--   withCterm term (poke p)
--   return p

-- deallocateCterm :: Ptr Cterm -> IO ()
-- deallocateCterm p = free p

-- allocateCterms :: NonEmpty OperatorTerm -> IO (Ptr Cterm)
-- allocateCterms terms = do
--   let !count = NonEmpty.length terms
--       !elemSize = let x = x in sizeOf (x :: Cterm)
--   p <- mallocBytes $ elemSize * count
--   forM_ (NonEmpty.zip terms (fromList [0 ..])) $ \(t, i) -> withCterm t (pokeElemOff p i)
--   return p

-- deallocateCterms :: Ptr Cterm -> IO ()
-- deallocateCterms p = free p

-- withCterm :: OperatorTerm -> (Cterm -> IO a) -> IO a
-- withCterm (OperatorTerm matrix sites) action =
--   withCsparse_matrix matrix $ \matrix' ->
--     S.unsafeWith (slData sites) $ \tuplesPtr ->
--       action $
--         Cterm
--           matrix'
--           (fromIntegral . slNumberTuples $ sites)
--           (fromIntegral . slTupleSize $ sites)
--           tuplesPtr
--           gatherPtr
--           scatterPtr
--   where
--     (gatherPtr, scatterPtr) = case slTupleSize sites of
--       1 -> (ls_internal_term_gather_1, ls_internal_term_scatter_1)
--       2 -> (ls_internal_term_gather_2, ls_internal_term_scatter_2)
--       3 -> (ls_internal_term_gather_3, ls_internal_term_scatter_3)
--       4 -> (ls_internal_term_gather_4, ls_internal_term_scatter_4)
--       _ -> error "Oops!"

-- toOperatorTerm' :: InteractionSpec -> Either Text OperatorTerm
-- toOperatorTerm' (InteractionSpec matrixSpec sitesSpec) = do
--   sites <- sitesListFromList sitesSpec
--   matrix <- denseToSparse =<< denseMatrixFromList matrixSpec
--   when (2 ^ (slTupleSize sites) /= ssmDimension matrix) $
--     Left $ "wrong matrix dimension"
--   Right $ OperatorTerm matrix sites

-- toOperatorTerm :: InteractionSpec -> OperatorTerm
-- toOperatorTerm spec = case toOperatorTerm' spec of
--   Right o -> o
--   Left e -> error e

-- foreign import capi "helpers.h ls_hs_apply_term"
--   ls_hs_apply_term :: Ptr Cterm -> Ptr Word64 -> Ptr Coutput_buffer -> IO ()

-- applyOperatorTerm' ::
--   OperatorTerm ->
--   [Word64] ->
--   IO (S.Vector Word64, S.Vector (Complex Double), Complex Double)
-- applyOperatorTerm' term bits =
--   withCterm term $ \term' -> with term' $ \termPtr ->
--     withArrayLen bits $ \numberWords bitsPtr -> do
--       outputSpins <- SM.new (bufferSize * numberWords)
--       outputCoeffs <- SM.new bufferSize
--       diagonal <- SM.unsafeWith outputSpins $ \spinsPtr ->
--         SM.unsafeWith outputCoeffs $ \coeffsPtr ->
--           alloca $ \diagonalPtr -> do
--             let (copyPtr, fillPtr) = case numberWords of
--                   1 -> (ls_internal_spin_copy_1, ls_internal_spin_fill_1)
--                 outputBuffer =
--                   Coutput_buffer
--                     spinsPtr
--                     coeffsPtr
--                     diagonalPtr
--                     (fromIntegral numberWords)
--                     copyPtr
--                     fillPtr
--             poke diagonalPtr 0
--             with outputBuffer $ \outPtr ->
--               ls_hs_apply_term termPtr bitsPtr outPtr
--             peek diagonalPtr
--       (,,) <$> S.freeze outputSpins
--         <*> S.freeze outputCoeffs
--         <*> pure diagonal
--   where
--     bufferSize = estimateBufferSizeForTerm term

-- applyOperatorTerm ::
--   OperatorTerm ->
--   Word64 ->
--   (S.Vector Word64, S.Vector (Complex Double), Complex Double)
-- applyOperatorTerm = undefined

-- toOperator :: FlatSpinBasis -> OperatorSpec -> SparseOperator
-- toOperator basis (OperatorSpec _ terms) = SparseOperator basis (toOperatorTerm <$> terms)

-- typedef struct ls_output_buffer {
--   uint64_t *spins;
--   _Complex double *coeffs;
--   _Complex double *const diagonal;
--   uint64_t const number_words;
--   ls_spin_copy_fn const spin_copy;
--   ls_spin_fill_fn const spin_fill;
-- } ls_output_buffer;
-- data {-# CTYPE "helpers.h" "ls_output_buffer" #-} Coutput_buffer
--   = Coutput_buffer
--       {-# UNPACK #-} !(Ptr Word64)
--       {-# UNPACK #-} !(Ptr (Complex Double))
--       {-# UNPACK #-} !(Ptr (Complex Double))
--       {-# UNPACK #-} !Word64
--       {-# UNPACK #-} !(FunPtr Cspin_copy_fn)
--       {-# UNPACK #-} !(FunPtr Cspin_fill_fn)

-- type Cspin_copy_fn = Ptr Word64 -> Ptr Word64 -> IO ()

-- type Cspin_fill_fn = Ptr Word64 -> Word64 -> Ptr Word64 -> IO ()

-- instance Storable Coutput_buffer where
--   sizeOf _ = 48
--   alignment _ = 8
--   peek p =
--     Coutput_buffer
--       <$> peekByteOff p 0
--       <*> peekByteOff p 8
--       <*> peekByteOff p 16
--       <*> peekByteOff p 24
--       <*> peekByteOff p 32
--       <*> peekByteOff p 40
--   poke p (Coutput_buffer spins coeffs diagonal number_words spin_copy spin_fill) = do
--     pokeByteOff p 0 spins
--     pokeByteOff p 8 coeffs
--     pokeByteOff p 16 diagonal
--     pokeByteOff p 24 number_words
--     pokeByteOff p 32 spin_copy
--     pokeByteOff p 40 spin_fill

-- trueCoutput_bufferSizeOf :: Int
-- trueCoutput_bufferSizeOf = fromIntegral [CU.pure| unsigned int { sizeof(ls_output_buffer) } |]

-- trueCoutput_bufferAlignment :: Int
-- trueCoutput_bufferAlignment = fromIntegral [CU.pure| unsigned int { __alignof__(ls_output_buffer) } |]

-- data {-# CTYPE "helpers.h" "ls_workspace" #-} Cworkspace
--   = Cworkspace
--       {-# UNPACK #-} !(Ptr Word64)
--       {-# UNPACK #-} !(Ptr (Complex Double))
--       {-# UNPACK #-} !(Ptr Double)

-- typedef struct ls_sparse_operator {
--   ls_flat_spin_basis const *basis;
--   unsigned number_terms;
--   ls_term *terms;
-- } ls_sparse_operator;
-- data {-# CTYPE "helpers.h" "ls_sparse_operator" #-} Csparse_operator
--   = Csparse_operator
--       {-# UNPACK #-} !(Ptr CFlatSpinBasis)
--       {-# UNPACK #-} !CUInt
--       {-# UNPACK #-} !(Ptr Cterm)

-- instance Storable Csparse_operator where
--   sizeOf _ = 24
--   alignment _ = 8
--   peek p =
--     Csparse_operator
--       <$> peekByteOff p 0
--       <*> peekByteOff p 8
--       <*> peekByteOff p 16
--   poke p (Csparse_operator basis number_terms terms) = do
--     pokeByteOff p 0 basis
--     pokeByteOff p 8 number_terms
--     pokeByteOff p 16 terms

-- trueCsparse_operatorSizeOf :: Int
-- trueCsparse_operatorSizeOf = fromIntegral [CU.pure| unsigned int { sizeof(ls_sparse_operator) } |]

-- trueCsparse_operatorAlignment :: Int
-- trueCsparse_operatorAlignment = fromIntegral [CU.pure| unsigned int { __alignof__(ls_sparse_operator) } |]

-- estimateBufferSizeForTerm :: OperatorTerm -> Int
-- estimateBufferSizeForTerm (OperatorTerm matrix (SitesList numberTuples _ _)) =
--   1 + maxNonZeroPerRow * numberTuples
--   where
--     maxNonZeroPerRow =
--       if ssmDimension matrix /= 0
--         then
--           fromIntegral . S.maximum $
--             S.generate
--               (ssmDimension matrix)
--               (\i -> scsrOffsets matrix ! (i + 1) - scsrOffsets matrix ! i)
--         else 0

-- data SparseOperator = SparseOperator FlatSpinBasis (NonEmpty OperatorTerm)

-- data SparseOperatorWrapper = SparseOperatorWrapper (Ptr Csparse_operator) (StablePtr SparseOperator)

-- mkSparseOperator :: FlatSpinBasis -> NonEmpty OperatorTerm -> SparseOperator
-- mkSparseOperator = SparseOperator

-- allocateCsparse_operator :: SparseOperator -> IO (Ptr Csparse_operator)
-- allocateCsparse_operator (SparseOperator (FlatSpinBasis basis) terms) = do
--   termsPtr <- allocateCterms terms
--   new $
--     Csparse_operator
--       (unsafeForeignPtrToPtr basis)
--       (fromIntegral $ NonEmpty.length terms)
--       termsPtr

-- deallocateCsparse_operator :: Ptr Csparse_operator -> IO ()
-- deallocateCsparse_operator p = do
--   (Csparse_operator _ _ termsPtr) <- peek p
--   deallocateCterms termsPtr
--   free p

-- allocateCsparse_operator ::

-- instance Storable SparseSquareMatrix where
--   sizeOf _ = 16
--   alignment _ = 4
--   peek p =
--     HalideDimension
--       <$> peekByteOff p 0
--       <*> peekByteOff p 4
--       <*> peekByteOff p 8
--       <*> peekByteOff p 12
--   poke p x = do
--     pokeByteOff p 0 (halideDimensionMin x)
--     pokeByteOff p 4 (halideDimensionExtent x)
--     pokeByteOff p 8 (halideDimensionStride x)
--     pokeByteOff p 12 (halideDimensionFlags x)

-- denseMatrixCountNonZero :: (Storable a, Eq a, Num a) => DenseMatrix v r c a -> Int
-- denseMatrixCountNonZero matrix = S.foldl' (\(!n) !x -> if x /= 0 then n + 1 else n) 0 (denseMatrixData matrix)

-- for (auto i = 0; i < numberRows; ++i) {
--   unsigned numberNonZero = 0;
--   for (auto j = 0; j < numberColumns; ++j) {
--     if (dense[i, j] != 0 && i != j) {
--       ++numberNonZero;
--       *(columns++) = j;
--       *(off_diag_elements++) = dense[i, j];
--     }
--   }
--   offsets[i + 1] = offsets[i] + numberNonZero;
-- }

-- denseToSparse :: DenseMatrix (Complex Double) -> Either Text SparseSquareMatrix
-- denseToSparse dense
--   | isDenseMatrixSquare dense = Right $
--     System.IO.Unsafe.unsafePerformIO $
--       do
--         let numberNonZero = denseMatrixCountNonZero dense
--             dimension = let (DenseMatrix (n, _) _) = dense in n
--         offsets <- SM.new (dimension + 1)
--         columns <- SM.new numberNonZero
--         offDiagElements <- SM.new numberNonZero
--         diagElements <- SM.new dimension
--         S.unsafeWith (denseMatrixData dense) $ \c_dense ->
--           SM.unsafeWith offsets $ \c_offsets ->
--             SM.unsafeWith columns $ \c_columns ->
--               SM.unsafeWith offDiagElements $ \c_off_diag_elements ->
--                 SM.unsafeWith diagElements $ \c_diag_elements ->
--                   ls_csr_matrix_from_dense
--                     (fromIntegral dimension)
--                     c_dense
--                     c_offsets
--                     c_columns
--                     c_off_diag_elements
--                     c_diag_elements
--         numberNonZero' <- fromIntegral <$> SM.read offsets dimension
--         SparseSquareMatrix dimension
--           <$> S.freeze offsets
--           <*> S.freeze (SM.take numberNonZero' columns)
--           <*> S.freeze (SM.take numberNonZero' offDiagElements)
--           <*> S.freeze diagElements
--   | otherwise = Left "expected a square matrix"

{- ORMOLU_DISABLE -}
-- foreign import capi unsafe "helpers.h ls_csr_matrix_from_dense"
--   ls_csr_matrix_from_dense :: CUInt -> Ptr (Complex Double) -> Ptr CUInt -> Ptr CUInt ->
--                               Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()

-- foreign import capi unsafe "helpers.h ls_dense_from_csr_matrix"
--   ls_dense_from_csr_matrix :: CUInt -> Ptr CUInt -> Ptr CUInt -> Ptr (Complex Double) ->
--                               Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()
{- ORMOLU_ENABLE -}

-- sparseToDense :: SparseSquareMatrix -> DenseMatrix (Complex Double)
-- sparseToDense sparse = System.IO.Unsafe.unsafePerformIO $ do
--   let dimension = ssmDimension sparse
--   elements <- SM.new (dimension * dimension)
--   S.unsafeWith (scsrOffsets sparse) $ \c_offsets ->
--     S.unsafeWith (ssmColumns sparse) $ \c_columns ->
--       S.unsafeWith (ssmOffDiagElements sparse) $ \c_off_diag_elements ->
--         S.unsafeWith (ssmDiagElements sparse) $ \c_diag_elements ->
--           SM.unsafeWith elements $ \c_dense ->
--             ls_dense_from_csr_matrix
--               (fromIntegral dimension)
--               c_offsets
--               c_columns
--               c_off_diag_elements
--               c_diag_elements
--               c_dense
--   DenseMatrix (dimension, dimension) <$> S.freeze elements

-- foreign import capi unsafe "helpers.h &ls_internal_term_gather_1"
--   ls_internal_term_gather_1 :: FunPtr (Ptr Word64 -> Ptr Cbit_index -> IO CUInt)

-- foreign import capi unsafe "helpers.h &ls_internal_term_gather_2"
--   ls_internal_term_gather_2 :: FunPtr (Ptr Word64 -> Ptr Cbit_index -> IO CUInt)

-- foreign import capi unsafe "helpers.h &ls_internal_term_gather_3"
--   ls_internal_term_gather_3 :: FunPtr (Ptr Word64 -> Ptr Cbit_index -> IO CUInt)

-- foreign import capi unsafe "helpers.h &ls_internal_term_gather_4"
--   ls_internal_term_gather_4 :: FunPtr (Ptr Word64 -> Ptr Cbit_index -> IO CUInt)

-- foreign import capi unsafe "helpers.h &ls_internal_term_scatter_1"
--   ls_internal_term_scatter_1 :: FunPtr (CUInt -> Ptr Cbit_index -> Ptr Word64 -> IO ())

-- foreign import capi unsafe "helpers.h &ls_internal_term_scatter_2"
--   ls_internal_term_scatter_2 :: FunPtr (CUInt -> Ptr Cbit_index -> Ptr Word64 -> IO ())

-- foreign import capi unsafe "helpers.h &ls_internal_term_scatter_3"
--   ls_internal_term_scatter_3 :: FunPtr (CUInt -> Ptr Cbit_index -> Ptr Word64 -> IO ())

-- foreign import capi unsafe "helpers.h &ls_internal_term_scatter_4"
--   ls_internal_term_scatter_4 :: FunPtr (CUInt -> Ptr Cbit_index -> Ptr Word64 -> IO ())

-- foreign import capi unsafe "helpers.h &ls_internal_spin_copy_1"
--   ls_internal_spin_copy_1 :: FunPtr (Ptr Word64 -> Ptr Word64 -> IO ())

-- foreign import capi unsafe "helpers.h &ls_internal_spin_fill_1"
--   ls_internal_spin_fill_1 :: FunPtr (Ptr Word64 -> Word64 -> Ptr Word64 -> IO ())
