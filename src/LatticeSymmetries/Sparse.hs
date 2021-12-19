{-# LANGUAGE CApiFFI #-}
{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE TemplateHaskell #-}

module LatticeSymmetries.Sparse where

-- import Data.Binary (Binary (..))

-- import Data.Vector.Binary

import Control.Exception.Safe (bracket, impureThrow, throwIO)
import qualified Control.Monad.Primitive as Primitive
import Control.Monad.ST
import qualified Control.Monad.ST.Unsafe (unsafeIOToST)
import Data.Aeson
import Data.Aeson.Types (typeMismatch)
import Data.Complex
import qualified Data.List
import qualified Data.List.NonEmpty as NonEmpty
import Data.Scientific (toRealFloat)
import qualified Data.Vector as B
import Data.Vector.Generic ((!))
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Storable as S
import qualified Data.Vector.Storable.Mutable as SM
import qualified Data.Vector.Unboxed as U
import Data.Yaml (decodeFileWithWarnings)
import Foreign.C.String (CString, peekCString)
import Foreign.C.Types (CInt (..), CUInt (..))
import Foreign.ForeignPtr
import Foreign.ForeignPtr.Unsafe (unsafeForeignPtrToPtr)
import Foreign.Marshal.Alloc (alloca, free, malloc, mallocBytes)
import Foreign.Marshal.Array (newArray, withArrayLen)
import Foreign.Marshal.Utils (new, with)
import Foreign.Ptr (FunPtr, Ptr, castPtr)
import Foreign.StablePtr
import Foreign.Storable (Storable (..))
import qualified GHC.Exts as GHC (IsList (..))
import qualified GHC.ForeignPtr as GHC (Finalizers (..), ForeignPtr (..), ForeignPtrContents (..))
import GHC.Generics
import qualified GHC.IORef as GHC (atomicSwapIORef)
import GHC.Prim
import qualified GHC.Ptr as GHC (Ptr (..))
import qualified Language.C.Inline as C
import qualified Language.C.Inline.Unsafe as CU
import LatticeSymmetries.Context
import LatticeSymmetries.IO
import LatticeSymmetries.Types
import qualified System.IO.Unsafe
import qualified System.Mem.Weak

C.context (C.baseCtx <> C.bsCtx <> C.funCtx <> lsCtx)
C.include "<lattice_symmetries/lattice_symmetries.h>"
C.include "helpers.h"

-- Row-major order (C layout)
data DenseMatrix a = DenseMatrix {denseMatrixShape :: (Int, Int), denseMatrixData :: S.Vector a}
  deriving stock (Show, Eq, Generic)

-- deriving anyclass (Binary)

data SparseMatrix a = SparseMatrix
  { smRows :: !Int,
    smCols :: !Int,
    smOffsets :: !(S.Vector CUInt),
    smIndices :: !(S.Vector CUInt),
    smData :: !(S.Vector a)
  }
  deriving stock (Show, Eq, Generic)

binarySearch :: (G.Vector v a, Ord a) => v a -> a -> Maybe Int
binarySearch v z = go 0 (G.length v)
  where
    go !l !u
      | l < u =
        -- NOTE: we assume that the vector is short enought such that @u + l@ does not overflow
        let !m = (u + l) `div` 2
         in case compare (v ! m) z of
              LT -> go (m + 1) u
              EQ -> Just m
              GT -> go l m
      | otherwise = Nothing

smIndex :: (Num a, Storable a) => SparseMatrix a -> (Int, Int) -> a
smIndex matrix (i, j) =
  case binarySearch (G.slice l (u - l) (smIndices matrix)) (fromIntegral j) of
    Just k -> smData matrix ! (l + k)
    Nothing -> 0
  where
    !l = fromIntegral $ smOffsets matrix ! i
    !u = fromIntegral $ smOffsets matrix ! (i + 1)
{-# INLINE smIndex #-}

preprocessCoo :: Num a => [(Int, Int, a)] -> [(Int, Int, a)]
preprocessCoo =
  fmap (Data.List.foldl1' (\(!i, !j, !x) (_, _, y) -> (i, j, x + y)))
    . Data.List.groupBy (\a b -> key a == key b)
    . sortOn key
  where
    key (i, j, _) = (i, j)

unsafeCooToCsr :: forall a. (Storable a, Num a) => [(Int, Int, a)] -> SparseMatrix a
unsafeCooToCsr coo = SparseMatrix n m offsets indices elements
  where
    coordinates = B.fromList coo
    indices = G.convert $ G.map (\(_, j, _) -> fromIntegral j) coordinates
    elements = G.convert $ G.map (\(_, _, x) -> x) coordinates
    n = (+ 1) . G.maximum . G.map (\(i, _, _) -> i) $ coordinates
    m = (+ 1) . fromIntegral . G.maximum $ indices
    offsets = runST $ do
      rs <- SM.replicate (n + 1) 0
      G.forM_ coordinates $ \(!i, _, _) ->
        SM.modify rs (+ 1) (i + 1)
      forM_ [0 .. n - 1] $ \(!i) -> do
        r <- SM.read rs i
        SM.modify rs (+ r) (i + 1)
      S.unsafeFreeze rs

cooToCsr :: forall a. (Storable a, Num a) => [(Int, Int, a)] -> SparseMatrix a
cooToCsr = unsafeCooToCsr . preprocessCoo

csrToCoo :: Storable a => SparseMatrix a -> [(Int, Int, a)]
csrToCoo csr = concat [row i | i <- [0 .. smRows csr - 1]]
  where
    row !i =
      let !l = fromIntegral $ smOffsets csr ! i
          !u = fromIntegral $ smOffsets csr ! (i + 1)
       in zipWith
            (\j x -> (i, fromIntegral j, x))
            (G.toList $ G.slice l (u - l) (smIndices csr))
            (G.toList $ G.slice l (u - l) (smData csr))

denseToCsr :: (Storable a, Eq a, Num a) => DenseMatrix a -> SparseMatrix a
denseToCsr dense =
  unsafeCooToCsr
    [ (i, j, x)
      | i <- [0 .. numberRows - 1],
        j <- [0 .. numberCols - 1],
        let x = indexDenseMatrix dense i j,
        x /= 0
    ]
  where
    (numberRows, numberCols) = denseMatrixShape dense

csrToDense :: Storable a => SparseMatrix a -> DenseMatrix a
csrToDense csr = runST $ do
  elements <- SM.new (n * m)
  forM_ (csrToCoo csr) $ \(i, j, x) ->
    SM.write elements (i * m + j) x
  DenseMatrix (n, m) <$> S.unsafeFreeze elements
  where
    !n = smRows csr
    !m = smCols csr

smKron :: (Storable a, Num a) => SparseMatrix a -> SparseMatrix a -> SparseMatrix a
smKron a b =
  cooToCsr
    [ (i₁ * n + i₂, j₁ * m + j₂, x₁ * x₂)
      | (i₁, j₁, x₁) <- csrToCoo a,
        (i₂, j₂, x₂) <- csrToCoo b
    ]
  where
    !n = smRows b
    !m = smCols b

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

isDenseMatrixSquare :: DenseMatrix a -> Bool
isDenseMatrixSquare (DenseMatrix (r, c) _) = r == c

isDenseMatrixEmpty :: DenseMatrix a -> Bool
isDenseMatrixEmpty (DenseMatrix (r, c) _) = r * c == 0

indexDenseMatrix :: Storable a => DenseMatrix a -> Int -> Int -> a
indexDenseMatrix (DenseMatrix (_, c) v) i j = v ! (c * i + j)

isDenseMatrixHermitian :: (Storable a, Num a, Eq a) => DenseMatrix (Complex a) -> Bool
isDenseMatrixHermitian matrix@(DenseMatrix (r, c) _) = isDenseMatrixSquare matrix && go 0 0
  where
    -- Iterate over upper triangle (including the diagonal) of the matrix
    go :: Int -> Int -> Bool
    go !i !j
      | j == c = let !i' = i + 1 in (i' >= r) || go i' i'
      | otherwise =
        let !p = indexDenseMatrix matrix i j == conjugate (indexDenseMatrix matrix j i)
         in p && go i (j + 1)

denseMatrixFromList :: Storable a => [[a]] -> Either Text (DenseMatrix a)
denseMatrixFromList [] = Right $ DenseMatrix (0, 0) S.empty
denseMatrixFromList (r : rs) =
  case all ((== nColumns) . length) rs of
    True -> Right $ DenseMatrix (length (r : rs), nColumns) v
    False -> Left "all rows of a matrix must have the same length"
  where
    nColumns = length r
    v = fromList $ mconcat (r : rs)

instance Storable a => GHC.IsList (DenseMatrix a) where
  type Item (DenseMatrix a) = [a]
  fromList rows = case denseMatrixFromList rows of
    Right m -> m
    Left msg -> error msg

data {-# CTYPE "ls_bit_index" #-} BitIndex = BitIndex !Word8 !Word8
  deriving stock (Show, Eq, Generic)

--  deriving anyclass (Binary)

-- instance Binary CUInt where
--   put (CUInt x) = Data.Binary.put x
--   get = CUInt <$> Data.Binary.get

data SparseSquareMatrix = SparseSquareMatrix
  { ssmDimension :: {-# UNPACK #-} !Int,
    ssmOffsets :: {-# UNPACK #-} !(S.Vector CUInt),
    ssmColumns :: {-# UNPACK #-} !(S.Vector CUInt),
    ssmOffDiagElements :: {-# UNPACK #-} !(S.Vector (Complex Double)),
    ssmDiagElements :: {-# UNPACK #-} !(S.Vector (Complex Double))
  }
  deriving stock (Show, Eq, Generic)

--   deriving anyclass (Binary)

isSparseMatrixHermitian :: SparseSquareMatrix -> Bool
isSparseMatrixHermitian = isDenseMatrixHermitian . sparseToDense

withCsparse_matrix :: SparseSquareMatrix -> (Csparse_matrix -> IO a) -> IO a
withCsparse_matrix matrix action =
  S.unsafeWith (ssmOffsets matrix) $ \offsetsPtr ->
    S.unsafeWith (ssmColumns matrix) $ \columnsPtr ->
      S.unsafeWith (ssmOffDiagElements matrix) $ \offDiagElementsPtr ->
        S.unsafeWith (ssmDiagElements matrix) $ \diagElementsPtr ->
          action $
            Csparse_matrix
              (fromIntegral . ssmDimension $ matrix)
              (fromIntegral . S.length . ssmColumns $ matrix)
              offsetsPtr
              columnsPtr
              offDiagElementsPtr
              diagElementsPtr

-- typedef struct ls_csr_matrix {
--     unsigned         dimension;
--     unsigned         number_nonzero;
--     unsigned*        offsets;
--     unsigned*        columns;
--     _Complex double* off_diag_elements;
--     _Complex double* diag_elements;
-- } ls_csr_matrix;
data {-# CTYPE "helpers.h" "ls_csr_matrix" #-} Csparse_matrix
  = Csparse_matrix
      {-# UNPACK #-} !CUInt
      {-# UNPACK #-} !CUInt
      {-# UNPACK #-} !(Ptr CUInt)
      {-# UNPACK #-} !(Ptr CUInt)
      {-# UNPACK #-} !(Ptr (Complex Double))
      {-# UNPACK #-} !(Ptr (Complex Double))

instance Storable Csparse_matrix where
  sizeOf _ = 40
  alignment _ = 8
  peek p =
    Csparse_matrix
      <$> peekByteOff p 0
      <*> peekByteOff p 4
      <*> peekByteOff p 8
      <*> peekByteOff p 16
      <*> peekByteOff p 24
      <*> peekByteOff p 32
  poke p (Csparse_matrix dimension number_nonzero offsets columns off_diag_elements diag_elements) = do
    pokeByteOff p 0 dimension
    pokeByteOff p 4 number_nonzero
    pokeByteOff p 8 offsets
    pokeByteOff p 16 columns
    pokeByteOff p 24 off_diag_elements
    pokeByteOff p 32 diag_elements

trueCsparse_matrixSizeOf :: Int
trueCsparse_matrixSizeOf = fromIntegral [CU.pure| unsigned int { sizeof(ls_csr_matrix) } |]

trueCsparse_matrixAlignment :: Int
trueCsparse_matrixAlignment = fromIntegral [CU.pure| unsigned int { __alignof__(ls_csr_matrix) } |]

-- typedef struct ls_bit_index {
--     uint8_t word;
--     uint8_t bit;
-- } ls_bit_index;
data {-# CTYPE "helpers.h" "ls_bit_index" #-} Cbit_index
  = Cbit_index {-# UNPACK #-} !Word8 {-# UNPACK #-} !Word8
  deriving (Show, Eq, Generic)

--  deriving anyclass (Binary)

instance Storable Cbit_index where
  sizeOf _ = 2
  alignment _ = 1
  peek p = Cbit_index <$> peekByteOff p 0 <*> peekByteOff p 1
  poke p (Cbit_index word bit) = pokeByteOff p 0 word >> pokeByteOff p 1 bit

trueCbit_indexSizeOf :: Int
trueCbit_indexSizeOf = fromIntegral [CU.pure| unsigned int { sizeof(ls_bit_index) } |]

trueCbit_indexAlignment :: Int
trueCbit_indexAlignment = fromIntegral [CU.pure| unsigned int { __alignof__(ls_bit_index) } |]

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
data {-# CTYPE "helpers.h" "ls_term" #-} Cterm
  = Cterm
      {-# UNPACK #-} !Csparse_matrix
      {-# UNPACK #-} !CUInt
      {-# UNPACK #-} !CUInt
      {-# UNPACK #-} !(Ptr Cbit_index)
      {-# UNPACK #-} !(FunPtr (Ptr Word64 -> Ptr Cbit_index -> IO CUInt))
      {-# UNPACK #-} !(FunPtr (CUInt -> Ptr Cbit_index -> Ptr Word64 -> IO ()))

instance Storable Cterm where
  sizeOf _ = 72
  alignment _ = 8
  peek p =
    Cterm
      <$> peekByteOff p 0
      <*> peekByteOff p 40
      <*> peekByteOff p 44
      <*> peekByteOff p 48
      <*> peekByteOff p 56
      <*> peekByteOff p 64
  poke p (Cterm matrix number_tuples tuple_size tuples gather_fn scatter_fn) = do
    pokeByteOff p 0 matrix
    pokeByteOff p 40 number_tuples
    pokeByteOff p 44 tuple_size
    pokeByteOff p 48 tuples
    pokeByteOff p 56 gather_fn
    pokeByteOff p 64 scatter_fn

trueCtermSizeOf :: Int
trueCtermSizeOf = fromIntegral [CU.pure| unsigned int { sizeof(ls_term) } |]

trueCtermAlignment :: Int
trueCtermAlignment = fromIntegral [CU.pure| unsigned int { __alignof__(ls_term) } |]

data SitesList = SitesList
  { slNumberTuples :: !Int,
    slTupleSize :: !Int,
    slData :: !(S.Vector Cbit_index)
  }
  deriving stock (Show, Eq, Generic)

--  deriving anyclass (Binary)

sitesListFromList :: [[Int]] -> Either Text SitesList
sitesListFromList rows = do
  (DenseMatrix (numberTuples, tupleSize) v) <- denseMatrixFromList rows
  when (S.any (< 0) v) $
    Left "sites list cannot have negative elements"
  when (S.any (> fromIntegral (maxBound :: Word16)) v) $
    Left "sites list cannot such large elements"
  Right $ SitesList numberTuples tupleSize (S.map toBitIndex v)
  where
    toBitIndex x = Cbit_index (fromIntegral (x `div` 64)) (fromIntegral (x `mod` 64))

-- newtype SitesList = SitesList [[Int]]

data OperatorTerm = OperatorTerm {otMatrix :: !SparseSquareMatrix, otSites :: !SitesList}
  deriving stock (Show, Eq, Generic)

--  deriving anyclass (Binary)

data {-# CTYPE "ls_hs_operator_term" #-} OperatorTermWrapper
  = OperatorTermWrapper
      {-# UNPACK #-} !(Ptr Cterm)
      {-# UNPACK #-} !(StablePtr OperatorTerm)

instance Storable OperatorTermWrapper where
  {-# INLINE sizeOf #-}
  sizeOf _ = 16
  {-# INLINE alignment #-}
  alignment _ = 8
  {-# INLINE peek #-}
  peek p = OperatorTermWrapper <$> peekByteOff p 0 <*> peekByteOff p 8
  {-# INLINE poke #-}
  poke p (OperatorTermWrapper term stable) =
    pokeByteOff p 0 term >> pokeByteOff p 8 stable

vectorFromPtr :: Storable a => Int -> Ptr a -> IO (S.Vector a)
vectorFromPtr n p = S.freeze =<< SM.unsafeFromForeignPtr0 <$> newForeignPtr_ p <*> pure n

denseMatrixFromPtr :: Int -> Int -> Ptr a -> IO (DenseMatrix a)
denseMatrixFromPtr = undefined

ls_hs_create_operator_term_from_dense ::
  CUInt ->
  Ptr (Complex Double) ->
  CUInt ->
  CUInt ->
  Ptr Word16 ->
  IO (Ptr OperatorTermWrapper)
ls_hs_create_operator_term_from_dense dimension matrixData numberTuples tupleSize tuplesData = do
  let dimension' = fromIntegral dimension
      numberTuples' = fromIntegral numberTuples
      tupleSize' = fromIntegral tupleSize
  matrixContents <- vectorFromPtr (dimension' * dimension') matrixData
  let matrix = case denseToSparse $ DenseMatrix (dimension', dimension') matrixContents of
        Right m -> m
        Left e -> error e
      toBitIndex x = Cbit_index (fromIntegral (x `div` 64)) (fromIntegral (x `mod` 64))
  sitesContents <- vectorFromPtr (numberTuples' * tupleSize') tuplesData
  let sites = SitesList numberTuples' tupleSize' (S.map toBitIndex sitesContents)
  when (2 ^ (slTupleSize sites) /= ssmDimension matrix) $
    error $ "wrong matrix dimension"
  let term = OperatorTerm matrix sites
  wrapper <- OperatorTermWrapper <$> allocateCterm term <*> newStablePtr term
  new wrapper

allocateCterm :: OperatorTerm -> IO (Ptr Cterm)
allocateCterm term = do
  p <- malloc
  withCterm term (poke p)
  return p

deallocateCterm :: Ptr Cterm -> IO ()
deallocateCterm p = free p

allocateCterms :: NonEmpty OperatorTerm -> IO (Ptr Cterm)
allocateCterms terms = do
  let !count = NonEmpty.length terms
      !elemSize = let x = x in sizeOf (x :: Cterm)
  p <- mallocBytes $ elemSize * count
  forM_ (NonEmpty.zip terms (fromList [0 ..])) $ \(t, i) -> withCterm t (pokeElemOff p i)
  return p

deallocateCterms :: Ptr Cterm -> IO ()
deallocateCterms p = free p

withCterm :: OperatorTerm -> (Cterm -> IO a) -> IO a
withCterm (OperatorTerm matrix sites) action =
  withCsparse_matrix matrix $ \matrix' ->
    S.unsafeWith (slData sites) $ \tuplesPtr ->
      action $
        Cterm
          matrix'
          (fromIntegral . slNumberTuples $ sites)
          (fromIntegral . slTupleSize $ sites)
          tuplesPtr
          gatherPtr
          scatterPtr
  where
    (gatherPtr, scatterPtr) = case slTupleSize sites of
      1 -> (ls_internal_term_gather_1, ls_internal_term_scatter_1)
      2 -> (ls_internal_term_gather_2, ls_internal_term_scatter_2)
      3 -> (ls_internal_term_gather_3, ls_internal_term_scatter_3)
      4 -> (ls_internal_term_gather_4, ls_internal_term_scatter_4)
      _ -> error "Oops!"

toOperatorTerm' :: InteractionSpec -> Either Text OperatorTerm
toOperatorTerm' (InteractionSpec matrixSpec sitesSpec) = do
  sites <- sitesListFromList sitesSpec
  matrix <- denseToSparse =<< denseMatrixFromList matrixSpec
  when (2 ^ (slTupleSize sites) /= ssmDimension matrix) $
    Left $ "wrong matrix dimension"
  Right $ OperatorTerm matrix sites

toOperatorTerm :: InteractionSpec -> OperatorTerm
toOperatorTerm spec = case toOperatorTerm' spec of
  Right o -> o
  Left e -> error e

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

applyOperatorTerm ::
  OperatorTerm ->
  Word64 ->
  (S.Vector Word64, S.Vector (Complex Double), Complex Double)
applyOperatorTerm = undefined

toOperator :: FlatSpinBasis -> OperatorSpec -> SparseOperator
toOperator basis (OperatorSpec _ terms) = SparseOperator basis (toOperatorTerm <$> terms)

-- typedef struct ls_output_buffer {
--   uint64_t *spins;
--   _Complex double *coeffs;
--   _Complex double *const diagonal;
--   uint64_t const number_words;
--   ls_spin_copy_fn const spin_copy;
--   ls_spin_fill_fn const spin_fill;
-- } ls_output_buffer;
data {-# CTYPE "helpers.h" "ls_output_buffer" #-} Coutput_buffer
  = Coutput_buffer
      {-# UNPACK #-} !(Ptr Word64)
      {-# UNPACK #-} !(Ptr (Complex Double))
      {-# UNPACK #-} !(Ptr (Complex Double))
      {-# UNPACK #-} !Word64
      {-# UNPACK #-} !(FunPtr Cspin_copy_fn)
      {-# UNPACK #-} !(FunPtr Cspin_fill_fn)

type Cspin_copy_fn = Ptr Word64 -> Ptr Word64 -> IO ()

type Cspin_fill_fn = Ptr Word64 -> Word64 -> Ptr Word64 -> IO ()

instance Storable Coutput_buffer where
  sizeOf _ = 48
  alignment _ = 8
  peek p =
    Coutput_buffer
      <$> peekByteOff p 0
      <*> peekByteOff p 8
      <*> peekByteOff p 16
      <*> peekByteOff p 24
      <*> peekByteOff p 32
      <*> peekByteOff p 40
  poke p (Coutput_buffer spins coeffs diagonal number_words spin_copy spin_fill) = do
    pokeByteOff p 0 spins
    pokeByteOff p 8 coeffs
    pokeByteOff p 16 diagonal
    pokeByteOff p 24 number_words
    pokeByteOff p 32 spin_copy
    pokeByteOff p 40 spin_fill

trueCoutput_bufferSizeOf :: Int
trueCoutput_bufferSizeOf = fromIntegral [CU.pure| unsigned int { sizeof(ls_output_buffer) } |]

trueCoutput_bufferAlignment :: Int
trueCoutput_bufferAlignment = fromIntegral [CU.pure| unsigned int { __alignof__(ls_output_buffer) } |]

data {-# CTYPE "helpers.h" "ls_workspace" #-} Cworkspace
  = Cworkspace
      {-# UNPACK #-} !(Ptr Word64)
      {-# UNPACK #-} !(Ptr (Complex Double))
      {-# UNPACK #-} !(Ptr Double)

-- typedef struct ls_sparse_operator {
--   ls_flat_spin_basis const *basis;
--   unsigned number_terms;
--   ls_term *terms;
-- } ls_sparse_operator;
data {-# CTYPE "helpers.h" "ls_sparse_operator" #-} Csparse_operator
  = Csparse_operator
      {-# UNPACK #-} !(Ptr CFlatSpinBasis)
      {-# UNPACK #-} !CUInt
      {-# UNPACK #-} !(Ptr Cterm)

instance Storable Csparse_operator where
  sizeOf _ = 24
  alignment _ = 8
  peek p =
    Csparse_operator
      <$> peekByteOff p 0
      <*> peekByteOff p 8
      <*> peekByteOff p 16
  poke p (Csparse_operator basis number_terms terms) = do
    pokeByteOff p 0 basis
    pokeByteOff p 8 number_terms
    pokeByteOff p 16 terms

trueCsparse_operatorSizeOf :: Int
trueCsparse_operatorSizeOf = fromIntegral [CU.pure| unsigned int { sizeof(ls_sparse_operator) } |]

trueCsparse_operatorAlignment :: Int
trueCsparse_operatorAlignment = fromIntegral [CU.pure| unsigned int { __alignof__(ls_sparse_operator) } |]

estimateBufferSizeForTerm :: OperatorTerm -> Int
estimateBufferSizeForTerm (OperatorTerm matrix (SitesList numberTuples _ _)) =
  1 + maxNonZeroPerRow * numberTuples
  where
    maxNonZeroPerRow =
      if ssmDimension matrix /= 0
        then
          fromIntegral . S.maximum $
            S.generate
              (ssmDimension matrix)
              (\i -> ssmOffsets matrix ! (i + 1) - ssmOffsets matrix ! i)
        else 0

data SparseOperator = SparseOperator FlatSpinBasis (NonEmpty OperatorTerm)

data SparseOperatorWrapper = SparseOperatorWrapper (Ptr Csparse_operator) (StablePtr SparseOperator)

mkSparseOperator :: FlatSpinBasis -> NonEmpty OperatorTerm -> SparseOperator
mkSparseOperator = SparseOperator

allocateCsparse_operator :: SparseOperator -> IO (Ptr Csparse_operator)
allocateCsparse_operator (SparseOperator (FlatSpinBasis basis) terms) = do
  termsPtr <- allocateCterms terms
  new $
    Csparse_operator
      (unsafeForeignPtrToPtr basis)
      (fromIntegral $ NonEmpty.length terms)
      termsPtr

deallocateCsparse_operator :: Ptr Csparse_operator -> IO ()
deallocateCsparse_operator p = do
  (Csparse_operator _ _ termsPtr) <- peek p
  deallocateCterms termsPtr
  free p

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

denseMatrixCountNonZero :: (Storable a, Eq a, Num a) => DenseMatrix a -> Int
denseMatrixCountNonZero matrix = S.foldl' (\(!n) !x -> if x /= 0 then n + 1 else n) 0 (denseMatrixData matrix)

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

denseToSparse :: DenseMatrix (Complex Double) -> Either Text SparseSquareMatrix
denseToSparse dense
  | isDenseMatrixSquare dense = Right $
    System.IO.Unsafe.unsafePerformIO $
      do
        let numberNonZero = denseMatrixCountNonZero dense
            dimension = let (DenseMatrix (n, _) _) = dense in n
        offsets <- SM.new (dimension + 1)
        columns <- SM.new numberNonZero
        offDiagElements <- SM.new numberNonZero
        diagElements <- SM.new dimension
        S.unsafeWith (denseMatrixData dense) $ \c_dense ->
          SM.unsafeWith offsets $ \c_offsets ->
            SM.unsafeWith columns $ \c_columns ->
              SM.unsafeWith offDiagElements $ \c_off_diag_elements ->
                SM.unsafeWith diagElements $ \c_diag_elements ->
                  ls_csr_matrix_from_dense
                    (fromIntegral dimension)
                    c_dense
                    c_offsets
                    c_columns
                    c_off_diag_elements
                    c_diag_elements
        numberNonZero' <- fromIntegral <$> SM.read offsets dimension
        SparseSquareMatrix dimension
          <$> S.freeze offsets
          <*> S.freeze (SM.take numberNonZero' columns)
          <*> S.freeze (SM.take numberNonZero' offDiagElements)
          <*> S.freeze diagElements
  | otherwise = Left "expected a square matrix"

{- ORMOLU_DISABLE -}
foreign import capi unsafe "helpers.h ls_csr_matrix_from_dense"
  ls_csr_matrix_from_dense :: CUInt -> Ptr (Complex Double) -> Ptr CUInt -> Ptr CUInt ->
                              Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()

foreign import capi unsafe "helpers.h ls_dense_from_csr_matrix"
  ls_dense_from_csr_matrix :: CUInt -> Ptr CUInt -> Ptr CUInt -> Ptr (Complex Double) ->
                              Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()
{- ORMOLU_ENABLE -}

sparseToDense :: SparseSquareMatrix -> DenseMatrix (Complex Double)
sparseToDense sparse = System.IO.Unsafe.unsafePerformIO $ do
  let dimension = ssmDimension sparse
  elements <- SM.new (dimension * dimension)
  S.unsafeWith (ssmOffsets sparse) $ \c_offsets ->
    S.unsafeWith (ssmColumns sparse) $ \c_columns ->
      S.unsafeWith (ssmOffDiagElements sparse) $ \c_off_diag_elements ->
        S.unsafeWith (ssmDiagElements sparse) $ \c_diag_elements ->
          SM.unsafeWith elements $ \c_dense ->
            ls_dense_from_csr_matrix
              (fromIntegral dimension)
              c_offsets
              c_columns
              c_off_diag_elements
              c_diag_elements
              c_dense
  DenseMatrix (dimension, dimension) <$> S.freeze elements

foreign import capi unsafe "helpers.h &ls_internal_term_gather_1"
  ls_internal_term_gather_1 :: FunPtr (Ptr Word64 -> Ptr Cbit_index -> IO CUInt)

foreign import capi unsafe "helpers.h &ls_internal_term_gather_2"
  ls_internal_term_gather_2 :: FunPtr (Ptr Word64 -> Ptr Cbit_index -> IO CUInt)

foreign import capi unsafe "helpers.h &ls_internal_term_gather_3"
  ls_internal_term_gather_3 :: FunPtr (Ptr Word64 -> Ptr Cbit_index -> IO CUInt)

foreign import capi unsafe "helpers.h &ls_internal_term_gather_4"
  ls_internal_term_gather_4 :: FunPtr (Ptr Word64 -> Ptr Cbit_index -> IO CUInt)

foreign import capi unsafe "helpers.h &ls_internal_term_scatter_1"
  ls_internal_term_scatter_1 :: FunPtr (CUInt -> Ptr Cbit_index -> Ptr Word64 -> IO ())

foreign import capi unsafe "helpers.h &ls_internal_term_scatter_2"
  ls_internal_term_scatter_2 :: FunPtr (CUInt -> Ptr Cbit_index -> Ptr Word64 -> IO ())

foreign import capi unsafe "helpers.h &ls_internal_term_scatter_3"
  ls_internal_term_scatter_3 :: FunPtr (CUInt -> Ptr Cbit_index -> Ptr Word64 -> IO ())

foreign import capi unsafe "helpers.h &ls_internal_term_scatter_4"
  ls_internal_term_scatter_4 :: FunPtr (CUInt -> Ptr Cbit_index -> Ptr Word64 -> IO ())

foreign import capi unsafe "helpers.h &ls_internal_spin_copy_1"
  ls_internal_spin_copy_1 :: FunPtr (Ptr Word64 -> Ptr Word64 -> IO ())

foreign import capi unsafe "helpers.h &ls_internal_spin_fill_1"
  ls_internal_spin_fill_1 :: FunPtr (Ptr Word64 -> Word64 -> Ptr Word64 -> IO ())
