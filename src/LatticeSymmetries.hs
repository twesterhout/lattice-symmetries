{-# LANGUAGE CApiFFI #-}
{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE TemplateHaskell #-}

-- {-# OPTIONS_GHC -fno-warn-orphans #-}

module LatticeSymmetries where

import qualified Control.Monad.Primitive as Primitive
import Control.Monad.ST
import qualified Control.Monad.ST.Unsafe (unsafeIOToST)
import Data.Aeson
import Data.Aeson.Types (typeMismatch)
import Data.Complex
import qualified Data.List.NonEmpty as NonEmpty
import Data.Scientific (toRealFloat)
import Data.Vector.Generic ((!))
import qualified Data.Vector.Storable as S
import qualified Data.Vector.Storable.Mutable as SM
import qualified Data.Vector.Unboxed as U
import Data.Yaml (decodeFileWithWarnings)
import Foreign.C.String (CString, peekCString)
import Foreign.C.Types (CInt (..), CUInt (..))
import Foreign.ForeignPtr
import Foreign.Marshal.Alloc (alloca)
import Foreign.Marshal.Array (withArrayLen)
import Foreign.Ptr (FunPtr, Ptr)
import Foreign.Storable (Storable (..))
import qualified GHC.Exts as GHC (IsList (..))
import qualified GHC.ForeignPtr as GHC (Finalizers (..), ForeignPtr (..), ForeignPtrContents (..))
import qualified GHC.IORef as GHC (atomicSwapIORef)
import qualified GHC.Ptr as GHC (Ptr (..))
import qualified Language.C.Inline as C
import qualified Language.C.Inline.Unsafe as CU
import LatticeSymmetries.Context
import qualified System.IO.Unsafe
import UnliftIO.Exception (bracket, impureThrow, throwIO)

C.context (C.baseCtx <> C.bsCtx <> C.funCtx <> lsCtx)
C.include "<lattice_symmetries/lattice_symmetries.h>"
C.include "helpers.h"

foo :: IO CInt
foo = do
  putStrLn "Hello world!"
  return 123

foreign export ccall "ls_hs_foo" foo :: IO CInt

-- | Exceptions thrown when an error occurs in @liblattice_symmetries@.
data LatticeSymmetriesException = LatticeSymmetriesException {eCode :: Int, eMessage :: Text}
  deriving stock (Show)

instance Exception LatticeSymmetriesException

data SpinEDException = SpinEDException Text
  deriving stock (Show)

instance Exception SpinEDException

-- | Retrieve textual representation of an error
getErrorMessage ::
  -- | Error code returned by @lattice_symmetries@ library
  Int ->
  -- | Explanation of the error
  IO Text
getErrorMessage c = bracket (ls_error_to_string (fromIntegral c)) ls_destroy_string $ \s ->
  toText <$> peekCString s

-- | Check the status code returned by @lattice_symmetries@ library. If it
-- indicates an error, 'LatticeSymmetriesException' is thrown.
checkStatus :: (MonadIO m, Integral a) => a -> m ()
checkStatus c
  | c == 0 = return ()
  | otherwise = do
    print =<< liftIO (getErrorMessage c')
    throwIO . LatticeSymmetriesException c' =<< liftIO (getErrorMessage c')
  where
    c' = fromIntegral c

data SymmetrySpec = SymmetrySpec !(NonEmpty Int) !Int
  deriving stock (Read, Show, Eq)

instance FromJSON SymmetrySpec where
  parseJSON = withObject "symmetry" $ \v ->
    SymmetrySpec
      <$> v .: "permutation"
      <*> v .: "sector"

instance ToJSON SymmetrySpec where
  toJSON (SymmetrySpec permutation sector) =
    object ["permutation" .= permutation, "sector" .= sector]

toSymmetry :: MonadIO m => SymmetrySpec -> m Symmetry
toSymmetry (SymmetrySpec p s) = mkSymmetry p s

data BasisSpec = BasisSpec !Int !(Maybe Int) !(Maybe Int) ![SymmetrySpec]
  deriving stock (Read, Show, Eq)

instance FromJSON BasisSpec where
  parseJSON = withObject "basis" $ \v ->
    BasisSpec
      <$> v .: "number_spins"
      <*> v .:? "hamming_weight"
      <*> v .:? "spin_inversion"
      <*> v .: "symmetries"

instance ToJSON BasisSpec where
  toJSON (BasisSpec numberSpins hammingWeight spinInversion symmetries) =
    object
      [ "number_spins" .= numberSpins,
        "hamming_weight" .= maybe Null toJSON hammingWeight,
        "spin_inversion" .= maybe Null toJSON spinInversion,
        "symmetries" .= symmetries
      ]

newtype WrappedBasisSpec = WrappedBasisSpec BasisSpec

instance FromJSON WrappedBasisSpec where
  parseJSON = withObject "config" $ \v -> WrappedBasisSpec <$> v .: "basis"

instance FromJSON (Complex Double) where
  parseJSON (Number x) = pure . fromReal . toRealFloat $ x
    where
      fromReal :: Num a => a -> Complex a
      fromReal x' = x' :+ 0
  parseJSON v@(Array xs) = case (toList xs) of
    [re, im] -> (:+) <$> parseJSON re <*> parseJSON im
    _ -> typeMismatch "Complex" v
  parseJSON v = typeMismatch "Complex" v

newtype DenseMatrixSpec = DenseMatrixSpec [[Complex Double]]
  deriving stock (Read, Show, Eq)
  deriving newtype (FromJSON)

-- Row-major order (C layout)
data DenseMatrix = DenseMatrix {denseMatrixShape :: (Int, Int), denseMatrixData :: S.Vector (Complex Double)}
  deriving stock (Show, Eq)

isDenseMatrixSquare :: DenseMatrix -> Bool
isDenseMatrixSquare (DenseMatrix (r, c) _) = r == c

isDenseMatrixEmpty :: DenseMatrix -> Bool
isDenseMatrixEmpty (DenseMatrix (r, c) _) = r * c == 0

indexDenseMatrix :: DenseMatrix -> Int -> Int -> Complex Double
indexDenseMatrix (DenseMatrix (_, c) v) i j = v ! (c * i + j)

isDenseMatrixHermitian :: DenseMatrix -> Bool
isDenseMatrixHermitian matrix@(DenseMatrix (r, c) _) = isDenseMatrixSquare matrix && go 0 0
  where
    -- Iterate over upper triangle (including the diagonal) of the matrix
    go :: Int -> Int -> Bool
    go !i !j
      | j == c = let !i' = i + 1 in (i' >= r) || go i' i'
      | otherwise =
        let !p = indexDenseMatrix matrix i j == conjugate (indexDenseMatrix matrix j i)
         in p && go i (j + 1)

denseMatrixFromList :: [[Complex Double]] -> Maybe DenseMatrix
denseMatrixFromList [] = Just $ DenseMatrix (0, 0) S.empty
denseMatrixFromList (r : rs) =
  case all ((== nColumns) . length) rs of
    True -> Just $ DenseMatrix (length (r : rs), nColumns) v
    False -> Nothing
  where
    nColumns = length r
    v = fromList $ mconcat (r : rs)

instance GHC.IsList DenseMatrix where
  type Item DenseMatrix = [Complex Double]
  fromList rows = case denseMatrixFromList rows of
    Just m -> m
    Nothing -> error $ "all rows must have the same length"

-- data InteractionSpec = InteractionSpec ![[Complex Double]] ![[Int]]
--   deriving stock (Read, Show, Eq)
--
-- instance FromJSON InteractionSpec where
--   parseJSON = withObject "interaction" $ \v ->
--     InteractionSpec
--       <$> v .: "matrix"
--       <*> v .: "sites"

loadRawBasis :: MonadIO m => Text -> m (Ptr CSpinBasis)
loadRawBasis path = do
  r <- liftIO $ decodeFileWithWarnings (toString path)
  case r of
    Left e -> throwIO e
    Right (warnings, (WrappedBasisSpec basisSpec)) -> do
      mapM_ print warnings
      toRawBasis basisSpec

ls_hs_load_basis_from_yaml :: CString -> IO (Ptr CSpinBasis)
ls_hs_load_basis_from_yaml c_str = loadRawBasis . toText =<< peekCString c_str

foreign export ccall "ls_hs_load_basis_from_yaml"
  ls_hs_load_basis_from_yaml :: CString -> IO (Ptr CSpinBasis)

{-
isPermutation :: Integral a => NonEmpty a -> Bool
isPermutation xs = and . toList $ NonEmpty.zipWith (==) (NonEmpty.sort xs) (fromList [0 ..])

permuteArray :: Integral a => NonEmpty a -> U.Vector Int -> U.Vector Int
permuteArray permutation xs = (xs U.!!) undefined

computePeriodicity :: Integral a => NonEmpty a -> Int
computePeriodicity = undefined

checkSymmetrySpec :: SymmetrySpec -> Either Text SymmetrySpec
checkSymmetrySpec (SymmetrySpec permutation sector) = do
  when (not . isPermutation $ permutation) $
    Left $ show permutation <> " is not a permutation of " <> show [0 .. length permutation]
  undefined
-}

data {-# CTYPE "lattice_symmetries/lattice_symmetries.h" "ls_symmetry" #-} CSymmetry

data {-# CTYPE "lattice_symmetries/lattice_symmetries.h" "ls_group" #-} CGroup

data {-# CTYPE "lattice_symmetries/lattice_symmetries.h" "ls_spin_basis" #-} CSpinBasis

data {-# CTYPE "lattice_symmetries/lattice_symmetries.h" "ls_flat_spin_basis" #-} CFlatSpinBasis

newtype Symmetry = Symmetry (ForeignPtr CSymmetry)

newtype SymmetryGroup = SymmetryGroup (ForeignPtr CGroup)

newtype SpinBasis = SpinBasis (ForeignPtr CSpinBasis)

newtype FlatSpinBasis = FlatSpinBasis (ForeignPtr CFlatSpinBasis)

mkSymmetry ::
  MonadIO m =>
  -- | Permutation
  NonEmpty Int ->
  -- | Symmetry sector
  Int ->
  m Symmetry
mkSymmetry !permutation !sector = do
  -- Make sure permutation and sector can be safely converted to unsigned
  -- representations. Everything else is checked by ls_create_symmetry
  when (any (< 0) . toList $ permutation) . throwIO . SpinEDException $
    "invalid permutation: " <> show permutation <> "; indices must be non-negative"
  when (sector < 0) . throwIO . SpinEDException $
    "invalid sector: " <> show sector <> "; expected a non-negative number"
  ptr <- liftIO . mkObject $ \ptrPtr ->
    withArrayLen (fromIntegral <$> toList permutation) $ \n permutationPtr ->
      ls_create_symmetry ptrPtr (fromIntegral n) permutationPtr (fromIntegral sector)
  fmap Symmetry . liftIO $ newForeignPtr ls_destroy_symmetry ptr

mkSymmetryGroup ::
  MonadIO m =>
  -- | Symmetry generators
  [Symmetry] ->
  m SymmetryGroup
mkSymmetryGroup !xs = do
  ptr <- liftIO . mkObject $ \ptrPtr ->
    withSymmetries xs $ \n xsPtr ->
      ls_create_group ptrPtr (fromIntegral n) xsPtr
  fmap SymmetryGroup . liftIO $ newForeignPtr ls_destroy_group ptr

mkRawBasis ::
  MonadIO m =>
  -- | Symmetry group
  SymmetryGroup ->
  -- | Number of spins
  Int ->
  -- | Hamming weight
  Maybe Int ->
  -- | Spin inversion
  Maybe Int ->
  m (Ptr CSpinBasis)
mkRawBasis !(SymmetryGroup g) !numberSpins !hammingWeight !spinInversion = do
  when (numberSpins <= 0) . throwIO . SpinEDException $
    "invalid number of spins: " <> show numberSpins <> "; expected a positive number"
  hammingWeight' <- case hammingWeight of
    Just x -> do
      when (x < 0) . throwIO . SpinEDException $
        "invalid Hamming weight: " <> show x <> "; expected a non-negative number"
      return $ fromIntegral x
    Nothing -> return (-1)
  spinInversion' <- case spinInversion of
    Just x -> do
      when (x /= 1 && x /= -1) . throwIO . SpinEDException $
        "invalid value for spin inversion: " <> show x <> "; expected either -1 or +1"
      return $ fromIntegral x
    Nothing -> return 0
  ptr <- liftIO . mkObject $ \ptrPtr ->
    withForeignPtr g $ \groupPtr ->
      ls_create_spin_basis ptrPtr groupPtr (fromIntegral numberSpins) hammingWeight' spinInversion'
  return ptr

mkBasis raw = fmap SpinBasis . liftIO $ newForeignPtr ls_destroy_spin_basis raw

toRawBasis :: MonadIO m => BasisSpec -> m (Ptr CSpinBasis)
toRawBasis (BasisSpec numberSpins hammingWeight spinInversion symmetries) = do
  -- We need to make sure we use "small" basis, otherwise we won't be able to
  -- build a list of representatives later
  when (numberSpins > 64) . throwIO . SpinEDException $
    "invalid number_spins: " <> show numberSpins <> "; exact diagonalization is not feasible "
      <> "for systems larger than 64 spins"
  -- logDebug "Building symmetry group..."
  symmetryGroup <- mkSymmetryGroup =<< mapM toSymmetry symmetries
  -- logInfo $ "Symmetry group contains " <> show (getGroupSize symmetryGroup) <> " elements"
  mkRawBasis symmetryGroup numberSpins hammingWeight spinInversion

-- | Extension of 'withForeignPtr' to lists of 'ForeignPtr's.
withManyForeignPtr :: [ForeignPtr a] -> (Int -> Ptr (Ptr a) -> IO b) -> IO b
withManyForeignPtr xs func = loop [] xs
  where
    loop acc (y : ys) = withForeignPtr y $ \y' -> loop (y' : acc) ys
    loop acc [] = withArrayLen (reverse acc) func

withSymmetries :: [Symmetry] -> (Int -> Ptr (Ptr CSymmetry) -> IO a) -> IO a
withSymmetries xs func = withManyForeignPtr pointers func
  where
    pointers = (\(Symmetry p) -> p) <$> xs

mkObject :: (Ptr (Ptr a) -> IO CInt) -> IO (Ptr a)
mkObject f = alloca $ \ptrPtr -> f ptrPtr >>= checkStatus >> peek ptrPtr

releaseForeignPtr :: ForeignPtr a -> IO (Ptr a)
releaseForeignPtr (GHC.ForeignPtr p foreignPtr) = do
  _ <- GHC.atomicSwapIORef refFinalizers GHC.NoFinalizers
  putStrLn "Shoot!"
  return (GHC.Ptr p)
  where
    refFinalizers = case foreignPtr of
      (GHC.PlainForeignPtr ref) -> ref
      _ -> error "releaseForeignPtr: only 'PlainForeignPtr's are supported"

data {-# CTYPE "ls_bit_index" #-} BitIndex = BitIndex !Word8 !Word8

data SparseSquareMatrix = SparseSquareMatrix
  { ssmDimension :: {-# UNPACK #-} !Int,
    ssmOffsets :: {-# UNPACK #-} !(S.Vector CUInt),
    ssmColumns :: {-# UNPACK #-} !(S.Vector CUInt),
    ssmOffDiagElements :: {-# UNPACK #-} !(S.Vector (Complex Double)),
    ssmDiagElements :: {-# UNPACK #-} !(S.Vector (Complex Double))
  }
  deriving stock (Show, Eq)

isSparseMatrixHermitian :: SparseSquareMatrix -> Bool
isSparseMatrixHermitian = isDenseMatrixHermitian . sparseToDense

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
data {-# CTYPE "helpers.h" "ls_bit_index" #-} Cbit_index = Cbit_index {-# UNPACK #-} !Word8 {-# UNPACK #-} !Word8

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

-- typedef struct ls_output_buffer {
--   uint64_t *spins;
--   _Complex double *coeffs;
--   _Complex double *const diagonal;
--   uint64_t number_words;
-- } ls_output_buffer;
data {-# CTYPE "helpers.h" "ls_output_buffer" #-} Coutput_buffer
  = Coutput_buffer
      {-# UNPACK #-} !(Ptr Word64)
      {-# UNPACK #-} !(Ptr (Complex Double))
      {-# UNPACK #-} !(Ptr (Complex Double))
      {-# UNPACK #-} !Word64

instance Storable Coutput_buffer where
  sizeOf _ = 32
  alignment _ = 8
  peek p =
    Coutput_buffer
      <$> peekByteOff p 0
      <*> peekByteOff p 8
      <*> peekByteOff p 16
      <*> peekByteOff p 24
  poke p (Coutput_buffer spins coeffs diagonal number_words) = do
    pokeByteOff p 0 spins
    pokeByteOff p 8 coeffs
    pokeByteOff p 16 diagonal
    pokeByteOff p 24 number_words

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

newtype SitesList = SitesList [[Int]]

data Term = Term {termMatrix :: !SparseSquareMatrix, termSites :: !SitesList}

estimateBufferSizeForTerm :: Term -> Int
estimateBufferSizeForTerm (Term matrix (SitesList sites)) = 1 + maxNonZeroPerRow * length sites
  where
    maxNonZeroPerRow =
      if ssmDimension matrix /= 0
        then
          fromIntegral . S.maximum $
            S.generate (ssmDimension matrix) (\i -> ssmOffsets matrix ! (i + 1) - ssmOffsets matrix ! i)
        else 0

data SparseOperator = SparseOperator FlatSpinBasis [Term]

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

denseMatrixCountNonZero :: DenseMatrix -> Int
denseMatrixCountNonZero matrix = S.foldl' (\n x -> if x /= 0 then n + 1 else n) 0 (denseMatrixData matrix)

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

denseToCSRSquare :: DenseMatrix -> Maybe SparseSquareMatrix
denseToCSRSquare dense
  | isDenseMatrixSquare dense = Just $
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
  | otherwise = Nothing

foreign import capi unsafe "helpers.h ls_csr_matrix_from_dense"
  ls_csr_matrix_from_dense :: CUInt -> Ptr (Complex Double) -> Ptr CUInt -> Ptr CUInt -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()

foreign import capi unsafe "helpers.h ls_dense_from_csr_matrix"
  ls_dense_from_csr_matrix :: CUInt -> Ptr CUInt -> Ptr CUInt -> Ptr (Complex Double) -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()

sparseToDense :: SparseSquareMatrix -> DenseMatrix
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

foreign import ccall unsafe "lattice_symmetries/lattice_symmetries.h ls_error_to_string"
  ls_error_to_string :: CInt -> IO CString

foreign import capi unsafe "lattice_symmetries/lattice_symmetries.h ls_destroy_string"
  ls_destroy_string :: CString -> IO ()

foreign import capi unsafe "lattice_symmetries/lattice_symmetries.h ls_enable_logging"
  ls_enable_logging :: IO ()

foreign import capi unsafe "lattice_symmetries/lattice_symmetries.h ls_disable_logging"
  ls_disable_logging :: IO ()

foreign import capi unsafe "lattice_symmetries/lattice_symmetries.h ls_create_symmetry"
  ls_create_symmetry :: Ptr (Ptr CSymmetry) -> CUInt -> Ptr CUInt -> CUInt -> IO CInt

foreign import capi unsafe "lattice_symmetries/lattice_symmetries.h &ls_destroy_symmetry"
  ls_destroy_symmetry :: FunPtr (Ptr CSymmetry -> IO ())

foreign import ccall unsafe "lattice_symmetries/lattice_symmetries.h ls_create_group"
  ls_create_group :: Ptr (Ptr CGroup) -> CUInt -> Ptr (Ptr CSymmetry) -> IO CInt

foreign import capi unsafe "lattice_symmetries/lattice_symmetries.h &ls_destroy_group"
  ls_destroy_group :: FunPtr (Ptr CGroup -> IO ())

foreign import capi unsafe "lattice_symmetries/lattice_symmetries.h ls_get_group_size"
  ls_get_group_size :: Ptr CGroup -> IO CUInt

foreign import ccall unsafe "lattice_symmetries/lattice_symmetries.h ls_create_spin_basis"
  ls_create_spin_basis :: Ptr (Ptr CSpinBasis) -> Ptr CGroup -> CUInt -> CInt -> CInt -> IO CInt

foreign import ccall unsafe "lattice_symmetries/lattice_symmetries.h &ls_destroy_spin_basis"
  ls_destroy_spin_basis :: FunPtr (Ptr CSpinBasis -> IO ())

foreign import ccall safe "lattice_symmetries/lattice_symmetries.h ls_build"
  ls_build :: Ptr CSpinBasis -> IO CInt

foreign import ccall safe "lattice_symmetries/lattice_symmetries.h ls_build_unsafe"
  ls_build_unsafe :: Ptr CSpinBasis -> Word64 -> Ptr Word64 -> IO CInt

foreign import ccall unsafe "lattice_symmetries/lattice_symmetries.h ls_get_number_states"
  ls_get_number_states :: Ptr CSpinBasis -> Ptr Word64 -> IO CInt

foreign import capi unsafe "lattice_symmetries/lattice_symmetries.h ls_flat_spin_basis_state_info"
  ls_flat_spin_basis_state_info :: Ptr CFlatSpinBasis -> Word64 -> Ptr () -> Ptr () -> Ptr (Complex Double) -> Ptr Double -> IO ()

-- foreign import ccall unsafe "lattice_symmetries/lattice_symmetries.h ls_get_states"
--   ls_get_states :: Ptr (Ptr ()) -> Ptr () -> IO CInt
--
-- foreign import ccall unsafe "lattice_symmetries/lattice_symmetries.h ls_states_get_data"
--   ls_states_get_data :: Ptr () -> Ptr Word64
--
-- foreign import ccall unsafe "lattice_symmetries/lattice_symmetries.h ls_states_get_size"
--   ls_states_get_size :: Ptr () -> Word64
--
-- foreign import ccall unsafe "lattice_symmetries/lattice_symmetries.h ls_destroy_states"
--   ls_destroy_states :: Ptr () -> IO ()

--
-- foreign import ccall unsafe "ls_get_sector" ls_get_sector :: Ptr () -> IO CUInt
--
-- foreign import ccall unsafe "ls_get_phase" ls_get_phase :: Ptr () -> IO CDouble
--
-- foreign import ccall unsafe "ls_get_periodicity" ls_get_periodicity :: Ptr () -> IO CUInt
