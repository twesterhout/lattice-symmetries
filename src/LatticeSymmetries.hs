{-# LANGUAGE CApiFFI #-}

-- {-# LANGUAGE QuasiQuotes #-}
-- {-# LANGUAGE TemplateHaskell #-}

module LatticeSymmetries where

import Control.Exception.Safe (MonadThrow, bracket, impureThrow, throwIO)
import Data.Complex
import qualified Data.HDF5 as H5
import qualified Data.HDF5.Types as H5
import Data.Vector.Storable (Vector)
import qualified Data.Vector.Storable as V
-- import Data.Yaml.Aeson
import Foreign.C.String (CString, peekCString)
import Foreign.C.Types (CBool (..), CInt (..), CUInt (..), CUShort (..))
import Foreign.ForeignPtr
import Foreign.Marshal.Alloc (alloca)
import Foreign.Marshal.Array (peekArray, pokeArray, withArrayLen)
import Foreign.Ptr (FunPtr, Ptr, nullPtr)
import Foreign.StablePtr
import Foreign.Storable (Storable (..))
-- import LatticeSymmetries.IO
-- import LatticeSymmetries.Types
import qualified System.IO.Unsafe

-- import UnliftIO.Exception (bracket, impureThrow, throwIO)

  {-
-- | Retrieve textual representation of an error
getErrorMessage ::
  -- | Error code returned by @lattice_symmetries@ library
  Int ->
  -- | Explanation of the error
  IO Text
getErrorMessage c = bracket (ls_error_to_string (fromIntegral c)) ls_destroy_string $ \s ->
  toText <$> peekCString s
-}

  {-
-- | Check the status code returned by @lattice_symmetries@ library. If it
-- indicates an error, 'LatticeSymmetriesException' is thrown.
checkStatus :: (MonadIO m, Integral a) => a -> m ()
checkStatus c
  | c == 0 = return ()
  | otherwise = do
    print =<< liftIO (getErrorMessage c')
    liftIO $ throwIO . LatticeSymmetriesException c' =<< getErrorMessage c'
  where
    c' = fromIntegral c
-}

-- loadRawBasis :: MonadIO m => Text -> m (Ptr Cspin_basis)
-- loadRawBasis path = do
--   r <- liftIO $ decodeFileWithWarnings (toString path)
--   case r of
--     Left e -> throwIO e
--     Right (warnings, (WrappedBasisSpec basisSpec)) -> do
--       mapM_ print warnings
--       toRawBasis basisSpec

-- ls_hs_internal_load_basis_from_yaml :: CString -> IO (Ptr Cspin_basis)
-- ls_hs_internal_load_basis_from_yaml c_str = loadRawBasis . toText =<< peekCString c_str
--
-- foreign export ccall "ls_hs_internal_load_basis_from_yaml"
--   ls_hs_internal_load_basis_from_yaml :: CString -> IO (Ptr Cspin_basis)

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

createDataset ::
  forall a.
  (H5.KnownDatatype a) =>
  CString ->
  CString ->
  CUInt ->
  Ptr Word64 ->
  IO ()
createDataset _fileName _datasetName _dim shapePtr = do
  fileName <- fromString <$> peekCString _fileName
  datasetName <- fromString <$> peekCString _datasetName
  shape <- fmap fromIntegral <$> peekArray (fromIntegral _dim) shapePtr
  H5.withFile fileName H5.WriteAppend $ \file -> do
    exists <- H5.exists file datasetName
    when exists $ H5.delete file datasetName
    void . join $ H5.createEmptyDataset file datasetName <$> H5.ofType @a <*> H5.ofShape shape

ls_hs_hdf5_create_dataset_u64 :: CString -> CString -> CUInt -> Ptr Word64 -> IO ()
ls_hs_hdf5_create_dataset_u64 = createDataset @Word64

ls_hs_hdf5_create_dataset_f32 :: CString -> CString -> CUInt -> Ptr Word64 -> IO ()
ls_hs_hdf5_create_dataset_f32 = createDataset @Double

ls_hs_hdf5_create_dataset_f64 :: CString -> CString -> CUInt -> Ptr Word64 -> IO ()
ls_hs_hdf5_create_dataset_f64 = createDataset @Double

ls_hs_hdf5_create_dataset_c64 :: CString -> CString -> CUInt -> Ptr Word64 -> IO ()
ls_hs_hdf5_create_dataset_c64 = createDataset @(Complex Float)

ls_hs_hdf5_create_dataset_c128 :: CString -> CString -> CUInt -> Ptr Word64 -> IO ()
ls_hs_hdf5_create_dataset_c128 = createDataset @(Complex Double)

-- foreign export ccall "ls_hs_hdf5_create_dataset_u64"
--   ls_hs_hdf5_create_dataset_u64 :: CString -> CString -> CUInt -> Ptr Word64 -> IO ()

-- foreign export ccall "ls_hs_hdf5_create_dataset_f32"
--   ls_hs_hdf5_create_dataset_f32 :: CString -> CString -> CUInt -> Ptr Word64 -> IO ()

-- foreign export ccall "ls_hs_hdf5_create_dataset_f64"
--   ls_hs_hdf5_create_dataset_f64 :: CString -> CString -> CUInt -> Ptr Word64 -> IO ()

-- foreign export ccall "ls_hs_hdf5_create_dataset_c64"
--   ls_hs_hdf5_create_dataset_c64 :: CString -> CString -> CUInt -> Ptr Word64 -> IO ()

-- foreign export ccall "ls_hs_hdf5_create_dataset_c128"
--   ls_hs_hdf5_create_dataset_c128 :: CString -> CString -> CUInt -> Ptr Word64 -> IO ()

writeDatasetChunk ::
  forall a.
  (Storable a, H5.KnownDatatype a) =>
  CString ->
  CString ->
  CUInt ->
  Ptr Word64 ->
  Ptr Word64 ->
  Ptr a ->
  IO ()
writeDatasetChunk _fileName _datasetName _dim _offsetPtr _shapePtr dataPtr = do
  fileName <- fromString <$> peekCString _fileName
  datasetName <- fromString <$> peekCString _datasetName
  let rank = fromIntegral _dim
  offset <- fromList . fmap fromIntegral <$> peekArray rank _offsetPtr
  shape <- fromList . fmap fromIntegral <$> peekArray rank _shapePtr
  let hyperslab = H5.Hyperslab offset (V.replicate rank 1) shape (V.replicate rank 1)
  fp <- newForeignPtr_ dataPtr
  let size = V.product shape
      buffer = V.unsafeFromForeignPtr0 fp size
  H5.withFile fileName H5.WriteAppend $ \file ->
    H5.open file datasetName
      >>= return
        . H5.sliceWithHyperslab hyperslab
      >>= H5.writeSelected buffer

ls_hs_hdf5_write_chunk_u64 ::
  CString -> CString -> CUInt -> Ptr Word64 -> Ptr Word64 -> Ptr Word64 -> IO ()
ls_hs_hdf5_write_chunk_u64 = writeDatasetChunk @Word64

ls_hs_hdf5_write_chunk_f32 ::
  CString -> CString -> CUInt -> Ptr Word64 -> Ptr Word64 -> Ptr Float -> IO ()
ls_hs_hdf5_write_chunk_f32 = writeDatasetChunk @Float

ls_hs_hdf5_write_chunk_f64 ::
  CString -> CString -> CUInt -> Ptr Word64 -> Ptr Word64 -> Ptr Double -> IO ()
ls_hs_hdf5_write_chunk_f64 = writeDatasetChunk @Double

ls_hs_hdf5_write_chunk_c64 ::
  CString -> CString -> CUInt -> Ptr Word64 -> Ptr Word64 -> Ptr (Complex Float) -> IO ()
ls_hs_hdf5_write_chunk_c64 = writeDatasetChunk @(Complex Float)

ls_hs_hdf5_write_chunk_c128 ::
  CString -> CString -> CUInt -> Ptr Word64 -> Ptr Word64 -> Ptr (Complex Double) -> IO ()
ls_hs_hdf5_write_chunk_c128 = writeDatasetChunk @(Complex Double)

-- foreign export ccall "ls_hs_hdf5_write_chunk_u64"
--   ls_hs_hdf5_write_chunk_u64 :: CString -> CString -> CUInt -> Ptr Word64 -> Ptr Word64 -> Ptr Word64 -> IO ()

-- foreign export ccall "ls_hs_hdf5_write_chunk_f64"
--   ls_hs_hdf5_write_chunk_f64 :: CString -> CString -> CUInt -> Ptr Word64 -> Ptr Word64 -> Ptr Double -> IO ()

readDatasetChunk ::
  forall a.
  (Storable a, H5.KnownDatatype a) =>
  CString ->
  CString ->
  CUInt ->
  Ptr Word64 ->
  Ptr Word64 ->
  Ptr a ->
  IO ()
readDatasetChunk _fileName _datasetName _dim _offsetPtr _shapePtr dataPtr = do
  fileName <- fromString <$> peekCString _fileName
  datasetName <- fromString <$> peekCString _datasetName
  let rank = fromIntegral _dim
  offset <- fromList . fmap fromIntegral <$> peekArray rank _offsetPtr
  shape <- fromList . fmap fromIntegral <$> peekArray rank _shapePtr
  let hyperslab = H5.Hyperslab offset (V.replicate rank 1) shape (V.replicate rank 1)
  fp <- newForeignPtr_ dataPtr
  let size = V.product shape
      buffer = V.unsafeFromForeignPtr0 fp size
  H5.withFile fileName H5.WriteAppend $ \file ->
    H5.open file datasetName
      >>= return
        . H5.sliceWithHyperslab hyperslab
      >>= H5.readSelectedInplace buffer

ls_hs_hdf5_read_chunk_u64 ::
  CString -> CString -> CUInt -> Ptr Word64 -> Ptr Word64 -> Ptr Word64 -> IO ()
ls_hs_hdf5_read_chunk_u64 = readDatasetChunk @Word64

ls_hs_hdf5_read_chunk_f32 ::
  CString -> CString -> CUInt -> Ptr Word64 -> Ptr Word64 -> Ptr Float -> IO ()
ls_hs_hdf5_read_chunk_f32 = readDatasetChunk @Float

ls_hs_hdf5_read_chunk_f64 ::
  CString -> CString -> CUInt -> Ptr Word64 -> Ptr Word64 -> Ptr Double -> IO ()
ls_hs_hdf5_read_chunk_f64 = readDatasetChunk @Double

ls_hs_hdf5_read_chunk_c64 ::
  CString -> CString -> CUInt -> Ptr Word64 -> Ptr Word64 -> Ptr (Complex Float) -> IO ()
ls_hs_hdf5_read_chunk_c64 = readDatasetChunk @(Complex Float)

ls_hs_hdf5_read_chunk_c128 ::
  CString -> CString -> CUInt -> Ptr Word64 -> Ptr Word64 -> Ptr (Complex Double) -> IO ()
ls_hs_hdf5_read_chunk_c128 = readDatasetChunk @(Complex Double)

-- foreign export ccall "ls_hs_hdf5_read_chunk_u64"
--   ls_hs_hdf5_read_chunk_u64 :: CString -> CString -> CUInt -> Ptr Word64 -> Ptr Word64 -> Ptr Word64 -> IO ()

-- foreign export ccall "ls_hs_hdf5_read_chunk_f64"
--   ls_hs_hdf5_read_chunk_f64 :: CString -> CString -> CUInt -> Ptr Word64 -> Ptr Word64 -> Ptr Double -> IO ()

ls_hs_hdf5_get_dataset_rank :: CString -> CString -> IO CUInt
ls_hs_hdf5_get_dataset_rank _fileName _datasetName = do
  fileName <- fromString <$> peekCString _fileName
  datasetName <- fromString <$> peekCString _datasetName
  rank <- H5.withFile fileName H5.ReadOnly $ \file ->
    H5.open file datasetName >>= return . H5.datasetRank
  pure (fromIntegral rank)

-- foreign export ccall "ls_hs_hdf5_get_dataset_rank"
--   ls_hs_hdf5_get_dataset_rank :: CString -> CString -> IO CUInt

ls_hs_hdf5_get_dataset_shape :: CString -> CString -> Ptr Word64 -> IO ()
ls_hs_hdf5_get_dataset_shape _fileName _datasetName sizePtr = do
  fileName <- fromString <$> peekCString _fileName
  datasetName <- fromString <$> peekCString _datasetName
  shape <- H5.withFile fileName H5.ReadOnly $ \file ->
    H5.open file datasetName >>= return . H5.datasetShape
  pokeArray sizePtr (fromIntegral <$> shape)

-- foreign export ccall "ls_hs_hdf5_get_dataset_shape"
--   ls_hs_hdf5_get_dataset_shape :: CString -> CString -> Ptr Word64 -> IO ()

-- writeDatasetChunk _fileName _datasetName offset size ptr = do
--   fileName <- fromString <$> peekCString _fileName
--   datasetName <- fromString <$> peekCString _datasetName
--   H5.withFile filename H5.AppendMode $ \file ->
--     H5.open

{-
loadBasisAndHamiltonianFromYAML :: Text -> IO (SpinBasis, Operator)
loadBasisAndHamiltonianFromYAML path = do
  r <- decodeFileWithWarnings (toString path)
  case r of
    Left e -> throwIO e
    Right (warnings, (ConfigSpec basisSpec operatorSpec)) -> do
      mapM_ print warnings
      let basis = mkBasis basisSpec
      return (basis, mkOperator basis operatorSpec)
-}

{-
ls_hs_basis_and_hamiltonian_from_yaml :: CString -> Ptr SpinBasisWrapper -> Ptr OperatorWrapper -> IO ()
ls_hs_basis_and_hamiltonian_from_yaml path basisPtr operatorPtr = do
  (basis, hamiltonian) <- loadBasisAndHamiltonianFromYAML . toText =<< peekCString path
  case basis of
    SpinBasis fp -> withForeignPtr fp $ \rawPtr ->
      poke basisPtr =<< SpinBasisWrapper rawPtr <$> newStablePtr basis
  case hamiltonian of
    Operator fp -> withForeignPtr fp $ \rawPtr ->
      poke operatorPtr =<< OperatorWrapper rawPtr <$> newStablePtr hamiltonian
-}

-- void ls_hs_basis_and_hamiltonian_from_yaml(char const *path,
--                                            ls_hs_spin_basis_v1 *basis,
--                                            ls_hs_operator_v1 *hamiltonian);
-- foreign export ccall "ls_hs_basis_and_hamiltonian_from_yaml"
--   ls_hs_basis_and_hamiltonian_from_yaml :: CString -> Ptr SpinBasisWrapper -> Ptr OperatorWrapper -> IO ()

{-
ls_hs_destroy_spin_basis :: Ptr SpinBasisWrapper -> IO ()
ls_hs_destroy_spin_basis ptr = do
  peek ptr >>= \(SpinBasisWrapper _ stablePtr) -> freeStablePtr stablePtr
  poke ptr $ SpinBasisWrapper nullPtr (castPtrToStablePtr nullPtr)
-}

-- void ls_hs_destroy_spin_basis(ls_hs_spin_basis_v1 *basis);
-- foreign export ccall "ls_hs_destroy_spin_basis"
--   ls_hs_destroy_spin_basis :: Ptr SpinBasisWrapper -> IO ()

{-
ls_hs_destroy_operator :: Ptr OperatorWrapper -> IO ()
ls_hs_destroy_operator ptr = do
  peek ptr >>= \(OperatorWrapper _ stablePtr) -> freeStablePtr stablePtr
  poke ptr $ OperatorWrapper nullPtr (castPtrToStablePtr nullPtr)
-}

-- void ls_hs_destroy_operator(ls_hs_operator_v1 *op);
-- foreign export ccall "ls_hs_destroy_operator"
--   ls_hs_destroy_operator :: Ptr OperatorWrapper -> IO ()

{-
toSymmetry :: SymmetrySpec -> Symmetry
toSymmetry (SymmetrySpec p s) = mkSymmetry p s

mkSymmetry ::
  -- | Permutation
  NonEmpty Int ->
  -- | Symmetry sector
  Int ->
  Symmetry
mkSymmetry !permutation !sector
  | any (< 0) permutation =
    error $
      "invalid permutation: " <> show permutation <> "; indices must be non-negative"
  | sector < 0 =
    error $
      "invalid sector: " <> show sector <> "; expected a non-negative number"
  | otherwise = System.IO.Unsafe.unsafePerformIO $ do
    ptr <- mkObject $ \ptrPtr ->
      withArrayLen (fromIntegral <$> toList permutation) $ \n permutationPtr ->
        ls_create_symmetry ptrPtr (fromIntegral n) permutationPtr (fromIntegral sector)
    Symmetry <$> newForeignPtr ls_destroy_symmetry ptr
-}

{-
mkSymmetryGroup ::
  -- | Symmetry generators
  [Symmetry] ->
  SymmetryGroup
mkSymmetryGroup symmetries = System.IO.Unsafe.unsafePerformIO $ do
  ptr <- liftIO . mkObject $ \ptrPtr ->
    withSymmetries symmetries $ \n symmetriesPtr ->
      ls_create_group ptrPtr (fromIntegral n) symmetriesPtr
  fmap SymmetryGroup . liftIO $ newForeignPtr ls_destroy_group ptr

mkRawBasis ::
  (MonadIO m, MonadThrow m) =>
  -- | Symmetry group
  SymmetryGroup ->
  -- | Number of spins
  Int ->
  -- | Hamming weight
  Maybe Int ->
  -- | Spin inversion
  Maybe Int ->
  m SpinBasis
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
  SpinBasis <$> liftIO (newForeignPtr ls_destroy_spin_basis ptr)
-}

{-
mkFlatBasis :: BasisSpec -> FlatSpinBasis
mkFlatBasis spec = System.IO.Unsafe.unsafePerformIO $ do
  flatBasisPtr <- alloca $ \basisPtrPtr ->
    withForeignPtr normalBasis $ \normalBasisPtr -> do
      checkStatus =<< ls_convert_to_flat_spin_basis basisPtrPtr normalBasisPtr
      peek basisPtrPtr
  FlatSpinBasis <$> newForeignPtr ls_destroy_flat_spin_basis flatBasisPtr
  where
    (SpinBasis normalBasis) = mkBasis spec

mkBasis :: BasisSpec -> SpinBasis
mkBasis (BasisSpec numberSpins hammingWeight spinInversion symmetries)
  | numberSpins > 64 =
    error $
      "invalid number_spins: " <> show numberSpins <> "; exact diagonalization is not feasible "
        <> "for systems larger than 64 spins"
  | otherwise =
    let symmetryGroup = mkSymmetryGroup (toSymmetry <$> symmetries)
     in System.IO.Unsafe.unsafePerformIO $
          mkRawBasis symmetryGroup numberSpins hammingWeight spinInversion
-}

{-
mkInteraction :: InteractionSpec -> Interaction
mkInteraction (InteractionSpec matrix sites) =
  System.IO.Unsafe.unsafePerformIO $
    mkInteraction' matrix sites

mkOperator :: SpinBasis -> OperatorSpec -> Operator
mkOperator basis (OperatorSpec _ interactions) =
  mkRawOperator basis $
    mkInteraction <$> interactions
-}

{-
mkRawOperator :: SpinBasis -> NonEmpty Interaction -> Operator
mkRawOperator (SpinBasis basis) terms = System.IO.Unsafe.unsafePerformIO $ do
  (code, ptr) <- liftIO $
    alloca $ \ptrPtr -> do
      c <- withInteractions (toList terms) $ \n interactionsPtr ->
        withForeignPtr basis $ \basisPtr ->
          ls_create_operator ptrPtr basisPtr (fromIntegral n) interactionsPtr
      if c == 0
        then (,) <$> pure c <*> peek ptrPtr
        else pure (c, nullPtr)
  checkStatus code
  fmap Operator . liftIO $ newForeignPtr ls_destroy_operator ptr

withInteractions :: [Interaction] -> (Int -> Ptr (Ptr Cinteraction) -> IO a) -> IO a
withInteractions xs func = withManyForeignPtr pointers func
  where
    pointers = (\(Interaction p) -> p) <$> xs
-}

{-
toMatrix ::
  (Monad m, Show r, Real r) =>
  -- | Expected dimension @n@ of the matrix
  Int ->
  -- | Square @matrix@ of dimension @n@
  [[Complex r]] ->
  -- | Row-major representation of @matrix@
  m (Vector (Complex Double))
toMatrix !dim !rows = do
  when (length rows /= dim) . impureThrow . SpinEDException $
    "invalid matrix: " <> show rows <> "; expected a square matrix of dimension " <> show dim
  when (any ((/= dim) . length) rows) . impureThrow . SpinEDException $
    "invalid matrix: " <> show rows <> "; expected a square matrix of dimension " <> show dim
  return . V.fromList . fmap (fmap (fromRational . toRational)) . concat $ rows

type CreateInteraction = Ptr (Ptr Cinteraction) -> Ptr (Complex Double) -> CUInt -> Ptr CUShort -> IO CInt

class MakeInteraction a where
  mkInteraction' :: (MonadIO m, MonadThrow m, Show r, Real r) => [[Complex r]] -> [a] -> m Interaction

instance MakeInteraction Int where
  mkInteraction' matrix sites =
    toMatrix 2 matrix >>= \matrix' ->
      unsafeMkInteraction ls_create_interaction1 (length sites) matrix' sites'
    where
      sites' = V.fromList sites

instance MakeInteraction (Int, Int) where
  mkInteraction' matrix sites =
    toMatrix 4 matrix >>= \matrix' ->
      unsafeMkInteraction ls_create_interaction2 (length sites) matrix' sites'
    where
      sites' = V.fromList $ sites >>= \(x₁, x₂) -> [x₁, x₂]

instance MakeInteraction (Int, Int, Int) where
  mkInteraction' matrix sites =
    toMatrix 8 matrix >>= \matrix' ->
      unsafeMkInteraction ls_create_interaction3 (length sites) matrix' sites'
    where
      sites' = V.fromList $ sites >>= \(x₁, x₂, x₃) -> [x₁, x₂, x₃]

instance MakeInteraction (Int, Int, Int, Int) where
  mkInteraction' matrix sites =
    toMatrix 16 matrix >>= \matrix' ->
      unsafeMkInteraction ls_create_interaction4 (length sites) matrix' sites'
    where
      sites' = V.fromList $ sites >>= \(x₁, x₂, x₃, x₄) -> [x₁, x₂, x₃, x₄]

instance MakeInteraction [Int] where
  mkInteraction' _ [] =
    throwIO . SpinEDException $
      "zero-point interactions (i.e. constant factors) are not supported"
  mkInteraction' matrix rows@(r : _) = case n of
    1 -> mkInteraction' matrix =<< mapM match1 rows
    2 -> mkInteraction' matrix =<< mapM match2 rows
    3 -> mkInteraction' matrix =<< mapM match3 rows
    4 -> mkInteraction' matrix =<< mapM match4 rows
    _ ->
      throwIO . SpinEDException $
        "currently only 1-, 2-, 3-, and 4-point interactions are supported, but received n="
          <> show n
    where
      n = length r
      match1 [x₁] = return x₁
      match1 _ = throwIO failure
      match2 [x₁, x₂] = return (x₁, x₂)
      match2 _ = throwIO failure
      match3 [x₁, x₂, x₃] = return (x₁, x₂, x₃)
      match3 _ = throwIO failure
      match4 [x₁, x₂, x₃, x₄] = return (x₁, x₂, x₃, x₄)
      match4 _ = throwIO failure
      failure =
        SpinEDException $
          "invalid sites: " <> show rows <> "; expected an array of length-" <> show n <> " tuples"
-}

{-
unsafeMkInteraction ::
  (MonadIO m, MonadThrow m) =>
  CreateInteraction ->
  Int ->
  Vector (Complex Double) ->
  Vector Int ->
  m Interaction
unsafeMkInteraction ffiCreate numberSites matrix sites = do
  when (V.any (< 0) sites) . throwIO . SpinEDException $
    "site indices must be all non-negative numbers"
  (code, ptr) <- liftIO $
    alloca $ \ptrPtr -> do
      c <- V.unsafeWith matrix $ \matrixPtr ->
        V.unsafeWith (V.map fromIntegral sites) $ \sitesPtr ->
          ffiCreate ptrPtr matrixPtr (fromIntegral numberSites) sitesPtr
      if c == 0
        then (,) <$> pure c <*> peek ptrPtr
        else pure (c, nullPtr)
  checkStatus code
  fmap Interaction . liftIO $ newForeignPtr ls_destroy_interaction ptr
-}

{-
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
-}

  {-
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
  ls_create_spin_basis :: Ptr (Ptr Cspin_basis) -> Ptr CGroup -> CUInt -> CInt -> CInt -> IO CInt

foreign import ccall unsafe "lattice_symmetries/lattice_symmetries.h &ls_destroy_spin_basis"
  ls_destroy_spin_basis :: FunPtr (Ptr Cspin_basis -> IO ())

foreign import ccall safe "lattice_symmetries/lattice_symmetries.h ls_build"
  ls_build :: Ptr Cspin_basis -> IO CInt

foreign import ccall safe "lattice_symmetries/lattice_symmetries.h ls_build_unsafe"
  ls_build_unsafe :: Ptr Cspin_basis -> Word64 -> Ptr Word64 -> IO CInt

foreign import ccall unsafe "lattice_symmetries/lattice_symmetries.h ls_get_number_states"
  ls_get_number_states :: Ptr Cspin_basis -> Ptr Word64 -> IO CInt

{- ORMOLU_DISABLE -}
foreign import capi unsafe "lattice_symmetries/lattice_symmetries.h ls_flat_spin_basis_state_info"
  ls_flat_spin_basis_state_info :: Ptr CFlatSpinBasis -> Word64 -> Ptr () -> Ptr () ->
                                   Ptr (Complex Double) -> Ptr Double -> IO ()
{- ORMOLU_ENABLE -}

foreign import capi unsafe "lattice_symmetries/lattice_symmetries.h ls_convert_to_flat_spin_basis"
  ls_convert_to_flat_spin_basis :: Ptr (Ptr CFlatSpinBasis) -> Ptr Cspin_basis -> IO CInt

foreign import capi unsafe "lattice_symmetries/lattice_symmetries.h &ls_destroy_flat_spin_basis"
  ls_destroy_flat_spin_basis :: FunPtr (Ptr CFlatSpinBasis -> IO ())

foreign import capi unsafe "lattice_symmetries/lattice_symmetries.h ls_create_interaction1"
  ls_create_interaction1 :: Ptr (Ptr Cinteraction) -> Ptr (Complex Double) -> CUInt -> Ptr CUShort -> IO CInt

foreign import capi unsafe "lattice_symmetries/lattice_symmetries.h ls_create_interaction2"
  ls_create_interaction2 :: Ptr (Ptr Cinteraction) -> Ptr (Complex Double) -> CUInt -> Ptr CUShort -> IO CInt

foreign import capi unsafe "lattice_symmetries/lattice_symmetries.h ls_create_interaction3"
  ls_create_interaction3 :: Ptr (Ptr Cinteraction) -> Ptr (Complex Double) -> CUInt -> Ptr CUShort -> IO CInt

foreign import capi unsafe "lattice_symmetries/lattice_symmetries.h ls_create_interaction4"
  ls_create_interaction4 :: Ptr (Ptr Cinteraction) -> Ptr (Complex Double) -> CUInt -> Ptr CUShort -> IO CInt

foreign import capi unsafe "lattice_symmetries/lattice_symmetries.h ls_interaction_is_real"
  ls_interaction_is_real :: Ptr () -> IO CBool

foreign import capi unsafe "lattice_symmetries/lattice_symmetries.h &ls_destroy_interaction"
  ls_destroy_interaction :: FunPtr (Ptr Cinteraction -> IO ())

foreign import ccall unsafe "lattice_symmetries/lattice_symmetries.h ls_create_operator"
  ls_create_operator :: Ptr (Ptr Coperator) -> Ptr Cspin_basis -> CUInt -> Ptr (Ptr Cinteraction) -> IO CInt

foreign import ccall unsafe "lattice_symmetries/lattice_symmetries.h &ls_destroy_operator"
  ls_destroy_operator :: FunPtr (Ptr Coperator -> IO ())
-}
