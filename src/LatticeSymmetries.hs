{-# LANGUAGE CApiFFI #-}

module LatticeSymmetries where

import Data.Aeson
import Data.Aeson.Types (typeMismatch)
import Data.Complex
import qualified Data.List.NonEmpty as NonEmpty
import Data.Scientific (toRealFloat)
import qualified Data.Vector.Unboxed as U
import Data.Yaml (decodeFileWithWarnings)
import Foreign.C.String (CString, peekCString)
import Foreign.C.Types (CInt (..), CUInt (..))
import Foreign.ForeignPtr
import Foreign.Marshal.Alloc (alloca)
import Foreign.Marshal.Array (withArrayLen)
import Foreign.Ptr (FunPtr, Ptr)
import Foreign.Storable (Storable (..))
import qualified GHC.ForeignPtr as GHC (Finalizers (..), ForeignPtr (..), ForeignPtrContents (..))
import qualified GHC.IORef as GHC (atomicSwapIORef)
import qualified GHC.Ptr as GHC (Ptr (..))
import UnliftIO.Exception (bracket, impureThrow, throwIO)

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

newtype Symmetry = Symmetry (ForeignPtr CSymmetry)

newtype SymmetryGroup = SymmetryGroup (ForeignPtr CGroup)

newtype SpinBasis = SpinBasis (ForeignPtr CSpinBasis)

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
