{-# LANGUAGE RankNTypes #-}

module ForeignLibrary () where

import Control.Exception.Safe (handleAny, handleAnyDeep)
import qualified Data.Aeson
import Data.List.Split (chunksOf)
import qualified Data.Text as Text
import Foreign.C.String (CString)
import Foreign.C.Types (CBool (..), CInt (..), CUInt (..))
import Foreign.Marshal.Alloc (alloca, free)
import Foreign.Marshal.Array (peekArray)
import Foreign.Marshal.Utils (fromBool, new)
import Foreign.Ptr (FunPtr, Ptr, nullPtr)
import Foreign.StablePtr
import Foreign.Storable (Storable (..))
import LatticeSymmetries
import LatticeSymmetries.Algebra
import LatticeSymmetries.Basis
import LatticeSymmetries.BitString
import LatticeSymmetries.ComplexRational (ComplexRational, fromComplexDouble)
import LatticeSymmetries.FFI
import LatticeSymmetries.Generator
import LatticeSymmetries.Operator
import LatticeSymmetries.Parser
import LatticeSymmetries.Utils
import Prettyprinter (Pretty (..))
import qualified Prettyprinter as Pretty
import Prettyprinter.Render.Text (renderStrict)
import Type.Reflection
import Prelude hiding (state)

foreign export ccall "ls_hs_hdf5_create_dataset_u64"
  ls_hs_hdf5_create_dataset_u64 :: CString -> CString -> CUInt -> Ptr Word64 -> IO ()

foreign export ccall "ls_hs_hdf5_create_dataset_f32"
  ls_hs_hdf5_create_dataset_f32 :: CString -> CString -> CUInt -> Ptr Word64 -> IO ()

foreign export ccall "ls_hs_hdf5_create_dataset_f64"
  ls_hs_hdf5_create_dataset_f64 :: CString -> CString -> CUInt -> Ptr Word64 -> IO ()

foreign export ccall "ls_hs_hdf5_create_dataset_c64"
  ls_hs_hdf5_create_dataset_c64 :: CString -> CString -> CUInt -> Ptr Word64 -> IO ()

foreign export ccall "ls_hs_hdf5_create_dataset_c128"
  ls_hs_hdf5_create_dataset_c128 :: CString -> CString -> CUInt -> Ptr Word64 -> IO ()

foreign export ccall "ls_hs_hdf5_write_chunk_u64"
  ls_hs_hdf5_write_chunk_u64 :: CString -> CString -> CUInt -> Ptr Word64 -> Ptr Word64 -> Ptr Word64 -> IO ()

foreign export ccall "ls_hs_hdf5_write_chunk_f64"
  ls_hs_hdf5_write_chunk_f64 :: CString -> CString -> CUInt -> Ptr Word64 -> Ptr Word64 -> Ptr Double -> IO ()

foreign export ccall "ls_hs_hdf5_read_chunk_u64"
  ls_hs_hdf5_read_chunk_u64 :: CString -> CString -> CUInt -> Ptr Word64 -> Ptr Word64 -> Ptr Word64 -> IO ()

foreign export ccall "ls_hs_hdf5_read_chunk_f64"
  ls_hs_hdf5_read_chunk_f64 :: CString -> CString -> CUInt -> Ptr Word64 -> Ptr Word64 -> Ptr Double -> IO ()

foreign export ccall "ls_hs_hdf5_get_dataset_rank"
  ls_hs_hdf5_get_dataset_rank :: CString -> CString -> IO CUInt

foreign export ccall "ls_hs_hdf5_get_dataset_shape"
  ls_hs_hdf5_get_dataset_shape :: CString -> CString -> Ptr Word64 -> IO ()

-- foreign export ccall "ls_hs_basis_and_hamiltonian_from_yaml"
--   ls_hs_basis_and_hamiltonian_from_yaml :: CString -> Ptr Types.SpinBasisWrapper -> Ptr Types.OperatorWrapper -> IO ()

-- foreign export ccall "ls_hs_destroy_spin_basis"
--   ls_hs_destroy_spin_basis :: Ptr Types.SpinBasisWrapper -> IO ()

-- foreign export ccall "ls_hs_destroy_operator"
--   ls_hs_destroy_operator :: Ptr Types.OperatorWrapper -> IO ()

-- foreign export ccall "ls_hs_create_basis"
--   ls_hs_create_basis :: Cparticle_type -> CInt -> CInt -> CInt -> IO (Ptr Cbasis)

foreign export ccall "ls_hs_clone_basis"
  ls_hs_clone_basis :: Ptr Cbasis -> IO (Ptr Cbasis)

ls_hs_clone_basis ptr = do
  _ <- basisIncRefCount ptr
  pure ptr

-- foreign export ccall "ls_hs_create_spin_basis_from_json"
--   ls_hs_create_spin_basis_from_json :: CString -> IO (Ptr Cbasis)

-- foreign export ccall "ls_hs_create_spinless_fermion_basis_from_json"
--   ls_hs_create_spinless_fermion_basis_from_json :: CString -> IO (Ptr Cbasis)

-- foreign export ccall "ls_hs_create_spinful_fermion_basis_from_json"
--   ls_hs_create_spinful_fermion_basis_from_json :: CString -> IO (Ptr Cbasis)

-- foreign export ccall "ls_hs_create_spin_basis_from_yaml"
--   ls_hs_create_spin_basis_from_yaml :: CString -> IO (Ptr Cbasis)

foreign export ccall "ls_hs_basis_to_json"
  ls_hs_basis_to_json :: Ptr Cbasis -> IO CString

ls_hs_basis_to_json cBasis =
  withReconstructedBasis cBasis $ \basis ->
    newCString $ toStrict (Data.Aeson.encode basis)

foreign export ccall "ls_hs_basis_from_json"
  ls_hs_basis_from_json :: CString -> IO (Ptr Cbasis)

ls_hs_basis_from_json cStr = handleAny (propagateErrorToC nullPtr) $ do
  (basis :: SomeBasis) <- decodeCString cStr
  foldSomeBasis borrowCbasis basis

foreign export ccall "ls_hs_destroy_string"
  ls_hs_destroy_string :: CString -> IO ()

ls_hs_destroy_string = free

-- foreign export ccall "ls_hs_spin_chain_10_basis"
--   ls_hs_spin_chain_10_basis :: IO (Ptr Cbasis)

-- foreign export ccall "ls_hs_spin_kagome_12_basis"
--   ls_hs_spin_kagome_12_basis :: IO (Ptr Cbasis)

-- foreign export ccall "ls_hs_spin_kagome_16_basis"
--   ls_hs_spin_kagome_16_basis :: IO (Ptr Cbasis)

-- foreign export ccall "ls_hs_spin_square_4x4_basis"
--   ls_hs_spin_square_4x4_basis :: IO (Ptr Cbasis)

foreign export ccall "ls_hs_min_state_estimate"
  ls_hs_min_state_estimate :: Ptr Cbasis -> IO Word64

ls_hs_min_state_estimate p =
  withReconstructedBasis p $ \basis ->
    let (BasisState n (BitString x)) = minStateEstimate (basisHeader basis)
     in if n > 64
          then throwC 0 "minimal state is not representable as a 64-bit integer"
          else pure $ fromIntegral x

foreign export ccall "ls_hs_max_state_estimate"
  ls_hs_max_state_estimate :: Ptr Cbasis -> IO Word64

ls_hs_max_state_estimate p =
  withReconstructedBasis p $ \basis ->
    let (BasisState n (BitString x)) = maxStateEstimate (basisHeader basis)
     in if n > 64
          then throwC (-1) "maximal state is not representable as a 64-bit integer"
          else pure $ fromIntegral x

ls_hs_basis_has_fixed_hamming_weight :: Ptr Cbasis -> IO CBool
ls_hs_basis_has_fixed_hamming_weight basis =
  fromBool
    <$> withReconstructedBasis basis (pure . hasFixedHammingWeight . basisHeader)

ls_hs_basis_build :: Ptr Cbasis -> IO ()
ls_hs_basis_build basis = withReconstructedBasis basis basisBuild

foreign export ccall "ls_hs_basis_build"
  ls_hs_basis_build :: Ptr Cbasis -> IO ()

foreign export ccall "ls_hs_basis_has_fixed_hamming_weight"
  ls_hs_basis_has_fixed_hamming_weight :: Ptr Cbasis -> IO CBool

foreign export ccall "ls_hs_basis_number_bits"
  ls_hs_basis_number_bits :: Ptr Cbasis -> IO CInt

ls_hs_basis_number_bits basisPtr =
  fromIntegral
    <$> withReconstructedBasis basisPtr (pure . getNumberBits . basisHeader)

foreign export ccall "ls_hs_basis_number_words"
  ls_hs_basis_number_words :: Ptr Cbasis -> IO CInt

ls_hs_basis_number_words basisPtr =
  fromIntegral
    <$> withReconstructedBasis basisPtr (pure . getNumberWords . basisHeader)

foreign export ccall "ls_hs_basis_state_to_string"
  ls_hs_basis_state_to_string :: Ptr Cbasis -> Ptr Word64 -> IO CString

ls_hs_basis_state_to_string basisPtr statePtr =
  withReconstructedBasis basisPtr $ \(basis :: Basis t) -> do
    let numberBits = getNumberBits (basisHeader basis)
        numberWords = getNumberWords (basisHeader basis)
    state <- BasisState @t numberBits <$> readBitString numberWords statePtr
    newCString . encodeUtf8 . toPrettyText $ state

newCexpr :: SomeExpr -> IO (Ptr Cexpr)
newCexpr !expr = do
  cExpr <- (new . Cexpr 1 . castStablePtrToPtr) =<< newStablePtr expr
  pure cExpr

withCexpr :: Ptr Cexpr -> (SomeExpr -> IO a) -> IO a
withCexpr !p f =
  f =<< deRefStablePtr =<< exprPeekPayload p

withCexpr2 :: Ptr Cexpr -> Ptr Cexpr -> (SomeExpr -> SomeExpr -> a) -> IO a
withCexpr2 !p1 !p2 f =
  withCexpr p1 $ \x1 ->
    withCexpr p2 $ \x2 ->
      pure (f x1 x2)

forCexpr :: (SomeExpr -> IO a) -> Ptr Cexpr -> IO a
forCexpr f p = withCexpr p f

foreign export ccall "ls_hs_create_expr"
  ls_hs_create_expr :: CString -> IO (Ptr Cexpr)

ls_hs_create_expr cStr =
  handleAnyDeep (propagateErrorToC nullPtr) $ do
    s <- peekUtf8 cStr
    let create :: forall t. ParticleTag t -> IO (Ptr Cexpr)
        create tag = newCexpr . SomeExpr tag $ mkExpr tag s
        result :: IO (Ptr Cexpr)
        result
          | isJust $ Text.find (\c -> c == 'σ' || c == 'S') s = create SpinTag
          | isJust $ Text.find (\c -> c == '↑' || c == '↓') s = create SpinfulFermionTag
          | otherwise = create SpinlessFermionTag
    result

foreign export ccall "ls_hs_destroy_expr"
  ls_hs_destroy_expr :: Ptr Cexpr -> IO ()

ls_hs_destroy_expr p = do
  refcount <- exprDecRefCount p
  when (refcount == 0) $ do
    freeStablePtr =<< exprPeekPayload p
    free p

foreign export ccall "ls_hs_expr_to_string"
  ls_hs_expr_to_string :: Ptr Cexpr -> IO CString

ls_hs_expr_to_string p =
  handleAnyDeep (propagateErrorToC nullPtr) $
    withCexpr p $
      newCString . encodeUtf8 . toPrettyText

foreign export ccall "ls_hs_expr_plus"
  ls_hs_expr_plus :: Ptr Cexpr -> Ptr Cexpr -> IO (Ptr Cexpr)

ls_hs_expr_plus a b = handleAnyDeep (propagateErrorToC nullPtr) $ withCexpr2 a b (+) >>= newCexpr

foreign export ccall "ls_hs_expr_minus"
  ls_hs_expr_minus :: Ptr Cexpr -> Ptr Cexpr -> IO (Ptr Cexpr)

ls_hs_expr_minus a b = handleAnyDeep (propagateErrorToC nullPtr) $ withCexpr2 a b (-) >>= newCexpr

foreign export ccall "ls_hs_expr_times"
  ls_hs_expr_times :: Ptr Cexpr -> Ptr Cexpr -> IO (Ptr Cexpr)

ls_hs_expr_times a b = handleAnyDeep (propagateErrorToC nullPtr) $ withCexpr2 a b (*) >>= newCexpr

foreign export ccall "ls_hs_expr_scale"
  ls_hs_expr_scale :: Ptr Cscalar -> Ptr Cexpr -> IO (Ptr Cexpr)

ls_hs_expr_scale c_z c_a =
  handleAnyDeep (propagateErrorToC nullPtr) $
    withCexpr c_a $ \a -> do
      z <- fromComplexDouble <$> peek c_z
      newCexpr $ scale (z :: ComplexRational) a

foreign export ccall "ls_hs_replace_indices"
  ls_hs_replace_indices :: Ptr Cexpr -> FunPtr Creplace_index -> IO (Ptr Cexpr)

ls_hs_replace_indices exprPtr fPtr =
  handleAnyDeep (propagateErrorToC nullPtr) $
    withCexpr exprPtr $ \expr -> do
      let f :: Int -> Int -> IO (Int, Int)
          f !s !i =
            alloca $ \spinPtr ->
              alloca $ \sitePtr -> do
                mkCreplace_index fPtr (fromIntegral s) (fromIntegral i) spinPtr sitePtr
                (,)
                  <$> (fromIntegral <$> peek spinPtr)
                  <*> (fromIntegral <$> peek sitePtr)
      newCexpr
        =<< case expr of
          SomeExpr SpinTag terms ->
            SomeExpr SpinTag <$> mapIndicesM (\i -> snd <$> f 0 i) terms
          SomeExpr SpinlessFermionTag terms ->
            SomeExpr SpinlessFermionTag <$> mapIndicesM (\i -> snd <$> f 0 i) terms
          SomeExpr SpinfulFermionTag terms ->
            let f' (s, i) = do
                  (s', i') <- f (fromEnum s) i
                  pure (toEnum s', i')
             in SomeExpr SpinfulFermionTag <$> mapIndicesM f' terms

foreign export ccall "ls_hs_expr_to_json"
  ls_hs_expr_to_json :: Ptr Cexpr -> IO CString

ls_hs_expr_to_json cExpr =
  withCexpr cExpr $ \expr -> do
    newCString $ toStrict (Data.Aeson.encode expr)

foreign export ccall "ls_hs_expr_from_json"
  ls_hs_expr_from_json :: CString -> IO (Ptr Cexpr)

ls_hs_expr_from_json cStr = handleAny (propagateErrorToC nullPtr) $ do
  !expr <- decodeCString cStr
  newCexpr expr

foreign export ccall "ls_hs_create_operator"
  ls_hs_create_operator :: Ptr Cbasis -> Ptr Cexpr -> IO (Ptr Coperator)

ls_hs_create_operator :: Ptr Cbasis -> Ptr Cexpr -> IO (Ptr Coperator)
ls_hs_create_operator basisPtr exprPtr =
  handleAny (propagateErrorToC nullPtr) $
    withReconstructedBasis basisPtr $ \basis ->
      withCexpr exprPtr $ \someExpr ->
        withSomeExpr @IsBasis someExpr $ \expr ->
          case matchParticleType2 basis expr of
            Just HRefl -> borrowCoperator (mkOperator basis expr)
            Nothing -> error "basis and expression have different particle types"

operatorBinaryFunction ::
  (forall t. IsBasis t => Operator t -> Operator t -> Operator t) ->
  Ptr Coperator ->
  Ptr Coperator ->
  IO (Ptr Coperator)
operatorBinaryFunction f aPtr bPtr =
  handleAny (propagateErrorToC nullPtr) $
    withReconstructedOperator aPtr $ \a ->
      withReconstructedOperator bPtr $ \b ->
        case matchParticleType2 a b of
          Just HRefl -> borrowCoperator (f a b)
          Nothing -> error "operators have different particle types"

ls_hs_operator_plus :: Ptr Coperator -> Ptr Coperator -> IO (Ptr Coperator)
ls_hs_operator_plus = operatorBinaryFunction (+)

ls_hs_operator_minus :: Ptr Coperator -> Ptr Coperator -> IO (Ptr Coperator)
ls_hs_operator_minus = operatorBinaryFunction (-)

ls_hs_operator_times :: Ptr Coperator -> Ptr Coperator -> IO (Ptr Coperator)
ls_hs_operator_times = operatorBinaryFunction (*)

ls_hs_operator_scale :: Ptr Cscalar -> Ptr Coperator -> IO (Ptr Coperator)
ls_hs_operator_scale cPtr opPtr =
  handleAny (propagateErrorToC nullPtr) $
    withReconstructedOperator opPtr $ \op -> do
      (c :: ComplexRational) <- fromComplexDouble <$> peek cPtr
      borrowCoperator $ scale c op

foreign export ccall "ls_hs_operator_plus"
  ls_hs_operator_plus :: Ptr Coperator -> Ptr Coperator -> IO (Ptr Coperator)

foreign export ccall "ls_hs_operator_minus"
  ls_hs_operator_minus :: Ptr Coperator -> Ptr Coperator -> IO (Ptr Coperator)

foreign export ccall "ls_hs_operator_times"
  ls_hs_operator_times :: Ptr Coperator -> Ptr Coperator -> IO (Ptr Coperator)

foreign export ccall "ls_hs_operator_scale"
  ls_hs_operator_scale :: Ptr Cscalar -> Ptr Coperator -> IO (Ptr Coperator)

foreign export ccall "ls_hs_operator_hermitian_conjugate"
  ls_hs_operator_hermitian_conjugate :: Ptr Coperator -> IO (Ptr Coperator)

ls_hs_operator_hermitian_conjugate opPtr =
  withReconstructedOperator opPtr $ \(Operator (OperatorHeader basis terms) _) ->
    borrowCoperator $
      operatorFromHeader (OperatorHeader basis (conjugateExpr terms))

foreign export ccall "ls_hs_operator_is_hermitian"
  ls_hs_operator_is_hermitian :: Ptr Coperator -> IO CBool

ls_hs_operator_is_hermitian opPtr =
  fromBool <$> withReconstructedOperator opPtr (pure . isHermitianExpr . opTerms . opHeader)

foreign export ccall "ls_hs_operator_is_identity"
  ls_hs_operator_is_identity :: Ptr Coperator -> IO CBool

ls_hs_operator_is_identity opPtr =
  fromBool <$> withReconstructedOperator opPtr (pure . isIdentityExpr . opTerms . opHeader)

foreign export ccall "ls_hs_operator_max_number_off_diag"
  ls_hs_operator_max_number_off_diag :: Ptr Coperator -> IO CInt

ls_hs_operator_max_number_off_diag opPtr =
  fromIntegral <$> withReconstructedOperator opPtr (pure . maxNumberOffDiag)

foreign export ccall "ls_hs_load_hamiltonian_from_yaml"
  ls_hs_load_hamiltonian_from_yaml :: CString -> IO (Ptr Coperator)

ls_hs_operator_pretty_terms :: Ptr Coperator -> IO CString
ls_hs_operator_pretty_terms p =
  withReconstructedOperator p $ \op ->
    newCString
      . encodeUtf8
      . renderStrict
      . Pretty.layoutPretty (Pretty.LayoutOptions Pretty.Unbounded)
      . pretty
      $ opTerms (opHeader op)

foreign export ccall "ls_hs_operator_pretty_terms"
  ls_hs_operator_pretty_terms :: Ptr Coperator -> IO CString

-- foreign export ccall "ls_hs_fatal_error"
--   ls_hs_fatal_error :: CString -> CString -> IO ()
