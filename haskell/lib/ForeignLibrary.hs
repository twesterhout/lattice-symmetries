{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE OverloadedLists #-}
{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE TemplateHaskell #-}
{-# OPTIONS_GHC -Wno-missed-extra-shared-lib #-}

module ForeignLibrary () where

import Data.Aeson qualified
import Data.Vector.Generic qualified as G
import Foreign.C.String (CString)
import Foreign.C.Types (CBool (..), CInt (..), CPtrdiff (..))
import Foreign.Marshal (alloca, free, fromBool, malloc, newArray)
import Foreign.Ptr (FunPtr, Ptr, nullPtr)
import Foreign.Storable (Storable (..))
import GHC.Exts (IsList (..))
import HeaderFileGeneration
import LatticeSymmetries.Algebra
import LatticeSymmetries.Basis
import LatticeSymmetries.BitString
import LatticeSymmetries.ComplexRational (ComplexRational, fromComplexDouble)
import LatticeSymmetries.Conversion
import LatticeSymmetries.Expr
import LatticeSymmetries.FFI
import LatticeSymmetries.Generator
import LatticeSymmetries.Group
import LatticeSymmetries.Operator
import LatticeSymmetries.Utils
import LatticeSymmetries.Yaml
import Type.Reflection
import Prelude hiding (state, toList)

ls_hs_destroy_string :: CString -> IO ()
ls_hs_destroy_string = free

newCencoded :: (HasCallStack, Data.Aeson.ToJSON a) => a -> IO CString
newCencoded = newCString . toStrict . Data.Aeson.encode

data SymmetryInfo = SymmetryInfo
  { periodicity :: !Int
  , phase :: !Double
  }
  deriving stock (Generic)
  deriving anyclass (Data.Aeson.ToJSON)

ls_hs_symmetry_more_info :: CString -> IO CString
ls_hs_symmetry_more_info jsonString = propagateErrorToC nullPtr $ do
  !(symm :: Symmetry) <- decodeCString jsonString
  newCencoded $
    SymmetryInfo
      { periodicity = getPeriodicity symm.permutation
      , phase = realToFrac symm.phase
      }

ls_hs_symmetries_from_generators :: CString -> IO CString
ls_hs_symmetries_from_generators jsonString = propagateErrorToC nullPtr $ do
  newCencoded . toList . either error id . groupRepresentationFromGenerators =<< decodeCString jsonString

-- {{{ Basis

ls_hs_clone_basis :: Ptr Cbasis -> IO (MutablePtr Cbasis)
ls_hs_clone_basis = cloneCbasis

ls_hs_destroy_basis :: MutablePtr Cbasis -> IO ()
ls_hs_destroy_basis = destroyCbasis

ls_hs_basis_to_json :: Ptr Cbasis -> IO CString
ls_hs_basis_to_json cBasis = withCbasis cBasis newCencoded

ls_hs_basis_from_json :: CString -> IO (MutablePtr Cbasis)
ls_hs_basis_from_json cStr =
  propagateErrorToC nullPtr $
    foldSomeBasis newCbasis =<< decodeCString cStr

ls_hs_min_state_estimate :: Ptr Cbasis -> IO Word64
ls_hs_min_state_estimate p = propagateErrorToC 0 $
  withCbasis p $ \someBasis ->
    withSomeBasis someBasis $ \basis ->
      let (BasisState n (BitString x)) = minStateEstimate basis
       in if n > 64
            then error "minimal state is not representable as a 64-bit integer"
            else pure $ fromIntegral x

ls_hs_max_state_estimate :: Ptr Cbasis -> IO Word64
ls_hs_max_state_estimate p = propagateErrorToC 0 $
  withCbasis p $ \someBasis ->
    withSomeBasis someBasis $ \basis ->
      let (BasisState n (BitString x)) = maxStateEstimate basis
       in if n > 64
            then error "maximal state is not representable as a 64-bit integer"
            else pure $ fromIntegral x

ls_hs_basis_has_fixed_hamming_weight :: Ptr Cbasis -> IO CBool
ls_hs_basis_has_fixed_hamming_weight basis =
  fromBool <$> withCbasis basis (foldSomeBasis (pure . hasFixedHammingWeight))

ls_hs_basis_has_spin_inversion_symmetry :: Ptr Cbasis -> IO CBool
ls_hs_basis_has_spin_inversion_symmetry basis =
  fromBool <$> withCbasis basis (foldSomeBasis (pure . hasSpinInversionSymmetry))

ls_hs_basis_has_permutation_symmetries :: Ptr Cbasis -> IO CBool
ls_hs_basis_has_permutation_symmetries basis =
  fromBool <$> withCbasis basis (foldSomeBasis (pure . hasPermutationSymmetries))

ls_hs_basis_requires_projection :: Ptr Cbasis -> IO CBool
ls_hs_basis_requires_projection basis =
  fromBool <$> withCbasis basis (foldSomeBasis (pure . requiresProjection))

foreign import ccall safe "ls_hs_build_representatives"
  ls_hs_build_representatives :: Ptr Cbasis -> Word64 -> Word64 -> IO ()

-- foreign export ccall "ls_hs_basis_build"
--   ls_hs_basis_build :: Ptr Cbasis -> IO ()

ls_hs_basis_build :: Ptr Cbasis -> IO ()
ls_hs_basis_build p = propagateErrorToC () $ do
  withCbasis p $ \someBasis ->
    withSomeBasis someBasis $ \basis ->
      if getNumberBits basis <= 64
        then do
          let (BasisState _ (BitString lower)) = minStateEstimate basis
              (BasisState _ (BitString upper)) = maxStateEstimate basis
          ls_hs_build_representatives p (fromIntegral lower) (fromIntegral upper)
        else error "too many bits"

-- foreign export ccall "ls_hs_basis_is_built"
--   ls_hs_basis_is_built :: Ptr Cbasis -> IO CBool

ls_hs_basis_is_built :: Ptr Cbasis -> IO CBool
ls_hs_basis_is_built =
  pure . fromBool . (/= nullPtr) . external_array_elts . cbasis_representatives <=< peek

-- foreign export ccall "ls_hs_basis_number_words"
--   ls_hs_basis_number_words :: Ptr Cbasis -> IO CInt

ls_hs_basis_number_bits :: Ptr Cbasis -> IO CInt
ls_hs_basis_number_bits basisPtr =
  fromIntegral <$> withCbasis basisPtr (foldSomeBasis (pure . getNumberBits))

ls_hs_basis_number_words :: Ptr Cbasis -> IO CInt
ls_hs_basis_number_words basisPtr =
  fromIntegral <$> withCbasis basisPtr (foldSomeBasis (pure . getNumberWords))

ls_hs_basis_is_real :: Ptr Cbasis -> IO CBool
ls_hs_basis_is_real basisPtr =
  withCbasis basisPtr (foldSomeBasis (pure . fromBool . isBasisReal))

-- foreign export ccall "ls_hs_basis_state_to_string"
--   ls_hs_basis_state_to_string :: Ptr Cbasis -> Ptr Word64 -> IO CString

ls_hs_basis_state_to_string :: Ptr Cbasis -> Ptr Word64 -> IO CString
ls_hs_basis_state_to_string basisPtr statePtr =
  withCbasis basisPtr $ \someBasis ->
    withSomeBasis someBasis $ \(basis :: Basis t) -> do
      let numberBits = getNumberBits basis
          numberWords = getNumberWords basis
      state <- BasisState @t numberBits <$> readBitString numberWords statePtr
      newCString . encodeUtf8 . toPrettyText $ state

-- foreign export ccall "ls_hs_fixed_hamming_state_to_index"
--   ls_hs_fixed_hamming_state_to_index :: Word64 -> CPtrdiff

ls_hs_fixed_hamming_state_to_index :: Word64 -> CPtrdiff
ls_hs_fixed_hamming_state_to_index = fromIntegral . fixedHammingStateToIndex

-- foreign export ccall "ls_hs_fixed_hamming_index_to_state"
--   ls_hs_fixed_hamming_index_to_state :: CPtrdiff -> CInt -> Word64

ls_hs_fixed_hamming_index_to_state :: CPtrdiff -> CInt -> Word64
ls_hs_fixed_hamming_index_to_state index hammingWeight =
  fixedHammingIndexToState (fromIntegral hammingWeight) (fromIntegral index)

-- }}}

-- {{{ Expr

-- foreign export ccall "ls_hs_expr_to_json"
--   ls_hs_expr_to_json :: Ptr Cexpr -> IO CString

ls_hs_expr_to_json :: Ptr Cexpr -> IO CString
ls_hs_expr_to_json cExpr =
  withCexpr cExpr $ \expr -> do
    newCString $ toStrict (Data.Aeson.encode expr)

-- foreign export ccall "ls_hs_expr_from_json"
--   ls_hs_expr_from_json :: CString -> IO (Ptr Cexpr)

ls_hs_expr_from_json :: CString -> IO (MutablePtr Cexpr)
ls_hs_expr_from_json cStr = propagateErrorToC nullPtr $ do
  !expr <- decodeCString cStr
  newCexpr expr

-- foreign export ccall "ls_hs_destroy_expr"
--   ls_hs_destroy_expr :: Ptr Cexpr -> IO ()

ls_hs_destroy_expr :: MutablePtr Cexpr -> IO ()
ls_hs_destroy_expr = destroyCexpr

-- foreign export ccall "ls_hs_expr_to_string"
--   ls_hs_expr_to_string :: Ptr Cexpr -> IO CString

ls_hs_expr_to_string :: Ptr Cexpr -> IO CString
ls_hs_expr_to_string p =
  propagateErrorToC nullPtr $
    withCexpr p $
      newCString . encodeUtf8 . toPrettyText

-- foreign export ccall "ls_hs_expr_plus"
--   ls_hs_expr_plus :: Ptr Cexpr -> Ptr Cexpr -> IO (Ptr Cexpr)

ls_hs_expr_plus :: Ptr Cexpr -> Ptr Cexpr -> IO (MutablePtr Cexpr)
ls_hs_expr_plus a b = propagateErrorToC nullPtr $ withCexpr2 a b (+) >>= newCexpr

-- foreign export ccall "ls_hs_expr_minus"
--   ls_hs_expr_minus :: Ptr Cexpr -> Ptr Cexpr -> IO (Ptr Cexpr)

ls_hs_expr_minus :: Ptr Cexpr -> Ptr Cexpr -> IO (MutablePtr Cexpr)
ls_hs_expr_minus a b = propagateErrorToC nullPtr $ withCexpr2 a b (-) >>= newCexpr

-- foreign export ccall "ls_hs_expr_times"
--   ls_hs_expr_times :: Ptr Cexpr -> Ptr Cexpr -> IO (Ptr Cexpr)

ls_hs_expr_times :: Ptr Cexpr -> Ptr Cexpr -> IO (MutablePtr Cexpr)
ls_hs_expr_times a b = propagateErrorToC nullPtr $ withCexpr2 a b (*) >>= newCexpr

-- foreign export ccall "ls_hs_expr_scale"
--   ls_hs_expr_scale :: Ptr Cscalar -> Ptr Cexpr -> IO (MutablePtr Cexpr)

ls_hs_expr_scale :: Ptr Cscalar -> Ptr Cexpr -> IO (MutablePtr Cexpr)
ls_hs_expr_scale c_z c_a =
  propagateErrorToC nullPtr $
    withCexpr c_a $ \a -> do
      z <- fromComplexDouble <$> peek c_z
      newCexpr $ scale (z :: ComplexRational) a

-- foreign export ccall "ls_hs_replace_indices"
--   ls_hs_replace_indices :: Ptr Cexpr -> FunPtr Creplace_index -> IO (Ptr Cexpr)

ls_hs_replace_indices :: Ptr Cexpr -> FunPtr Creplace_index -> IO (MutablePtr Cexpr)
ls_hs_replace_indices exprPtr fPtr =
  propagateErrorToC nullPtr $
    withCexpr exprPtr $ \expr -> do
      let f :: Int -> Int -> IO (Int, Int)
          f !s !i =
            alloca $ \spinPtr ->
              alloca $ \sitePtr -> do
                mkCreplace_index fPtr (fromIntegral s) (fromIntegral i) spinPtr sitePtr
                (,)
                  <$> (fromIntegral <$> peek spinPtr)
                  <*> (fromIntegral <$> peek sitePtr)
      newCexpr =<< case expr of
        SomeExpr SpinTag terms ->
          SomeExpr SpinTag . simplifyExpr <$> mapIndicesM (fmap snd . f 0) terms
        SomeExpr SpinlessFermionTag terms ->
          SomeExpr SpinlessFermionTag . simplifyExpr <$> mapIndicesM (fmap snd . f 0) terms
        SomeExpr SpinfulFermionTag terms ->
          let f' (s, i) = do
                (s', i') <- f (fromEnum s) i
                pure (toEnum s', i')
           in SomeExpr SpinfulFermionTag . simplifyExpr <$> mapIndicesM f' terms

-- foreign export ccall "ls_hs_expr_equal"
--   ls_hs_expr_equal :: Ptr Cexpr -> Ptr Cexpr -> IO CBool

ls_hs_expr_equal :: Ptr Cexpr -> Ptr Cexpr -> IO CBool
ls_hs_expr_equal aPtr bPtr =
  withCexpr aPtr $ \a ->
    withCexpr bPtr $ \b ->
      pure $ fromBool (a == b)

-- foreign export ccall "ls_hs_expr_adjoint"
--   ls_hs_expr_adjoint :: Ptr Cexpr -> IO (Ptr Cexpr)

ls_hs_expr_adjoint :: Ptr Cexpr -> IO (MutablePtr Cexpr)
ls_hs_expr_adjoint = flip withCexpr $ \(SomeExpr tag expr) ->
  newCexpr $ SomeExpr tag (conjugateExpr expr)

-- foreign export ccall "ls_hs_expr_is_hermitian"
--   ls_hs_expr_is_hermitian :: Ptr Cexpr -> IO CBool

ls_hs_expr_is_hermitian :: Ptr Cexpr -> IO CBool
ls_hs_expr_is_hermitian = flip withCexpr $ foldSomeExpr (pure . fromBool . isHermitianExpr)

-- foreign export ccall "ls_hs_expr_is_real"
--   ls_hs_expr_is_real :: Ptr Cexpr -> IO CBool

ls_hs_expr_is_real :: Ptr Cexpr -> IO CBool
ls_hs_expr_is_real = flip withCexpr $ foldSomeExpr (pure . fromBool . isRealExpr)

-- foreign export ccall "ls_hs_expr_is_identity"
--   ls_hs_expr_is_identity :: Ptr Cexpr -> IO CBool

ls_hs_expr_is_identity :: Ptr Cexpr -> IO CBool
ls_hs_expr_is_identity = flip withCexpr $ foldSomeExpr (pure . fromBool . isIdentityExpr)

-- }}}

-- {{{ Operator

-- foreign export ccall "ls_hs_create_operator"
--   ls_hs_create_operator :: Ptr Cbasis -> Ptr Cexpr -> IO (Ptr Coperator)

ls_hs_create_operator :: Ptr Cbasis -> Ptr Cexpr -> IO (MutablePtr Coperator)
ls_hs_create_operator basisPtr exprPtr = propagateErrorToC nullPtr $ do
  withCbasis basisPtr $ \someBasis ->
    withSomeBasis someBasis $ \basis ->
      withCexpr exprPtr $ \someExpr ->
        withSomeExpr someExpr $ \expr ->
          case matchParticleType2 basis expr of
            Just HRefl -> newCoperator (Just basisPtr) (mkOperator basis expr)
            Nothing -> error "basis and expression have different particle types"

-- foreign export ccall "ls_hs_clone_operator"
--   ls_hs_clone_operator :: Ptr Coperator -> IO (Ptr Coperator)

ls_hs_clone_operator :: Ptr Coperator -> IO (MutablePtr Coperator)
ls_hs_clone_operator = cloneCoperator

-- foreign export ccall "ls_hs_destroy_operator"
--   destroyCoperator :: Ptr Coperator -> IO ()

ls_hs_destroy_operator :: MutablePtr Coperator -> IO ()
ls_hs_destroy_operator = destroyCoperator

-- foreign export ccall "ls_hs_operator_max_number_off_diag"
--   ls_hs_operator_max_number_off_diag :: Ptr Coperator -> IO CInt

ls_hs_operator_max_number_off_diag :: MutablePtr Coperator -> IO CInt
ls_hs_operator_max_number_off_diag opPtr =
  fmap fromIntegral $
    withCoperator opPtr $ \someOp ->
      withSomeOperator someOp $ \op ->
        pure $ maxNumberOffDiag op

-- foreign export ccall "ls_hs_operator_get_expr"
--   ls_hs_operator_get_expr :: Ptr Coperator -> IO (Ptr Cexpr)

ls_hs_operator_get_expr :: Ptr Coperator -> IO (MutablePtr Cexpr)
ls_hs_operator_get_expr opPtr =
  withCoperator opPtr $ \someOp ->
    withSomeOperator someOp $ \op ->
      newCexpr $
        SomeExpr
          (getParticleTag . opBasis $ op)
          (opTerms op)

-- foreign export ccall "ls_hs_operator_get_basis"
--   ls_hs_operator_get_basis :: Ptr Coperator -> IO (Ptr Cbasis)

ls_hs_operator_get_basis :: Ptr Coperator -> IO (MutablePtr Cbasis)
ls_hs_operator_get_basis = ls_hs_clone_basis . coperator_basis <=< peek

ls_hs_operator_abelian_representations :: Ptr Coperator -> IO CString
ls_hs_operator_abelian_representations opPtr = do
  withCoperator opPtr $ \someOp ->
    withSomeOperator someOp $ \op ->
      newCString . toStrict . Data.Aeson.encode $
        operatorAbelianRepresentations op

-- foreign export ccall "ls_hs_load_hamiltonian_from_yaml"
--   ls_hs_load_hamiltonian_from_yaml :: CString -> IO (Ptr Coperator)
--
-- ls_hs_load_hamiltonian_from_yaml cFilename =
--   foldSomeOperator borrowCoperator =<< hamiltonianFromYAML =<< peekUtf8 cFilename

-- foreign import ccall "ls_hs_destroy_basis_v2"
--   ls_hs_destroy_basis_v2 :: Ptr Cbasis -> IO ()

-- foreign import ccall "ls_hs_destroy_operator_v2"
--   ls_hs_destroy_operator_v2 :: Ptr Coperator -> IO ()

ls_hs_prepare_hphi :: Ptr Coperator -> CString -> IO ()
ls_hs_prepare_hphi opPtr pathPtr = propagateErrorToC () $ do
  path <- peekUtf8 pathPtr
  withCoperator opPtr $ \op -> convertedToInteractions op (prepareHPhi path)

ls_hs_prepare_mvmc :: Ptr Coperator -> CString -> IO ()
ls_hs_prepare_mvmc opPtr pathPtr = propagateErrorToC () $ do
  path <- peekUtf8 pathPtr
  withCoperator opPtr $ \op -> convertedToInteractions op (prepareVMC path)

toCyaml_config :: ConfigSpec -> IO (Ptr Cyaml_config)
toCyaml_config (ConfigSpec basis maybeHamiltonian observables) = do
  p <- malloc
  -- print basis
  basisPtr <- withSomeBasis basis newCbasis
  -- withSomeBasis basis $ \b ->
  --   withForeignPtr (basisContents b) $ \ptr -> do
  --     _ <- basisIncRefCount ptr
  --     pure ptr
  -- print "1)"
  -- borrowCbasis
  hamiltonianPtr <- case maybeHamiltonian of
    Just h -> withSomeOperator h (newCoperator (Just basisPtr))
    Nothing -> pure nullPtr
  -- print "2)"
  observablesPtr <-
    (newArray =<<) $
      G.toList <$> G.mapM (foldSomeOperator (newCoperator (Just basisPtr))) observables
  -- print "3)"
  poke p $
    Cyaml_config basisPtr hamiltonianPtr (fromIntegral (G.length observables)) observablesPtr
  pure p

-- foreign export ccall "ls_hs_load_yaml_config"
--   ls_hs_load_yaml_config :: CString -> IO (Ptr Cyaml_config)

ls_hs_load_yaml_config :: CString -> IO (MutablePtr Cyaml_config)
ls_hs_load_yaml_config cFilename =
  propagateErrorToC nullPtr $
    toCyaml_config =<< configFromYAML =<< peekUtf8 cFilename

-- foreign export ccall "ls_hs_destroy_yaml_config"
--   ls_hs_destroy_yaml_config :: Ptr Cyaml_config -> IO ()

ls_hs_destroy_yaml_config :: MutablePtr Cyaml_config -> IO ()
ls_hs_destroy_yaml_config p
  | p == nullPtr = pure ()
  | otherwise = do
      (Cyaml_config basisPtr hamiltonianPtr numberObservables observablesPtr) <- peek p
      -- logDebug' "ls_hs_destroy_yaml_config 1) ..."
      forM_ @[] [0 .. fromIntegral numberObservables - 1] $
        destroyCoperator <=< peekElemOff observablesPtr
      -- logDebug' "ls_hs_destroy_yaml_config 2) ..."
      when (observablesPtr /= nullPtr) $ free observablesPtr
      -- logDebug' "ls_hs_destroy_yaml_config 3) ..."
      when (hamiltonianPtr /= nullPtr) $ destroyCoperator hamiltonianPtr
      -- logDebug' "ls_hs_destroy_yaml_config 4) ..."
      destroyCbasis basisPtr
      -- logDebug' "ls_hs_destroy_yaml_config 5) ..."
      free p

-- ls_hs_operator_pretty_terms :: Ptr Coperator -> IO CString
-- ls_hs_operator_pretty_terms p =
--   withReconstructedOperator p $ \op ->
--     newCString
--       . encodeUtf8
--       . renderStrict
--       . Pretty.layoutPretty (Pretty.LayoutOptions Pretty.Unbounded)
--       . pretty
--       $ opTerms (opHeader op)

-- foreign export ccall "ls_hs_operator_pretty_terms"
--   ls_hs_operator_pretty_terms :: Ptr Coperator -> IO CString

-- foreign export ccall "ls_hs_fatal_error"
--   ls_hs_fatal_error :: CString -> CString -> IO ()

$(pure [])

typesTable
  [ ([t|()|], "void")
  , ([t|CBool|], "bool")
  , ([t|CInt|], "int")
  , ([t|Word64|], "uint64_t")
  , ([t|CPtrdiff|], "ptrdiff_t")
  , ([t|Double|], "double")
  , ([t|CString|], "char const *")
  , ([t|Cbasis|], "ls_hs_basis")
  , ([t|Cexpr|], "ls_hs_expr")
  , ([t|Creplace_index|], "ls_hs_index_replacement_type")
  , ([t|Cscalar|], "ls_hs_scalar")
  , ([t|Coperator|], "ls_hs_operator")
  , ([t|Cyaml_config|], "ls_hs_yaml_config")
  ]

headerFile "lattice_symmetries_functions.h"

addVerbatimPrefix
  [ "#include \"lattice_symmetries_types.h\""
  , "#include <stdint.h>"
  , ""
  , "#if defined(__cplusplus)"
  , "extern \"C\" {"
  , "#endif"
  , ""
  , "/* python-cffi: START */"
  ]

addVerbatimSuffix
  [ "/* python-cffi: STOP */"
  , ""
  , "#if defined(__cplusplus)"
  , "} // extern \"C\""
  , "#endif"
  ]

addDeclarations
  [ "ls_hs_symmetry_more_info"
  , -- , "ls_hs_destroy_symmetry"
    -- , "ls_hs_symmetry_sector"
    -- , "ls_hs_symmetry_phase"
    -- , "ls_hs_symmetry_periodicity"
    -- , "ls_hs_symmetry_length"
    -- , "ls_hs_symmetry_permutation"
    -- "ls_hs_destroy_permutation"
    "ls_hs_symmetries_from_generators"
  , "ls_hs_clone_basis"
  , "ls_hs_destroy_basis"
  , "ls_hs_basis_from_json"
  , "ls_hs_basis_to_json"
  , "ls_hs_destroy_string"
  , "ls_hs_min_state_estimate"
  , "ls_hs_max_state_estimate"
  , "ls_hs_basis_has_fixed_hamming_weight"
  , "ls_hs_basis_has_spin_inversion_symmetry"
  , "ls_hs_basis_has_permutation_symmetries"
  , "ls_hs_basis_requires_projection"
  , "ls_hs_basis_build"
  , "ls_hs_basis_is_built"
  , "ls_hs_basis_number_bits"
  , "ls_hs_basis_number_words"
  , "ls_hs_basis_state_to_string"
  , "ls_hs_basis_is_real"
  , "ls_hs_fixed_hamming_state_to_index"
  , "ls_hs_fixed_hamming_index_to_state"
  , "ls_hs_expr_to_json"
  , "ls_hs_expr_from_json"
  , "ls_hs_destroy_expr"
  , "ls_hs_expr_to_string"
  , "ls_hs_expr_plus"
  , "ls_hs_expr_minus"
  , "ls_hs_expr_times"
  , "ls_hs_expr_scale"
  , "ls_hs_replace_indices"
  , "ls_hs_expr_equal"
  , "ls_hs_expr_adjoint"
  , "ls_hs_expr_is_hermitian"
  , "ls_hs_expr_is_real"
  , "ls_hs_expr_is_identity"
  , "ls_hs_create_operator"
  , "ls_hs_clone_operator"
  , "ls_hs_destroy_operator"
  , "ls_hs_operator_max_number_off_diag"
  , "ls_hs_operator_get_expr"
  , "ls_hs_operator_get_basis"
  , "ls_hs_operator_abelian_representations"
  , "ls_hs_prepare_hphi"
  , "ls_hs_prepare_mvmc"
  , "ls_hs_load_yaml_config"
  , "ls_hs_destroy_yaml_config"
  ]

-- $( do
-- info <- reify ''Csymmetry
-- runIO $ print info
-- Just nm <- lookupValueName "ls_hs_symmetry_from_json"
-- (VarI (Name name flavour) tp _) <- reify nm
-- runIO $ print (name, flavour, tp)
-- s <- genCType =<< [t|CInt|]
-- runIO $ print s
-- s <- genCType =<< [t|Ptr CInt|]
-- runIO $ print s
-- s <- genCType =<< [t|MutablePtr CInt|]
-- runIO $ print s
-- s <- getArgTypes <$> [t|Int -> Float -> IO Int|]
-- runIO $ print s
-- s <- genCSignature "ls_hs_symmetry_from_json"
-- runIO $ print s
-- printDeclarations
-- [d||]
-- )
