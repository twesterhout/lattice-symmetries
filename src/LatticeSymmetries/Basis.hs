-- {-# LANGUAGE CApiFFI #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE OverloadedLists #-}
{-# LANGUAGE QuantifiedConstraints #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeOperators #-}

module LatticeSymmetries.Basis
  ( -- * High-level interface
    ParticleTy (..),
    BasisState (..),
    Basis (..),
    -- BasisHeader (..),
    SpinfulOccupation (..),
    SomeBasis (..),
    -- basisBuild,
    IndexType,
    GeneratorType,
    IsBasis,
    withSomeBasis,
    foldSomeBasis,

    -- ** Creating bases

    -- mkSpinBasis,
    -- mkSpinlessFermionicBasis,
    -- mkSpinfulFermionicBasis,

    -- ** Querying information
    getNumberBits,
    getNumberWords,
    minStateEstimate,
    maxStateEstimate,
    isStateIndexIdentity,
    hasFixedHammingWeight,
    hasSpinInversionSymmetry,
    hasPermutationSymmetries,
    requiresProjection,
    getParticleTag,
    isBasisReal,
    fixedHammingStateToIndex,
    fixedHammingIndexToState,

    -- ** Low-level interface
    Cbasis (..),
    -- basisFromYAML,
    -- objectFromYAML,
    -- Factor,
    -- borrowCbasis,
    newCbasis,
    cloneCbasis,
    destroyCbasis,
    withCbasis,
    -- Cparticle_type (..),
    -- Cbasis_kernels (..),
    -- createCbasis_kernels,
    -- destroyCbasis_kernels,
    -- stateIndex,
    flattenIndex,
    matchParticleType2,
    -- withParticleType,
    -- withReconstructedBasis,
  )
where

import Control.Exception.Safe (bracket, handleAny)
import Data.Aeson
import Data.Aeson.Types (Pair)
import Data.Bits
import Data.ByteString (packCString, useAsCString)
import Data.ByteString.Internal (ByteString (..))
import qualified Data.Maybe (fromJust)
import Data.Yaml.Aeson
import Foreign.C.String (CString)
import Foreign.C.Types
import Foreign.ForeignPtr
import Foreign.Marshal.Alloc (alloca, free, mallocBytes)
import Foreign.Marshal.Utils (copyBytes, fromBool, new, with)
import Foreign.Ptr
import Foreign.StablePtr
import Foreign.Storable
import GHC.ForeignPtr
import LatticeSymmetries.Algebra (Algebra (..))
import LatticeSymmetries.BitString
import LatticeSymmetries.FFI
import LatticeSymmetries.Generator
import LatticeSymmetries.Group
-- import LatticeSymmetries.IO
import LatticeSymmetries.NonbranchingTerm
import LatticeSymmetries.Utils
import Prettyprinter (Doc, Pretty (..))
import qualified Prettyprinter as Pretty
import Prettyprinter.Render.Text (renderStrict)
import System.IO.Unsafe (unsafePerformIO)
import Type.Reflection

-- | Hilbert space basis vector parametrized by the particle type.
data BasisState (t :: ParticleTy) = BasisState {-# UNPACK #-} !Int !BitString
  deriving stock (Show, Eq, Ord)

unsafeCastBasisState :: BasisState t1 -> BasisState t2
unsafeCastBasisState (BasisState n bits) = BasisState n bits

prettyBitString :: Int -> Integer -> Doc ann
prettyBitString n bits = mconcat $ (prettyBool . testBit bits) <$> reverse [0 .. n - 1]
  where
    prettyBool True = "1"
    prettyBool False = "0"

matchParticleType2 ::
  forall (t1 :: ParticleTy) (t2 :: ParticleTy) proxy1 proxy2.
  (Typeable t1, Typeable t2) =>
  proxy1 t1 ->
  proxy2 t2 ->
  Maybe (t1 :~~: t2)
matchParticleType2 _ _ = case eqTypeRep (typeRep @t1) (typeRep @t2) of
  Just HRefl -> Just HRefl
  Nothing -> Nothing

instance Typeable t => Pretty (BasisState t) where
  pretty x =
    case particleDispatch @t of
      SpinTag -> prettySpin x
      SpinfulFermionTag -> prettyFermion x
      SpinlessFermionTag ->
        let (BasisState n bits) = x
         in prettySpin (BasisState n bits :: BasisState 'SpinTy)
    where
      prettySpin (BasisState n bits) = "|" <> prettyBitString n (unBitString bits) <> "⟩"
      prettyFermion (BasisState n bits) =
        let up = unBitString bits `shiftR` (n `div` 2)
         in mconcat
              [ "|",
                prettyBitString (n `div` 2) up,
                "⟩",
                "|",
                prettyBitString (n `div` 2) (unBitString bits),
                "⟩"
              ]

-- instance Pretty (BasisState 'SpinTy) where
--   pretty (BasisState n bits) = "|" <> prettyBitString n (unBitString bits) <> "⟩"
--
-- instance Pretty (BasisState 'SpinfulFermionTy)
--
-- instance Pretty (BasisState 'SpinlessFermionTy) where
--   pretty (BasisState n bits) = pretty (BasisState n bits :: BasisState 'SpinTy)

-- instance Typeable t => Pretty (BasisState t) where
--   pretty (BasisState n bits)
--     | Just HRefl <- eqTypeRep (typeRep @t) (typeRep @'SpinTy) =
--         "|" <> prettyBitString n (unBitString bits) <> "⟩"
--     | Just HRefl <- eqTypeRep (typeRep @t) (typeRep @'SpinfulFermionTy) =
--         let up = unBitString bits `shiftR` (n `div` 2)
--          in mconcat
--               [ "|",
--                 prettyBitString (n `div` 2) up,
--                 "⟩",
--                 "|",
--                 prettyBitString (n `div` 2) (unBitString bits),
--                 "⟩"
--               ]
--     | Just HRefl <- eqTypeRep (typeRep @t) (typeRep @'SpinlessFermionTy) =
--         pretty (BasisState n bits :: BasisState 'SpinTy)
--     | otherwise = error "this cannot happen by construction"

class
  ( -- Enum (GeneratorType t),
    -- Bounded (GeneratorType t),
    -- HasMatrixRepresentation (GeneratorType t),
    -- HasNonbranchingRepresentation (Generator Int (GeneratorType t)),
    Typeable t,
    Algebra (GeneratorType t),
    Pretty (Generator (IndexType t) (GeneratorType t)),
    HasProperGeneratorType t,
    HasProperIndexType t
    -- Ord (IndexType t),
    -- HasSiteIndex (IndexType t),
    -- Pretty (IndexType t),
    -- Pretty (GeneratorType t),
    -- Pretty (BasisState t),
    -- HasSiteIndex (IndexType t),
  ) =>
  IsBasis t

instance IsBasis 'SpinTy

instance IsBasis 'SpinfulFermionTy

instance IsBasis 'SpinlessFermionTy

-- isBasisDict :: proxy (t :: ParticleTy) -> Dict (IsBasis t)
-- isBasisDict _
--   | Just HRefl <- eqTypeRep (typeRep @t) (typeRep @'SpinTy) =
--     "|" <> prettyBitString n (unBitString bits) <> "⟩"
--   | Just HRefl <- eqTypeRep (typeRep @t) (typeRep @'SpinfulFermionTy) =
--     let up = unBitString bits `shiftR` (n `div` 2)
--      in mconcat
--           [ "|",
--             prettyBitString (n `div` 2) up,
--             "⟩",
--             "|",
--             prettyBitString (n `div` 2) (unBitString bits),
--             "⟩"
--           ]
--   | Just HRefl <- eqTypeRep (typeRep @t) (typeRep @'SpinlessFermionTy) =
--     pretty (BasisState n bits :: BasisState 'SpinTy)

data Basis (t :: ParticleTy) where
  SpinBasis :: !Int -> !(Maybe Int) -> !(Maybe Int) -> !Symmetries -> Basis 'SpinTy
  SpinfulFermionBasis :: !Int -> !SpinfulOccupation -> Basis 'SpinfulFermionTy
  SpinlessFermionBasis :: !Int -> !(Maybe Int) -> Basis 'SpinlessFermionTy

deriving stock instance Show (Basis t)

deriving stock instance Eq (Basis t)

-- data Basis t = Basis
--   { basisHeader :: !(BasisHeader t),
--     basisContents :: !(ForeignPtr Cbasis)
--   }
--   deriving stock (Show)

data SomeBasis where
  SomeBasis :: IsBasis t => Basis t -> SomeBasis

-- data SomeBasis where
--   SomeBasis :: IsBasis t => Basis t -> SomeBasis

deriving stock instance Show SomeBasis

instance Eq SomeBasis where
  (==) (SomeBasis a@(SpinBasis _ _ _ _)) (SomeBasis b@(SpinBasis _ _ _ _)) = a == b
  (==) (SomeBasis a@(SpinBasis _ _ _ _)) (SomeBasis _) = False
  (==) (SomeBasis a@(SpinfulFermionBasis _ _)) (SomeBasis b@(SpinfulFermionBasis _ _)) = a == b
  (==) (SomeBasis a@(SpinfulFermionBasis _ _)) (SomeBasis _) = False
  (==) (SomeBasis a@(SpinlessFermionBasis _ _)) (SomeBasis b@(SpinlessFermionBasis _ _)) = a == b
  (==) (SomeBasis a@(SpinlessFermionBasis _ _)) (SomeBasis _) = False

withSomeBasis :: SomeBasis -> (forall t. IsBasis t => Basis t -> a) -> a
withSomeBasis x f = case x of SomeBasis basis -> f basis
{-# INLINE withSomeBasis #-}

foldSomeBasis :: (forall t. IsBasis t => Basis t -> a) -> SomeBasis -> a
foldSomeBasis f x = withSomeBasis x f
{-# INLINE foldSomeBasis #-}

-- withSomeBasisHeader :: SomeBasisHeader -> (forall t. IsBasis t => BasisHeader t -> a) -> a
-- withSomeBasisHeader x f = case x of
--   SomeBasisHeader header@(SpinHeader _ _ _ _) -> f header
--   SomeBasisHeader header@(SpinfulFermionHeader _ _) -> f header
--   SomeBasisHeader header@(SpinlessFermionHeader _ _) -> f header
-- {-# INLINE withSomeBasisHeader #-}

-- particleTyDispatch :: Text -> (forall (t :: ParticleTy). Proxy t -> a) -> a
-- particleTyDispatch t f
--   | t == "spin" || t == "spin-1/2" = f (Proxy @SpinTy)
--   | t == "spinless" || t == "spinless-fermion" || t == "spinless fermion" = f (Proxy @SpinlessFermionTy)
--   | t == "spinful" || t == "spinful-fermion" || t == "spinful fermion" = f (Proxy @SpinfulFermionTy)
--   | otherwise = error "invalid particle type"

instance FromJSON SomeBasis where
  parseJSON = basisHeaderFromJSON

instance ToJSON SomeBasis where
  toJSON = foldSomeBasis toJSON

instance ToJSON (Basis t) where
  toJSON = basisHeaderToJSON

parseSpinfulOccupation :: Object -> Parser SpinfulOccupation
parseSpinfulOccupation v = do
  (r :: Maybe Value) <- v .:? "number_particles"
  case r of
    Nothing -> pure SpinfulNoOccupation
    Just n -> case n of
      Number _ -> fmap SpinfulTotalParticles $ parseJSON n
      Array _ -> fmap (\(up, down) -> SpinfulPerSector up down) $ parseJSON n
      _ -> mzero

encodeSpinfulOccupation :: SpinfulOccupation -> [Pair]
encodeSpinfulOccupation SpinfulNoOccupation = []
encodeSpinfulOccupation (SpinfulTotalParticles n) = [("number_particles", toJSON n)]
encodeSpinfulOccupation (SpinfulPerSector up down) =
  [("number_particles", toJSON ([up, down] :: [Int]))]

basisHeaderToJSON :: Basis t -> Value
basisHeaderToJSON x
  | (SpinBasis n h i g) <- x =
    object
      [ "particle" .= SpinTy,
        "number_spins" .= n,
        "hamming_weight" .= h,
        "spin_inversion" .= i,
        "symmetries" .= g
      ]
  | (SpinfulFermionBasis n o) <- x =
    object $
      ["particle" .= SpinfulFermionTy, "number_sites" .= n] <> encodeSpinfulOccupation o
  | (SpinlessFermionBasis n o) <- x =
    object $
      ["particle" .= SpinlessFermionTy, "number_sites" .= n, "number_particles" .= o]

basisHeaderFromJSON :: Value -> Parser SomeBasis
basisHeaderFromJSON = withObject "Basis" $ \v -> do
  tp <- v .:! "particle" .!= SpinTy
  case tp of
    SpinTy ->
      fmap SomeBasis . join $
        mkSpinHeader
          <$> (v .: "number_spins")
          <*> (v .:? "hamming_weight")
          <*> (v .:? "spin_inversion")
          <*> (v .:! "symmetries" .!= emptySymmetries)
    SpinlessFermionTy ->
      fmap SomeBasis $
        SpinlessFermionBasis <$> (v .: "number_sites") <*> (v .:? "number_particles")
    SpinfulFermionTy ->
      fmap SomeBasis $
        SpinfulFermionBasis <$> (v .: "number_sites") <*> parseSpinfulOccupation v

mkSpinHeader ::
  MonadFail m =>
  -- | Number of sites
  Int ->
  -- | Hamming weight
  Maybe Int ->
  -- | Spin inversion
  Maybe Int ->
  -- | Lattice symmetries
  Symmetries ->
  m (Basis 'SpinTy)
mkSpinHeader n _ _ _
  | n <= 0 = fail $ "invalid number of spins: " <> show n
mkSpinHeader n (Just h) _ _
  | h < 0 || h > n = fail $ "invalid Hamming weight: " <> show h
mkSpinHeader _ _ (Just i) _
  | i /= (-1) && i /= 1 = fail $ "invalid spin inversion: " <> show i
mkSpinHeader n (Just h) (Just _) _
  | n /= 2 * h = fail $ "invalid spin inversion: " <> show n <> " spins, but Hamming weight is " <> show h
mkSpinHeader n _ _ g
  | not (nullSymmetries g) && symmetriesGetNumberBits g /= n = fail "invalid symmetries"
mkSpinHeader n h i g = pure $ SpinBasis n h i g

-- instance IsBasis t => FromJSON (Basis t) where
--   parseJSON x = fmap basisFromHeader (parseJSON x)

-- instance FromJSON SomeBasis where
--   parseJSON x = do
--     header <- parseJSON x
--     pure $ withSomeBasisHeader header (SomeBasis . basisFromHeader)

-- instance ToJSON (Basis t) where
--   toJSON = toJSON . basisHeader

-- instance Eq (Basis t) where
--   (==) a b = basisHeader a == basisHeader b

-- type Factor t = Generator (IndexType t) (GeneratorType t)

data SpinfulOccupation
  = SpinfulNoOccupation
  | SpinfulTotalParticles !Int
  | -- | @SpinfulPerSector numberDown numberUp@
    SpinfulPerSector !Int !Int
  deriving stock (Show, Eq)

-- mkSpinBasis ::
--   HasCallStack =>
--   -- | Number of sites
--   Int ->
--   -- | Hamming weight
--   Maybe Int ->
--   -- | Spin inversion
--   Maybe Int ->
--   -- | Lattice symmetries
--   Symmetries ->
--   Basis 'SpinTy
-- mkSpinBasis n _ _ _
--   | n <= 0 = error $ "invalid number of spins: " <> show n
-- mkSpinBasis n (Just h) _ _
--   | h < 0 || h > n = error $ "invalid Hamming weight: " <> show h
-- mkSpinBasis _ _ (Just i) _
--   | i /= (-1) && i /= 1 = error $ "invalid spin inversion: " <> show i
-- mkSpinBasis n (Just h) (Just _) _
--   | n /= 2 * h = error $ "invalid spin inversion: " <> show n <> " spins, but Hamming weight is " <> show h
-- mkSpinBasis n _ _ g
--   | not (nullSymmetries g) && symmetriesGetNumberBits g /= n = error "invalid symmetries"
-- mkSpinBasis n h i g = basisFromHeader $ SpinHeader n h i g

-- mkSpinfulFermionicBasis ::
--   HasCallStack =>
--   -- | Number of sites
--   Int ->
--   -- | Number of particles
--   SpinfulOccupation ->
--   Basis 'SpinfulFermionTy
-- mkSpinfulFermionicBasis n _
--   | n <= 0 = error $ "invalid number of sites: " <> show n
-- mkSpinfulFermionicBasis n (SpinfulTotalParticles p)
--   | p < 0 || p > n = error $ "invalid number of particles: " <> show p
-- mkSpinfulFermionicBasis n (SpinfulPerSector u d)
--   | u < 0 || u > n || d < 0 || d > n =
--       error $
--         "invalid number of particles: " <> show u <> " up, " <> show d <> " down"
-- mkSpinfulFermionicBasis n occupation = basisFromHeader $ SpinfulFermionHeader n occupation

-- mkSpinlessFermionicBasis ::
--   HasCallStack =>
--   -- | Number of sites
--   Int ->
--   -- | Number of particles
--   Maybe Int ->
--   Basis 'SpinlessFermionTy
-- mkSpinlessFermionicBasis n _
--   | n <= 0 = error $ "invalid number of sites: " <> show n
-- mkSpinlessFermionicBasis n (Just h)
--   | h < 0 || h > n = error $ "invalid number of particles: " <> show h
-- mkSpinlessFermionicBasis n h = basisFromHeader $ SpinlessFermionHeader n h

-- ls_hs_destroy_string :: Ptr CChar -> IO ()
-- ls_hs_destroy_string = free

-- ls_hs_basis_to_json :: Ptr Cbasis -> IO CString
-- ls_hs_basis_to_json cBasis =
--   withReconstructedBasis cBasis $ \basis ->
--     newCString $ toStrict (Data.Aeson.encode basis)

-- ls_hs_basis_from_json :: CString -> IO (Ptr Cbasis)
-- ls_hs_basis_from_json cStr = handleAny (propagateErrorToC nullPtr) $ do
--   (basis :: SomeBasis) <- decodeCString cStr
--   foldSomeBasis borrowCbasis basis

-- ls_hs_create_spin_basis_from_json :: CString -> IO (Ptr Cbasis)
-- ls_hs_create_spin_basis_from_json cStr = do
--   (basis :: Basis 'SpinTy) <- decodeCString cStr
--   borrowCbasis basis

-- ls_hs_create_spin_basis_from_yaml :: CString -> IO (Ptr Cbasis)
-- ls_hs_create_spin_basis_from_yaml cFilename =
--   foldSomeBasis borrowCbasis =<< basisFromYAML =<< peekUtf8 cFilename

-- loadRawBasis :: MonadIO m => Text -> m (Ptr Cspin_basis)
-- loadRawBasis path = do
--   r <- liftIO $ decodeFileWithWarnings (toString path)
--   case r of
--     Left e -> throwIO e
--     Right (warnings, (WrappedBasisSpec basisSpec)) -> do
--       mapM_ print warnings
--       toRawBasis basisSpec

-- ls_hs_create_spinful_fermion_basis_from_json :: CString -> IO (Ptr Cbasis)
-- ls_hs_create_spinful_fermion_basis_from_json cStr = do
--   (basis :: Basis 'SpinfulFermionTy) <- decodeCString cStr
--   borrowCbasis basis

-- ls_hs_create_spinless_fermion_basis_from_json :: CString -> IO (Ptr Cbasis)
-- ls_hs_create_spinless_fermion_basis_from_json cStr = do
--   (basis :: Basis 'SpinlessFermionTy) <- decodeCString cStr
--   borrowCbasis basis

-- ls_hs_create_basis :: HasCallStack => Cparticle_type -> CInt -> CInt -> CInt -> IO (Ptr Cbasis)
-- ls_hs_create_basis particleType numberSites numberParticles numberUp
--   | particleType == c_LS_HS_SPIN =
--     let h = if numberUp == -1 then Nothing else Just (fromIntegral numberUp)
--         basis = mkSpinBasis (fromIntegral numberSites) h Nothing emptySymmetries
--      in borrowCbasis basis
--   | particleType == c_LS_HS_SPINFUL_FERMION =
--     let p
--           | numberParticles == -1 = SpinfulNoOccupation
--           | numberUp == -1 = SpinfulTotalParticles (fromIntegral numberParticles)
--           | otherwise = SpinfulPerSector (fromIntegral numberUp) (fromIntegral (numberParticles - numberUp))
--         basis = mkSpinfulFermionicBasis (fromIntegral numberSites) p
--      in borrowCbasis basis
--   | particleType == c_LS_HS_SPINLESS_FERMION =
--     let p = if numberParticles == -1 then Nothing else Just (fromIntegral numberParticles)
--         basis = mkSpinlessFermionicBasis (fromIntegral numberSites) p
--      in borrowCbasis basis
--   | otherwise = error $ "invalid particle type: " <> show particleType

-- newtype ChapelArray = ChapelArray (ForeignPtr Cexternal_array)

-- data Representatives = Representatives {rStates :: !ChapelArray}

binomial :: Int -> Int -> Maybe Int
binomial n k
  | n <= 0 || k > n = Just 0
  | otherwise = toIntegralSized $ factorial n `div` factorial (n - k) `div` factorial k
  where
    factorial :: Int -> Integer
    factorial x = product ([1 .. fromIntegral x] :: [Integer])

fixedHammingStateToIndex :: Word64 -> Int
fixedHammingStateToIndex = go 0 1
  where
    go !i !k !α
      | α /= 0 =
        let c = countTrailingZeros α
            α' = α .&. (α - 1)
            i' = i + Data.Maybe.fromJust (binomial c k)
         in go i' (k + 1) α'
      | otherwise = i

fixedHammingIndexToState :: Int -> Int -> Word64
fixedHammingIndexToState hammingWeight = go hammingWeight 0
  where
    go !i !state !index
      | i <= 0 = state
      | otherwise = go (i - 1) state' index'
      where
        (c, contribution) = inner index i (i - 1) 0
        -- (Data.Maybe.fromJust (binomial c i))
        state' = state .|. bit c
        index' = index - contribution
    inner !index !i !c !contribution
      | c >= numberBits = (c, contribution)
      | otherwise =
        if contribution' > index
          then (c, contribution)
          else inner index i c' contribution'
      where
        numberBits = 64
        c' = c + 1
        contribution' = Data.Maybe.fromJust (binomial c' i)

-- uint64_t
-- ls_hs_combinadics_index_to_state(ptrdiff_t index, int const hamming_weight,
--                                  ls_hs_combinadics_kernel_data const *cache) {
--   int const number_bits = cache->dimension - 1;
--   uint64_t state = 0;
--   for (int i = hamming_weight; i > 0; --i) {
--     // We are searching for the largest c such that
--     // binomial(c, i, cache) <= index
--     int c = i - 1;
--     ptrdiff_t contribution = (ptrdiff_t)binomial(c, i, cache);
--     while (c < number_bits) {
--       int const new_c = c + 1;
--       ptrdiff_t const new_contribution = (ptrdiff_t)binomial(new_c, i, cache);
--       if (new_contribution > index) {
--         break;
--       }
--       c = new_c;
--       contribution = new_contribution;
--     }
--
--     state |= ((uint64_t)1) << c;
--     index -= contribution;
--   }
--   return state;
-- }

-- ptrdiff_t
-- ls_hs_combinadics_state_to_index(uint64_t alpha,
--                                  ls_hs_combinadics_kernel_data const *cache) {
--   // fprintf(stderr, "rank_via_combinadics(%zu) = ", alpha);
--   ptrdiff_t i = 0;
--   for (int k = 1; alpha != 0; ++k) {
--     int c = __builtin_ctzl(alpha);
--     alpha &= alpha - 1;
--     // fprintf(stderr, "binomial(%i, %i) = %zi\n", c, k, binomial(c, k, cache));
--     i += (ptrdiff_t)binomial(c, k, cache);
--   }
--   // fprintf(stderr, "%zi\n", i);
--   return i;
-- }

-- foreign import ccall safe "ls_hs_build_representatives"
--   ls_hs_build_representatives :: Ptr Cbasis -> Word64 -> Word64 -> IO ()

-- basisBuild :: HasCallStack => Basis t -> IO ()
-- basisBuild basis
--   | getNumberBits (basisHeader basis) <= 64 =
--       withForeignPtr (basisContents basis) $ \basisPtr -> do
--         let (BasisState _ (BitString lower)) = minStateEstimate (basisHeader basis)
--             (BasisState _ (BitString upper)) = maxStateEstimate (basisHeader basis)
--         ls_hs_build_representatives basisPtr (fromIntegral lower) (fromIntegral upper)
--   | otherwise = withFrozenCallStack $ error "Too many bits"

-- stateIndex :: Basis t -> BasisState t -> Maybe Int
-- stateIndex basis (BasisState _ (BitString α))
--   | α <= fromIntegral (maxBound :: Word64) = unsafePerformIO . withCbasis basis $ \basisPtr ->
--       with (fromIntegral α :: Word64) $ \spinsPtr ->
--         alloca $ \indicesPtr -> do
--           (f, env) <- basisPeekStateIndexKernel basisPtr
--           if f == nullFunPtr
--             then pure Nothing
--             else do
--               mkCindex_kernel f 1 spinsPtr 1 indicesPtr 1 env
--               i <- peek indicesPtr
--               if i < 0
--                 then pure Nothing
--                 else pure $ Just (fromIntegral i)
--   | otherwise = Nothing

flattenIndex :: Basis t -> IndexType t -> Int
flattenIndex basis i = case basis of
  SpinBasis _ _ _ _ -> i
  SpinfulFermionBasis n _ ->
    case i of
      (SpinUp, k) -> k
      (SpinDown, k) -> n + k
  SpinlessFermionBasis _ _ -> i

unFlattenIndex :: Basis t -> Int -> IndexType t
unFlattenIndex basis i = case basis of
  SpinBasis _ _ _ _ -> i
  SpinlessFermionBasis _ _ -> i
  SpinfulFermionBasis n _
    | i > n -> (SpinDown, i - n)
    | otherwise -> (SpinUp, i)

getNumberSites :: Basis t -> Int
getNumberSites x = case x of
  SpinBasis n _ _ _ -> n
  SpinfulFermionBasis n _ -> n
  SpinlessFermionBasis n _ -> n
{-# INLINE getNumberSites #-}

getNumberBits :: Basis t -> Int
getNumberBits x = case x of
  SpinBasis n _ _ _ -> n
  SpinfulFermionBasis n _ -> 2 * n
  SpinlessFermionBasis n _ -> n
{-# INLINE getNumberBits #-}

getNumberWords :: Basis t -> Int
getNumberWords x = (n + 63) `div` 64
  where
    n = getNumberBits x
{-# INLINE getNumberWords #-}

getNumberParticles :: Basis t -> Maybe Int
getNumberParticles x = case x of
  SpinBasis n _ _ _ -> Just n
  SpinfulFermionBasis _ occupation ->
    case occupation of
      SpinfulNoOccupation -> Nothing
      SpinfulTotalParticles p -> Just p
      SpinfulPerSector u d -> Just (u + d)
  SpinlessFermionBasis _ p -> p
{-# INLINE getNumberParticles #-}

getNumberUp :: Basis t -> Maybe Int
getNumberUp x = case x of
  SpinBasis _ h _ _ -> h
  SpinfulFermionBasis _ occupation ->
    case occupation of
      SpinfulNoOccupation -> Nothing
      SpinfulTotalParticles _ -> Nothing
      SpinfulPerSector u _ -> Just u
  SpinlessFermionBasis _ _ -> Nothing
{-# INLINE getNumberUp #-}

getParticleType :: Basis t -> Cparticle_type
getParticleType x = case x of
  SpinBasis _ _ _ _ -> c_LS_HS_SPIN
  SpinfulFermionBasis _ _ -> c_LS_HS_SPINFUL_FERMION
  SpinlessFermionBasis _ _ -> c_LS_HS_SPINLESS_FERMION
{-# INLINE getParticleType #-}

getSpinInversion :: Basis t -> Maybe Int
getSpinInversion x = case x of
  SpinBasis _ _ i _ -> i
  SpinfulFermionBasis _ _ -> Nothing
  SpinlessFermionBasis _ _ -> Nothing
{-# INLINE getSpinInversion #-}

getParticleTag :: Basis t -> ParticleTag t
getParticleTag x = case x of
  SpinBasis _ _ _ _ -> SpinTag
  SpinfulFermionBasis _ _ -> SpinfulFermionTag
  SpinlessFermionBasis _ _ -> SpinlessFermionTag
{-# INLINE getParticleTag #-}

isStateIndexIdentity :: Basis t -> Bool
isStateIndexIdentity x = case x of
  SpinBasis _ Nothing Nothing g -> nullSymmetries g
  SpinfulFermionBasis _ SpinfulNoOccupation -> True
  SpinlessFermionBasis _ Nothing -> True
  _ -> False
{-# INLINE isStateIndexIdentity #-}

requiresProjection :: Basis t -> Bool
requiresProjection x = hasPermutationSymmetries x || hasSpinInversionSymmetry x
{-# INLINE requiresProjection #-}

hasSpinInversionSymmetry :: Basis t -> Bool
hasSpinInversionSymmetry x = case x of
  SpinBasis _ _ i _ -> isJust i
  SpinfulFermionBasis _ _ -> False
  SpinlessFermionBasis _ _ -> False
{-# INLINE hasSpinInversionSymmetry #-}

hasPermutationSymmetries :: Basis t -> Bool
hasPermutationSymmetries x = case x of
  SpinBasis _ _ _ g -> not (nullSymmetries g)
  SpinfulFermionBasis _ _ -> False
  SpinlessFermionBasis _ _ -> False
{-# INLINE hasPermutationSymmetries #-}

hasFixedHammingWeight :: Basis t -> Bool
hasFixedHammingWeight x = case x of
  SpinBasis _ (Just _) _ _ -> True
  SpinfulFermionBasis _ (SpinfulTotalParticles _) -> True
  SpinlessFermionBasis _ (Just _) -> True
  _ -> False
{-# INLINE hasFixedHammingWeight #-}

maxStateEstimate :: Basis t -> BasisState t
maxStateEstimate x = case x of
  SpinBasis n (Just h) i _ ->
    BasisState n . BitString $
      if isNothing i
        then (bit h - 1) `shiftL` (n - h)
        else (bit h - 1) `shiftL` (n - h - 1)
  -- TODO: improve the bound for fixed spin inversion
  SpinBasis n Nothing _ _ -> BasisState n . BitString $ bit n - 1
  SpinfulFermionBasis n SpinfulNoOccupation ->
    unsafeCastBasisState $
      maxStateEstimate (SpinlessFermionBasis (2 * n) Nothing)
  SpinfulFermionBasis n (SpinfulTotalParticles p) ->
    unsafeCastBasisState $
      maxStateEstimate (SpinlessFermionBasis (2 * n) (Just p))
  SpinfulFermionBasis n (SpinfulPerSector down up) ->
    let (BasisState _ minDown) = maxStateEstimate (SpinlessFermionBasis n (Just down))
        (BasisState _ minUp) = maxStateEstimate (SpinlessFermionBasis n (Just up))
     in BasisState (2 * n) $ (minUp `shiftL` n) .|. minDown
  SpinlessFermionBasis n p ->
    unsafeCastBasisState $
      maxStateEstimate (SpinBasis n p Nothing emptySymmetries)

minStateEstimate :: Basis t -> BasisState t
minStateEstimate x = case x of
  SpinBasis n (Just h) _ _ -> BasisState n . BitString $ bit h - 1
  -- TODO: improve the bound for fixed spin inversion
  SpinBasis n Nothing _ _ -> BasisState n . BitString $ zeroBits
  SpinfulFermionBasis n SpinfulNoOccupation ->
    unsafeCastBasisState $
      minStateEstimate (SpinlessFermionBasis (2 * n) Nothing)
  SpinfulFermionBasis n (SpinfulTotalParticles p) ->
    unsafeCastBasisState $
      minStateEstimate (SpinlessFermionBasis (2 * n) (Just p))
  SpinfulFermionBasis n (SpinfulPerSector down up) ->
    let (BasisState _ minDown) = minStateEstimate (SpinlessFermionBasis n (Just down))
        (BasisState _ minUp) = minStateEstimate (SpinlessFermionBasis n (Just up))
     in BasisState (2 * n) $ (minUp `shiftL` n) .|. minDown
  SpinlessFermionBasis n p ->
    unsafeCastBasisState $
      minStateEstimate (SpinBasis n p Nothing emptySymmetries)

isBasisReal :: Basis t -> Bool
isBasisReal x = case x of
  SpinBasis _ _ _ g -> areSymmetriesReal g
  SpinfulFermionBasis _ _ -> True
  SpinlessFermionBasis _ _ -> True

optionalNatural :: Maybe Int -> CInt
optionalNatural x = case x of
  (Just n) -> fromIntegral n
  Nothing -> (-1)
{-# INLINE optionalNatural #-}

newCbasis :: HasCallStack => IsBasis t => Basis t -> IO (Ptr Cbasis)
newCbasis x = do
  -- logDebug' $ "basisFromHeader" <> show x
  -- fp <- mallocForeignPtr
  kernels <- createCbasis_kernels x
  payload <- newStablePtr (SomeBasis x)
  -- addForeignPtrConcFinalizer fp (destroyCbasis_kernels kernels)
  -- withForeignPtr fp $ \ptr ->
  --   poke ptr $
  new $
    Cbasis
      { cbasis_refcount = 1,
        cbasis_number_sites = fromIntegral (getNumberSites x),
        cbasis_number_particles = optionalNatural (getNumberParticles x),
        cbasis_number_up = optionalNatural (getNumberUp x),
        cbasis_particle_type = getParticleType x,
        cbasis_spin_inversion = maybe 0 fromIntegral (getSpinInversion x),
        cbasis_state_index_is_identity = fromBool (isStateIndexIdentity x),
        cbasis_requires_projection = fromBool (requiresProjection x),
        cbasis_kernels = kernels,
        cbasis_representatives = emptyExternalArray,
        cbasis_haskell_payload = castStablePtrToPtr payload
      }

foreign import ccall unsafe "ls_hs_destroy_external_array"
  ls_hs_destroy_external_array :: Ptr Cexternal_array -> IO ()

destroyCbasis :: HasCallStack => Ptr Cbasis -> IO ()
destroyCbasis p = do
  logDebug' $ "ls_hs_destroy_basis " <> show p
  refcount <- basisDecRefCount p
  when (refcount == 1) $ do
    logDebug' $ "fully destroying " <> show p
    x <- peek p
    with (cbasis_representatives x) ls_hs_destroy_external_array
    destroyCbasis_kernels (cbasis_kernels x)
    freeStablePtr . castPtrToStablePtr $ cbasis_haskell_payload x

-- pure $ Basis x fp

-- basisFromHeader :: BasisHeader t -> Basis t
-- basisFromHeader x = unsafePerformIO $ do
--   -- logDebug' $ "basisFromHeader" <> show x
--   fp <- mallocForeignPtr
--   kernels <- createCbasis_kernels x
--   addForeignPtrConcFinalizer fp (destroyCbasis_kernels kernels)
--   withForeignPtr fp $ \ptr ->
--     poke ptr $
--       Cbasis
--         { cbasis_refcount = 0,
--           cbasis_number_sites = fromIntegral (getNumberSites x),
--           cbasis_number_particles = optionalNatural (getNumberParticles x),
--           cbasis_number_up = optionalNatural (getNumberUp x),
--           cbasis_particle_type = getParticleType x,
--           cbasis_state_index_is_identity = fromBool (isStateIndexIdentity x),
--           cbasis_requires_projection = fromBool (requiresProjection x),
--           cbasis_kernels = kernels,
--           cbasis_representatives = emptyExternalArray,
--           cbasis_haskell_payload = nullPtr
--         }
--   pure $ Basis x fp
-- {-# NOINLINE basisFromHeader #-}

-- borrowCbasis :: HasCallStack => Basis t -> IO (Ptr Cbasis)
-- borrowCbasis basis = do
--   withForeignPtr (basisContents basis) $ \ptr -> do
--     logDebug' $ "borrowCbasis " <> show ptr
--     _ <- basisIncRefCount ptr
--     basisPokePayload ptr =<< newStablePtr basis
--   pure $ unsafeForeignPtrToPtr (basisContents basis)

cloneCbasis :: HasCallStack => Ptr Cbasis -> IO (Ptr Cbasis)
cloneCbasis ptr = do
  logDebug' $ "ls_hs_clone_basis " <> show ptr
  _ <- basisIncRefCount ptr
  pure ptr

withCbasis :: Ptr Cbasis -> (SomeBasis -> IO a) -> IO a
withCbasis p f = f =<< deRefStablePtr . castPtrToStablePtr . cbasis_haskell_payload =<< peek p

-- withForeignPtr (basisContents x) action
-- withCbasis :: Basis t -> (Ptr Cbasis -> IO a) -> IO a
-- withCbasis x action = withForeignPtr (basisContents x) action

-- withParticleType ::
--   HasCallStack =>
--   Cparticle_type ->
--   (forall (t :: ParticleTy). IsBasis t => Proxy t -> a) ->
--   a
-- withParticleType t action
--   | t == c_LS_HS_SPIN = action (Proxy :: Proxy 'SpinTy)
--   | t == c_LS_HS_SPINFUL_FERMION = action (Proxy :: Proxy 'SpinfulFermionTy)
--   | t == c_LS_HS_SPINLESS_FERMION = action (Proxy :: Proxy 'SpinlessFermionTy)
--   | otherwise = error $ "invalid particle type: " <> show t

-- withReconstructedBasis :: forall a. HasCallStack => Ptr Cbasis -> (forall (t :: ParticleTy). IsBasis t => Basis t -> IO a) -> IO a
-- withReconstructedBasis p action = do
--   logDebug' $ "withReconstructedBasis " <> show p
--
--   let run :: forall (t :: ParticleTy). IsBasis t => Proxy t -> IO a
--       run _ = do
--         (x :: Basis t) <- deRefStablePtr =<< basisPeekPayload p
--         action x
--   tp <- basisPeekParticleType p
--   withParticleType tp run

-- data Ccombinadics_kernel_data

data Cbinary_search_kernel_data

-- foreign import ccall safe "lattice_symmetries_haskell.h ls_hs_internal_create_combinadics_kernel_data"
--   ls_hs_internal_create_combinadics_kernel_data :: CInt -> CBool -> IO (Ptr Ccombinadics_kernel_data)

-- foreign import ccall safe "lattice_symmetries_haskell.h ls_hs_internal_destroy_combinadics_kernel_data"
--   ls_hs_internal_destroy_combinadics_kernel_data :: Ptr Ccombinadics_kernel_data -> IO ()

-- foreign import ccall safe "lattice_symmetries_haskell.h &ls_hs_state_index_combinadics_kernel"
--   ls_hs_state_index_combinadics_kernel :: FunPtr Cindex_kernel

foreign import ccall safe "lattice_symmetries_haskell.h ls_hs_destroy_state_index_binary_search_kernel_data"
  ls_hs_destroy_state_index_binary_search_kernel_data :: Ptr Cbinary_search_kernel_data -> IO ()

foreign import ccall safe "lattice_symmetries_haskell.h &ls_hs_state_index_binary_search_kernel"
  ls_hs_state_index_binary_search_kernel :: FunPtr Cindex_kernel

-- foreign import ccall safe "lattice_symmetries_haskell.h &ls_hs_state_index_identity_kernel"
--   ls_hs_state_index_identity_kernel :: FunPtr Cindex_kernel

data Chalide_kernel_data

foreign import ccall safe "lattice_symmetries_haskell.h ls_internal_create_halide_kernel_data"
  ls_internal_create_halide_kernel_data :: Ptr Cpermutation_group -> CInt -> IO (Ptr Chalide_kernel_data)

foreign import ccall safe "lattice_symmetries_haskell.h ls_internal_destroy_halide_kernel_data"
  ls_internal_destroy_halide_kernel_data :: Ptr Chalide_kernel_data -> IO ()

foreign import ccall safe "lattice_symmetries_haskell.h &ls_hs_is_representative_halide_kernel"
  ls_hs_is_representative_halide_kernel :: FunPtr Cis_representative_kernel

foreign import ccall safe "lattice_symmetries_haskell.h &ls_hs_state_info_halide_kernel"
  ls_hs_state_info_halide_kernel :: FunPtr Cstate_info_kernel

setStateInfoKernel :: Basis t -> Cbasis_kernels -> IO Cbasis_kernels
setStateInfoKernel (SpinBasis n _ i g) k
  | n <= 64 = do
    kernelData <-
      bracket (newCpermutation_group g) destroyCpermutation_group $ \gPtr ->
        -- withForeignPtr (symmetriesContents g) $ \gPtr ->
        ls_internal_create_halide_kernel_data gPtr (maybe 0 fromIntegral i)
    pure $
      k
        { cbasis_state_info_kernel = ls_hs_state_info_halide_kernel,
          cbasis_state_info_data = castPtr kernelData,
          cbasis_is_representative_kernel = ls_hs_is_representative_halide_kernel,
          cbasis_is_representative_data = castPtr kernelData
        }
setStateInfoKernel _ k =
  pure $
    k
      { cbasis_state_info_kernel = nullFunPtr,
        cbasis_state_info_data = nullPtr,
        cbasis_is_representative_kernel = nullFunPtr,
        cbasis_is_representative_data = nullPtr
      }

setStateIndexKernel :: Basis t -> Cbasis_kernels -> IO Cbasis_kernels
setStateIndexKernel x k =
  pure $
    k
      { cbasis_state_index_kernel = nullFunPtr,
        cbasis_state_index_data = nullPtr
      }

-- \| isStateIndexIdentity x =
--     pure $
--       k
--         { cbasis_state_index_kernel = ls_hs_state_index_identity_kernel,
--           cbasis_state_index_data = nullPtr
--         }
-- \| hasFixedHammingWeight x && not (requiresProjection x) = do
--     kernelData <-
--       ls_hs_internal_create_combinadics_kernel_data
--         (fromIntegral (getNumberBits x))
--         (fromBool False)
--     pure $
--       k
--         { cbasis_state_index_kernel = ls_hs_state_index_combinadics_kernel,
--           cbasis_state_index_data = castPtr kernelData
--         }
-- \| otherwise = case x of
--     (SpinfulFermionHeader n (SpinfulPerSector _ _)) -> do
--       kernelData <-
--         ls_hs_internal_create_combinadics_kernel_data
--           (fromIntegral n) -- NOTE: not 2 * n !
--           (fromBool True)
--       pure $
--         k
--           { cbasis_state_index_kernel = ls_hs_state_index_combinadics_kernel,
--             cbasis_state_index_data = castPtr kernelData
--           }
--     _ -> error "not implemented"

createCbasis_kernels :: Basis t -> IO (Ptr Cbasis_kernels)
createCbasis_kernels x =
  new <=< setStateInfoKernel x <=< setStateIndexKernel x $
    Cbasis_kernels
      { cbasis_state_info_kernel = nullFunPtr,
        cbasis_state_info_data = nullPtr,
        cbasis_is_representative_kernel = nullFunPtr,
        cbasis_is_representative_data = nullPtr,
        cbasis_state_index_kernel = nullFunPtr,
        cbasis_state_index_data = nullPtr
      }

destroyCindex_kernel :: HasCallStack => Cbasis_kernels -> IO ()
destroyCindex_kernel p
  -- \| kernel == ls_hs_state_index_combinadics_kernel = ls_hs_internal_destroy_combinadics_kernel_data (castPtr env)
  | kernel == ls_hs_state_index_binary_search_kernel = ls_hs_destroy_state_index_binary_search_kernel_data (castPtr env)
  -- \| kernel == ls_hs_state_index_identity_kernel && env == nullPtr = pure ()
  | kernel == nullFunPtr && env == nullPtr = pure ()
  | otherwise = error "failed to automatically deallocate state_index kernel"
  where
    kernel = cbasis_state_index_kernel p
    env = cbasis_state_index_data p

destroyCstate_info_kernel :: HasCallStack => Cbasis_kernels -> IO ()
destroyCstate_info_kernel p
  | cbasis_state_info_kernel p == ls_hs_state_info_halide_kernel
      && cbasis_is_representative_kernel p == ls_hs_is_representative_halide_kernel =
    do
      ls_internal_destroy_halide_kernel_data (castPtr (cbasis_state_info_data p))
      when (cbasis_state_info_data p /= cbasis_is_representative_data p) $
        ls_internal_destroy_halide_kernel_data (castPtr (cbasis_is_representative_data p))
  | cbasis_state_info_kernel p == nullFunPtr
      && cbasis_state_info_data p == nullPtr
      && cbasis_is_representative_kernel p == nullFunPtr
      && cbasis_is_representative_data p == nullPtr =
    pure ()
  | otherwise = error "failed to automatically deallocate state_info and is_representative kernels"

destroyCbasis_kernels :: HasCallStack => Ptr Cbasis_kernels -> IO ()
destroyCbasis_kernels p = do
  kernels <- peek p
  destroyCindex_kernel kernels
  destroyCstate_info_kernel kernels
  free p

-- data BasisSpec = BasisSpec !Int !(Maybe Int) !(Maybe Int) ![SymmetrySpec]
--   deriving stock (Read, Show, Eq)
--
-- instance FromJSON BasisSpec where
--   parseJSON = withObject "basis" $ \v ->
--     BasisSpec
--       <$> v .: "number_spins"
--       <*> v .:? "hamming_weight"
--       <*> v .:? "spin_inversion"
--       <*> v .:! "symmetries" .!= []
--
-- instance ToJSON BasisSpec where
--   toJSON (BasisSpec numberSpins hammingWeight spinInversion symmetries) =
--     object
--       [ "number_spins" .= numberSpins,
--         "hamming_weight" .= maybe Null toJSON hammingWeight,
--         "spin_inversion" .= maybe Null toJSON spinInversion,
--         "symmetries" .= symmetries
--       ]

-- basisFromSpec :: BasisSpec -> Basis 'SpinTy
-- basisFromSpec (BasisSpec n h i gs) = mkSpinBasis n h i (symmetriesFromSpec gs)

foreign export ccall "ls_hs_basis_number_bits"
  ls_hs_basis_number_bits :: Ptr Cbasis -> IO CInt

ls_hs_basis_number_bits basisPtr =
  fromIntegral <$> withCbasis basisPtr (foldSomeBasis (pure . getNumberBits))
