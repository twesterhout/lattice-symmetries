{-# LANGUAGE DuplicateRecordFields #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE OverloadedLists #-}
{-# LANGUAGE OverloadedRecordDot #-}
{-# LANGUAGE QuantifiedConstraints #-}
{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE TypeOperators #-}

module LatticeSymmetries.Basis
  ( -- * High-level interface
    BasisState (..)
  , Basis (..)
  -- BasisHeader (..),
  , SpinfulOccupation (..)
  , SomeBasis (..)
  -- basisBuild,
  , IndexType
  , GeneratorType
  , IsBasis
  , withSomeBasis
  , foldSomeBasis

    -- ** Creating bases

  -- mkSpinBasis,
  -- mkSpinlessFermionicBasis,
  -- mkSpinfulFermionicBasis,

    -- ** Querying information
  , getNumberBits
  , getNumberWords
  , getNumberSites
  , getHammingWeight
  , minStateEstimate
  , maxStateEstimate
  , isStateIndexIdentity
  , hasFixedHammingWeight
  , hasSpinInversionSymmetry
  , hasPermutationSymmetries
  , requiresProjection
  , getParticleTag
  , isBasisReal

    -- ** Low-level interface
  , newCbasis
  , ls_hs_basis_from_json
  , ls_hs_basis_to_json
  , ls_hs_destroy_basis
  , ls_hs_init_basis_info
  , ls_hs_init_is_representative_kernel
  , ls_hs_init_state_info_kernel
  , ls_hs_init_state_to_index_kernel
  -- , Cbasis (..)
  -- basisFromYAML,
  -- objectFromYAML,
  -- Factor,
  -- borrowCbasis,
  -- , newCbasis
  -- , cloneCbasis
  -- , destroyCbasis
  , withCbasis
  -- Cparticle_type (..),
  -- Cbasis_kernels (..),
  -- createCbasis_kernels,
  -- destroyCbasis_kernels,
  -- stateIndex,
  -- , flattenIndex
  -- , unFlattenIndex
  , matchParticleType2
  -- withParticleType,
  -- withReconstructedBasis,
  )
where

import Data.Aeson
import Data.Aeson.Types (Pair, Parser)
import Data.Bits
import Data.Maybe (fromJust)
import Data.Vector.Generic qualified as G
import Data.Vector.Storable qualified as S
import Foreign (fromBool)
import Foreign.C (CString)
import Foreign.C.Types
import Foreign.Marshal.Alloc
import Foreign.Ptr
import Foreign.StablePtr
import Language.C.Inline.Unsafe qualified as CU
import LatticeSymmetries.Algebra
import LatticeSymmetries.BitString
import LatticeSymmetries.Context
import LatticeSymmetries.FFI
import LatticeSymmetries.Generator
import LatticeSymmetries.Group
import LatticeSymmetries.Lowering
import LatticeSymmetries.Permutation (Permutation)
import LatticeSymmetries.Utils
import Text.Read (read)
import Prelude hiding (state)

importLS

data Basis (t :: ParticleTy) where
  SpinBasis :: !Int -> !(Maybe Int) -> !(Maybe Int) -> !(Representation Permutation) -> Basis 'SpinTy
  SpinfulFermionBasis :: !Int -> !SpinfulOccupation -> Basis 'SpinfulFermionTy
  SpinlessFermionBasis :: !Int -> !(Maybe Int) -> Basis 'SpinlessFermionTy

deriving stock instance Show (Basis t)

deriving stock instance Eq (Basis t)

data SomeBasis where
  SomeBasis :: IsBasis t => Basis t -> SomeBasis

deriving stock instance Show SomeBasis

instance Eq SomeBasis where
  (==) (SomeBasis a@(SpinBasis {})) (SomeBasis b@(SpinBasis {})) = a == b
  (==) (SomeBasis a@(SpinfulFermionBasis {})) (SomeBasis b@(SpinfulFermionBasis {})) = a == b
  (==) (SomeBasis a@(SpinlessFermionBasis {})) (SomeBasis b@(SpinlessFermionBasis {})) = a == b
  (==) _ _ = False

withSomeBasis :: SomeBasis -> (forall t. IsBasis t => Basis t -> a) -> a
withSomeBasis x f = case x of SomeBasis basis -> f basis
{-# INLINE withSomeBasis #-}

foldSomeBasis :: (forall t. IsBasis t => Basis t -> a) -> SomeBasis -> a
foldSomeBasis f x = withSomeBasis x f
{-# INLINE foldSomeBasis #-}

-- mapSomeBasisM :: Monad m => (forall t. IsBasis t => Basis t -> m (Basis t)) -> SomeBasis -> m SomeBasis
-- mapSomeBasisM f x = case x of SomeBasis basis -> SomeBasis <$> f basis
-- {-# INLINE mapSomeBasisM #-}

-- mapSomeBasis :: (forall t. IsBasis t => Basis t -> Basis t) -> SomeBasis -> SomeBasis
-- mapSomeBasis f x = case x of SomeBasis basis -> SomeBasis (f basis)
-- {-# INLINE mapSomeBasis #-}

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
      Number _ -> SpinfulTotalParticles <$> parseJSON n
      Array _ -> uncurry SpinfulPerSector <$> parseJSON n
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
        [ "particle" .= SpinTy
        , "number_spins" .= n
        , "hamming_weight" .= h
        , "spin_inversion" .= i
        , "symmetries" .= unRepresentation g
        ]
  | (SpinfulFermionBasis n o) <- x =
      object $
        ["particle" .= SpinfulFermionTy, "number_sites" .= n]
          <> encodeSpinfulOccupation o
  | (SpinlessFermionBasis n o) <- x =
      object
        ["particle" .= SpinlessFermionTy, "number_sites" .= n, "number_particles" .= o]

basisHeaderFromJSON :: Value -> Parser SomeBasis
basisHeaderFromJSON = withObject "Basis" $ \v -> do
  tp <- v .:! "particle" .!= SpinTy
  case tp of
    SpinTy -> do
      g <- eitherToParser . fromGenerators =<< (v .:! "symmetries" .!= [])
      fmap SomeBasis
        . join
        $ mkSpinHeader
          <$> (v .: "number_spins")
          <*> (v .:? "hamming_weight")
          <*> (v .:? "spin_inversion")
          <*> pure g
    SpinlessFermionTy ->
      fmap SomeBasis $
        SpinlessFermionBasis
          <$> (v .: "number_sites")
          <*> (v .:? "number_particles")
    SpinfulFermionTy ->
      fmap SomeBasis $
        SpinfulFermionBasis
          <$> (v .: "number_sites")
          <*> parseSpinfulOccupation v

mkSpinHeader
  :: MonadFail m
  => Int
  -- ^ Number of sites
  -> Maybe Int
  -- ^ Hamming weight
  -> Maybe Int
  -- ^ Spin inversion
  -> Representation Permutation
  -- ^ Lattice symmetries
  -> m (Basis 'SpinTy)
mkSpinHeader n _ _ _
  | n <= 0 = fail $ "invalid number of spins: " <> show n
mkSpinHeader n (Just h) _ _
  | h < 0 || h > n = fail $ "invalid Hamming weight: " <> show h
mkSpinHeader _ _ (Just i) _
  | i /= (-1) && i /= 1 = fail $ "invalid spin inversion: " <> show i
mkSpinHeader n (Just h) (Just k) _
  | n /= 2 * h = fail $ "invalid spin inversion: " <> show k <> "; " <> show n <> " spins, but the Hamming weight is " <> show h
mkSpinHeader n _ _ g
  | maybe False (/= n) g.numberBits = fail "invalid symmetries"
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

-- flattenIndex :: Basis t -> IndexType t -> Int
-- flattenIndex basis i = case basis of
--   SpinBasis {} -> i
--   SpinfulFermionBasis n _ ->
--     case i of
--       (SpinUp, k) -> k
--       (SpinDown, k) -> n + k
--   SpinlessFermionBasis _ _ -> i

-- unFlattenIndex :: Basis t -> Int -> IndexType t
-- unFlattenIndex basis i = case basis of
--   SpinBasis {} -> i
--   SpinlessFermionBasis {} -> i
--   SpinfulFermionBasis n _
--     | i >= n -> (SpinDown, i - n)
--     | otherwise -> (SpinUp, i)

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

-- getParticleType :: Basis t -> Cparticle_type
-- getParticleType x = case x of
--   SpinBasis {} -> c_LS_HS_SPIN
--   SpinfulFermionBasis {} -> c_LS_HS_SPINFUL_FERMION
--   SpinlessFermionBasis {} -> c_LS_HS_SPINLESS_FERMION
-- {-# INLINE getParticleType #-}

getSpinInversion :: Basis t -> Maybe Int
getSpinInversion x = case x of
  SpinBasis _ _ i _ -> i
  SpinfulFermionBasis _ _ -> Nothing
  SpinlessFermionBasis _ _ -> Nothing
{-# INLINE getSpinInversion #-}

getParticleTag :: Basis t -> ParticleTag t
getParticleTag x = case x of
  SpinBasis {} -> SpinTag
  SpinfulFermionBasis {} -> SpinfulFermionTag
  SpinlessFermionBasis {} -> SpinlessFermionTag
{-# INLINE getParticleTag #-}

isStateIndexIdentity :: Basis t -> Bool
isStateIndexIdentity x = case x of
  -- TODO: only fixed Hamming weight or lattice symmetries break the trivial index<->state mapping
  -- Setting spin inversion limits the domain to representatives, but keeps the index<->state identity.
  SpinBasis _ Nothing _ g -> G.null (unRepresentation g)
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
  _ -> False
{-# INLINE hasSpinInversionSymmetry #-}

hasPermutationSymmetries :: Basis t -> Bool
hasPermutationSymmetries x = case getPermutationSymmetries x of
  Nothing -> False
  Just g -> not $ G.null (unRepresentation g)
{-# INLINE hasPermutationSymmetries #-}

getPermutationSymmetries :: Basis t -> Maybe (Representation Permutation)
getPermutationSymmetries = \case
  SpinBasis _ _ _ g -> Just g
  _ -> Nothing
{-# INLINE getPermutationSymmetries #-}

hasFixedHammingWeight :: Basis t -> Bool
hasFixedHammingWeight x = case x of
  SpinBasis _ (Just _) _ _ -> True
  SpinfulFermionBasis _ (SpinfulTotalParticles _) -> True
  SpinfulFermionBasis _ (SpinfulPerSector _ _) -> True
  SpinlessFermionBasis _ (Just _) -> True
  _ -> False
{-# INLINE hasFixedHammingWeight #-}

getHammingWeight :: Basis t -> Maybe Int
getHammingWeight x = case x of
  SpinBasis _ h _ _ -> h
  SpinfulFermionBasis _ (SpinfulTotalParticles n) -> Just n
  SpinfulFermionBasis _ (SpinfulPerSector u d) -> Just (u + d)
  SpinlessFermionBasis _ (Just n) -> Just n
  _ -> Nothing

maxStateEstimate :: Basis t -> BasisState t
maxStateEstimate x = case x of
  SpinBasis n (Just h) i _ ->
    BasisState n
      . BitString
      $ if isNothing i
        then (bit h - 1) `shiftL` (n - h)
        else (bit h - 1) `shiftL` (n - h - 1)
  SpinBasis n Nothing (Just _) _ -> BasisState n . BitString $ bit (n - 1) - 1
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
      maxStateEstimate (SpinBasis n p Nothing (Representation G.empty))

minStateEstimate :: Basis t -> BasisState t
minStateEstimate x = case x of
  SpinBasis n (Just h) _ _ -> BasisState n . BitString $ bit h - 1
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
      minStateEstimate (SpinBasis n p Nothing (Representation G.empty))

isBasisReal :: Basis t -> Bool
isBasisReal x = case x of
  SpinBasis _ _ _ g -> isRepresentationReal g
  SpinfulFermionBasis _ _ -> True
  SpinlessFermionBasis _ _ -> True

-- optionalNatural :: Maybe Int -> CInt
-- optionalNatural x = case x of
--   (Just n) -> fromIntegral n
--   Nothing -> (-1)
-- {-# INLINE optionalNatural #-}

-- newCbasis :: HasCallStack => IsBasis t => Basis t -> IO (Ptr Cbasis)
-- newCbasis x = do
--   -- logDebug' $ "basisFromHeader" <> show x
--   -- fp <- mallocForeignPtr
--   kernels <- createCbasis_kernels x
--   payload <- newStablePtr (SomeBasis x)
--   -- addForeignPtrConcFinalizer fp (destroyCbasis_kernels kernels)
--   -- withForeignPtr fp $ \ptr ->
--   --   poke ptr $
--   new
--     $ Cbasis
--       { cbasis_refcount = AtomicCInt 1
--       , cbasis_number_sites = fromIntegral (getNumberSites x)
--       , cbasis_number_particles = optionalNatural (getNumberParticles x)
--       , cbasis_number_up = optionalNatural (getNumberUp x)
--       , cbasis_particle_type = getParticleType x
--       , cbasis_spin_inversion = maybe 0 fromIntegral (getSpinInversion x)
--       , cbasis_state_index_is_identity = fromBool (isStateIndexIdentity x)
--       , cbasis_requires_projection = fromBool (requiresProjection x)
--       , cbasis_kernels = kernels
--       , cbasis_representatives = emptyExternalArray
--       , cbasis_haskell_payload = castStablePtrToPtr payload
--       }

-- foreign import ccall unsafe "ls_hs_internal_destroy_external_array"
--   ls_hs_internal_destroy_external_array :: Ptr Cexternal_array -> IO ()

-- void ls_hs_destroy_basis(ls_hs_basis *basis) {
--     LS_CHECK(basis != NULL, "trying to destroy a NULL basis");
--     if (ls_hs_internal_object_dec_ref_count(&basis->base) == 1) {
--         ls_hs_internal_destroy_external_array(&basis->local_representatives);
--         ls_hs_internal_destroy_external_array(&basis->norms);
--         ls_hs_destroy_is_representative_kernel(basis->is_representative_kernel);
--         ls_hs_destroy_state_info_kernel(basis->is_representative_kernel);
--         hs_free_stable_ptr(basis->haskell_payload);
--     }
-- }

ls_hs_destroy_basis :: Ptr Cbasis -> IO ()
ls_hs_destroy_basis = ls_hs_destroy_object $ \p -> do
  [CU.block| void {
    ls_hs_basis* p = $(ls_hs_basis* p);

    LS_HS_INTERNAL_DESTROY_KERNEL(ls_hs_is_representative_kernel_type_v2,
                                  p->is_representative_kernel,
                                  p->is_representative_destructor);
    LS_HS_INTERNAL_DESTROY_KERNEL(ls_hs_state_info_kernel_type_v2,
                                  p->state_info_kernel,
                                  p->state_info_destructor);
    LS_HS_INTERNAL_DESTROY_KERNEL(ls_hs_state_to_index_kernel_type,
                                  p->state_to_index_kernel,
                                  p->state_to_index_destructor);

    ls_hs_internal_destroy_external_array(&p->local_representatives);
    ls_hs_internal_destroy_external_array(&p->local_norms);

    if (p->info != NULL) {
      if (p->info->characters != NULL) {
        // remove const qualifiers
        free((void *)p->info->characters);
      }
      free(p->info);
    }
  } |]

ls_hs_init_is_representative_kernel :: Ptr Cbasis -> IO ()
ls_hs_init_is_representative_kernel basisPtr = withCbasis basisPtr $ foldSomeBasis $ \basis ->
  if hasPermutationSymmetries basis
    then do
      kernel <- createIsRepresentativeKernel_v2 (fromJust $ getPermutationSymmetries basis) (getSpinInversion basis)
      [CU.block| void {
        ls_hs_basis* const basisPtr = $(ls_hs_basis* basisPtr);
        LS_HS_INTERNAL_INIT_KERNEL(ls_hs_is_representative_kernel_type_v2,
                                   &basisPtr->is_representative_kernel,
                                   &basisPtr->is_representative_destructor,
                                   $(ls_hs_is_representative_kernel_type_v2 kernel),
                                   &ls_hs_internal_destroy_is_representative_kernel);
      } |]
    else error "do not call ls_hs_init_is_representative_kernel on a basis without permutation symmetries"

ls_hs_init_state_info_kernel :: Ptr Cbasis -> IO ()
ls_hs_init_state_info_kernel basisPtr = withCbasis basisPtr $ foldSomeBasis $ \basis ->
  if hasPermutationSymmetries basis
    then do
      kernel <- createStateInfoKernel_v2 (fromJust $ getPermutationSymmetries basis) (getSpinInversion basis)
      [CU.block| void {
        ls_hs_basis* const basisPtr = $(ls_hs_basis* basisPtr);
        LS_HS_INTERNAL_INIT_KERNEL(ls_hs_state_info_kernel_type_v2,
                                   &basisPtr->state_info_kernel,
                                   &basisPtr->state_info_destructor,
                                   $(ls_hs_state_info_kernel_type_v2 kernel),
                                   &ls_hs_internal_destroy_state_info_kernel);
      } |]
    else error "do not call ls_hs_init_state_info_kernel on a basis without permutation symmetries"

isBasisBuilt :: Ptr Cbasis -> IO Bool
isBasisBuilt basis = toEnum . fromIntegral <$> [CU.exp| bool { $(ls_hs_basis const* basis)->local_representatives.elts != NULL } |]

ls_hs_init_state_to_index_kernel :: Ptr Cbasis -> IO ()
ls_hs_init_state_to_index_kernel basisPtr = withCbasis basisPtr $ foldSomeBasis $ \basis ->
  if isJust (getHammingWeight basis) && not (hasPermutationSymmetries basis)
    then do
      kernel <- createFixedHammingStateToIndexKernel (getNumberBits basis) (fromJust $ getHammingWeight basis)
      [CU.block| void {
        ls_hs_basis* const basisPtr = $(ls_hs_basis* basisPtr);
        LS_HS_INTERNAL_INIT_KERNEL(ls_hs_state_to_index_kernel_type,
                                   &basisPtr->state_to_index_kernel,
                                   &basisPtr->state_to_index_destructor,
                                   $(ls_hs_state_to_index_kernel_type kernel),
                                   &ls_hs_internal_destroy_fixed_hamming_state_to_index_kernel);
      } |]
    else
      isBasisBuilt basisPtr >>= \case
        True -> do
          let numberBits = fromIntegral (getNumberBits basis)
          prefixBits <- maybe 26 read <$> lookupEnv "LS_HS_CACHE_NUM_BITS"
          [CU.block| void {
            ls_hs_basis* const basisPtr = $(ls_hs_basis* basisPtr);
            ls_hs_state_to_index_kernel_type const kernel =
              ls_hs_internal_mk_binary_search_state_to_index_kernel((int64_t)basisPtr->local_representatives.num_elts,
                                                                    (uint64_t const*)basisPtr->local_representatives.elts,
                                                                    $(unsigned numberBits),
                                                                    $(unsigned prefixBits));
            LS_HS_INTERNAL_INIT_KERNEL(ls_hs_state_to_index_kernel_type,
                                       &basisPtr->state_to_index_kernel,
                                       &basisPtr->state_to_index_destructor,
                                       kernel,
                                       &ls_hs_internal_destroy_binary_search_state_to_index_kernel);
          } |]
        False -> error "cannot initialize the binary_search_state_to_index_kernel because the basis has not been built"

newCbasis :: IsBasis t => Basis t -> IO (Ptr Cbasis)
newCbasis x = do
  payload <- castStablePtrToPtr <$> newStablePtr (SomeBasis x)
  -- NOTE: important to initialize memory to 0 such that we don't have to manually initialize fields
  p <- callocBytes $ fromIntegral [CU.pure| size_t { sizeof(ls_hs_basis) } |]
  [CU.block| void { ls_hs_internal_object_init(&$(ls_hs_basis* p)->base, 1, $(void* payload)); } |]
  pure p

ls_hs_init_basis_info :: Ptr Cbasis -> IO ()
ls_hs_init_basis_info p =
  withCbasis p $ foldSomeBasis $ \basis -> do
    let has_permutation_symmetries = fromBool (hasPermutationSymmetries basis)
        requires_projection = fromBool (requiresProjection basis)
        is_state_index_identity = fromBool (isStateIndexIdentity basis)
        is_real = fromBool (isBasisReal basis)
        number_bits = fromIntegral (getNumberBits basis)
        number_words = fromIntegral (getNumberWords basis)
        number_sites = fromIntegral (getNumberSites basis)
        number_particles = maybe (-1) fromIntegral (getNumberParticles basis)
        hamming_weight = maybe (-1) fromIntegral (getHammingWeight basis)
        number_up = maybe (-1) fromIntegral (getNumberUp basis)
        spin_inversion = maybe 0 fromIntegral (getSpinInversion basis)
        min_state_estimate = case minStateEstimate basis of
          BasisState k (BitString n) -> if k <= 64 then fromIntegral n else minBound
        max_state_estimate = case maxStateEstimate basis of
          BasisState k (BitString n) -> if k <= 64 then fromIntegral n else maxBound
    characters <-
      case getPermutationSymmetries basis of
        Just (Representation gs) ->
          S.unsafeWith (G.convert (coerce . (.character) <$> gs)) $ \src -> do
            let n = fromIntegral (G.length gs)
            [CU.block| ls_hs_scalar const* {
              size_t const num_bytes = $(size_t n) * sizeof(ls_hs_scalar);
              ls_hs_scalar* dest = aligned_alloc(/*alignment=*/64U, /*size=*/num_bytes);
              LS_CHECK(dest != NULL, "ls_hs_init_basis_info: aligned_alloc failed");
              memcpy(dest, $(ls_hs_scalar const* src), num_bytes);
              return dest;
            } |]
        Nothing -> pure nullPtr
    [CU.block| void {
      ls_hs_basis* const p = $(ls_hs_basis* p);
      if (atomic_load(&p->info) == NULL) {
        ls_hs_basis_info* info = calloc(1, sizeof(ls_hs_basis_info));
        // This should really never happen, but let's check anyway
        if (info == NULL) {
          fprintf(stderr, "ls_hs_init_basis_info: failed to allocate memory, aborting...\n");
          abort();
        }
        *info = (ls_hs_basis_info){
          .has_permutation_symmetries = $(bool has_permutation_symmetries),
          .requires_projection = $(bool requires_projection),
          .is_state_index_identity = $(bool is_state_index_identity),
          .is_real = $(bool is_real),
          .number_bits = $(int number_bits),
          .number_words = $(int number_words),
          .number_sites = $(int number_sites),
          .number_particles = $(int number_particles),
          .number_up = $(int number_up),
          .hamming_weight = $(int hamming_weight),
          .spin_inversion = $(int spin_inversion),
          .min_state_estimate = $(uint64_t min_state_estimate),
          .max_state_estimate = $(uint64_t max_state_estimate),
          .characters = $(ls_hs_scalar const* characters)
        };

        ls_hs_basis_info* _Atomic expected = NULL;
        ls_hs_basis_info* _Atomic desired = info;
        // atomic_compare_exchange_strong atomically does the following:
        // if (p->info == expected) {
        //   p->info = desired;
        //   return true;
        // }
        // else {
        //   expected = p->info;
        //   return false;
        // }
        //
        // In our case the return value false means that another thread beat us to constructing
        // info, so we don't need to set p->info anymore.
        if (!atomic_compare_exchange_strong(&p->info, &expected, desired)) {
          free(info);
        }
      }
    } |]

ls_hs_basis_from_json :: CString -> IO CString
ls_hs_basis_from_json = newCencoded <=< rightM (foldSomeBasis newCbasis) <=< decodeCString @SomeBasis

ls_hs_basis_to_json :: Ptr Cbasis -> IO CString
ls_hs_basis_to_json = foldCbasis newCencoded

--   -- logDebug' $ "basisFromHeader" <> show x
--   -- fp <- mallocForeignPtr
--   kernels <- createCbasis_kernels x
--   payload <- newStablePtr (SomeBasis x)
--   -- addForeignPtrConcFinalizer fp (destroyCbasis_kernels kernels)
--   -- withForeignPtr fp $ \ptr ->
--   --   poke ptr $
--         ls_hs_destroy_state_info_kernel(basis->is_representative_kernel);
--         hs_free_stable_ptr(basis->haskell_payload);

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

-- cloneCbasis :: HasCallStack => Ptr Cbasis -> IO (Ptr Cbasis)
-- cloneCbasis ptr = do
--   _ <- basisIncRefCount ptr
--   pure ptr

withCbasis :: Ptr Cbasis -> (SomeBasis -> IO a) -> IO a
withCbasis p f =
  f
    =<< (deRefStablePtr . castPtrToStablePtr)
    =<< [CU.exp| void* { $(ls_hs_basis* p)->base.haskell_payload } |]

foldCbasis :: (SomeBasis -> IO a) -> Ptr Cbasis -> IO a
foldCbasis f p = withCbasis p f

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

-- data Cbinary_search_kernel_data

-- foreign import ccall safe "lattice_symmetries_haskell.h ls_hs_internal_create_combinadics_kernel_data"
--   ls_hs_internal_create_combinadics_kernel_data :: CInt -> CBool -> IO (Ptr Ccombinadics_kernel_data)

-- foreign import ccall safe "lattice_symmetries_haskell.h ls_hs_internal_destroy_combinadics_kernel_data"
--   ls_hs_internal_destroy_combinadics_kernel_data :: Ptr Ccombinadics_kernel_data -> IO ()

-- foreign import ccall safe "lattice_symmetries_haskell.h &ls_hs_state_index_combinadics_kernel"
--   ls_hs_state_index_combinadics_kernel :: FunPtr Cindex_kernel

-- foreign import ccall safe "lattice_symmetries_haskell.h ls_hs_destroy_state_index_binary_search_kernel_data"
--   ls_hs_destroy_state_index_binary_search_kernel_data :: Ptr Cbinary_search_kernel_data -> IO ()
--
-- foreign import ccall safe "lattice_symmetries_haskell.h &ls_hs_state_index_binary_search_kernel"
--   ls_hs_state_index_binary_search_kernel :: FunPtr Cindex_kernel

-- foreign import ccall safe "lattice_symmetries_haskell.h &ls_hs_state_index_identity_kernel"
--   ls_hs_state_index_identity_kernel :: FunPtr Cindex_kernel

-- data Chalide_kernel_data

-- foreign import ccall safe "lattice_symmetries_haskell.h ls_internal_create_halide_kernel_data"
--   ls_internal_create_halide_kernel_data :: Ptr Cpermutation_group -> CInt -> IO (Ptr Chalide_kernel_data)
--
-- foreign import ccall safe "lattice_symmetries_haskell.h ls_internal_destroy_halide_kernel_data"
--   ls_internal_destroy_halide_kernel_data :: Ptr Chalide_kernel_data -> IO ()
--
-- foreign import ccall safe "lattice_symmetries_haskell.h &ls_hs_is_representative_halide_kernel"
--   ls_hs_is_representative_halide_kernel :: FunPtr Cis_representative_kernel
--
-- foreign import ccall safe "lattice_symmetries_haskell.h &ls_hs_state_info_halide_kernel"
--   ls_hs_state_info_halide_kernel :: FunPtr Cstate_info_kernel

-- setStateInfoKernel :: Basis t -> Cbasis_kernels -> IO Cbasis_kernels
-- setStateInfoKernel (SpinBasis n _ i g) k
--   | n <= 64 && not (nullSymmetries g) = do
--       kernelData <-
--         bracket (newCpermutation_group g) destroyCpermutation_group $ \gPtr ->
--           ls_internal_create_halide_kernel_data gPtr (maybe 0 fromIntegral i)
--       pure
--         $ k
--           { cbasis_state_info_kernel = ls_hs_state_info_halide_kernel
--           , cbasis_state_info_data = castPtr kernelData
--           , cbasis_is_representative_kernel = ls_hs_is_representative_halide_kernel
--           , cbasis_is_representative_data = castPtr kernelData
--           }
-- setStateInfoKernel _ k =
--   pure
--     $ k
--       { cbasis_state_info_kernel = nullFunPtr
--       , cbasis_state_info_data = nullPtr
--       , cbasis_is_representative_kernel = nullFunPtr
--       , cbasis_is_representative_data = nullPtr
--       }

-- setStateIndexKernel :: Basis t -> Cbasis_kernels -> IO Cbasis_kernels
-- setStateIndexKernel _ k =
--   pure
--     $ k
--       { cbasis_state_index_kernel = nullFunPtr
--       , cbasis_state_index_data = nullPtr
--       }

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

-- createCbasis_kernels :: Basis t -> IO (Ptr Cbasis_kernels)
-- createCbasis_kernels x =
--   new
--     <=< setStateInfoKernel x
--     <=< setStateIndexKernel x
--     $ Cbasis_kernels
--       { cbasis_state_info_kernel = nullFunPtr
--       , cbasis_state_info_data = nullPtr
--       , cbasis_is_representative_kernel = nullFunPtr
--       , cbasis_is_representative_data = nullPtr
--       , cbasis_state_index_kernel = nullFunPtr
--       , cbasis_state_index_data = nullPtr
--       }
--
-- destroyCindex_kernel :: HasCallStack => Cbasis_kernels -> IO ()
-- destroyCindex_kernel p
--   -- \| kernel == ls_hs_state_index_combinadics_kernel = ls_hs_internal_destroy_combinadics_kernel_data (castPtr env)
--   | kernel == ls_hs_state_index_binary_search_kernel = ls_hs_destroy_state_index_binary_search_kernel_data (castPtr env)
--   -- \| kernel == ls_hs_state_index_identity_kernel && env == nullPtr = pure ()
--   | kernel == nullFunPtr && env == nullPtr = pure ()
--   | otherwise = error "failed to automatically deallocate state_index kernel"
--   where
--     kernel = cbasis_state_index_kernel p
--     env = cbasis_state_index_data p
--
-- destroyCstate_info_kernel :: HasCallStack => Cbasis_kernels -> IO ()
-- destroyCstate_info_kernel p
--   | cbasis_state_info_kernel p
--       == ls_hs_state_info_halide_kernel
--       && cbasis_is_representative_kernel p
--       == ls_hs_is_representative_halide_kernel =
--       do
--         ls_internal_destroy_halide_kernel_data (castPtr (cbasis_state_info_data p))
--         when (cbasis_state_info_data p /= cbasis_is_representative_data p)
--           $ ls_internal_destroy_halide_kernel_data (castPtr (cbasis_is_representative_data p))
--   | cbasis_state_info_kernel p
--       == nullFunPtr
--       && cbasis_state_info_data p
--       == nullPtr
--       && cbasis_is_representative_kernel p
--       == nullFunPtr
--       && cbasis_is_representative_data p
--       == nullPtr =
--       pure ()
--   | otherwise = error "failed to automatically deallocate state_info and is_representative kernels"
--
-- destroyCbasis_kernels :: HasCallStack => Ptr Cbasis_kernels -> IO ()
-- destroyCbasis_kernels p = do
--   kernels <- peek p
--   destroyCindex_kernel kernels
--   destroyCstate_info_kernel kernels
--   free p

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

-- foreign export ccall "ls_hs_basis_number_bits"
--   ls_hs_basis_number_bits :: Ptr Cbasis -> IO CInt
--
-- ls_hs_basis_number_bits basisPtr =
--   fromIntegral <$> withCbasis basisPtr (foldSomeBasis (pure . getNumberBits))
