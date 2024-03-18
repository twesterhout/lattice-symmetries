{-# LANGUAGE GADTs #-}
{-# LANGUAGE OverloadedLists #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE UndecidableInstances #-}

module LatticeSymmetries.Some
  ( SomeExpr (..)
  , mkSomeExpr
  , withSomeExpr
  , foldSomeExpr
  , mapSomeExpr

    -- * FFI
  , newCexpr
  , withCexpr
  , foldCexpr
  , ls_hs_destroy_expr
  , ls_hs_expr_from_json
  , ls_hs_expr_to_json
  , ls_hs_expr_to_string
  , ls_hs_expr_plus
  , ls_hs_expr_minus
  , ls_hs_expr_times
  , ls_hs_expr_negate
  , ls_hs_expr_scale
  , ls_hs_expr_adjoint
  , ls_hs_expr_equal
  , ls_hs_expr_is_hermitian
  , ls_hs_expr_is_real
  , ls_hs_expr_is_identity
  , ls_hs_replace_indices
  , ls_hs_expr_permutation_group
  , ls_hs_expr_abelian_permutation_group
  , ls_hs_expr_spin_inversion_invariant
  , ls_hs_expr_conserves_number_particles
  , ls_hs_expr_particle_type
  , ls_hs_expr_number_sites
  )
where

import Data.Aeson
import Data.Aeson qualified as Aeson
import Data.Aeson.Types (JSONPathElement (..), Parser, parserThrowError)
import Data.ByteString (packCString)
import Data.IntMap qualified as IntMap
import Data.Map.Strict qualified as Map
import Data.Stream.Monadic qualified
import Data.Vector qualified as B
import Data.Vector.Fusion.Bundle.Monadic qualified as Bundle
import Data.Vector.Fusion.Util (Id (unId))
import Data.Vector.Generic qualified as G
import Foreign (fromBool)
import Foreign.C (CBool, CInt, CString)
import Foreign.Ptr (Ptr, castPtr)
import GHC.Show (Show (showsPrec))
import LatticeSymmetries.Algebra
import LatticeSymmetries.ComplexRational
import LatticeSymmetries.Context
import LatticeSymmetries.Expr
import LatticeSymmetries.FFI
import LatticeSymmetries.Generator
import LatticeSymmetries.Group (abelianSubgroup)
import LatticeSymmetries.Parser hiding (Parser)
import LatticeSymmetries.Utils
import Prettyprinter (Pretty (..))
import Prelude hiding (Product, Sum, identity, toList)

data SomeExpr where
  SomeExpr :: IsBasis t => !(ParticleTag t) -> !(Expr t) -> SomeExpr

instance Pretty SomeExpr where
  pretty (SomeExpr SpinTag expr) = pretty expr
  pretty (SomeExpr SpinlessFermionTag expr) = pretty expr
  pretty (SomeExpr SpinfulFermionTag expr) = pretty expr

instance Show SomeExpr where
  showsPrec prec (SomeExpr SpinTag expr) = showsPrec prec expr
  showsPrec prec (SomeExpr SpinlessFermionTag expr) = showsPrec prec expr
  showsPrec prec (SomeExpr SpinfulFermionTag expr) = showsPrec prec expr

withSomeExpr
  :: SomeExpr
  -> (forall t. IsBasis t => Expr t -> a)
  -> a
withSomeExpr (SomeExpr SpinTag a) f = f a
withSomeExpr (SomeExpr SpinfulFermionTag a) f = f a
withSomeExpr (SomeExpr SpinlessFermionTag a) f = f a

foldSomeExpr
  :: (forall t. IsBasis t => Expr t -> a)
  -> SomeExpr
  -> a
foldSomeExpr f expr = withSomeExpr expr f

mapSomeExpr
  :: (forall t. IsBasis t => Expr t -> Expr t)
  -> SomeExpr
  -> SomeExpr
mapSomeExpr f expr = case expr of
  SomeExpr SpinTag x -> SomeExpr SpinTag (f x)
  SomeExpr SpinfulFermionTag x -> SomeExpr SpinfulFermionTag (f x)
  SomeExpr SpinlessFermionTag x -> SomeExpr SpinlessFermionTag (f x)

instance Eq SomeExpr where
  (==) (SomeExpr SpinTag a) (SomeExpr SpinTag b) = a == b
  (==) (SomeExpr SpinlessFermionTag a) (SomeExpr SpinlessFermionTag b) = a == b
  (==) (SomeExpr SpinfulFermionTag a) (SomeExpr SpinfulFermionTag b) = a == b
  (==) _ _ = False

exprBinaryOp
  :: HasCallStack
  => (forall t. IsBasis t => Expr t -> Expr t -> Expr t)
  -> SomeExpr
  -> SomeExpr
  -> Either Text SomeExpr
exprBinaryOp op = binaryOp (\(a :: Expr t) b -> SomeExpr (particleDispatch @t) (op a b))

binaryOp
  :: HasCallStack
  => (forall t. IsBasis t => Expr t -> Expr t -> a)
  -> SomeExpr
  -> SomeExpr
  -> Either Text a
binaryOp op (SomeExpr SpinTag a) (SomeExpr SpinTag b) = pure $ op a b
binaryOp op (SomeExpr SpinlessFermionTag a) (SomeExpr SpinlessFermionTag b) = pure $ op a b
binaryOp op (SomeExpr SpinfulFermionTag a) (SomeExpr SpinfulFermionTag b) = pure $ op a b
binaryOp _ (SomeExpr t1 _) (SomeExpr t2 _) =
  Left $
    "Expressions are defined for different particle types: "
      <> show (particleTagToType t1)
      <> " and "
      <> show (particleTagToType t2)

-- instance Num SomeExpr where
--   (+) = binaryOp (+)
--   (-) = binaryOp (-)
--   (*) = binaryOp (*)
--   negate = mapSomeExpr negate
--   abs = mapSomeExpr abs
--   signum = mapSomeExpr signum
--   fromInteger _ =
--     error
--       $ "Num instance of SomeExpr does not implement fromInteger; "
--       <> "consider constructing an explicit identity 𝟙₀ and then scaling it"

instance CanScale ComplexRational SomeExpr where
  scale z = mapSomeExpr (scale z)

instance ToJSON SomeExpr where
  toJSON x@(SomeExpr tag _) =
    object
      [ "particle" .= particleTagToType tag
      , "expression" .= withSomeExpr x toPrettyText
      ]

determineParticleTag :: SExpr -> Maybe ParticleTy
determineParticleTag = go
  where
    go (SSum ts) = foldl' (<|>) Nothing (fmap go ts)
    go (SScaled _ p) = go p
    go (SProduct (t :| ts)) = foldl' (<|>) (go t) (fmap go ts)
    go (SPrimitive (SSpinOp {})) = Just SpinTy
    go (SPrimitive (SFermionOp _ Nothing _)) = Just SpinlessFermionTy
    go (SPrimitive (SFermionOp _ (Just _) _)) = Just SpinfulFermionTy
    go (SPrimitive SIdentity) = Nothing

mkSomeExpr :: Maybe ParticleTy -> Text -> Either Text SomeExpr
mkSomeExpr tp s = do
  sexpr <- parseExprEither s
  case maybe (determineParticleTag sexpr) pure tp of
    Just SpinTy -> SomeExpr SpinTag <$> fromSExpr SpinTag sexpr
    Just SpinlessFermionTy -> SomeExpr SpinlessFermionTag <$> fromSExpr SpinlessFermionTag sexpr
    Just SpinfulFermionTy -> SomeExpr SpinfulFermionTag <$> fromSExpr SpinfulFermionTag sexpr
    Nothing -> fail "could not determine the particle type of the expression; try specifying it explicitly"

exprFromJSON :: (Maybe ParticleTy -> Text -> Either Text a) -> (B.Vector (B.Vector Int) -> a -> Either Text a) -> Value -> Parser a
exprFromJSON f expand = withObject "Expr" $ \v -> do
  tp <- v .:? "particle"
  s <- v .: "expression"
  sites <- v .:? "sites"
  let expand' x = maybe (pure x) (`expand` x) sites
  case f tp s >>= expand' of
    Left e -> parserThrowError [Key "expression"] (toString e)
    Right x -> pure x

expandExpr :: forall t m v. (IsBasis t, MonadFail m, G.Vector v Int) => B.Vector (v Int) -> Expr t -> m (Expr t)
expandExpr indices expr =
  fmap simplifyExpr $
    Data.Stream.Monadic.foldl1M' (\a b -> pure $ Expr (unExpr a + unExpr b)) $
      Data.Stream.Monadic.mapM replaceSiteIndicesList $
        (Bundle.elements . Bundle.trans (pure . unId)) $
          G.stream indices
  where
    !indices0 = getSiteIndex @t <$> collectIndices expr

    replaceSiteIndicesList :: v Int -> m (Expr t)
    replaceSiteIndicesList is
      | G.length is == G.length indices0 = pure $ replaceSiteIndices mapping expr
      | otherwise = fail $ "cannot replace indices " <> show indices0 <> " with " <> show (G.toList is) <> "; size mismatch"
      where
        mapping = IntMap.fromAscList (zip (G.toList indices0) (G.toList is))

instance FromJSON SomeExpr where
  parseJSON = exprFromJSON mkSomeExpr expandSomeExpr
    where
      expandSomeExpr :: B.Vector (B.Vector Int) -> SomeExpr -> Either Text SomeExpr
      expandSomeExpr indices (SomeExpr tag expr)
        | G.null indices = pure . SomeExpr tag . Expr . Sum . G.singleton . Scaled 0 $ Product G.empty
        | otherwise = SomeExpr tag <$> expandExpr indices expr

ls_hs_destroy_expr :: MutablePtr Cexpr -> IO ()
ls_hs_destroy_expr = ls_hs_destroy_object (const (pure ())) . castPtr @Cexpr @Cobject

newCexpr :: SomeExpr -> IO (MutablePtr Cexpr)
newCexpr = fmap (castPtr @Cobject @Cexpr) . newCobject

withCexpr :: Ptr Cexpr -> (SomeExpr -> IO a) -> IO a
withCexpr p = withCobject (castPtr @Cexpr @Cobject p)

foldCexpr :: (SomeExpr -> IO a) -> Ptr Cexpr -> IO a
foldCexpr = flip withCexpr

withCexpr2 :: Ptr Cexpr -> Ptr Cexpr -> (SomeExpr -> SomeExpr -> a) -> IO a
withCexpr2 p1 p2 f = withCexpr p1 $ \x1 -> withCexpr p2 $ \x2 -> pure (f x1 x2)

ls_hs_expr_to_json :: Ptr Cexpr -> IO CString
ls_hs_expr_to_json = foldCexpr newCencoded

-- | Text -> Either Text (Ptr Cexpr)
ls_hs_expr_from_json :: CString -> IO CString
ls_hs_expr_from_json = newCencoded <=< rightM newCexpr <=< decodeCString

ls_hs_expr_to_string :: Ptr Cexpr -> IO CString
ls_hs_expr_to_string = foldCexpr $ newCString . encodeUtf8 . toPrettyText

ls_hs_expr_plus :: Ptr Cexpr -> Ptr Cexpr -> IO CString
ls_hs_expr_plus a b = withCexpr2 a b (exprBinaryOp (+)) >>= rightM newCexpr >>= newCencoded

ls_hs_expr_minus :: Ptr Cexpr -> Ptr Cexpr -> IO CString
ls_hs_expr_minus a b = withCexpr2 a b (exprBinaryOp (-)) >>= rightM newCexpr >>= newCencoded

ls_hs_expr_times :: Ptr Cexpr -> Ptr Cexpr -> IO CString
ls_hs_expr_times a b = withCexpr2 a b (exprBinaryOp (*)) >>= rightM newCexpr >>= newCencoded

ls_hs_expr_negate :: Ptr Cexpr -> IO (MutablePtr Cexpr)
ls_hs_expr_negate = foldCexpr (newCexpr . mapSomeExpr negate)

ls_hs_expr_scale :: Double -> Double -> Ptr Cexpr -> IO (MutablePtr Cexpr)
ls_hs_expr_scale re im = foldCexpr (newCexpr . mapSomeExpr (ComplexRational (realToFrac re) (realToFrac im) `scale`))

ls_hs_expr_equal :: Ptr Cexpr -> Ptr Cexpr -> IO CBool
ls_hs_expr_equal a b = fromBool . fromRight False <$> withCexpr2 a b (binaryOp (==))

ls_hs_expr_adjoint :: Ptr Cexpr -> IO (MutablePtr Cexpr)
ls_hs_expr_adjoint = foldCexpr (newCexpr . mapSomeExpr conjugateExpr)

ls_hs_expr_is_hermitian :: Ptr Cexpr -> IO CBool
ls_hs_expr_is_hermitian = foldCexpr $ foldSomeExpr (pure . fromBool . isHermitianExpr)

ls_hs_expr_is_real :: Ptr Cexpr -> IO CBool
ls_hs_expr_is_real = foldCexpr $ foldSomeExpr (pure . fromBool . isRealExpr)

ls_hs_expr_is_identity :: Ptr Cexpr -> IO CBool
ls_hs_expr_is_identity = foldCexpr $ foldSomeExpr (pure . fromBool . isIdentityExpr)

replaceSiteIndicesFn :: forall t. IsBasis t => IntMap Int -> IndexType t -> IndexType t
replaceSiteIndicesFn mapping = case particleDispatch @t of
  SpinTag -> \i -> fromMaybe i (IntMap.lookup i mapping)
  SpinlessFermionTag -> \i -> fromMaybe i (IntMap.lookup i mapping)
  SpinfulFermionTag -> \(σ, i) -> (σ, fromMaybe i (IntMap.lookup i mapping))

replaceSpinIndicesFn :: forall t. IsBasis t => Map SpinIndex SpinIndex -> IndexType t -> IndexType t
replaceSpinIndicesFn mapping = case particleDispatch @t of
  SpinTag -> id
  SpinlessFermionTag -> id
  SpinfulFermionTag -> \(σ, i) -> (fromMaybe σ (Map.lookup σ mapping), i)

getSiteIndex :: forall t. IsBasis t => IndexType t -> Int
getSiteIndex = case particleDispatch @t of
  SpinTag -> id
  SpinlessFermionTag -> id
  SpinfulFermionTag -> snd

replaceSiteIndices :: forall t. IsBasis t => IntMap Int -> Expr t -> Expr t
replaceSiteIndices mapping expr = simplifyExpr $ mapIndices (replaceSiteIndicesFn @t mapping) expr

replaceSpinIndices :: Map SpinIndex SpinIndex -> SomeExpr -> SomeExpr
replaceSpinIndices mapping = mapSomeExpr $ \(expr :: Expr t) ->
  simplifyExpr $ mapIndices (replaceSpinIndicesFn @t mapping) expr

ls_hs_replace_indices :: Ptr Cexpr -> CString -> IO CString
ls_hs_replace_indices expr jsonString = do
  s <- toLazy <$> packCString jsonString
  newCencoded
    =<< rightM newCexpr
    =<< case Aeson.decode @[(Int, Int)] s of
      Just m -> withCexpr expr (pure . Right . mapSomeExpr (replaceSiteIndices (IntMap.fromList m)))
      Nothing -> case Aeson.decode @[(SpinIndex, SpinIndex)] s of
        Just m -> withCexpr expr (pure . Right . replaceSpinIndices (Map.fromList m))
        Nothing -> pure . Left @Text $ "ls_hs_replace_indices: invalid argument: expected a list of index replacements"

ls_hs_expr_permutation_group :: Ptr Cexpr -> IO CString
ls_hs_expr_permutation_group = foldCexpr $ foldSomeExpr $ newCencoded . exprPermutationGroup Nothing

ls_hs_expr_abelian_permutation_group :: Ptr Cexpr -> IO CString
ls_hs_expr_abelian_permutation_group =
  foldCexpr $ foldSomeExpr $ newCencoded . abelianSubgroup . exprPermutationGroup Nothing

ls_hs_expr_spin_inversion_invariant :: Ptr Cexpr -> IO CBool
ls_hs_expr_spin_inversion_invariant = foldCexpr $ foldSomeExpr $ pure . fromBool . isInvariantUponSpinInversion

ls_hs_expr_conserves_number_particles :: Ptr Cexpr -> IO CBool
ls_hs_expr_conserves_number_particles = foldCexpr $ foldSomeExpr $ pure . fromBool . conservesNumberParticles

ls_hs_expr_particle_type :: Ptr Cexpr -> IO CString
ls_hs_expr_particle_type = foldCexpr $ \(SomeExpr t _) -> newCString . encodeUtf8 . toPrettyText $ particleTagToType t

ls_hs_expr_number_sites :: Ptr Cexpr -> IO CInt
ls_hs_expr_number_sites = foldCexpr $ foldSomeExpr $ pure . fromIntegral . estimateNumberSites
