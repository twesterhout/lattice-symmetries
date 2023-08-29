{-# LANGUAGE GADTs #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE OverloadedLists #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE NoMonomorphismRestriction #-}

module LatticeSymmetries.Expr
  ( Expr (..)
  , mkExpr
  , mkExprEither
  , SomeExpr (..)
  , mkSomeExpr
  , mkSomeExprEither
  , withSomeExpr
  , foldSomeExpr
  , mapSomeExpr
  , conjugateExpr
  , isHermitianExpr
  , isRealExpr
  , isIdentityExpr
  , mapGenerators
  , mapGeneratorsM
  , mapIndices
  , mapIndicesM
  -- mapCoeffs,
  , simplifyExpr
  , replicateSiteIndices

    -- ** FFI helpers
  , Cexpr (..)
  , newCexpr
  , destroyCexpr
  , withCexpr
  , withCexpr2
  )
where

import Data.Aeson
import Data.Aeson.Types (JSONPathElement (..), Parser, parserThrowError)
import Data.List qualified as List
import Data.Map.Strict qualified as Map
import Data.Set qualified as Set
import Data.Vector.Generic qualified as G
import Foreign.Marshal (free, new)
import Foreign.Ptr (Ptr)
import Foreign.StablePtr
import Foreign.Storable
import LatticeSymmetries.Algebra
import LatticeSymmetries.Basis
import LatticeSymmetries.ComplexRational
import LatticeSymmetries.Generator
import LatticeSymmetries.Parser
import LatticeSymmetries.Utils
import Prettyprinter (Pretty (..))
import Prettyprinter qualified as Pretty
import Prelude hiding (Product, Sum, identity, toList)

newtype Expr t = Expr
  { unExpr :: Polynomial ComplexRational (Generator (IndexType t) (GeneratorType t))
  }
  deriving stock (Generic)

deriving instance (Eq (IndexType t), Eq (GeneratorType t)) => Eq (Expr t)

deriving instance (Show (IndexType t), Show (GeneratorType t)) => Show (Expr t)

mapGenerators :: (g ~ Generator (IndexType t) (GeneratorType t)) => (g -> g) -> Expr t -> Expr t
mapGenerators f = Expr . fmap (fmap (fmap f)) . unExpr

mapGeneratorsM :: (g ~ Generator (IndexType t) (GeneratorType t), Monad m) => (g -> m g) -> Expr t -> m (Expr t)
mapGeneratorsM f =
  fmap (Expr . Sum)
    . G.mapM (\(Scaled c p) -> Scaled c . Product <$> G.mapM f (unProduct p))
    . unSum
    . unExpr

mapIndices :: (IndexType t -> IndexType t) -> Expr t -> Expr t
mapIndices f = mapGenerators (\(Generator i g) -> Generator (f i) g)

mapIndicesM :: Monad m => (IndexType t -> m (IndexType t)) -> Expr t -> m (Expr t)
mapIndicesM f = mapGeneratorsM (\(Generator i g) -> Generator <$> f i <*> pure g)

collectIndices :: Ord (IndexType t) => Expr t -> [IndexType t]
collectIndices = List.nub . List.sort . collectSum . unExpr
  where
    collectSum (Sum v) = concatMap collectScaled (G.toList v)
    collectScaled (Scaled _ p) = collectProduct p
    collectProduct (Product v) = (\(Generator i _) -> i) <$> G.toList v

foldlCoeffs' :: (a -> ComplexRational -> a) -> a -> Expr t -> a
foldlCoeffs' combine x₀ (Expr (Sum s)) =
  G.foldl' (\ !x (Scaled c _) -> combine x c) x₀ s

-- mapCoeffs :: (ComplexRational -> ComplexRational) -> Expr t -> Expr t
-- mapCoeffs f = Expr . fmap (\(Scaled c p) -> Scaled (f c) p) . unExpr

simplifyExpr :: (Algebra (GeneratorType t), Ord (IndexType t)) => Expr t -> Expr t
simplifyExpr = Expr . simplifyPolynomial . unExpr

conjugateExpr :: (Algebra (GeneratorType t), Ord (IndexType t)) => Expr t -> Expr t
conjugateExpr = simplifyExpr . Expr . conjugateSum . unExpr
  where
    conjugateSum = fmap conjugateScaled
    conjugateScaled (Scaled c p) = Scaled (conjugate c) (conjugateProduct p)
    conjugateProduct (Product gs) = Product . G.reverse $ conjugateGenerator <$> gs

isIdentityExpr :: Algebra (GeneratorType t) => Expr t -> Bool
isIdentityExpr = isIdentitySum . unExpr
  where
    isIdentitySum (Sum terms) = G.length terms == 1 && isIdentityScaled (G.head terms)
    isIdentityScaled (Scaled _ (Product p)) = G.length p == 1 && isIdentityGenerator (G.head p)
    isIdentityGenerator (Generator _ g) = isIdentity g

isHermitianExpr :: (Ord (IndexType t), Eq (GeneratorType t), Algebra (GeneratorType t)) => Expr t -> Bool
isHermitianExpr terms = terms == conjugateExpr terms

isRealExpr :: Expr t -> Bool
isRealExpr = foldlCoeffs' (\f c -> f && imagPart c == 0) True

instance (Algebra (GeneratorType t), Ord (IndexType t)) => Num (Expr t) where
  (+) a b = simplifyExpr . Expr $ unExpr a + unExpr b
  (-) a b = simplifyExpr . Expr $ unExpr a - unExpr b
  (*) a b = simplifyExpr . Expr $ unExpr a * unExpr b
  negate = Expr . negate . unExpr
  abs = Expr . abs . unExpr
  signum = Expr . signum . unExpr
  fromInteger _ =
    error $
      "Num instance of Expr does not implement fromInteger; "
        <> "consider constructing an explicit identity 𝟙₀ and then scaling it"

instance CanScale ComplexRational (Expr t) where
  scale c a
    | c == 0 = Expr []
    | otherwise = Expr $ c `scale` unExpr a

instance Pretty (Generator (IndexType t) (GeneratorType t)) => Pretty (Expr t) where
  pretty (Expr (Sum terms)) = prettySum (G.toList terms)
    where
      prettyScaled (Scaled c g)
        | c == 1 = prettyProduct g
        | c == -1 = "-" <> prettyProduct g
        | otherwise = pretty c <> " " <> prettyProduct g
      prettyProduct (Product v) =
        Pretty.encloseSep mempty mempty " " (G.toList $ fmap pretty v)
      prettySum [] = "0"
      prettySum (term : others) = prettyScaled term <> prettySumRest others
      prettySumRest [] = ""
      prettySumRest (Scaled c p : others)
        | (realPart c == 0 && imagPart c < 0) || (imagPart c == 0 && realPart c < 0) =
            " - " <> prettyScaled (Scaled (-c) p) <> prettySumRest others
        | otherwise = " + " <> prettyScaled (Scaled c p) <> prettySumRest others

-- instance Pretty g => Pretty (Product g) where
--   pretty (Product v) =
--     Pretty.encloseSep mempty mempty " " (G.toList $ fmap pretty v)
--
-- instance Pretty g => Pretty (Sum g) where
--   pretty (Sum v)
--     | not (G.null v) = Pretty.encloseSep mempty mempty " + " (G.toList $ fmap pretty v)
--     | otherwise = pretty (0 :: Int)

instance Pretty SomeExpr where
  pretty (SomeExpr SpinTag expr) = pretty expr
  pretty (SomeExpr SpinlessFermionTag expr) = pretty expr
  pretty (SomeExpr SpinfulFermionTag expr) = pretty expr

data SomeExpr where
  SomeExpr :: IsBasis t => !(ParticleTag t) -> !(Expr t) -> SomeExpr

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

binaryOp
  :: HasCallStack
  => (forall t. (Algebra (GeneratorType t), Ord (IndexType t)) => Expr t -> Expr t -> Expr t)
  -> SomeExpr
  -> SomeExpr
  -> SomeExpr
binaryOp op (SomeExpr tag@SpinTag a) (SomeExpr SpinTag b) = SomeExpr tag $ op a b
binaryOp op (SomeExpr tag@SpinlessFermionTag a) (SomeExpr SpinlessFermionTag b) = SomeExpr tag $ op a b
binaryOp op (SomeExpr tag@SpinfulFermionTag a) (SomeExpr SpinfulFermionTag b) = SomeExpr tag $ op a b
binaryOp _ (SomeExpr t1 _) (SomeExpr t2 _) =
  error $
    "Expressions are defined for different particle types: "
      <> show (particleTagToType t1)
      <> " and "
      <> show (particleTagToType t2)

instance Num SomeExpr where
  (+) = binaryOp (+)
  (-) = binaryOp (-)
  (*) = binaryOp (*)
  negate = mapSomeExpr negate
  abs = mapSomeExpr abs
  signum = mapSomeExpr signum
  fromInteger _ =
    error $
      "Num instance of SomeExpr does not implement fromInteger; "
        <> "consider constructing an explicit identity 𝟙₀ and then scaling it"

instance CanScale ComplexRational SomeExpr where
  scale z = mapSomeExpr (scale z)

instance ToJSON SomeExpr where
  toJSON x@(SomeExpr tag _) =
    object
      [ "particle" .= particleTagToType tag
      , "expression" .= withSomeExpr x toPrettyText
      ]

-- tableFromLowLevelMapping
--   :: (i ~ IndexType t, Ord i)
--   => ParticleTag t
--   -> Int
--   -> Ptr CInt
--   -> Ptr CInt
--   -> IO (Map i i)
-- tableFromLowLevelMapping tag count fromPtr toPtr =
--   Map.fromList
--     <$> case tag of
--       SpinTag ->
--         zip
--           <$> (fmap fromIntegral <$> peekArray count fromPtr)
--           <*> (fmap fromIntegral <$> peekArray count toPtr)
--       SpinlessFermionTag ->
--         zip
--           <$> (fmap fromIntegral <$> peekArray count fromPtr)
--           <*> (fmap fromIntegral <$> peekArray count toPtr)
--       SpinfulFermionTag ->
--         zip
--           <$> (toSpinfulIndex <$> peekArray (2 * count) fromPtr)
--           <*> (toSpinfulIndex <$> peekArray (2 * count) toPtr)
--   where
--     toSpinfulIndex :: [CInt] -> [(SpinIndex, Int)]
--     toSpinfulIndex [] = []
--     toSpinfulIndex (s : i : rest) = (toEnum (fromIntegral s), fromIntegral i) : toSpinfulIndex rest
--     toSpinfulIndex _ = error "this cannot happen by construction"

collectSiteIndices :: forall t. HasProperIndexType t => Expr t -> [Int]
collectSiteIndices = Set.toList . Set.fromList . fmap getSiteIndex . collectIndices

replicateSiteIndices
  :: forall t
   . (HasProperIndexType t, Algebra (GeneratorType t))
  => [[Int]]
  -> Expr t
  -> Expr t
replicateSiteIndices newIndices expr =
  case newIndices of
    [] -> Expr []
    (i : is) -> foldl' (+) (replace i) $ fmap replace is
  where
    oldSiteIndices = collectSiteIndices expr
    k = length oldSiteIndices
    replace siteIndices
      | Map.size mapping == k = mapIndices (mapSiteIndex (mapping Map.!)) expr
      | otherwise =
          error $
            "wrong number of site indices: " <> show (length siteIndices) <> "; expected " <> show k
      where
        mapping = Map.fromList (zip oldSiteIndices siteIndices)

newtype Cexpr = Cexpr
  { unCexpr :: StablePtr SomeExpr
  }
  deriving stock (Eq)
  deriving newtype (Storable)

newCexpr :: SomeExpr -> IO (Ptr Cexpr)
newCexpr expr = (new . Cexpr) =<< newStablePtr expr

destroyCexpr :: Ptr Cexpr -> IO ()
destroyCexpr p = do
  freeStablePtr . unCexpr =<< peek p
  free p

withCexpr :: Ptr Cexpr -> (SomeExpr -> IO a) -> IO a
withCexpr p f = f =<< deRefStablePtr . unCexpr =<< peek p

withCexpr2 :: Ptr Cexpr -> Ptr Cexpr -> (SomeExpr -> SomeExpr -> a) -> IO a
withCexpr2 p1 p2 f =
  withCexpr p1 $ \x1 ->
    withCexpr p2 $ \x2 ->
      pure (f x1 x2)

fromSExpr :: (HasCallStack, IsBasis t) => ParticleTag t -> SExpr -> Either Text (Expr t)
fromSExpr t (SSum terms) = do
  exprs <- mapM (fromSExpr t) terms
  pure $ foldl' (+) (Expr []) exprs
fromSExpr t (SScaled c term) = scale c <$> fromSExpr t term
fromSExpr t (SProduct (term :| terms)) = do
  expr <- fromSExpr t term
  exprs <- mapM (fromSExpr t) terms
  pure $ foldl' (*) expr exprs
fromSExpr SpinTag (SPrimitive (SSpinOp c t i)) = pure $
  case t of
    SSpinPlus -> Expr [scaled [Generator i SpinPlus]]
    SSpinMinus -> Expr [scaled [Generator i SpinMinus]]
    SSpinZ -> Expr [scaled [Generator i SpinZ]]
    SSpinX -> Expr [scaled [Generator i SpinPlus], scaled [Generator i SpinMinus]]
    SSpinY ->
      scale (ComplexRational 0 (-1)) $
        Expr [scaled [Generator i SpinPlus], scale (-1 :: ℂ) (scaled [Generator i SpinMinus])]
  where
    scaled = if c == 'S' then Scaled 0.5 else Scaled 1
fromSExpr SpinfulFermionTag (SPrimitive (SFermionOp t (Just s) i)) = pure $
  case t of
    SFermionCreate -> Expr [Scaled 1 [Generator (s, i) FermionCreate]]
    SFermionAnnihilate -> Expr [Scaled 1 [Generator (s, i) FermionAnnihilate]]
    SFermionNumber -> Expr [Scaled 1 [Generator (s, i) FermionCount]]
fromSExpr SpinlessFermionTag (SPrimitive (SFermionOp t Nothing i)) = pure $
  case t of
    SFermionCreate -> Expr [Scaled 1 [Generator i FermionCreate]]
    SFermionAnnihilate -> Expr [Scaled 1 [Generator i FermionAnnihilate]]
    SFermionNumber -> Expr [Scaled 1 [Generator i FermionCount]]
fromSExpr SpinTag (SPrimitive _) = fail "expected an expression for spin-1/2 particles"
fromSExpr SpinfulFermionTag (SPrimitive _) = fail "expected an expression for spinful fermions"
fromSExpr SpinlessFermionTag (SPrimitive _) = fail "expected an expression for spinless fermions"

mkExprEither :: (HasCallStack, IsBasis t) => ParticleTag t -> Text -> Either Text (Expr t)
mkExprEither tag s = parseExprEither s >>= fromSExpr tag

mkExpr :: (HasCallStack, IsBasis t) => ParticleTag t -> Text -> Expr t
mkExpr tag s = either error id (mkExprEither tag s)

determineParticleTag :: SExpr -> ParticleTy
determineParticleTag expr = fromMaybe SpinTy (go expr)
  where
    go (SSum ts) = foldl' (<|>) Nothing (fmap go ts)
    go (SScaled _ p) = go p
    go (SProduct (t :| ts)) = foldl' (<|>) (go t) (fmap go ts)
    go (SPrimitive (SSpinOp {})) = Just SpinTy
    go (SPrimitive (SFermionOp _ Nothing _)) = Just SpinlessFermionTy
    go (SPrimitive (SFermionOp _ (Just _) _)) = Just SpinfulFermionTy

mkSomeExprEither :: Maybe ParticleTy -> Text -> Either Text SomeExpr
mkSomeExprEither tp s = do
  sexpr <- parseExprEither s
  case fromMaybe (determineParticleTag sexpr) tp of
    SpinTy -> SomeExpr SpinTag <$> fromSExpr SpinTag sexpr
    SpinlessFermionTy -> SomeExpr SpinlessFermionTag <$> fromSExpr SpinlessFermionTag sexpr
    SpinfulFermionTy -> SomeExpr SpinfulFermionTag <$> fromSExpr SpinfulFermionTag sexpr

mkSomeExpr :: HasCallStack => Maybe ParticleTy -> Text -> SomeExpr
mkSomeExpr tp s = either error id $ mkSomeExprEither tp s

exprFromJSON :: (Maybe ParticleTy -> Text -> Either Text a) -> ([[Int]] -> a -> a) -> Value -> Parser a
exprFromJSON f expand = withObject "Expr" $ \v -> do
  tp <- v .:? "particle"
  s <- v .: "expression"
  sites <- v .:? "sites"
  let expand' x = case sites of
        Just is -> expand is x
        Nothing -> x
  case fmap expand' (f tp s) of
    Left e -> parserThrowError [Key "expression"] (toString e)
    Right x -> pure x

instance IsBasis t => FromJSON (Expr t) where
  parseJSON = exprFromJSON (exprParserConcrete (particleDispatch @t)) replicateSiteIndices
    where
      exprParserConcrete :: ParticleTag t -> Maybe ParticleTy -> Text -> Either Text (Expr t)
      exprParserConcrete tag Nothing s = mkExprEither tag s
      exprParserConcrete tag (Just tp) s
        | particleTagToType tag == tp = mkExprEither tag s
        | otherwise = Left $ "invalid particle type: " <> show tp <> "; expected " <> show (particleTagToType tag)

instance FromJSON SomeExpr where
  parseJSON = exprFromJSON mkSomeExprEither expandSomeExpr
    where
      expandSomeExpr :: [[Int]] -> SomeExpr -> SomeExpr
      expandSomeExpr indices = mapSomeExpr (replicateSiteIndices indices)
