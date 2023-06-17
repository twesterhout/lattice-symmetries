{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE OverloadedLists #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE NoMonomorphismRestriction #-}

module LatticeSymmetries.Expr
  ( Expr (..)
  , SomeExpr (..)
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
import LatticeSymmetries.Utils
import Prettyprinter (Pretty (..))
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

mapCoeffs :: (ComplexRational -> ComplexRational) -> Expr t -> Expr t
mapCoeffs f = Expr . fmap (\(Scaled c p) -> Scaled (f c) p) . unExpr

simplifyExpr :: (Algebra (GeneratorType t), Ord (IndexType t)) => Expr t -> Expr t
simplifyExpr = Expr . simplifyPolynomial . unExpr

conjugateExpr :: (Algebra (GeneratorType t), Ord (IndexType t)) => Expr t -> Expr t
conjugateExpr = simplifyExpr . mapGenerators conjugateGenerator . mapCoeffs conjugate

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
  pretty (Expr terms) = pretty terms

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
