{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE DefaultSignatures #-}
{-# LANGUAGE FunctionalDependencies #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE InstanceSigs #-}
{-# LANGUAGE OverloadedLists #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE NoMonomorphismRestriction #-}

module LatticeSymmetries.Expr
  ( Expr (..),
    SomeExpr (..),
    withSomeExpr,
    foldSomeExpr,
    mapSomeExpr,
    conjugateExpr,
    isHermitianExpr,
    isRealExpr,
    isIdentityExpr,
    mapGenerators,
    mapGeneratorsM,
    mapIndices,
    mapIndicesM,
    -- mapCoeffs,
    simplifyExpr,
    replicateSiteIndices,

    -- ** FFI helpers
    Cexpr (..),
    newCexpr,
    destroyCexpr,
    withCexpr,
    withCexpr2,
  )
where

import Control.Exception (assert)
import Control.Monad.ST
import Data.Aeson
import Data.Aeson.Types (parserThrowError)
import qualified Data.List as List
import qualified Data.Map.Strict as Map
import qualified Data.Set as Set
import qualified Data.Text as Text
import Data.Vector (Vector)
import qualified Data.Vector.Algorithms.Intro
import qualified Data.Vector.Fusion.Bundle as Bundle (inplace)
import Data.Vector.Fusion.Bundle.Size (toMax)
import Data.Vector.Fusion.Stream.Monadic (Step (..), Stream (..))
import Data.Vector.Generic ((!))
import qualified Data.Vector.Generic as G
import Foreign.C.Types (CInt (..))
import Foreign.Marshal (free, new, peekArray)
import Foreign.Ptr (Ptr)
import Foreign.StablePtr
import Foreign.Storable
import GHC.Exts (IsList (..))
import LatticeSymmetries.Algebra
import LatticeSymmetries.Basis
import LatticeSymmetries.ComplexRational
import LatticeSymmetries.Dense
import LatticeSymmetries.Generator
import LatticeSymmetries.NonbranchingTerm
import LatticeSymmetries.Utils
import Prettyprinter (Doc, Pretty (..))
import qualified Prettyprinter as Pretty
import Prettyprinter.Render.Text (renderStrict)
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
  fmap Expr
    . fmap Sum
    . G.mapM (\(Scaled c p) -> Scaled c . Product <$> G.mapM f (unProduct p))
    . unSum
    . unExpr

mapIndices :: (IndexType t -> IndexType t) -> Expr t -> Expr t
mapIndices f = mapGenerators (\(Generator i g) -> Generator (f i) g)

mapIndicesM :: Monad m => (IndexType t -> m (IndexType t)) -> Expr t -> m (Expr t)
mapIndicesM f = mapGeneratorsM (\(Generator i g) -> Generator <$> (f i) <*> pure g)

collectIndices :: Ord (IndexType t) => Expr t -> [IndexType t]
collectIndices = List.nub . List.sort . collectSum . unExpr
  where
    collectSum (Sum v) = concatMap collectScaled (G.toList v)
    collectScaled (Scaled _ p) = collectProduct p
    collectProduct (Product v) = (\(Generator i _) -> i) <$> (G.toList v)

-- mapGenerators' :: (g -> g) -> Polynomial c g -> Polynomial c g
-- mapGenerators' f = fmap (fmap (fmap f))

-- foldlGenerators' :: (a -> g -> a) -> a -> Polynomial c g -> a
-- foldlGenerators' combine x₀ (Sum s) =
--   G.foldl' (\ !x (Scaled _ (Product p)) -> G.foldl' combine x p) x₀ s

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
    | otherwise = Expr $ c `scale` (unExpr a)

instance Pretty (Generator (IndexType t) (GeneratorType t)) => Pretty (Expr t) where
  pretty (Expr terms) = pretty terms

instance Pretty SomeExpr where
  pretty (SomeExpr SpinTag expr) = pretty expr
  pretty (SomeExpr SpinlessFermionTag expr) = pretty expr
  pretty (SomeExpr SpinfulFermionTag expr) = pretty expr

data SomeExpr where
  SomeExpr :: IsBasis t => !(ParticleTag t) -> !(Expr t) -> SomeExpr

withSomeExpr ::
  SomeExpr ->
  (forall t. IsBasis t => Expr t -> a) ->
  a
withSomeExpr (SomeExpr SpinTag a) f = f a
withSomeExpr (SomeExpr SpinfulFermionTag a) f = f a
withSomeExpr (SomeExpr SpinlessFermionTag a) f = f a

foldSomeExpr ::
  (forall t. IsBasis t => Expr t -> a) ->
  SomeExpr ->
  a
foldSomeExpr f expr = withSomeExpr expr f

mapSomeExpr ::
  (forall t. IsBasis t => Expr t -> Expr t) ->
  SomeExpr ->
  SomeExpr
mapSomeExpr f expr = case expr of
  SomeExpr SpinTag x -> SomeExpr SpinTag (f x)
  SomeExpr SpinfulFermionTag x -> SomeExpr SpinfulFermionTag (f x)
  SomeExpr SpinlessFermionTag x -> SomeExpr SpinlessFermionTag (f x)

binaryOp ::
  (forall t. (Algebra (GeneratorType t), Ord (IndexType t)) => Expr t -> Expr t -> Expr t) ->
  SomeExpr ->
  SomeExpr ->
  SomeExpr
binaryOp op (SomeExpr tag@SpinTag a) (SomeExpr SpinTag b) = SomeExpr tag $ op a b
binaryOp op (SomeExpr tag@SpinlessFermionTag a) (SomeExpr SpinlessFermionTag b) = SomeExpr tag $ op a b
binaryOp op (SomeExpr tag@SpinfulFermionTag a) (SomeExpr SpinfulFermionTag b) = SomeExpr tag $ op a b
binaryOp op (SomeExpr t₁ _) (SomeExpr t₂ _) =
  error $
    "Expressions are defined for different particle types: "
      <> show (particleTagToType t₁)
      <> " and "
      <> show (particleTagToType t₂)

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
    object $
      [ "particle" .= particleTagToType tag,
        "expression" .= withSomeExpr x toPrettyText
      ]

tableFromLowLevelMapping ::
  (i ~ IndexType t, Ord i) =>
  ParticleTag t ->
  Int ->
  Ptr CInt ->
  Ptr CInt ->
  IO (Map i i)
tableFromLowLevelMapping tag count fromPtr toPtr =
  Map.fromList
    <$> case tag of
      SpinTag ->
        zip
          <$> (fmap fromIntegral <$> peekArray count fromPtr)
          <*> (fmap fromIntegral <$> peekArray count toPtr)
      SpinlessFermionTag ->
        zip
          <$> (fmap fromIntegral <$> peekArray count fromPtr)
          <*> (fmap fromIntegral <$> peekArray count toPtr)
      SpinfulFermionTag ->
        zip
          <$> (toSpinfulIndex <$> peekArray (2 * count) fromPtr)
          <*> (toSpinfulIndex <$> peekArray (2 * count) toPtr)
  where
    toSpinfulIndex :: [CInt] -> [(SpinIndex, Int)]
    toSpinfulIndex [] = []
    toSpinfulIndex (s : i : rest) = (toEnum (fromIntegral s), fromIntegral i) : toSpinfulIndex rest
    toSpinfulIndex _ = error "this cannot happen by construction"

collectSiteIndices :: forall t. HasProperIndexType t => Expr t -> [Int]
collectSiteIndices = Set.toList . Set.fromList . fmap getSiteIndex . collectIndices

replicateSiteIndices ::
  forall t.
  (HasProperIndexType t, Algebra (GeneratorType t)) =>
  [[Int]] ->
  Expr t ->
  Expr t
replicateSiteIndices newIndices expr@(Expr poly) =
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

newtype {-# CTYPE "lattice_symmetries_haskell.h" "ls_hs_expr" #-} Cexpr = Cexpr
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
