{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE DefaultSignatures #-}
{-# LANGUAGE FunctionalDependencies #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE InstanceSigs #-}
{-# LANGUAGE OverloadedLists #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE NoMonomorphismRestriction #-}

module LatticeSymmetries.Algebra
  ( -- * Polynomials
    Scaled (..),
    Product (..),
    Sum (..),
    Polynomial,
    Expr (..),
    SomeExpr (..),
    withSomeExpr,

    -- * Fermionic and spin (bosonic) algebra
    CommutatorType (..),
    Algebra (..),

    -- * Simplification
    expandProduct,
    productToCanonical,
    simplifyProductNoIndices,
    simplifyProduct,
    simplifyPolynomial,
    conjugatePolynomial,
    CanScale (..),
    forSiteIndices,
    HasSiteIndex,

    -- * Expressions
    conjugateExpr,
    isHermitianExpr,
    isIdentityExpr,
    mapGenerators,
    mapGeneratorsM,
    mapIndices,
    mapIndicesM,
    mapCoeffs,
    simplifyExpr,
  )
where

import Control.Exception (assert)
import Control.Monad.ST
import Data.Aeson
import Data.Aeson.Types (parserThrowError)
import qualified Data.List as List
import qualified Data.Map.Strict as Map
import qualified Data.Text as Text
import Data.Vector (Vector)
import qualified Data.Vector.Algorithms.Intro
import qualified Data.Vector.Fusion.Bundle as Bundle (inplace)
import Data.Vector.Fusion.Bundle.Size (toMax)
import Data.Vector.Fusion.Stream.Monadic (Step (..), Stream (..))
import Data.Vector.Generic ((!))
import qualified Data.Vector.Generic as G
import Foreign.C.Types (CInt (..))
import Foreign.Marshal.Array (peekArray)
import Foreign.Ptr (Ptr)
import GHC.Exts (IsList (..))
import LatticeSymmetries.ComplexRational
import LatticeSymmetries.Dense
import LatticeSymmetries.Generator
import LatticeSymmetries.NonbranchingTerm
import LatticeSymmetries.Utils
import Prettyprinter (Doc, Pretty (..))
import qualified Prettyprinter as Pretty
import Prettyprinter.Render.Text (renderStrict)
import Prelude hiding (Product, Sum, identity, toList)

-- | Represents a term of the form @c × g@ where @c@ is typically a scalar and @g@ is some
-- expression.
data Scaled c g = Scaled !c !g
  deriving stock (Eq, Ord, Show, Generic)

-- | A product of @g@s.
newtype Product g = Product {unProduct :: Vector g}
  deriving stock (Eq, Ord, Show, Generic)
  deriving newtype (IsList, Functor, Semigroup, Foldable)

-- | A sum of @g@s.
newtype Sum g = Sum {unSum :: Vector g}
  deriving stock (Eq, Show, Generic)
  deriving newtype (IsList, Functor, Semigroup, Monoid, Foldable)

-- | A polynomial in variable @g@ over a field @c@.
type Polynomial c g = Sum (Scaled c (Product g))

newtype Expr t = Expr
  { unExpr :: Polynomial ComplexRational (Generator (IndexType t) (GeneratorType t))
  }
  deriving stock (Generic)

deriving instance (Eq (IndexType t), Eq (GeneratorType t)) => Eq (Expr t)

deriving instance (Show (IndexType t), Show (GeneratorType t)) => Show (Expr t)

-- | Specifies the type of commutator to use for an algebra.
data CommutatorType
  = -- | Commutator \([a, b] = ab - ba\)
    Commutator
  | -- | Anticommutator \(\{a, b\} = ab + ba\)
    Anticommutator
  deriving stock (Show, Eq)

class Ord g => Algebra g where
  -- | The type of commutator this algebra uses.
  --
  -- We have 'Commutator' for spin (or bosonic) systems and 'Anticommutator' for fermionic systems.
  nonDiagonalCommutatorType :: CommutatorType

  -- | Check whether a given generator is an identity.
  isIdentity :: g -> Bool

  -- | Check whether a given generator is diagonal.
  isDiagonal :: g -> Bool

  -- | Compute the Hermitian conjugate of a generator
  conjugateGenerator :: g -> g

  -- | Compute the commutator between two generators.
  --
  -- Depending on the algebra, this operation will compute either the commutator or the
  -- anticommutator. In other words, this function defines the following:
  -- \(g_a g_b = \pm g_b g_a + \sum_i c_i g_i\).
  commute :: Fractional c => g -> g -> (CommutatorType, Sum (Scaled c g))

  -- | Simplify a product of multiple generators.
  simplifyOrderedProduct :: Fractional c => g -> g -> Sum (Scaled c g)

class Num c => CanScale c a where
  scale :: c -> a -> a

data CombineNeighborsHelper a
  = CombineNeighborsFirst
  | CombineNeighborsPrevious !a
  | CombineNeighborsDone

combineNeighborsImpl :: Monad m => (a -> a -> Bool) -> (a -> a -> a) -> Stream m a -> Stream m a
{-# INLINE combineNeighborsImpl #-}
combineNeighborsImpl equal combine (Stream step s₀) = Stream step' (CombineNeighborsFirst, s₀)
  where
    {-# INLINE step' #-}
    step' (CombineNeighborsFirst, s) = do
      r <- step s
      case r of
        Yield a s' -> pure $ Skip (CombineNeighborsPrevious a, s')
        Skip s' -> pure $ Skip (CombineNeighborsFirst, s')
        Done -> pure $ Done
    step' (CombineNeighborsPrevious a, s) = do
      r <- step s
      case r of
        Yield b s' ->
          if equal a b
            then pure $ Skip (CombineNeighborsPrevious (combine a b), s')
            else pure $ Yield a (CombineNeighborsPrevious b, s')
        Skip s' -> pure $ Skip (CombineNeighborsPrevious a, s')
        Done -> pure $ Yield a (CombineNeighborsDone, s)
    step' (CombineNeighborsDone, _) = pure $ Done

combineNeighbors :: G.Vector v a => (a -> a -> Bool) -> (a -> a -> a) -> v a -> v a
combineNeighbors equal combine =
  G.unstream . Bundle.inplace (combineNeighborsImpl equal combine) toMax . G.stream

instance Functor (Scaled c) where
  fmap f (Scaled c g) = Scaled c (f g)

instance (Pretty c, Pretty g) => Pretty (Scaled c g) where
  pretty (Scaled c g) = pretty c <> " × " <> pretty g

instance
  HasNonbranchingRepresentation g =>
  HasNonbranchingRepresentation (Scaled ComplexRational g)
  where
  nonbranchingRepresentation (Scaled c g) = t {nbtV = c * nbtV t}
    where
      t = nonbranchingRepresentation g

instance
  HasNonbranchingRepresentation g =>
  HasNonbranchingRepresentation (Product g)
  where
  nonbranchingRepresentation (Product v)
    | not (G.null v) = G.foldl1' (<>) . G.map nonbranchingRepresentation $ v
    | otherwise = error "empty products do not have a nonbranching representation"

instance Num c => CanScale c (CommutatorType, Sum (Scaled c g)) where
  scale z (tp, terms) = (tp, scale z terms)

instance Algebra SpinGeneratorType where
  nonDiagonalCommutatorType = Commutator
  isIdentity g = g == SpinIdentity
  isDiagonal g = case g of
    SpinIdentity -> True
    SpinZ -> True
    SpinPlus -> False
    SpinMinus -> False
  conjugateGenerator g = case g of
    SpinIdentity -> SpinIdentity
    SpinZ -> SpinZ
    SpinPlus -> SpinMinus
    SpinMinus -> SpinPlus
  commute ::
    forall c.
    Fractional c =>
    SpinGeneratorType ->
    SpinGeneratorType ->
    (CommutatorType, Sum (Scaled c SpinGeneratorType))
  commute a b
    | a == b || isIdentity a || isIdentity b = (Commutator, [])
    | a > b = (-1 :: c) `scale` commute b a
  commute SpinZ SpinPlus = (Commutator, [Scaled 2 SpinPlus])
  commute SpinZ SpinMinus = (Commutator, [Scaled (-2) SpinMinus])
  commute SpinPlus SpinMinus = (Commutator, [Scaled 1 SpinZ])
  commute _ _ = error "should never happen"
  simplifyOrderedProduct ::
    forall c.
    Fractional c =>
    SpinGeneratorType ->
    SpinGeneratorType ->
    Sum (Scaled c SpinGeneratorType)
  simplifyOrderedProduct SpinIdentity b = [Scaled 1 b]
  simplifyOrderedProduct SpinZ SpinZ = [Scaled 1 SpinIdentity]
  simplifyOrderedProduct SpinZ SpinPlus = [Scaled 1 SpinPlus]
  simplifyOrderedProduct SpinZ SpinMinus = [Scaled (-1) SpinMinus]
  simplifyOrderedProduct SpinPlus SpinPlus = []
  simplifyOrderedProduct SpinPlus SpinMinus =
    [ Scaled ((1 :: c) / 2) SpinIdentity,
      Scaled ((1 :: c) / 2) SpinZ
    ]
  simplifyOrderedProduct SpinMinus SpinMinus = []
  simplifyOrderedProduct _ _ = error "should never happened because the product is ordered"

instance Algebra FermionGeneratorType where
  nonDiagonalCommutatorType = Anticommutator
  isIdentity g = g == FermionIdentity
  isDiagonal g = case g of
    FermionIdentity -> True
    FermionCount -> True
    FermionCreate -> False
    FermionAnnihilate -> False
  conjugateGenerator g = case g of
    FermionIdentity -> FermionIdentity
    FermionCount -> FermionCount
    FermionCreate -> FermionAnnihilate
    FermionAnnihilate -> FermionCreate
  commute ::
    forall c.
    Fractional c =>
    FermionGeneratorType ->
    FermionGeneratorType ->
    (CommutatorType, Sum (Scaled c FermionGeneratorType))
  commute a b
    | a == b || isIdentity a || isIdentity b = (Commutator, [])
    | a > b = commute b a
  commute FermionCount FermionCreate = (Anticommutator, [Scaled 1 FermionCreate])
  commute FermionCount FermionAnnihilate = (Anticommutator, [Scaled 1 FermionAnnihilate])
  commute FermionCreate FermionAnnihilate = (Anticommutator, [Scaled 1 FermionIdentity])
  commute _ _ = error "should never happen"
  simplifyOrderedProduct ::
    forall c.
    Fractional c =>
    FermionGeneratorType ->
    FermionGeneratorType ->
    Sum (Scaled c FermionGeneratorType)
  simplifyOrderedProduct FermionIdentity b = [Scaled 1 b]
  simplifyOrderedProduct FermionCount FermionCount = [Scaled 1 FermionCount]
  simplifyOrderedProduct FermionCount FermionCreate = [Scaled 1 FermionCreate]
  simplifyOrderedProduct FermionCount FermionAnnihilate = []
  simplifyOrderedProduct FermionCreate FermionCreate = []
  simplifyOrderedProduct FermionCreate FermionAnnihilate = [Scaled 1 FermionCount]
  simplifyOrderedProduct FermionAnnihilate FermionAnnihilate = []
  simplifyOrderedProduct _ _ = error "should never happened because the product is ordered"

instance (Algebra g, Ord i) => Algebra (Generator i g) where
  nonDiagonalCommutatorType = nonDiagonalCommutatorType @g
  isIdentity (Generator _ g) = isIdentity g
  isDiagonal (Generator _ g) = isDiagonal g
  conjugateGenerator (Generator i g) = Generator i (conjugateGenerator g)
  commute (Generator i a) (Generator j b)
    | i == j = case commute a b of
        (tp, ts) -> (tp, pack <$> ts)
    -- TODO: the following is probably unsafe, but it does work for both spins and fermions
    | isDiagonal a || isDiagonal b = (Commutator, [])
    | otherwise = (nonDiagonalCommutatorType @g, [])
    where
      pack = fmap (Generator i)
  simplifyOrderedProduct (Generator i a) (Generator j b)
    | i == j = fmap (fmap (Generator i)) $ simplifyOrderedProduct a b
    | otherwise = error "cannot simplify product of operators on different sites"

instance Num c => CanScale c (Scaled c g) where
  scale c (Scaled c' g) = Scaled (c * c') g

instance Traversable Product where
  traverse f (Product v) = Product <$> traverse f v

instance Pretty g => Pretty (Product g) where
  pretty (Product v) =
    Pretty.encloseSep mempty mempty " " (G.toList $ fmap pretty v)

instance Pretty g => Pretty (Sum g) where
  pretty (Sum v) =
    Pretty.encloseSep mempty mempty " + " (G.toList $ fmap pretty v)

instance CanScale c g => CanScale c (Sum g) where
  scale c (Sum v) = Sum $ G.map (scale c) v

instance Traversable Sum where
  traverse f (Sum v) = Sum <$> traverse f v

instance Num c => Num (Sum (Scaled c (Product g))) where
  (+) = (<>)
  (Sum a) * (Sum b) = Sum $ multiply <$> a <*> b
    where
      multiply (Scaled c g) (Scaled c' g') = Scaled (c * c') (g <> g')
  negate = scale (-1 :: c)
  abs = fmap (\(Scaled c g) -> Scaled (abs c) g)
  signum _ = error "Num instance of Sum does not implement signum"
  fromInteger _ = error "Num instance of Sum does not implement fromInteger"

sortVectorBy :: G.Vector v a => (a -> a -> Ordering) -> v a -> v a
sortVectorBy comp v = runST $ do
  buffer <- G.thaw v
  Data.Vector.Algorithms.Intro.sortBy comp buffer
  G.unsafeFreeze buffer

-- | Swaps generators at positions @i@ and @i+1@
swapGenerators ::
  (Fractional c, Algebra g) =>
  Int ->
  Product g ->
  Polynomial c g
swapGenerators i (Product v) = newTerms
  where
    !before = G.take i v
    !a = v ! i
    !b = v ! (i + 1)
    !after = G.drop (i + 2) v
    combined c terms =
      [Scaled c (Product $ before <> [b, a] <> after)]
        <> fmap (\(Scaled zᵢ gᵢ) -> Scaled zᵢ (Product $ before <> [gᵢ] <> after)) terms
    newTerms = case commute a b of
      (Commutator, terms) -> combined 1 terms
      (Anticommutator, terms) -> combined (-1) terms

-- | Reorder terms in the product.
--
-- Since the (anti)commutator is not always zero, the result is a polynomial rather than a monomial.
productToCanonical ::
  forall c g.
  (Fractional c, Algebra g) =>
  Product g ->
  Polynomial c g
productToCanonical t₀ = go (Scaled 1 t₀) 0 False
  where
    go :: Scaled c (Product g) -> Int -> Bool -> Sum (Scaled c (Product g))
    go x@(Scaled c t@(Product v)) !i !keepGoing
      | i < G.length v - 1 =
          case compare (v ! i) (v ! (i + 1)) of
            GT -> c `scale` foldMap (\x' -> go x' i True) (swapGenerators i t)
            _ -> go x (i + 1) keepGoing
      | keepGoing = go x 0 False
      | otherwise = Sum (G.singleton x)

-- | Reduce the degree of a product of terms that belong to the same site.
simplifyProductNoIndices :: forall c g. (Fractional c, Algebra g) => Product g -> Sum (Scaled c g)
simplifyProductNoIndices (Product v) = case G.toList v of
  [] -> []
  (g : gs) -> Sum $ go (Scaled 1 g) gs
  where
    go :: Scaled c g -> [g] -> Vector (Scaled c g)
    go r [] = [r]
    go (Scaled c g₁) (g₂ : gs) = do
      -- NOTE: this will fail if g₁ and g₂ belong to different sites or have different spin index.
      (Scaled c' g') <- unSum $ simplifyOrderedProduct g₁ g₂
      go (Scaled (c * c') g') gs

-- | Simplify a product of terms
simplifyProduct ::
  forall c g i.
  (Fractional c, Algebra g, Ord i) =>
  Product (Generator i g) ->
  Polynomial c (Generator i g)
simplifyProduct =
  fmap (fmap dropRedundantIdentities)
    . expandProduct
    . fromList
    . fmap (simplifyProductNoIndices . fromList)
    . List.groupBy (\(Generator i _) (Generator j _) -> i == j)
    . toList
  where
    dropRedundantIdentities (Product v)
      | not (G.null useful) = Product useful
      | otherwise = Product $ G.take 1 identities
      where
        (identities, useful) = G.partition (\(Generator _ g) -> isIdentity g) v

collectTerms :: (Num c, Eq c, Eq g) => Sum (Scaled c g) -> Sum (Scaled c g)
collectTerms = Sum . dropZeros . combine . unSum
  where
    dropZeros = G.filter (\(Scaled c _) -> c /= 0)
    combine =
      combineNeighbors
        (\(Scaled _ g) (Scaled _ g') -> g == g')
        (\(Scaled c g) (Scaled c' _) -> Scaled (c + c') g)

foldScaled ::
  CanScale c b =>
  (a -> b) ->
  Scaled c a ->
  b
foldScaled f (Scaled c p) = scale c (f p)

reorderTerms :: Ord g => Sum (Scaled c g) -> Sum (Scaled c g)
reorderTerms (Sum v) = Sum $ sortVectorBy ordering v
  where
    ordering (Scaled _ a) (Scaled _ b) = compare a b

collectIdentities ::
  forall c g i.
  (Eq c, Fractional c, Algebra g, Ord i) =>
  Polynomial c (Generator i g) ->
  Polynomial c (Generator i g)
collectIdentities (Sum terms)
  | G.null identities = Sum $ otherTerms
  | otherwise = Sum $ otherTerms `G.snoc` megaIdentity
  where
    isId (Scaled _ (Product p)) = G.length p == 1 && isIdentity (G.head p)
    (identities, otherTerms) = G.partition isId terms
    megaIdentity = G.foldl1' (\(Scaled c g) (Scaled c' _) -> (Scaled (c + c') g)) identities

simplifyPolynomial ::
  forall c g i.
  (Eq c, Fractional c, Algebra g, Ord i) =>
  Polynomial c (Generator i g) ->
  Polynomial c (Generator i g)
simplifyPolynomial =
  collectTerms
    . reorderTerms
    . collectIdentities
    . termsToCanonical

termsToCanonical ::
  forall c g i.
  (Eq c, Fractional c, Algebra g, Ord i) =>
  Polynomial c (Generator i g) ->
  Polynomial c (Generator i g)
termsToCanonical =
  foldMap (foldScaled simplifyProduct)
    . foldMap (foldScaled (productToCanonical @c))

expandProduct ::
  forall c g.
  Num c =>
  Product (Sum (Scaled c g)) ->
  Sum (Scaled c (Product g))
expandProduct (Product v)
  | G.null v = Sum G.empty
  | otherwise = G.foldl1' (*) $ fmap asProducts v
  where
    -- asProducts :: Sum (Scaled c g) -> Sum (Scaled c (Product g))
    asProducts = fmap (fmap (Product . G.singleton)) -- \(Scaled c g) -> (Scaled c (Product (G.singleton g))))

conjugatePolynomial ::
  forall c g i.
  (Eq c, Fractional c, Algebra g, Ord i) =>
  Polynomial c (Generator i g) ->
  Polynomial c (Generator i g)
conjugatePolynomial = simplifyPolynomial . fmap conjugateScaled
  where
    conjugateProduct = Product . G.reverse . G.map conjugateGenerator . unProduct
    conjugateScaled = fmap conjugateProduct

collectIndices :: Ord i => Polynomial c (Generator i g) -> [i]
collectIndices = List.nub . List.sort . collectSum
  where
    collectSum (Sum v) = concatMap collectScaled (G.toList v)
    collectScaled (Scaled _ p) = collectProduct p
    collectProduct (Product v) = (\(Generator i _) -> i) <$> (G.toList v)

mapGenerators :: (g -> g) -> Polynomial c g -> Polynomial c g
mapGenerators f = fmap (fmap (fmap f))

mapGeneratorsM :: Monad m => (g -> m g) -> Polynomial c g -> m (Polynomial c g)
mapGeneratorsM f =
  fmap Sum
    . G.mapM (\(Scaled c p) -> Scaled c . Product <$> G.mapM f (unProduct p))
    . unSum

foldlGenerators' :: (a -> g -> a) -> a -> Polynomial c g -> a
foldlGenerators' combine x₀ (Sum s) =
  G.foldl' (\ !x (Scaled _ (Product p)) -> G.foldl' combine x p) x₀ s

mapIndices :: (IndexType t -> IndexType t) -> Expr t -> Expr t
mapIndices f = Expr . mapGenerators (\(Generator i g) -> Generator (f i) g) . unExpr

mapIndicesM :: Monad m => (IndexType t -> m (IndexType t)) -> Expr t -> m (Expr t)
mapIndicesM f = fmap Expr . mapGeneratorsM (\(Generator i g) -> Generator <$> (f i) <*> pure g) . unExpr

mapCoeffs :: (ComplexRational -> ComplexRational) -> Expr t -> Expr t
mapCoeffs f = Expr . fmap (\(Scaled c p) -> Scaled (f c) p) . unExpr

foldlCoeffs' :: (a -> c -> a) -> a -> Polynomial c g -> a
foldlCoeffs' combine x₀ (Sum s) =
  G.foldl' (\ !x (Scaled c _) -> combine x c) x₀ s

simplifyExpr :: (Algebra (GeneratorType t), Ord (IndexType t)) => Expr t -> Expr t
simplifyExpr = Expr . simplifyPolynomial . unExpr

conjugateExpr :: (Algebra (GeneratorType t), Ord (IndexType t)) => Expr t -> Expr t
conjugateExpr = simplifyExpr . Expr . mapGenerators conjugateGenerator . unExpr . mapCoeffs conjugate

isIdentityExpr :: Algebra (GeneratorType t) => Expr t -> Bool
isIdentityExpr = isIdentitySum . unExpr
  where
    isIdentitySum (Sum terms) = G.length terms == 1 && isIdentityScaled (G.head terms)
    isIdentityScaled (Scaled _ (Product p)) = G.length p == 1 && isIdentityGenerator (G.head p)
    isIdentityGenerator (Generator _ g) = isIdentity g

isHermitianExpr :: (Ord (IndexType t), Eq (GeneratorType t), Algebra (GeneratorType t)) => Expr t -> Bool
isHermitianExpr terms = terms == conjugateExpr terms

isRealExpr :: Expr t -> Bool
isRealExpr = foldlCoeffs' (\f c -> f && imagPart c == 0) True . unExpr

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
  SomeExpr :: !(ParticleTag t) -> !(Expr t) -> SomeExpr

withSomeExpr ::
  forall c a.
  (c 'SpinTy, c 'SpinlessFermionTy, c 'SpinfulFermionTy) =>
  SomeExpr ->
  (forall t. c t => Expr t -> a) ->
  a
withSomeExpr (SomeExpr SpinTag a) f = f a
withSomeExpr (SomeExpr SpinfulFermionTag a) f = f a
withSomeExpr (SomeExpr SpinlessFermionTag a) f = f a

forSomeExpr ::
  forall c a.
  (c 'SpinTy, c 'SpinlessFermionTy, c 'SpinfulFermionTy) =>
  (forall t. c t => Expr t -> a) ->
  SomeExpr ->
  a
forSomeExpr f expr = withSomeExpr @c expr f

mapSomeExpr ::
  forall c a.
  (c 'SpinTy, c 'SpinlessFermionTy, c 'SpinfulFermionTy) =>
  (forall t. c t => Expr t -> Expr t) ->
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

class Num (Expr t) => ProvidesNum t

instance Num (Expr t) => ProvidesNum t

instance Num SomeExpr where
  (+) = binaryOp (+)
  (-) = binaryOp (-)
  (*) = binaryOp (*)
  negate = mapSomeExpr @ProvidesNum negate
  abs = mapSomeExpr @ProvidesNum abs
  signum = mapSomeExpr @ProvidesNum signum
  fromInteger _ =
    error $
      "Num instance of SomeExpr does not implement fromInteger; "
        <> "consider constructing an explicit identity 𝟙₀ and then scaling it"

class CanScale ComplexRational (Expr t) => ProvidesCanScale t

instance CanScale ComplexRational (Expr t) => ProvidesCanScale t

instance CanScale ComplexRational SomeExpr where
  scale z = mapSomeExpr @ProvidesCanScale (scale z)

class Pretty (Expr t) => ProvidesPretty t

instance Pretty (Expr t) => ProvidesPretty t

instance ToJSON SomeExpr where
  toJSON x@(SomeExpr tag _) =
    object $
      [ "particle" .= particleTagToType tag,
        "expression" .= withSomeExpr @ProvidesPretty x toPrettyText
      ]

tableFromLowLevelMapping :: (i ~ IndexType t, Ord i) => ParticleTag t -> Int -> Ptr CInt -> Ptr CInt -> IO (Map i i)
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

{-
replaceIndices :: Ord i => Polynomial c (Generator i g) -> [(i, i)] -> Polynomial c (Generator i g)
replaceIndices poly map = replaceSum poly
  where
    replaceSum = fmap replaceScaled
    replaceScaled = fmap replaceProduct
    replaceProduct = fmap replaceGenerator
    replaceGenerator (Generator i g) = Generator i' g
      where
        i' = case [to | (from, to) <- map, from == i] of
          [to] -> to
          [] -> error "index missing in mapping"
          _ -> error "multiple indices found in mapping"

forIndices :: Ord i => Polynomial c (Generator i g) -> [[i]] -> Polynomial c (Generator i g)
forIndices poly indices = mconcat (fmap processOne indices)
  where
    processOne newIndices = replaceIndices poly (zipWith (,) symbols newIndices)
    symbols = collectIndices poly
-}

replaceSiteIndices :: HasSiteIndex i => Polynomial c (Generator i g) -> [(Int, Int)] -> Polynomial c (Generator i g)
replaceSiteIndices poly indexMap = replaceSum poly
  where
    replaceSum = fmap replaceScaled
    replaceScaled = fmap replaceProduct
    replaceProduct = fmap replaceGenerator
    replaceGenerator (Generator i g) = Generator (mapSiteIndex f i) g
      where
        f i' = case [to | (from, to) <- indexMap, from == getSiteIndex i'] of
          [to] -> to
          [] -> error "index missing in mapping"
          _ -> error "multiple indices found in mapping"

forSiteIndices :: (Ord i, HasSiteIndex i) => Polynomial c (Generator i g) -> [[Int]] -> Polynomial c (Generator i g)
forSiteIndices poly indices = mconcat (fmap processOne indices)
  where
    processOne newIndices = replaceSiteIndices poly (zipWith (,) symbols newIndices)
    symbols = fmap getSiteIndex $ collectIndices poly

{-
groupTerms ::
  forall c i g.
  (HasMatrixRepresentation g, Algebra g, Ord i, Num c) =>
  i ->
  Polynomial c (Generator i g) ->
  Sum (LoweredOperator i (Polynomial c g))
groupTerms i₀ = step4 . step3 . step2 . step1
  where
    step1 ::
      Sum (Scaled c (Product (Generator i g))) ->
      Sum (Scaled c (LoweredOperator i (Product g)))
    step1 = fmap (\(Scaled c x) -> Scaled c (lower0 i₀ x))
    step2 ::
      Sum (Scaled c (LoweredOperator i (Product g))) ->
      Sum (Scaled c (LoweredOperator i (Product g)))
    step2 (Sum v) = Sum $ sortVectorBy withoutOrderingCoefficients v
    step3 ::
      Sum (Scaled c (LoweredOperator i (Product g))) ->
      Sum (LoweredOperator i (Polynomial c g))
    step3 = fmap (\(Scaled c (LoweredOperator i s g)) -> LoweredOperator i s [Scaled c g])
    step4 ::
      Sum (LoweredOperator i (Polynomial c g)) ->
      Sum (LoweredOperator i (Polynomial c g))
    step4 (Sum v) =
      Sum $
        combineNeighbors
          (\(LoweredOperator i1 s1 _) (LoweredOperator i2 s2 _) -> i1 == i2 && s1 == s2)
          (\(LoweredOperator i1 s1 g1) (LoweredOperator _ _ g2) -> (LoweredOperator i1 s1 (g1 + g2)))
          v
-}

-- lowerToMatrix ::
--   forall c g.
--   (HasMatrixRepresentation g, Algebra g, ComplexFloating c) =>
--   Polynomial c g ->
--   CsrMatrix
-- lowerToMatrix = sumToMatrix
--   where
--     productToMatrix (Product v) = csrKronMany $ csrMatrixFromDense @Vector . matrixRepresentation <$> G.toList v
--     scaledToMatrix (Scaled c p) = csrScale (toComplexDouble c) (productToMatrix p)
--     sumToMatrix (Sum v) = G.foldl1' (+) (G.map scaledToMatrix v)
