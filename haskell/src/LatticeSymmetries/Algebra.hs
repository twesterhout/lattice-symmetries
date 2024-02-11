{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE InstanceSigs #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE OverloadedLists #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE NoMonomorphismRestriction #-}

module LatticeSymmetries.Algebra
  ( -- * Polynomials
    Scaled (..)
  , Product (..)
  , Sum (..)
  , Polynomial
  , CanScale (..)
  , CommutatorType (..)
  , Algebra (..)
  , simplifyPolynomial
  , swapGenerators
  , HasProperGeneratorType
  , HasProperIndexType
  , IsBasis
  )
where

import Data.List qualified as List
import Data.Vector (Vector)
import Data.Vector.Fusion.Bundle qualified as Bundle (inplace)
import Data.Vector.Fusion.Bundle.Size (toMax)
import Data.Vector.Fusion.Stream.Monadic (Step (..), Stream (..))
import Data.Vector.Generic ((!))
import Data.Vector.Generic qualified as G
import GHC.Exts (IsList (..))
import LatticeSymmetries.ComplexRational
import LatticeSymmetries.Generator
import LatticeSymmetries.NonbranchingTerm
import Prettyprinter (Pretty (..))
import Prelude hiding (Product, Sum, identity, toList)
import LatticeSymmetries.Utils (sortVectorBy)

-- | Represents a term of the form @c × g@ where @c@ is typically a scalar and @g@ is some
-- expression.
data Scaled c g = Scaled !c !g
  deriving stock (Eq, Ord, Show, Generic)

instance Functor (Scaled c) where
  fmap f (Scaled c g) = Scaled c (f g)

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

-- | Specifies the type of commutator to use for an algebra.
data CommutatorType
  = -- | Commutator \([a, b] = ab - ba\)
    Commutator
  | -- | Anticommutator \(\{a, b\} = ab + ba\)
    Anticommutator
  deriving stock (Show, Eq)

class (Ord g, HasIdentity g) => Algebra g where
  -- | The type of commutator this algebra uses.
  --
  -- We have 'Commutator' for spin (or bosonic) systems and 'Anticommutator' for fermionic systems.
  nonDiagonalCommutatorType :: CommutatorType

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
        Done -> pure Done
    step' (CombineNeighborsPrevious a, s) = do
      r <- step s
      case r of
        Yield b s' ->
          if equal a b
            then pure $ Skip (CombineNeighborsPrevious (combine a b), s')
            else pure $ Yield a (CombineNeighborsPrevious b, s')
        Skip s' -> pure $ Skip (CombineNeighborsPrevious a, s')
        Done -> pure $ Yield a (CombineNeighborsDone, s)
    step' (CombineNeighborsDone, _) = pure Done

combineNeighbors :: G.Vector v a => (a -> a -> Bool) -> (a -> a -> a) -> v a -> v a
combineNeighbors equal combine =
  G.unstream . Bundle.inplace (combineNeighborsImpl equal combine) toMax . G.stream

instance HasNonbranchingRepresentation g => HasNonbranchingRepresentation (Scaled ComplexRational g) where
  nonbranchingRepresentation (Scaled c g) = t {nbtV = c * nbtV t}
    where
      t = nonbranchingRepresentation g

instance HasNonbranchingRepresentation g => HasNonbranchingRepresentation (Product g) where
  nonbranchingRepresentation (Product v)
    | not (G.null v) = G.foldl1' (<>) . G.map nonbranchingRepresentation $ v
    | otherwise =
        -- an empty Product is equivalent to an identity, and identities for all particle types are the same
        nonbranchingRepresentation (Generator (0 :: Int) SpinIdentity)

instance Num c => CanScale c (CommutatorType, Sum (Scaled c g)) where
  scale z (tp, terms) = (tp, scale z terms)

instance Algebra SpinGeneratorType where
  nonDiagonalCommutatorType = Commutator
  isDiagonal = \case
    SpinIdentity -> True
    SpinZ -> True
    SpinPlus -> False
    SpinMinus -> False
  conjugateGenerator = \case
    SpinIdentity -> SpinIdentity
    SpinZ -> SpinZ
    SpinPlus -> SpinMinus
    SpinMinus -> SpinPlus
  commute
    :: forall c
     . Fractional c
    => SpinGeneratorType
    -> SpinGeneratorType
    -> (CommutatorType, Sum (Scaled c SpinGeneratorType))
  commute a b
    | a == b || isIdentity a || isIdentity b = (Commutator, [])
    | a > b = (-1 :: c) `scale` commute b a
  commute SpinZ SpinPlus = (Commutator, [Scaled 2 SpinPlus])
  commute SpinZ SpinMinus = (Commutator, [Scaled (-2) SpinMinus])
  commute SpinPlus SpinMinus = (Commutator, [Scaled 1 SpinZ])
  commute _ _ = error "should never happen"
  simplifyOrderedProduct
    :: forall c
     . Fractional c
    => SpinGeneratorType
    -> SpinGeneratorType
    -> Sum (Scaled c SpinGeneratorType)
  simplifyOrderedProduct SpinIdentity b = [Scaled 1 b]
  simplifyOrderedProduct SpinZ SpinZ = [Scaled 1 SpinIdentity]
  simplifyOrderedProduct SpinZ SpinPlus = [Scaled 1 SpinPlus]
  simplifyOrderedProduct SpinZ SpinMinus = [Scaled (-1) SpinMinus]
  simplifyOrderedProduct SpinPlus SpinPlus = []
  simplifyOrderedProduct SpinPlus SpinMinus = [Scaled ((1 :: c) / 2) SpinIdentity, Scaled ((1 :: c) / 2) SpinZ]
  simplifyOrderedProduct SpinMinus SpinMinus = []
  simplifyOrderedProduct _ _ = error "should never happened because the product is ordered"

instance Algebra FermionGeneratorType where
  nonDiagonalCommutatorType = Anticommutator
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
  commute
    :: forall c
     . Fractional c
    => FermionGeneratorType
    -> FermionGeneratorType
    -> (CommutatorType, Sum (Scaled c FermionGeneratorType))
  commute a b
    | a == b || isIdentity a || isIdentity b = (Commutator, [])
    | a > b = commute b a
  commute FermionCount FermionCreate = (Anticommutator, [Scaled 1 FermionCreate])
  commute FermionCount FermionAnnihilate = (Anticommutator, [Scaled 1 FermionAnnihilate])
  commute FermionCreate FermionAnnihilate = (Anticommutator, [Scaled 1 FermionIdentity])
  commute _ _ = error "should never happen"
  simplifyOrderedProduct
    :: forall c
     . Fractional c
    => FermionGeneratorType
    -> FermionGeneratorType
    -> Sum (Scaled c FermionGeneratorType)
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
  isDiagonal (Generator _ g) = isDiagonal g
  conjugateGenerator (Generator i g) = Generator i (conjugateGenerator g)
  commute (Generator i a) (Generator j b)
    | i == j = let (tp, ts) = commute a b in (tp, pack <$> ts)
    -- TODO: the following is probably unsafe, but it does work for both spins and fermions
    -- Since i != j in this case, we are essentially computing [1⊗a, b⊗1],
    | isDiagonal a || isDiagonal b = (Commutator, [])
    | otherwise = (nonDiagonalCommutatorType @g, [])
    where
      pack = fmap (Generator i)
  simplifyOrderedProduct (Generator i a) (Generator j b)
    | i == j = fmap (Generator i) <$> simplifyOrderedProduct a b
    | otherwise = error "cannot simplify product of operators on different sites"

instance Num c => CanScale c (Scaled c g) where
  scale c (Scaled c' g) = Scaled (c * c') g

instance Traversable Product where
  traverse f (Product v) = Product <$> traverse f v

-- instance Pretty g => Pretty (Product g) where
--   pretty (Product v) =
--     Pretty.encloseSep mempty mempty " " (G.toList $ fmap pretty v)
--
-- instance Pretty g => Pretty (Sum g) where
--   pretty (Sum v)
--     | not (G.null v) = Pretty.encloseSep mempty mempty " + " (G.toList $ fmap pretty v)
--     | otherwise = pretty (0 :: Int)

instance CanScale c g => CanScale c (Sum g) where
  scale c (Sum v) = Sum $ G.map (scale c) v

instance Traversable Sum where
  traverse f (Sum v) = Sum <$> traverse f v

instance Num c => Num (Polynomial c g) where
  (+) = (<>)
  (Sum a) * (Sum b) = Sum $ multiply <$> a <*> b
    where
      multiply (Scaled c g) (Scaled c' g') = Scaled (c * c') (g <> g')
  negate = scale (-1 :: c)
  abs = fmap (\(Scaled c g) -> Scaled (abs c) g)
  signum _ = error "Num instance of Sum does not implement signum"
  fromInteger _ = error "Num instance of Sum does not implement fromInteger"

-- | Swaps generators at positions @i@ and @i+1@
swapGenerators :: (Fractional c, Algebra g) => Int -> Product g -> Polynomial c g
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
productToCanonical
  :: forall c g
   . (Fractional c, Algebra g)
  => Product g
  -> Polynomial c g
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
  [] -> error "simplifyProductNoIndices does not work on empty Products"
  (g : gs) -> Sum $ go (Scaled 1 g) gs
  where
    go :: Scaled c g -> [g] -> Vector (Scaled c g)
    go r [] = [r]
    go (Scaled c g₁) (g₂ : gs) = do
      -- NOTE: this will fail if g₁ and g₂ belong to different sites or have different spin index.
      (Scaled c' g') <- unSum $ simplifyOrderedProduct g₁ g₂
      go (Scaled (c * c') g') gs

-- | Simplify a product of terms
simplifyProduct
  :: forall c g i
   . (Fractional c, Algebra g, Ord i)
  => Product (Generator i g)
  -> Polynomial c (Generator i g)
simplifyProduct =
  fmap (fmap dropIdentities)
    . expandProduct
    . fromList
    . fmap (simplifyProductNoIndices . fromList)
    . List.groupBy (\(Generator i _) (Generator j _) -> i == j)
    . toList
  where
    dropIdentities (Product v) = Product $ G.filter (not . isIdentity) v

collectTerms :: (Num c, Eq c, Eq g) => Sum (Scaled c g) -> Sum (Scaled c g)
collectTerms = Sum . dropZeros . combine . unSum
  where
    dropZeros = G.filter (\(Scaled c _) -> c /= 0)
    combine =
      combineNeighbors
        (\(Scaled _ g) (Scaled _ g') -> g == g')
        (\(Scaled c g) (Scaled c' _) -> Scaled (c + c') g)

foldScaled
  :: CanScale c b
  => (a -> b)
  -> Scaled c a
  -> b
foldScaled f (Scaled c p) = scale c (f p)

reorderTerms :: Ord g => Sum (Scaled c g) -> Sum (Scaled c g)
reorderTerms (Sum v) = Sum $ sortVectorBy ordering v
  where
    ordering (Scaled _ a) (Scaled _ b) = compare a b

-- collectIdentities ::
--   forall c g i.
--   (Eq c, Fractional c, Algebra g, Ord i) =>
--   Polynomial c (Generator i g) ->
--   Polynomial c (Generator i g)
-- collectIdentities (Sum terms)
--   | G.null identities = Sum $ otherTerms
--   | otherwise = Sum $ otherTerms `G.snoc` megaIdentity
--   where
--     isId (Scaled _ (Product p)) = G.length p == 1 && isIdentity (G.head p)
--     (identities, otherTerms) = G.partition isId terms
--     megaIdentity = G.foldl1' (\(Scaled c g) (Scaled c' _) -> (Scaled (c + c') g)) identities

simplifyPolynomial
  :: forall c g i
   . (Eq c, Fractional c, Algebra g, Ord i)
  => Polynomial c (Generator i g)
  -> Polynomial c (Generator i g)
simplifyPolynomial =
  collectTerms
    . reorderTerms
    -- . collectIdentities
    . termsToCanonical

termsToCanonical
  :: forall c g i
   . (Eq c, Fractional c, Algebra g, Ord i)
  => Polynomial c (Generator i g)
  -> Polynomial c (Generator i g)
termsToCanonical =
  foldMap (foldScaled simplifyProduct)
    . foldMap (foldScaled (productToCanonical @c))

expandProduct
  :: forall c g
   . Num c
  => Product (Sum (Scaled c g))
  -> Sum (Scaled c (Product g))
expandProduct (Product v)
  | G.null v = Sum [Scaled 1 (Product [])]
  | otherwise = G.foldl1' (*) $ fmap asProducts v
  where
    -- asProducts :: Sum (Scaled c g) -> Sum (Scaled c (Product g))
    asProducts = fmap (fmap (Product . G.singleton)) -- \(Scaled c g) -> (Scaled c (Product (G.singleton g))))

class IsGeneratorType (GeneratorType t) => HasProperGeneratorType t

instance IsGeneratorType (GeneratorType t) => HasProperGeneratorType t

class IsIndexType (IndexType t) => HasProperIndexType t

instance IsIndexType (IndexType t) => HasProperIndexType t

type IsGeneratorType g = (Eq g, Ord g, Pretty g, Algebra g, HasNonbranchingRepresentation (Generator Int g))

type IsIndexType i = (Eq i, Ord i {-HasSiteIndex i,-}, Pretty i)

class
  ( Typeable t
  , HasProperGeneratorType t
  , HasProperIndexType t
  , Pretty (Generator (IndexType t) (GeneratorType t))
  ) =>
  IsBasis t

instance IsBasis 'SpinTy

instance IsBasis 'SpinfulFermionTy

instance IsBasis 'SpinlessFermionTy
