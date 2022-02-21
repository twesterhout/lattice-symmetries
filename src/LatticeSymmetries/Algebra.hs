{-# LANGUAGE DefaultSignatures #-}
{-# LANGUAGE FunctionalDependencies #-}
{-# LANGUAGE OverloadedLists #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE NoMonomorphismRestriction #-}

module LatticeSymmetries.Algebra
  ( SpinGeneratorType (..),
    FermionGeneratorType (..),
    HasMatrixRepresentation (..),
    Generator (..),
    Scaled (..),
    Sum (..),
    Product (..),
    scale,
    -- sumToCanonical,
    expandProduct,
    simplify,
  )
where

import Control.Exception (assert)
import Control.Monad.ST
-- import Data.Complex
import qualified Data.List as List
-- import Data.Ratio ((%))
import Data.Vector (Vector)
-- import qualified Data.Vector
import qualified Data.Vector.Algorithms.Intro
import Data.Vector.Generic ((!))
import qualified Data.Vector.Generic as G
-- import qualified Data.Vector.Generic.Mutable as GM
import GHC.Exts (IsList (..))
-- import qualified GHC.Show as GHC
import LatticeSymmetries.Sparse (DenseMatrix, KnownDenseMatrix, combineNeighbors, denseDot, denseMatMul)
import qualified LatticeSymmetries.Sparse as Sparse
import Prelude hiding (Product, Sum, identity, toList)

data SpinGeneratorType = SpinIdentity | SpinZ | SpinPlus | SpinMinus
  deriving stock (Eq, Ord, Show, Enum, Bounded, Generic)

data FermionGeneratorType = FermionIdentity | FermionCount | FermionCreate | FermionAnnihilate
  deriving stock (Eq, Ord, Show, Enum, Bounded, Generic)

data Generator i g = Generator !i !g
  deriving stock (Eq, Ord, Show, Generic)

data Scaled c g = Scaled !c !g
  deriving stock (Eq, Ord, Show, Generic)

withoutOrderingCoefficients :: Ord g => Scaled c g -> Scaled c g -> Ordering
withoutOrderingCoefficients (Scaled _ a) (Scaled _ b) = compare a b

-- class HasIdentity g where
--   identity :: g

-- instance HasIdentity SpinGeneratorType where
--   identity = SpinIdentity

-- instance HasIdentity FermionGeneratorType where
--   identity = FermionIdentity

class KnownNat r => HasMatrixRepresentation g r | g -> r where
  matrixRepresentation :: (G.Vector v c, Num c) => g -> DenseMatrix v r r c

instance HasMatrixRepresentation SpinGeneratorType 2 where
  matrixRepresentation = spinMatrixRepresentation

instance HasMatrixRepresentation FermionGeneratorType 2 where
  matrixRepresentation = fermionMatrixRepresentation

spinMatrixRepresentation :: (G.Vector v c, Num c) => SpinGeneratorType -> DenseMatrix v 2 2 c
spinMatrixRepresentation g = fromList $ case g of
  SpinIdentity -> [[1, 0], [0, 1]]
  SpinZ -> [[1, 0], [0, -1]]
  SpinPlus -> [[0, 1], [0, 0]]
  SpinMinus -> [[0, 0], [1, 0]]

fermionMatrixRepresentation :: (G.Vector v c, Num c) => FermionGeneratorType -> DenseMatrix v 2 2 c
fermionMatrixRepresentation g = fromList $ case g of
  FermionIdentity -> [[1, 0], [0, 1]]
  FermionCount -> [[1, 0], [0, 0]]
  FermionCreate -> [[0, 1], [0, 0]]
  FermionAnnihilate -> [[0, 0], [1, 0]]

getValues :: (Enum g, Bounded g) => [g]
getValues = enumFromTo minBound maxBound

getBasisExpansion ::
  (Enum g, Bounded g, HasMatrixRepresentation g r, KnownDenseMatrix v r r c, Eq c, Fractional c) =>
  DenseMatrix v r r c ->
  [(c, g)]
getBasisExpansion m = filter (\(c, _) -> c /= 0) $ getFactor <$> getValues
  where
    getFactor g =
      let !a = matrixRepresentation g
          !dot = denseDot a m
          !norm = denseDot a a
       in (dot / norm, g)

getBasisExpansion' ::
  (Enum g, Bounded g, HasMatrixRepresentation g r, KnownDenseMatrix v r r c, Eq c, Fractional c) =>
  DenseMatrix v r r c ->
  Sum (Scaled c g)
getBasisExpansion' m =
  Sum
    . G.map (\(c, g) -> Scaled c g)
    . G.filter (\(!c, _) -> c /= 0)
    . G.map getFactor
    $ G.fromList getValues
  where
    getFactor !g =
      let !a = matrixRepresentation g
          !dot = denseDot a m
          !norm = denseDot a a
       in (dot / norm, g)

isIdentity :: forall g r. HasMatrixRepresentation g r => g -> Bool
isIdentity g = matrixRepresentation g == (Sparse.denseEye :: DenseMatrix Vector r r Rational)

data CommutatorType = Commutator | Anticommutator
  deriving stock (Show, Eq)

data CommutationResult c g = CommutationResult [(c, g)] | AnticommutationResult [(c, g)]
  deriving stock (Show, Eq)

-- deriving anyclass (Functor)

class Ord g => Algebra g where
  nonDiagonalCommutatorType :: CommutatorType

  isDiagonal :: g -> Bool
  default isDiagonal :: forall r. HasMatrixRepresentation g r => g -> Bool
  isDiagonal = Sparse.isDenseMatrixDiagonal @Vector @r @Rational . matrixRepresentation

  commute :: (Eq c, Fractional c) => g -> g -> CommutationResult c g
  default commute ::
    forall c r.
    (Enum g, Bounded g, HasMatrixRepresentation g r, Eq c, Fractional c) =>
    g ->
    g ->
    CommutationResult c g
  commute a b = case nonDiagonalCommutatorType @g of
    Commutator -> CommutationResult commutator
    Anticommutator -> AnticommutationResult anticommutator
    where
      -- isDiagonal a || isDiagonal b = CommutationResult []
      aMatrix :: DenseMatrix Vector r r c
      aMatrix = matrixRepresentation a
      bMatrix = matrixRepresentation b
      ab = aMatrix `denseMatMul` bMatrix
      ba = bMatrix `denseMatMul` aMatrix
      commutator = getBasisExpansion $ ab - ba
      anticommutator = getBasisExpansion $ ab + ba

instance Algebra SpinGeneratorType where
  nonDiagonalCommutatorType = Commutator

instance Algebra FermionGeneratorType where
  nonDiagonalCommutatorType = Anticommutator

instance (Algebra g, HasMatrixRepresentation g r, Ord i) => Algebra (Generator i g) where
  isDiagonal (Generator _ g) = isDiagonal g
  nonDiagonalCommutatorType = nonDiagonalCommutatorType @g
  commute (Generator i a) (Generator j b)
    | i /= j && (isDiagonal a || isDiagonal b) = CommutationResult []
    | i /= j = case nonDiagonalCommutatorType @g of
      Commutator -> CommutationResult []
      Anticommutator -> AnticommutationResult []
    | otherwise = case commute a b of
      CommutationResult ts -> CommutationResult (pack <$> ts)
      AnticommutationResult ts -> AnticommutationResult (pack <$> ts)
    where
      pack (c, g) = (c, Generator i g)

-- infix 9 :*:

class Num c => CanScale c a where
  scale :: c -> a -> a

instance Num c => CanScale c (Scaled c g) where
  scale c (Scaled c' g) = Scaled (c * c') g

newtype Product g = Product (Vector g)
  deriving stock (Eq, Ord, Show, Generic)
  deriving newtype (IsList, Functor, Semigroup, Foldable)

instance Traversable Product where
  traverse f (Product v) = Product <$> traverse f v

newtype Sum g = Sum (Vector g)
  deriving stock (Eq, Show, Generic)
  deriving newtype (IsList, Functor, Semigroup, Monoid, Foldable)

instance CanScale c g => CanScale c (Sum g) where
  scale c (Sum v) = Sum $ G.map (scale c) v

instance Traversable Sum where
  traverse f (Sum v) = Sum <$> traverse f v

instance Num c => Num (Sum (Scaled c (Product g))) where
  (+) = (<>)
  (Sum a) * (Sum b) = Sum . G.fromList $ multiply <$> G.toList a <*> G.toList b
    where
      multiply (Scaled c g) (Scaled c' g') = Scaled (c * c') (g <> g')
  negate = scale (-1 :: c)
  abs = fmap (\(Scaled c g) -> Scaled (abs c) g)
  signum _ = error "Num instance of Sum does not implement signum"
  fromInteger _ = error "Num instance of Sum does not implement fromInteger"

-- data Monomial c g = !c :*: !(Vector g)
--   deriving stock (Eq, Ord, Show, Generic)

-- data SignedMonomial c g = SignedMonomial !c !(Vector g) ![(i, i)]
--   deriving stock (Eq, Ord, Show, Generic)

-- instance (Ord c, Ord g) => Ord (Monomial c g) where
--   compare (c₁ :*: t₁) (c₂ :*: t₂) = compare (t₁, c₁) (t₂, c₂)

-- instance Functor (Monomial c) where
--   fmap f (c :*: v) = c :*: (fmap f v)

-- instance Num c => IsList (Monomial c g) where
--   type Item (Monomial c g) = g
--   fromList gs = 1 :*: (G.fromList gs)
--   toList = error "IsList instance of Monomial does not implement toList"

-- monomialMapCoeff :: (c -> d) -> Monomial c g -> Monomial d g
-- monomialMapCoeff f (c :*: g) = f c :*: g

-- instance (Num c, HasIdentity g) => HasIdentity (Monomial c g) where
--   identity = 1 :*: G.singleton identity

-- newtype Polynomial c g = Polynomial (Vector (Monomial c g))
--   deriving stock (Eq, Show, Generic)

-- instance Num c => IsList (Polynomial c g) where
--   type Item (Polynomial c g) = (Monomial c g)
--   fromList = Polynomial . G.fromList
--   toList = error "IsList instance of Polynomial does not implement toList"

-- instance Functor (Polynomial c) where
--   fmap f (Polynomial v) = Polynomial $ fmap (fmap f) v

-- polynomialMapCoeff :: (c -> d) -> Polynomial c g -> Polynomial d g
-- polynomialMapCoeff f (Polynomial v) = Polynomial $ G.map (monomialMapCoeff f) v

-- instance (Num c, HasIdentity g) => HasIdentity (Polynomial c g) where
--   identity = Polynomial . G.singleton $ identity

-- identityPolynomial :: (HasMatrixRepresentation g r, Fractional c) => Polynomial c g
-- identityPolynomial = undefined

-- scaleMonomial :: Num c => c -> Monomial c g -> Monomial c g
-- scaleMonomial c = monomialMapCoeff (c *)

-- scalePolynomial :: Num c => c -> Polynomial c g -> Polynomial c g
-- scalePolynomial c = polynomialMapCoeff (c *)

-- instance (Num c) => Num (Polynomial c g) where
--   (+) (Polynomial a) (Polynomial b) = Polynomial (a G.++ b)
--   (*) (Polynomial a) (Polynomial b) = Polynomial . G.fromList $ (<>) <$> G.toList a <*> G.toList b
--   negate a = (-1) `scalePolynomial` a
--   abs = mapPolynomial (\(z :*: g) -> abs z :*: g)
--   signum = error "Num instance of Polynomial does not implement signum"
--   fromInteger _ = error "Num instance of Polynomial does not implement fromInteger"

-- monomialDropIdentities :: (Num c, HasIdentity g, Eq g) => Monomial c g -> Monomial c g
-- monomialDropIdentities (c :*: v)
--   | G.null v' = (c :*: G.singleton identity)
--   | otherwise = (c :*: v')
--   where
--     v' = G.filter (/= identity) v

-- dropIdentities :: (Num c, HasIdentity g, Eq g) => Polynomial c g -> Polynomial c g
-- dropIdentities (Polynomial v) = Polynomial $ G.map monomialDropIdentities v

-- dropZeros :: (Num c, Eq c) => Polynomial c g -> Polynomial c g
-- dropZeros (Polynomial v) = Polynomial $ G.filter (\(c :*: g) -> c /= 0) v

sortVectorBy :: G.Vector v a => (a -> a -> Ordering) -> v a -> v a
sortVectorBy comp v = runST $ do
  buffer <- G.thaw v
  Data.Vector.Algorithms.Intro.sortBy comp buffer
  G.unsafeFreeze buffer

-- sortMonomials :: (Ord c, Ord g) => Polynomial c g -> Polynomial c g
-- sortMonomials (Polynomial v) = Polynomial $
--   runST $ do
--     buffer <- G.thaw v
--     Data.Vector.Algorithms.Intro.sort buffer
--     G.unsafeFreeze buffer

-- groupMonomials :: (Num c, Eq c, Eq g) => Polynomial c g -> Polynomial c g
-- groupMonomials (Polynomial v) =
--   dropZeros . Polynomial $
--     combineNeighbors
--       (\(_ :*: a) (_ :*: b) -> a == b)
--       (\(a :*: t) (b :*: _) -> (a + b) :*: t)
--       v

-- commute ::
--   forall g c r.
--   (Enum g, Bounded g, HasMatrixRepresentation g r, Algebra g, Eq c, Fractional c) =>
--   g ->
--   g ->
--   [(c, g)]
-- commute a b =
--   getBasisExpansion $
--     (aMatrix `denseMatMul` bMatrix) - c * (bMatrix `denseMatMul` aMatrix)
--   where
--     (aMatrix :: DenseMatrix Vector r r c) = matrixRepresentation a
--     bMatrix = matrixRepresentation b
--     c = case commutatorTypeOf (Proxy @g) of
--       Commutator -> 1
--       Anticommutator -> -1

-- | Swaps generators at positions i and i + 1
swapGenerators' ::
  (Fractional c, Eq c, Algebra g) =>
  Int ->
  Product g ->
  Sum (Scaled c (Product g))
swapGenerators' i (Product v) = Sum (G.fromList newTerms)
  where
    !before = G.take i v
    !a = v ! i
    !b = v ! (i + 1)
    !after = G.drop (i + 2) v
    combined c terms =
      Scaled c (Product $ before <> [b, a] <> after) :
        [Scaled zᵢ (Product $ before <> [gᵢ] <> after) | (zᵢ, gᵢ) <- terms]
    newTerms = case commute a b of
      CommutationResult terms -> combined 1 terms
      AnticommutationResult terms -> combined (-1) terms

productToCanonical ::
  forall c g.
  (Fractional c, Eq c, Algebra g) =>
  Product g ->
  Sum (Scaled c (Product g))
productToCanonical t₀ = go (Scaled 1 t₀) 0 False
  where
    go :: Scaled c (Product g) -> Int -> Bool -> Sum (Scaled c (Product g))
    go x@(Scaled c t@(Product v)) i keepGoing
      | i < G.length v - 1 =
        case compare (v ! i) (v ! (i + 1)) of
          GT ->
            let newTerms = swapGenerators' i t
             in scale c $ foldMap (\x' -> go x' i True) newTerms
          _ -> go x (i + 1) keepGoing
      | keepGoing = go x 0 False
      | otherwise = Sum (G.singleton x)

-- sumToCanonical ::
--   (Fractional c, Eq c, Algebra g) =>
--   Sum (Scaled c (Product g)) ->
--   Sum (Scaled c (Product g))
-- sumToCanonical = order . foldMap single
--   where
--     single (Scaled c p) = c `scale` productToCanonical p
--     order (Sum v) = Sum $ sortVectorBy withoutOrderingCoefficients v

simplifyPrimitiveProduct' ::
  forall c g r.
  (HasCallStack, Enum g, Bounded g, HasMatrixRepresentation g r, Fractional c, Eq c) =>
  Product g ->
  Sum (Scaled c g)
simplifyPrimitiveProduct' (Product !v)
  | G.null v = []
  | otherwise = getBasisExpansion' $ G.foldl' combine acc₀ v
  where
    acc₀ :: DenseMatrix Vector r r c
    acc₀ = Sparse.denseEye
    combine !acc !g = acc `denseMatMul` (matrixRepresentation g)

-- in trace ("combining " <> show acc <> " with " <> show g <> " --> " <> show r) r
-- basis x = let y = getBasisExpansion' x
--    in trace ("basis expansion of " <> show x <> " is " <> show y) y

dropIdentities' :: HasMatrixRepresentation g r => Product (Generator i g) -> Product (Generator i g)
dropIdentities' (Product v)
  | not (G.null useful) = Product useful
  | otherwise = Product $ G.take 1 identities
  where
    (identities, useful) = G.partition (\(Generator _ g) -> isIdentity g) v

simplifyProduct ::
  forall c g i r.
  (HasCallStack, Enum g, Bounded g, HasMatrixRepresentation g r, Fractional c, Eq c, Eq i) =>
  Product (Generator i g) ->
  Sum (Scaled c (Product (Generator i g)))
simplifyProduct =
  fmap (\(Scaled c p) -> Scaled c (dropIdentities' p))
    . expandProduct
    . fromList
    . fmap (simplifyOneSite . fromList)
    . List.groupBy (\(Generator i _) (Generator j _) -> i == j)
    . toList
  where
    unpack (Generator _ g) = g
    pack i (Scaled c g) = Scaled c (Generator i g)
    simplifyOneSite :: HasCallStack => Product (Generator i g) -> Sum (Scaled c (Generator i g))
    simplifyOneSite p@(Product !v)
      | G.null v = assert False []
      | otherwise =
        let (Generator i _) = G.head v
         in fmap (pack i) . simplifyPrimitiveProduct' . fmap unpack $ p

dropZeros' :: (Num c, Eq c) => Sum (Scaled c g) -> Sum (Scaled c g)
dropZeros' (Sum v) = Sum $ G.filter (\(Scaled c _) -> c /= 0) v

combineFactors' :: (Num c, Eq g) => Sum (Scaled c g) -> Sum (Scaled c g)
combineFactors' (Sum v) =
  Sum $
    combineNeighbors
      (\(Scaled _ g) (Scaled _ g') -> g == g')
      (\(Scaled c g) (Scaled c' _) -> Scaled (c + c') g)
      v

simplify ::
  forall c g r i.
  (HasCallStack, Fractional c, Eq c, Enum g, Bounded g, HasMatrixRepresentation g r, Algebra g, Ord i) =>
  Sum (Scaled c (Product (Generator i g))) ->
  Sum (Scaled c (Product (Generator i g)))
simplify =
  dropZeros'
    . combineFactors'
    . reorder
    . foldMap (withScaled simplifyProduct)
    . foldMap (withScaled productToCanonical)
  where
    withScaled ::
      (Product (Generator i g) -> Sum (Scaled c (Product (Generator i g)))) ->
      Scaled c (Product (Generator i g)) ->
      Sum (Scaled c (Product (Generator i g)))
    withScaled f (Scaled c p) = scale c (f p)
    reorder (Sum v) = Sum $ sortVectorBy withoutOrderingCoefficients v

---
--- monomialDropIdentities :: HasMatrixRepresentation g r => Monomial c (Generator i g) -> Monomial c (Generator i g)
--- monomialDropIdentities (c :*: v)
---   | not (G.null useful) = c :*: useful
---   | otherwise = c :*: (G.take 1 identities)
---   where
---     (identities, useful) = G.partition (\(Generator _ g) -> isIdentity g) v
---
-- simplifyMonomial ::
--   (Fractional c, Ord c, Ord i, Enum g, Bounded g, Ord g, HasMatrixRepresentation g r) =>
--   Monomial c (Generator i g) ->
--   Polynomial c (Generator i g)
-- simplifyMonomial (z :*: v) =
--   -- groupMonomials
--   -- . sortMonomials
--   scalePolynomial z
--     . polynomialMapCoeff fromRational
--     . mapPolynomial monomialDropIdentities
--     . product
--     . fmap simplifyOneSite
--     . List.groupBy (\(Generator i _) (Generator j _) -> i == j)
--     . G.toList
--     $ v

expandProduct ::
  forall c g.
  Num c =>
  Product (Sum (Scaled c g)) ->
  Sum (Scaled c (Product g))
expandProduct (Product v)
  | G.null v = Sum G.empty
  | otherwise = G.foldl1' (*) $ fmap asProducts v
  where
    asProducts :: Sum (Scaled c g) -> Sum (Scaled c (Product g))
    asProducts = fmap (\(Scaled c g) -> (Scaled c (Product (G.singleton g))))

-- | Swaps generators at positions i and i + 1
-- swapGenerators ::
--   forall c g.
--   (Fractional c, Eq c, Algebra g) =>
--   Int ->
--   Monomial c g ->
--   Polynomial c g
-- swapGenerators i (z :*: v) =
--   dropZeros
--     . scalePolynomial z
--     . Polynomial
--     . G.fromList
--     $ newTerms
--   where
--     before = G.take i v
--     !a = v ! i
--     !b = v ! (i + 1)
--     after = G.drop (i + 2) v
--     combined c terms =
--       c :*: (G.concat [before, G.singleton b, G.singleton a, after]) :
--         [zᵢ :*: (G.concat [before, G.singleton gᵢ, after]) | (zᵢ, gᵢ) <- terms]
--     newTerms = case commute a b of
--       CommutationResult terms -> combined 1 terms
--       AnticommutationResult terms -> combined (-1) terms

-- monomialToCanonical ::
--   forall c g r.
--   (Fractional c, Eq c, Ord c, Algebra g) =>
--   Monomial c g ->
--   Polynomial c g
-- monomialToCanonical t₀ = Polynomial $ go t₀ 0 False
--   where
--     go :: Monomial c g -> Int -> Bool -> Vector (Monomial c g)
--     go term@(_ :*: v) i keepGoing
--       | i < G.length v - 1 =
--         case compare (v ! i) (v ! (i + 1)) of
--           GT ->
--             let (Polynomial newTerms) = swapGenerators i term
--              in G.concatMap (\x -> go x i True) newTerms
--           _ -> go term (i + 1) keepGoing
--       | keepGoing = go term 0 False
--       | otherwise = G.singleton term

-- polynomialToCanonical ::
--   forall c g r.
--   (Fractional c, Eq c, Ord c, Algebra g) =>
--   Polynomial c g ->
--   Polynomial c g
-- polynomialToCanonical = sortMonomials . concatMapPolynomial monomialToCanonical

--
-- data SpinGenerator = SpinGenerator {sgId :: !Int, sgType :: !SpinGeneratorType}
--   deriving stock (Eq, Ord, Generic)
--
-- data FermionGenerator = FermionGenerator {fgId :: !Int, fgType :: !FermionGeneratorType}
--   deriving stock (Eq, Ord, Generic)

-- instance GHC.Show FermionGenerator where
--   show (FermionGenerator i t) = letter <> dagger <> index
--     where
--       letter = "c"
--       dagger = if t == FermionCreate then "†" else ""
--       index = toSubscript <$> (show i)
--       toSubscript c = case c of
--         '0' -> '₀'
--         '1' -> '₁'
--         '2' -> '₂'
--         '3' -> '₃'
--         '4' -> '₄'
--         '5' -> '₅'
--         '6' -> '₆'
--         '7' -> '₇'
--         '8' -> '₈'
--         '9' -> '₉'

-- instance Num c => AlgebraGenerator c FermionGenerator where
--   commute (FermionGenerator i t₁) (FermionGenerator j t₂)
--     | t₁ == t₂ = (-1, [], 0)
--     | i /= j = (-1, [], 0)
--     | otherwise = (-1, [], 1)

-- instance Num c => AlgebraGenerator c SpinGenerator where
--   commute (SpinGenerator i tᵢ) (SpinGenerator j tⱼ)
--     | i /= j = (1, [], 0)
--     | tᵢ == tⱼ = (1, [], 0)
--   commute (SpinGenerator i SpinPlus) (SpinGenerator _ SpinMinus) = (1, [(1, (SpinGenerator i SpinZ))], 0)
--   commute (SpinGenerator i SpinMinus) (SpinGenerator _ SpinPlus) = (1, [(-1, (SpinGenerator i SpinZ))], 0)
--   commute (SpinGenerator i SpinPlus) (SpinGenerator _ SpinZ) = (1, [(-2, (SpinGenerator i SpinPlus))], 0)
--   commute (SpinGenerator i SpinZ) (SpinGenerator _ SpinPlus) = (1, [(2, (SpinGenerator i SpinPlus))], 0)
--   commute (SpinGenerator i SpinMinus) (SpinGenerator _ SpinZ) = (1, [(2, (SpinGenerator i SpinMinus))], 0)
--   commute (SpinGenerator i SpinZ) (SpinGenerator _ SpinMinus) = (1, [(-2, (SpinGenerator i SpinMinus))], 0)

-- viaOrdered :: (Fractional c, Ord g) => (g -> g -> (c, [(c, g)])) -> g -> g -> (c, [(c, g)])
-- viaOrdered f a b = case compare a b of
--   EQ -> (1, [])
--   LT -> f a b
--   -- f a b computes ab - c*ba = X
--   -- then ba - 1/c*ab = -1/c*X
--   GT ->
--     let (c, terms) = f b a
--      in (-1 / c, fmap (\(x, g) -> (-x / c, g)) terms)

-- instance Algebra SpinGeneratorType where
--   commute = viaOrdered $ \a b -> case (a, b) of
--     (SpinIdentity, _) -> (1, [])
--     (SpinZ, SpinPlus) -> (1, [(2, SpinPlus)])
--     (SpinZ, SpinMinus) -> (1, [(-2, SpinMinus)])
--     (SpinPlus, SpinMinus) -> (1, [(1, SpinZ)])
--     _ -> error "this should not have happened"
--   needsSignFlip _ = False
--
-- instance Algebra FermionGeneratorType where
--   commute = viaOrdered $ \a b -> case (a, b) of
--     -- [1, g] = 0 for any g
--     (FermionIdentity, _) -> (1, [])
--     -- {n, a†} = na† + a†n = (a†a)a† + a†(a†a) = (1 - aa†)a† + 0 = a†
--     (FermionCount, FermionCreate) -> (-1, [(1, FermionCreate)])
--     -- {n, a} = na + an = (a†a)a + a(a†a) = 0 + a(1 - aa†) = a
--     (FermionCount, FermionAnnihilate) -> (-1, [(1, FermionAnnihilate)])
--     -- {a†, a} = 1
--     (FermionCreate, FermionAnnihilate) -> (-1, [(1, FermionIdentity)])
--     _ -> error "this should not have happened"
--   needsSignFlip g = case g of
--     FermionIdentity -> False
--     FermionCount -> False
--     FermionCreate -> True
--     FermionAnnihilate -> True

-- commute (Generator ia a) (Generator ib b)
--   | ia == ib = let (c, terms) = commute a b in (c, liftToIndex ia <$> terms)
--   | otherwise = (1, [])
--   where
--     liftToIndex i (z, g) = (z, Generator i g)
-- needsSignFlip (Generator _ g) = needsSignFlip g

-- instance Ord i => Algebra (Generator i FermionGeneratorType) where
--   commute (Generator ia a) (Generator ib b)
--     | ia == ib = let (c, terms) = commute a b in (c, liftToIndex ia <$> terms)
--     -- NOTE: This is tricky. In general, we have that for fermionic operators on different sites:
--     -- {a, b} = ab + ba = 0.
--     -- However, if now either a or b is identity 1, then we want to use the commutator:
--     -- [a, 1] = [1, b] = 0
--     | a == FermionIdentity || b == FermionIdentity = (1, [])
--     | otherwise = (-1, [])
--     where
--       liftToIndex i (z, g) = (z, Generator i g)
--   needsSignFlip (Generator _ g) = needsSignFlip g

-- instance Num c => Semigroup (Monomial c g) where
--   (<>) (ca :*: ga) (cb :*: gb) = (ca * cb) :*: (G.concat [ga, gb])

-- mapPolynomial :: (Monomial c1 g1 -> Monomial c2 g2) -> Polynomial c1 g1 -> Polynomial c2 g2
-- mapPolynomial f (Polynomial v) = Polynomial $ G.map f v

-- concatMapPolynomial :: (Monomial c1 g1 -> Polynomial c2 g2) -> Polynomial c1 g1 -> Polynomial c2 g2
-- concatMapPolynomial f (Polynomial v) = Polynomial $ G.concatMap f' v
--   where
--     f' x = let (Polynomial y) = f x in y

-- simplifyPrimitiveProduct ::
--   forall g r.
--   (Enum g, Bounded g, HasMatrixRepresentation g r) =>
--   [g] ->
--   [(Rational, g)]
-- simplifyPrimitiveProduct v
--   | null v = []
--   | otherwise = getBasisExpansion $ foldl' combine acc₀ v
--   where
--     acc₀ :: DenseMatrix Vector r r Rational
--     acc₀ = Sparse.denseEye
--     combine !acc !g = acc `denseMatMul` (matrixRepresentation g)

-- simplifyOneSite ::
--   forall g i r.
--   (Enum g, Bounded g, HasMatrixRepresentation g r) =>
--   [Generator i g] ->
--   Polynomial Rational (Generator i g)
-- simplifyOneSite gs@((Generator i _) : _) =
--   Polynomial
--     . G.fromList
--     . fmap pack
--     . simplifyPrimitiveProduct
--     . fmap unpack
--     $ gs
--   where
--     unpack (Generator _ g) = g
--     pack (c, g) = c :*: G.singleton (Generator i g)
-- simplifyOneSite [] = []

-- monomialDropIdentities :: HasMatrixRepresentation g r => Monomial c (Generator i g) -> Monomial c (Generator i g)
-- monomialDropIdentities (c :*: v)
--   | not (G.null useful) = c :*: useful
--   | otherwise = c :*: (G.take 1 identities)
--   where
--     (identities, useful) = G.partition (\(Generator _ g) -> isIdentity g) v

-- simplifyMonomial ::
--   (Fractional c, Ord c, Ord i, Enum g, Bounded g, Ord g, HasMatrixRepresentation g r) =>
--   Monomial c (Generator i g) ->
--   Polynomial c (Generator i g)
-- simplifyMonomial (z :*: v) =
--   -- groupMonomials
--   -- . sortMonomials
--   scalePolynomial z
--     . polynomialMapCoeff fromRational
--     . mapPolynomial monomialDropIdentities
--     . product
--     . fmap simplifyOneSite
--     . List.groupBy (\(Generator i _) (Generator j _) -> i == j)
--     . G.toList
--     $ v

-- simplifyPolynomial ::
--   (Fractional c, Ord c, Ord g, Enum g, Bounded g, HasMatrixRepresentation g r, Ord i) =>
--   Polynomial c (Generator i g) ->
--   Polynomial c (Generator i g)
-- simplifyPolynomial (Polynomial v) = groupMonomials . sortMonomials . sum $ simplifyMonomial <$> (G.toList v)

data SignMask i = SignMask ![(i, i)]
  deriving stock (Show, Eq, Ord)

getSignMask ::
  (HasMatrixRepresentation g r, Algebra g, Ord i) =>
  i ->
  Product (Generator i g) ->
  SignMask i
getSignMask iₘᵢₙ (Product gs)
  | G.null gs' = SignMask []
  | otherwise = SignMask . reverse $ go [] (G.length edges `mod` 2 == 1) (G.toList edges)
  where
    needsSignMask = not . isDiagonal
    gs' = G.map (\(Generator i _) -> i) . G.filter needsSignMask $ gs
    edges = G.cons iₘᵢₙ gs'
    go acc _ [] = acc
    go acc !flag (i : j : rest) =
      let acc' = if flag then (i, j) : acc else acc
       in go acc' (not flag) (j : rest)
    go _ _ [_] = error "this should not have happened"

data LoweredOperator i g = LoweredOperator !(Vector i) !(SignMask i) !g
  deriving stock (Show, Eq, Ord)

isStrictlyOrdered :: Ord i => Vector i -> Bool
isStrictlyOrdered v = G.all id $ G.zipWith (\i j -> i < j) v (G.drop 1 v)

extractIndices :: Ord i => Product (Generator i g) -> Vector i
extractIndices (Product v)
  | isStrictlyOrdered is = is
  | otherwise = error "indices are not strictly ordered"
  where
    is = G.map (\(Generator i _) -> i) v

extractGenerators :: Product (Generator i g) -> Product g
extractGenerators (Product v) = Product $ G.map (\(Generator _ g) -> g) v

lower0 ::
  (HasMatrixRepresentation g r, Algebra g, Ord i) =>
  i ->
  Product (Generator i g) ->
  LoweredOperator i (Product g)
lower0 i₀ p = LoweredOperator (extractIndices p) (getSignMask i₀ p) (extractGenerators p)

type Polynomial c g = Sum (Scaled c (Product g))

groupTerms ::
  forall c i g r.
  (HasMatrixRepresentation g r, Algebra g, Ord i, Num c) =>
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

-- lower = undefined

-- spinSimplifyProduct ::
--   (HasCallStack, Fractional c) =>
--   SpinGeneratorType ->
--   SpinGeneratorType ->
--   Polynomial c SpinGeneratorType
-- spinSimplifyProduct a b = case (a, b) of
--   -- We simply drop identities: g⋅1 = 1⋅g = g
--   (SpinIdentity, _) -> [1 :*: [b]]
--   (_, SpinIdentity) -> [1 :*: [a]]
--   -- σᶻ⋅σ⁺ = σ⁺
--   (SpinZ, SpinPlus) -> [1 :*: [SpinPlus]]
--   (SpinPlus, SpinZ) -> [(-1) :*: [SpinPlus]]
--   -- σᶻ⋅σ⁻ = -σ⁻
--   (SpinZ, SpinMinus) -> [(-1) :*: [SpinMinus]]
--   (SpinMinus, SpinZ) -> [1 :*: [SpinMinus]]
--   -- σ⁺⋅σ⁻ = 1/2 (1 + σᶻ)
--   (SpinPlus, SpinMinus) -> [half :*: [SpinIdentity], half :*: [SpinZ]]
--   (SpinMinus, SpinPlus) -> [half :*: [SpinIdentity], (-half) :*: [SpinZ]]
--   -- σ⁺⋅σ⁺ = 0
--   (SpinPlus, SpinPlus) -> []
--   -- σ⁻⋅σ⁻ = 0
--   (SpinMinus, SpinMinus) -> []
--   _ -> error "this should not have happened"
--   where
--     half = fromRational (1 % 2)

-- fermionSimplifyProduct ::
--   (HasCallStack, Num c) =>
--   FermionGeneratorType ->
--   FermionGeneratorType ->
--   Polynomial c FermionGeneratorType
-- fermionSimplifyProduct a b = case (a, b) of
--   -- We simply drop identities: g⋅1 = 1⋅g = g
--   (FermionIdentity, _) -> [1 :*: [b]]
--   (_, FermionIdentity) -> [1 :*: [a]]
--   -- n⋅a† = (a†⋅a)⋅a† = (1 - a⋅a†)⋅a† = a†
--   (FermionCount, FermionCreate) -> [1 :*: [FermionCreate]]
--   -- n⋅a = (a†⋅a)⋅a = 0
--   (FermionCount, FermionAnnihilate) -> []
--   -- a†⋅n = a†⋅(a†⋅a) = 0
--   (FermionCreate, FermionCount) -> []
--   -- a⋅n = a⋅(a†⋅a) = a⋅(1 - a⋅a†) = a
--   (FermionAnnihilate, FermionCount) -> [1 :*: [FermionAnnihilate]]
--   -- a†⋅a = n
--   (FermionCreate, FermionAnnihilate) -> [1 :*: [FermionCount]]
--   -- a⋅a† = 1 - n
--   (FermionAnnihilate, FermionCreate) -> [identity, (-1) :*: [FermionCount]]
--   -- a†⋅a† = 0
--   (FermionCreate, FermionCreate) -> []
--   -- a⋅a = 0
--   (FermionAnnihilate, FermionAnnihilate) -> []
--   _ -> error "this should not have happened"

-- go (c :*: v)
--   | G.length v == 0 = []
--   | G.length v == 1 = [c :*: v]
--   | otherwise =
--     let a = v ! 0
--         b = v ! 1
--         suffix = c :*: G.drop 2 v
--      in concatMapPolynomial go $
--           mapPolynomial (<> suffix) (simplifyProduct a b)

-- spinSimplifyProduct ::
--   Fractional c =>
--   SpinGeneratorType ->
--   SpinGeneratorType ->
--   Expr c SpinGeneratorType
-- spinSimplifyProduct (SpinGenerator ia ta) (SpinGenerator ib tb)
--   | ia == ib && ta <= tb = case (ta, tb) of
--     (SpinPlus, SpinPlus) -> zeroExpr
--     (SpinPlus, SpinMinus) -> [half :*: IdentityTerm, half :*: [SpinGenerator ia SpinZ]]
--     (SpinPlus, SpinZ) -> [(-1) :*: [SpinGenerator ia SpinPlus]]
--     -- (SpinMinus, SpinPlus) -> [half :*: IdentityTerm, (-half) :*: [SpinGenerator ia SpinZ]]
--     (SpinMinus, SpinMinus) -> zeroExpr
--     (SpinMinus, SpinZ) -> [1 :*: [SpinGenerator ia SpinMinus]]
--     -- (SpinZ, SpinPlus) -> [1 :*: [SpinGenerator ia SpinPlus]]
--     -- (SpinZ, SpinMinus) -> [(-1) :*: [SpinGenerator ia SpinMinus]]
--     (SpinZ, SpinZ) -> [1 :*: IdentityTerm]
--     _ -> error "this should not have happened"
--   | ta > tb = error "expected generators in canonical order"
--   | otherwise = error "expected generators on the same site"
--   where
--     half = fromRational (1 % 2)

-- data Term g = Term !(Vector g) | IdentityTerm
--   deriving stock (Eq, Ord, Show, Generic)
--
-- instance IsList (Term g) where
--   type Item (Term g) = g
--   fromList = Term . G.fromList
--   toList = error "IsList instance of Term does not implement toList"
--
-- instance IsList (Expr c g) where
--   type Item (Expr c g) = ScaledTerm c g
--   fromList = Expr . G.fromList
--   toList = error "IsList instance of Expr does not implement toList"

-- tOrder :: Term g -> Int
-- tOrder (Term v) = G.length v
--
-- data Expr c g = Expr !(Vector (ScaledTerm c g))
--   deriving stock (Eq, Show, Generic)
--
-- zeroExpr :: Expr c g
-- zeroExpr = Expr G.empty
--
-- scaleTerm :: Num c => c -> ScaledTerm c g -> ScaledTerm c g
-- scaleTerm z (c :*: g) = (z * c) :*: g
--
-- scaleExpr :: Num c => c -> Expr c g -> Expr c g
-- scaleExpr z (Expr v) = Expr $ G.map (scaleTerm z) v
--
-- fromTerms :: Num c => [Term g] -> Expr c g
-- fromTerms ts = Expr . G.fromList $ (1 :*:) <$> ts

-- type ℂ = Complex Double

-- | Specifies the algebra
--
-- g_α g_β - c g_β g_α = F_αβ + \sum_γ f^γ_ab g_γ
--
-- @commute ga gb@ returns a tuple of @(c, [(f^γ_ab, g_γ)], F_αβ)@
-- class Ord g => AlgebraGenerator c g where
--   commute :: g -> g -> (c, [(c, g)], c)
