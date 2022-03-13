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
    Algebra (..),
    Polynomial,
    scale,
    -- sumToCanonical,
    expandProduct,
    simplify,
    groupTerms,
    -- lowerToMatrix,
    forIndices,
  )
where

import Control.Exception (assert)
import Control.Monad.ST
import Data.Bits
import qualified Data.List as List
import qualified Data.Text as Text
import Data.Vector (Vector)
import qualified Data.Vector.Algorithms.Intro
import Data.Vector.Generic ((!))
import qualified Data.Vector.Generic as G
import GHC.Exts (IsList (..))
import LatticeSymmetries.BitString
import LatticeSymmetries.CSR
import LatticeSymmetries.ComplexRational
import LatticeSymmetries.Dense
import LatticeSymmetries.NonbranchingTerm
import LatticeSymmetries.Sparse (combineNeighbors)
import Numeric.Natural
import Text.PrettyPrint.ANSI.Leijen (Pretty (..))
import qualified Text.PrettyPrint.ANSI.Leijen as Pretty
import Prelude hiding (Product, Sum, identity, toList)

data SpinGeneratorType = SpinIdentity | SpinZ | SpinPlus | SpinMinus
  deriving stock (Eq, Ord, Show, Enum, Bounded, Generic)

instance Pretty SpinGeneratorType where
  pretty x = case x of
    SpinIdentity -> "1"
    SpinZ -> "σᶻ"
    SpinPlus -> "σ⁺"
    SpinMinus -> "σ⁻"

data FermionGeneratorType = FermionIdentity | FermionCount | FermionCreate | FermionAnnihilate
  deriving stock (Eq, Ord, Show, Enum, Bounded, Generic)

instance Pretty FermionGeneratorType where
  pretty x = case x of
    FermionIdentity -> "1"
    FermionCount -> "n"
    FermionCreate -> "c†"
    FermionAnnihilate -> "c"

data Generator i g = Generator !i !g
  deriving stock (Eq, Ord, Show, Generic)

toSubscript :: Int -> Text
toSubscript n = Text.map h (show n)
  where
    h '0' = '₀'
    h '1' = '₁'
    h '2' = '₂'
    h '3' = '₃'
    h '4' = '₄'
    h '5' = '₅'
    h '6' = '₆'
    h '7' = '₇'
    h '8' = '₈'
    h '9' = '₉'

instance Pretty g => Pretty (Generator Int g) where
  pretty (Generator i g) = pretty g <> Pretty.text (Text.unpack (toSubscript i))

data Scaled c g = Scaled !c !g
  deriving stock (Eq, Ord, Show, Generic)

instance Functor (Scaled c) where
  fmap f (Scaled c g) = Scaled c (f g)

instance (Pretty c, Pretty g) => Pretty (Scaled c g) where
  pretty (Scaled c g) = pretty c <> Pretty.text " × " <> pretty g

withoutOrderingCoefficients :: Ord g => Scaled c g -> Scaled c g -> Ordering
withoutOrderingCoefficients (Scaled _ a) (Scaled _ b) = compare a b

-- class HasIdentity g where
--   identity :: g

-- instance HasIdentity SpinGeneratorType where
--   identity = SpinIdentity

-- instance HasIdentity FermionGeneratorType where
--   identity = FermionIdentity

class HasMatrixRepresentation g where
  matrixRepresentation :: (G.Vector v c, Num c) => g -> DenseMatrix v c

instance HasMatrixRepresentation SpinGeneratorType where
  matrixRepresentation = spinMatrixRepresentation

instance HasMatrixRepresentation FermionGeneratorType where
  matrixRepresentation = fermionMatrixRepresentation

spinMatrixRepresentation :: (G.Vector v c, Num c) => SpinGeneratorType -> DenseMatrix v c
spinMatrixRepresentation g = fromList $ case g of
  SpinIdentity -> [[1, 0], [0, 1]]
  SpinZ -> [[1, 0], [0, -1]]
  SpinPlus -> [[0, 1], [0, 0]]
  SpinMinus -> [[0, 0], [1, 0]]

fermionMatrixRepresentation :: (G.Vector v c, Num c) => FermionGeneratorType -> DenseMatrix v c
fermionMatrixRepresentation g = fromList $ case g of
  FermionIdentity -> [[1, 0], [0, 1]]
  FermionCount -> [[1, 0], [0, 0]]
  FermionCreate -> [[0, 1], [0, 0]]
  FermionAnnihilate -> [[0, 0], [1, 0]]

instance HasNonbranchingRepresentation (Generator Int SpinGeneratorType) where
  nonbranchingRepresentation (Generator _ SpinIdentity) =
    NonbranchingTerm 1 zeroBits zeroBits zeroBits zeroBits zeroBits
  nonbranchingRepresentation (Generator i SpinZ) =
    NonbranchingTerm (-1) zeroBits zeroBits zeroBits zeroBits (bit i)
  nonbranchingRepresentation (Generator i SpinPlus) =
    NonbranchingTerm 1 (bit i) (bit i) zeroBits (bit i) zeroBits
  nonbranchingRepresentation (Generator i SpinMinus) =
    NonbranchingTerm 1 (bit i) zeroBits (bit i) (bit i) zeroBits

instance HasNonbranchingRepresentation (Generator Int FermionGeneratorType) where
  nonbranchingRepresentation (Generator _ FermionIdentity) =
    NonbranchingTerm 1 zeroBits zeroBits zeroBits zeroBits zeroBits
  nonbranchingRepresentation (Generator i FermionCount) =
    NonbranchingTerm 1 (bit i) (bit i) (bit i) zeroBits zeroBits
  nonbranchingRepresentation (Generator i FermionCreate) =
    NonbranchingTerm 1 (bit i) (bit i) zeroBits (bit i) (BitString (bit i - 1))
  nonbranchingRepresentation (Generator i FermionAnnihilate) =
    NonbranchingTerm 1 (bit i) zeroBits (bit i) (bit i) (BitString (bit i - 1))

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

getValues :: (Enum g, Bounded g) => [g]
getValues = enumFromTo minBound maxBound

getBasisExpansion ::
  (Enum g, Bounded g, HasMatrixRepresentation g, G.Vector v c, Eq c, Fractional c) =>
  DenseMatrix v c ->
  [(c, g)]
getBasisExpansion m = filter (\(c, _) -> c /= 0) $ getFactor <$> getValues
  where
    getFactor g =
      let !a = matrixRepresentation g
          !dot = denseDot a m
          !norm = denseDot a a
       in (dot / norm, g)

getBasisExpansion' ::
  (Enum g, Bounded g, HasMatrixRepresentation g, G.Vector v c, Eq c, Fractional c) =>
  DenseMatrix v c ->
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

isIdentity :: forall g. HasMatrixRepresentation g => g -> Bool
isIdentity g = m == (denseEye (dmRows m) :: DenseMatrix Vector Rational)
  where
    m = matrixRepresentation g

data CommutatorType = Commutator | Anticommutator
  deriving stock (Show, Eq)

data CommutationResult c g = CommutationResult [(c, g)] | AnticommutationResult [(c, g)]
  deriving stock (Show, Eq)

-- deriving anyclass (Functor)

class Ord g => Algebra g where
  nonDiagonalCommutatorType :: CommutatorType

  isDiagonal :: g -> Bool
  default isDiagonal :: HasMatrixRepresentation g => g -> Bool
  isDiagonal = denseIsDiagonal @Vector @Rational . matrixRepresentation

  commute :: (Eq c, Fractional c) => g -> g -> CommutationResult c g
  default commute ::
    forall c.
    (Enum g, Bounded g, HasMatrixRepresentation g, Eq c, Fractional c) =>
    g ->
    g ->
    CommutationResult c g
  commute a b = case nonDiagonalCommutatorType @g of
    Commutator -> CommutationResult commutator
    Anticommutator -> AnticommutationResult anticommutator
    where
      -- isDiagonal a || isDiagonal b = CommutationResult []
      aMatrix :: DenseMatrix Vector c
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

instance (Algebra g, HasMatrixRepresentation g, Ord i) => Algebra (Generator i g) where
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

instance Pretty g => Pretty (Product g) where
  pretty (Product v) =
    Pretty.encloseSep mempty mempty (Pretty.char ' ') (G.toList $ fmap pretty v)

newtype Sum g = Sum (Vector g)
  deriving stock (Eq, Show, Generic)
  deriving newtype (IsList, Functor, Semigroup, Monoid, Foldable)

instance Pretty g => Pretty (Sum g) where
  pretty (Sum v) =
    Pretty.encloseSep mempty mempty (Pretty.text " + ") (G.toList $ fmap pretty v)

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
  forall c g.
  (HasCallStack, Enum g, Bounded g, HasMatrixRepresentation g, Fractional c, Eq c) =>
  Product g ->
  Sum (Scaled c g)
simplifyPrimitiveProduct' (Product !v)
  | G.null v = []
  | otherwise = getBasisExpansion' $ G.foldl' combine acc₀ v
  where
    acc₀ :: DenseMatrix Vector c
    acc₀ = denseEye $ dmRows (matrixRepresentation (G.head v) :: DenseMatrix Vector c)
    combine !acc !g = acc `denseMatMul` (matrixRepresentation g)

-- in trace ("combining " <> show acc <> " with " <> show g <> " --> " <> show r) r
-- basis x = let y = getBasisExpansion' x
--    in trace ("basis expansion of " <> show x <> " is " <> show y) y

dropIdentities' :: HasMatrixRepresentation g => Product (Generator i g) -> Product (Generator i g)
dropIdentities' (Product v)
  | not (G.null useful) = Product useful
  | otherwise = Product $ G.take 1 identities
  where
    (identities, useful) = G.partition (\(Generator _ g) -> isIdentity g) v

simplifyProduct ::
  forall c g i.
  (HasCallStack, Enum g, Bounded g, HasMatrixRepresentation g, Fractional c, Eq c, Eq i) =>
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
  forall c g i.
  (HasCallStack, Fractional c, Eq c, Enum g, Bounded g, HasMatrixRepresentation g, Algebra g, Ord i) =>
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

data SignMask i = SignMask ![(i, i)]
  deriving stock (Show, Eq, Ord)

getSignMask ::
  (HasMatrixRepresentation g, Algebra g, Ord i) =>
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

data LoweredTerm i g = LoweredTerm !(Vector i) !(SignMask i) !g
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
  (HasMatrixRepresentation g, Algebra g, Ord i) =>
  i ->
  Product (Generator i g) ->
  LoweredOperator i (Product g)
lower0 i₀ p = LoweredOperator (extractIndices p) (getSignMask i₀ p) (extractGenerators p)

type Polynomial c g = Sum (Scaled c (Product g))

collectIndices :: Ord i => Polynomial c (Generator i g) -> [i]
collectIndices = List.nub . List.sort . collectSum
  where
    collectSum (Sum v) = concatMap collectScaled (G.toList v)
    collectScaled (Scaled _ p) = collectProduct p
    collectProduct (Product v) = (\(Generator i _) -> i) <$> (G.toList v)

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
