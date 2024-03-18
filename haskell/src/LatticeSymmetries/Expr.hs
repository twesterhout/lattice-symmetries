{-# LANGUAGE GADTs #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE OverloadedLists #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE NoMonomorphismRestriction #-}

module LatticeSymmetries.Expr
  ( Expr (..)
  , mkExpr
  , conjugateExpr
  , isHermitianExpr
  , isRealExpr
  , isIdentityExpr
  , simplifyExpr
  , mapGenerators
  , mapGeneratorsM
  , mapIndices
  , mapIndicesM
  , collectIndices
  , exprToHypergraph
  , exprPermutationGroup
  , fromSExpr
  , estimateNumberSites
  , isInvariantUponSpinInversion
  , conservesNumberParticles
  )
where

import Data.Set qualified as Set
import Data.Vector.Generic qualified as G
import LatticeSymmetries.Algebra
import LatticeSymmetries.Automorphisms
import LatticeSymmetries.ComplexRational
import LatticeSymmetries.Generator
import LatticeSymmetries.Parser
import LatticeSymmetries.Permutation (Permutation (..))
import Prettyprinter (Pretty (..))
import Prettyprinter qualified as Pretty
import Prelude hiding (Product, Sum, identity, toList)
import LatticeSymmetries.Utils (unique, sortVectorBy)
import Data.Vector (Vector)

newtype Expr t = Expr
  { unExpr :: Polynomial ComplexRational (Generator (IndexType t) (GeneratorType t))
  }
  deriving stock (Generic)

deriving instance (Eq (IndexType t), Eq (GeneratorType t)) => Eq (Expr t)

deriving instance (Show (IndexType t), Show (GeneratorType t)) => Show (Expr t)

mapGenerators
  :: (Generator (IndexType t1) (GeneratorType t1) -> Generator (IndexType t2) (GeneratorType t2))
  -> Expr t1
  -> Expr t2
mapGenerators f = Expr . fmap (fmap (fmap f)) . unExpr

mapGeneratorsM
  :: Monad m
  => (Generator (IndexType t1) (GeneratorType t1) -> m (Generator (IndexType t2) (GeneratorType t2)))
  -> Expr t1
  -> m (Expr t2)
mapGeneratorsM f =
  fmap (Expr . Sum)
    . G.mapM (\(Scaled c p) -> Scaled c . Product <$> G.mapM f (unProduct p))
    . unSum
    . unExpr

mapIndices :: (IndexType t -> IndexType t) -> Expr t -> Expr t
mapIndices f = mapGenerators (\(Generator i g) -> Generator (f i) g)

mapIndicesM :: Monad m => (IndexType t -> m (IndexType t)) -> Expr t -> m (Expr t)
mapIndicesM f = mapGeneratorsM (\(Generator i g) -> Generator <$> f i <*> pure g)

collectIndices :: Ord (IndexType t) => Expr t -> Vector (IndexType t)
collectIndices = unique . sortVectorBy compare . collectSum . unExpr
  where
    collectSum (Sum v) = G.concatMap collectScaled v
    collectScaled (Scaled _ p) = collectProduct p
    collectProduct (Product v) = (\(Generator i _) -> i) <$> v

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
      prettyProduct (Product v)
        | G.null v = "I"
        | otherwise = Pretty.encloseSep mempty mempty " " (G.toList $ fmap pretty v)
      prettySum [] = pretty (0 :: ComplexRational) <> " I"
      prettySum (term : others) = prettyScaled term <> prettySumRest others
      prettySumRest [] = ""
      prettySumRest (Scaled c p : others)
        | (realPart c == 0 && imagPart c < 0) || (imagPart c == 0 && realPart c < 0) =
            " - " <> prettyScaled (Scaled (-c) p) <> prettySumRest others
        | otherwise = " + " <> prettyScaled (Scaled c p) <> prettySumRest others

-- collectSiteIndices :: forall t. HasProperIndexType t => Expr t -> [Int]
-- collectSiteIndices = Set.toList . Set.fromList . fmap getSiteIndex . collectIndices
--
-- replicateSiteIndices
--   :: forall t
--    . (HasProperIndexType t, Algebra (GeneratorType t))
--   => [[Int]]
--   -> Expr t
--   -> Expr t
-- replicateSiteIndices newIndices expr =
--   case newIndices of
--     [] -> Expr []
--     (i : is) -> foldl' (+) (replace i) $ fmap replace is
--   where
--     oldSiteIndices = collectSiteIndices expr
--     k = length oldSiteIndices
--     replace siteIndices
--       | Map.size mapping == k = mapIndices (mapSiteIndex (mapping Map.!)) expr
--       | otherwise =
--           error
--             $ "wrong number of site indices: "
--             <> show (length siteIndices)
--             <> "; expected "
--             <> show k
--       where
--         mapping = Map.fromList (zip oldSiteIndices siteIndices)

fromSExpr :: forall t. (HasCallStack, Algebra (GeneratorType t), Ord (IndexType t)) =>
  ParticleTag t -> SExpr -> Either Text (Expr t)
fromSExpr t0 expr0 = simplifyExpr <$> go t0 expr0
  where
    go :: ParticleTag t -> SExpr -> Either Text (Expr t)
    go t (SSum terms) = do
      exprs <- mapM (go t) terms
      pure $ foldl' (\a b -> Expr $ unExpr a + unExpr b) (Expr []) exprs
    go t (SScaled c term) = scale c <$> go t term
    go t (SProduct (term :| terms)) = do
      expr <- go t term
      exprs <- mapM (go t) terms
      pure $ foldl' (\a b -> Expr $ unExpr a * unExpr b) expr exprs
    go _ (SPrimitive SIdentity) = pure $ Expr [Scaled 1 []]
    go SpinTag (SPrimitive (SSpinOp c t i)) = pure $
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
    go SpinfulFermionTag (SPrimitive (SFermionOp t (Just s) i)) = pure $
      case t of
        SFermionCreate -> Expr [Scaled 1 [Generator (s, i) FermionCreate]]
        SFermionAnnihilate -> Expr [Scaled 1 [Generator (s, i) FermionAnnihilate]]
        SFermionNumber -> Expr [Scaled 1 [Generator (s, i) FermionCount]]
    go SpinlessFermionTag (SPrimitive (SFermionOp t Nothing i)) = pure $
      case t of
        SFermionCreate -> Expr [Scaled 1 [Generator i FermionCreate]]
        SFermionAnnihilate -> Expr [Scaled 1 [Generator i FermionAnnihilate]]
        SFermionNumber -> Expr [Scaled 1 [Generator i FermionCount]]
    go SpinTag (SPrimitive _) = fail "expected an expression for spin-1/2 particles"
    go SpinfulFermionTag (SPrimitive _) = fail "expected an expression for spinful fermions"
    go SpinlessFermionTag (SPrimitive _) = fail "expected an expression for spinless fermions"

mkExpr :: (HasCallStack, Algebra (GeneratorType t), Ord (IndexType t)) => ParticleTag t -> Text -> Either Text (Expr t)
mkExpr tag s = fromSExpr tag =<< parseExprEither s

permuteExprG :: IsBasis t => (IndexType t -> Int) -> (Int -> IndexType t) -> Permutation -> Expr t -> Expr t
permuteExprG to from (unPermutation -> p) expr = simplifyExpr $ mapIndices remap expr
  where
    remap i = from (p G.! to i)

isInvariantG :: IsBasis t => (IndexType t -> Int) -> (Int -> IndexType t) -> Expr t -> Permutation -> Bool
isInvariantG to from expr permutation = permuteExprG to from permutation expr == expr

estimateNumberSites :: forall t. IsBasis t => Expr t -> Int
estimateNumberSites expr =
  1 + case particleDispatch @t of
    SpinTag -> m (collectIndices expr)
    SpinlessFermionTag -> m (collectIndices expr)
    SpinfulFermionTag -> m (snd <$> collectIndices expr)
  where
    m v | G.null v = 0
        | otherwise = G.maximum v

mapHypergraph :: (Ord a, Ord b) => (a -> b) -> Hypergraph a -> Hypergraph b
mapHypergraph f (Hypergraph vertices edges) =
  Hypergraph (Set.map f vertices) (Set.map (Set.map f) edges)

exprToHypergraph :: forall t. IsBasis t => Expr t -> Hypergraph (IndexType t)
exprToHypergraph e@(Expr (Sum terms)) = Hypergraph vertices (Set.fromList . fmap Set.fromList $ hyperedges)
  where
    n = estimateNumberSites e
    vertices :: Set (IndexType t)
    vertices = case particleDispatch @t of
      SpinTag -> Set.fromAscList [0 .. n]
      SpinlessFermionTag -> Set.fromAscList [0 .. n]
      SpinfulFermionTag -> Set.fromList [(s, i) | i <- [0 .. n], s <- [SpinUp, SpinDown]]
    hyperedges = G.toList $ fmap scaledToSet terms
    scaledToSet (Scaled _ p) = productToSet p
    productToSet (Product v) = G.toList $ (\(Generator i _) -> i) <$> v

exprPermutationGroup :: forall t. IsBasis t => Maybe Int -> Expr t -> PermutationGroup
exprPermutationGroup numberSites expr = hypergraphAutomorphisms (isInvariantG to from expr) hypergraph
  where
    !k = fromMaybe (estimateNumberSites expr) numberSites
    !to = flattenIndex @t k
    !from = unFlattenIndex @t k
    vertices = Set.fromAscList
      . enumFromTo 0
      $ case particleDispatch @t of
        SpinTag -> k - 1
        SpinlessFermionTag -> k - 1
        SpinfulFermionTag -> 2 * k - 1
    hypergraph = (mapHypergraph to (exprToHypergraph expr)) {vertices = vertices}

isInvariantUponSpinInversion :: forall t. IsBasis t => Expr t -> Bool
isInvariantUponSpinInversion expr = case particleDispatch @t of
  SpinTag ->
    let prim = \case
          SpinIdentity -> (1 :: Int, SpinIdentity)
          SpinZ -> (-1, SpinZ)
          SpinPlus -> (1, SpinMinus)
          SpinMinus -> (1, SpinPlus)
        gen (Generator i g) = let (c, g') = prim g in (c, Generator i g')
        prod (Scaled z (Product gs)) =
          let (cs', gs') = G.unzip (gen <$> gs) in Scaled (fromIntegral (G.product cs') * z) (Product gs')
        expr' = simplifyExpr . Expr @t . fmap prod . unExpr $ expr
     in expr' == expr
  SpinlessFermionTag -> True
  SpinfulFermionTag -> undefined

conservesNumberParticles :: forall t. IsBasis t => Expr t -> Bool
conservesNumberParticles expr = case particleDispatch @t of
  SpinTag ->
    let n = estimateNumberSites expr
        m = Expr . Sum $ G.generate n (\i -> Scaled 1 (Product (G.singleton (Generator i SpinZ))))
     in expr * m == m * expr
  SpinlessFermionTag -> undefined
  SpinfulFermionTag -> undefined

