module LatticeSymmetries.Algebra where

import Control.Monad.ST
import Data.Complex
import Data.Vector (Vector)
import qualified Data.Vector
import qualified Data.Vector.Algorithms.Intro
import Data.Vector.Generic ((!), (//))
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Generic.Mutable as GM
import qualified GHC.Show as GHC
import LatticeSymmetries.Sparse (combineNeighbors)

data SpinGeneratorType = SpinPlus | SpinMinus | SpinZ

data FermionGeneratorType = FermionAnnihilate | FermionCreate
  deriving stock (Eq, Ord, Show, Generic)

data SpinGenerator = SpinGenerator {sgType :: !SpinGeneratorType, sgId :: !Int}

data FermionGenerator = FermionGenerator {fgId :: !Int, fgType :: !FermionGeneratorType}
  deriving stock (Eq, Ord, Generic)

instance GHC.Show FermionGenerator where
  show (FermionGenerator i t) = letter <> dagger <> index
    where
      letter = "c"
      dagger = if t == FermionCreate then "†" else ""
      index = toSubscript <$> (show i)
      toSubscript c = case c of
        '0' -> '₀'
        '1' -> '₁'
        '2' -> '₂'
        '3' -> '₃'
        '4' -> '₄'
        '5' -> '₅'
        '6' -> '₆'
        '7' -> '₇'
        '8' -> '₈'
        '9' -> '₉'

instance Num c => AlgebraGenerator c FermionGenerator where
  commute (FermionGenerator i t₁) (FermionGenerator j t₂)
    | t₁ == t₂ = (-1, [], 0)
    | i /= j = (-1, [], 0)
    | otherwise = (-1, [], 1)

data Term g = Term !(Vector g)
  deriving stock (Eq, Ord, Show, Generic)

tOrder :: Term g -> Int
tOrder (Term v) = G.length v

data ScaledTerm c g = ScaledTerm {stCoeff :: !c, stTerm :: !(Term g)}
  deriving stock (Eq, Show, Generic)

instance (Ord g, Ord c) => Ord (ScaledTerm c g) where
  compare (ScaledTerm c₁ t₁) (ScaledTerm c₂ t₂) =
    case compare t₁ t₂ of
      EQ -> compare c₁ c₂
      x -> x

data Expr c g = Expr !(Vector (ScaledTerm c g))
  deriving stock (Eq, Show, Generic)

scaleTerm :: Num c => c -> ScaledTerm c g -> ScaledTerm c g
scaleTerm z (ScaledTerm c g) = ScaledTerm (z * c) g

scaleExpr :: Num c => c -> Expr c g -> Expr c g
scaleExpr z (Expr v) = Expr $ G.map (scaleTerm z) v

type ℂ = Complex Double

class Ord g => AlgebraGenerator c g where
  commute :: g -> g -> (c, [(c, g)], c)

swapGenerators :: (Num c, Eq c, AlgebraGenerator c g) => Int -> Term g -> Expr c g
swapGenerators i (Term v) =
  Expr . fromList . filter ((/= 0) . stCoeff) . mconcat $
    [ [ScaledTerm c (Term $ G.concat [before, G.singleton b, G.singleton a, after])],
      [ScaledTerm z (Term $ G.concat [before, G.singleton g, after]) | (z, g) <- newTerms],
      [ScaledTerm f (Term $ G.concat [before, after])]
    ]
  where
    before = G.take i v
    a = v ! i
    b = v ! (i + 1)
    after = G.drop (i + 2) v
    (c, newTerms, f) = commute a b

toCanonical ::
  forall c g.
  (Num c, Eq c, Ord c, AlgebraGenerator c g) =>
  Term g ->
  Expr c g
toCanonical t₀ = Expr . groupTerms . sortTerms $ go (ScaledTerm (1 :: c) t₀) 0 False
  where
    go st@(ScaledTerm z t@(Term v)) i keepGoing
      | i < tOrder t - 1 =
        case compare (v ! i) (v ! (i + 1)) of
          GT ->
            let (Expr exprs) = z `scaleExpr` swapGenerators i t
             in G.concatMap (\x -> go x i True) exprs
          _ -> go st (i + 1) keepGoing
      | keepGoing = go st 0 False
      | otherwise = G.singleton st
    sortTerms v = runST $ do
      buffer <- G.thaw v
      Data.Vector.Algorithms.Intro.sort buffer
      G.unsafeFreeze buffer
    groupTerms =
      G.filter ((/= 0) . stCoeff)
        . combineNeighbors
          (\a b -> stTerm a == stTerm b)
          (\a b -> ScaledTerm (stCoeff a + stCoeff b) (stTerm a))
