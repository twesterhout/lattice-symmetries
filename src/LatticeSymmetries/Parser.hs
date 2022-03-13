{-# LANGUAGE OverloadedLists #-}

module LatticeSymmetries.Parser
  ( pFermionicOperator,
    pSpinOperator,
    pOperatorString,
    SpinIndex (..),
    mkSpinOperator,
    operatorFromString,
  )
where

import Data.Bits
import qualified Data.List.NonEmpty as NonEmpty
import Data.Ratio
import Data.String (IsString (..))
import qualified Data.Vector.Generic as G
import LatticeSymmetries.Algebra
import LatticeSymmetries.Basis
import LatticeSymmetries.ComplexRational
import LatticeSymmetries.Operator
import Text.Parsec
import qualified Text.Read (read)
import Prelude hiding (Product, Sum, (<|>))

type ℝ = Rational

type ℂ = ComplexRational

pSubscriptDigit :: Stream s m Char => ParsecT s u m Char
pSubscriptDigit =
  choice
    [ char '₀' *> pure '0',
      char '₁' *> pure '1',
      char '₂' *> pure '2',
      char '₃' *> pure '3',
      char '₄' *> pure '4',
      char '₅' *> pure '5',
      char '₆' *> pure '6',
      char '₇' *> pure '7',
      char '₈' *> pure '8',
      char '₉' *> pure '9'
    ]
    <?> "index subscript (one of ₀, ₁, ₂, ₃, ₄, ₅, ₆, ₇, ₈, ₉)"

pUnicodeSubscriptNumber :: Stream s m Char => ParsecT s u m Int
pUnicodeSubscriptNumber = readInteger <$> many1 pSubscriptDigit
  where
    readInteger s = case readEither s of
      Right x -> x
      Left _ -> error "should not have happened"

pSpin :: Stream s m Char => ParsecT s u m SpinIndex
pSpin = choice [char '↑' *> pure SpinUp, char '↓' *> pure SpinDown] <?> "spin label (one of ↑, ↓)"

data FermionicOperatorType
  = FermionicCreationOperator
  | FermionicAnnihilationOperator
  | FermionicNumberCountingOperator
  deriving (Show, Eq)

data SpinOperatorType
  = SpinPlusOperator
  | SpinMinusOperator
  | SpinXOperator
  | SpinYOperator
  | SpinZOperator
  deriving (Show, Eq)

data PrimitiveOperator
  = SpinfulFermionicOperator !FermionicOperatorType !Char !SpinIndex !Int
  | SpinlessFermionicOperator !FermionicOperatorType !Char !Int
  | SpinOperator !SpinOperatorType !Char !Int
  deriving (Show, Eq)

getSiteIndex :: PrimitiveOperator -> Int
getSiteIndex (SpinfulFermionicOperator _ _ _ i) = i
getSiteIndex (SpinlessFermionicOperator _ _ i) = i
getSiteIndex (SpinOperator _ _ i) = i

class KnownIndex i where
  pIndex :: Stream s m Char => ParsecT s u m i

instance KnownIndex Int where
  pIndex = pUnicodeSubscriptNumber

instance KnownIndex (SpinIndex, Int) where
  pIndex = flip (,) <$> pUnicodeSubscriptNumber <*> pSpin

pFermionicOperatorType :: Stream s m Char => ParsecT s u m (Char, FermionGeneratorType)
pFermionicOperatorType =
  choice
    [ (,) <$> char 'n' <*> pure FermionCount,
      (,) <$> char 'c' <*> pType
    ]
  where
    pType =
      isJust <$> optionMaybe (char '†') >>= \x -> case x of
        True -> pure FermionCreate
        False -> pure FermionAnnihilate

pFermionicOperator ::
  (KnownIndex i, Num c, Stream s m Char) =>
  ParsecT s u m (Sum (Scaled c (Generator i FermionGeneratorType)))
pFermionicOperator = do
  (_, t) <- pFermionicOperatorType
  i <- pIndex
  pure $ [Scaled 1 (Generator i t)]

pSpinOperator ::
  Stream s m Char =>
  ParsecT s u m (Sum (Scaled ComplexRational (Generator Int SpinGeneratorType)))
pSpinOperator = do
  c <- oneOf "σS" <?> "one of σ, S"
  superscript <- oneOf "⁺⁻ˣʸᶻ" <?> "one of ⁺, ⁻, ˣ, ʸ, ᶻ"
  i <- pUnicodeSubscriptNumber
  let pre :: ComplexRational
      pre = if c == 'S' then fromRational (1 % 2) else 1
      t = case superscript of
        '⁺' -> [Scaled 1 (Generator i SpinPlus)]
        '⁻' -> [Scaled 1 (Generator i SpinMinus)]
        'ˣ' -> [Scaled 1 (Generator i SpinPlus), Scaled 1 (Generator i SpinMinus)]
        'ʸ' ->
          (ComplexRational 0 (-1))
            `scale` [Scaled 1 (Generator i SpinPlus), Scaled (-1) (Generator i SpinMinus)]
        'ᶻ' -> [Scaled 1 (Generator i SpinZ)]
  pure $ pre `scale` t

-- pPrimitiveOperator :: Stream s m Char => ParsecT s u m PrimitiveOperator
-- pPrimitiveOperator = (pSpinOperator <|> pFermionicOperator) <?> "primitive operator"

pOperatorString ::
  (Stream s m Char, Num c) =>
  ParsecT s u m (Sum (Scaled c (Generator i g))) ->
  ParsecT s u m (Sum (Scaled c (Product (Generator i g))))
pOperatorString pPrimitive =
  expandProduct . fromList
    <$> (pPrimitive `sepBy1` spaces <?> "operator string")

mkSpinOperator ::
  HasCallStack =>
  Text ->
  [[Int]] ->
  Polynomial ℂ (Generator Int SpinGeneratorType)
mkSpinOperator s indices = case parse (pOperatorString pSpinOperator) "" s of
  Left e -> error (show e)
  Right x -> simplify $ forIndices x indices

class IsBasis basis => KnownBasis basis where
  getPrimitiveParser ::
    Stream s m Char =>
    basis ->
    ParsecT s u m (Sum (Scaled ComplexRational (Generator (IndexType basis) (GeneratorType basis))))

instance KnownBasis SpinBasis where
  getPrimitiveParser _ = pSpinOperator

instance KnownBasis SpinfulFermionicBasis where
  getPrimitiveParser _ = pFermionicOperator

instance KnownBasis SpinlessFermionicBasis where
  getPrimitiveParser _ = pFermionicOperator

operatorFromString ::
  ( KnownBasis basis,
    Ord (IndexType basis),
    Bounded (GeneratorType basis),
    Enum (GeneratorType basis),
    HasMatrixRepresentation (GeneratorType basis),
    Algebra (GeneratorType basis)
  ) =>
  basis ->
  Text ->
  [[IndexType basis]] ->
  Operator ComplexRational basis
operatorFromString basis s indices = case parse (pOperatorString (getPrimitiveParser basis)) "" s of
  Left e -> error (show e)
  Right x ->
    let terms = simplify $ forIndices x indices
     in Operator basis terms

pBasisState :: Stream s m Char => ParsecT s u m BasisState
pBasisState = do
  _ <- char '|'
  s <- many1 (char '0' <|> char '1')
  _ <- char '⟩'
  let go !n !count [] = BasisState count (BitString n)
      go !n !count (c : cs) = go ((n `shiftL` 1) .|. x) (count + 1) cs
        where
          x = case c of
            '0' -> 0
            '1' -> 1
            _ -> error "should never happen"
  pure $ go 0 0 s

instance IsString BasisState where
  fromString s = case parse pBasisState "" s of
    Left e -> error (show e)
    Right x -> x

-- isSpinful :: FermionicOperator -> Bool
-- isSpinful x = case x of
--   SpinfulFermionicOperator _ _ _ _ -> True
--   SpinlessFermionicOperator _ _ _ -> False

-- isSpinOperator :: PrimitiveOperator -> Bool
-- isSpinOperator x = case x of
--   SpinOperator _ _ _ -> True
--   _ -> False

-- isSpinfulFermionicOperator :: PrimitiveOperator -> Bool
-- isSpinfulFermionicOperator x = case x of
--   SpinfulFermionicOperator _ _ _ _ -> True
--   _ -> False

-- isSpinlessFermionicOperator :: PrimitiveOperator -> Bool
-- isSpinlessFermionicOperator x = case x of
--   SpinlessFermionicOperator _ _ _ -> True
--   _ -> False

-- isConsistent :: NonEmpty PrimitiveOperator -> Bool
-- isConsistent xs
--   | all isSpinOperator xs = True
--   | all isSpinfulFermionicOperator xs = True
--   | all isSpinlessFermionicOperator xs = True
--   | otherwise = False

-- sortIndices :: [PrimitiveOperator] -> (Parity, [PrimitiveOperator])
-- sortIndices xs = (p, fromList ys)
--   where
--     (p, ys) = sortByWithParity (comparing getSiteIndex) (toList xs)

-- groupIndices :: [PrimitiveOperator] -> [NonEmpty PrimitiveOperator]
-- groupIndices xs = NonEmpty.groupWith getSiteIndex xs

-- Left "operators are are of different type"

-- data Parity = Even | Odd
--   deriving (Show, Eq)

-- sortByWithParity :: forall a. (a -> a -> Ordering) -> [a] -> (Parity, [a])
-- sortByWithParity cmp = go 0
--   where
--     toParity :: Int -> Parity
--     toParity k
--       | k `mod` 2 == 0 = Even
--       | otherwise = Odd
--     go :: Int -> [a] -> (Parity, [a])
--     go !n xs =
--       let (n', xs') = bubble n xs
--        in if n' > n then go n' xs' else (trace (show n') (toParity n'), xs')
--     bubble :: Int -> [a] -> (Int, [a])
--     bubble !n (x : y : ys) =
--       case cmp x y of
--         GT -> let (n', ys') = bubble (n + 1) (x : ys) in (n', y : ys')
--         _ -> let (n', ys') = bubble n (y : ys) in (n', x : ys')
--     bubble !n xs = (n, xs)

-- sortIndices :: OperatorString -> (Parity, OperatorString)
-- sortIndices (SpinOperatorString s) =
--   let (p, s') = sortByWithParity (comparing getSiteIndex) (toList s)
--    in (p, SpinOperatorString (fromList s'))
-- sortIndices (FermionicOperatorString s) =
--   let (p, s') = sortByWithParity (comparing getSiteIndex) (toList s)
--    in (p, FermionicOperatorString (fromList s'))

-- normalizeIndices :: NonEmpty PrimitiveOperator -> Either Text (NonEmpty PrimitiveOperator)
-- normalizeIndices = undefined
