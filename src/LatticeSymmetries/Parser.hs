module LatticeSymmetries.Parser
  ( pPrimitiveOperator,
    pOperatorString,
    PrimitiveOperator (..),
    SpinOperatorType (..),
    SpinIndex (..),
    FermionicOperatorType (..),
    sortByWithParity,
    Parity (..),
  )
where

import qualified Data.List.NonEmpty as NonEmpty
import Text.Parsec
import Prelude hiding ((<|>))

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

data SpinIndex = SpinUp | SpinDown
  deriving (Show, Eq)

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

pFermionicOperatorType :: Stream s m Char => ParsecT s u m (Char, FermionicOperatorType)
pFermionicOperatorType =
  choice
    [ (,) <$> char 'n' <*> pure FermionicNumberCountingOperator,
      (,) <$> letter <*> pType
    ]
  where
    pType =
      isJust <$> optionMaybe (char '†') >>= \x -> case x of
        True -> pure FermionicCreationOperator
        False -> pure FermionicAnnihilationOperator

pFermionicOperator :: Stream s m Char => ParsecT s u m PrimitiveOperator
pFermionicOperator = do
  (c, t) <- pFermionicOperatorType
  s <- optionMaybe pSpin
  i <- pUnicodeSubscriptNumber
  pure $ case s of
    Just s' -> SpinfulFermionicOperator t c s' i
    Nothing -> SpinlessFermionicOperator t c i

pSpinOperator :: Stream s m Char => ParsecT s u m PrimitiveOperator
pSpinOperator = do
  c <- oneOf "σS" <?> "one of σ, S"
  superscript <- oneOf "⁺⁻ˣʸᶻ" <?> "one of ⁺, ⁻, ˣ, ʸ, ᶻ"
  i <- pUnicodeSubscriptNumber
  let t = case superscript of
        '⁺' -> SpinPlusOperator
        '⁻' -> SpinMinusOperator
        'ˣ' -> SpinXOperator
        'ʸ' -> SpinYOperator
        'ᶻ' -> SpinZOperator
  pure $ SpinOperator t c i

pPrimitiveOperator :: Stream s m Char => ParsecT s u m PrimitiveOperator
pPrimitiveOperator = (pSpinOperator <|> pFermionicOperator) <?> "primitive operator"

pOperatorString :: Stream s m Char => ParsecT s u m (NonEmpty PrimitiveOperator)
pOperatorString = fromList <$> pPrimitiveOperator `sepBy1` spaces <?> "operator string"

-- isSpinful :: FermionicOperator -> Bool
-- isSpinful x = case x of
--   SpinfulFermionicOperator _ _ _ _ -> True
--   SpinlessFermionicOperator _ _ _ -> False

isSpinOperator :: PrimitiveOperator -> Bool
isSpinOperator x = case x of
  SpinOperator _ _ _ -> True
  _ -> False

isSpinfulFermionicOperator :: PrimitiveOperator -> Bool
isSpinfulFermionicOperator x = case x of
  SpinfulFermionicOperator _ _ _ _ -> True
  _ -> False

isSpinlessFermionicOperator :: PrimitiveOperator -> Bool
isSpinlessFermionicOperator x = case x of
  SpinlessFermionicOperator _ _ _ -> True
  _ -> False

isConsistent :: NonEmpty PrimitiveOperator -> Bool
isConsistent xs
  | all isSpinOperator xs = True
  | all isSpinfulFermionicOperator xs = True
  | all isSpinlessFermionicOperator xs = True
  | otherwise = False

sortIndices :: [PrimitiveOperator] -> (Parity, [PrimitiveOperator])
sortIndices xs = (p, fromList ys)
  where
    (p, ys) = sortByWithParity (comparing getSiteIndex) (toList xs)

groupIndices :: [PrimitiveOperator] -> [NonEmpty PrimitiveOperator]
groupIndices xs = NonEmpty.groupWith getSiteIndex xs

-- Left "operators are are of different type"

data Parity = Even | Odd
  deriving (Show, Eq)

sortByWithParity :: forall a. (a -> a -> Ordering) -> [a] -> (Parity, [a])
sortByWithParity cmp = go 0
  where
    toParity :: Int -> Parity
    toParity k
      | k `mod` 2 == 0 = Even
      | otherwise = Odd
    go :: Int -> [a] -> (Parity, [a])
    go !n xs =
      let (n', xs') = bubble n xs
       in if n' > n then go n' xs' else (trace (show n') (toParity n'), xs')
    bubble :: Int -> [a] -> (Int, [a])
    bubble !n (x : y : ys) =
      case cmp x y of
        GT -> let (n', ys') = bubble (n + 1) (x : ys) in (n', y : ys')
        _ -> let (n', ys') = bubble n (y : ys) in (n', x : ys')
    bubble !n xs = (n, xs)

-- sortIndices :: OperatorString -> (Parity, OperatorString)
-- sortIndices (SpinOperatorString s) =
--   let (p, s') = sortByWithParity (comparing getSiteIndex) (toList s)
--    in (p, SpinOperatorString (fromList s'))
-- sortIndices (FermionicOperatorString s) =
--   let (p, s') = sortByWithParity (comparing getSiteIndex) (toList s)
--    in (p, FermionicOperatorString (fromList s'))

normalizeIndices :: NonEmpty PrimitiveOperator -> Either Text (NonEmpty PrimitiveOperator)
normalizeIndices = undefined
