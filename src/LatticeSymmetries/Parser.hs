module LatticeSymmetries.Parser
  ( pFermionicOperator,
    pFermionicString,
    FermionicOperator (..),
    FermionicOperatorType (..),
    SpinIndex (..),
  )
where

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
    <?> "index subscript (one of ₀₁₂₃₄₅₆₇₈₉)"

pUnicodeSubscriptNumber :: Stream s m Char => ParsecT s u m Int
pUnicodeSubscriptNumber = readInteger <$> many1 pSubscriptDigit
  where
    readInteger s = case readEither s of
      Right x -> x
      Left _ -> error "should not have happened"

data SpinIndex = SpinUp | SpinDown
  deriving (Show, Eq)

pSpin :: Stream s m Char => ParsecT s u m SpinIndex
pSpin = choice [char '↑' *> pure SpinUp, char '↓' *> pure SpinDown] <?> "spin label (one of ↑↓)"

data FermionicOperatorType
  = FermionicCreationOperator
  | FermionicAnnihilationOperator
  | FermionicNumberCountingOperator
  deriving (Show, Eq)

data FermionicOperator
  = SpinfulFermionicOperator !FermionicOperatorType !Char !SpinIndex !Int
  | SpinlessFermionicOperator !FermionicOperatorType !Char !Int
  deriving (Show, Eq)

data SpinOperatorType
  = SpinPlusOperator
  | SpinMinusOperator
  | SpinXOperator
  | SpinYOperator
  | SpinZOperator
  deriving (Show, Eq)

data SpinOperator = SpinOperator !SpinOperatorType !Char !Int
  deriving (Show, Eq)

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

pFermionicOperator :: Stream s m Char => ParsecT s u m FermionicOperator
pFermionicOperator = do
  (c, t) <- pFermionicOperatorType
  s <- optionMaybe pSpin
  i <- pUnicodeSubscriptNumber
  pure $ case s of
    Just s' -> SpinfulFermionicOperator t c s' i
    Nothing -> SpinlessFermionicOperator t c i

pFermionicString :: Stream s m Char => ParsecT s u m (NonEmpty FermionicOperator)
pFermionicString = fromList <$> pFermionicOperator `sepBy1` spaces

pSpinOperator :: Stream s m Char => ParsecT s u m SpinOperator
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

pSpinString :: Stream s m Char => ParsecT s u m (NonEmpty SpinOperator)
pSpinString = fromList <$> pSpinOperator `sepBy1` spaces <?> "string of spin operators"

data OperatorString
  = SpinOperatorString (NonEmpty SpinOperator)
  | FermionicOperatorString (NonEmpty FermionicOperator)
  deriving (Show, Eq)

pOperatorString :: Stream s m Char => ParsecT s u m OperatorString
pOperatorString = (SpinOperatorString <$> pSpinString) <|> (FermionicOperatorString <$> pFermionicString)
