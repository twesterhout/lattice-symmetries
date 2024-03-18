{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE DuplicateRecordFields #-}
{-# LANGUAGE OverloadedStrings #-}

module LatticeSymmetries.Parser
  ( parseExprEither
  , parsePhaseEither
  , SExpr (..)
  , SOp (..)
  , SSpin (..)
  , SFermion (..)
  , Parser
  -- * Exports for testing
  , pReal
  , pCoeff
  )
where

import Control.Monad.Combinators
import Control.Monad.Combinators.NonEmpty qualified as NonEmpty
import Data.Ratio
import Data.Text qualified as Text
import LatticeSymmetries.ComplexRational
import LatticeSymmetries.Generator
import Text.Megaparsec
import Text.Megaparsec.Char
import Text.Megaparsec.Char.Lexer (decimal, float, lexeme, signed)
import Prelude hiding (Product, Sum, many, optional, some, (<|>))

data SFermion
  = SFermionCreate
  | SFermionAnnihilate
  | SFermionNumber
  deriving stock (Show, Eq)

data SSpin
  = SSpinPlus
  | SSpinMinus
  | SSpinX
  | SSpinY
  | SSpinZ
  deriving stock (Show, Eq)

data SOp
  = SFermionOp !SFermion !(Maybe SpinIndex) !Int
  | SSpinOp !Char !SSpin !Int
  | SIdentity
  deriving stock (Show, Eq)

data SExpr
  = SSum ![SExpr]
  | SScaled !ℂ !SExpr
  | SProduct !(NonEmpty SExpr)
  | SPrimitive !SOp
  deriving stock (Show, Eq)

type Parser = Parsec Void Text

-- | Parse the phase
parsePhaseEither :: forall a. Integral a => Text -> Either Text (Ratio a)
parsePhaseEither s = case parse p "" s of
  Left e -> Left $ "Error parsing '" <> s <> "': " <> Text.pack (errorBundlePretty e)
  Right x -> Right x
  where
    p :: Parser (Ratio a)
    p = (%) <$> decimal <*> (char '/' *> decimal)

parseExprEither :: Text -> Either Text SExpr
parseExprEither s = case parse pSum "" s of
  Left e -> Left $ "Error parsing '" <> s <> "': " <> Text.pack (errorBundlePretty e)
  Right x -> Right x

pSum :: Parser SExpr
pSum = fmap SSum . many $ do
  sign <- fromMaybe '+' <$> optional (lexeme space (oneOf ['-', '+']))
  (SScaled c expr) <- lexeme space pScaled
  case sign of
    '-' -> pure $ SScaled (-c) expr
    '+' -> pure $ SScaled c expr
    _ -> fail "can never happen"

bracketted :: Parser a -> Parser a
bracketted = between (char '(') (char ')')

pScaled :: Parser SExpr
pScaled = do
  c <- fromMaybe 1 <$> optional (lexeme space (try pCoeff))
  p <- lexeme space pProduct
  pure $ SScaled c p

pProduct :: Parser SExpr
pProduct = fmap SProduct . NonEmpty.some $ do
  choice
    [ try (lexeme space (bracketted pSum))
    , lexeme space pPrimitive
    ]


pReal :: Parser Rational
pReal = try (realToFrac @Double <$> float) <|> try ((%) <$> decimal <*> (char '/' *> decimal)) <|> decimal

pCoeff :: Parser ℂ
pCoeff = try ((0 `ComplexRational`) <$> signed space pImaginary) <|> try ((`ComplexRational` 0) <$> signed space pReal) <|> bothParts
  where
    pImaginary = (pImaginaryUnit $> 1) <|> (pReal <* pImaginaryUnit)
    pImaginaryUnit = char 'j' <|> char 'ⅈ' <|> (char 'i' *> char 'm')
    pPlusOrMinus = (char '+' <|> char '-') >>= \case
      '+' -> pure id
      '-' -> pure negate
      _ -> error "cannot happen"
    bothParts = between (lexeme space (char '(')) (char ')') $ do
      r <- lexeme space (signed space pReal)
      f <- lexeme space pPlusOrMinus
      i <- lexeme space pImaginary
      pure $ ComplexRational r (f i)

pPrimitive :: Parser SExpr
pPrimitive = SPrimitive <$> (pPrimitiveIdentity <|> pPrimitiveSpin <|> pPrimitiveFermion)

pPrimitiveIdentity :: Parser SOp
pPrimitiveIdentity = char 'I' $> SIdentity

pPrimitiveSpin :: Parser SOp
pPrimitiveSpin = SSpinOp <$> prefix <*> superscript <*> subscript
  where
    prefix =
      choice
        [ char 'σ'
        , char 'S'
        , try (string "\\sigma" $> 'σ')
        ]
    superscript =
      choice
        [ char 'ˣ' $> SSpinX
        , char 'ʸ' $> SSpinY
        , char 'ᶻ' $> SSpinZ
        , char '⁺' $> SSpinPlus
        , char '⁻' $> SSpinMinus
        , optional (char '^')
            >> choice
              [ char 'x' $> SSpinX
              , char 'y' $> SSpinY
              , char 'z' $> SSpinZ
              , char '+' $> SSpinPlus
              , char '-' $> SSpinMinus
              ]
        ]
    subscript = pUnicodeSubscriptNumber <|> (optional (char '_') >> decimal)

pPrimitiveFermion :: Parser SOp
pPrimitiveFermion = do
  op <-
    choice
      [ try (string "c†") $> SFermionCreate
      , try (char 'c' >> optional (char '^') >> string "\\dagger") $> SFermionCreate
      , char 'c' $> SFermionAnnihilate
      , char 'n' $> SFermionNumber
      ]
  (i, s) <-
    choice
      [ (,) <$> pUnicodeSubscriptNumber <*> pSpin
      , (,) <$> (optional (char '_') >> decimal) <*> pSpin
      ]
  pure $ SFermionOp op s i
  where
    pSpin =
      optional $
        choice
          [ char '↑' $> SpinUp
          , char '↓' $> SpinDown
          , try (optional (char '_') >> string "\\up" $> SpinUp)
          , try (optional (char '_') >> string "\\down" $> SpinDown)
          ]

pSubscriptDigit :: Parser Char
pSubscriptDigit =
  choice
    [ char '₀' $> '0'
    , char '₁' $> '1'
    , char '₂' $> '2'
    , char '₃' $> '3'
    , char '₄' $> '4'
    , char '₅' $> '5'
    , char '₆' $> '6'
    , char '₇' $> '7'
    , char '₈' $> '8'
    , char '₉' $> '9'
    ]
    <?> "index subscript (one of ₀, ₁, ₂, ₃, ₄, ₅, ₆, ₇, ₈, ₉)"

pUnicodeSubscriptNumber :: Parser Int
pUnicodeSubscriptNumber = readInteger <$> some pSubscriptDigit
  where
    readInteger s = case readEither s of
      Right x -> x
      Left _ -> error "should not have happened"
