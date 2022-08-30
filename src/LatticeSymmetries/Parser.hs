{-# LANGUAGE OverloadedLists #-}

module LatticeSymmetries.Parser
  ( pSpinlessFermionicOperator,
    pSpinfulFermionicOperator,
    pSpinOperator,
    pOperatorString,
    pBasisState,
    pNumber,
    pExpr,
    SpinIndex (..),
    -- mkSpinOperator,
    termsFromText,
    operatorFromString,
    hamiltonianFromYAML,
    ls_hs_load_hamiltonian_from_yaml,
  )
where

import Data.Aeson (FromJSON (..), withObject, (.:))
import qualified Data.Aeson.Key
import qualified Data.Aeson.KeyMap
import Data.Aeson.Types (Parser, Value)
import Data.Bits
-- import qualified Data.List.NonEmpty as NonEmpty
import Data.Ratio
import Foreign.C.String (CString)
import Foreign.Ptr (Ptr)
-- import Data.String (IsString (..))
-- import qualified Data.Vector.Generic as G

import LatticeSymmetries.Algebra
import LatticeSymmetries.Basis
import LatticeSymmetries.BitString
import LatticeSymmetries.ComplexRational
import LatticeSymmetries.FFI
import LatticeSymmetries.Generator
import LatticeSymmetries.Operator
import LatticeSymmetries.Utils
import Text.Parsec
import Text.Parsec.Language (emptyDef)
import Text.Parsec.Token
import Type.Reflection
-- import qualified Text.Read (read)
import Prelude hiding (Product, Sum, many, optional, (<|>))

-- type ℝ = Rational

-- type ℂ = ComplexRational

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
  deriving stock (Show, Eq)

data SpinOperatorType
  = SpinPlusOperator
  | SpinMinusOperator
  | SpinXOperator
  | SpinYOperator
  | SpinZOperator
  deriving stock (Show, Eq)

data PrimitiveOperator
  = SpinfulFermionicOperator !FermionicOperatorType !Char !SpinIndex !Int
  | SpinlessFermionicOperator !FermionicOperatorType !Char !Int
  | SpinOperator !SpinOperatorType !Char !Int
  deriving stock (Show, Eq)

-- getSiteIndex :: PrimitiveOperator -> Int
-- getSiteIndex (SpinfulFermionicOperator _ _ _ i) = i
-- getSiteIndex (SpinlessFermionicOperator _ _ i) = i
-- getSiteIndex (SpinOperator _ _ i) = i

class KnownIndex i where
  pIndex :: Stream s m Char => ParsecT s u m i

instance KnownIndex Int where
  pIndex = pUnicodeSubscriptNumber

instance KnownIndex (SpinIndex, Int) where
  pIndex = (,) <$> pSpin <*> pUnicodeSubscriptNumber

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

pSpinfulFermionicOperator ::
  (Num c, Stream s m Char) =>
  ParsecT s u m (Expr 'SpinfulFermionTy)
pSpinfulFermionicOperator = do
  (_, t) <- pFermionicOperatorType
  i <- pIndex
  pure $ Expr [Scaled 1 [Generator i t]]

pSpinlessFermionicOperator ::
  (Num c, Stream s m Char) =>
  ParsecT s u m (Expr 'SpinlessFermionTy)
pSpinlessFermionicOperator = do
  (_, t) <- pFermionicOperatorType
  i <- pIndex
  pure $ Expr [Scaled 1 [Generator i t]]

pSpinOperator ::
  Stream s m Char =>
  ParsecT s u m (Expr 'SpinTy)
pSpinOperator = do
  c <- oneOf "σS" <?> "one of σ, S"
  superscript <- oneOf "⁺⁻ˣʸᶻ" <?> "one of ⁺, ⁻, ˣ, ʸ, ᶻ"
  i <- pUnicodeSubscriptNumber
  let pre :: ComplexRational
      pre = if c == 'S' then fromRational (1 % 2) else 1
      t = case superscript of
        '⁺' -> [Scaled 1 [Generator i SpinPlus]]
        '⁻' -> [Scaled 1 [Generator i SpinMinus]]
        'ˣ' -> [Scaled 1 [Generator i SpinPlus], Scaled 1 [Generator i SpinMinus]]
        'ʸ' ->
          (ComplexRational 0 (-1))
            `scale` [Scaled 1 [Generator i SpinPlus], Scaled (-1) [Generator i SpinMinus]]
        'ᶻ' -> [Scaled 1 [Generator i SpinZ]]
        _ -> error "should not have happened"
  pure $ Expr $ pre `scale` t

-- pPrimitiveOperator :: Stream s m Char => ParsecT s u m PrimitiveOperator
-- pPrimitiveOperator = (pSpinOperator <|> pFermionicOperator) <?> "primitive operator"

myDef :: Stream s m Char => GenLanguageDef s u m
myDef =
  LanguageDef
    { commentStart = "",
      commentEnd = "",
      commentLine = "",
      nestedComments = True,
      identStart = letter <|> char '_',
      identLetter = alphaNum <|> oneOf "_'",
      opStart = oneOf ":!#$%&*+./<=>?@\\^|-~",
      opLetter = oneOf ":!#$%&*+./<=>?@\\^|-~",
      reservedOpNames = [],
      reservedNames = [],
      caseSensitive = True
    }

skipSpaces :: Stream s m Char => ParsecT s u m ()
skipSpaces = do
  _ <- many space
  pure ()

-- skipSpaces1 :: Stream s m Char => ParsecT s u m ()
-- skipSpaces1 = do
--   _ <- Parsec.many1 Parsec.space
--   pure ()

pNumber :: Stream s m Char => ParsecT s u m ComplexRational
pNumber = (try pComplex <|> (ComplexRational <$> pReal <*> pure 0)) <?> "real or complex number"
  where
    pReal = (<?> "real number") $ do
      r <- naturalOrFloat (Text.Parsec.Token.makeTokenParser myDef)
      pure $ either toRational toRational r
    pImag = (<?> "imaginary part") $ do
      c <- char '+' <|> char '-'
      skipSpaces
      i <- pReal
      (char 'i' >> char 'm') <|> (char 'I') <|> (char 'ⅈ')
      if c == '+'
        then pure i
        else pure (-i)
    pComplex = do
      char '(' >> skipSpaces
      r <- pReal
      i <- option 0 (try pImag)
      skipSpaces >> char ')' >> skipSpaces
      pure $ ComplexRational r i

pOperatorString ::
  (Stream s m Char, Ord (IndexType t), Algebra (GeneratorType t)) =>
  ParsecT s u m (Expr t) ->
  ParsecT s u m (Expr t)
pOperatorString pPrimitive = do
  c <- option 1 $ do z <- pNumber; optional (char '×' >> skipSpaces); pure z
  ops <- pPrimitive `sepBy1` spaces <?> "operator string"
  let expr = case (unExpr <$> ops) of
        (o : os) -> scale c . simplifyExpr . Expr $ foldr (*) o os
        [] -> error "should never happen"
  pure expr

pExpr ::
  forall s u m t.
  (Stream s m Char, Ord (IndexType t), Algebra (GeneratorType t)) =>
  ParsecT s u m (Expr t) ->
  ParsecT s u m (Expr t)
pExpr pPrimitive = do
  t <- pTerm
  skipSpaces
  ts <- many $ do
    t' <- pPlusTerm <|> pMinusTerm
    skipSpaces
    pure t'
  pure $ foldl' (+) t ts
  where
    pProduct p = do
      terms <- many1 $ do t <- p; skipSpaces; pure t
      case terms of
        (t : ts) -> pure $ foldl' (*) t ts
        [] -> error "should never happen"
    pScaledExpr = do
      z <- pNumber
      optional (char '×' >> skipSpaces)
      expr <- pParenExpr <|> pPrimitive
      pure $ scale z expr
    pParenExpr =
      between (char '(' >> skipSpaces) (skipSpaces >> char ')' >> skipSpaces) $
        pExpr pPrimitive
    -- pNormalExpr = pProduct pPrimitive
    pTerm = pProduct (try pScaledExpr <|> pParenExpr <|> pPrimitive)
    pPlusTerm = char '+' >> skipSpaces >> pTerm
    pMinusTerm = char '-' >> skipSpaces >> (fmap negate pTerm)

-- mkSpinOperator ::
--   HasCallStack =>
--   Text ->
--   [[Int]] ->
--   Polynomial ℂ (Generator Int SpinGeneratorType)
-- mkSpinOperator s indices = case parse (pOperatorString pSpinOperator) "" s of
--   Left e -> error (show e)
--   Right x -> simplify $ forIndices x indices

-- getPrimitiveParser :: Stream s m Char => Basis t -> ParsecT s u m (Sum (Scaled ComplexRational (Factor t)))
-- getPrimitiveParser b = case basisHeader b of
--   SpinHeader _ _ _ _ -> pSpinOperator
--   SpinfulFermionHeader _ _ -> pFermionicOperator
--   SpinlessFermionHeader _ _ -> pFermionicOperator

getPrimitiveParser ::
  forall (t :: ParticleTy) s m proxy u.
  (Stream s m Char, Typeable t) =>
  proxy t ->
  ParsecT s u m (Expr t)
getPrimitiveParser _
  | Just HRefl <- eqTypeRep (typeRep @t) (typeRep @'SpinTy) = pSpinOperator
  | Just HRefl <- eqTypeRep (typeRep @t) (typeRep @'SpinlessFermionTy) = pSpinlessFermionicOperator
  | Just HRefl <- eqTypeRep (typeRep @t) (typeRep @'SpinfulFermionTy) = pSpinfulFermionicOperator
  | otherwise = error "this should never happen by construction"

termsFromText ::
  forall t.
  IsBasis t =>
  Text ->
  [[Int]] ->
  Polynomial ComplexRational (Factor t)
termsFromText s indices = case parse (pOperatorString (getPrimitiveParser (Proxy @t))) "" s of
  Left e -> error $ "failed to parse " <> show s <> ": " <> show e
  Right x -> forSiteIndices (unExpr x) indices

operatorFromString ::
  forall t.
  IsBasis t =>
  Basis t ->
  Text ->
  [[Int]] ->
  Operator t
operatorFromString basis s indices =
  operatorFromHeader . OperatorHeader basis $
    Expr $
      termsFromText @t s indices

pBasisState :: Stream s m Char => ParsecT s u m (BasisState t)
pBasisState = do
  _ <- char '|'
  s <- many1 (char '0' <|> char '1')
  _ <- char '⟩'
  let go !n !size [] = BasisState size (BitString n)
      go !n !size (c : cs) = go ((n `shiftL` 1) .|. x) (size + 1) cs
        where
          x = case c of
            '0' -> 0
            '1' -> 1
            _ -> error "should never happen"
  pure $ go 0 0 s

instance IsString (BasisState t) where
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

termsFromJSON :: forall t. IsBasis t => Value -> Parser (Polynomial ComplexRational (Factor t))
termsFromJSON = withObject "Term" $ \v -> do
  expr <- v .: "expression"
  sites <- v .: "sites"
  pure $ termsFromText @t expr sites

data OperatorTermSpec = OperatorTermSpec !Text ![[Int]]

instance FromJSON OperatorTermSpec where
  parseJSON = withObject "Term" $ \v ->
    OperatorTermSpec <$> v .: "expression" <*> v .: "sites"

operatorFromJSON :: IsBasis t => Basis t -> Value -> Parser (Operator t)
operatorFromJSON basis = withObject "Operator" $ \v -> do
  (terms :: NonEmpty OperatorTermSpec) <- v .: "terms"
  let (o :| os) = fmap (\(OperatorTermSpec expr sites) -> operatorFromString basis expr sites) terms
  pure $ foldl' (+) o os

newtype BasisAndHamiltonianConfig = BasisAndHamiltonianConfig SomeOperator

instance FromJSON BasisAndHamiltonianConfig where
  parseJSON = withObject "Config" $ \v -> do
    basis <- v .: "basis"
    withSomeBasis basis $ \basis' ->
      case Data.Aeson.KeyMap.lookup (Data.Aeson.Key.fromString "hamiltonian") v of
        Just h -> BasisAndHamiltonianConfig . SomeOperator <$> operatorFromJSON basis' h
        Nothing -> fail "missing 'hamiltonian' key"

hamiltonianFromYAML :: HasCallStack => Text -> IO SomeOperator
hamiltonianFromYAML path = (\(BasisAndHamiltonianConfig x) -> x) <$> objectFromYAML "Hamiltonian" path

ls_hs_load_hamiltonian_from_yaml :: CString -> IO (Ptr Coperator)
ls_hs_load_hamiltonian_from_yaml cFilename =
  foldSomeOperator borrowCoperator =<< hamiltonianFromYAML =<< peekUtf8 cFilename
