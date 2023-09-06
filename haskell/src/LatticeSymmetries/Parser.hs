{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE DuplicateRecordFields #-}
{-# LANGUAGE OverloadedStrings #-}

module LatticeSymmetries.Parser
  ( parseExprEither
  , SExpr (..)
  , SOp (..)
  , SSpin (..)
  , SFermion (..)
  )
where

import Control.Monad.Combinators
import Control.Monad.Combinators.NonEmpty qualified as NonEmpty
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
  deriving stock (Show, Eq)

data SExpr
  = SSum ![SExpr]
  | SScaled !ℂ !SExpr
  | SProduct !(NonEmpty SExpr)
  | SPrimitive !SOp
  deriving stock (Show, Eq)

type Parser = Parsec Void Text

parseExprEither :: Text -> Either Text SExpr
parseExprEither s = case parse pSum "" s of
  Left e -> Left $ Text.pack (errorBundlePretty e)
  Right x -> Right x

pSum :: Parser SExpr
pSum = fmap SSum . many $ do
  sign <- fromMaybe '+' <$> optional (lexeme space (oneOf ['-', '+']))
  (SScaled c expr) <- lexeme space pScaled
  case sign of
    '-' -> pure $ SScaled (-c) expr
    '+' -> pure $ SScaled c expr
    _ -> fail "can never happen"

pScaled :: Parser SExpr
pScaled = do
  c <- fromMaybe 1 <$> optional (lexeme space (try pCoeff))
  p <- lexeme space pProduct
  pure $ SScaled c p

pProduct :: Parser SExpr
pProduct = fmap SProduct . NonEmpty.some $ do
  choice
    [ try (lexeme space (between (char '(') (char ')') pSum))
    , lexeme space pPrimitive
    ]

pCoeff :: Parser ℂ
pCoeff =
  choice
    [ between (lexeme space (char '(')) (char ')') $
        -- lexeme space $
        choice
          [ try bothP
          , try ((0 `ComplexRational`) <$> imagP)
          , (`ComplexRational` 0) <$> realP
          ]
    , try $ (0 `ComplexRational`) <$> imagP
    , (`ComplexRational` 0) <$> realP
    ]
  where
    imaginaryUnit =
      choice
        [ string "j"
        , string "ⅈ"
        , string "im"
        ]
    realP = realToFrac @Double <$> (try float <|> decimal)
    imagP = do
      c <- fromMaybe 1 <$> optional (try realP)
      _ <- imaginaryUnit
      pure c
    bothP = do
      !r <- lexeme space (signed space realP)
      !i <- lexeme space (signed space imagP)
      pure (ComplexRational r i)

pPrimitive :: Parser SExpr
pPrimitive = SPrimitive <$> (pPrimitiveSpin <|> pPrimitiveFermion)

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

-- pSum :: Parser SyntaxTree
-- pSum = do
--   x <- lexeme space pProduct
--   ys <- many $ do
--     sign <- lexeme space (oneOf ['-', '+'])
--     p@(SProduct terms) <- lexeme space pProduct
--     pure $ case sign of
--       '-' -> SProduct $ NonEmpty.cons (SNumber (-1)) terms
--       '+' -> p
--       _ -> error "can never happen"
--   pure $ SSum (x :| ys)

{-
pSubscriptDigit :: Stream s m Char => ParsecT s u m Char
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

pUnicodeSubscriptNumber :: Stream s m Char => ParsecT s u m Int
pUnicodeSubscriptNumber = readInteger <$> many1 pSubscriptDigit
  where
    readInteger s = case readEither s of
      Right x -> x
      Left _ -> error "should not have happened"

pSpin :: Stream s m Char => ParsecT s u m SpinIndex
pSpin = choice [char '↑' *> pure SpinUp, char '↓' *> pure SpinDown] <?> "spin label (one of ↑, ↓)"

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
    [ (,) <$> char 'n' <*> pure FermionCount
    , (,) <$> char 'c' <*> pType
    ]
  where
    pType =
      isJust <$> optionMaybe (char '†') >>= \x -> case x of
        True -> pure FermionCreate
        False -> pure FermionAnnihilate

pSpinfulFermionicOperator
  :: Stream s m Char
  => ParsecT s u m (Expr 'SpinfulFermionTy)
pSpinfulFermionicOperator = do
  (_, t) <- pFermionicOperatorType
  i <- pIndex
  pure $ Expr [Scaled 1 [Generator i t]]

pSpinlessFermionicOperator
  :: Stream s m Char
  => ParsecT s u m (Expr 'SpinlessFermionTy)
pSpinlessFermionicOperator = do
  (_, t) <- pFermionicOperatorType
  i <- pIndex
  pure $ Expr [Scaled 1 [Generator i t]]

pSpinOperator
  :: Stream s m Char
  => ParsecT s u m (Expr 'SpinTy)
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
    { commentStart = ""
    , commentEnd = ""
    , commentLine = ""
    , nestedComments = True
    , identStart = letter <|> char '_'
    , identLetter = alphaNum <|> oneOf "_'"
    , opStart = oneOf ":!#$%&*+./<=>?@\\^|-~"
    , opLetter = oneOf ":!#$%&*+./<=>?@\\^|-~"
    , reservedOpNames = []
    , reservedNames = []
    , caseSensitive = True
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
      _ <- (char 'i' >> char 'm') <|> (char 'I') <|> (char 'ⅈ')
      if c == '+'
        then pure i
        else pure (-i)
    pComplex = do
      char '(' >> skipSpaces
      r <- pReal
      i <- option 0 (try pImag)
      skipSpaces >> char ')' >> skipSpaces
      pure $ ComplexRational r i

pOperatorString
  :: (Stream s m Char, Ord (IndexType t), Algebra (GeneratorType t))
  => ParsecT s u m (Expr t)
  -> ParsecT s u m (Expr t)
pOperatorString pPrimitive = do
  c <- option 1 $ do z <- pNumber; optional (char '×' >> skipSpaces); pure z
  ops <- pPrimitive `sepBy1` spaces <?> "operator string"
  let expr = case (unExpr <$> ops) of
        (o : os) -> scale c . simplifyExpr . Expr $ foldr (*) o os
        [] -> error "should never happen"
  pure expr

pExpr
  :: forall s u m t
   . (Stream s m Char, Ord (IndexType t), Algebra (GeneratorType t))
  => ParsecT s u m (Expr t)
  -> ParsecT s u m (Expr t)
pExpr pPrimitive = do
  isMinus <- isJust <$> optionMaybe (char '-' >> skipSpaces)
  t <- pTerm
  skipSpaces
  ts <- many $ do
    t' <- pPlusTerm <|> pMinusTerm
    skipSpaces
    pure t'
  pure $
    (if isMinus then negate else id) $
      foldl' (+) t ts
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

instance IsString (Expr 'SpinTy) where
  fromString = mkExpr SpinTag . toText

instance IsString (Expr 'SpinlessFermionTy) where
  fromString = mkExpr SpinlessFermionTag . toText

instance IsString (Expr 'SpinfulFermionTy) where
  fromString = mkExpr SpinfulFermionTag . toText

mkExpr' :: ParticleTag t -> Text -> Either Text (Expr t)
mkExpr' t s = case run t of
  Left e -> Left (show e)
  Right x -> Right x
  where
    run :: ParticleTag t -> Either ParseError (Expr t)
    run SpinTag = parse (pExpr pSpinOperator) "" s
    run SpinlessFermionTag = parse (pExpr pSpinlessFermionicOperator) "" s
    run SpinfulFermionTag = parse (pExpr pSpinfulFermionicOperator) "" s

mkExpr :: HasCallStack => ParticleTag t -> Text -> Expr t
mkExpr tag s = either error id $ mkExpr' tag s

mkSomeExpr' :: Maybe ParticleTy -> Text -> Either Text SomeExpr
mkSomeExpr' tp s
  | Just SpinTy <- tp = make SpinTag
  | Just SpinlessFermionTy <- tp = make SpinlessFermionTag
  | Just SpinfulFermionTy <- tp = make SpinfulFermionTag
  | isJust $ Text.find (\c -> c == 'σ' || c == 'S') s = make SpinTag
  | isJust $ Text.find (\c -> c == '↑' || c == '↓') s = make SpinfulFermionTag
  | otherwise = make SpinlessFermionTag
  where
    make :: IsBasis t => ParticleTag t -> Either Text SomeExpr
    make tag = fmap (SomeExpr tag) $ mkExpr' tag s

mkSomeExpr :: HasCallStack => Maybe ParticleTy -> Text -> SomeExpr
mkSomeExpr tp s = either error id $ mkSomeExpr' tp s

-- instance FromJSON SomeExpr where
--   parseJSON = withObject "SomeExpr" $ \v -> do
--     tp <- v .:? "particle"
--     s <- v .: "expression"
--     sites <- v .:? "sites"
--     let expand x = case sites of
--           Just is -> Expr $ forSiteIndices (unExpr x) is
--           Nothing -> x
--     case fmap expand of
--       Left e -> parserThrowError [Key "expression"] (show e)
--       Right x -> pure x

-- (pure . SomeExpr tag) $ mkExprSafe tag s
-- let make :: ParticleTag t -> Parser SomeExpr
--     make tag = either (parserThrowError [Key "expression"] . show) (pure . SomeExpr tag) $ mkExprSafe tag s
--     dispatch
--       | Just SpinTy <- tp = make SpinTag
--       | Just SpinlessFermionTy <- tp = make SpinlessFermionTag
--       | Just SpinfulFermionTy <- tp = make SpinfulFermionTag
--       | isJust $ Text.find (\c -> c == 'σ' || c == 'S') s = make SpinTag
--       | isJust $ Text.find (\c -> c == '↑' || c == '↓') s = make SpinfulFermionTag
--       | otherwise = make SpinlessFermionTag
-- dispatch

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

{-
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
-}

{-
operatorFromString ::
  forall t.
  IsBasis t =>
  Basis t ->
  Text ->
  Maybe [[Int]] ->
  Operator t
operatorFromString basis s maybeIndices =
  operatorFromHeader . OperatorHeader basis $ expand expr
  where
    expr = mkExpr (getParticleTag (basisHeader basis)) s
    expand x = case maybeIndices of
      Just indices -> Expr $ forSiteIndices (unExpr x) indices
      Nothing -> x
-}

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
-- termsFromJSON :: forall t. IsBasis t => Value -> Parser (Polynomial ComplexRational (Factor t))
-- termsFromJSON = withObject "Term" $ \v -> do
--   expr <- v .: "expression"
--   sites <- v .: "sites"
--   pure $ termsFromText @t expr sites

objectFromYAML :: (HasCallStack, FromJSON a) => Text -> Text -> IO a
objectFromYAML _ filename = do
  r <- decodeFileWithWarnings (toString filename)
  case r of
    Left e -> error $ toText $ prettyPrintParseException e
    Right (_, x) -> pure x

instance FromJSON ConfigSpec where
  parseJSON = withObject "Config" $ \v -> do
    someBasis <- v .: "basis"
    maybeHamiltonian <-
      withSomeBasis someBasis $ \basis ->
        case Data.Aeson.KeyMap.lookup (Data.Aeson.Key.fromString "hamiltonian") v of
          Just h ->
            Just . SomeOperator (getParticleTag basis)
              <$> operatorFromJSON basis h
          Nothing -> pure Nothing
    observables <-
      withSomeBasis someBasis $ \basis ->
        case Data.Aeson.KeyMap.lookup (Data.Aeson.Key.fromString "observables") v of
          Just h ->
            flip (withArray "Observables") h $
              G.mapM (fmap (SomeOperator (getParticleTag basis)) . operatorFromJSON basis)
          Nothing -> pure G.empty
    pure $ ConfigSpec someBasis maybeHamiltonian observables

configFromYAML :: HasCallStack => Text -> IO ConfigSpec
configFromYAML path = objectFromYAML "config" path

-- basisFromYAML :: HasCallStack => Text -> IO SomeBasis
-- basisFromYAML = fmap (\(ConfigSpec basis _ _) -> basis) . configFromYAML

-- hamiltonianFromYAML :: HasCallStack => Text -> IO SomeOperator
-- hamiltonianFromYAML path = do
--   (ConfigSpec basis maybeHamiltonian _) <- configFromYAML path
--   case maybeHamiltonian of
--     Just h -> pure h
--     Nothing -> error $ "missing 'hamiltonian' key in '" <> path <> "'"

-- observablesFromYAML :: HasCallStack => Text -> IO (Vector SomeOperator)
-- observablesFromYAML = fmap (\(ConfigSpec _ _ observables) -> observables) . configFromYAML

-- basisFromYAML :: HasCallStack => Text -> IO SomeBasis
-- basisFromYAML path = (\(BasisOnlyConfig x) -> x) <$> objectFromYAML "Basis" path

-- data OperatorTermSpec = OperatorTermSpec !Text !(Maybe [[Int]])

-- instance FromJSON OperatorTermSpec where
--   parseJSON = withObject "Term" $ \v ->
--     OperatorTermSpec <$> v .: "expression" <*> v .:? "sites"

operatorFromJSON :: IsBasis t => Basis t -> Value -> Parser (Operator t)
operatorFromJSON basis = withObject "Operator" $ \v -> do
  ((t :| ts) :: NonEmpty (Expr t)) <- v .: "terms"
  pure $ mkOperator basis (foldl' (+) t ts)
-}

-- newtype BasisAndHamiltonianConfig = BasisAndHamiltonianConfig SomeOperator
--
-- instance FromJSON BasisAndHamiltonianConfig where
--   parseJSON = withObject "Config" $ \v -> do
--     basis <- v .: "basis"
--     withSomeBasis basis $ \basis' ->
--       case Data.Aeson.KeyMap.lookup (Data.Aeson.Key.fromString "hamiltonian") v of
--         Just h -> BasisAndHamiltonianConfig . SomeOperator <$> operatorFromJSON basis' h
--         Nothing -> fail "missing 'hamiltonian' key"
--
-- hamiltonianFromYAML :: HasCallStack => Text -> IO SomeOperator
-- hamiltonianFromYAML path = (\(BasisAndHamiltonianConfig x) -> x) <$> objectFromYAML "Hamiltonian" path
--
-- newtype ObservablesSpec = ObservablesSpec (Vector SomeOperator)
--
-- instance FromJSON ObservablesSpec where
--   parseJSON = withObject "Config" $ \v -> do
--     basis <- v .: "basis"
--     withSomeBasis basis $ \basis' ->
--       case Data.Aeson.KeyMap.lookup (Data.Aeson.Key.fromString "observables") v of
--         Just h ->
--           flip (withArray "Observables") h $
--             fmap ObservablesSpec
--               . G.mapM (fmap SomeOperator . operatorFromJSON basis')
--         Nothing -> pure $ ObservablesSpec G.empty
