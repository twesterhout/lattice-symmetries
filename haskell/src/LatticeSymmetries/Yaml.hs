module LatticeSymmetries.Yaml where

import Data.Aeson
import Data.Aeson.Key qualified
import Data.Aeson.KeyMap qualified
import Data.Aeson.Types (JSONPathElement (..), Parser, parserThrowError)
import Data.Vector (Vector)
import Data.Vector.Generic qualified as G
import Data.Yaml.Aeson (decodeFileWithWarnings, prettyPrintParseException)
import LatticeSymmetries.Basis
import LatticeSymmetries.Expr
import LatticeSymmetries.Generator
import LatticeSymmetries.Operator

data ConfigSpec = ConfigSpec !SomeBasis !(Maybe SomeOperator) !(Vector SomeOperator)

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
configFromYAML = objectFromYAML "config"

operatorFromJSON :: IsBasis t => Basis t -> Value -> Parser (Operator t)
operatorFromJSON basis = withObject "Operator" $ \v -> do
  ((t :| ts) :: NonEmpty (Expr t)) <- v .: "terms"
  pure $ mkOperator basis (foldl' (+) t ts)

objectFromYAML :: (HasCallStack, FromJSON a) => Text -> Text -> IO a
objectFromYAML _ filename = do
  r <- decodeFileWithWarnings (toString filename)
  case r of
    Left e -> error $ toText $ prettyPrintParseException e
    Right (_, x) -> pure x
