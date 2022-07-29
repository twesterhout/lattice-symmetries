module LatticeSymmetries.IO where

import Data.Aeson
import Data.Aeson.Types (typeMismatch)
import Data.Complex
import Data.Scientific

data SymmetrySpec = SymmetrySpec !(NonEmpty Int) !Int
  deriving stock (Read, Show, Eq)

instance FromJSON SymmetrySpec where
  parseJSON = withObject "symmetry" $ \v ->
    SymmetrySpec
      <$> v .: "permutation"
      <*> v .: "sector"

instance ToJSON SymmetrySpec where
  toJSON (SymmetrySpec permutation sector) =
    object ["permutation" .= permutation, "sector" .= sector]

data BasisSpec = BasisSpec !Int !(Maybe Int) !(Maybe Int) ![SymmetrySpec]
  deriving stock (Read, Show, Eq)

instance FromJSON BasisSpec where
  parseJSON = withObject "basis" $ \v ->
    BasisSpec
      <$> v .: "number_spins"
      <*> v .:? "hamming_weight"
      <*> v .:? "spin_inversion"
      <*> v .:! "symmetries" .!= []

instance ToJSON BasisSpec where
  toJSON (BasisSpec numberSpins hammingWeight spinInversion symmetries) =
    object
      [ "number_spins" .= numberSpins,
        "hamming_weight" .= maybe Null toJSON hammingWeight,
        "spin_inversion" .= maybe Null toJSON spinInversion,
        "symmetries" .= symmetries
      ]

newtype WrappedBasisSpec = WrappedBasisSpec BasisSpec

instance FromJSON WrappedBasisSpec where
  parseJSON = withObject "config" $ \v -> WrappedBasisSpec <$> v .: "basis"

instance FromJSON (Complex Double) where
  parseJSON (Number x) = pure . fromReal . toRealFloat $ x
    where
      fromReal :: Num a => a -> Complex a
      fromReal x' = x' :+ 0
  parseJSON v@(Array xs) = case (toList xs) of
    [re, im] -> (:+) <$> parseJSON re <*> parseJSON im
    _ -> typeMismatch "Complex" v
  parseJSON v = typeMismatch "Complex" v

newtype DenseMatrixSpec = DenseMatrixSpec [[Complex Double]]
  deriving stock (Read, Show, Eq)
  deriving newtype (FromJSON)

data InteractionSpec = InteractionSpec ![[Complex Double]] ![[Int]]
  deriving stock (Read, Show, Eq)

instance FromJSON InteractionSpec where
  parseJSON = withObject "interaction" $ \v ->
    InteractionSpec
      <$> v .: "matrix"
      <*> v .: "sites"

data OperatorSpec = OperatorSpec !Text !(NonEmpty InteractionSpec)
  deriving stock (Read, Show)

instance FromJSON OperatorSpec where
  parseJSON = withObject "operator" $ \v ->
    OperatorSpec
      <$> v .: "name"
      <*> v .: "terms"

data ConfigSpec = ConfigSpec !BasisSpec !OperatorSpec
  deriving stock (Read, Show)

instance FromJSON ConfigSpec where
  parseJSON = withObject "config" $ \v ->
    ConfigSpec
      <$> v .: "basis"
      <*> v .: "hamiltonian"
