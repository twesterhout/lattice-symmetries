module LatticeSymmetries.Generator
  ( SpinIndex (..),
    HasSiteIndex (..),
    SpinGeneratorType (..),
    FermionGeneratorType (..),
    Generator (..),
    -- HasMatrixRepresentation (..),
  )
where

import Data.Bits
import qualified Data.Text as Text
import qualified Data.Vector.Generic as G
import LatticeSymmetries.BitString
import LatticeSymmetries.Dense
import LatticeSymmetries.NonbranchingTerm
import Prettyprinter (Doc, Pretty (..))
import qualified Prettyprinter as Pretty
import Prelude hiding (Product, Sum, identity, toList)

data SpinIndex = SpinUp | SpinDown
  deriving stock (Show, Eq, Ord)

instance Pretty SpinIndex where
  pretty SpinUp = "↑"
  pretty SpinDown = "↓"

class HasSiteIndex i where
  getSiteIndex :: i -> Int
  mapSiteIndex :: (Int -> Int) -> i -> i

instance HasSiteIndex Int where
  getSiteIndex = id
  mapSiteIndex f i = f i

instance HasSiteIndex (SpinIndex, Int) where
  getSiteIndex (_, i) = i
  mapSiteIndex f (σ, i) = (σ, f i)

data SpinGeneratorType = SpinIdentity | SpinZ | SpinPlus | SpinMinus
  deriving stock (Eq, Ord, Show, Enum, Bounded, Generic)

instance Pretty SpinGeneratorType where
  pretty x = case x of
    SpinIdentity -> "1"
    SpinZ -> "σᶻ"
    SpinPlus -> "σ⁺"
    SpinMinus -> "σ⁻"

data FermionGeneratorType = FermionIdentity | FermionCount | FermionCreate | FermionAnnihilate
  deriving stock (Eq, Ord, Show, Enum, Bounded, Generic)

instance Pretty FermionGeneratorType where
  pretty x = case x of
    FermionIdentity -> "1"
    FermionCount -> "n"
    FermionCreate -> "c†"
    FermionAnnihilate -> "c"

data Generator i g = Generator !i !g
  deriving stock (Eq, Ord, Show, Generic)

toSubscript :: HasCallStack => Int -> Text
toSubscript n = Text.map h (show n)
  where
    h '0' = '₀'
    h '1' = '₁'
    h '2' = '₂'
    h '3' = '₃'
    h '4' = '₄'
    h '5' = '₅'
    h '6' = '₆'
    h '7' = '₇'
    h '8' = '₈'
    h '9' = '₉'
    h _ = error "invalid character"

instance Pretty g => Pretty (Generator Int g) where
  pretty (Generator i g) = pretty g <> pretty (Text.unpack (toSubscript i))

instance Pretty g => Pretty (Generator (SpinIndex, Int) g) where
  pretty (Generator (σ, i) g) = pretty g <> pretty σ <> pretty (Text.unpack (toSubscript i))

-- class HasMatrixRepresentation g where
--   matrixRepresentation :: (G.Vector v c, Num c) => g -> DenseMatrix v c
--
-- instance HasMatrixRepresentation SpinGeneratorType where
--   matrixRepresentation = spinMatrixRepresentation
--
-- instance HasMatrixRepresentation FermionGeneratorType where
--   matrixRepresentation = fermionMatrixRepresentation
--
-- spinMatrixRepresentation :: (G.Vector v c, Num c) => SpinGeneratorType -> DenseMatrix v c
-- spinMatrixRepresentation g = fromList $ case g of
--   SpinIdentity -> [[1, 0], [0, 1]]
--   SpinZ -> [[1, 0], [0, -1]]
--   SpinPlus -> [[0, 1], [0, 0]]
--   SpinMinus -> [[0, 0], [1, 0]]
--
-- fermionMatrixRepresentation :: (G.Vector v c, Num c) => FermionGeneratorType -> DenseMatrix v c
-- fermionMatrixRepresentation g = fromList $ case g of
--   FermionIdentity -> [[1, 0], [0, 1]]
--   FermionCount -> [[1, 0], [0, 0]]
--   FermionCreate -> [[0, 1], [0, 0]]
--   FermionAnnihilate -> [[0, 0], [1, 0]]

instance HasNonbranchingRepresentation (Generator Int SpinGeneratorType) where
  nonbranchingRepresentation (Generator _ SpinIdentity) =
    NonbranchingTerm 1 zeroBits zeroBits zeroBits zeroBits zeroBits
  nonbranchingRepresentation (Generator i SpinZ) =
    NonbranchingTerm (-1) zeroBits zeroBits zeroBits zeroBits (bit i)
  nonbranchingRepresentation (Generator i SpinPlus) =
    NonbranchingTerm 1 (bit i) (bit i) zeroBits (bit i) zeroBits
  nonbranchingRepresentation (Generator i SpinMinus) =
    NonbranchingTerm 1 (bit i) zeroBits (bit i) (bit i) zeroBits

instance HasNonbranchingRepresentation (Generator Int FermionGeneratorType) where
  nonbranchingRepresentation (Generator _ FermionIdentity) =
    NonbranchingTerm 1 zeroBits zeroBits zeroBits zeroBits zeroBits
  nonbranchingRepresentation (Generator i FermionCount) =
    NonbranchingTerm 1 (bit i) (bit i) (bit i) zeroBits zeroBits
  nonbranchingRepresentation (Generator i FermionCreate) =
    NonbranchingTerm 1 (bit i) (bit i) zeroBits (bit i) (BitString (bit i - 1))
  nonbranchingRepresentation (Generator i FermionAnnihilate) =
    NonbranchingTerm 1 (bit i) zeroBits (bit i) (bit i) (BitString (bit i - 1))
