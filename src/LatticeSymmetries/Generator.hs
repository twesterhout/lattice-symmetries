{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE GADTs #-}

-- |
-- Module      : LatticeSymmetries.Generator
-- Description : Bosonic and fermionic algebra generators
-- Copyright   : (c) Tom Westerhout, 2022
-- Stability   : experimental
module LatticeSymmetries.Generator
  ( SpinIndex (..),
    SpinGeneratorType (..),
    FermionGeneratorType (..),
    Generator (..),
    ParticleTy (..),
    ParticleTag (..),
    particleTagToType,
    particleDispatch,
    IndexType (..),
    GeneratorType (..),
    HasSiteIndex (..),
    -- HasMatrixRepresentation (..),
  )
where

import Data.Aeson
import Data.Aeson.Types (Pair)
import Data.Bits
import qualified Data.Text as Text
import qualified Data.Vector.Generic as G
import LatticeSymmetries.BitString
import LatticeSymmetries.Dense
import LatticeSymmetries.NonbranchingTerm
import Prettyprinter (Doc, Pretty (..))
import qualified Prettyprinter as Pretty
import Prettyprinter.Render.Text (renderStrict)
import Type.Reflection
import Prelude hiding (Product, Sum, identity, toList)

-- | Particle type
data ParticleTy
  = -- | A localized spin-1/2 particle.
    SpinTy
  | -- | An electron (i.e. a fermion with spin-1/2)
    SpinfulFermionTy
  | -- | A spinless fermion
    SpinlessFermionTy
  deriving stock (Show, Eq, Typeable)

instance Pretty ParticleTy where
  pretty SpinTy = "spin-1/2"
  pretty SpinfulFermionTy = "spinful-fermion"
  pretty SpinlessFermionTy = "spinless-fermion"

instance FromJSON ParticleTy where
  parseJSON = withText "ParticleTy" f
    where
      f t
        | t == "spin" || t == "spin-1/2" = pure SpinTy
        | t == "spinful" || t == "spinful-fermion" || t == "spinful fermion" = pure SpinfulFermionTy
        | t == "spinless" || t == "spinless-fermion" || t == "spinless fermion" =
            pure SpinlessFermionTy
        | otherwise = fail "invalid particle type"

instance ToJSON ParticleTy where
  toJSON = String . renderStrict . Pretty.layoutCompact . pretty

data ParticleTag (t :: ParticleTy) where
  SpinTag :: ParticleTag 'SpinTy
  SpinfulFermionTag :: ParticleTag 'SpinfulFermionTy
  SpinlessFermionTag :: ParticleTag 'SpinlessFermionTy

particleTagToType :: ParticleTag t -> ParticleTy
particleTagToType x = case x of
  SpinTag -> SpinTy
  SpinfulFermionTag -> SpinfulFermionTy
  SpinlessFermionTag -> SpinlessFermionTy

particleDispatch :: forall (t :: ParticleTy). (HasCallStack, Typeable t) => ParticleTag t
particleDispatch
  | Just HRefl <- eqTypeRep (typeRep @t) (typeRep @'SpinTy) = SpinTag
  | Just HRefl <- eqTypeRep (typeRep @t) (typeRep @'SpinfulFermionTy) = SpinfulFermionTag
  | Just HRefl <- eqTypeRep (typeRep @t) (typeRep @'SpinlessFermionTy) = SpinlessFermionTag
  | otherwise = error "this should never happen by construction"

-- | Index for the spin sector.
--
-- __Note the ordering__: spin up appears before spin down.
data SpinIndex
  = -- | ↑
    SpinUp
  | -- | ↓
    SpinDown
  deriving stock (Show, Eq, Ord)

instance Enum SpinIndex where
  toEnum 0 = SpinUp
  toEnum 1 = SpinDown
  toEnum i = error $ "invalid spin index: " <> show i <> "; expected either 0 or 1"
  fromEnum SpinUp = 0
  fromEnum SpinDown = 1

instance Pretty SpinIndex where
  pretty SpinUp = "↑"
  pretty SpinDown = "↓"

-- | Generators for the algebra of spin-1/2 particles.
data SpinGeneratorType
  = -- | Identity \( 1 = \begin{pmatrix} 1 & 0\\ 0 & 1\end{pmatrix} \)
    SpinIdentity
  | -- | Pauli matrix \( \sigma^z = \begin{pmatrix} 1 & 0\\ 0 & -1\end{pmatrix} \)
    SpinZ
  | -- | \( \sigma^{+} = \sigma^x + 𝕚\sigma^y = \begin{pmatrix} 0 & 1\\ 0 & 0 \end{pmatrix} \)
    SpinPlus
  | -- | \( \sigma^{-} = \sigma^x - 𝕚\sigma^y = \begin{pmatrix} 0 & 0\\ 1 & 0 \end{pmatrix} \)
    SpinMinus
  deriving stock (Eq, Ord, Show, Enum, Bounded, Generic)

instance Pretty SpinGeneratorType where
  pretty x = case x of
    SpinIdentity -> "1"
    SpinZ -> "σᶻ"
    SpinPlus -> "σ⁺"
    SpinMinus -> "σ⁻"

-- | Generators for the fermionic algebra.
data FermionGeneratorType
  = -- | Identity 𝟙
    FermionIdentity
  | -- | Number counting operator
    FermionCount
  | -- | Creation operator \( c^\dagger \)
    FermionCreate
  | -- | Annihilation operator \( c \)
    FermionAnnihilate
  deriving stock (Eq, Ord, Show, Enum, Bounded, Generic)

instance Pretty FermionGeneratorType where
  pretty x = case x of
    FermionIdentity -> "1"
    FermionCount -> "n"
    FermionCreate -> "c†"
    FermionAnnihilate -> "c"

type family IndexType (t :: ParticleTy) where
  IndexType 'SpinTy = Int
  IndexType 'SpinfulFermionTy = (SpinIndex, Int)
  IndexType 'SpinlessFermionTy = Int

type family GeneratorType (t :: ParticleTy) where
  GeneratorType 'SpinTy = SpinGeneratorType
  GeneratorType 'SpinfulFermionTy = FermionGeneratorType
  GeneratorType 'SpinlessFermionTy = FermionGeneratorType

type IsGeneratorType g = (Eq g, Ord g, Pretty g)

type IsIndexType i = (Eq i, Ord i, Pretty i)

-- | A generator (either spin or fermionic) which is not associated with an index @i@. The index
-- could be the site index or a tuple of spin and site indices.
data Generator i g = Generator !i !g
  deriving stock (Eq, Ord, Show, Generic)

class HasSiteIndex i where
  getSiteIndex :: i -> Int
  mapSiteIndex :: (Int -> Int) -> i -> i

instance HasSiteIndex Int where
  getSiteIndex = id
  mapSiteIndex f i = f i

instance HasSiteIndex (SpinIndex, Int) where
  getSiteIndex (_, i) = i
  mapSiteIndex f (σ, i) = (σ, f i)

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
  pretty (Generator i g) = pretty g <> pretty (toSubscript i)

instance Pretty g => Pretty (Generator (SpinIndex, Int) g) where
  pretty (Generator (σ, i) g) = pretty g <> pretty σ <> pretty (toSubscript i)

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
