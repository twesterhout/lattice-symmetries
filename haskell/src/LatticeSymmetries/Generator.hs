{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE UndecidableInstances #-}

-- |
-- Module      : LatticeSymmetries.Generator
-- Description : Bosonic and fermionic algebra generators
-- Copyright   : (c) Tom Westerhout, 2022
-- Stability   : experimental
module LatticeSymmetries.Generator
  ( SpinIndex (..)
  , SpinGeneratorType (..)
  , FermionGeneratorType (..)
  , Generator (..)
  , ParticleTy (..)
  , BasisState (..)
  , unsafeCastBasisState
  , invertBasisState
  , prettyBitString
  , ParticleTag (..)
  , particleTagToType
  , particleDispatch
  , matchParticleType2
  , IndexType
  , GeneratorType
  , withConstraint
  -- HasSiteIndex (..),
  -- , getSiteIndex
  , mapSiteIndex
  , HasIdentity (..)
  , flattenIndex
  , unFlattenIndex
  -- HasMatrixRepresentation (..),
  )
where

import Data.Aeson
import Data.Bits
import Data.Constraint
import Data.Constraint.Deferrable
import Data.Text qualified as Text
import LatticeSymmetries.BitString
import LatticeSymmetries.NonbranchingTerm
import Prettyprinter (Pretty (..))
import Prettyprinter qualified as Pretty
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

-- | Hilbert space basis vector parametrized by the particle type.
data BasisState (t :: ParticleTy) = BasisState {-# UNPACK #-} !Int !BitString
  deriving stock (Show, Eq, Ord)

instance Typeable t => Pretty (BasisState t) where
  pretty x = case particleDispatch @t of
    SpinTag -> prettySpin x
    SpinfulFermionTag -> prettyFermion x
    SpinlessFermionTag ->
      let (BasisState n bits) = x
       in prettySpin (BasisState n bits :: BasisState 'SpinTy)
    where
      prettySpin (BasisState n bits) = "|" <> prettyBitString n (unBitString bits) <> "⟩"
      prettyFermion (BasisState n bits) =
        let up = unBitString bits `shiftR` (n `div` 2)
         in mconcat
              [ "|"
              , prettyBitString (n `div` 2) up
              , "⟩"
              , "|"
              , prettyBitString (n `div` 2) (unBitString bits)
              , "⟩"
              ]

unsafeCastBasisState :: BasisState t1 -> BasisState t2
unsafeCastBasisState (BasisState n bits) = BasisState n bits

-- | Invert all bits in a BasisState
invertBasisState :: BasisState t -> BasisState t
invertBasisState (BasisState n bits) = BasisState n (bits `xor` mask)
  where
    mask = BitString (bit n - 1)

prettyBitString :: Int -> Integer -> Pretty.Doc ann
prettyBitString n bits = mconcat $ prettyBool . testBit bits <$> reverse [0 .. n - 1]
  where
    prettyBool True = "1"
    prettyBool False = "0"

-- | Type-level analog of 'ParticleTy'.
data ParticleTag (t :: ParticleTy) where
  SpinTag :: ParticleTag 'SpinTy
  SpinfulFermionTag :: ParticleTag 'SpinfulFermionTy
  SpinlessFermionTag :: ParticleTag 'SpinlessFermionTy

deriving instance Show (ParticleTag t)

deriving instance Eq (ParticleTag t)

-- | Get the runtime representation of the particle type.
particleTagToType :: ParticleTag t -> ParticleTy
particleTagToType x = case x of
  SpinTag -> SpinTy
  SpinfulFermionTag -> SpinfulFermionTy
  SpinlessFermionTag -> SpinlessFermionTy

-- | Get the type-level representation of the particle type.
particleDispatch :: forall (t :: ParticleTy). (HasCallStack, Typeable t) => ParticleTag t
particleDispatch
  | Just HRefl <- eqTypeRep (typeRep @t) (typeRep @'SpinTy) = SpinTag
  | Just HRefl <- eqTypeRep (typeRep @t) (typeRep @'SpinfulFermionTy) = SpinfulFermionTag
  | Just HRefl <- eqTypeRep (typeRep @t) (typeRep @'SpinlessFermionTy) = SpinlessFermionTag
  | otherwise = error "this should never happen by construction"

matchParticleType2
  :: forall (t1 :: ParticleTy) (t2 :: ParticleTy) proxy1 proxy2
   . (Typeable t1, Typeable t2)
  => proxy1 t1
  -> proxy2 t2
  -> Maybe (t1 :~~: t2)
matchParticleType2 _ _ = case eqTypeRep (typeRep @t1) (typeRep @t2) of
  Just HRefl -> Just HRefl
  Nothing -> Nothing

withConstraint
  :: forall (c :: ParticleTy -> Constraint) (t :: ParticleTy) a
   . (HasCallStack, Typeable t, c 'SpinTy, c 'SpinfulFermionTy, c 'SpinlessFermionTy)
  => (c t => a)
  -> a
withConstraint f = case particleDispatch @t of
  SpinTag -> f
  SpinfulFermionTag -> f
  SpinlessFermionTag -> f

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

instance ToJSON SpinIndex where
  toJSON = String . renderStrict . Pretty.layoutCompact . pretty

instance FromJSON SpinIndex where
  parseJSON (Data.Aeson.String t)
    | t == "↑" || Text.toLower t == "up" = pure SpinUp
    | t == "↓" || Text.toLower t == "down" = pure SpinDown
  parseJSON x = do
    (k :: Int) <- parseJSON x
    case k of
      0 -> pure SpinUp
      1 -> pure SpinDown
      _ -> fail $ "invalid spin index: " <> show k

class HasIdentity g where
  isIdentity :: g -> Bool

-- | Generators for the algebra of spin-1/2 particles.
data SpinGeneratorType
  = -- | Identity \( 1 = \begin{pmatrix} 1 & 0\\ 0 & 1\end{pmatrix} \)
    SpinIdentity
  | -- | Pauli matrix \( \sigma^z = \begin{pmatrix} 1 & 0\\ 0 & -1\end{pmatrix} \)
    SpinZ
  | -- | \( \sigma^{+} = \sigma^x + 𝕚\sigma^y = \begin{pmatrix} 0 & 1\\ 0 & 0 \end{pmatrix} \)
    SpinPlus
  | -- | \( \sigma^{\-\} = \sigma^x - 𝕚\sigma^y = \begin{pmatrix} 0 & 0\\ 1 & 0 \end{pmatrix} \))
    SpinMinus
  deriving stock (Eq, Ord, Show, Enum, Bounded, Generic)

instance Pretty SpinGeneratorType where
  pretty x = case x of
    SpinIdentity -> "I"
    SpinZ -> "σᶻ"
    SpinPlus -> "σ⁺"
    SpinMinus -> "σ⁻"

instance HasIdentity SpinGeneratorType where
  isIdentity = (== SpinIdentity)

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
    FermionIdentity -> "I"
    FermionCount -> "n"
    FermionCreate -> "c†"
    FermionAnnihilate -> "c"

instance HasIdentity FermionGeneratorType where
  isIdentity = (== FermionIdentity)

type family IndexType (t :: ParticleTy) where
  IndexType 'SpinTy = Int
  IndexType 'SpinfulFermionTy = (SpinIndex, Int)
  IndexType 'SpinlessFermionTy = Int

type family GeneratorType (t :: ParticleTy) where
  GeneratorType 'SpinTy = SpinGeneratorType
  GeneratorType 'SpinfulFermionTy = FermionGeneratorType
  GeneratorType 'SpinlessFermionTy = FermionGeneratorType

-- | A generator (either spin or fermionic) which is not associated with an index @i@. The index
-- could be the site index or a tuple of spin and site indices.
data Generator i g = Generator !i !g
  deriving stock (Eq, Ord, Show, Generic)

instance HasIdentity g => HasIdentity (Generator i g) where
  isIdentity (Generator _ g) = isIdentity g

class HasSiteIndex i where
  getSiteIndex :: i -> Int
  mapSiteIndex :: (Int -> Int) -> i -> i

instance HasSiteIndex Int where
  getSiteIndex = id
  mapSiteIndex f = f

instance HasSiteIndex (SpinIndex, Int) where
  getSiteIndex (_, i) = i
  mapSiteIndex f (σ, i) = (σ, f i)

flattenIndex :: forall t. Typeable t => Int -> IndexType t -> Int
flattenIndex numberSites = case particleDispatch @t of
  SpinTag -> id
  SpinlessFermionTag -> id
  SpinfulFermionTag -> \case
    (σ, k) | k < numberSites -> fromEnum σ * numberSites + k
    (_, k) -> error $ "index exceeds numberSites: " <> show k <> " > " <> show numberSites

unFlattenIndex :: forall t. Typeable t => Int -> Int -> IndexType t
unFlattenIndex numberSites = case particleDispatch @t of
  SpinTag -> id
  SpinlessFermionTag -> id
  SpinfulFermionTag -> \i ->
    if i >= numberSites then (SpinDown, i - numberSites) else (SpinUp, i)

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
    h c = error $ "invalid character: '" <> Text.singleton c <> "'"

instance (Pretty g, HasIdentity g) => Pretty (Generator Int g) where
  pretty (Generator i g) = pretty g <> (if isIdentity g then "" else pretty (toSubscript i))

instance (Pretty g, HasIdentity g) => Pretty (Generator (SpinIndex, Int) g) where
  pretty (Generator (σ, i) g) = pretty g <> (if isIdentity g then "" else pretty (toSubscript i) <> pretty σ)

instance HasNonbranchingRepresentation (Generator Int SpinGeneratorType) where
  nonbranchingRepresentation :: Generator Int SpinGeneratorType -> NonbranchingTerm
  nonbranchingRepresentation (Generator _ SpinIdentity) =
    NonbranchingTerm 1 zeroBits zeroBits zeroBits zeroBits zeroBits
  nonbranchingRepresentation (Generator i SpinZ) =
    NonbranchingTerm 1 zeroBits zeroBits zeroBits zeroBits (bit i)
  nonbranchingRepresentation (Generator i SpinPlus) =
    NonbranchingTerm 1 (bit i) (bit i) zeroBits (bit i) zeroBits
  nonbranchingRepresentation (Generator i SpinMinus) =
    NonbranchingTerm 1 (bit i) zeroBits (bit i) (bit i) zeroBits

instance HasNonbranchingRepresentation (Generator Int FermionGeneratorType) where
  nonbranchingRepresentation :: Generator Int FermionGeneratorType -> NonbranchingTerm
  nonbranchingRepresentation (Generator _ FermionIdentity) =
    NonbranchingTerm 1 zeroBits zeroBits zeroBits zeroBits zeroBits
  nonbranchingRepresentation (Generator i FermionCount) =
    NonbranchingTerm 1 (bit i) (bit i) (bit i) zeroBits zeroBits
  nonbranchingRepresentation (Generator i FermionCreate) =
    NonbranchingTerm 1 (bit i) (bit i) zeroBits (bit i) (BitString (bit i - 1))
  nonbranchingRepresentation (Generator i FermionAnnihilate) =
    NonbranchingTerm 1 (bit i) zeroBits (bit i) (bit i) (BitString (bit i - 1))
