-- |
-- Module      : LatticeSymmetries.ComplexRational
-- Description : Arbitrary-precision complex numbers
-- Copyright   : (c) Tom Westerhout, 2022
-- Stability   : experimental
module LatticeSymmetries.ComplexRational
  ( ComplexRational (..)
  , ℂ
  , ConvertibleToComplexDouble (..)
  , realPart
  , imagPart
  , conjugate
  , magnitudeSquared
  , Cscalar
  )
where

import Data.Complex (Complex (..))
import Foreign.C.Types (CDouble (..))
import Prettyprinter (Doc, Pretty (..))
import Prettyprinter qualified as Pretty

-- | Arbitrary precision complex number \(\mathbb{C}\) built from two 'Rational's.
data ComplexRational = ComplexRational {-# UNPACK #-} !Rational {-# UNPACK #-} !Rational
  deriving stock (Eq, Show)

type ℂ = ComplexRational

prettyRational :: Rational -> Doc ann
prettyRational x
  | realToFrac (fromRational x :: Double) == x = pretty (fromRational x :: Double)
  | otherwise = Pretty.parens $ pretty (numerator x) <> "/" <> pretty (denominator x)

instance Pretty ComplexRational where
  pretty (ComplexRational r i)
    | i == 0 = prettyRational r
    | r == 0 && i == 1 = "ⅈ"
    | r == 0 && i == -1 = "-ⅈ"
    | r == 0 = prettyRational i <> "ⅈ"
    | otherwise = Pretty.parens $ prettyRational r <> (if i > 0 then " + " else " - ") <> prettyRational (abs i) <> "ⅈ"

realPart :: ComplexRational -> Rational
realPart (ComplexRational r _) = r
{-# INLINE realPart #-}

imagPart :: ComplexRational -> Rational
imagPart (ComplexRational _ i) = i
{-# INLINE imagPart #-}

conjugate :: ComplexRational -> ComplexRational
conjugate (ComplexRational r i) = ComplexRational r (-i)
{-# INLINE conjugate #-}

magnitudeSquared :: ComplexRational -> Rational
magnitudeSquared (ComplexRational r i) = r * r + i * i
{-# INLINE magnitudeSquared #-}

instance Num ComplexRational where
  {-# INLINE (+) #-}
  {-# INLINE (-) #-}
  {-# INLINE (*) #-}
  {-# INLINE fromInteger #-}
  (ComplexRational r i) + (ComplexRational r' i') = ComplexRational (r + r') (i + i')
  (ComplexRational r i) - (ComplexRational r' i') = ComplexRational (r - r') (i - i')
  (ComplexRational r i) * (ComplexRational r' i') = ComplexRational (r * r' - i * i') (r * i' + i * r')
  negate (ComplexRational r i) = ComplexRational (-r) (-i)
  abs _ = error "Num instance of ComplexRational does not implement abs"
  signum _ = error "Num instance of ComplexRational does not implement signum"
  fromInteger n = ComplexRational (fromInteger n) 0

instance Fractional ComplexRational where
  {-# INLINE (/) #-}
  {-# INLINE fromRational #-}
  (ComplexRational r i) / z'@(ComplexRational r' i') =
    ComplexRational ((r * r' + i * i') / d) ((-r * i' + i * r') / d)
    where
      d = magnitudeSquared z'
  fromRational a = ComplexRational a 0

type Cscalar = Complex CDouble

class Fractional a => ConvertibleToComplexDouble a where
  toComplexDouble :: a -> Cscalar
  fromComplexDouble :: Cscalar -> a

instance ConvertibleToComplexDouble ComplexRational where
  toComplexDouble (ComplexRational r i) = (fromRational r) :+ (fromRational i)
  fromComplexDouble (r :+ i) = ComplexRational (toRational r) (toRational i)

instance ConvertibleToComplexDouble (Complex CDouble) where
  toComplexDouble = id
  fromComplexDouble = id

instance ConvertibleToComplexDouble (Complex Double) where
  toComplexDouble = coerce
  fromComplexDouble = coerce
