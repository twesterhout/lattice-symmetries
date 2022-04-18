module LatticeSymmetries.ComplexRational
  ( ComplexRational (..),
    ConvertibleToComplexDouble (..),
    realPart,
    imagPart,
    conjugate,
    magnitudeSquared,
    Cscalar,
  )
where

import Data.Complex (Complex (..))
import Foreign.C.Types (CDouble (..))
import Text.PrettyPrint.ANSI.Leijen (Pretty (..))
import qualified Text.PrettyPrint.ANSI.Leijen as Pretty

data ComplexRational = ComplexRational {-# UNPACK #-} !Rational {-# UNPACK #-} !Rational
  deriving stock (Eq, Show)

prettyRational :: Rational -> Pretty.Doc
prettyRational x
  | realToFrac (fromRational x :: Double) == x = Pretty.double (fromRational x :: Double)
  | otherwise = Pretty.rational x

instance Pretty ComplexRational where
  pretty (ComplexRational r i)
    | i == 0 = prettyRational r
    | otherwise = Pretty.parens $ prettyRational r <> Pretty.text " + " <> prettyRational i <> "𝕚"

realPart :: ComplexRational -> Rational
realPart (ComplexRational r _) = r
{-# INLINE realPart #-}

imagPart :: ComplexRational -> Rational
imagPart (ComplexRational _ i) = i
{-# INLINE imagPart #-}

conjugate :: ComplexRational -> ComplexRational
conjugate (ComplexRational r i) = (ComplexRational r (-i))
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
  (ComplexRational r i) / (ComplexRational r' i') =
    ComplexRational ((r * r' + i * i') / d) ((-r * i' + i * r') / d)
    where
      d = r' * r' + i' * i'
  fromRational a = ComplexRational a 0

-- class Fractional a => ComplexFloating a where
--   toComplexDouble :: a -> Complex Double
--   fromComplexDouble :: Complex Double -> a

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
