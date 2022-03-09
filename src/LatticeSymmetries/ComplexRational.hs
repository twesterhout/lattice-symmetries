module LatticeSymmetries.ComplexRational
  ( ComplexRational (..),
    ComplexFloating (..),
    realPart,
    imagPart,
    conjugate,
    magnitudeSquared,
  )
where

import Data.Complex (Complex (..))

data ComplexRational = ComplexRational {-# UNPACK #-} !Rational {-# UNPACK #-} !Rational
  deriving (Eq, Show)

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

class Fractional a => ComplexFloating a where
  toComplexDouble :: a -> Complex Double
  fromComplexDouble :: Complex Double -> a

instance ComplexFloating ComplexRational where
  toComplexDouble (ComplexRational r i) = (fromRational r) :+ (fromRational i)
  fromComplexDouble (r :+ i) = ComplexRational (toRational r) (toRational i)
