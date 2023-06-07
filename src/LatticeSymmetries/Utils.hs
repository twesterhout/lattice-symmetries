-- |
-- Module      : LatticeSymmetries.Utils
-- Description : Random utilities
-- Copyright   : (c) Tom Westerhout, 2022
-- Stability   : experimental
module LatticeSymmetries.Utils
  ( -- ** Looping constructs
    loopM
  , iFoldM

    -- ** Error handling
  , throwC
  , propagateErrorToC

    -- ** String handling
  , peekUtf8
  , newCString
  , decodeCString
  , toPrettyText

    -- ** Testing utilities
  , ApproxEq (..)
  )
where

import Data.Aeson
import Data.ByteString (packCString, useAsCString)
import Data.ByteString.Internal (ByteString (..))
import Data.Complex
import Foreign.C.String (CString)
import Foreign.C.Types (CChar, CDouble (..))
import Foreign.ForeignPtr (withForeignPtr)
import Foreign.Marshal.Alloc (mallocBytes)
import Foreign.Marshal.Utils (copyBytes)
import Foreign.Ptr (Ptr, castPtr)
import Foreign.Storable (Storable (..))
import Prettyprinter (Pretty (..))
import Prettyprinter qualified as Pretty
import Prettyprinter.Render.Text (renderStrict)

loopM :: Monad m => i -> (i -> Bool) -> (i -> i) -> (i -> m ()) -> m ()
loopM i₀ cond inc action = go i₀
  where
    go !i
      | cond i = do () <- action i; go (inc i)
      | otherwise = pure ()
{-# INLINE loopM #-}

iFoldM :: Monad m => i -> (i -> Bool) -> (i -> i) -> a -> (a -> i -> m a) -> m a
iFoldM i₀ cond inc x₀ action = go x₀ i₀
  where
    go !x !i
      | cond i = do !x' <- action x i; go x' (inc i)
      | otherwise = pure x
{-# INLINE iFoldM #-}

foreign import ccall unsafe "ls_hs_error"
  ls_hs_error :: CString -> IO ()

-- | Invoke the 'ls_hs_error' error handling function with the given message.
throwC :: HasCallStack => Text -> IO a
throwC msg = do
  useAsCString (encodeUtf8 msg) ls_hs_error
  error "this should never happen, because ls_hs_error should not return"

propagateErrorToC :: HasCallStack => SomeException -> IO a
propagateErrorToC = throwC . show

infix 4 ≈

class ApproxEq a where
  (≈) :: a -> a -> Bool
  approx :: Double -> Double -> a -> a -> Bool

instance ApproxEq Double where
  approx rtol atol a b = abs (a - b) <= max atol (rtol * max (abs a) (abs b))
  (≈) = approx 1.4901161193847656e-8 8.881784197001252e-16

instance ApproxEq (Complex Double) where
  approx rtol atol (ra :+ ia) (rb :+ ib) =
    approx rtol atol ra rb && approx rtol atol ia ib
  (≈) = approx 1.4901161193847656e-8 8.881784197001252e-16

instance ApproxEq (Complex CDouble) where
  approx = coerce (approx :: Double -> Double -> Complex Double -> Complex Double -> Bool)
  (≈) = coerce ((≈) :: Complex Double -> Complex Double -> Bool)

peekUtf8 :: CString -> IO Text
peekUtf8 c_str = decodeUtf8 <$> packCString c_str

newCString :: ByteString -> IO CString
newCString (PS fp _ l) = do
  (buf :: Ptr CChar) <- mallocBytes (l + 1)
  withForeignPtr fp $ \(p :: Ptr Word8) -> do
    copyBytes buf (castPtr p) l
    pokeByteOff buf l (0 :: CChar)
  pure buf

-- | Read JSON from a 'CString'.
decodeCString :: (HasCallStack, FromJSON a) => CString -> IO a
decodeCString cStr = do
  s <- packCString cStr
  case eitherDecode (fromStrict s) of
    Right x -> pure x
    Left msg -> error (toText msg)

toPrettyText :: Pretty a => a -> Text
toPrettyText = renderStrict . Pretty.layoutPretty (Pretty.LayoutOptions Pretty.Unbounded) . pretty
