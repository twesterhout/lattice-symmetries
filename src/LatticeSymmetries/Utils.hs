{-# LANGUAGE ImplicitParams #-}
{-# LANGUAGE RankNTypes #-}

module LatticeSymmetries.Utils
  ( loopM,
    iFoldM,
    logDebug',
    logInfo',
    logWarning',
    logError',
    ApproxEq (..),
    peekUtf8,
    newCString,
    decodeCString,
    propagateErrorToC,
    throwC,
    toPrettyText,
  )
where

import Colog
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
import GHC.Stack (freezeCallStack, popCallStack)
import Prettyprinter (Doc, Pretty (..))
import qualified Prettyprinter as Pretty
import Prettyprinter.Render.Text (renderStrict)
import System.IO.Unsafe (unsafePerformIO)

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

defaultLogAction :: LogAction IO Message
defaultLogAction = filterBySeverity Info msgSeverity $ cmap fmtMessage logTextStderr

currentLogAction :: IORef (LogAction IO Message)
currentLogAction = unsafePerformIO $ newIORef defaultLogAction
{-# NOINLINE currentLogAction #-}

withDefaultLogger :: HasCallStack => LoggerT Message IO () -> IO ()
withDefaultLogger f = withFrozenCallStack $ do
  logAction <- readIORef currentLogAction
  usingLoggerT logAction f

withFrozenPoppedCallStack :: HasCallStack => (HasCallStack => a) -> a
withFrozenPoppedCallStack do_this =
  let ?callStack = freezeCallStack . popCallStack . popCallStack $ callStack
   in do_this

logDebug' :: HasCallStack => Text -> IO ()
logDebug' t = withFrozenPoppedCallStack $ withDefaultLogger (logDebug t)

logInfo' :: HasCallStack => Text -> IO ()
logInfo' t = withFrozenPoppedCallStack $ withDefaultLogger (logInfo t)

logWarning' :: HasCallStack => Text -> IO ()
logWarning' t = withFrozenPoppedCallStack $ withDefaultLogger (logWarning t)

logError' :: HasCallStack => Text -> IO ()
logError' t = withFrozenPoppedCallStack $ withDefaultLogger (logError t)

foreign import ccall unsafe "ls_hs_error"
  ls_hs_error :: CString -> IO ()

throwC :: HasCallStack => a -> Text -> IO a
throwC x₀ msg = do
  logError' "Throwing error to C ..."
  useAsCString (encodeUtf8 msg) ls_hs_error
  pure x₀

propagateErrorToC :: HasCallStack => a -> SomeException -> IO a
propagateErrorToC x₀ = \e -> throwC x₀ (show e)

-- ls_hs_fatal_error :: HasCallStack => CString -> CString -> IO ()
-- ls_hs_fatal_error c_func c_msg = withFrozenCallStack $ do
--   func <- peekUtf8 c_func
--   msg <- peekUtf8 c_msg
--   logError' $ "[" <> func <> "] " <> msg
--   exitFailure

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

decodeCString :: (HasCallStack, FromJSON a) => CString -> IO a
decodeCString cStr = do
  s <- packCString cStr
  case eitherDecode (fromStrict s) of
    Right x -> pure x
    Left msg -> error (toText msg)

toPrettyText :: Pretty a => a -> Text
toPrettyText = renderStrict . Pretty.layoutPretty (Pretty.LayoutOptions Pretty.Unbounded) . pretty
