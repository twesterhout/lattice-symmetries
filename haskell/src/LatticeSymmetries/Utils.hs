{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE UndecidableInstances #-}
{-# OPTIONS_GHC -Wno-orphans #-}

-- |
-- Module      : LatticeSymmetries.Utils
-- Description : Random utilities
-- Copyright   : (c) Tom Westerhout, 2022
-- Stability   : experimental
module LatticeSymmetries.Utils
  ( -- ** Looping constructs
    loopM
  , iFoldM
  , rightM
  , sortVectorBy
  , unique

    -- ** Error handling
  -- , throwC
  -- , propagateErrorToC

    -- ** String handling
  , peekUtf8
  , newCString
  , newCencoded
  , decodeCString
  , toPrettyText
  , renderDoc
  , eitherToParser
  , viaAeson
  , viaAesonIO
  , ls_hs_destroy_string
  , MutablePtr

    -- ** Testing utilities
  , ApproxEq (..)
  )
where

import Control.Monad.ST (runST)
import Data.Aeson
import Data.Aeson qualified as Aeson
import Data.Aeson.Types (Parser)
import Data.Aeson.Types qualified as Aeson
import Data.ByteString (packCString)
import Data.ByteString.Internal (ByteString (..))
import Data.Complex
import Data.Containers.ListUtils (nubOrd)
import Data.Vector.Algorithms.Intro qualified
import Data.Vector.Generic qualified as G
import Foreign.C.String (CString)
import Foreign.C.Types (CChar, CDouble (..))
import Foreign.ForeignPtr (withForeignPtr)
import Foreign.Marshal.Alloc (free, mallocBytes)
import Foreign.Marshal.Utils (copyBytes)
import Foreign.Ptr (Ptr, castPtr)
import Foreign.Storable (Storable (..))
import Language.Halide (Arguments (..), UnCurry (..))
import Prettyprinter (Pretty (..))
import Prettyprinter qualified as Pretty
import Prettyprinter.Render.Text (renderStrict)

type MutablePtr a = Ptr a

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

sortVectorBy :: G.Vector v a => (a -> a -> Ordering) -> v a -> v a
sortVectorBy comp v = runST $ do
  buffer <- G.thaw v
  Data.Vector.Algorithms.Intro.sortBy comp buffer
  G.unsafeFreeze buffer

-- foreign import ccall unsafe "ls_hs_error"
--   ls_hs_error :: CString -> IO ()
-- ls_hs_error :: CString -> IO ()
-- ls_hs_error = undefined

-- | Invoke the 'ls_hs_error' error handling function with the given message.
-- throwC :: HasCallStack => a -> Text -> IO a
-- throwC def msg = do
--   useAsCString (encodeUtf8 msg) ls_hs_error
--   pure def

-- propagateErrorToC :: (HasCallStack, NFData a) => a -> IO a -> IO a
-- propagateErrorToC def = handleAnyDeep (throwC def . show)

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

newCencoded :: (HasCallStack, Data.Aeson.ToJSON a) => a -> IO CString
newCencoded = newCString . {-(\x -> trace (decodeUtf8 x) x) .-} toStrict . Data.Aeson.encode

-- | Read JSON from a 'CString'.
decodeCString :: FromJSON a => CString -> IO (Either Text a)
decodeCString cStr = do
  s <- packCString cStr
  case eitherDecode (fromStrict s) of
    Right x -> pure $ Right x
    Left msg -> pure $ Left (toText msg)

eitherToParser :: Either Text a -> Parser a
eitherToParser = either (fail . toString) pure

renderDoc :: Pretty.Doc ann -> Text
renderDoc = renderStrict . Pretty.layoutPretty (Pretty.LayoutOptions Pretty.Unbounded)

toPrettyText :: Pretty a => a -> Text
toPrettyText = renderDoc . pretty

instance FromJSON a => FromJSON (Complex a) where
  parseJSON (Aeson.Object v) = (:+) <$> v .: "real" <*> v .: "imag"
  parseJSON (Aeson.Array (G.toList -> [r, i])) = (:+) <$> parseJSON r <*> parseJSON i
  parseJSON invalid = Aeson.prependFailure "parsing Complex failed, " (Aeson.typeMismatch "Object or Array" invalid)

instance ToJSON a => ToJSON (Complex a) where
  toJSON (a :+ b) = object ["real" .= a, "imag" .= b]

class ArgumentsFromJSON (ts :: [Type]) where
  argumentsFromJSON :: [Aeson.Value] -> Aeson.Parser (Arguments ts)
  expectedNumberArgs :: Proxy ts -> Int

instance ArgumentsFromJSON '[] where
  argumentsFromJSON [] = pure Nil
  argumentsFromJSON _ = fail "parseJSON: expected an empty array"
  expectedNumberArgs _ = 0

instance (FromJSON t, ArgumentsFromJSON ts) => ArgumentsFromJSON (t ': ts) where
  argumentsFromJSON (x : xs) = (:::) <$> parseJSON x <*> argumentsFromJSON xs
  argumentsFromJSON [] = fail "parseJSON: expected an non-empty array"
  expectedNumberArgs _ = 1 + expectedNumberArgs (Proxy @ts)

instance ArgumentsFromJSON ts => FromJSON (Arguments ts) where
  parseJSON = Aeson.withArray "Arguments" $ \v ->
    if G.length v == n
      then argumentsFromJSON (G.toList v)
      else fail $ "expected " <> show n <> " arguments, but " <> show (G.length v) <> " were provided"
    where
      n = expectedNumberArgs (Proxy @ts)

viaAeson
  :: forall args r f
   . ( UnCurry f args (Either Text r)
     , FromJSON (Arguments args)
     , ToJSON r
     )
  => f
  -> CString
  -> IO CString
viaAeson f cStr =
  decodeCString @(Arguments args) cStr
    >>= either (newCencoded @(Either Text r) . Left) (newCencoded . uncurryG f)

viaAesonIO
  :: forall args r f
   . (UnCurry f args (IO (Either Text r)), FromJSON (Arguments args), ToJSON r)
  => f
  -> CString
  -> IO CString
viaAesonIO f cStr =
  decodeCString @(Arguments args) cStr
    >>= either (newCencoded @(Either Text r) . Left) (newCencoded <=< uncurryG f)

ls_hs_destroy_string :: CString -> IO ()
ls_hs_destroy_string = free

rightM :: (b -> IO c) -> Either a b -> IO (Either a c)
rightM _ (Left a) = pure (Left a)
rightM f (Right b) = Right <$> f b

unique :: (Ord a, G.Vector v a) => v a -> v a
unique = G.fromList . nubOrd . G.toList
