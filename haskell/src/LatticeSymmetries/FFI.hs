{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE TemplateHaskell #-}

module LatticeSymmetries.FFI where

import Data.Typeable (typeOf)
import Foreign
import Foreign.C.Types
import Language.C.Inline.Unsafe qualified as CU
import LatticeSymmetries.Context

importLS

ls_hs_destroy_object :: forall a. IsCobject a => (Ptr a -> IO ()) -> Ptr a -> IO ()
ls_hs_destroy_object custom p = do
  let p' = castPtr @a @Cobject p
  refcount <- [CU.exp| int { ls_hs_internal_object_dec_ref_count($(ls_hs_object* p')) } |]
  case compare refcount 1 of
    GT -> do
      -- don't have to do anything yet
      pure ()
    EQ -> do
      -- call the destructor
      custom p
      freeStablePtr . castPtrToStablePtr =<< [CU.exp| void* { $(ls_hs_object* p')->haskell_payload } |]
    LT ->
      -- should never happen
      error $ "ls_hs_destroy_object: refcount < 1 when destroying " <> show (typeOf p)

newCobject :: IsCobject b => a -> IO (Ptr b)
newCobject = fmap castPtr . newCobjectSized 0

newCobjectSized :: Int -> a -> IO (Ptr Cobject)
newCobjectSized sizeBytes x = do
  payload <- castStablePtrToPtr <$> newStablePtr x
  p <- callocBytes $ max (fromIntegral [CU.pure| size_t { sizeof(ls_hs_object) } |]) sizeBytes
  [CU.block| void { ls_hs_internal_object_init($(ls_hs_object* p), 1, $(void* payload)); } |]
  pure p

withCobject :: forall b a r. IsCobject b => Ptr b -> (a -> IO r) -> IO r
withCobject (castPtr @b @Cobject -> p) f =
  f
    =<< (deRefStablePtr . castPtrToStablePtr)
    =<< [CU.exp| void* { $(ls_hs_object* p)->haskell_payload } |]

foldCobject :: forall b a r. IsCobject b => (a -> IO r) -> Ptr b -> IO r
foldCobject f p = withCobject p f
