{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskellQuotes #-}
{-# OPTIONS_GHC -Wno-orphans #-}

module LatticeSymmetries.Context
  ( importLS
  , RawIsRepresentativeKernel
  , RawStateInfoKernel
  , RawStateToIndexKernel
  , Cnonbranching_terms
  , Cobject
  , Cbasis
  , Cexpr
  , Coperator
  , Cscalar
  , IsCobject
  )
where

import Data.Aeson (FromJSON (parseJSON), ToJSON (toJSON))
import Data.Complex (Complex)
import Data.Map.Strict qualified as Map
import Foreign (FunPtr, Ptr, WordPtr (WordPtr), ptrToWordPtr, wordPtrToPtr)
import Foreign.C (CDouble)
import Language.C.Inline qualified as C
import Language.C.Inline.Context (Context (..))
import Language.C.Types (CIdentifier, TypeSpecifier (..))
import Language.Halide (RawHalideBuffer)
import Language.Haskell.TH (DecsQ, Q, TypeQ, lookupTypeName)
import Language.Haskell.TH qualified as TH
import Prelude hiding (optional)

importLS :: DecsQ
importLS =
  concat
    <$> sequence
      [ C.context =<< lsCxt
      , C.include "<stdatomic.h>"
      , C.include "<stdint.h>" -- uint64_t, int64_t
      , C.include "<stdlib.h>" -- abort
      , C.include "<stdio.h>" -- fprintf, stderr
      , C.include "<string.h>" -- memset
      , C.include "<lattice_symmetries_types.h>"
      , C.include "<HalideRuntime.h>" -- halide_buffer_t
      ]

type RawIsRepresentativeKernel = Ptr RawHalideBuffer -> Ptr RawHalideBuffer -> IO ()
type RawStateInfoKernel = Ptr RawHalideBuffer -> Ptr RawHalideBuffer -> Ptr RawHalideBuffer -> IO ()
type RawStateToIndexKernel = Ptr RawHalideBuffer -> Ptr RawHalideBuffer -> IO ()

type Cscalar = Complex CDouble

data Cnonbranching_terms

data Cobject

data Cbasis

data Cexpr

data Coperator

class Typeable a => IsCobject a

instance IsCobject Cobject

instance IsCobject Cbasis

instance IsCobject Cexpr

instance IsCobject Coperator

lsCxt :: Q C.Context
lsCxt = do
  typePairs <- Map.fromList <$> lsTypePairs

  pure $ C.fptrCtx <> C.bsCtx <> C.baseCtx <> mempty {ctxTypesTable = typePairs}

lsTypePairs :: Q [(TypeSpecifier, TypeQ)]
lsTypePairs = (mandatory <>) <$> optionals [("halide_buffer_t", "RawHalideBuffer")]
  where
    optional :: (CIdentifier, String) -> Q [(TypeSpecifier, TypeQ)]
    optional (cName, hsName) = do
      hsType <- lookupTypeName hsName
      pure $ maybe [] (\x -> [(TypeName cName, pure (TH.ConT x))]) hsType
    optionals :: [(CIdentifier, String)] -> Q [(TypeSpecifier, TypeQ)]
    optionals pairs = concat <$> mapM optional pairs
    mandatory =
      [ (TypeName "ls_hs_is_representative_kernel_type_v2", [t|FunPtr RawIsRepresentativeKernel|])
      , (TypeName "ls_hs_state_info_kernel_type_v2", [t|FunPtr RawStateInfoKernel|])
      , (TypeName "ls_hs_state_to_index_kernel_type", [t|FunPtr RawStateToIndexKernel|])
      , (TypeName "ls_hs_basis", [t|Cbasis|])
      , (TypeName "ls_hs_expr", [t|Cexpr|])
      , (TypeName "ls_hs_nonbranching_terms", [t|Cnonbranching_terms|])
      , (TypeName "ls_hs_object", [t|Cobject|])
      , (TypeName "ls_hs_operator", [t|Coperator|])
      , (TypeName "ls_hs_scalar", [t|Cscalar|])
      ]

instance IsCobject a => ToJSON (Ptr a) where
  toJSON = toJSON . (\(WordPtr w) -> w) . ptrToWordPtr

instance IsCobject a => FromJSON (Ptr a) where
  parseJSON = fmap (wordPtrToPtr . WordPtr) . parseJSON
