{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE TemplateHaskell #-}

-- |
-- Copyright: (c) 2021 Tom Westerhout
-- SPDX-License-Identifier: BSD-3-Clause
-- Maintainer: Tom Westerhout <14264576+twesterhout@users.noreply.github.com>
module LatticeSymmetries.Context
  ( lsCtx,
  )
where

import qualified Data.Map as Map
import Foreign.C.Types (CInt)
import Language.C.Inline.Context (Context (..))
import qualified Language.C.Types as Types
import qualified Language.Haskell.TH as TH

newtype ErrorCode = ErrorCode CInt

lsTypesTable :: Map.Map Types.TypeSpecifier TH.TypeQ
lsTypesTable =
  Map.fromList
    [ (Types.TypeName "ls_error_code", [t|ErrorCode|])
    ]

-- | Provides type mappings for better interoperability with "Language.C.Inline".
lsCtx :: Context
lsCtx = mempty {ctxTypesTable = lsTypesTable}
