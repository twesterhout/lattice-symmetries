{-# LANGUAGE OverloadedStrings #-}

module LatticeSymmetries.Context
  ( importLS
  )
where

import Data.Map.Strict qualified as Map
import Language.C.Inline qualified as C
import Language.C.Inline.Context (Context (..))
import Language.C.Types (CIdentifier, TypeSpecifier (..))
import Language.Haskell.TH (DecsQ, Q, TypeQ, lookupTypeName)
import Language.Haskell.TH qualified as TH
import Prelude hiding (optional)

importLS :: DecsQ
importLS =
  concat
    <$> sequence
      [ C.context =<< lsCxt
      , C.include "<stdatomic.h>"
      , C.include "<lattice_symmetries_types.h>"
      ]

lsCxt :: Q C.Context
lsCxt = do
  typePairs <- Map.fromList <$> lsTypePairs
  pure $ C.fptrCtx <> C.bsCtx <> C.baseCtx <> mempty {ctxTypesTable = typePairs}

lsTypePairs :: Q [(TypeSpecifier, TypeQ)]
lsTypePairs =
  optionals
    [ ("ls_hs_nonbranching_terms", "Cnonbranching_terms")
    , ("ls_hs_basis", "Cbasis")
    , ("ls_hs_operator", "Coperator")
    , ("Halide::Internal::Dimension", "CxxDimension")
    , ("Halide::Internal::FusedPair", "FusedPair")
    , ("Halide::Internal::PrefetchDirective", "PrefetchDirective")
    , ("Halide::Internal::ReductionVariable", "ReductionVariable")
    , ("Halide::Internal::Split", "Split")
    , ("Halide::Internal::StageSchedule", "CxxStageSchedule")
    , ("Halide::Argument", "CxxArgument")
    , ("Halide::Buffer", "CxxBuffer")
    , ("Halide::LoopLevel", "CxxLoopLevel")
    , ("Halide::Stage", "CxxStage")
    , ("Halide::Range", "CxxRange")
    , ("Halide::RDom", "CxxRDom")
    , ("halide_buffer_t", "Language.Halide.Buffer.RawHalideBuffer")
    , ("halide_device_interface_t", "HalideDeviceInterface")
    , ("halide_dimension_t", "HalideDimension")
    , ("halide_trace_event_t", "TraceEvent")
    ]
  where
    optional :: (CIdentifier, String) -> Q [(TypeSpecifier, TypeQ)]
    optional (cName, hsName) = do
      hsType <- lookupTypeName hsName
      pure $ maybe [] (\x -> [(TypeName cName, pure (TH.ConT x))]) hsType
    optionals :: [(CIdentifier, String)] -> Q [(TypeSpecifier, TypeQ)]
    optionals pairs = concat <$> mapM optional pairs
