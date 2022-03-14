module LatticeSymmetries.Indexing where

import qualified Data.Vector.Generic as G
import Data.Vector.Storable (Vector)
import qualified Data.Vector.Storable as S

data LookupRange = LookupRange {-# UNPACK #-} !CPtrdiff {-# UNPACK #-} !CPtrdiff
  deriving (Show, Eq)

data BasisLookupCache = BasisLookupCache
  { blcRepresentatives :: !(Vector Word64),
    blcRanges :: !(Vector LookupRange)
  }
