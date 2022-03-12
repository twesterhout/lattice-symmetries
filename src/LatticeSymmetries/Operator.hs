module LatticeSymmetries.Operator where

import Data.Vector (Vector)
import LatticeSymmetries.Algebra
import LatticeSymmetries.Basis

data Operator c basis = Operator
  { opTerms :: !(Polynomial c (Generator (IndexType basis) (GeneratorType basis)))
  }

-- typedef struct ls_hs_nonbranching_terms {
--   int number_terms;
--   int number_bits;
--   // number_words = ceil(number_bits / 64)
--   ls_hs_scalar const *v; // array of shape [number_terms]
--   uint64_t const *m;     // array of shape [number_terms, number_words]
--   uint64_t const *l;     // array of shape [number_terms, number_words]
--   uint64_t const *r;     // array of shape [number_terms, number_words]
--   uint64_t const *x;     // array of shape [number_terms, number_words]
--   uint64_t const *s;     // array of shape [number_terms, number_words]
--   // all arrays are contiguous in row-major order
-- } ls_hs_nonbranching_terms;

type Cscalar = Complex CDouble

data Cnonbranching_terms = Cnonbranching_terms
  { cnonbranching_terms_number_terms :: {-# UNPACK #-} !CInt,
    cnonbranching_terms_number_bits :: {-# UNPACK #-} !CInt,
    cnonbranching_terms_v :: {-# UNPACK #-} !(Ptr Cscalar),
    cnonbranching_terms_m :: {-# UNPACK #-} !(Ptr Word64),
    cnonbranching_terms_l :: {-# UNPACK #-} !(Ptr Word64),
    cnonbranching_terms_r :: {-# UNPACK #-} !(Ptr Word64),
    cnonbranching_terms_x :: {-# UNPACK #-} !(Ptr Word64),
    cnonbranching_terms_s :: {-# UNPACK #-} !(Ptr Word64)
  }

--
-- typedef struct ls_hs_operator {
--   ls_hs_basis const *basis;
--   ls_hs_nonbranching_terms const *off_diag_terms;
--   ls_hs_nonbranching_terms const *diag_terms;
--   bool needs_projection;
-- } ls_hs_operator;

data Coperator = Coperator
  { coperator_basis :: {-# UNPACK #-} !(Ptr Cbasis),
    coperator_off_diag_terms :: {-# UNPACK #-} !(Ptr Cnonbranching_terms),
    coperator_diag_terms :: {-# UNPACK #-} !(Ptr Cnonbranching_terms),
    coperator_needs_projection :: {-# UNPACK #-} !CBool
  }

-- applyOperator :: Operator c basis -> BasisState -> (Vector c, Vector BasisState)
-- applyOperator = undefined
