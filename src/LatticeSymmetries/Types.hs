module LatticeSymmetries.Types where

import Foreign.ForeignPtr

-- | Exceptions thrown when an error occurs in @liblattice_symmetries@.
data LatticeSymmetriesException = LatticeSymmetriesException {eCode :: Int, eMessage :: Text}
  deriving stock (Show)

instance Exception LatticeSymmetriesException

data SpinEDException = SpinEDException Text
  deriving stock (Show)

instance Exception SpinEDException

data {-# CTYPE "lattice_symmetries/lattice_symmetries.h" "ls_symmetry" #-} CSymmetry

data {-# CTYPE "lattice_symmetries/lattice_symmetries.h" "ls_group" #-} CGroup

data {-# CTYPE "lattice_symmetries/lattice_symmetries.h" "ls_spin_basis" #-} Cspin_basis

data {-# CTYPE "lattice_symmetries/lattice_symmetries.h" "ls_flat_spin_basis" #-} CFlatSpinBasis

data {-# CTYPE "lattice_symmetries/lattice_symmetries.h" "ls_operator" #-} Coperator

data {-# CTYPE "lattice_symmetries/lattice_symmetries.h" "ls_interaction" #-} Cinteraction

newtype Symmetry = Symmetry (ForeignPtr CSymmetry)

newtype SymmetryGroup = SymmetryGroup (ForeignPtr CGroup)

newtype SpinBasis = SpinBasis (ForeignPtr Cspin_basis)

newtype FlatSpinBasis = FlatSpinBasis (ForeignPtr CFlatSpinBasis)

newtype Operator = Operator (ForeignPtr Coperator)

newtype Interaction = Interaction (ForeignPtr Cinteraction)
