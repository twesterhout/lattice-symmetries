module ForeignLibrary () where

import Control.Concurrent
import Data.Complex
import Data.List.Split (chunksOf)
import qualified Data.Text as T
import qualified Data.Vector.Generic as G
import Foreign.C.String (CString, peekCString)
import Foreign.C.Types (CBool (..), CInt (..), CPtrdiff (..), CSize (..), CUInt (..), CUShort (..))
import Foreign.ForeignPtr
import Foreign.Marshal.Array (peekArray)
import Foreign.Marshal.Utils (fromBool)
import Foreign.Ptr (FunPtr, Ptr)
import Foreign.StablePtr
import LatticeSymmetries
import LatticeSymmetries.Algebra
import LatticeSymmetries.Basis
import LatticeSymmetries.ComplexRational
import LatticeSymmetries.FFI
import LatticeSymmetries.Generator
import LatticeSymmetries.Operator
import LatticeSymmetries.Parser
-- import qualified LatticeSymmetries.Types as Types
import LatticeSymmetries.Utils
import System.Mem
import Text.PrettyPrint.ANSI.Leijen (hardline, pretty, putDoc)

foreign export ccall "ls_hs_hdf5_create_dataset_u64"
  ls_hs_hdf5_create_dataset_u64 :: CString -> CString -> CUInt -> Ptr Word64 -> IO ()

foreign export ccall "ls_hs_hdf5_create_dataset_f32"
  ls_hs_hdf5_create_dataset_f32 :: CString -> CString -> CUInt -> Ptr Word64 -> IO ()

foreign export ccall "ls_hs_hdf5_create_dataset_f64"
  ls_hs_hdf5_create_dataset_f64 :: CString -> CString -> CUInt -> Ptr Word64 -> IO ()

foreign export ccall "ls_hs_hdf5_create_dataset_c64"
  ls_hs_hdf5_create_dataset_c64 :: CString -> CString -> CUInt -> Ptr Word64 -> IO ()

foreign export ccall "ls_hs_hdf5_create_dataset_c128"
  ls_hs_hdf5_create_dataset_c128 :: CString -> CString -> CUInt -> Ptr Word64 -> IO ()

foreign export ccall "ls_hs_hdf5_write_chunk_u64"
  ls_hs_hdf5_write_chunk_u64 :: CString -> CString -> CUInt -> Ptr Word64 -> Ptr Word64 -> Ptr Word64 -> IO ()

foreign export ccall "ls_hs_hdf5_write_chunk_f64"
  ls_hs_hdf5_write_chunk_f64 :: CString -> CString -> CUInt -> Ptr Word64 -> Ptr Word64 -> Ptr Double -> IO ()

foreign export ccall "ls_hs_hdf5_read_chunk_u64"
  ls_hs_hdf5_read_chunk_u64 :: CString -> CString -> CUInt -> Ptr Word64 -> Ptr Word64 -> Ptr Word64 -> IO ()

foreign export ccall "ls_hs_hdf5_read_chunk_f64"
  ls_hs_hdf5_read_chunk_f64 :: CString -> CString -> CUInt -> Ptr Word64 -> Ptr Word64 -> Ptr Double -> IO ()

foreign export ccall "ls_hs_hdf5_get_dataset_rank"
  ls_hs_hdf5_get_dataset_rank :: CString -> CString -> IO CUInt

foreign export ccall "ls_hs_hdf5_get_dataset_shape"
  ls_hs_hdf5_get_dataset_shape :: CString -> CString -> Ptr Word64 -> IO ()

-- foreign export ccall "ls_hs_basis_and_hamiltonian_from_yaml"
--   ls_hs_basis_and_hamiltonian_from_yaml :: CString -> Ptr Types.SpinBasisWrapper -> Ptr Types.OperatorWrapper -> IO ()

-- foreign export ccall "ls_hs_destroy_spin_basis"
--   ls_hs_destroy_spin_basis :: Ptr Types.SpinBasisWrapper -> IO ()

-- foreign export ccall "ls_hs_destroy_operator"
--   ls_hs_destroy_operator :: Ptr Types.OperatorWrapper -> IO ()

-- foreign export ccall "ls_hs_create_basis"
--   ls_hs_create_basis :: Cparticle_type -> CInt -> CInt -> CInt -> IO (Ptr Cbasis)

foreign export ccall "ls_hs_clone_basis"
  ls_hs_clone_basis :: Ptr Cbasis -> IO (Ptr Cbasis)

foreign export ccall "ls_hs_create_spin_basis_from_json"
  ls_hs_create_spin_basis_from_json :: CString -> IO (Ptr Cbasis)

foreign export ccall "ls_hs_create_spinless_fermion_basis_from_json"
  ls_hs_create_spinless_fermion_basis_from_json :: CString -> IO (Ptr Cbasis)

foreign export ccall "ls_hs_create_spinful_fermion_basis_from_json"
  ls_hs_create_spinful_fermion_basis_from_json :: CString -> IO (Ptr Cbasis)

-- foreign export ccall "ls_hs_spin_chain_10_basis"
--   ls_hs_spin_chain_10_basis :: IO (Ptr Cbasis)

-- foreign export ccall "ls_hs_spin_kagome_12_basis"
--   ls_hs_spin_kagome_12_basis :: IO (Ptr Cbasis)

-- foreign export ccall "ls_hs_spin_kagome_16_basis"
--   ls_hs_spin_kagome_16_basis :: IO (Ptr Cbasis)

-- foreign export ccall "ls_hs_spin_square_4x4_basis"
--   ls_hs_spin_square_4x4_basis :: IO (Ptr Cbasis)

foreign export ccall "ls_hs_min_state_estimate"
  ls_hs_min_state_estimate :: Ptr Cbasis -> IO Word64

foreign export ccall "ls_hs_max_state_estimate"
  ls_hs_max_state_estimate :: Ptr Cbasis -> IO Word64

ls_hs_basis_has_fixed_hamming_weight :: Ptr Cbasis -> IO CBool
ls_hs_basis_has_fixed_hamming_weight basis =
  fromBool
    <$> withReconstructedBasis basis (pure . hasFixedHammingWeight . basisHeader)

ls_hs_basis_build :: Ptr Cbasis -> IO ()
ls_hs_basis_build basis = withReconstructedBasis basis basisBuild

foreign export ccall "ls_hs_basis_build"
  ls_hs_basis_build :: Ptr Cbasis -> IO ()

foreign export ccall "ls_hs_basis_has_fixed_hamming_weight"
  ls_hs_basis_has_fixed_hamming_weight :: Ptr Cbasis -> IO CBool

ls_hs_create_operator :: Ptr Cbasis -> CString -> CInt -> CInt -> Ptr CInt -> IO (Ptr Coperator)
ls_hs_create_operator basisPtr cStr numberTuples tupleSize tuplesPtr =
  withReconstructedBasis basisPtr $ \basis -> do
    indices <-
      chunksOf (fromIntegral tupleSize)
        <$> fmap fromIntegral
        <$> peekArray (fromIntegral (numberTuples * tupleSize)) tuplesPtr
    s <- peekCString cStr
    let !operator = operatorFromString basis (T.pack s) indices
    borrowCoperator operator

foreign export ccall "ls_hs_create_operator"
  ls_hs_create_operator :: Ptr Cbasis -> CString -> CInt -> CInt -> Ptr CInt -> IO (Ptr Coperator)

ls_hs_operator_plus :: Ptr Coperator -> Ptr Coperator -> IO (Ptr Coperator)
ls_hs_operator_plus aPtr bPtr =
  withReconstructedOperator aPtr $ \a ->
    withSameTypeAs a (withReconstructedOperator bPtr) $ \b ->
      borrowCoperator (a + b)

foreign export ccall "ls_hs_operator_plus"
  ls_hs_operator_plus :: Ptr Coperator -> Ptr Coperator -> IO (Ptr Coperator)

ls_hs_print_terms :: Ptr Coperator -> IO ()
ls_hs_print_terms p =
  withReconstructedOperator p $ \operator ->
    putDoc (pretty (opTerms (opHeader operator)) <> hardline)

foreign export ccall "ls_hs_print_terms"
  ls_hs_print_terms :: Ptr Coperator -> IO ()

-- foreign export ccall "ls_hs_fatal_error"
--   ls_hs_fatal_error :: CString -> CString -> IO ()
