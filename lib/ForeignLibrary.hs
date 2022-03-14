module ForeignLibrary () where

import Data.Complex
import Data.List.Split (chunksOf)
import qualified Data.Text as T
import Foreign.C.String (CString, peekCString)
import Foreign.C.Types (CBool (..), CInt (..), CUInt (..), CUShort (..))
import Foreign.Marshal.Array (peekArray)
import Foreign.Ptr (FunPtr, Ptr)
import LatticeSymmetries
import LatticeSymmetries.Algebra
import LatticeSymmetries.Basis
import LatticeSymmetries.Generator
import LatticeSymmetries.Operator
import LatticeSymmetries.Parser
import qualified LatticeSymmetries.Types as Types
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

foreign export ccall "ls_hs_basis_and_hamiltonian_from_yaml"
  ls_hs_basis_and_hamiltonian_from_yaml :: CString -> Ptr Types.SpinBasisWrapper -> Ptr Types.OperatorWrapper -> IO ()

foreign export ccall "ls_hs_destroy_spin_basis"
  ls_hs_destroy_spin_basis :: Ptr Types.SpinBasisWrapper -> IO ()

foreign export ccall "ls_hs_destroy_operator"
  ls_hs_destroy_operator :: Ptr Types.OperatorWrapper -> IO ()

ls_hs_create_spin_basis :: CInt -> CInt -> IO (Ptr Cbasis)
ls_hs_create_spin_basis numberSites hammingWeight = do
  putStrLn "Creating SpinBasis ..."
  createCbasis $ SpinBasis numberSites' hammingWeight'
  where
    numberSites' = fromIntegral numberSites
    hammingWeight'
      | hammingWeight == -1 = Nothing
      | otherwise = Just (fromIntegral hammingWeight)

ls_hs_create_spinful_fermionic_basis :: CInt -> CInt -> CInt -> IO (Ptr Cbasis)
ls_hs_create_spinful_fermionic_basis numberSites numberUp numberDown = do
  putStrLn "Creating SpinfulFermionicBasis ..."
  createCbasis $ SpinfulFermionicBasis numberSites' occupation
  where
    numberSites' = fromIntegral numberSites
    occupation
      | numberUp == -1 && numberDown == -1 = SpinfulNoOccupation
      | numberUp >= 0 && numberDown == -1 = SpinfulTotalParticles (fromIntegral numberUp)
      | otherwise = SpinfulPerSector (fromIntegral numberUp) (fromIntegral numberDown)

ls_hs_create_basis_kernels :: Ptr Cbasis -> IO (Ptr Cbasis_kernels)
ls_hs_create_basis_kernels p = withReconstructedBasis p createCbasis_kernels

foreign export ccall "ls_hs_create_spin_basis"
  ls_hs_create_spin_basis :: CInt -> CInt -> IO (Ptr Cbasis)

foreign export ccall "ls_hs_create_spinful_fermionic_basis"
  ls_hs_create_spinful_fermionic_basis :: CInt -> CInt -> CInt -> IO (Ptr Cbasis)

foreign export ccall "ls_hs_destroy_basis_v2"
  destroyCbasis :: Ptr Cbasis -> IO ()

foreign export ccall "ls_hs_create_basis_kernels"
  ls_hs_create_basis_kernels :: Ptr Cbasis -> IO (Ptr Cbasis_kernels)

foreign export ccall "ls_hs_destroy_basis_kernels"
  destroyCbasis_kernels :: Ptr Cbasis_kernels -> IO ()

ls_hs_create_operator :: Ptr Cbasis -> CString -> CInt -> CInt -> Ptr CInt -> IO (Ptr Coperator)
ls_hs_create_operator basisPtr cStr numberTuples tupleSize tuplesPtr =
  withReconstructedBasis basisPtr $ \basis -> do
    indices <-
      chunksOf (fromIntegral tupleSize)
        <$> fmap fromIntegral
        <$> peekArray (fromIntegral (numberTuples * tupleSize)) tuplesPtr
    s <- peekCString cStr
    let operator = operatorFromString basis (T.pack s) indices
    putDoc (pretty (opTerms operator) <> hardline)
    createCoperator operator

foreign export ccall "ls_hs_create_operator"
  ls_hs_create_operator :: Ptr Cbasis -> CString -> CInt -> CInt -> Ptr CInt -> IO (Ptr Coperator)

foreign export ccall "ls_hs_destroy_operator_v2"
  destroyCoperator :: Ptr Coperator -> IO ()
