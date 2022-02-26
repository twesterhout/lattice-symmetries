module ForeignLibrary where

import Data.Complex
import Foreign.C.String (CString)
import Foreign.C.Types (CBool (..), CInt (..), CUInt (..), CUShort (..))
import Foreign.Ptr (FunPtr, Ptr)
import LatticeSymmetries
import LatticeSymmetries.Types

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
  ls_hs_basis_and_hamiltonian_from_yaml :: CString -> Ptr SpinBasisWrapper -> Ptr OperatorWrapper -> IO ()

foreign export ccall "ls_hs_destroy_spin_basis"
  ls_hs_destroy_spin_basis :: Ptr SpinBasisWrapper -> IO ()

foreign export ccall "ls_hs_destroy_operator"
  ls_hs_destroy_operator :: Ptr OperatorWrapper -> IO ()

main :: IO ()
main = putStrLn "Hello world!"
