module LatticeSymmetries
  ( foo,
  )
where

import Foreign.C.Types (CInt (..))

foo :: IO CInt
foo = do
  putStrLn "Hello world!"
  return 123

foreign export ccall "ls_hs_foo" foo :: IO CInt
