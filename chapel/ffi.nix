{ stdenv
, chapel
, lattice-symmetries
, version
}:

stdenv.mkDerivation {
  pname = "lattice-symmetries-chapel-ffi";
  inherit version;
  dontUnpack = true;

  buildPhase = ''
    c2chapel \
      ${lattice-symmetries.haskell}/include/lattice_symmetries_functions.h \
      -DLS_C2CHAPEL \
      -I${chapel}/runtime/include \
      -I${lattice-symmetries.kernels}/include \
      >FFI.chpl

    # Remove the declaration of chpl_external_array since it's already
    # present in the ExternalArray module
    sed -i -e '/extern record chpl_external_array/,+5d' FFI.chpl
    sed -i -E '/chpl_make_external_array(_ptr)?(_free)?\(/d' FFI.chpl
    sed -i -e '/cleanupOpaqueArray(/d' FFI.chpl
    sed -i -e '/chpl_free_external_array(/d' FFI.chpl
    sed -i -e '/chpl_call_free_func(/d' FFI.chpl
    sed -i 's/extern type ls_hs_scalar = _Complex double/extern type ls_hs_scalar = complex(128)/' FFI.chpl
  '';

  installPhase = "install -m 644 FFI.chpl $out";
  buildInputs = [ lattice-symmetries.haskell lattice-symmetries.kernels ];
  nativeBuildInputs = [ chapel ];
}
