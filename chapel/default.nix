{ stdenv
, chapel
, chapelFixupBinary
, lattice-symmetries
, version
}:

stdenv.mkDerivation {
  pname = "lattice-symmetries-chapel";
  inherit version;
  src = ./.;

  configurePhase = "ln --symbolic ${lattice-symmetries.ffi} src/FFI.chpl";
  makeFlags = [
    "PREFIX=$(out)"
    "OPTIMIZATION=--fast"
    "CHPL_CFLAGS='-I${lattice-symmetries.kernels}/include'"
    "CHPL_LDFLAGS='-L${lattice-symmetries.haskell.lib}/lib'"
  ];
  preInstall = ''
    for f in $(ls lib); do
      chapelFixupBinary lib/$f
    done
    for f in $(ls bin); do
      chapelFixupBinary bin/$f
    done
  '';

  buildInputs = [ lattice-symmetries.kernels lattice-symmetries.haskell.lib ];
  nativeBuildInputs = [ chapel chapelFixupBinary ];
}
