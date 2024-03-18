{ chapel
, halide
, lattice-symmetries
, removeReferencesTo
, stdenv
, version
}:

stdenv.mkDerivation {
  pname = "lattice-symmetries-chapel";
  inherit version;
  src = ./.;

  configurePhase = "ln --symbolic ${lattice-symmetries.ffi} src/FFI.chpl";
  preBuild = ''
    makeFlagsArray+=(
      PREFIX="$out"
      CHPL_CFLAGS="-I${lattice-symmetries.kernels_v2}/include -I${lattice-symmetries.haskell}/include --no-ieee-float --local --no-debug --fast"
      CHPL_LDFLAGS="-L${lattice-symmetries.haskell}/lib"
      HALIDE_PATH="${halide}"
    )
  '';
  preInstall = ''
    for f in $(ls lib); do
      chapelFixupBinary lib/$f
      # These are only needed during compilation
      remove-references-to -t ${lattice-symmetries.kernels_v2} lib/$f
    done
  '';

  buildInputs = [
    halide
    lattice-symmetries.kernels_v2
    lattice-symmetries.haskell
  ];
  nativeBuildInputs = [
    chapel
    removeReferencesTo
  ];
}
