{ stdenv
, chapel
, chapelFixupBinary
, halide
, lib
, lattice-symmetries
, version
, sse4_2Support ? stdenv.targetPlatform.sse4_2Support
}:

stdenv.mkDerivation {
  pname = "lattice-symmetries-chapel";
  inherit version;
  src = ./.;

  configurePhase = "ln --symbolic ${lattice-symmetries.ffi} src/FFI.chpl";
  preBuild = ''
    makeFlagsArray+=(
      PREFIX="$out"
      OPTIMIZATION="--debug -g --preserve-inlined-line-numbers --task-tracking --savec tmp"
      CHPL_CFLAGS="--no-ieee-float -I${lattice-symmetries.kernels}/include"
      CHPL_LDFLAGS="-L${lattice-symmetries.haskell.lib}/lib"
      HALIDE_PATH="${halide}"
    )
  ''; # + lib.optionalString sse4_2Support "makeFlagsArray+=( CHPL_TARGET_CPU=\"nehalem\" )";
  preInstall = ''
    # for f in $(ls lib); do
    #   chapelFixupBinary lib/$f
    # done
    # for f in $(ls bin); do
    #   chapelFixupBinary bin/$f
    # done
  '';
  dontStrip = true;

  buildInputs = [ lattice-symmetries.kernels2 lattice-symmetries.haskell.lib ];
  nativeBuildInputs = [
    (chapel.override {
      compiler = "llvm";
      settings = {
        CHPL_LIB_PIC = "pic";
        CHPL_UNWIND = "system";
      };
    })
    chapelFixupBinary
  ];
}
