{ version
, target
, chapel
, chapelFixupBinary
, stdenv
, lattice-symmetries
, hdf5
, halide
, singularity-tools
}:

let
  drv = stdenv.mkDerivation {
    pname = target;
    inherit version;
    src = ./.;
    configurePhase = "ln --symbolic ${lattice-symmetries.ffi} src/FFI.chpl";

    buildPhase = ''
      make \
        CHPL_CFLAGS='-I${lattice-symmetries.kernels_v2}/include --fast --no-ieee-float' \
        CHPL_LDFLAGS='-L${lattice-symmetries.haskell.lib}/lib' \
        HALIDE_PATH='${halide}' \
        bin/${target}

      for f in $(ls bin); do
        chapelFixupBinary bin/$f
      done
    '';

    installPhase = ''
      mkdir -p $out/bin
      install -Dm 755 bin/${target} $out/bin
    '';

    buildInputs = [
      lattice-symmetries.kernels_v2
      lattice-symmetries.haskell
    ];
    nativeBuildInputs = [ chapel chapelFixupBinary ];
  };
in
singularity-tools.buildImage {
  name = drv.pname;
  contents = [ drv ];
  diskSize = 10240;
  memSize = 5120;
}
