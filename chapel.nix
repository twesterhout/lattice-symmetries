{ version
}:

final: prev: {
  lattice-symmetries-chapel-ffi = final.stdenv.mkDerivation
    {
      pname = "lattice-symmetries-chapel-ffi";
      inherit version;
      src = ./include;

      buildPhase = ''
        c2chapel \
          lattice_symmetries.h \
          -I${final.chapel}/runtime/include \
          >> FFI.chpl

        # sed -i 's/extern type ls_hs_scalar = _Complex double/extern type ls_hs_scalar = complex(128)/' FFI.chpl
      '';

      installPhase = "install -m 644 FFI.chpl $out";
      buildInputs = [ ];
      nativeBuildInputs = [ final.chapel ];
    };

  lattice-symmetries-chapel = final.stdenv.mkDerivation {
    pname = "lattice-symmetries-chapel";
    inherit version;
    src = final.lib.fileset.toSource {
      root = ./.;
      fileset = final.lib.fileset.unions [
        ./include
        ./chapel
      ];
    };

    configurePhase = ''
      ln --symbolic ${final.lattice-symmetries-chapel-ffi} chapel/src/FFI.chpl;
      cat chapel/src/FFI.chpl
    '';

    preBuild = ''
      cd chapel/

      makeFlagsArray+=(
        PREFIX="$out"
        CHPL_CFLAGS="-I../include --no-ieee-float --local"
      )
    '';

    preInstall = ''
      patchelf --debug --remove-needed libstdc++.so.6 ./lib/liblattice_symmetries_chapel.so
      patchelf --debug --remove-needed libgcc_s.so.1 ./lib/liblattice_symmetries_chapel.so
      # for f in $(ls lib); do
      #   chapelFixupBinary lib/$f
      #   # These are only needed during compilation
      #   remove-references-to -t $${lattice-symmetries.kernels_v2} lib/$f
      # done
    '';

    buildInputs = [
      final.halide
    ];

    nativeBuildInputs = [
      final.patchelf
      final.chapel

      # removeReferencesTo
    ];
  };

}
