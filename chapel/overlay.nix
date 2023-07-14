{ version
}:

final: prev: {
  lattice-symmetries = (prev.lattice-symmetries or { }) // {
    ffi = final.callPackage ./ffi.nix { inherit version; };
    chapel = final.callPackage ./. { inherit version; };
    distributed =
      let
        chapelBuild = target: final.stdenv.mkDerivation {
          pname = target;
          inherit version;
          src = ./.;
          configurePhase = "ln --symbolic ${final.lattice-symmetries.ffi} src/FFI.chpl";
          buildPhase = ''
            # CHPL_HOST_MEM=cstdlib CHPL_TARGET_MEM=cstdlib
            make \
              CHPL_COMM=gasnet \
              CHPL_COMM_SUBSTRATE=ibv \
              CHPL_LAUNCHER=none \
              OPTIMIZATION=--fast \
              CHPL_CFLAGS='-I${final.lattice-symmetries.kernels}/include' \
              CHPL_LDFLAGS='-L${final.lattice-symmetries.haskell.lib}/lib' \
              bin/${target}
  
            for f in $(ls bin); do
              chapelFixupBinary bin/$f
            done
          '';
          installPhase = ''
            mkdir -p $out/bin
            install -Dm 755 bin/* $out/bin
          '';
          nativeBuildInputs = with final; [ chapel chapelFixupBinary ];
        };

        toContainer = drv: final.singularity-tools.buildImage {
          name = drv.pname;
          contents = [ drv ];
          diskSize = 10240;
          memSize = 5120;
        };
      in
        {
          test-matrix-vector-product = toContainer (chapelBuild "TestMatrixVectorProduct");
          benchmark-states-enumeration = toContainer (chapelBuild "BenchmarkStatesEnumeration");
          benchmark-matrix-vector-product = toContainer (chapelBuild "BenchmarkMatrixVectorProduct");
        };
  };
}
