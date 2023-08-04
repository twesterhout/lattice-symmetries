{ version
}:

final: prev: {
  atomic_queue = final.callPackage ./atomic_queue.nix { };
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
            #  CHPL_GASNET_SEGMENT=everything
            #  CHPL_HOST_MEM=cstdlib CHPL_TARGET_MEM=cstdlib

            make \
              CHPL_COMM=gasnet \
              CHPL_COMM_SUBSTRATE=ibv \
              CHPL_GASNET_SEGMENT=fast \
              CHPL_HOST_MEM=jemalloc CHPL_TARGET_MEM=jemalloc \
              CHPL_LAUNCHER=none \
              OPTIMIZATION=--fast \
              CHPL_CFLAGS='-I${final.lattice-symmetries.kernels}/include' \
              CHPL_LDFLAGS='-L${final.lattice-symmetries.haskell.lib}/lib' \
              HDF5_CFLAGS='-I${final.hdf5.dev}/include' \
              HDF5_LDFLAGS='-L${final.hdf5}/lib -lhdf5_hl -lhdf5 -lrt' \
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
        benchmark-block-hashed = toContainer (chapelBuild "BenchmarkBlockHashed");
      };
  };
}
