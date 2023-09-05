{ version
}:

final: prev: {
  atomic_queue = final.callPackage ./atomic_queue.nix { };
  lattice-symmetries = (prev.lattice-symmetries or { }) // {
    ffi = final.callPackage ./ffi.nix { inherit version; };
    chapel = final.callPackage ./. { inherit version; };
    test-data = final.stdenv.mkDerivation {
      pname = "lattice-symmetries-test-data";
      inherit version;
      src = final.fetchurl {
        url = "https://surfdrive.surf.nl/files/index.php/s/OK5527Awfgl1hT2/download?path=%2Fdata%2Fmatvec ";
        hash = "sha256-D+SsAaqB/0rH/ZChYBixg+FXcGcp1cIypzVhOJZx5iI=";
      };
      unpackPhase = "unzip $src";
      dontConfigure = true;
      dontBuild = true;
      installPhase = ''
        mkdir -p $out/share/data
        cp -r matvec $out/share/data/
      '';
      nativeBuildInputs = with final; [ unzip ];
    };
    distributed =
      let
        chapelBuild = target: makeFlags:
          let
            finalMakeFlags =
              final.lib.concatStringsSep " " (makeFlags ++
                [
                  "CHPL_HOST_MEM=jemalloc"
                  "CHPL_TARGET_MEM=jemalloc"
                  "CHPL_LAUNCHER=none"
                  "OPTIMIZATION=--fast"
                  "CHPL_CFLAGS='-I${final.lattice-symmetries.kernels}/include'"
                  "CHPL_LDFLAGS='-L${final.lattice-symmetries.haskell.lib}/lib'"
                  "HDF5_CFLAGS='-I${final.hdf5.dev}/include'"
                  "HDF5_LDFLAGS='-L${final.hdf5}/lib -lhdf5_hl -lhdf5 -lrt'"
                ]);
          in
          final.stdenv.mkDerivation {
            pname = target;
            inherit version;
            src = ./.;
            configurePhase = "ln --symbolic ${final.lattice-symmetries.ffi} src/FFI.chpl";

            passthru.finalMakeFlags = finalMakeFlags;
            buildPhase = ''
              make ${finalMakeFlags} \
                   bin/${target}

              for f in $(ls bin); do
                chapelFixupBinary bin/$f
              done
            '';
            doCheck = target == "TestMatrixVectorProduct";
            checkPhase = ''
              make ${finalMakeFlags} \
                   TEST_DATA=${final.lattice-symmetries.test-data} \
                   CHPL_ARGS='--numLocales=1' \
                   check-matrix-vector-product
            '';

            installPhase = ''
              mkdir -p $out/bin
              install -Dm 755 bin/${target} $out/bin
            '';
            nativeBuildInputs = with final; [ chapel chapelFixupBinary ];
          };

        toContainer = drv: final.singularity-tools.buildImage {
          name = drv.pname;
          contents = [ drv ];
          diskSize = 10240;
          memSize = 5120;
        };

        buildContainers = makeFlags:
          builtins.foldl'
            (acc: s: acc // { "${s}" = toContainer (chapelBuild s makeFlags); })
            { }
            [
              "TestMatrixVectorProduct"
              "BenchmarkStatesEnumeration"
              "BenchmarkMatrixVectorProduct"
              "BenchmarkBlockHashed"
            ];
      in
      {
        smp = buildContainers [
          "CHPL_COMM=gasnet"
          "CHPL_COMM_SUBSTRATE=smp"
        ];
        ibv = buildContainers [
          "CHPL_COMM=gasnet"
          "CHPL_COMM_SUBSTRATE=ibv"
          "CHPL_GASNET_SEGMENT=fast"
        ];
      };
  };
}
