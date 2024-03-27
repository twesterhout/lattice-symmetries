{ version
, enableDebugging ? false
, enableSanitizers ? false
}:

final: prev: {
  lattice-symmetries = (prev.lattice-symmetries or { }) // {

    ffi = final.callPackage ./ffi.nix { inherit version; };

    chapel = final.callPackage ./. {
      inherit version enableDebugging enableSanitizers;
      chapel = final.chapel.override {
        compiler = if enableSanitizers then "gnu" else "llvm";
        settings = {
          CHPL_GMP = "none";
          CHPL_RE2 = "none";
          CHPL_UNWIND = "none";
          CHPL_LIB_PIC = "pic";
          CHPL_TARGET_CPU = "none";
        } // final.lib.optionalAttrs enableSanitizers {
          CHPL_TARGET_MEM = "cstdlib";
          CHPL_HOST_MEM = "cstdlib";
          CHPL_TASKS = "fifo";
          CHPL_SANITIZE_EXE = "address";
        };
      };
    };

    # BenchmarkSingleLocaleMatrixVector = final.callPackage ./matrix-vector-product.nix {
    #   inherit version;
    #   target = "BenchmarkSingleLocaleMatrixVector";
    #   chapel = final.chapel.override {
    #     settings = { "CHPL_TARGET_CPU" = "znver3"; };
    #   };
    # };

    # test-data = final.stdenv.mkDerivation {
    #   pname = "lattice-symmetries-test-data";
    #   inherit version;
    #   src = final.fetchurl {
    #     url = "https://surfdrive.surf.nl/files/index.php/s/jyaDzuXyTPI1HUs/download";
    #     hash = "sha256-BUGgD8PYOW5s039PkuATeyRIvSGEIn3PGRcDNE1j5hs=";
    #     # url = "https://github.com/twesterhout/lattice-symmetries/releases/download/test-v1/matvec.zip";
    #     # hash = "sha256-DR+HpYZkgxigXooXCkluGZz7ChRy8xqCKeVd4B3nGDQ=";
    #     # https://surfdrive.surf.nl/files/index.php/s/OK5527Awfgl1hT2/download?path=%2Fdata%2Fv3%2Fmatvec";
    #     # sha256-laaL7WemYccjGU9B0fPN/SzyBv+vOEdIy7BgJ2JzKRw=";
    #   };
    #   unpackPhase = "unzip $src";
    #   dontConfigure = true;
    #   dontBuild = true;
    #   installPhase = ''
    #     mkdir -p $out/share/data/matvec
    #     cp v4/*.{yaml,h5} $out/share/data/matvec/
    #   '';
    #   nativeBuildInputs = with final; [ unzip ];
    # };
    # distributed =
    #   let
    #     chapelBuild = target: makeFlags:
    #       let
    #         finalMakeFlags =
    #           final.lib.concatStringsSep " " (makeFlags ++
    #             [
    #               "CHPL_HOST_MEM=jemalloc"
    #               "CHPL_TARGET_MEM=jemalloc"
    #               "CHPL_LAUNCHER=none"
    #               "OPTIMIZATION='--debug -g'"
    #               "CHPL_CFLAGS='-I${final.lattice-symmetries.kernels}/include'"
    #               "CHPL_LDFLAGS='-L${final.lattice-symmetries.haskell.lib}/lib'"
    #               "HDF5_CFLAGS='-I${final.hdf5.dev}/include'"
    #               "HDF5_LDFLAGS='-L${final.hdf5}/lib -lhdf5_hl -lhdf5 -lrt'"
    #             ]);
    #       in
    #       final.stdenv.mkDerivation {
    #         pname = target;
    #         inherit version;
    #         src = ./.;
    #         configurePhase = "ln --symbolic ${final.lattice-symmetries.ffi} src/FFI.chpl";

    #         passthru.finalMakeFlags = finalMakeFlags;
    #         buildPhase = ''
    #           make ${finalMakeFlags} \
    #                bin/${target}

    #           for f in $(ls bin); do
    #             chapelFixupBinary bin/$f
    #           done
    #         '';
    #         doCheck = target == "TestMatrixVectorProduct";
    #         checkPhase = ''
    #           export GASNET_PSHM_NODES=1
    #           make ${finalMakeFlags} \
    #                TEST_DATA=${final.lattice-symmetries.test-data}/share \
    #                CHPL_ARGS='--numLocales=1' \
    #                check-matrix-vector-product
    #         '';

    #         installPhase = ''
    #           mkdir -p $out/bin
    #           install -Dm 755 bin/${target} $out/bin
    #         '';
    #         nativeBuildInputs = with final; [ chapel chapelFixupBinary ];
    #       };

    #     toContainer = drv: final.singularity-tools.buildImage {
    #       name = drv.pname;
    #       contents = [ drv ];
    #       diskSize = 10240;
    #       memSize = 5120;
    #     };

    #     buildContainers = makeFlags:
    #       builtins.foldl'
    #         (acc: s: acc // { "${s}" = toContainer (chapelBuild s makeFlags); })
    #         { }
    #         [
    #           "TestMatrixVectorProduct"
    #           "BenchmarkStatesEnumeration"
    #           "BenchmarkMatrixVectorProduct"
    #           "BenchmarkBlockHashed"
    #         ];
    #   in
    #   {
    #     smp = buildContainers [
    #       "CHPL_COMM=gasnet"
    #       "CHPL_COMM_SUBSTRATE=smp"
    #       "CHPL_TARGET_CPU=nehalem"
    #     ];
    #     ibv = buildContainers [
    #       "CHPL_COMM=gasnet"
    #       "CHPL_COMM_SUBSTRATE=ibv"
    #       "CHPL_GASNET_SEGMENT=fast"
    #       "CHPL_TARGET_CPU=nehalem"
    #     ];
    #     test = chapelBuild "TestMatrixVectorProduct" [
    #       "CHPL_COMM=gasnet"
    #       "CHPL_COMM_SUBSTRATE=smp"
    #       "CHPL_TARGET_CPU=nehalem"
    #     ];
    #   };
  };
}
