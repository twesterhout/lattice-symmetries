{ version
}:

final: prev: {
  chapel = prev.chapel.overrideAttrs (attrs: {
    src = prev.fetchFromGitHub {
      owner = "chapel-lang";
      repo = "chapel";
      # BROKEN
      # rev = "e3a9c913516ac9abf48c9a8b86c199953f12030f";
      # hash = "sha256-MzCIzJdFAjK/BNx6C6gaF/3Y9lmw08CauVJfu6N+YrE=";
      # BROKEN
      # rev = "17beb276d68be9a373550c54c148d370f4e84109";
      # hash = "sha256-Ds1ff1U/D+1ihhOUr5GMdPlHuL7TARgoysAb+Fyxyqw=";

      # BROKEN Aug 9
      # rev = "1bb9a6146ed204843e1fcc381ed0d27bb272b7b0";
      # hash = "sha256-I/rdz4R7fI2kgSw7qlIxr/pvw+M+lXFuJie6U71tXCk=";

      # BROKEN Aug 8
      rev = "57e4a64aa1de4b8ee6648ab39abba9492783be24";
      hash = "sha256-Ql22LRIWgfJua6piJztjTiGRBS/SQr1otBHQePjhqtU=";

      # Fails to compile :/
      # rev = "d0f526b9685d3c651d52f87723b7b23d83c9817a";
      # hash = "sha256-bI8JQMxwI9S41AMpF72wzuAyFwxQ2/DKlfhr9tiZ0z4=";

      # WORKS Aug 9
      # rev = "038cdee419658caffd906cbcf755a230cd8702e1";
      # hash = "sha256-l/TGfsAC6ujF/0gD5g1FInu1+ItUs2mhuOsZ/ooUxlw=";

      # WORKS Aug 8
      # rev = "ad9b35af938fb4ee88149563a1c1f9853a6fbd8d";
      # hash = "sha256-OnH8Ob6IQo6nCf+5TemU3vWzCyoDs6XfArUPv8ixlNY=";

      # WORKS Aug 8
      # rev = "127dc6302e1d3e0b5900a2706101c08d61f889d1";
      # hash = "sha256-aoRXLfGLEbtaD4cdrYl3zOb36zLNQ9KW882cgN5oSD0=";

      # WORKS Aug 7
      # rev = "a4f05385b6e25fd1c24614a8376aa5f449bc62a6";
      # hash = "sha256-Muy14+5rOAoMCwyDIphY/QYF1fOyXl/hIhApjW3QIDQ=";

      # WORKS Aug 4
      # rev = "22f138fe25e66586f7b35699ff0699f7900c6a39";
      # hash = "sha256-phIws0/fN8xjjpgLggGZb/8OQLETjMdzI2RTpBbZGoE=";

      # WORKS July 17
      # rev = "f361e5dd2bda02c2dd1826a05b8b791a4d8909b2";
      # hash = "sha256-Lnj4umnx/Tv2+FD69rPj+NcTw4IXaDfjtF9LskwSuAg=";

      # WORKS
      # rev = "4585257c03c11e4d9aff16ff395f7217f5c162b7";
      # hash = "sha256-cDsdypLlh2YFym6rb0fiaX2ZW16By00HYrow2jDpKH0=";

      # WORKS
      # owner = "bradcray";
      # repo = "chapel";
      # rev = "15fb5d5b2e7d6de146b8c39aa69fb00ce11e5f17";
      # hash = "sha256-FNV2INfzpSsa1a3C+Su3Xi0e8Pe9TCaUhwHoQnqt9XE=";
    };
  });

  atomic_queue = final.callPackage ./atomic_queue.nix { };
  lattice-symmetries = (prev.lattice-symmetries or { }) // {
    ffi = final.callPackage ./ffi.nix { inherit version; };
    chapel = final.callPackage ./. { inherit version; };
    test-data = final.stdenv.mkDerivation {
      pname = "lattice-symmetries-test-data";
      inherit version;
      src = final.fetchzip {
        url = "https://github.com/twesterhout/lattice-symmetries/releases/download/test-v1/matvec.zip";
        # https://surfdrive.surf.nl/files/index.php/s/OK5527Awfgl1hT2/download?path=%2Fdata%2Fv3%2Fmatvec";
        hash = "sha256-DR+HpYZkgxigXooXCkluGZz7ChRy8xqCKeVd4B3nGDQ=";
        # sha256-laaL7WemYccjGU9B0fPN/SzyBv+vOEdIy7BgJ2JzKRw=";
      };
      # unpackPhase = "unzip $src";
      dontConfigure = true;
      dontBuild = true;
      installPhase = ''
        mkdir -p $out/share/data/matvec
        cp *.h5 $out/share/data/matvec/
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
              export GASNET_PSHM_NODES=1
              make ${finalMakeFlags} \
                   TEST_DATA=${final.lattice-symmetries.test-data}/share \
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
        test = chapelBuild "TestMatrixVectorProduct" [ "CHPL_COMM=gasnet" "CHPL_COMM_SUBSTRATE=smp" ];
      };
  };
}
