{
  description = "twesterhout/lattice-symmetries-haskell";

  nixConfig = {
    extra-experimental-features = "nix-command flakes";
    extra-substituters = "https://twesterhout-chapel.cachix.org";
    extra-trusted-public-keys = "twesterhout-chapel.cachix.org-1:bs5PQPqy21+rP2KJl+O40/eFVzdsTe6m7ZTiOEE7PaI=";
  };

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
    nix-filter.url = "github:numtide/nix-filter";
    flake-compat = {
      url = "github:edolstra/flake-compat";
      flake = false;
    };
    nix-chapel = {
      url = "github:twesterhout/nix-chapel";
      inputs.nixpkgs.follows = "nixpkgs";
      inputs.flake-utils.follows = "flake-utils";
    };
  };

  outputs = inputs: inputs.flake-utils.lib.eachDefaultSystem (system:
    with builtins;
    let
      version = "2.1.0";

      inherit (inputs.nixpkgs) lib;
      pkgs = import inputs.nixpkgs {
        inherit system;
        overlays = [
          # An overlay to replace ghc961 with a custom one that has
          # the static RTS libraries compiled with -fPIC. This lets us use
          # these static libraries to build a self-contained shared library.
          (final: prev:
            let
              ourGhc = prev.haskell.compiler.ghc961.override {
                enableRelocatedStaticLibs = true;
              };
            in
            lib.recursiveUpdate prev {
              haskell.packages.ghc961.ghc = ourGhc;
              haskell.compiler.ghc961 = ourGhc;
            })
        ];
      };

      haskell-sources = root: inputs.nix-filter.lib {
        inherit root;
        include = [
          "app"
          "cbits"
          "lib"
          "src"
          "test"
          (inputs.nix-filter.lib.matchExt "cabal")
          "cabal.project"
          "LICENSE"
          "README.md"
        ];
      };

      overrideHaskellPackages = hp: withPic:
        hp.override {
          overrides = self: super: {
            # Our project
            lattice-symmetries-haskell =
              (self.callCabal2nix "lattice-symmetries-haskell" ./haskell {
                kernels = lattice-symmetries-kernels;
              }).overrideAttrs
                (attrs: {
                  outputs = (attrs.outputs or [ ]) ++ [ "lib" ];
                  postInstall = ''
                    ${attrs.postInstall or ""}

                    echo "Installing foreign library to $lib ..."
                    mkdir -p $lib/lib
                    install -v -Dm 755 \
                      $out/lib/ghc-*/lib/liblattice_symmetries_haskell.so \
                      $lib/lib/

                    mkdir -p $out/include
                    install -v -Dm 644 \
                      lattice_symmetries_functions.h \
                      $out/include/
                  '';
                });
          } // lib.optionalAttrs withPic {
            # Ensure that all Haskell packages are built with -fPIC
            mkDerivation = args: (super.mkDerivation args).overrideAttrs (attrs: {
              configureFlags = (attrs.configureFlags or [ ]) ++ [
                "--ghc-option=-fPIC"
                "--ghc-option=-fexternal-dynamic-refs"
              ];
            });
          } // lib.optionalAttrs (lib.versionAtLeast hp.ghc.version "9.6.1") {
            # Loosen constraints to make them build with 9.6.1
            tagged = pkgs.haskell.lib.doJailbreak super.tagged;
            zigzag = pkgs.haskell.lib.doJailbreak super.zigzag;
            typerep-map = pkgs.haskell.lib.doJailbreak super.typerep-map;
            relude = pkgs.haskell.lib.dontCheck (pkgs.haskell.lib.doJailbreak super.relude);
            bytebuild = pkgs.haskell.lib.doJailbreak super.bytebuild;
            chronos = pkgs.haskell.lib.doJailbreak super.chronos;
          };
        };

      chapel = inputs.nix-chapel.packages.${system}.chapel;
      chapelFixupBinary = inputs.nix-chapel.packages.${system}.chapelFixupBinary;

      lattice-symmetries-kernels = pkgs.stdenv.mkDerivation {
        pname = "lattice-symmetries-kernels";
        inherit version;
        src = ./kernels;

        dontConfigure = true;
        makeFlags = [
          "HALIDE_PATH=${pkgs.halide}"
          "PREFIX=$(out)"
        ];
      };

      lattice-symmetries-haskell =
        let
          hp = overrideHaskellPackages pkgs.haskell.packages.ghc961 true;
        in
        hp.lattice-symmetries-haskell;

      lattice-symmetries-chapel-ffi = pkgs.stdenv.mkDerivation {
        pname = "lattice-symmetries-chapel-ffi";
        inherit version;
        unpackPhase = "true";
        buildPhase = ''
          c2chapel \
            ${lattice-symmetries-haskell}/include/lattice_symmetries_functions.h \
            -DLS_C2CHAPEL \
            -I${chapel}/runtime/include \
            -I${lattice-symmetries-kernels}/include \
            >FFI.chpl

          # Remove the declaration of chpl_external_array since it's already
          # present in the ExternalArray module
          sed -i -e '/extern record chpl_external_array/,+5d' FFI.chpl
          sed -i -E '/chpl_make_external_array(_ptr)?(_free)?\(/d' FFI.chpl
          sed -i -e '/cleanupOpaqueArray(/d' FFI.chpl
          sed -i -e '/chpl_free_external_array(/d' FFI.chpl
          sed -i -e '/chpl_call_free_func(/d' FFI.chpl
          sed -i 's/extern type ls_hs_scalar = _Complex double/extern type ls_hs_scalar = complex(128)/' FFI.chpl
        '';
        installPhase = ''
          install -m 644 FFI.chpl $out
        '';
        buildInputs = [
          lattice-symmetries-haskell
          lattice-symmetries-kernels
        ];
        nativeBuildInputs = [
          chapel
        ];
      };


      lattice-symmetries-chapel = pkgs.stdenv.mkDerivation {
        pname = "lattice-symmetries-chapel";
        inherit version;
        src = ./chapel;

        configurePhase = ''
          ln --symbolic ${lattice-symmetries-chapel-ffi} src/FFI.chpl
        '';

        makeFlags = [
          "PREFIX=$(out)"
          # "OPTIMIZATION=--fast"
          "CHPL_CFLAGS='-I${lattice-symmetries-kernels}/include'"
          "CHPL_LDFLAGS='-L${lattice-symmetries-haskell.lib}/lib'"
        ];

        buildInputs = [ lattice-symmetries-kernels lattice-symmetries-haskell.lib ];
        nativeBuildInputs = [ chapel chapelFixupBinary ];
      };

      test-matrix-vector = pkgs.stdenv.mkDerivation {
        pname = "test-matrix-vector";
        inherit version;
        src = ./chapel;

        configurePhase = ''
          ln --symbolic ${lattice-symmetries-chapel-ffi} src/FFI.chpl
        '';

        buildPhase = ''
          make \
            CHPL_COMM=gasnet \
            CHPL_COMM_SUBSTRATE=smp \
            OPTIMIZATION=--fast \
            CHPL_CFLAGS='-I${lattice-symmetries-kernels}/include' \
            CHPL_LDFLAGS='-L${lattice-symmetries-haskell.lib}/lib' \
            HDF5_CFLAGS='-I${pkgs.hdf5.dev}/include' \
            HDF5_LDFLAGS='-L${pkgs.hdf5}/lib -lhdf5_hl -lhdf5 -lrt' \
            bin/TestMatrixVectorProduct

          chapelFixupBinary bin/TestMatrixVectorProduct
          chapelFixupBinary bin/TestMatrixVectorProduct_real
        '';

        installPhase = ''
          mkdir -p $out/bin
          install -Dm 755 bin/TestMatrixVectorProduct* $out/bin
        '';

        nativeBuildInputs = [ chapel chapelFixupBinary ];
      };

      lattice-symmetries-python = pkgs.python3Packages.buildPythonPackage {
        pname = "lattice-symmetries";
        inherit version;
        src = ./python;

        buildInputs = [
          lattice-symmetries-kernels
          lattice-symmetries-haskell
          lattice-symmetries-chapel
        ];
        propagatedBuildInputs = with pkgs.python3Packages; [
          cffi
          loguru
          numpy
          scipy
        ];

        postPatch = ''
          awk '/python-cffi: START/{flag=1;next}/python-cffi: STOP/{flag=0}flag' \
            ${lattice-symmetries-kernels}/include/lattice_symmetries_types.h \
            >lattice_symmetries/extracted_declarations.h
          awk '/python-cffi: START/{flag=1;next}/python-cffi: STOP/{flag=0}flag' \
            ${lattice-symmetries-haskell}/include/lattice_symmetries_functions.h \
            >>lattice_symmetries/extracted_declarations.h
        '';

        installCheckPhase = ''
          # we want to import the installed module that also contains the compiled library
          rm -rf lattice_symmetries
          runHook pytestCheckPhase
        '';
        nativeCheckInputs = with pkgs.python3Packages; [
          pytestCheckHook
        ];
      };
    in
    {
      packages = {
        kernels = lattice-symmetries-kernels;
        haskell = lattice-symmetries-haskell;
        chapel = lattice-symmetries-chapel;
        test-matrix-vector = test-matrix-vector;
        python = lattice-symmetries-python;
        ghc = ghc;
      };
      devShells.kernels = pkgs.mkShell {
        buildInputs = with pkgs; [
          halide
        ];
        nativeBuildInputs = with pkgs; [
          gnumake
          gcc
        ];
        shellHook = ''
          export HALIDE_PATH=${pkgs.halide}
        '';
      };
      devShells.chapel = with pkgs; mkShell {
        # packages = [ lattice-symmetries-chapel ];
        buildInputs = [
          lattice-symmetries-kernels
          lattice-symmetries-haskell
          lattice-symmetries-haskell.lib
          hdf5
          hdf5.dev
        ];
        nativeBuildInputs = [
          chapel
          gcc
          pkg-config
        ];
        shellHook = ''
          export LS_KERNELS="${lattice-symmetries-kernels}";
          export LS_HASKELL="${lattice-symmetries-haskell}";
          export LS_HASKELL_LIB="${lattice-symmetries-haskell.lib}";
          export CHPL_CFLAGS="-I${lattice-symmetries-kernels}/include"
          export CHPL_LDFLAGS="-L${lattice-symmetries-haskell.lib}/lib"
          export HDF5_CFLAGS="-I${hdf5.dev}/include"
          export HDF5_LDFLAGS="-L${hdf5}/lib -lhdf5_hl -lhdf5 -lrt"
        '';
      };
      devShells.default =
        let
          hp = overrideHaskellPackages pkgs.haskellPackages false;
        in
        hp.shellFor {
          packages = ps: [ ps.lattice-symmetries-haskell ];
          withHoogle = true;
          nativeBuildInputs = with pkgs; with hp; [
            cabal-fmt
            cabal-install
            fourmolu
            haskell-language-server
            hsc2hs
            nil
            nixpkgs-fmt
          ];
          shellHook = ''
            export LD_LIBRARY_PATH=${lattice-symmetries-kernels}/lib:$LD_LIBRARY_PATH;
          '';
        };
    }
  );
}
