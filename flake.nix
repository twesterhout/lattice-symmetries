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
    nix-chapel.url = "github:twesterhout/nix-chapel";
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
      haskellPackages =
        pkgs.haskell.packages.ghc961.override {
          overrides = self: super: {
            # Ensure that all Haskell packages are built with -fPIC
            mkDerivation = args: (super.mkDerivation args).overrideAttrs (attrs: {
              configureFlags = (attrs.configureFlags or [ ]) ++ [
                "--ghc-option=-fPIC"
                "--ghc-option=-fexternal-dynamic-refs"
              ];
            });
            # Loosen constraints to make them build with 9.6.1
            tagged = pkgs.haskell.lib.doJailbreak super.tagged;
            zigzag = pkgs.haskell.lib.doJailbreak super.zigzag;
            typerep-map = pkgs.haskell.lib.doJailbreak super.typerep-map;
            relude = pkgs.haskell.lib.dontCheck (pkgs.haskell.lib.doJailbreak super.relude);
            bytebuild = pkgs.haskell.lib.doJailbreak super.bytebuild;
            chronos = pkgs.haskell.lib.doJailbreak super.chronos;
            # Our project
            lattice-symmetries-haskell =
              (self.callCabal2nix "lattice-symmetries-haskell" (haskell-sources ./.) {
                kernels = lattice-symmetries-kernels;
              }).overrideAttrs
                (attrs: {
                  postInstall = ''
                    ${attrs.postInstall or ""}
                    ln --symbolic \
                      $out/lib/ghc-*/lib/liblattice_symmetries_haskell.so \
                      $out/lib/
                  '';
                });
          };
        };
      # haskellPackages = ghc.override # pkgs.haskell.packages.ghc961.override
      #   {
      #     overrides = self: super: {
      #       lattice-symmetries-haskell =
      #         (self.callCabal2nix "lattice-symmetries-haskell" (haskell-sources ./.) {
      #           kernels = lattice-symmetries-kernels;
      #         }).overrideAttrs
      #           (attrs: {
      #             postInstall = ''
      #               ${attrs.postInstall or ""}
      #               ln --symbolic \
      #                 $out/lib/ghc-*/liblattice_symmetries_haskell.so \
      #                 $out/lib/
      #             '';
      #           });
      #     };
      #   };
      # (mapAttrs (name: value:
      #   if (value ? overrideAttrs)
      #   then
      #     value.overrideAttrs
      #       (attrs: {
      #         #
      #       })
      #   else value))
      # ];
      chapel = inputs.nix-chapel.packages.${system}.default;

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

      lattice-symmetries-chapel = pkgs.stdenv.mkDerivation {
        pname = "lattice-symmetries-chapel";
        inherit version;
        src = ./chapel;

        dontConfigure = true;
        makeFlags = [
          "PREFIX=$(out)"
          "OPTIMIZATION=--fast"
          "CHPL_CFLAGS='-I${lattice-symmetries-kernels}/include'"
          "CHPL_LDFLAGS='-L${haskellPackages.lattice-symmetries-haskell}/lib'"
        ];

        buildInputs = [
          lattice-symmetries-kernels
          haskellPackages.lattice-symmetries-haskell
        ];
        nativeBuildInputs = [
          chapel
        ];
      };

      lattice-symmetries-python = pkgs.python3Packages.buildPythonPackage {
        pname = "lattice-symmetries";
        inherit version;
        src = ./python;

        buildInputs = [
          lattice-symmetries-kernels
          haskellPackages.lattice-symmetries-haskell
          lattice-symmetries-chapel
        ];
        propagatedBuildInputs = with pkgs.python3Packages; [
          cffi
          loguru
          numpy
          scipy
        ];

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
        haskell = haskellPackages.lattice-symmetries-haskell;
        chapel = lattice-symmetries-chapel;
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
      devShells.chapel = pkgs.mkShell {
        buildInputs = [
          lattice-symmetries-kernels
          haskellPackages.lattice-symmetries-haskell
        ];
        nativeBuildInputs = [
          chapel
        ];
        shellHook = ''
          export CHPL_CFLAGS="-I ${lattice-symmetries-kernels}/include";
          export CHPL_LDFLAGS="-L ${haskellPackages.lattice-symmetries-haskell}/lib/ghc-${haskellPackages.ghc.version}";
        '';
      };
      devShells.default = haskellPackages.shellFor {
        packages = ps: [ ps.lattice-symmetries-haskell ];
        # withHoogle = true;
        nativeBuildInputs = with pkgs; with pkgs.haskell.packages.ghc92; [
          # cabal-fmt
          cabal-install
          # fourmolu
          # haskell-language-server
          hsc2hs
          nil
          nixpkgs-fmt
        ];
      };
    }
  );
}
