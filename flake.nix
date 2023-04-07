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
      pkgs = import inputs.nixpkgs { inherit system; };
      haskellPackages = pkgs.haskellPackages;
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

      haskell-sources = inputs.nix-filter.lib {
        root = ./.;
        include = [
          "app"
          "cbits"
          "lib"
          "src"
          "test"
          "lattice-symmetries-haskell.cabal"
          "cabal.project"
          "LICENSE"
          "README.md"
        ];
      };

      lattice-symmetries-haskell =
        (haskellPackages.callCabal2nix "lattice-symmetries-haskell" haskell-sources {
          kernels = lattice-symmetries-kernels;
        }).overrideAttrs (attrs: {
          postInstall = ''
            ln --symbolic \
              $out/lib/ghc-${haskellPackages.ghc.version}/liblattice_symmetries_haskell.so \
              $out/lib/
          '';
        });

      lattice-symmetries-chapel = pkgs.stdenv.mkDerivation {
        pname = "lattice-symmetries-chapel";
        inherit version;
        src = ./chapel;

        dontConfigure = true;
        makeFlags = [
          "PREFIX=$(out)"
          "OPTIMIZATION=--fast"
          "CHPL_CFLAGS='-I${lattice-symmetries-kernels}/include'"
          "CHPL_LDFLAGS='-L${lattice-symmetries-haskell}/lib'"
        ];

        buildInputs = [
          lattice-symmetries-kernels
          lattice-symmetries-haskell
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
          lattice-symmetries-haskell
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
        haskell = lattice-symmetries-haskell;
        chapel = lattice-symmetries-chapel;
        python = lattice-symmetries-python;
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
          lattice-symmetries-haskell
        ];
        nativeBuildInputs = [
          chapel
        ];
        shellHook = ''
          export CHPL_CFLAGS="-I ${lattice-symmetries-kernels}/include";
          export CHPL_LDFLAGS="-L ${lattice-symmetries-haskell}/lib/ghc-${haskellPackages.ghc.version}";
        '';
      };
      devShells.default = haskellPackages.shellFor {
        packages = ps: [ lattice-symmetries-haskell ];
        withHoogle = true;
        nativeBuildInputs = with pkgs; with haskellPackages; [
          cabal-install
          hsc2hs
          haskell-language-server
          nil
          fourmolu
          cabal-fmt
          nixpkgs-fmt
        ];
      };
    }
  );
}
