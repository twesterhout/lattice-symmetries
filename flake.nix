{
  description = "twesterhout/lattice-symmetries";

  nixConfig = {
    extra-substituters = "https://twesterhout-chapel.cachix.org";
    extra-trusted-public-keys = "twesterhout-chapel.cachix.org-1:bs5PQPqy21+rP2KJl+O40/eFVzdsTe6m7ZTiOEE7PaI=";
  };

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
    nix-chapel = {
      url = "github:twesterhout/nix-chapel";
      inputs.nixpkgs.follows = "nixpkgs";
      inputs.flake-utils.follows = "flake-utils";
    };
  };

  outputs = { self, nixpkgs, flake-utils, nix-chapel }:
    let
      inherit (nixpkgs) lib;
      version = "2.2.0";

      kernels-overlay = import ./kernels/overlay.nix { inherit version; };
      haskell-overlay = { withPic }: import ./haskell/overlay.nix { inherit lib withPic; };
      chapel-overlay = import ./chapel/overlay.nix { inherit version; };
      python-overlay = import ./python/overlay.nix { inherit version; };

      composed-overlay = { withPic }: lib.foldl' lib.composeExtensions (_: _: { }) ([
        nix-chapel.overlays.default
        kernels-overlay
        (haskell-overlay { inherit withPic; })
        chapel-overlay
        python-overlay
      ]
      ++ lib.optionals withPic [
        # An overlay to replace ghc96 with a custom one that has
        # the static RTS libraries compiled with -fPIC. This lets us use
        # these static libraries to build a self-contained shared library.
        (final: prev:
          let
            ourGhc =
              if prev.stdenv.isLinux then
                prev.haskell.compiler.ghc962.override { enableRelocatedStaticLibs = true; }
              else
                prev.haskell.compiler.ghc962;
          in
          lib.recursiveUpdate prev {
            haskell.packages.ghc962.ghc = ourGhc;
            haskell.compiler.ghc962 = ourGhc;
          })
      ]);

      pkgs-for = args: system: import nixpkgs {
        inherit system;
        overlays = [ (composed-overlay args) ];
      };

    in
    {
      overlays.default = composed-overlay { withPic = true; };

      templates.default = {
        path = builtins.toPath "${./.}/template";
        description = "Python project template that uses lattice-symmetries";
      };

      packages = flake-utils.lib.eachDefaultSystemMap (system:
        with (pkgs-for { withPic = true; } system); {
          inherit (lattice-symmetries) kernels haskell chapel python distributed;
          inherit atomic_queue;
        });

      devShells = flake-utils.lib.eachDefaultSystemMap (system:
        let
          pkgs = pkgs-for { withPic = true; } system;
          pkgsNoPic = pkgs-for { withPic = false; } system;
        in
        {
          default = self.outputs.devShells.${system}.haskell;
          kernels = with pkgs; mkShell {
            buildInputs = [ halide ];
            nativeBuildInputs = [ gnumake gcc ];
            shellHook = "export HALIDE_PATH=${halide}";
          };
          haskell = with pkgsNoPic; haskellPackages.shellFor {
            packages = ps: [ ps.lattice-symmetries-haskell ];
            withHoogle = true;
            nativeBuildInputs = with haskellPackages; [
              cabal-fmt
              cabal-install
              fourmolu
              haskell-language-server
              hsc2hs
              nil
              nixpkgs-fmt
              (python3Packages.grip.overrideAttrs (attrs: {
                src = fetchFromGitHub {
                  owner = "Antonio-R1";
                  repo = "grip";
                  rev = "d2efd3c6a896c01cfd7624b6504107e7b3b4b20f";
                  hash = "sha256-0wgIM7Ll5WELvAOiu1TLyoNSrhJ22Y1SRbWqa3BDF3k=";
                };
                checkPhase = "true";
                installCheckPhase = "true";
              }))
            ];
            shellHook = ''
              if [ ! -f libkernels.so ]; then
                gcc -shared -o libkernels.so ${lattice-symmetries.kernels}/lib/libkernels.a
              fi
              export LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH;
            '';
          };
          chapel = with pkgs; mkShell {
            buildInputs = [
              lattice-symmetries.kernels
              lattice-symmetries.haskell
              hdf5
              hdf5.dev
            ];
            nativeBuildInputs = [
              pkgs.chapel
              chapelFixupBinary
              gcc
              pkg-config
              prettierd
            ];
            shellHook = ''
              export LS_KERNELS="${lattice-symmetries.kernels}";
              export LS_HASKELL="${lattice-symmetries.haskell}";
              export LS_HASKELL_LIB="${lattice-symmetries.haskell.lib}";
              export CHPL_CFLAGS="-I${lattice-symmetries.kernels}/include"
              export CHPL_LDFLAGS="-L${lattice-symmetries.haskell.lib}/lib"
              export HDF5_CFLAGS="-I${hdf5.dev}/include"
              export HDF5_LDFLAGS="-L${hdf5}/lib -lhdf5_hl -lhdf5 -lrt"
            '';
          };
          python = with pkgs; lattice-symmetries.python.overrideAttrs (attrs: {
            nativeBuildInputs = (attrs.nativeBuildInputs or [ ]) ++ [
              python3Packages.black
              nodePackages.pyright
            ];
          });
        });
    };
}
