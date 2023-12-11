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
    haskell-python-tools = {
      url = "github:twesterhout/haskell-python-tools.nix";
      inputs.nixpkgs.follows = "nixpkgs";
    };
  };

  outputs = { self, nixpkgs, flake-utils, nix-chapel, haskell-python-tools }:
    let
      inherit (nixpkgs) lib;
      inherit (haskell-python-tools.lib)
        doInstallForeignLibs
        doEnableRelocatedStaticLibs;
      version = "2.2.0";

      kernels-overlay = import ./kernels/overlay.nix { inherit version; };
      haskell-overlay = import ./haskell/overlay.nix {
        inherit lib doInstallForeignLibs doEnableRelocatedStaticLibs;
      };
      chapel-overlay = import ./chapel/overlay.nix { inherit version; };
      python-overlay = import ./python/overlay.nix { inherit version; };

      composed-overlay = lib.composeManyExtensions [
        nix-chapel.overlays.default
        kernels-overlay
        haskell-overlay
        chapel-overlay
        python-overlay
      ];

      pkgs-for = system: import nixpkgs {
        inherit system;
        overlays = [ composed-overlay ];
      };

    in
    {
      overlays.default = composed-overlay;

      templates.default = {
        path = builtins.toPath "${./.}/template";
        description = "Python project template that uses lattice-symmetries";
      };

      packages = flake-utils.lib.eachDefaultSystemMap (system:
        with pkgs-for system; {
          inherit (lattice-symmetries) kernels haskell chapel python distributed test-data;
          inherit atomic_queue;
          inherit haskellPackages;
        });

      devShells = flake-utils.lib.eachDefaultSystemMap (system:
        let
          pkgs = pkgs-for system;
        in
        {
          default = self.outputs.devShells.${system}.haskell;
          kernels = with pkgs; mkShell {
            buildInputs = [ halide ];
            nativeBuildInputs = [ gnumake cmake clang clang-tools ];
            shellHook = "export HALIDE_PATH=${halide}";
          };
          haskell = with pkgs; haskellPackages.shellFor {
            packages = ps: [ ps.lattice-symmetries-haskell ];
            withHoogle = true;
            nativeBuildInputs = with haskellPackages; [
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
              export CHPL_COMM=gasnet
              export CHPL_COMM_SUBSTRATE=smp
              export CHPL_RT_OVERSUBSCRIBED=yes
              export CHPL_HOST_MEM=jemalloc
              export CHPL_TARGET_MEM=jemalloc
              export CHPL_LAUNCHER=none
              export OPTIMIZATION=--fast
              export CHPL_CFLAGS='-I${pkgs.lattice-symmetries.kernels}/include'
              export CHPL_LDFLAGS='-L${pkgs.lattice-symmetries.haskell.lib}/lib'
              export HDF5_CFLAGS='-I${pkgs.hdf5.dev}/include'
              export HDF5_LDFLAGS='-L${pkgs.hdf5}/lib -lhdf5_hl -lhdf5 -lrt'
              export TEST_DATA='${pkgs.lattice-symmetries.test-data}/share'

              rm -f src/FFI.chpl
              ln --symbolic ${pkgs.lattice-symmetries.ffi} src/FFI.chpl
            '';
          };
          python = with pkgs; lattice-symmetries.python.overrideAttrs (attrs: {
            nativeBuildInputs = (attrs.nativeBuildInputs or [ ]) ++ [
              python3Packages.black
              nodePackages.pyright
              gdb
            ];
          });
        });
    };
}
