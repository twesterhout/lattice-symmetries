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
    halide-haskell = {
      url = "github:twesterhout/halide-haskell/dynamic-dimensions";
      inputs.nixpkgs.follows = "nixpkgs";
      inputs.flake-utils.follows = "flake-utils";
    };
    haskell-python-tools = {
      url = "github:twesterhout/haskell-python-tools.nix";
      inputs.nixpkgs.follows = "nixpkgs";
    };
  };

  outputs = { self, nixpkgs, flake-utils, nix-chapel, halide-haskell, haskell-python-tools }:
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
        (final: prev: {
          halide = prev.halide.overrideAttrs
            (attrs: {
              patches = (prev.patches or [ ]) ++ [
                (final.fetchpatch {
                  name = "strict-prototypes-fix.patch";
                  url = "https://github.com/twesterhout/Halide/commit/24831b77f51f8def7fe850ba4a921e746b7a3725.patch";
                  hash = "sha256-uR3jn88UzfqpMrLFVZBrVhw4orTAXXqwNiV6qlXdjdA=";
                })
              ];
            });

          petsc = (prev.petsc.override {
            mpiSupport = true;
            petsc-optimized = true;
            petsc-scalar-type = "complex";
          }).overrideAttrs
            (attrs: rec {
              version = "3.20.5";
              name = "${attrs.pname}-${version}";
              # NOTE: We cannot use fetchFromGitLab because then PETSc tries to download SOWING...
              src = final.fetchurl {
                url = "https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-${version}.tar.gz";
                hash = "sha256-+05jd1hzevkQsF8wp4UkVjORbNCpKbe2RHrRAo2k6lo=";
              };
            });

          slepc = final.callPackage ./nix/slepc.nix { };
        })
        nix-chapel.overlays.default
        halide-haskell.overlays.default
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
          inherit (lattice-symmetries)
            kernels_v2
            haskell
            chapel
            python
            apptainer-python-minimal;
          inherit haskellPackages;
          inherit python3Packages;
          inherit petsc slepc;
        });

      devShells = flake-utils.lib.eachDefaultSystemMap (system:
        let
          pkgs = pkgs-for system;
        in
        {
          default = self.outputs.devShells.${system}.haskell;
          kernels = with pkgs; mkShell {
            buildInputs = [ halide libffcall.dev libffcall.out ];
            nativeBuildInputs = [ cmake gnumake ninja clang clang-tools ];
            shellHook = ''
              export HALIDE_PATH=${halide}
              # export CMAKE_LIBRARY_PATH=${libffcall.out}/lib:$CMAKE_LIBRARY_PATH
              # export CMAKE_INCLUDE_PATH=${libffcall.dev}/include:$CMAKE_INCLUDE_PATH
            '';
          };
          haskell = with pkgs; haskellPackages.shellFor {
            packages = ps: [ ps.lattice-symmetries-haskell ];
            withHoogle = true;
            nativeBuildInputs = with haskellPackages; [
              cabal-install
              cabal-fmt
              fourmolu
              haskell-language-server
              python3Packages.grip
            ];
            shellHook = ''
              if [ ! -f libkernels_v2.so ]; then
                gcc -shared -o libkernels_v2.so ${lattice-symmetries.kernels_v2}/lib/libkernels_v2.a
              fi
              export LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH;
            '';
          };
          chapel = with pkgs; mkShell {
            buildInputs = [
              lattice-symmetries.kernels_v2
              lattice-symmetries.haskell
              hdf5
              hdf5.dev
              halide
            ];
            nativeBuildInputs = [
              # (chapel.override {
              #   compiler = "gnu";
              #   settings = { CHPL_TARGET_MEM = "cstdlib"; CHPL_HOST_MEM = "cstdlib"; CHPL_UNWIND = "none"; CHPL_TASKS = "fifo"; CHPL_SANITIZE_EXE = "address"; CHPL_LIB_PIC = "none"; };
              # })
              # (chapel.override {
              #   llvmPackages = llvmPackages_16;
              #   compiler = "llvm";
              #   # settings = { CHPL_LIB_PIC = "pic"; CHPL_UNWIND = "system"; };
              # })
              # chapelFixupBinary
              # (chapel.override {
              #   compiler = "gnu";
              #   settings = { CHPL_COMM = "gasnet"; CHPL_COMM_SUBSTRATE = "smp"; CHPL_UNWIND = "none"; };
              # })
              (chapel.override {
                settings = { "CHPL_LIB_PIC" = "pic"; "CHPL_TARGET_CPU" = "haswell"; };
              })
              # (pr_XXX.override { compiler = "gnu"; })
              gcc
              gdb
              valgrind
              pkg-config
              prettierd
              parallel
            ];
            shellHook = ''
              # export CHPL_COMM=gasnet
              # export CHPL_COMM_SUBSTRATE=smp
              # export CHPL_RT_OVERSUBSCRIBED=yes
              # export CHPL_HOST_MEM=jemalloc
              # export CHPL_TARGET_MEM=jemalloc
              # export CHPL_LAUNCHER=none
              export CHPL_CFLAGS='-I${pkgs.lattice-symmetries.kernels_v2}/include --no-ieee-float'
              export CHPL_LDFLAGS='-L${pkgs.lattice-symmetries.haskell.lib}/lib'
              export HDF5_CFLAGS='-I${pkgs.hdf5.dev}/include'
              export HDF5_LDFLAGS='-L${pkgs.hdf5}/lib -lhdf5_hl -lhdf5 -lrt'
              export HALIDE_PATH='${pkgs.halide}'

              rm -f src/FFI.chpl
              ln --symbolic ${pkgs.lattice-symmetries.ffi} src/FFI.chpl
            '';
          };
          python = with pkgs; lattice-symmetries.python.overrideAttrs (attrs: {
            nativeBuildInputs = (attrs.nativeBuildInputs or [ ]) ++ [
              python3Packages.black
              nodePackages.pyright
              gdb
              valgrind
            ];
          });
        });
    };
}
