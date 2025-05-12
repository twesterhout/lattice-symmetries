{
  description = "twesterhout/lattice-symmetries";
  inputs = {
    # nixpkgs.url = "github:nixos/nixpkgs/nixos-unstable";
    nixpkgs.follows = "quantum-nix/nixpkgs";
    nix-gl-host.url = "github:numtide/nix-gl-host";
    nix-gl-host.inputs.nixpkgs.follows = "nixpkgs";
    quantum-nix.url = "github:twesterhout/quantum-nix";
    # quantum-nix.inputs.nixpkgs.follows = "nixpkgs";
    quantum-nix.inputs.nix-gl-host.follows = "nix-gl-host";
  };
  nixConfig = {
    extra-substituters = [
      "https://twesterhout.cachix.org"
      "https://nix-community.cachix.org"
      "https://cuda-maintainers.cachix.org"
    ];
    extra-trusted-public-keys = [
      "twesterhout.cachix.org-1:AtBrVtHiRtg7piQOwT9IWx3N/+q+lM6RrpxzpeT3zAE="
      "nix-community.cachix.org-1:mB9FSh9qf2dCimDSUo8Zy7bkq5CX+/rkCWyvRCYg3Fs="
      "cuda-maintainers.cachix.org-1:0dq3bujKpuEPMCX6U4WylrUDZ9JyUG0VpVZa7CNfq5E="
    ];
  };

  outputs = inputs:
    let
      version = "3.0.0";

      inherit (inputs.nixpkgs) lib;
      forEachSystem = f: lib.mapAttrs f inputs.nixpkgs.legacyPackages;

      overlay = lib.composeManyExtensions [
        inputs.nix-gl-host.overlays.default
        inputs.quantum-nix.overlays.default
        (import ./python.nix { inherit version; })

        # igraph's test suite pulls in a ton of dependencies...
        (final: prev: {
          pythonPackagesExtensions = prev.pythonPackagesExtensions ++ [
            (python-final: python-prev: lib.optionalAttrs (python-prev.python.pythonOlder "3.12") {
              igraph = python-prev.igraph.overridePythonAttrs (attrs: { doCheck = false; });
            })];
        })
      ];
      pkgs-for-cpu = system: import inputs.nixpkgs { inherit system; overlays = [ overlay ]; };
      pkgs-for-cuda = system: import inputs.nixpkgs {
        inherit system;
        config = { allowUnfree = true; cudaSupport = true; cudaCapabilities = [ "7.0" ]; cudaForwardCompat = true; };
        overlays = [ overlay ];
      };
    in
    {
      overlays.default = overlay;
      packages = forEachSystem (system: _:
        let pkgs = pkgs-for-cpu system; in {
          inherit (pkgs) python3Packages python311Packages python312Packages;
          apptainer = pkgs.singularity-tools.buildImage {
            name = "lattice-symmetries";
            contents = [
              pkgs.coreutils
              pkgs.python3.stdenv.cc
              (pkgs.python3.withPackages (ps: [ ps.lattice-symmetries ]))
            ];
            diskSize = 10240;
            memSize = 5120;
          };
          docker = pkgs.dockerTools.buildImage {
            name = "lattice-symmetries";
            tag = "latest";
            copyToRoot = pkgs.buildEnv {
              name = "image-root";
              paths = [
                (pkgs.python3.withPackages (ps: [ ps.lattice-symmetries ]))
                pkgs.zig
              ];
              pathsToLink = [ "/bin" ];
            };
          };
        });
      devShells = forEachSystem (system: _:
        let pkgs = pkgs-for-cpu system;
        in
        {
          python = pkgs.python3Packages.lattice-symmetries.overridePythonAttrs (attrs: {
            nativeBuildInputs = with pkgs; (attrs.nativeBuildInputs or []) ++ [ pkgs.zig pkgs.nix-tree pkgs.python3Packages.ipython ]; # pkgs.nix-gl-host ];
          });
          testing = with pkgs; mkShell { nativeBuildInputs = [ (python3.withPackages (ps: with ps; [ lattice-symmetries ])) ]; };
        });
      formatter = forEachSystem (system: pkgs: pkgs.nixpkgs-fmt);
    };
}
