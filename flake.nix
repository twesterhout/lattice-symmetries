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
      ];
      pkgs-for = system: import inputs.nixpkgs { inherit system; overlays = [ overlay ]; };
    in
    {
      overlays.default = overlay;
      packages = forEachSystem (system: _:
        let pkgs = pkgs-for system; in {
          inherit (pkgs) python3Packages python311Packages python312Packages;
        });
      devShells = forEachSystem (system: _:
        let pkgs = pkgs-for system;
        in
        {
          python = pkgs.python3Packages.lattice-symmetries;
          testing = with pkgs; mkShell {
            nativeBuildInputs = [
              (python3.withPackages (ps: with ps; [ lattice-symmetries ]))
            ];
          };
        });
      formatter = forEachSystem (system: pkgs: pkgs.nixpkgs-fmt);
    };
}
