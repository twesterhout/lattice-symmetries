{
  description = "A project using lattice-symmetries";

  # We specify additional cache locations to avoid a very lengthly compilation of
  # dependencies. Since Nix ensures reproducibility, it is safe to fetch pre-built
  # packages from cache.
  nixConfig = {
    extra-substituters = "https://twesterhout-chapel.cachix.org";
    extra-trusted-public-keys = "twesterhout-chapel.cachix.org-1:bs5PQPqy21+rP2KJl+O40/eFVzdsTe6m7ZTiOEE7PaI=";
  };

  inputs = {
    # We follow nixpkgs and flake-utils from lattice-symmetries to have a higher
    # probability of dependencies being cached.
    nixpkgs.follows = "lattice-symmetries/nixpkgs";
    flake-utils.follows = "lattice-symmetries/flake-utils";
    lattice-symmetries.url = "github:twesterhout/lattice-symmetries";
    flake-compat = {
      url = "github:edolstra/flake-compat";
      flake = false;
    };
  };

  outputs = inputs:
    let
      pkgs-for = system: import inputs.nixpkgs {
        inherit system;
        overlays = [ inputs.lattice-symmetries.overlays.default ];
      };

      # Our Python dependencies
      my-python-packages = ps: with ps; [
        # jupyter
        lattice-symmetries
        ipykernel
        loguru
        matplotlib
        numpy
        pandas
        scipy
        seaborn
      ];
    in
    {
      packages = inputs.flake-utils.lib.eachDefaultSystemMap (system:
        with (pkgs-for system); {
          default = singularity-tools.buildImage {
            name = "my-project";
            contents = [
              (python3.withPackages my-python-packages)
              coreutils
              less # for displaying docs in the Python interpreter
            ];
            diskSize = 10240;
            memSize = 5120;
          };
        });

      devShells = inputs.flake-utils.lib.eachDefaultSystemMap (system:
        with (pkgs-for system); {
          default = mkShell {
            nativeBuildInputs = [
              (python3.withPackages my-python-packages)
              # LSP support for Python
              python3Packages.black
              nodePackages.pyright
            ];
            shellHook = ''
              export PROMPT_COMMAND=""
              export PS1='🐍 Python ${python3.version} \w $ '
              export LS_PATH=${lattice-symmetries.python}
            '';
          };
        });
    };
}
