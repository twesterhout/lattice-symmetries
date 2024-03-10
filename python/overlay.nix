{ version
}:

final: prev: {
  pythonPackagesExtensions = prev.pythonPackagesExtensions ++ [
    (python-final: python-prev: {

      grip = python-prev.grip.overrideAttrs (attrs: {
        src = final.fetchFromGitHub {
          owner = "Antonio-R1";
          repo = "grip";
          rev = "d2efd3c6a896c01cfd7624b6504107e7b3b4b20f";
          hash = "sha256-0wgIM7Ll5WELvAOiu1TLyoNSrhJ22Y1SRbWqa3BDF3k=";
        };
        checkPhase = "true";
        installCheckPhase = "true";
      });

      lattice-symmetries = python-final.buildPythonPackage rec {
        pname = "lattice-symmetries";
        inherit version;
        src = ./.;

        buildInputs = with final; [
          lattice-symmetries.kernels_v2
          lattice-symmetries.haskell
          lattice-symmetries.chapel
        ];
        propagatedBuildInputs = with python-final; [
          cffi
          loguru
          numpy
          scipy
          sympy
        ];

        nativeCheckInputs = with python-final; [
          pip
          pytestCheckHook
          igraph
        ];

        postPatch = ''
          # - Let CFFI see chunks between `python-cffi: START` and `python-cffi: STOP`
          # - Hide `LS_HS_ATOMIC` from CFFI
          awk '/python-cffi: START/{flag=1;next}/python-cffi: STOP/{flag=0}flag' \
            ${final.lattice-symmetries.kernels_v2}/include/lattice_symmetries_types.h \
            | sed -E 's/LS_HS_ATOMIC\(([^)]+)\)/\1/' \
            >lattice_symmetries/extracted_declarations.h
          awk '/python-cffi: START/{flag=1;next}/python-cffi: STOP/{flag=0}flag' \
            ${final.lattice-symmetries.haskell}/include/lattice_symmetries_functions.h \
            >>lattice_symmetries/extracted_declarations.h
          awk '/python-cffi: START/{flag=1;next}/python-cffi: STOP/{flag=0}flag' \
            ${final.lattice-symmetries.chapel}/include/lattice_symmetries_chapel.h \
            >>lattice_symmetries/extracted_declarations.h
        '';

        preCheck = "rm -rf lattice_symmetries";

        checkPhase = ''
          runHook preCheck
          python3 -m pytest --color=yes --capture=no test/test_api.py | tee output.txt
          grep -q -E '(FAILURES|failed)' output.txt && exit 1
          runHook postCheck
        '';

        preShellHook = ''
          if test -e setup.py; then
            rm -rf build/ lattice_symmetries/*.so
            ${postPatch}
          fi
        '';
      };

    })
  ];
  lattice-symmetries = (prev.lattice-symmetries or { }) // {
    python = final.python3Packages.lattice-symmetries;
    apptainer-python-minimal = final.singularity-tools.buildImage {
      name = "my-project";
      contents = with final; [
        (python3.withPackages (ps: with ps; [
          lattice-symmetries
          loguru
          scipy
          numpy
        ]))
        coreutils
        less # for displaying docs in the Python interpreter
      ];
      diskSize = 10240;
      memSize = 5120;
    };
  };
}
