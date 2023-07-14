{ version
}:

final: prev: {
  pythonPackagesExtensions = prev.pythonPackagesExtensions ++ [
    (python-final: python-prev: {
      lattice-symmetries = python-final.buildPythonPackage {
        pname = "lattice-symmetries";
        inherit version;
        src = ./.;

        buildInputs = with prev; [
          lattice-symmetries.kernels
          lattice-symmetries.haskell
          lattice-symmetries.chapel
        ];
        propagatedBuildInputs = with python-final; [
          cffi
          loguru
          numpy
          scipy
        ];

        postPatch = ''
          awk '/python-cffi: START/{flag=1;next}/python-cffi: STOP/{flag=0}flag' \
            ${prev.lattice-symmetries.kernels}/include/lattice_symmetries_types.h \
            >lattice_symmetries/extracted_declarations.h
          awk '/python-cffi: START/{flag=1;next}/python-cffi: STOP/{flag=0}flag' \
            ${prev.lattice-symmetries.haskell}/include/lattice_symmetries_functions.h \
            >>lattice_symmetries/extracted_declarations.h
        '';

        preCheck = "rm -rf lattice_symmetries";

        checkPhase = ''
          runHook preCheck
          python3 -m pytest --color=yes --capture=no test/test_api.py | tee output.txt
          grep -q -E '(FAILURES|failed)' output.txt && exit 1
          runHook postCheck
        '';

        nativeCheckInputs = with python-final; [ pytestCheckHook ];
      };
    })
  ];
  lattice-symmetries = (prev.lattice-symmetries or { }) // {
    python = prev.python3Packages.lattice-symmetries;
  };
}
