from setuptools import setup

setup(
    packages=["lattice_symmetries"],
    cffi_modules=["lattice_symmetries/_build_extension.py:ffibuilder"],
)
