from setuptools import setup, Extension
import os
import re


def get_version(package):
    pwd = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(pwd, package, "__init__.py"), "r") as input:
        result = re.search(r'__version__\s*=\s*[\'"]([^\'"]*)[\'"]', input.read())
    if not result:
        raise ValueError("failed to determine {} version".format(package))
    return result.group(1)


setup(
    name="lattice-symmetries",
    version=get_version("lattice_symmetries"),
    description="See README.md",
    url="http://github.com/twesterhout/lattice-symmetries-haskell",
    author="Tom Westerhout",
    author_email="14264576+twesterhout@users.noreply.github.com",
    license="BSD3",
    packages=["lattice_symmetries"],
    setup_requires=["cffi>=1.15.0"],
    cffi_modules=["lattice_symmetries/build_extension.py:ffibuilder"],
    install_requires=[
        "cffi>=1.15.0",
        "numpy>=1.23.0",
        "scipy>=1.8.0",
    ],
    include_package_data=True,
    # package_data={"lattice-symmetries": ["lattice-symmetries-chapel-*"]}
    zip_safe=False,
)

