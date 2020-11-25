from setuptools import setup

setup(
    name="lattice-symmetries",
    version="0.2.0",
    description="Easily work with quantum many-body bases taking all lattice symmetries into account",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3 :: Only",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    url="http://github.com/twesterhout/lattice-symmetries",
    author="Tom Westerhout",
    author_email="14264576+twesterhout@users.noreply.github.com",
    license="BSD3",
    packages=["lattice_symmetries"],
    install_requires=["numpy"],
    zip_safe=False,
)
