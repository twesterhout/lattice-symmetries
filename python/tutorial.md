# `lattice_symmetries` tutorial

## Contents

- [Installing](#Installing)
- [Introduction, basic concepts and functions](#Intro)
    - [Basis](#Basis)
    - [Expressions](#Expressions)
    - [Operators](#Operators)
    - [Symmetries](#Symmetries)
- [Examples](#Examples)
    - [Simple ED](#Simple_ED)
    - [More complicated ED](#Complex_ED)
    - [Time evolution](#Time_evolution)

## Installing

The first step before we can apply `lattice_symmetries` is to install Nix. The installation process slightly depends on your system and can be found in the [official documentaion](https://nix.dev/install-nix#install-nix).
In the case of WSL (Windows Subsystem for Linux), it is preferred to install a single-user version:

```sh
curl -L https://nixos.org/nix/install | sh -s -- --no-daemon
```

You can check that Nix is installed by opening a new terminal and typing:

```sh
nix --version
```

After that, it is necessary to add the support of flakes and experimental features into the configuration Nix file `nix.conf`. The file can be located in one of the two paths (starting from the root directory): `~/.config/nix/nix.conf` or `\etc\nix\nix.conf`. If the file is absent, you should create it in one of the mentioned directories. The file should contain the following line:

```sh
experimental-features = nix-command flakes
```
Now, the Nix is ready.

The next step is to actually install `lattice_symmetries`. For that, you need to clone the `lattice_symmetries` from github:

```sh
git clone https://github.com/twesterhout/lattice-symmetries.git
cd lattice-symmetries
git submodule update --init --recursive
```

Then you should move to the `lattice-symmetries/python` folder and prepare Nix files:
```sh
cd python
patchPhase || nix
```

Now you can build everything using Nix:

```sh
nix develop .#python
```

Voila, and you have the working package.
If you open a new terminal, the last step should be repeated ,i.e., moving to  `lattice-symmetries/python` and typing `nix develop .#python`.

## Introduction, basic concepts and functions

The `lattice_symmetries` is a powerful package that nicely works with many-body quantum spin systems
and takes into account system symmetries.
The `lattice_symmetries` implements fast matrix-vector multiplication that can be applied to various problems, 
such as time evolution or exact diagonalization of many-body Hamiltonians.

The basic objects upon which the machinery is built are Basis, Expression, Operator, and Symmetry.
- Basis objects allow to construct basis of many-body hamiltonians (spin and fermion) with appropriate symmetries.
The Basis objects als
- Expression objects is a nice way to work with symbolic representations of operators. It is possible to add

### Basis

### Expressions

### Operators

### Symmetries

The full power of `lattice_symmetries` manifests if one use symmetries when constructing 
symmetry adapted basis and linear operators acting on teh corresponding Hilbert space.

## Examples

Here we will take a look at different examples of lattice_symmetries applications.

### Simple ED

Now we will consider the simplest example of exact diagonalization.

### More complicated ED

Now let's consider a more complicated example of ED.

### Time evolution

Another example of capabilities of `lattice_symmetries` is time evolution.
