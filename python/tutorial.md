# `lattice_symmetries` tutorial

## Contents

- [Installing](#Installing)
- [Introduction, basic concepts and functions](#Intro)
    - [Basis](#Basis)
        - [Spin basis](#Spin_basis)
        - [Fermionic basis](#Fermionic basis)
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

The basic objects upon which the machinery is built are `Basis`, `Expression`, `Operator`, and `Symmetry`.
- `Basis` objects allow to construct basis of many-body Hilbert space consisting of spins or fermions with appropriate symmetries.
Each basis state is represented as a sequence of $0$s and $1$s (i.e. bits), which can be interpreted as a binary number.
- `Expression` objects is a nice way to work with symbolic representations of operators. It is possible to sum expressions, multiply them to numbers and each other.
`Expressions` allow not to think about an explicit matrix representation of an operator, and user can work directly with analytical formulae as if they are written in paper. 
- `Operator` object is an actual operator that can act on individual basis states and their linear combinations. 
- `Symmetry` is a method to work with symmetry adapted basises.
If an operator has symmmetries, it is useful to work in symmetry-adapted basis, since it can drastically decrease the dimension of the Hilbert space.

Now we will take a look at basic functions and methods for these objects. 

### Basis
#### Spin basis
Let's take a look at simple examples, and at first we will not consider additional symmetries.
The simplest example would be a spin basis:

```sh
import lattice_symmetries as ls
 
basis = ls.SpinBasis(3) #We create basis consisting of three $\frac{1}{2}$-spins (each spin can be |0⟩ or |1⟩) 
basis.build() #build basis
print(basis.index(basis.states)) #Print array of indices of basis states
print(basis.number_states) #Print total number of states in the basis
for i in range(basis.number_states): #Here we print all basis states as they stored in memory. Without symmetries, basis states equal their index
    print(basis.states[i], basis.state_to_string(i)) #The last function represents basis states as a sequence of $0$s and $1$s  
```

The result looks like:
```sh
[0 1 2 3 4 5 6 7]
8
0 |000⟩
1 |001⟩
2 |010⟩
3 |011⟩
4 |100⟩
5 |101⟩
6 |110⟩
7 |111⟩
```
The basis states are equal to their indices and binary represenations as they should.

#### Fermionic basis
We can also consider fermionic basis:


### Expressions

Let's take a look at a few examples.


### Operators

Based on an expression and a basis, we can build the corresponding operator acting on the chosen basis. Let's consider an example.



### Symmetries

The full power of `lattice_symmetries` manifests if one use symmetries when constructing 
symmetry adapted basis and linear operators acting on teh corresponding Hilbert space.

## Examples

Here we will take a look at different examples of `lattice_symmetries` applications.

### Simple ED

Now we will consider the simplest example of exact diagonalization.

### More complicated ED

Now let's consider a more complicated example of ED.

### Time evolution

Another example of capabilities of `lattice_symmetries` is time evolution.
