# `lattice_symmetries` tutorial

## Contents

- [Installing](#Installing)
- [Basic concepts and functions](#Basic-concepts-and-functions)
    - [Basis](#Basis)
        - [Spin basis](#Spin-basis)
        - [Spinless Fermionic basis](#Spinless-fermionic-basis)
        - [Spinful Fermionic basis](#Spinful-fermionic-basis)
    - [Expressions](#Expressions)
    - [Operators](#Operators)
    - [Symmetries](#Symmetries)
- [Examples](#Examples)
    - [Simple ED](#Simple-ED)
    - [More complicated ED](#More-complicated-ED)
    - [Time evolution](#Time-evolution)

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

## Basic concepts and functions

The `lattice_symmetries` is a powerful package that nicely works with many-body quantum spin systems
and takes into account system symmetries.
The `lattice_symmetries` implements fast matrix-vector multiplication that can be applied to various problems, 
such as time evolution or exact diagonalization of many-body Hamiltonians.

The basic objects upon which the machinery is built are `Basis`, `Expression`, `Operator`, and `Symmetry`.
- `Basis` object stores a basis of a many-body Hilbert space consisting of spins or fermions with appropriate symmetries.
Each basis state is represented as a sequence of 0s and 1s (i.e. a sequence of bits), which can be interpreted as a binary number.
- `Expression` object is a nice way to work with symbolic representations of operators. It is possible to sum expressions, and multiply them by numbers and each other.
`Expression` allows not to think about an explicit matrix representation of an operator, and the user can work directly with analytical formulae. 
- `Operator` object is an actual operator that can act on individual basis states and their linear combinations. 
- `Symmetry` is a method to work with symmetry adapted basises.
If an operator has symmmetries, it is useful to work in symmetry-adapted basis, since it can drastically decrease the dimension of the Hilbert space.

Now we will take a look at basic functions and methods for these objects. 

### Basis
#### Spin basis
Let's look at simple examples; at first we will not consider additional symmetries.

The simplest example would be a spin basis:

```sh
import lattice_symmetries as ls
 
basis = ls.SpinBasis(3) #We create basis consisting of three $\frac{1}{2}$-spins (each spin can be |0⟩ or |1⟩) 
basis.build() #build basis
print(basis.index(basis.states)) #Print array of indices of basis states
print(basis.number_states) #Print total number of states in the basis
for i in range(basis.number_states): #Here we print all basis states as they stored in memory. Without symmetries, basis states equal their index
    print(basis.states[i], basis.state_to_string(basis.states[i]))  
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
The basis states are equal to their indices and binary representations, as they should be.

We can also consider only a part of the whole Hilbert space with a given number of spin ups (i.e. hamming weight of binary representations).
We can specify it as follows:
```sh
basis = ls.SpinBasis(3, hamming_weight=2) #We want the subspace with only 2 spins up 
basis.build()
for i in range(basis.number_states): #Print the states in the basis
    print(basis.states[i], basis.state_to_string(basis.states[i]))  
```
The output is:
```sh
3 |011⟩
5 |101⟩
6 |110⟩
```
We see that basis states include only spins with 2 spins up. It is also interesting to note that in this case a basis state is not equal to its index.

Sometimes, our system has spin inversion symmetry, and we can additionally shorten our basis. In this case, we can specify it as follows:
```sh
basis = ls.SpinBasis(4, spin_inversion=-1) #Spin inversion is present. It is not nessecary to specify hamming weight here, since it is fixed by spin inversion symmetry. 
basis.build()
for i in range(basis.number_states):
    print(basis.states[i], basis.state_to_string(basis.states[i]))  
```
We have the basis states:
```sh
3 |011⟩
5 |101⟩
6 |110⟩
```

#### Spinless Fermionic basis
We can also consider the basis of fermions without spins. The basis states are stored as integers as for spin basis. However the binary representation has a second quantization interpretation.
Each basis state is given by the sequence of 0s and 1s, where 1 means a fermion on the corresponding site, and 0 means that the site is vacant.
Let's consider the simplest example of fermions on two sites:
```sh
basis = ls.SpinlessFermionBasis(2) #We create fermionic basis on 2 sites
basis.build()
for i in range(basis.number_states):
    print(basis.states[i], basis.state_to_string(basis.states[i])) 
```
which gives:
```sh
0 |00⟩
1 |01⟩
2 |10⟩
3 |11⟩111
```
as one would expect.

We can specify the number of particles as well:
```sh
basis = ls.SpinlessFermionBasis(4, number_particles=3) #We create fermionic basis on 4 sites with only 3 fermions
basis.build()
for i in range(basis.number_states):
    print(basis.states[i], basis.state_to_string(basis.states[i])) 
```
which gives:
```sh
7 |0111⟩
11 |1011⟩
13 |1101⟩
14 |1110⟩
```
We can see that the basis consists of states with three fermions.

#### Spinful Fermionic basis

The last case includes fermions with spin. The binary representations of basis states can be read as a pair of (fermions with spin up on a lattice, fermions with spin down on a lattice).
We can create a basis of spinful fermions as follows:
```sh
basis = ls.SpinfulFermionBasis(2) #We create basis of spinful fermions on 2 sites
basis.build()
for i in range(basis.number_states):
    print(basis.states[i], basis.state_to_string(basis.states[i])) 
```
which gives:
```sh
0 |00⟩|00⟩
1 |00⟩|01⟩
2 |00⟩|10⟩
3 |00⟩|11⟩
4 |01⟩|00⟩
5 |01⟩|01⟩
6 |01⟩|10⟩
7 |01⟩|11⟩
8 |10⟩|00⟩
9 |10⟩|01⟩
10 |10⟩|10⟩
11 |10⟩|11⟩
12 |11⟩|00⟩
13 |11⟩|01⟩
14 |11⟩|10⟩
15 |11⟩|11⟩
```
We see that binary representation now means the second quantization of fermions with two spins.

As before, we can specify a sector of the total Hilbert space with a given number of fermions with spin down and spin up:

```sh
basis = ls.SpinfulFermionBasis(2, number_particles=(2,1)) #We specify the numbers of fermions with spins down and up (N_down, N_up)=(2,1)
basis.build()
for i in range(basis.number_states):
    print(basis.states[i], basis.state_to_string(basis.states[i])) 
```
which gives:
```sh
7 |01⟩|11⟩
11 |10⟩|11⟩
```
as expected.

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
