# `lattice_symmetries` tutorial

## Contents

- [Installing](#Installing)
- [Basic concepts and functions](#Basic-concepts-and-functions)
    - [Basis](#Basis)
        - [Spin basis](#Spin-basis)
        - [Spinless Fermionic basis](#Spinless-fermionic-basis)
        - [Spinful Fermionic basis](#Spinful-fermionic-basis)
    - [Expressions](#Expressions)
        - [Primitive operators](#Primitive-operators)
        - [Operator algebra](#Operator-algebra)
        - [Complex expressions](#Complex-expressions)
        - [Properties](#Properties)
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
3 |11⟩
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
basis = ls.SpinfulFermionBasis(2, number_particles=(2,1)) #We specify the numbers of fermions with spins up and down (N_up, N_down)=(2,1)
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

Expressions are an easy way to work with operators (for example, Hamiltonians) on a symbolic level using second quantization formalism.
This means that you can use primitive operators such as $\sigma^x$, $\sigma^y$, and $\sigma^z$ to build expressions for your Hamiltonian and observables.
It is also possible to sum different expressions and multiply them to each other to compose more complicated expressions.
Let's consider at first the simplest examples. 

#### Primitive operators

At first we need to import `Expr` from lattice-symmetries:
```sh
from lattice_symmetries import Expr
```
Now we will consider primitives operators defined on site with number 0 as an example, but you can also construct primitive operators residing on other lattice sites:

 - $\sigma^x$ and $S^x$:
   ```pycon
   >>> Expr("σˣ₀") == Expr("\\sigma^x_0") #It is possible to use different notations for Expr for primitive operators.
   #Here we check that they agree. Index 0 means the index of the corresponding site.
   True
   >>> Expr("Sˣ₀") == 0.5 * Expr("σˣ₀") #We check that Sˣ₀ is the same as 0.5*σˣ₀.
   True
   >>> Expr("σˣ₀").to_dense() #It is also possible to visualize the expressions as matrices.
   [[0, 1],
    [1, 0]]

   ```
   Now, we will take a look at other Pauli and spin matrices.

 - $\sigma^y$ and $S^y$:
   ```pycon
   >>> Expr("σʸ₀") == Expr("\\sigma^y_0")
   True
   >>> Expr("Sʸ₀") == 0.5 * ls.Expr("σʸ₀")
   True
   >>> Expr("σʸ₀").to_dense()
   [[0, -1j],
    [1j, 0]]
   ```

 - $\sigma^z$ and $S^z$:
   ```pycon
   >>> Expr("σᶻ₀") == Expr("\\sigma^z_0")
   True
   >>> Expr("Sᶻ₀") == 0.5 * Expr("σᶻ₀")
   True
   >>> Expr("σᶻ₀").to_dense()
   [[1, 0],
    [0, -1]]
   ```
    We see that everything works as one would expect.

 - Identity operator $\mathbb{1}$:
   ```pycon
   >>> Expr("I", particle="spin-1/2").to_dense()
   [[1, 0],
    [0, 1]]
   ```
   (*Note:* in this case, we need to explicitly specify the particle type because it cannot be deduced from the expression)



#### Operator algebra

Primitives can be combined using the `+`, `-`, and `*` operations to build more complex operators.
Furthermore, expressions can also be multiplied by scalars from the left using the `*` operator.

For example, here are a few ways to write down the Heisenberg interaction between sites 0 and 1

$$
\mathbf{S}_0 \cdot \mathbf{S}_1 = S^x_0 S^x_1 + S^y_0 S^y_1 + S^z_0 S^z_1
$$

```pycon
>>> Expr("Sˣ₀ Sˣ₁ + Sʸ₀ Sʸ₁ + Sᶻ₀ Sᶻ₁") == Expr("Sx0 Sx1 + Sy0 Sy1 + Sz0 Sz1")
True
>>> Expr("Sˣ₀ Sˣ₁ + Sʸ₀ Sʸ₁ + Sᶻ₀ Sᶻ₁") == \
             Expr("Sˣ₀") * Expr("Sˣ₁") + Expr("Sʸ₀") * Expr("Sʸ₁") + Expr("Sᶻ₀") * Expr("Sᶻ₁")
True
>>> Expr("Sˣ₀ Sˣ₁ + Sʸ₀ Sʸ₁ + Sᶻ₀ Sᶻ₁") == Expr("0.25 (σˣ₀ σˣ₁ + σʸ₀ σʸ₁ + σᶻ₀ σᶻ₁)")
True
>>> Expr("Sˣ₀ Sˣ₁ + Sʸ₀ Sʸ₁ + Sᶻ₀ Sᶻ₁") == Expr("0.5 (σ⁺₀ σ⁻₁ + σ⁺₁ σ⁻₀) + 0.25 σᶻ₀ σᶻ₁")
True

#Above we checked that different definitions of this interaction agree. Let's take a look at the corresponding matrix:
>>> Expr("Sx0 Sx1 + Sy0 Sy1 + Sz0 Sz1").to_dense()

[[ 0.25  0.    0.    0.  ]
 [ 0.   -0.25  0.5   0.  ]
 [ 0.    0.5  -0.25  0.  ]
 [ 0.    0.    0.    0.25]]
```

Under the hood, lattice-symmetries rewrites all the expressions into the canonical form, simplifying stuff along the way:

```pycon
>>> str(Expr("Sˣ₀ Sˣ₁ + Sʸ₀ Sʸ₁ + Sᶻ₀ Sᶻ₁"))
"0.25 σᶻ₀ σᶻ₁ + 0.5 σ⁺₀ σ⁻₁ + 0.5 σ⁻₀ σ⁺₁"
>>> str(Expr("0.5 (σˣ₁ + 1im σʸ₁) - σ⁺₁"))
"0.0 I"
>>> str(Expr("σ⁺₁ σ⁺₁"))
"0.0 I"
```
#### Complex expressions

So far, we defined expressions only on a few number of sites, which can be indeed easily written explicitly. Here we will consider more complicated expressions, which require other techniques.
One of the ways is to use the function `replace_indices`. For example, we can construct the sum of $\sigma^z$ operators.

```pycon
import operator #import relevant methods
from functools import reduce

expr1 = ls.Expr("σᶻ₁") #define an elementary expression
many_exprs = [expr1.replace_indices({1: i}) for i in range(4)] #we apply the permutation of vertices and make the corresponding array
expr2 = reduce(operator.add, many_exprs) #sum all elementary operations
print(expr2)
>>>
σᶻ₀ + σᶻ₁ + σᶻ₂ + σᶻ₃
```

Another way is to define an expression on a graph defined by edges:

```pycon
edges = [(i,) for i in range(4)] #define graph of 4 vertices 
expr = ls.Expr("σᶻ₁", sites=edges) #the elementary expression is applied for all vertices
>>>
σᶻ₀ + σᶻ₁ + σᶻ₂ + σᶻ₃
```
One can see that the results are the same.

#### Properties

One can make standart operations on expressions, such make an adjoint, as well

```pycon
a = ls.Expr("c†₀ c₁")
b=a.adjoint() #Makes an adjoint expression
```

It is also possible to check various properties of expressions:
```pycon
>>> a=Expr("\\sigma^z_0")
>>> a.is_real #If the corresponding matrix real
True
>>> a.is_hermitian #If the corresponding matrix hermitian
True
>>> a.is_identity #If the corresponding matrix identity
False
>>> a.number_sites #Number of sites in the expression
1
```

### Operators

Based on an expression and a basis, we can build the corresponding operator acting on the chosen basis. Let's at first consider the Hubbard model on two sites without spins.
One of the ways to contruct an operator is to create an expression and then create the operator:

```pycon
op_expr = ls.Expr("-c†₀ c₁-c†₁ c₀+2.0 n₀ n₀+2.0 n₁ n₁")
basis=ls.SpinlessFermionBasis(2)
basis.build() #It is nessecary to build the basis before creating an operator
opr=ls.Operator(basis, op_expr) #Here we create the operator
```
Another way is to create the expression for each small expression and then add them to each other:
```pycon
op_expr = ls.Expr("-c†₀ c₁")
op_expr -= ls.Expr("c†₁ c₀")
op_expr += 2*ls.Expr("n₀ n₀")
op_expr += 2*ls.Expr("n₁ n₁")
opr=ls.Operator(basis, op_expr)
opr.shape() #We can check the shape of the operator
>>>
(4,4)
```

We can make a multiplication of the operator on a vector using standart python notation.

```pycon
vec1=np.array([1.0,0,0,0]) #One should work with floats
opr@vec1
>>>
[0. 0. 0. 0.]

vec2=np.array([0,1.0,0,0])
opr@vec2
>>>
[ 0.  2. -1.  0.]

vec3=np.array([1.0,2.0,0,0])
opr@vec3
>>>
[ 0.  4. -2.  0.]
```
From this operations one can build the matrix form of an operator.

### Symmetries

The full power of `lattice_symmetries` manifests if one use symmetries when constructing 
symmetry adapted basis and linear operators acting on teh corresponding Hilbert space. 
The symmetries can be constructed with the help of expressions as well. 

## Examples

Here we will take a look at different examples of `lattice_symmetries` applications.

### Simple ED

We will strat with the simplest example of exact diagonalization. We will consider Heisenberg chain on 10 sites and diagonalization with the help of symmetries.
For that, we combine all the considerd methods.
At first we create the expression for Heisenberg chain on 10 sites using the methods from the section [Complex expressions](#Complex-expressions).

```pycon
# Constructing the basis
    basis = ls.SpinBasis(
        number_spins=number_spins,
        # NOTE: we don't actually need to specify hamming_weight when spin_inversion
        # is set. The library will guess that hamming_weight = number_spins / 2.
        spin_inversion=-1,
        symmetries=symmetries,
    )
    basis.build()  # Build the list of representatives, we need it since we're doing ED

# Constructinf the expression
edges = [(i, (i + 1) % number_spins) for i in range(number_spins)]
expr = ls.Expr("2 (σ⁺₀ σ⁻₁ + σ⁺₁ σ⁻₀) + σᶻ₀ σᶻ₁", sites=edges)
print("Expression:", expr)

# Alternatively, we can create expr in an algebraic way:
two_site_expr = ls.Expr("2 (σ⁺₀ σ⁻₁ + σ⁺₁ σ⁻₀) + σᶻ₀ σᶻ₁")
many_exprs = [two_site_expr.replace_indices({0: i, 1: j}) for (i, j) in edges]
expr2 = reduce(operator.add, many_exprs)
assert expr == expr2 #we check that the expressions are equal

# Construct the Hamiltonian
hamiltonian = ls.Operator(basis, expr)


# Diagonalize the Hamiltonian using ARPACK
eigenvalues, eigenstates = scipy.sparse.linalg.eigsh(hamiltonian, k=1, which="SA")
print("Ground state energy is {}".format(eigenvalues[0]))
assert np.isclose(eigenvalues[0], -18.06178542)
```

### More complicated ED

Now let's consider a more complicated example of ED.

### Time evolution

Another example of capabilities of `lattice_symmetries` is time evolution.
