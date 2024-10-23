# `lattice_symmetries` tutorial

## Contents

- [Installing](#Installing)
- [Basic concepts and functions](#Basic-concepts-and-functions)
    - [Expressions](#Expressions)
        - [Primitive expressions](#Primitive-expressions)
        - [Algebra of expressions](#Algebra-of-expressions)
        - [Complex expressions](#Complex-expressions)
        - [Properties](#Properties)
	- [Basis](#Basis)
        - [Spin basis](#Spin-basis)
        - [Spinless Fermionic basis](#Spinless-fermionic-basis)
        - [Spinful Fermionic basis](#Spinful-fermionic-basis)
		- [Basis from Expressions](#Basis-from-expressions)
    - [Operators](#Operators)
    - [Symmetry](#Symmetry)
		- [Symmetry as Permutations](#Symmetry-as-permutations)
		- [Symmetry from Expressions](#Symmetry-from-expressions)
		- [Symmetry-adapted basis](#Symmetry-adapted-basis)
			- [Manually generated representations](#Manually-generated-representations)
			- [Functions `hilbert_space_sectors` and `ground_state_sectors`](#Sectors-section)
- [Examples](#Examples)
    - [ED of 1D chain](#Ed-of-1d-chain)
    - [ED of Star graph](#Ed-of-star-graph)
    - [Time evolution](#Time-evolution)
	- [Spin structure factor](#Spin-structure-factor)

## Installing

The first step before we can apply `lattice_symmetries` is to install Nix. The installation process slightly depends on your system and can be found in the [official documentation](https://nix.dev/install-nix#install-nix).
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

Voila and you have the working package.
If you open a new terminal, the last step should be repeated i.e., moving to  `lattice-symmetries/python` and typing `nix develop .#python`.

## Basic concepts and functions

The `lattice_symmetries` is a powerful package that nicely works with many-body quantum spin systems
and takes into account system symmetries.
The `lattice_symmetries` implements fast matrix-vector multiplication that can be applied to various problems, 
such as time evolution or exact diagonalization of many-body Hamiltonians.

The basic objects upon which the machinery is built are `Basis`, `Expression`, `Operator`, and `Symmetry`.
The `Symmetries` is not a separate class, however, it lies in the core of `lattice_symmetries`, so we will consider it separately.
- `Expression` object is a nice way to work with symbolic representations of operators. It is possible to sum expressions and multiply them by numbers and each other.
`Expression` does not require the user to consider an explicit matrix representation of an operator and allows the user to work directly with analytical formulae. 
- `Basis` object stores a basis of a many-body Hilbert space consisting of spins or fermions with appropriate symmetries.
Each basis state is represented as a sequence of 0s and 1s (i.e. a sequence of bits), which can be interpreted as a binary number.
- `Operator` object is an actual operator that can act on individual basis states and their linear combinations. 
- `Symmetry` is a method to work with symmetry-adapted basis.
If an operator has symmetries, it is useful to work with a symmetry-adapted basis since it can drastically decrease the dimension of the Hilbert space.

Now, we will take a look at basic functions and methods for these entities. 

### Expressions

Expressions are an easy way to work with operators (such as Hamiltonians) on a symbolic level using the second quantization formalism.
This means that you can use primitive symbols of well-known operators such as $\sigma^x$, $\sigma^y$, and $\sigma^z$ to build expressions for your Hamiltonian and observables.
It is also possible to sum different expressions and multiply them to each other to compose more complicated expressions.
Let's consider the simplest examples first. 

#### Primitive expressions

At first, we need to import `Expr` from lattice-symmetries:
```pycon
from lattice_symmetries import Expr
```
Now we will consider primitive symbols defined on a site with index 0 as an example, but you can also construct them residing on other lattice sites:

 - $\sigma^x$ and $S^x$:
   ```pycon
   >>> Expr("σˣ₀") == Expr("\\sigma^x_0") #It is possible to use different notations for Expr for primitives.
   #Here, we check that they agree. Index 0 means the index of the corresponding site.
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



#### Algebra of expressions

Primitives can be combined using the `+`, `-`, and `*` operations to build more complex expressions.
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

#Above, we checked that different definitions of this interaction agree. Let's take a look at the corresponding matrix:
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

So far, we have defined expressions only for a system with a few sites, which can be easily written explicitly. Here, we will consider more complicated expressions, which require other techniques.
One of the ways to do this is by using the function `replace_indices`. For example, we can construct the sum of $\sigma^z$s defined on different sites: 

```pycon
import operator #Import relevant methods
from functools import reduce

expr1 = Expr("σᶻ₁") #Define an elementary expression
many_exprs = [expr1.replace_indices({1: i}) for i in range(4)] #We apply the permutation of vertices and make the corresponding array
expr2 = reduce(operator.add, many_exprs) #Sum all elementary operations
print(expr2)
>>>
σᶻ₀ + σᶻ₁ + σᶻ₂ + σᶻ₃
```

Another way is to define an expression on a (hyper)graph defined by edges:

```pycon
edges = [(i,) for i in range(4)] #Define graph of 4 vertices. Edges consist only of one element, because the elementary expression is linear
expr = Expr("σᶻ₁", sites=edges) #The elementary expression is applied for all vertices
>>>
σᶻ₀ + σᶻ₁ + σᶻ₂ + σᶻ₃
```

One can use the function `on` for the same effect:
```pycon
e=Expr("σᶻ₁")
expt=e.on(edges)
>>>
σᶻ₀ + σᶻ₁ + σᶻ₂ + σᶻ₃
```

It is also possible to use `iGraph` to define an underlying (hyper)graph:
```pycon
import igraph as ig

e=Expr("σ⁺₀ σ⁻₁ + σ⁺₁ σ⁻₀") #Here, we have two-site interaction, so we need an actual graph with edges consisting of two vertices
expr=e.on(ig.Graph.Lattice(dim=[2, 2])) # The graph is square 2x2
>>>
σ⁺₀ σ⁻₁ + σ⁺₀ σ⁻₂ + σ⁻₀ σ⁺₁ + σ⁻₀ σ⁺₂ + σ⁺₁ σ⁻₃ + σ⁻₁ σ⁺₃ + σ⁺₂ σ⁻₃ + σ⁻₂ σ⁺₃
```
However, the method of defining a Hamiltonian on a graph works correctly only if the elementary expression (defined on one site) is homogeneous.
It is also important that the edges of a (hyper)graph should have the same number of vertices as the degree of the elementary expression.

#### Properties

One can make standard operations on expressions, such as making an adjoint, as well:

```pycon
a = Expr("c†₀ c₁")
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

### Basis

In this section, we will consider the generation of bases of various types.
We will not consider symmetry-adapted bases yet; we will do so in the section [Symmetry-adapted basis](#Symmetry-adapted-basis).


#### Spin basis
Let's look at simple examples; at first we will not consider additional symmetries.

The simplest example would be a spin basis:

```pycon
import lattice_symmetries as ls
 
basis = ls.SpinBasis(3) #We create basis consisting of three $\frac{1}{2}$-spins (each spin can be |0⟩ or |1⟩) 
basis.build() #build basis
print(basis.index(basis.states)) #Print array of indices of basis states
print(basis.number_states) #Print total number of states in the basis

def show_basis(basis): #Here we print all basis states as they are stored in memory 
	for i in range(basis.number_states): 
		print(basis.states[i], basis.state_to_string(basis.states[i])) 
show_basis(basis)
>>>
[0 1 2 3 4 5 6 7]
8
0 |000⟩ #Without symmetries, basis states are equal to their index
1 |001⟩
2 |010⟩
3 |011⟩
4 |100⟩
5 |101⟩
6 |110⟩
7 |111⟩
```
The basis states are equal to their indices and binary representations, as they should be.

We can also consider only a part of the whole Hilbert space with a given number of spins up (i.e. the hamming weight of binary representations).
We can specify it as follows:
```pycon
basis = ls.SpinBasis(3, hamming_weight=2) #We want the subspace with only 2 spins up 
basis.build()
show_basis(basis)
>>>
3 |011⟩
5 |101⟩
6 |110⟩
```
We see that basis states include only spins with 2 spins up. It is also interesting to note that, in this case, a basis state is not equal to its index.

Sometimes, our system has spin inversion symmetry, and we can additionally shorten our basis. In this case, we can specify it as follows:
```pycon
basis = ls.SpinBasis(4, spin_inversion=-1) #Spin inversion is present. It is not nessecary to specify hamming weight here, since it is fixed by spin inversion symmetry. 
basis.build()
show_basis(basis) 
>>>
3 |011⟩
5 |101⟩
6 |110⟩
```

#### Spinless Fermionic basis
We can also consider the basis of fermions without spins. The basis states are stored as integers as for spin basis. However, the binary representation has a second quantization interpretation.
Each basis state is given by the sequence of 0s and 1s, where 1 means a fermion on the corresponding site, and 0 means that the site is vacant.
Let's consider the simplest example of fermions on two sites:
```pycon
basis = ls.SpinlessFermionBasis(2) #We create a fermionic basis on 2 sites
basis.build()
show_basis(basis)
>>>
0 |00⟩
1 |01⟩
2 |10⟩
3 |11⟩
```
as one would expect.

We can specify the number of particles as well:
```pycon
basis = ls.SpinlessFermionBasis(4, number_particles=3) #We create a fermionic basis on 4 sites with only 3 fermions
basis.build()
show_basis(basis)
>>>
7 |0111⟩
11 |1011⟩
13 |1101⟩
14 |1110⟩
```
We can see that the basis consists of states with three fermions.

#### Spinful Fermionic basis

The last case includes fermions with spin. The binary representations of basis states can be read as a pair of (fermions with spin up on a lattice, fermions with spin down on a lattice).
We can create a basis of spinful fermions as follows:
```pycon
basis = ls.SpinfulFermionBasis(2) #We create a basis of spinful fermions on 2 sites
basis.build()
show_basis(basis)
>>>
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

```pycon
basis = ls.SpinfulFermionBasis(2, number_particles=(2,1)) #We specify the numbers of fermions with spins up and down (N_up, N_down)=(2,1)
basis.build()
show_basis(basis)
>>>
7 |01⟩|11⟩
11 |10⟩|11⟩
```

#### Basis from Expressions

Before we created a basis explicitly, however, there was a way to construct a basis directly from expressions.
However, in this case, one can obtain only a full Fock basis without restriction on the number of particles.

```pycon
expr=ls.Expr("Sˣ₀ Sˣ₁ + Sʸ₀ Sʸ₁ + Sᶻ₀ Sᶻ₁")
basis=expr.full_basis()
basis.build()
```


### Operators

Based on an expression and a basis, we can build the corresponding operator acting on the chosen basis. Let's at first consider the Hubbard model on two sites without spins.
One of the ways to construct an operator is to create an expression and then create the operator:

```pycon
import lattice_symmetries as ls

op_expr = ls.Expr("-c†₀ c₁-c†₁ c₀+2.0 n₀ n₀+2.0 n₁ n₁")
basis=ls.SpinlessFermionBasis(2)
basis.build() #It is necessary  to build the basis before creating an operator
opr=ls.Operator(op_expr, basis) #Here we create the operator
```
Another way is to create the expression for each small expression and then add them to each other:
```pycon
op_expr = ls.Expr("-c†₀ c₁")
op_expr -= ls.Expr("c†₁ c₀")
op_expr += 2*ls.Expr("n₀ n₀")
op_expr += 2*ls.Expr("n₁ n₁")
opr=ls.Operator(op_expr, basis)
opr.shape() #We can check the shape of the operator
>>>
(4,4)
```

We can multiply the operator on a vector using standard Python notation.

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

The operators for the symmetry-adapted basis are created and manipulated in the same way.
We discuss this in more detail in the section [Symmetry-adapted basis](#Symmetry-adapted-basis).

### Symmetry

#### Symmetry as Permutations

The full power of `lattice_symmetries` manifests if one uses symmetries when constructing 
symmetry-adapted basis and linear operators acting on the corresponding Hilbert space. 
The symmetries are constructed with the help of expressions, and are represented as a permutation group of indices;
`lattice_symmetries` uses sympy to represent permutations, therefore one can take a look at [`sympy` documentation](https://docs.sympy.org/latest/modules/combinatorics/permutations.html) for more details.

Let's take a look at a couple of examples:
```pycon
p = ls.Permutation([1,2,3,0]) #This permutation shifts the indices, so that 0->1, 1->2, 2->3, 3->0
>>> (0 1 2 3) #The same permutation written in the cycle representation
```
```pycon
p = ls.Permutation([0,1,2,3]) #This is the identity permutation
>>> (3) #No cycles, therefore identity. 
#If the largest index doesn't move, it is shown in brackets; the format is used by sympy
```
```pycon
p=ls.Permutation([1,0,3,2]) #Exchange indices 0<->1 and 2<->3
>>> (0 1)(2 3) #Two cycles
```

We can use some of the `sympy` functions to analyze permutations.
For example, if one wants to find the order of a permutation:
```pycon
p = ls.Permutation([1,2,3,0])
p.order()
>>>
4
>>>
p**(p.order()) #check that p^4==1
>>>
Permutation(3)
```

#### Symmetry from Expressions

In the previous section, we constructed the symmetries by hand, however, this can only be done for relatively small and simple systems.
In most cases, the more powerful and convenient way is to ask `lattice_symmetries` to calculate the symmetries of a given expression. 

Since we need the characters to construct symmetry-adapted basis, there are two possibilities.

- One can find all permutation symmetries of an expression.
This option can be used to study specific sectors and does not cover the whole Hilbert space.

As a simple test case, we can consider a chain of 3 spins with the symmetry group $D_3=S_3$. This group is already non-abelian.

```pycon
import igraph as ig

e=ls.Expr("σ^x_0 σ^x_1") #A simple expression defined on an edge
expr=e.on(ig.Graph.Lattice(dim=[3], circular=True)) # The periodic chain with 3 sites
sym=expr.permutation_group()
>>>
PermutationGroup([
    (1 2),
    (2)(0 1),
    (0 1 2),
    (0 2 1),
    (0 2)]) 
#The identity element is always present, so it is not shown.
#The group has 6 elements, as it should
```

- Another option is to find the maximum abelian subgroup of the symmetry group (usually it coincides with the translation subgroup).
In this case, the sectors cover the whole Hilbert space.

For the same 3-site chain, we have:
```pycon
ab_sym=expr.abelian_permutation_group()
>>>
PermutationGroup([
    (0 1 2), #translation by 1
    (0 2 1)]) #translation by 2
#The maximal abelian sugroup consists of 2 translations and identity
```

These functions return a sympy permutation group. This allows us to apply all the existing `sympy` functions and analyze the symmetries of expressions.
For more details, one can look [at `sympy` documentation](https://docs.sympy.org/latest/modules/combinatorics/perm_groups.html).

#### Symmetry-adapted basis

To make calculations with the help of symmetries, we need to construct a symmetry-adapted basis.
The basis is constructed with the help of one-dimensional representations of the operator symmetry group,
or, in other words, the symmetry (permutation) group of the corresponding expression.

We describe the characters of the symmetry group as a list of tuples:
```pycon
[(permutation, rational number)]
```
where a permutation is a symmetry of the expression, and a rational number is its phase, defining the Hilbert space sector.
To be more precise, the phase is described as follows:

$$
\hat P |\psi\rangle=e^{2\pi i \phi}|\psi\rangle
$$

where $\hat P$ is the operator representing the permutation, $\psi$ is a vector from the chosen sector and $\phi$ is the rational number.

One can think about this list of tuples as a list of symmetry generators with corresponding phases.
However, the generators are not necessarily independent, and one can pass the whole symmetry group with the appropriate phases.
If one considers only abelian symmetries and independent generators, then the generator phases can be any rational number of the form $k/n$, 
where $n$ is the order of the generator, and $k$ is a natural number.
However, if the symmetry group is non-abelian, the phases should be chosen in a proper way,
so that characters of group elements organize a one-dimension representation of the symmetry group.

Now, we are ready to build a symmetry-adapted basis. 
There are two main ways to do that. The first way is to generate the representations by hand, and the second is to use the built-in functions of `lattice-symmetries`.

##### Manually generated representations

The simplest example would be:
```pycon
p = ls.Permutation([1,2,0]) #translation by 1
b0 = ls.SpinBasis(3, symmetries=[(p, ls.Rational(0, 1))]) #We specify the phase to be zero. It defines the trivial character
b0.build()
```

The order of the permutation is 3, wo we can have 2 other possible sectors:
```pycon
b1 = ls.SpinBasis(3, symmetries=[(p, ls.Rational(1, 3))]) #The phase is 1/3, the total phase is 2/3π
b1.build()

b2 = ls.SpinBasis(3, symmetries=[(p, ls.Rational(2, 3))]) #The phase is 2/3, the total phase is 4/3π
b2.build()
```
When making calculations, one needs to choose which sectors are suitable for the given task.

Let's consider another example, where we will work with the symmetries of expressions:
```pycon
e=ls.Expr("σ^x_0 σ^x_1")
expr=e.on(ig.Graph.Full(5)) # We define our lattice to be a complete graph with 5 vertices
sym=expr.permutation_group()
rep=[(p,ls.Rational(0,1)) for p in sym] #The trivial one-dimensional representation 
#Notice that we can use the whole symmetry group for constructing 1D-representation
basis = ls.SpinBasis(5, symmetries=rep) 
basis.build()
```

It should be noticed, that one can define a representation by any of its subset, and `lattice_symmetry` will generate the whole representation by the given generators.
Of course, the generators should be consistent with each other.

<a name="Sectors-section"></a>
##### Functions `hilbert_space_sectors` and `ground_state_sectors`

There is a possibility to generate the symmetry-adapted basis in `lattice_symmetries automatically.` One can find all sectors of the maximal abelian subgroup using `hilbert_space_sectors`,
or one can find all one-dimensional representations of the whole symmetry group by `ground_state_sectors`.
It is important to note that the function also calculate bases for different number of particles/hamming weight and spin inversion, so that the whole Hilbert space is covered. 
Let's consider an example of the most straightforward example of 3-sites chain:
```pycon
import igraph as ig

e=ls.Expr("σ^x_0 σ^x_1") #A simple expression defined on an edge
expr=e.on(ig.Graph.Lattice(dim=[3], circular=True))
```

We can find the sectors (symmetry adapted bases) using an `Expression`, as we did before.

```pycon
basis_list=expr.hilbert_space_sectors()
```

The result is the array of symmetrized bases for each of the possible representations of the maximal Abelian subgroup, hamming weight and spin inversion.
It should be noted that this list covers the whole Hilbert space and can be used for the exhaustive search of the relevant sectors.
One can check the corresponding representations as follows:
```
for basis in basis_list:
	#Print all possible sectors depending on translation sector, hamming weight and spin inversion
	print(basis.symmetries, basis.hamming_weight, basis.spin_inversion) 
>>>
[(Permutation(2), 0), (Permutation(0, 1, 2), 1/3), (Permutation(0, 2, 1), 2/3)] None -1
[(Permutation(2), 0), (Permutation(0, 1, 2), 2/3), (Permutation(0, 2, 1), 1/3)] None -1
[(Permutation(2), 0), (Permutation(0, 1, 2), 0), (Permutation(0, 2, 1), 0)] None -1
[(Permutation(2), 0), (Permutation(0, 1, 2), 1/3), (Permutation(0, 2, 1), 2/3)] None 1
[(Permutation(2), 0), (Permutation(0, 1, 2), 2/3), (Permutation(0, 2, 1), 1/3)] None 1
[(Permutation(2), 0), (Permutation(0, 1, 2), 0), (Permutation(0, 2, 1), 0)] None 1
```
We see that the hamming weight is not defined. However, the spin inversion is possible, and for each value of spin inversion we have three one-dimensional representations of the translation group.
We can also check the expression indeed does not conserve the number of particles, so it is impossible to define basis for the chosen hamming weight:

```pycon
e.conserves_number_particles
>>> False
```

Another useful function is `ground_state_sectors`:

```pycon
basis_list=expr.ground_state_sectors()
```
The result is the array of symmetrized bases for every one-dimensional representation of the whole symmetry group and possible values of particles/hamming weight and spin inversion.
The bases do not cover the whole Hilbert space, but they can be used for specific purposes. For example, if the ground state is unique, it for sure lies in one of the given sectors.

We can check, that the result is expected by printing the corresponding representations:
```
for basis in basis_list:
	#Print all possible sectors depending on translation sector, hamming weight and spin inversion
	print(basis.symmetries, basis.hamming_weight, basis.spin_inversion)
>>>
[(Permutation(2), 0), (Permutation(1, 2), 1/2), (Permutation(2)(0, 1), 1/2), (Permutation(0, 1, 2), 0), (Permutation(0, 2, 1), 0), (Permutation(0, 2), 1/2)] None -1
[(Permutation(2), 0), (Permutation(1, 2), 0), (Permutation(2)(0, 1), 0), (Permutation(0, 1, 2), 0), (Permutation(0, 2, 1), 0), (Permutation(0, 2), 0)] None -1
[(Permutation(2), 0), (Permutation(1, 2), 1/2), (Permutation(2)(0, 1), 1/2), (Permutation(0, 1, 2), 0), (Permutation(0, 2, 1), 0), (Permutation(0, 2), 1/2)] None 1
[(Permutation(2), 0), (Permutation(1, 2), 0), (Permutation(2)(0, 1), 0), (Permutation(0, 1, 2), 0), (Permutation(0, 2, 1), 0), (Permutation(0, 2), 0)] None 1
```
We have 2 representations with commuting parity and translations for each values of spin inversion.

It is important to notice that bases are not built yet. Therefore after choosing one of the bases, one should build it using `basis.build()` to make calculations.	

The last step is to make calculations in symmetry adapted basis. For that we need to create an `Operator` object.
When making an operator, one needs to be sure that the basis symmetries are consistent with the symmetries of the correspondning expression.
An operator for symmetry-adapted basis is built in the same way as without symmetries. Schematically, the code would look as follows:

```pycon
Make an expression
Calculate symmetries/choose symmetry-adapted basis
Build symmetry-adapted basis

opr=ls.Operator(expression,basis) # Make the operator in symmetry-adapted basis
```

## Examples

Here, we will take a look at different examples of `lattice_symmetries` applications.

### ED of 1D chain

We will start with the simplest example of exact diagonalization. We will consider Heisenberg chain on 10 sites and diagonalization with the help of symmetries.
For that, we will combine methods described in the previous sections.

First, we generate the expression for the Hamiltonian.

```pycon
import lattice_symmetries as ls
import numpy as np
import scipy

number_spins = 10  # System size

# Constructing the expression of the Hamiltonian
edges = [(i, (i + 1) % number_spins) for i in range(number_spins)]
expr = ls.Expr("2 (σ⁺₀ σ⁻₁ + σ⁺₁ σ⁻₀) + σᶻ₀ σᶻ₁", sites=edges)
print("Expression:", expr)

# Alternatively, we can create expr in an algebraic way:
two_site_expr = ls.Expr("2 (σ⁺₀ σ⁻₁ + σ⁺₁ σ⁻₀) + σᶻ₀ σᶻ₁")
many_exprs = [two_site_expr.replace_indices({0: i, 1: j}) for (i, j) in edges]
expr2 = reduce(operator.add, many_exprs)
assert expr == expr2 #we check that the expressions are equal
```

Since the lattice is simple, we can construct the symmetries by hand and use one-dimensional representation of the whole symmetry group.
It is worth noting that it is only possible, because we choose the proper sector for parity and translations permutations, which don't commute.
```pycon
# Constructing symmetries
sites = np.arange(number_spins)
# Momentum in x direction with phase φ=1/2 (i.e., with the eigenvalue λ=exp(-2πiφ)=-1)
T = (ls.Permutation((sites + 1) % number_spins), ls.Rational(1, 2))
# Parity with eigenvalue λ=-1
P = (ls.Permutation(sites[::-1]), ls.Rational(1, 2))
```

Now, we construct the spin basis with zero magnetization.
```pycon
# Constructing the basis
basis = ls.SpinBasis(
	number_spins=number_spins,
	# NOTE: we don't actually need to specify hamming_weight when spin_inversion
	# is set. The library will guess that hamming_weight = number_spins // 2.
	spin_inversion=-1,
	symmetries=[T, P], #Here, we use the characters of the whole symmetry group 	
)
basis.build()  # Build the list of representatives. We need it since we're doing ED
```
 The last step is to construct the actual Hamiltonian matrix and diagonalize it.
```pycon
# Construct the Hamiltonian
hamiltonian = ls.Operator(expr, basis)

# Diagonalize the Hamiltonian using ARPACK. The fast matrix-vector multiplication of `lattice-symmetries` will be used,
#because hamiltonian is a lattice-symmetry operator
eigenvalues, eigenstates = scipy.sparse.linalg.eigsh(hamiltonian, k=1, which="SA")
print("Ground state energy is {}".format(eigenvalues[0]))
assert np.isclose(eigenvalues[0], -18.06178542)
```
Above, we used the symmetry-adapted basis, and therefore we were restricted to the specific sector of the total Hilbert space.
Since the considered system is relatively small, we can make the diagonalization in the full basis and check that we chose the right sector.

```pycon
new_basis = ls.SpinBasis(
	number_spins=number_spins,
	spin_inversion=-1,
)
new_basis.build()
new_hamiltonian = ls.Operator(expr, new_basis)
eigenvalues, eigenstates = scipy.sparse.linalg.eigsh(hamiltonian, k=1, which="SA")
print("Ground state energy is {}".format(eigenvalues[0]))
assert np.isclose(eigenvalues[0], -18.06178542)
```
Indeed, the outcome is the same within machine precision. 


### ED of Star graph

In the previous example, we considered a 1D chain, where we knew the graph's symmetries and could code them by hand.
Sometimes, that task requires much time or is even practically impossible.
Nevertheless, symmetries can still help, especially when the required sector lies in the one-dimensional representation of the whole symmetry group.
As an example, we will consider the Heisenberg model on a star graph, where every vertex connects only to the center:
 
```pycon
import lattice_symmetries as ls
import numpy as np
import scipy
import igraph as ig

number_spins=7
expr = ls.Expr("2 (σ⁺₀ σ⁻₁ + σ⁺₁ σ⁻₀) + σᶻ₀ σᶻ₁").on(ig.Graph.Star(number_spins))
sym=expr.permutation_group()
rep=[(p,ls.Rational(0,1)) for p in sym] #The trivial one-dimensional representation 
basis = ls.SpinBasis(number_spins, symmetries=rep) #We don't choose specific magnetization
basis.build()
hamiltonian = ls.Operator(expr, basis)
eigenvalues, eigenstates = scipy.sparse.linalg.eigsh(hamiltonian, k=1, which="SA")
print("Ground state energy is {}".format(eigenvalues[0]))
assert np.isclose(eigenvalues[0], -8)
```
Since the system is small, we can check that we chose the right sector:

```pycon
basis = ls.SpinBasis(number_spins) #No symmetries
basis.build()
hamiltonian = ls.Operator(expr, basis)
eigenvalues, eigenstates = scipy.sparse.linalg.eigsh(hamiltonian, k=1, which="SA")
print("Ground state energy is {}".format(eigenvalues[0]))
assert np.isclose(eigenvalues[0], -8)
```

Of course, the considered example is somewhat artificial, nevertheless, it nicely demonstrates one of the possible applications of `lattice_symmetries`.

### Time evolution

Another example of the capabilities of `lattice_symmetries` is time evolution.
To apply time evolution, we will use the Chebyshev expansion of the matrix exponent:

$$
e^{iHt}=e^{-i\frac{bt}{a}}[J_0(at)+\sum\limits^{\infty} 2(-i)^n J_n(at)T_n(H')]
$$

where we rescale the original Hamiltonian $H$ with bandwidth $[E_g, E_s]$ to the Hamiltonian $H'$ with bandwidth $[-1+\epsilon,1-\epsilon]$, so that the series converges.
The rescaling takes the form:

$$
H'=\frac{H-b}{a}
$$

where $a=1/(2-\epsilon)$, $b=(E_g+E_s)/2$, and $\epsilon$ is a safety factor to make the spectrum of $H'$ be within $[-1,1]$.

We will consider the basic version of the algorithm, where instead of usual `numpy` arrays we will use `lattice_symmetries` operators.

### Spin structure factor

We can use `lattice-symmetries` to calculate static spin structure factor:

$$
S({\bf q})=\frac{1}{N}\sum e^{i{\bf q}\cdot ({\bf r_i}-{\bf r_j})}\langle S_i S_j \rangle
$$

The idea is to use the symmetries of the system to make calculations faster.