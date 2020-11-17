# Lattice Symmetries
![Linux](https://github.com/twesterhout/lattice-symmetries/workflows/Ubuntu/badge.svg)
![OS X](https://github.com/twesterhout/lattice-symmetries/workflows/MacOS/badge.svg)

A package to simplify working with symmetry-adapted quantum many-body bases
(think spin systems). This package is written with two main applications in
mind:

  * Exact diagonalization;
  * Experiments with neural quantum states and symmetries.

`lattice_symmetries` provides a relatively low-level (and high performance)
interface to working with symmetries and Hilbert space bases and operators. If
all you want to do it to diagonalize a spin Hamiltonian, have a loop at
[`SpinED`](...) application which does all the heavy-lifting for you.

Have a look at [online documentation](...) for installation instructions and
usage examples.


## Why don't you just use QuSpin?

[QuSpin](...) is a very general Python package for working with many-body
systems using exact-diagonalization-like methods. I am a big fan! However, it is
very tightly coupled to Python and could not be easily incorporated into
variational Monte Carlo code I was working on. This was my main reason for
partially reimplementing QuSpin's functionality. Hence, this package provides a
portable C interface with minimal dependencies.

It turns out, however, that this package is also considerably faster than QuSpin
and can thus handle larger systems. For systems larger than 32 spins, storing
the Hamiltonian as a matrix (even sparse) is impractical. However, computing
matrix elements on the fly requires that symmetry applications are very
efficient, and the implementation in QuSpin is not quite there yet.


## Citing

If you are using this package in your research, please, consider citing the
following paper (**WIP**):
```
@unpublished{citekey,
  author = "",
  title  = "",
  note   = "",
  month  = "",
  year   = "",
  annote = "".
}
```


## Okay, but how does it work?

Doxygen documentation focuses on the API (Application Programming Interface) and
the algorithms which make this code work fast enough might seem like black magic.
The best way to learn about them is to read the source code. To simplify the
process we provide a short guide.

The code can be roughly divided into a few layers. At the very bottom we have
error handling (`error_handling.*` files) and utilities (`bits512.hpp`,
`intrusive_ptr.hpp` files). Next, we have the permutation-to-Benes-network
compiler (`permutation.*` files). Intermediate representation if these networks
are then converted into more compact form optimized for forward propagation
(`network.*` files define the representation and `kernel.*` files implement the
 forward propagation). Using Benes networks we build symmetries (`symmetry.*`).
They keep track of the characters of the symmetry group, allow to find
representative vectors, etc. Given a bunch of symmetry operators, we construct
the symmetry group (`group.*` files). Finally, this symmetry group is used to
construct a Hilbert space basis (`basis.*` files). For small systems, one can
explicitly construct a list of all representative vectors of the basis. We call
it "cache" (`cache.*` files). It allows one to determine the Hilbert space
dimension, assigns indices to representative vectors etc.

Having a Hilbert space basis is typically not enough. We need to define a
Hamiltonian. All operators in `lattice_symmetries` are assumed to be Hermitian
and are defined as sums of "interaction" terms (`operator.*` files). These are
just small operators defined on 1-, 2-, 3-, or 4-spin subsystems. Operators can
then be applied to explicit vectors (e.g. for exact diagonalization) of in a CPS
(continuation-passing-style) for more general applications (e.g. neural quantum
states).
