**WORK IN PROGRESS**


# Lattice Symmetries

A package to simplify working with symmetrized quantum many-body bases (think
spin systems). This package is written with two main applications in mind:

  * Exact diagonalization;
  * Experiments with neural quantum states and symmetries.


## Reading the source code

Doxygen documentation focuses on the API (Application Programming Interface) and
the algorithms which made this code work fast enough may seem like black magic.
The best way to learn about them is to read the source code. To simplify the
process we provide a short description for each file:

* `error_handling.hpp`/`error_handling.cpp` teaches C++11 `<system_error>` (and
hence Outcome) about our custom error code type (`ls_error_code`). It also
defines `LATTICE_SYMMETRIES_CHECK` and `LATTICE_SYMMETRIES_ASSERT` macros which
are used to sanity checking (i.e. *not* for validating user input!).

* `bits.hpp` defines `bits512` type which is just a contiguous array of 512
bits (one cache line on most modern processors). It is used to represent spin
configurations when system size exceeds 64.

* `permutation.hpp`/`permutation.cpp` implements the compilation of a
permutation to a Benes network. The Benes network itself is not optimized yet at
this point (that happens in `network.hpp`).

* `network.hpp`/`network.cpp` takes a network produced by `permutation.hpp` and
depending on the system size converts it to a more compact representation. For
this representation we implement high-performance kernels to propagate spin
configurations efficiently through the network.
