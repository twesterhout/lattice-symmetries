# lattice_symmetries tutorial

## Contents

-[Installing] (#installing)
-[Introduction] (#introduction)
-[Examples] (#examples)
-[ED] (#ed)

## Installing

The installation of lattice_symmetries requires Nix. 
Therefore, at first you need to install Nix. Then you need to fork the last version from guthub.
The last step is to configure Nix, so that every depencence works nicely. Voila, and you have a working package.

## Introduction

The `lattice_symmetries` is a powerful package that nicely works with many-body quantum spin systems
and takes into account system symmetries.
The `lattice_symmetries` implements fast matrix-vector multiplication that can be applied to various problems, 
such as time evolution or exact diagonalization of many-body Hamiltonians. 

## Examples

Here we will take a look at different examples of lattice_symmetries applications.

### ED

Now we will consider the simplest example of exact diagonalization.
