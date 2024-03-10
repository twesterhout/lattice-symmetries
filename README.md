> ⚠️ **INFO**
>
> This is a Haskell rewrite of the original
> [lattice-symmetries](https://github.com/twesterhout/lattice-symmetries/v1). At
> some point, this package will completely replace the first version of
> lattice-symmetries.

# lattice-symmetries [![Build](https://github.com/twesterhout/lattice-symmetries-haskell/actions/workflows/ci.yml/badge.svg)](https://github.com/twesterhout/lattice-symmetries-haskell/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

A package to simplify working with symmetry-adapted quantum many-body bases.

## Permutations




## Expressions

In lattice-symmetries, we work with in the second quantization formalism. This
means that you can use primitive operators such as $\sigma^x$, $\sigma^y$, and
$\sigma^z$ to build expressions for your Hamiltonian and observables.

### Primitive operators

 - $\sigma^x$ and $S^x$:
   ```pycon
   >>> Expr("σˣ₀") == Expr("\\sigma^x_0")
   True
   >>> Expr("Sˣ₀") == 0.5 * Expr("σˣ₀")
   True
   >>> Expr("σˣ₀").to_dense()
   [[0, 1],
    [1, 0]]
   ```

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

 - $\mathbb{1}$:
   ```pycon
   >>> Expr("I", particle="spin-1/2").to_dense()
   [[1, 0],
    [0, 1]]
   ```
   (*Note:* we need to explicitly specify the particle type because it cannot be deduced from the expression)

> ⚠️ **Note:** In all expressions above, we used the index 0 as en example, but
> you can also construct primitive operators residing on other lattice sites.


### Operator algebra

Primitive can be combined using the `+`, `-`, and `*` operations to build more complex operators. Furthermore, expressions can also be multiplied by scalars from the left using the `*` operator.

For example, here are a few ways to write down the Heisenberg interaction

$$
\mathbf{S}_i \cdot \mathbf{S}_j = S^x_i S^x_j + S^y_i S^y_j + S^z_i S^z_j
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


### More sites

At some point


