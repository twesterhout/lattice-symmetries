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




## Hamiltonians

#### Spins

<table>
<tr><th>Maths</th><th>Code</th></tr>
<tr>
<td>

$$
\mathbf{S}_i \cdot \mathbf{S}_j = S^x_i S^x_j + S^y_i S^y_j + S^z_i S^z_j
$$

</td>
<td>

`"Sˣ₀ Sˣ₁ + Sʸ₀ Sʸ₁ + Sᶻ₀ Sᶻ₁"`
or
`"Sx0 Sx0 + Sy1 Sy1 + Sz0 Sz1"`

</td>
</tr>
<tr>
<td>

$$
\mathbf{S}_i \cdot \mathbf{S}_j = \frac{1}{4} \left( \sigma^x_i \sigma^x_j + \sigma^y_i \sigma^y_j + \sigma^z_i \sigma^z_j \right)
$$

</td>
<td>

`"0.25 (σˣ₀ σˣ₁ + σʸ₀ σʸ₁ + σᶻ₀ σᶻ₁)"`

</td>
</tr>
<tr>
<td>

$$
\sigma^{+}_i \sigma^{-}_j
$$

</td>
<td>

`"σ⁺₀ σ⁻₁"` or
`"\sigma^+_0 \sigma^-_1"` or
`"\sigma+0 \sigma-1"`

</td>
</tr>
</table>

#### Electrons

<table>
<tr><th>Maths</th><th>Code</th></tr>
<tr><td>

$$
c^\dagger_{i\uparrow}c_{j\uparrow} + c^\dagger_{i\downarrow}c_{j\downarrow}
$$

</td><td>

`"c†₀↑ c₁↑ + c†₀↓ c₁↓"`

</td></tr>
<tr><td>

$$
n_{i\uparrow} n_{i\downarrow}
$$

</td><td>

`"n₀↑ n₀↓"` or
`"n0up n0down"`

</td></tr>
</table>
