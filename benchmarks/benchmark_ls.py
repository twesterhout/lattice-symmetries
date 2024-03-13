import igraph as ig
import lattice_symmetries as ls
import numpy as np
import timeit


def main():
    n = 25
    b = ls.SpinBasis(n, hamming_weight=n // 2)
    b.build()
    e = ls.Expr(
        "2 (σ⁺₀ σ⁻₁ + σ⁺₁ σ⁻₀) + σᶻ₀ σᶻ₁", sites=ig.Graph.Ring(n=n, circular=True)
    )
    h = ls.Operator(e, b)
    x = np.random.rand(b.number_states)
    y = np.empty_like(x)
    print(timeit.repeat(lambda: h.apply_to_state_vector(x, out=y), repeat=4, number=1))


if __name__ == "__main__":
    main()
