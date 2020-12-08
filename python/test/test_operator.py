import numpy as np
import lattice_symmetries as ls


def test_apply():
    basis = ls.SpinBasis(ls.Group([]), number_spins=10, hamming_weight=5)
    basis.build()
    # fmt: off
    matrix = np.array([[1,  0,  0, 0],
                       [0, -1,  2, 0],
                       [0,  2, -1, 0],
                       [0,  0,  0, 1]])
    # fmt: on
    edges = [(i, (i + 1) % basis.number_spins) for i in range(basis.number_spins)]
    operator = ls.Operator(basis, [ls.Interaction(matrix, edges)])

    for (i, spin) in enumerate(basis.states):
        spins, coeffs = operator.apply(spin)
        v = np.zeros(basis.number_states, dtype=np.float64)
        v[i] = 1.0
        v = operator(v)
        for (s, c) in zip(spins, coeffs):
            assert v[basis.index(s[0])] == c

test_apply()
