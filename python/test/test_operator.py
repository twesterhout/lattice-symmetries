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


def test_non_hermitian_matvec():
    basis = ls.SpinBasis(ls.Group([]), number_spins=2)
    basis.build()

    rng = np.random.default_rng()
    matrix = rng.random((4, 4)) + rng.random((4, 4)) * 1j - (0.5 + 0.5j)
    operator = ls.Operator(basis, [ls.Interaction(matrix, [(0, 1)])])

    for (i, spin) in enumerate(basis.states):
        v = np.zeros(basis.number_states, dtype=np.complex128)
        v[i] = 1.0
        predicted = operator(v)
        expected = matrix[[0, 2, 1, 3], :][:, [0, 2, 1, 3]] @ v
        # print(predicted, expected)
        assert np.allclose(predicted, expected)


def test_index():
    basis = ls.SpinBasis(ls.Group([]), number_spins=4, hamming_weight=2)
    basis.build()

    states = basis.states
    print(states)
    print(ls.batched_index(basis, states))
    # print(ls.batched_index(basis, np.array([3, 9, 10, 12, 8], dtype=np.uint64)))
    # for i in range(len(states)):
    #     assert i == ls._numba_index(basis._payload.value, states[i])

    bits = np.zeros((len(states), 8), dtype=np.uint64)
    bits[:, 0] = states
    print(ls.batched_state_info(basis, bits))
    print(basis.state_info(states[1]))

test_non_hermitian_matvec()
