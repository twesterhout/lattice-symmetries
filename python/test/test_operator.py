import numpy as np
import lattice_symmetries as ls

ls.enable_logging()


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


def test_batched_apply():
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

    for batch_size in [1, 2, 10, 20, 50]:
        x = basis.states[11 : 11 + batch_size]  # just a random range
        spins, coeffs, counts = operator.batched_apply(x)
        assert len(counts) == x.shape[0]
        offset = 0
        for i in range(batch_size):
            expected_spins, expected_coeffs = operator.apply(x[i])
            assert np.all(spins[offset : offset + counts[i]] == expected_spins)
            assert np.all(coeffs[offset : offset + counts[i]] == expected_coeffs)
            offset += counts[i]
        assert offset == spins.shape[0]
        assert offset == coeffs.shape[0]


def test_non_hermitian_matvec():
    basis = ls.SpinBasis(ls.Group([]), number_spins=2)
    basis.build()

    rng = np.random.default_rng()
    matrix = rng.random((4, 4)) + rng.random((4, 4)) * 1j - (0.5 + 0.5j)
    operator = ls.Operator(basis, [ls.Interaction(matrix, [(0, 1)])])

    # print(operator.to_csr().toarray())
    # print(matrix[[0, 2, 1, 3], :][:, [0, 2, 1, 3]])
    assert np.all(operator.to_csr().toarray() == matrix[[0, 2, 1, 3], :][:, [0, 2, 1, 3]])

    for (i, spin) in enumerate(basis.states):
        v = np.zeros(basis.number_states, dtype=np.complex128)
        v[i] = 1.0
        predicted = operator(v)
        expected = matrix[[0, 2, 1, 3], :][:, [0, 2, 1, 3]] @ v
        # print(predicted, expected)
        assert np.allclose(predicted, expected)

    v = rng.random(4) + rng.random(4) * 1j - (0.5 + 0.5j)
    predicted = operator(v)
    expected = matrix[[0, 2, 1, 3], :][:, [0, 2, 1, 3]] @ v
    # print(predicted)
    # print(expected)
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


# test_batched_apply()
test_non_hermitian_matvec()
