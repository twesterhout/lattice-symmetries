import numpy as np
import lattice_symmetries as ls

ls.enable_logging()

import systems


def test_4_spins():
    # fmt: off
    matrix = np.array([[1,  0,  0, 0],
                       [0, -1,  2, 0],
                       [0,  2, -1, 0],
                       [0,  0,  0, 1]])
    # fmt: on
    number_spins = 4
    edges = [(i, (i + 1) % number_spins) for i in range(number_spins)]

    basis = ls.SpinBasis(ls.Group([]), number_spins=4, hamming_weight=2)
    basis.build()
    assert basis.number_states == 6
    operator = ls.Operator(basis, [ls.Interaction(matrix, edges)])
    assert np.isclose(ls.diagonalize(operator, k=1)[0], -8)

    basis = ls.SpinBasis(ls.Group([]), number_spins=4, hamming_weight=2, spin_inversion=1)
    basis.build()
    assert basis.number_states == 3
    operator = ls.Operator(basis, [ls.Interaction(matrix, edges)])
    assert np.isclose(ls.diagonalize(operator, k=1)[0], -8)

    T = ls.Symmetry([1, 2, 3, 0], sector=0)
    basis = ls.SpinBasis(ls.Group([T]), number_spins=4, hamming_weight=2, spin_inversion=1)
    basis.build()
    assert basis.number_states == 2
    operator = ls.Operator(basis, [ls.Interaction(matrix, edges)])
    assert np.isclose(ls.diagonalize(operator, k=1)[0], -8)


def test_index():
    L_x, L_y = (4, 6)
    backend = "ls"
    symmetries = systems.square_lattice_symmetries(L_x, L_y)
    nearest, _ = systems.square_lattice_edges(L_x, L_y)
    basis = systems.make_basis(
        symmetries,
        backend=backend,
        number_spins=L_x * L_y,
        hamming_weight=(L_x * L_y) // 2,
    )
    # print(basis.number_states)
    hamiltonian = systems.make_heisenberg(basis, nearest, backend=backend)

    # indices = ls.batched_index(basis, basis.states)
    # assert np.all(indices == np.arange(basis.number_states, dtype=np.uint64))
    for i in range(basis.number_states):
        index = basis.index(basis.states[i])
        assert index == i

    assert np.all(basis.batched_index(basis.states) == np.arange(basis.number_states))

    spins = np.zeros((10000, 8), dtype=np.uint64)
    spins[:, 0] = basis.states[:10000]
    ls.batched_state_info(basis, spins)
    # evals, evecs = hamiltonian.eigsh(k=1, which='SA')
    # evals, evecs = ls.diagonalize(hamiltonian)
    # print(evals)


test_index()
# test_4_spins()
