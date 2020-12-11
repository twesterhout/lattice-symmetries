import lattice_symmetries as ls
import numpy as np
import numba


def reference_cast_to_basis(dest_basis, src_basis, src_state):
    dest_state = np.zeros_like(src_state, shape=(dest_basis.number_states,))
    for (i, x) in enumerate(dest_basis.states):
        _, _, dest_norm = dest_basis.state_info(x)
        r, c, src_norm = src_basis.state_info(x)
        if src_norm != 0:
            index = src_basis.index(r)
            dest_state[i] = src_state[index] * src_norm * c.real / dest_norm
    assert np.isclose(np.linalg.norm(dest_state), 1.0)
    return dest_state


def cast_to_basis(dest_basis, src_basis, src_state, batch_size=100):
    dest_state = np.zeros_like(src_state, shape=(dest_basis.number_states,))
    use_real = dest_state.dtype == np.float32 or dest_state.dtype == np.float64

    def task(i, j):
        x = dest_basis.states[i:j]
        _, _, dest_norm = ls.batched_state_info(dest_basis, x)
        r, c, src_norm = ls.batched_state_info(src_basis, x)
        mask = src_norm != 0
        dest_norm = dest_norm[mask]
        r = r[:, 0][mask]
        c = c[mask]
        src_norm = src_norm[mask]

        if use_real:
            assert np.allclose(c.imag, 0.0)
            c = c.real
        index = ls.batched_index(src_basis, r)
        dest_state[i:j][mask] = src_state[index] * src_norm * c.real / dest_norm

    i = 0
    while i + batch_size <= dest_basis.number_states:
        task(i, i + batch_size)
        i += batch_size
    if i != dest_basis.number_states:
        task(i, dest_basis.number_states)
    return dest_state


def test_cast_to_basis():
    T = ls.Symmetry([1, 2, 3, 4, 5, 6, 7, 8, 9, 0], sector=5)
    P = ls.Symmetry([9, 8, 7, 6, 5, 4, 3, 2, 1, 0], sector=1)
    basis1 = ls.SpinBasis(ls.Group([T]), number_spins=10, hamming_weight=5, spin_inversion=-1)
    basis1.build()

    matrix = np.array([[1, 0, 0, 0], [0, -1, 2, 0], [0, 2, -1, 0], [0, 0, 0, 1]])
    edges = [(i, (i + 1) % basis1.number_spins) for i in range(basis1.number_spins)]
    operator1 = ls.Operator(basis1, [ls.Interaction(matrix, edges)])
    E1, x1 = ls.diagonalize(operator1)
    x1 = x1.squeeze()

    basis2 = ls.SpinBasis(ls.Group([P]), number_spins=10, hamming_weight=5)
    basis2.build()
    operator2 = ls.Operator(basis2, [ls.Interaction(matrix, edges)])

    E2, x2 = ls.diagonalize(operator2)
    x2 = x2.squeeze()
    assert np.isclose(E1, E2)

    y = reference_cast_to_basis(basis1, basis2, x2)
    y2 = cast_to_basis(basis1, basis2, x2)
    assert np.allclose(y, y2)
    assert np.isclose(np.abs(np.dot(x1, y)), 1.0)


test_cast_to_basis()
