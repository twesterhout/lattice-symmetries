import lattice_symmetries as ls
from loguru import logger
import numpy as np
import timeit
from systems import square_lattice_symmetries, make_basis


def reference_cast_to_basis(dest_basis, src_basis, src_state):
    dest_state = np.zeros_like(src_state, shape=(dest_basis.number_states,))
    use_real = dest_state.dtype == np.float32 or dest_state.dtype == np.float64
    for (i, x) in enumerate(dest_basis.states):
        _, _, dest_norm = dest_basis.state_info(x)
        r, c, src_norm = src_basis.state_info(x)
        if use_real:
            assert np.isclose(c.imag, 0.0)
            c = c.real
        if src_norm != 0:
            index = src_basis.index(r)
            dest_state[i] = src_state[index] * src_norm * c / dest_norm
    return dest_state


def cast_to_basis(dest_basis, src_basis, src_state, batch_size=10000):
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


def benchmark_cast_to_basis():
    rng = np.random.default_rng()
    for L_y, L_x in [(4, 5)]:
        logger.info("Benchmarking for {}x{}...", L_y, L_x)
        symmetries = square_lattice_symmetries(L_x, L_y)

        basis1 = make_basis(
            symmetries, backend="ls", number_spins=L_x * L_y, hamming_weight=(L_x * L_y) // 2
        )
        x1 = rng.standard_normal(size=basis1.number_states, dtype=np.float64)

        basis2 = make_basis(
            [], backend="ls", number_spins=L_x * L_y, hamming_weight=(L_x * L_y) // 2
        )

        logger.info("Benchmarking reference_cast_to_basis...")
        reference_times = timeit.repeat(
            lambda: reference_cast_to_basis(basis2, basis1, x1), repeat=3, number=1
        )
        logger.info("Reference: {}", reference_times)
        logger.info("Benchmarking cast_to_basis...")
        batched_times = timeit.repeat(
            lambda: cast_to_basis(basis2, basis1, x1), repeat=3, number=1
        )
        logger.info("Numba:     {}", batched_times)


if __name__ == "__main__":
    benchmark_cast_to_basis()
