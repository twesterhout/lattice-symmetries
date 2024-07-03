from typing import Literal, overload

import numba
import numpy as np
import numpy.typing as npt
import scipy.special
from numba import float32, float64, uint64

from . import Basis


def index_to_spin(index: npt.NDArray[np.uint64], number_spins: int) -> npt.NDArray[np.bool_]:
    return (
        (index.reshape(-1, 1).astype(np.int64) & (1 << np.arange(number_spins).astype(np.int64)))
    ) > 0


def spin_to_index(
    spin: npt.NDArray[np.bool_ | np.int_ | np.uint], number_spins: int
) -> npt.NDArray[np.uint64]:
    a = 2 ** np.arange(number_spins, dtype=np.int64)
    return spin.dot(a)


@numba.jit("uint64(uint64, uint64)", nogil=True, nopython=True)
def _binom(n, k):
    r"""Compute the number of ways to choose k elements out of a pile of n.

    :param n: the size of the pile of elements
    :param k: the number of elements to take from the pile
    :return: the number of ways to choose k elements out of a pile of n
    """
    assert 0 <= n and n < 40
    assert 0 <= k and k <= n

    if k == 0 or k == n:
        return 1
    total_ways = 1
    for i in range(min(k, n - k)):
        total_ways = total_ways * (n - i) // (i + 1)
    return total_ways


@numba.jit("uint64[:, :](uint64, uint64)", nogil=True, nopython=True)
def make_binom_cache(max_n, max_k):
    assert 0 < max_k and max_k <= max_n
    assert max_n < 40

    cache = np.zeros((max_n + 1, max_k + 2), dtype=np.uint64)
    for i in range(max_n + 1):
        for j in range(min(i, max_k) + 2):
            cache[i, j] = _binom(i, j) if j <= i else 0
    return cache


@numba.jit("uint64(uint64)", nogil=True, nopython=True)
def _hamming_weight(n: int):
    # See https://stackoverflow.com/a/9830282
    n = (n & 0x5555555555555555) + ((n & 0xAAAAAAAAAAAAAAAA) >> 1)
    n = (n & 0x3333333333333333) + ((n & 0xCCCCCCCCCCCCCCCC) >> 2)
    n = (n & 0x0F0F0F0F0F0F0F0F) + ((n & 0xF0F0F0F0F0F0F0F0) >> 4)
    n = (n & 0x00FF00FF00FF00FF) + ((n & 0xFF00FF00FF00FF00) >> 8)
    n = (n & 0x0000FFFF0000FFFF) + ((n & 0xFFFF0000FFFF0000) >> 16)
    n = (n & 0x00000000FFFFFFFF) + ((n & 0xFFFFFFFF00000000) >> 32)
    return n


@numba.jit("uint64[:](uint64, uint64)", nogil=True, nopython=True)
def generate_binaries(length: int, hamming: int) -> npt.NDArray[np.uint64]:
    assert length > 0
    assert hamming >= 0 and hamming <= length
    if hamming == 0:
        return np.zeros(1, dtype=np.uint64)

    size = _binom(length, hamming)
    set_of_vectors = np.empty(size, dtype=np.uint64)
    val = (1 << hamming) - 1
    for i in range(size):
        set_of_vectors[i] = val
        c = val & -val
        r = val + c
        val = (((r ^ val) >> 2) // c) | r
    return set_of_vectors


@numba.jit("uint64(uint64, uint64, uint64)", nogil=True, nopython=True)
def _merge(vec1: int, vec2: int, dim2: int) -> int:
    """
    (vec1 << dim2) | vec2
    """
    return (vec1 << dim2) | vec2


@numba.jit("uint64(uint64, uint64, uint64[:, :])", nogil=True, nopython=True)
def _get_index(x: int, dim: int, binom_cache: npt.NDArray[np.uint64]) -> int:
    r"""Compute the lexicographic index of a binary number among those with
    the same Hamming weight.

    :param x: the binary number to find the index for
    :param dim: the number of bits in the binary number
    :param binom_cache: a 2D array of precomputed binomial coefficients
    :return: the lexicographic index of the input binary number

    This function uses a precomputed binomial coefficient cache for efficiency.
    It is optimized with Numba's just-in-time compilation and assumes x is a
    valid binary number with dim bits.
    """

    n = 0
    h = 0
    if x & 1 == 1:
        h += 1
    x >>= 1
    for i in range(1, dim):
        if x & 1 == 1:
            h += 1
            n += binom_cache[i, h]
        x >>= 1
    return n


@numba.jit(
    "complex128(uint64, uint64, uint64, uint64, uint64, complex128[::1], uint64[:, :])",
    nogil=True,
    nopython=True,
)
def _density_matrix_element(
    sector_dim: int,
    dim: int,
    hamming: int,
    vec1: int,
    vec2: int,
    amplitudes: npt.NDArray[np.complex128],
    binom_cache: npt.NDArray[np.uint64],
) -> complex:
    sector_hamming = _hamming_weight(vec1)
    assert sector_hamming == _hamming_weight(vec2)

    smallest_number = 2 ** (hamming - sector_hamming) - 1
    k = dim - sector_dim
    index1 = _get_index(_merge(vec1, smallest_number, k), dim, binom_cache)
    index2 = _get_index(_merge(vec2, smallest_number, k), dim, binom_cache)
    size_of_traced = _binom(k, hamming - sector_hamming)
    matrix_element = np.dot(
        amplitudes[index1 : index1 + size_of_traced].conj(),
        amplitudes[index2 : index2 + size_of_traced],
    )

    return matrix_element


@numba.jit(
    "Tuple((complex128[:, :], uint64[:]))(uint64, uint64, uint64, uint64, complex128[::1])",
    nogil=True,
    nopython=True,
    parallel=False,
)
def _sector_density_matrix(
    sector_dim: int,
    dim: int,
    sector_hamming: int,
    hamming: int,
    amplitudes: npt.NDArray[np.complex128],
) -> tuple[npt.NDArray[np.complex128], npt.NDArray[np.uint64]]:
    assert sector_hamming <= hamming and sector_dim <= dim
    assert hamming - sector_hamming <= dim - sector_dim

    sector_basis = generate_binaries(sector_dim, sector_hamming)
    binom_cache = make_binom_cache(dim, hamming)
    n = len(sector_basis)
    matrix = np.empty((n, n), dtype=np.complex128)
    for i in range(len(sector_basis)):
        for j in range(i, len(sector_basis)):
            matrix[i, j] = _density_matrix_element(
                sector_dim,
                dim,
                hamming,
                sector_basis[i],
                sector_basis[j],
                amplitudes,
                binom_cache,
            )
            matrix[j, i] = np.conj(matrix[i, j])

    return matrix, sector_basis


def _density_matrix_blocks(
    sub_dim: int, dim: int, hamming: int, amplitudes: npt.NDArray[np.complex128]
) -> tuple[list[npt.NDArray[np.complex128]], list[npt.NDArray[np.uint64]]]:
    matrices = []
    sector_basis_states = []
    for sub_hamming in range(max(0, hamming - (dim - sub_dim)), min(hamming, sub_dim) + 1):
        matrix, sector_basis = _sector_density_matrix(
            sub_dim, dim, sub_hamming, hamming, amplitudes
        )
        matrices.append(matrix)
        sector_basis_states.append(sector_basis)

    return matrices, sector_basis_states


def _move_sites_to_back(
    wavefunction: npt.NDArray[np.complex128], basis: ls.Basis, sites: list[int]
) -> npt.NDArray[np.complex128]:
    unpacked_states = index_to_spin(basis.states, basis.number_sites)
    permuted_sites = sorted(set(range(basis.number_sites)) - set(sites)) + sites
    assert set(permuted_sites) == set(range(basis.number_sites))
    assert len(permuted_sites) == basis.number_sites
    permuted_unpacked_states = unpacked_states[:, np.argsort(permuted_sites)]
    permuted_states = spin_to_index(permuted_unpacked_states, basis.number_sites)
    return wavefunction[basis.index(permuted_states)]


def _reduced_density_matrix_blocks(
    wavefunction: npt.NDArray[np.complex128], basis: ls.Basis, sites: list[int]
) -> tuple[list[npt.NDArray[np.complex128]], list[npt.NDArray[np.uint64]]]:
    if basis.hamming_weight is None:
        raise ValueError("Only bases with a fixed hamming weight are supported.")
    if basis.symmetries != []:
        raise ValueError("Symmetries are not supported (yet)")
    if basis.spin_inversion is not None:
        raise ValueError("Spin inversion is not supported (yet)")

    hamming_weight = basis.hamming_weight
    n_sites = basis.number_sites
    wavefunction = _move_sites_to_back(wavefunction, basis, sites)
    return _density_matrix_blocks(len(sites), n_sites, hamming_weight, wavefunction)


@overload
def reduced_density_matrix(
    wavefunction: npt.NDArray[np.complex128],
    basis: Basis,
    sites: list[int],
    return_states: Literal[True] = True,
) -> tuple[npt.NDArray[np.complex128], npt.NDArray[np.uint64]]: ...


@overload
def reduced_density_matrix(
    wavefunction: npt.NDArray[np.complex128],
    basis: Basis,
    sites: list[int],
    return_states: Literal[False] = True,
) -> npt.NDArray[np.complex128]: ...


def reduced_density_matrix(
    wavefunction: npt.NDArray[np.complex128],
    basis: Basis,
    sites: list[int],
    return_states: bool = True,
) -> npt.NDArray[np.complex128] | tuple[npt.NDArray[np.complex128], npt.NDArray[np.uint64]]:
    matrices, sector_basis_states = _reduced_density_matrix_blocks(wavefunction, basis, sites)
    if return_states:
        return scipy.linalg.block_diag(*matrices), np.concatenate(sector_basis_states)

    all_states = np.arange(2 ** len(sites), dtype=np.uint64)
    full_matrix = np.zeros((2 ** len(sites), 2 ** len(sites)), dtype=np.complex128)
    for matrix, states in zip(matrices, sector_basis_states):
        state_idxs = np.searchsorted(all_states, states)
        full_matrix[np.ix_(state_idxs, state_idxs)] = matrix
    return full_matrix
