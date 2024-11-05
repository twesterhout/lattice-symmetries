import lattice_symmetries as ls
import lattice_symmetries._kernels
import sympy
import time
import math
from loguru import logger
from sympy import Rational
from sympy.combinatorics import Permutation, PermutationGroup
from sympy.physics.quantum import pauli, represent
from scipy.special import comb
from sympy.core.singleton import S
from functools import reduce
import random
import operator
import importlib

# from lattice_symmetries import Expr
# import math
import igraph as ig
# import glob
# import json
import numpy as np
from numpy.testing import assert_equal

# import os
import scipy.sparse.linalg
import pytest
from pytest import raises, approx
import hypothesis
import hypothesis.strategies as st
import hypothesis.extra.numpy

has_quspin = importlib.util.find_spec("quspin") is not None

our_phases = (hypothesis.Phase.explicit, hypothesis.Phase.reuse, hypothesis.Phase.generate)

rng = np.random.default_rng(seed=123)


def test_trailing_zeros():
    from lattice_symmetries.basis import trailing_zeros
    assert trailing_zeros(0) == 0
    assert trailing_zeros(1) == 0
    assert trailing_zeros(2) == 1
    assert trailing_zeros(3) == 0
    assert trailing_zeros(4) == 2
    assert trailing_zeros(7) == 0
    assert trailing_zeros(8) == 3
    assert trailing_zeros(40) == 3
    assert trailing_zeros(128) == 7

def test_generate_offset_ranges():
    from lattice_symmetries._kernels import generate_offset_ranges

    representatives = np.array([0, 1, 2, 4, 8, 9, 16, 17], dtype=np.uint64)
    number_bits = 3
    shift = 61
    offsets, range_size = generate_offset_ranges(representatives, number_bits, shift)
    np.testing.assert_equal(offsets, [0, 0, 0, 0, 0, 0, 0, 0, 8])
    assert range_size == 8
    
    number_bits = 3
    shift = 2
    offsets, range_size = generate_offset_ranges(representatives, number_bits, shift)
    np.testing.assert_equal(offsets, [0, 3, 4, 5, 5, 5, 5, 5, 8])
    assert range_size == 3

    representatives = np.array([0, 1, 2, 3], dtype=np.uint64)
    number_bits = 0
    for shift in range(64):
        offsets, range_size = generate_offset_ranges(representatives, number_bits, shift)
        assert len(offsets) == 2
        assert offsets[0] == 0
        assert offsets[1] == len(representatives)
        assert range_size == len(representatives)

    representatives = np.array([42], dtype=np.uint64)
    number_bits = 2
    shift = 4
    offsets, range_size = generate_offset_ranges(representatives, number_bits, shift)
    np.testing.assert_equal(offsets, [0, 0, 0, 0, 1])
    assert range_size == 1

    # Test with numbers that have large gaps between them
    representatives = np.array([0, (1 << 32), (2 << 32), (3 << 32)], dtype=np.uint64)
    number_bits = 2
    shift = 32
    offsets, range_size = generate_offset_ranges(representatives, number_bits, shift)
    assert len(offsets) == 5  # 2^2 + 1
    assert offsets[0] == 0
    assert offsets[1] == 1
    assert offsets[2] == 2
    assert offsets[3] == 3
    assert offsets[4] == 4
    assert range_size == 1

def test_state_to_index_binary_search():
    from lattice_symmetries._kernels import state_to_index_kernel

    # Test with small numbers first
    representatives = np.array([0, 1, 2, 4, 8, 9, 16, 17], dtype=np.int64)
    kernel = state_to_index_kernel(representatives, prefix_bits=1)
    alpha = np.array([0, 1, 4, 5], dtype=np.int64)
    out = np.zeros(alpha.size, dtype=np.int32)
    kernel.callable(alpha, representatives, out)
    np.testing.assert_equal(out, [0, 1, 3, -1])

    # Test with larger numbers and more complex patterns
    representatives = np.array([0, 3, 7, 15, 31, 63, 127, 255, 511, 1023], dtype=np.int64)
    kernel = state_to_index_kernel(representatives, prefix_bits=4)
    alpha = np.array([31, 0, 255, 7, 1024, 15, 63, 512, 3, 127, 8, 511, 32, 16, 256, 64, 1023], dtype=np.int64)
    out = np.zeros(alpha.size, dtype=np.int32)
    kernel.callable(alpha, representatives, out)
    np.testing.assert_equal(out, [4, 0, 7, 2, -1, 3, 5, -1, 1, 6, -1, 8, -1, -1, -1, -1, 9])


def test_SpinBasis():
    basis = ls.SpinBasis(3)
    basis.build()
    np.testing.assert_equal(basis.index(basis.states), np.arange(2**3))
    assert basis.number_states == 2**3
    assert [basis.state_to_string(basis.states[i]) for i in range(basis.number_states)] == [
        "|000⟩",
        "|001⟩",
        "|010⟩",
        "|011⟩",
        "|100⟩",
        "|101⟩",
        "|110⟩",
        "|111⟩",
    ]

    basis = ls.SpinBasis(3, hamming_weight=2)  # We want the subspace with only 2 spins up
    basis.build()
    assert [basis.state_to_string(basis.states[i]) for i in range(basis.number_states)] == [
        "|011⟩",
        "|101⟩",
        "|110⟩",
    ]

    basis = ls.SpinBasis(4, hamming_weight=2, spin_inversion=-1)
    basis.build()
    assert [basis.state_to_string(basis.states[i]) for i in range(basis.number_states)] == [
        "|0011⟩",
        "|0101⟩",
        "|0110⟩",
    ]


def reverse_bits(x, n_bits):
    x = np.array(x)
    x_reversed = np.zeros_like(x)
    for _ in range(n_bits):
        x_reversed = (x_reversed << 1) | x & 1
        x >>= 1
    return x_reversed


def permute_bits(x, permutation):
    x = np.array(x)
    x_permuted = np.zeros_like(x)
    for i in range(permutation.size):
        x_permuted |= ((x >> permutation[i]) & 1) << i
    return x_permuted


def from_quspin_states(basis_states, basis):
    # QuSpin orders basis states dirrefently:
    #   - It stores the spin 0 in the most significant bit, and lattice-symmetries in the least significant bit
    basis_states = reverse_bits(basis_states, basis.number_bits)
    #   - It uses 1 to represent ↑, and 0 to represent ↓, but lattice symmetries does the inverse 😭
    basis_states = basis_states ^ ((1 << basis.number_bits) - 1)

    representatives = basis_states
    characters = np.ones(basis_states.size, dtype=np.complex128)
    if basis.spin_inversion is not None:
        mask = (1 << basis.number_bits) - 1
        inverted = representatives ^ mask
        characters[inverted < representatives] = basis.spin_inversion
        representatives = np.minimum(representatives, inverted)
    return representatives, characters


@pytest.mark.skipif(not has_quspin, reason="QuSpin not available")
@pytest.mark.parametrize("number_bits", list(range(1, 17)))
def test_SpinBasis_quspin_no_symmetries(number_bits):
    import quspin

    basis = ls.SpinBasis(number_spins=number_bits)
    basis.build()

    quspin_basis = quspin.basis.spin_basis_1d(number_bits)
    ref, _ = from_quspin_states(quspin_basis.states, basis)
    ref.sort()

    np.testing.assert_equal(basis.states, ref)


@st.composite
def random_number_bits_and_hamming_weight(draw, max_number_states=2**16):
    number_bits = draw(st.integers(min_value=1, max_value=64))
    predicate = lambda h: sympy.binomial(number_bits, h) < max_number_states
    hamming_weights = list(filter(predicate, range(number_bits + 1)))
    hamming_weight = draw(st.sampled_from(hamming_weights))
    return number_bits, hamming_weight


@st.composite
def random_basis_info_1d(draw, max_number_states=2**16):
    number_bits = draw(st.integers(min_value=1, max_value=64))

    predicate = lambda h: sympy.binomial(number_bits, h) / number_bits < max_number_states
    choices = (st.sampled_from(list(filter(predicate, range(number_bits + 1)))),)
    if 2**number_bits / number_bits < max_number_states:
        choices = choices + (st.just(None),)
    if number_bits % 2 == 0:
        choices = choices + (st.just(number_bits // 2),)
    hamming_weight = draw(st.one_of(*choices))

    if hamming_weight is None:
        choices = (None, -1, 1)
    else:
        choices = (None, -1, 1) if 2 * hamming_weight == number_bits else (None,)
    spin_inversion = draw(st.sampled_from(choices))

    if number_bits > 1:
        translation = draw(st.integers(min_value=1, max_value=number_bits - 1))
        p = Permutation(list(range(1, number_bits)) + [0])
        symmetries = [(p, Rational(translation, number_bits))]
    else:
        symmetries = []

    return ls.BasisInfo(number_bits, hamming_weight, spin_inversion, symmetries)


@pytest.mark.skipif(not has_quspin, reason="QuSpin not available")
@hypothesis.given(random_number_bits_and_hamming_weight(max_number_states=2**16))
@hypothesis.settings(max_examples=10, deadline=None, phases=our_phases)
def test_SpinBasis_quspin_U1(args):
    import quspin

    number_bits, hamming_weight = args
    basis = ls.SpinBasis(number_spins=number_bits, hamming_weight=hamming_weight)
    basis.build()

    quspin_basis = quspin.basis.spin_basis_1d(number_bits, Nup=number_bits - hamming_weight)
    ref, _ = from_quspin_states(quspin_basis.states, basis)
    ref.sort()

    np.testing.assert_equal(basis.states, ref)


@pytest.mark.skipif(not has_quspin, reason="QuSpin not available")
@hypothesis.given(st.integers(min_value=1, max_value=16), st.sampled_from([-1, 1]))
@hypothesis.settings(max_examples=10, deadline=None, phases=our_phases)
def test_SpinBasis_quspin_Z2(number_bits, spin_inversion):
    import quspin

    basis = ls.SpinBasis(number_spins=number_bits, spin_inversion=spin_inversion)
    basis.build()

    quspin_basis = quspin.basis.spin_basis_1d(number_bits, zblock=spin_inversion)
    ref, _ = from_quspin_states(quspin_basis.states, basis)
    ref.sort()

    np.testing.assert_equal(basis.states, ref)

    if number_bits % 2 == 0:
        basis = ls.SpinBasis(number_bits, number_bits // 2, spin_inversion=spin_inversion)
        basis.build()

        quspin_basis = quspin.basis.spin_basis_1d(
            number_bits, Nup=number_bits // 2, zblock=spin_inversion
        )
        ref, _ = from_quspin_states(quspin_basis.states, basis)
        ref.sort()

        np.testing.assert_equal(basis.states, ref)


@pytest.mark.skipif(not has_quspin, reason="QuSpin not available")
@hypothesis.example(ls.BasisInfo(4, 2))
@hypothesis.example(ls.BasisInfo(10, 0))
@hypothesis.example(ls.BasisInfo(10, 6))
@hypothesis.example(ls.BasisInfo(20, 17))
@hypothesis.given(random_basis_info_1d(max_number_states=2**10))
@hypothesis.settings(max_examples=5, deadline=None, phases=our_phases)
def test_SpinBasis_quspin_1d(info):
    from quspin.basis import spin_basis_1d

    basis = ls.SpinBasis(
        number_spins=info.number_bits,
        hamming_weight=info.hamming_weight,
        spin_inversion=info.spin_inversion,
        symmetries=info.symmetries,
    )
    logger.debug(info)
    tick = time.perf_counter()
    basis.build()
    tock = time.perf_counter()
    logger.debug(f"build() took {tock - tick}")

    Nup = info.number_bits - info.hamming_weight if info.hamming_weight is not None else None
    if len(info.symmetries) > 0:
        kblock = int(info.number_bits * info.symmetries[0][1])
        if kblock != 0:
            kblock = info.number_bits - kblock
    else:
        kblock = None
    zblock = info.spin_inversion

    tick = time.perf_counter()
    quspin_basis = spin_basis_1d(info.number_bits, Nup=Nup, kblock=kblock, zblock=zblock)
    tock = time.perf_counter()
    logger.debug(f"spin_basis_1d() took {tock - tick}")
    ref = quspin_basis.states ^ (2**info.number_bits - 1)

    np.testing.assert_equal(basis.states, ref)


def test_is_representative_examples():
    symmetries = [(Permutation([1, 2, 3, 0]), Rational(0, 4))]
    basis = ls.SpinBasis(4, symmetries=symmetries)
    kernel = lattice_symmetries._kernels.is_representative_kernel(basis.info, verbose=True)

    states = np.arange(2**basis.info.number_bits).astype(np.int64)
    out = np.zeros(states.size, dtype=np.uint16)
    kernel.callable(states, out)

    for a, n in zip(states, out):
        print(basis.state_to_string(a), n)

    basis.build()
    print(basis.states)

    symmetries = [(Permutation(list(range(1, 20)) + [0]), Rational(0, 20))]
    basis = ls.SpinBasis(20, hamming_weight=17, symmetries=symmetries)
    kernel = lattice_symmetries._kernels.is_representative_kernel(basis.info, verbose=True)
    states = np.array(
        [0b00011111111111111111, 0b00111111111111111110, 0b01111111111111111100]
    ).astype(np.int64)
    out = np.zeros(states.size, dtype=np.uint16)
    kernel.callable(states, out)
    for a, n in zip(states, out):
        print(basis.state_to_string(a), n)


@st.composite
def random_pauli_term(draw, number_sites, max_order=3):
    operators = [pauli.SigmaX, pauli.SigmaY, pauli.SigmaZ, pauli.SigmaPlus, pauli.SigmaMinus]
    p = draw(st.sampled_from(operators))
    c = draw(st.complex_numbers(min_magnitude=1e-3, max_magnitude=10, allow_subnormal=False))
    i = st.integers(min_value=0, max_value=number_sites - 1)
    indices = st.lists(i, min_size=0, max_size=max_order).map(sorted)
    return c * reduce(operator.mul, map(p, draw(indices)), S.One)


@st.composite
def random_pauli_expression(draw, max_number_sites=8, max_terms=10, max_order=3):
    number_sites = draw(st.integers(min_value=1, max_value=max_number_sites - 1))
    elements = st.lists(
        random_pauli_term(number_sites=number_sites, max_order=max_order),
        min_size=1,
        max_size=max_terms,
    )
    return reduce(operator.add, draw(elements))


@hypothesis.given(random_pauli_term(number_sites=4, max_order=3))
@hypothesis.settings(max_examples=10, deadline=None, phases=our_phases)
def test_foo(term):
    print(term)


def test_fixed_hamming_compilation():
    import halide as hl

    target = hl.get_jit_target_from_environment()
    target = target.with_feature(hl.TargetFeature.NoAsserts)
    target = target.with_feature(hl.TargetFeature.NoBoundsQuery)
    target = target.with_feature(hl.TargetFeature.AVX512_Zen4)
    _ = lattice_symmetries._kernels.xored_state_to_index_kernel(
        ls.BasisInfo(10, hamming_weight=3), target=target
    )


def test_fixed_hamming_state_to_index_examples():
    kernel = lattice_symmetries._kernels.xored_state_to_index_kernel(ls.BasisInfo(10, 2))

    # NOTE: uncomment the following line to dlclose the shared library early.
    # If all goes well, the program should crash (likely, with a segmentation fault)
    # lattice_symmetries._kernels.COMPILER.ffi.release(kernel.ffi_fun_ptr)

    states = np.arange(2**10, dtype=np.uint64)
    states = states[[x.bit_count() == 2 for x in states]]

    # Invoking using Halide::Callable
    np.random.seed(42)
    x = np.random.choice(states, size=10).astype(np.int64)
    mask = np.zeros(1, dtype=np.int64)
    out = np.zeros((1, len(x)), dtype=np.int64)
    kernel.callable(x, mask, states.view(np.int64), out)
    np.testing.assert_equal(states[out[0]], x)

    # Invoking using raw function pointers
    x_buf, x_buf_keep_alive = lattice_symmetries._kernels.create_halide_buffer_view(x)
    mask_buf, mask_buf_keep_alive = lattice_symmetries._kernels.create_halide_buffer_view(mask)
    states_buf, states_buf_keep_alive = lattice_symmetries._kernels.create_halide_buffer_view(
        states.view(np.int64)
    )
    out = np.zeros((1, len(x)), dtype=np.int64)
    out_buf, out_buf_keep_alive = lattice_symmetries._kernels.create_halide_buffer_view(out)
    kernel.ffi_fun_ptr(x_buf, mask_buf, states_buf, out_buf)
    np.testing.assert_equal(states[out[0]], x)


@hypothesis.given(
    st.builds(
        lambda a, b: (max(a, b), min(a, b)),
        st.integers(min_value=1, max_value=20),
        st.integers(min_value=1, max_value=20),
    ),
    st.integers(min_value=1, max_value=1000),
)
@hypothesis.settings(max_examples=10, deadline=None, phases=our_phases)
def test_fixed_hamming_state_to_index(args, batch_size):
    number_sites, hamming_weight = args
    info = ls.BasisInfo(number_sites, hamming_weight)
    kernel = lattice_symmetries._kernels.xored_state_to_index_kernel(info)
    states = np.arange(2**number_sites, dtype=np.uint64)
    states = states[[x.bit_count() == hamming_weight for x in states]]
    # Invoking using Halide::Callable
    np.random.seed(42)
    x = np.random.choice(states, size=batch_size).astype(np.int64)
    mask = np.full((7,), fill_value=0, dtype=np.int64)
    out = np.zeros((mask.size, x.size), dtype=np.int64)
    kernel.callable(x, mask, states.view(np.int64), out)
    np.testing.assert_equal(states[out], x.reshape(1, -1) ^ mask.reshape(-1, 1))
    # Invoking using raw function pointers
    out = np.zeros((mask.size, x.size), dtype=np.int64)
    x_buf, x_buf_keep_alive = lattice_symmetries._kernels.create_halide_buffer_view(x)
    mask_buf, mask_buf_keep_alive = lattice_symmetries._kernels.create_halide_buffer_view(mask)
    states_buf, states_buf_keep_alive = lattice_symmetries._kernels.create_halide_buffer_view(
        states.view(np.int64)
    )
    out_buf, out_buf_keep_alive = lattice_symmetries._kernels.create_halide_buffer_view(out)
    kernel.ffi_fun_ptr(x_buf, mask_buf, states_buf, out_buf)
    np.testing.assert_equal(states[out], x.reshape(1, -1) ^ mask.reshape(-1, 1))


def notest_is_representative_examples():
    info = ls.BasisInfo(number_bits=10)
    kernel = ls._kernels.is_representative_kernel(info)

    # Invoking using Halide::Callable
    x = np.random.choice(2**info.number_bits, size=10).astype(np.int64)
    out = np.zeros(len(x), dtype=np.uint16)
    kernel.callable(x, out)
    np.testing.assert_equal(out, np.ones(len(x), dtype=np.uint16))

    # Invoking using raw function pointers
    # x_buf, x_buf_keep_alive = lattice_symmetries._kernels.create_halide_buffer_view(x)
    # out = np.zeros(len(x), dtype=np.int64)
    # out_buf, out_buf_keep_alive = lattice_symmetries._kernels.create_halide_buffer_view(out)
    # kernel.ffi_fun_ptr(x_buf, out_buf)
    # np.testing.assert_equal(states[np.asarray(out)], x)


def test_enumerate_basis_states():
    states, norms = ls.enumerate_basis_states(ls.BasisInfo(number_bits=5))
    assert_equal(states, np.arange(2**5, dtype=np.uint64))
    assert_equal(norms, np.ones(2**5, dtype=np.uint16))

    states, norms = ls.enumerate_basis_states(ls.BasisInfo(number_bits=4, hamming_weight=2))
    assert_equal(states, np.array([0b0011, 0b0101, 0b0110, 0b1001, 0b1010, 0b1100]))
    assert_equal(norms, np.ones(6, dtype=np.uint16))


@hypothesis.given(
    st.builds(
        lambda a, b: (max(a, b), min(a, b)),
        st.integers(min_value=1, max_value=20),
        st.integers(min_value=1, max_value=20),
    ),
)
@hypothesis.settings(max_examples=10, deadline=None, phases=our_phases)
def test_fixed_hamming_state_to_index_scalar(args):
    number_sites, hamming_weight = args
    batch_size = 16
    states = np.arange(2**number_sites, dtype=np.uint64)
    states = states[[x.bit_count() == hamming_weight for x in states]]

    np.random.seed(42)
    xs = np.random.choice(states, size=batch_size)

    for i, x in enumerate(xs):
        assert states[ls.fixed_hamming_state_to_index(int(x))] == x


@hypothesis.given(
    st.builds(
        lambda a, b: (max(a, b), min(a, b)),
        st.integers(min_value=1, max_value=20),
        st.integers(min_value=1, max_value=20),
    ),
)
@hypothesis.settings(max_examples=10, deadline=None, phases=our_phases)
def test_fixed_hamming_index_to_state_scalar(args):
    number_sites, hamming_weight = args
    batch_size = 16
    states = np.arange(2**number_sites, dtype=np.uint64)
    states = states[[x.bit_count() == hamming_weight for x in states]]

    np.random.seed(42)
    indices = np.random.choice(len(states), size=batch_size)

    for i in indices:
        x = ls.fixed_hamming_index_to_state(
            i, number_sites=number_sites, hamming_weight=hamming_weight
        )
        assert x == states[i]


def test_permutation_to_benes_network_examples():
    benes = ls.permutation_to_benes_network(Permutation([0, 1, 2, 3]))
    assert benes(0b0100) == 0b0100
    benes = ls.permutation_to_benes_network(Permutation([1, 2, 0, 3]))
    assert benes(0b0100) == 0b0010
    assert benes(0b0101) == 0b0110
    benes = ls.permutation_to_benes_network(Permutation([1, 2, 3, 0, 5, 6, 7, 4, 9, 10, 11, 8]))
    assert benes(0b100111010010) == 0b110011100001
    benes = ls.permutation_to_benes_network(Permutation([0]))
    assert benes(0b0) == 0b0
    assert benes(0b1) == 0b1
    benes = ls.permutation_to_benes_network(Permutation([]))
    assert benes(0b0) == 0b0
    assert benes(0b10100110) == 0b10100110


@st.composite
def permutations(draw, min_size=1, max_size=1000):
    size = draw(st.integers(min_value=min_size, max_value=max_size))
    arr = draw(st.permutations(list(range(size))))
    return Permutation(arr)


def reference_permute_bits(p: Permutation, bits: int):
    assert bits < 2**p.size
    s = ("{:0" + str(p.size) + "b}").format(bits)[::-1]
    s = "".join(s[int(i)] for i in p.array_form)[::-1]
    return int(s, base=2)


@hypothesis.given(permutations())
@hypothesis.settings(max_examples=10, deadline=None, phases=our_phases)
def test_permutation_to_benes_network(permutation):
    random.seed(42)
    benes = ls.permutation_to_benes_network(permutation)
    for i in range(10):
        bits = random.randint(0, 2**permutation.size - 1)
        assert benes(bits) == reference_permute_bits(permutation, bits)


normal_complex = st.complex_numbers(
    max_magnitude=10, allow_infinity=False, allow_nan=False, allow_subnormal=False
)
two_same_length_vectors = st.integers(min_value=0, max_value=1000).flatmap(
    lambda n: st.lists(
        hypothesis.extra.numpy.arrays(np.complex128, n, elements=normal_complex),
        min_size=2,
        max_size=2,
    )
)


@hypothesis.given(normal_complex, two_same_length_vectors)
@hypothesis.settings(deadline=None, phases=our_phases)
def test_axpy(alpha, arrays):
    x, y = arrays
    x[np.isnan(x) | np.isinf(x)] = 0
    y[np.isnan(y) | np.isinf(y)] = 0

    ref = alpha * x + y
    ls.axpy(alpha, x, y)  # NOTE: modifies y inplace
    np.testing.assert_allclose(y, ref)


# def test_diag_matrix_kernel():
#     from lattice_symmetries.expression import (
#         pauli_expression_to_nonbranching_terms,
#         lower_nonbranching_terms,
#         pauli_expression_to_matrix,
#     )
#
#     info = ls.BasisInfo(2)
#     e = sum(rng.uniform(0, 1) * pauli.SigmaZ(i) for i in range(info.number_bits))
#     print(e)
#
#     op = ls._kernels.LoweredOperator(pauli_expression_to_nonbranching_terms(e))
#
#     states = np.arange(2**info.number_bits)
#     out_re = np.zeros(states.size, dtype=np.float64)
#     out_im = np.zeros(states.size, dtype=np.float64)
#     op.diag_kernel.callable(states, op.terms.v_re_diag, op.terms.v_im_diag, out_re, out_im)
#     ref = np.diag(np.asarray(pauli_expression_to_matrix(e), dtype=np.complex128))
#     np.testing.assert_allclose(out_re + 1j * out_im, ref)
#
#     out = np.zeros(states.size, dtype=np.complex128)
#     x = rng.uniform(0, 1, size=states.size) + rng.uniform(0, 1, size=states.size) * 1j
#     op.apply_diag(states, x.reshape(-1, 1), out.reshape(-1, 1))
#     np.testing.assert_allclose(out, ref * x)

# def random_operator(max_order: int=3):
#     for order
#     pass


@hypothesis.given(random_pauli_expression(max_number_sites=5))
@hypothesis.example(sympy.Float(0.001))
@hypothesis.example(-9.999 - 2.2250738585072e-309 * sympy.I - 10.0 * sympy.I * pauli.SigmaX(1))
@hypothesis.example(pauli.SigmaPlus(0))
@hypothesis.example(
    (0.5 + 5.96046447753906e-8 * sympy.I) * pauli.SigmaPlus(1)
    + (-1.40129846432482e-45 - 5.80934524614774 * sympy.I) * pauli.SigmaZ(2)
)
@hypothesis.example(
    (-0.32052911031739 + 9.99486173438328 * sympy.I) * pauli.SigmaPlus(0) * pauli.SigmaPlus(1)
    + 10.0 * pauli.SigmaX(0)
)
@hypothesis.settings(max_examples=10, deadline=None, phases=our_phases)
def test_matrix_apply(e):
    from lattice_symmetries.expression import (
        pauli_expression_to_nonbranching_terms,
        pauli_expression_to_matrix,
        _collect_indices,
    )

    number_bits = 1 + max(_collect_indices(e), default=0)
    info = ls.BasisInfo(number_bits)
    terms = pauli_expression_to_nonbranching_terms(e)
    states = np.arange(2**info.number_bits).astype(np.int64)
    op = ls._kernels.LoweredOperator(info=info, terms=terms, symm=None, state_to_index_info=None)
    ref = np.asarray(pauli_expression_to_matrix(e), dtype=np.complex128)

    # x = np.zeros(states.size, dtype=np.complex128)
    # x[0] = 1
    x = rng.uniform(0, 1, size=states.size).astype(
        np.float32
    )  # + rng.uniform(0, 1, size=states.size) * 1j
    out = np.zeros(states.size, dtype=np.float32)

    op.apply(states, x, out)
    np.testing.assert_allclose(out, (ref @ x).real, rtol=1e-6, atol=1e-7)

@pytest.mark.skipif(not has_quspin, reason="QuSpin not available")
def test_matrix_apply_quspin_U1():
    np.random.seed(42)
    import quspin
    from lattice_symmetries.expression import (
        pauli_expression_to_nonbranching_terms,
        pauli_expression_to_matrix,
        _collect_indices,
    )

    number_bits = 10
    quspin_basis = quspin.basis.spin_basis_1d(number_bits, Nup=number_bits // 2)
    quspin_J_nn = [[1, i, (i + 1) % number_bits] for i in range(number_bits)]
    quspin_static = [["xx", quspin_J_nn], ["yy", quspin_J_nn], ["zz", quspin_J_nn]]
    quspin_dynamic = []
    hamiltonian = quspin.operators.hamiltonian(quspin_static, quspin_dynamic, dtype=np.float64, basis=quspin_basis)

    number_states = quspin_basis.states.size
    x = np.random.rand(number_states)
    # x = np.zeros(number_states)
    # x[1] = 1
    ref = hamiltonian.dot(x)

    mk = lambda n: sum(
        pauli.SigmaX(i) * pauli.SigmaX((i + 1) % n)
        + pauli.SigmaY(i) * pauli.SigmaY((i + 1) % n)
        + pauli.SigmaZ(i) * pauli.SigmaZ((i + 1) % n)
        for i in range(n)
    )

    e = mk(number_bits)

    basis = ls.SpinBasis(number_spins=number_bits, hamming_weight=number_bits // 2)
    basis.build()

    terms = pauli_expression_to_nonbranching_terms(e)
    op = ls._kernels.LoweredOperator(info=basis.info, terms=terms, symm=None, state_to_index_info=None)

    out = np.zeros(number_states, dtype=np.float32)
    op.apply(basis.states, x.astype(np.float32), out)

    np.testing.assert_allclose(out, ref, rtol=1e-6, atol=1e-7)

def test_matrix_apply_example():
    mk = lambda n: sum(
        pauli.SigmaX(i) * pauli.SigmaX((i + 1) % n)
        + pauli.SigmaY(i) * pauli.SigmaY((i + 1) % n)
        + pauli.SigmaZ(i) * pauli.SigmaZ((i + 1) % n)
        for i in range(n)
    )

    from lattice_symmetries.expression import (
        pauli_expression_to_nonbranching_terms,
        pauli_expression_to_matrix,
        _collect_indices,
    )

    e = mk(5)
    number_bits = 1 + max(_collect_indices(e), default=0)
    info = ls.BasisInfo(number_bits)
    terms = pauli_expression_to_nonbranching_terms(e)
    states = np.arange(2**info.number_bits).astype(np.int64)
    op = ls._kernels.LoweredOperator(info=info, terms=terms, symm=None, state_to_index_info=None)
    ref = np.asarray(pauli_expression_to_matrix(e), dtype=np.complex128)

    x = rng.uniform(0, 1, size=states.size).astype(
        np.float32
    )  # + rng.uniform(0, 1, size=states.size) * 1j
    out = np.zeros(states.size, dtype=np.float32)

    # op.off_diag_kernel.callable(states, op.terms.v_re_2d, op.terms.v_im_2d, x, out)
    op.apply(states, x, out)

    # op.apply(states, x.reshape(-1, 1), out.reshape(-1, 1))

    np.testing.assert_allclose(out, (ref @ x).real, rtol=1e-6, atol=1e-7)

    # info = ls.BasisInfo(2)
    # e = mk(info.number_bits)
    # terms = pauli_expression_to_nonbranching_terms(e)

    # states = np.arange(2**info.number_bits)
    # # for alpha in states:
    # #     print([t.act_on_ket(alpha) for t in terms])

    # op = ls._kernels.LoweredOperator(terms)

    # ref = np.asarray(pauli_expression_to_matrix(e), dtype=np.complex128)
    # # ref -= np.diag(np.diag(ref))

    # out = np.zeros(states.size, dtype=np.complex128)
    # # x = np.zeros(states.size, dtype=np.complex128)
    # # x[0] = 1
    # x = rng.uniform(0, 1, size=states.size) + rng.uniform(0, 1, size=states.size) * 1j
    # op.apply(states, x.reshape(-1, 1), out.reshape(-1, 1))
    # # print(ref)

    # np.testing.assert_allclose(out, ref @ x)

    # target = hl.get_jit_target_from_environment()
    # target = target.with_feature(hl.TargetFeature.NoAsserts)
    # target = target.with_feature(hl.TargetFeature.NoBoundsQuery)
    # target = target.with_feature(hl.TargetFeature.AVX512_Zen4)

def test_matrix_apply_example_symmetries():
    mk = lambda n: sum(
        pauli.SigmaX(i) * pauli.SigmaX((i + 1) % n)
        + pauli.SigmaY(i) * pauli.SigmaY((i + 1) % n)
        + pauli.SigmaZ(i) * pauli.SigmaZ((i + 1) % n)
        for i in range(n)
    )

    from lattice_symmetries.expression import (
        pauli_expression_to_nonbranching_terms,
        pauli_expression_to_matrix,
        _collect_indices,
    )

    number_bits = 4
    e = mk(number_bits)
    symmetries = [(Permutation(list(range(1, number_bits)) + [0]), Rational(0, 1))]
    basis = ls.SpinBasis(number_spins=number_bits, symmetries=symmetries)
    basis.build()

    terms = pauli_expression_to_nonbranching_terms(e)
    symm = ls._kernels.LoweredSymmetries(basis.info.symmetries)


    total_bits = int(basis.states.max()).bit_length()
    prefix_bits = 1
    shift = max(0, total_bits - prefix_bits)
    offsets, range_size = ls._kernels.generate_offset_ranges(basis.states, prefix_bits, shift)
    state_to_index_info = ls._kernels.StateToIndexInfo(
        offsets=offsets,
        shift=shift,
        prefix_bits=prefix_bits,
        range_size=range_size,
    )

    op = ls._kernels.LoweredOperator(info=basis.info, terms=terms, symm=symm, state_to_index_info=state_to_index_info)
    ref = np.asarray(pauli_expression_to_matrix(e), dtype=np.complex128)

    x = rng.uniform(0, 1, size=basis.states.size).astype(
        np.float32
    )  # + rng.uniform(0, 1, size=basis.states.size) * 1j
    out = np.zeros(basis.states.size, dtype=np.float32)

    # op.off_diag_kernel.callable(states, op.terms.v_re_2d, op.terms.v_im_2d, x, out)
    op.apply(basis.states, x, out)

    # op.apply(states, x.reshape(-1, 1), out.reshape(-1, 1))
    print(out)


def test_matrix_apply_example_2():
    from lattice_symmetries.expression import (
        pauli_expression_to_nonbranching_terms,
        pauli_expression_to_matrix,
        _collect_indices,
    )

    e = - 0.22463139456015124 * pauli.SigmaZ(1) * pauli.SigmaPlus(0) * pauli.SigmaMinus(1) * pauli.SigmaZ(0) * pauli.SigmaZ(1) * pauli.SigmaZ(1) \
        + 0.4892602908941869 * pauli.SigmaZ(1) \
        + 0.437521418673926 * pauli.SigmaPlus(0) * pauli.SigmaMinus(0) * pauli.SigmaPlus(1) * pauli.SigmaMinus(0) * pauli.SigmaPlus(0) * pauli.SigmaMinus(1) \
        - 0.19991708088098303 * pauli.SigmaPlus(1) * pauli.SigmaMinus(1) * pauli.SigmaZ(0) * pauli.SigmaPlus(1) * pauli.SigmaMinus(1) * pauli.SigmaZ(1) \
        + 0.481139511164027 * pauli.SigmaZ(1) * pauli.SigmaZ(1) * pauli.SigmaZ(0) * pauli.SigmaZ(0) * pauli.SigmaZ(0) * pauli.SigmaZ(1) \
        + 0.377901909342394 * pauli.SigmaPlus(0) * pauli.SigmaMinus(0) * pauli.SigmaZ(0) \
        - 0.13375235591083356 * pauli.SigmaZ(1) * pauli.SigmaZ(1) * pauli.SigmaPlus(0) * pauli.SigmaMinus(1) * pauli.SigmaZ(0) \
        + 0.26327045984176356 * pauli.SigmaZ(1) * pauli.SigmaZ(0) * pauli.SigmaZ(0) \
        - 0.15872758255707742 * pauli.SigmaPlus(0) * pauli.SigmaMinus(0) * pauli.SigmaPlus(1) * pauli.SigmaMinus(0) * pauli.SigmaPlus(1) * pauli.SigmaMinus(1) \
        + 0.23851168527648248 * pauli.SigmaZ(1) * pauli.SigmaZ(0) * pauli.SigmaZ(0)

    number_bits = 2
    hamming_weight = 0
    symmetries = [(Permutation([0, 1]), Rational(0, 1))]
    basis = ls.SpinBasis(number_spins=number_bits, hamming_weight=hamming_weight, symmetries=symmetries)
    basis.build()

    terms = pauli_expression_to_nonbranching_terms(e)
    symm = ls._kernels.LoweredSymmetries(basis.info.symmetries)
    total_bits = 2
    prefix_bits = 0
    shift = max(0, total_bits - prefix_bits)
    offsets, range_size = ls._kernels.generate_offset_ranges(basis.states, prefix_bits, shift)
    state_to_index_info = ls._kernels.StateToIndexInfo(
        offsets=offsets,
        shift=shift,
        prefix_bits=prefix_bits,
        range_size=range_size,
    )
    op = ls._kernels.LoweredOperator(info=basis.info, terms=terms, symm=symm, state_to_index_info=state_to_index_info)
    x = np.asarray([0.920056717860089], dtype=np.float32)
    out = np.zeros(basis.states.size, dtype=np.float32)
    op.apply(basis.states, x, out)

    np.testing.assert_allclose(out, [1.5182470275151454], rtol=1e-6, atol=1e-7)


def test_simplify_pauli_expression():
    from lattice_symmetries.expression import simplify_pauli_expression
    e1 = pauli.SigmaX(0) * pauli.SigmaX(1) + pauli.SigmaY(0) * pauli.SigmaY(1)
    e2 = 2 * (pauli.SigmaPlus(0) * pauli.SigmaMinus(1) + pauli.SigmaPlus(1) * pauli.SigmaMinus(0))
    assert simplify_pauli_expression(e1) == simplify_pauli_expression(e2)

def test_matrix_apply_example_3():
    from lattice_symmetries.expression import (
        pauli_expression_to_nonbranching_terms,
        pauli_expression_to_matrix,
        _collect_indices,
    )

    number_bits = 3
    hamming_weight = None
    symmetries = [(Permutation([1, 2, 0]), Rational(0, 1))]
    e = sum(pauli.SigmaPlus(i) * pauli.SigmaMinus((i + 1) % number_bits) + pauli.SigmaPlus((i + 1) % number_bits) * pauli.SigmaMinus(i) for i in range(number_bits))

    basis = ls.SpinBasis(number_spins=number_bits, hamming_weight=hamming_weight, symmetries=symmetries)
    basis.build()
    np.testing.assert_array_equal(basis.states, np.array([0, 1, 3, 7], dtype=np.uint64))
    np.testing.assert_array_equal(basis.norms, np.array([3, 1, 1, 3], dtype=np.uint16))

    terms = pauli_expression_to_nonbranching_terms(e)
    print(terms)
    symm = ls._kernels.LoweredSymmetries(basis.info.symmetries)
    total_bits = number_bits
    prefix_bits = 10
    shift = max(0, total_bits - prefix_bits)
    offsets, range_size = ls._kernels.generate_offset_ranges(basis.states, prefix_bits, shift)
    state_to_index_info = ls._kernels.StateToIndexInfo(
        offsets=offsets,
        shift=shift,
        prefix_bits=prefix_bits,
        range_size=range_size,
    )
    op = ls._kernels.LoweredOperator(info=basis.info, terms=terms, symm=symm, state_to_index_info=state_to_index_info)

    out = np.zeros(basis.states.size, dtype=np.float32)
    x = np.asarray([1, 0, 0, 0], dtype=np.float32)
    op.apply(basis.states, basis.norms, x, out)
    np.testing.assert_allclose(out, [0, 0, 0, 0], rtol=1e-6, atol=1e-7)

    out = np.zeros(basis.states.size, dtype=np.float32)
    x = np.asarray([0, 1, 0, 0], dtype=np.float32)
    op.apply(basis.states, basis.norms, x, out)
    np.testing.assert_allclose(out, [0, 2, 0, 0], rtol=1e-6, atol=1e-7)

    out = np.zeros(basis.states.size, dtype=np.float32)
    x = np.asarray([0, 0, 1, 0], dtype=np.float32)
    op.apply(basis.states, basis.norms, x, out)
    np.testing.assert_allclose(out, [0, 0, 2, 0], rtol=1e-6, atol=1e-7)

    out = np.zeros(basis.states.size, dtype=np.float32)
    x = np.asarray([0, 0, 0, 1], dtype=np.float32)
    op.apply(basis.states, basis.norms, x, out)
    np.testing.assert_allclose(out, [0, 0, 0, 0], rtol=1e-6, atol=1e-7)

def test_matrix_apply_example_4():
    from lattice_symmetries.expression import (
        pauli_expression_to_nonbranching_terms,
        pauli_expression_to_matrix,
        _collect_indices,
    )

    number_bits = 6
    hamming_weight = 3
    symmetries = [(Permutation([1, 2, 3, 4, 5, 0]), Rational(1, 2))]
    e = sum(pauli.SigmaPlus(i) * pauli.SigmaMinus((i + 1) % number_bits) + pauli.SigmaPlus((i + 1) % number_bits) * pauli.SigmaMinus(i) for i in range(number_bits))

    basis = ls.SpinBasis(number_spins=number_bits, hamming_weight=hamming_weight, symmetries=symmetries)
    basis.build()
    np.testing.assert_array_equal(basis.states, np.array([7, 11, 13, 21], dtype=np.uint64))
    np.testing.assert_array_equal(basis.norms, np.array([1, 1, 1, 3], dtype=np.uint16))

    op = ls.Operator(ls.Expr(e), basis)
    out = np.zeros(basis.states.size, dtype=np.float32)
    x = np.asarray([0.7148606065054317, 0.5972164020369808, 0.7639741193866884, 0.01096456606490337], dtype=np.float32)
    op.apply_to_state_vector(x, out=out)
    np.testing.assert_allclose(out, [-0.16675771734970757, 2.2238176597714503, 0.49856338307588804, 0.2888328390039035], rtol=1e-6, atol=1e-7)

    # terms = pauli_expression_to_nonbranching_terms(e)
    # print(terms)
    # symm = ls._kernels.LoweredSymmetries(basis.info.symmetries)
    # total_bits = number_bits
    # prefix_bits = 0
    # shift = max(0, total_bits - prefix_bits)
    # offsets, range_size = ls._kernels.generate_offset_ranges(basis.states, prefix_bits, shift)
    # state_to_index_info = ls._kernels.StateToIndexInfo(
    #     offsets=offsets,
    #     shift=shift,
    #     prefix_bits=prefix_bits,
    #     range_size=range_size,
    # )
    # op = ls._kernels.LoweredOperator(info=basis.info, terms=terms, symm=symm, state_to_index_info=state_to_index_info)

    # out = np.zeros(basis.states.size, dtype=np.float32)
    # x = np.asarray([0.7148606065054317, 0.5972164020369808, 0.7639741193866884, 0.01096456606490337], dtype=np.float32)
    # op.apply(basis.states, basis.norms, x, out)
    # np.testing.assert_allclose(out, [-0.16675771734970757, 2.2238176597714503, 0.49856338307588804, 0.2888328390039035], rtol=1e-6, atol=1e-7)

def test_readme():
    from lattice_symmetries.expression import Expr

    assert Expr("σˣ₀") == Expr("\\sigma^x_0")
    assert np.array_equal(Expr("σˣ₀").to_dense(), np.array([[0, 1], [1, 0]]))
    assert Expr("Sˣ₀") == 0.5 * Expr("σˣ₀")
    assert Expr("σʸ₀") == Expr("\\sigma^y_0")
    assert np.array_equal(Expr("σʸ₀").to_dense(), np.array([[0, -1j], [1j, 0]]))
    assert Expr("Sʸ₀") == 0.5 * Expr("σʸ₀")
    assert Expr("σᶻ₀") == Expr("\\sigma^z_0")
    assert np.array_equal(Expr("σᶻ₀").to_dense(), np.array([[1, 0], [0, -1]]))
    assert Expr("Sᶻ₀") == 0.5 * Expr("σᶻ₀")
    assert np.array_equal(Expr("I", particle="spin-1/2").to_dense(), np.array([[1, 0], [0, 1]]))

    assert Expr("Sˣ₀ Sˣ₁ + Sʸ₀ Sʸ₁ + Sᶻ₀ Sᶻ₁") == Expr("Sx0 Sx1 + Sy0 Sy1 + Sz0 Sz1")
    # fmt: off
    assert Expr("Sˣ₀ Sˣ₁ + Sʸ₀ Sʸ₁ + Sᶻ₀ Sᶻ₁") == \
        Expr("Sˣ₀") * Expr("Sˣ₁") + Expr("Sʸ₀") * Expr("Sʸ₁") + Expr("Sᶻ₀") * Expr("Sᶻ₁")
    # fmt: on
    assert Expr("Sˣ₀ Sˣ₁ + Sʸ₀ Sʸ₁ + Sᶻ₀ Sᶻ₁") == Expr("0.5 (σ⁺₀ σ⁻₁ + σ⁺₁ σ⁻₀) + 0.25 σᶻ₀ σᶻ₁")

    # assert str(Expr("Sˣ₀ Sˣ₁ + Sʸ₀ Sʸ₁ + Sᶻ₀ Sᶻ₁")) == "0.25 σᶻ₀ σᶻ₁ + 0.5 σ⁺₀ σ⁻₁ + 0.5 σ⁻₀ σ⁺₁"
    # assert str(Expr("0.5 (σˣ₁ + 1im σʸ₁) - σ⁺₁")) == "0.0 I"
    # assert str(Expr("σ⁺₁ σ⁺₁")) == "0.0 I"


def test_expr_construction():
    _ = ls.Expr("2 I", particle="spin-1/2")
    _ = ls.Expr("2 I + S+3")
    _ = ls.Expr(expression="2 I + S+3", sites=[[2], [4]])
    _ = ls.Expr("5 σ⁺₀ σ⁻₁ + (8 + 3im) σ⁻₁")
    # _ = ls.Expr("-2 (c†₀ c₁ + c†₁ c₀)")

    # with raises(ValueError, match=r".*particle type.*"):
    #     ls.Expr(expression="2 I")
    # with raises(ValueError, match=r".*cannot replace.*"):
    #     ls.Expr(expression="S+3", sites=[[2], [3, 4]])
    # Ilya's expressions
    _ = ls.Expr("2 (σ⁺₀ σ⁻₁ + σ⁺₁ σ⁻₀) + σᶻ₀ σᶻ₁")
    _ = ls.Expr(
        "0.5 (σˣ₀ + σᶻ₀)(σˣ₁ + σᶻ₁)(σ⁺₀ σ⁻₁ + σ⁺₁ σ⁻₀)(σˣ₁ + σᶻ₁)(σˣ₀ + σᶻ₀) + 0.25 (σˣ₁ + σᶻ₁)(σˣ₀ + σᶻ₀)σᶻ₀ σᶻ₁(σˣ₀ + σᶻ₀)(σˣ₁ + σᶻ₁)"
    )

    def heisenberg_expr_rot(phi):
        C = f"({np.sin(phi)} σˣ₀ + {np.cos(phi)} σᶻ₀)({np.sin(phi)} σˣ₁ + {np.cos(phi)} σᶻ₁)"
        return f"2 {C}(σ⁺₀ σ⁻₁ + σ⁺₁ σ⁻₀) {C} + {C} σᶻ₀ σᶻ₁{C}"

    _ = ls.Expr(heisenberg_expr_rot(math.pi / 4))


def test_expr_replace_indices():
    a = ls.Expr("2 I + S+3")
    b = ls.Expr("2 I + S+57")
    assert a.replace_indices({3: 57}) == b


def test_expr_arithmetic():
    a = ls.Expr("2 I + S+3")
    b = ls.Expr("2 (σ⁺₀ σ⁻₁ + σ⁺₁ σ⁻₀) + σᶻ₀ σᶻ₁")
    assert a + b == b + a
    assert a + a == 2 * a
    assert a - b == -(b - a)
    assert a.adjoint() == ls.Expr("2 I + S-3")
    assert b.adjoint() == b


def test_expr_properties():
    a = ls.Expr("2 I + S+3")
    b = ls.Expr("2 (σ⁺₀ σ⁻₁ + σ⁺₁ σ⁻₀) + σᶻ₀ σᶻ₁")
    # assert a.is_real
    assert not a.is_hermitian
    assert not a.is_identity
    # assert b.is_real
    assert b.is_hermitian
    assert not b.is_identity
    # assert a.number_sites == 4
    # assert b.number_sites == 2

def test_sz_conserved_1d_heisenberg():
    number_spins = 8
    hamming_weight = number_spins // 2
    b = ls.SpinBasis(number_spins=number_spins, hamming_weight=hamming_weight)
    e = ls.Expr("2 (σ⁺₀ σ⁻₁ + σ⁺₁ σ⁻₀) + σᶻ₀ σᶻ₁").on(ig.Graph.Ring(n=number_spins, circular=True))
    h = ls.Operator(e, b)

    b.build()
    evals, evecs = scipy.sparse.linalg.eigsh(h, k=5)
    np.testing.assert_allclose(evals, [-14.60437363574869, -12.513676255378364, -10.798512593101313, -9.834954035579337, -9.83495403557929], rtol=1e-6, atol=1e-7)
#
# def test_basis_state_to_index():
#     # fmt: off
#     g = [
#         (ls.Permutation([3, 5, 6, 7, 9, 10, 11, 12, 13, 8, 14, 15, 16, 17, 18, 19, 1, 2, 0, 4]), ls.Rational(0)),
#         (ls.Permutation([11, 8, 14, 15, 12, 13, 18, 19, 1, 16, 17, 0, 4, 5, 2, 3, 9, 10, 6, 7]), ls.Rational(0))
#     ]
#     # fmt: on
#     b = ls.SpinBasis(number_spins=20, hamming_weight=2, symmetries=g)
#     b.build()
#     assert np.array_equal(b.index(b.states), np.arange(b.number_states))
#
#     bitstrings = np.arange(1 << 20, dtype=np.uint64)
#     indices = np.zeros_like(bitstrings, dtype=np.int64)
#     indices[:] = -1
#     indices[b.states] = np.arange(b.number_states)
#
#     for k in [1, 2, 3, 4, 8, 16, 128, 155, 1099]:
#         for _ in range(10):
#             batch = rng.choice(bitstrings, size=k)
#             assert np.array_equal(b.index(bitstrings[batch]), indices[batch])
#
#
# def test_basis_construction():
#     _ = ls.SpinBasis(4)
#     _ = ls.SpinBasis(4, 2)
#     _ = ls.SpinBasis(number_spins=4, spin_inversion=-1)
#
#     b = ls.SpinBasis(number_spins=64, hamming_weight=1)
#     b.build()
#     assert np.array_equal(b.states, 1 << np.arange(64, dtype=np.uint64))
#     assert (b.is_representative(b.states) > 0).all()
#
#     with raises(ValueError, match=r".*invalid spin inversion.*"):
#         ls.SpinBasis(number_spins=4, spin_inversion=5)
#     with raises(ValueError, match=r".*invalid Hamming weight.*"):
#         ls.SpinBasis(number_spins=2, hamming_weight=3)
#
#     _ = ls.Expr("2 I + S+3")
#
#
# def test_csr_matvec():
#     for dtype in [np.dtype("float64"), np.dtype("complex128")]:
#         for _ in range(10):
#             [n, m] = rng.integers(low=1, high=1000, size=2).tolist()
#             matrix = scipy.sparse.random(
#                 n, m, density=0.25, format="csr", dtype=dtype, random_state=rng
#             )
#             x = rng.random(matrix.shape[1])
#             if np.iscomplexobj(matrix.data):
#                 x = x + 1j * rng.random(matrix.shape[1])
#             assert np.allclose(ls._csr_matvec(matrix, x), matrix @ x)
#
#
# def test_build_matrix():
#     basis = ls.SpinBasis(2)
#     basis.build()
#     matrix = ls.Operator(ls.Expr("2 (σ⁺₀ σ⁻₁ + σ⁺₁ σ⁻₀) + σᶻ₀ σᶻ₁"), basis)
#     matrix.build_matrix()
#
#     h = np.array([[1, 0, 0, 0], [0, -1, 2, 0], [0, 2, -1, 0], [0, 0, 0, 1]])
#     x = np.random.rand(basis.number_states) + 1j * np.random.rand(basis.number_states)
#     assert np.allclose(matrix @ x, h @ x)
#
#
# def test_number_off_diag():
#     b = ls.SpinBasis(number_spins=10, hamming_weight=5)
#     b.build()
#     e = ls.Expr("2 (σ⁺₀ σ⁻₁ + σ⁺₁ σ⁻₀) + σᶻ₀ σᶻ₁", sites=ig.Graph.Ring(n=10, circular=True))
#     h = ls.Operator(e, b)
#
#     h.build_matrix()
#     real_max_number_off_diag = int(np.max(np.diff(h._off_diag_csr.indptr)))
#     # dynamite estimates the number of nonzero elements per row as 10
#     # we check that lattice-symmetries computes a similar estimate
#     assert real_max_number_off_diag <= h.estimate_max_number_off_diag() <= 10
#     assert real_max_number_off_diag <= h._payload.max_number_off_diag <= 10
#
#
# def test_off_diag_partial():
#     b = ls.SpinBasis(number_spins=6, hamming_weight=3)
#     b.build()
#     e = ls.Expr(
#         "2 (σ⁺₀ σ⁻₁ + σ⁺₁ σ⁻₀) + σᶻ₀ σᶻ₁",
#         sites=ig.Graph.Lattice(dim=[2, 3], circular=True),
#     )
#
#     for b in e.hilbert_space_sectors():
#         if b.hamming_weight != e.number_sites // 2:
#             continue
#         b.build()
#         h = ls.Operator(e, b)
#
#         x = rng.random(size=b.number_states)
#         y = h @ x
#
#         assert np.allclose(h.off_diag_to_csr() @ x + h.diag_to_array() * x, y)
#
#         for k in [1, 2, 10, 100]:
#             indices = rng.integers(b.number_states, size=k)
#             states = b.states[indices]
#             m = h.off_diag_to_triple(states=states, convert_to_index=True)
#             m = scipy.sparse.csr_matrix(m, shape=(states.size, b.number_states))
#             assert np.allclose(m @ x + h.diag_to_array()[indices] * x[indices], y[indices])
#
#
# def test_off_diag_ilya():
#     b = ls.SpinBasis(number_spins=6, hamming_weight=3)
#     e = ls.Expr(
#         "2 (σ⁺₀ σ⁻₁ + σ⁺₁ σ⁻₀) + σᶻ₀ σᶻ₁",
#         sites=ig.Graph.Lattice(dim=[2, 3], circular=True),
#     )
#
#     for b in e.hilbert_space_sectors():
#         if b.hamming_weight != e.number_sites // 2:
#             continue
#         b.build()
#         h = ls.Operator(e, b)
#
#         x = rng.random(size=b.number_states)
#         y = h @ x
#
#         for k in [1, 2, 10, 100]:
#             indices = rng.integers(b.number_states, size=k)
#             states = b.states[indices]
#
#             matrix, other_states = h.to_partial_csr(states)
#             assert np.allclose(matrix @ x[b.index(other_states)], y[b.index(states)])
#
#
# def test_basis_properties():
#     a = ls.SpinBasis(4)
#     assert a.spin_inversion is None
#     assert a.hamming_weight is None
#     assert not a.is_built
#     assert a.min_state_estimate == 0
#     assert a.max_state_estimate == 15
#     assert not a.has_permutation_symmetries
#     assert not a.requires_projection
#     assert a.is_state_index_identity
#
#     b = ls.SpinBasis(3, hamming_weight=2)
#     assert b.spin_inversion is None
#     assert b.hamming_weight == 2
#     assert not b.is_built
#     assert b.min_state_estimate == 0b11
#     assert b.max_state_estimate == 0b110
#     assert not b.has_permutation_symmetries
#     assert not b.requires_projection
#     assert not b.is_state_index_identity
#
#     c = ls.SpinBasis(3, spin_inversion=-1)
#     assert c.spin_inversion == -1
#     assert c.hamming_weight is None
#     assert not c.is_built
#     assert c.min_state_estimate == 0b0
#     assert c.max_state_estimate == 0b011
#     assert not c.has_permutation_symmetries
#     assert c.requires_projection
#     assert c.is_state_index_identity
#
#
# def test_basis_build():
#     a = ls.SpinBasis(2)
#     a.build()
#     assert np.array_equal(a.states, np.array([0, 1, 2, 3]))
#
#
# def test_operator_construction():
#     _ = ls.Operator(ls.Expr("2 (σ⁺₀ σ⁻₁ + σ⁺₁ σ⁻₀) + σᶻ₀ σᶻ₁"))
#
#
# def test_operator_matvec():
#     basis = ls.SpinBasis(2)
#     basis.build()
#     matrix = ls.Operator(ls.Expr("2 (σ⁺₀ σ⁻₁ + σ⁺₁ σ⁻₀) + σᶻ₀ σᶻ₁"), basis)
#
#     x = np.random.rand(matrix.shape[0])
#     y = matrix @ x
#     z = np.array([[1, 0, 0, 0], [0, -1, 2, 0], [0, 2, -1, 0], [0, 0, 0, 1]]) @ x
#     assert np.allclose(y, z)
#
#     x = np.random.rand(matrix.shape[0]) + np.random.rand(matrix.shape[0]) * 1j
#     y = matrix @ x
#     z = np.array([[1, 0, 0, 0], [0, -1, 2, 0], [0, 2, -1, 0], [0, 0, 0, 1]]) @ x
#     assert np.allclose(y, z)
#
#
# def test_operator_to_csr():
#     basis = ls.SpinBasis(2)
#     basis.build()
#     matrix = ls.Operator(ls.Expr("2 (σ⁺₀ σ⁻₁ + σ⁺₁ σ⁻₀) + σᶻ₀ σᶻ₁"), basis)
#
#     h = np.array([[1, 0, 0, 0], [0, -1, 2, 0], [0, 2, -1, 0], [0, 0, 0, 1]])
#     assert np.array_equal(matrix.diag_to_array(), h.diagonal())
#     assert np.array_equal(matrix.off_diag_to_csr().todense(), h - np.diag(h.diagonal()))
#
#
# def test_randomized_matvec():
#     test_data = os.environ.get("TEST_DATA")
#     if test_data is not None:
#         for filename in glob.glob(f"{test_data}/*_expr.json"):
#             print(filename)
#             with open(filename, "r") as f:
#                 matrix = ls.Operator.from_json(f.read())
#             with open(filename.replace("expr", "arrays"), "r") as f:
#                 o = json.load(f)
#                 r = np.asarray(o["states"], dtype=np.uint64)
#                 x = np.asarray(o["x_real"]) + np.asarray(o["x_imag"]) * 1j
#                 y = np.asarray(o["y_real"]) + np.asarray(o["y_imag"]) * 1j
#
#             if not matrix.basis.permutation_group.is_trivial:
#                 for p in matrix.basis.permutation_group.elements:
#                     assert (
#                         matrix.expression.replace_indices(
#                             dict((i, p(i)) for i in range(matrix.basis.number_sites))
#                         )
#                         == matrix.expression
#                     )
#
#             matrix.basis.build()
#             assert np.array_equal(matrix.basis.states, r)
#             # print(matrix.basis.states.shape)
#             for k in range(10):
#                 # z = np.zeros_like(x)
#                 z = matrix @ x
#                 # matrix.apply_to_state_vector(x, out=z)
#                 if not np.allclose(z, y):
#                     print(k)
#                     for i in range(len(r)):
#                         if not np.isclose(z[i], y[i]):
#                             print(i, z[i], y[i])
#                 assert np.allclose(matrix @ x, y)
#
#             diag = matrix.diag_to_array()
#             off_diag = matrix.off_diag_to_csr()
#             assert np.allclose(off_diag @ x + diag * x, y)
#
#
# # def test_kagome_ground_state():
# #     expr = ls.Expr(
# #         "1.0 σᶻ₀ σᶻ₁ + 1.0 σᶻ₀ σᶻ₃ + 1.0 σᶻ₀ σᶻ₈ + 1.0 σᶻ₀ σᶻ₁₀ + 2.0 σ⁺₀ σ⁻₁ + 2.0 σ⁺₀ σ⁻₃ + 2.0 σ⁺₀ σ⁻₈ + 2.0 σ⁺₀ σ⁻₁₀ + 2.0 σ⁻₀ σ⁺₁ + 2.0 σ⁻₀ σ⁺₃ + 2.0 σ⁻₀ σ⁺₈ + 2.0 σ⁻₀ σ⁺₁₀ + 1.0 σᶻ₁ σᶻ₂ + 0.8 σᶻ₁ σᶻ₃ + 0.8 σᶻ₁ σᶻ₉ + 2.0 σ⁺₁ σ⁻₂ + 1.6 σ⁺₁ σ⁻₃ + 1.6 σ⁺₁ σ⁻₉ + 2.0 σ⁻₁ σ⁺₂ + 1.6 σ⁻₁ σ⁺₃ + 1.6 σ⁻₁ σ⁺₉ + 1.0 σᶻ₂ σᶻ₄ + 1.0 σᶻ₂ σᶻ₉ + 1.0 σᶻ₂ σᶻ₁₀ + 2.0 σ⁺₂ σ⁻₄ + 2.0 σ⁺₂ σ⁻₉ + 2.0 σ⁺₂ σ⁻₁₀ + 2.0 σ⁻₂ σ⁺₄ + 2.0 σ⁻₂ σ⁺₉ + 2.0 σ⁻₂ σ⁺₁₀ + 1.0 σᶻ₃ σᶻ₅ + 0.8 σᶻ₃ σᶻ₁₁ + 2.0 σ⁺₃ σ⁻₅ + 1.6 σ⁺₃ σ⁻₁₁ + 2.0 σ⁻₃ σ⁺₅ + 1.6 σ⁻₃ σ⁺₁₁ + 0.8 σᶻ₄ σᶻ₆ + 1.0 σᶻ₄ σᶻ₇ + 0.8 σᶻ₄ σᶻ₁₀ + 1.6 σ⁺₄ σ⁻₆ + 2.0 σ⁺₄ σ⁻₇ + 1.6 σ⁺₄ σ⁻₁₀ + 1.6 σ⁻₄ σ⁺₆ + 2.0 σ⁻₄ σ⁺₇ + 1.6 σ⁻₄ σ⁺₁₀ + 1.0 σᶻ₅ σᶻ₆ + 1.0 σᶻ₅ σᶻ₈ + 1.0 σᶻ₅ σᶻ₁₁ + 2.0 σ⁺₅ σ⁻₆ + 2.0 σ⁺₅ σ⁻₈ + 2.0 σ⁺₅ σ⁻₁₁ + 2.0 σ⁻₅ σ⁺₆ + 2.0 σ⁻₅ σ⁺₈ + 2.0 σ⁻₅ σ⁺₁₁ + 1.0 σᶻ₆ σᶻ₇ + 0.8 σᶻ₆ σᶻ₈ + 2.0 σ⁺₆ σ⁻₇ + 1.6 σ⁺₆ σ⁻₈ + 2.0 σ⁻₆ σ⁺₇ + 1.6 σ⁻₆ σ⁺₈ + 1.0 σᶻ₇ σᶻ₉ + 1.0 σᶻ₇ σᶻ₁₁ + 2.0 σ⁺₇ σ⁻₉ + 2.0 σ⁺₇ σ⁻₁₁ + 2.0 σ⁻₇ σ⁺₉ + 2.0 σ⁻₇ σ⁺₁₁ + 0.8 σᶻ₈ σᶻ₁₀ + 1.6 σ⁺₈ σ⁻₁₀ + 1.6 σ⁻₈ σ⁺₁₀ + 0.8 σᶻ₉ σᶻ₁₁ + 1.6 σ⁺₉ σ⁻₁₁ + 1.6 σ⁻₉ σ⁺₁₁"
# #     )
# #     assert expr.is_real
# #     basis = ls.SpinBasis(12)
# #     basis.build()
# #     assert basis.is_real
# #
# #     hamiltonian = ls.Operator(expr, basis)
# #     # energy, _ = scipy.sparse.linalg.eigsh(hamiltonian, k=1, which="SA")
# #     # assert energy == approx(-19.95338528)
# #
# #     hamiltonian.build_matrix()
# #     energy, _ = scipy.sparse.linalg.eigsh(hamiltonian, k=1, which="SA")
# #     assert energy == approx(-19.95338528)
#
#
# def test_axpy():
#     for k in range(10):
#         n = np.random.randint(1, 1000)
#         alpha = np.random.rand() + np.random.rand() * 1j
#         x = np.random.rand(n) + np.random.rand(n) * 1j
#         y = np.random.rand(n) + np.random.rand(n) * 1j
#         z = alpha * x + y
#         ls._axpy(alpha, x, y)
#         assert np.allclose(z, y)
#
#
# def test_expr_permutation_group():
#     def heisenberg_on_graph(g):
#         return ls.Expr("2 (σ⁺₀ σ⁻₁ + σ⁺₁ σ⁻₀) + σᶻ₀ σᶻ₁", sites=g)
#
#     def check(g):
#         p = ls.PermutationGroup(list(map(ls.Permutation, g.get_automorphisms_vf2())))
#         assert heisenberg_on_graph(g).permutation_group() == p
#
#     # Chains
#     for n in [3, 4, 5, 6, 10, 25]:
#         check(ig.Graph.Ring(n=n, circular=True))
#         check(ig.Graph.Ring(n=n, circular=False))
#
#     # Square lattice
#     check(ig.Graph.Lattice(dim=[3, 3], circular=True))
#     check(ig.Graph.Lattice(dim=[3, 3], circular=False))
#
#     # Disconnected graph
#     check(ig.Graph(n=4, edges=[[0, 1], [2, 3]]))
#
#     # Complete graphs
#     for n in [3, 4, 5]:
#         check(ig.Graph.Full(n=n, directed=False, loops=False))
#
#     # Trees
#     check(ig.Graph.Tree(n=7, children=2))
#     check(ig.Graph.Tree(n=5, children=3))
#
#
# # long test
# def test_symmetries_kagome():
#     def find_ground_state_energy(basis_json: str, expression_str: str):
#         basis = ls.Basis(basis_json)
#         basis.build()
#         expression = ls.Expr(expression_str, particle="spin-1/2")
#         hamiltonian = ls.Operator(expression, basis)
#         hamiltonian.build_matrix()
#         energies, _ = scipy.sparse.linalg.eigsh(hamiltonian, k=1, which="SA")
#         return energies[0]
#
#     # test non-trivial sectors
#
#     basis_json_sym = '{"particle": "spin-1/2", "number_sites": 24, "hamming_weight": 12, "spin_inversion": null, "symmetries": [[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23], "0/1"], [[2, 3, 0, 1, 5, 4, 8, 9, 6, 7, 11, 10, 14, 15, 12, 13, 17, 16, 20, 21, 18, 19, 23, 22], "0/1"], [[6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5], "1/2"], [[8, 9, 6, 7, 11, 10, 14, 15, 12, 13, 17, 16, 20, 21, 18, 19, 23, 22, 2, 3, 0, 1, 5, 4], "1/2"], [[12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], "0/1"], [[14, 15, 12, 13, 17, 16, 20, 21, 18, 19, 23, 22, 2, 3, 0, 1, 5, 4, 8, 9, 6, 7, 11, 10], "0/1"], [[18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17], "1/2"], [[20, 21, 18, 19, 23, 22, 2, 3, 0, 1, 5, 4, 8, 9, 6, 7, 11, 10, 14, 15, 12, 13, 17, 16], "1/2"]]}'
#     expression_str_sym = "σᶻ₀ σᶻ₁ + σᶻ₀ σᶻ₃ + σᶻ₀ σᶻ₄ + σᶻ₀ σᶻ₂₂ + 2.0 σ⁺₀ σ⁻₁ + 2.0 σ⁺₀ σ⁻₃ + 2.0 σ⁺₀ σ⁻₄ + 2.0 σ⁺₀ σ⁻₂₂ + 2.0 σ⁻₀ σ⁺₁ + 2.0 σ⁻₀ σ⁺₃ + 2.0 σ⁻₀ σ⁺₄ + 2.0 σ⁻₀ σ⁺₂₂ + σᶻ₁ σᶻ₂ + σᶻ₁ σᶻ₄ + σᶻ₁ σᶻ₂₃ + 2.0 σ⁺₁ σ⁻₂ + 2.0 σ⁺₁ σ⁻₄ + 2.0 σ⁺₁ σ⁻₂₃ + 2.0 σ⁻₁ σ⁺₂ + 2.0 σ⁻₁ σ⁺₄ + 2.0 σ⁻₁ σ⁺₂₃ + σᶻ₂ σᶻ₃ + σᶻ₂ σᶻ₅ + σᶻ₂ σᶻ₂₃ + 2.0 σ⁺₂ σ⁻₃ + 2.0 σ⁺₂ σ⁻₅ + 2.0 σ⁺₂ σ⁻₂₃ + 2.0 σ⁻₂ σ⁺₃ + 2.0 σ⁻₂ σ⁺₅ + 2.0 σ⁻₂ σ⁺₂₃ + σᶻ₃ σᶻ₅ + σᶻ₃ σᶻ₂₂ + 2.0 σ⁺₃ σ⁻₅ + 2.0 σ⁺₃ σ⁻₂₂ + 2.0 σ⁻₃ σ⁺₅ + 2.0 σ⁻₃ σ⁺₂₂ + σᶻ₄ σᶻ₆ + σᶻ₄ σᶻ₉ + 2.0 σ⁺₄ σ⁻₆ + 2.0 σ⁺₄ σ⁻₉ + 2.0 σ⁻₄ σ⁺₆ + 2.0 σ⁻₄ σ⁺₉ + σᶻ₅ σᶻ₇ + σᶻ₅ σᶻ₈ + 2.0 σ⁺₅ σ⁻₇ + 2.0 σ⁺₅ σ⁻₈ + 2.0 σ⁻₅ σ⁺₇ + 2.0 σ⁻₅ σ⁺₈ + σᶻ₆ σᶻ₇ + σᶻ₆ σᶻ₉ + σᶻ₆ σᶻ₁₀ + 2.0 σ⁺₆ σ⁻₇ + 2.0 σ⁺₆ σ⁻₉ + 2.0 σ⁺₆ σ⁻₁₀ + 2.0 σ⁻₆ σ⁺₇ + 2.0 σ⁻₆ σ⁺₉ + 2.0 σ⁻₆ σ⁺₁₀ + σᶻ₇ σᶻ₈ + σᶻ₇ σᶻ₁₀ + 2.0 σ⁺₇ σ⁻₈ + 2.0 σ⁺₇ σ⁻₁₀ + 2.0 σ⁻₇ σ⁺₈ + 2.0 σ⁻₇ σ⁺₁₀ + σᶻ₈ σᶻ₉ + σᶻ₈ σᶻ₁₁ + 2.0 σ⁺₈ σ⁻₉ + 2.0 σ⁺₈ σ⁻₁₁ + 2.0 σ⁻₈ σ⁺₉ + 2.0 σ⁻₈ σ⁺₁₁ + σᶻ₉ σᶻ₁₁ + 2.0 σ⁺₉ σ⁻₁₁ + 2.0 σ⁻₉ σ⁺₁₁ + σᶻ₁₀ σᶻ₁₂ + σᶻ₁₀ σᶻ₁₅ + 2.0 σ⁺₁₀ σ⁻₁₂ + 2.0 σ⁺₁₀ σ⁻₁₅ + 2.0 σ⁻₁₀ σ⁺₁₂ + 2.0 σ⁻₁₀ σ⁺₁₅ + σᶻ₁₁ σᶻ₁₃ + σᶻ₁₁ σᶻ₁₄ + 2.0 σ⁺₁₁ σ⁻₁₃ + 2.0 σ⁺₁₁ σ⁻₁₄ + 2.0 σ⁻₁₁ σ⁺₁₃ + 2.0 σ⁻₁₁ σ⁺₁₄ + σᶻ₁₂ σᶻ₁₃ + σᶻ₁₂ σᶻ₁₅ + σᶻ₁₂ σᶻ₁₆ + 2.0 σ⁺₁₂ σ⁻₁₃ + 2.0 σ⁺₁₂ σ⁻₁₅ + 2.0 σ⁺₁₂ σ⁻₁₆ + 2.0 σ⁻₁₂ σ⁺₁₃ + 2.0 σ⁻₁₂ σ⁺₁₅ + 2.0 σ⁻₁₂ σ⁺₁₆ + σᶻ₁₃ σᶻ₁₄ + σᶻ₁₃ σᶻ₁₆ + 2.0 σ⁺₁₃ σ⁻₁₄ + 2.0 σ⁺₁₃ σ⁻₁₆ + 2.0 σ⁻₁₃ σ⁺₁₄ + 2.0 σ⁻₁₃ σ⁺₁₆ + σᶻ₁₄ σᶻ₁₅ + σᶻ₁₄ σᶻ₁₇ + 2.0 σ⁺₁₄ σ⁻₁₅ + 2.0 σ⁺₁₄ σ⁻₁₇ + 2.0 σ⁻₁₄ σ⁺₁₅ + 2.0 σ⁻₁₄ σ⁺₁₇ + σᶻ₁₅ σᶻ₁₇ + 2.0 σ⁺₁₅ σ⁻₁₇ + 2.0 σ⁻₁₅ σ⁺₁₇ + σᶻ₁₆ σᶻ₁₈ + σᶻ₁₆ σᶻ₂₁ + 2.0 σ⁺₁₆ σ⁻₁₈ + 2.0 σ⁺₁₆ σ⁻₂₁ + 2.0 σ⁻₁₆ σ⁺₁₈ + 2.0 σ⁻₁₆ σ⁺₂₁ + σᶻ₁₇ σᶻ₁₉ + σᶻ₁₇ σᶻ₂₀ + 2.0 σ⁺₁₇ σ⁻₁₉ + 2.0 σ⁺₁₇ σ⁻₂₀ + 2.0 σ⁻₁₇ σ⁺₁₉ + 2.0 σ⁻₁₇ σ⁺₂₀ + σᶻ₁₈ σᶻ₁₉ + σᶻ₁₈ σᶻ₂₁ + σᶻ₁₈ σᶻ₂₂ + 2.0 σ⁺₁₈ σ⁻₁₉ + 2.0 σ⁺₁₈ σ⁻₂₁ + 2.0 σ⁺₁₈ σ⁻₂₂ + 2.0 σ⁻₁₈ σ⁺₁₉ + 2.0 σ⁻₁₈ σ⁺₂₁ + 2.0 σ⁻₁₈ σ⁺₂₂ + σᶻ₁₉ σᶻ₂₀ + σᶻ₁₉ σᶻ₂₂ + 2.0 σ⁺₁₉ σ⁻₂₀ + 2.0 σ⁺₁₉ σ⁻₂₂ + 2.0 σ⁻₁₉ σ⁺₂₀ + 2.0 σ⁻₁₉ σ⁺₂₂ + σᶻ₂₀ σᶻ₂₁ + σᶻ₂₀ σᶻ₂₃ + 2.0 σ⁺₂₀ σ⁻₂₁ + 2.0 σ⁺₂₀ σ⁻₂₃ + 2.0 σ⁻₂₀ σ⁺₂₁ + 2.0 σ⁻₂₀ σ⁺₂₃ + σᶻ₂₁ σᶻ₂₃ + 2.0 σ⁺₂₁ σ⁻₂₃ + 2.0 σ⁻₂₁ σ⁺₂₃"
#     energy_sym = find_ground_state_energy(basis_json_sym, expression_str_sym)
#     basis_json_nosym = '{"particle": "spin-1/2", "number_sites": 24, "hamming_weight": 12, "spin_inversion": null, "symmetries": []}'
#     expression_str_nosym = "σᶻ₀ σᶻ₁ + σᶻ₀ σᶻ₃ + σᶻ₀ σᶻ₄ + σᶻ₀ σᶻ₂₂ + 2.0 σ⁺₀ σ⁻₁ + 2.0 σ⁺₀ σ⁻₃ + 2.0 σ⁺₀ σ⁻₄ + 2.0 σ⁺₀ σ⁻₂₂ + 2.0 σ⁻₀ σ⁺₁ + 2.0 σ⁻₀ σ⁺₃ + 2.0 σ⁻₀ σ⁺₄ + 2.0 σ⁻₀ σ⁺₂₂ + σᶻ₁ σᶻ₂ + σᶻ₁ σᶻ₄ + σᶻ₁ σᶻ₂₃ + 2.0 σ⁺₁ σ⁻₂ + 2.0 σ⁺₁ σ⁻₄ + 2.0 σ⁺₁ σ⁻₂₃ + 2.0 σ⁻₁ σ⁺₂ + 2.0 σ⁻₁ σ⁺₄ + 2.0 σ⁻₁ σ⁺₂₃ + σᶻ₂ σᶻ₃ + σᶻ₂ σᶻ₅ + σᶻ₂ σᶻ₂₃ + 2.0 σ⁺₂ σ⁻₃ + 2.0 σ⁺₂ σ⁻₅ + 2.0 σ⁺₂ σ⁻₂₃ + 2.0 σ⁻₂ σ⁺₃ + 2.0 σ⁻₂ σ⁺₅ + 2.0 σ⁻₂ σ⁺₂₃ + σᶻ₃ σᶻ₅ + σᶻ₃ σᶻ₂₂ + 2.0 σ⁺₃ σ⁻₅ + 2.0 σ⁺₃ σ⁻₂₂ + 2.0 σ⁻₃ σ⁺₅ + 2.0 σ⁻₃ σ⁺₂₂ + σᶻ₄ σᶻ₆ + σᶻ₄ σᶻ₉ + 2.0 σ⁺₄ σ⁻₆ + 2.0 σ⁺₄ σ⁻₉ + 2.0 σ⁻₄ σ⁺₆ + 2.0 σ⁻₄ σ⁺₉ + σᶻ₅ σᶻ₇ + σᶻ₅ σᶻ₈ + 2.0 σ⁺₅ σ⁻₇ + 2.0 σ⁺₅ σ⁻₈ + 2.0 σ⁻₅ σ⁺₇ + 2.0 σ⁻₅ σ⁺₈ + σᶻ₆ σᶻ₇ + σᶻ₆ σᶻ₉ + σᶻ₆ σᶻ₁₀ + 2.0 σ⁺₆ σ⁻₇ + 2.0 σ⁺₆ σ⁻₉ + 2.0 σ⁺₆ σ⁻₁₀ + 2.0 σ⁻₆ σ⁺₇ + 2.0 σ⁻₆ σ⁺₉ + 2.0 σ⁻₆ σ⁺₁₀ + σᶻ₇ σᶻ₈ + σᶻ₇ σᶻ₁₀ + 2.0 σ⁺₇ σ⁻₈ + 2.0 σ⁺₇ σ⁻₁₀ + 2.0 σ⁻₇ σ⁺₈ + 2.0 σ⁻₇ σ⁺₁₀ + σᶻ₈ σᶻ₉ + σᶻ₈ σᶻ₁₁ + 2.0 σ⁺₈ σ⁻₉ + 2.0 σ⁺₈ σ⁻₁₁ + 2.0 σ⁻₈ σ⁺₉ + 2.0 σ⁻₈ σ⁺₁₁ + σᶻ₉ σᶻ₁₁ + 2.0 σ⁺₉ σ⁻₁₁ + 2.0 σ⁻₉ σ⁺₁₁ + σᶻ₁₀ σᶻ₁₂ + σᶻ₁₀ σᶻ₁₅ + 2.0 σ⁺₁₀ σ⁻₁₂ + 2.0 σ⁺₁₀ σ⁻₁₅ + 2.0 σ⁻₁₀ σ⁺₁₂ + 2.0 σ⁻₁₀ σ⁺₁₅ + σᶻ₁₁ σᶻ₁₃ + σᶻ₁₁ σᶻ₁₄ + 2.0 σ⁺₁₁ σ⁻₁₃ + 2.0 σ⁺₁₁ σ⁻₁₄ + 2.0 σ⁻₁₁ σ⁺₁₃ + 2.0 σ⁻₁₁ σ⁺₁₄ + σᶻ₁₂ σᶻ₁₃ + σᶻ₁₂ σᶻ₁₅ + σᶻ₁₂ σᶻ₁₆ + 2.0 σ⁺₁₂ σ⁻₁₃ + 2.0 σ⁺₁₂ σ⁻₁₅ + 2.0 σ⁺₁₂ σ⁻₁₆ + 2.0 σ⁻₁₂ σ⁺₁₃ + 2.0 σ⁻₁₂ σ⁺₁₅ + 2.0 σ⁻₁₂ σ⁺₁₆ + σᶻ₁₃ σᶻ₁₄ + σᶻ₁₃ σᶻ₁₆ + 2.0 σ⁺₁₃ σ⁻₁₄ + 2.0 σ⁺₁₃ σ⁻₁₆ + 2.0 σ⁻₁₃ σ⁺₁₄ + 2.0 σ⁻₁₃ σ⁺₁₆ + σᶻ₁₄ σᶻ₁₅ + σᶻ₁₄ σᶻ₁₇ + 2.0 σ⁺₁₄ σ⁻₁₅ + 2.0 σ⁺₁₄ σ⁻₁₇ + 2.0 σ⁻₁₄ σ⁺₁₅ + 2.0 σ⁻₁₄ σ⁺₁₇ + σᶻ₁₅ σᶻ₁₇ + 2.0 σ⁺₁₅ σ⁻₁₇ + 2.0 σ⁻₁₅ σ⁺₁₇ + σᶻ₁₆ σᶻ₁₈ + σᶻ₁₆ σᶻ₂₁ + 2.0 σ⁺₁₆ σ⁻₁₈ + 2.0 σ⁺₁₆ σ⁻₂₁ + 2.0 σ⁻₁₆ σ⁺₁₈ + 2.0 σ⁻₁₆ σ⁺₂₁ + σᶻ₁₇ σᶻ₁₉ + σᶻ₁₇ σᶻ₂₀ + 2.0 σ⁺₁₇ σ⁻₁₉ + 2.0 σ⁺₁₇ σ⁻₂₀ + 2.0 σ⁻₁₇ σ⁺₁₉ + 2.0 σ⁻₁₇ σ⁺₂₀ + σᶻ₁₈ σᶻ₁₉ + σᶻ₁₈ σᶻ₂₁ + σᶻ₁₈ σᶻ₂₂ + 2.0 σ⁺₁₈ σ⁻₁₉ + 2.0 σ⁺₁₈ σ⁻₂₁ + 2.0 σ⁺₁₈ σ⁻₂₂ + 2.0 σ⁻₁₈ σ⁺₁₉ + 2.0 σ⁻₁₈ σ⁺₂₁ + 2.0 σ⁻₁₈ σ⁺₂₂ + σᶻ₁₉ σᶻ₂₀ + σᶻ₁₉ σᶻ₂₂ + 2.0 σ⁺₁₉ σ⁻₂₀ + 2.0 σ⁺₁₉ σ⁻₂₂ + 2.0 σ⁻₁₉ σ⁺₂₀ + 2.0 σ⁻₁₉ σ⁺₂₂ + σᶻ₂₀ σᶻ₂₁ + σᶻ₂₀ σᶻ₂₃ + 2.0 σ⁺₂₀ σ⁻₂₁ + 2.0 σ⁺₂₀ σ⁻₂₃ + 2.0 σ⁻₂₀ σ⁺₂₁ + 2.0 σ⁻₂₀ σ⁺₂₃ + σᶻ₂₁ σᶻ₂₃ + 2.0 σ⁺₂₁ σ⁻₂₃ + 2.0 σ⁻₂₁ σ⁺₂₃"
#
#     energy_nosym = find_ground_state_energy(basis_json_nosym, expression_str_nosym)
#     assert energy_nosym == approx(energy_sym)
#
#
# def notest_issue_pim_1():
#     n = 10
#     expr = ls.Expr("Sx0 Sx1 + Sy0 Sy1 + Sz0 Sz1", sites=ig.Graph.Ring(n=n, circular=True))
#     translation = ls.Permutation([(1 + i) % n for i in range(n)])
#     for k in range(8):
#         b = ls.SpinBasis(number_spins=n, symmetries=[(translation, ls.Rational(k, n))])
#         b.build()
#         h = ls.Operator(expr, b)
#         energies, _ = scipy.sparse.linalg.eigsh(h, k=4)
#         # print(energies)
#     # print("done")
#
#
# def notest_issue_pim_2():
#     n = 3
#     tx = ls.Permutation([n * ((i + 1) % n) + j for i in range(n) for j in range(n)])
#     ty = ls.Permutation([n * i + ((j + 1) % n) for i in range(n) for j in range(n)])
#
#     sites = list(range(n * n))
#     edges = []
#     for i in range(n * n):
#         edges.append((sites[i], sites[tx(i)]))
#         edges.append((sites[i], sites[ty(i)]))
#
#     expr = ls.Expr("Sx0 Sx1 + Sy0 Sy1 + Sz0 Sz1", sites=edges)
#     for kx in range(n):
#         for ky in range(1, n):
#             g = [(tx, ls.Rational(kx, n)), (ty, ls.Rational(ky, n))]
#             b = ls.SpinBasis(number_spins=n * n, symmetries=g)
#             b.build()
#             h = ls.Operator(expr, b)
#             energies, _ = scipy.sparse.linalg.eigsh(h, k=4)
#             # print(energies)
#     # print("done")


# test_issue_pim_2()

# def test_permutation_construction():
#     p = ls.Permutation([0, 1, 2])
#     assert p.periodicity == 1
#     assert p.permutation == [0, 1, 2]
#     with raises(ValueError, match=r".*permutation.*"):
#         ls.Permutation([1, 2, 4])


# def sum1(xs):
#     s = None
#     for x in xs:
#         if s is not None:
#             s += x
#         else:
#             s = x
#     return s
#
#
# def test_symmetry():
#     a = ls.Symmetry([0, 1, 2], sector=0)
#     assert a.sector == 0
#     assert len(a) == 3
#     assert a.permutation.tolist() == [0, 1, 2]
#     del a
#
#     a = ls.Symmetry(list(range(10)), sector=123)
#     assert a.phase == 0
#
#     a = ls.Symmetry(list(range(10)), sector=-1)
#     assert a.phase == 0
#
#
# def test_symmetries():
#     a = ls.Symmetry([1, 2, 3, 0], sector=0)
#     b = ls.Symmetry([3, 2, 1, 0], sector=0)
#     c = ls.Symmetries([a, b])
#     assert len(c) == 8
#
#     # with raises(SystemError):
#     #     ls.Symmetries([ls.Symmetry([1, 2, 0], sector=1), ls.Symmetry([1, 2, 0], sector=2)])
#
#     print(c.compile())
#
#
# def test_index():
#     basis = ls.SpinBasis(4)
#     basis.build()
#     assert np.array_equal(basis.index(basis.states), basis.states)
#     assert np.array_equal(basis.index(basis.states[2]), 2)
#
#
# def test_kagome_symmetries():
#     expr = ls.Expr(
#         "1.0 σᶻ₀ σᶻ₁ + 1.0 σᶻ₀ σᶻ₃ + 1.0 σᶻ₀ σᶻ₈ + 1.0 σᶻ₀ σᶻ₁₀ + 2.0 σ⁺₀ σ⁻₁ + 2.0 σ⁺₀ σ⁻₃ + 2.0 σ⁺₀ σ⁻₈ + 2.0 σ⁺₀ σ⁻₁₀ + 2.0 σ⁻₀ σ⁺₁ + 2.0 σ⁻₀ σ⁺₃ + 2.0 σ⁻₀ σ⁺₈ + 2.0 σ⁻₀ σ⁺₁₀ + 1.0 σᶻ₁ σᶻ₂ + 0.8 σᶻ₁ σᶻ₃ + 0.8 σᶻ₁ σᶻ₉ + 2.0 σ⁺₁ σ⁻₂ + 1.6 σ⁺₁ σ⁻₃ + 1.6 σ⁺₁ σ⁻₉ + 2.0 σ⁻₁ σ⁺₂ + 1.6 σ⁻₁ σ⁺₃ + 1.6 σ⁻₁ σ⁺₉ + 1.0 σᶻ₂ σᶻ₄ + 1.0 σᶻ₂ σᶻ₉ + 1.0 σᶻ₂ σᶻ₁₀ + 2.0 σ⁺₂ σ⁻₄ + 2.0 σ⁺₂ σ⁻₉ + 2.0 σ⁺₂ σ⁻₁₀ + 2.0 σ⁻₂ σ⁺₄ + 2.0 σ⁻₂ σ⁺₉ + 2.0 σ⁻₂ σ⁺₁₀ + 1.0 σᶻ₃ σᶻ₅ + 0.8 σᶻ₃ σᶻ₁₁ + 2.0 σ⁺₃ σ⁻₅ + 1.6 σ⁺₃ σ⁻₁₁ + 2.0 σ⁻₃ σ⁺₅ + 1.6 σ⁻₃ σ⁺₁₁ + 0.8 σᶻ₄ σᶻ₆ + 1.0 σᶻ₄ σᶻ₇ + 0.8 σᶻ₄ σᶻ₁₀ + 1.6 σ⁺₄ σ⁻₆ + 2.0 σ⁺₄ σ⁻₇ + 1.6 σ⁺₄ σ⁻₁₀ + 1.6 σ⁻₄ σ⁺₆ + 2.0 σ⁻₄ σ⁺₇ + 1.6 σ⁻₄ σ⁺₁₀ + 1.0 σᶻ₅ σᶻ₆ + 1.0 σᶻ₅ σᶻ₈ + 1.0 σᶻ₅ σᶻ₁₁ + 2.0 σ⁺₅ σ⁻₆ + 2.0 σ⁺₅ σ⁻₈ + 2.0 σ⁺₅ σ⁻₁₁ + 2.0 σ⁻₅ σ⁺₆ + 2.0 σ⁻₅ σ⁺₈ + 2.0 σ⁻₅ σ⁺₁₁ + 1.0 σᶻ₆ σᶻ₇ + 0.8 σᶻ₆ σᶻ₈ + 2.0 σ⁺₆ σ⁻₇ + 1.6 σ⁺₆ σ⁻₈ + 2.0 σ⁻₆ σ⁺₇ + 1.6 σ⁻₆ σ⁺₈ + 1.0 σᶻ₇ σᶻ₉ + 1.0 σᶻ₇ σᶻ₁₁ + 2.0 σ⁺₇ σ⁻₉ + 2.0 σ⁺₇ σ⁻₁₁ + 2.0 σ⁻₇ σ⁺₉ + 2.0 σ⁻₇ σ⁺₁₁ + 0.8 σᶻ₈ σᶻ₁₀ + 1.6 σ⁺₈ σ⁻₁₀ + 1.6 σ⁻₈ σ⁺₁₀ + 0.8 σᶻ₉ σᶻ₁₁ + 1.6 σ⁺₉ σ⁻₁₁ + 1.6 σ⁻₉ σ⁺₁₁"
#     )
#     # top_shift = ls.Symmetry([5, 6, 7, 8, 9, 0, 1, 2, 3, 4, 11, 10], sector=0)
#     right_shift = ls.Symmetry([2, 10, 0, 4, 3, 7, 11, 5, 9, 8, 1, 6], sector=1)
#     assert right_shift.periodicity == 2
#     assert right_shift.phase == 0.5
#     assert expr == expr.replace_indices(dict(zip(range(12), right_shift.permutation)))
#     symmetries = ls.Symmetries([right_shift])
#     basis = ls.SpinBasis(
#         symmetries=symmetries, number_spins=12, hamming_weight=6, spin_inversion=None
#     )
#     basis.build()
#     # print(basis.states)
#     # print(basis.state_info(basis.states))
#     hamiltonian = ls.Operator(basis, expr)
#     assert basis.is_real
#     assert expr.is_real
#     assert hamiltonian.is_real
#     energy, _ = scipy.sparse.linalg.eigsh(hamiltonian, k=1, which="SA")
#     assert energy == approx(-19.95338528)
#
#
# def test_operator_apply():
#     basis = ls.SpinBasis(2)
#     basis.build()
#     expr = ls.Expr("1.0 σᶻ₀ σᶻ₁ + 2.0 σ⁺₀ σ⁻₁ + 2.0 σ⁻₀ σ⁺₁")
#     hamiltonian = ls.Operator(basis, expr)
#     assert len(hamiltonian.apply_off_diag_to_basis_state(basis.states[1])) == 1
#
#
# def test_simple_spin_expr():
#     sp = np.array([[0, 1], [0, 0]])
#     sm = np.array([[0, 0], [1, 0]])
#     sz = np.diag([1, -1])
#     s0 = np.eye(2)
#
#     def check(number_spins, expression, matrix_ref):
#         basis = ls.SpinBasis(number_spins)
#         basis.build()
#         operator = ls.Operator(basis, ls.Expr(expression))
#         cols = []
#         for i in range(basis.number_states):
#             v = np.zeros(basis.number_states)
#             v[i] = 1
#             cols.append(operator @ v)
#         matrix = np.vstack(cols).T
#         assert matrix.tolist() == matrix_ref.tolist()
#         matrix = operator.to_csr().todense()
#         assert matrix.tolist() == matrix_ref.tolist()
#
#     check(1, "σ⁻₀", sm)
#     check(1, "σ⁺₀", sp)
#     check(1, "σᶻ₀", sz)
#     check(2, "σ⁺₀ σ⁻₁", np.kron(sm, sp))
#     check(2, "σ⁺₀ σᶻ₁", np.kron(sz, sp))
#     check(2, "σ⁻₁", np.kron(sm, s0))
#     check(3, "σ⁺₀ σᶻ₁ σ⁻₂", np.kron(sm, np.kron(sz, sp)))
#     check(3, "σ⁺₀ σ⁻₂", np.kron(sm, np.kron(s0, sp)))
#
#     check(1, "σʸ₀", -1j * sp + 1j * sm)
#     check(3, "3im σ⁺₀ σᶻ₁ σʸ₂", 3j * np.kron((1j * sm - 1j * sp), np.kron(sz, sp)))
#
#
# def test_prepare_hphi():
#     basis = ls.SpinBasis(2)
#     expr = ls.Expr("σᶻ₀ σᶻ₁ + 2 σ⁺₀ σ⁻₁ + 2 σ⁻₀ σ⁺₁")
#     hamiltonian = ls.Operator(basis, expr)
#     hamiltonian.prepare_inputs_for_hphi("/tmp/lattice-symmetries-python/hphi")
#
#
# def test_prepare_mvmc():
#     # basis = ls.SpinBasis(4)
#     # expr = ls.Expr("σᶻ₀ σᶻ₁ + 2 σ⁺₀ σ⁻₁ + 2 σ⁻₀ σ⁺₁")
#     # expr = sum((expr.replace_indices({0: i, 1: j}) for (i, j) in [(0, 1), (1, 2), (2, 3), (3, 0)]), ls.Expr(""))
#     basis = ls.SpinfulFermionBasis(4, number_particles=4)
#
#     hopping = ls.Expr("- (c†₁↑ c₀↑ + c†₀↑ c₁↑ + c†₁↓ c₀↓ + c†₀↓ c₁↓)")
#     coulomb = ls.Expr("4.0 n₀↑ n₀↓")
#     expr = hopping + coulomb
#     for i, j in [(1, 2), (2, 3), (3, 0)]:
#         expr += hopping.replace_indices({0: i, 1: j})
#     for i in [1, 2, 3]:
#         expr += coulomb.replace_indices({0: i})
#     print(expr)
#     hamiltonian = ls.Operator(basis, expr)
#     hamiltonian.prepare_inputs_for_mvmc("/tmp/lattice-symmetries-python/mvmc")
#
#
# def test_anisotropic_kagome_9():
#     # fmt: off
#     nearest = [
#         (0, 1), (1, 2), (2, 0),
#         (3, 4), (4, 5), (5, 3),
#         (6, 7), (7, 8), (8, 6),
#     ]
#     next_nearest = [
#         (0, 5), (0, 7),
#         (1, 3), (1, 8),
#         (2, 4), (2, 6),
#         (3, 8),
#         (4, 6),
#         (5, 7),
#     ]
#     # fmt: on
#
#     basis = ls.SpinfulFermionBasis(number_sites=9, number_particles=3)
#     basis.build()
#     hopping = lambda i, j: ls.Expr("c†₁↑ c₀↑ + c†₀↑ c₁↑ + c†₁↓ c₀↓ + c†₀↓ c₁↓").replace_indices(
#         {0: i, 1: j}
#     )
#     coulomb = lambda i: ls.Expr("n₀↑ n₀↓").replace_indices({0: i})
#
#     t1 = -0.3251
#     t2 = 0.0845
#     U = 2.8
#     expr = (
#         t1 * sum1((hopping(i, j) for i, j in nearest))
#         + t2 * sum1((hopping(i, j) for i, j in next_nearest))
#         + U * sum1((coulomb(i) for i in range(9)))
#     )
#     print(expr)
#     hamiltonian = ls.Operator(basis, expr)
#     hamiltonian.prepare_inputs_for_hphi("/tmp/kagome")
#     energy, state = scipy.sparse.linalg.eigsh(hamiltonian, k=1, which="SA", tol=1e-6)
#     print(energy)
#
#
# def notest_vs_hphi():
#     prefix = "../../test"
#     for folder in os.listdir(prefix):
#         print(folder)
#         config = ls.load_yaml_config(os.path.join(prefix, folder, "hamiltonian.yaml"))
#         config.basis.build()
#         energy, state = scipy.sparse.linalg.eigsh(config.hamiltonian, k=1, which="SA", tol=1e-6)
#         with open(os.path.join(prefix, folder, "HPhi", "output", "zvo_energy.dat")) as f:
#             for line in f.readlines():
#                 if line.startswith("Energy"):
#                     ref_energy = float(line.strip().split(" ")[-1])
#         print(energy, ref_energy)
#         assert ref_energy is not None
#         assert energy == approx(ref_energy)
#
#
# def test_apply_off_diag_projected():
#     basis1 = ls.SpinBasis(4)
#     basis1.build()
#     expr = ls.Expr(
#         "σᶻ₀ σᶻ₁ + 2 σ⁺₀ σ⁻₁ + 2 σ⁻₀ σ⁺₁",
#         sites=[(0, 1), (1, 2), (2, 3), (3, 0)],
#     )
#     hamiltonian1 = ls.Operator(basis1, expr)
#     assert hamiltonian1.apply_off_diag_to_basis_state(int("0101", base=2)) == [
#         ((2 + 0j), 6),
#         ((2 + 0j), 3),
#         ((2 + 0j), 12),
#         ((2 + 0j), 9),
#     ]
#
#     group = ls.Symmetries([ls.Symmetry([1, 2, 3, 0], sector=0)])
#     basis = ls.SpinBasis(4, symmetries=group)
#     basis.build()
#     hamiltonian = ls.Operator(basis, expr)
#     for c, x in hamiltonian.apply_off_diag_to_basis_state(int("0101", base=2)):
#         assert x == 3
#         assert c == approx(math.sqrt(2))
#
#     group = ls.Symmetries([ls.Symmetry([1, 2, 3, 0], sector=2)])
#     basis = ls.SpinBasis(4, symmetries=group)
#     basis.build()
#     hamiltonian = ls.Operator(basis, expr)
#     assert hamiltonian.apply_off_diag_to_basis_state(int("0101", base=2)) == [
#         (approx(-math.sqrt(2)), 3),
#         (approx(math.sqrt(2)), 3),
#         (approx(math.sqrt(2)), 3),
#         (approx(-math.sqrt(2)), 3),
#     ]
#
#
# # def test_matvec():
# #     if os.en
#
#
# def get_csr_hamiltonian(hamiltonian):
#     basis = hamiltonian.basis
#     states, coeffs, row_idxs = hamiltonian.apply_off_diag_to_basis_state(basis.states)
#     columns = basis.index(states)
#     off_diagonal_matrix = scipy.sparse.csr_matrix(
#         (coeffs, columns, row_idxs),
#         shape=(basis.number_states, basis.number_states),
#     )
#     diagonal = scipy.sparse.diags(hamiltonian.apply_diag_to_basis_state(basis.states))
#     return off_diagonal_matrix + diagonal
#
#
# def test_correct_imag_part():
#     expr = ls.Expr(
#         # "σᶻ₀ σᶻ₁ +
#         "2 σ⁺₀ σ⁻₁"
#         # + 2 σ⁻₀ σ⁺₁"
#         ,
#         sites=[(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)],
#     )
#     basis = ls.SpinBasis(5, hamming_weight=2)
#     basis.build()
#     hamiltonian = ls.Operator(basis, expr)
#
#     x = np.zeros(basis.number_states, dtype=complex)
#     #                  43210
#     x[basis.index(int("00011", base=2))] += 1.0
#     x[basis.index(int("00110", base=2))] += np.exp(+2j * np.pi * 1 / 5)
#     x[basis.index(int("01100", base=2))] += np.exp(+2j * np.pi * 2 / 5)
#     x[basis.index(int("11000", base=2))] += np.exp(+2j * np.pi * 3 / 5)
#     x[basis.index(int("10001", base=2))] += np.exp(+2j * np.pi * 4 / 5)
#     x /= np.linalg.norm(x)
#
#     y = np.zeros(basis.number_states, dtype=complex)
#     y[basis.index(int("00101", base=2))] += 1.0
#     y[basis.index(int("01010", base=2))] += np.exp(+2j * np.pi * 1 / 5)
#     y[basis.index(int("10100", base=2))] += np.exp(+2j * np.pi * 2 / 5)
#     y[basis.index(int("01001", base=2))] += np.exp(+2j * np.pi * 3 / 5)
#     y[basis.index(int("10010", base=2))] += np.exp(+2j * np.pi * 4 / 5)
#     y /= np.linalg.norm(y)
#
#     matrix = np.array(
#         [
#             [np.vdot(x, hamiltonian @ x), np.vdot(x, hamiltonian @ y)],
#             [np.vdot(y, hamiltonian @ x), np.vdot(y, hamiltonian @ y)],
#         ]
#     )
#
#     group = ls.Symmetries([ls.Symmetry([4, 0, 1, 2, 3], sector=1)])
#     basis = ls.SpinBasis(5, hamming_weight=2, symmetries=group)
#     hamiltonian = ls.Operator(basis, expr)
#     basis.build()
#
#     cols = []
#     for i in range(basis.number_states):
#         v = np.zeros(basis.number_states)
#         v[i] = 1
#         cols.append(hamiltonian @ v)
#     matrix2 = np.vstack(cols).T
#     matrix3 = hamiltonian.to_csr().todense()
#     print(matrix)
#     print(matrix2)
#     print(matrix3)
#     assert np.allclose(matrix, matrix2)
#     assert np.allclose(matrix, matrix3)
#
#
# def test_convert_to_csr():
#     def check(basis, expr):
#         assert expr == expr.adjoint()
#         hamiltonian = ls.Operator(basis, expr)
#         basis.build()
#         matrix = hamiltonian.to_csr()
#         assert matrix.has_canonical_format
#         assert (matrix != get_csr_hamiltonian(hamiltonian)).nnz == 0
#         np.random.seed(42)
#         for i in range(5):
#             x = np.random.rand(basis.number_states)  # + 1j * np.random.rand(basis.number_states)
#             y1 = hamiltonian @ x
#             y2 = matrix @ x
#             if not np.allclose(y1, y2):
#                 print(y1.tolist()[:10])
#                 print(y2.tolist()[:10])
#             assert np.allclose(y1, y2)
#
#     basis = ls.SpinBasis(2)
#     expr = ls.Expr("σᶻ₀ σᶻ₁ + 2 σ⁺₀ σ⁻₁ + 2 σ⁻₀ σ⁺₁")
#     check(basis, expr)
#
#     group = ls.Symmetries([ls.Symmetry([1, 2, 3, 0], sector=0)])
#     basis = ls.SpinBasis(4, symmetries=group)
#     expr = ls.Expr(
#         "σᶻ₀ σᶻ₁ + 2 σ⁺₀ σ⁻₁ + 2 σ⁻₀ σ⁺₁",
#         sites=[(0, 1), (1, 2), (2, 3), (3, 0)],
#     )
#     check(basis, expr)
#
#     basis = ls.Basis.from_json(
#         '{"hamming_weight":null,"number_spins":9,"particle":"spin-1/2","spin_inversion":null,"symmetries":[{"permutation":[0,1,2,3,4,5,6,7,8],"sector":0},{"permutation":[1,2,0,4,5,3,7,8,6],"sector":0},{"permutation":[2,0,1,5,3,4,8,6,7],"sector":0},{"permutation":[3,4,5,6,7,8,0,1,2],"sector":0},{"permutation":[4,5,3,7,8,6,1,2,0],"sector":0},{"permutation":[5,3,4,8,6,7,2,0,1],"sector":0},{"permutation":[6,7,8,0,1,2,3,4,5],"sector":0},{"permutation":[7,8,6,1,2,0,4,5,3],"sector":0},{"permutation":[8,6,7,2,0,1,5,3,4],"sector":0}]}'
#     )
#     expr = ls.Expr(
#         "-σᶻ₀ σᶻ₁ - σᶻ₀ σᶻ₂ - σᶻ₀ σᶻ₃ - σᶻ₀ σᶻ₆ - 2.0 σ⁺₀ σ⁻₁ - 2.0 σ⁺₀ σ⁻₂ - 2.0 σ⁺₀ σ⁻₃ - 2.0 σ⁺₀ σ⁻₆ - 2.0 σ⁻₀ σ⁺₁ - 2.0 σ⁻₀ σ⁺₂ - 2.0 σ⁻₀ σ⁺₃ - 2.0 σ⁻₀ σ⁺₆ - σᶻ₁ σᶻ₂ - σᶻ₁ σᶻ₄ - σᶻ₁ σᶻ₇ - 2.0 σ⁺₁ σ⁻₂ - 2.0 σ⁺₁ σ⁻₄ - 2.0 σ⁺₁ σ⁻₇ - 2.0 σ⁻₁ σ⁺₂ - 2.0 σ⁻₁ σ⁺₄ - 2.0 σ⁻₁ σ⁺₇ - σᶻ₂ σᶻ₅ - σᶻ₂ σᶻ₈ - 2.0 σ⁺₂ σ⁻₅ - 2.0 σ⁺₂ σ⁻₈ - 2.0 σ⁻₂ σ⁺₅ - 2.0 σ⁻₂ σ⁺₈ - σᶻ₃ σᶻ₄ - σᶻ₃ σᶻ₅ - σᶻ₃ σᶻ₆ - 2.0 σ⁺₃ σ⁻₄ - 2.0 σ⁺₃ σ⁻₅ - 2.0 σ⁺₃ σ⁻₆ - 2.0 σ⁻₃ σ⁺₄ - 2.0 σ⁻₃ σ⁺₅ - 2.0 σ⁻₃ σ⁺₆ - σᶻ₄ σᶻ₅ - σᶻ₄ σᶻ₇ - 2.0 σ⁺₄ σ⁻₅ - 2.0 σ⁺₄ σ⁻₇ - 2.0 σ⁻₄ σ⁺₅ - 2.0 σ⁻₄ σ⁺₇ - σᶻ₅ σᶻ₈ - 2.0 σ⁺₅ σ⁻₈ - 2.0 σ⁻₅ σ⁺₈ - σᶻ₆ σᶻ₇ - σᶻ₆ σᶻ₈ - 2.0 σ⁺₆ σ⁻₇ - 2.0 σ⁺₆ σ⁻₈ - 2.0 σ⁻₆ σ⁺₇ - 2.0 σ⁻₆ σ⁺₈ - σᶻ₇ σᶻ₈ - 2.0 σ⁺₇ σ⁻₈ - 2.0 σ⁻₇ σ⁺₈"
#     )
#     check(basis, expr)
#
#     basis = ls.Basis.from_json(
#         '{"hamming_weight":null,"number_spins":14,"particle":"spin-1/2","spin_inversion":null,"symmetries":[{"permutation":[0,1,2,3,4,5,6,7,8,9,10,11,12,13],"sector":0},{"permutation":[1,2,3,4,5,6,7,8,9,10,11,12,13,0],"sector":11},{"permutation":[2,3,4,5,6,7,8,9,10,11,12,13,0,1],"sector":4},{"permutation":[3,4,5,6,7,8,9,10,11,12,13,0,1,2],"sector":5},{"permutation":[4,5,6,7,8,9,10,11,12,13,0,1,2,3],"sector":1},{"permutation":[5,6,7,8,9,10,11,12,13,0,1,2,3,4],"sector":13},{"permutation":[6,7,8,9,10,11,12,13,0,1,2,3,4,5],"sector":5},{"permutation":[7,8,9,10,11,12,13,0,1,2,3,4,5,6],"sector":1},{"permutation":[8,9,10,11,12,13,0,1,2,3,4,5,6,7],"sector":2},{"permutation":[9,10,11,12,13,0,1,2,3,4,5,6,7,8],"sector":1},{"permutation":[10,11,12,13,0,1,2,3,4,5,6,7,8,9],"sector":6},{"permutation":[11,12,13,0,1,2,3,4,5,6,7,8,9,10],"sector":9},{"permutation":[12,13,0,1,2,3,4,5,6,7,8,9,10,11],"sector":3},{"permutation":[13,0,1,2,3,4,5,6,7,8,9,10,11,12],"sector":3}]}'
#     )
#     expr = ls.Expr(
#         "-σᶻ₀ σᶻ₁ - σᶻ₀ σᶻ₁₃ - 2.0 σ⁺₀ σ⁻₁ - 2.0 σ⁺₀ σ⁻₁₃ - 2.0 σ⁻₀ σ⁺₁ - 2.0 σ⁻₀ σ⁺₁₃ - σᶻ₁ σᶻ₂ - 2.0 σ⁺₁ σ⁻₂ - 2.0 σ⁻₁ σ⁺₂ - σᶻ₂ σᶻ₃ - 2.0 σ⁺₂ σ⁻₃ - 2.0 σ⁻₂ σ⁺₃ - σᶻ₃ σᶻ₄ - 2.0 σ⁺₃ σ⁻₄ - 2.0 σ⁻₃ σ⁺₄ - σᶻ₄ σᶻ₅ - 2.0 σ⁺₄ σ⁻₅ - 2.0 σ⁻₄ σ⁺₅ - σᶻ₅ σᶻ₆ - 2.0 σ⁺₅ σ⁻₆ - 2.0 σ⁻₅ σ⁺₆ - σᶻ₆ σᶻ₇ - 2.0 σ⁺₆ σ⁻₇ - 2.0 σ⁻₆ σ⁺₇ - σᶻ₇ σᶻ₈ - 2.0 σ⁺₇ σ⁻₈ - 2.0 σ⁻₇ σ⁺₈ - σᶻ₈ σᶻ₉ - 2.0 σ⁺₈ σ⁻₉ - 2.0 σ⁻₈ σ⁺₉ - σᶻ₉ σᶻ₁₀ - 2.0 σ⁺₉ σ⁻₁₀ - 2.0 σ⁻₉ σ⁺₁₀ - σᶻ₁₀ σᶻ₁₁ - 2.0 σ⁺₁₀ σ⁻₁₁ - 2.0 σ⁻₁₀ σ⁺₁₁ - σᶻ₁₁ σᶻ₁₂ - 2.0 σ⁺₁₁ σ⁻₁₂ - 2.0 σ⁻₁₁ σ⁺₁₂ - σᶻ₁₂ σᶻ₁₃ - 2.0 σ⁺₁₂ σ⁻₁₃ - 2.0 σ⁻₁₂ σ⁺₁₃"
#     )
#     check(basis, expr)
#
#
# def test_csr_matvec():
#     basis = ls.SpinBasis(2)
#     basis.build()
#     expr = ls.Expr("σᶻ₀ σᶻ₁ + 2 σ⁺₀ σ⁻₁ + 2 σ⁻₀ σ⁺₁")
#     hamiltonian = ls.Operator(basis, expr)
#     matrix = hamiltonian.to_csr()
#     x = np.random.rand(basis.number_states).astype(np.complex128)
#     assert np.allclose(hamiltonian @ x, ls.matrix_vector_product_csr(matrix, x))
#
#
# def test_abelian_representations():
#     basis = ls.SpinBasis(3)
#     basis.build()
#     expr = ls.Expr("σᶻ₀ σᶻ₁ + 2 σ⁺₀ σ⁻₁ + 2 σ⁻₀ σ⁺₁", sites=[[0, 1], [1, 2], [2, 0]])
#     hamiltonian = ls.Operator(basis, expr)
#     # print(hamiltonian.to_csr().todense().real)
#     for r in hamiltonian.abelian_representations():
#         print(r)
#
#
# def test_weird_segfault_1():
#     basis = ls.SpinBasis(8, symmetries=ls.Symmetries([ls.Symmetry([7, 0, 1, 2, 3, 4, 5, 6], 0)]))
#     basis.build()
#     expr = ls.Expr(
#         # "σᶻ₀ σᶻ₁ + ",
#         "2 σ⁺₀ σ⁻₁ + 2 σ⁻₀ σ⁺₁",
#         sites=[(i, (i + 1) % basis.number_bits) for i in range(basis.number_bits)],
#     )
#     hamiltonian = ls.Operator(basis, expr)
#     # print(ls.SpinBasis.from_json(basis.to_json()))
#     print("starting ...")
#     print(basis.number_bits)
#     print(basis.number_states)
#     # x = np.random.rand(basis.number_states)
#     # y = hamiltonian @ x
#     energy, _ = scipy.sparse.linalg.eigsh(hamiltonian, k=1, which="SA")
#     print(energy)
#
#
# def test_ground_state_in_abelian_representations():
#     for k in [8]:  # , 3, 5, 6, 8]:
#         basis = ls.SpinBasis(k)
#         basis.build()
#         expr = ls.Expr(
#             "σᶻ₀ σᶻ₁ + 2 σ⁺₀ σ⁻₁ + 2 σ⁻₀ σ⁺₁",
#             sites=[(i, (i + 1) % basis.number_bits) for i in range(basis.number_bits)],
#         )
#         hamiltonian = ls.Operator(basis, expr)
#         real_energy, _ = scipy.sparse.linalg.eigsh(hamiltonian, k=1, which="SA")
#
#         energies = []
#         for r in hamiltonian.abelian_representations():
#             symm_basis = ls.SpinBasis(basis.number_bits, symmetries=r)
#             symm_basis.build()
#             symm_hamiltonian = ls.Operator(symm_basis, expr)
#             # if symm_basis.number_states > 2:
#             #     print(symm_basis.number_states)
#             #     e = scipy.sparse.linalg.eigsh(symm_hamiltonian, k=1, which="SA")[0][0]
#             # else:
#             assert symm_basis.number_states > 0
#             e = scipy.linalg.eigvalsh(symm_hamiltonian.to_csr().todense())
#             energies.append(e[0])
#
#         assert np.any(np.isclose(real_energy, energies))


# test_kagome_symmetries()
# test_weird_segfault_1()
# test_randomized_matvec()
