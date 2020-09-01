import ctypes
from ctypes import *
import os
import sys
import math
import subprocess
import weakref
from typing import List, Optional, Tuple
import numpy as np


def __load_shared_library():
    # Find library location
    result = subprocess.run(
        ["pkg-config", "--variable=libdir", "lattice_symmetries"], capture_output=True, text=True
    )
    if result.returncode != 0:
        raise ImportError("Failed to load lattice_symmetries C library")
    prefix = result.stdout.strip()
    # Determine shared library extension
    if sys.platform == "linux":
        extension = ".so"
    elif sys.platform == "darwin":
        extension = ".dylib"
    else:
        raise ImportError("Unsupported platform: {}".format(sys.platform))
    # Load the library
    lib = ctypes.CDLL(os.path.join(prefix, "liblattice_symmetries{}".format(extension)))
    return lib


_lib = __load_shared_library()


def __preprocess_library():
    # fmt: off
    info = [
        # Error messages
        ("ls_error_to_string", [c_int], POINTER(c_char)),
        ("ls_destroy_string", [POINTER(c_char)], None),
        # Symmetry
        ("ls_create_symmetry", [POINTER(c_void_p), c_uint, POINTER(c_uint), c_bool, c_uint], c_int),
        ("ls_destroy_symmetry", [c_void_p], None),
        ("ls_get_sector", [c_void_p], c_uint),
        ("ls_get_flip", [c_void_p], c_bool),
        ("ls_get_phase", [c_void_p], c_double),
        ("ls_get_eigenvalue", [c_void_p, c_double * 2], None),
        ("ls_get_periodicity", [c_void_p], c_uint),
        ("ls_symmetry_get_number_spins", [c_void_p], c_uint),
        # Group
        ("ls_create_group", [POINTER(c_void_p), c_uint, POINTER(c_void_p)], c_int),
        ("ls_destroy_group", [c_void_p], None),
        ("ls_get_group_size", [c_void_p], c_uint),
        # Basis
        ("ls_create_spin_basis", [POINTER(c_void_p), c_void_p, c_uint, c_int], c_int),
        ("ls_destroy_spin_basis", [c_void_p], None),
        ("ls_get_number_spins", [c_void_p], c_uint),
        ("ls_get_number_bits", [c_void_p], c_uint),
        ("ls_get_hamming_weight", [c_void_p], c_int),
        ("ls_has_symmetries", [c_void_p], c_bool),
        ("ls_get_number_states", [c_void_p, POINTER(c_uint64)], c_int),
        ("ls_build", [c_void_p], c_int),
        ("ls_get_state_info", [c_void_p, POINTER(c_uint64), POINTER(c_uint64), c_double * 2, POINTER(c_double)], None),
        ("ls_get_index", [c_void_p, POINTER(c_uint64), POINTER(c_uint64)], c_int),
        ("ls_get_states", [POINTER(c_void_p), c_void_p], c_int),
        ("ls_destroy_states", [c_void_p], None),
        ("ls_states_get_data", [c_void_p], POINTER(c_uint64)),
        ("ls_states_get_size", [c_void_p], c_uint64),
        ("ls_save_cache", [c_void_p, c_char_p], c_int),
        ("ls_load_cache", [c_void_p, c_char_p], c_int),
        # Interaction
        ("ls_create_interaction1", [POINTER(c_void_p), c_void_p, c_uint, POINTER(c_uint16)], c_int),
        ("ls_create_interaction2", [POINTER(c_void_p), c_void_p, c_uint, POINTER(c_uint16 * 2)], c_int),
        ("ls_create_interaction3", [POINTER(c_void_p), c_void_p, c_uint, POINTER(c_uint16 * 3)], c_int),
        ("ls_create_interaction4", [POINTER(c_void_p), c_void_p, c_uint, POINTER(c_uint16 * 4)], c_int),
        ("ls_destroy_interaction", [c_void_p], None),
        # Operator
        ("ls_create_operator", [POINTER(c_void_p), c_void_p, c_uint, POINTER(c_void_p)], c_int),
        ("ls_destroy_operator", [c_void_p], None),
        ("ls_operator_matvec_f64", [c_void_p, c_uint64, POINTER(c_double), POINTER(c_double)], c_int),
    ]
    # fmt: on
    for (name, argtypes, restype) in info:
        f = getattr(_lib, name)
        f.argtypes = argtypes
        f.restype = restype


__preprocess_library()


def _get_error_message(status: int) -> str:
    """Convert `ls_error_code` produced by lattice_symmetries C library into a
    human-readable string.
    """
    raw = _lib.ls_error_to_string(status)
    msg = ctypes.string_at(raw).decode()
    _lib.ls_destroy_string(raw)
    return msg


class LatticeSymmetriesException(Exception):
    """Used to report errors from lattice_symmetries C library.
    """

    def __init__(self, error_code: int):
        """Constructs the exception. `error_code` is the status code obtained
        from the C library.
        """
        self.status = error_code
        self.message = _get_error_message(error_code)
        super().__init__(self.message + " (error code: {})".format(self.status))


def _check_error(status: int):
    """Check `status` and raise a ``LatticeSymmetriesException`` in case of an
    error.
    """
    if status != 0:
        raise LatticeSymmetriesException(status)


def _create_symmetry(permutation, flip, sector) -> c_void_p:
    permutation = np.asarray(permutation, dtype=np.uint32)
    symmetry = c_void_p()
    _check_error(
        _lib.ls_create_symmetry(
            ctypes.byref(symmetry),
            permutation.size,
            permutation.ctypes.data_as(POINTER(c_uint)),
            flip,
            sector,
        )
    )
    return symmetry


class Symmetry:
    """Symmetry operator.

    See documentation of ls_symmetry for more information.
    """

    def __init__(self, permutation: List[int], sector: int, flip: bool = False):
        self._payload = _create_symmetry(permutation, flip, sector)
        self._finalizer = weakref.finalize(self, _lib.ls_destroy_symmetry, self._payload)

    @property
    def sector(self) -> int:
        """Return symmetry sector"""
        return _lib.ls_get_sector(self._payload)

    @property
    def flip(self) -> bool:
        """Return whether symmetry applies spin inversion"""
        return _lib.ls_get_flip(self._payload)

    @property
    def phase(self) -> float:
        """Return phase of symmetry eigenvalue"""
        return _lib.ls_get_phase(self._payload)

    @property
    def eigenvalue(self) -> complex:
        """Return symmetry eigenvalue"""
        out = (c_double * 2)()
        _lib.ls_get_eigenvalue(self._payload, out)
        return complex(out[0], out[1])

    @property
    def periodicity(self) -> int:
        """Return periodicity of the operator"""
        return _lib.ls_get_periodicity(self._payload)

    @property
    def number_spins(self) -> int:
        """Return number of spins on which the operator acts"""
        return _lib.ls_symmetry_get_number_spins(self._payload)


def _create_group(generators) -> c_void_p:
    # Things will break really badly if an element of the generators list
    # happens to be a Group or SpinBasis. They also have _payload attribute
    # which will also return a c_void_p, but C code will not be happy... :/
    if not all(map(lambda x: isinstance(x, Symmetry), generators)):
        raise TypeError("expected List[Symmetry]")
    view = (c_void_p * len(generators))()
    for i in range(len(generators)):
        view[i] = generators[i]._payload
    group = c_void_p()
    _check_error(_lib.ls_create_group(ctypes.byref(group), len(generators), view))
    return group


class Group:
    def __init__(self, generators: List[Symmetry]):
        self._payload = _create_group(generators)
        self._finalizer = weakref.finalize(self, _lib.ls_destroy_group, self._payload)

    def __len__(self):
        return _lib.ls_get_group_size(self._payload)


def _create_spin_basis(group, number_spins, hamming_weight) -> c_void_p:
    if not isinstance(group, Group):
        raise TypeError("expected Group, but got {}".format(type(group)))
    if hamming_weight is None:
        hamming_weight = -1
    basis = c_void_p()
    _check_error(
        _lib.ls_create_spin_basis(ctypes.byref(basis), group._payload, number_spins, hamming_weight)
    )
    return basis


def _int_to_bits(x: int, is_big: bool) -> ctypes.Array:
    x = int(x)
    if is_big:
        bits = (c_uint64 * 8)()
        for i in range(8):
            bits[i] = x & 0xFFFFFFFFFFFFFFFF
            x >>= 64
    else:
        bits = (c_uint64 * 1)(x)
    return bits


def _bits_to_int(bits: ctypes.Array) -> int:
    if len(bits) > 1:
        x = int(bits[7])
        for i in range(6, -1, -1):
            x <<= 64
            x |= int(bits[i])
    else:
        x = int(bits[0])
    return x


class SpinBasis:
    def __init__(self, group, number_spins, hamming_weight):
        self._payload = _create_spin_basis(group, number_spins, hamming_weight)
        self._finalizer = weakref.finalize(self, _lib.ls_destroy_spin_basis, self._payload)

    @property
    def number_spins(self) -> int:
        return _lib.ls_get_number_spins(self._payload)

    @property
    def number_bits(self) -> int:
        return _lib.ls_get_number_bits(self._payload)

    @property
    def hamming_weight(self) -> Optional[int]:
        r = _lib.ls_get_hamming_weight(self._payload)
        return None if r == -1 else r

    @property
    def has_symmetries(self) -> bool:
        return _lib.ls_has_symmetries(self._payload)

    @property
    def number_states(self) -> int:
        r = c_uint64()
        _check_error(_lib.ls_get_number_states(self._payload, ctypes.byref(r)))
        return r.value

    def build(self):
        _check_error(_lib.ls_build(self._payload))

    def state_info(self, bits: int) -> Tuple[int, complex, float]:
        is_big = self.number_bits > 64
        bits = _int_to_bits(bits, is_big)
        representative = (c_uint64 * 8)() if is_big else (c_uint64 * 1)()
        character = (c_double * 2)()
        norm = c_double()
        _lib.ls_get_state_info(self._payload, bits, representative, character, byref(norm))
        return _bits_to_int(representative), complex(character[0], character[1]), norm.value

    def index(self, bits: int) -> int:
        i = c_uint64()
        _check_error(_lib.ls_get_index(self._payload, _int_to_bits(bits, False), byref(i)))
        return i.value

    @property
    def states(self) -> np.ndarray:
        states = c_void_p()
        _check_error(_lib.ls_get_states(byref(states), self._payload))
        Array = c_uint64 * _lib.ls_states_get_size(states)
        array = Array.from_address(cast(_lib.ls_states_get_data(states), c_void_p).value)
        weakref.finalize(array, _lib.ls_destroy_states, states)
        return np.frombuffer(array, dtype=np.uint64)

    def save_cache(self, filename: str):
        _check_error(_lib.ls_save_cache(self._payload, bytes(filename, "utf-8")))

    def load_cache(self, filename: str):
        _check_error(_lib.ls_load_cache(self._payload, bytes(filename, "utf-8")))


# def _create_spin_basis(group, number_spins, hamming_weight) -> c_void_p:
#     if not isinstance(group, Group):
#         raise TypeError("expected Group, but got {}".format(type(group)))
#     if hamming_weight is None:
#         hamming_weight = -1
#     basis = c_void_p()
#     _check_error(
#         _lib.ls_create_spin_basis(ctypes.byref(basis), group._payload, number_spins, hamming_weight)
#     )
#     return basis


def _deduce_number_spins(matrix) -> int:
    if matrix.ndim != 2:
        ndim = matrix.ndim
        raise ValueError("'matrix' must be a matrix, but got a {}-dimensional array".format(ndim))
    n = matrix.shape[0]
    if matrix.shape != (n, n):
        shape = matrix.shape
        raise ValueError("'matrix' must be square, but got an array of shape {}".format(shape))

    error = ValueError("'matrix' must have shape 2ⁿ x 2ⁿ where n > 0 is the number of spins")
    if n < 2:
        raise error
    number_spins = round(math.log2(n))
    if 1 << number_spins != n:
        raise error
    if number_spins not in {1, 2, 3, 4}:
        msg = "'Interaction' currently only supports interactions between 1, 2, 3 or 4 spins"
        raise ValueError(msg)
    return number_spins


def _create_interaction(matrix, sites) -> c_void_p:
    matrix = np.asarray(matrix, dtype=np.complex128, order="C")
    number_spins = _deduce_number_spins(matrix)
    sites = np.asarray(sites, dtype=np.uint16, order="C")
    if sites.ndim == 1:
        sites = sites.reshape(-1, 1)
    if sites.ndim != 2 or sites.shape[1] != number_spins:
        raise ValueError(
            "'sites' must be a list of tuples and each tuple must have length {}"
            "".format(number_spins)
        )
    f = {
        1: _lib.ls_create_interaction1,
        2: _lib.ls_create_interaction2,
        3: _lib.ls_create_interaction3,
        4: _lib.ls_create_interaction4,
    }[number_spins]

    interaction = c_void_p()
    matrix_ptr = matrix.ctypes.data_as(c_void_p)
    number_sites = sites.shape[0]
    sites_ptr = sites.ctypes.data_as(POINTER(c_uint16 * number_spins))
    _check_error(f(byref(interaction), matrix_ptr, number_sites, sites_ptr))
    return interaction


class Interaction:
    def __init__(self, matrix, sites):
        self._payload = _create_interaction(matrix, sites)
        self._finalizer = weakref.finalize(self, _lib.ls_destroy_interaction, self._payload)


def _create_operator(basis: SpinBasis, terms: List[Interaction]) -> c_void_p:
    if not isinstance(basis, SpinBasis):
        raise TypeError("expected SpinBasis, but got {}".format(type(basis)))
    if not all(map(lambda x: isinstance(x, Interaction), terms)):
        raise TypeError("expected List[Interaction]")
    view = (c_void_p * len(terms))()
    for i in range(len(terms)):
        view[i] = terms[i]._payload
    op = c_void_p()
    _check_error(_lib.ls_create_operator(ctypes.byref(op), basis._payload, len(terms), view))
    return op


class Operator:
    def __init__(self, basis, terms):
        self._payload = _create_operator(basis, terms)
        self._finalizer = weakref.finalize(self, _lib.ls_destroy_operator, self._payload)
        self.basis = basis

    def __call__(self, x, out=None):
        if x.ndim != 1:
            raise ValueError("'x' must be a vector, but got a {}-dimensional array".format(x.ndim))
        assert x.dtype == np.float64
        if out is None:
            out = np.empty_like(x)
        _check_error(
            _lib.ls_operator_matvec_f64(
                self._payload,
                len(x),
                x.ctypes.data_as(POINTER(c_double)),
                out.ctypes.data_as(POINTER(c_double)),
            )
        )
        return out


def diagonalize(hamiltonian: Operator, k: int = 1, dtype=None, tol=0):
    import gc
    import numpy as np
    import scipy.sparse.linalg

    hamiltonian.basis.build()
    n = hamiltonian.basis.number_states
    dtype = np.float64
    # if dtype is not None:
    #     if dtype not in {np.float32, np.float64, np.complex64, np.complex128}:
    #         raise ValueError(
    #             "invalid dtype: {}; expected float32, float64, complex64 or complex128"
    #             "".format(dtype)
    #         )
    #     if not hamiltonian.is_real and dtype in {np.float32, np.float64}:
    #         raise ValueError(
    #             "invalid dtype: {}; Hamiltonian is complex -- expected either complex64 "
    #             "or complex128".format(dtype)
    #         )
    # else:
    #     dtype = np.float64 if hamiltonian.is_real else np.complex128

    def matvec(x):
        gc.collect()
        return hamiltonian(x)

    def number_lanczos_vectors():
        # free = psutil.virtual_memory().free
        # usage = np.dtype(dtype).itemsize * n
        # need = 20 * usage
        # if need > free:
        #     import warnings

        #     count = (2 * free // 3) // usage
        #     warnings.warn(
        #         "Not enough memory to store the default=20 Lanczos vectors. "
        #         "Need ~{:.1f}GB, but have only ~{:.1f}GB. Will use {} Lanczos "
        #         "vectors instead.".format(need / 1024 ** 3, free / 1024 ** 3, count)
        #     )
        #     return count
        return None

    op = scipy.sparse.linalg.LinearOperator(shape=(n, n), matvec=matvec, dtype=dtype)
    return scipy.sparse.linalg.eigsh(op, k=k, ncv=number_lanczos_vectors(), which="SA", tol=tol)
