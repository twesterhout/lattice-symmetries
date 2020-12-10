# Copyright (c) 2019-2020, Tom Westerhout
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# * Neither the name of the copyright holder nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import ctypes
from ctypes import *
import os
import sys
import math
import subprocess
import warnings
import weakref
from typing import List, Optional, Tuple
import numpy as np

__version__ = "0.2.2"

# Enable import warnings
warnings.filterwarnings("default", category=ImportWarning)


def __library_name() -> str:
    """Get lattice_symmetries C library file name with correct extension."""
    if sys.platform == "linux":
        extension = ".so"
    elif sys.platform == "darwin":
        extension = ".dylib"
    else:
        raise ImportError("Unsupported platform: {}".format(sys.platform))
    return "liblattice_symmetries{}".format(extension)


def __package_path() -> str:
    """Get current package installation path"""
    return os.path.dirname(os.path.realpath(__file__))


def __load_shared_library():
    """Load lattice_symmetries C library"""
    libname = __library_name()
    # First, try the current directory.
    prefix = __package_path()
    if os.path.exists(os.path.join(prefix, libname)):
        return ctypes.CDLL(os.path.join(prefix, libname))
    # Next, try using conda
    if os.path.exists(os.path.join(sys.prefix, "conda-meta")):
        prefix = os.path.join(sys.prefix, "lib")
        try:
            return ctypes.CDLL(os.path.join(prefix, libname))
        except:
            warnings.warn(
                "Using python from Conda, but '{}' library was not found in "
                "the current environment. Will try pkg-config now...".format(libname),
                ImportWarning,
            )
    # Finally, try to determine the prefix using pkg-config
    result = subprocess.run(
        ["pkg-config", "--variable=libdir", "lattice_symmetries"], capture_output=True, text=True
    )
    if result.returncode != 0:
        raise ImportError("Failed to load lattice_symmetries C library")
    prefix = result.stdout.strip()
    return ctypes.CDLL(os.path.join(prefix, __library_name()))


_lib = __load_shared_library()

ls_bits512 = c_uint64 * 8
ls_callback = CFUNCTYPE(c_int, POINTER(ls_bits512), POINTER(c_double * 2), c_void_p)


def __preprocess_library():
    # fmt: off
    info = [
        # Error messages
        ("ls_error_to_string", [c_int], POINTER(c_char)),
        ("ls_destroy_string", [POINTER(c_char)], None),
        # Symmetry
        ("ls_create_symmetry", [POINTER(c_void_p), c_uint, POINTER(c_uint), c_uint], c_int),
        ("ls_destroy_symmetry", [c_void_p], None),
        ("ls_get_sector", [c_void_p], c_uint),
        ("ls_get_phase", [c_void_p], c_double),
        ("ls_get_eigenvalue", [c_void_p, c_double * 2], None),
        ("ls_get_periodicity", [c_void_p], c_uint),
        ("ls_symmetry_get_number_spins", [c_void_p], c_uint),
        # Group
        ("ls_create_group", [POINTER(c_void_p), c_uint, POINTER(c_void_p)], c_int),
        ("ls_destroy_group", [c_void_p], None),
        ("ls_get_group_size", [c_void_p], c_uint),
        # Basis
        ("ls_create_spin_basis", [POINTER(c_void_p), c_void_p, c_uint, c_int, c_int], c_int),
        ("ls_destroy_spin_basis", [c_void_p], None),
        ("ls_get_number_spins", [c_void_p], c_uint),
        ("ls_get_number_bits", [c_void_p], c_uint),
        ("ls_get_hamming_weight", [c_void_p], c_int),
        ("ls_has_symmetries", [c_void_p], c_bool),
        ("ls_get_number_states", [c_void_p, POINTER(c_uint64)], c_int),
        ("ls_build", [c_void_p], c_int),
        # ("ls_get_state_info", [c_void_p, POINTER(ls_bits512), POINTER(ls_bits512), c_double * 2, POINTER(c_double)], None),
        ("ls_get_state_info", [c_void_p, POINTER(c_uint64), POINTER(c_uint64), c_void_p, POINTER(c_double)], None),
        ("ls_get_index", [c_void_p, c_uint64, POINTER(c_uint64)], c_int),
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
        ("ls_operator_max_buffer_size", [c_void_p], c_uint64),
        ("ls_operator_apply", [c_void_p, POINTER(ls_bits512), ls_callback, c_void_p], c_int),
        ("ls_operator_matmat", [c_void_p, c_int, c_uint64, c_uint64, c_void_p, c_uint64, c_void_p, c_uint64], c_int),
        ("ls_operator_expectation", [c_void_p, c_int, c_uint64, c_uint64, c_void_p, c_uint64, c_void_p], c_int),
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
    """Used to report errors from lattice_symmetries C library."""

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


def _get_dtype(dtype):
    """Convert NumPy datatype to ls_datatype enum"""
    if dtype == np.float32:
        return 0
    if dtype == np.float64:
        return 1
    if dtype == np.complex64:
        return 2
    if dtype == np.complex128:
        return 3
    raise ValueError(
        "unexpected datatype: {}; currently only float32, float64, complex64, and "
        "complex128 are supported".format(dtype)
    )


def _create_symmetry(permutation, sector) -> c_void_p:
    permutation = np.asarray(permutation, dtype=np.uint32)
    symmetry = c_void_p()
    _check_error(
        _lib.ls_create_symmetry(
            ctypes.byref(symmetry),
            permutation.size,
            permutation.ctypes.data_as(POINTER(c_uint)),
            sector,
        )
    )
    return symmetry


class Symmetry:
    """Symmetry operator."""

    def __init__(self, permutation: List[int], sector: int):
        """Create a symmetry given a `permutation` of sites, `flip` indicating
        whether to apply global spin inversion, and `sector` specifying the
        eigenvalue
        """
        self._payload = _create_symmetry(permutation, sector)
        self._finalizer = weakref.finalize(self, _lib.ls_destroy_symmetry, self._payload)

    @property
    def sector(self) -> int:
        """Return symmetry sector"""
        return _lib.ls_get_sector(self._payload)

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
    """Symmetry group"""

    def __init__(self, generators: List[Symmetry]):
        """Construct a symmetry group from a list of generators"""
        self._payload = _create_group(generators)
        self._finalizer = weakref.finalize(self, _lib.ls_destroy_group, self._payload)

    def __len__(self):
        return _lib.ls_get_group_size(self._payload)


def _create_spin_basis(group, number_spins, hamming_weight, spin_inversion) -> c_void_p:
    if not isinstance(group, Group):
        raise TypeError("expected Group, but got {}".format(type(group)))
    if hamming_weight is None:
        hamming_weight = -1
    if spin_inversion is None:
        spin_inversion = 0
    basis = c_void_p()
    _check_error(
        _lib.ls_create_spin_basis(
            ctypes.byref(basis), group._payload, number_spins, hamming_weight, spin_inversion
        )
    )
    return basis


def _int_to_ls_bits512(x: int) -> ls_bits512:
    x = int(x)
    bits = (c_uint64 * 8)()
    for i in range(8):
        bits[i] = x & 0xFFFFFFFFFFFFFFFF
        x >>= 64
    return bits


def _ls_bits512_to_int(bits: ls_bits512) -> int:
    x = int(bits[7])
    for i in range(6, -1, -1):
        x <<= 64
        x |= int(bits[i])
    return x


class SpinBasis:
    """Hilbert space basis for a spin system"""

    def __init__(
        self,
        group: Group,
        number_spins: int,
        hamming_weight: Optional[int],
        spin_inversion: Optional[int] = None,
    ):
        """Construct a spin basis given a symmetry group, number of spins in
        the system, and (optionally) the Hamming weight to which to restrict
        the Hilbert space.
        """
        self._payload = _create_spin_basis(group, number_spins, hamming_weight, spin_inversion)
        self._finalizer = weakref.finalize(self, _lib.ls_destroy_spin_basis, self._payload)

    @property
    def number_spins(self) -> int:
        """Number of spins in the system"""
        return _lib.ls_get_number_spins(self._payload)

    @property
    def number_bits(self) -> int:
        """Number of bits used to represent the spin configuration"""
        return _lib.ls_get_number_bits(self._payload)

    @property
    def hamming_weight(self) -> Optional[int]:
        """Hamming weight of all spin configurations, `None` if it varies"""
        r = _lib.ls_get_hamming_weight(self._payload)
        return None if r == -1 else r

    @property
    def has_symmetries(self) -> bool:
        """Whether lattice symmetries were used to construct the basis"""
        return _lib.ls_has_symmetries(self._payload)

    @property
    def number_states(self) -> int:
        """Number of states in the basis (i.e. dimension of the Hilbert space)"""
        r = c_uint64()
        _check_error(_lib.ls_get_number_states(self._payload, ctypes.byref(r)))
        return r.value

    def build(self):
        """Build internal cache"""
        _check_error(_lib.ls_build(self._payload))

    def state_info(self, bits: int) -> Tuple[int, complex, float]:
        """For a spin configuration `bits` obtain its representative, corresponding
        group character, and orbit norm
        """
        bits = _int_to_ls_bits512(bits)
        representative = ls_bits512()
        character = (c_double * 2)()
        norm = c_double()
        _lib.ls_get_state_info(
            self._payload, byref(bits), byref(representative), character, byref(norm)
        )
        return _ls_bits512_to_int(representative), complex(character[0], character[1]), norm.value

    def index(self, bits: int) -> int:
        """Obtain index of a representative in `self.states` array"""
        i = c_uint64()
        _check_error(_lib.ls_get_index(self._payload, bits, byref(i)))
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


import numba

_ls_get_index = _lib.ls_get_index
_ls_get_state_info = _lib.ls_get_state_info


def _int_to_ptr_generator(pointer_type):
    @numba.extending.intrinsic
    def _int_to_ptr(typingctx, src):
        from numba import types

        # Check for accepted types
        if isinstance(src, types.Integer):
            # Custom code generation
            def codegen(context, builder, signature, args):
                [src] = args
                llrtype = context.get_value_type(signature.return_type)
                return builder.inttoptr(src, llrtype)

            # Create expected type signature
            _signature = pointer_type(types.intp)
            return _signature, codegen

    return _int_to_ptr


_int_to_uint64_ptr = _int_to_ptr_generator(numba.types.CPointer(numba.types.uint64))
_int_to_void_ptr = _int_to_ptr_generator(numba.types.voidptr)
_int_to_float64_ptr = _int_to_ptr_generator(numba.types.CPointer(numba.types.float64))

# @numba.extending.intrinsic
# def _int_to_uint64_ptr(typingctx, src):
#     from numba import types
#
#     # check for accepted types
#     if isinstance(src, types.Integer):
#         # defines the custom code generation
#         def codegen(context, builder, signature, args):
#             [src] = args
#             llrtype = context.get_value_type(signature.return_type)
#             return builder.inttoptr(src, llrtype)
#
#         # create the expected type signature
#         _signature = types.CPointer(types.uint64)(types.intp)
#         return _signature, codegen


@numba.jit(nopython=True, nogil=True, parallel=True)
def _batched_index_helper(basis, spins):
    basis_ptr = _int_to_void_ptr(basis)
    batch_size = spins.shape[0]
    indices = np.empty((batch_size,), dtype=np.uint64)
    stride = indices.strides[0]
    status = 0
    for i in numba.prange(batch_size):
        index_ptr = _int_to_uint64_ptr(indices.ctypes.data + i * stride)
        local_status = _ls_get_index(basis_ptr, spins[i], index_ptr)
        if local_status != 0:
            status = max(status, local_status)
    return status, indices


@numba.jit(nopython=True)
def _batched_state_info_helper(basis, spins):
    basis_ptr = _int_to_void_ptr(basis)
    batch_size = spins.shape[0]
    representative = np.zeros((batch_size, 8), dtype=np.uint64)
    eigenvalue = np.empty((batch_size,), dtype=np.complex128)
    norm = np.empty((batch_size,), dtype=np.float64)
    for i in numba.prange(batch_size):
        _ls_get_state_info(
            basis_ptr,
            _int_to_uint64_ptr(spins.ctypes.data + i * spins.strides[0]),
            _int_to_uint64_ptr(representative.ctypes.data + i * representative.strides[0]),
            _int_to_void_ptr(eigenvalue.ctypes.data + i * eigenvalue.strides[0]),
            _int_to_float64_ptr(norm.ctypes.data + i * norm.strides[0]),
        )
    return representative, eigenvalue, norm


def batched_index(basis: SpinBasis, spins: np.ndarray) -> np.ndarray:
    if not isinstance(basis, SpinBasis):
        raise TypeError("basis must be a SpinBasis, but got {}".format(type(basis)))
    if not isinstance(spins, np.ndarray) or spins.dtype != np.uint64:
        raise TypeError("spins must be a 1D NumPy array of uint64")
    if spins.ndim != 1:
        raise ValueError("spins has wrong shape: {}; expected a 1D array".format(spins.shape))
    status, indices = _batched_index_helper(basis._payload.value, spins)
    _check_error(status)
    return indices


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
        if x.ndim != 1 and x.ndim != 2:
            raise ValueError(
                "'x' must either a vector or a matrix, but got a {}-dimensional array"
                "".format(x.ndim)
            )
        x_was_a_vector = False
        if x.ndim == 1:
            x_was_a_vector = True
            x = x.reshape(-1, 1)
        if not x.flags["F_CONTIGUOUS"]:
            warnings.warn(
                "Operator.__call__ works with Fortran-contiguous (i.e. column-major), "
                "but 'x' is not. A copy of 'x' will be created with proper memory order, "
                "but note that this will incur performance and memory (!) overhead..."
            )
            x = np.asfortranarray(x)
        if out is None:
            out = np.empty_like(x, order="F")
        else:
            if not out.flags["F_CONTIGUOUS"]:
                warnings.warn(
                    "Operator.__call__ works with Fortran-contiguous (i.e. column-major), "
                    "but 'out' is not. A copy of 'out' will be created with proper memory order, "
                    "but note that this will incur performance and memory (!) overhead..."
                )
                out = np.asfortranarray(out)
            if x.dtype != out.dtype:
                raise ValueError(
                    "datatypes of 'x' and 'out' do not match: {} vs {}".format(x.dtype, out.dtype)
                )
        _check_error(
            _lib.ls_operator_matmat(
                self._payload,
                _get_dtype(x.dtype),
                x.shape[0],
                x.shape[1],
                x.ctypes.data_as(c_void_p),
                x.strides[1] // x.itemsize,
                out.ctypes.data_as(c_void_p),
                out.strides[1] // out.itemsize,
            )
        )
        if x_was_a_vector:
            out = np.squeeze(out)
        return out

    def expectation(self, x):
        if x.ndim != 1 and x.ndim != 2:
            raise ValueError(
                "'x' must either a vector or a matrix, but got a {}-dimensional array"
                "".format(x.ndim)
            )
        x_was_a_vector = False
        if x.ndim == 1:
            x_was_a_vector = True
            x = x.reshape(-1, 1)
        if not x.flags["F_CONTIGUOUS"]:
            warnings.warn(
                "Operator.expectation works with Fortran-contiguous (i.e.  column-major), "
                "but 'x' is not. A copy of 'x' will be created with proper memory order, "
                "but note that this will incur performance and memory (!) overhead..."
            )
            x = np.asfortranarray(x)
        out = np.empty(x.shape[1], dtype=np.complex128)
        _check_error(
            _lib.ls_operator_expectation(
                self._payload,
                _get_dtype(x.dtype),
                x.shape[0],
                x.shape[1],
                x.ctypes.data_as(c_void_p),
                x.strides[1] // x.itemsize,
                out.ctypes.data_as(c_void_p),
            )
        )
        if x_was_a_vector:
            out = complex(out)
        return out

    @property
    def max_buffer_size(self):
        return int(_lib.ls_operator_max_buffer_size(self._payload))

    def apply(self, x: int):
        max_size = self.max_buffer_size
        spins = np.empty((max_size, 8), dtype=np.uint64)
        coeffs = np.empty((max_size, 2), dtype=np.float64)
        i = 0
        e = None

        def callback(spin, coeff, cxt):
            nonlocal i, e
            try:
                spins[i] = spin.contents
                coeffs[i] = coeff.contents
                i += 1
                return 0
            except Exception as _e:
                e = _e
                return -1

        status = _lib.ls_operator_apply(
            self._payload, byref(_int_to_ls_bits512(x)), ls_callback(callback), None
        )
        if status == -1:
            assert e is not None
            raise e
        _check_error(status)
        coeffs = coeffs.view(np.complex128).reshape(-1)
        return spins[:i], coeffs[:i]


def diagonalize(hamiltonian: Operator, k: int = 1, dtype=None, **kwargs):
    import gc
    import scipy.sparse.linalg

    hamiltonian.basis.build()
    n = hamiltonian.basis.number_states
    if dtype is None:
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
    return scipy.sparse.linalg.eigsh(op, k=k, which="SA", **kwargs)
