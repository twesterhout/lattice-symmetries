# Copyright (c) 2019-2021, Tom Westerhout
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
"""Wrapper around lattice_symmetries C library providing some handy functions for constructing and
working with quantum many-body bases.

See <https://github.com/twesterhout/lattice-symmetries> for more info.
"""

__version__ = "0.8.3"
__author__ = "Tom Westerhout <14264576+twesterhout@users.noreply.github.com>"

__all__ = [
    "LatticeSymmetriesException",
    "Symmetry",
    "Group",
    "SpinBasis",
    "Interaction",
    "Operator",
    "diagonalize",
    "enable_logging",
    "disable_logging",
    "is_logging_enabled",
]

import ctypes
from ctypes import (
    CFUNCTYPE,
    POINTER,
    byref,
    c_void_p,
    c_bool,
    c_char,
    c_char_p,
    c_int,
    c_uint,
    c_uint8,
    c_uint16,
    c_uint32,
    c_uint64,
    c_double,
)
import inspect
import math
import numpy as np
import os
import subprocess
import sys
import time
from typing import List, Optional, Tuple, Union
import warnings
import weakref


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
    """Get current package installation path."""
    return os.path.dirname(os.path.realpath(__file__))


def __load_shared_library():
    """Load lattice_symmetries C library."""
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
        # Debug logging
        ("ls_enable_logging", [], None),
        ("ls_disable_logging", [], None),
        ("ls_is_logging_enabled", [], c_bool),
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
        ("ls_symmetry_get_network_depth", [c_void_p], c_uint),
        ("ls_symmetry_get_network_masks", [c_void_p, c_void_p, c_uint64], c_int),
        ("ls_symmetry_get_permutation", [c_void_p, POINTER(c_uint32)], None),
        ("ls_batched_apply_symmetry", [c_void_p, c_uint64, POINTER(c_uint64), c_uint64], None),
        ("ls_symmetry_sizeof", [], c_uint64),
        # Group
        ("ls_create_group", [POINTER(c_void_p), c_uint, POINTER(c_void_p)], c_int),
        ("ls_destroy_group", [c_void_p], None),
        ("ls_get_group_size", [c_void_p], c_uint),
        ("ls_group_get_number_spins", [c_void_p], c_int),
        ("ls_group_get_network_depth", [c_void_p], c_int),
        ("ls_group_dump_symmetry_info", [c_void_p, c_void_p, POINTER(c_double)], c_int),
        ("ls_group_get_symmetries", [c_void_p], c_void_p),
        # Basis
        ("ls_create_spin_basis", [POINTER(c_void_p), c_void_p, c_uint, c_int, c_int], c_int),
        ("ls_destroy_spin_basis", [c_void_p], None),
        ("ls_get_number_spins", [c_void_p], c_uint),
        ("ls_get_number_bits", [c_void_p], c_uint),
        ("ls_get_hamming_weight", [c_void_p], c_int),
        ("ls_has_symmetries", [c_void_p], c_bool),
        ("ls_get_number_states", [c_void_p, POINTER(c_uint64)], c_int),
        ("ls_build", [c_void_p], c_int),
        ("ls_build_unsafe", [c_void_p, c_uint64, POINTER(c_uint64)], c_int),
        # ("ls_get_state_info", [c_void_p, POINTER(ls_bits512), POINTER(ls_bits512), c_double * 2, POINTER(c_double)], None),
        ("ls_get_state_info", [c_void_p, POINTER(c_uint64), POINTER(c_uint64), c_void_p, POINTER(c_double)], None),
        ("ls_batched_get_state_info", [c_void_p, c_uint64, POINTER(c_uint64), c_uint64,
                                       POINTER(c_uint64), c_uint64,
                                       c_void_p, c_uint64,
                                       POINTER(c_double), c_uint64], None),
        ("ls_get_index", [c_void_p, c_uint64, POINTER(c_uint64)], c_int),
        ("ls_batched_get_index", [c_void_p, c_uint64, POINTER(c_uint64), c_uint64, POINTER(c_uint64), c_uint64], c_int),
        ("ls_get_states", [POINTER(c_void_p), c_void_p], c_int),
        ("ls_destroy_states", [c_void_p], None),
        ("ls_states_get_data", [c_void_p], POINTER(c_uint64)),
        ("ls_states_get_size", [c_void_p], c_uint64),
        ("ls_save_cache", [c_void_p, c_char_p], c_int),
        ("ls_load_cache", [c_void_p, c_char_p], c_int),
        # Flat basis
        ("ls_convert_to_flat_spin_basis", [POINTER(c_void_p), c_void_p], c_int),
        ("ls_destroy_flat_spin_basis", [c_void_p], None),
        ("ls_get_buffer_size_for_flat_spin_basis", [c_void_p], c_uint64),
        ("ls_serialize_flat_spin_basis", [c_void_p, POINTER(c_char), c_uint64], c_int),
        ("ls_deserialize_flat_spin_basis", [POINTER(c_void_p), POINTER(c_char), c_uint64], c_int),
        ("ls_flat_spin_basis_number_spins", [c_void_p], c_uint),
        ("ls_flat_spin_basis_hamming_weight", [c_void_p], c_int),
        ("ls_flat_spin_basis_spin_inversion", [c_void_p], c_int),
        ("ls_flat_spin_basis_state_info", [c_void_p, c_uint64, c_void_p,
                                           c_void_p, POINTER(c_double), POINTER(c_double)], None),
        ("ls_flat_spin_basis_is_representative", [c_void_p, c_uint64, c_void_p,
                                                  POINTER(c_uint8), POINTER(c_double)], None),
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
        ("ls_batched_operator_apply", [c_void_p, c_uint64, POINTER(c_uint64),
                                       POINTER(c_uint64), c_void_p, POINTER(c_uint64)], c_uint64),
        ("ls_operator_matmat", [c_void_p, c_int, c_uint64, c_uint64, c_void_p, c_uint64, c_void_p, c_uint64], c_int),
        ("ls_operator_expectation", [c_void_p, c_int, c_uint64, c_uint64, c_void_p, c_uint64, c_void_p], c_int),
    ]
    # fmt: on
    for (name, argtypes, restype) in info:
        f = getattr(_lib, name)
        f.argtypes = argtypes
        f.restype = restype


__preprocess_library()


def enable_logging() -> None:
    """Turn on debug logging in lattice_symmetries C library."""
    _lib.ls_enable_logging()


def disable_logging() -> None:
    """Turn off debug logging in lattice_symmetries C library."""
    _lib.ls_disable_logging()


def is_logging_enabled() -> bool:
    """Return whether debug logging is currently enabled."""
    return _lib.ls_is_logging_enabled()


def debug_log(msg: str, end: str = "\n") -> None:
    if is_logging_enabled():
        current_frame = inspect.currentframe()
        parent_frame = inspect.getouterframes(current_frame)[1]
        try:
            filename = parent_frame.filename
            line = parent_frame.lineno
            function = parent_frame.function
        finally:
            del parent_frame
            del current_frame

        if len(filename) > 40:
            filename = "..." + filename[-37:]
        current_time = time.time()
        millisec = int(round(1000 * (current_time - int(current_time))))
        time_str = time.strftime("%H:%M:%S", time.localtime(int(current_time)))
        sys.stderr.write(
            "\x1b[1m\x1b[97m[Debug]\x1b[0m [{}.{:03d}] [{}:{}:{}] {}{}".format(
                time_str, millisec, filename, line, function, msg, end
            )
        )


def _get_error_message(status: int) -> str:
    """Convert `ls_error_code` produced by lattice_symmetries C library into a
    human-readable string.
    """
    raw = _lib.ls_error_to_string(status)
    msg = ctypes.string_at(raw).decode()
    _lib.ls_destroy_string(raw)
    return msg


class LatticeSymmetriesException(Exception):
    """Exception type which is used to report errors from lattice_symmetries C library."""

    def __init__(self, error_code: int):
        """Constructs the exception. `error_code` is the status code obtained from the C library."""
        self.status = error_code
        self.message = _get_error_message(error_code)
        super().__init__(self.message + " (error code: {})".format(self.status))


def _check_error(status: int) -> None:
    """Check `status` and raise a `LatticeSymmetriesException` in case of an error."""
    if status != 0:
        raise LatticeSymmetriesException(status)


def _get_dtype(dtype: np.dtype) -> int:
    """Convert NumPy datatype to `ls_datatype` enum"""
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


def _create_symmetry(permutation: List[int], sector: int) -> c_void_p:
    assert isinstance(sector, int)
    permutation = np.asarray(permutation, dtype=np.uint32)
    symmetry = c_void_p()
    _check_error(
        _lib.ls_create_symmetry(
            byref(symmetry),
            permutation.size,
            permutation.ctypes.data_as(POINTER(c_uint)),
            sector,
        )
    )
    return symmetry


def _destroy(fn):
    known_destructors = [
        (_lib.ls_destroy_symmetry, "Symmetry"),
        (_lib.ls_destroy_group, "Group"),
        (_lib.ls_destroy_spin_basis, "SpinBasis"),
        (_lib.ls_destroy_flat_spin_basis, "FlatSpinBasis"),
        (_lib.ls_destroy_states, "states array"),
        (_lib.ls_destroy_interaction, "Interaction"),
        (_lib.ls_destroy_operator, "Operator"),
        (_lib.ls_destroy_string, "C-string"),
    ]
    name = None
    for (k, v) in known_destructors:
        if k == fn:
            name = v
            break
    if name is None:
        raise ValueError("Unknown destructor: {}".format(fn))

    def wrapper(*args, **kwargs):
        debug_log("Destroying {} on Python side...".format(name))
        return fn(*args, **kwargs)

    return wrapper


class Symmetry:
    """Symmetry operator (wrapper around `ls_symmetry` C type).

    >>> # Lattice momentum with eigenvalue -ⅈ for a chain of 4 spins.
    >>> p = lattice_symmetries.Symmetry([1, 2, 3, 0], sector=1)
    >>> p.sector
    1
    >>> p.periodicity
    4
    >>> p.eigenvalue
    -1j
    """

    def __init__(self, permutation: List[int], sector: int):
        """Create a symmetry given a `permutation` of sites and `sector` specifying the eigenvalue."""
        self._payload = _create_symmetry(permutation, sector)
        self._finalizer = weakref.finalize(self, _destroy(_lib.ls_destroy_symmetry), self._payload)

    @staticmethod
    def _view_pointer(p: c_void_p, parent=None):
        s = Symmetry([], 0)
        s._payload = p
        s._finalizer = None
        return s

    @property
    def sector(self) -> int:
        """Symmetry sector."""
        return _lib.ls_get_sector(self._payload)

    @property
    def phase(self) -> float:
        """Phase of the eigenvalue."""
        return _lib.ls_get_phase(self._payload)

    @property
    def eigenvalue(self) -> complex:
        """Symmetry eigenvalue."""
        out = (c_double * 2)()
        _lib.ls_get_eigenvalue(self._payload, out)
        return complex(out[0], out[1])

    @property
    def periodicity(self) -> int:
        """Periodicity of the symmetry operator."""
        return _lib.ls_get_periodicity(self._payload)

    @property
    def number_spins(self) -> int:
        """Number of spins on which the symmetry operator acts."""
        return _lib.ls_symmetry_get_number_spins(self._payload)

    @property
    def network_depth(self) -> int:
        """Depth of the underlying Benes network."""
        return _lib.ls_symmetry_get_network_depth(self._payload)

    @property
    def network_masks(self) -> np.ndarray:
        """Masks of the underlying Benes network."""
        width = 8 if self.number_spins > 64 else 1
        masks = np.empty((self.network_depth, width), dtype=np.uint64)
        _lib.ls_symmetry_get_network_masks(
            self._payload,
            masks.ctypes.data_as(c_void_p),
            1,
        )
        return masks

    @property
    def permutation(self) -> np.ndarray:
        """Underlying permutation."""
        out = np.empty((self.number_spins,), dtype=np.uint32)
        _lib.ls_symmetry_get_permutation(
            self._payload,
            out.ctypes.data_as(POINTER(c_uint32)),
        )
        return out

    @staticmethod
    def load_from_yaml(src):
        """Load Symmetry from a parsed YAML document."""
        return Symmetry(src["permutation"], src["sector"])

    def __call__(self, spins: np.ndarray) -> None:
        if not isinstance(spins, np.ndarray) or spins.dtype != np.uint64 or spins.ndim != 2:
            raise TypeError("'spins' must be a 2D NumPy array of uint64")
        if (self.number_spins + 63) // 64 != spins.shape[1]:
            raise ValueError(
                "expected 'spins' to have {} columns, but it has {}"
                "".format((self.number_spins + 63) // 64, spins.shape[1])
            )
        if not spins.flags["C_CONTIGUOUS"]:
            spins = np.ascontiguousarray(spins)

        batch_size = spins.shape[0]
        _lib.ls_batched_apply_symmetry(
            self._payload,
            spins.shape[0],
            spins.ctypes.data_as(POINTER(c_uint64)),
            spins.strides[0] // spins.itemsize,
        )


def _create_group(generators: List[Symmetry]) -> c_void_p:
    # Things will break really badly if an element of the generators list
    # happens to be a Group or SpinBasis. They also have _payload attribute
    # which will also return a c_void_p, but C code will not be happy... :/
    if not all(map(lambda x: isinstance(x, Symmetry), generators)):
        raise TypeError("'generators' must be a List[Symmetry]")
    view = (c_void_p * len(generators))()
    for i in range(len(generators)):
        view[i] = generators[i]._payload
    group = c_void_p()
    _check_error(_lib.ls_create_group(byref(group), len(generators), view))
    return group


class Group:
    """Symmetry group (wrapper around `ls_group` C type).

    >>> T = lattice_symmetries.Symmetry([1, 2, 3, 0], sector=0) # translation
    >>> P = lattice_symmetries.Symmetry([3, 2, 1, 0], sector=0) # parity
    >>> group = lattice_symmetries.Group([T, P])
    >>> len(group)
    8
    """

    def __init__(self, generators: List[Symmetry]):
        """Construct a symmetry group from a list of generators."""
        self._payload = _create_group(generators)
        self._finalizer = weakref.finalize(self, _destroy(_lib.ls_destroy_group), self._payload)

    def __len__(self):
        return _lib.ls_get_group_size(self._payload)

    @property
    def network_depth(self):
        depth = _lib.ls_group_get_network_depth(self._payload)
        if depth < 0:
            return None
        return depth

    @property
    def number_spins(self):
        n = _lib.ls_group_get_number_spins(self._payload)
        if n < 0:
            return None
        return n

    def dump_symmetry_info(self):
        if len(self) == 0:
            raise ValueError("expected a non-empty group")
        depth = self.network_depth
        number_masks = len(self)
        mask_size = 8 if self.number_spins > 64 else 1
        masks = np.empty((depth, number_masks, mask_size), dtype=np.uint64)
        eigenvalues = np.empty((number_masks,), dtype=np.complex128)
        _check_error(
            _lib.ls_group_dump_symmetry_info(
                self._payload,
                masks.ctypes.data_as(c_void_p),
                eigenvalues.ctypes.data_as(POINTER(c_double)),
            )
        )
        return masks, eigenvalues

    @property
    def symmetries(self):
        """Symmetries of this group."""
        symmetries = []
        n = len(self)
        p = _lib.ls_group_get_symmetries(self._payload)
        for i in range(n):
            s = Symmetry._view_pointer(p + i * _lib.ls_symmetry_sizeof())
            symmetries.append(Symmetry(s.permutation, s.sector))
        return symmetries


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
            byref(basis), group._payload, number_spins, hamming_weight, spin_inversion
        )
    )
    return basis


def _int_to_ls_bits512(x: int) -> ls_bits512:
    x = int(x)
    bits = ls_bits512()
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
    """Hilbert space basis for a spin system (wrapper around `ls_spin_basis` C type)."""

    def __init__(
        self,
        group: Group,
        number_spins: int,
        hamming_weight: Optional[int] = None,
        spin_inversion: Optional[int] = None,
    ):
        """Construct a spin basis given a symmetry group, number of spins in the system,
        (optionally) the Hamming weight to which to restrict the Hilbert space, and (optionally) the
        phase the system acquires upon global spin inversion.
        """
        self._payload = _create_spin_basis(group, number_spins, hamming_weight, spin_inversion)
        self._finalizer = weakref.finalize(
            self, _destroy(_lib.ls_destroy_spin_basis), self._payload
        )

    @property
    def number_spins(self) -> int:
        """Number of spins in the system."""
        return _lib.ls_get_number_spins(self._payload)

    @property
    def number_bits(self) -> int:
        """Number of bits used to represent the spin configuration."""
        return _lib.ls_get_number_bits(self._payload)

    @property
    def hamming_weight(self) -> Optional[int]:
        """Hamming weight of all spin configurations, `None` if it varies."""
        r = _lib.ls_get_hamming_weight(self._payload)
        return None if r == -1 else r

    @property
    def has_symmetries(self) -> bool:
        """Whether lattice symmetries were used to construct the basis."""
        return _lib.ls_has_symmetries(self._payload)

    @property
    def number_states(self) -> int:
        """Number of states in the basis (i.e. dimension of the Hilbert space). This attribute is
        available only after a call to `build`."""
        r = c_uint64()
        _check_error(_lib.ls_get_number_states(self._payload, byref(r)))
        return r.value

    def build(self, representatives: Optional[np.ndarray] = None) -> None:
        """Build internal cache."""
        if representatives is None:
            _check_error(_lib.ls_build(self._payload))
        else:
            if not isinstance(representatives, np.ndarray) or representatives.dtype != np.uint64:
                raise TypeError(
                    "representatives must be a 1D NumPy array of uint64, but got {}"
                    "".format(type(representatives))
                )
            if not representatives.flags["C_CONTIGUOUS"]:
                warnings.warn(
                    "SpinBasis.build expects 'representatives' to be C-contiguous. A copy of "
                    "'representatives' will be created with proper memory order, but note that "
                    "this will uncur memory (!) overhead..."
                )
                representatives = np.ascontiguousarray(representatives)
            _check_error(
                _lib.ls_build_unsafe(
                    self._payload,
                    len(representatives),
                    representatives.ctypes.data_as(POINTER(c_uint64)),
                )
            )

    def state_info(self, bits: Union[int, np.ndarray]) -> Tuple[int, complex, float]:
        """For a spin configuration `bits` obtain its representative, corresponding
        group character, and orbit norm.
        """
        if isinstance(bits, np.ndarray):
            if bits.dtype != np.uint64 or bits.shape != (8,):
                raise TypeError(
                    "'bits' must be an 8-element 1D NumPy array of uint64, but got {}; did you mean"
                    "to call batched_state_info instead?".format(bits)
                )
            spin = ls_bits512()
            spin[:] = bits
        else:
            spin = _int_to_ls_bits512(bits)
        representative = ls_bits512()
        character = (c_double * 2)()
        norm = c_double()
        _lib.ls_get_state_info(
            self._payload,
            ctypes.cast(byref(spin), POINTER(c_uint64)),
            ctypes.cast(byref(representative), POINTER(c_uint64)),
            character,
            byref(norm),
        )
        return _ls_bits512_to_int(representative), complex(character[0], character[1]), norm.value

    def batched_state_info(self, spins: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Batched version of `self.state_info`. `batched_state_info` is equivalent to looping over `spins` and calling `self.state_info` for each element, but is much faster.
        """
        if (
            not isinstance(spins, np.ndarray)
            or spins.dtype != np.uint64
            or spins.ndim not in [1, 2]
            or (spins.ndim == 2 and spins.shape[1] != 8)
        ):
            raise TypeError(
                "'spins' must be a either 2D NumPy array of uint64 of shape (batch_size, 8) or 1D NumpyArray of uint64 of shape (batch_size,)"
            )
        one_column = spins.ndim == 1
        if one_column:
            spins = np.pad(
                spins.reshape(-1, 1),
                [(0, 0), (0, 7)],
                "constant",
                constant_values=0,
            )
        if not spins.flags["C_CONTIGUOUS"]:
            spins = np.ascontiguousarray(spins)
        batch_size = spins.shape[0]
        representative = np.zeros((batch_size, 8), dtype=np.uint64)
        eigenvalue = np.empty((batch_size,), dtype=np.complex128)
        norm = np.empty((batch_size,), dtype=np.float64)
        _lib.ls_batched_get_state_info(
            self._payload,
            spins.shape[0],
            spins.ctypes.data_as(POINTER(c_uint64)),
            spins.strides[0] // (8 * spins.itemsize),
            representative.ctypes.data_as(POINTER(c_uint64)),
            representative.strides[0] // (8 * representative.itemsize),
            eigenvalue.ctypes.data_as(POINTER(c_double)),
            eigenvalue.strides[0] // eigenvalue.itemsize,
            norm.ctypes.data_as(POINTER(c_double)),
            norm.strides[0] // norm.itemsize,
        )

        if one_column:
            representative = representative[:, 0]

        return representative, eigenvalue, norm

    def index(self, bits: int) -> int:
        """Obtain index of a representative in `self.states` array. This function is available only
        after a call to `self.build`."""
        bits = int(bits)
        i = c_uint64()
        _check_error(_lib.ls_get_index(self._payload, bits, byref(i)))
        return i.value

    def batched_index(self, spins: np.ndarray) -> np.ndarray:
        """Batched version of `self.index`. `batched_index` is equivalent to looping over `spins`
        and calling `self.index` for each element, but is much faster.
        """
        if not isinstance(spins, np.ndarray) or spins.dtype != np.uint64 or spins.ndim != 1:
            raise TypeError("'spins' must be a 1D NumPy array of uint64")
        out = np.empty(spins.shape, dtype=np.uint64)
        _check_error(
            _lib.ls_batched_get_index(
                self._payload,
                spins.shape[0],
                spins.ctypes.data_as(POINTER(c_uint64)),
                spins.strides[0] // spins.itemsize,
                out.ctypes.data_as(POINTER(c_uint64)),
                out.strides[0] // out.itemsize,
            )
        )
        return out

    @property
    def states(self) -> np.ndarray:
        """Array of representatives. This attribute is available only after a call to `self.build`."""
        states = c_void_p()
        _check_error(_lib.ls_get_states(byref(states), self._payload))
        Array = c_uint64 * _lib.ls_states_get_size(states)
        array = Array.from_address(ctypes.cast(_lib.ls_states_get_data(states), c_void_p).value)
        weakref.finalize(array, _lib.ls_destroy_states, states)
        return np.frombuffer(array, dtype=np.uint64)

    @staticmethod
    def load_from_yaml(src):
        """Load SpinBasis from a parsed YAML document."""
        number_spins = src["number_spins"]
        hamming_weight = src.get("hamming_weight")
        spin_inversion = src.get("spin_inversion")
        group = Group(list(map(Symmetry.load_from_yaml, src["symmetries"])))
        return SpinBasis(group, number_spins, hamming_weight, spin_inversion)


def _create_flat_spin_basis(basis: SpinBasis) -> c_void_p:
    if not isinstance(basis, SpinBasis):
        raise TypeError("expected SpinBasis, but got {}".format(type(group)))
    flat_basis = c_void_p()
    _check_error(_lib.ls_convert_to_flat_spin_basis(byref(flat_basis), basis._payload))
    return flat_basis


class FlatSpinBasis:
    def __init__(
        self,
        basis: SpinBasis,
    ):
        self._payload = _create_flat_spin_basis(basis)
        self._finalizer = weakref.finalize(
            self, _destroy(_lib.ls_destroy_flat_spin_basis), self._payload
        )

    @property
    def number_spins(self) -> int:
        return _lib.ls_flat_spin_basis_number_spins(self._payload)

    @property
    def hamming_weight(self) -> Optional[int]:
        r = _lib.ls_flat_spin_basis_hamming_weight(self._payload)
        return None if r == -1 else r

    @property
    def spin_inversion(self) -> Optional[int]:
        r = _lib.ls_flat_spin_basis_spin_inversion(self._payload)
        return None if r == 0 else r

    def state_info(self, spins: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        if (
            not isinstance(spins, np.ndarray)
            or spins.dtype != np.uint64
            or spins.ndim != 1
        ):
            raise TypeError("'spins' must be a 1D NumPy array of uint64 of shape (batch_size,)")
        if not spins.flags["C_CONTIGUOUS"]:
            spins = np.ascontiguousarray(spins)
        batch_size = spins.shape[0]
        representative = np.zeros((batch_size,), dtype=np.uint64)
        eigenvalue = np.empty((batch_size,), dtype=np.complex128)
        norm = np.empty((batch_size,), dtype=np.float64)
        _lib.ls_flat_spin_basis_state_info(
            self._payload,
            spins.shape[0],
            spins.ctypes.data_as(c_void_p),
            representative.ctypes.data_as(c_void_p),
            eigenvalue.ctypes.data_as(POINTER(c_double)),
            norm.ctypes.data_as(POINTER(c_double)),
        )
        return representative, eigenvalue, norm

    def is_representative(self, spins: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        if (
            not isinstance(spins, np.ndarray)
            or spins.dtype != np.uint64
            or spins.ndim != 1
        ):
            raise TypeError("'spins' must be a 1D NumPy array of uint64 of shape (batch_size,)")
        if not spins.flags["C_CONTIGUOUS"]:
            spins = np.ascontiguousarray(spins)
        batch_size = spins.shape[0]
        is_repr = np.zeros((batch_size,), dtype=np.uint8)
        norm = np.empty((batch_size,), dtype=np.float64)
        _lib.ls_flat_spin_basis_is_representative(
            self._payload,
            spins.shape[0],
            spins.ctypes.data_as(c_void_p),
            is_repr.ctypes.data_as(POINTER(c_uint8)),
            norm.ctypes.data_as(POINTER(c_double)),
        )
        return is_repr, norm

    def serialize(self) -> np.ndarray:
        n = _lib.ls_get_buffer_size_for_flat_spin_basis(self._payload)
        buf = np.zeros((n,), dtype=np.uint8)
        ("ls_serialize_flat_spin_basis", [c_void_p, POINTER(c_char), c_uint64], c_int),
        ("ls_deserialize_flat_spin_basis", [POINTER(c_void_p), POINTER(c_char), c_uint64], c_int),
        _check_error(
            _lib.ls_serialize_flat_spin_basis(self._payload, buf.ctypes.data_as(POINTER(c_char)), n)
        )
        return buf

    @staticmethod
    def deserialize(buf: np.ndarray):
        if buf.dtype != np.uint8:
            raise TypeError("'buf' has wrong dtype: {}; expected uint8".format(buf.dtype))
        buf = np.ascontiguousarray(buf)
        payload = c_void_p()
        _check_error(
            _lib.ls_deserialize_flat_spin_basis(
                byref(payload), buf.ctypes.data_as(POINTER(c_char)), buf.size
            )
        )
        basis = FlatSpinBasis.__new__(FlatSpinBasis)
        basis._payload = payload
        basis._finalizer = weakref.finalize(
            basis, _destroy(_lib.ls_destroy_flat_spin_basis), basis._payload
        )
        return basis


# import numba
#
# _ls_get_state_info = _lib.ls_get_state_info
#
#
# def _int_to_ptr_generator(pointer_type):
#     @numba.extending.intrinsic
#     def _int_to_ptr(typingctx, src):
#         from numba import types
#
#         # Check for accepted types
#         if isinstance(src, types.Integer):
#             # Custom code generation
#             def codegen(context, builder, signature, args):
#                 [src] = args
#                 llrtype = context.get_value_type(signature.return_type)
#                 return builder.inttoptr(src, llrtype)
#
#             # Create expected type signature
#             _signature = pointer_type(types.intp)
#             return _signature, codegen
#
#     return _int_to_ptr
#
#
# _int_to_uint64_ptr = _int_to_ptr_generator(numba.types.CPointer(numba.types.uint64))
# _int_to_void_ptr = _int_to_ptr_generator(numba.types.voidptr)
# _int_to_float64_ptr = _int_to_ptr_generator(numba.types.CPointer(numba.types.float64))

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


# @numba.jit(nopython=True, nogil=True, parallel=True)
# def _batched_index_helper(basis, spins):
#     basis_ptr = _int_to_void_ptr(basis)
#     batch_size = spins.shape[0]
#     indices = np.empty((batch_size,), dtype=np.uint64)
#     stride = indices.strides[0]
#     status = 0
#     for i in numba.prange(batch_size):
#         if status == 0:
#             index_ptr = _int_to_uint64_ptr(indices.ctypes.data + i * stride)
#             local_status = _ls_get_index(basis_ptr, spins[i], index_ptr)
#             if local_status != 0:
#                 status = max(status, local_status)
#     return status, indices


def batched_index(basis: SpinBasis, spins: np.ndarray) -> np.ndarray:
    warnings.warn(
        "Freestanding `batched_index(basis, spins)` function is deprecated. "
        "Please, use `basis.batched_index(spins)` instead.",
        DeprecationWarning,
    )
    return basis.batched_index(spins)


# @numba.jit(nopython=True, nogil=True, parallel=True)
# def _batched_state_info_helper(basis, spins):
#     basis_ptr = _int_to_void_ptr(basis)
#     batch_size = spins.shape[0]
#     representative = np.zeros((batch_size, 8), dtype=np.uint64)
#     eigenvalue = np.empty((batch_size,), dtype=np.complex128)
#     norm = np.empty((batch_size,), dtype=np.float64)
#     for i in numba.prange(batch_size):
#         _ls_get_state_info(
#             basis_ptr,
#             _int_to_uint64_ptr(spins.ctypes.data + i * spins.strides[0]),
#             _int_to_uint64_ptr(representative.ctypes.data + i * representative.strides[0]),
#             _int_to_void_ptr(eigenvalue.ctypes.data + i * eigenvalue.strides[0]),
#             _int_to_float64_ptr(norm.ctypes.data + i * norm.strides[0]),
#         )
#     return representative, eigenvalue, norm


def batched_state_info(
    basis: SpinBasis, spins: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    warnings.warn(
        "Freestanding `batched_state_info(basis, spins)` function is deprecated. "
        "Please, use `basis.batched_state_info(spins)` instead.",
        DeprecationWarning,
    )
    r = basis.batched_state_info(spins)
    # For testing purposes only:
    # old = _batched_state_info_helper(basis._payload.value, spins)
    # assert all(np.all(x == y) for (x, y) in zip(r, old))
    return r


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
    sites_ptr = sites.ctypes.data_as(
        POINTER(c_uint16 * number_spins) if number_spins > 1 else POINTER(c_uint16)
    )
    _check_error(f(byref(interaction), matrix_ptr, number_sites, sites_ptr))
    return interaction


def _list_to_complex(x):
    if isinstance(x, (list, tuple)) and len(x) == 2:
        return complex(x[0], x[1])
    return x


class Interaction:
    """1-, 2-, 3-, or 4-point interaction term (wrapper around `ls_interaction` C type)."""

    def __init__(self, matrix: np.ndarray, sites):
        """Create Interaction term given a matrix which specifies the interaction and a list of
        sites on which to act.
        """
        self._payload = _create_interaction(matrix, sites)
        self._finalizer = weakref.finalize(
            self, _destroy(_lib.ls_destroy_interaction), self._payload
        )

    @staticmethod
    def load_from_yaml(src):
        """Load Interaction from a parsed YAML document."""
        matrix = []
        for row in src["matrix"]:
            matrix.append([_list_to_complex(element) for element in row])
        return Interaction(matrix, src["sites"])


def _create_operator(basis: SpinBasis, terms: List[Interaction]) -> c_void_p:
    if not isinstance(basis, SpinBasis):
        raise TypeError("expected SpinBasis, but got {}".format(type(basis)))
    if not all(map(lambda x: isinstance(x, Interaction), terms)):
        raise TypeError("expected List[Interaction]")
    view = (c_void_p * len(terms))()
    for i in range(len(terms)):
        view[i] = terms[i]._payload
    op = c_void_p()
    _check_error(_lib.ls_create_operator(byref(op), basis._payload, len(terms), view))
    return op


class Operator:
    def __init__(self, basis, terms):
        self._payload = _create_operator(basis, terms)
        self._finalizer = weakref.finalize(self, _destroy(_lib.ls_destroy_operator), self._payload)
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

    def batched_apply(self, x: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        x = np.asarray(x, dtype=np.uint64)
        if x.ndim == 1:
            x = np.hstack([x.reshape(-1, 1), np.zeros((x.shape[0], 7), dtype=np.uint64)])
        elif x.ndim == 2:
            if x.shape[1] != 8:
                raise ValueError("'x' has wrong shape: {}; expected (?, 8)".format(x.shape))
            x = np.ascontiguousarray(x)
        else:
            raise ValueError("'x' has wrong shape: {}; expected a 2D array".format(x.shape))

        max_size = x.shape[0] * self.max_buffer_size
        spins = np.empty((max_size, 8), dtype=np.uint64)
        coeffs = np.empty(max_size, dtype=np.complex128)
        counts = np.empty(x.shape[0], dtype=np.uint64)
        written = _lib.ls_batched_operator_apply(
            self._payload,
            x.shape[0],
            x.ctypes.data_as(POINTER(c_uint64)),
            spins.ctypes.data_as(POINTER(c_uint64)),
            coeffs.ctypes.data_as(c_void_p),
            counts.ctypes.data_as(POINTER(c_uint64)),
        )
        return spins[:written], coeffs[:written], counts.astype(np.int64)

    def to_csr(self):
        import scipy.sparse

        self.basis.build()
        spins, coeffs, counts = self.batched_apply(self.basis.states)
        indices = self.basis.batched_index(spins[:, 0])
        row_indices = np.empty((self.basis.number_states + 1,), dtype=np.int64)
        row_indices[0] = 0
        row_indices[1:] = np.cumsum(counts)
        col_indices = indices.astype(np.int64)
        if np.all(coeffs.imag == 0):
            coeffs = np.ascontiguousarray(coeffs.real)
        return scipy.sparse.csr_matrix(
            (coeffs, col_indices, row_indices),
            shape=(self.basis.number_states, self.basis.number_states),
        )

    @staticmethod
    def load_from_yaml(src, basis: SpinBasis):
        """Load Operator from a parsed YAML document."""
        terms = list(map(Interaction.load_from_yaml, src["terms"]))
        return Operator(basis, terms)


def diagonalize(hamiltonian: Operator, k: int = 1, dtype=None, **kwargs):
    import gc
    import scipy.sparse.linalg

    hamiltonian.basis.build()
    n = hamiltonian.basis.number_states
    if dtype is None:
        dtype = np.float64

    def matvec(x):
        gc.collect()
        return hamiltonian(x)

    op = scipy.sparse.linalg.LinearOperator(shape=(n, n), matvec=matvec, dtype=dtype)
    return scipy.sparse.linalg.eigsh(op, k=k, which="SA", **kwargs)
