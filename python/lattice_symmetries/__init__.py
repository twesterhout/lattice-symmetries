# Copyright (c) 2022, Tom Westerhout
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

import json
import weakref
import typing
from typing import Any, List, Dict, Optional, Tuple, Union, overload, Callable

import numpy as np
from numpy.typing import NDArray
from loguru import logger
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import LinearOperator
from scipy.special import binom

from lattice_symmetries._ls_hs import ffi, lib

try:
    import igraph as ig
except:
    ig = None

__version__ = "2.2.0"


class _RuntimeInitializer:
    def __init__(self):
        logger.trace("Initializing Haskell runtime...")
        lib.ls_hs_init()
        logger.trace("Initializing Chapel runtime...")
        lib.ls_chpl_init()

    def __del__(self):
        # NOTE: The order of these should actually be reversed, but ls_chpl_finalize calls exit(0) :/
        logger.trace("Deinitializing Haskell runtime...")
        lib.ls_hs_exit()
        logger.trace("Deinitializing Chapel runtime...")
        # sys.stdout.flush()
        # sys.stderr.flush()
        lib.ls_chpl_finalize()


_runtime_init = _RuntimeInitializer()


# class Symmetry:
#     """Lattice symmetry.
#
#     >>> # Lattice momentum with eigenvalue -ⅈ for a chain of 4 spins.
#     >>> p = lattice_symmetries.Symmetry([1, 2, 3, 0], sector=1)
#     >>> p.sector
#     1
#     >>> p.permutation
#     [1, 2, 3, 0]
#     >>> p.phase
#     0.25
#     """
#
#     permutation: NDArray[np.int32]
#     sector: int
#     phase: float
#     periodicity: int
#     character: complex
#
#     def __init__(self, permutation: ArrayLike, sector: int):
#         self.permutation = np.asarray(permutation, dtype=np.int32)
#         self.sector = int(sector)
#         for name, value in _from_json(lib.ls_hs_symmetry_more_info(_to_json(self))).items():
#             setattr(self, name, value)
#
#     def __repr__(self) -> str:
#         return "Symmetry(permutation={}, sector={})".format(self.permutation, self.sector)
#
#     def __len__(self) -> int:
#         return len(self.permutation)


# class Symmetries:
#     """Lattice symmetry group representation."""
#
#     elements: List[Symmetry]
#
#     def __init__(self, generators: List[Symmetry] = []):
#         # TODO: add a Haskell function that expects a list of generators rather
#         # than a JSON representation of Symmetries to avoid recursive calls
#         self.elements = _from_json(lib.ls_hs_symmetries_from_generators(_to_json(generators)))
#
#     def __repr__(self):
#         return "Symmetries({})".format(self.elements)
#
#     def __len__(self):
#         return len(self.elements)
#
#     def __iter__(self):
#         return iter(self.elements)
#
#     def compile(self):
#         return _from_json(lib.ls_hs_compile_symmetries(_to_json(self.elements)))


class LatticeSymmetriesEncoder(json.JSONEncoder):
    def default(self, o):
        name = type(o).__name__
        try:
            encoder = getattr(self, f"encode_{name}")
        except AttributeError:
            super().default(o)
        else:
            encoded = encoder(o)
            if isinstance(encoded, dict):
                encoded["__type__"] = name
            return encoded

    # def encode_Symmetry(self, obj: Symmetry):
    #     return {"permutation": obj.permutation.tolist(), "sector": obj.sector}

    # def encode_Symmetries(self, obj: Symmetries):
    #     return obj.elements

    def encode_complex(self, obj: complex):
        return {"real": obj.real, "imag": obj.imag}


class LatticeSymmetriesDecoder(object):
    def object_hook(self, obj):
        if isinstance(obj, dict):
            if "Left" in obj or "Right" in obj:
                return self.decode_Result(obj)
            if "__type__" in obj:
                name = obj["__type__"]
                decoder = getattr(self, f"decode_{name}", lambda x: x)
                return decoder(obj)
        return obj

    def decode_Result(self, obj):
        if "Left" in obj:
            raise ValueError(obj["Left"])
        elif "Right" in obj:
            return self.object_hook(obj["Right"])
        else:
            return obj

    # def decode_Symmetry(self, obj):
    #     return Symmetry(permutation=obj["permutation"], sector=obj["sector"])

    def decode_Complex(self, obj):
        return complex(obj["real"], obj["imag"])


def _to_json(obj) -> bytes:
    return json.dumps(obj, cls=LatticeSymmetriesEncoder).encode("utf-8")


def _from_haskell_string(c_str) -> str:
    s = ffi.string(c_str).decode("utf-8")
    lib.ls_hs_destroy_string(c_str)
    return s


def _from_json(obj) -> Any:
    if isinstance(obj, ffi.CData):
        obj = _from_haskell_string(obj)
    return json.loads(obj, object_hook=LatticeSymmetriesDecoder().object_hook)


# def _basis_state_to_array(state: int, number_words: int) -> ffi.CData:
#     state = int(state)
#     arr = ffi.new("uint64_t[]", number_words)
#     for i in range(number_words):
#         arr[i] = state & 0xFFFFFFFFFFFFFFFF
#         state >>= 64
#     return arr


# class Basis:
#     _payload: ffi.CData
#     _finalizer: weakref.finalize
#
#     @overload
#     def __init__(self, payload: ffi.CData):
#         ...
#
#     @overload
#     def __init__(self, json_string: bytes):
#         ...
#
#     def _init_from_json(self, json_string: bytes):
#         print("_init_from_json:", json_string.decode("utf-8"))
#         self._payload = lib.ls_hs_basis_from_json(json_string)
#         assert self._payload != ffi.NULL
#         self._finalizer = weakref.finalize(self, lib.ls_hs_destroy_basis, self._payload)
#
#     def _init_from_payload(self, payload: ffi.CData):
#         self._payload = payload
#         assert self._payload != ffi.NULL
#         self._finalizer = weakref.finalize(self, lib.ls_hs_destroy_basis, self._payload)
#
#     def __init__(self, *args, **kwargs):
#         if "json_string" in kwargs:
#             self._init_from_json(*args, **kwargs)
#         elif "payload" in kwargs:
#             self._init_from_payload(*args, **kwargs)
#         else:
#             raise NotImplementedError("Invalid arguments: {}, {}".format(args, kwargs))
#
#     @property
#     def is_built(self) -> bool:
#         return self._payload.representatives.elts != ffi.NULL
#
#     def check_is_built(self):
#         if not self.is_built:
#             raise ValueError(
#                 "basis states have not been built yet; "
#                 "if you wish to do so, use the basis.build() function"
#             )
#
#     def build(self) -> None:
#         """Generate a list of representatives.
#
#         These can later be accessed using the `number_states` and `states` attributes.
#         """
#         if not self.is_built:
#             lib.ls_hs_basis_build(self._payload)
#         assert self.is_built
#
#     @property
#     def number_states(self) -> int:
#         self.check_is_built()
#         return int(self._payload.representatives.num_elts)
#
#     @property
#     def states(self) -> NDArray[np.uint64]:
#         n = self.number_states
#         p = int(ffi.cast("uintptr_t", self._payload.representatives.elts))
#         # int(ffi.cast("uintptr_t", lib.ls_hs_basis_states(self._payload)))
#
#         class StatesWrapper(object):
#             def __init__(self, wrapped):
#                 self.wrapped = wrapped
#                 self.__array_interface__ = {
#                     "version": 3,
#                     "typestr": "u8",
#                     "data": (p, True),
#                     "shape": (n,),
#                 }
#
#         return np.array(StatesWrapper(self), copy=False)
#
#     @property
#     def min_state_estimate(self) -> int:
#         return lib.ls_hs_min_state_estimate(self._payload)
#
#     @property
#     def max_state_estimate(self) -> int:
#         return lib.ls_hs_max_state_estimate(self._payload)
#
#     @property
#     def number_bits(self) -> int:
#         return lib.ls_hs_basis_number_bits(self._payload)
#
#     @property
#     def number_words(self) -> int:
#         return lib.ls_hs_basis_number_words(self._payload)
#
#     @property
#     def has_spin_inversion_symmetry(self) -> bool:
#         return lib.ls_hs_basis_has_spin_inversion_symmetry(self._payload) != 0
#
#     @property
#     def has_permutation_symmetries(self) -> bool:
#         return lib.ls_hs_basis_has_permutation_symmetries(self._payload) != 0
#
#     @property
#     def is_real(self) -> bool:
#         return lib.ls_hs_basis_is_real(self._payload) != 0
#
#     @property
#     def requires_projection(self) -> bool:
#         return lib.ls_hs_basis_requires_projection(self._payload) != 0
#
#     def state_to_string(self, state: int) -> str:
#         """Pretty-print a basis state."""
#         arr = _basis_state_to_array(state, self.number_words)
#         c_str = lib.ls_hs_basis_state_to_string(self._payload, arr)
#         s = ffi.string(c_str).decode("utf-8")
#         lib.ls_hs_destroy_string(c_str)
#         return s
#
#     def state_info(self, x):
#         assert self.number_bits <= 64
#         is_scalar = isinstance(x, int)
#         if is_scalar:
#             x = np.array([x], dtype=np.uint64)
#         else:
#             x = np.asarray(x, dtype=np.uint64, order="C")
#
#         count = x.shape[0]
#         betas: NDArray[np.uint64]
#         characters: NDArray[np.complex128]
#         norms: NDArray[np.float64]
#         if self.has_permutation_symmetries:
#             betas = np.zeros_like(x)
#             characters = np.zeros(count, dtype=np.complex128)
#             norms = np.zeros(count, dtype=np.float64)
#
#             x_ptr = ffi.from_buffer("uint64_t[]", x, require_writable=False)
#             betas_ptr = ffi.from_buffer("uint64_t[]", betas, require_writable=True)
#             characters_ptr = ffi.from_buffer("ls_hs_scalar[]", characters, require_writable=True)
#             norms_ptr = ffi.from_buffer("double[]", norms, require_writable=True)
#             lib.ls_hs_state_info(
#                 self._payload, count, x_ptr, 1, betas_ptr, 1, characters_ptr, norms_ptr
#             )
#         elif self.has_spin_inversion_symmetry:
#             mask = (1 << self.number_bits) - 1
#             betas = np.bitwise_xor(x, np.uint64(mask))
#             when = betas < x
#
#             betas = np.where(when, betas, x)
#             characters = np.where(when, float(self.spin_inversion), 1.0)
#             norms = np.ones(count, dtype=np.float64)
#         else:
#             assert not self.requires_projection
#             betas = x
#             characters = np.ones(count, dtype=np.complex128)
#             norms = np.ones(count, dtype=np.float64)
#
#         if is_scalar:
#             return (int(betas[0]), complex(characters[0]), float(norms[0]))
#         else:
#             return (betas, characters, norms)
#
#     def index(self, x: int | NDArray[np.uint64]) -> int | NDArray[np.int64]:
#         """Return the index of a basis state."""
#         assert self.number_bits <= 64
#         is_scalar = False
#         x = np.asarray(x, dtype=np.uint64, order="C")
#         if x.ndim == 0:
#             is_scalar = True
#             x = np.expand_dims(x, axis=0)
#
#         count = x.shape[0]
#         indices = np.zeros(count, dtype=np.int64)
#
#         x_ptr = ffi.from_buffer("uint64_t const*", x, require_writable=False)
#         indices_ptr = ffi.from_buffer("ptrdiff_t *", indices, require_writable=True)
#         lib.ls_hs_state_index(self._payload, count, x_ptr, 1, indices_ptr, 1)
#
#         if is_scalar:
#             return int(indices[0])
#         else:
#             return indices
#
#     @staticmethod
#     def from_json(json_string: str) -> "Basis":
#         _assert_subtype(json_string, str)
#         return Basis(json_string=json_string.encode("utf-8"))
#
#     def to_json(self) -> str:
#         c_str = lib.ls_hs_basis_to_json(self._payload)
#         s = ffi.string(c_str).decode("utf-8")
#         lib.ls_hs_destroy_string(c_str)
#         return s
#
#
# class SpinBasis(Basis):
#     def __init__(
#         self,
#         number_spins: int,
#         hamming_weight: Optional[int] = None,
#         spin_inversion: Optional[int] = None,
#         symmetries: Optional[Symmetries] = None,
#     ):
#         """Create a Hilbert space basis for `number_spins` spin-1/2 particles."""
#         super().__init__(
#             json_string=_to_json(
#                 {
#                     "particle": "spin-1/2",
#                     "number_spins": number_spins,
#                     "hamming_weight": hamming_weight,
#                     "spin_inversion": spin_inversion,
#                     "symmetries": symmetries if symmetries is not None else Symmetries(),
#                 }
#             )
#         )
#
#     @property
#     def has_fixed_hamming_weight(self) -> bool:
#         return lib.ls_hs_basis_has_fixed_hamming_weight(self._payload)
#
#     @property
#     def spin_inversion(self) -> Optional[int]:
#         i = self._payload.spin_inversion
#         return int(i) if i != 0 else None
#
#
# class SpinlessFermionBasis(Basis):
#     def __init__(
#         self,
#         number_sites: int,
#         number_particles: Optional[int] = None,
#     ):
#         """Create a Hilbert space basis for spinless fermions living on a lattice with
#         `number_sites` sites. The number of fermions may be optionally specified by the
#         `number_particles` argument."""
#         super().__init__(
#             json_string=_to_json(
#                 {
#                     "particle": "spinless-fermion",
#                     "number_sites": number_sites,
#                     "number_particles": number_particles,
#                 }
#             )
#         )
#
#
# class SpinfulFermionBasis(Basis):
#     def __init__(
#         self,
#         number_sites: int,
#         number_particles: Union[None, int, Tuple[int, int]] = None,
#     ):
#         """Create a Hilbert space basis for spinful fermions living on a lattice with `number_sites`
#         sites. The total number of fermions may be optionally specified with an integer
#         `number_particles` argument. Alternatively, the number of fermions with spin up and spin
#         down can be specified separately by setting `number_particles` to a tuple
#         `(N_up, N_down)`."""
#         super().__init__(
#             json_string=_to_json(
#                 {
#                     "particle": "spinful-fermion",
#                     "number_sites": number_sites,
#                     "number_particles": number_particles,
#                 }
#             )
#         )


# def _normalize_site_indices(sites):
#     sites_arr = np.asarray(sites, dtype=np.int32, order="C")
#     if sites_arr.ndim == 0 or sites_arr.ndim > 2:
#         raise ValueError("invalid array of site indices: {}".format(sites))
#     if sites_arr.ndim == 1:
#         sites_arr = sites_arr.reshape(-1, 1)
#     return sites_arr


def _assert_subtype(variable, required_type):
    if not isinstance(variable, required_type):
        raise TypeError("expected a '{}', but got '{}'".format(required_type, type(variable)))


# def _chpl_external_array_as_ndarray(arr: ffi.CData, dtype) -> NDArray[Any]:
#     if not isinstance(dtype, np.dtype):
#         dtype = np.dtype(dtype)
#     buf = ffi.buffer(arr.elts, arr.num_elts * dtype.itemsize)
#     weakref.finalize(buf, lambda: lib.ls_hs_internal_destroy_external_array(arr))
#     return np.frombuffer(buf, dtype=dtype)


# def _to_spin_index(i) -> str:
#     if isinstance(i, int):
#         if i == 0:
#             return "↑"
#         elif i == 1:
#             return "↓"
#         else:
#             raise ValueError("invalid spin index: {}; expected either 0 or 1".format(i))
#     elif isinstance(i, str):
#         if i == "↑" or i == "↓":
#             return i
#         else:
#             raise ValueError("invalid spin index: {}; expected either ↑ or ↓".format(i))
#     else:
#         raise TypeError("invalid spin index: {}".format(i))


# def _from_spin_index(i) -> int:
#     if isinstance(i, int):
#         if i == 0:
#             return 0
#         elif i == 1:
#             return 1
#         else:
#             raise ValueError("invalid spin index: {}; expected either 0 or 1".format(i))
#     elif isinstance(i, str):
#         if i == "↑":
#             return 0
#         elif i == "↓":
#             return 1
#         else:
#             raise ValueError("invalid spin index: {}; expected either ↑ or ↓".format(i))
#     else:
#         raise TypeError("invalid spin index: {}".format(i))


# @ffi.def_extern()
# def python_replace_indices(s, i, new_s_ptr, new_i_ptr):
#     assert replace_indices_impl is not None
#     (new_s, new_i) = replace_indices_impl(s, i)
#     new_s_ptr[0] = new_s
#     new_i_ptr[0] = new_i


class ExternalArrayWrapper:
    payload: ffi.CData
    keep_alive: Optional[Any]
    finalizer: Optional[weakref.finalize]

    def __init__(
        self,
        arr: ffi.CData,
        typestr: str,
        keep_alive: Optional[Any] = None,
        finalizer=lib.ls_hs_internal_destroy_external_array,
    ):
        self.payload = arr
        self.keep_alive = keep_alive

        if finalizer is not None:
            self.finalizer = weakref.finalize(self, finalizer, self.payload)
        else:
            self.finalizer = None

        n = arr.num_elts
        p = int(ffi.cast("uintptr_t", arr.elts))
        self.__array_interface__ = {
            "version": 3,
            "typestr": typestr,
            "data": (p, True),
            "shape": (n,),
        }


class HsWrapper(object):
    _payload: ffi.CData
    _finalizer: Optional[weakref.finalize]

    def __init__(self, payload: ffi.CData, finalizer: Optional[Callable[[ffi.CData], None]]):
        assert isinstance(payload, ffi.CData)
        assert payload != 0
        self._payload = payload
        self._finalizer = (
            weakref.finalize(self, finalizer, self._payload) if finalizer is not None else None
        )


class Permutation(HsWrapper):
    _info: Optional[Dict[str, Any]]

    def __init__(
        self,
        permutation: List[int] = [],
        payload: ffi.CData | int = 0,
        finalizer: Callable[[ffi.CData], None] = lib.ls_hs_destroy_permutation,
    ):
        if len(permutation) > 0:
            permutation = [int(i) for i in permutation]
            payload = _from_json(
                lib.ls_hs_create_permutation(json.dumps(permutation).encode("utf-8"))
            )
            self.__init__(payload=payload)
        else:
            if isinstance(payload, int):
                payload = ffi.cast("ls_hs_permutation *", payload)
            super().__init__(payload=payload, finalizer=finalizer)
            self._info = None

    def _init_info(self) -> None:
        if self._info is None:
            self._info = _from_json(lib.ls_hs_permutation_info(self._payload))

    @property
    def periodicity(self) -> int:
        self._init_info()
        assert self._info is not None
        return self._info["periodicity"]

    @property
    def permutation(self) -> List[int]:
        self._init_info()
        assert self._info is not None
        return self._info["permutation"]


class Basis(HsWrapper):
    def __init__(
        self,
        json_string: str = "",
        payload: ffi.CData | int = 0,
        finalizer: Callable[[ffi.CData], None] = lib.ls_hs_destroy_basis,
    ):
        if len(json_string) > 0:
            payload = _from_json(lib.ls_hs_basis_from_json(json_string.encode("utf-8")))
        elif payload == 0:
            raise ValueError(
                "incompatible arguments: either 'payload' and 'finalizer' "
                "or 'json_string' may be specified"
            )

        if isinstance(payload, int):
            payload = ffi.cast("ls_hs_basis *", payload)

        super().__init__(payload=payload, finalizer=finalizer)
        # NOTE: Important!! Pre-initialize the basis info
        lib.ls_hs_init_basis_info(self._payload)

    def _get_info(self):
        return self._payload.info

    @property
    def hamming_weight(self) -> Optional[int]:
        h = self._get_info().hamming_weight
        return int(h) if h != -1 else None

    @property
    def spin_inversion(self) -> Optional[int]:
        i = self._get_info().spin_inversion
        return int(i) if i != 0 else None

    @property
    def min_state_estimate(self) -> int:
        return int(self._get_info().min_state_estimate)

    @property
    def max_state_estimate(self) -> int:
        return int(self._get_info().max_state_estimate)

    @property
    def number_bits(self) -> int:
        return int(self._get_info().number_bits)

    @property
    def number_words(self) -> int:
        return int(self._get_info().number_words)

    @property
    def has_permutation_symmetries(self) -> bool:
        return bool(self._get_info().has_permutation_symmetries)

    @property
    def requires_projection(self) -> bool:
        return bool(self._get_info().requires_projection)

    @property
    def is_state_index_identity(self) -> bool:
        return bool(self._get_info().is_state_index_identity)

    @property
    def is_real(self) -> bool:
        return bool(self._get_info().is_real)

    @property
    def is_built(self) -> bool:
        return self._payload.local_representatives.elts != ffi.NULL

    def check_is_built(self):
        if not self.is_built:
            raise ValueError(
                "basis states have not been built yet; "
                "if you wish to do so, use the basis.build() function"
            )

    def build(self) -> None:
        """Generate a list of representatives.

        These can later be accessed using the `number_states` and `states` attributes.
        """
        if not self.is_built:
            lib.ls_chpl_local_enumerate_states(
                self._payload, self.min_state_estimate, self.max_state_estimate
            )
        assert self.is_built

    @property
    def number_states(self) -> int:
        self.check_is_built()
        return self._payload.local_representatives.num_elts

    @property
    def states(self) -> NDArray[np.uint64]:
        self.check_is_built()
        arr = ExternalArrayWrapper(
            arr=self._payload.local_representatives,
            typestr="u8",
            keep_alive=self,
            finalizer=None,
        )
        return np.array(arr, copy=False)

    def index(self, x: int | NDArray[np.uint64]) -> int | NDArray[np.int64]:
        """Return the index of a basis state."""
        if self.number_bits > 64:
            raise ValueError(
                "it is impractical to compute indices of states with more than 64 bits"
            )

        is_scalar = False
        x = np.asarray(x, dtype=np.uint64, order="C")
        if x.ndim == 0:
            is_scalar = True
            x = np.expand_dims(x, axis=0)

        if self.is_state_index_identity:
            indices = x.astype(np.int64)
        elif not self.has_permutation_symmetries:
            if self._get_info().state_to_index_kernel == ffi.NULL:
                lib.ls_hs_init_state_to_index_kernel(self._payload)

            count = x.shape[0]
            indices = np.zeros(count, dtype=np.int64)

            pass
        else:
            raise NotImplementedError()

        # count = x.shape[0]
        # indices = np.zeros(count, dtype=np.int64)

        # x_ptr = ffi.from_buffer("uint64_t const*", x, require_writable=False)
        # indices_ptr = ffi.from_buffer("ptrdiff_t *", indices, require_writable=True)
        # lib.ls_hs_state_index(self._payload, count, x_ptr, 1, indices_ptr, 1)

        if is_scalar:
            return int(indices[0])
        else:
            return indices


#
#     def build(self) -> None:
#         """Generate a list of representatives.
#
#         These can later be accessed using the `number_states` and `states` attributes.
#         """
#         if not self.is_built:
#             lib.ls_hs_basis_build(self._payload)
#         assert self.is_built
#
#     @property
#     def number_states(self) -> int:
#         self.check_is_built()
#         return int(self._payload.representatives.num_elts)
#
#     @property
#     def states(self) -> NDArray[np.uint64]:
#         n = self.number_states
#         p = int(ffi.cast("uintptr_t", self._payload.representatives.elts))
#         # int(ffi.cast("uintptr_t", lib.ls_hs_basis_states(self._payload)))
#
#         class StatesWrapper(object):
#             def __init__(self, wrapped):
#                 self.wrapped = wrapped
#                 self.__array_interface__ = {
#                     "version": 3,
#                     "typestr": "u8",
#                     "data": (p, True),
#                     "shape": (n,),
#                 }
#
#         return np.array(StatesWrapper(self), copy=False)
#


class SpinBasis(Basis):
    def __init__(
        self,
        number_spins: int,
        hamming_weight: Optional[int] = None,
        spin_inversion: Optional[int] = None,
    ):
        """Create a Hilbert space basis for `number_spins` spin-1/2 particles."""
        super().__init__(
            json_string=json.dumps(
                {
                    "particle": "spin-1/2",
                    "number_spins": number_spins,
                    "hamming_weight": hamming_weight,
                    "spin_inversion": spin_inversion,
                }
            )
        )


class Expr(HsWrapper):
    def __init__(
        self,
        expression: str = "",
        sites: Any = None,
        particle: Optional[str] = None,
        payload: ffi.CData | int = 0,
        finalizer: Callable[[ffi.CData], None] = lib.ls_hs_destroy_expr,
    ):
        if len(expression) > 0:
            json_object: Dict[str, Any] = {"expression": expression}
            if sites is not None:
                if ig is not None and isinstance(sites, ig.Graph):
                    json_object["sites"] = [list(edge.tuple) for edge in sites.es]
                else:
                    json_object["sites"] = sites
            if particle is not None:
                json_object["particle"] = particle
            payload = _from_json(lib.ls_hs_expr_from_json(json.dumps(json_object).encode("utf-8")))
        else:
            if not (len(expression) == 0 and sites is None and particle is None):
                raise ValueError(
                    "incompatible arguments: either 'payload' and 'finalizer' "
                    "or 'expression', 'sites', and 'particle' may be specified"
                )
        if isinstance(payload, int):
            payload = ffi.cast("ls_hs_expr *", payload)
        super().__init__(payload=payload, finalizer=finalizer)

    def permutation_group(self) -> List[List[int]]:
        return _from_json(lib.ls_hs_expr_permutation_group(self._payload))

    # def to_json(self) -> str:
    #     c_str = lib.ls_hs_expr_to_json(self._payload)
    #     s = ffi.string(c_str).decode("utf-8")
    #     lib.ls_hs_destroy_string(c_str)
    #     return s

    def __str__(self) -> str:
        """Get the string representation of the underlying expression."""
        return _from_haskell_string(lib.ls_hs_expr_to_string(self._payload))

    @property
    def is_real(self) -> bool:
        return lib.ls_hs_expr_is_real(self._payload) != 0

    @property
    def is_identity(self) -> bool:
        return lib.ls_hs_expr_is_identity(self._payload) != 0

    @property
    def is_hermitian(self) -> bool:
        return lib.ls_hs_expr_is_hermitian(self._payload) != 0

    def replace_indices(self, mapping: Dict[int, int]) -> "Expr":
        return Expr(
            payload=_from_json(
                lib.ls_hs_replace_indices(
                    self._payload, _to_json([(int(k), int(v)) for k, v in mapping.items()])
                )
            )
        )

    def on(self, graph) -> "Expr":
        if ig is not None and isinstance(graph, ig.Graph):
            return Expr(expression=str(self), sites=graph)
        return Expr(expression=str(self), sites=[list(edge) for edge in graph])

    def adjoint(self) -> "Expr":
        return Expr(payload=lib.ls_hs_expr_adjoint(self._payload))

    def scale(self, coeff: complex) -> "Expr":
        coeff = complex(coeff)
        return Expr(payload=lib.ls_hs_expr_scale(coeff.real, coeff.imag, self._payload))

    def __eq__(self, other: object) -> bool:
        _assert_subtype(other, Expr)
        other = typing.cast(Expr, other)
        return lib.ls_hs_expr_equal(self._payload, other._payload) != 0

    def __add__(self, other: "Expr") -> "Expr":
        _assert_subtype(other, Expr)
        return Expr(payload=_from_json(lib.ls_hs_expr_plus(self._payload, other._payload)))

    def __sub__(self, other: "Expr") -> "Expr":
        _assert_subtype(other, Expr)
        return Expr(payload=_from_json(lib.ls_hs_expr_minus(self._payload, other._payload)))

    def __mul__(self, other: "Expr") -> "Expr":
        _assert_subtype(other, Expr)
        return Expr(payload=_from_json(lib.ls_hs_expr_times(self._payload, other._payload)))

    def __neg__(self) -> "Expr":
        return Expr(payload=lib.ls_hs_expr_negate(self._payload))

    def __rmul__(self, other: complex) -> "Expr":
        if np.isscalar(other):
            return self.scale(typing.cast(complex, other))
        else:
            return NotImplemented

    @property
    def particle_type(self) -> str:
        return _from_haskell_string(lib.ls_hs_expr_particle_type(self._payload))

    @property
    def conserves_number_particles(self) -> bool:
        return bool(lib.ls_hs_expr_conserves_number_particles(self._payload))

    @property
    def spin_inversion_invariant(self) -> bool:
        return bool(lib.ls_hs_expr_spin_inversion_invariant(self._payload))

    @property
    def number_sites(self) -> int:
        return int(lib.ls_hs_expr_number_sites(self._payload))

    def get_possible_spin_inversion(self, desired: bool | int) -> List[Optional[int]]:
        if self.spin_inversion_invariant:
            if isinstance(desired, bool):
                return [-1, 1] if bool(desired) else [None]
            else:
                return [desired]
        else:
            if isinstance(desired, bool):
                return [None]
            else:
                raise ValueError(
                    f"specified spin_inversion={desired}, but the expression "
                    "is not spin inversion-invariant"
                )

    def get_possible_hamming_weights(self, desired: bool | int) -> List[Optional[int]]:
        if self.conserves_number_particles:
            if isinstance(desired, bool):
                return list(range(0, self.number_sites + 1)) if desired else [None]
            else:
                return [int(desired)]
        else:
            if isinstance(desired, bool):
                return [None]
            else:
                raise ValueError(
                    f"specified particle_conservation={desired}, but the expression "
                    "does not conserve the number of particles"
                )

    def symmetric_basis(
        self,
        particle_conservation: bool | int = True,
        spin_inversion: bool | int = True,
        symmetries: bool = False,
    ):
        if symmetries:
            return NotImplemented

        tp = self.particle_type
        if tp == "spin-1/2":
            n = self.number_sites
            configs = []
            for h in self.get_possible_hamming_weights(desired=particle_conservation):
                if h is None or 2 * h == n:
                    for i in self.get_possible_spin_inversion(desired=spin_inversion):
                        configs.append((h, i))
                else:
                    configs.append((h, None))

            def estimate(args):
                h, i = args
                if h is not None:
                    d = binom(n, h)
                    h = -h
                else:
                    d = 1 << n
                if i is not None:
                    d /= 2
                return (d, h, i)

            configs = sorted(configs, key=estimate, reverse=True)
            for h, i in configs:
                yield SpinBasis(number_spins=n, hamming_weight=h, spin_inversion=i)

        else:
            return NotImplemented


class Operator(HsWrapper, LinearOperator):
    _basis: Basis
    _expression: Expr

    def __init__(
        self,
        expression: Optional[Expr] = None,
        basis: Optional[Basis] = None,
        payload: ffi.CData | int = 0,
        finalizer: Callable[[ffi.CData], None] = lib.ls_hs_destroy_operator,
    ):
        if expression is not None:
            _assert_subtype(expression, Expr)
            if basis is None:
                basis = next(expression.symmetric_basis())
            _assert_subtype(basis, Basis)

            payload = _from_json(lib.ls_hs_create_operator(basis._payload, expression._payload))
        else:
            if basis is not None or payload == 0:
                raise ValueError(
                    "incompatible arguments: either 'payload' and 'finalizer' "
                    "or 'expression' and 'basis' may be specified"
                )

        if isinstance(payload, int):
            payload = ffi.cast("ls_hs_operator *", payload)

        # Manually increase refcounts, because self.basis and self.expression
        # Python objects will automatically with decrease them on destruction
        lib.ls_hs_internal_object_inc_ref_count(ffi.addressof(payload.basis, "base"))  # type: ignore
        lib.ls_hs_internal_object_inc_ref_count(ffi.addressof(payload.expr, "base"))  # type: ignore
        self._basis = Basis(payload=payload.basis)  # type: ignore
        self._expression = Expr(payload=payload.expr)  # type: ignore

        super().__init__(payload=payload, finalizer=finalizer)

    @staticmethod
    def from_json(json_string: str) -> "Operator":
        return Operator(
            payload=_from_json(lib.ls_hs_operator_from_json(json_string.encode("utf-8")))
        )

    @property
    def shape(self) -> Tuple[int, int]:
        n = self.basis.number_states
        return (n, n)

    @property
    def basis(self):
        return self._basis

    @property
    def expression(self):
        return self._expression

    def apply_to_state_vector(self, vector: NDArray, out: Optional[NDArray] = None) -> NDArray:
        self.basis.check_is_built()

        if vector.ndim > 2 or vector.shape[0] != self.basis.number_states:
            raise ValueError(
                "'vector' has invalid shape: {}; expected {}"
                "".format(vector.shape, (self.basis.number_states,))
            )

        operator_is_real = self.basis.is_real and self.expression.is_real
        vector_is_real = np.isrealobj(vector)

        if vector_is_real and operator_is_real:
            vector = np.asfortranarray(vector, dtype=np.float64)
            c_type_str = "double[]"
            fn = lib.ls_chpl_matrix_vector_product_f64
        else:
            vector = np.asfortranarray(vector, dtype=np.complex128)
            c_type_str = "ls_hs_scalar[]"
            fn = lib.ls_chpl_matrix_vector_product_c128

        if out is None:
            out = np.zeros_like(vector)
        else:
            if out.shape != vector.shape:
                raise ValueError(
                    "'out' has invalid shape: {}; expected {}" "".format(out.shape, vector.shape)
                )
            if not out.flags.f_contiguous:
                raise ValueError("'out' must be Fortran contiguous")

        x_ptr = ffi.from_buffer(c_type_str, vector, require_writable=False)
        y_ptr = ffi.from_buffer(c_type_str, out, require_writable=True)
        num_vectors = out.shape[1] if out.ndim > 1 else 1
        fn(self._payload, num_vectors, x_ptr, y_ptr)
        return out

    def _matvec(self, x):
        return self.apply_to_state_vector(x)

    #
    #     # def adjoint(self):
    #     #     return Operator(self.basis, lib.ls_hs_operator_hermitian_conjugate(self._payload))
    #
    #     # @property
    #     # def is_identity(self) -> bool:
    #     #     return lib.ls_hs_operator_is_identity(self._payload)
    #
    #     # @property
    #     # def is_hermitian(self) -> bool:
    #     #     return lib.ls_hs_operator_is_hermitian(self._payload)
    #     @property
    #     def is_real(self) -> bool:
    #         return self.basis.is_real and self.expression.is_real
    #
    def __add__(self, other):
        """Add two operators."""
        if isinstance(other, Operator):
            return Operator(self.expression + other.expression, self.basis)
        else:
            return NotImplemented

    def __sub__(self, other):
        """Subtract `other` Operator from `self`."""
        if isinstance(other, Operator):
            return Operator(self.expression - other.expression, self.basis)
        else:
            return NotImplemented

    def scale(self, coeff: complex) -> "Operator":
        return Operator(self.expression.scale(coeff), self.basis)

    def __mul__(self, other):
        """Multiply `self` by `other` Operator."""
        if isinstance(other, Operator):
            return Operator(self.expression * other.expression, self.basis)
        else:
            return NotImplemented

    def __rmul__(self, other: complex):
        """Scale `self` by a scalar `other`."""
        if np.isscalar(other):
            return self.scale(typing.cast(complex, other))
        else:
            return NotImplemented

    def __matmul__(self, other: NDArray[Any]) -> NDArray[Any]:
        return self.apply_to_state_vector(other)

    def __repr__(self):
        return "<Operator defined on {}>".format(self.basis.__class__.__name__)

    def diag_to_array(self) -> NDArray[Any]:
        self.basis.check_is_built()
        diag = ffi.new("chpl_external_array *")
        lib.ls_chpl_extract_diag_c128(self._payload, diag)
        arr = ExternalArrayWrapper(arr=diag, typestr="c16")
        return np.array(arr, copy=False)

    def off_diag_to_csr(self) -> csr_matrix:
        self.basis.check_is_built()
        c_coeffs = ffi.new("chpl_external_array *")
        c_offsets = ffi.new("chpl_external_array *")
        c_indices = ffi.new("chpl_external_array *")

        lib.ls_chpl_off_diag_to_csr_c128(self._payload, c_coeffs, c_offsets, c_indices)
        coeffs = ExternalArrayWrapper(c_coeffs, "c16")
        offsets = ExternalArrayWrapper(c_offsets, "i8")
        indices = ExternalArrayWrapper(c_indices, "i8")
        return csr_matrix((coeffs, indices, offsets), shape=self.shape, dtype=np.complex128)


#
#     # def apply_off_diag_to_basis_states(
#     #     op: ls.Operator, alphas: npt.NDArray[np.uint64]
#     # ) -> tuple[npt.NDArray[np.uint64], npt.NDArray[np.complex128], npt.NDArray[np.int64]]:
#     #     alphas = np.asarray(alphas, order="C", dtype=np.uint64)
#     #     alphas_ptr = ls.ffi.from_buffer("uint64_t[]", alphas, require_writable=False)
#     #
#     #     betas = ls.ffi.new("chpl_external_array *")
#     #     coeffs = ls.ffi.new("chpl_external_array *")
#     #     offsets = ls.ffi.new("chpl_external_array *")
#     #     kernels = ls.lib.ls_hs_internal_get_chpl_kernels()
#     #     kernels.operator_apply_off_diag(
#     #         op._payload, alphas.size, alphas_ptr, betas, coeffs, offsets, 0
#     #     )
#     #
#     #     offsets_arr = ls._chpl_external_array_as_ndarray(offsets, np.int64)
#     #     betas_arr = ls._chpl_external_array_as_ndarray(betas, np.uint64)[: offsets_arr[-1]]
#     #     coeffs_arr = ls._chpl_external_array_as_ndarray(coeffs, np.complex128)[: offsets_arr[-1]]
#     #     return betas_arr, coeffs_arr, offsets_arr
#     #
#     #
#     # def apply_diag_to_basis_states(
#     #     op: ls.Operator, alphas: npt.NDArray[np.uint64]
#     # ) -> npt.NDArray[np.float64]:
#     #     alphas = np.asarray(alphas, order="C", dtype=np.uint64)
#     #     alphas_ptr = ls.ffi.from_buffer("uint64_t[]", alphas, require_writable=False)
#     #
#     #     coeffs = ls.ffi.new("chpl_external_array *")
#     #     kernels = ls.lib.ls_hs_internal_get_chpl_kernels()
#     #     kernels.operator_apply_diag(op._payload, alphas.size, alphas_ptr, coeffs, 0)
#     #
#     #     coeffs_arr = ls._chpl_external_array_as_ndarray(coeffs, np.float64)
#     #     return coeffs_arr
#
#     @overload
#     def apply_diag_to_basis_state(self, alphas: int) -> float:
#         ...
#
#     @overload
#     def apply_diag_to_basis_state(self, alphas: NDArray[np.uint64]) -> NDArray[np.float64]:
#         ...
#
#     def apply_diag_to_basis_state(
#         self, alphas: int | NDArray[np.uint64]
#     ) -> float | NDArray[np.float64]:
#         if isinstance(alphas, (int, np.integer)):
#             is_scalar = True
#             count = 1
#             alphas_ptr = _basis_state_to_array(int(alphas), self.basis.number_words)
#         else:
#             is_scalar = False
#             alphas = np.asarray(alphas, order="C", dtype=np.uint64)
#             alphas_ptr = ffi.from_buffer("uint64_t[]", alphas, require_writable=False)
#             count = alphas.size
#
#         coeffs = ffi.new("chpl_external_array *")
#         kernels = lib.ls_hs_internal_get_chpl_kernels()
#         kernels.operator_apply_diag(self._payload, count, alphas_ptr, coeffs, 0)
#
#         coeffs_arr = _chpl_external_array_as_ndarray(coeffs, np.float64)
#         if is_scalar:
#             return float(coeffs_arr[0])
#         else:
#             return coeffs_arr
#
#     @overload
#     def apply_off_diag_to_basis_state(self, alphas: int) -> List[Tuple[complex, int]]:
#         ...
#
#     @overload
#     def apply_off_diag_to_basis_state(
#         self, alphas: NDArray[np.uint64]
#     ) -> Tuple[NDArray[np.uint64], NDArray[np.complex128], NDArray[np.int64]]:
#         ...
#
#     def apply_off_diag_to_basis_state(
#         self, alphas: int | NDArray[np.uint64]
#     ) -> (
#         List[Tuple[complex, int]]
#         | Tuple[NDArray[np.uint64], NDArray[np.complex128], NDArray[np.int64]]
#     ):
#         if isinstance(alphas, (int, np.integer)):
#             is_scalar = True
#             count = 1
#             alphas_ptr = _basis_state_to_array(int(alphas), self.basis.number_words)
#         else:
#             is_scalar = False
#             alphas = np.asarray(alphas, order="C", dtype=np.uint64)
#             alphas_ptr = ffi.from_buffer("uint64_t[]", alphas, require_writable=False)
#             count = alphas.size
#
#         betas = ffi.new("chpl_external_array *")
#         coeffs = ffi.new("chpl_external_array *")
#         offsets = ffi.new("chpl_external_array *")
#         kernels = lib.ls_hs_internal_get_chpl_kernels()
#         kernels.operator_apply_off_diag(self._payload, count, alphas_ptr, betas, coeffs, offsets, 0)
#
#         offsets_arr = _chpl_external_array_as_ndarray(offsets, np.int64)
#         betas_arr = _chpl_external_array_as_ndarray(betas, np.uint64)[: offsets_arr[-1]]
#         coeffs_arr = _chpl_external_array_as_ndarray(coeffs, np.complex128)[: offsets_arr[-1]]
#         if is_scalar:
#             return list(zip(coeffs_arr, betas_arr))
#         else:
#             return betas_arr, coeffs_arr, offsets_arr
#
#     def to_csr(self):
#         row_offsets = ffi.new("chpl_external_array *")
#         col_indices = ffi.new("chpl_external_array *")
#         matrix_elements = ffi.new("chpl_external_array *")
#         kernels = lib.ls_hs_internal_get_chpl_kernels()
#         kernels.operator_to_csr(self._payload, row_offsets, col_indices, matrix_elements, 0)
#
#         offsets_arr = _chpl_external_array_as_ndarray(row_offsets, np.int32)
#         indices_arr = _chpl_external_array_as_ndarray(col_indices, np.int32)
#         coeffs_arr = _chpl_external_array_as_ndarray(matrix_elements, np.complex128)
#         size = offsets_arr[-1]
#         dim = self.basis.number_states
#
#         return scipy.sparse.csr_matrix(
#             (coeffs_arr[:size], indices_arr[:size], offsets_arr), shape=(dim, dim)
#         )
#
#     def _check_basis_is_built(self, attribute):
#         if not self.basis.is_built:
#             raise AttributeError(
#                 "'Operator' object has no attribute '{}' (did you forget to build the basis?)"
#                 "".format(attribute),
#                 name=attribute,
#                 obj=self,
#             )
#
#     @property
#     def dtype(self):
#         self._check_basis_is_built("dtype")
#         return np.dtype("float64") if self.is_real else np.dtype("complex128")
#
#     @property
#     def shape(self):
#         self._check_basis_is_built("shape")
#         n = self.basis.number_states
#         return (n, n)
#
#     def _matvec(self, x):
#         return self.apply_to_state_vector(x)
#
#     def prepare_inputs_for_hphi(self, folder: str) -> None:
#         lib.ls_hs_prepare_hphi(self._payload, folder.encode("utf-8"))
#
#     def prepare_inputs_for_mvmc(self, folder: str) -> None:
#         lib.ls_hs_prepare_mvmc(self._payload, folder.encode("utf-8"))
#
#     def abelian_representations(self):
#         # with process_symmetries_impl_lock:
#         #     global process_symmetries_impl
#         #     representations = []
#         #     process_symmetries_impl = lambda p: representations.append(Symmetries(p))
#         return _from_json(lib.ls_hs_operator_abelian_representations(self._payload))
#
#
# def load_yaml_config(filename: str):
#     config = lib.ls_hs_load_yaml_config(filename.encode("utf-8"))
#     basis = Basis(payload=lib.ls_hs_clone_basis(config.basis))
#     hamiltonian = None
#     if config.hamiltonian != 0:
#         hamiltonian = Operator(lib.ls_hs_clone_operator(config.hamiltonian))
#     if config.number_observables != 0:
#         raise NotImplementedError
#     lib.ls_hs_destroy_yaml_config(config)
#     Config = namedtuple("Config", ["basis", "hamiltonian", "observables"], defaults=[None, None])
#     return Config(basis, hamiltonian)
#
#
# def matrix_vector_product_csr(
#     matrix: scipy.sparse.csr_matrix,
#     x: NDArray[np.complex128],
#     out: None | NDArray[np.complex128] = None,
# ):
#     data = np.require(matrix.data, dtype=np.complex128, requirements=["C"])
#     indptr = np.require(matrix.indptr, dtype=np.int32, requirements=["C"])
#     indices = np.require(matrix.indices, dtype=np.int32, requirements=["C"])
#     matrix_elements = ffi.from_buffer("ls_hs_scalar[]", matrix.data, require_writable=False)
#     row_offsets = ffi.from_buffer("int32_t[]", matrix.indptr, require_writable=False)
#     col_indices = ffi.from_buffer("int32_t[]", matrix.indices, require_writable=False)
#
#     x = np.asarray(x, order="C", dtype=np.complex128)
#     x_ptr = ffi.from_buffer("ls_hs_scalar[]", x, require_writable=False)
#     if out is None:
#         out = np.empty_like(x)
#     else:
#         out = np.asarray(out, order="C", dtype=np.complex128)
#     y_ptr = ffi.from_buffer("ls_hs_scalar[]", out, require_writable=True)
#
#     kernels = lib.ls_hs_internal_get_chpl_kernels()
#     kernels.matrix_vector_product_csr_i32_c128(
#         matrix.shape[0],
#         matrix.shape[1],
#         matrix.nnz,
#         matrix_elements,
#         row_offsets,
#         col_indices,
#         x_ptr,
#         y_ptr,
#         -1,
#     )
#     return out
