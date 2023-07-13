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
import threading
import weakref
from collections import namedtuple
from functools import singledispatchmethod
from typing import Any, List, Optional, Tuple, Union, overload

import numpy as np
import scipy.sparse.linalg
from loguru import logger
from numpy.typing import ArrayLike, NDArray
from scipy.sparse.linalg import LinearOperator

import lattice_symmetries
from lattice_symmetries._ls_hs import ffi, lib

__version__ = "2.1.0"


class _RuntimeInitializer:
    def __init__(self):
        logger.debug("Initializing Haskell runtime...")
        lib.ls_hs_init()
        logger.debug("Initializing Chapel runtime...")
        lib.ls_chpl_init()
        logger.debug("Setting Python exception handler...")
        lib.set_python_exception_handler()

    def __del__(self):
        # NOTE: The order of these should actually be reversed, but ls_chpl_finalize calls exit(0) :/
        # logger.debug("Deinitializing Haskell runtime...")
        lib.ls_hs_exit()
        # logger.debug("Deinitializing Chapel runtime...")
        lib.ls_chpl_finalize()


_runtime_init = _RuntimeInitializer()


class Symmetry:
    """Symmetry operator.

    >>> # Lattice momentum with eigenvalue -ⅈ for a chain of 4 spins.
    >>> p = lattice_symmetries.Symmetry([1, 2, 3, 0], sector=1)
    >>> p.sector
    1
    >>> p.permutation
    [1, 2, 3, 0]
    """

    _payload: ffi.CData
    _finalizer: weakref.finalize

    @singledispatchmethod
    def __init__(self, permutation: ArrayLike, sector: int):
        permutation = np.asarray(permutation, dtype=int).tolist()
        sector = int(sector)
        json_object = {"permutation": permutation, "sector": sector}
        self._payload = lib.ls_hs_symmetry_from_json(json.dumps(json_object).encode("utf-8"))
        assert self._payload != 0
        self._finalizer = weakref.finalize(self, lib.ls_hs_destroy_symmetry, self._payload)

    @__init__.register
    def _(self, payload: ffi.CData, finalizer: weakref.finalize):
        self._payload = payload
        self._finalizer = finalizer

    @property
    def sector(self) -> int:
        return lib.ls_hs_symmetry_sector(self._payload)

    def __len__(self) -> int:
        return lib.ls_hs_symmetry_length(self._payload)

    @property
    def permutation(self) -> NDArray[np.int32]:
        p = lib.ls_hs_symmetry_permutation(self._payload)
        buffer = ffi.buffer(p, len(self) * ffi.sizeof("int"))
        array = np.copy(np.frombuffer(buffer, dtype=np.int32))
        lib.ls_hs_destroy_permutation(p)
        return array

    def json_object(self):
        return {"permutation": self.permutation.tolist(), "sector": self.sector}


class Symmetries:
    _payload: ffi.CData
    _finalizer: weakref.finalize
    _generators: List[Symmetry]

    def __init__(self, generators: List[Symmetry] = []):
        json_string = json.dumps([g.json_object() for g in generators]).encode("utf-8")
        self._payload = lib.ls_hs_symmetries_from_json(json_string)
        assert self._payload != 0
        self._finalizer = weakref.finalize(self, lib.ls_hs_destroy_symmetries, self._payload)
        self._generators = generators

    def __len__(self) -> int:
        return len(self.generators)

    @property
    def generators(self) -> List[Symmetry]:
        return self._generators

    def json_object(self):
        return [g.json_object() for g in self.generators]


def _basis_state_to_array(state: int, number_words: int) -> ffi.CData:
    state = int(state)
    arr = ffi.new("uint64_t[]", number_words)
    for i in range(number_words):
        arr[i] = state & 0xFFFFFFFFFFFFFFFF
        state >>= 64
    return arr


class Basis:
    _payload: ffi.CData
    _finalizer: Optional[weakref.finalize]

    @singledispatchmethod
    def __init__(self, arg, *args, **kwargs):
        raise NotImplementedError("Invalid argument type: {}".format(type(arg)))

    @__init__.register
    def _(self, payload: ffi.CData, owner: bool = True):
        self._payload = payload
        if owner:
            self._finalizer = weakref.finalize(self, lib.ls_hs_destroy_basis, self._payload)
        else:
            self._finalizer = None

    @__init__.register
    def _(self, json_string: str):
        self._payload = lib.ls_hs_basis_from_json(json_string.encode("utf-8"))
        assert self._payload != 0
        self._finalizer = weakref.finalize(self, lib.ls_hs_destroy_basis, self._payload)

    # @__init__.register
    # def _(self, json_object: dict):
    #     self.__init__(json.dumps(json_object))

    # def __init__(self, **kwargs):
    #     json_string = json.dumps(kwargs).encode("utf-8")
    #     self._payload = lib.ls_hs_basis_from_json(json_string)
    #     self._finalizer = weakref.finalize(self, lib.ls_hs_destroy_basis, self._payload)
    #     self._representatives = None

    @property
    def is_built(self) -> bool:
        return lib.ls_hs_basis_is_built(self._payload)

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
            lib.ls_hs_basis_build(self._payload)
        assert self.is_built

    @property
    def number_states(self) -> int:
        self.check_is_built()
        return lib.ls_hs_basis_number_states(self._payload)

    @property
    def states(self) -> NDArray[np.uint64]:
        n = self.number_states
        p = lib.ls_hs_basis_states(self._payload)
        buffer = ffi.buffer(p, n * ffi.sizeof("uint64_t"))
        return np.frombuffer(buffer, dtype=np.uint64)

    @property
    def min_state_estimate(self) -> int:
        return lib.ls_hs_min_state_estimate(self._payload)

    @property
    def max_state_estimate(self) -> int:
        return lib.ls_hs_max_state_estimate(self._payload)

    @property
    def number_bits(self) -> int:
        return lib.ls_hs_basis_number_bits(self._payload)

    @property
    def number_words(self) -> int:
        return lib.ls_hs_basis_number_words(self._payload)

    @property
    def has_spin_inversion_symmetry(self) -> bool:
        return lib.ls_hs_basis_has_spin_inversion_symmetry(self._payload) != 0

    @property
    def has_permutation_symmetries(self) -> bool:
        return lib.ls_hs_basis_has_permutation_symmetries(self._payload) != 0

    @property
    def requires_projection(self) -> bool:
        return lib.ls_hs_basis_requires_projection(self._payload) != 0

    def state_to_string(self, state: int) -> str:
        """Pretty-print a basis state."""
        arr = _basis_state_to_array(state, self.number_words)
        c_str = lib.ls_hs_basis_state_to_string(self._payload, arr)
        s = ffi.string(c_str).decode("utf-8")
        lib.ls_hs_destroy_string(c_str)
        return s

    def state_info(self, x):
        assert self.number_bits <= 64
        is_scalar = isinstance(x, int)
        if is_scalar:
            x = np.array([x], dtype=np.uint64)
        else:
            x = np.asarray(x, dtype=np.uint64, order="C")

        count = x.shape[0]
        betas: NDArray[np.uint64]
        characters: NDArray[np.complex128]
        norms: NDArray[np.float64]
        if self.has_permutation_symmetries:
            betas = np.zeros_like(x)
            characters = np.zeros(count, dtype=np.complex128)
            norms = np.zeros(count, dtype=np.float64)

            x_ptr = ffi.from_buffer("uint64_t[]", x, require_writable=False)
            betas_ptr = ffi.from_buffer("uint64_t[]", betas, require_writable=True)
            characters_ptr = ffi.from_buffer("ls_hs_scalar[]", characters, require_writable=True)
            norms_ptr = ffi.from_buffer("double[]", norms, require_writable=True)
            lib.ls_hs_state_info(
                self._payload, count, x_ptr, 1, betas_ptr, 1, characters_ptr, norms_ptr
            )
        elif self.has_spin_inversion_symmetry:
            mask = (1 << self.number_bits) - 1
            betas = np.bitwise_xor(x, np.uint64(mask))
            when = betas < x

            betas = np.where(when, betas, x)
            characters = np.where(when, float(self.spin_inversion), 1.0)
            norms = np.ones(count, dtype=np.float64)
        else:
            assert not self.requires_projection
            betas = x
            characters = np.ones(count, dtype=np.complex128)
            norms = np.ones(count, dtype=np.float64)

        if is_scalar:
            return (int(betas[0]), complex(characters[0]), float(norms[0]))
        else:
            return (betas, characters, norms)

    def index(self, x: Union[int, NDArray[np.uint64]]) -> Union[int, NDArray[np.int64]]:
        """Return the index of a basis state."""
        assert self.number_bits <= 64
        is_scalar = False
        x = np.asarray(x, dtype=np.uint64, order="C")
        if x.ndim == 0:
            is_scalar = True
            x = np.expand_dims(x, axis=0)

        count = x.shape[0]
        indices = np.zeros(count, dtype=np.int64)

        x_ptr = ffi.from_buffer("uint64_t const*", x, require_writable=False)
        indices_ptr = ffi.from_buffer("ptrdiff_t *", indices, require_writable=True)
        lib.ls_hs_state_index(self._payload, count, x_ptr, 1, indices_ptr, 1)

        if is_scalar:
            return int(indices[0])
        else:
            return indices

    @staticmethod
    def from_json(json_string: str) -> "Basis":
        _assert_subtype(json_string, str)
        return Basis(**json.loads(json_string))

    def to_json(self) -> str:
        c_str = lib.ls_hs_basis_to_json(self._payload)
        s = ffi.string(c_str).decode("utf-8")
        lib.ls_hs_destroy_string(c_str)
        return s


class SpinBasis(Basis):
    def __init__(
        self,
        number_spins: int,
        hamming_weight: Optional[int] = None,
        spin_inversion: Optional[int] = None,
        symmetries: Optional[Symmetries] = None,
    ):
        """Create a Hilbert space basis for `number_spins` spin-1/2 particles."""
        super().__init__(
            json.dumps(
                {
                    "particle": "spin-1/2",
                    "number_spins": number_spins,
                    "hamming_weight": hamming_weight,
                    "spin_inversion": spin_inversion,
                    "symmetries": symmetries.json_object() if symmetries is not None else [],
                }
            )
        )

    @property
    def has_fixed_hamming_weight(self) -> bool:
        return lib.ls_hs_basis_has_fixed_hamming_weight(self._payload)

    @property
    def spin_inversion(self) -> Optional[int]:
        i = self._payload.spin_inversion
        return int(i) if i != 0 else None


class SpinlessFermionBasis(Basis):
    def __init__(
        self,
        number_sites: int,
        number_particles: Optional[int] = None,
    ):
        """Create a Hilbert space basis for spinless fermions living on a lattice with
        `number_sites` sites. The number of fermions may be optionally specified by the
        `number_particles` argument."""
        super().__init__(
            json.dumps(
                {
                    "particle": "spinless-fermion",
                    "number_sites": number_sites,
                    "number_particles": number_particles,
                }
            )
        )


class SpinfulFermionBasis(Basis):
    def __init__(
        self,
        number_sites: int,
        number_particles: Union[None, int, Tuple[int, int]] = None,
    ):
        """Create a Hilbert space basis for spinful fermions living on a lattice with `number_sites`
        sites. The total number of fermions may be optionally specified with an integer
        `number_particles` argument. Alternatively, the number of fermions with spin up and spin
        down can be specified separately by setting `number_particles` to a tuple
        `(N_up, N_down)`."""
        super().__init__(
            json.dumps(
                {
                    "particle": "spinful-fermion",
                    "number_sites": number_sites,
                    "number_particles": number_particles,
                }
            )
        )


def _normalize_site_indices(sites):
    sites_arr = np.asarray(sites, dtype=np.int32, order="C")
    if sites_arr.ndim == 0 or sites_arr.ndim > 2:
        raise ValueError("invalid array of site indices: {}".format(sites))
    if sites_arr.ndim == 1:
        sites_arr = sites_arr.reshape(-1, 1)
    return sites_arr


def _assert_subtype(variable, required_type):
    if not isinstance(variable, required_type):
        raise TypeError("expected a '{}', but got '{}'".format(required_type, type(variable)))


def _chpl_external_array_as_ndarray(arr: ffi.CData, dtype) -> NDArray[Any]:
    if not isinstance(dtype, np.dtype):
        dtype = np.dtype(dtype)
    buf = ffi.buffer(arr.elts, arr.num_elts * dtype.itemsize)
    weakref.finalize(buf, lambda: lib.ls_hs_internal_destroy_external_array(arr))
    return np.frombuffer(buf, dtype=dtype)


def _to_spin_index(i) -> str:
    if isinstance(i, int):
        if i == 0:
            return "↑"
        elif i == 1:
            return "↓"
        else:
            raise ValueError("invalid spin index: {}; expected either 0 or 1".format(i))
    elif isinstance(i, str):
        if i == "↑" or i == "↓":
            return i
        else:
            raise ValueError("invalid spin index: {}; expected either ↑ or ↓".format(i))
    else:
        raise TypeError("invalid spin index: {}".format(i))


def _from_spin_index(i) -> int:
    if isinstance(i, int):
        if i == 0:
            return 0
        elif i == 1:
            return 1
        else:
            raise ValueError("invalid spin index: {}; expected either 0 or 1".format(i))
    elif isinstance(i, str):
        if i == "↑":
            return 0
        elif i == "↓":
            return 1
        else:
            raise ValueError("invalid spin index: {}; expected either ↑ or ↓".format(i))
    else:
        raise TypeError("invalid spin index: {}".format(i))


replace_indices_impl = None
replace_indices_impl_lock = threading.Lock()


@ffi.def_extern()
def python_replace_indices(s, i, new_s_ptr, new_i_ptr):
    assert replace_indices_impl is not None
    (new_s, new_i) = replace_indices_impl(s, i)
    new_s_ptr[0] = new_s
    new_i_ptr[0] = new_i


class Expr(object):
    _payload: ffi.CData
    _finalizer: Optional[weakref.finalize]

    @singledispatchmethod
    def __init__(self, expression: str, sites: Optional[List[List[int]]] = None):
        json_object = {"expression": expression}
        json_object |= {"sites": sites} if sites is not None else dict()
        self._payload = lib.ls_hs_expr_from_json(json.dumps(json_object).encode("utf-8"))
        assert self._payload != 0
        self._finalizer = weakref.finalize(self, lib.ls_hs_destroy_expr, self._payload)

    @__init__.register
    def _(self, payload: ffi.CData, owner: bool = True):
        self._payload = payload
        if owner:
            self._finalizer = weakref.finalize(self, lib.ls_hs_destroy_expr, self._payload)
        else:
            self._finalizer = None

    def to_json(self) -> str:
        c_str = lib.ls_hs_expr_to_json(self._payload)
        s = ffi.string(c_str).decode("utf-8")
        lib.ls_hs_destroy_string(c_str)
        return s

    def __str__(self):
        """Get the string representation of the underlying expression."""
        c_str = lib.ls_hs_expr_to_string(self._payload)
        s = ffi.string(c_str).decode("utf-8")
        lib.ls_hs_destroy_string(c_str)
        return s

    def replace_indices(self, mapping):
        if isinstance(mapping, dict):
            if len(mapping) == 0:
                return self
            element = next(iter(mapping))
            if isinstance(element, int):

                # Mapping over site indices
                def f(s, i):
                    if i in mapping:
                        return (s, mapping[i])
                    else:
                        return (s, i)

            elif isinstance(element, str):

                # Mapping over spin indices
                def f(s, i):
                    k = _to_spin_index(s)
                    if k in mapping:
                        return (_from_spin_index(mapping[k]), i)
                    else:
                        return (s, i)

            elif isinstance(element, tuple):

                # Mapping over both
                def f(s, i):
                    k = (_to_spin_index(s), i)
                    if k in mapping:
                        (new_s, new_i) = mapping[k]
                        return (_from_spin_index(new_s), new_i)
                    else:
                        return (s, i)

            else:
                raise ValueError("invalid mapping: {}".format(mapping))

            with replace_indices_impl_lock:
                global replace_indices_impl
                replace_indices_impl = f
                r = Expr(lib.ls_hs_replace_indices(self._payload, lib.python_replace_indices))
                replace_indices_impl = None
                return r

        else:
            raise NotImplementedError

    def adjoint(self) -> "Expr":
        return Expr(lib.ls_hs_expr_adjoint(self._payload))

    def __eq__(self, other: object) -> bool:
        _assert_subtype(other, Expr)
        return lib.ls_hs_expr_equal(self._payload, other._payload) != 0

    def __add__(self, other):
        _assert_subtype(other, Expr)
        return Expr(lib.ls_hs_expr_plus(self._payload, other._payload))

    def __sub__(self, other):
        _assert_subtype(other, Expr)
        return Expr(lib.ls_hs_expr_minus(self._payload, other._payload))

    def scale(self, coeff: complex) -> "Expr":
        coeff = complex(coeff)
        c_coeff = ffi.new("ls_hs_scalar const*", [coeff.real, coeff.imag])
        return Expr(lib.ls_hs_expr_scale(c_coeff, self._payload))

    def __mul__(self, other):
        _assert_subtype(other, Expr)
        return Expr(lib.ls_hs_expr_times(self._payload, other._payload))

    def __rmul__(self, other: complex):
        if np.isscalar(other):
            return self.scale(other)
        else:
            return NotImplemented


# str(ls.Expr("5 σ⁺₀ σ⁻₁ + (8 + 3im) σ⁻₁"))
# str(ls.Expr("-2 (c†₀ c₁ + c†₁ c₀)"))


class Operator(LinearOperator):
    _payload: ffi.CData
    _finalizer: weakref.finalize

    # def _make_from_payload(self, basis: Basis, expression: Expr, payload: ffi.CData):
    #     """
    #     !!! warning
    #         This is an internal function. Do not use it directly unless you know what you're doing.

    #     Create a quantum operator from its C representation.
    #     """
    #     self._payload = payload
    #     self._finalizer = weakref.finalize(self, lib.ls_hs_destroy_operator_v2, self._payload)
    #     self.basis = basis
    #     self.expression = expression

    # def _make_from_expression(self, basis: Basis, expression: str, sites: ArrayLike):
    #     """Create a quantum operator from a mathematical expression.

    #     `basis` specifies the Hilbert space on which the operator will be defined. `expression` is a
    #     mathematical expression specifying the interaction, e.g. `"σ⁺₀ σ⁻₁"` or `"n↑₀ n↓₁"`.
    #     """
    #     _assert_subtype(expression, str)
    #     c_sites = _normalize_site_indices(sites)
    #     c_expression = expression.encode("utf-8")
    #     self._payload = lib.ls_hs_create_operator(
    #         basis._payload,
    #         c_expression,
    #         c_sites.shape[0],
    #         c_sites.shape[1],
    #         ffi.from_buffer("int[]", c_sites),
    #     )
    #     self._finalizer = weakref.finalize(self, lib.ls_hs_destroy_operator_v2, self._payload)
    #     self.basis = basis

    @singledispatchmethod
    def __init__(self, basis: Basis, expression: Expr):
        _assert_subtype(basis, Basis)
        _assert_subtype(expression, Expr)
        self._payload = lib.ls_hs_create_operator(basis._payload, expression._payload)
        self._finalizer = weakref.finalize(self, lib.ls_hs_destroy_operator, self._payload)

    @__init__.register
    def _(self, payload: ffi.CData, owner: bool = True):
        self._payload = payload
        if owner:
            self._finalizer = weakref.finalize(self, lib.ls_hs_destroy_operator, self._payload)
        else:
            self._finalizer = None

    @property
    def basis(self):
        return Basis(lib.ls_hs_operator_get_basis(self._payload))

    @property
    def expression(self):
        return Expr(lib.ls_hs_operator_get_expr(self._payload))

    # def adjoint(self):
    #     return Operator(self.basis, lib.ls_hs_operator_hermitian_conjugate(self._payload))

    # @property
    # def is_identity(self) -> bool:
    #     return lib.ls_hs_operator_is_identity(self._payload)

    # @property
    # def is_hermitian(self) -> bool:
    #     return lib.ls_hs_operator_is_hermitian(self._payload)

    def __add__(self, other):
        """Add two operators."""
        _assert_subtype(other, Operator)
        return Operator(self.basis, self.expression + other.expression)

    def __sub__(self, other):
        """Subtract `other` Operator from `self`."""
        _assert_subtype(other, Operator)
        return Operator(self.basis, self.expression - other.expression)

    def scale(self, coeff: complex) -> "Operator":
        return Operator(self.basis, self.expression.scale(coeff))

    def __mul__(self, other):
        """Multiply `self` by `other` Operator."""
        if isinstance(other, Operator):
            return Operator(self.basis, self.expression * other.expression)
        else:
            return NotImplemented

    def __rmul__(self, other: complex):
        """Scale `self` by a scalar `other`."""
        if np.isscalar(other):
            return self.scale(other)
        else:
            return NotImplemented

    def __matmul__(self, other: NDArray[Any]) -> NDArray[Any]:
        return self.apply_to_state_vector(other)

    def __repr__(self):
        return "<Operator defined on {}>".format(self.basis.__class__.__name__)

    def apply_diag_to_basis_state(self, state: int) -> float:
        arr = _basis_state_to_array(state, self.basis.number_words)
        coeffs = ffi.new("chpl_external_array *")
        kernels = lib.ls_hs_internal_get_chpl_kernels()
        kernels.operator_apply_diag(self._payload, 1, arr, coeffs, 0)

        coeffs_arr = _chpl_external_array_as_ndarray(coeffs, np.float64)
        return float(coeffs_arr[0])

    def apply_off_diag_to_basis_state(self, state: int) -> List[Tuple[complex, int]]:
        arr = _basis_state_to_array(state, self.basis.number_words)
        betas = ffi.new("chpl_external_array *")
        coeffs = ffi.new("chpl_external_array *")
        offsets = ffi.new("chpl_external_array *")
        kernels = lib.ls_hs_internal_get_chpl_kernels()
        kernels.operator_apply_off_diag(self._payload, 1, arr, betas, coeffs, offsets, 0)

        offsets_arr = _chpl_external_array_as_ndarray(offsets, np.int64)
        betas_arr = _chpl_external_array_as_ndarray(betas, np.uint64)[: offsets_arr[1]]
        coeffs_arr = _chpl_external_array_as_ndarray(coeffs, np.complex128)[: offsets_arr[1]]
        return list(zip(coeffs_arr, betas_arr))

    def apply_to_state_vector(self, vector: NDArray[np.float64]) -> NDArray[np.float64]:
        self._check_basis_is_built("apply_to_state_vector")
        if vector.dtype != np.float64:
            raise TypeError(
                "expected a NDArray[float64], but got {}[{}]"
                "".format(type(vector).__name__, type(vector.dtype))
            )

        out = np.empty_like(vector)
        x_ptr = ffi.from_buffer("double[]", vector, require_writable=False)
        y_ptr = ffi.from_buffer("double[]", out, require_writable=True)

        kernels = lib.ls_hs_internal_get_chpl_kernels()
        kernels.matrix_vector_product(self._payload, 1, x_ptr, y_ptr)
        return out

    def _check_basis_is_built(self, attribute):
        if not self.basis.is_built:
            raise AttributeError(
                "'Operator' object has no attribute '{}' (did you forget to build the basis?)"
                "".format(attribute),
                name=attribute,
                obj=self,
            )

    @property
    def dtype(self):
        self._check_basis_is_built("dtype")
        return np.dtype("float64")

    @property
    def shape(self):
        self._check_basis_is_built("shape")
        n = self.basis.number_states
        return (n, n)

    def _matvec(self, x):
        return self.apply_to_state_vector(x)


def load_yaml_config(filename: str):
    config = lib.ls_hs_load_yaml_config(filename.encode("utf-8"))
    basis = Basis(lib.ls_hs_clone_basis(config.basis))
    hamiltonian = None
    if config.hamiltonian != 0:
        hamiltonian = Operator(lib.ls_hs_clone_operator(config.hamiltonian))
    if config.observables != 0:
        raise NotImplementedError
    lib.ls_hs_destroy_yaml_config(config)
    Config = namedtuple("Config", ["basis", "hamiltonian", "observables"], defaults=[None, None])
    return Config(basis, hamiltonian)


def test_01():
    basis = SpinBasis(2)
    basis.build()
    a = Expr("2 (σ⁺₀ σ⁻₁ + σ⁺₁ σ⁻₀) + σᶻ₀ σᶻ₁")
    h = Operator(basis, a)
    (e, v) = scipy.sparse.linalg.eigsh(h, k=3, which="SA")
    return h, e, v
