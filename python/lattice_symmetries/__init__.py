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

from ._ls_hs import ffi, lib

import json
import numpy as np
from numpy.typing import ArrayLike, NDArray
import scipy.sparse.linalg
from scipy.sparse.linalg import LinearOperator
import threading
from typing import overload, Any, List, Optional, Tuple, Union
import weakref

__version__ = "2.0.0"

class _RuntimeInitializer:
    def __init__(self):
        print("Initializing Haskell runtime...")
        lib.ls_hs_init()
        print("Initializing Chapel runtime...")
        lib.ls_chpl_init()

    def __del__(self):
        print("Deinitializing Chapel runtime...")
        lib.ls_chpl_finalize()
        print("Deinitializing Haskell runtime...")
        lib.ls_hs_exit()


_runtime_init = _RuntimeInitializer()
lib.set_python_exception_handler()


def _basis_state_to_array(state: int, number_words: int) -> ffi.CData:
    state = int(state)
    arr = ffi.new("uint64_t[]", number_words)
    for i in range(number_words):
        arr[i] = state & 0xFFFFFFFFFFFFFFFF
        state >>= 64
    return arr


class Basis:
    _payload: ffi.CData
    _finalizer: weakref.finalize
    _representatives: Optional[NDArray[np.uint64]]

    def __init__(self, **kwargs):
        json_string = json.dumps(kwargs).encode("utf-8")
        self._payload = lib.ls_hs_basis_from_json(json_string)
        self._finalizer = weakref.finalize(self, lib.ls_hs_destroy_basis_v2, self._payload)
        self._representatives = None

    @property
    def is_built(self) -> bool:
        return self._representatives is not None

    def build(self) -> None:
        """Generate a list of representatives.

        These can later be accessed using the `number_states` and `states` attributes.
        """
        if not self.is_built:
            # The actual compute-intensive part
            lib.ls_hs_basis_build(self._payload)
            # Create a buffer
            n = self.number_states
            p = lib.ls_hs_basis_states(self._payload)
            buffer = ffi.buffer(p, n * ffi.sizeof("uint64_t"))
            # Create an array
            self._representatives = np.frombuffer(buffer, dtype=np.uint64)
        assert self.is_built

    @property
    def number_states(self) -> int:
        n = lib.ls_hs_basis_number_states(self._payload)
        if n == -1:
            raise ValueError(
                "basis states have not been built yet; "
                "if you wish to do so, use the basis.build() function"
            )
        return n

    @property
    def states(self) -> NDArray[np.uint64]:
        if self._representatives is None:
            raise ValueError(
                "basis states have not been built yet; "
                "if you wish to do so, use the basis.build() function"
            )
        return self._representatives

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

    def state_to_string(self, state: int) -> str:
        arr = _basis_state_to_array(state, self.number_words)
        c_str = lib.ls_hs_basis_state_to_string(self._payload, arr)
        s = ffi.string(c_str).decode("utf-8")
        lib.ls_hs_destroy_string(c_str)
        return s


class SpinBasis(Basis):
    def __init__(
        self,
        number_spins: int,
        hamming_weight: Optional[int] = None,
        spin_inversion: Optional[int] = None,
    ):
        """Create a Hilbert space basis for `number_spins` spin-1/2 particles."""
        super().__init__(
            particle="spin-1/2",
            number_spins=number_spins,
            hamming_weight=hamming_weight,
            spin_inversion=spin_inversion,
        )

    @property
    def has_fixed_hamming_weight(self) -> bool:
        return lib.ls_hs_basis_has_fixed_hamming_weight(self._payload)


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
            particle="spinless-fermion",
            number_sites=number_sites,
            number_particles=number_particles,
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
            particle="spinful-fermion",
            number_sites=number_sites,
            number_particles=number_particles,
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
    weakref.finalize(buf, lambda: lib.ls_hs_destroy_external_array(arr))
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
    _finalizer: weakref.finalize

    def _make_from_payload(self, payload: ffi.CData):
        self._payload = payload
        assert self._payload != 0
        self._finalizer = weakref.finalize(self, lib.ls_hs_destroy_expr, self._payload)

    def _make_from_expression(self, expression: str):
        _assert_subtype(expression, str)
        c_expression = expression.encode("utf-8")
        self._payload = lib.ls_hs_create_expr(c_expression)
        assert self._payload != 0
        self._finalizer = weakref.finalize(self, lib.ls_hs_destroy_expr, self._payload)

    def __init__(self, *args, **kwargs):
        if len(args) == 0 and len(kwargs) == 1 and "payload" in kwargs:
            self._make_from_payload(*args, **kwargs)
        else:
            self._make_from_expression(*args, **kwargs)

    @staticmethod
    def from_json(json_string: str) -> "Expr":
        _assert_subtype(json_string, str)
        c_str = json_string.encode("utf-8")
        return Expr(payload=lib.ls_hs_expr_from_json(c_str))

    def to_json(self) -> str:
        c_str = lib.ls_hs_expr_to_json(self._payload)
        s = ffi.string(c_str).decode("utf-8")
        lib.ls_hs_destroy_string(c_str)
        return s

    def __str__(self):
        """Get the string representation of the underlying expression"""
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
                r = Expr(
                    payload=lib.ls_hs_replace_indices(self._payload, lib.python_replace_indices)
                )
                replace_indices_impl = None
                return r

        else:
            NotImplemented

    def __add__(self, other):
        _assert_subtype(other, Expr)
        return Expr(payload=lib.ls_hs_expr_plus(self._payload, other._payload))

    def __sub__(self, other):
        _assert_subtype(other, Expr)
        return Expr(payload=lib.ls_hs_expr_minus(self._payload, other._payload))

    def scale(self, coeff: complex) -> "Expr":
        coeff = complex(coeff)
        c_coeff = ffi.new("ls_hs_scalar const*", [coeff.real, coeff.imag])
        return Expr(payload=lib.ls_hs_expr_scale(c_coeff, self._payload))

    def __mul__(self, other):
        _assert_subtype(other, Expr)
        return Expr(payload=lib.ls_hs_expr_times(self._payload, other._payload))

    def __rmul__(self, other: complex):
        if np.isscalar(other):
            return self.scale(other)
        else:
            NotImplemented


# str(ls.Expr("5 σ⁺₀ σ⁻₁ + (8 + 3im) σ⁻₁"))
# str(ls.Expr("-2 (c†₀ c₁ + c†₁ c₀)"))


class Operator(LinearOperator):
    _payload: ffi.CData
    _finalizer: weakref.finalize
    basis: Basis
    expression: Expr

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

    def __init__(self, basis: Basis, expression: Expr):
        _assert_subtype(basis, Basis)
        _assert_subtype(expression, Expr)
        self._payload = lib.ls_hs_create_operator(basis._payload, expression._payload)
        self._finalizer = weakref.finalize(self, lib.ls_hs_destroy_operator_v2, self._payload)
        self.basis = basis
        self.expression = expression

    # def __init__(self, basis: Basis, *args, **kwargs):
    #     _assert_subtype(basis, Basis)

    #     if len(args) == 1 and len(kwargs) == 0:
    #         self._make_from_payload(basis, *args, **kwargs)
    #     elif len(args) == 0 and len(kwargs) == 1 and "payload" in kwargs:
    #         self._make_from_payload(basis, *args, **kwargs)
    #     else:
    #         self._make_from_expression(basis, *args, **kwargs)

    def adjoint(self):
        return Operator(self.basis, lib.ls_hs_operator_hermitian_conjugate(self._payload))

    @property
    def is_identity(self) -> bool:
        return lib.ls_hs_operator_is_identity(self._payload)

    @property
    def is_hermitian(self) -> bool:
        return lib.ls_hs_operator_is_hermitian(self._payload)

    def __add__(self, other):
        """Add two operators."""
        _assert_subtype(other, Operator)
        return Operator(self.basis, lib.ls_hs_operator_plus(self._payload, other._payload))

    def __sub__(self, other):
        """Subtract `other` Operator from `self`."""
        _assert_subtype(other, Operator)
        return Operator(self.basis, lib.ls_hs_operator_minus(self._payload, other._payload))

    def scale(self, coeff: complex) -> "Operator":
        coeff = complex(coeff)
        c_coeff = ffi.new("ls_hs_scalar const*", [coeff.real, coeff.imag])
        return Operator(self.basis, lib.ls_hs_operator_scale(c_coeff, self._payload))

    def __mul__(self, other):
        """Multiply `self` by `other` Operator."""
        if isinstance(other, Operator):
            return Operator(self.basis, lib.ls_hs_operator_times(self._payload, other._payload))
        else:
            NotImplemented

    def __rmul__(self, other: complex):
        """Scale `self` by a scalar `other`."""
        if np.isscalar(other):
            return self.scale(other)
        else:
            NotImplemented

    def __matmul__(self, other: NDArray[Any]) -> NDArray[Any]:
        return self.apply_to_state_vector(other)

    def __str__(self):
        """Get the string representation of the underlying expression"""
        c_str = lib.ls_hs_operator_pretty_terms(self._payload)
        s = ffi.string(c_str).decode("utf-8")
        lib.ls_hs_destroy_string(c_str)
        return s

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


def test_01():
    basis = SpinBasis(2)
    basis.build()
    a = Expr("2 (σ⁺₀ σ⁻₁ + σ⁺₁ σ⁻₀) + σᶻ₀ σᶻ₁")
    h = Operator(basis, a)
    (e, v) = scipy.sparse.linalg.eigsh(h, k=3, which="SA")
    return h, e, v
