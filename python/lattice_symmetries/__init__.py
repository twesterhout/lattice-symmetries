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
from typing import overload, Optional, Tuple, Union
import weakref


class _RuntimeInitializer:
    def __init__(self):
        print("Initializing Haskell runtime...")
        lib.ls_hs_init()
        print("Initializing Chapel runtime...")
        lib.chpl_library_init_wrapper()

    def __del__(self):
        print("Deinitializing Haskell runtime...")
        lib.ls_hs_exit()
        print("Deinitializing Chapel runtime...")
        lib.chpl_library_finalize()


_runtime_init = _RuntimeInitializer()
lib.set_python_exception_handler()


class Basis:
    _payload: ffi.CData
    _finalizer: weakref.finalize
    _representatives: Optional[NDArray[np.uint64]]

    def __init__(self, **kwargs):
        json_string = json.dumps(kwargs).encode("utf-8")
        self._payload = lib.ls_hs_basis_from_json(json_string)
        self._finalizer = weakref.finalize(self, lib.ls_hs_destroy_basis_v2, self._payload)
        self._representatives = None

    def build(self) -> None:
        """Generate a list of representatives.

        These can later be accessed using the `number_states` and `states` attributes.
        """
        if self._representatives is None:
            # The actual compute-intensive part
            lib.ls_hs_basis_build(self._payload)
            # Create a buffer
            n = self.number_states
            p = lib.ls_hs_basis_states(self._payload)
            buffer = ffi.buffer(p, n * ffi.sizeof("uint64_t"))
            # Create an array
            self._representatives = np.frombuffer(buffer, dtype=np.uint64)

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
        state = int(state)
        n = self.number_words
        arr = ffi.new("uint64_t[]", n)
        for i in range(n):
            arr[i] = state & 0xFFFFFFFFFFFFFFFF
            state >>= 64
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


class Operator:
    _payload: ffi.CData
    _finalizer: weakref.finalize
    basis: Basis

    def _make_from_payload(self, basis: Basis, payload: ffi.CData):
        """
        !!! warning
            This is an internal function. Do not use it directly unless you know what you're doing.

        Create a quantum operator from its C representation.
        """
        self._payload = payload
        self._finalizer = weakref.finalize(self, lib.ls_hs_destroy_operator_v2, self._payload)
        self.basis = basis

    def _make_from_expression(self, basis: Basis, expression: str, sites: ArrayLike):
        """Create a quantum operator from a mathematical expression.

        `basis` specifies the Hilbert space on which the operator will be defined. `expression` is a
        mathematical expression specifying the interaction, e.g. `"σ⁺₀ σ⁻₁"` or `"n↑₀ n↓₁"`.
        """
        _assert_subtype(expression, str)
        c_sites = _normalize_site_indices(sites)
        c_expression = expression.encode("utf-8")
        self._payload = lib.ls_hs_create_operator(
            basis._payload,
            c_expression,
            c_sites.shape[0],
            c_sites.shape[1],
            ffi.from_buffer("int[]", c_sites),
        )
        self._finalizer = weakref.finalize(self, lib.ls_hs_destroy_operator_v2, self._payload)
        self.basis = basis

    def __init__(self, basis: Basis, *args, **kwargs):
        _assert_subtype(basis, Basis)

        if len(args) == 1 and len(kwargs) == 0:
            self._make_from_payload(basis, *args, **kwargs)
        elif len(args) == 0 and len(kwargs) == 1 and "payload" in kwargs:
            self._make_from_payload(basis, *args, **kwargs)
        else:
            self._make_from_expression(basis, *args, **kwargs)

    def conj(self):
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

    def __mul__(self, other):
        """Multiply `self` by `other` Operator."""
        _assert_subtype(other, Operator)
        return Operator(self.basis, lib.ls_hs_operator_times(self._payload, other._payload))

    def __rmul__(self, other: complex):
        """Scale `self` by a scalar `other`."""
        coeff = complex(other)
        c_coeff = ffi.new("ls_hs_scalar const*", [coeff.real, coeff.imag])
        return Operator(self.basis, lib.ls_hs_operator_scale(c_coeff, self._payload))

    def __str__(self):
        """Get the string representation of the underlying expression"""
        c_str = lib.ls_hs_operator_pretty_terms(self._payload)
        s = ffi.string(c_str).decode("utf-8")
        lib.ls_hs_destroy_string(c_str)
        return s
