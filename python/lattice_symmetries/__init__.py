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
from typing import Optional, Tuple, Union
import weakref


class _RuntimeInitializer:
    def __init__(self):
        print("Initializing Haskell runtime...")
        lib.ls_hs_init()

    def __del__(self):
        print("Deinitializing Haskell runtime...")
        lib.ls_hs_exit()


_runtime_init = _RuntimeInitializer()
lib.set_python_exception_handler()


class Basis:
    _payload: ffi.CData
    _finalizer: weakref.finalize

    def __init__(self, **kwargs):
        json_string = json.dumps(kwargs).encode("utf-8")
        self._payload = lib.ls_hs_basis_from_json(json_string)
        self._finalizer = weakref.finalize(self, lib.ls_hs_destroy_basis_v2, self._payload)


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















