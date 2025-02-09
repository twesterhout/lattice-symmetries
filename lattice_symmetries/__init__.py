# Copyright (c) 2022-2024, Tom Westerhout
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

__version__ = "3.0.0"

from sympy.physics.quantum.pauli import SigmaX, SigmaY, SigmaZ

# Fix printing of pauli expressions
def _proper_print_contents(t):
    def f(self, printer, *args): return t.__name__ + "(" + printer._print(self.name) + ")"
    return f
SigmaX._print_contents = _proper_print_contents(SigmaX)
SigmaY._print_contents = _proper_print_contents(SigmaY)
SigmaZ._print_contents = _proper_print_contents(SigmaZ)


# from . import _ls, _numpy_helper, _parser, _axpy, _benes, _representation, _kernels, basis, expression, matrix
from . import _parser, expression, _benes

# from ._axpy import axpy
# from ._kernels import BasisInfo
from ._benes import BenesNetwork, perm2benes
# from ._representation import generate_representation
# from .basis import (
#     Basis,
#     SpinBasis,
#     fixed_hamming_state_to_index,
#     fixed_hamming_index_to_state,
#     enumerate_basis_states,
# )
from .expression import Expr, heisenberg, ising
# from .matrix import Operator

# _ls.lib.ls_chpl_init()


# result = _ls.ffi.new("ls_numpy_array_1d *")
# _ls.lib.the_ultimate_solution(_ls.lib.ls_alloc_numpy_array_1d, result)
# print(result.handle)
# print(_ls.ffi.from_handle(result.handle))
# _ls.lib.ls_PyObject_decref(result.handle)
