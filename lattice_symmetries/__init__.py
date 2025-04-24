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

import contextlib, time
from sympy import Rational
from sympy.physics.quantum.pauli import SigmaX, SigmaY, SigmaZ, SigmaPlus, SigmaMinus
from sympy.combinatorics import Permutation, PermutationGroup

# Fix printing of pauli expressions
def _fix_print():
    def _proper(t):
        def f(self, p, *args): return t.__name__ + "(" + p._print(self.name) + ")"
        return f
    for o in [SigmaX, SigmaY, SigmaZ, SigmaPlus, SigmaMinus]: o._print_contents = _proper(o)
_fix_print()

@contextlib.contextmanager
def measure_time(): tick = tock = time.perf_counter(); yield lambda: tock - tick; tock = time.perf_counter() 

from . import _benes, _parser, _representation, expression, compiler, basis, matrix
from ._benes import BenesNetwork, perm2benes
from ._representation import generate_representation
from .expression import Expr, pauli2nbts, heisenberg, ising
from .compiler import COMPILER, KERNELS, BasisInfo
from .basis import B, Basis, SpinBasis
from .matrix import O, Operator
