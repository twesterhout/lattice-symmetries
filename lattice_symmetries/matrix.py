import numpy as np, sympy, lattice_symmetries as ls
from numpy.typing import NDArray
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import LinearOperator

from lattice_symmetries.basis import Basis, SpinBasis
from lattice_symmetries.expression import Expr, pauli2nbts


class O(LinearOperator):
    b: Basis; e: Expr; dtype: np.dtype; diag_ctx: any; off_diag_ctx: any
    shape = property(lambda self: (self.b.number_states,) * 2)
    def __init__(self, expr: Expr, basis=None, dtype=None):
        self.e = expr
        self.b = ls.Basis(ls.BasisInfo(expr.number_sites)) if basis is None else basis
        ts = pauli2nbts(expr.raw)
        need_cplx = any(sympy.im(t.v) != sympy.S.Zero for t in ts)
        self.diag_ctx, self.off_diag_ctx = ls.compiler.oc_ctx_t(ts)
        self.dtype = dtype if dtype is not None else \
            np.dtype("complex128") if need_cplx else np.dtype("float64")
    def _matvec(self, x): return self.apply_to_state_vector(x)
    def apply_to_state_vector(self, vector: NDArray, out=None):
        if np.issubdtype(self.dtype, np.complexfloating):
            vector = vector.astype(self.dtype, copy=False)
        return self._prepare_Matvec()(self.b.states, self.b.norms, vector, out=out)
    def _prepare_Matvec(self):
        self.b._check_is_built(); self.b._prepare_bs_ctx(); self.b._prepare_search_ctx()
        return ls.compiler.Matvec(self.diag_ctx, self.off_diag_ctx, self.b.bs_ctx, self.b.search_ctx)
        
Operator = O
