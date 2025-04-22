import numpy as np, sympy, lattice_symmetries as ls
from numpy.typing import NDArray
from scipy.sparse.linalg import LinearOperator

from lattice_symmetries.basis import Basis, SpinBasis
from lattice_symmetries.expression import Expr, pauli2nbts


class O(LinearOperator):
    b: Basis; e: Expr; dtype: np.dtype; diag_ctx: any; off_diag_ctx: any
    shape = property(lambda self: (self.b.number_states,) * 2)
    basis = property(lambda self: self.b)
    expression = property(lambda self: self.e)
    def __init__(self, expr: Expr, basis=None, dtype=None):
        self.e = expr
        self.b = ls.Basis(ls.BasisInfo(expr.number_sites)) if basis is None else basis
        ts = pauli2nbts(expr.raw)
        need_cplx = any(sympy.im(t.v) != sympy.S.Zero for t in ts)
        self.diag_ctx, self.off_diag_ctx = ls.compiler.oc_ctx_t(ts)
        self.dtype = dtype if dtype is not None else \
            np.dtype("complex128") if need_cplx else np.dtype("float64")
    def _matvec(self, x): return self.apply_to_state_vector(x)
    def apply_to_state_vector(self, x: NDArray, out=None, b=0, n=None):
        if np.issubdtype(self.dtype, np.complexfloating): x = x.astype(self.dtype, copy=False)
        if n is None: n = x.size - b
        alpha0, norm0, x0 = self.b.states[b:b + n], self.b.norms[b:b + n], x[b:b + n]
        return self._prepare_Matvec()(alpha0, norm0, x0, x, out=out)
    def apply_to_state_vector_parallel(self, x: NDArray, out: NDArray, batch_size=64):
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        sz, rk = comm.Get_size(), comm.Get_rank()
        rs, o, n = [], 0, x.size
        t = 0
        while o < n:
            my_n = min(batch_size, (n - o) // sz)
            if my_n >= 64:
                my_o = o + rk * my_n
                with ls.measure_time() as dt:
                    self.apply_to_state_vector(x, out=out[my_o:my_o + my_n], b=my_o, n=my_n)
                t += dt()
                rs.append(comm.Iallgather(MPI.IN_PLACE, out[o:o + sz * my_n]))
                o += sz * my_n
            else:
                with ls.measure_time() as dt:
                    out[o:n] = self.apply_to_state_vector(x, b=o, n=n - o)
                t += dt()
                o += n - o
        MPI.Request.Waitall(rs)
        logger.trace(f"{rk}/{sz}: time in apply_to_state_vector {t}")
    def _prepare_Matvec(self):
        self.b._check_is_built(); self.b._prepare_bs_ctx(); self.b._prepare_search_ctx()
        return ls.compiler.Matvec(self.diag_ctx, self.off_diag_ctx, self.b.bs_ctx, self.b.search_ctx)
        
Operator = O
