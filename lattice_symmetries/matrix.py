import numpy as np, sympy, lattice_symmetries as ls
from loguru import logger
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
        """
        Initialize a new quantum operator.

        Creates a linear operator representing a quantum operator described by the
        given expression in the provided basis.

        Parameters
        ----------
        expr : Expr
            The symbolic expression defining the quantum operator, typically
            containing Pauli matrices.
        basis : Basis, optional
            The basis in which the operator will act. If None, a new basis is
            automatically created using the number of sites in the expression,
            without any symmetries.
        dtype : numpy.dtype, optional
            The data type to use for operator coefficients. If None, the type is
            automatically determined: complex128 if the expression contains complex
            coefficients, otherwise float64.

        Notes
        -----
        When basis=None, a new basis is created that doesn't use any symmetries,
        which can be much less efficient for large systems. Consider creating an
        appropriate basis with symmetries for better performance.

        The operator is internally represented using non-branching terms that
        are optimized for efficient matrix-vector multiplication.
        """
        self.e = expr
        self.b = ls.Basis(ls.BasisInfo(expr.number_sites)) if basis is None else basis
        ts = pauli2nbts(expr.raw)
        need_cplx = any(sympy.im(t.v) != sympy.S.Zero for t in ts)
        self.diag_ctx, self.off_diag_ctx = ls.compiler.oc_ctx_t(ts)
        self.dtype = dtype if dtype is not None else \
            np.dtype("complex128") if need_cplx else np.dtype("float64")
    def _matvec(self, x):
        """
        Matrix-vector multiplication implementation for LinearOperator interface.

        This is an internal method that implements the abstract _matvec method
        required by scipy.sparse.linalg.LinearOperator. It simply delegates to
        the more feature-rich apply_to_state_vector method with default parameters.

        Parameters
        ----------
        x : NDArray
            The input state vector to which the operator will be applied.

        Returns
        -------
        NDArray
            The resulting vector after applying the operator to the input state.

        See Also
        --------
        apply_to_state_vector : The full implementation with additional parameters
                                for more fine-grained control over the operation.
        """
        return self.apply_to_state_vector(x)
    def apply_to_state_vector(self, x: NDArray, out=None, b=0, n=None):
        """
        Apply this operator to a state vector.

        This method applies the quantum operator to a given state vector, which can be
        thought of as multiplying the operator's matrix representation by the vector.
        The implementation leverages optimized matrix-vector product kernels.

        Parameters
        ----------
        x : NDArray
            The input state vector to which the operator will be applied.
        out : NDArray, optional
            Output array for the result. If provided, it must have appropriate shape
            and dtype. If None, a new array is allocated.
        b : int, default=0
            Starting index in the input vector. Useful for processing only a 
            segment of the full state vector.
        n : int, optional
            Number of elements to process starting from index `b`. 
            If None, processes all elements from `b` to the end of the vector.

        Returns
        -------
        NDArray
            The resulting vector after applying the operator to the input state.

        Notes
        -----
        For complex operators, the input vector is automatically cast to the
        appropriate complex datatype if needed. The function uses the basis states
        and norms stored in the operator's basis for the computation.
        """
        if np.issubdtype(self.dtype, np.complexfloating): x = x.astype(self.dtype, copy=False)
        if n is None: n = x.size - b
        alpha0, norm0, x0 = self.b.states[b:b + n], self.b.norms[b:b + n], x[b:b + n]
        return self._prepare_Matvec()(alpha0, norm0, x0, x, out=out)
    def apply_to_state_vector_parallel(self, x: NDArray, out: NDArray, batch_size=2**20):
        """
        Apply this operator to a state vector in parallel using MPI.

        This method applies the quantum operator to a given state vector in parallel
        using MPI for distributed computing, which can significantly improve performance
        for large state vectors. The work is divided across MPI processes, with each
        process handling a portion of the computation.

        Parameters
        ----------
        x : NDArray
            The full input state vector to which the operator will be applied.
            Each MPI process must have a complete copy of this vector.
        out : NDArray
            The full output array for the result, representing the entire Hilbert space.
            Must be pre-allocated with the appropriate shape and dtype. Each MPI process
            will compute its contribution to this full result, and the results will be
            automatically gathered and combined.
        batch_size : int, default=2**20
            Minimum number of elements to process in each batch. Controls the 
            granularity of parallelization.

        Notes
        -----
        WARNING: THIS FUNCTION IS EXPERIMENTAL. The API is likely to change in future
        releases.

        Unlike apply_to_state_vector, this method requires a pre-allocated output array
        and uses MPI for parallel execution. All MPI processes must have access to the
        complete input vector, and the final result will be the complete output vector 
        across the entire Hilbert space, not just a local chunk.

        See Also
        --------
        apply_to_state_vector : The sequential implementation that provides more
                               fine-grained control over the operation.
        """
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
