from lattice_symmetries.basis import Basis
from lattice_symmetries.expression import Expr, pauli_expression_to_nonbranching_terms
from lattice_symmetries._kernels import LoweredOperator, LoweredSymmetries
import numpy as np
from numpy.typing import NDArray
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import LinearOperator


class Operator(LinearOperator):
    basis: Basis
    expression: Expr
    dtype: np.dtype
    lowered_operator: LoweredOperator | None

    def __init__(
        self,
        expression: Expr,
        basis: Basis | None = None,
        dtype: np.dtype = np.dtype("float32"),
    ):
        self.expression = expression
        self.basis = basis
        self.dtype = dtype
        self.lowered_operator = None

    @property
    def shape(self) -> tuple[int, int]:
        n = self.basis.number_states
        return (n, n)

    def apply_to_state_vector(self, vector: NDArray, out: NDArray | None = None, verbose: bool = False) -> NDArray:
        self.basis.check_is_built()

        if self.lowered_operator is None:
            terms = pauli_expression_to_nonbranching_terms(self.expression.raw)
            symm = LoweredSymmetries(self.basis.info.symmetries) if self.basis.info.has_permutation_symmetries else None
            self.lowered_operator = LoweredOperator(info=self.basis.info, terms=terms, symm=symm, state_to_index_info=self.basis.state_to_index_info, verbose=verbose)

        vector = np.asarray(vector, dtype=self.dtype, order="F")
        if out is None:
            out = np.zeros(self.basis.states.size, dtype=self.dtype)
        else:
            out = np.asarray(out, dtype=self.dtype, order="F")

        self.lowered_operator.apply(self.basis.states, self.basis.norms, vector, out)
        return out

    def _matvec(self, x):
        return self.apply_to_state_vector(x)