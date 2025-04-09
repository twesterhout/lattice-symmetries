import igraph as ig, lattice_symmetries as ls, numpy as np, timeit, quspin
from sympy import simplify
from quspin.operators import quantum_LinearOperator, hamiltonian
from quspin.basis import spin_basis_general


def measure_quspin(n, matrix=True, dtype=np.float64):
    graph = ig.Graph.Ring(n=n, circular=n > 2)
    basis = spin_basis_general(N=n, pauli=-1)
    nearest = [e.tuple for e in graph.es]
    static = [
        ["+-", [[2, i, j] for (i, j) in nearest]],
        ["-+", [[2, i, j] for (i, j) in nearest]],
        ["zz", [[1, i, j] for (i, j) in nearest]],
    ]

    kwargs = dict(basis=basis, dtype=dtype, check_symm=False, check_herm=False)
    h = hamiltonian(static, [], **kwargs) if matrix \
        else quantum_LinearOperator(static, **kwargs)

    rng = np.random.default_rng(5)
    x = rng.random(basis.Ns, dtype=np.float64)
    out = np.zeros(basis.Ns, dtype=np.float64)
    r = timeit.repeat(lambda: h.dot(x), repeat=8, number=1)
    return np.min(r).item(), np.std(r).item()

def measure_ls(n):
    rng = np.random.default_rng(5)
    h = ls.heisenberg(ig.Graph.Ring(n, circular=n > 2))
    terms = ls.expression.pauli2nbts(simplify(h.raw))
    kernels = ls.compiler.build_kernels()
    p_diag, p_off_diag, _keep_alive = ls.compiler.get_ctxs(terms)
    alpha = np.arange(2**n, dtype=np.uint64)
    x = rng.random(2**n, dtype=np.float64)
    out = np.zeros(2**n, dtype=np.float64)
    assert x.size >= 64
    ffi = ls.COMPILER.ffi
    p_diag.alpha0 = ffi.from_buffer("const uint64_t*", alpha, require_writable=False)
    p_diag.X = ffi.from_buffer("const double*", x, require_writable=False)
    p_off_diag.alpha0 = ffi.from_buffer("const uint64_t*", alpha, require_writable=False)
    p_off_diag.X = ffi.from_buffer("const double*", x, require_writable=False)
    f = lambda: kernels.matvec(out.size, p_diag, p_off_diag, ffi.from_buffer("double*", out, require_writable=True))
    r = timeit.repeat(f, repeat=8, number=1)
    return np.min(r), np.std(r)

if __name__ == "__main__":
    print(",".join(map(str, measure_quspin(20, matrix=True))))
    # print(",".join(map(str, measure_quspin(20, matrix=False))))
    # print(",".join(map(str, measure_ls(20))))

