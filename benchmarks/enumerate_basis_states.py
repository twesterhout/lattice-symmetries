import igraph as ig, lattice_symmetries as ls, numpy as np, time, timeit, quspin
from sympy import simplify, Rational
from sympy.combinatorics import Permutation
from functools import reduce
from quspin.operators import quantum_LinearOperator, hamiltonian
from quspin.basis import spin_basis_general


def measure_quspin(dim):
    graph = ig.Graph.Lattice(dim, circular=True)
    sites = np.arange(np.prod(dim)).reshape(dim)
    Tx = np.roll(sites, shift=-1, axis=1).ravel()
    Ty = np.roll(sites, shift=-1, axis=0).ravel()
    Px = np.flip(sites, axis=1).ravel()
    Py = np.flip(sites, axis=0).ravel()
    Z = -(sites + 1).ravel()
    ts = []
    for _ in range(1):
        basis = spin_basis_general(N=np.prod(dim), pauli=-1, make_basis=False, 
            kxblock=(Tx, 0),
            kyblock=(Ty, 0), pxblock=(Px, 0), pyblock=(Py, 0),
            zblock=(Z, 0),
            )
        tick = time.perf_counter()
        basis.make()
        tock = time.perf_counter()
        ts.append(tock - tick)
    print(basis.Ns)
    # print(basis.states)
    return np.min(ts), np.std(ts)

def measure_ls(dim):
    graph = ig.Graph.Lattice(dim, circular=True)
    sites = np.arange(np.prod(dim)).reshape(dim)
    Tx = Permutation(np.roll(sites, shift=-1, axis=1).ravel())
    Ty = Permutation(np.roll(sites, shift=-1, axis=0).ravel())
    Px = Permutation(np.flip(sites, axis=1).ravel())
    Py = Permutation(np.flip(sites, axis=0).ravel())
    info = ls.BasisInfo(bits=np.prod(dim), inversion=1,
        symmetries=[(Tx, Rational(0)), (Ty, Rational(0)), (Px, Rational(0)), (Py, Rational(0))
                    ])
    kernels = ls.compiler.build_kernels()
    f = ls.compiler.EnumerateStates(info, kernels)
    
    ts = []
    for _ in range(1):
        tick = time.perf_counter()
        states, norms = f()
        tock = time.perf_counter()
        ts.append(tock - tick)
    print(states.size)
    # print(states)
    return np.min(ts), np.std(ts)


    
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
    # print(",".join(map(str, measure_quspin((6, 6)))))
    print(",".join(map(str, measure_ls((7, 6)))))
    # print(",".join(map(str, measure_quspin(20, matrix=False))))
    # print(",".join(map(str, measure_ls(20))))

