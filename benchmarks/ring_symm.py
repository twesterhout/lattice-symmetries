import igraph as ig, lattice_symmetries as ls, numpy as np, timeit, quspin
from sympy import simplify
from quspin.operators import quantum_LinearOperator, hamiltonian
from quspin.basis import spin_basis_general


def setup(shape):
    # n = shape
    # graph = ig.Graph.Ring(n=n, circular=n > 2)
    # T = np.roll(np.arange(n), shift=-1)
    # P = np.arange(n)[::-1]
    # return [(ls.Permutation(T), ls.Rational(0)), (ls.Permutation(P), ls.Rational(0))]
    graph = ig.Graph.Lattice(tuple(reversed(shape)), circular=True)
    sites = np.arange(np.prod(shape)).reshape(shape)
    Tx = np.roll(sites, shift=-1, axis=1).ravel()
    Ty = np.roll(sites, shift=-1, axis=0).ravel()
    Px = np.flip(sites, axis=1).ravel()
    Py = np.flip(sites, axis=0).ravel()
    # R = np.rot90(sites).ravel()
    return graph, [(ls.Permutation(p), ls.Rational(0)) for p in [Tx, Ty, Px, Py]]


def measure_quspin(shape, matrix=True, dtype=np.float64, **kwargs):
    graph, symms = setup(shape)
    blocks = {
        chr(ord("a") + i) + "block": (p.array_form, int(k * p.order()))
        for i, (p, k) in enumerate(symms)
    }
    n = len(graph.vs)
    basis = spin_basis_general(N=n, pauli=-1, **blocks, **kwargs)
    nearest = [e.tuple for e in graph.es]
    static = [
        ["+-", [[2, i, j] for (i, j) in nearest]],
        ["-+", [[2, i, j] for (i, j) in nearest]],
        ["zz", [[1, i, j] for (i, j) in nearest]],
    ]
    print(basis.Ns)

    kwargs = dict(basis=basis, dtype=dtype, check_herm=False) #check_symm=False,
    with ls.measure_time() as dt:
        h = hamiltonian(static, [], **kwargs) if matrix \
            else quantum_LinearOperator(static, **kwargs)
    print(f"Operator construction took {dt()}")

    rng = np.random.default_rng(5)
    if dtype == np.complex128:
        x = rng.random(2*basis.Ns, dtype=np.float64).view(np.complex128)
        out = np.zeros(basis.Ns, dtype=np.complex128)
    else:
        assert dtype == np.float64
        x = rng.random(basis.Ns, dtype=np.float64)
        out = np.zeros(basis.Ns, dtype=np.float64)

    print("Measuring ...")
    r = timeit.repeat(lambda: h.dot(x, out=out), repeat=1, number=1)
    return np.min(r).item(), np.std(r).item()

def measure_ls(shape,cplx=False):
    graph, symms = setup(shape)
    e = ls.heisenberg(graph)
    b = ls.SpinBasis(len(graph.vs), symmetries=symms)
    o = ls.O(e, b)
    o.b.build()
    print(o.b.number_states)

    rng = np.random.default_rng(5)
    if cplx:
        x = rng.random(2*o.b.number_states, dtype=np.float64).view(np.complex128)
        out = np.zeros(o.b.number_states, dtype=np.complex128)
    else:
        x = rng.random(o.b.number_states, dtype=np.float64)
        out = np.zeros(o.b.number_states, dtype=np.float64)
    assert x.size >= 64
    f = lambda: o.apply_to_state_vector(x, out=out)
    r = timeit.repeat(f, repeat=1, number=1)
    return np.min(r), np.std(r)

if __name__ == "__main__":
    # print(",".join(map(str, measure_quspin(30, matrix=True))))
    print(",".join(map(str, measure_ls((5, 6),cplx=False))))
    print(",".join(map(str, measure_ls((5, 6),cplx=True))))
    # print(",".join(map(str, measure_quspin((5, 6),dtype=np.float64,matrix=False)))) # , Ns_block_est=239123150))))
    # print(",".join(map(str, measure_quspin((5, 6),dtype=np.complex128,matrix=False)))) # , Ns_block_est=239123150))))

