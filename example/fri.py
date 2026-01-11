import numpy as np, scipy.sparse.linalg, sympy, igraph as ig, quspin, uniplot
from sympy import Rational
from sympy.combinatorics import Permutation
from numba import jit

def quspin_heisenberg(graph, symmetries=None):
    import quspin.basis, quspin.operators
    bits = len(graph.vs)
    if symmetries is None: symmetries = info.symmetries
    blocks = {
        chr(ord("a") + i) + "block": (p.array_form, int(k * p.order()))
        for i, (p, k) in enumerate(symmetries)
    }
    # print(blocks)
    basis = quspin.basis.spin_basis_general(bits, Nup=bits // 2, pauli=-1, **blocks)
    nearest = [e.tuple for e in graph.es]
    static = [
        ["+-", [[2, i, j] for (i, j) in nearest]],
        ["-+", [[2, i, j] for (i, j) in nearest]],
        ["zz", [[1, i, j] for (i, j) in nearest]],
        # ["z", [[-0.125, i] for i in range(info.bits)]]
    ]
    # QuSpin doesn't detect that the Hamiltonian is real if some characters are -1...
    hamiltonian = quspin.operators.hamiltonian(static, [], basis=basis, dtype=np.dtype("complex128"), check_symm=False, check_herm=False)
    return basis, hamiltonian

def naive_compress(x, m):
    order = np.argsort(np.abs(x))
    y = x.copy()
    y[order[:len(x) - m]] = 0
    return y

@jit
def pivotal(x):
    y = np.copy(x)
    i, j, k = 0, 1, 2
    d = y.size
    a, b = y[i], y[j]
    while k < d:
        if k < d and (a < 1e-12 or a > 1 - 1e-12): a = y[k]; i = k; k += 1
        if k < d and (b < 1e-12 or b > 1 - 1e-12): b = y[k]; j = k; k += 1
        u = np.random.rand()
        add = a + b
        if (add > 1) and (add < 2):
            if u < (1 - b)/(2 - add): b = add - 1; a = 1
            else: a = add - 1; b = 1
        elif (add > 0) and (add <= 1):
            if u < b / add: b = add; a = 0
            else: a = add; b = 0
        y[i], y[j] = a, b
    return y

def pivotal_compress(x, m):
    if m >= x.size: return x
    x_abs = np.abs(x)
    order = np.argsort(x_abs)[::-1] # decreasing order

    # print("x", x)
    # print("x_abs", x_abs)
    # print("order", order)

    x_abs_ordered = x_abs[order]
    cumsums_ordered = np.cumsum(x_abs_ordered[::-1])[::-1]
    condition = x_abs_ordered[:m] >= (1 / (m - np.arange(m))) * cumsums_ordered[:m]
    q = np.argmin(condition)

    D = order[:q]

    # Computing p
    p = x_abs.copy()
    p[D] = 0
    pnorm = cumsums_ordered[q]
    p *= (m - q) / pnorm

    # print(np.sum(p))
    # print(p)
    S = np.rint(pivotal(p)).nonzero()[0]
    # print(S)

    psi = np.zeros_like(x)
    psi[D] = x[D]
    psi[S] = x[S] / p[S]
    # print(psi)
    return psi

    
n = 28
p = Permutation(np.roll(np.arange(n), shift=-1))
p2 = Permutation(np.arange(n)[::-1])
g = ig.Graph.Ring(n, circular=True)
basis, hamiltonian = quspin_heisenberg(g, symmetries=[(p, Rational(0, 2)), (p2, Rational(0))])
m = hamiltonian.tocsr()
m.data = np.ascontiguousarray(m.data.real)
# print(scipy.sparse.linalg.eigsh(m, k=1, return_eigenvectors=False, which="LA"))
es, vs = scipy.sparse.linalg.eigsh(m, k=1, which="SA")
print(es)

print(np.abs(vs[:, 0] - np.mean(np.stack([pivotal_compress(vs[:, 0], round(0.5 * basis.Ns)) for _ in range(200)]), axis=0)).sum())


rng = np.random.default_rng(seed=123)
x = rng.random(basis.Ns, dtype=np.float64)
x /= np.linalg.norm(x)

rs1 = []
y = x
y = y / np.linalg.norm(y, ord=1)
u = None
# u = rng.random(basis.Ns)
# rs1.append(np.vdot(u, m @ y) / np.vdot(u, y))
for i in range(300):
    # yn = naive_compress(m @ y, m=round(0.25 * basis.Ns))
    yn = pivotal_compress(m @ y, m=round(0.01 * basis.Ns))
    if u is not None:
        rs1.append(np.vdot(u, yn) / np.vdot(u, y))
    yn = yn / np.linalg.norm(yn, ord=1)
    # rs1.append(np.vdot(u, yn) / np.vdot(u, y))
    y = yn
    if i == 100:
        u = y
        

print(es.item())
print(np.mean(rs1[-100:]), np.std(rs1[-100:]), np.abs(es.item() - np.mean(rs1[-100:])) / np.abs(es.item()))
uniplot.plot(np.abs(rs1 - es), y_as_log=True, character_set="braille")
uniplot.histogram(rs1[-100:] - es)


    
    



rng = np.random.default_rng(seed=123)
x = rng.random(basis.Ns, dtype=np.float64)
x /= np.linalg.norm(x)

rs2 = []
y = x
rs.append(np.vdot(y, m @ y))
for i in range(100):
    y = m @ y
    # y /= np.linalg.norm(y)
    # dumb compress
    # print(np.sum(np.abs(y)[np.argsort(np.abs(y))[:5000]]**2))
    # y[np.argsort(np.abs(y))[:5000]] = 0.0
    # y = naive_compress(y, basis.Ns // 10)
    y = pivotal_compress(y, basis.Ns // 5) # naive_compress(y, basis.Ns // 10)
    y /= np.linalg.norm(y)
    rs2.append(np.vdot(y, m @ y))

rs3 = []
y = x
rs.append(np.vdot(y, m @ y))
for i in range(100):
    y = m @ y
    # y /= np.linalg.norm(y)
    # dumb compress
    # print(np.sum(np.abs(y)[np.argsort(np.abs(y))[:5000]]**2))
    # y[np.argsort(np.abs(y))[:5000]] = 0.0
    y = naive_compress(y, basis.Ns // 5)
    # y = pivotal_compress(y, basis.Ns // 10) # naive_compress(y, basis.Ns // 10)
    y /= np.linalg.norm(y)
    rs3.append(np.vdot(y, m @ y))

uniplot.plot([rs1 - es, rs2 - es, rs3 - es], x_min=20, y_as_log=True, legend_labels=["1", "2", "3"], character_set="braille")
