import numpy as np
import lattice_symmetries as ls
import pytest
import math

ls.enable_logging()

import systems


def test_empty():
    with pytest.raises(ls.LatticeSymmetriesException):
        ls.SpinBasis(ls.Group([]), number_spins=0)


def test_huge():
    with pytest.raises(ls.LatticeSymmetriesException):
        ls.SpinBasis(ls.Group([]), number_spins=1000)


def test_1_spin():
    basis = ls.SpinBasis(ls.Group([]), number_spins=1)
    basis.build()
    assert basis.states.tolist() == [0, 1]
    assert basis.state_info(0) == (0, 1.0, 1.0)
    assert basis.state_info(1) == (1, 1.0, 1.0)
    basis = ls.SpinBasis(ls.Group([]), number_spins=1, hamming_weight=0)
    basis.build()
    assert basis.states.tolist() == [0]
    assert basis.state_info(0) == (0, 1.0, 1.0)
    basis = ls.SpinBasis(ls.Group([]), number_spins=1, spin_inversion=-1)
    basis.build()
    assert basis.states.tolist() == [0]
    assert basis.state_info(0) == (0, 1.0, pytest.approx(1 / math.sqrt(2)))
    assert basis.state_info(1) == (0, -1.0, pytest.approx(1 / math.sqrt(2)))


def test_2_spins():
    basis = ls.SpinBasis(ls.Group([]), number_spins=2)
    basis.build()
    assert basis.states.tolist() == [0, 1, 2, 3]
    with pytest.raises(ls.LatticeSymmetriesException):
        ls.SpinBasis(ls.Group([]), number_spins=2, hamming_weight=2, spin_inversion=1)
    with pytest.raises(ls.LatticeSymmetriesException):
        ls.SpinBasis(ls.Group([]), number_spins=2, hamming_weight=2, spin_inversion=-1)
    basis = ls.SpinBasis(ls.Group([ls.Symmetry([1, 0], sector=1)]), number_spins=2)
    basis.build()
    assert basis.states.tolist() == [1]
    assert basis.state_info(0) == (0, 1.0, 0.0)
    assert basis.state_info(1) == (1, 1.0, pytest.approx(1 / math.sqrt(2)))
    assert basis.state_info(2) == (1, -1.0, pytest.approx(1 / math.sqrt(2)))
    assert basis.state_info(3) == (3, 1.0, 0.0)


def test_4_spins():
    # fmt: off
    matrix = np.array([[1,  0,  0, 0],
                       [0, -1,  2, 0],
                       [0,  2, -1, 0],
                       [0,  0,  0, 1]])
    # fmt: on
    number_spins = 4
    edges = [(i, (i + 1) % number_spins) for i in range(number_spins)]

    basis = ls.SpinBasis(ls.Group([]), number_spins=4, hamming_weight=2)
    basis.build()
    assert basis.number_states == 6
    operator = ls.Operator(basis, [ls.Interaction(matrix, edges)])
    assert np.isclose(ls.diagonalize(operator, k=1)[0], -8)

    basis = ls.SpinBasis(ls.Group([]), number_spins=4, hamming_weight=2, spin_inversion=1)
    basis.build()
    assert basis.number_states == 3
    operator = ls.Operator(basis, [ls.Interaction(matrix, edges)])
    assert np.isclose(ls.diagonalize(operator, k=1)[0], -8)

    T = ls.Symmetry([1, 2, 3, 0], sector=0)
    basis = ls.SpinBasis(ls.Group([T]), number_spins=4, hamming_weight=2, spin_inversion=1)
    basis.build()
    assert basis.number_states == 2
    operator = ls.Operator(basis, [ls.Interaction(matrix, edges)])
    assert np.isclose(ls.diagonalize(operator, k=1)[0], -8)


def test_index():
    L_x, L_y = (4, 6)
    backend = "ls"
    symmetries = systems.square_lattice_symmetries(L_x, L_y)
    nearest, _ = systems.square_lattice_edges(L_x, L_y)
    basis = systems.make_basis(
        symmetries,
        backend=backend,
        number_spins=L_x * L_y,
        hamming_weight=(L_x * L_y) // 2,
    )
    # print(basis.number_states)
    hamiltonian = systems.make_heisenberg(basis, nearest, backend=backend)

    # indices = ls.batched_index(basis, basis.states)
    # assert np.all(indices == np.arange(basis.number_states, dtype=np.uint64))
    for i in range(basis.number_states):
        index = basis.index(basis.states[i])
        assert index == i

    assert np.all(basis.batched_index(basis.states) == np.arange(basis.number_states))

    spins = np.zeros((10000, 8), dtype=np.uint64)
    spins[:, 0] = basis.states[:10000]
    basis.batched_state_info(spins)
    # evals, evecs = hamiltonian.eigsh(k=1, which='SA')
    # evals, evecs = ls.diagonalize(hamiltonian)
    # print(evals)


def notest_construction():
    symmetries = [
        # ls.Symmetry([18, 19,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17], sector=5),
        # ls.Symmetry([19, 18, 17, 16, 15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1, 0], sector=0)
    ]
    basis = ls.SpinBasis(
        ls.Group(symmetries), number_spins=20, hamming_weight=10, spin_inversion=None
    )
    basis.build()

    # fmt: off
    interactions = [
        ls.Interaction([[0.25, 0, 0, 0], [0, -0.25, 0.5, 0], [0,  0.5 , -0.25,  0.  ], [ 0.  ,  0.
            ,  0.  ,  0.25]], [[0, 1], [2, 3], [4, 5], [6, 7], [8, 9], [10, 11], [12, 13], [14, 15],
                [16, 17], [18, 19]]),

        ls.Interaction([[ 0.25,  0.  ,  0.  ,  0.  ],
               [ 0.  , -0.25,  0.5 ,  0.  ],
               [ 0.  ,  0.5 , -0.25,  0.  ],
               [ 0.  ,  0.  ,  0.  ,  0.25]],
            [[0, 2], [1, 3], [2, 4], [3, 5], [4, 6], [5, 7], [6, 8], [7, 9], [8, 10], [9, 11], [10, 12], [11, 13], [12, 14], [13, 15], [14, 16], [15, 17], [16, 18], [17, 19], [18, 0], [19, 1]]
        ),
        ls.Interaction([[-0.0625, -0.    , -0.    , -0.    ],
               [-0.    ,  0.0625, -0.125 , -0.    ],
               [-0.    , -0.125 ,  0.0625, -0.    ],
               [-0.    , -0.    , -0.    , -0.0625]],
            [[0, 4], [1, 5], [2, 6], [3, 7], [4, 8], [5, 9], [6, 10], [7, 11], [8, 12], [9, 13], [10, 14], [11, 15], [12, 16], [13, 17], [14, 18], [15, 19], [16, 0], [17, 1], [18, 2], [19, 3]]
        ),
        ls.Interaction([[ 0.02777778,  0.        ,  0.        ,  0.        ],
               [ 0.        , -0.02777778,  0.05555556,  0.        ],
               [ 0.        ,  0.05555556, -0.02777778,  0.        ],
               [ 0.        ,  0.        ,  0.        ,  0.02777778]],
            [[0, 6], [1, 7], [2, 8], [3, 9], [4, 10], [5, 11], [6, 12], [7, 13], [8, 14], [9, 15],
                [10, 16], [11, 17], [12, 18], [13, 19], [14, 0], [15, 1], [16, 2], [17, 3], [18, 4],
                [19, 5]]),
        ls.Interaction([[-0.015625, -0.      , -0.      , -0.      ],
               [-0.      ,  0.015625, -0.03125 , -0.      ],
               [-0.      , -0.03125 ,  0.015625, -0.      ],
               [-0.      , -0.      , -0.      , -0.015625]],
            [[0, 8], [1, 9], [2, 10], [3, 11], [4, 12], [5, 13], [6, 14], [7, 15], [8, 16], [9, 17],
                [10, 18], [11, 19], [12, 0], [13, 1], [14, 2], [15, 3], [16, 4], [17, 5], [18, 6],
                [19, 7]]),
        ls.Interaction([[ 0.01,  0.  ,  0.  ,  0.  ],
               [ 0.  , -0.01,  0.02,  0.  ],
               [ 0.  ,  0.02, -0.01,  0.  ],
               [ 0.  ,  0.  ,  0.  ,  0.01]],
            [[0, 10], [1, 11], [2, 12], [3, 13], [4, 14], [5, 15], [6, 16], [7, 17], [8, 18], [9,
                19], [10, 0], [11, 1], [12, 2], [13, 3], [14, 4], [15, 5], [16, 6], [17, 7], [18,
                    8], [19, 9]]),
        ls.Interaction([[-0.00694444, -0.        , -0.        , -0.        ],
                   [-0.        ,  0.00694444, -0.01388889, -0.        ],
                   [-0.        , -0.01388889,  0.00694444, -0.        ],
                   [-0.        , -0.        , -0.        , -0.00694444]],
            [[0, 12], [1, 13], [2, 14], [3, 15], [4, 16], [5, 17], [6, 18], [7, 19], [8, 0], [9, 1],
                [10, 2], [11, 3], [12, 4], [13, 5], [14, 6], [15, 7], [16, 8], [17, 9], [18, 10],
                [19, 11]]),
        ls.Interaction([[ 0.00510204,  0.        ,  0.        ,  0.        ],
               [ 0.        , -0.00510204,  0.01020408,  0.        ],
               [ 0.        ,  0.01020408, -0.00510204,  0.        ],
               [ 0.        ,  0.        ,  0.        ,  0.00510204]],
            [[0, 14], [1, 15], [2, 16], [3, 17], [4, 18], [5, 19], [6, 0], [7, 1], [8, 2], [9, 3],
                [10, 4], [11, 5], [12, 6], [13, 7], [14, 8], [15, 9], [16, 10], [17, 11], [18, 12],
                [19, 13]]),
        ls.Interaction([[-0.00390625, -0.        , -0.        , -0.        ],
               [-0.        ,  0.00390625, -0.0078125 , -0.        ],
               [-0.        , -0.0078125 ,  0.00390625, -0.        ],
               [-0.        , -0.        , -0.        , -0.00390625]],
            [[0, 16], [1, 17], [2, 18], [3, 19], [4, 0], [5, 1], [6, 2], [7, 3], [8, 4], [9, 5],
                [10, 6], [11, 7], [12, 8], [13, 9], [14, 10], [15, 11], [16, 12], [17, 13], [18,
                    14], [19, 15]]),
        ls.Interaction([[ 0.00308642,  0.        ,  0.        ,  0.        ],
               [ 0.        , -0.00308642,  0.00617284,  0.        ],
               [ 0.        ,  0.00617284, -0.00308642,  0.        ],
               [ 0.        ,  0.        ,  0.        ,  0.00308642]],
            [[0, 18], [1, 19], [2, 0], [3, 1], [4, 2], [5, 3], [6, 4], [7, 5], [8, 6], [9, 7], [10,
                8], [11, 9], [12, 10], [13, 11], [14, 12], [15, 13], [16, 14], [17, 15], [18, 16],
                [19, 17]])
    ]
    # fmt: on
    operator = ls.Operator(basis, interactions)

    e, _ = ls.diagonalize(operator, k=5)
    print(e)


def test_construct_flat_basis():
    basis = ls.SpinBasis(ls.Group([]), number_spins=4, hamming_weight=2)
    flat_basis = ls.FlatSpinBasis(basis)
    assert flat_basis.number_spins == 4
    assert flat_basis.hamming_weight == 2
    assert flat_basis.spin_inversion is None

    basis = ls.SpinBasis(ls.Group([ls.Symmetry([1, 2, 3, 0], sector=1)]), number_spins=4, hamming_weight=2)
    flat_basis = ls.FlatSpinBasis(basis)
    assert flat_basis.number_spins == 4
    assert flat_basis.hamming_weight == 2
    assert flat_basis.spin_inversion is None
    # print(flat_basis.serialize())
    buf = flat_basis.serialize()
    other_basis = ls.FlatSpinBasis.deserialize(buf)
    assert other_basis.number_spins == 4
    assert other_basis.hamming_weight == 2
    assert other_basis.spin_inversion is None
    assert np.all(other_basis.serialize() == buf)

def test_state_info_flat_basis():
    basis = ls.SpinBasis(ls.Group([ls.Symmetry([1, 2, 3, 0], sector=1)]), number_spins=4)
    basis.build()
    flat = ls.FlatSpinBasis(basis)
    full = ls.SpinBasis(ls.Group([]), number_spins=4)
    full.build()
    r1, e1, n1 = basis.batched_state_info(
        np.hstack((full.states.reshape(-1, 1), np.zeros((full.number_states, 7), dtype=np.uint64)))
    )
    r2, e2, n2 = flat.state_info(full.states)
    assert np.all(r1[:, 0] == r2)
    assert np.all(n1 == n2)
    assert np.all(e1 == e2)

    is_r2, n2 = flat.is_representative(full.states)
    assert np.all(basis.states == full.states[is_r2.view(np.bool_)])
    # assert np.app(n1 == n2)
    # if not np.all(e1 == e2):
    #     for i in range(e1.shape[0]):
    #         if e1[i] != e2[i]:
    #             print(i, full.states[i], r1[i], n1[i], ":", e1[i], e2[i])
    #             # assert False




# test_construct_flat_basis()
test_state_info_flat_basis()
# test_index()
# test_4_spins()
# test_construction()
