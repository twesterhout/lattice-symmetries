import numpy as np, lattice_symmetries as ls

def test_SpinBasis():
    basis = ls.SpinBasis(3)
    basis.build()
    np.testing.assert_equal(basis.index(basis.states), np.arange(2**3))
    assert basis.number_states == 2**3
    assert [basis.state_to_string(basis.states[i]) for i in range(basis.number_states)] == [
        "|000⟩",
        "|001⟩",
        "|010⟩",
        "|011⟩",
        "|100⟩",
        "|101⟩",
        "|110⟩",
        "|111⟩",
    ]

    # basis = ls.SpinBasis(3, hamming_weight=2)  # We want the subspace with only 2 spins up
    # basis.build()
    # assert [basis.state_to_string(basis.states[i]) for i in range(basis.number_states)] == [
    #     "|011⟩",
    #     "|101⟩",
    #     "|110⟩",
    # ]

    # basis = ls.SpinBasis(4, hamming_weight=2, spin_inversion=-1)
    # basis.build()
    # assert [basis.state_to_string(basis.states[i]) for i in range(basis.number_states)] == [
    #     "|0011⟩",
    #     "|0101⟩",
    #     "|0110⟩",
    # ]
