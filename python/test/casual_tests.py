import glob
import os
import lattice_symmetries as ls
import numpy as np
import scipy.sparse.linalg
#from pytest import approx


def test_index():
    basis = ls.SpinBasis(4)
    basis.build()
    print(basis.index(basis.states))
    print(basis.number_states)
    for i in range(basis.number_states):
        print(i, basis.states[i], basis.state_to_string(i))
    #assert np.array_equal(basis.index(basis.states), basis.states)
    #assert np.array_equal(basis.index(basis.states[2]), 2)

def test_operator_apply():
    basis = ls.SpinBasis(2)
    basis.build()
    expr = ls.Expr("1.0 σᶻ₀ σᶻ₁ + 2.0 σ⁺₀ σ⁻₁ + 2.0 σ⁻₀ σ⁺₁")
    hamiltonian = ls.Operator(basis, expr)
    print(hamiltonian.expression)
    print(hamiltonian.basis.state_to_string(0))
    #len(hamiltonian.apply_off_diag_to_basis_state(basis.states[1]))


def test_prepare_hphi():
    basis = ls.SpinBasis(2)
    expr = ls.Expr("σᶻ₀ σᶻ₁ + 2 σ⁺₀ σ⁻₁ + 2 σ⁻₀ σ⁺₁")
    hamiltonian = ls.Operator(basis, expr)
    hamiltonian.prepare_inputs_for_hphi("/tmp/lattice-symmetries-python/hphi")


def main():
    #test_index()
    test_operator_apply()

main()
