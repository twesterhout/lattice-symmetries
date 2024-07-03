import itertools

import lattice_symmetries as ls
import numpy as np
import numpy.typing as npt
import scipy


def test_reduced_density_matrix_hamming():
    def rdms_reference(
        ground_state: npt.NDArray[np.float64], sites: list[int], n_spins: int
    ) -> npt.NDArray[np.float64]:
        ground_state_reshaped = ground_state.reshape((2,) * n_spins).transpose(
            range(n_spins - 1, -1, -1)
        )
        idxs = tuple(sorted(set(range(n_spins)) - set(sites)))
        ground_state_transposed = ground_state_reshaped.transpose(tuple(reversed(sites)) + idxs)
        ground_state_transposed = ground_state_transposed.reshape(
            (2 ** len(sites), 2 ** (n_spins - len(sites)))
        )
        return ground_state_transposed @ ground_state_transposed.T

    expr_str = "σᶻ₀ σᶻ₁ + σᶻ₀ σᶻ₃ + σᶻ₀ σᶻ₄ + σᶻ₀ σᶻ₁₂ + 2.0 σ⁺₀ σ⁻₁ + 2.0 σ⁺₀ σ⁻₃ + 2.0 σ⁺₀ σ⁻₄ + 2.0 σ⁺₀ σ⁻₁₂ + 2.0 σ⁻₀ σ⁺₁ + 2.0 σ⁻₀ σ⁺₃ + 2.0 σ⁻₀ σ⁺₄ + 2.0 σ⁻₀ σ⁺₁₂ + σᶻ₁ σᶻ₂ + σᶻ₁ σᶻ₅ + σᶻ₁ σᶻ₁₃ + 2.0 σ⁺₁ σ⁻₂ + 2.0 σ⁺₁ σ⁻₅ + 2.0 σ⁺₁ σ⁻₁₃ + 2.0 σ⁻₁ σ⁺₂ + 2.0 σ⁻₁ σ⁺₅ + 2.0 σ⁻₁ σ⁺₁₃ + σᶻ₂ σᶻ₃ + σᶻ₂ σᶻ₆ + σᶻ₂ σᶻ₁₄ + 2.0 σ⁺₂ σ⁻₃ + 2.0 σ⁺₂ σ⁻₆ + 2.0 σ⁺₂ σ⁻₁₄ + 2.0 σ⁻₂ σ⁺₃ + 2.0 σ⁻₂ σ⁺₆ + 2.0 σ⁻₂ σ⁺₁₄ + σᶻ₃ σᶻ₇ + σᶻ₃ σᶻ₁₅ + 2.0 σ⁺₃ σ⁻₇ + 2.0 σ⁺₃ σ⁻₁₅ + 2.0 σ⁻₃ σ⁺₇ + 2.0 σ⁻₃ σ⁺₁₅ + σᶻ₄ σᶻ₅ + σᶻ₄ σᶻ₇ + σᶻ₄ σᶻ₈ + 2.0 σ⁺₄ σ⁻₅ + 2.0 σ⁺₄ σ⁻₇ + 2.0 σ⁺₄ σ⁻₈ + 2.0 σ⁻₄ σ⁺₅ + 2.0 σ⁻₄ σ⁺₇ + 2.0 σ⁻₄ σ⁺₈ + σᶻ₅ σᶻ₆ + σᶻ₅ σᶻ₉ + 2.0 σ⁺₅ σ⁻₆ + 2.0 σ⁺₅ σ⁻₉ + 2.0 σ⁻₅ σ⁺₆ + 2.0 σ⁻₅ σ⁺₉ + σᶻ₆ σᶻ₇ + σᶻ₆ σᶻ₁₀ + 2.0 σ⁺₆ σ⁻₇ + 2.0 σ⁺₆ σ⁻₁₀ + 2.0 σ⁻₆ σ⁺₇ + 2.0 σ⁻₆ σ⁺₁₀ + σᶻ₇ σᶻ₁₁ + 2.0 σ⁺₇ σ⁻₁₁ + 2.0 σ⁻₇ σ⁺₁₁ + σᶻ₈ σᶻ₉ + σᶻ₈ σᶻ₁₁ + σᶻ₈ σᶻ₁₂ + 2.0 σ⁺₈ σ⁻₉ + 2.0 σ⁺₈ σ⁻₁₁ + 2.0 σ⁺₈ σ⁻₁₂ + 2.0 σ⁻₈ σ⁺₉ + 2.0 σ⁻₈ σ⁺₁₁ + 2.0 σ⁻₈ σ⁺₁₂ + σᶻ₉ σᶻ₁₀ + σᶻ₉ σᶻ₁₃ + 2.0 σ⁺₉ σ⁻₁₀ + 2.0 σ⁺₉ σ⁻₁₃ + 2.0 σ⁻₉ σ⁺₁₀ + 2.0 σ⁻₉ σ⁺₁₃ + σᶻ₁₀ σᶻ₁₁ + σᶻ₁₀ σᶻ₁₄ + 2.0 σ⁺₁₀ σ⁻₁₁ + 2.0 σ⁺₁₀ σ⁻₁₄ + 2.0 σ⁻₁₀ σ⁺₁₁ + 2.0 σ⁻₁₀ σ⁺₁₄ + σᶻ₁₁ σᶻ₁₅ + 2.0 σ⁺₁₁ σ⁻₁₅ + 2.0 σ⁻₁₁ σ⁺₁₅ + σᶻ₁₂ σᶻ₁₃ + σᶻ₁₂ σᶻ₁₅ + 2.0 σ⁺₁₂ σ⁻₁₃ + 2.0 σ⁺₁₂ σ⁻₁₅ + 2.0 σ⁻₁₂ σ⁺₁₃ + 2.0 σ⁻₁₂ σ⁺₁₅ + σᶻ₁₃ σᶻ₁₄ + 2.0 σ⁺₁₃ σ⁻₁₄ + 2.0 σ⁻₁₃ σ⁺₁₄ + σᶻ₁₄ σᶻ₁₅ + 2.0 σ⁺₁₄ σ⁻₁₅ + 2.0 σ⁻₁₄ σ⁺₁₅"
    number_sites = 16
    basis = ls.SpinBasis(
        number_spins=number_sites,
        hamming_weight=number_sites // 2,
    )
    basis_no_hamming = ls.SpinBasis(number_spins=number_sites, hamming_weight=None)
    basis.build()
    basis_no_hamming.build()
    expression = ls.Expr(expr_str, particle="spin-1/2")
    hamiltonian = ls.Operator(expression, basis)
    hamiltonian_no_hamming = ls.Operator(expression, basis_no_hamming)
    ground_state = scipy.sparse.linalg.eigsh(hamiltonian, k=1, which="SA")[1].ravel()
    ground_state_no_hamming = scipy.sparse.linalg.eigsh(hamiltonian_no_hamming, k=1, which="SA")[
        1
    ].ravel()

    for length in [1, 2, 3]:
        for sites in itertools.combinations(range(basis.number_sites), length):
            reference = rdms_reference(ground_state_no_hamming, list(sites), basis.number_sites)
            rdm = ls.rdms.reduced_density_matrix(
                ground_state.astype(np.complex128), basis, list(sites), return_states=False
            )
            assert np.allclose(reference, rdm)
