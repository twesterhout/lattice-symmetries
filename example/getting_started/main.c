#include <lattice_symmetries/lattice_symmetries.h>

#include <complex.h>
#include <stdio.h>
#include <stdlib.h>

void print_error_message_and_exit(ls_error_code const status)
{
    char const* msg = ls_error_to_string(status);
    fprintf(stderr, "Error (%i): %s\n", status, msg);
    ls_destroy_string(msg);
    exit(1);
}

int main(void)
{
    unsigned const number_spins   = 10U;
    int const      hamming_weight = 5;

    unsigned* permutation = malloc(number_spins * sizeof(unsigned));
    if (permutation == NULL) {
        fprintf(stderr, "Error: failed to allocate memory\n");
        return 1;
    }

    // Momentum in x direction with eigenvalue π
    for (unsigned i = 0U; i < number_spins; ++i) {
        permutation[i] = (i + 1U) % number_spins;
    }
    ls_symmetry*  momentum;
    ls_error_code status = ls_create_symmetry(&momentum, number_spins, permutation, /*flip=*/false,
                                              /*sector=*/number_spins / 2);
    if (status != LS_SUCCESS) { goto fail1; }

    // Parity with eigenvalue π
    for (unsigned i = 0U; i < number_spins; ++i) {
        permutation[i] = number_spins - 1U - i;
    }
    ls_symmetry* parity;
    status = ls_create_symmetry(&parity, number_spins, permutation, /*flip=*/false, /*sector=*/1);
    if (status != LS_SUCCESS) { goto fail2; }

    // Global spin inversion with eigenvalue π
    for (unsigned i = 0U; i < number_spins; ++i) {
        permutation[i] = i;
    }
    ls_symmetry* inversion;
    status = ls_create_symmetry(&inversion, number_spins, permutation, /*flip=*/true, /*sector=*/1);
    if (status != LS_SUCCESS) { goto fail3; }

    // Constructing the group
    ls_symmetry const* generators[] = {momentum, parity, inversion};
    ls_group*          group;
    status = ls_create_group(&group, 3, generators);
    if (status != LS_SUCCESS) { goto fail4; }
    printf("Symmetry group contains %u elements\n", ls_get_group_size(group));

    // Constructing the basis
    ls_spin_basis* basis;
    status = ls_create_spin_basis(&basis, group, number_spins, hamming_weight);
    if (status != LS_SUCCESS) { goto fail5; }
    status = ls_build(basis);
    if (status != LS_SUCCESS) { goto fail6; }

    uint64_t number_states;
    status = ls_get_number_states(basis, &number_states);
    if (status != LS_SUCCESS) { goto fail6; }
    printf("Hilbert space dimension is %zu\n", number_states);

    // Heisenberg Hamiltonian
    // It is easier to compute kronecker product of Pauli matrices in Mathematica or Python and just
    // hardcode the result here.
    // clang-format off
    _Complex double const matrix[] = {1.0,  0.0,  0.0, 0.0,
                                      0.0, -1.0,  2.0, 0.0,
                                      0.0,  2.0, -1.0, 0.0,
                                      0.0,  0.0,  0.0, 1.0};
    // clang-format on
    uint16_t(*edges)[2] = malloc(number_spins * sizeof(uint16_t[2]));
    if (edges == NULL) {
        status = LS_OUT_OF_MEMORY;
        goto fail6;
    }
    for (unsigned i = 0U; i < number_spins; ++i) {
        edges[i][0] = i;
        edges[i][1] = (i + 1) % number_spins;
    }
    ls_interaction* term;
    status = ls_create_interaction2(&term, matrix, number_spins, edges);
    if (status != LS_SUCCESS) { goto fail7; }
    ls_operator* hamiltonian;
    status = ls_create_operator(&hamiltonian, basis, 1, (ls_interaction const* const*)&term);
    if (status != LS_SUCCESS) { goto fail8; }

    // Unfortunately, we cannot diagonalize the operator without relying on external libraries. We
    // can, however, compute its expectation value
    if (number_states == 13) {
        double const ground_state[13] = {
            1.522647886368464539e-03,  -1.295335714623147005e-02, 4.027096491793154959e-02,
            -1.084815641037164685e-01, 1.919478508837144104e-01,  -5.008341285296182693e-02,
            8.000569172532195905e-02,  -2.027799615731489813e-01, -1.531995338521515981e-01,
            5.985187387435567663e-01,  -4.627911117494983295e-01, -2.884870225547441214e-01,
            4.695442358221240675e-01};

        _Complex double energy;
        status = ls_operator_expectation(hamiltonian, LS_FLOAT64, number_states, 1, ground_state, 1,
                                         &energy);
        if (status != LS_SUCCESS) { goto fail9; }
        printf("Ground state energy is %f\n", creal(energy));
    }

    // Cleaning up
fail9:
    ls_destroy_operator(hamiltonian);
fail8:
    ls_destroy_interaction(term);
fail7:
    free(edges);
fail6:
    ls_destroy_spin_basis(basis);
fail5:
    ls_destroy_group(group);
fail4:
    ls_destroy_symmetry(inversion);
fail3:
    ls_destroy_symmetry(parity);
fail2:
    ls_destroy_symmetry(momentum);
fail1:
    free(permutation);
    if (status != LS_SUCCESS) { print_error_message_and_exit(status); }
    return 0;
}
