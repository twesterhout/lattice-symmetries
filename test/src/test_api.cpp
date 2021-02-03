#include "cpu/search_sorted.hpp"
#include "lattice_symmetries/lattice_symmetries.h"
#include <catch2/catch.hpp>
#include <complex>
#include <cstdio>
#include <numeric>

TEST_CASE("toggles debug logging", "[api]")
{
    REQUIRE(ls_is_logging_enabled() == false);
    ls_enable_logging();
    REQUIRE(ls_is_logging_enabled() == true);
    ls_disable_logging();
    REQUIRE(ls_is_logging_enabled() == false);
}

TEST_CASE("constructs symmetries", "[api]")
{
    {
        unsigned const permutation[] = {1, 2, 3, 4, 5, 6, 0};
        ls_symmetry*   self          = nullptr;
        ls_error_code  status = ls_create_symmetry(&self, std::size(permutation), permutation, 0);
        REQUIRE(status == LS_SUCCESS);
        REQUIRE(self != nullptr);
        REQUIRE(ls_get_sector(self) == 0);
        REQUIRE(ls_get_periodicity(self) == 7);
        REQUIRE(ls_get_phase(self) == 0.0);
        REQUIRE(ls_symmetry_get_number_spins(self) == 7);
        std::complex<double> eigenvalue;
        ls_get_eigenvalue(self, &eigenvalue);
        REQUIRE(eigenvalue == 1.0);
        ls_destroy_symmetry(self);
    }

    {
        unsigned const permutation[] = {1, 2, 3, 0};
        ls_symmetry*   self          = nullptr;
        ls_error_code  status = ls_create_symmetry(&self, std::size(permutation), permutation, 3);
        REQUIRE(status == LS_SUCCESS);
        REQUIRE(self != nullptr);
        REQUIRE(ls_get_sector(self) == 3);
        REQUIRE(ls_get_periodicity(self) == 4);
        REQUIRE(ls_get_phase(self) == 0.75);
        REQUIRE(ls_symmetry_get_number_spins(self) == 4);
        ls_destroy_symmetry(self);
    }

    {
        unsigned const permutation[] = {3, 2, 1, 0};
        ls_symmetry*   self          = nullptr;
        ls_error_code  status = ls_create_symmetry(&self, std::size(permutation), permutation, 1);
        REQUIRE(status == LS_SUCCESS);
        REQUIRE(self != nullptr);
        REQUIRE(ls_get_sector(self) == 1);
        REQUIRE(ls_get_periodicity(self) == 2);
        REQUIRE(ls_get_phase(self) == 0.5);
        REQUIRE(ls_symmetry_get_number_spins(self) == 4);
        std::complex<double> eigenvalue;
        ls_get_eigenvalue(self, &eigenvalue);
        REQUIRE(eigenvalue == -1.0);
        ls_destroy_symmetry(self);
    }

    {
        unsigned      permutation[] = {1, 2, 2, 0};
        ls_symmetry*  self          = nullptr;
        ls_error_code status = ls_create_symmetry(&self, std::size(permutation), permutation, 3);
        REQUIRE(status == LS_INVALID_PERMUTATION);
    }

    {
        std::vector<unsigned> permutation(1000);
        std::iota(std::begin(permutation), std::end(permutation), 0U);
        ls_symmetry*  self = nullptr;
        ls_error_code status =
            ls_create_symmetry(&self, std::size(permutation), permutation.data(), 0);
        REQUIRE(status == LS_PERMUTATION_TOO_LONG);
    }
}

TEST_CASE("constructs symmetry groups", "[api]")
{
    {
        unsigned const permutation[] = {1, 2, 3, 0};
        ls_symmetry*   symmetry      = nullptr;
        ls_group*      group         = nullptr;
        ls_error_code  status =
            ls_create_symmetry(&symmetry, std::size(permutation), permutation, 0);
        REQUIRE(status == LS_SUCCESS);

        ls_symmetry const* generators[] = {symmetry};
        status = ls_create_group(&group, std::size(generators), generators);
        REQUIRE(status == LS_SUCCESS);
        REQUIRE(ls_get_group_size(group) == 4);
        ls_destroy_group(group);
        ls_destroy_symmetry(symmetry);
    }

    {
        unsigned              number_spins = 70U;
        std::vector<unsigned> T;
        std::vector<unsigned> P;
        for (auto i = 0U; i < number_spins; ++i) {
            T.push_back((i + 1U) % number_spins);
            P.push_back(number_spins - 1U - i);
        }

        ls_symmetry*  translation = nullptr;
        ls_symmetry*  parity      = nullptr;
        ls_group*     group       = nullptr;
        ls_error_code status =
            ls_create_symmetry(&translation, T.size(), T.data(), number_spins / 2);
        REQUIRE(status == LS_SUCCESS);
        status = ls_create_symmetry(&parity, P.size(), P.data(), 1);
        REQUIRE(status == LS_SUCCESS);

        ls_symmetry const* generators[] = {translation, parity};
        status = ls_create_group(&group, std::size(generators), generators);
        REQUIRE(status == LS_SUCCESS);
        REQUIRE(ls_get_group_size(group) == 2 * number_spins);
        ls_destroy_group(group);
        ls_destroy_symmetry(parity);
        ls_destroy_symmetry(translation);
    }

    {
        ls_group*     group  = nullptr;
        ls_error_code status = ls_create_group(&group, 0, nullptr);
        REQUIRE(status == LS_SUCCESS);
        REQUIRE(group != nullptr);
        REQUIRE(ls_get_group_size(group) == 0);
        ls_destroy_group(group);
    }

    {
        ls_group*     group  = nullptr;
        ls_error_code status = ls_create_trivial_group(&group, 25U);
        REQUIRE(status == LS_SUCCESS);
        REQUIRE(group != nullptr);
        REQUIRE(ls_get_group_size(group) == 1);
        ls_destroy_group(group);
    }
}

TEST_CASE("constructs basis", "[api]")
{
    {
        unsigned const permutation[] = {1, 2, 3, 0};
        ls_symmetry*   symmetry      = nullptr;
        ls_error_code  status =
            ls_create_symmetry(&symmetry, std::size(permutation), permutation, 0);
        REQUIRE(status == LS_SUCCESS);

        ls_symmetry const* generators[] = {symmetry};
        ls_group*          group        = nullptr;
        status = ls_create_group(&group, std::size(generators), generators);
        REQUIRE(status == LS_SUCCESS);
        ls_destroy_symmetry(symmetry);

        ls_spin_basis* basis = nullptr;
        status               = ls_create_spin_basis(&basis, group, 4, -1, 0);
        REQUIRE(status == LS_SUCCESS);
        ls_destroy_group(group);
        REQUIRE(ls_get_number_spins(basis) == 4);
        REQUIRE(ls_get_hamming_weight(basis) == -1);
        REQUIRE(ls_get_spin_inversion(basis) == 0);
        REQUIRE(ls_has_symmetries(basis) == true);

        status = ls_build(basis);
        REQUIRE(status == LS_SUCCESS);

        uint64_t count;
        status = ls_get_number_states(basis, &count);
        REQUIRE(status == LS_SUCCESS);
        REQUIRE(count == 6); // 1 + 1 + 1 + 2 + 1

        ls_states* states = nullptr;
        status            = ls_get_states(&states, basis);
        REQUIRE(status == LS_SUCCESS);

        auto const begin = ls_states_get_data(states);
        auto const end   = begin + count;
        for (auto i = 0U; i < count; ++i) {
            uint64_t index;
            status = ls_get_index(basis, begin[i], &index);
            REQUIRE(status == LS_SUCCESS);
            REQUIRE(index == i);
        }

        ls_destroy_states(states);
        ls_destroy_spin_basis(basis);
    }

    {
        ls_group*     group  = nullptr;
        ls_error_code status = ls_create_group(&group, 0, nullptr);
        REQUIRE(status == LS_SUCCESS);

        ls_spin_basis* basis = nullptr;
        status               = ls_create_spin_basis(&basis, group, 4, 2, 1);
        REQUIRE(status == LS_SUCCESS);

        ls_destroy_group(group);

        status = ls_build(basis);
        REQUIRE(status == LS_SUCCESS);

        uint64_t count;
        status = ls_get_number_states(basis, &count);
        REQUIRE(status == LS_SUCCESS);
        REQUIRE(count == 3);

        ls_destroy_spin_basis(basis);
    }
}

TEST_CASE("finds correct states", "[api]")
{
    {
        ls_error_code status;
        // clang-format off
        unsigned const sites[] =
            { 0,  1,  2,  3,  4,  5,
              6,  7,  8,  9, 10, 11,
             12, 13, 14, 15, 16, 17,
             18, 19, 20, 21, 22, 23};
        unsigned const T_x[] =
            { 1,  2,  3,  4,  5,  0,
              7,  8,  9, 10, 11,  6,
             13, 14, 15, 16, 17, 12,
             19, 20, 21, 22, 23, 18};
        unsigned const T_y[] =
            { 6,  7,  8,  9, 10, 11,
             12, 13, 14, 15, 16, 17,
             18, 19, 20, 21, 22, 23,
              0,  1,  2,  3,  4,  5};
        unsigned const P_x[] =
            {  5,  4,  3,  2,  1,  0,
              11, 10,  9,  8,  7,  6,
              17, 16, 15, 14, 13, 12,
              23, 22, 21, 20, 19, 18};
        unsigned const P_y[] =
            {18, 19, 20, 21, 22, 23,
             12, 13, 14, 15, 16, 17,
              6,  7,  8,  9, 10, 11,
              0,  1,  2,  3,  4,  5};
        // clang-format on
        ls_symmetry* T_x_symmetry = nullptr;
        status                    = ls_create_symmetry(&T_x_symmetry, std::size(T_x), T_x, 0);
        REQUIRE(status == LS_SUCCESS);

        ls_symmetry* T_y_symmetry = nullptr;
        status                    = ls_create_symmetry(&T_y_symmetry, std::size(T_y), T_y, 0);
        REQUIRE(status == LS_SUCCESS);

        ls_symmetry* P_x_symmetry = nullptr;
        status                    = ls_create_symmetry(&P_x_symmetry, std::size(P_x), P_x, 0);
        REQUIRE(status == LS_SUCCESS);

        ls_symmetry* P_y_symmetry = nullptr;
        status                    = ls_create_symmetry(&P_y_symmetry, std::size(P_y), P_y, 0);
        REQUIRE(status == LS_SUCCESS);

        ls_symmetry const* generators[] = {T_x_symmetry, T_y_symmetry, P_x_symmetry, P_y_symmetry};
        ls_group*          group        = nullptr;
        status = ls_create_group(&group, std::size(generators), generators);
        REQUIRE(status == LS_SUCCESS);
        ls_destroy_symmetry(T_x_symmetry);
        ls_destroy_symmetry(T_y_symmetry);
        ls_destroy_symmetry(P_x_symmetry);
        ls_destroy_symmetry(P_y_symmetry);

        ls_spin_basis* basis = nullptr;
        status               = ls_create_spin_basis(&basis, group, 24, 12, 1);
        REQUIRE(status == LS_SUCCESS);
        ls_destroy_group(group);

        status = ls_build(basis);
        REQUIRE(status == LS_SUCCESS);

        uint64_t count;
        status = ls_get_number_states(basis, &count);
        REQUIRE(status == LS_SUCCESS);
        REQUIRE(count == 15578U); // Computed using QuSpin with the same symmetries

        ls_states* states = nullptr;
        status            = ls_get_states(&states, basis);
        REQUIRE(status == LS_SUCCESS);

        auto const begin = ls_states_get_data(states);
        auto const end   = begin + count;
        for (auto i = 0UL; i < count; ++i) {
            uint64_t index;
            status = ls_get_index(basis, begin[i], &index);
            REQUIRE(status == LS_SUCCESS);
            REQUIRE(index == i);

            uint64_t const other_index = lattice_symmetries::search_sorted(begin, count, begin[i]);
            REQUIRE(other_index == i);
        }

        ls_destroy_states(states);
        ls_destroy_spin_basis(basis);
    }

    {
        ls_enable_logging();
        LATTICE_SYMMETRIES_LOG_DEBUG("avx2: %i\navx:  %i\nsse4: %i\n", ls_has_avx2(), ls_has_avx(),
                                     ls_has_sse4());
    }
}
