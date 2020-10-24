#include "lattice_symmetries/lattice_symmetries.h"
#include <catch2/catch.hpp>
#include <complex>
#include <cstdio>
#include <numeric>

TEST_CASE("constructs symmetries", "[api]")
{
    {
        unsigned const permutation[] = {1, 2, 3, 4, 5, 6, 0};
        ls_symmetry*   self          = nullptr;
        ls_error_code  status =
            ls_create_symmetry(&self, std::size(permutation), permutation, false, 0);
        REQUIRE(status == LS_SUCCESS);
        REQUIRE(self != nullptr);
        REQUIRE(ls_get_sector(self) == 0);
        REQUIRE(ls_get_flip(self) == false);
        REQUIRE(ls_get_periodicity(self) == 7);
        REQUIRE(ls_get_phase(self) == 0.0);
        std::complex<double> eigenvalue;
        ls_get_eigenvalue(self, &eigenvalue);
        REQUIRE(eigenvalue == 1.0);
        ls_destroy_symmetry(self);
    }

    {
        unsigned const permutation[] = {1, 2, 3, 0};
        ls_symmetry*   self          = nullptr;
        ls_error_code  status =
            ls_create_symmetry(&self, std::size(permutation), permutation, true, 3);
        REQUIRE(status == LS_SUCCESS);
        REQUIRE(self != nullptr);
        REQUIRE(ls_get_sector(self) == 3);
        REQUIRE(ls_get_flip(self) == true);
        REQUIRE(ls_get_periodicity(self) == 4);
        REQUIRE(ls_get_phase(self) == 0.75);
        ls_destroy_symmetry(self);
    }

    {
        unsigned const permutation[] = {1, 2, 0};
        ls_symmetry*   self          = nullptr;
        ls_error_code  status =
            ls_create_symmetry(&self, std::size(permutation), permutation, true, 3);
        REQUIRE(status == LS_SUCCESS);
        REQUIRE(self != nullptr);
        REQUIRE(ls_get_sector(self) == 3);
        REQUIRE(ls_get_flip(self) == true);
        REQUIRE(ls_get_periodicity(self) == 6);
        REQUIRE(ls_get_phase(self) == 0.5);
        ls_destroy_symmetry(self);
    }

    {
        unsigned      permutation[] = {1, 2, 2, 0};
        ls_symmetry*  self          = nullptr;
        ls_error_code status =
            ls_create_symmetry(&self, std::size(permutation), permutation, true, 3);
        REQUIRE(status == LS_INVALID_PERMUTATION);
    }

    {
        std::vector<unsigned> permutation(1000);
        std::iota(std::begin(permutation), std::end(permutation), 0U);
        ls_symmetry*  self = nullptr;
        ls_error_code status =
            ls_create_symmetry(&self, std::size(permutation), permutation.data(), false, 0);
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
            ls_create_symmetry(&symmetry, std::size(permutation), permutation, false, 0);
        REQUIRE(status == LS_SUCCESS);

        ls_symmetry const* generators[] = {symmetry};
        status = ls_create_group(&group, std::size(generators), generators);
        REQUIRE(status == LS_SUCCESS);
        REQUIRE(ls_get_group_size(group) == 4);
        ls_destroy_group(group);
        ls_destroy_symmetry(symmetry);
    }
}

TEST_CASE("constructs basis", "[api]")
{
    {
        unsigned const permutation[] = {1, 2, 3, 0};
        ls_symmetry*   symmetry      = nullptr;
        ls_error_code  status =
            ls_create_symmetry(&symmetry, std::size(permutation), permutation, false, 0);
        REQUIRE(status == LS_SUCCESS);

        ls_symmetry const* generators[] = {symmetry};
        ls_group*          group        = nullptr;
        status = ls_create_group(&group, std::size(generators), generators);
        REQUIRE(status == LS_SUCCESS);
        ls_destroy_symmetry(symmetry);

        ls_spin_basis* basis = nullptr;
        status               = ls_create_spin_basis(&basis, group, 4, -1);
        REQUIRE(status == LS_SUCCESS);

        ls_destroy_group(group);

        status = ls_build(basis);
        REQUIRE(status == LS_SUCCESS);

        // uint64_t count;
        // status = ls_get_number_states(basis, &count);
        // REQUIRE(status == LS_SUCCESS);
        // printf("%zu\n", count);

        // ls_states* states = nullptr;
        // status            = ls_get_states(&states, basis);
        // REQUIRE(status == LS_SUCCESS);

        // auto const begin = ls_states_get_data(states);
        // auto const end   = begin + count;
        // for (auto i = begin; i != end; ++i) {
        //     printf("%zu, ", *i);
        // }
        // printf("\n");

        // ls_destroy_states(states);
        ls_destroy_spin_basis(basis);
    }

    {
        ls_group*     group  = nullptr;
        ls_error_code status = ls_create_group(&group, 0, nullptr);
        REQUIRE(status == LS_SUCCESS);

        ls_spin_basis* basis = nullptr;
        status               = ls_create_spin_basis(&basis, group, 4, 2);
        REQUIRE(status == LS_SUCCESS);

        ls_destroy_group(group);

        status = ls_build(basis);
        REQUIRE(status == LS_SUCCESS);

        uint64_t count;
        status = ls_get_number_states(basis, &count);
        REQUIRE(status == LS_SUCCESS);
        REQUIRE(count == 6);

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
        status = ls_create_symmetry(&T_x_symmetry, std::size(T_x), T_x, false, 0);
        REQUIRE(status == LS_SUCCESS);

        ls_symmetry* T_y_symmetry = nullptr;
        status = ls_create_symmetry(&T_y_symmetry, std::size(T_y), T_y, false, 0);
        REQUIRE(status == LS_SUCCESS);

        ls_symmetry* P_x_symmetry = nullptr;
        status = ls_create_symmetry(&P_x_symmetry, std::size(P_x), P_x, false, 0);
        REQUIRE(status == LS_SUCCESS);

        ls_symmetry* P_y_symmetry = nullptr;
        status = ls_create_symmetry(&P_y_symmetry, std::size(P_y), P_y, false, 0);
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
        status               = ls_create_spin_basis(&basis, group, 24, 12);
        REQUIRE(status == LS_SUCCESS);

        ls_destroy_group(group);

        status = ls_build(basis);
        REQUIRE(status == LS_SUCCESS);

        uint64_t count;
        status = ls_get_number_states(basis, &count);
        REQUIRE(status == LS_SUCCESS);
        REQUIRE(count == 29606U);

        ls_states* states = nullptr;
        status            = ls_get_states(&states, basis);
        REQUIRE(status == LS_SUCCESS);

        auto const begin = ls_states_get_data(states);
        auto const end   = begin + count;
        for (auto i = 0UL; i < count; ++i) {
            uint64_t index;
            status = ls_get_index(basis, &begin[i], &index);
            REQUIRE(status == LS_SUCCESS);
            REQUIRE(index == i);
        }

        ls_destroy_states(states);
        ls_destroy_spin_basis(basis);
    }
}
