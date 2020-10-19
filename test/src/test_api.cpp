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
