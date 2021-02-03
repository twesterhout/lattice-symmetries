#include <lattice_symmetries/lattice_symmetries.h>
#include <stdlib.h>

#define sizeof_array(xs) (sizeof(xs) / sizeof((xs)[0]))

static ls_spin_basis* make_basis_6x4(void)
{
    //   0,  1,  2,  3,  4,  5,
    //   6,  7,  8,  9, 10, 11,
    //  12, 13, 14, 15, 16, 17,
    //  18, 19, 20, 21, 22, 23,
    static unsigned const width       = 6U;
    static unsigned const height      = 4U;
    static unsigned const system_size = width * height;

    unsigned identity[height][width];
    unsigned offset = 0;
    for (unsigned i = 0; i < height; ++i) {
        for (unsigned j = 0; j < width; ++j, ++offset) {
            identity[i][j] = offset;
        }
    }

    unsigned permutation[system_size];
    offset = 0;
    for (unsigned i = 0; i < height; ++i) {
        for (unsigned j = 0; j < width; ++j, ++offset) {
            permutation[offset] = identity[i][(j + 1) % width];
        }
    }
    ls_symmetry*  T_x    = NULL;
    ls_error_code status = ls_create_symmetry(&T_x, system_size, permutation, 0);
    LATTICE_SYMMETRIES_CHECK(status == LS_SUCCESS, "");

    ls_symmetry const* generators[] = {T_x};
    ls_group*          group        = NULL;
    status                          = ls_create_group(&group, 1, generators);
    LATTICE_SYMMETRIES_CHECK(status == LS_SUCCESS, "");

    ls_spin_basis* basis = NULL;
    status               = ls_create_spin_basis(&basis, group, system_size, system_size / 2, 1);
    LATTICE_SYMMETRIES_CHECK(status == LS_SUCCESS, "");

    ls_destroy_group(group);
    ls_destroy_symmetry(T_x);
    return basis;
}

static ls_operator* make_operator_6x4(ls_spin_basis const* basis)
{
    static unsigned const width       = 6U;
    static unsigned const height      = 4U;
    static unsigned const system_size = width * height;

    unsigned identity[height][width];
    unsigned offset = 0;
    for (unsigned i = 0; i < height; ++i) {
        for (unsigned j = 0; j < width; ++j, ++offset) {
            identity[i][j] = offset;
        }
    }

    // clang-format off
    _Complex double matrix[4][4] = {
        {1.0,  0.0,  0.0, 0.0},
        {0.0, -1.0,  2.0, 0.0},
        {0.0,  2.0, -1.0, 0.0},
        {0.0,  0.0,  0.0, 1.0}};
    // clang-format on
    uint16_t edges[2 * system_size][2];
    offset = 0;
    for (unsigned i = 0; i < height; ++i) {
        for (unsigned j = 0; j < width; ++j) {
            edges[offset][0] = identity[i][j];
            edges[offset][1] = identity[(i + 1) % height][j];
            ++offset;
            edges[offset][0] = identity[i][j];
            edges[offset][1] = identity[i][(j + 1) % width];
            ++offset;
        }
    }

    ls_interaction* interaction = NULL;
    ls_error_code status = ls_create_interaction2(&interaction, matrix, sizeof_array(edges), edges);
    LATTICE_SYMMETRIES_CHECK(status == LS_SUCCESS, "");

    ls_operator*          op      = NULL;
    ls_interaction const* terms[] = {interaction};
    status                        = ls_create_operator(&op, basis, 1, terms);
    LATTICE_SYMMETRIES_CHECK(status == LS_SUCCESS, "");

    ls_destroy_interaction(interaction);
    return op;
}

void do_measure(ls_operator const* op)
{

}

int main(int argc, char** argv)
{
    ls_enable_logging();
    ls_spin_basis* basis = make_basis_6x4();
    ls_build(basis);
    ls_operator* op = make_operator_6x4(basis);

    ls_destroy_operator(op);
    ls_destroy_spin_basis(basis);
    return 0;
}
