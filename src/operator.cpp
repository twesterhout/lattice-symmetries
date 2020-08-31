#include "operator.hpp"
#include "basis.hpp"
#include "bits.hpp"
#include "error_handling.hpp"
#include "lattice_symmetries/lattice_symmetries.h"
#include <omp.h>
#include <algorithm>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <span.hpp>
#include <variant>
#include <vector>

namespace lattice_symmetries {

namespace {
    template <unsigned N> auto transpose(std::complex<double> (&matrix)[N][N]) noexcept -> void
    {
        for (auto i = 0U; i < N; ++i) {
            for (auto j = 0U; j < i; ++j) {
                std::swap(matrix[i][j], matrix[j][i]);
            }
        }
    }
} // namespace

template <unsigned NumberSpins> struct interaction_t {
    struct free_fn_t {
        auto operator()(void* p) const noexcept -> void { std::free(p); }
    };

    struct alignas(64) matrix_t {
        std::complex<double> payload[NumberSpins][NumberSpins];

        explicit matrix_t(std::complex<double> const* data) noexcept
        {
            std::memcpy(payload, data, NumberSpins * NumberSpins * sizeof(std::complex<double>));
        }
    };
    using sites_t = std::vector<std::array<uint16_t, NumberSpins>>;

    std::unique_ptr<matrix_t> matrix;
    sites_t                   sites;

    interaction_t(std::complex<double> const*                        _matrix,
                  tcb::span<std::array<uint16_t, NumberSpins> const> _sites)
        : matrix{std::make_unique<matrix_t>(_matrix)}, sites{std::begin(_sites), std::end(_sites)}
    {
        transpose(matrix->payload);
    }
};

namespace {
    template <class Bits, std::size_t N>
    constexpr auto gather_bits(Bits const& bits, std::array<uint16_t, N> const& indices) noexcept
        -> unsigned
    {
        // ============================= IMPORTANT ==============================
        // The order is REALLY important here. This is done to adhere to the
        // definition of cronecker product, i.e. that
        //
        //     kron(A, B) =  A00 B     A01 B     A02 B  ...
        //                   A10 B     A11 B     A12 B  ...
        //                    .
        //                    .
        //                    .
        //
        // In other words, if you change it to
        //    r |= test_bit(bits, i);
        //    r <<= 1U;
        // shit will break in really difficult to track ways...
        auto r = 0U;
        for (auto const i : indices) {
            LATTICE_SYMMETRIES_ASSERT(i < 64, "index out of bounds");
            r <<= 1U;
            r |= test_bit(bits, i);
        }
        return r;
    }

    template <class Bits, std::size_t N>
    constexpr auto scatter_bits(Bits& bits, unsigned r, std::array<uint16_t, N> const& indices)
        -> void
    {
        for (auto i = N; i-- > 0;) {
            LATTICE_SYMMETRIES_ASSERT(indices[i] < 64, "index out of bounds");
            set_bit_to(bits, indices[i], r & 1U);
            r >>= 1U;
        }
    }

    template <class Bits, class OffDiag> struct interaction_apply_fn_t {
        Bits const&           x;
        std::complex<double>& diagonal;
        OffDiag               off_diag;

        static constexpr bool is_noexcept = noexcept(std::declval<OffDiag const&>()(
            std::declval<Bits const&>(), std::declval<std::complex<double> const&>()));

        template <unsigned N>
        auto operator()(interaction_t<N> const& self) const noexcept(is_noexcept) -> void
        {
            for (auto const edge : self.sites) {
                auto const  k    = gather_bits(x, edge);
                auto const& data = self.matrix->payload[k];
                for (auto n = 0U; n < std::size(data); ++n) {
                    if (data[n] != 0.0) {
                        if (n == k) { diagonal += data[n]; }
                        else {
                            auto y = x;
                            scatter_bits(y, n, edge);
                            off_diag(y, data[n]);
                        }
                    }
                }
            }
        }
    };

} // namespace

} // namespace lattice_symmetries

using namespace lattice_symmetries;

struct ls_interaction {
    std::variant<interaction_t<1>, interaction_t<2>, interaction_t<3>, interaction_t<4>> payload;

    template <class T, class Arg, class... Args>
    ls_interaction(std::in_place_type_t<T> tag, Arg&& arg, Args&&... args)
        : payload{tag, std::forward<Arg>(arg), std::forward<Args>(args)...}
    {}
};

namespace lattice_symmetries {
namespace {
    template <class Bits, class OffDiag>
    auto apply(ls_interaction const& interaction, Bits const& spin, std::complex<double>& diagonal,
               OffDiag off_diag) noexcept(interaction_apply_fn_t<Bits, OffDiag>::is_noexcept)
        -> void
    {
        std::visit(interaction_apply_fn_t<Bits, OffDiag>{spin, diagonal, std::move(off_diag)},
                   interaction.payload);
    }
} // namespace
} // namespace lattice_symmetries

struct ls_operator {
    struct basis_deleter_fn_t {
        auto operator()(ls_spin_basis* p) const noexcept -> void { ls_destroy_spin_basis(p); }
    };
    using basis_ptr_t = std::unique_ptr<ls_spin_basis, basis_deleter_fn_t>;

    basis_ptr_t                 basis;
    std::vector<ls_interaction> terms;
    bool                        is_real;
};

extern "C" ls_error_code ls_create_interaction1(ls_interaction** ptr, void const* matrix_2x2,
                                                unsigned const number_nodes, uint16_t const* nodes)
{
    auto p = std::make_unique<ls_interaction>(
        std::in_place_type_t<interaction_t<1>>{},
        reinterpret_cast<std::complex<double> const*>(matrix_2x2),
        tcb::span<std::array<uint16_t, 1> const>{
            reinterpret_cast<std::array<uint16_t, 1> const*>(nodes), number_nodes});
    *ptr = p.release();
    return LS_SUCCESS;
}

extern "C" ls_error_code ls_create_interaction2(ls_interaction** ptr, void const* matrix_4x4,
                                                unsigned const number_edges,
                                                uint16_t const (*edges)[2])
{
    auto p = std::make_unique<ls_interaction>(
        std::in_place_type_t<interaction_t<2>>{},
        reinterpret_cast<std::complex<double> const*>(matrix_4x4),
        tcb::span<std::array<uint16_t, 2> const>{
            reinterpret_cast<std::array<uint16_t, 2> const*>(edges), number_edges});
    *ptr = p.release();
    return LS_SUCCESS;
}

extern "C" ls_error_code ls_create_interaction3(ls_interaction** ptr, void const* matrix_8x8,
                                                unsigned const number_triangles,
                                                uint16_t const (*triangles)[3])
{
    auto p = std::make_unique<ls_interaction>(
        std::in_place_type_t<interaction_t<3>>{},
        reinterpret_cast<std::complex<double> const*>(matrix_8x8),
        tcb::span<std::array<uint16_t, 3> const>{
            reinterpret_cast<std::array<uint16_t, 3> const*>(triangles), number_triangles});
    *ptr = p.release();
    return LS_SUCCESS;
}

extern "C" ls_error_code ls_create_interaction4(ls_interaction** ptr, void const* matrix_16x16,
                                                unsigned const number_plaquettes,
                                                uint16_t const (*plaquettes)[4])
{
    auto p = std::make_unique<ls_interaction>(
        std::in_place_type_t<interaction_t<4>>{},
        reinterpret_cast<std::complex<double> const*>(matrix_16x16),
        tcb::span<std::array<uint16_t, 4> const>{
            reinterpret_cast<std::array<uint16_t, 4> const*>(plaquettes), number_plaquettes});
    *ptr = p.release();
    return LS_SUCCESS;
}

extern "C" void ls_destroy_interaction(ls_interaction* interaction)
{
    std::default_delete<ls_interaction>{}(interaction);
}

extern "C" bool ls_interaction_is_real(ls_interaction const* interaction)
{
    return std::visit(
        [](auto const& x) noexcept {
            for (auto const& column : x.matrix->payload) {
                for (auto const element : column) {
                    if (element.imag() != 0.0) { return false; }
                }
            }
            return true;
        },
        interaction->payload);
}

#if 0
extern "C" bool ls_operator_is_real(ls_operator const* op)
{
    return is_real(*op->basis)
           && std::all_of(std::begin(op->terms), std::end(op->terms),
                          [](auto const& x) { return ls_interaction_is_real(&x); });
}
#else
extern "C" bool ls_operator_is_real(ls_operator const* op) { return op->is_real; }
#endif

namespace lattice_symmetries {

template <class T, class = void> struct is_complex : std::false_type {};
template <class T>
struct is_complex<std::complex<T>, std::enable_if_t<std::is_floating_point<T>::value>>
    : std::true_type {};
template <class T> inline constexpr bool is_complex_v = is_complex<T>::value;

constexpr auto to_bits(bits512 const& x) noexcept -> uint64_t const* { return x.words; }
constexpr auto to_bits(bits512& x) noexcept -> uint64_t* { return x.words; }
constexpr auto to_bits(bits64 const& x) noexcept -> uint64_t const* { return &x; }
constexpr auto to_bits(bits64& x) noexcept -> uint64_t* { return &x; }

template <class Bits, class Callback>
auto apply_helper(ls_operator const& op, Bits const& spin, Callback callback) noexcept(
    noexcept(std::declval<Callback&>()(std::declval<Bits const&>(),
                                       std::declval<std::complex<double> const&>())))
    -> outcome::result<void>
{
    auto                 repr = spin;
    std::complex<double> eigenvalue;
    double               norm;

    ls_get_state_info(op.basis.get(), to_bits(spin), to_bits(repr), &eigenvalue, &norm);
    if (norm == 0.0) { return LS_INVALID_STATE; }
    auto const old_norm = norm;
    auto       diagonal = std::complex<double>{0.0, 0.0};
    auto const off_diag = [&](Bits const& x, std::complex<double> const& c) {
        ls_get_state_info(op.basis.get(), to_bits(x), to_bits(repr), &eigenvalue, &norm);
        if (norm > 0.0) { callback(repr, c * norm / old_norm * eigenvalue); }
    };
    for (auto const& term : op.terms) {
        apply(term, spin, diagonal, std::cref(off_diag));
    }
    callback(spin, diagonal);
    return LS_SUCCESS;
}

template <class T>
auto apply_helper(ls_operator const& op, T const* x, T* y) noexcept -> outcome::result<void>
{
    OUTCOME_TRY(states, [&op]() -> outcome::result<tcb::span<uint64_t const>> {
        ls_states* states = nullptr;
        auto const status = ls_get_states(&states, op.basis.get());
        if (status != LS_SUCCESS) { return status; }
        return tcb::span<uint64_t const>{ls_states_get_data(states), ls_states_get_size(states)};
    }());
    if constexpr (!is_complex_v<T>) {
        if (!op.is_real) { return LS_OPERATOR_IS_COMPLEX; }
    }

    auto const chunk_size =
        std::max(500UL, states.size() / (100UL * static_cast<unsigned>(omp_get_max_threads())));
    using acc_t = std::conditional_t<is_complex_v<T>, std::complex<double>, double>;
#pragma omp parallel for schedule(dynamic, chunk_size) default(none)                               \
    firstprivate(x, y, chunk_size, states) shared(op)
    for (auto i = uint64_t{0}; i < states.size(); ++i) {
        auto acc    = acc_t{0.0};
        auto status = apply_helper(
            op, states[i], [&acc, &op, x](auto const spin, auto const& coeff) noexcept {
                uint64_t   index;
                auto const _status = ls_get_index(op.basis.get(), to_bits(spin), &index);
                LATTICE_SYMMETRIES_ASSERT(_status == LS_SUCCESS, "");
                if constexpr (is_complex_v<T>) {
                    acc += std::conj(coeff) * static_cast<acc_t>(x[index]);
                }
                else {
                    acc += coeff.real() * x[index];
                }
            });
        LATTICE_SYMMETRIES_ASSERT(status, "");
        y[i] = static_cast<T>(acc);
    }
    return LS_SUCCESS;
}

template <class T> auto apply(ls_operator const* op, T const* x, T* y) noexcept -> ls_error_code
{
    auto r = apply_helper(*op, x, y);
    if (!r) {
        if (r.error().category() == get_error_category()) {
            return static_cast<ls_error_code>(r.error().value());
        }
        return LS_SYSTEM_ERROR;
    }
    return LS_SUCCESS;
}

} // namespace lattice_symmetries

extern "C" ls_error_code ls_operator_matvec_f32(ls_operator const* op, float const* x, float* y)
{
    return apply(op, x, y);
}

extern "C" ls_error_code ls_operator_matvec_f64(ls_operator const* op, double const* x, double* y)
{
    return apply(op, x, y);
}

extern "C" ls_error_code ls_operator_matvec_c64(ls_operator const* op, void const* x, void* y)
{
    using C = std::complex<float>;
    return apply(op, static_cast<C const*>(x), static_cast<C*>(y));
}

extern "C" ls_error_code ls_operator_matvec_c128(ls_operator const* op, void const* x, void* y)
{
    using C = std::complex<double>;
    return apply(op, static_cast<C const*>(x), static_cast<C*>(y));
}
