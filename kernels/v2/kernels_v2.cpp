#include "lattice_symmetries_types.h"
#include <Halide.h>
#include <callback.h>
#include <cstdint>

using namespace Halide;

namespace lattice_symmetries {

auto get_inversion_mask(unsigned number_bits) -> uint64_t {
    if (number_bits == 0) return 0;
    if (number_bits == 64) return ~uint64_t{0};
    return (uint64_t{1} << number_bits) - 1;
}

template <class T> auto from_env(std::string const &name, T const default_value) {
    auto *const value = getenv(name.c_str());
    return (value != nullptr) ? static_cast<T>(std::stoi(value)) : default_value;
}

template <class T>
auto mk_axpby_kernel()
    -> std::function<int(T, T, halide_buffer_t *, halide_buffer_t *, halide_buffer_t *)> {
    Param<T> _alpha_re{"alpha_re"};
    Param<T> _alpha_im{"alpha_im"};
    ImageParam _x{halide_type_of<T>(), 2, "x"};
    ImageParam _y{halide_type_of<T>(), 2, "y"};

    _x.dim(0).set_min(0).set_extent(2).set_stride(1);
    _x.dim(1).set_min(0).set_stride(2);
    _y.dim(0).set_min(0).set_extent(2).set_stride(1);
    _y.dim(1).set_min(0).set_stride(2);

    Func out{"out"};
    Var c{"c"};
    Var i{"i"};

    Func re{"re"};
    Func im{"im"};
    re(i) = _alpha_re * _x(0, i) - _alpha_im * _x(1, i) + _y(0, i);
    im(i) = _alpha_re * _x(1, i) + _alpha_im * _x(0, i) + _y(1, i);

    out(c, i) = select(c == 0, re(i), im(i));
    out.bound(c, 0, 2);

    Var inner{"inner"};
    Var outer{"outer"};

    auto const block_size = static_cast<int>(from_env("LS_HS_AXPBY_BLOCK_SIZE", 8));
    constexpr auto strategy = TailStrategy::GuardWithIf;
    out.split(i, outer, inner, block_size, strategy)
        .reorder(inner, c, outer)
        .vectorize(inner)
        .unroll(c);
    if (static_cast<bool>(from_env("LS_HS_AXPBY_UNROLL", 0))) {
        out.unroll(outer, 2, strategy);
    }

    out.specialize(_alpha_re == 0);
    out.specialize(_alpha_im == 0);
    // re.split(i, outer, inner, block_size, strategy)
    //     .vectorize(inner)
    //     .compute_at(out, outer)
    //     .store_in(MemoryType::Register);
    // im.split(i, outer, inner, block_size, strategy)
    //     .vectorize(inner)
    //     .compute_at(out, outer)
    //     .store_in(MemoryType::Register);

    out.output_buffer().dim(0).set_min(0).set_extent(2).set_stride(1);
    out.output_buffer().dim(1).set_min(0).set_extent(_x.dim(1).extent()).set_stride(2);
    // out.print_loop_nest();

    auto target = get_jit_target_from_environment();
    target.set_features({Target::Feature::NoAsserts, Target::Feature::NoBoundsQuery});
    if (static_cast<bool>(from_env("LS_HS_DEBUG_KERNELS", false))) {
        out.compile_to_lowered_stmt("axpy.html", {_alpha_re, _alpha_im, _x, _y}, Halide::HTML,
                                    target);
    }
    auto callable = out.compile_to_callable({_alpha_re, _alpha_im, _x, _y}, target);
    return callable.template make_std_function<T, T, halide_buffer_t *, halide_buffer_t *,
                                               halide_buffer_t *>();
}

std::once_flag get_axpby_kernel_flag;

auto const &get_axpby_kernel() {
    static std::function<int(double, double, halide_buffer_t *, halide_buffer_t *,
                             halide_buffer_t *)>
        kernel;
    std::call_once(get_axpby_kernel_flag, [&]() { kernel = mk_axpby_kernel<double>(); });
    return kernel;
}

auto invoke_axpby_kernel(void *data, va_alist alist) noexcept -> void {
    auto *fn_ptr = static_cast<std::function<int(double, double, halide_buffer_t *,
                                                 halide_buffer_t *, halide_buffer_t *)> *>(data);

    va_start_void(alist);
    auto const alpha_re = va_arg_double(alist);
    auto const alpha_im = va_arg_double(alist);
    auto *const x = va_arg_ptr(alist, halide_buffer_t *);
    auto *const y = va_arg_ptr(alist, halide_buffer_t *);
    auto *const out = va_arg_ptr(alist, halide_buffer_t *);

    try {
        (*fn_ptr)(alpha_re, alpha_im, x, y, out);
    } catch (RuntimeError &e) {
        fprintf(stderr, "Halide::RuntimeError: %s\n", e.what());
        abort();
    } catch (InternalError &e) {
        fprintf(stderr, "Halide::InternalError: %s\n", e.what());
        abort();
    }
    va_return_void(alist);
}

using FixedHammingStateToIndex = std::function<int(halide_buffer_t *, halide_buffer_t *)>;

auto mk_fixed_hamming_state_to_index_kernel(int const number_sites, int const hamming_weight,
                                            halide_buffer_t const *_binomials)
    -> FixedHammingStateToIndex {
    ImageParam _x{halide_type_of<uint64_t>(), 1, "x"};
    Func out{"out"};

    if (hamming_weight == 0 || hamming_weight == number_sites) {
        Var batch_idx{"batch_idx"};
        out(batch_idx) = cast<int64_t>(0);
    } else {
        auto binomials = Buffer{Runtime::Buffer<uint64_t, 2>{*_binomials}.copy(), "binomials"};

        Var batch_idx{"batch_idx"};
        RDom k{0, min(Expr{hamming_weight}, binomials.dim(0).extent() - 1), "k"};

        Func temp{"temp"};

        temp(batch_idx) = Tuple{cast<int64_t>(0), _x(batch_idx)};
        auto index = temp(batch_idx)[0];
        auto state = temp(batch_idx)[1];
        auto n = cast<int32_t>(count_trailing_zeros(state));
        // min(cast<int32_t>(count_trailing_zeros(state)), binomials.dim(1).extent() - 1);
        temp(batch_idx) = {index + binomials(k + 1, n), state & (state - 1)};

        out(batch_idx) = temp(batch_idx)[0];

        Var i_inner{"i_inner"};
        Var i_outer{"i_outer"};

        auto *const vector_size_str = getenv("LS_S2I_BLOCK_SIZE");
        auto vector_size = 8;
        if (vector_size_str != nullptr) {
            vector_size = std::stoi(vector_size_str);
        }

        auto *const unroll_str = getenv("LS_S2I_UNROLL");
        auto unroll = true;
        if (unroll_str != nullptr) {
            unroll = std::stoi(unroll_str) != 0;
        }

        out.split(batch_idx, i_outer, i_inner, vector_size, TailStrategy::GuardWithIf)
            .vectorize(i_inner);
        temp.compute_at(out, i_inner).store_in(MemoryType::Register);
        if (unroll) {
            temp.update(0).unroll(k);
        }
    }

    _x.dim(0).set_min(0).set_stride(1);
    out.output_buffer().dim(0).set_min(0).set_stride(1).set_extent(_x.dim(0).extent());

    // out.print_loop_nest();
    auto target = get_jit_target_from_environment();
    target.set_features({Target::Feature::NoAsserts, Target::Feature::NoBoundsQuery});
    // out.compile_to_lowered_stmt("state_to_index.html", {_x}, Halide::HTML, target);
    // out.compile_to_assembly("state_to_index.S", {_x}, "fixed_hamming_state_to_index", target);
    auto callable = out.compile_to_callable({_x}, target);
    return callable.make_std_function<halide_buffer_t *, halide_buffer_t *>();
}

auto invoke_fixed_hamming_state_to_index_kernel(void *data, va_alist alist) noexcept -> void {
    auto *fn_ptr = static_cast<FixedHammingStateToIndex *>(data);
    va_start_void(alist);
    auto *basis_states = va_arg_ptr(alist, halide_buffer_t *);
    auto *indices = va_arg_ptr(alist, halide_buffer_t *);
    try {
        (*fn_ptr)(basis_states, indices);
    } catch (RuntimeError &e) {
        fprintf(stderr, "Halide::RuntimeError: %s\n", e.what());
        abort();
    } catch (InternalError &e) {
        fprintf(stderr, "Halide::InternalError: %s\n", e.what());
        abort();
    }
    va_return_void(alist);
}

auto bit_permute_step_64(Expr x, Expr m, Expr d) -> Expr {
    auto y = ((x >> d) ^ x) & m;
    return (x ^ y) ^ (y << d);
}

template <class X, class Masks>
auto mk_permuted(X const& x, Masks const& masks, halide_buffer_t const * _shifts) {
    Var batch_idx{"batch_idx"};
    Var group_idx{"group_idx"};
    Func y{"y"};

    auto const depth = masks.dim(1).extent();
    auto const *shifts_raw = reinterpret_cast<uint64_t const *>(_shifts->host);
    Expr p = x(batch_idx);
    for (auto i = 0; i < depth; ++i) {
        p = bit_permute_step_64(p, masks(group_idx, i), cast<uint64_t>(Expr{shifts_raw[i]}));
    }
    y(batch_idx, group_idx) = p;

    return y;
}

auto mk_is_representative_kernel(halide_buffer_t const *_masks, halide_buffer_t const *_eigvals_re,
                                 halide_buffer_t const *_shifts, unsigned number_bits, int spin_inversion)
    -> Halide::Callable {
    auto masks = Halide::Buffer{Halide::Runtime::Buffer<uint64_t, 2>{*_masks}.copy(), "masks"};
    auto eigvals_re =
        Halide::Buffer{Halide::Runtime::Buffer<double, 1>{*_eigvals_re}.copy(), "eigvals_re"};
    auto shifts = Halide::Buffer{Halide::Runtime::Buffer<uint64_t, 1>{*_shifts}.copy(), "shifts"};
    ImageParam _x{halide_type_of<uint64_t>(), 1, "x"};

    auto const depth = masks.dim(1).extent();
    auto const number_masks = masks.dim(0).extent();
    auto const count = _x.dim(0).extent();
    auto const inversion_mask = Expr{get_inversion_mask(number_bits)};

    Var batch_idx{"batch_idx"};
    Var group_idx{"group_idx"};
    auto y = mk_permuted(_x, masks, _shifts);

    Func temp{"temp"};
    RDom rgroup_idx{1, number_masks - 1, "rgroup_idx"};

    auto init_n = cast<uint16_t>(1);
    if (spin_inversion != 0) {
        auto const inverted = _x(batch_idx) ^ inversion_mask;
        auto const is_greater = inverted > _x(batch_idx);
        init_n = cast<uint16_t>(is_greater);
    }
    temp(batch_idx) = init_n;

    auto const current_n = temp(batch_idx);
    auto const is_greater = y(batch_idx, rgroup_idx) > _x(batch_idx);
    auto const is_equal = y(batch_idx, rgroup_idx) == _x(batch_idx);
    auto const is_trivial = eigvals_re(rgroup_idx) == cast<double>(1);
    auto next_n = 
        cast<uint16_t>(is_greater || (is_equal && is_trivial)) * (current_n + cast<uint16_t>(is_equal));
    if (spin_inversion != 0) {
        auto inverted = y(batch_idx, rgroup_idx) ^ inversion_mask;
        auto const is_greater = inverted > _x(batch_idx);
        auto const is_equal = inverted == _x(batch_idx);
        auto const is_trivial = eigvals_re(rgroup_idx) == cast<double>(spin_inversion);
        next_n = cast<uint16_t>(is_greater || (is_equal && is_trivial)) * (next_n + cast<uint16_t>(is_equal));
    }

    rgroup_idx.where(current_n > 0);
    temp(batch_idx) = next_n;

    Func norm{"norm"};
    norm(batch_idx) = temp(batch_idx);

    auto const block_size = from_env("LS_HS_IS_REPRESENTATIVE_BLOCK_SIZE", 1);
    LS_CHECK(block_size <= LS_HS_MAX_BLOCK_SIZE, "block_size too big");
    // Shapes & strides
    _x.dim(0).set_min(0).set_stride(1).set_extent(block_size * (_x.dim(0).extent() / block_size));
    norm.output_buffer().dim(0).set_min(0).set_stride(1).set_extent(_x.dim(0).extent());

    if (block_size > 1) {
        Var outer{"outer"};
        Var inner{"inner"};
        norm.split(batch_idx, outer, inner, block_size).vectorize(inner);

        temp.split(batch_idx, outer, inner, block_size).vectorize(inner);
        temp.update(0).split(batch_idx, outer, inner, block_size).reorder(inner, rgroup_idx, outer).vectorize(inner);
        temp.compute_at(norm, outer).store_in(MemoryType::Register);
    }
    else {
        temp.compute_at(norm, batch_idx).store_in(MemoryType::Register);
    }

    auto target = get_jit_target_from_environment();
    // target.set_features({Target::Feature::NoAsserts, Target::Feature::NoBoundsQuery});
    if (static_cast<bool>(from_env("LS_HS_DEBUG_KERNELS", false))) {
        norm.print_loop_nest();
        norm.compile_to_lowered_stmt("is_representative.html", {_x}, Halide::HTML, target);
    }
    auto callable = norm.compile_to_callable({_x}, target);
    return callable;
}

auto invoke_is_representative_kernel(void *data, va_alist alist) noexcept -> void {
    auto *ctx = static_cast<std::tuple<Halide::Callable, void *> *>(data);
    auto const &fn = std::get<0>(*ctx);
    va_start_void(alist);
    auto *basis_states = va_arg_ptr(alist, halide_buffer_t *);
    auto *norms = va_arg_ptr(alist, halide_buffer_t *);
    try {
        fn(basis_states, norms);
    } catch (RuntimeError &e) {
        fprintf(stderr, "Halide::RuntimeError: %s\n", e.what());
        abort();
    } catch (InternalError &e) {
        fprintf(stderr, "Halide::InternalError: %s\n", e.what());
        abort();
    }
    va_return_void(alist);
}

auto mk_state_info_kernel_v3(halide_buffer_t const *_masks, halide_buffer_t const *_shifts,
                             unsigned number_bits, bool spin_inversion) -> Halide::Callable {
    auto masks = Halide::Buffer{Halide::Runtime::Buffer<uint64_t, 2>{*_masks}.copy(), "masks"};
    auto shifts = Halide::Buffer{Halide::Runtime::Buffer<uint64_t, 1>{*_shifts}.copy(), "shifts"};
    ImageParam _x{halide_type_of<uint64_t>(), 1, "basis_states"};

    auto depth = masks.dim(1).extent();
    auto number_masks = masks.dim(0).extent();

    auto const inversion_mask = Expr{get_inversion_mask(number_bits)};

    Var batch_idx{"batch_idx"};
    Var group_idx{"group_idx"};
    RDom depth_idx{0, depth, "depth_idx"};
    Func y{"y"};

    auto const unroll = static_cast<bool>(from_env("LS_HS_STATE_INFO_UNROLL", 1));
    if (unroll) {
        auto const *shifts_raw = reinterpret_cast<uint64_t const *>(_shifts->host);
        Expr p = _x(batch_idx);
        for (auto i = 0; i < depth; ++i) {
            p = bit_permute_step_64(p, masks(group_idx, i), cast<uint64_t>(Expr{shifts_raw[i]}));
        }
        y(batch_idx, group_idx) = p;
    } else {
        y(batch_idx, group_idx) = _x(batch_idx);
        y(batch_idx, group_idx) = bit_permute_step_64(
            y(batch_idx, group_idx), masks(group_idx, depth_idx), shifts(depth_idx));
    }

    Func temp{"temp"};
    RDom rgroup_idx{1, number_masks - 1, "rgroup_idx"};
    auto init_r = _x(batch_idx);
    auto init_j = cast<int32_t>(0);
    if (spin_inversion) {
        auto const inverted = _x(batch_idx) ^ inversion_mask;
        auto const is_less = inverted < init_r;
        init_r = select(is_less, inverted, init_r);
        init_j = select(is_less, number_masks, init_j);
    }
    temp(batch_idx) = Tuple{init_r, init_j};
    auto const current_r = temp(batch_idx)[0];
    auto const current_j = temp(batch_idx)[1];
    auto const is_less = y(batch_idx, rgroup_idx) < current_r;
    auto next_r = select(is_less, y(batch_idx, rgroup_idx), current_r);
    auto next_j = select(is_less, rgroup_idx, current_j);
    if (spin_inversion) {
        auto const inverted = y(batch_idx, rgroup_idx) ^ inversion_mask;
        auto const is_less2 = inverted < next_r;
        next_r = select(is_less2, inverted, next_r);
        next_j = select(is_less2, number_masks + rgroup_idx, next_j);
    }
    temp(batch_idx) = Tuple{next_r, next_j};

    auto const block_size = from_env("LS_HS_STATE_INFO_BLOCK_SIZE", 16);
    LS_CHECK(block_size <= LS_HS_MAX_BLOCK_SIZE, "block_size too big");
    Func out{"out"};
    out(batch_idx) = temp(batch_idx);

    _x.dim(0).set_min(0).set_stride(1).set_extent(block_size * (_x.dim(0).extent() / block_size));
    out.output_buffers()[0].dim(0).set_min(0).set_stride(1).set_extent(_x.dim(0).extent());

    Var outer{"outer"};
    Var inner{"inner"};

    out.split(batch_idx, outer, inner, block_size).vectorize(inner);

    temp.split(batch_idx, outer, inner, block_size).vectorize(inner);
    temp.update(0)
        .split(batch_idx, outer, inner, block_size)
        .reorder(inner, rgroup_idx, outer)
        .vectorize(inner);
    temp.compute_at(out, outer).store_in(MemoryType::Register);

    if (unroll) {
        y.split(batch_idx, outer, inner, block_size)
            .reorder(inner, group_idx, outer)
            .vectorize(inner);
    } else {
        y.split(batch_idx, outer, inner, block_size).vectorize(inner);
        y.update(0)
            .split(batch_idx, outer, inner, block_size)
            .reorder(inner, depth_idx, group_idx, outer)
            .vectorize(inner);
    }
    y.compute_at(temp, rgroup_idx).store_in(MemoryType::Register);

    // out.print_loop_nest();
    auto target = get_jit_target_from_environment();
    target.set_features({Target::Feature::NoAsserts, Target::Feature::NoBoundsQuery});
    if (static_cast<bool>(from_env("LS_HS_DEBUG_KERNELS", false))) {
        out.print_loop_nest();
        out.compile_to_lowered_stmt("state_info.html", {_x}, Halide::HTML, target);
    }
    auto callable = out.compile_to_callable({_x}, target);
    return callable;
}

auto mk_state_info_kernel(halide_buffer_t const *_masks, halide_buffer_t const *_shifts,
                          int spin_inversion) -> Halide::Callable {

    auto masks = Halide::Buffer{Halide::Runtime::Buffer<uint64_t, 2>{*_masks}.copy(), "masks"};
    auto shifts = Halide::Buffer{Halide::Runtime::Buffer<uint64_t, 1>{*_shifts}.copy(), "shifts"};
    // auto const flip_mask = 0;

    ImageParam _x{halide_type_of<uint64_t>(), 1, "basis_states"};
    // ImageParam _repr{halide_type_of<uint64_t>(), 1, "representatives"};
    // ImageParam _indices{halide_type_of<int32_t>(), 1, "indices"};

    auto depth = masks.dim(0).extent();
    auto number_masks = masks.dim(1).extent();
    auto const chunk_size = 2 * 4;
    auto number_chunks = number_masks / chunk_size;
    auto number_rest = number_masks - number_chunks * chunk_size;

    auto const reduction_step_impl = [](Expr /*i*/, Tuple current, Expr y, Expr j) -> Tuple {
        auto current_r = current[0];
        auto current_j = current[1];

        auto is_less = y < current_r;
        auto next_r = select(is_less, y, current_r);
        auto next_j = select(is_less, j, current_j);

        return {next_r, next_j};
    };

    auto const reduction_step = [spin_inversion, reduction_step_impl](Expr i, Tuple current, Expr y,
                                                                      Expr j) -> Tuple {
        if (spin_inversion == 0) {
            return reduction_step_impl(i, current, y, j);
        }
        // if (spin_inversion == 1) {
        //     auto temp = reduction_step_impl(i, current, y, j);
        //     return reduction_step_impl(i, temp, y ^ flip_mask, c_re, c_im);
        // }
        // if (spin_inversion == -1) {
        //     auto temp = reduction_step_impl(i, current, y, c_re, c_im);
        //     return reduction_step_impl(i, temp, y ^ flip_mask, -c_re, -c_im);
        // }
        throw std::runtime_error{"invalid spin_inversion"};
    };

    auto const combine_step = [](Tuple current, Tuple other) -> Tuple {
        auto current_r = current[0];
        auto current_j = current[1];
        auto r = other[0];
        auto j = other[1];

        auto is_less = r < current_r;
        auto next_r = select(is_less, r, current_r);
        auto next_j = select(is_less, j, current_j);

        return {next_r, next_j};
    };

    Var i{"i"};

    // Apply masks to generate transformed spin configurations
    RDom k{0, depth, "k"};
    Func y_batched{"y_batched"};
    Var j_outer{"j_outer"};
    Var j_inner{"j_inner"};
    y_batched(i, j_outer, j_inner) = _x(i);
    y_batched(i, j_outer, j_inner) = bit_permute_step_64(
        y_batched(i, j_outer, j_inner), masks(k, j_outer * chunk_size + j_inner), shifts(k));

    Func y_scalar{"y_scalar"};
    Var j_tail{"j_tail"};
    y_scalar(i, j_tail) = _x(i);
    y_scalar(i, j_tail) = bit_permute_step_64(y_scalar(i, j_tail), masks(k, j_tail), shifts(k));

    // Compute vectorized reduction
    RDom m_main{0, number_chunks, "m_main"};
    Func custom_reduction{"custom_reduction"};
    custom_reduction(i, j_inner) = Tuple{_x(i), cast<int32_t>(-1)};
    custom_reduction(i, j_inner) =
        reduction_step(i, custom_reduction(i, j_inner), y_batched(i, m_main, j_inner),
                       m_main * chunk_size + j_inner);

    // Reduce across one chunk
    RDom m_lane{1, chunk_size - 1, "m_lane"};
    Func reduction_lane{"reduction_lane"};
    reduction_lane(i) = Tuple{custom_reduction(i, 0)[0], custom_reduction(i, 0)[1]};
    reduction_lane(i) = combine_step(reduction_lane(i), custom_reduction(i, m_lane));

    // Compute reduction over remaining elements
    Func reduction_scalar{"reduction_scalar"};
    RDom m_tail{number_chunks * chunk_size, number_rest, "m_tail"};
    reduction_scalar(i) = Tuple{reduction_lane(i)[0], reduction_lane(i)[1]};
    reduction_scalar(i) = reduction_step(i, reduction_scalar(i), y_scalar(i, m_tail), m_tail);

    // Schedule
    y_batched.compute_at(custom_reduction, m_main);

    custom_reduction.compute_at(reduction_scalar, i);
    custom_reduction.vectorize(j_inner);
    custom_reduction.update(0).vectorize(j_inner);

    return reduction_scalar.compile_to_callable({_x});
}

auto invoke_state_info_kernel(void *data, va_alist alist) noexcept -> void {
    auto *ctx = static_cast<std::tuple<Halide::Callable, void *> *>(data);
    auto const &fn = std::get<0>(*ctx);
    va_start_void(alist);
    auto *basis_states = va_arg_ptr(alist, halide_buffer_t *);
    auto *representatives = va_arg_ptr(alist, halide_buffer_t *);
    auto *indices = va_arg_ptr(alist, halide_buffer_t *);
    try {
        fn(basis_states, representatives, indices);
    } catch (RuntimeError &e) {
        fprintf(stderr, "Halide::RuntimeError: %s\n", e.what());
        abort();
    } catch (InternalError &e) {
        fprintf(stderr, "Halide::InternalError: %s\n", e.what());
        abort();
    }
    va_return_void(alist);
}

} // namespace lattice_symmetries

extern "C" ls_hs_state_to_index_kernel_type ls_hs_internal_mk_fixed_hamming_state_to_index_kernel(
    int const number_sites, int const hamming_weight, halide_buffer_t const *binomials) {
    using namespace lattice_symmetries;
    try {
        // fprintf(stderr, "mk_fixed_hamming_state_to_index_kernel\n");
        auto fn = mk_fixed_hamming_state_to_index_kernel(number_sites, hamming_weight, binomials);
        auto closure = alloc_callback(&invoke_fixed_hamming_state_to_index_kernel,
                                      reinterpret_cast<void *>(new FixedHammingStateToIndex{fn}));
        return reinterpret_cast<ls_hs_state_to_index_kernel_type>(closure);
    } catch (CompileError &e) {
        fprintf(stderr, "Halide::CompileError: %s\n", e.what());
        abort();
    } catch (InternalError &e) {
        fprintf(stderr, "Halide::InternalError: %s\n", e.what());
        abort();
    } catch (RuntimeError &e) {
        fprintf(stderr, "Halide::RuntimeError: %s\n", e.what());
        abort();
    }
}
extern "C" void ls_hs_internal_destroy_fixed_hamming_state_to_index_kernel(
    ls_hs_state_to_index_kernel_type closure) {
    using namespace lattice_symmetries;

    if (closure == nullptr) {
        return;
    }
    LS_CHECK(is_callback(reinterpret_cast<void *>(closure)),
             "trying to destroy a normal function pointer. This "
             "should never happen. Please, submit a bug report.");

    FixedHammingStateToIndex *fn_ptr = reinterpret_cast<FixedHammingStateToIndex *>(
        callback_data(reinterpret_cast<callback_t>(closure)));
    LS_CHECK(fn_ptr != nullptr,
             "callback_data is NULL. This should never happen. Please, submit a bug report.");

    delete fn_ptr;
    free_callback(reinterpret_cast<callback_t>(closure));
}

extern "C" ls_hs_is_representative_kernel_type_v2
ls_hs_internal_mk_is_representative_kernel(halide_buffer_t const *masks,
                                           halide_buffer_t const *eigvals_re,
                                           halide_buffer_t const *shifts,
                                           unsigned number_bits,
                                           int spin_inversion) {
    try {
        auto fn = lattice_symmetries::mk_is_representative_kernel(masks, eigvals_re, shifts,
                                                                  number_bits, spin_inversion);
        auto closure = alloc_callback(&lattice_symmetries::invoke_is_representative_kernel,
                                      new std::tuple<Halide::Callable, void *>{fn, nullptr});
        return reinterpret_cast<ls_hs_is_representative_kernel_type_v2>(closure);
    } catch (CompileError &e) {
        fprintf(stderr, "Halide::CompileError: %s\n", e.what());
        abort();
    } catch (InternalError &e) {
        fprintf(stderr, "Halide::InternalError: %s\n", e.what());
        abort();
    } catch (RuntimeError &e) {
        fprintf(stderr, "Halide::RuntimeError: %s\n", e.what());
        abort();
    }
}
extern "C" void
ls_hs_internal_destroy_is_representative_kernel(ls_hs_is_representative_kernel_type_v2 closure) {
    using namespace lattice_symmetries;

    if (closure == nullptr) {
        return;
    }
    LS_CHECK(is_callback(reinterpret_cast<void *>(closure)),
             "trying to destroy a normal function pointer. This "
             "should never happen. Please, submit a bug report.");

    auto *data = reinterpret_cast<std::tuple<Halide::Callable, void *> *>(
        callback_data(reinterpret_cast<callback_t>(closure)));
    LS_CHECK(data != nullptr,
             "callback_data is NULL. This should never happen. Please, submit a bug report.");
    delete data;
    free_callback(reinterpret_cast<callback_t>(closure));
}

extern "C" ls_hs_state_info_kernel_type_v2
ls_hs_internal_mk_state_info_kernel_v3(halide_buffer_t const *masks, halide_buffer_t const *shifts,
                                       unsigned number_bits, bool spin_inversion) {
    try {
        auto fn =
            lattice_symmetries::mk_state_info_kernel_v3(masks, shifts, number_bits, spin_inversion);
        auto closure = alloc_callback(&lattice_symmetries::invoke_state_info_kernel,
                                      new std::tuple<Halide::Callable, void *>{fn, nullptr});
        return reinterpret_cast<ls_hs_state_info_kernel_type_v2>(closure);
    } catch (CompileError &e) {
        fprintf(stderr, "Halide::CompileError: %s\n", e.what());
        abort();
    } catch (InternalError &e) {
        fprintf(stderr, "Halide::InternalError: %s\n", e.what());
        abort();
    } catch (RuntimeError &e) {
        fprintf(stderr, "Halide::RuntimeError: %s\n", e.what());
        abort();
    }
}
extern "C" void
ls_hs_internal_destroy_state_info_kernel_v3(ls_hs_state_info_kernel_type_v2 closure) {
    using namespace lattice_symmetries;

    if (closure == nullptr) {
        return;
    }
    LS_CHECK(is_callback(reinterpret_cast<void *>(closure)),
             "trying to destroy a normal function pointer. This "
             "should never happen. Please, submit a bug report.");

    auto *data = reinterpret_cast<std::tuple<Halide::Callable, void *> *>(
        callback_data(reinterpret_cast<callback_t>(closure)));
    LS_CHECK(data != nullptr,
             "callback_data is NULL. This should never happen. Please, submit a bug report.");
    delete data;
    free_callback(reinterpret_cast<callback_t>(closure));
}

extern "C" ls_hs_state_info_kernel_type_v2
ls_hs_internal_mk_state_info_kernel(halide_buffer_t const *masks, halide_buffer_t const *shifts,
                                    int spin_inversion) {
    try {
        auto fn = lattice_symmetries::mk_state_info_kernel(masks, shifts, spin_inversion);
        auto closure = alloc_callback(&lattice_symmetries::invoke_state_info_kernel,
                                      new std::tuple<Halide::Callable, void *>{fn, nullptr});
        return reinterpret_cast<ls_hs_state_info_kernel_type_v2>(closure);
    } catch (CompileError &e) {
        fprintf(stderr, "Halide::CompileError: %s\n", e.what());
        abort();
    } catch (InternalError &e) {
        fprintf(stderr, "Halide::InternalError: %s\n", e.what());
        abort();
    } catch (RuntimeError &e) {
        fprintf(stderr, "Halide::RuntimeError: %s\n", e.what());
        abort();
    }
}
extern "C" void ls_hs_internal_destroy_state_info_kernel(ls_hs_state_info_kernel_type_v2 closure) {
    using namespace lattice_symmetries;

    if (closure == nullptr) {
        return;
    }
    LS_CHECK(is_callback(reinterpret_cast<void *>(closure)),
             "trying to destroy a normal function pointer. This "
             "should never happen. Please, submit a bug report.");

    auto *data = reinterpret_cast<std::tuple<Halide::Callable, void *> *>(
        callback_data(reinterpret_cast<callback_t>(closure)));
    LS_CHECK(data != nullptr,
             "callback_data is NULL. This should never happen. Please, submit a bug report.");
    delete data;
    free_callback(reinterpret_cast<callback_t>(closure));
}

extern "C" void ls_hs_internal_axpy(int64_t const size, double alpha_re, double alpha_im,
                                    ls_hs_scalar const *x, ls_hs_scalar const *y,
                                    ls_hs_scalar *out) {
    try {
        LS_CHECK(size < std::numeric_limits<int32_t>::max(), "integer overflow");
        auto const &fn = lattice_symmetries::get_axpby_kernel();
        halide_dimension_t dims[2] = {halide_dimension_t{0, (int32_t)2, 1},
                                      halide_dimension_t{0, (int32_t)size, 2}};
        auto x_buf = halide_buffer_t{.device = 0,
                                     .device_interface = 0,
                                     .host = (uint8_t *)x,
                                     .flags = 0,
                                     .type = halide_type_of<double>(),
                                     .dimensions = 2,
                                     .dim = dims,
                                     .padding = 0};
        auto y_buf = halide_buffer_t{.device = 0,
                                     .device_interface = 0,
                                     .host = (uint8_t *)y,
                                     .flags = 0,
                                     .type = halide_type_of<double>(),
                                     .dimensions = 2,
                                     .dim = dims,
                                     .padding = 0};
        auto out_buf = halide_buffer_t{.device = 0,
                                       .device_interface = 0,
                                       .host = (uint8_t *)out,
                                       .flags = 0,
                                       .type = halide_type_of<double>(),
                                       .dimensions = 2,
                                       .dim = dims,
                                       .padding = 0};
        fn(alpha_re, alpha_im, &x_buf, &y_buf, &out_buf);
    } catch (CompileError &e) {
        fprintf(stderr, "Halide::CompileError: %s\n", e.what());
        abort();
    } catch (InternalError &e) {
        fprintf(stderr, "Halide::InternalError: %s\n", e.what());
        abort();
    } catch (RuntimeError &e) {
        fprintf(stderr, "Halide::RuntimeError: %s\n", e.what());
        abort();
    }
}
