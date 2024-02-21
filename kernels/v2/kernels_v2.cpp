#include "lattice_symmetries_types.h"
#include <Halide.h>
#include <callback.h>
#include <cstdint>

using namespace Halide;

namespace lattice_symmetries {

auto from_env(std::string const &name, int const default_value) {
    auto *const value = getenv(name.c_str());
    return (value != nullptr) ? std::stoi(value) : default_value;
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

auto mk_is_representative_kernel(halide_buffer_t const *_masks, halide_buffer_t const *_eigvals_re,
                                 halide_buffer_t const *_shifts, int spin_inversion)
    -> Halide::Callable {
    auto masks = Halide::Buffer{Halide::Runtime::Buffer<uint64_t, 2>{*_masks}.copy(), "masks"};
    auto eigvals_re =
        Halide::Buffer{Halide::Runtime::Buffer<double, 1>{*_eigvals_re}.copy(), "eigvals_re"};
    auto shifts = Halide::Buffer{Halide::Runtime::Buffer<uint64_t, 1>{*_shifts}.copy(), "shifts"};
    // auto const flip_mask = 0;

    ImageParam _x{halide_type_of<uint64_t>(), 1, "x"};
    // ImageParam _is_representative{"is_representative", 1};
    // ImageParam _norm{"norm", 1};

    auto count = _x.dim(0).extent();
    auto depth = masks.dim(0).extent();
    auto number_masks = masks.dim(1).extent();
    auto const chunk_size = 4; // natural_vector_size(type_of<uint64_t>());
    auto number_chunks = number_masks / chunk_size;
    auto number_rest = number_masks - number_chunks * chunk_size;

    Var i{"i"};
    Var j_outer{"j_outer"};
    Var j_inner{"j_inner"};

    // Transform input
    RDom k{0, depth, "k"};
    Func y{"y"};
    y(i, j_outer, j_inner) = _x(i);
    y(i, j_outer, j_inner) = bit_permute_step_64(
        y(i, j_outer, j_inner), masks(k, j_outer * chunk_size + j_inner), shifts(k));
    Func temp{"temp"};
    if (spin_inversion == 0) {
        temp(i, j_outer, j_inner) = Tuple{cast<double>(y(i, j_outer, j_inner) == _x(i)) *
                                              eigvals_re(j_outer * chunk_size + j_inner),
                                          y(i, j_outer, j_inner) >= _x(i)};
        // } else if (spin_inversion == 1) {
        //   auto y_flipped = y(i, j_outer, j_inner) ^ flip_mask;
        //   temp(i, j_outer, j_inner) =
        //       Tuple{(cast<double>(y(i, j_outer, j_inner) == _x(i)) +
        //              cast<double>(y_flipped == _x(i))) *
        //                 eigvals_re(j_outer * chunk_size + j_inner),
        //             y(i, j_outer, j_inner) >= _x(i) && y_flipped >= _x(i)};
        // } else if (spin_inversion == -1) {
        //   auto y_flipped = y(i, j_outer, j_inner) ^ flip_mask;
        //   temp(i, j_outer, j_inner) =
        //       Tuple{(cast<double>(y(i, j_outer, j_inner) == _x(i)) -
        //              cast<double>(y_flipped == _x(i))) *
        //                 eigvals_re(j_outer * chunk_size + j_inner),
        //             y(i, j_outer, j_inner) >= _x(i) && y_flipped >= _x(i)};
    } else {
        throw std::runtime_error{"invalid spin_inversion"};
    }

    // Horizontal reduction
    RDom m_lane{0, chunk_size, "m_lane"};
    Func lane_reduction{"lane_reduction"};
    lane_reduction(i, j_outer) = Tuple{cast<double>(0), cast<bool>(true)};
    lane_reduction(i, j_outer) =
        Tuple{lane_reduction(i, j_outer)[0] + temp(i, j_outer, m_lane)[0],
              lane_reduction(i, j_outer)[1] && temp(i, j_outer, m_lane)[1]};

    // Vertical reduction
    RDom m_batch{0, number_chunks, "m_batch"};
    Func batch_reduction{"batch_reduction"};
    batch_reduction(i) = Tuple{cast<double>(0), cast<bool>(true)};
    m_batch.where(batch_reduction(i)[1]);
    batch_reduction(i) = Tuple{batch_reduction(i)[0] + lane_reduction(i, m_batch)[0],
                               batch_reduction(i)[1] && lane_reduction(i, m_batch)[1]};

    // Compute reduction over remaining elements
    Func scalar_reduction{"scalar_reduction"};
    RDom m_tail{0, number_rest, "m_tail"};
    scalar_reduction(i) = Tuple{batch_reduction(i)[0], batch_reduction(i)[1]};
    // m_tail.where(batch_reduction(i)[1]);
    scalar_reduction(i) = Tuple{scalar_reduction(i)[0] + temp(i, number_chunks, m_tail)[0],
                                scalar_reduction(i)[1] && temp(i, number_chunks, m_tail)[1]};

    // Save results
    auto const cutoff = Expr{1024 * std::numeric_limits<double>::epsilon()};
    auto const n = scalar_reduction(i)[0];
    auto const rounded = select(n < cutoff, cast<double>(0), n);

    Func norm{"norm"};
    norm(i) = select(n > cutoff && scalar_reduction(i)[1], n, 0);

    // Shapes & strides
    _x.dim(0).set_min(0).set_stride(1);

    // Schedule
    temp.in(lane_reduction).vectorize(j_inner);
    temp.in(lane_reduction).compute_at(lane_reduction, j_outer);
    temp.in(scalar_reduction).compute_at(scalar_reduction, i);
    // lane_reduction.compute_at(batch_reduction, m_batch);
    scalar_reduction.compute_at(norm, i);

    // norm.print_loop_nest();
    return norm.compile_to_callable({_x});
}

auto invoke_is_representative_kernel(void *data, va_alist alist) noexcept -> void {
    auto *ctx = static_cast<std::tuple<Halide::Callable, void *> *>(data);
    auto const &fn = std::get<0>(*ctx);
    va_start_void(alist);
    auto *basis_states = va_arg_ptr(alist, halide_buffer_t *);
    auto *norms = va_arg_ptr(alist, halide_buffer_t *);
    try {
        fn(basis_states, norms);
    } catch (InternalError &e) {
        fprintf(stderr, "Halide::InternalError: %s\n", e.what());
        abort();
    }
    va_return_void(alist);
}

auto mk_state_info_kernel_v3(halide_buffer_t const *_masks, halide_buffer_t const *_shifts,
                             int spin_inversion) -> Halide::Callable {
    auto masks = Halide::Buffer{Halide::Runtime::Buffer<uint64_t, 2>{*_masks}.copy(), "masks"};
    auto shifts = Halide::Buffer{Halide::Runtime::Buffer<uint64_t, 1>{*_shifts}.copy(), "shifts"};
    ImageParam _x{halide_type_of<uint64_t>(), 1, "basis_states"};

    auto depth = masks.dim(1).extent();
    auto number_masks = masks.dim(0).extent();

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
    RDom rgroup_idx{0, number_masks, "rgroup_idx"};
    temp(batch_idx) = Tuple{_x(batch_idx), cast<int32_t>(-1)};
    auto const current_r = temp(batch_idx)[0];
    auto const current_j = temp(batch_idx)[1];
    auto const is_less = y(batch_idx, rgroup_idx) < current_r;
    auto const next_r = select(is_less, y(batch_idx, rgroup_idx), current_r);
    auto const next_j = select(is_less, rgroup_idx, current_j);
    temp(batch_idx) = Tuple{next_r, next_j};

    Func out{"out"};
    out(batch_idx) = temp(batch_idx);

    _x.dim(0).set_min(0).set_stride(1);
    out.output_buffers()[0].dim(0).set_min(0).set_stride(1).set_extent(_x.dim(0).extent());

    Var outer{"outer"};
    Var inner{"inner"};
    auto const block_size = from_env("LS_HS_STATE_INFO_BLOCK_SIZE", 16);
    constexpr auto strategy = TailStrategy::RoundUp;

    out.split(batch_idx, outer, inner, block_size, strategy).vectorize(inner);

    temp.split(batch_idx, outer, inner, block_size, strategy).vectorize(inner);
    temp.update(0)
        .split(batch_idx, outer, inner, block_size, strategy)
        .reorder(inner, rgroup_idx, outer)
        .vectorize(inner);
    temp.compute_at(out, outer).store_in(MemoryType::Register);

    if (unroll) {
        y.split(batch_idx, outer, inner, block_size, strategy)
            .reorder(inner, group_idx, outer)
            .vectorize(inner);
    } else {
        y.split(batch_idx, outer, inner, block_size, strategy).vectorize(inner);
        y.update(0)
            .split(batch_idx, outer, inner, block_size, strategy)
            .reorder(inner, depth_idx, group_idx, outer)
            .vectorize(inner);
    }
    y.compute_at(temp, rgroup_idx).store_in(MemoryType::Register);

    // out.print_loop_nest();
    auto target = get_jit_target_from_environment();
    target.set_features({Target::Feature::NoAsserts, Target::Feature::NoBoundsQuery});
    // out.compile_to_lowered_stmt("state_info.html", {_x}, Halide::HTML, target);
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
                                           halide_buffer_t const *shifts, int spin_inversion) {
    try {
        auto fn = lattice_symmetries::mk_is_representative_kernel(masks, eigvals_re, shifts,
                                                                  spin_inversion);
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
                                       int spin_inversion) {
    try {
        auto fn = lattice_symmetries::mk_state_info_kernel_v3(masks, shifts, spin_inversion);
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
