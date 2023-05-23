#include "Halide.h"

using namespace Halide;

auto bit_permute_step_64(Expr x, Expr m, Expr d) -> Expr {
  auto y = ((x >> d) ^ x) & m;
  return (x ^ y) ^ (y << d);
}

class state_info_generator : public Halide::Generator<state_info_generator> {
public:
  GeneratorParam<int> _spin_inversion{"spin_inversion", /*default value*/ 0};

  Input<Buffer<uint64_t>> _x{"x", 1};
  Input<uint64_t> _flip_mask{"flip_mask"};
  Input<Buffer<uint64_t>> _masks{"masks", 2};
  Input<Buffer<double>> _eigvals_re{"eigvals_re", 1};
  Input<Buffer<double>> _eigvals_im{"eigvals_im", 1};
  Input<Buffer<uint64_t>> _shifts{"shifts", 1};
  Output<Buffer<uint64_t>> _repr{"representative", 1};
  Output<Buffer<double>> _character{"character", 2};
  Output<Buffer<double>> _norm{"norm"};

  auto reduction_step_impl(Expr i, Tuple current, Expr y, Expr c_re,
                           Expr c_im) const -> Tuple {
    auto current_r = current[0];
    auto current_c_re = current[1];
    auto current_c_im = current[2];
    auto current_n = current[3];

    auto is_less = y < current_r;
    auto next_r = select(is_less, y, current_r);
    auto next_c_re = select(is_less, c_re, current_c_re);
    auto next_c_im = select(is_less, c_im, current_c_im);

    auto is_equal = y == _x(i);
    auto next_n = select(is_equal, current_n + c_re, current_n);

    return {next_r, next_c_re, next_c_im, next_n};
  }

  auto reduction_step(Expr i, Tuple current, Expr y, Expr c_re, Expr c_im) const
      -> Tuple {
    if (_spin_inversion == 0) {
      return reduction_step_impl(i, current, y, c_re, c_im);
    }
    if (_spin_inversion == 1) {
      auto temp = reduction_step_impl(i, current, y, c_re, c_im);
      return reduction_step_impl(i, temp, y ^ _flip_mask, c_re, c_im);
    }
    if (_spin_inversion == -1) {
      auto temp = reduction_step_impl(i, current, y, c_re, c_im);
      return reduction_step_impl(i, temp, y ^ _flip_mask, -c_re, -c_im);
    }
    throw std::runtime_error{"invalid spin_inversion"};
  }

  auto reduction_step(Tuple current, Tuple other) const -> Tuple {
    auto current_r = current[0];
    auto current_c_re = current[1];
    auto current_c_im = current[2];
    auto current_n = current[3];
    auto r = other[0];
    auto c_re = other[1];
    auto c_im = other[2];
    auto n = other[3];

    auto is_less = r < current_r;
    auto next_r = select(is_less, r, current_r);
    auto next_c_re = select(is_less, c_re, current_c_re);
    auto next_c_im = select(is_less, c_im, current_c_im);
    auto next_n = current_n + n;

    return {next_r, next_c_re, next_c_im, next_n};
  }

  void generate() {
    auto depth = _masks.dim(0).extent();
    auto number_masks = _masks.dim(1).extent();
    auto const chunk_size = 2 * natural_vector_size(type_of<uint64_t>());
    auto number_chunks = number_masks / chunk_size;
    auto number_rest = number_masks - number_chunks * chunk_size;

    Var i{"i"};

    // Apply masks to generate transformed spin configurations
    RDom k{0, depth, "k"};
    Func y_batched{"y_batched"};
    Var j_outer{"j_outer"};
    Var j_inner{"j_inner"};
    y_batched(i, j_outer, j_inner) = _x(i);
    y_batched(i, j_outer, j_inner) = bit_permute_step_64(
        y_batched(i, j_outer, j_inner),
        _masks(k, j_outer * chunk_size + j_inner), _shifts(k));

    Func y_scalar{"y_scalar"};
    Var j_tail{"j_tail"};
    y_scalar(i, j_tail) = _x(i);
    y_scalar(i, j_tail) =
        bit_permute_step_64(y_scalar(i, j_tail), _masks(k, j_tail), _shifts(k));

    // Compute vectorized reduction
    RDom m_main{0, number_chunks, "m_main"};
    Func custom_reduction{"custom_reduction"};
    custom_reduction(i, j_inner) =
        Tuple{_x(i), cast<double>(1), cast<double>(0), cast<double>(0)};
    custom_reduction(i, j_inner) = reduction_step(
        i, custom_reduction(i, j_inner), y_batched(i, m_main, j_inner),
        _eigvals_re(m_main * chunk_size + j_inner),
        _eigvals_im(m_main * chunk_size + j_inner));

    // Reduce across one chunk
    RDom m_lane{1, chunk_size - 1, "m_lane"};
    Func reduction_lane{"reduction_lane"};
    reduction_lane(i) =
        Tuple{custom_reduction(i, 0)[0], custom_reduction(i, 0)[1],
              custom_reduction(i, 0)[2], custom_reduction(i, 0)[3]};
    reduction_lane(i) =
        reduction_step(reduction_lane(i), custom_reduction(i, m_lane));

    // Compute reduction over remaining elements
    Func reduction_scalar{"reduction_scalar"};
    RDom m_tail{number_chunks * chunk_size, number_rest, "m_tail"};
    reduction_scalar(i) = Tuple{reduction_lane(i)[0], reduction_lane(i)[1],
                                reduction_lane(i)[2], reduction_lane(i)[3]};
    reduction_scalar(i) =
        reduction_step(i, reduction_scalar(i), y_scalar(i, m_tail),
                       _eigvals_re(m_tail), _eigvals_im(m_tail));

    // Store results
    Var q;
    _repr(i) = reduction_scalar(i)[0];
    _character(i, q) = undef<double>();
    _character(i, 0) = reduction_scalar(i)[1];
    _character(i, 1) = reduction_scalar(i)[2];
    if (_spin_inversion == 0) {
      _norm(i) = sqrt(reduction_scalar(i)[3] / number_masks);
    } else {
      _norm(i) = sqrt(reduction_scalar(i)[3] / (2 * number_masks));
    }

    auto batch_size = _x.dim(0).extent();
    _x.dim(0).set_min(0).set_stride(1);
    _masks.dim(0).set_min(0).set_stride(number_masks);
    _masks.dim(1).set_min(0).set_stride(1);
    _eigvals_re.dim(0).set_min(0).set_stride(1).set_extent(number_masks);
    _eigvals_im.dim(0).set_min(0).set_stride(1).set_extent(number_masks);
    _shifts.dim(0).set_min(0).set_stride(1).set_extent(depth);
    _repr.dim(0).set_min(0).set_stride(1).set_extent(batch_size);
    _character.dim(0).set_min(0).set_stride(2).set_extent(batch_size);
    _character.dim(1).set_min(0).set_stride(1).set_extent(2);
    _norm.dim(0).set_min(0).set_stride(1).set_extent(batch_size);

    // Schedule
    y_batched.compute_at(custom_reduction, m_main);

    custom_reduction.compute_at(reduction_scalar, i);
    custom_reduction.vectorize(j_inner);
    custom_reduction.update(0).vectorize(j_inner);
    // reduction_lane.compute_at(reduction_scalar, i);
    reduction_scalar.compute_at(_repr, i);
    _character.update(0).compute_with(_repr, i);
    _character.update(1).compute_with(_repr, i);
    _norm.compute_with(_repr, i);
  }
};

class is_representative_generator
    : public Halide::Generator<is_representative_generator> {
public:
  GeneratorParam<int> _spin_inversion{"spin_inversion", /*default value*/ 0};

  Input<Buffer<uint64_t>> _x{"x", 1};
  Input<uint64_t> _flip_mask{"flip_mask"};
  Input<Buffer<uint64_t>> _masks{"masks", 2};
  Input<Buffer<double>> _eigvals_re{"eigvals_re", 1};
  Input<Buffer<uint64_t>> _shifts{"shifts", 1};
  Output<Buffer<uint8_t>> _is_representative{"is_representative", 1};
  Output<Buffer<double>> _norm{"norm", 1};

  void generate() {
    auto count = _x.dim(0).extent();
    auto depth = _masks.dim(0).extent();
    auto number_masks = _masks.dim(1).extent();
    auto const chunk_size = natural_vector_size(type_of<uint64_t>());
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
        y(i, j_outer, j_inner), _masks(k, j_outer * chunk_size + j_inner),
        _shifts(k));
    Func temp{"temp"};
    if (_spin_inversion == 0) {
      temp(i, j_outer, j_inner) =
          Tuple{cast<double>(y(i, j_outer, j_inner) == _x(i)) *
                    _eigvals_re(j_outer * chunk_size + j_inner),
                y(i, j_outer, j_inner) >= _x(i)};
    } else if (_spin_inversion == 1) {
      auto y_flipped = y(i, j_outer, j_inner) ^ _flip_mask;
      temp(i, j_outer, j_inner) =
          Tuple{(cast<double>(y(i, j_outer, j_inner) == _x(i)) +
                 cast<double>(y_flipped == _x(i))) *
                    _eigvals_re(j_outer * chunk_size + j_inner),
                y(i, j_outer, j_inner) >= _x(i) && y_flipped >= _x(i)};
    } else if (_spin_inversion == -1) {
      auto y_flipped = y(i, j_outer, j_inner) ^ _flip_mask;
      temp(i, j_outer, j_inner) =
          Tuple{(cast<double>(y(i, j_outer, j_inner) == _x(i)) -
                 cast<double>(y_flipped == _x(i))) *
                    _eigvals_re(j_outer * chunk_size + j_inner),
                y(i, j_outer, j_inner) >= _x(i) && y_flipped >= _x(i)};
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
    batch_reduction(i) =
        Tuple{batch_reduction(i)[0] + lane_reduction(i, m_batch)[0],
              batch_reduction(i)[1] && lane_reduction(i, m_batch)[1]};

    // Compute reduction over remaining elements
    Func scalar_reduction{"scalar_reduction"};
    RDom m_tail{0, number_rest, "m_tail"};
    scalar_reduction(i) = Tuple{batch_reduction(i)[0], batch_reduction(i)[1]};
    // m_tail.where(batch_reduction(i)[1]);
    scalar_reduction(i) =
        Tuple{scalar_reduction(i)[0] + temp(i, number_chunks, m_tail)[0],
              scalar_reduction(i)[1] && temp(i, number_chunks, m_tail)[1]};

    // Save results
    _norm(i) = scalar_reduction(i)[0];
    _is_representative(i) = select(scalar_reduction(i)[0] > 0,
                                   cast<uint8_t>(scalar_reduction(i)[1]), 0);

    // Shapes & strides
    _x.dim(0).set_min(0).set_stride(1);
    _masks.dim(0).set_min(0).set_stride(number_masks);
    _masks.dim(1).set_min(0).set_stride(1);
    _shifts.dim(0).set_min(0).set_stride(1).set_extent(depth);
    _eigvals_re.dim(0).set_min(0).set_stride(1).set_extent(number_masks);
    _is_representative.dim(0).set_min(0).set_stride(1).set_extent(count);
    _norm.dim(0).set_min(0).set_stride(1).set_extent(count);

    // Schedule
    temp.in(lane_reduction).vectorize(j_inner);
    temp.in(lane_reduction).compute_at(lane_reduction, j_outer);
    temp.in(scalar_reduction).compute_at(scalar_reduction, i);
    // lane_reduction.compute_at(batch_reduction, m_batch);
    scalar_reduction.compute_at(_is_representative, i);
    _norm.compute_with(_is_representative, i);
  }
};

#if 0
class apply_off_diagonal_generator
    : public Halide::Generator<apply_off_diagonal_generator> {
public:
  GeneratorParam<int> _number_words{"number_words", /*default value*/ 1};

  Input<Buffer<double>> _v{"v", 2};
  Input<Buffer<uint64_t>> _m{"m", 2};
  Input<Buffer<uint64_t>> _r{"r", 2};
  Input<Buffer<uint64_t>> _x{"x", 2};
  Input<Buffer<uint64_t>> _s{"s", 2};

  Input<Buffer<uint64_t>> _alpha{"alpha", 2};
  Output<Buffer<uint64_t>> _beta{"beta", 3};
  Output<Buffer<double>> _coeff{"coeff", 3};

  void _set_shapes_and_strides() {
    auto const number_terms = _v.dim(1).extent();
    auto const batch_size = _alpha.dim(1).extent();
    _v.dim(0).set_min(0).set_stride(1).set_extent(2);
    _v.dim(1).set_min(0).set_stride(2);
    _m.dim(0).set_min(0).set_stride(1).set_extent(_number_words);
    _m.dim(1).set_min(0).set_stride(_number_words).set_extent(number_terms);
    _r.dim(0).set_min(0).set_stride(1).set_extent(_number_words);
    _r.dim(1).set_min(0).set_stride(_number_words).set_extent(number_terms);
    _x.dim(0).set_min(0).set_stride(1).set_extent(_number_words);
    _x.dim(1).set_min(0).set_stride(_number_words).set_extent(number_terms);
    _s.dim(0).set_min(0).set_stride(1).set_extent(_number_words);
    _s.dim(1).set_min(0).set_stride(_number_words).set_extent(number_terms);
    _alpha.dim(0).set_min(0).set_stride(1).set_extent(_number_words);
    _alpha.dim(1).set_min(0).set_stride(_number_words);
    _beta.dim(0).set_min(0).set_stride(1).set_extent(_number_words);
    _beta.dim(1).set_min(0).set_stride(_number_words).set_extent(number_terms);
    _beta.dim(2)
        .set_min(0)
        .set_stride(_number_words * number_terms)
        .set_extent(batch_size);
    _coeff.dim(0).set_min(0).set_stride(1).set_extent(2);
    _coeff.dim(1).set_min(0).set_stride(2).set_extent(number_terms);
    _coeff.dim(2)
        .set_min(0)
        .set_stride(2 * number_terms)
        .set_extent(batch_size);
  }

  void generate() {
    Var batch_idx{"batch_idx"};
    Var term_idx{"term_idx"};
    Var word_idx{"word_idx"};
    Var re_im_idx{"re_im_idx"};
    RDom k{0, _number_words, "k"};

    Func delta{"delta"};
    if (_number_words == 1) {
      delta(term_idx, batch_idx) = cast<double>(
          (_alpha(0, batch_idx) & _m(0, term_idx)) == _r(0, term_idx));
    } else {
      delta(term_idx, batch_idx) = cast<bool>(true);
      delta(term_idx, batch_idx) =
          delta(term_idx, batch_idx) &&
          (_alpha(k, batch_idx) & _m(k, term_idx)) == _r(k, term_idx);
    }

    Func count{"count"};
    if (_number_words == 1) {
      count(term_idx, batch_idx) =
          popcount(_alpha(0, batch_idx) & _s(0, term_idx));
    } else {
      count(term_idx, batch_idx) = cast<int64_t>(0);
      count(term_idx, batch_idx) +=
          popcount(_alpha(k, batch_idx) & _s(k, term_idx));
    }

    Func sign{"sign"};
    sign(term_idx, batch_idx) = 1 - 2 * (count(term_idx, batch_idx) % 2);

    _coeff(re_im_idx, term_idx, batch_idx) = _v(re_im_idx, term_idx) *
                                             delta(term_idx, batch_idx) *
                                             sign(term_idx, batch_idx);
    _beta(word_idx, term_idx, batch_idx) =
        _alpha(word_idx, batch_idx) ^ _x(word_idx, term_idx);

    // Shapes & strides
    _coeff.bound(re_im_idx, 0, 2);
    _beta.bound(word_idx, 0, _number_words);
    _set_shapes_and_strides();

    // Scheduling
    if (!auto_schedule) {
      Var batch_idxi("batch_idxi");
      Var re_im_idxi("re_im_idxi");
      Var term_idxi("term_idxi");
      Var word_idxi("word_idxi");

      if (_number_words == 1) {
        auto const number_terms = _v.dim(1).extent();
        auto const batch_size = _alpha.dim(1).extent();
        _beta.compute_root()
            .specialize(number_terms >= 32 && batch_size >= 32)
            .split(term_idx, term_idx, term_idxi, 32,
                   TailStrategy::ShiftInwards)
            .split(batch_idx, batch_idx, batch_idxi, 32,
                   TailStrategy::ShiftInwards)
            .reorder({word_idx, term_idxi, batch_idxi, term_idx, batch_idx})
            .fuse(word_idx, term_idxi, word_idx)
            .vectorize(word_idx, 8);

        _coeff.vectorize(re_im_idx);
        _coeff.specialize(number_terms >= 8 && batch_size >= 8)
            .split(term_idx, term_idx, term_idxi, 8, TailStrategy::ShiftInwards)
            .split(batch_idx, batch_idx, batch_idxi, 8,
                   TailStrategy::ShiftInwards)
            .reorder({re_im_idx, term_idxi, batch_idxi, batch_idx, term_idx});
        sign.store_in(MemoryType::Stack)
            .split(term_idx, term_idx, term_idxi, 8, TailStrategy::GuardWithIf)
            .split(batch_idx, batch_idx, batch_idxi, 8,
                   TailStrategy::GuardWithIf)
            .vectorize(batch_idxi)
            .compute_at(_coeff, batch_idx)
            .reorder({batch_idxi, term_idxi, batch_idx, term_idx})
            .reorder_storage(batch_idx, term_idx);
        delta.store_in(MemoryType::Stack)
            .split(term_idx, term_idx, term_idxi, 8, TailStrategy::GuardWithIf)
            .split(batch_idx, batch_idx, batch_idxi, 8,
                   TailStrategy::GuardWithIf)
            .vectorize(batch_idxi)
            .compute_at(_coeff, batch_idx)
            .reorder({batch_idxi, term_idxi, batch_idx, term_idx})
            .reorder_storage(batch_idx, term_idx);
      }
    }
  }

  void schedule() {
    if (auto_schedule) {
      auto const number_terms = 100;
      auto const batch_size = 1000;
      _v.set_estimates({{0, 2}, {0, number_terms}});
      _m.set_estimates({{0, _number_words}, {0, number_terms}});
      _r.set_estimates({{0, _number_words}, {0, number_terms}});
      _x.set_estimates({{0, _number_words}, {0, number_terms}});
      _s.set_estimates({{0, _number_words}, {0, number_terms}});
      _alpha.set_estimates({{0, _number_words}, {0, batch_size}});
      _beta.set_estimates(
          {{0, _number_words}, {0, number_terms}, {0, batch_size}});
      _coeff.set_estimates({{0, 2}, {0, number_terms}, {0, batch_size}});
    }
  }
};
#endif

// void ls_hs_operator_apply_off_diag_kernel(
//     ls_hs_operator const *op, ptrdiff_t batch_size, uint64_t const *alphas,
//     ptrdiff_t alphas_stride, uint64_t *betas, ptrdiff_t betas_stride,
//     ls_hs_scalar *coeffs) {
//   // fprintf(stderr, "ls_hs_operator_apply_off_diag_kernel ...\n");
//   if (op->off_diag_terms == NULL ||
//       op->off_diag_terms->number_terms == 0) { // nothing to apply
//     // fprintf(stderr, "Nothing to do...\n");
//     return;
//   }
//
//   int const number_words = (op->off_diag_terms->number_bits + 63) / 64;
//   uint64_t *restrict const temp =
//       malloc((size_t)number_words * sizeof(uint64_t));
//   if (temp == NULL) {
//     // fprintf(stderr, "%s\n", "failed to allocate memory");
//     abort();
//   }
//
//   ls_hs_nonbranching_terms const *restrict terms = op->off_diag_terms;
//   ptrdiff_t const number_terms = terms->number_terms;
//   for (ptrdiff_t batch_idx = 0; batch_idx < batch_size; ++batch_idx) {
//     uint64_t const *restrict const alpha = alphas + batch_idx *
//     alphas_stride; for (ptrdiff_t term_idx = 0; term_idx < number_terms;
//     ++term_idx) {
//       ls_hs_scalar const v = terms->v[term_idx];
//       uint64_t const *restrict const m = terms->m + term_idx * number_words;
//       uint64_t const *restrict const r = terms->r + term_idx * number_words;
//       uint64_t const *restrict const x = terms->x + term_idx * number_words;
//       uint64_t const *restrict const s = terms->s + term_idx * number_words;
//       uint64_t *restrict const beta =
//           betas + (batch_idx * number_terms + term_idx) * betas_stride;
//
//       bitstring_and(number_words, alpha, m, temp);
//       int const delta = bitstring_equal(number_words, temp, r);
//       bitstring_and(number_words, alpha, s, temp);
//       int const sign = 1 - 2 * (bitstring_popcount(number_words, temp) % 2);
//       coeffs[batch_idx * number_terms + term_idx] = v * (delta * sign);
//       bitstring_xor(number_words, alpha, x, beta);
//     }
//   }
//   free(temp);
// }

HALIDE_REGISTER_GENERATOR(state_info_generator, state_info)
HALIDE_REGISTER_GENERATOR(is_representative_generator, is_representative)
// HALIDE_REGISTER_GENERATOR(apply_off_diagonal_generator, apply_off_diagonal)
