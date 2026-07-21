#include "corvus/corvus.h"
#include "src/erf_data.h"
#include "src/erfc_tail_data.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "src/erfc.cpp"
#include "hwy/foreach_target.h"

#include "src/erf_core-inl.h"
#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

// Select one of three per-interval constants by the tail-interval masks.
template <class D, class M>
HWY_INLINE op::V<D> Sel3(D d, M m1, M m2, double v0, double v1, double v2) {
  return op::IfThenElse(m1, op::Set(d, v0),
                        op::IfThenElse(m2, op::Set(d, v1), op::Set(d, v2)));
}

// erfc in two regions.
//
// Core, |x| <= 6: from the shared erf table, assembled compensated so the
// cancellation in 1 - erf(x) never loses the low bits:
//   sE = sign(x)*E_hi; hi = 1 - sE (Fast2Sum head, |E_hi| <= 1);
//   lo = (1 - hi) - sE (exact residual); erfc = hi + (lo - sign(x)*small).
// Near the erf saturation region E_hi == 1 and E_lo is exactly -erfc(r) to
// full precision, so relative accuracy survives all the way to 6.
//
// Tail, |x| > 6: erfc(a) = e^{-a^2} * G(1/a) * (1/a) with G fitted per
// interval (see tools/gen_erfc_tail_poly.py). The squared argument is split
// exactly (ssq + sl = a*a via FMA) because a 1/2-ULP error in a^2 alone
// would be amplified to ~a^2 * 2^-53 relative error by the exponential --
// ~360 ULP at a = 27. e^{-ssq-sl} = Exp(-ssq)*(1 - sl) to O(sl^2).
// Accuracy here is gated by the backend Exp; measured by test_erfc_ulp.
// Clamping a to 28 (past the erfc underflow point ~27.3) keeps inf lanes
// out of the inf*0 = NaN trap in the exact split.
//
// x < 0 mirrors via erfc(x) = 2 - erfc(|x|); in the tail erfc(|x|) < 2^-52
// so the subtraction rounds to exactly 2, matching erfc's saturation.
template <class D, class M>
HWY_INLINE op::V<D> ErfcCoreVec(D d, op::V<D> x, op::V<D> ax, M nan) {
  const auto one = op::Set(d, 1.0);
  const auto sgn = op::CopySign(one, x);  // +/-1

  auto ac = op::Min(ax, op::Set(d, 6.0));
  ac = op::IfThenElse(nan, op::Zero(d), ac);  // safe table index for NaN lanes
  const auto parts = ErfTableCore(d, ac);

  const auto sE = op::Mul(sgn, parts.e_hi);
  const auto hi = op::Sub(one, sE);
  const auto lo = op::Sub(op::Sub(one, hi), sE);  // exact Fast2Sum residual
  return op::Add(hi, op::Sub(lo, op::Mul(sgn, parts.small)));
}

template <class D>
HWY_INLINE op::V<D> ErfcTailVec(D d, op::V<D> x, op::V<D> ax) {
  const auto one = op::Set(d, 1.0);

  const auto at = op::Min(ax, op::Set(d, 28.0));
  const auto ssq = op::Mul(at, at);
  const auto sl = op::SquareLow(d, at, ssq);  // exact: at^2 = ssq + sl
  const auto u = op::Div(one, at);

  const auto m1 = op::Lt(at, op::Set(d, detail::kErfcTailBound1));
  const auto m2 = op::Lt(at, op::Set(d, detail::kErfcTailBound2));
  const auto scale = Sel3(d, m1, m2, detail::kErfcTailScale[0],
                          detail::kErfcTailScale[1], detail::kErfcTailScale[2]);
  const auto shift = Sel3(d, m1, m2, detail::kErfcTailShift[0],
                          detail::kErfcTailShift[1], detail::kErfcTailShift[2]);
  const auto s = op::MulAdd(u, scale, shift);

  const auto* c = detail::kErfcTailCoef;
  auto poly = Sel3(d, m1, m2, c[0][detail::kErfcTailNCoef - 1],
                   c[1][detail::kErfcTailNCoef - 1],
                   c[2][detail::kErfcTailNCoef - 1]);
  for (int k = detail::kErfcTailNCoef - 2; k >= 0; --k) {
    poly = op::MulAdd(poly, s, Sel3(d, m1, m2, c[0][k], c[1][k], c[2][k]));
  }

  auto ex = op::Exp(d, op::Neg(ssq));
  ex = op::Mul(ex, op::Sub(one, sl));
  const auto tail_pos = op::Mul(op::Mul(ex, poly), u);
  return op::IfThenElse(op::Lt(x, op::Zero(d)),
                        op::Sub(op::Set(d, 2.0), tail_pos), tail_pos);
}

template <class D>
HWY_INLINE op::V<D> ErfcVec(D d, op::V<D> x) {
  const auto ax = op::Abs(x);
  const auto nan = op::IsNaN(x);
  // NaN lanes have tail_m false and are handled in the core branch.
  const auto tail_m = op::Gt(ax, op::Set(d, 6.0));

  // Real workloads are usually single-region per vector (Gaussian CDFs
  // essentially never leave |x| <= 6), so skip the unused path: the tail's
  // Div+Exp+Horner roughly halves core-only throughput if always computed.
  op::V<D> res;
  if (op::AllFalse(d, tail_m)) {
    res = ErfcCoreVec(d, x, ax, nan);
  } else if (op::AllTrue(d, tail_m)) {
    res = ErfcTailVec(d, x, ax);
  } else {
    res = op::IfThenElse(tail_m, ErfcTailVec(d, x, ax),
                         ErfcCoreVec(d, x, ax, nan));
  }
  res = op::IfThenElse(nan, x, res);  // propagate NaN (payload preserved)
  return res;
}

void ErfcImpl(const double* in, double* out, size_t n) {
  const op::ScalableTag<double> d;
  const size_t N = op::Lanes(d);
  size_t i = 0;
  for (; i + N <= n; i += N) {
    op::Store(ErfcVec(d, op::Load(d, in + i)), d, out + i);
  }
  if (i < n) {
    // Masked tail: same code path as the full lanes, no scalar fallback.
    op::StoreN(ErfcVec(d, op::LoadN(d, in + i, n - i)), d, out + i, n - i);
  }
}

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace corvus {

HWY_EXPORT(ErfcImpl);

void Erfc(std::span<const double> in, std::span<double> out) {
  HWY_DYNAMIC_DISPATCH(ErfcImpl)(in.data(), out.data(), in.size());
}

}  // namespace corvus
#endif  // HWY_ONCE
