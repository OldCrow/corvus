// erf kernel assembly (#8: hoisted from the dispatch TU so it is
// reachable by a corvus_kernel_test_target() test and reusable from
// another TU, like every other family; the method itself lives in
// src/erf_core-inl.h). Per-target include guard (Highway -inl.h idiom).
#if defined(CORVUS_ERF_INL_H_) == defined(HWY_TARGET_TOGGLE)
#ifdef CORVUS_ERF_INL_H_
#undef CORVUS_ERF_INL_H_
#else
#define CORVUS_ERF_INL_H_
#endif

#include "src/erf_core-inl.h"
#include "src/erf_data.h"
#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {
namespace op = ops;

// Table cell j = 0 (ac <= 1/512, where Round(ac*256) is 0 with ties to
// even), #80. There E = E_lo = 0 and the core returns RN(S*RN(p)), with S
// = 2/sqrt(pi) rounded: two roundings plus S's own 2^-54. That left ~26%
// of the cell not correctly rounded. Here 2/sqrt(pi) = Shi + Slo (dd),
//   erf(ac) = Shi*ac + (Shi*q + Slo*ac),  q = ac^3 (-1/3 + w/10 - w^2/42),
// w = ac^2; the x^9/216 truncation is <= 2^-79.7 relative. Shi*ac is
// split exactly (hi + lo) and the result rounds once: the correction's
// own error is ~2^-70 relative, so only lanes that close to a midpoint
// can miss. Two scalings:
//   * ac >= 2^-900 (main): unscaled. Shi*ac >= 2^-900 sits above the
//     ProdLow floor and Slo*ac >= 2^-956 is normal.
//   * ac < 2^-900 (tiny): q is below 2^-1800 relative, so drop it. Scale
//     to k = ac*2^1074, an exact power-of-two shift (k < 2^174). The
//     result's grid is then the integers while S*k < 2^52 (subnormal
//     outputs), so round S*k to an integer from its dd form, ties
//     impossible (S is irrational). Above that, round to 53 bits; the
//     scale back by 2^-1074 is exact both ways.
// Outlined (rare path behind ErfVec's AllFalse skip).
template <class D>
HWY_NOINLINE op::V<D> ErfCellZeroVec(D d, op::V<D> ac) {
  const auto shi = op::Set(d, 0x1.20dd750429b6dp+0);
  const auto slo = op::Set(d, 0x1.1ae3a914fed80p-56);
  const auto half = op::Set(d, 0.5);
  const auto one = op::Set(d, 1.0);

  // Main: ac in [2^-900, 2^-9] (tiny lanes ride along and are discarded).
  const auto w = op::Mul(ac, ac);
  const auto poly = op::MulAdd(
      w, op::MulAdd(w, op::Set(d, -1.0 / 42.0), op::Set(d, 0.1)),
      op::Set(d, -1.0 / 3.0));
  const auto q = op::Mul(op::Mul(ac, w), poly);
  const auto hi = op::Mul(shi, ac);
  // FLOOR: main lanes have Shi*ac >= 2^-900; tiny lanes are discarded.
  const auto lo = op::ProdLow(d, shi, ac, hi);
  const auto c = op::MulAdd(shi, q, op::Mul(slo, ac));
  const auto res = op::Add(hi, op::Add(lo, c));

  // Tiny: k = ac*2^1074 in [0, 2^174), in two exact steps since 2^1074
  // overflows (the Min keeps discarded lanes finite).
  const auto k = op::Mul(op::Mul(op::Min(ac, op::Set(d, 0x1p-900)), op::Set(d, 0x1p+600)),
                         op::Set(d, 0x1p+474));
  const auto khi = op::Mul(shi, k);
  // FLOOR: Shi*k is 0 (exact) or >= Shi, since k is a nonzero integer.
  const auto klo = op::Add(op::ProdLow(d, shi, k, khi), op::Mul(slo, k));
  // Integer grid (khi < 2^53): n = Round(khi), and khi - n is exact (|.|
  // <= 1/2 on a common grid). |klo| <= 1/2 + 0.13 (Slo*k), so the
  // remainder r = (khi - n) + klo lies in (-1.13, 1.13), and one step of
  // +-1 finishes the rounding. n <= 2^53,
  // so n*2^-1074 is exact.
  const auto n = op::Round(khi);
  const auto r = op::Add(op::Sub(khi, n), klo);
  auto ni = op::IfThenElse(op::Gt(r, half), op::Add(n, one), n);
  ni = op::IfThenElse(op::Lt(r, op::Neg(half)), op::Sub(n, one), ni);
  // 53-bit grid (khi >= 2^53): khi + klo rounds once to >= 2^52, a normal
  // result after the scale back.
  const auto kint = op::IfThenElse(op::Lt(khi, op::Set(d, 0x1p+53)), ni,
                                   op::Add(khi, klo));
  const auto tiny_res = op::Mul(kint, op::Set(d, 0x1p-1074));

  return op::IfThenElse(op::Lt(ac, op::Set(d, 0x1p-900)), tiny_res, res);
}

// See ErfTableCore (src/erf_core-inl.h) for the method; erf runs the
// correction to d^6 (the d^7/d^8 terms exist for erfc's relative
// precision, ~28% of erf's NEON cost, below 0.012 ulp here). Max 1 ULP on
// validated tiers. |x| >= kErfAMax saturates to +/-1 (the bound is the
// smallest double whose erf rounds to 1, so this is correctly rounded);
// NaN is masked explicitly rather than relying on per-target Min/convert
// NaN semantics.
template <class D>
static CORVUS_OUTLINE op::V<D> ErfVec(D d, op::V<D> x) {
  const auto one = op::Set(d, 1.0);
  const auto amax = op::Set(d, detail::kErfAMax);

  const auto ax = op::Abs(x);
  const auto nan = op::IsNaN(x);
  auto ac = op::Min(ax, amax);
  ac = op::IfThenElse(nan, op::Zero(d), ac);  // safe table index for NaN lanes

  const auto core = ErfTableCore<6>(d, ac);
  auto res = op::Add(core.e_hi, core.small);
  const auto j0 = op::Ge(op::Set(d, 0x1p-9), ac);  // cell j = 0 (and NaN lanes)
  CORVUS_ROUTE(d, "erf/cell0", j0);
  if (!op::AllFalse(d, j0)) {
    CORVUS_ROUTE_SCOPE(d, j0);
    res = op::IfThenElse(j0, ErfCellZeroVec(d, ac), res);
  }

  CORVUS_ROUTE(d, "erf/core", op::Lt(ax, amax));
  CORVUS_ROUTE(d, "erf/saturated", op::Ge(ax, amax));
  CORVUS_ROUTE(d, "erf/special/nan", nan);
  res = op::IfThenElse(op::Ge(ax, amax), one, res);
  res = op::IfThenElse(nan, x, res);      // propagate NaN (payload preserved)
  return op::CopySign(res, x);            // erf is odd; also -0 -> -0
}

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#endif  // include guard
