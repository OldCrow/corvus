#include "corvus/corvus.h"
#include "src/erf_data.h"
#include "src/erfc_tail_data.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "src/erfc.cpp"
#include "hwy/foreach_target.h"

#include "src/dd-inl.h"
#include "src/erf_core-inl.h"
#include "src/erfc_core-inl.h"
#include "src/exp_dd-inl.h"
#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

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
// exactly (ssq + sl = a*a via ops::SquareLow) because a 1/2-ULP error in a^2
// alone would be amplified to ~a^2 * 2^-53 relative error by the exponential
// -- ~360 ULP at a = 27. The split pair is fed to corvus's own exp_dd, which
// consumes both halves inside its argument reduction; the earlier
// Exp(-ssq)*(1 - sl) first-order correction is gone with it, as is the
// backend Exp that used to set this region's error bound (5 ULP, ~59% not
// correctly rounded).
//
// Everything after the exponential is assembled in double-double and rounded
// exactly once, at the end:
//   1/a as a dd (DdRecip) -- a rounded 1/a would put its own 1/2 ULP straight
//     into the result, since erfc is proportional to it. Its high word alone
//     is still fine as the POLYNOMIAL's argument: G is nearly flat there
//     (u*G'/G = O(u^2) <= 0.03), so u's rounding is attenuated ~30x.
//   exp_dd in mantissa/exponent form, so the power-of-two scaling happens
//     after the multiplications -- that is what keeps the subnormal band
//     (a > ~26.6) to a single rounding rather than one per factor.
// Clamping a to 28 (past the erfc underflow point ~27.3) keeps inf lanes
// out of the inf*0 = NaN trap in the exact split.
//
// x < 0 mirrors via erfc(x) = 2 - erfc(|x|); in the tail erfc(|x|) < 2^-52
// so the subtraction rounds to exactly 2, matching erfc's saturation.
template <class D, class M>
static HWY_INLINE op::V<D> ErfcCoreVec(D d, op::V<D> x, op::V<D> ax, M nan) {
  // Safe table index for NaN lanes; ErfcCoreDd does its own <= 6 clamp.
  const auto ax_safe = op::IfThenElse(nan, op::Zero(d), ax);
  const auto pair = ErfcCoreDd(d, x, ax_safe);
  return op::Add(pair.hi, pair.lo);
}

template <class D>
static HWY_INLINE op::V<D> ErfcTailVec(D d, op::V<D> x, op::V<D> ax) {
  const auto at = op::Min(ax, op::Set(d, 28.0));
  const auto ssq = op::Mul(at, at);
  const auto sl = op::SquareLow(d, at, ssq);  // exact: at^2 = ssq + sl
  const auto ur = DdRecip(d, at);             // 1/at to ~2^-105
  const auto u = ur.hi;
  const auto poly = ErfcTailGFromU(d, at, u);

  const auto ex = ExpDdFrac(d, op::Neg(ssq), op::Neg(sl));
  const auto m = DdMul(d, ex.m, DdMulD(d, ur, poly));  // e^{-a^2} mantissa * G/a
  const auto tail_pos = ScaleTwo(d, DdToDouble(m), ex.e);
  return op::IfThenElse(op::Lt(x, op::Zero(d)),
                        op::Sub(op::Set(d, 2.0), tail_pos), tail_pos);
}

template <class D>
static HWY_INLINE op::V<D> ErfcVec(D d, op::V<D> x) {
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

static void ErfcImpl(const double* in, double* out, size_t n) {
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

// The dispatch calls stay INSIDE namespace corvus: with a single compiled
// target (the SSE2 cap) Highway collapses HWY_DYNAMIC_DISPATCH to
// N_SSE2::FUNC, and a globally qualified call would then name a namespace
// that does not exist. It compiles at every other tier, so only the cap
// sweep catches it.
void erfc(std::span<const double> in, std::span<double> out) noexcept {
  HWY_DYNAMIC_DISPATCH(ErfcImpl)(in.data(), out.data(), in.size());
}

}  // namespace corvus
#endif  // HWY_ONCE
