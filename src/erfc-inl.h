// erfc kernel assembly (#8: hoisted from the dispatch TU so it is
// reachable by a corvus_kernel_test_target() test and reusable from
// another TU, like every other family; the compensated core lives in
// src/erfc_core-inl.h). Per-target include guard (Highway -inl.h idiom).
#if defined(CORVUS_ERFC_INL_H_) == defined(HWY_TARGET_TOGGLE)
#ifdef CORVUS_ERFC_INL_H_
#undef CORVUS_ERFC_INL_H_
#else
#define CORVUS_ERFC_INL_H_
#endif

#include "src/dd-inl.h"
#include "src/erf_core-inl.h"
#include "src/erf_data.h"
#include "src/erfc_core-inl.h"
#include "src/erfc_tail_data.h"
#include "src/erfc_tail_dd-inl.h"  // GammaTailG, the dd G (#74)
#include "src/exp_dd-inl.h"
#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {
namespace op = ops;

// erfc in two regions.
//
// Core, |x| <= 6: from the shared erf table, assembled compensated so the
// cancellation in 1 - erf(x) never loses the low bits:
//   sE = sign(x)*E_hi; hi = 1 - sE (Fast2Sum head, |E_hi| <= 1);
//   lo = (1 - hi) - sE (exact residual); erfc = hi + (lo - sign(x)*small).
// Near the erf saturation region E_hi == 1 and E_lo is exactly -erfc(r) to
// full precision, so relative accuracy survives all the way to 6.
//
// Tail, |x| > 6: erfc(a) = e^{-a^2} * G(a) / a, G(a) = a e^{a^2} erfc(a).
// The squared argument is split exactly (ssq + sl = a*a via ops::SquareLow)
// because a 1/2-ULP error in a^2 alone would be amplified to ~a^2 * 2^-53
// relative error by the exponential -- ~360 ULP at a = 27. The split pair is
// fed to corvus's own exp_dd, which consumes both halves inside its argument
// reduction; the earlier Exp(-ssq)*(1 - sl) first-order correction is gone
// with it, as is the backend Exp that used to set this region's error bound
// (5 ULP, ~59% not correctly rounded).
//
// G IS THE SHARED dd FIT (#74). It used to be erfc's own plain-double G
// (ErfcTailGFromU, ~1.7 x 2^-53 relative: fit 0.81 + Horner rounding 0.89),
// and that error set the tail's accuracy: up to 2.07 ulp real on normal
// results (the 2-ULP bound) and 1.45 in the top subnormal binade, against a
// 1-ULP bound there. GammaTailG (src/erfc_tail_dd-inl.h) is the same
// function as a dd to 2^-64 on z in [5, 28.5], built for the Temme tails
// (#64, #69); erfc's live tail, a in (6, 28], lies inside it. With it every
// sampled tail row is correctly rounded or a near-tie (max 0.50 real, three
// tiers). Cost: the dd Horner and products take the tail-only bench from
// ~5.5 to ~9.8 ns/el (Zen 4, indicative); the core is untouched.
// ErfcTailGFromU stays for erfinv's far Halley step, which needs only
// double accuracy.
//
// Everything after the exponential is assembled in double-double and rounded
// exactly once, at the end:
//   1/a as a dd (DdRecip) -- a rounded 1/a would put its own 1/2 ULP straight
//     into the result, since erfc is proportional to it; r = 1/(2a), G's
//     argument, is its exact half.
//   exp_dd in mantissa/exponent form, and the product m times 2^e rounded by
//     DdScaleToDouble. Results are subnormal for a in (26.5433, 27.22602)
//     (erfc(a) < 2^-1022 from 26.5433; below half of 2^-1074 past 27.22602).
//     Before #74 the dd was rounded to double and THEN scaled, which rounds
//     twice on that band; DdScaleToDouble rounds once.
// Clamping a to 28 (past the erfc underflow point) keeps inf lanes out of the
// inf*0 = NaN trap in the exact split, and the lower clamp at 5 keeps core
// lanes that share the vector inside G's domain (their result is discarded).
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
  CORVUS_ROUTE_SCOPE(d, op::Gt(ax, op::Set(d, detail::kErfcTailLo)));
  CORVUS_ROUTE(d, "erfc/tail_clamp", op::Gt(ax, op::Set(d, detail::kErfcTailHi)));
  const auto at = op::Max(op::Min(ax, op::Set(d, detail::kErfcTailHi)),
                          op::Set(d, detail::kGammaTailZLo));
  const auto ssq = op::Mul(at, at);
  const auto sl = op::SquareLow(d, at, ssq);  // exact: at^2 = ssq + sl
  const auto ur = DdRecip(d, at);             // 1/at to ~2^-105
  // G as a dd from the shared Temme fit at r = 1/(2a), exactly half of ur.
  const Dd<D> r2{op::Mul(ur.hi, op::Set(d, 0.5)), op::Mul(ur.lo, op::Set(d, 0.5))};
  const auto gdd = GammaTailG(d, r2);

  const auto ex = ExpDdFrac(d, op::Neg(ssq), op::Neg(sl));
  const auto m = DdMul(d, ex.m, DdMul(d, ur, gdd));  // e^{-a^2} mantissa * G/a
  const auto tail_pos = DdScaleToDouble(d, m, ex.e);
  CORVUS_ROUTE(d, "erfc/tail_underflow", op::Eq(tail_pos, op::Zero(d)));
  {
    CORVUS_ROUTE_SCOPE(d, op::Gt(tail_pos, op::Zero(d)));
    CORVUS_ROUTE(d, "erfc/tail_subnormal",
                 op::Lt(tail_pos, op::Set(d, 0x1p-1022)));
  }
  CORVUS_ROUTE(d, "erfc/tail_negative", op::Lt(x, op::Zero(d)));
  return op::IfThenElse(op::Lt(x, op::Zero(d)),
                        op::Sub(op::Set(d, 2.0), tail_pos), tail_pos);
}

template <class D>
static HWY_INLINE op::V<D> ErfcVec(D d, op::V<D> x) {
  const auto ax = op::Abs(x);
  const auto nan = op::IsNaN(x);
  // NaN lanes have tail_m false and are handled in the core branch.
  const auto tail_m = op::Gt(ax, op::Set(d, detail::kErfcTailLo));

  // Real workloads are usually single-region per vector (Gaussian CDFs
  // essentially never leave |x| <= 6), so skip the unused path: the tail's
  // Div+Exp+Horner roughly halves core-only throughput if always computed.
  CORVUS_ROUTE(d, "erfc/core", op::Ge(op::Set(d, detail::kErfcTailLo), ax));
  CORVUS_ROUTE(d, "erfc/tail", tail_m);
  CORVUS_ROUTE(d, "erfc/special/nan", nan);

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

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#endif  // include guard
