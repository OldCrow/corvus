// Public exp kernel (#32): a thin assembly over the corvus-owned exp_dd
// core. Per-target include guard (Highway -inl.h idiom).
//
// METHOD
//   ExpDdFrac gives exp(x) in mantissa/exponent form (m.hi + m.lo, e) with
//   ~2^-70 relative error (budget certified by tools/gen_exp_table.py's
//   self-check); the public kernel rounds m to one double and lets ScaleTwo
//   apply 2^e in two exact-then-rounding steps. Total error on a normal
//   result: 0.5 ulp (final rounding) + ~2^-17 ulp (core), i.e. correctly
//   rounded except within ~2^-17 ulp of a tie. On a subnormal result the
//   fl(m.hi + m.lo) rounding (2^-53 relative) composes with ScaleTwo's
//   single subnormal rounding for <= ~0.51 ulp in the output's own ulp.
//
// WHY THERE ARE NO THRESHOLD BLENDS
//   The core's kExpXMax = 1100 clamp bounds |e| <= 1600 (ScaleTwo's
//   precondition) and cannot change any representable result -- and the
//   clamp makes every non-NaN special come out right BY CONSTRUCTION:
//     +inf -> clamped to 1100 -> m*2^2031-ish -> RN overflow -> +inf
//     -inf -> clamped to -1100 -> RN(m*2^-2031-ish) -> +0
//     overflow/underflow at the exact correctly-rounded thresholds:
//       ScaleTwo's second multiply IS the result's rounding, so
//       RN(m*2^e) crosses to inf (or to 0 through the subnormals)
//       exactly where the correctly rounded exp does. The reference set
//       brackets both boundaries with +-160-ulp bit ladders; the ULP
//       gate holds those rows to exact inf / exact +0.
//   Only NaN needs a blend (payload-preserving IfThenElse at the end; the
//   clamp's Min/Max on a NaN lane is unspecified, the core's table gather
//   is index-masked so nothing reads out of bounds).
#if defined(CORVUS_EXP_INL_H_) == defined(HWY_TARGET_TOGGLE)
#ifdef CORVUS_EXP_INL_H_
#undef CORVUS_EXP_INL_H_
#else
#define CORVUS_EXP_INL_H_
#endif

#include "src/exp_dd-inl.h"
#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

// OUTLINED core call [MSVC BUILD-TIME GATE, AGENTS.md]: the driver inlines
// the kernel twice per export (full-vector + masked tail), and ExpDdFrac is
// a table gather plus a polynomial. One consumer, so hosted here rather
// than beside the shared OutlinedExpDd pair.
template <class D>
HWY_NOINLINE ExpDdParts<D> OutlinedExpDdParts(D d, op::V<D> x) {
  return ExpDdFrac(d, x, op::Zero(d));
}

template <class D>
HWY_INLINE op::V<D> ExpVec(D d, op::V<D> x) {
  const ExpDdParts<D> parts = OutlinedExpDdParts(d, x);
  const auto v = op::Add(parts.m.hi, parts.m.lo);
  const auto res = ScaleTwo(d, v, parts.e);
  return op::IfThenElse(op::IsNaN(x), x, res);  // payload-preserving
}

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#endif  // include guard
