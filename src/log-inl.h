// Public log and log1p kernels (#32): thin assemblies over the corvus-owned
// log_dd core. Per-target include guard (Highway -inl.h idiom).
//
// METHOD
//   log:   LogDdAny(x) -- the centred-mantissa table core (~2^-70 relative,
//          budget certified by tools/gen_log_table.py's self-check) with its
//          2^600 subnormal prescale -- then one final rounding fl(hi + lo).
//          Total <= 0.5 + ~2^-17 ulp.
//   log1p: s = TwoSum(1, x) captures 1 + x EXACTLY for every double x
//          (Sterbenz makes 1 + x itself exact on [-1, -1/2]; elsewhere the
//          residual carries what the sum rounded away), then LogDdAny's dd
//          overload: log(s.hi) + t*(1 - t/2) with t = s.lo/s.hi. A TwoSum
//          pair is normalized, so |t| <= 2^-52 always and the kept
//          quadratic term suffices (dropped cubic ~2^-106 relative); for
//          tiny x the path degenerates to x - x^2/2 with full relative
//          accuracy. NOTE: no representable double x > -1 has a subnormal
//          1 + x (the grid bottoms out at 1 + x = 2^-53 exactly, at
//          nextafter(-1, 0)), so LogDdAny's prescale is never exercised
//          from log1p -- routed through it anyway for one shared shape.
//
// SPECIALS (explicit blends; discarded lanes still execute the core, whose
// slot index is bit-masked into table range, so any input is gather-safe):
//   log:   +-0 -> -inf; x < 0 -> qNaN; +inf -> +inf; NaN propagates.
//   log1p: x == -1 -> -inf; x < -1 -> qNaN; +inf -> +inf; NaN propagates;
//          +-0 -> x (the core would return +0 for -0: RN(+0 + -0) = +0
//          drops the sign, so the sign-of-zero contract is a blend here,
//          unlike sin's by-construction sign).
#if defined(CORVUS_LOG_INL_H_) == defined(HWY_TARGET_TOGGLE)
#ifdef CORVUS_LOG_INL_H_
#undef CORVUS_LOG_INL_H_
#else
#define CORVUS_LOG_INL_H_
#endif

#include <limits>

#include "src/dd-inl.h"
#include "src/log_dd-inl.h"
#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

template <class D>
HWY_INLINE op::V<D> LogVec(D d, op::V<D> x) {
  const Dd<D> r = OutlinedLogDd(d, x);  // shared HWY_NOINLINE wrapper
  auto res = op::Add(r.hi, r.lo);

  const auto zero = op::Zero(d);
  const auto inf = op::Set(d, std::numeric_limits<double>::infinity());
  const auto qnan = op::Set(d, std::numeric_limits<double>::quiet_NaN());
  res = op::IfThenElse(op::Eq(x, zero), op::Neg(inf), res);  // both zeros
  res = op::IfThenElse(op::Lt(x, zero), qnan, res);
  res = op::IfThenElse(op::Eq(x, inf), inf, res);
  return op::IfThenElse(op::IsNaN(x), x, res);  // payload-preserving
}

template <class D>
HWY_INLINE op::V<D> Log1pVec(D d, op::V<D> x) {
  const auto one = op::Set(d, 1.0);
  const auto s = TwoSum(d, one, x);  // exact 1 + x for every double
  const Dd<D> r = OutlinedLogDd(d, Dd<D>{s.hi, s.lo});
  auto res = op::Add(r.hi, r.lo);

  const auto zero = op::Zero(d);
  const auto inf = op::Set(d, std::numeric_limits<double>::infinity());
  const auto qnan = op::Set(d, std::numeric_limits<double>::quiet_NaN());
  const auto neg1 = op::Set(d, -1.0);
  res = op::IfThenElse(op::Eq(x, neg1), op::Neg(inf), res);
  res = op::IfThenElse(op::Lt(x, neg1), qnan, res);
  res = op::IfThenElse(op::Eq(x, inf), inf, res);
  res = op::IfThenElse(op::Eq(x, zero), x, res);  // log1p(+-0) = +-0
  return op::IfThenElse(op::IsNaN(x), x, res);  // payload-preserving
}

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#endif  // include guard
