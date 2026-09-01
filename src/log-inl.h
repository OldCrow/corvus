// Public log and log1p kernels (#32): thin assemblies over the corvus-owned
// log_dd core. Per-target include guard (Highway -inl.h idiom).
//
// METHOD
//   log:   LogDdAny(x) -- the centred-mantissa table core (~2^-70 relative,
//          budget certified by tools/gen_log_table.py's self-check) with its
//          2^600 subnormal prescale -- then one final rounding fl(hi + lo).
//          Total <= 0.5 + ~2^-17 ulp.
//   log1p: two regimes, split at |x| <= 2^-30 (kLog1pTinyCut).
//          MAIN: s = TwoSum(1, x) captures 1 + x EXACTLY for every double x
//          (Sterbenz makes 1 + x itself exact on [-1, -1/2]; elsewhere the
//          residual carries what the sum rounded away), then LogDdAny's dd
//          overload: log(s.hi) + t*(1 - t/2) with t = s.lo/s.hi.
//          TINY: the dd overload's correction c = t*fl(1 - t/2) is plain
//          double, and for |x| ~ 2^-53 the result IS c: LogDd(1) is exactly
//          zero, so c's own rounding (the factor quantizes to one part in
//          2^53, and the t^2/2 term drops entirely for |t| <= 2^-53) lands
//          UNATTENUATED in a result of the same size -- up to ~0.66 ulp
//          from truth, misrounding up to 21% of the 2^-53 binade (#35 H1;
//          the v0.8.0 header claimed "full relative accuracy" here, a
//          relative-vs-absolute frame slip). Below the cut, compute the
//          series directly in exact-product arithmetic instead:
//          x - x^2/2 + x^3/3, with x^2 exact by TwoProd, the x^2/2
//          subtraction captured by Fast2Sum (|x^2/2| <= 2^-61|x| << |x|),
//          and the residual + x^2 tail + cubic folded into one low word.
//          Error budget at the cut (worst case both regimes): neglected
//          x^4/4 <= 2^-122 relative; series-side rounding ~2^-23 ulp;
//          main-path absolute error ~2^-105 against ulp(2^-30) = 2^-82
//          is ~2^-23 ulp. Total <= 0.5 + ~2^-23 ulp everywhere -- inside
//          the library's CR-except-near-tie standard with margin.
//          NOTE: no representable double x > -1 has a subnormal
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

  // Tiny regime (see header): series in exact products. sq = x^2 exactly;
  // Fast2Sum's precondition |x| >= |sq.hi/2| holds for every |x| <= 1/2
  // (equality only at 0, where every term is a zero). The x^3 coefficient's
  // sq.lo*(-0.5) is an exact halving and x*sq.hi rounds ~2^-53 relative on
  // a term that is itself <= 2^-60 of the result, so both MulAdds are
  // fusion-indifferent: FMA and mul+add tiers agree bit for bit.
  // |x| <= 2^-30, via Ge (the facade carries no Le); NaN lanes compare
  // false and ride the main path into the payload-preserving blend below.
  const auto tiny = op::Ge(op::Set(d, 0x1p-30), op::Abs(x));
  if (!op::AllFalse(d, tiny)) {
    const auto sq = TwoProd(d, x, x);
    const auto mhalf = op::Set(d, -0.5);
    const auto s1 = Fast2Sum(d, x, op::Mul(sq.hi, mhalf));
    auto lo = op::MulAdd(sq.lo, mhalf, s1.lo);
    // fl(1/3)'s own rounding enters ~2^-60-relative -- far below budget.
    lo = op::MulAdd(op::Mul(x, sq.hi), op::Set(d, 1.0 / 3.0), lo);
    res = op::IfThenElse(tiny, op::Add(s1.hi, lo), res);
  }

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
