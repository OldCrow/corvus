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
//   fl(m.hi + m.lo) rounding composes with ScaleTwo's single subnormal
//   rounding: in the TOP subnormal binade the first rounding is worth up
//   to 0.25 ulp of the output (m in [1,2) against e = -1023), so the
//   construction bound is ~0.75 ulp in the output's own ulp -- 0.7495
//   measured (#35 L4 corrected the earlier ~0.51 figure, which held only
//   for deeper subnormals). Still within the 1-ULP gate: total < 1
//   guarantees at most 1 ULP from correctly rounded, which is what the
//   gate deliberately admits for the two-rounding path.
//
// FAST PATH (#43)
//   The dd core costs ~4x a plain table exp, and almost every lane does not
//   need it: its result only has to be correct to well inside half an ulp
//   to round correctly. ExpFast reuses the core's own tables and reduction
//   constants, carries the result as an unevaluated y + rho, and applies
//   the rounding test of CRlibm (de Dinechin, Lauter, Muller): y is the
//   correctly rounded result if y == fl(y + fl(rho*c)), c = 1 + 2^-11.
//   Why: passing means |rho|*c*(1 - 2^-53) <= g, the half-gap to y's
//   neighbour on rho's side, and g >= |y|*2^-54; so the true value, within
//   eps of y + rho, stays inside y's rounding interval whenever
//   eps <= |y|*2^-54*(c*(1 - 2^-53) - 1), i.e. eps/|y| <= ~2^-65. Lanes that fail the test, and lanes
//   with |x| > 708 (subnormal or saturating results, +-inf, NaN), take
//   the dd path below -- rare enough that exp costs the fast path. Bits are unchanged from the dd path:
//   both are correctly rounded wherever the dd path is (0-ULP gate).
//   Error budget of y + rho against exp(x)/2^e, in absolute units where
//   y is in [0.997, 2.006] (ulp(y) >= 2^-53), on EVERY tier -- ProdLow is
//   exact with or without FMA, and every other rounding is on a term of
//   magnitude <= 2^-17:
//     q(r) evaluation and r^2 (~3 roundings on a 2^-17 term)  2^-67.6
//     tail MulAdds without FMA (products <= 2^-17)           3 x 2^-71
//     polynomial truncation r^7/5040, |r| <= 0.00271        2^-71
//     dropped rlo*r, tl*r^2 cross terms                      2^-70
//     k*L3, the part of ln2/N the reduction cannot hold      2^-77
//     table representation (dd split of 2^(j/N))              2^-106
//   Total < 2^-66.8 absolute, so eps/|y| < 2^-66.7 (|y| >= 0.997): a 3.2x
//   margin under the test's 2^-65. Lanes failing it: ~2^-11 of the
//   uniform ones.
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

#include <limits>

#include "src/dd_special-inl.h"  // Expm1DdSeries (#45)
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

// The dd path: every input, specials included. Outlined on every compiler
// (rare path behind ExpVec's AllTrue skip).
template <class D>
HWY_NOINLINE op::V<D> ExpDdVec(D d, op::V<D> x) {
  const ExpDdParts<D> parts = OutlinedExpDdParts(d, x);
  const auto v = op::Add(parts.m.hi, parts.m.lo);
  const auto res = ScaleTwo(d, v, parts.e);
  return op::IfThenElse(op::IsNaN(x), x, res);  // payload-preserving
}

// Fast path: exp(x) for |x| <= 708 (normal, finite results), with a
// mask of the lanes whose result is certified correctly rounded. Other
// lanes hold unspecified values. See the file header for the budget.
template <class D>
struct ExpFastResult {
  op::V<D> v;
  op::M<D> ok;
};

template <class D>
CORVUS_OUTLINE ExpFastResult<D> ExpFast(D d, op::V<D> x) {
  const op::SignedTag<D> di;

  // Reduction, as in ExpDdFrac: k*L1 is exact (34-bit L1, |k| < 2^18) and
  // x - k*L1 is exact by Sterbenz (k*L1 lies within a factor 1 +- 1/2k of
  // x). r rounds once; rlo recovers that rounding: r1 - r is exact by
  // Sterbenz whenever |r| >= 2|k*L2| (~2^-25), and below that the whole
  // error is <= ulp(r)/2 <= 2^-77.
  const auto kf = op::Round(op::Mul(x, op::Set(d, detail::kExpInvL)));
  const auto ki = op::ConvertToInt(di, kf);  // garbage out of range: masked
  const auto nkf = op::Neg(kf);
  const auto r1 = op::Sub(x, op::Mul(kf, op::Set(d, detail::kExpL1)));
  const auto l2 = op::Set(d, detail::kExpL2);
  const auto r = op::MulAdd(nkf, l2, r1);
  const auto rlo = op::MulAdd(nkf, l2, op::Sub(r1, r));

  // e^r - 1 = r + r^2*q(r), q = 1/2 + r/6 + r^2/24 + r^3/120 + r^4/720.
  auto q = op::Set(d, 1.0 / 720.0);
  q = op::MulAdd(q, r, op::Set(d, 1.0 / 120.0));
  q = op::MulAdd(q, r, op::Set(d, 1.0 / 24.0));
  q = op::MulAdd(q, r, op::Set(d, 1.0 / 6.0));
  q = op::MulAdd(q, r, op::Set(d, 0.5));
  const auto g = op::Mul(op::Mul(r, r), q);  // <= 2^-17.9

  // Index masked into range for every input, as in ExpDdFrac.
  const auto j = op::And(ki, op::Set(di, detail::kExpN - 1));
  const auto th = op::GatherIndex(d, detail::kExpTableHi, j);
  const auto tl = op::GatherIndex(d, detail::kExpTableLo, j);

  // T*(1 + (e^r - 1))*(1 + rlo) with T = th + tl:
  //   th + th*r exactly as y0 + e0 + pl; then the small terms.
  const auto p = op::Mul(th, r);
  // FLOOR: th*r < 2^-968 only for k = 0, |x| < 2^-967 (for k != 0, |r| stays
  // far above it): y = th = 1 and |t| < 2^-966, so a 2^-1073 error in pl
  // moves neither y + t (CR value 1) nor the rounding test.
  const auto pl = op::ProdLow(d, th, r, p);  // exact on every tier
  const auto y0 = op::Add(th, p);
  const auto e0 = op::Sub(p, op::Sub(y0, th));  // Fast2Sum: |th| >= 1 > |p|
  auto t = op::MulAdd(th, rlo, tl);
  t = op::MulAdd(tl, r, t);
  t = op::Add(op::Add(t, pl), e0);
  t = op::MulAdd(th, g, t);  // <= 2^-17
  const auto y = op::Add(y0, t);
  const auto rho = op::Sub(t, op::Sub(y, y0));  // Fast2Sum: exact

  // Rounding test (file header), and the range the fast path covers
  // (|x| <= 708 keeps y*2^e normal and finite; NaN compares false).
  const auto ok = op::And(
      op::Ge(op::Set(d, 708.0), op::Abs(x)),
      op::Eq(op::Add(y, op::Mul(rho, op::Set(d, 1.0 + 0x1p-11))), y));

  // y*2^e, one exact multiply (normal result).
  const auto e = op::ShiftRight<detail::kExpLog2N>(ki);
  const auto scale =
      op::BitCast(d, op::ShiftLeft<52>(op::Add(e, op::Set(di, int64_t{1023}))));
  return {op::Mul(y, scale), ok};
}

template <class D>
HWY_INLINE op::V<D> ExpVec(D d, op::V<D> x) {
  const ExpFastResult<D> f = ExpFast(d, x);
  // Routes (census only): the fast path, the dd fallback by cause (rounding
  // test failed inside |x| <= 708; result-range bands outside it; NaN).
  CORVUS_ROUTE(d, "exp/fast", f.ok);
  CORVUS_ROUTE(d, "exp/dd/round",
               op::And(op::Ge(op::Set(d, 708.0), op::Abs(x)), op::MaskNot(f.ok)));
  CORVUS_ROUTE(d, "exp/dd/edge",
               op::And(op::Gt(op::Abs(x), op::Set(d, 708.0)),
                       op::And(op::Ge(op::Set(d, 0x1.62e42fefa39efp9), x),
                               op::Ge(x, op::Set(d, -0x1.6232bdd7abcd2p9)))));
  CORVUS_ROUTE(d, "exp/dd/overflow", op::Gt(x, op::Set(d, 0x1.62e42fefa39efp9)));
  CORVUS_ROUTE(d, "exp/dd/subnormal",
               op::And(op::Gt(op::Set(d, -0x1.6232bdd7abcd2p9), x),
                       op::Ge(x, op::Set(d, -0x1.74910d52d3051p9))));
  CORVUS_ROUTE(d, "exp/dd/underflow", op::Gt(op::Set(d, -0x1.74910d52d3051p9), x));
  CORVUS_ROUTE(d, "exp/special/nan", op::IsNaN(x));
  if (op::AllTrue(d, f.ok)) return f.v;
  return op::IfThenElse(f.ok, f.v, ExpDdVec(d, x));
}

// ------------------------------------------------------------------------
// Public expm1 (#45): e^x - 1, relative to the result, rounded once.
//
// SERIES, |x| <= kExpm1Cut = 2^-6: Expm1DdSeries (src/dd_special-inl.h),
// ~2^-67.6 relative, then one rounding. Subnormal x return x: the w^2/2
// term is below half an ulp of x (the series' FLOOR note).
//
// EXPONENTIAL, |x| > 2^-6: ExpDdFrac's e^x = m * 2^e (~2^-70 relative),
// and the 1 is subtracted IN m's SCALE, m (-) 2^-e, before anything is
// rounded: one rounding of the dd difference, then an exact power-of-two
// scaling. e^x - 1 and e^x differ by a factor 1 - e^-x, so the core's
// error is amplified by e^x/|e^x - 1| <= ~65 at the cut: ~2^-64 relative,
// so correctly rounded except within ~2^-11 ulp of a tie. 2^-e is exact
// for every e this branch sees (e >= -56 below; ScaleTwo's two steps take
// 2^-e to zero for e > 1074, where it is negligible anyway). Overflow is
// exp's construction: the rounded mantissa times 2^e overflows exactly
// when the correctly rounded result does (|e| <= 1600 by the core's clamp).
//
// x <= -38: e^x < 2^-54.8, below half an ulp of -1, so the result is -1
// exactly; those lanes never reach the core.
//
// Specials: expm1(+-0) = +-0 (the series would turn -0 into +0), expm1(-inf)
// = -1, expm1(+inf) = +inf, NaN propagates (payload preserved).
template <class D>
HWY_INLINE op::V<D> Expm1Vec(D d, op::V<D> x) {
  const auto zero = op::Zero(d);
  const auto one = op::Set(d, 1.0);
  const auto inf = op::Set(d, std::numeric_limits<double>::infinity());
  const auto ax = op::Abs(x);
  const auto m_ser = op::Ge(op::Set(d, kExpm1Cut), ax);  // NaN: false
  const auto m_neg1 = op::Ge(op::Set(d, -38.0), x);       // -inf included
  const auto m_exp = op::And(op::Gt(ax, op::Set(d, kExpm1Cut)),
                             op::Gt(x, op::Set(d, -38.0)));  // +inf included
  CORVUS_ROUTE(d, "expm1/series", m_ser);
  CORVUS_ROUTE(d, "expm1/exp", m_exp);
  CORVUS_ROUTE(d, "expm1/neg1", m_neg1);

  auto res = op::Neg(one);  // x <= -38
  if (!op::AllFalse(d, m_ser)) {
    CORVUS_ROUTE_SCOPE(d, m_ser);
    const auto xs = op::IfThenElse(m_ser, x, zero);
    res = op::IfThenElse(m_ser, DdToDouble(Expm1DdSeries(d, Dd<D>{xs, zero})),
                         res);
  }
  if (!op::AllFalse(d, m_exp)) {
    CORVUS_ROUTE_SCOPE(d, m_exp);
    // +inf is clamped by ExpDdFrac like any large x and overflows below.
    const auto xe = op::IfThenElse(m_exp, x, one);
    const auto parts = ExpDdFrac(d, xe, zero);
    const auto t = ScaleTwo(d, one, op::Neg(parts.e));  // 2^-e, exact
    const auto mm = DdAddD(d, parts.m, op::Neg(t));
    res = op::IfThenElse(m_exp, ScaleTwo(d, DdToDouble(mm), parts.e), res);
  }
  res = op::IfThenElse(op::Eq(x, zero), x, res);  // +-0 keeps its sign
  res = op::IfThenElse(op::Eq(x, inf), inf, res);
  CORVUS_ROUTE(d, "expm1/special/nan", op::IsNaN(x));
  return op::IfThenElse(op::IsNaN(x), x, res);  // payload-preserving
}

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#endif  // include guard
