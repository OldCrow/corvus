// psi_1(x) = d^2/dx^2 log Gamma(x), the trigamma function, over the whole real
// axis. Per-target include guard (Highway -inl.h idiom).
//
// EVERYTHING BELOW IS ACCUMULATED IN DOUBLE-DOUBLE AND ROUNDED EXACTLY ONCE.
// Only two assemblies subtract at all -- the [2, 8) down-walk and the negative
// axis -- and both are mild (1.3 and 0.15 bits respectively); the dd carriage
// is here mostly so that the *sum* of a fitted zone value and an exactly
// reciprocated pole term keeps its last bits, not to rescue a cancellation.
//
// WHY THE METRIC IS RELATIVE EVERYWHERE, WITH NO ABSOLUTE BAND
//   psi_1(x) = sum_{n>=0} 1/(x + n)^2 is a sum of squares, so it is strictly
//   POSITIVE wherever it is finite -- on the negative axis too, where the
//   reflection turns it into pi^2/sin^2(pi x) - psi_1(1 - x) and the per-
//   interval minima rise monotonically from 8.933 (the global minimum, at
//   x ~ -0.4957) to pi^2. psi_1 has no zero anywhere. That is the whole
//   difference in accuracy doctrine from digamma and lgamma, both of which
//   carry an absolute band around zeros they cannot reproduce exactly: here a
//   single relative bound holds on the entire real line, and the ULP test
//   gates one metric.
//
// REGIONS (x > 0; the negative axis reflects onto these -- see TrigammaVec)
//   (0, 1)     one up-step, psi_1(x) = psi_1(1 + x) + 1/x^2, with 1 + x NEVER
//              formed: the zone polynomial is evaluated at the exact shift
//              t1 = x - (c - 1) instead, c = 1.5 being the fit's own centre.
//              Both terms are positive, so this branch has no cancellation at
//              all; what it does have is a pole, and 1/x^2 is carried in dd
//              end to end because a naive (1/x)^2 or 1/(x*x) is not reliably
//              correctly rounded (24-46% 1-ULP misses, measured by the probe
//              -- see src/trigamma_data.h's warning at the deep-tiny guard).
//   [1, 2)     the zone: psi_1 = P(t), t = x - c. A PLAIN VALUE fit, not
//              digamma's product form -- there is no root to divide out. c is
//              exactly representable, so t needs no dd centre; TwoSum is still
//              used because on (0, 1) the shift by c - 1 is not always exact.
//   [2, X0)    masked fixed-step down-walk to [1, 2) by
//              psi_1(z + 1) = psi_1(z) - 1/z^2, at most
//              kTrigammaWalkDepth = 6 steps.
//   [X0, cut)  X0 = 8. Bernoulli asymptotic in the DIRECT (unfactored) form
//              psi_1 = 1/x + 1/(2x^2) + x^-3 S(x^-2), K = 11 terms. No
//              logarithm appears anywhere: unlike digamma, this family does
//              not consume log_dd at all.
//   [cut, inf) cut = 2^89. fl(1/x) alone; see kTrigammaAsymCut below.
//
// WHY THE RECURRENCE GOES DOWN, NOT UP
//   x - k is EXACTLY representable for every integer k < x < 2^52 (Sterbenz),
//   while x + k need not be, so walking down keeps every shifted argument
//   exact and lands them all in [1, 2). The one step that must go UP, from
//   (0, 1), is the reason that branch carries its own shifted centre rather
//   than a shifted argument.
//
// THE NEGATIVE AXIS
//   psi_1(x) = pi^2/sin^2(pi x) - psi_1(1 - x). Both halves are dd:
//   * y = 1 - x is EXACT (TwoSum). The positive pipeline dispatches on y.hi
//     and the residual folds back as y.lo * psi_1'(y.hi) = y.lo * psi_2(y.hi);
//     a ~2^-30 tetragamma suffices there, because y.lo is already 2^-53
//     relative and psi_1 >= 8.93 on this axis, which puts the whole
//     correction's own error near 2^-56 relative of the result.
//   * pi^2/sin^2(pi x) is built from the EXACT reduction u = x - round(x)
//     (|u| <= 1/2, exact for every double). With sinc(u) = sin(pi u)/(pi u),
//         u * sinc(u) = sin(pi u)/pi   exactly by construction,
//     so 1/(u sinc u)^2 = pi^2/sin^2(pi u) = pi^2/sin^2(pi x). NO cos table
//     and no parity handling: sin^2 has period pi, so the (-1)^n from
//     x = n + u disappears under the square. This is where trigamma is
//     strictly simpler than digamma's cot ratio, not merely different.
//   The 1/u^2 pole appears explicitly rather than as an overflow of
//   pi^2/sin^2, which is what makes the pole values fall out of the
//   arithmetic instead of needing to be manufactured.
//
// POLE DOCTRINE (scipy parity, and deliberately NOT digamma's)
//   Every pole of psi_1 is a DOUBLE pole, so its sign is unambiguous: the
//   answer is +inf at +0 AND at -0, at every negative integer (which includes
//   -inf and every double <= -2^53, all of which are integers), and wherever
//   the reciprocal-square overflows. psi_1(+inf) = +0. digamma returns NaN at
//   its negative-integer poles because a simple pole has no sign; that
//   reasoning does not apply here.
#if defined(CORVUS_TRIGAMMA_INL_H_) == defined(HWY_TARGET_TOGGLE)
#ifdef CORVUS_TRIGAMMA_INL_H_
#undef CORVUS_TRIGAMMA_INL_H_
#else
#define CORVUS_TRIGAMMA_INL_H_
#endif

#include <limits>

#include "src/dd-inl.h"
#include "src/ops-inl.h"
#include "src/trigamma_data.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

namespace op = ops;

// Scale used by the deep-tiny reciprocal-square, chosen so that the SQUARE of
// the scaled reciprocal cannot overflow: the branch runs only below
// kTrigammaDeepTinyGuard = 2^-480, where a*2^512 lands in (0, 2^32) and the
// clamp below pins it at 1 for everything at or under 2^-512. Exact, being a
// power of two.
constexpr double kTrigammaTinyScale = 0x1p512;

// P(t) = L0 + t*(L1 + ... + t*S(t)) with dd LEAD coefficients and a plain
// double tail S, for a t that is itself a dd. Shared by the zone (3 leads) and
// the reflection sinc (3 leads).
//
// THE SEAM IS NOT ROUNDED. S and the single product t*S are the only parts in
// plain double, but that product is formed as dd (DdMulD) rather than as one
// rounded double: t is carried at dd precision precisely because the caller
// built it with TwoSum/TwoProd, and collapsing the seam would throw away the
// precision that construction exists to protect. tools/gen_trigamma_data.py
// replays this exact shape (eval_lead_tail_dd) and refuses to emit a table
// that misses budget.
//
// This and TrigammaLeadTailScalar below duplicate digamma's two evaluators
// rather than sharing them. That is the same call PLAN.md's design already
// made for the reflection constants: hoisting a helper out of a SHIPPED family
// invokes the byte-identity protocol on digamma's whole gate set, and ~30
// duplicated lines are cheaper than a revalidation pass. If a third consumer
// appears, hoist all of it at once.
template <class D>
HWY_INLINE Dd<D> TrigammaLeadTailDd(D d, Dd<D> t, const double* lead_hi,
                                    const double* lead_lo, int n_lead,
                                    const double* coef, int n_coef) {
  auto s = op::Set(d, coef[n_coef - 1]);
  for (int k = n_coef - 2; k >= 0; --k) {
    s = op::MulAdd(s, t.hi, op::Set(d, coef[k]));
  }
  Dd<D> acc = DdMulD(d, t, s);
  acc = DdAdd(d, acc,
              Dd<D>{op::Set(d, lead_hi[n_lead - 1]),
                    op::Set(d, lead_lo[n_lead - 1])});
  for (int k = n_lead - 2; k >= 0; --k) {
    acc = DdAdd(d, DdMul(d, acc, t),
                Dd<D>{op::Set(d, lead_hi[k]), op::Set(d, lead_lo[k])});
  }
  return acc;
}

// Same shape for a plain-double argument (the asymptotic's w = x^-2): the
// lead/tail seam IS one rounding, because w carries no more than working
// precision to begin with.
template <class D>
HWY_INLINE Dd<D> TrigammaLeadTailScalar(D d, op::V<D> t, const double* lead_hi,
                                        const double* lead_lo, int n_lead,
                                        const double* coef, int n_coef) {
  auto s = op::Set(d, coef[n_coef - 1]);
  for (int k = n_coef - 2; k >= 0; --k) {
    s = op::MulAdd(s, t, op::Set(d, coef[k]));
  }
  Dd<D> acc{op::Set(d, lead_hi[n_lead - 1]), op::Set(d, lead_lo[n_lead - 1])};
  acc = DdAddD(d, acc, op::Mul(s, t));
  for (int k = n_lead - 2; k >= 0; --k) {
    acc = DdAdd(d, DdMulD(d, acc, t),
                Dd<D>{op::Set(d, lead_hi[k]), op::Set(d, lead_lo[k])});
  }
  return acc;
}

// 1/a^2 as a dd, ~2^-105 relative. Every reciprocal-square in this kernel goes
// through here, and all of them take the reciprocal FIRST.
//
// WHY NOT 1/TwoProd(a, a). Squaring first is one dd operation cheaper, but it
// puts a^2 -- not a -- through Dekker's split on non-FMA targets, and a^2 is
// the quantity that leaves the safe range at BOTH ends: it overflows the
// split's 2^996 ceiling for large a, and for a near the deep-tiny guard
// (a ~ 2^-480) the split's low limb squares to ~2^-1014, one binade from
// subnormal, where the residual quietly loses bits. Reciprocating first keeps
// every operand of every ProdLow inside [2^-480, 2^512] for every caller here.
template <class D>
HWY_INLINE Dd<D> TrigammaRecipSqDd(D d, op::V<D> a) {
  const auto r = DdRecip(d, a);
  return DdMul(d, r, r);
}

// psi_1 on [1, 2) as P(t), t = x - c in dd. A VALUE fit: unlike digamma's
// zone there is no factored-out root, because psi_1 has no zero to reproduce.
// Outlined like every other region core here -- fully inlined, the single
// export becomes one enormous function per target and MSVC's optimizer is
// superlinear in function size (AGENTS.md; the 2026-07-29 CI timeouts).
// Contraction is off, so outlining cannot change FP semantics.
template <class D>
HWY_NOINLINE Dd<D> TrigammaZoneDd(D d, Dd<D> t) {
  return TrigammaLeadTailDd(d, t, detail::kTrigammaZoneLeadHi,
                            detail::kTrigammaZoneLeadLo,
                            detail::kTrigammaZoneLead,
                            detail::kTrigammaZoneCoef,
                            detail::kTrigammaZoneNCoef);
}

// psi_1 on (0, X0): the (0, 1) up-step, the zone, and the down-walk between
// them, as one pass.
//
// THE WALK. Step j fires when x >= j + 1, so the number of steps is
// floor(x) - 1 and every lane lands in [1, 2). x - j is EXACT (Sterbenz,
// x < X0 = 8) and each weight is the dd reciprocal-square of that exact
// double, so the only error in the sum is dd assembly at ~2^-105. The
// accumulator is frozen by SELECT, never by adding a zero: DdAdd
// renormalizes, and a renormalization is value-preserving but not
// bit-preserving, which would make a lane's answer depend on how many steps
// its NEIGHBOURS needed. test_trigamma_smoke's lane-mix check polices that.
//
// CANCELLATION, and why the walk is the one bucket allowed 2 ULP. The
// subtracted sum is exact to dd, so the walk's ABSOLUTE error is the zone
// fit's own absolute error at the landing point -- but the output shrinks
// monotonically from psi_1(2) = 0.645 to psi_1(8) = 0.133 as the walk
// deepens, while the landing value can be as large as psi_1(1) = 1.645. The
// same absolute error against a smaller output is a larger relative one, by up
// to psi_1(1)/psi_1(8) = 12.36x -- the SECOND CORRECTION recorded in
// tools/gen_trigamma_data.py, whose predicted worst point (x just above 7,
// landing nearest 1) is where the generator's replay actually finds it.
//
// DOMAIN CLAMPS. This core also runs on the lanes routed to the asymptotic
// branch, so x is clamped to X0 on entry, and clamped from BELOW at the
// deep-tiny guard so that a discarded subnormal lane cannot overflow 1/x^2
// (those lanes are answered by the driver's scaled reciprocal-square instead).
// The walk's step is floored at 1 for the same reason. All three clamps are
// no-ops on every live lane, so the mixed and all-in-region paths are
// bit-identical.
template <class D>
HWY_NOINLINE Dd<D> TrigammaLowDd(D d, op::V<D> x_in) {
  const auto one = op::Set(d, 1.0);
  const auto x =
      op::Min(op::Max(x_in, op::Set(d, detail::kTrigammaDeepTinyGuard)),
              op::Set(d, detail::kTrigammaX0));
  const auto lo1 = op::Lt(x, one);

  Dd<D> s{op::Zero(d), op::Zero(d)};
  auto y = x;
  for (int j = 1; j <= detail::kTrigammaWalkDepth; ++j) {
    const auto fire = op::Ge(x, op::Set(d, static_cast<double>(j + 1)));
    // The fire masks shrink monotonically in j, so once no lane fires every
    // remaining step would select its own old value back -- breaking out
    // returns the identical result.
    if (op::AllFalse(d, fire)) break;
    const auto step =
        op::Max(op::Sub(x, op::Set(d, static_cast<double>(j))), one);
    const auto sn = DdAdd(d, s, TrigammaRecipSqDd(d, step));
    s = Dd<D>{op::IfThenElse(fire, sn.hi, s.hi),
              op::IfThenElse(fire, sn.lo, s.lo)};
    y = op::IfThenElse(fire, step, y);
  }

  // The shifted argument. c = 1.5 and c - 1 = 0.5 are both exact doubles, so
  // unlike digamma's irrational root there is no dd centre to correct for.
  // TwoSum rather than a bare subtraction all the same: on [1, 2) the shift is
  // exact by Sterbenz, but on (0, 1) x - 0.5 rounds once x drops below ~2^-53
  // of the centre, and that residual is the whole reason the branch shifts the
  // CENTRE instead of forming 1 + x.
  const auto centre =
      op::IfThenElse(lo1, op::Set(d, detail::kTrigammaZoneCentreM1),
                     op::Set(d, detail::kTrigammaZoneCentre));
  const Dd<D> t = TwoSum(d, op::IfThenElse(lo1, x, y), op::Neg(centre));

  Dd<D> res = TrigammaZoneDd(d, t);

  const auto walked = op::Ge(x, op::Set(d, 2.0));
  const auto rw = DdAdd(d, res, Dd<D>{op::Neg(s.hi), op::Neg(s.lo)});
  res = Dd<D>{op::IfThenElse(walked, rw.hi, res.hi),
              op::IfThenElse(walked, rw.lo, res.lo)};

  // (0, 1): ADD 1/x^2. Both terms are positive -- the up-step is the one
  // recurrence in this kernel with no cancellation whatsoever -- and the pole
  // is carried by the reciprocal-square rather than by any fit.
  const auto r01 = DdAdd(d, res, TrigammaRecipSqDd(d, x));
  return Dd<D>{op::IfThenElse(lo1, r01.hi, res.hi),
               op::IfThenElse(lo1, r01.lo, res.lo)};
}

// psi_1 for x >= X0: 1/x + 1/(2x^2) + x^-3 S(x^-2), K = 11 Bernoulli terms
// with one dd head, in the DIRECT (unfactored) form -- the probe measured it
// better conditioned than the 1/x-factored alternative, and it is also the
// form the generator replays.
//
// THE CUT, kTrigammaAsymCut = 2^89. Above it the answer is fl(1/x) alone. The
// dropped part is 1/(2x^2) + x^-3 S, dominated by 1/(2x^2), whose size
// relative to 1/x is 1/(2x) = 2^-90 at the cut and falls from there -- 2^-37
// of a half ulp, so it cannot move the rounded double at all, not merely
// rarely. The cut is not only an economy: it retires every large-operand dd
// op in this kernel. Below it, DdRecip's quotient stays under 2^-89 and
// ops::ProdLow's non-FMA Dekker split (which overflows for operands past
// 2^996, AGENTS.md) is never approached; above it there is nothing but one
// correctly-rounded division.
//
// Note also w = (1/x)^2 and NOT 1/(x*x): x*x overflows past ~1.3e154, which
// would be inside the domain if the corrections were evaluated there. They are
// not -- the argument is clamped to the cut first -- but the clamp is a
// routing decision and the squaring order is a correctness one, so both are
// enforced independently.
template <class D>
HWY_NOINLINE Dd<D> TrigammaAsymDd(D d, op::V<D> x) {
  const auto one = op::Set(d, 1.0);
  const auto cut = op::Set(d, detail::kTrigammaAsymCut);

  // Corrections evaluated at the clamped argument; at and above the cut they
  // are selected away, and the clamp is what keeps every dd operand small.
  const auto xc = op::Min(x, cut);
  const auto q = op::Div(one, xc);
  const auto w = op::Mul(q, q);
  const auto r = DdRecip(d, xc);  // dd 1/x, the leading term itself

  // 1/(2x^2) from the dd reciprocal, not from a rounded 2*x*x: that term is
  // ~6% of the result at X0, so a half ulp of its own would land at 2^-57 of
  // the answer -- inside the budget rather than under it.
  const auto half = DdMulD(d, DdMul(d, r, r), op::Set(d, 0.5));

  const auto s = TrigammaLeadTailScalar(
      d, w, detail::kTrigammaAsymHeadHi, detail::kTrigammaAsymHeadLo,
      detail::kTrigammaAsymHead, detail::kTrigammaAsymCoef,
      detail::kTrigammaAsymNCoef);
  const auto term = DdMul(d, DdMulD(d, s, w), r);  // x^-3 S(x^-2)

  // Every term is positive and strictly decreasing, so summing smallest-first
  // costs nothing and there is no cancellation anywhere on this branch.
  const Dd<D> res = DdAdd(d, r, DdAdd(d, half, term));

  const auto big = op::Ge(x, cut);
  return Dd<D>{op::IfThenElse(big, op::Div(one, x), res.hi),
               op::IfThenElse(big, op::Zero(d), res.lo)};
}

// psi_1 for x > 0, as a dd. The caller rounds -- the reflection subtracts this
// from a term of comparable size, so it needs the full precision.
template <class D>
HWY_INLINE Dd<D> TrigammaPosDd(D d, op::V<D> x) {
  const auto x0 = op::Set(d, detail::kTrigammaX0);
  const auto big = op::Ge(x, x0);

  // Both branches are domain-clamped (TrigammaLowDd clamps its own input), so
  // whichever one a lane is not taking cannot fault or manufacture a NaN that
  // a select would then have to launder.
  if (op::AllTrue(d, big)) return TrigammaAsymDd(d, x);
  if (op::AllFalse(d, big)) return TrigammaLowDd(d, x);
  const auto a = TrigammaAsymDd(d, op::Max(x, x0));
  const auto b = TrigammaLowDd(d, x);
  return Dd<D>{op::IfThenElse(big, a.hi, b.hi),
               op::IfThenElse(big, a.lo, b.lo)};
}

// tetragamma to ~2^-30 relative, plain double throughout. ONLY used for the
// y.lo * psi_2(y.hi) correction on the reflection's dd argument, where y.lo is
// already 2^-53 relative and psi_1 >= 8.93, which puts the whole correction's
// error near 2^-56 of the result.
//
// WALK FORM: mirrors src/trigamma_data.h's own note, which pins it. Step up to
// y >= kTrigammaRoughTetraFloor accumulating -2/y^3 (the recurrence is
// psi_2(y) = psi_2(y + 1) - 2/y^3), then the K = 6 Bernoulli asymptotic
//     psi_2(y) ~= -(w + w/y + w^2 Q(w)),  w = 1/y^2.
// Five steps suffice: the argument is 1 - x > 1 for every live lane.
template <class D>
HWY_NOINLINE op::V<D> TrigammaRoughTetragamma(D d, op::V<D> y_in) {
  const auto one = op::Set(d, 1.0);
  const auto two = op::Set(d, 2.0);
  // Live lanes already satisfy y > 1; the floor keeps discarded ones (which
  // may be tiny, where 1/y^3 overflows) inside the walk's domain.
  auto y = op::Max(y_in, one);
  auto s = op::Zero(d);
  const auto floorv = op::Set(d, detail::kTrigammaRoughTetraFloor);
  for (int k = 0; k < 5; ++k) {
    const auto fire = op::Lt(y, floorv);
    if (op::AllFalse(d, fire)) break;
    const auto dec = op::Div(two, op::Mul(op::Mul(y, y), y));
    s = op::IfThenElse(fire, op::Sub(s, dec), s);
    y = op::IfThenElse(fire, op::Add(y, one), y);
  }

  const auto w = op::Div(one, op::Mul(y, y));  // 0 for a huge discarded lane
  auto p = op::Set(
      d, detail::kTrigammaRoughTetraCoef[detail::kTrigammaRoughTetraN - 1]);
  for (int k = detail::kTrigammaRoughTetraN - 2; k >= 0; --k) {
    p = op::MulAdd(p, w, op::Set(d, detail::kTrigammaRoughTetraCoef[k]));
  }
  const auto asym = op::Neg(op::Add(op::Add(w, op::Div(w, y)),
                                    op::Mul(op::Mul(w, w), p)));
  return op::Add(s, asym);
}

// pi^2/sin^2(pi x) as a dd, from the exact reduction u already computed by the
// caller. See the file header for why neither pi nor a cos table appears:
// u * sinc(u) IS sin(pi u)/pi, and the square kills the (-1)^n parity.
//
// v = u^2 goes through TwoProd, i.e. the capability-guarded ops::ProdLow --
// never a bare MulSub, which is silently zero on non-FMA targets (AGENTS.md).
//
// RECIPROCATE FIRST, THEN SQUARE. Both orders carry the same ~2^-104 relative
// error (two dd operations either way), so the choice is made on RANGE. The
// caller floors |u| at kTrigammaDeepTinyGuard = 2^-480, so denom >= ~2^-480
// and this order puts operands of at most 2^480 through ops::ProdLow's Dekker
// split. Squaring first would instead hand DdRecipDd a denominator near
// 2^-960 and a quotient near 2^960, whose split intermediate reaches 2^987 --
// nine binades from the 2^996 ceiling, which is not a margin worth taking for
// no accuracy gain.
template <class D>
HWY_NOINLINE Dd<D> TrigammaReflDd(D d, op::V<D> u) {
  const auto v = TwoProd(d, u, u);  // exact
  const auto sinc = TrigammaLeadTailDd(
      d, v, detail::kTrigammaSincLeadHi, detail::kTrigammaSincLeadLo,
      detail::kTrigammaSincLead, detail::kTrigammaSincCoef,
      detail::kTrigammaSincNCoef);
  const auto denom = DdMulD(d, sinc, u);  // = sin(pi u)/pi, either sign
  const auto rec = DdRecipDd(d, denom);
  return DdMul(d, rec, rec);
}

// Correctly-rounded 1/a^2 for a tiny positive a, including the overflow to
// +inf, as ONE rounding.
//
// This is the deep-tiny lane of both poles -- x -> 0+ on (0, 1) and u -> 0 at
// a negative-integer pole -- where src/trigamma_data.h's guard says the zone
// term (at most pi^2/6) is below 2^-950 relative of 1/a^2 and can be dropped
// outright. It cannot be computed in place, because the ANSWER overflows below
// a = 2^-512 = 1/sqrt(DBL_MAX) while dd arithmetic cannot deliver an infinity:
// TwoProd of an overflowing product returns inf - inf = NaN in its residual,
// and Fast2Sum launders that into a NaN result. So the whole computation is
// done scaled, exactly.
//
//   as = clamp(a * 2^512, 1, 2^60)          exact (power of two)
//   V  = dd 1/as^2                          in (2^-120, 1], no overflow
//   out = fl(V * 2^512) * 2^512
//
// The first scaling is exact and lands fl(V * 2^512) in [2^392, 2^512], so the
// dd's single rounding happens there; the second is exact too, or overflows to
// +inf. Rounding commutes with exact power-of-two scaling, so the result is
// fl(V * 2^1024) = fl(1/a^2) -- one rounding, correctly rounded, with the
// infinity arriving on its own.
//
// The lower clamp is not hygiene, it is the pole: every a <= 2^-512 pins
// as = 1, V = 1, and 2^512 * 2^512 = +inf, which is exactly the answer there.
// The upper clamp only exists because this runs on lanes the driver will throw
// away, where a * 2^512 may itself be inf.
template <class D>
HWY_NOINLINE op::V<D> TrigammaTinyRecipSq(D d, op::V<D> a) {
  const auto s = op::Set(d, kTrigammaTinyScale);
  const auto as = op::Min(op::Max(op::Mul(a, s), op::Set(d, 1.0)),
                          op::Set(d, 0x1p60));
  const auto q = TrigammaRecipSqDd(d, as);
  return op::Mul(op::Add(op::Mul(q.hi, s), op::Mul(q.lo, s)), s);
}

// HWY_NOINLINE like the region cores, and for the same MSVC-codegen reason:
// the driver is inlined TWICE per export (full-vector and masked-tail call
// sites), so leaving it inline doubles the largest function the optimizer
// sees.
template <class D>
HWY_NOINLINE op::V<D> TrigammaVec(D d, op::V<D> x) {
  const auto zero = op::Zero(d);
  const auto one = op::Set(d, 1.0);
  const auto inf = op::Set(d, std::numeric_limits<double>::infinity());
  const auto maxf = op::Set(d, (std::numeric_limits<double>::max)());
  const auto guard = op::Set(d, detail::kTrigammaDeepTinyGuard);

  const auto neg = op::Lt(x, zero);           // false for -0.0, as required
  const auto isint = op::Eq(x, op::Round(x));  // true for +-inf as well

  // SCRUB. Masked-off lanes still execute every op (AGENTS.md), and although
  // this kernel gathers nothing, an unscrubbed lane would push infinities and
  // NaNs through the dd assemblies for no purpose -- and -inf in particular
  // would make u = -inf - (-inf) = NaN. Every lane the pipeline must not
  // evaluate is replaced by a benign in-domain point of the RIGHT SIGN: +1.5
  // for the specials on the positive side, -1.5 for the poles, so a scrubbed
  // negative lane still takes a well-behaved reflection.
  auto xs = op::IfThenElse(neg, op::IfThenElse(isint, op::Set(d, -1.5), x), x);
  xs = op::IfThenElse(op::Eq(x, zero), op::Set(d, 1.5), xs);
  xs = op::IfThenElse(op::Gt(x, maxf), op::Set(d, 1.5), xs);  // +inf
  xs = op::IfThenElse(op::IsNaN(x), op::Set(d, 1.5), xs);

  // y = 1 - x, EXACT. Positive lanes compute it too and discard it.
  const auto ydd = TwoSum(d, one, op::Neg(xs));
  const auto arg = op::IfThenElse(neg, ydd.hi, xs);
  Dd<D> g = TrigammaPosDd(d, arg);

  const auto u = op::Sub(xs, op::Round(xs));  // exact; |u| <= 1/2
  if (!op::AllFalse(d, neg)) {
    // psi_1(1 - x) at the exact dd argument: the pipeline ran on y.hi, so the
    // residual enters as y.lo * psi_2(y.hi). Selected in rather than added
    // unconditionally -- y.lo is not zero on positive lanes (1 - x is inexact
    // for large x) and adding it there would be wrong, not merely wasteful.
    const auto gc =
        DdAddD(d, g, op::Mul(ydd.lo, TrigammaRoughTetragamma(d, arg)));

    // Discarded lanes get a fixed benign u; live ones are floored in
    // magnitude at the point where the scaled reciprocal-square takes over.
    const auto ur = op::IfThenElse(
        neg, op::CopySign(op::Max(op::Abs(u), guard), u), op::Set(d, 0.25));
    const auto rf = TrigammaReflDd(d, ur);
    const auto r = DdAdd(d, rf, Dd<D>{op::Neg(gc.hi), op::Neg(gc.lo)});
    g = Dd<D>{op::IfThenElse(neg, r.hi, g.hi), op::IfThenElse(neg, r.lo, g.lo)};
  }

  auto out = DdToDouble(g);

  // The deep-tiny pole lane, where 1/a^2 alone IS the answer to 2^-950
  // relative: a = x approaching 0 on the positive side, a = |u| approaching a
  // negative-integer pole on the other. One helper serves both because the
  // reflection's sinc factor differs from 1 by less than 2^-959 there.
  const auto ta = op::IfThenElse(neg, op::Abs(u), x);
  out = op::IfThenElse(op::Lt(ta, guard), TrigammaTinyRecipSq(d, ta), out);

  // Specials, outermost last. Every pole is a DOUBLE pole and so unsigned:
  // +inf at +-0, at every negative integer (which covers -inf and every
  // double <= -2^53), and psi_1(+inf) = +0.
  out = op::IfThenElse(op::Eq(x, zero), inf, out);
  out = op::IfThenElse(neg, op::IfThenElse(isint, inf, out), out);
  out = op::IfThenElse(op::Gt(x, maxf), zero, out);
  return op::IfThenElse(op::IsNaN(x), x, out);  // payload preserved
}

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#endif  // include guard
