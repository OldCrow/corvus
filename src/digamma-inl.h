// psi(x) = d/dx log Gamma(x), the digamma function, over the whole real axis.
// Per-target include guard (Highway -inl.h idiom).
//
// EVERYTHING BELOW IS ACCUMULATED IN DOUBLE-DOUBLE AND ROUNDED EXACTLY ONCE.
// Three of the four assemblies end in a subtraction of two same-signed
// quantities -- the (0,1) up-step, the asymptotic corrections and the whole
// negative axis -- and carrying ~106 bits is what makes those cancellations
// cost bits instead of digits.
//
// REGIONS (x > 0; the negative axis reflects onto these -- see DigammaVec)
//   (0, 1)     one up-step, psi(x) = psi(1 + x) - 1/x, with 1 + x NEVER
//              formed: the zone polynomial is evaluated at the exact shift
//              t1 = x - (x0 - 1) instead. fl(1 + x) would round, and near
//              x -> 1- that rounding is ~2^-52.5 relative of psi(1 + x).
//   [1, 2)     the zone: psi = t (*) P(t), t = x - x0 as a dd, x0 the unique
//              positive root. This is the full width-1 recurrence landing
//              interval -- under integer steps nothing narrower can be a
//              landing target -- so one polynomial serves every walked lane.
//   [2, X0)    masked fixed-step down-walk to [1, 2) by
//              psi(z) = psi(z - 1) + 1/(z - 1), at most kDigammaWalkDepth = 6
//              steps.
//   [X0, inf)  X0 = 8. Bernoulli asymptotic
//              psi = log x - 1/(2x) - x^-2 S(x^-2), K = 9 terms.
//
// WHY THE RECURRENCE GOES DOWN, NOT UP
//   x - k is EXACTLY representable for every integer k < x < 2^52 (x is a
//   multiple of ulp(x) <= 1 and |x - k| < |x|, so the difference needs no
//   finer grid than x already sits on -- Sterbenz). x + k is on the same grid
//   but LARGER, so it can need a bit the format does not have. Walking down
//   keeps every shifted argument exact and lands them all in [1, 2). The one
//   step that must go UP, from (0, 1), is the reason that branch carries its
//   own shifted centre rather than an argument.
//
// THE ROOT, AND WHY THE ZONE IS A PRODUCT
//   psi has exactly one positive zero, x0 ~ 1.4616. The zone evaluates
//   t (*) P(t) with t = x - x0 carried as a dd, so a zero t gives a zero
//   result whatever P is and relative accuracy survives arbitrarily close to
//   the root -- the same reason lgamma fits lgamma(c + t)/t rather than
//   lgamma itself. x0 is irrational, so t is formed by TwoSum against the dd
//   pair (kDigammaRootHi, kDigammaRootLo) and never by a bare subtraction.
//
// THE NEGATIVE AXIS
//   psi(x) = psi(1 - x) - pi*cot(pi x). Both halves are assembled in dd:
//   * y = 1 - x is EXACT (TwoSum). The positive pipeline dispatches on y.hi
//     and the residual is folded back as y.lo * psi'(y.hi); a ~2^-40
//     trigamma suffices for that correction, since y.lo is already 2^-53
//     relative. Threading a dd argument through all four regions instead
//     would double the width of every shift for no measurable gain.
//   * pi*cot(pi x) is built from the EXACT reduction u = x - round(x)
//     (|u| <= 1/2, exact for every double: |x| < 1/2 gives u = x, and
//     otherwise |round(x) - x| <= 1/2 <= |x|/2 puts Sterbenz in range).
//     pi cancels analytically: with sinc(u) = sin(pi u)/(pi u),
//         cos(pi u) / (u * sinc(u)) = pi*cos(pi u)/sin(pi u) = pi*cot(pi x),
//     the (-1)^n from x = n + u cancelling in the ratio. So no pi ever
//     multiplies anything here, and the 1/u pole -- which is exactly the
//     term that should diverge at a pole -- appears explicitly instead of as
//     an overflow of pi/sin(pi x).
//   ACCURACY DOCTRINE (lgamma analog):
//   psi has zeros at 20+ points on the negative axis with no closed form, so
//   near them a fixed absolute error is an unbounded relative one. The bound
//   is relative where |psi| >= 1 and absolute (2^-53 class) inside the zero
//   bands; see docs/ACCURACY.md for the measured split.
#if defined(CORVUS_DIGAMMA_INL_H_) == defined(HWY_TARGET_TOGGLE)
#ifdef CORVUS_DIGAMMA_INL_H_
#undef CORVUS_DIGAMMA_INL_H_
#else
#define CORVUS_DIGAMMA_INL_H_
#endif

#include <limits>

#include "src/dd-inl.h"
#include "src/digamma_data.h"
#include "src/log_dd-inl.h"
#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

namespace op = ops;

// Above this the asymptotic form returns log(x) ALONE.
//
// DERIVATION. The dropped terms are 1/(2x) + x^-2 S(x^-2), dominated by
// 1/(2x); against psi ~ log x the relative size is 1/(2x log x). The design's
// own criterion (1/(2x) < 2^-60 log x) is first met at x ~ 2^53.6. This cut is
// placed much higher, at 2^85, where the dropped part is
//     1/(2*2^85*log 2^85) = 2.19e-28 = 2^-92.0
// relative -- 2^-39 of a half ulp, so it cannot move the rounded double at
// all, not merely rarely.
//
// WHY A CUT EXISTS AT ALL. Two large-argument hazards live above it and both
// vanish with the corrections:
//   (a) x^-2 must never be formed as 1/(x*x): x*x overflows past ~1.3e154,
//       which is inside the domain. It is built as r = 1/x then w = r*r, and
//       the extra rounding is irrelevant -- w enters a term that is itself
//       only ~6e-4 of psi at X0 and falls as x^-2 from there.
//   (b) ops::ProdLow's non-FMA fallback is Dekker's split, whose intermediate
//       a*(2^27+1) overflows for |a| > 2^996 (AGENTS.md; lgamma's Stirling
//       product is the precedent). The only dd operand that scales with x is
//       2x inside DdRecip(2x). Clamping the corrections to x <= 2^85 keeps
//       every dd operand under 2^86, so there is no non-FMA special case and
//       no scale/unscale dance -- above the cut those ops simply do not run.
// LogDd itself is safe at any positive normal x: its slot index comes from
// the exponent bits, bounded by construction (AGENTS.md gather rule).
constexpr double kDigammaAsymCut = 0x1p85;

// Below this |argument| the reciprocal limb IS the answer: on (0, 1) psi(x) =
// -1/x + psi(1 + x) with |psi(1 + x)| <= 0.58, and near a negative pole
// psi(x) = -1/u + O(log|x|). At 2^-960 the neglected part is under
// 710 * 2^-960 = 2^-950 relative, far below a half ulp, so a plain
// correctly-rounded division reproduces psi exactly -- including the overflow
// to +-inf that is the documented answer at the poles and for tiny arguments.
//
// The threshold is not merely convenient: it is also what keeps DdRecip's and
// DdRecipDd's quotient under 2^996 and out of the Dekker ceiling of (b) above.
constexpr double kDigammaDirectRecip = 0x1p-960;

// P(t) = L0 + t*(L1 + ... + t*S(t)) with dd LEAD coefficients and a plain
// double tail S, for a t that is itself a dd. Shared by the zone (2 leads) and
// the reflection sinc/cos pair (3 leads each).
//
// THE SEAM IS NOT ROUNDED. S and the single product t*S are the only parts in
// plain double, but that product is formed as dd (DdMulD) rather than as one
// rounded double: t is carried at dd precision precisely because the caller
// built it with TwoSum/TwoProd, and collapsing the seam would throw away the
// precision that construction exists to protect. Compare LeadTailScalar, whose
// t is an ordinary double and whose seam is therefore one honest rounding.
// The generator replays this exact shape and refuses to emit a table that
// misses budget -- see tools/gen_digamma_data.py.
template <class D>
HWY_INLINE Dd<D> LeadTailDd(D d, Dd<D> t, const double* lead_hi,
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
HWY_INLINE Dd<D> LeadTailScalar(D d, op::V<D> t, const double* lead_hi,
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

// psi on [1, 2) as t (*) P(t), t = x - x0 in dd. Outlined like every other
// region core here: fully inlined, the single export becomes one enormous
// function per target and MSVC's optimizer is superlinear in function size
// (AGENTS.md). Contraction is off, so outlining cannot change FP semantics.
template <class D>
HWY_NOINLINE Dd<D> DigammaZoneDd(D d, Dd<D> t) {
  return DdMul(d, t,
               LeadTailDd(d, t, detail::kDigammaZoneLeadHi,
                          detail::kDigammaZoneLeadLo, detail::kDigammaZoneLead,
                          detail::kDigammaZoneCoef,
                          detail::kDigammaZoneNCoef));
}

// psi on (0, X0): the (0, 1) up-step, the zone, and the down-walk between
// them, as one pass.
//
// THE WALK. Step j fires when x >= j + 1, so the number of steps is
// floor(x) - 1 and every lane lands in [1, 2). x - j is EXACT (Sterbenz,
// x < X0 = 8) and each weight is DdRecip of that exact double, so the only
// error in the sum is dd assembly at ~2^-105. The accumulator is frozen by
// SELECT, never by adding a zero: DdAdd renormalizes, and a renormalization
// is value-preserving but not bit-preserving, which would make a lane's
// answer depend on how many steps its NEIGHBOURS needed. test_digamma_smoke's
// lane-mix check is what polices that.
//
// Cancellation is mild: psi(land) is in [-0.578, 0.424) and the sum is
// positive, so for x in [2, 3) -- the worst case -- the result 0.42 comes
// from 1.00 - 0.58, about 1.3 bits.
//
// DOMAIN CLAMPS. This core also runs on the lanes routed to the asymptotic
// branch, so x is clamped to X0 on entry and the reciprocal's argument is
// floored: a discarded lane at 1e308 would otherwise walk to a step of
// -1e308 and a (0,1) reciprocal of zero. Both clamps are no-ops on every live
// lane, so the mixed and all-in-region paths are bit-identical.
template <class D>
HWY_NOINLINE Dd<D> DigammaLowDd(D d, op::V<D> x_in) {
  const auto one = op::Set(d, 1.0);
  const auto x = op::Min(x_in, op::Set(d, detail::kDigammaX0));
  const auto lo1 = op::Lt(x, one);

  Dd<D> s{op::Zero(d), op::Zero(d)};
  auto y = x;
  for (int j = 1; j <= detail::kDigammaWalkDepth; ++j) {
    const auto fire = op::Ge(x, op::Set(d, static_cast<double>(j + 1)));
    // The fire masks shrink monotonically in j, so once no lane fires every
    // remaining step would select its own old value back -- breaking out
    // returns the identical result.
    if (op::AllFalse(d, fire)) break;
    // Exact by Sterbenz; the floor only touches lanes whose fire is false
    // (a firing lane has x - j >= 1 by the mask itself).
    const auto step =
        op::Max(op::Sub(x, op::Set(d, static_cast<double>(j))), one);
    const auto sn = DdAdd(d, s, DdRecip(d, step));
    s = Dd<D>{op::IfThenElse(fire, sn.hi, s.hi),
              op::IfThenElse(fire, sn.lo, s.lo)};
    y = op::IfThenElse(fire, step, y);
  }

  // The shifted argument, exact in both branches: y - x0 for [1, X0) and
  // x - (x0 - 1) for (0, 1). x0.hi - 1 is exact by Sterbenz (x0.hi is in
  // [1, 2)) and reuses kDigammaRootLo unchanged, which is the whole point of
  // shifting the CENTRE rather than the argument. TwoSum rather than a bare
  // subtraction because on (0, 1) the difference is not always exact.
  const auto centre = op::IfThenElse(lo1, op::Set(d, detail::kDigammaRootM1Hi),
                                     op::Set(d, detail::kDigammaRootHi));
  Dd<D> t = TwoSum(d, op::IfThenElse(lo1, x, y), op::Neg(centre));
  t = DdAddD(d, t, op::Set(d, -detail::kDigammaRootLo));

  Dd<D> res = DigammaZoneDd(d, t);

  const auto walked = op::Ge(x, op::Set(d, 2.0));
  const auto rw = DdAdd(d, res, s);
  res = Dd<D>{op::IfThenElse(walked, rw.hi, res.hi),
              op::IfThenElse(walked, rw.lo, res.lo)};

  // (0, 1): subtract 1/x. The pole is not a cancellation -- as x -> 0 the
  // -1/x term simply dominates -- and the probe-measured cancellation ratio
  // over the rest of the interval is under 1.74, about one bit. Lanes below
  // kDigammaDirectRecip are floored here and answered by the driver's direct
  // division instead.
  const auto rc = DdRecip(d, op::Max(x, op::Set(d, kDigammaDirectRecip)));
  const auto r01 = DdAdd(d, res, Dd<D>{op::Neg(rc.hi), op::Neg(rc.lo)});
  return Dd<D>{op::IfThenElse(lo1, r01.hi, res.hi),
               op::IfThenElse(lo1, r01.lo, res.lo)};
}

// psi for x >= X0: log x - 1/(2x) - x^-2 S(x^-2), K = 9 Bernoulli terms with
// one dd head. X0 = 8 is accuracy-forced: the probe table puts K = 9 at
// 2^-56.5 relative there, while X0 <= 5 cannot reach 2^-55 at any K.
//
// The two large-argument hazards and the cut that removes them are derived at
// kDigammaAsymCut. Note the grouping: log x is the whole answer to within
// 6e-4 at X0 and improves from there, so the corrections carry no
// cancellation of their own.
template <class D>
HWY_NOINLINE Dd<D> DigammaAsymDd(D d, op::V<D> x) {
  const auto one = op::Set(d, 1.0);
  const auto cut = op::Set(d, kDigammaAsymCut);
  const auto l = LogDd(d, x);

  // Corrections evaluated at the clamped argument; above the cut they are
  // selected away, and the clamp is what keeps every dd operand small.
  const auto xc = op::Min(x, cut);
  const auto r = op::Div(one, xc);  // never 1/(x*x): x*x overflows past 1.3e154
  const auto w = op::Mul(r, r);

  const auto term = DdMulD(
      d,
      LeadTailScalar(d, w, detail::kDigammaAsymHeadHi,
                     detail::kDigammaAsymHeadLo, detail::kDigammaAsymHead,
                     detail::kDigammaAsymCoef, detail::kDigammaAsymNCoef),
      w);
  const auto hr = DdRecip(d, op::Add(xc, xc));  // 1/(2x); the doubling is exact

  Dd<D> res = DdAdd(d, l, Dd<D>{op::Neg(hr.hi), op::Neg(hr.lo)});
  res = DdAdd(d, res, Dd<D>{op::Neg(term.hi), op::Neg(term.lo)});

  const auto big = op::Gt(x, cut);
  return Dd<D>{op::IfThenElse(big, l.hi, res.hi),
               op::IfThenElse(big, l.lo, res.lo)};
}

// psi for x > 0, as a dd. The caller rounds -- the reflection subtracts this
// from a term of comparable size, so it needs the full precision.
template <class D>
HWY_INLINE Dd<D> DigammaPosDd(D d, op::V<D> x) {
  const auto x0 = op::Set(d, detail::kDigammaX0);
  const auto big = op::Ge(x, x0);

  // Both branches are domain-clamped (DigammaLowDd clamps its own input), so
  // whichever one a lane is not taking cannot fault or manufacture a NaN that
  // a select would then have to launder.
  if (op::AllTrue(d, big)) return DigammaAsymDd(d, x);
  if (op::AllFalse(d, big)) return DigammaLowDd(d, x);
  const auto a = DigammaAsymDd(d, op::Max(x, x0));
  const auto b = DigammaLowDd(d, x);
  return Dd<D>{op::IfThenElse(big, a.hi, b.hi),
               op::IfThenElse(big, a.lo, b.lo)};
}

// trigamma to ~2^-40 relative, plain double throughout. ONLY used for the
// y.lo * psi'(y.hi) correction on the reflection's dd argument, where y.lo is
// already 2^-53 relative -- 2^-40 there lands at 2^-93 of the result.
//
// WALK FORM (mirrors src/digamma_data.h's own note, which pins it): step up
// to y >= kDigammaRoughTrigammaFloor = 6 accumulating 1/y^2, then the K = 8
// Bernoulli asymptotic. The floor is what buys the margin, not K: at y = 2
// the asymptotic measures ~2.4e-5 at its best, nowhere near 2^-40. Five steps
// suffice, since the argument is 1 - x > 1 for every live lane.
template <class D>
HWY_NOINLINE op::V<D> DigammaRoughTrigamma(D d, op::V<D> y_in) {
  const auto one = op::Set(d, 1.0);
  // Live lanes already satisfy y > 1; the floor keeps discarded ones (which
  // may be tiny, where 1/y^2 overflows) inside the walk's domain.
  auto y = op::Max(y_in, one);
  auto s = op::Zero(d);
  const auto floorv = op::Set(d, detail::kDigammaRoughTrigammaFloor);
  for (int k = 0; k < 5; ++k) {
    const auto fire = op::Lt(y, floorv);
    if (op::AllFalse(d, fire)) break;
    const auto inc = op::Div(one, op::Mul(y, y));
    s = op::IfThenElse(fire, op::Add(s, inc), s);
    y = op::IfThenElse(fire, op::Add(y, one), y);
  }

  const auto w = op::Div(one, op::Mul(y, y));  // 0 for a huge discarded lane
  auto p = op::Set(
      d, detail::kDigammaRoughTrigammaCoef[detail::kDigammaRoughTrigammaN - 1]);
  for (int k = detail::kDigammaRoughTrigammaN - 2; k >= 0; --k) {
    p = op::MulAdd(p, w, op::Set(d, detail::kDigammaRoughTrigammaCoef[k]));
  }
  const auto asym =
      op::Add(op::Add(op::Div(one, y), op::Mul(op::Set(d, 0.5), w)),
              op::Div(op::Mul(w, p), y));
  return op::Add(s, asym);
}

// pi*cot(pi x) as a dd, from the exact reduction u already computed by the
// caller. See the file header for why pi never appears: the ratio
// cos(pi u) / (u * sinc(u)) IS pi*cot(pi x), parity included.
//
// v = u^2 goes through TwoProd, i.e. the capability-guarded ops::ProdLow --
// never a bare MulSub, which is silently zero on non-FMA targets (AGENTS.md).
// The leading term of cos is 1 - (pi^2/2) v, large enough that v's own
// rounding would show. The cos fit is targeted at ABSOLUTE 2^-58 rather than
// relative: cos(pi u) has a zero at u = 1/2, where the cot assembly inherits
// the same zero-band doctrine as the result itself.
//
// The caller guarantees |u| >= kDigammaDirectRecip, which keeps DdRecipDd's
// quotient under 2^960 and so under the Dekker ceiling of ops::ProdLow's
// non-FMA path; below that threshold psi IS -1/u and the driver says so
// directly.
template <class D>
HWY_NOINLINE Dd<D> DigammaCotDd(D d, op::V<D> u) {
  const auto v = TwoProd(d, u, u);  // exact
  const auto sinc =
      LeadTailDd(d, v, detail::kDigammaSincLeadHi, detail::kDigammaSincLeadLo,
                 detail::kDigammaSincLead, detail::kDigammaSincCoef,
                 detail::kDigammaSincNCoef);
  const auto cs =
      LeadTailDd(d, v, detail::kDigammaCosLeadHi, detail::kDigammaCosLeadLo,
                 detail::kDigammaCosLead, detail::kDigammaCosCoef,
                 detail::kDigammaCosNCoef);
  // u * sinc(u) = sin(pi u)/pi, so the ratio below already carries the pi.
  return DdMul(d, cs, DdRecipDd(d, DdMulD(d, sinc, u)));
}

// HWY_NOINLINE like the region cores, and for the same MSVC-codegen reason:
// the driver is inlined TWICE per export (full-vector and masked-tail call
// sites), so leaving it inline doubles the largest function the optimizer
// sees.
template <class D>
HWY_NOINLINE op::V<D> DigammaVec(D d, op::V<D> x) {
  const auto zero = op::Zero(d);
  const auto one = op::Set(d, 1.0);
  const auto inf = op::Set(d, std::numeric_limits<double>::infinity());
  const auto qnan = op::Set(d, std::numeric_limits<double>::quiet_NaN());
  const auto maxf = op::Set(d, (std::numeric_limits<double>::max)());
  const auto cut = op::Set(d, kDigammaDirectRecip);

  const auto neg = op::Lt(x, zero);
  const auto isint = op::Eq(x, op::Round(x));  // true for +-inf as well

  // SCRUB. Masked-off lanes still execute every op (AGENTS.md), and although
  // this kernel gathers nothing, an unscrubbed lane would push infinities and
  // NaNs through the dd assemblies for no purpose. Every lane the pipeline
  // must not evaluate is replaced by a benign in-domain point of the RIGHT
  // SIGN -- +1.5 for the specials on the positive side, -1.5 for the poles,
  // so a scrubbed negative lane still takes a well-behaved reflection.
  auto xs = op::IfThenElse(neg, op::IfThenElse(isint, op::Set(d, -1.5), x), x);
  xs = op::IfThenElse(op::Eq(x, zero), op::Set(d, 1.5), xs);
  xs = op::IfThenElse(op::Gt(x, maxf), op::Set(d, 1.5), xs);  // +inf
  xs = op::IfThenElse(op::IsNaN(x), op::Set(d, 1.5), xs);

  // y = 1 - x, EXACT. Positive lanes compute it too and discard it.
  const auto ydd = TwoSum(d, one, op::Neg(xs));
  const auto arg = op::IfThenElse(neg, ydd.hi, xs);
  Dd<D> g = DigammaPosDd(d, arg);

  const auto u = op::Sub(xs, op::Round(xs));  // exact; |u| <= 1/2
  if (!op::AllFalse(d, neg)) {
    // psi(1 - x) at the exact dd argument: the pipeline ran on y.hi, so the
    // residual enters as y.lo * psi'(y.hi). Selected in rather than added
    // unconditionally -- y.lo is not zero on positive lanes (1 - x is inexact
    // for large x) and adding it there would be wrong, not merely wasteful.
    const auto gc = DdAddD(d, g, op::Mul(ydd.lo, DigammaRoughTrigamma(d, arg)));

    // Discarded lanes get a fixed benign u; live ones are floored in
    // magnitude at the point where the direct division below takes over.
    const auto ucot = op::IfThenElse(
        neg, op::CopySign(op::Max(op::Abs(u), cut), u), op::Set(d, 0.25));
    const auto cot = DigammaCotDd(d, ucot);
    const auto r = DdAdd(d, gc, Dd<D>{op::Neg(cot.hi), op::Neg(cot.lo)});
    g = Dd<D>{op::IfThenElse(neg, r.hi, g.hi), op::IfThenElse(neg, r.lo, g.lo)};
  }

  auto out = DdToDouble(g);

  // The reciprocal limb alone, where it is the whole answer to 2^-950
  // relative (see kDigammaDirectRecip). One correctly-rounded division, which
  // also delivers the documented +-inf when 1/x or 1/u overflows -- and note
  // x = -0 lands here rather than in the branch below, since -0 < 0 is false,
  // giving +inf as required.
  out = op::IfThenElse(
      neg,
      op::IfThenElse(op::Lt(op::Abs(u), cut), op::Neg(op::Div(one, u)), out),
      op::IfThenElse(op::Lt(x, cut), op::Neg(op::Div(one, x)), out));

  // Specials, outermost last. psi(+-0) = -+inf is the signed-zero pole
  // convention (scipy parity); every negative integer is a pole and gives
  // NaN, which also covers -inf and every double <= -2^53 (all integers).
  out = op::IfThenElse(op::Eq(x, zero), op::Neg(op::CopySign(inf, x)), out);
  out = op::IfThenElse(neg, op::IfThenElse(isint, qnan, out), out);
  out = op::IfThenElse(op::Gt(x, maxf), inf, out);  // psi(+inf) = +inf
  return op::IfThenElse(op::IsNaN(x), x, out);      // payload preserved
}

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#endif  // include guard
