// Regularized incomplete gamma P(a,x) and Q(a,x) = 1 - P(a,x).
// Per-target include guard (Highway -inl.h idiom).
//
// Both public functions share one kernel: four region cores plus a router
// that differs between them in exactly two places -- which core owns the
// R1/R4 overlap, and which side of the pair gets complemented. Everything
// here follows PLAN.md, "Phase C part 2 -- regularized incomplete gamma P/Q
// detail design", whose region map and fixed lengths are probe-validated;
// the tables and every threshold live in src/gamma_data.h.
//
// THE ONE RULE THE WHOLE DESIGN TURNS ON: always compute the SMALLER of the
// pair directly, and get the other one as 1 (-) it, as a double-double
// subtraction rounded once. A directly computed side keeps full relative
// accuracy however tiny it gets; the complement is by construction >= ~0.4,
// where 1 (-) small costs nothing. Computing the large side directly and
// subtracting would be the opposite trade and is never done.
//
// REGIONS (lambda = x/a, a_T = kGammaAT = 20)
//   R1  series, P-direct   {a < 20, x <= a+1} u {a >= 20, lambda <= 1/2}
//       P = e^E' * sum_{n>=0} x^n / ((a+1)...(a+n)),  E' = E - log a.
//       Fixed cap kGammaSeriesN with a per-lane freeze mask; terms are
//       plain double (each is a factor of the next, so their rounding is
//       attenuated by the tail they head), the sum is dd.
//   R2  backward CF, Q-direct  {a < 20, x > a+1} u {a >= 20, lambda >= 2}
//       Legendre's continued fraction run BACKWARD from a fixed depth
//       kGammaCfN with NO convergence test. Backward is not a stylistic
//       choice: forward Lentz measured 24 ULP against backward's 4, and
//       worse, it FALSE-CONVERGES -- (a,x) = (30,133) satisfied a 1e-16
//       relative stopping test at N = 8 with 70 ULP of error. A fixed depth
//       has no stopping test to fool and contracts rounding as it goes.
//   R3  Temme uniform asymptotic  {a >= 20, 1/2 < lambda < 2}
//       The ridge, where both P and Q are O(1) and neither series nor CF
//       converges usefully. See GammaTemme.
//   R4  small-a Q-direct   {0 < a <= 3/2, 0 < x <= 4}
//       The corner that closes small-a complement accuracy: there Q ~ a*E1(x)
//       is genuinely tiny while P -> 1, so R2's CF (which owns this corner
//       by the R2 rule) is fine for Q but R1's P is what gamma_p wants.
//       See GammaSmallQ.
//
// ROUTING DIFFERENCE. R4 overlaps R1 (for x <= a+1) and R2 (for x > a+1).
// gamma_q takes R4 for the whole {a <= 3/2, x <= 4} box, because R4 computes
// Q directly there. gamma_p keeps R1 wherever R1 applies: R4's Q is O(1)
// there and 1 (-) Q would destroy the relative accuracy of a tiny P.
//
// PREFACTOR. E = ln(x^a e^{-x} / Gamma(a)) = a*log(x) (-) x (-) lgamma(a),
// all in dd, with lgamma from LgammaPosDd (src/lgamma-inl.h), which covers
// every positive a including its own 0 < a < 1/2 log-shift band -- so the
// argument fed to it is always exact and no separate large-a Stirling form
// is needed. R1 additionally folds (-) log a, which is what turns Gamma(a)
// into Gamma(a+1) WITHOUT ever forming the inexact argument 1 + a and
// without a 1/a that would overflow for subnormal a. R3 needs no lgamma at
// all: the extracted Temme coefficients absorb the Stirling remainder and
// the 1/sqrt(2*pi*a) factor is formed directly.
//
// SATURATION. e^E underflows to zero for E below about -745, and every
// region that ends in an exponential clamps its argument at
// kGammaExpFloor = -800 and reports the clamp (GammaClampE). Two things
// hang on that mask: the result is forced to an exact 0/1 pair, and the
// (a,x) fed to the series and the CF are scrubbed to a benign (1,3) --
// without which the CF's j*(j-a) overflows for a above ~4e306, which is
// squarely inside the domain (a = 1e307, x = 2e307 is an ordinary R2
// point whose Q is zero). The same clamp catches the NaN that a*log(x)
// produces by inf - inf once a exceeds ~2.5e305.
#if defined(CORVUS_GAMMA_INL_H_) == defined(HWY_TARGET_TOGGLE)
#ifdef CORVUS_GAMMA_INL_H_
#undef CORVUS_GAMMA_INL_H_
#else
#define CORVUS_GAMMA_INL_H_
#endif

#include <limits>

#include "src/dd-inl.h"
#include "src/erfc_core-inl.h"
#include "src/exp_dd-inl.h"
#include "src/gamma_data.h"
#include "src/lgamma-inl.h"
#include "src/log_dd-inl.h"
#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

namespace op = ops;

// --- kernel-internal constants ------------------------------------------
// These are structural (crossover points and loop shapes chosen from the
// error budget in the comments at their use sites), not fitted table data,
// so they live here rather than in the generated src/gamma_data.h -- the
// same split erfinv makes for kErfinvTwoOverSqrtPi.

// Expm1Dd series/exponential crossover. Below it the k <= 6 series truncates
// at w^7/5040 <= 2^-72 relative; above it e^w (-) 1 loses at most
// log2(1/|w|) = 10 bits of the dd's ~2^-105, i.e. lands at ~2^-95.
inline constexpr double kGammaExpm1Cut = 0x1.0p-10;

// Log1pmxDd's series: how many leading terms of T get dd treatment. A
// rounding introduced at Horner step k reaches T attenuated by u^k, so with
// |u| <= 1/16 the plain-double tail from k = 6 contributes
// 2^-53 * u^6 / (T * 8) ~ 2^-79 relative -- comfortably under the 2^-74
// series truncation the generator's self-check (g) certifies.
inline constexpr int kGammaPhiLead = 6;

// The leading six coefficients are carried as dd pairs
// (detail::kGammaPhiCoefLo holds the low words), and this is not
// belt-and-braces: rounding 1/3 to a double is 2^-55.9 absolute, which is
// 2^-58.5 RELATIVE on phi at the cut -- fine for phi, not for a*phi, the
// argument of an exponential (a relative eps on phi is an a*phi*eps
// relative error on the result, and a*phi reaches ~740 with |u| still at
// the cut when a ~ 4e5; the plain-double table measured 12 ULP at
// a = 3.8e5, lambda = 1.062). The dd leads put the coefficient
// contribution at 2^-79, leaving the 2^-75 truncation dominant -- the same
// lead/tail split as lgamma's zone polynomials, for the same reason. The
// generator self-checks each pair against the exact rational.

// R4's box in x. The alternating series' length kGammaR4N was sized by the
// generator's self-check (f) for exactly this bound (the 4^n/n! tail at
// x = 4 is what sets it), so the two must move together.
inline constexpr double kGammaR4XMax = 0x1.0p+2;

// R4's box in a: 1 + a must stay inside lgamma's centre-2 zone so that
// lgamma(1+a) can be evaluated at an EXACT argument (see GammaSmallQ).
inline constexpr double kGammaR4AMax = detail::kLgammaZoneHi - 1.0;  // 3/2

// Per-lane freeze threshold for the two summed series: a term below this
// times the running sum cannot move the dd result.
inline constexpr double kGammaFreezeEps = 0x1.0p-60;

// Temme core/tail split, exact: a*phi = 36 <=> z = sqrt(a*phi) = 6, the same
// threshold erfc's own core/tail split uses, so the two halves of this
// kernel consume exactly the two halves of erfc's machinery.
inline constexpr double kGammaTemmeZ2Split = 36.0;

// Ceiling on a when forming 1/sqrt(2*pi*a): 2*pi*a overflows above ~2.9e307
// and DdSqrt would then return NaN. Above this clamp the only non-saturated
// R3 lane possible is x == a exactly (any other x differs from a by at least
// ulp(a), which already puts a*phi far past kGammaExpFloor), and there the
// whole S/sqrt(2*pi*a) term is below 1e-150 against a result of 1/2 -- so
// the clamped value and the true one round identically.
inline constexpr double kGammaTwoPiAClamp = 0x1.0p+1000;

// --- indicator helpers ---------------------------------------------------
// The ops facade deliberately exposes no mask AND/OR (erfcinv's comment says
// the same), but the region map is a boolean expression over five
// comparisons. Carrying the predicates as 1.0/0.0 vectors makes AND a
// multiply, OR a max and NOT a subtract from one, all exact, and converts
// back to a mask with a single compare.
template <class D, class M>
HWY_INLINE op::V<D> Ind(D d, M m) {
  return op::IfThenElse(m, op::Set(d, 1.0), op::Zero(d));
}
template <class D>
HWY_INLINE op::V<D> IndNot(D d, op::V<D> v) {
  return op::Sub(op::Set(d, 1.0), v);
}
template <class D>
HWY_INLINE op::M<D> IndMask(D d, op::V<D> v) {
  return op::Gt(v, op::Set(d, 0.5));
}

// A region's contribution: the direct-side value rounded exactly once, and
// the same quantity as a dd so the other side can be formed as 1 (-) it
// before any rounding happens. For the regions that end in a power-of-two
// scaling the two are NOT the same computation -- see GammaScale.
template <class D>
struct GammaVal {
  op::V<D> v;
  Dd<D> dd;
};

// Apply exp_dd's power-of-two scaling LAST (the erfc pattern), twice over.
//   v  = ScaleTwo(round(m), e): the direct answer, ONE rounding even when the
//        result lands in the subnormal band.
//   dd = the scaled pair, for the complement. Where the complement is taken
//        the result is >= ~0.4 and the scaling is exact, so the pair is a
//        faithful dd there; where it is not, the complement rounds to 1
//        regardless.
template <class D>
HWY_INLINE GammaVal<D> GammaScale(D d, Dd<D> m, op::V<op::SignedTag<D>> e) {
  return {ScaleTwo(d, DdToDouble(m), e),
          Dd<D>{ScaleTwo(d, m.hi, e), ScaleTwo(d, m.lo, e)}};
}

// Exponential argument after the underflow guard, plus the mask saying the
// guard fired. See the file header's SATURATION note.
template <class D>
struct GammaExpArg {
  op::V<D> hi;
  op::V<D> lo;
  op::M<D> sat;
};

template <class D>
HWY_INLINE GammaExpArg<D> GammaClampE(D d, op::V<D> eh, op::V<D> el) {
  const auto zero = op::Zero(d);
  const auto floorv = op::Set(d, detail::kGammaExpFloor);
  // a*log(x) overflows to +/-inf for a above ~2.5e305 and its dd residual to
  // NaN with it. Every such lane is many hundreds of e-foldings past the
  // underflow point, so mapping NaN onto the floor is not a fudge: it selects
  // the answer (an exact 0) that the finite arithmetic was going to give.
  const auto bad = op::IsNaN(eh);
  eh = op::IfThenElse(bad, floorv, eh);
  el = op::IfThenElse(bad, zero, el);
  const auto sat = op::Ge(floorv, eh);
  // Clamping hi without clearing lo would hand ExpDdFrac an unnormalized
  // pair, which its argument reduction is not contracted to accept.
  return {op::Max(eh, floorv), op::IfThenElse(sat, zero, el), sat};
}

// ------------------------------------------------------------------------
// phi(u) = u - log1p(u), to dd precision RELATIVE to phi.
//
// This is the shared primitive the whole ridge rests on (the incomplete beta
// will reuse it). Note what it is NOT: u (-) LogDd(fl(1+u)). For small u,
// phi ~ u^2/2 while log1p(u) ~ u, so that spelling amplifies log_dd's ~2^-68
// by 2/u -- at u = 2^-20 that is 2^-49 relative, and multiplied by an a of
// 1e6 it is a wrong answer, not a rounding.
//
// SMALL |u| (<= kGammaPhiCut = 1/16): phi = u^2 * T(u) with
// T(u) = sum_k (-1)^k u^k / (k+2), the series in src/gamma_data.h. No
// cancellation anywhere -- u^2 is formed exactly (TwoProd) and T is O(1).
// Truncation at 18 terms is u^18/10 <= 2^-75 relative at the cut.
//
// The u.lo handling is the subtle part. Evaluating the series at u.hi and
// stopping there costs (2/3)*|u.lo| <= 2^-58 relative at the cut -- fine for
// phi itself, fatal once multiplied by an a large enough to make a*phi ~ 700.
// The fix is one first-order correction, and it is written on phi rather
// than on T because phi' = u/(1+u) is closed form:
//     phi(u.hi + u.lo) = phi(u.hi) + u.lo * u.hi/(1 + u.hi) + O(u.lo^2).
// The neglected term is (u.lo^2/2)/(1+u)^2, which against phi ~ u^2/2 is
// u.lo^2/u^2 <= 2^-106 relative -- uniformly, at every u.
//
// LARGE |u|: no cancellation to protect against (phi >= 2^-9 while u <= 1),
// so the direct difference is used, with 1 + u carried EXACTLY as a dd
// (TwoSum, plus u.lo folded into the low word) into LogDdAny. At the cut
// this inherits log_dd's error amplified by 2/u; log_dd near argument 1 is
// far better than its ~2^-68 global bound (the table's L_j is exact to
// 2^-107 and log1p(r) is relative to itself), which is what keeps the seam
// continuous in accuracy as well as in value.
template <class D>
HWY_INLINE Dd<D> Log1pmxDd(D d, Dd<D> u) {
  const auto one = op::Set(d, 1.0);
  const auto uh = u.hi;
  constexpr int kN =
      static_cast<int>(sizeof(detail::kGammaPhiCoef) / sizeof(double));
  const auto* c = detail::kGammaPhiCoef;

  // --- |u| <= 1/16 --------------------------------------------------------
  auto s = op::Set(d, c[kN - 1]);
  for (int k = kN - 2; k >= kGammaPhiLead; --k) {
    s = op::MulAdd(s, uh, op::Set(d, c[k]));
  }
  Dd<D> t{op::Set(d, c[kGammaPhiLead - 1]),
          op::Set(d, detail::kGammaPhiCoefLo[kGammaPhiLead - 1])};
  t = DdAddD(d, t, op::Mul(s, uh));
  for (int k = kGammaPhiLead - 2; k >= 0; --k) {
    t = DdAdd(d, DdMulD(d, t, uh),
              Dd<D>{op::Set(d, c[k]), op::Set(d, detail::kGammaPhiCoefLo[k])});
  }
  auto ser = DdMul(d, TwoProd(d, uh, uh), t);
  ser = DdAddD(d, ser, op::Mul(u.lo, op::Div(uh, op::Add(one, uh))));

  // --- |u| > 1/16 ---------------------------------------------------------
  auto w = TwoSum(d, one, uh);  // exact
  w.lo = op::Add(w.lo, u.lo);
  const auto lg = LogDdAny(d, w);
  const auto big = DdAdd(d, u, Dd<D>{op::Neg(lg.hi), op::Neg(lg.lo)});

  const auto m = op::Ge(op::Set(d, detail::kGammaPhiCut), op::Abs(uh));
  return Dd<D>{op::IfThenElse(m, ser.hi, big.hi),
               op::IfThenElse(m, ser.lo, big.lo)};
}

// e^w - 1 for a dd w, to dd precision relative to the result.
//
// Small |w| is the only interesting case: e^w (-) 1 there is a subtraction of
// two numbers within 2^-10 of each other, which throws away ten bits of the
// exponential's hundred. The series keeps them. Structured like Log1pDd:
// w and w^2/2 in dd (w^2/2 is 2^-11 relative to w, so a rounded w^2 would
// already show), everything from w^3 on in plain double, where |w^3/6| is
// 2^-23 relative and its own rounding lands at 2^-76.
template <class D>
HWY_INLINE Dd<D> Expm1Dd(D d, Dd<D> w) {
  const auto one = op::Set(d, 1.0);
  const auto wh = w.hi;

  auto q = op::Set(d, 1.0 / 720.0);
  q = op::MulAdd(q, wh, op::Set(d, 1.0 / 120.0));
  q = op::MulAdd(q, wh, op::Set(d, 1.0 / 24.0));
  q = op::MulAdd(q, wh, op::Set(d, 1.0 / 6.0));
  const auto sq = TwoProd(d, wh, wh);
  const auto half = op::Set(d, 0.5);
  // w^2/2 as a dd; halving is exact and the 2*w.hi*w.lo cross term is kept
  // for the same reason Log1pDd keeps it.
  const Dd<D> half_sq{op::Mul(sq.hi, half),
                      op::MulAdd(wh, w.lo, op::Mul(sq.lo, half))};
  auto ser = DdAdd(d, w, half_sq);
  ser = DdAddD(d, ser, op::Mul(op::Mul(sq.hi, wh), q));

  const auto ex = ExpDd(d, wh, w.lo);
  const auto big = DdAddD(d, ex, op::Neg(one));

  const auto m = op::Ge(op::Set(d, kGammaExpm1Cut), op::Abs(wh));
  return Dd<D>{op::IfThenElse(m, ser.hi, big.hi),
               op::IfThenElse(m, ser.lo, big.lo)};
}

// ------------------------------------------------------------------------
// R1: sum_{n>=0} x^n / ((a+1)(a+2)...(a+n)), the series factor of
//     P(a,x) = x^a e^{-x} / Gamma(a+1) * sum.
//
// Every term is positive and t_n / t_{n-1} = x/(a+n) <= 1 on the whole
// region (x <= a+1 in the small-a band, x <= a/2 in the large-a band), so
// the terms are non-increasing and the freeze mask below is monotone by
// construction. The cap kGammaSeriesN is a worst-case bound from the
// generator's self-check (d), not a guess: the deepest point in the region
// (a -> 20-, x = a+1) needs 52 terms.
//
// The freeze is per lane and STICKY -- once a lane's live flag drops it
// never comes back. A frozen lane's partial sum is then left literally
// untouched (the accumulator is SELECTED, not multiplied by a zero term):
// DdAddD of an exact zero is value-preserving but renormalizes, and a lane
// must not be able to tell how many renormalizations its neighbours' terms
// dragged it through. With the select, the vector-level break only skips
// iterations that would have been no-ops for every lane, so a point's answer
// is bit-identical evaluated alone, in a different lane position, or mixed
// with points from other regions -- which is what test_gamma's lane-mix
// check verifies.
//
// x is clamped to a+1, the region's own upper bound: a no-op on every live
// lane, and it keeps the term ratio at or below 1 for the DISCARDED lanes
// this core is also evaluated on (an x of 1e300 from an R2 lane would
// otherwise run the terms up to inf and, worse, never freeze).
template <class D>
HWY_INLINE Dd<D> GammaSeriesSum(D d, op::V<D> a, op::V<D> x) {
  const auto one = op::Set(d, 1.0);
  const auto zero = op::Zero(d);
  const auto half = op::Set(d, 0.5);
  const auto eps = op::Set(d, kGammaFreezeEps);
  const auto xc = op::Min(x, op::Add(a, one));

  auto t = one;
  Dd<D> s{one, zero};
  auto live = one;
  for (int n = 1; n <= detail::kGammaSeriesN; ++n) {
    t = op::Mul(t, op::Div(xc, op::Add(a, op::Set(d, static_cast<double>(n)))));
    const auto lm = op::Gt(live, half);
    const auto sn = DdAddD(d, s, t);
    s = Dd<D>{op::IfThenElse(lm, sn.hi, s.hi), op::IfThenElse(lm, sn.lo, s.lo)};
    live = op::IfThenElse(op::Lt(t, op::Mul(s.hi, eps)), zero, live);
    if (op::AllTrue(d, op::Eq(live, zero))) break;
  }
  return s;
}

// R2: 1/K, K the backward evaluation of Legendre's continued fraction
//     Q(a,x) = x^a e^{-x} / Gamma(a) * 1/K,
//     K = (x+1-a) - 1*(1-a)/((x+3-a) - 2*(2-a)/((x+5-a) - ...)).
//
// Started from the bare tail K_N = x + 2N + 1 - a and walked down to j = 1
// at a FIXED depth, with no convergence test at all -- see the file header
// for why a test is worse than useless here. 1/K is taken as a dd because Q
// is directly proportional to it, so a rounded reciprocal would put its own
// half ulp straight into the answer.
template <class D>
HWY_INLINE Dd<D> GammaCfRecip(D d, op::V<D> a, op::V<D> x) {
  auto k = op::Sub(
      op::Add(x, op::Set(d, static_cast<double>(2 * detail::kGammaCfN + 1))), a);
  for (int j = detail::kGammaCfN; j >= 1; --j) {
    const auto jv = op::Set(d, static_cast<double>(j));
    const auto num = op::Mul(jv, op::Sub(jv, a));
    const auto den =
        op::Sub(op::Add(x, op::Set(d, static_cast<double>(2 * j - 1))), a);
    k = op::Sub(den, op::Div(num, k));
  }
  return DdRecip(d, k);
}

// R4: Q(a,x) for 0 < a <= 3/2, 0 < x <= 4, computed DIRECTLY.
//
// From gamma(a,x) = x^a * sum_{n>=0} (-x)^n / (n! (a+n)) and Gamma(a) =
// Gamma(1+a)/a,
//     a*Gamma(a,x) = [Gamma(1+a) - x^a] - a*x^a*S,   S = sum_{n>=1} ...
//     Q(a,x)       = a*Gamma(a,x) / Gamma(1+a).
// Multiplying through by a (rather than forming Gamma(a,x) and dividing by
// Gamma(a)) removes the only 1/a in sight, which matters: a is allowed to be
// subnormal, where 1/a is +inf.
//
// THE CANCELLATION. As a -> 0 both Gamma(1+a) and x^a approach 1, and their
// difference is the entire answer; at x = e^-gamma they additionally cancel
// against each other to leading order in a. So neither is ever formed as
// "something minus one": Gamma(1+a) - 1 comes from Expm1Dd of lgamma(1+a),
// and x^a - 1 from Expm1Dd(a*log x), and it is those two dd expm1 values
// that get subtracted. The +1s are restored only afterwards, where they are
// harmless.
//
// THE EXACT ARGUMENT. lgamma(1+a) is NOT taken from LgammaPosDd(1+a): fl(1+a)
// rounds, and for a below 2^-53 it rounds to exactly 1, which would return a
// flat zero and lose the -gamma*a that is the whole signal. Instead lgamma's
// OWN zone polynomials are evaluated at the exact shifted argument -- t = a
// against centre 1 for a <= 1/2, t = a-1 (exact by Sterbenz) against centre 2
// for the rest -- which is exactly what lgamma-inl.h does internally for
// x in (0, 1/2) and [3/2, 5/2].
template <class D>
HWY_INLINE Dd<D> GammaSmallQ(D d, op::V<D> a, op::V<D> x) {
  const auto one = op::Set(d, 1.0);
  const auto zero = op::Zero(d);
  const auto half = op::Set(d, 0.5);

  // Domain clamp to R4's own box, exactly as LgammaPosDd clamps the argument
  // of the branch it is not taking: a no-op on every live lane, and it keeps
  // the discarded ones (this core runs on the whole vector) out of x^a = inf
  // and out of the zone polynomials' extrapolation range.
  a = op::Min(a, op::Set(d, kGammaR4AMax));
  x = op::Min(x, op::Set(d, kGammaR4XMax));

  const auto c1 = op::Ge(op::Set(d, detail::kLgammaZoneLo), a);  // a <= 1/2
  const auto t = op::IfThenElse(c1, a, op::Sub(a, one));
  const auto lg1a = DdMulD(d, ZoneBracket(d, t, c1), t);  // lgamma(1+a)
  const auto g1a_m1 = Expm1Dd(d, lg1a);                   // Gamma(1+a) - 1
  const auto g1a = DdAddD(d, g1a_m1, one);

  const auto xa_m1 = Expm1Dd(d, DdMulD(d, LogDdAny(d, x), a));  // x^a - 1
  const auto xa = DdAddD(d, xa_m1, one);

  // S = sum_{n>=1} (-x)^n / (n! (a+n)). Alternating with |x| <= 4, so terms
  // rise until n ~ x and then fall; the sticky freeze below is what makes
  // the early rise harmless. The 1/n weights are the dd table in
  // gamma_data.h (a rounded 1/n would be a 2^-53 relative error on a term
  // that can be the leading one), and 1/(a+n) is DdRecipDd of the EXACT
  // a+n (TwoSum -- a+n is not representable for a general a).
  const auto negx = op::Neg(x);
  const auto eps = op::Set(d, kGammaFreezeEps);
  Dd<D> term{one, zero};
  Dd<D> s{zero, zero};
  auto live = one;
  for (int n = 1; n <= detail::kGammaR4N; ++n) {
    term = DdMul(d, DdMulD(d, term, negx),
                 Dd<D>{op::Set(d, detail::kGammaRecipNHi[n - 1]),
                       op::Set(d, detail::kGammaRecipNLo[n - 1])});
    const auto wgt =
        DdRecipDd(d, TwoSum(d, a, op::Set(d, static_cast<double>(n))));
    const auto contrib = DdMul(d, term, wgt);
    const auto lm = op::Gt(live, half);  // sticky; see GammaSeriesSum
    const auto sn = DdAdd(d, s, contrib);
    s = Dd<D>{op::IfThenElse(lm, sn.hi, s.hi), op::IfThenElse(lm, sn.lo, s.lo)};
    live = op::IfThenElse(
        op::Lt(op::Abs(contrib.hi), op::Mul(op::Abs(s.hi), eps)), zero, live);
    if (op::AllTrue(d, op::Eq(live, zero))) break;
  }

  auto num = DdAdd(d, g1a_m1, Dd<D>{op::Neg(xa_m1.hi), op::Neg(xa_m1.lo)});
  const auto axs = DdMulD(d, DdMul(d, xa, s), a);
  num = DdAdd(d, num, Dd<D>{op::Neg(axs.hi), op::Neg(axs.lo)});
  return DdMul(d, num, DdRecipDd(d, g1a));
}

// R3: Temme's uniform asymptotic expansion on the ridge.
//
// With lambda = x/a and phi(lambda) = lambda - 1 - log(lambda) >= 0, the
// smaller of the pair is
//     small = 1/2*erfc(z) + sgn * e^{-a*phi}/sqrt(2*pi*a) * S(eta, 1/a),
//     z = sqrt(a*phi),  eta = sgn*sqrt(2*phi),  sgn = sign(x - a),
// with sgn = +1 selecting Q (x >= a) and sgn = -1 selecting P. S is the
// Chebyshev table in gamma_data.h: Clenshaw down each of the 11 rows in the
// mapped variable, then Horner in 1/a across the rows. The coefficients were
// extracted clean-room from the oracle by a Vandermonde solve in 1/a; they
// absorb the Stirling remainder, which is why this region touches no lgamma.
//
// FIVE THINGS THAT WOULD OTHERWISE COST THE REGION ITS ACCURACY:
//  1. lambda - 1 is formed as (x - a)/a with x - a EXACT (Sterbenz holds on
//     the whole region, since lambda in (1/2, 2)) and 1/a as a dd. Near the
//     ridge that difference IS the signal; fl(x/a) - 1 would round it away.
//  2. a*phi comes from Log1pmxDd, never from a naive lambda - 1 - log lambda.
//  3. e^{-z^2} is e^{-a*phi} taken from the dd a*phi -- NEVER from squaring
//     the rounded z. Near z = 6 a half ulp of z is 2^-48 in z^2, which is
//     2^-48 RELATIVE in the exponential; the dd a*phi is 2^-64 there.
//  4. z is a dd and its low word is the first-order correction to erfc:
//     1/2*erfc(z.hi + z.lo) = 1/2*erfc(z.hi) - z.lo/sqrt(pi) * e^{-z^2}, so
//     the correction folds into the same bracket S already multiplies.
//  5. The core/tail split is on a*phi against an exact 36.0, i.e. z against
//     6 -- erfc's own split. Below it ErfcCoreDd supplies 1/2*erfc(z) as a
//     dd and nothing underflows. Above it the whole expression carries the
//     factor e^{-a*phi}, so erfc is taken in its tail form
//     e^{-z^2} G(1/z)/z (REUSING erfc_tail_data's G -- z <= sqrt(800) ~ 28.3
//     stays inside its fitted range) and the power-of-two scaling is applied
//     last, which is what keeps the subnormal band a*phi in [700, 760] at a
//     single rounding.
template <class D>
struct GammaTemmeOut {
  GammaVal<D> val;
  op::V<D> is_p;  // 1.0 where the directly computed side is P
  op::M<D> sat;
};

template <class D>
HWY_INLINE GammaTemmeOut<D> GammaTemme(D d, op::V<D> a, op::V<D> x) {
  const auto one = op::Set(d, 1.0);
  const auto zero = op::Zero(d);

  // Domain clamp (no-op on every live lane): keeps 1/a finite for the
  // discarded lanes this core also runs on, a subnormal a among them.
  a = op::Max(a, op::Set(d, detail::kGammaAT));

  const auto xma = op::Sub(x, a);  // exact by Sterbenz on this region
  const auto sgn = op::CopySign(one, xma);
  const auto u = DdMul(d, TwoSum(d, x, op::Neg(a)), DdRecip(d, a));
  const auto phi = Log1pmxDd(d, u);
  const auto aphi = DdMulD(d, phi, a);

  const Dd<D> two_phi{op::Add(phi.hi, phi.hi), op::Add(phi.lo, phi.lo)};
  const auto eta = op::CopySign(DdSqrt(d, two_phi).hi, xma);
  const auto z = DdSqrt(d, aphi);

  // S(eta, 1/a). The Chebyshev variable maps the fitted eta band onto
  // [-1, 1]; the reciprocal of the half-width is a compile-time constant.
  const auto tc = op::Mul(op::Sub(eta, op::Set(d, detail::kGammaEtaMid)),
                          op::Set(d, 1.0 / detail::kGammaEtaHalf));
  const auto tc2 = op::Add(tc, tc);
  const auto r = op::Div(one, a);
  auto sum = zero;
  for (int k = detail::kGammaTemmeK - 1; k >= 0; --k) {
    const double* row = detail::kGammaTemmeCheb[k];
    auto b1 = zero;
    auto b2 = zero;
    for (int j = detail::kGammaTemmeNCoef - 1; j >= 1; --j) {
      const auto nb = op::Sub(op::MulAdd(tc2, b1, op::Set(d, row[j])), b2);
      b2 = b1;
      b1 = nb;
    }
    sum = op::MulAdd(sum, r,
                     op::Sub(op::MulAdd(tc, b1, op::Set(d, row[0])), b2));
  }

  const Dd<D> twopi{op::Set(d, detail::kGammaTwoPiHi),
                    op::Set(d, detail::kGammaTwoPiLo)};
  const auto rv = DdRecipDd(
      d, DdSqrt(d, DdMulD(d, twopi,
                          op::Min(a, op::Set(d, kGammaTwoPiAClamp)))));
  const auto s_rv = DdMulD(d, rv, op::Mul(sgn, sum));

  const auto ea = GammaClampE(d, op::Neg(aphi.hi), op::Neg(aphi.lo));

  // --- core, a*phi <= 36 --------------------------------------------------
  // ErfcCoreDd's table index is round(min(z,6+1/1024)*256) and is NOT masked,
  // so a NaN z from a discarded lane must be scrubbed before it -- the same
  // one-op guard erfinv's HalleyMid carries (AGENTS.md, value-derived
  // gathers). Only discarded lanes are affected.
  const auto zs = op::IfThenElse(op::IsNaN(z.hi), zero, z.hi);
  const auto ec = ErfcCoreDd(d, zs, zs);
  const auto half = op::Set(d, 0.5);
  const Dd<D> half_erfc{op::Mul(ec.hi, half), op::Mul(ec.lo, half)};
  const auto exd = ExpDd(d, ea.hi, ea.lo);
  const auto brk_core = DdAddD(
      d, s_rv, op::Neg(op::Mul(z.lo, op::Set(d, detail::kGammaInvSqrtPiHi))));
  const auto core = DdAdd(d, half_erfc, DdMul(d, brk_core, exd));

  // --- tail, a*phi > 36 ---------------------------------------------------
  const auto uz = op::Div(one, zs);
  const auto gt = op::Div(ErfcTailGFromU(d, zs, uz), op::Add(zs, zs));
  const auto brk_tail = DdAddD(d, s_rv, gt);
  const auto exf = ExpDdFrac(d, ea.hi, ea.lo);
  const auto tail = GammaScale(d, DdMul(d, exf.m, brk_tail), exf.e);

  const auto tm = op::Gt(aphi.hi, op::Set(d, kGammaTemmeZ2Split));
  GammaVal<D> val{op::IfThenElse(tm, tail.v, DdToDouble(core)),
                  Dd<D>{op::IfThenElse(tm, tail.dd.hi, core.hi),
                        op::IfThenElse(tm, tail.dd.lo, core.lo)}};
  return {val, Ind(d, op::Lt(xma, zero)), ea.sat};
}

// ------------------------------------------------------------------------
// The driver. kP selects gamma_p (true) or gamma_q (false); the two differ
// only in who owns the R1/R4 overlap and in which side is complemented.
template <bool kP, class D>
HWY_INLINE op::V<D> GammaVec(D d, op::V<D> a_in, op::V<D> x_in) {
  const auto one = op::Set(d, 1.0);
  const auto zero = op::Zero(d);
  const auto inf = op::Set(d, std::numeric_limits<double>::infinity());
  const auto qnan = op::Set(d, std::numeric_limits<double>::quiet_NaN());

  // --- scrub every lane whose result is decided by the specials table -----
  // Masked-off lanes still execute every op (AGENTS.md), and this kernel
  // divides by a, takes logs of x and multiplies j*(j-a): a NaN, a zero or
  // an infinity left in place would propagate through gathers and dd
  // residuals rather than sit quietly. (1, 3) is an ordinary interior point
  // of every region's domain.
  const auto safe_a = one;
  const auto safe_x = op::Set(d, 3.0);
  auto a = a_in;
  auto x = x_in;
  {
    // NaN in either operand (a + x also catches inf + -inf).
    const auto m = op::IsNaN(op::Add(a, x));
    a = op::IfThenElse(m, safe_a, a);
    x = op::IfThenElse(m, safe_x, x);
  }
  {
    const auto m = op::Ge(zero, op::Min(a, x));  // a <= 0 or x <= 0
    a = op::IfThenElse(m, safe_a, a);
    x = op::IfThenElse(m, safe_x, x);
  }
  {
    const auto m = op::Eq(op::Max(a, x), inf);
    a = op::IfThenElse(m, safe_a, a);
    x = op::IfThenElse(m, safe_x, x);
  }

  // --- region map ---------------------------------------------------------
  // lambda <= 1/2 and lambda >= 2 are tested as 2x <= a and x >= 2a: the
  // doublings are exact, and where one overflows to inf the comparison is
  // still right (no double is >= 2a for a above 2^1023).
  const auto i_small = Ind(d, op::Lt(a, op::Set(d, detail::kGammaAT)));
  const auto i_xle = Ind(d, op::Ge(op::Add(a, one), x));  // x <= a+1
  const auto i_lo = Ind(d, op::Ge(a, op::Add(x, x)));
  const auto i_hi = Ind(d, op::Ge(x, op::Add(a, a)));
  const auto i_box = op::Mul(Ind(d, op::Ge(op::Set(d, kGammaR4AMax), a)),
                             Ind(d, op::Ge(op::Set(d, kGammaR4XMax), x)));

  const auto i_big = IndNot(d, i_small);
  auto i_r1 = op::MulAdd(i_small, i_xle, op::Mul(i_big, i_lo));
  auto i_r2 = op::MulAdd(i_small, IndNot(d, i_xle), op::Mul(i_big, i_hi));
  const auto i_r3 = op::Mul(i_big, op::Mul(IndNot(d, i_lo), IndNot(d, i_hi)));
  // R4 is a subset of R1 u R2 (a <= 3/2 < 20), so whichever of the two it
  // takes over from loses exactly those lanes.
  auto i_r4 = i_box;
  if (kP) {
    i_r4 = op::Mul(i_box, IndNot(d, i_xle));  // R1 keeps x <= a+1
    i_r2 = op::Mul(i_r2, IndNot(d, i_r4));
  } else {
    i_r1 = op::Mul(i_r1, IndNot(d, i_r4));
    i_r2 = op::Mul(i_r2, IndNot(d, i_r4));
  }

  const auto m_r1 = IndMask(d, i_r1);
  const auto m_r2 = IndMask(d, i_r2);
  const auto m_r3 = IndMask(d, i_r3);
  const auto m_r4 = IndMask(d, i_r4);

  // --- assembled per-lane direct side -------------------------------------
  GammaVal<D> val{zero, Dd<D>{zero, zero}};
  auto is_p = zero;

  const bool need12 = !op::AllFalse(d, m_r1) || !op::AllFalse(d, m_r2);
  if (need12) {
    // E = a*log(x) - x - lgamma(a), the log of the common prefactor.
    auto eb = DdMulD(d, LogDdAny(d, x), a);
    eb = DdAddD(d, eb, op::Neg(x));
    const auto lg = LgammaPosDd(d, a);
    eb = DdAdd(d, eb, Dd<D>{op::Neg(lg.hi), op::Neg(lg.lo)});

    // Saturated lanes are scrubbed before the series and the CF -- the CF's
    // j*(j-a) overflows to -inf for a above ~4e306, which is an ordinary
    // (if utterly saturated) R2 argument. Each region scrubs on ITS OWN
    // saturation mask: R1's exponent is E - log a, which for tiny a is
    // hundreds of e-foldings ABOVE E, so sharing R2's mask would scrub live
    // R1 lanes.
    if (!op::AllFalse(d, m_r1)) {
      // (-) log a turns Gamma(a) into Gamma(a+1) without forming 1 + a and
      // without a 1/a that would be +inf for subnormal a. The two logs
      // cancel to ~0 for tiny a, but only their ABSOLUTE error matters here
      // (it is the argument of an exponential), and that is ~2^-89.
      const auto la = LogDdAny(d, a);
      const auto e1 =
          DdAdd(d, eb, Dd<D>{op::Neg(la.hi), op::Neg(la.lo)});
      const auto ea = GammaClampE(d, e1.hi, e1.lo);
      const auto as = op::IfThenElse(ea.sat, safe_a, a);
      const auto xs = op::IfThenElse(ea.sat, safe_x, x);
      const auto ex = ExpDdFrac(d, ea.hi, ea.lo);
      auto r1 = GammaScale(d, DdMul(d, ex.m, GammaSeriesSum(d, as, xs)), ex.e);
      r1.v = op::IfThenElse(ea.sat, zero, r1.v);
      r1.dd.hi = op::IfThenElse(ea.sat, zero, r1.dd.hi);
      r1.dd.lo = op::IfThenElse(ea.sat, zero, r1.dd.lo);
      val.v = op::IfThenElse(m_r1, r1.v, val.v);
      val.dd.hi = op::IfThenElse(m_r1, r1.dd.hi, val.dd.hi);
      val.dd.lo = op::IfThenElse(m_r1, r1.dd.lo, val.dd.lo);
      is_p = op::IfThenElse(m_r1, one, is_p);
    }
    if (!op::AllFalse(d, m_r2)) {
      const auto ea = GammaClampE(d, eb.hi, eb.lo);
      const auto as = op::IfThenElse(ea.sat, safe_a, a);
      const auto xs = op::IfThenElse(ea.sat, safe_x, x);
      const auto ex = ExpDdFrac(d, ea.hi, ea.lo);
      auto r2 = GammaScale(d, DdMul(d, ex.m, GammaCfRecip(d, as, xs)), ex.e);
      r2.v = op::IfThenElse(ea.sat, zero, r2.v);
      r2.dd.hi = op::IfThenElse(ea.sat, zero, r2.dd.hi);
      r2.dd.lo = op::IfThenElse(ea.sat, zero, r2.dd.lo);
      val.v = op::IfThenElse(m_r2, r2.v, val.v);
      val.dd.hi = op::IfThenElse(m_r2, r2.dd.hi, val.dd.hi);
      val.dd.lo = op::IfThenElse(m_r2, r2.dd.lo, val.dd.lo);
    }
  }
  if (!op::AllFalse(d, m_r3)) {
    auto t = GammaTemme(d, a, x);
    t.val.v = op::IfThenElse(t.sat, zero, t.val.v);
    t.val.dd.hi = op::IfThenElse(t.sat, zero, t.val.dd.hi);
    t.val.dd.lo = op::IfThenElse(t.sat, zero, t.val.dd.lo);
    val.v = op::IfThenElse(m_r3, t.val.v, val.v);
    val.dd.hi = op::IfThenElse(m_r3, t.val.dd.hi, val.dd.hi);
    val.dd.lo = op::IfThenElse(m_r3, t.val.dd.lo, val.dd.lo);
    is_p = op::IfThenElse(m_r3, t.is_p, is_p);
  }
  if (!op::AllFalse(d, m_r4)) {
    const auto q4 = GammaSmallQ(d, a, x);
    val.v = op::IfThenElse(m_r4, DdToDouble(q4), val.v);
    val.dd.hi = op::IfThenElse(m_r4, q4.hi, val.dd.hi);
    val.dd.lo = op::IfThenElse(m_r4, q4.lo, val.dd.lo);
  }

  // The complement is formed from the dd BEFORE any rounding, so the single
  // rounding of 1 - direct is the only one it carries.
  const auto comp = DdToDouble(
      DdAdd(d, Dd<D>{one, zero}, Dd<D>{op::Neg(val.dd.hi), op::Neg(val.dd.lo)}));
  const auto m_p = IndMask(d, is_p);
  auto res = kP ? op::IfThenElse(m_p, val.v, comp)
                : op::IfThenElse(m_p, comp, val.v);

  // --- specials (SciPy limits; see the public header) ---------------------
  // Applied last, in increasing priority. The two conjunctions are written
  // arithmetically because a + x is zero only when both are (both are known
  // non-negative by the time it decides anything) and min(a,x) is infinite
  // only when both are.
  const auto p_at = kP ? one : zero;   // value at "x exhausted": P=1, Q=0
  const auto p_at0 = kP ? zero : one;  // value at "no mass yet":  P=0, Q=1
  res = op::IfThenElse(op::Eq(x_in, zero), p_at0, res);
  res = op::IfThenElse(op::Eq(a_in, zero), p_at, res);
  res = op::IfThenElse(op::Eq(x_in, inf), p_at, res);
  res = op::IfThenElse(op::Eq(a_in, inf), p_at0, res);
  res = op::IfThenElse(op::Eq(op::Add(a_in, x_in), zero), qnan, res);
  res = op::IfThenElse(op::Eq(op::Min(a_in, x_in), inf), qnan, res);
  res = op::IfThenElse(op::Lt(a_in, zero), qnan, res);
  res = op::IfThenElse(op::Lt(x_in, zero), qnan, res);
  res = op::IfThenElse(op::IsNaN(a_in), a_in, res);  // payload preserved
  res = op::IfThenElse(op::IsNaN(x_in), x_in, res);
  return res;
}

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#endif  // include guard
