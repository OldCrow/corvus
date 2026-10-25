// Regularized incomplete gamma P(a,x) and Q(a,x) = 1 - P(a,x).
// Per-target include guard (Highway -inl.h idiom).
//
// Both public functions share one kernel: four region cores plus a router
// that differs between them in exactly two places -- which core owns the
// R1/R4 overlap, and which side of the pair gets complemented. The region
// map and the fixed series/CF lengths are probe-validated and reproduced
// here, not re-derived;
// the tables and every threshold live in src/gamma_data.h.
//
// THE ONE RULE THE WHOLE DESIGN TURNS ON: always compute the SMALLER of the
// pair directly, and get the other one as 1 (-) it, as a double-double
// subtraction rounded once. A directly computed side keeps full relative
// accuracy however tiny it gets. The complement is NOT always large: in R1
// with a just above 3/2 and x near a+1, Q is ~0.13 (#56), which amplifies
// P's error ~5x. What keeps it correctly rounded is that the complement is
// taken from the direct side's DD, before any rounding -- so the direct
// side's dd must be accurate to well below a double ulp, which is why R1's
// terms and R2's last CF steps are carried in dd. Computing the large side
// directly and subtracting would be the opposite trade and is never done.
//
// REGIONS (lambda = x/a, a_T = kGammaAT = 20)
//   R1  series, P-direct   {a < 20, x <= a+1} u {a >= 20, lambda <= 1/2}
//       P = e^E' * sum_{n>=0} x^n / ((a+1)...(a+n)),  E' = E - log a.
//       Fixed cap kGammaSeriesN with a per-lane freeze mask; terms AND sum
//       are dd (#56). A plain-double term's rounding is NOT attenuated:
//       each term is a factor of the next, so a rounding at step k rides
//       into every later term, and near the ridge (x/(a+n) ~ 1) those
//       terms do not decay -- the plain form measured 3.11 ulp on P and 6.08
//       on its complement.
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
#include "src/dd_special-inl.h"
#include "src/erfc_core-inl.h"
#include "src/erfc_tail_dd-inl.h"  // GammaTailG, shared with beta (#69)
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
// GammaCfRecip's double-double tail (#56): the last kGammaCfDdSteps steps of
// the backward CF run in dd, the deeper ones in double. The recurrence
// contracts the deep steps' rounding, so only the last few reach the
// answer: 2 steps left R2 at 0.66 ulp on a 6,000-point ridge probe, 4
// steps are correctly rounded there and bit-identical to 8 on all but 9 of
// 20,000 R2 probe points (all half-ulp ties, both within 0.51 ulp), at ~3%
// cost against 8's ~25% (Zen 4, indicative).
inline constexpr int kGammaCfDdSteps = 4;

// Temme core/tail split, exact: a*phi = 36 <=> z = sqrt(a*phi) = 6, the same
// threshold erfc's own core/tail split uses, so the two halves of this
// kernel consume exactly the two halves of erfc's machinery.
inline constexpr double kGammaTemmeZ2Split = 36.0;

// Ceiling on a when forming 1/sqrt(2*pi*a): 2*pi*a overflows above ~2.9e307
// and DdSqrt would then return NaN. Above this clamp the only non-saturated
// R3 lane possible is x == a exactly (any other x differs from a by at least
// ulp(a), which already puts a*phi far past kGammaExpFloor), and there the
// whole S/sqrt(2*pi*a) term is below 2^-450 against a result of 1/2 -- so
// the clamped value and the true one round identically. Must also sit under
// ops::ProdLow's 2^996 non-FMA Dekker ceiling minus the 2^27 split factor:
// anything above that overflows the Dekker split of the clamped a on
// non-FMA tiers (the u and a*phi sites clear the same ceiling with the
// prescale below).
inline constexpr double kGammaTwoPiAClamp = 0x1.0p+900;

// Exact prescale for GammaTemme's Dekker-split operands (BetaPsiCore's
// pattern, beta-inl.h), shared by the gammainv forward's twin site.
// u = (x - a)/a is 0-homogeneous in (x, a) and a*phi is 1-homogeneous, so
// scaling both operands down by an exact power of two and scaling a*phi
// back up is bit-identical wherever the unscaled arithmetic was finite:
// power-of-two scaling commutes with every rounding in TwoSum, DdRecip and
// DdMulD, and no intermediate over/underflows (scaled a lands in
// (2^700, 2^824]). Without it, ops::ProdLow's non-FMA Dekker split of a
// overflows for a > 2^996 and GammaClampE laundered the NaN into a
// saturated (1, 0) pair on the diagonal x = a -- the one unsaturated point
// of that band, whose true value is ~ 1/2 (#12).
inline constexpr double kGammaScaleAbove = 0x1.0p+900;
inline constexpr double kGammaScaleDown = 0x1.0p-200;
inline constexpr double kGammaScaleUp = 0x1.0p+200;

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
//        the direct side is at most ~0.87 (Q ~ 0.13 at worst, #56) and far
//        from the subnormal band, so the scaling is exact and the pair is a
//        faithful dd; where the direct side is subnormal the complement
//        rounds to 1 regardless.
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
// R1: sum_{n>=0} x^n / ((a+1)(a+2)...(a+n)), the series factor of
//     P(a,x) = x^a e^{-x} / Gamma(a+1) * sum.
//
// Every term is positive and t_n / t_{n-1} = x/(a+n) <= 1 on the whole
// region (x <= a+1 in the small-a band, x <= a/2 in the large-a band), so
// the terms are non-increasing and the freeze mask below is monotone by
// construction. Terms and ratio are dd (#56): a+n is an exact TwoSum and
// x/(a+n) = x * DdRecipDd(a+n), so no step adds a double's rounding to the
// chain every later term inherits (the cost: ~40% on R1, Zen 4,
// indicative). The cap kGammaSeriesN is a worst-case bound from the
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
//
// All four region cores are HWY_NOINLINE: fully inlined, the
// two exports each become one enormous function per target and MSVC's
// optimizer is superlinear in function size -- large enough to run a CI
// job past its timeout inside code generation.
// Outlining the cores caps the largest
// function the optimizer sees. Results are bit-identical (contraction is
// off, so inlining cannot change FP semantics) and the call cost is a few
// cycles per vector against the cores' hundreds of flops.
template <class D>
HWY_NOINLINE Dd<D> GammaSeriesSum(D d, op::V<D> a, op::V<D> x) {
  const auto one = op::Set(d, 1.0);
  const auto zero = op::Zero(d);
  const auto half = op::Set(d, 0.5);
  const auto eps = op::Set(d, kGammaFreezeEps);
  const auto xc = op::Min(x, op::Add(a, one));

  Dd<D> t{one, zero};
  Dd<D> s{one, zero};
  auto live = one;
  for (int n = 1; n <= detail::kGammaSeriesN; ++n) {
    const auto apn = TwoSum(d, a, op::Set(d, static_cast<double>(n)));  // exact
    // FLOOR: tiny x (or a decayed term) puts these below it; each step loses
    // <= 2^-1071 absolute against s >= 1 (terms >= 0, ratio <= 1), <= 2^-1065
    // relative over 64 steps; every caller (gamma, gammainv, beta, betainv)
    // scales s or takes its log, so it stays relative.
    t = DdMul(d, t, DdMulD(d, DdRecipDd(d, apn), xc));
    const auto lm = op::Gt(live, half);
    const auto sn = DdAdd(d, s, t);
    s = Dd<D>{op::IfThenElse(lm, sn.hi, s.hi), op::IfThenElse(lm, sn.lo, s.lo)};
    live = op::IfThenElse(op::Lt(t.hi, op::Mul(s.hi, eps)), zero, live);
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
// for why a test is worse than useless here. The last kGammaCfDdSteps steps
// run in dd (#56): the recurrence contracts what the deep steps round, but
// the last few each put a plain double's half ulp into K, which is worth
// ~2 ulp on Q near the ridge. j - a and x + 2j - 1 - a are formed as exact
// TwoSums for the same reason. 1/K is a dd (DdRecipDd) because Q is
// directly proportional to it.
template <class D>
HWY_NOINLINE Dd<D> GammaCfRecip(D d, op::V<D> a, op::V<D> x) {
  const auto ma = op::Neg(a);
  auto kd = op::Sub(
      op::Add(x, op::Set(d, static_cast<double>(2 * detail::kGammaCfN + 1))), a);
  for (int j = detail::kGammaCfN; j > kGammaCfDdSteps; --j) {
    const auto jv = op::Set(d, static_cast<double>(j));
    const auto num = op::Mul(jv, op::Sub(jv, a));
    const auto den =
        op::Sub(op::Add(x, op::Set(d, static_cast<double>(2 * j - 1))), a);
    kd = op::Sub(den, op::Div(num, kd));
  }
  Dd<D> k{kd, op::Zero(d)};
  for (int j = kGammaCfDdSteps; j >= 1; --j) {
    const auto jv = op::Set(d, static_cast<double>(j));
    const auto num = DdMulD(d, TwoSum(d, jv, ma), jv);  // j*(j-a)
    const auto den = DdAddD(
        d, TwoSum(d, x, op::Set(d, static_cast<double>(2 * j - 1))), ma);
    // FLOOR: |num| = j|j-a| is 0 (exact) or >= 2^-52 (j - a exact), so
    // q < 2^-968 needs |k| > 2^916, i.e. x > ~2^915 (bracket iterates); K =
    // den - q is then as large, and q's <= 2^-1072 loss is ~2^-1987 of K.
    const auto q = DdMul(d, num, DdRecipDd(d, k));
    k = DdAdd(d, den, Dd<D>{op::Neg(q.hi), op::Neg(q.lo)});
  }
  return DdRecipDd(d, k);
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
HWY_NOINLINE Dd<D> GammaSmallQ(D d, op::V<D> a, op::V<D> x) {
  const auto one = op::Set(d, 1.0);
  const auto zero = op::Zero(d);
  const auto half = op::Set(d, 0.5);

  // Domain clamp to R4's own box, exactly as LgammaPosDd clamps the argument
  // of the branch it is not taking: a no-op on every live lane, and it keeps
  // the discarded ones (this core runs on the whole vector) out of x^a = inf
  // and out of the zone polynomials' extrapolation range.
  a = op::Min(a, op::Set(d, kGammaR4AMax));
  x = op::Min(x, op::Set(d, kGammaR4XMax));

  // lgamma(1+a), always from the CENTRE-2 zone (#67): t = a - 1 (Sterbenz)
  // for a > 1/2, and for a <= 1/2 t = a with lgamma(1+a) = lgamma(2+a) -
  // log(1+a). The centre-1 zone is ~2^-57 relative; Q's error is lg1a's
  // ABSOLUTE error over Q, and Q reaches ~3e-3 in this box, so centre-1
  // cost up to 4 ulp at a ~ 0.43, x ~ 4. Centre-2 is 50-100x tighter there.
  // log(1+a) must be RELATIVE-accurate to dd: for tiny a the O(a) terms of
  // num cancel to Q ~ a*E1(x), so a 2^-53 relative error in log(1+a) costs
  // ~100 ulp. For a <= 2^-7 it is Log1pDd(a) (dd, relative to itself). Above,
  // 1 + a is exact as TwoSum(1, a) and LogDd of that dd is accurate: its
  // plain-double l/h correction is <= 2^-53 against a result >= 2^-7. Then
  // lg1a = a*(B2(a) - 1) keeps its relative accuracy down to subnormal a
  // (the cancellation is ~1.7x).
  const auto c1 = op::Ge(op::Set(d, detail::kLgammaZoneLo), a);  // a <= 1/2
  CORVUS_ROUTE(d, "gamma/R4/a_le_half", c1);
  CORVUS_ROUTE(d, "gamma/R4/a_gt_half",
               op::Gt(a, op::Set(d, detail::kLgammaZoneLo)));
  const auto t = op::IfThenElse(c1, a, op::Sub(a, one));
  const auto centre2 = op::Gt(zero, one);  // all false: centre-2 coefficients
  // FLOOR: (R4 note, for every site in this function) |B*t| < 2^-968 needs
  // a < ~2^-967 (|B| >= 0.42); a*log x, axs and q reach it for tiny a too,
  // each losing <= 2^-1072 absolute (~2^-1070 in all) against Q >= a*E1(4)
  // ~ a*2^-8. Tiny x puts the S terms below it; that loss (<= 2^-1065) is
  // scaled by a*xa, so it stays ~2^-1053 relative to Q. Harmless for a >=
  // 2^-1000 (measured: no tier split on 300 rows a in [2^-1012, 2^-940]),
  // and for gamma_p (1 - Q). Below a = 2^-1000 it is up to ~3 ulp of a
  // subnormal Q on no-FMA tiers -- tracked by #71, and #83 for gamma_q_inv.
  auto lg1a = DdMulD(d, ZoneBracket(d, t, centre2), t);  // lgamma(2+t)
  if (!op::AllFalse(d, c1)) {
    CORVUS_ROUTE_SCOPE(d, c1);
    const auto r_cut = op::Set(d, 0x1.0p-7);  // Log1pDd's domain
    const auto tiny = op::Ge(r_cut, a);
    CORVUS_ROUTE(d, "gamma/R4/a_le_half/log", op::Gt(a, r_cut));
    Dd<D> l1a = LogDd(d, TwoSum(d, one, op::Max(a, r_cut)));  // 1+a exact
    if (!op::AllFalse(d, tiny)) {
      CORVUS_ROUTE_SCOPE(d, tiny);
      CORVUS_ROUTE(d, "gamma/R4/a_le_half/log1p", tiny);
      const auto ls = Log1pDd(d, Dd<D>{op::Min(a, r_cut), zero});
      l1a = Dd<D>{op::IfThenElse(tiny, ls.hi, l1a.hi),
                  op::IfThenElse(tiny, ls.lo, l1a.lo)};
    }
    const auto sub = DdAdd(d, lg1a, Dd<D>{op::Neg(l1a.hi), op::Neg(l1a.lo)});
    lg1a = Dd<D>{op::IfThenElse(c1, sub.hi, lg1a.hi),
                 op::IfThenElse(c1, sub.lo, lg1a.lo)};
  }
  const auto g1a_m1 = Expm1Dd(d, lg1a);                   // Gamma(1+a) - 1
  const auto g1a = DdAddD(d, g1a_m1, one);

  // FLOOR: |a*log x| < 2^-968 needs a < 2^-915 (log x is 0 exactly or in
  // [2^-53, 745]); see the R4 note at lg1a (#71 below a = 2^-1000).
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
    // FLOOR: live terms fall below it only for x < ~2^-900 (contrib >=
    // 2^-60 |S|, S ~ x); see the R4 note at lg1a.
    term = DdMul(d, DdMulD(d, term, negx),
                 Dd<D>{op::Set(d, detail::kGammaRecipNHi[n - 1]),
                       op::Set(d, detail::kGammaRecipNLo[n - 1])});
    const auto wgt =
        DdRecipDd(d, TwoSum(d, a, op::Set(d, static_cast<double>(n))));
    // FLOOR: as the term above (wgt <= 1); see the R4 note at lg1a.
    const auto contrib = DdMul(d, term, wgt);
    const auto lm = op::Gt(live, half);  // sticky; see GammaSeriesSum
    const auto sn = DdAdd(d, s, contrib);
    s = Dd<D>{op::IfThenElse(lm, sn.hi, s.hi), op::IfThenElse(lm, sn.lo, s.lo)};
    live = op::IfThenElse(
        op::Lt(op::Abs(contrib.hi), op::Mul(op::Abs(s.hi), eps)), zero, live);
    if (op::AllTrue(d, op::Eq(live, zero))) break;
  }

  auto num = DdAdd(d, g1a_m1, Dd<D>{op::Neg(xa_m1.hi), op::Neg(xa_m1.lo)});
  // FLOOR: below it for tiny S (x < ~2^-900) or a; see the R4 note at lg1a.
  const auto axs = DdMulD(d, DdMul(d, xa, s), a);
  num = DdAdd(d, num, Dd<D>{op::Neg(axs.hi), op::Neg(axs.lo)});
  // FLOOR: num ~ Q < 2^-968 needs a < ~2^-958: the loss is <= 2^-1072 on Q,
  // up to ~3 ulp of a subnormal Q on no-FMA -- tracked by #71 (R4 note).
  const auto q = DdMul(d, num, DdRecipDd(d, g1a));
  // SIGN CLAMP (#51). Q > 0 always, but for a subnormal a the three O(a)
  // terms of num cancel to a true Q far below denorm_min (gamma_q(2^-1074, 2)
  // = 2.4e-325), and their rounding residue can come out as -denorm_min or
  // -0. Clamping q.hi <= 0 to +0 is value-safe -- it only fires where the
  // true Q rounds to +0 anyway -- and keeps the contract that a probability
  // is never negative, signed zero included.
  const auto nonpos = op::Ge(zero, q.hi);
  CORVUS_ROUTE(d, "gamma/R4/sign_clamp", nonpos);
  return Dd<D>{op::IfThenElse(nonpos, zero, q.hi),
               op::IfThenElse(nonpos, zero, q.lo)};
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
//     e^{-z^2} G(1/z)/z (G from gamma's own dd fit, GammaTailG -- erfc's
//     plain-double G cost up to 2.1 ulp here, #64) and the power-of-two
//     scaling is applied
//     last, which is what keeps the subnormal band a*phi in [700, 760] at a
//     single rounding.
template <class D>
struct GammaTemmeOut {
  GammaVal<D> val;
  op::V<D> is_p;  // 1.0 where the directly computed side is P
  op::M<D> sat;
};

template <class D>
HWY_NOINLINE GammaTemmeOut<D> GammaTemme(D d, op::V<D> a, op::V<D> x) {
  const auto one = op::Set(d, 1.0);
  const auto zero = op::Zero(d);

  // Domain clamp (no-op on every live lane): keeps 1/a finite for the
  // discarded lanes this core also runs on, a subnormal a among them.
  a = op::Max(a, op::Set(d, detail::kGammaAT));

  const auto xma = op::Sub(x, a);  // exact by Sterbenz on this region
  const auto sgn = op::CopySign(one, xma);
  CORVUS_ROUTE(d, "gamma/R3/p_direct", op::Lt(xma, zero));
  CORVUS_ROUTE(d, "gamma/R3/q_direct", op::Ge(xma, zero));
  // TWIN: gammainv-inl.h repeats this prescale for its own forward (search
  // kGammaScaleAbove there).
  // Dekker-ceiling prescale (kGammaScaleAbove): u and phi are computed in
  // the scaled frame, a*phi scaled back with two exact multiplies. A
  // DISCARDED tiny-x lane's x*dn can underflow inexactly; such a lane fed
  // these sites NaN before the prescale existed, and its z is scrubbed at
  // the gather below either way.
  const auto big = op::Gt(a, op::Set(d, kGammaScaleAbove));
  CORVUS_ROUTE(d, "gamma/R3/prescale", big);
  const auto dn = op::IfThenElse(big, op::Set(d, kGammaScaleDown), one);
  const auto up = op::IfThenElse(big, op::Set(d, kGammaScaleUp), one);
  const auto as = op::Mul(a, dn);
  const auto xs = op::Mul(x, dn);
  const auto u = DdMul(d, TwoSum(d, xs, op::Neg(as)), DdRecip(d, as));
  const auto phi = Log1pmxDd(d, u);
  const auto aphi_s = DdMulD(d, phi, as);
  const auto aphi = Dd<D>{op::Mul(aphi_s.hi, up), op::Mul(aphi_s.lo, up)};

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

  // GammaClampE's NaN -> floor rule assumes every NaN lane is hundreds of
  // e-foldings past underflow. The diagonal x = a is the one unsaturated
  // point of the huge-a band, so it must never arrive as NaN -- the
  // prescale above is what guarantees it (#12). Debug-only premise check.
  HWY_DASSERT(op::AllFalse(
      d,
      IndMask(d, op::Mul(Ind(d, op::IsNaN(aphi.hi)), Ind(d, op::Eq(x, a))))));
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
  // FLOOR: brk_core*exd < 2^-968 on tail lanes (a*phi > 36, exd tiny), whose
  // core the tm select below discards; on core lanes (z <= 6) the sum is
  // ~1/2 erfc(z) >= 2^-57, so the <= 2^-1072 loss is invisible.
  const auto core = DdAdd(d, half_erfc, DdMul(d, brk_core, exd));

  // TWIN: beta-inl.h BetaR3Temme tail (cpsi > 36), a copy of this block;
  // it matches since #69 (same GammaTailG, dd G/(2z), guarded z.lo fold).
  // --- tail, a*phi > 36 ---------------------------------------------------
  // gt = G/(2z) = 1/2 erfc(z) e^{z^2} at z.hi, as a dd (#57: it was a
  // rounded double, and z.lo was dropped; #64: G itself is now a dd, and
  // GammaTailG's comment in src/erfc_tail_dd-inl.h covers the z clamp
  // below). z.lo enters through the
  // slope:
  // with f(z) = 1/2 erfc(z) e^{z^2}, f' = 2z f - 1/sqrt(pi), so
  // f(z.hi + z.lo) = gt - z.lo (1/sqrt(pi) - 2 z gt) + O(z.lo^2) -- the same
  // first-order fold the core applies to its own erfc term.
  //
  // ONLY ON UNSATURATED LANES (a*phi <= -kGammaExpFloor = 800, i.e.
  // z <= 28.3). The bracket 1/sqrt(pi) - 2z*gt is ~1/(2 sqrt(pi) z^2) but
  // is formed by a subtraction carrying ~2^-53 absolute, so the computed
  // correction's error relative to gt is ~2^-105 z^2: 2^-94 at z = 28.3,
  // and garbage for the huge z of a saturated lane (z ~ 1e50 at a ~ 1e172),
  // where it kept brk_tail from underflowing to zero. The forward forces
  // those lanes to an exact 0/1 anyway, but gammainv takes the LOG of this
  // dd on every R3 lane, and a nonzero garbage tail there returned a wrong
  // root (beyond-resolution bucket, 1.9e15 ULP).
  const auto zt = op::Min(op::Max(zs, op::Set(d, detail::kGammaTailZLo)),
                          op::Set(d, detail::kGammaTailZHi));
  const auto two_z = op::Add(zt, zt);
  const auto r2z = DdRecip(d, two_z);
  const auto gtd = DdMul(d, GammaTailG(d, r2z), r2z);
  const auto live_tail =
      op::Ge(op::Set(d, -detail::kGammaExpFloor), aphi.hi);
  const auto zlo_corr = op::IfThenElse(
      live_tail,
      op::Mul(z.lo, op::Sub(op::Set(d, detail::kGammaInvSqrtPiHi),
                            op::Mul(two_z, gtd.hi))),
      zero);
  const auto brk_tail = DdAdd(d, s_rv, DdAddD(d, gtd, op::Neg(zlo_corr)));
  const auto exf = ExpDdFrac(d, ea.hi, ea.lo);
  const auto tail = GammaScale(d, DdMul(d, exf.m, brk_tail), exf.e);

  const auto tm = op::Gt(aphi.hi, op::Set(d, kGammaTemmeZ2Split));
  CORVUS_ROUTE(d, "gamma/R3/core",
               op::Ge(op::Set(d, kGammaTemmeZ2Split), aphi.hi));
  CORVUS_ROUTE(d, "gamma/R3/tail", tm);
  const GammaVal<D> val{op::IfThenElse(tm, tail.v, DdToDouble(core)),
                        Dd<D>{op::IfThenElse(tm, tail.dd.hi, core.hi),
                              op::IfThenElse(tm, tail.dd.lo, core.lo)}};
  return {val, Ind(d, op::Lt(xma, zero)), ea.sat};
}

// ------------------------------------------------------------------------
// The driver. kP selects gamma_p (true) or gamma_q (false); the two differ
// only in who owns the R1/R4 overlap and in which side is complemented.
// HWY_NOINLINE like the region cores, and for the same MSVC-codegen
// reason: even with the cores outlined this function still inlines
// LgammaPosDd, Log1pmxDd, ExpDdFrac and the erfc core, and inlining it
// into the two export loops re-created a function big enough to stall
// cl.exe's back end past the CI timeout.
template <bool kP, class D>
HWY_NOINLINE op::V<D> GammaVec(D d, op::V<D> a_in, op::V<D> x_in) {
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

  // TWIN: gammainv-inl.h repeats this region map ("region map: GammaVec's");
  // a routing change here must be mirrored there.
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
  CORVUS_ROUTE(d, "gamma/R1", m_r1);
  CORVUS_ROUTE(d, "gamma/R2", m_r2);
  CORVUS_ROUTE(d, "gamma/R3", m_r3);
  CORVUS_ROUTE(d, "gamma/R4", m_r4);

  // --- assembled per-lane direct side -------------------------------------
  GammaVal<D> val{zero, Dd<D>{zero, zero}};
  auto is_p = zero;

  const bool need12 = !op::AllFalse(d, m_r1) || !op::AllFalse(d, m_r2);
  if (need12) {
    // E = a*log(x) - x - lgamma(a), the log of the common prefactor.
    // FLOOR: |a*log x| < 2^-968 needs a < 2^-915 (|log x| is 0 exactly or in
    // [2^-53, 745]); then lgamma(a) >= 634 dominates E, whose own dd rounding
    // (~2^-97) dwarfs the <= 2^-1072 loss, and e^E makes it relative error.
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
      CORVUS_ROUTE_SCOPE(d, m_r1);
      // (-) log a turns Gamma(a) into Gamma(a+1) without forming 1 + a and
      // without a 1/a that would be +inf for subnormal a. The two logs
      // cancel to ~0 for tiny a, but only their ABSOLUTE error matters here
      // (it is the argument of an exponential), and that is ~2^-89.
      const auto la = LogDdAny(d, a);
      const auto e1 =
          DdAdd(d, eb, Dd<D>{op::Neg(la.hi), op::Neg(la.lo)});
      const auto ea = GammaClampE(d, e1.hi, e1.lo);
      CORVUS_ROUTE(d, "gamma/R1/saturated", ea.sat);
      CORVUS_ROUTE(d, "gamma/R1/nan_exponent", op::IsNaN(e1.hi));
      const auto as = op::IfThenElse(ea.sat, safe_a, a);
      const auto xs = op::IfThenElse(ea.sat, safe_x, x);
      const auto ex = ExpDdFrac(d, ea.hi, ea.lo);
      auto r1 = GammaScale(d, DdMul(d, ex.m, GammaSeriesSum(d, as, xs)), ex.e);
      r1.v = op::IfThenElse(ea.sat, zero, r1.v);
      r1.dd.hi = op::IfThenElse(ea.sat, zero, r1.dd.hi);
      r1.dd.lo = op::IfThenElse(ea.sat, zero, r1.dd.lo);
      CORVUS_ROUTE(d, "gamma/R1/subnormal",
                   op::And(op::Lt(r1.v, op::Set(d, 0x1.0p-1022)),
                           op::MaskNot(ea.sat)));
      val.v = op::IfThenElse(m_r1, r1.v, val.v);
      val.dd.hi = op::IfThenElse(m_r1, r1.dd.hi, val.dd.hi);
      val.dd.lo = op::IfThenElse(m_r1, r1.dd.lo, val.dd.lo);
      is_p = op::IfThenElse(m_r1, one, is_p);
    }
    if (!op::AllFalse(d, m_r2)) {
      CORVUS_ROUTE_SCOPE(d, m_r2);
      const auto ea = GammaClampE(d, eb.hi, eb.lo);
      CORVUS_ROUTE(d, "gamma/R2/saturated", ea.sat);
      CORVUS_ROUTE(d, "gamma/R2/nan_exponent", op::IsNaN(eb.hi));
      const auto as = op::IfThenElse(ea.sat, safe_a, a);
      const auto xs = op::IfThenElse(ea.sat, safe_x, x);
      const auto ex = ExpDdFrac(d, ea.hi, ea.lo);
      auto r2 = GammaScale(d, DdMul(d, ex.m, GammaCfRecip(d, as, xs)), ex.e);
      r2.v = op::IfThenElse(ea.sat, zero, r2.v);
      r2.dd.hi = op::IfThenElse(ea.sat, zero, r2.dd.hi);
      r2.dd.lo = op::IfThenElse(ea.sat, zero, r2.dd.lo);
      CORVUS_ROUTE(d, "gamma/R2/subnormal",
                   op::And(op::Lt(r2.v, op::Set(d, 0x1.0p-1022)),
                           op::MaskNot(ea.sat)));
      val.v = op::IfThenElse(m_r2, r2.v, val.v);
      val.dd.hi = op::IfThenElse(m_r2, r2.dd.hi, val.dd.hi);
      val.dd.lo = op::IfThenElse(m_r2, r2.dd.lo, val.dd.lo);
    }
  }
  if (!op::AllFalse(d, m_r3)) {
    CORVUS_ROUTE_SCOPE(d, m_r3);
    auto t = GammaTemme(d, a, x);
    CORVUS_ROUTE(d, "gamma/R3/saturated", t.sat);
    t.val.v = op::IfThenElse(t.sat, zero, t.val.v);
    t.val.dd.hi = op::IfThenElse(t.sat, zero, t.val.dd.hi);
    t.val.dd.lo = op::IfThenElse(t.sat, zero, t.val.dd.lo);
    CORVUS_ROUTE(d, "gamma/R3/subnormal",
                 op::And(op::Lt(t.val.v, op::Set(d, 0x1.0p-1022)),
                         op::MaskNot(t.sat)));
    val.v = op::IfThenElse(m_r3, t.val.v, val.v);
    val.dd.hi = op::IfThenElse(m_r3, t.val.dd.hi, val.dd.hi);
    val.dd.lo = op::IfThenElse(m_r3, t.val.dd.lo, val.dd.lo);
    is_p = op::IfThenElse(m_r3, t.is_p, is_p);
  }
  if (!op::AllFalse(d, m_r4)) {
    CORVUS_ROUTE_SCOPE(d, m_r4);
    const auto q4 = GammaSmallQ(d, a, x);
    CORVUS_ROUTE(d, "gamma/R4/subnormal",
                 op::Lt(DdToDouble(q4), op::Set(d, 0x1.0p-1022)));
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
  CORVUS_ROUTE(d, "gamma/special/x_zero", op::Eq(x_in, zero));
  res = op::IfThenElse(op::Eq(x_in, zero), p_at0, res);
  CORVUS_ROUTE(d, "gamma/special/a_zero", op::Eq(a_in, zero));
  res = op::IfThenElse(op::Eq(a_in, zero), p_at, res);
  CORVUS_ROUTE(d, "gamma/special/x_inf", op::Eq(x_in, inf));
  res = op::IfThenElse(op::Eq(x_in, inf), p_at, res);
  CORVUS_ROUTE(d, "gamma/special/a_inf", op::Eq(a_in, inf));
  res = op::IfThenElse(op::Eq(a_in, inf), p_at0, res);
  CORVUS_ROUTE(d, "gamma/special/sum_zero", op::Eq(op::Add(a_in, x_in), zero));
  res = op::IfThenElse(op::Eq(op::Add(a_in, x_in), zero), qnan, res);
  CORVUS_ROUTE(d, "gamma/special/both_inf", op::Eq(op::Min(a_in, x_in), inf));
  res = op::IfThenElse(op::Eq(op::Min(a_in, x_in), inf), qnan, res);
  CORVUS_ROUTE(d, "gamma/special/a_neg", op::Lt(a_in, zero));
  res = op::IfThenElse(op::Lt(a_in, zero), qnan, res);
  CORVUS_ROUTE(d, "gamma/special/x_neg", op::Lt(x_in, zero));
  res = op::IfThenElse(op::Lt(x_in, zero), qnan, res);
  CORVUS_ROUTE(d, "gamma/special/a_nan", op::IsNaN(a_in));
  res = op::IfThenElse(op::IsNaN(a_in), a_in, res);  // payload preserved
  CORVUS_ROUTE(d, "gamma/special/x_nan", op::IsNaN(x_in));
  res = op::IfThenElse(op::IsNaN(x_in), x_in, res);
  return res;
}

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#endif  // include guard
