// Regularized incomplete beta I_x(a,b) -- beta_p and beta_q.
// Per-target include guard (Highway -inl.h idiom).
//
// Both public functions share one kernel: four region cores, two prefactor
// paths, and a router. Everything here follows PLAN.md, "Regularized
// incomplete beta -- detail design" together with its G1a, G1b and G1c
// corrections (the routing ORDER, the R2 orientation rule, the R3 ratio-band
// membership and the R4 box are all probe-validated there and are reproduced,
// not re-derived); the tables and every threshold live in src/beta_data.h,
// emitted and self-checked by tools/gen_beta_data.py.
//
// THE RULE THE WHOLE DESIGN TURNS ON is gamma's: compute ONE side of the pair
// directly, in double-double, and get the other as 1 (-) it rounded once. The
// difference from gamma is that WHICH side is computed is decided by the
// REGION, not by a single predicate: the router picks an ORIENTATION -- the
// argument triple (alpha, beta, xi) actually handed to a core -- and the
// direct/complement handout is a final select. The G1a escalation was exactly
// the failure to separate those two things.
//
// ORIENTATION AND ROUTING (route_final() in tools/gen_beta_data.py, verbatim;
// c = a+b, nu = ab/c, p = a/c, q = b/c, y = 1-x):
//   0. min(a,b) <= eps_R4, evaluated TINY-FIRST as (tau, B, xi_tau):
//      if tau*|ln xi_tau| <= ln2 AND xi_tau <= xi1 AND B*xi_tau <= B1 -> R4.
//      Otherwise fall through (this hoist above R1 is the G1b correction: R1's
//      box has no alpha floor and would otherwise steal (1e-20, 1, 0.4) and
//      evaluate a flat 1.0).
//   1. R1 if EITHER orientation has xi <= xi1 AND beta*xi <= B1.
//   2. R3 if nu >= T_ridge AND x/p in [1/2, 2] AND y/q in [1/2, 2] -- the
//      ridge RATIO band (the G1c correction; the earlier "cpsi <= 800 strip"
//      implied a zeta domain no 32 KiB tensor can fit at dd level).
//      Orientation by the mean predicate, i.e. x/p <= 1.
//   3. Else R2, orientation by the pinned rule x < (a+1)/(c+2).
// The evaluated side is <= 1 - 2^-12 for R1/R2/R3 (generator self-check (e)),
// which is all the complement needs; R4 satisfies the same requirement by
// CONSTRUCTION rather than by routing -- it assembles the small side
// analytically and never complements a rounded near-one value.
//
// REGIONS (canonical direct side (alpha, beta, xi) after routing)
//   R1  power series, BPSER analog. I = xi^alpha/B(alpha,beta) * S,
//       S = sum_{n>=0} t_n/(alpha+n), t_n = t_{n-1}(n-beta)xi/n. Written as
//       exp(E1)*(1 + alpha*sum_{n>=1}) so no 1/alpha is ever formed.
//   R2  backward continued fraction, DLMF 8.17.22, fixed depth kBetaN2 with
//       no convergence test (gamma's GammaCfRecip lesson: a stopping test
//       false-converges, a fixed depth cannot be fooled).
//   R3  Temme erf-form on the ridge; gamma's GammaTemme with cpsi replacing
//       a*phi and a two-dimensional (zeta, p) correction table.
//   R4  tiny-min analytic assembly, the APSER form:
//       Qtilde = -expm1(w + log1p(tau*Sigma)) -- the SMALL side directly, in
//       log space, so no cancellation of near-one quantities ever happens.
//
// PREFACTOR. E = alpha*ln xi + beta*ln y - ln B(alpha,beta), with y carried as
// the EXACT dd 1 (-) xi (a rounded 1-xi costs beta*2^-53 absolute in E, fatal
// at large beta), and never as three independently rounded lgammas at large c.
// Two paths (see BetaPrefactor):
//   PB  c > C_lg AND min >= Z0: the analytic Stirling difference
//       E = 1/2 ln(nu/2pi) (-) cpsi (-) Delta, whose derivation removes the
//       cancellation on paper.
//   PA  otherwise: -ln B = LgammaDiffDd(max, min) (-) LgammaPosDd(min), the
//       analytic lgamma difference. This is the design's P3 form used for P1's
//       range as well -- see the note on LgammaDiffDd for why, and NOTE that
//       it means the rounded c is NEVER an lgamma argument anywhere in this
//       kernel, which is a stronger statement than the design's own P1
//       (LgammaPosDd(c.hi) plus a c.lo*psi(c.hi) correction) makes.
//
// SATURATION. e^E underflows below about -745; every region that ends in an
// exponential clamps at kBetaExpFloor = -800 and reports the clamp
// (BetaClampE). Saturated lanes are forced to an exact 0/1 pair AND their
// (alpha, beta, xi) is scrubbed to a benign interior point before the series
// and the continued fraction run -- without which the CF's (c+m) term and the
// series' terms would carry infinities from perfectly ordinary (if utterly
// saturated) arguments.
//
// ACCURACY. Per-tier gates pinned at G4 (2026-08-05/06) against the
// harness-certified reference set; docs/ACCURACY.md carries the audited
// per-region table.
#if defined(CORVUS_BETA_INL_H_) == defined(HWY_TARGET_TOGGLE)
#ifdef CORVUS_BETA_INL_H_
#undef CORVUS_BETA_INL_H_
#else
#define CORVUS_BETA_INL_H_
#endif

#include <limits>

#include "src/beta_data.h"
#include "src/dd-inl.h"
#include "src/dd_special-inl.h"
#include "src/erfc_core-inl.h"
#include "src/exp_dd-inl.h"
// gamma-inl.h: ONLY for the (C) gamma-limit slice's two template cores
// (GammaSeriesSum, GammaCfRecip) -- templates instantiate what is called,
// so this does not pull gamma's other cores into beta.o. The slice exists
// because the beta CF is structurally degenerate above kBetaGammaLim
// (PLAN.md escalation (C)); its dd argument-sensitivity lives entirely in
// e^-t, which beta's own dd prefactor machinery absorbs.
#include "src/gamma-inl.h"
#include "src/lgamma-inl.h"
#include "src/log_dd-inl.h"
#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

namespace op = ops;

// --- kernel-internal constants -------------------------------------------
// Structural (crossovers and loop shapes chosen from the error budget at the
// use sites), not fitted table data, so they live here rather than in the
// generated src/beta_data.h -- the split erfinv and gamma both make.

// Per-lane freeze threshold for the two summed series: a term below this times
// the running sum cannot move the dd result.
//
// 2^-105, NOT gamma's 2^-60. gamma's series is a factor of its answer, so a
// 2^-60 truncation of the sum is a 2^-60 relative error in the result. R4's is
// not: its Sigma reaches the answer through log1p(tau*Sigma) which then
// CANCELS against w to ~2^-15 of itself, so a 2^-60 truncation of Sigma
// arrives amplified by |tau*Sigma|/|y|. Measured at (a, b, x) =
// (70.17, 7.4438e-6, 0.9), where the amplification is 3e4: the answer came
// back 16 ULP off, and the freeze was the whole of it. 2^-105 is the dd
// representation's own resolution, so the freeze now only ever skips terms the
// accumulator could not have held anyway; the cost is that most lanes run the
// full fixed cap (64 for R1, 48 for R4) instead of stopping early.
inline constexpr double kBetaFreezeEps = 0x1.0p-105;

// R3 core/tail split, exact: cpsi = 36 <=> z = sqrt(cpsi) = 6, the same
// threshold erfc's own core/tail split uses, so the two halves of this region
// consume exactly the two halves of erfc's machinery.
inline constexpr double kBetaZ2Split = 36.0;

// Ceiling on nu when forming 1/sqrt(2*pi*nu): 2*pi*nu overflows above ~2.9e307
// and DdSqrt would then return NaN. Above this clamp the whole S/sqrt(2*pi*nu)
// term is below 2^-450 against a result of order 1/2, so the clamped value and
// the true one round identically (kGammaTwoPiAClamp, same argument). The clamp
// must ALSO sit under ops::ProdLow's 2^996 non-FMA Dekker ceiling minus the
// 2^27 split factor: the original 2^1000 value overflowed the split of nu
// inside DdMulD(twopi, nu) on SSE tiers -- NaN on the a = b >= 2^998 diagonal,
// caught by the G4 capped sweep (FMA targets never see it).
inline constexpr double kBetaTwoPiNuClamp = 0x1.0p+900;

// R4 subnormal-tau reframing [ELEVENTH correction; see BetaR4Tiny]. All
// five are exact powers of two: the two thresholds partition the hazard
// (tau below kBetaTauSubn puts tau-scale products at the subnormal
// boundary, where non-FMA Dekker residuals collapse), kBetaTauUp/Down are
// the linear-lane reframing pair, and kBetaTinyPairUp is the joint scale
// for the both-tiny closed form's division.
inline constexpr double kBetaTauSubn = 0x1.0p-950;
inline constexpr double kBetaBTinyCut = 0x1.0p-140;
inline constexpr double kBetaTauUp = 0x1.0p+700;
inline constexpr double kBetaTauDown = 0x1.0p-700;
inline constexpr double kBetaTinyPairUp = 0x1.0p+900;

// Exact power-of-two prescale for the cpsi machinery. See BetaPsiCore: it is
// what makes c = a+b representable, and what keeps every Dekker-split operand
// under ops::ProdLow's 2^996 non-FMA ceiling (AGENTS.md).
inline constexpr double kBetaScaleAbove = 0x1.0p+900;
inline constexpr double kBetaScaleDown = 0x1.0p-200;
inline constexpr double kBetaScaleUp = 0x1.0p+200;

// Ceiling applied to z before 1/z inside BinetDiffDd, purely to keep the
// Dekker split of DdRecip's argument in range; Binet's difference is already
// below 2^-1000 there, so the clamp cannot move any budgeted quantity.
inline constexpr double kBetaBinetZClamp = 0x1.0p+500;

// Steps of the Gamma(z+1) = z*Gamma(z) up-walk in LgammaDiffDd. Ten unit steps
// carry any z > 0 into [Z0, Z0+1) with Z0 = 10.
inline constexpr int kBetaWalkSteps = 10;

// 2*pi, 1/sqrt(pi) and ln(2*pi) as dd pairs. Same bit patterns as
// src/gamma_data.h's kGammaTwoPi*/kGammaInvSqrtPi*; repeated here rather than
// included because beta deliberately does not depend on gamma's TU.
inline constexpr double kBetaTwoPiHi = 0x1.921fb54442d18p+2;
inline constexpr double kBetaTwoPiLo = 0x1.1a62633145c07p-52;
inline constexpr double kBetaInvSqrtPiHi = 0x1.20dd750429b6dp-1;
inline constexpr double kBetaLnTwoPiHi = 0x1.d67f1c864beb5p+0;
inline constexpr double kBetaLnTwoPiLo = -0x1.65b5a1b7ff5dfp-54;

// Benign interior point every core is scrubbed to on its inactive lanes.
// R4's box does not contain it (its first parameter must be <= eps_R4), so
// that core scrubs to its own interior point -- see BetaR4Tiny's call site.
inline constexpr double kBetaSafeA = 2.0;
inline constexpr double kBetaSafeB = 3.0;
inline constexpr double kBetaSafeX = 0.25;

// --- small dd conveniences ------------------------------------------------
template <class D>
HWY_INLINE Dd<D> DdNeg(D, Dd<D> a) {
  return Dd<D>{op::Neg(a.hi), op::Neg(a.lo)};
}
template <class D>
HWY_INLINE Dd<D> DdSub(D d, Dd<D> a, Dd<D> b) {
  return DdAdd(d, a, DdNeg(d, b));
}
// a/b for dd operands, WITHOUT forming 1/b: b.hi may be subnormal (the walk in
// LgammaDiffDd divides by a subnormal max-parameter) where DdRecipDd's 1/b.hi
// would be infinite.
//
// A huge divisor is taken down by an exact power of two first -- the quotient
// is unchanged, and it keeps b.hi under ops::ProdLow's 2^996 non-FMA Dekker
// ceiling (AGENTS.md). Numerators that underflow under that scaling had
// quotients below 2^-1722 and were going to be zero anyway. Precondition:
// |a/b| stays inside the Dekker range too, which every call site here
// satisfies by construction (all quotients are O(1) or smaller except R2's
// bounded (beta-m)xi/(alpha+2m), which saturation caps near 800).
template <class D>
HWY_INLINE Dd<D> DdDivDd(D d, Dd<D> a, Dd<D> b) {
  const auto one = op::Set(d, 1.0);
  const auto s = op::IfThenElse(op::Gt(op::Abs(b.hi),
                                       op::Set(d, kBetaScaleAbove)),
                                op::Set(d, kBetaScaleDown), one);
  const Dd<D> an{op::Mul(a.hi, s), op::Mul(a.lo, s)};
  const Dd<D> bn{op::Mul(b.hi, s), op::Mul(b.lo, s)};
  const auto q = op::Div(an.hi, bn.hi);
  const auto rem = DdSub(d, an, DdMulD(d, bn, q));
  return Fast2Sum(d, q, op::Div(rem.hi, bn.hi));
}
// log(1+w) for a dd w > -1, to dd precision relative to the result. Written on
// Log1pmxDd because that is the primitive that is accurate relative to
// phi = w - log1p(w); w (-) phi then has no cancellation at either end (phi is
// O(w^2) for small w, and both are positive and well separated for w ~ 1).
template <class D>
HWY_INLINE Dd<D> Log1pDdWide(D d, Dd<D> w) {
  return DdSub(d, w, Log1pmxDd(d, w));
}

// --- indicator helpers ----------------------------------------------------
// The ops facade deliberately exposes no mask AND/OR, and the region map is a
// boolean expression over a dozen comparisons. Carrying predicates as 1.0/0.0
// vectors makes AND a multiply, OR a max and NOT a subtract from one, all
// exact, and converts back to a mask with a single compare. (Same device as
// gamma-inl.h; distinct names so a TU may include both headers.)
template <class D, class M>
HWY_INLINE op::V<D> BetaInd(D d, M m) {
  return op::IfThenElse(m, op::Set(d, 1.0), op::Zero(d));
}
template <class D>
HWY_INLINE op::V<D> BetaIndNot(D d, op::V<D> v) {
  return op::Sub(op::Set(d, 1.0), v);
}
template <class D>
HWY_INLINE op::M<D> BetaIndMask(D d, op::V<D> v) {
  return op::Gt(v, op::Set(d, 0.5));
}

// A region's contribution: the direct-side value rounded exactly once, and the
// same quantity as a dd so the other side can be formed as 1 (-) it before any
// rounding. For the regions that end in a power-of-two scaling the two are NOT
// the same computation -- see BetaScale.
template <class D>
struct BetaVal {
  op::V<D> v;
  Dd<D> dd;
};

// Apply exp_dd's power-of-two scaling LAST (the erfc/gamma pattern), twice
// over: v is the direct answer with ONE rounding even in the subnormal band,
// dd is the scaled pair for the complement (which is only ever taken where the
// result is >= ~0.4 and the scaling is therefore exact).
template <class D>
HWY_INLINE BetaVal<D> BetaScale(D d, Dd<D> m, op::V<op::SignedTag<D>> e) {
  return {ScaleTwo(d, DdToDouble(m), e),
          Dd<D>{ScaleTwo(d, m.hi, e), ScaleTwo(d, m.lo, e)}};
}

// Exponential argument after the underflow guard, plus the mask saying it
// fired. gamma's GammaClampE with kBetaExpFloor.
template <class D>
struct BetaExpArg {
  op::V<D> hi;
  op::V<D> lo;
  op::M<D> sat;
};

template <class D>
HWY_INLINE BetaExpArg<D> BetaClampE(D d, op::V<D> eh, op::V<D> el) {
  const auto zero = op::Zero(d);
  const auto floorv = op::Set(d, detail::kBetaExpFloor);
  // alpha*log(xi) overflows to -inf for huge alpha and its dd residual to NaN
  // with it; so does an infinite cpsi. Every such lane is many hundreds of
  // e-foldings past the underflow point, so mapping NaN onto the floor selects
  // the answer (an exact 0) the finite arithmetic was going to give.
  const auto bad = op::IsNaN(eh);
  eh = op::IfThenElse(bad, floorv, eh);
  el = op::IfThenElse(bad, zero, el);
  const auto sat = op::Ge(floorv, eh);
  // Clamping hi without clearing lo would hand ExpDdFrac an unnormalized pair,
  // which its argument reduction is not contracted to accept.
  return {op::Max(eh, floorv), op::IfThenElse(sat, zero, el), sat};
}

// ------------------------------------------------------------------------
// OUTLINED log and exp [MSVC BUILD-TIME GATE, AGENTS.md]. These are thin
// wrappers whose only purpose is the HWY_NOINLINE: log_dd (via LogDdAny) and
// exp_dd (via ExpDd) are each large, and this file reaches them from about
// eleven call sites -- inside LgammaDiffDd, the driver's routing and
// prefactor assembly, and the gamma-limit slice. Inlined, each of those
// becomes its own copy of a table gather plus a polynomial IN EVERY ONE OF
// THE COMPILED TARGETS, and cl.exe's optimizer is superlinear in function
// size: this is the same fix that took betainv.cpp (src/betainv-inl.h) from
// past 45 minutes and 7 GB on one MSVC invocation to 127 s, retro-applied
// here per PLAN.md's "MSVC build-time headroom" item. Bit-identity is
// guaranteed by contraction-off and verified by byte-comparing the ULP
// tables across the change.
template <class D>
HWY_NOINLINE Dd<D> BetaLog(D d, Dd<D> x) {
  return LogDdAny(d, x);
}
template <class D>
HWY_NOINLINE Dd<D> BetaLog(D d, op::V<D> x) {
  return BetaLog(d, Dd<D>{x, op::Zero(d)});
}
template <class D>
HWY_NOINLINE Dd<D> BetaExpDd(D d, op::V<D> xh, op::V<D> xl) {
  return ExpDd(d, xh, xl);
}

// ------------------------------------------------------------------------
// Binet's function phi(z) = lgamma(z) - [(z-1/2)ln z - z + 1/2 ln 2pi], the
// Stirling tail, as sum_k B_2k/(2k(2k-1) z^(2k-1)) for z >= kBetaZ0 = 10.
//
// HOME OF THIS FUNCTION [G3 decision]. The design left it open whether to
// expose lgamma-inl.h's Stirling tail or write a fresh one. Fresh, here:
// lgamma's tail is the same series but truncated for ITS OWN X0 = 8 into a
// nine-entry table, and is welded into LgammaStirling's 2^-200 scaling block;
// beta_data.h already carries kBetaBinetCoef[16] emitted and self-checked
// (check (g), 2^-74.6) at Z0 = 10 precisely so that this file needs nothing
// from lgamma's internals. src/lgamma-inl.h is therefore UNTOUCHED and the
// dd_special-hoist byte-identity protocol does not apply.
//
// Plain double throughout: |phi| <= 1/(12*10) = 2^-6.9 here, so a half-ulp of
// it is 2^-60 absolute, and phi enters E additively against a 2^-56 budget.
// z = +inf gives r = 0 and an exact 0, which is the right limit and is a real
// input (c = a+b overflows for a, b both near DBL_MAX).
template <class D>
HWY_INLINE op::V<D> BinetVal(D d, op::V<D> z) {
  const auto r = op::Div(op::Set(d, 1.0), z);
  const auto w = op::Mul(r, r);
  auto p = op::Set(d, detail::kBetaBinetCoef[detail::kBetaKB - 1]);
  for (int k = detail::kBetaKB - 2; k >= 0; --k) {
    p = op::MulAdd(p, w, op::Set(d, detail::kBetaBinetCoef[k]));
  }
  return op::Mul(p, r);
}

// phi(z*(1+w)) - phi(z) for z >= kBetaZ0 and 0 <= w <= 1, to dd precision
// RELATIVE to the difference.
//
// NOT computed as BinetVal(z+m) - BinetVal(z). That difference is
// ~w/(12z^2) against two values of ~1/(12z), i.e. it cancels by a factor w --
// and w is m/z with m free to be subnormal, so the cancellation is unbounded
// while the accuracy demanded of the result is not. (The difference feeds
// LgammaDiffDd, where it carries ~2^-11 of a quantity that must be relatively
// accurate; a naive difference is already short of budget at m ~ 1e-12.)
//
// The exact factorization instead: with rho = z/(z(1+w)) = 1/(1+w),
//   phi(z(1+w)) - phi(z) = sum_k c_k z^-(2k-1) (rho^(2k-1) - 1)
//                        = (rho-1) * sum_k c_k z^-(2k-1) G_k,
//   G_k = sum_{j=0}^{2k-2} rho^j,  G_1 = 1, G_{k+1} = G_k + rho^(2k-1)(1+rho).
// rho - 1 = -w/(1+w) is formed in dd and is relatively accurate however small
// w is; the G_k recurrence and the coefficient sum have no cancellation at all
// (every term has one sign), so plain double suffices for them -- except the
// k = 1 term, which is ~1000x the rest at z = 10 and is carried as a dd
// (c_1 = 1/12 exactly, so its dd pair is kBetaRecipN[11]).
template <class D>
HWY_NOINLINE Dd<D> BinetDiffDd(D d, Dd<D> z, Dd<D> w) {
  const auto one = op::Set(d, 1.0);
  const auto g = DdMul(d, DdNeg(d, w), DdRecipDd(d, DdAddD(d, w, one)));
  const auto rho = op::Add(one, g.hi);
  const auto rho2 = op::Mul(rho, rho);

  const auto zc = op::Min(z.hi, op::Set(d, kBetaBinetZClamp));
  const auto rd = DdRecip(d, zc);
  const auto r = rd.hi;
  const auto r2 = op::Mul(r, r);

  const Dd<D> c1{op::Set(d, detail::kBetaRecipNHi[11]),
                 op::Set(d, detail::kBetaRecipNLo[11])};
  const auto lead = DdMul(d, rd, c1);

  auto gk = one;   // G_k
  auto pw = rho;   // rho^(2k-1)
  auto rp = r;     // r^(2k-1)
  auto rest = op::Zero(d);
  for (int k = 1; k < detail::kBetaKB; ++k) {
    gk = op::Add(gk, op::MulAdd(pw, rho, pw));
    pw = op::Mul(pw, rho2);
    rp = op::Mul(rp, r2);
    rest = op::MulAdd(op::Set(d, detail::kBetaBinetCoef[k]), op::Mul(rp, gk),
                      rest);
  }
  return DdMul(d, g, DdAddD(d, lead, rest));
}

// ------------------------------------------------------------------------
// lgamma(M + m) - lgamma(M) for 0 < m <= M, as a dd, RELATIVE to itself.
//
// This is the design's LgammaDiffDd, and it is what lets the prefactor avoid
// ever handing a rounded c = a+b to an lgamma: -ln B(alpha,beta) is
// lgamma(c) - lgamma(alpha) - lgamma(beta) = LgammaDiffDd(M, m) - lgamma(m)
// with M = max, m = min. The design reserved that spelling for P3 (c > C_lg,
// small first parameter) and gave P1 three LgammaPosDd calls plus a
// c.lo*DigammaRough(c.hi) correction; this kernel uses the difference form for
// P1's range too, for two reasons:
//   * it removes the rounded lgamma argument outright rather than correcting
//     for it, which is the design's own exact-c doctrine taken one step
//     further; and
//   * beta_data.h's DigammaRough table is fitted and self-checked (check (i))
//     only on (0, 2*Z0] = (0, 20], while P1's c reaches C_lg = 256, so the
//     shipped table cannot serve the correction it was emitted for.
// kBetaDigammaCoef/kBetaDigammaDeg are consequently unused by this kernel.
//
// TWO PIECES.
//  (a) An up-walk. For M < Z0 the Stirling tail is out of its table's range,
//      so Gamma(z+1) = z*Gamma(z) is used k <= 10 times:
//          lgamma(M+m) - lgamma(M)
//            = [lgamma(M+k+m) - lgamma(M+k)] - sum_{j<k} log1p(m/(M+j)).
//      The sum is NOT k separate log1p calls. Every factor (1 + m/(M+j)) is
//      >= 1, so their product is accumulated as 1 + s with
//      s <- s + eps + s*eps -- no cancellation anywhere, s relatively accurate
//      to ~2^-105 -- and the whole sum is one log1p(s) at the end.
//      M + j is carried as the EXACT TwoSum pair: a rounded M+j would put its
//      own 2^-53 into every walk term, and for tiny M the j = 0 term IS the
//      answer.
//  (b) The Stirling difference at z = M+k >= Z0, w = m/z:
//          m*ln z - z*phi(w) + (m-1/2)*log1p(w) + dBinet,
//      phi = Log1pmxDd, log1p(w) = w (-) phi(w). THE GROUPING IS THE POINT.
//      The textbook form is m*ln z + (z+m-1/2)*log1p(w) - m, whose leading
//      z*log1p(w) cancels the -m to many digits; using phi removes that
//      cancellation on paper, since z*log1p(w) - m == -z*phi(w) identically
//      when w = m/z. It also makes the result far less sensitive to w's own
//      rounding: d/dw of z*log1p(w) is z, but of z*phi(w) is only z*w/(1+w)
//      ~ m. That matters because w = m/z can land in or near the SUBNORMAL
//      range (m = 1e-300 over z = 1e8 gives 1e-308), where a dd cannot hold a
//      low word and w carries only 53 bits. Measured: at
//      (a, b, x) = (1e-300, 1e8, 8e-8) an explicitly compensated
//      (z*w - m) - z*phi(w) form -- algebraically identical, and what an
//      earlier version of this function used -- put the whole z*dw through and
//      landed 2.8e4 ULP out; the phi form is exact to the last bit there.
//      m - 1/2 is an exact TwoSum: m reaches 256 in P1's range and
//      (m-1/2)*log1p(w) is then ~177, whose last bits are budgeted.
//
// SCALING. The only operand that can leave ops::ProdLow's non-FMA Dekker range
// is z itself (M is unbounded above). w and the two products that carry a z
// are computed on an exactly down-scaled pair and scaled back; m is bounded by
// 256 at every live call site (PA has c <= C_lg or min < Z0; R4 hands
// tau <= kBetaPrTauMax = 2.5 for its lgdiff term and m <= 3/2 for its
// lgamma(1+tau) form [NINTH correction]), so no other operand is at risk.
// That audit is the AGENTS.md 2^996 rule discharged for this function.
template <class D>
HWY_NOINLINE Dd<D> LgammaDiffDd(D d, op::V<D> bigm, op::V<D> smallm) {
  const auto one = op::Set(d, 1.0);
  const auto zero = op::Zero(d);
  const auto z0 = op::Set(d, detail::kBetaZ0);

  // --- (a) up-walk ---
  Dd<D> s{zero, zero};
  auto kf = zero;
  for (int j = 0; j < kBetaWalkSteps; ++j) {
    const auto mj = TwoSum(d, bigm, op::Set(d, static_cast<double>(j)));
    // The fire masks shrink monotonically in j, so once no lane fires every
    // remaining step is a no-op and breaking out returns the same value.
    const auto fire = op::Gt(z0, mj.hi);
    if (op::AllFalse(d, fire)) break;
    const auto eps = DdDivDd(d, Dd<D>{smallm, zero}, mj);
    const auto sn = DdAdd(d, s, DdAdd(d, eps, DdMul(d, s, eps)));
    s = Dd<D>{op::IfThenElse(fire, sn.hi, s.hi),
              op::IfThenElse(fire, sn.lo, s.lo)};
    kf = op::Add(kf, op::IfThenElse(fire, one, zero));
  }
  const auto walk = Log1pDdWide(d, s);  // exactly 0 when no step fired

  // --- (b) Stirling difference at z = M + k ---
  const auto z = TwoSum(d, bigm, kf);  // exact
  const auto big = op::Gt(z.hi, op::Set(d, kBetaScaleAbove));
  const auto dn = op::IfThenElse(big, op::Set(d, kBetaScaleDown), one);
  const auto up = op::IfThenElse(big, op::Set(d, kBetaScaleUp), one);
  const Dd<D> zs{op::Mul(z.hi, dn), op::Mul(z.lo, dn)};
  const auto ms = op::Mul(smallm, dn);

  const auto w = DdDivDd(d, Dd<D>{ms, zero}, zs);  // = m/z, scale-invariant
  const auto phi = Log1pmxDd(d, w);
  const auto l1p = DdSub(d, w, phi);

  const auto lz = BetaLog(d, z);
  auto out = DdMulD(d, lz, smallm);  // m*ln z

  // z*phi(w), formed on the down-scaled z and brought back exactly.
  const auto zp = DdMul(d, zs, phi);
  const Dd<D> zp_u{op::Mul(zp.hi, up), op::Mul(zp.lo, up)};
  out = DdSub(d, out, zp_u);

  out = DdAdd(d, out, DdMul(d, TwoSum(d, smallm, op::Set(d, -0.5)), l1p));
  out = DdAdd(d, out, BinetDiffDd(d, z, w));
  return DdSub(d, out, walk);
}

// ------------------------------------------------------------------------
// The lambda / u / v / cpsi / nu machinery, shared by R3 and by the PB
// prefactor exactly as the design requires (ONE implementation).
//
//   lambda = alpha - c*xi = alpha*y - beta*xi,   u = -lambda/alpha,
//   v = +lambda/beta,      cpsi = alpha*phi(u) + beta*phi(v) >= 0,
//   nu = alpha*beta/c,     zeta = sign(lambda)*sqrt(cpsi/nu).
// alpha*u + beta*v = 0 identically, which is the cancellation the PB
// derivation removes on paper.
//
// LAMBDA IS THE DELICATE PART, and it is NOT alpha (-) c*xi. Three reasons:
//   * c*xi overflows exactly where c does;
//   * alpha*y and beta*xi are each Theta(nu) -- alpha*y = nu(1+v),
//     beta*xi = nu(1+u) -- not Theta(c), so the cancellation is measured
//     against nu and not against the much larger c; and
//   * both products are EXACT as TwoProds, so taking the leading difference
//     with an exact TwoSum first and folding the four residual words in
//     afterwards leaves an absolute error ~nu*2^-157 instead of ~nu*2^-105.
//     That matters near the ridge: an absolute delta in lambda moves the
//     result by delta/sqrt(2*nu*pi), so the coarser assembly would start
//     costing half-ulps at nu ~ 1e32, which is inside the domain. (Beyond
//     nu ~ 1e35 the only unsaturated lanes left are lambda exactly zero --
//     one ulp of xi moves lambda by c*2^-53, far past the saturation
//     threshold -- and an exactly zero lambda stays exactly zero here, which
//     is what makes I_{1/2}(a,a) reproduce 1/2.)
//
// EXACT PRESCALE. c = a+b is not representable for a, b both near DBL_MAX, and
// alpha itself can exceed ops::ProdLow's 2^996 non-FMA Dekker ceiling. Both
// are fixed by one exact power of two: u, v, p, q and zeta are 0-homogeneous
// in (alpha, beta) so they are unchanged, while cpsi and nu are 1-homogeneous
// and come back with a single exact multiply. This is the design's c-overflow
// guard, sited here [G3 decision]. It is safe because BOTH consumers of this
// function have min(alpha,beta) >= Z0 = 10 (R3 needs nu >= T_ridge = 32 <= min,
// PB needs min >= Z0), so the down-scaled small parameter cannot underflow.
template <class D>
struct BetaPsi {
  Dd<D> cpsi;  // unscaled
  Dd<D> nu;    // unscaled
  op::V<D> zeta;
  op::V<D> p;
  op::V<D> q;
  op::V<D> lam;  // scaled-frame lambda hi word; only its SIGN is consumed
};

template <class D>
HWY_NOINLINE BetaPsi<D> BetaPsiCore(D d, op::V<D> a, op::V<D> b, Dd<D> xi,
                                    Dd<D> y) {
  const auto one = op::Set(d, 1.0);

  const auto big = op::Gt(op::Max(a, b), op::Set(d, kBetaScaleAbove));
  const auto dn = op::IfThenElse(big, op::Set(d, kBetaScaleDown), one);
  const auto up = op::IfThenElse(big, op::Set(d, kBetaScaleUp), one);
  const auto as = op::Mul(a, dn);
  const auto bs = op::Mul(b, dn);
  const auto c = TwoSum(d, as, bs);  // exact and finite by construction

  const auto a0 = TwoProd(d, as, y.hi);
  const auto a1 = TwoProd(d, as, y.lo);
  const auto b0 = TwoProd(d, bs, xi.hi);
  const auto b1 = TwoProd(d, bs, xi.lo);
  const auto lead = TwoSum(d, a0.hi, op::Neg(b0.hi));  // exact
  Dd<D> lam{lead.hi, lead.lo};
  lam = DdAddD(d, lam, a0.lo);
  lam = DdAddD(d, lam, op::Neg(b0.lo));
  lam = DdAddD(d, lam, a1.hi);
  lam = DdAddD(d, lam, op::Neg(b1.hi));
  lam = DdAddD(d, lam, a1.lo);
  lam = DdAddD(d, lam, op::Neg(b1.lo));

  const auto u = DdMul(d, DdNeg(d, lam), DdRecip(d, as));
  const auto v = DdMul(d, lam, DdRecip(d, bs));

  // 1 + u AND 1 + v BY CLOSED FORM ON THE CORNER LANES [u -> -1 fix,
  // 2026-08-12; the shipped-since-v0.1.0 beta_p/beta_q defect pair]. The
  // identities
  //     1 + u = (alpha - lambda)/alpha = c*xi/alpha,
  //     1 + v = (beta + lambda)/beta  = c*y/beta
  // have NO subtraction anywhere, so the corner w they produce is
  // dd-relative at every size, where the generic spelling inside 1-arg
  // Log1pmxDd (TwoSum(1, u.hi) + u.lo) degenerates for u near -1: at
  // 1 + u < 2^-53, u.hi rounds to exactly -1 and LogDdAny receives a
  // zero-high pair (NaN -> the exact-0 return at (19, 1e5, 5.2e-21)); at
  // 1 + u merely small, the folded u.lo rivals the high word and LogDdAny
  // drops the cubic in lo/hi (the 1.4e-4 E error at (19, 1e5, 1.73e-19)).
  //
  // Corner select at u.hi < -1/2 (v mirrored; betainv's fixed-frame call
  // reaches the v corner even though beta's routed frame cannot):
  //  * Non-corner lanes keep the generic spelling BIT-IDENTICALLY -- for
  //    u in (-1/2, -1/16] Sterbenz makes TwoSum(1, u.hi) exact with zero
  //    low word, and |u.lo| <= 2^-53 stays far under w.hi >= 1/2, so the
  //    pair is a normalized dd and nothing needed fixing there.
  //  * Ordering (c (x) xi) (/) alpha, NOT (c/alpha) (x) xi: c/alpha
  //    reaches 2^1021 (alpha = Z0 scaled 2^-200 against c_s ~ 2^824) and
  //    would overflow the non-FMA Dekker split (2^996 ceiling,
  //    docs/NUMERICAL-DOCTRINE.md). c (x) xi is split-safe (c_s <= 2^901
  //    unscaled-frame, >= 2^700 whenever the prescale fired, so no live
  //    lane's product goes subnormal either), and the quotient is the
  //    corner w < ~1/2, in range by construction.
  //  * The numerator is scrubbed to alpha (quotient exactly 1) on
  //    non-corner lanes so the division never sees the huge-w lanes
  //    (w = c*xi/alpha is unbounded when xi ~ 1 with alpha << c; those
  //    lanes take the generic spelling anyway).
  // Corner accuracy: w carries ~2^-103 relative (one DdMul + one DdDivDd),
  // log(w) takes that to 2^-103 ABSOLUTE in phi, and alpha*phi <= ~750
  // (saturation) bounds the E contribution at 2^-93-class -- invisible
  // against the 2^-60s budget.
  const auto mzero = op::Zero(d);
  const auto mhalf = op::Set(d, -0.5);
  const auto cor_u = op::Lt(u.hi, mhalf);
  const auto cor_v = op::Lt(v.hi, mhalf);
  const auto cxi = DdMul(d, c, xi);
  const auto cyv = DdMul(d, c, y);
  const Dd<D> num_u{op::IfThenElse(cor_u, cxi.hi, as),
                    op::IfThenElse(cor_u, cxi.lo, mzero)};
  const Dd<D> num_v{op::IfThenElse(cor_v, cyv.hi, bs),
                    op::IfThenElse(cor_v, cyv.lo, mzero)};
  const auto wuc = DdDivDd(d, num_u, Dd<D>{as, mzero});
  const auto wvc = DdDivDd(d, num_v, Dd<D>{bs, mzero});
  auto wug = TwoSum(d, one, u.hi);  // exact (generic spelling)
  wug.lo = op::Add(wug.lo, u.lo);
  auto wvg = TwoSum(d, one, v.hi);  // exact
  wvg.lo = op::Add(wvg.lo, v.lo);
  const Dd<D> wu{op::IfThenElse(cor_u, wuc.hi, wug.hi),
                 op::IfThenElse(cor_u, wuc.lo, wug.lo)};
  const Dd<D> wv{op::IfThenElse(cor_v, wvc.hi, wvg.hi),
                 op::IfThenElse(cor_v, wvc.lo, wvg.lo)};

  // Both phi are >= 0, so this sum never cancels.
  const auto cs = DdAdd(d, DdMulD(d, Log1pmxDd(d, u, wu), as),
                        DdMulD(d, Log1pmxDd(d, v, wv), bs));

  const auto rc = DdRecipDd(d, c);
  // nu = alpha*(beta/c): beta/c <= 1 and alpha*(beta/c) <= min(alpha,beta), so
  // this ordering cannot overflow where alpha*beta would.
  const auto nus = DdMulD(d, DdMulD(d, rc, bs), as);
  const auto z2 = DdMul(d, cs, DdRecipDd(d, nus));

  // Aggregate-initialized, never default-constructed: an implicitly
  // generated default constructor is instantiated OUTSIDE the per-target
  // attribute region on Apple Clang, and the NEON_BF16 slice of a newer
  // system Highway then refuses to inline the vector members'
  // always_inline ctors into it (caught by the first macOS CI run of
  // this TU; gamma/erf never default-construct and were unaffected).
  const BetaPsi<D> out{Dd<D>{op::Mul(cs.hi, up), op::Mul(cs.lo, up)},
                       Dd<D>{op::Mul(nus.hi, up), op::Mul(nus.lo, up)},
                       op::CopySign(DdSqrt(d, z2).hi, lam.hi),
                       op::Div(as, c.hi),
                       op::Div(bs, c.hi),
                       lam.hi};
  return out;
}

// ------------------------------------------------------------------------
// R1: 1 + alpha * sum_{n>=1} t_n/(alpha+n), the series factor of
//     I = xi^alpha/B(alpha,beta) * sum_{n>=0} t_n/(alpha+n),
//     t_n = t_{n-1}*(n-beta)*xi/n.
//
// Multiplying the n >= 1 tail by alpha (and folding -log alpha into the
// exponent instead of dividing) is what keeps a subnormal alpha finite: the
// n = 0 term is 1/alpha, which is +inf there, and the whole series is
// 1 + O(alpha) as alpha -> 0.
//
// TERMS ARE dd, unlike gamma's R1. gamma's terms decrease monotonically from
// t_0 = 1, so a rounded recurrence puts its error only into terms that are
// already small; here the terms RISE until n ~ beta*xi (<= B1 = 8 on the whole
// region) so the dominant term is the eighth or so, and a plain-double
// recurrence would hand it ~8 roundings, i.e. 2^-50. n - beta is an exact
// TwoSum and 1/n is the generator's dd table for the same reason.
//
// 1/(alpha+n) is DdRecipDd of the EXACT TwoSum(alpha, n): alpha+n is not
// representable for a general alpha, and this weight is a factor of the term
// (gamma R4 precedent).
//
// The freeze is per lane and STICKY, and a frozen lane's accumulator is
// SELECTED rather than having a zero added to it: DdAdd of an exact zero is
// value-preserving but RENORMALIZES, and a lane must not be able to tell how
// many renormalizations its neighbours dragged it through. test_beta_smoke's
// lane-mix check is what polices this. The one place the freeze could fire
// early is n = ceil(beta), where n - beta can be ~0 -- but every later term
// carries that same factor, so all of them are equally small and nothing is
// lost. There is exactly one such n, since n - beta has one zero.
//
// HWY_NOINLINE from day one on every core AND on the driver (AGENTS.md): fully
// inlined, each export becomes one enormous function per target and MSVC's
// optimizer is superlinear in function size -- gamma's TU hit the CI Windows
// job's 25-minute timeout inside code generation before its cores were
// outlined, and this TU is heavier still.
template <class D>
HWY_NOINLINE Dd<D> BetaR1Series(D d, op::V<D> a, op::V<D> b, op::V<D> xi) {
  const auto one = op::Set(d, 1.0);
  const auto zero = op::Zero(d);
  const auto half = op::Set(d, 0.5);
  const auto eps = op::Set(d, kBetaFreezeEps);

  // HUGE-beta EXACT PRESCALE [2026-08-12, third defect found by the u -> -1
  // corner reference rows; never sampled before -- R1 with beta within
  // ~100x of DBL_MAX needs beta*xi <= B1, i.e. xi ~ 8e-307-class]. Two
  // independent hazards, one cure:
  //  * The recurrence's grouping t (x) ((n-beta)/n) peaks at |t|*beta/n
  //    BEFORE xi rescales it -- overflow at (19, 1e307, 7.6e-307), where
  //    |t_3|*beta/4 = 1.83e308 > DBL_MAX (witness: NaN out of a healthy
  //    3.55e-4 row). The grouping itself is load-bearing elsewhere and is
  //    NOT changed.
  //  * beta > 2^996 breaks ops::ProdLow's non-FMA Dekker split
  //    (docs/NUMERICAL-DOCTRINE.md) inside the same DdMul.
  // Scaling beta down and xi up by the same exact power of two fixes both:
  // (n - beta)*2^-200 is an exact TwoSum of exactly-scaled operands, the
  // intermediate peak drops to ~2^836, and the xi*2^200 factor restores
  // the true term exactly (xi <= XI1 < 1/2 keeps the upscale finite;
  // subnormal xi upscales exactly). Non-big lanes multiply by 1.0 -- an
  // IEEE identity -- so every previously-validated lane is BIT-IDENTICAL.
  const auto rbig = op::Gt(b, op::Set(d, kBetaScaleAbove));
  const auto rdn = op::IfThenElse(rbig, op::Set(d, kBetaScaleDown), one);
  const auto rup = op::IfThenElse(rbig, op::Set(d, kBetaScaleUp), one);
  const auto bsd = op::Mul(b, rdn);
  const auto xu = op::Mul(xi, rup);

  Dd<D> t{one, zero};
  Dd<D> s{zero, zero};
  auto live = one;
  for (int n = 1; n <= detail::kBetaN1; ++n) {
    const auto nv = op::Set(d, static_cast<double>(n));
    // (n - beta) * 2^-s, EXACT (both scalings exact powers of two)
    const auto nb = TwoSum(d, op::Mul(nv, rdn), op::Neg(bsd));
    const Dd<D> rn{op::Set(d, detail::kBetaRecipNHi[n - 1]),
                   op::Set(d, detail::kBetaRecipNLo[n - 1])};
    t = DdMulD(d, DdMul(d, t, DdMul(d, nb, rn)), xu);
    const auto wgt = DdRecipDd(d, TwoSum(d, a, nv));
    const auto contrib = DdMul(d, t, wgt);
    const auto lm = op::Gt(live, half);
    const auto sn = DdAdd(d, s, contrib);
    s = Dd<D>{op::IfThenElse(lm, sn.hi, s.hi), op::IfThenElse(lm, sn.lo, s.lo)};
    live = op::IfThenElse(
        op::Lt(op::Abs(contrib.hi), op::Mul(op::Abs(s.hi), eps)), zero, live);
    if (op::AllTrue(d, op::Eq(live, zero))) break;
  }
  return DdAddD(d, DdMulD(d, s, a), one);
}

// ------------------------------------------------------------------------
// R2: F, the backward evaluation of DLMF 8.17.22's continued fraction,
//     I = xi^alpha y^beta/(alpha B(alpha,beta)) * F,
//     F = 1/f,  f = 1 + d_1/(1 + d_2/(1 + ...)),
//     d_{2m}   = m(beta-m)xi / ((alpha+2m-1)(alpha+2m)),
//     d_{2m+1} = -(alpha+m)(c+m)xi / ((alpha+2m)(alpha+2m+1)).
// Fixed depth kBetaN2 with NO convergence test (gamma's GammaCfRecip lesson).
//
// THE RECURRENCE AND EVERY d ARE DOUBLE-DOUBLE, which gamma's are not
// [G3 deviation, and the review checklist's "d-term order of operations" item
// is discharged here rather than by a cheaper spelling]. Two independent
// reasons, both measured:
//
//  1. CONDITIONING. A relative perturbation eps in d_1 moves the returned F by
//     (F - 1)*eps -- the recurrence's last step is f_1 = 1 + d_1/f_2 and F is
//     1/f_1, so the closer f_1 sits to zero the more it magnifies. The pinned
//     orientation rule xi < (alpha+1)/(c+2) puts d_1 = -c*xi/(alpha+1) just
//     above -1, which is precisely where f_1 is small: F reaches 200 at
//     (a, b, x) = (7896.5, 7.4438e-6, 0.995) and 1.25e7 at (1, 1e10, 8e-8).
//     A double-precision d there is 2^-53*F, i.e. 90 and 5e6 ULP. Everything
//     that enters d therefore carries dd: alpha+m and the denominators are
//     EXACT TwoSums, and the single division is a dd division.
//  2. THE INEXACT xi. In the swapped orientation xi = 1 - x is inexact
//     whenever x < 1/2, and F is sensitive to that residual in a way the
//     prefactor is not: at (1, 1e10, 8e-8), F is exactly 1/y, so
//     dlnF/dln xi = xi/y = 1.25e7 and a 2^-54 absolute error in xi is 2^-30
//     relative in the answer. (R1 and R4 need no such treatment: their
//     membership forces xi <= 0.45, hence x >= 0.55 in the swapped case,
//     hence a Sterbenz-exact 1 - x.)
//
// ORDER OF OPERATIONS AND OVERFLOW, the design's own rule and it is load
// bearing: each d is TWO factors, a ratio bounded by ~1 times a bounded
// quotient, and never the raw (alpha+m)(c+m) product. alpha is NOT bounded on
// unsaturated lanes -- xi is carried as a dd, so xi can sit 1 - 5e-249 below
// one and alpha*ln xi stays finite for alpha ~ 1e250, at which point
// (alpha+m)*(c+m) is 1e500. (An earlier form here folded the two ratios into a
// single division for speed and produced -nan at (1e-6, 7.7e249, 5.1e-249);
// the reference set's huge-parameter rows caught it.) The two bounded pieces
// are (c+m)*xi <= alpha + 33 -- which is exactly what the orientation rule
// buys -- and (beta-m)*xi, capped near 800 by saturation.
//
// c and beta themselves still reach 1e308, past ops::ProdLow's 2^996 non-FMA
// Dekker ceiling, so the pair (c, xi) is carried on an exact power-of-two
// prescale: c*xi = (c*s)*(xi/s) identically, and with s = 2^-200 both factors
// are in range (xi <= 1 bounds xi/s by 2^200, and the branch only fires when
// c > 2^900 forces xi < 1e-270 anyway). The divisors alpha+2m are handled by
// DdDivDd's own scaling guard.
template <class D>
HWY_NOINLINE Dd<D> BetaR2Cf(D d, op::V<D> a, op::V<D> b, Dd<D> c, Dd<D> xi) {
  const auto one = op::Set(d, 1.0);
  const auto zero = op::Zero(d);

  const auto big = op::Gt(c.hi, op::Set(d, kBetaScaleAbove));
  const auto sdn = op::IfThenElse(big, op::Set(d, kBetaScaleDown), one);
  const auto sup = op::IfThenElse(big, op::Set(d, kBetaScaleUp), one);
  const Dd<D> cs{op::Mul(c.hi, sdn), op::Mul(c.lo, sdn)};
  // beta*s can flush a subnormal beta to zero, but only when c > 2^900 forces
  // alpha > 2^900 as well, where beta - m is m to 300 digits: the flush IS the
  // correct value there.
  const auto bs = op::Mul(b, sdn);
  const Dd<D> xis{op::Mul(xi.hi, sup), op::Mul(xi.lo, sup)};

  Dd<D> f{one, zero};
  for (int k = detail::kBetaN2; k >= 1; --k) {
    // m = floor(k/2) BY CONSTRUCTION (d_{2m} at even k, d_{2m+1} at odd);
    // the integer division is the point.
    // NOLINTNEXTLINE(bugprone-integer-division)
    const double mv = static_cast<double>(k / 2);
    const auto m = op::Set(d, mv);
    Dd<D> dk{zero, zero};  // brace-init: see BetaPsiCore's ctor note
    if ((k & 1) == 0) {
      // d_{2m} = [m/(alpha+2m-1)] * [(beta-m)xi/(alpha+2m)]
      const auto bx = DdSub(d, DdMulD(d, xis, bs), DdMulD(d, xi, m));
      const auto r1 = DdDivDd(d, Dd<D>{m, zero},
                              TwoSum(d, a, op::Set(d, 2.0 * mv - 1.0)));
      const auto r2 =
          DdDivDd(d, bx, TwoSum(d, a, op::Set(d, 2.0 * mv)));
      dk = DdMul(d, r1, r2);
    } else {
      // d_{2m+1} = -[(alpha+m)/(alpha+2m)] * [(c+m)xi/(alpha+2m+1)]
      const auto cx = DdAdd(d, DdMul(d, xis, cs), DdMulD(d, xi, m));
      const auto r1 = DdDivDd(d, TwoSum(d, a, m),
                              TwoSum(d, a, op::Set(d, 2.0 * mv)));
      const auto r2 =
          DdDivDd(d, cx, TwoSum(d, a, op::Set(d, 2.0 * mv + 1.0)));
      dk = DdNeg(d, DdMul(d, r1, r2));
    }
    f = DdAddD(d, DdMul(d, dk, DdRecipDd(d, f)), one);
  }
  return DdRecipDd(d, f);
}

// ------------------------------------------------------------------------
// R3: Temme's erf-form on the ridge. Gamma's GammaTemme with cpsi in place of
// a*phi; z^2 = cpsi parallels z^2 = a*phi exactly.
//
//   direct = 1/2 erfc(z) + sg * e^{-cpsi}/sqrt(2 pi nu) * S(zeta, p, 1/nu),
//   z = sqrt(cpsi),  zeta = sign(lambda) sqrt(cpsi/nu),  sg = -sign(lambda),
//   S = sum_k e_k(zeta,p)/nu^k.
// The sign convention is the generator's (r3_R_at): zeta carries
// sign(lambda) = sign(alpha - c*xi), which is OPPOSITE gamma's eta -- the
// gamma-limit mapping is eta = -zeta*sqrt(2) (cross-check (d)).
//
// S is the tensor table in beta_data.h, evaluated exactly as that header's
// recipe states: nested Clenshaw, in the p direction inside each of the 25
// zeta rows and in the zeta direction across them, then Horner in 1/nu across
// the 10 orders. The stored half is p <= 1/2; for p > 1/2 the table is
// evaluated at (-zeta, q) and NEGATED (self-check (h),
// e_k(zeta,p) = -e_k(-zeta,1-p)).
//
// FIVE THINGS THAT WOULD OTHERWISE COST THE REGION ITS ACCURACY -- gamma's
// list, carried over:
//  1. lambda is assembled exactly; see BetaPsiCore.
//  2. cpsi comes from Log1pmxDd, never from a naive u - log(1+u).
//  3. e^{-z^2} is e^{-cpsi} taken from the DD cpsi, NEVER from squaring the
//     rounded z: near z = 6 a half ulp of z is 2^-48 in z^2 and therefore
//     2^-48 RELATIVE in the exponential.
//  4. z is a dd and its low word is the first-order correction to erfc,
//     1/2 erfc(z.hi + z.lo) = 1/2 erfc(z.hi) - z.lo/sqrt(pi) e^{-z^2}, which
//     folds into the same bracket S already multiplies.
//  5. The core/tail split is on cpsi against an exact 36, i.e. z against 6 --
//     erfc's own split. Above it the whole expression carries e^{-cpsi}, so
//     erfc is taken in its tail form e^{-z^2}G(1/z)/z (reusing
//     erfc_tail_data's G; z <= sqrt(800) ~ 28.3 stays inside its fitted range)
//     and the power-of-two scaling is applied LAST.
//
// WHICH SIDE R3 COMPUTES IS DECIDED BY lambda, NOT BY THE ROUTER. The router
// orients on xi/p <= 1, which is the same predicate as lambda >= 0 in exact
// arithmetic but is evaluated on a ratio carrying three roundings; when the
// two disagree (|xi/p - 1| below ~2^-52, i.e. a point sitting on the ridge
// crest) this core still returns the genuinely smaller side, and the driver
// must be told which one that is. Hence lam_pos, and hence the XNOR against
// the router's own orientation at the call site. Gamma's GammaTemme reports
// is_p for the same reason; the difference is only that gamma's predicate and
// its sign come from the same subtraction, so the two can never disagree
// there. Without this the crest points come out as the OTHER member of the
// pair -- 71.111, 640, 0.1 (where xi/p rounds to exactly 1.0 while lambda is
// negative) returned 0.4867 for a true P of 0.5133.
template <class D>
struct BetaR3Out {
  BetaVal<D> val;
  op::V<D> lam_pos;  // 1.0 where lambda >= 0, i.e. where val is I_xi(alpha,beta)
  op::M<D> sat;
  // The TAIL branch's bracket, so a consumer that wants ln(val) rather than
  // val can take it in LOG space: val = e^-cpsi * brk there, so
  // ln val = -cpsi (+) ln brk exactly, with no exponential in between. This
  // file's own driver does not use it -- it wants the value -- but the
  // INVERSE does, and it must: past cpsi ~ 708 the assembled value is
  // subnormal and its log is quantized to a few bits, which for an inverse
  // is not a small error in the answer but a flat spot in the residual (two
  // y's 1e-8 apart returning the identical logit). Added for that consumer;
  // no arithmetic here changed, this field is a value already computed.
  Dd<D> brk;
};

template <class D>
HWY_INLINE op::V<D> BetaClenshawP(D d, const double* row, op::V<D> u,
                                  op::V<D> u2) {
  auto b1 = op::Set(d, row[detail::kBetaR3NP - 1]);
  auto b2 = op::Zero(d);
  for (int m = detail::kBetaR3NP - 2; m >= 1; --m) {
    const auto nb = op::Sub(op::MulAdd(u2, b1, op::Set(d, row[m])), b2);
    b2 = b1;
    b1 = nb;
  }
  return op::Sub(op::MulAdd(u, b1, op::Set(d, row[0])), b2);
}

// i_gl: 1.0 on gamma-limit lanes (max(alpha,beta) >= kBetaGammaLim), which
// additionally evaluate the p->0-edge depth-extension rows -- see below.
template <class D>
HWY_NOINLINE BetaR3Out<D> BetaR3Temme(D d, const BetaPsi<D>& ps,
                                      op::V<D> i_gl) {
  const auto one = op::Set(d, 1.0);
  const auto zero = op::Zero(d);
  const auto half = op::Set(d, 0.5);

  // Table symmetry: only p <= 1/2 is stored.
  const auto pswap = op::Gt(ps.p, half);
  const auto pe = op::IfThenElse(pswap, ps.q, ps.p);
  const auto ze = op::IfThenElse(pswap, op::Neg(ps.zeta), ps.zeta);
  const auto tsign = op::IfThenElse(pswap, op::Set(d, -1.0), one);

  const auto t = op::Mul(ze, op::Set(d, 1.0 / detail::kBetaZetaMax));
  const auto t2 = op::Add(t, t);
  const auto u = op::Mul(op::Sub(pe, op::Set(d, detail::kBetaR3PMid)),
                         op::Set(d, 1.0 / detail::kBetaR3PHalf));
  const auto u2 = op::Add(u, u);
  const auto r = op::Div(one, ps.nu.hi);

  auto sum = zero;
  for (int k = detail::kBetaR3K - 1; k >= 0; --k) {
    auto b1 = zero;
    auto b2 = zero;
    for (int n = detail::kBetaR3NZ - 1; n >= 1; --n) {
      const auto rn = BetaClenshawP(d, detail::kBetaR3Cheb[k][n], u, u2);
      const auto nb = op::Sub(op::MulAdd(t2, b1, rn), b2);
      b2 = b1;
      b1 = nb;
    }
    const auto r0 = BetaClenshawP(d, detail::kBetaR3Cheb[k][0], u, u2);
    const auto row = op::Sub(op::MulAdd(t, b1, r0), b2);
    sum = op::MulAdd(sum, r, op::Mul(tsign, row));
  }

  // [(C) slice DEPTH EXTENSION]: gamma-limit ridge lanes add the p->0-edge
  // rows k = 10..12 (kBetaR3GlExt, 1D Chebyshev in the same t) -- the main
  // K=10 truncation alone is 2^-50-class at nu = kBetaGlRidgeMin, 2^-60
  // with these (generator check (c)'s extension lattice). The rows are
  // p-edge values, applied at EVERY slice nu: whenever they are
  // non-negligible (small nu) the lane's own p <= nu/kBetaGammaLim is
  // tiny, so the edge coefficients are the right ones by construction;
  // at large nu the r^10 weight erases them. tsign carries the p > 1/2
  // symmetry exactly as for the main rows.
  {
    auto es = zero;
    for (int k = detail::kBetaR3GlK - 1; k >= 0; --k) {
      auto b1 = zero;
      auto b2 = zero;
      for (int n = detail::kBetaR3NZ - 1; n >= 1; --n) {
        const auto nb = op::Sub(
            op::MulAdd(t2, b1, op::Set(d, detail::kBetaR3GlExt[k][n])), b2);
        b2 = b1;
        b1 = nb;
      }
      const auto row = op::Sub(
          op::MulAdd(t, b1, op::Set(d, detail::kBetaR3GlExt[k][0])), b2);
      es = op::MulAdd(es, r, row);
    }
    const auto r2 = op::Mul(r, r);
    const auto r4 = op::Mul(r2, r2);
    const auto r10 = op::Mul(op::Mul(r4, r4), r2);
    sum = op::MulAdd(op::Mul(i_gl, r10), op::Mul(tsign, es), sum);
  }

  const auto z = DdSqrt(d, ps.cpsi);
  const Dd<D> twopi{op::Set(d, kBetaTwoPiHi), op::Set(d, kBetaTwoPiLo)};
  const auto rv = DdRecipDd(
      d, DdSqrt(d, DdMulD(d, twopi,
                          op::Min(ps.nu.hi, op::Set(d, kBetaTwoPiNuClamp)))));
  const auto sg = op::Neg(op::CopySign(one, ps.lam));
  const auto s_rv = DdMulD(d, rv, op::Mul(sg, sum));

  const auto ea = BetaClampE(d, op::Neg(ps.cpsi.hi), op::Neg(ps.cpsi.lo));

  // --- core, cpsi <= 36 ---------------------------------------------------
  // ErfcCoreDd's table index is round(min(z, 6+1/1024)*256) and is NOT masked,
  // so a NaN z from a discarded lane must be scrubbed before it -- the same
  // one-op guard erfinv's HalleyMid carries (AGENTS.md, value-derived
  // gathers). Only discarded lanes are affected.
  const auto zs = op::IfThenElse(op::IsNaN(z.hi), zero, z.hi);
  const auto ec = ErfcCoreDd(d, zs, zs);
  const Dd<D> half_erfc{op::Mul(ec.hi, half), op::Mul(ec.lo, half)};
  const auto exd = BetaExpDd(d, ea.hi, ea.lo);
  const auto brk_core = DdAddD(
      d, s_rv, op::Neg(op::Mul(z.lo, op::Set(d, kBetaInvSqrtPiHi))));
  const auto core = DdAdd(d, half_erfc, DdMul(d, brk_core, exd));

  // --- tail, cpsi > 36 ----------------------------------------------------
  const auto uz = op::Div(one, zs);
  const auto gt = op::Div(ErfcTailGFromU(d, zs, uz), op::Add(zs, zs));
  const auto brk_tail = DdAddD(d, s_rv, gt);
  const auto exf = ExpDdFrac(d, ea.hi, ea.lo);
  const auto tail = BetaScale(d, DdMul(d, exf.m, brk_tail), exf.e);

  const auto tm = op::Gt(ps.cpsi.hi, op::Set(d, kBetaZ2Split));
  const BetaVal<D> val{op::IfThenElse(tm, tail.v, DdToDouble(core)),
                       Dd<D>{op::IfThenElse(tm, tail.dd.hi, core.hi),
                             op::IfThenElse(tm, tail.dd.lo, core.lo)}};
  return {val, BetaInd(d, op::Ge(ps.lam, zero)), ea.sat, brk_tail};
}

// ------------------------------------------------------------------------
// R4: the tiny-min corner, evaluated TINY-FIRST as (tau, B, xi_tau) and
// producing the SMALL side directly:
//
//   Qtilde = 1 - I_{xi_tau}(tau, B) = -expm1(w + log1p(tau*Sigma)),
//   w      = tau*ln xi_tau + [lgamma(B+tau) - lgamma(B)] - lgamma(1+tau),
//   Sigma  = sum_{n>=1} t_n/(tau+n),   t_n = t_{n-1}(n-B)xi_tau/n.
//
// This is the design's "gamma-R4 verbatim in beta clothing" requirement met in
// LOG SPACE, which is strictly stronger than expanding products of (1+E_i):
// I_{xi_tau}(tau,B) tends to 1 as tau -> 0, and nothing here ever forms it.
// Every piece of w is individually O(tau) and individually relative-accurate:
//   * the lgamma difference is analytic (LgammaDiffDd) -- lgamma(B+tau) and
//     lgamma(B) are equal to the last bit once tau is small, so their
//     difference cannot be taken numerically;
//   * lgamma(1+tau) is the analytic difference LgammaDiffDd(1, tau) [NINTH
//     correction; see the block comment at its call below], never
//     LgammaPosDd(1+tau) -- fl(1+tau) rounds to exactly 1 below 2^-53 and
//     would return a flat zero, losing the -gamma*tau that is the whole
//     signal (the GammaSmallQ mechanism);
//   * log1p(tau*Sigma) is a log1p, not a log of 1 + something.
// The identity Sigma*e^w telescoping to the BPSER value is the reference
// generator's own small_tau_oracle, validated there at 2^-188 against the
// continued fraction on their overlap band.
//
// tau*|ln xi_tau| <= ln 2 on R4's own box, so |w| stays below ~1 there.
// POST-ROUTED lanes [SEVENTH correction; EIGHTH widens the gate] have no
// ln-2 cap and tau up to
// kBetaPrTauMax = 2.5, but their defining property (R1 value > kBetaNearOne)
// means w and log1p(tau*Sigma) CANCEL to below ~2^-10 -- the dd addition
// carries that cancellation at ~2^-105 absolute, so the relative error of
// the small result stays ~2^-94-class, and Expm1Dd receives the already-
// combined tiny argument in its series branch. The series itself is
// tau-benign (tau appears only in the 1/(tau+n) weights and the exact
// assembly); the generator's check (f) post-route lattice proves the
// N = kBetaR4N truncation over the widened domain.
template <class D>
HWY_NOINLINE Dd<D> BetaR4Tiny(D d, op::V<D> tau, op::V<D> bb, op::V<D> xi,
                              Dd<D> lxi) {
  const auto one = op::Set(d, 1.0);
  const auto zero = op::Zero(d);
  const auto half = op::Set(d, 0.5);
  const auto eps = op::Set(d, kBetaFreezeEps);

  // --- subnormal-tau guard [ELEVENTH correction] --------------------------
  // For tau < 2^-950 every term of w is tau-scale and the intermediate
  // PRODUCTS land at or below the subnormal boundary, where the non-FMA
  // Dekker residual in ops::ProdLow collapses (sub-products round at
  // 2^-1074 granularity): measured 2709 ULP at
  // (1.57e-311, 4.6e-210, 0.44) on the SSE4 capped sweep, clean on FMA
  // targets. Two exact power-of-two reframings, selected per lane:
  //  * tau < 2^-950, bb > 2^-140: the assembly is linear in tau to
  //    relative O(tau_w/bb) <= 2^-110, so run it ENTIRELY at
  //    tau_w = tau * 2^700 (in [2^-374, 2^-250]: every product normal)
  //    and unscale the result once at the end -- the final multiply is
  //    the single subnormal rounding the contract allows.
  //  * bb <= 2^-140 (both-tiny corner; tau <= bb since tau is the min):
  //    w = -log1p(r) + tau*(ln xi + Sigma + ...) with r = tau/bb, and the
  //    tau-linear remainder is relatively bounded by bb*|ln xi| <=
  //    2^-140 * 745 < 2^-130, so Qtilde = r/(1+r) exactly, with r from a
  //    JOINTLY 2^900-scaled division (both operands exact, all normal;
  //    r <= 1 because tau <= bb). Handled after the main assembly below.
  const auto tsub = op::Gt(op::Set(d, kBetaTauSubn), tau);
  const auto btiny = op::Ge(op::Set(d, kBetaBTinyCut), bb);
  const auto tau0 = tau;  // the shortcut divides the ORIGINAL pair
  tau = op::IfThenElse(tsub, op::Mul(tau, op::Set(d, kBetaTauUp)),
                       tau);

  Dd<D> t{one, zero};
  Dd<D> s{zero, zero};
  auto live = one;
  for (int n = 1; n <= detail::kBetaR4N; ++n) {
    const auto nv = op::Set(d, static_cast<double>(n));
    const auto nb = TwoSum(d, nv, op::Neg(bb));  // n - B, EXACT
    const Dd<D> rn{op::Set(d, detail::kBetaRecipNHi[n - 1]),
                   op::Set(d, detail::kBetaRecipNLo[n - 1])};
    t = DdMulD(d, DdMul(d, t, DdMul(d, nb, rn)), xi);
    const auto wgt = DdRecipDd(d, TwoSum(d, tau, nv));
    const auto contrib = DdMul(d, t, wgt);
    const auto lm = op::Gt(live, half);
    const auto sn = DdAdd(d, s, contrib);
    s = Dd<D>{op::IfThenElse(lm, sn.hi, s.hi), op::IfThenElse(lm, sn.lo, s.lo)};
    live = op::IfThenElse(
        op::Lt(op::Abs(contrib.hi), op::Mul(op::Abs(s.hi), eps)), zero, live);
    if (op::AllTrue(d, op::Eq(live, zero))) break;
  }

  // lgamma(1+tau) via the SAME analytic difference machinery as the
  // lgamma(B+tau) - lgamma(B) term below [NINTH correction]:
  //   lgamma(1+tau) = LgammaDiffDd(1, tau)      (tau <= 1)
  //                 = LgammaDiffDd(2, tau - 1)  (1 < tau <= kBetaPrTauMax)
  // exact identities, since lgamma(1) = lgamma(2) = 0 and tau - 1 is an
  // exact subtraction for tau in (1, 2.5] (integer off a finer grid).
  // m = tau or tau-1 stays in (0, 3/2] <= M in {1, 2}: inside
  // LgammaDiffDd's precondition, its walk range (M < Z0 fires <= 7 of
  // the 10 unit steps), and its already-exercised Log1pmxDd argument
  // range (the PA call sites reach comparable s).
  // WHY NOT the previous three-zone ZoneBracket form: that reused
  // lgamma's zone polynomial, whose replayed-evaluation budget
  // (ZONE_TARGET = 3e-17, floored by the tail Horner's own DOUBLE
  // rounding) is relative to lg1 ITSELF. Post-routed lanes cancel w
  // down to Qtilde ~ 2^-11..2^-17, amplifying that budget by
  // |lg1|/Qtilde ~ 2^13..2^17: measured 55/209 ULP against the
  // certified reference set, implied dw = (1..3)e-18 x lg1 -- squarely
  // the polynomial's class, twenty orders above the dd class of every
  // other component in w. Same disease shape as the reference oracle's
  // round-6 lg_diff defect: a truncation/rounding budget proven
  // relative to a component's own scale is void once the assembly
  // cancels below it. Original R4 lanes (tau <= eps_R4) were immune --
  // there |lg1| ~ gamma*tau and Qtilde ~ tau*O(10), no amplification --
  // but take the analytic form too: one code path, and the walk is
  // cheap next to the N=48 series above.
  const auto gt1 = op::Gt(tau, one);
  const auto lg1 = LgammaDiffDd(
      d, op::IfThenElse(gt1, op::Add(one, one), one),
      op::IfThenElse(gt1, op::Sub(tau, one), tau));

  auto w = DdMulD(d, lxi, tau);
  w = DdAdd(d, w, LgammaDiffDd(d, bb, tau));
  w = DdSub(d, w, lg1);
  w = DdAdd(d, w, Log1pDdWide(d, DdMulD(d, s, tau)));
  auto q = DdNeg(d, Expm1Dd(d, w));

  // [ELEVENTH correction] unscale the tau_w lanes (exact except the one
  // allowed final subnormal rounding), then override the both-tiny lanes
  // with the r/(1+r) closed form. Dead-lane execution of the shortcut can
  // produce inf/NaN from the 2^900 up-scale of a large bb; those lanes are
  // discarded by the select and feed no gather (masked-execution rule).
  const auto dnw = op::IfThenElse(tsub, op::Set(d, kBetaTauDown), one);
  q = Dd<D>{op::Mul(q.hi, dnw), op::Mul(q.lo, dnw)};
  const auto sc = op::Set(d, kBetaTinyPairUp);
  const auto r = DdDivDd(d, Dd<D>{op::Mul(tau0, sc), zero},
                         Dd<D>{op::Mul(bb, sc), zero});
  const auto qs = DdMul(d, r, DdRecipDd(d, DdAddD(d, r, one)));
  const auto mbt =
      BetaIndMask(d, op::Mul(BetaInd(d, tsub), BetaInd(d, btiny)));
  return Dd<D>{op::IfThenElse(mbt, qs.hi, q.hi),
               op::IfThenElse(mbt, qs.lo, q.lo)};
}

// ------------------------------------------------------------------------
// The driver. kP selects beta_p (true) or beta_q (false); the two differ only
// in the final handout, since the direct side is decided by the region.
//
// HWY_NOINLINE like the region cores, and for the same MSVC-codegen reason:
// even with the log/exp calls routed through BetaLog/BetaExpDd above, this
// function still inlines LgammaPosDd, ExpDdFrac twice and the dd assembly,
// and inlining it into the two export loops would re-create a function big
// enough to stall cl.exe's back end past the CI timeout.
template <bool kP, class D>
HWY_NOINLINE op::V<D> BetaVec(D d, op::V<D> a_in, op::V<D> b_in,
                              op::V<D> x_in) {
  const auto one = op::Set(d, 1.0);
  const auto zero = op::Zero(d);
  const auto half = op::Set(d, 0.5);
  const auto inf = op::Set(d, std::numeric_limits<double>::infinity());
  const auto qnan = op::Set(d, std::numeric_limits<double>::quiet_NaN());

  // --- scrub every lane whose result is decided by the specials table -----
  // Masked-off lanes still execute every op (AGENTS.md), and this kernel takes
  // logs of xi and 1-xi, divides by alpha and beta, and multiplies (c+m): a
  // NaN, a zero, a one or an infinity left in place would propagate through
  // gathers and dd residuals rather than sit quietly.
  const auto safe_a = op::Set(d, kBetaSafeA);
  const auto safe_b = op::Set(d, kBetaSafeB);
  const auto safe_x = op::Set(d, kBetaSafeX);
  auto a = a_in;
  auto b = b_in;
  auto x = x_in;
  {
    const auto m = op::IsNaN(op::Add(op::Add(a, b), x));
    a = op::IfThenElse(m, safe_a, a);
    b = op::IfThenElse(m, safe_b, b);
    x = op::IfThenElse(m, safe_x, x);
  }
  {
    // a <= 0 or b <= 0 or a = inf or b = inf or x outside (0, 1). The ops
    // facade exposes no mask OR, so the disjunction rides on the same 1.0/0.0
    // indicator device the region map uses.
    auto bad = BetaInd(d, op::Ge(zero, op::Min(a, b)));
    bad = op::Max(bad, BetaInd(d, op::Eq(op::Max(a, b), inf)));
    bad = op::Max(bad, BetaInd(d, op::Ge(zero, x)));
    bad = op::Max(bad, BetaInd(d, op::Ge(x, one)));
    const auto mm = BetaIndMask(d, bad);
    a = op::IfThenElse(mm, safe_a, a);
    b = op::IfThenElse(mm, safe_b, b);
    x = op::IfThenElse(mm, safe_x, x);
  }

  const auto y0 = TwoSum(d, one, op::Neg(x));  // 1 - x, EXACT
  const auto c_raw = op::Add(a, b);            // may be +inf; only compared

  // --- routing (route_final, in order) ------------------------------------
  // Every predicate below is written so that an overflowing intermediate still
  // decides correctly: b*x and B*xi_tau go to +inf and fail their <= B1 test,
  // nu is the harmonic form 1/(1/a+1/b), the ratio band is xi*(1 + b/a), and
  // the R2 threshold is 1/(1 + (b+1)/(a+1)).
  const auto tau = op::Min(a, b);
  const auto bmax = op::Max(a, b);
  const auto sw4 = op::Lt(b, a);  // tiny-first
  const Dd<D> xt{op::IfThenElse(sw4, y0.hi, x),
                 op::IfThenElse(sw4, y0.lo, zero)};
  const auto lxt = BetaLog(d, xt);

  const auto xi1 = op::Set(d, detail::kBetaXi1);
  const auto b1v = op::Set(d, detail::kBetaB1);
  // R2's orientation threshold, evaluated in the TINY-FIRST frame.
  const auto thr_t = op::Div(
      one, op::Add(one, op::Div(op::Add(bmax, one), op::Add(tau, one))));
  // R4's xi cap, WIDENED [G3 correction, fourth routing fix -- see the report
  // and the escalation note below]. The design's box caps xi_tau at xi1 = 0.45
  // and asserts that everything failing the box "lands correctly in R1/R2
  // below". That is true of the other two caps but NOT of this one: for
  // xi_tau in (xi1, thr_t) neither R1 orientation fires (both xi_tau and
  // 1 - xi_tau exceed xi1) and R2 then evaluates the TINY-FIRST side, which
  // for tau <= eps_R4 is 1 - O(tau) -- the near-one member of the pair, so the
  // complement it hands back is pure noise. Witness from the shipped reference
  // set: (a, b, x) = (1.7614e-86, 5.8567e-3, 0.45744) wants Q = 3.0105e-84 and
  // got 2.0143e-23. The window is narrow -- thr_t = (tau+1)/(tau+B+2) exceeds
  // xi1 only for B < ~0.24, and is always below 1/2 -- and inside it the
  // series is untroubled: measured worst N = 48 truncation over the whole
  // extension is 2^-57.3 (against 2^-71.9 inside the design's own box and a
  // generator target of 2^-58), because B < 0.24 keeps every term ratio below
  // xi_tau. The other two caps are unchanged.
  auto i_r4 = op::Mul(
      BetaInd(d, op::Ge(op::Set(d, detail::kBetaEpsR4), tau)),
      op::Mul(BetaInd(d, op::Ge(op::Set(d, detail::kBetaLn2),
                                op::Mul(tau, op::Abs(lxt.hi)))),
              op::Mul(op::Max(BetaInd(d, op::Ge(xi1, xt.hi)),
                              BetaInd(d, op::Lt(xt.hi, thr_t))),
                      BetaInd(d, op::Ge(b1v, op::Mul(bmax, xt.hi))))));

  const auto r1n = op::Mul(BetaInd(d, op::Ge(xi1, x)),
                           BetaInd(d, op::Ge(b1v, op::Mul(b, x))));
  const auto r1s = op::Mul(BetaInd(d, op::Ge(xi1, y0.hi)),
                           BetaInd(d, op::Ge(b1v, op::Mul(a, y0.hi))));
  auto i_r1 = op::Mul(op::Max(r1n, r1s), BetaIndNot(d, i_r4));

  const auto nu_r = op::Div(one, op::Add(op::Div(one, a), op::Div(one, b)));
  const auto rat1 = op::Mul(x, op::Add(one, op::Div(b, a)));       // xi/p
  const auto rat2 = op::Mul(y0.hi, op::Add(one, op::Div(a, b)));   // y/q
  const auto lo_r = op::Set(d, detail::kBetaXiRatioLo);
  const auto hi_r = op::Set(d, detail::kBetaXiRatioHi);
  const auto band = op::Mul(
      op::Mul(BetaInd(d, op::Ge(rat1, lo_r)), BetaInd(d, op::Ge(hi_r, rat1))),
      op::Mul(BetaInd(d, op::Ge(rat2, lo_r)), BetaInd(d, op::Ge(hi_r, rat2))));
  // [(C) gamma-limit slice, ridge part]: above kBetaGammaLim the CF is
  // structurally degenerate, so the in-band ridge floor drops from T_ridge
  // to kBetaGlRidgeMin = 20 -- exactly gamma's own kGammaAT, and exactly
  // the p -> 0 edge where this table's e_k are anchored to gamma's
  // validated c_k (generator check (c)'s extension lattice proves the
  // 1/nu extrapolation there; gamma's own table is likewise applied down
  // to a = 20 from a ladder extracted far higher).
  const auto i_glr = op::Mul(
      BetaInd(d, op::Ge(bmax, op::Set(d, detail::kBetaGammaLim))),
      BetaInd(d, op::Ge(nu_r, op::Set(d, detail::kBetaGlRidgeMin))));
  auto i_r3 = op::Mul(
      op::Mul(op::Max(BetaInd(d, op::Ge(nu_r, op::Set(d, detail::kBetaTRidge))),
                      i_glr),
              band),
      op::Mul(BetaIndNot(d, i_r4), BetaIndNot(d, i_r1)));

  const auto thr = op::Div(
      one, op::Add(one, op::Div(op::Add(b, one), op::Add(a, one))));
  auto i_r2 = op::Mul(BetaIndNot(d, i_r4),
                      op::Mul(BetaIndNot(d, i_r1), BetaIndNot(d, i_r3)));

  // Orientation, one indicator over disjoint regions.
  const auto i_sw =
      op::MulAdd(i_r4, BetaInd(d, sw4),
                 op::MulAdd(i_r1, op::Mul(BetaIndNot(d, r1n), r1s),
                            op::MulAdd(i_r3, BetaInd(d, op::Gt(rat1, one)),
                                       op::Mul(i_r2, BetaIndNot(
                                           d, BetaInd(d, op::Lt(x, thr)))))));
  const auto m_sw = BetaIndMask(d, i_sw);

  const auto alpha = op::IfThenElse(m_sw, b, a);
  const auto beta = op::IfThenElse(m_sw, a, b);
  const Dd<D> xi{op::IfThenElse(m_sw, y0.hi, x),
                 op::IfThenElse(m_sw, y0.lo, zero)};
  const Dd<D> yv{op::IfThenElse(m_sw, x, y0.hi),
                 op::IfThenElse(m_sw, zero, y0.lo)};

  const auto m_r1 = BetaIndMask(d, i_r1);
  const auto m_r2 = BetaIndMask(d, i_r2);
  const auto m_r3 = BetaIndMask(d, i_r3);
  const auto m_r4 = BetaIndMask(d, i_r4);

  // The direct side is P where the native orientation was evaluated -- except
  // in R4, which computes the COMPLEMENT of its own (tiny-first) triple. R3
  // overrides this below from lambda's exact sign (see BetaR3Out).
  const auto i_nat = BetaIndNot(d, i_sw);
  auto is_p = op::IfThenElse(m_r4, i_sw, i_nat);

  // --- prefactor path split -----------------------------------------------
  const auto i_e = op::Max(i_r1, i_r2);
  const auto i_pb = op::Mul(
      i_e, op::Mul(BetaInd(d, op::Gt(c_raw, op::Set(d, detail::kBetaClg))),
                   BetaInd(d, op::Ge(op::Min(alpha, beta),
                                     op::Set(d, detail::kBetaZ0)))));
  const auto i_pa = op::Mul(i_e, BetaIndNot(d, i_pb));
  const auto m_pa = BetaIndMask(d, i_pa);
  const auto m_pb = BetaIndMask(d, i_pb);
  const auto m_psi = BetaIndMask(d, op::Max(i_r3, i_pb));

  // --- shared cpsi machinery (R3 and PB) ----------------------------------
  // Scrubbed outside its own lanes: the prescale's no-underflow argument rests
  // on min(alpha,beta) >= Z0, which only holds where this is live.
  // Aggregate-initialized (see BetaPsiCore's ctor note); the scrub values
  // are unchanged: cpsi 0, nu 1, zeta 0, p = q = 1/2, lam 0.
  BetaPsi<D> ps{Dd<D>{zero, zero}, Dd<D>{one, zero}, zero, half, half, zero};
  const bool need_psi = !op::AllFalse(d, m_psi);
  if (need_psi) {
    const auto pa = op::IfThenElse(m_psi, alpha, safe_a);
    const auto pb = op::IfThenElse(m_psi, beta, safe_b);
    const Dd<D> px{op::IfThenElse(m_psi, xi.hi, safe_x),
                   op::IfThenElse(m_psi, xi.lo, zero)};
    const Dd<D> py{op::IfThenElse(m_psi, yv.hi, op::Set(d, 1.0 - kBetaSafeX)),
                   op::IfThenElse(m_psi, yv.lo, zero)};
    ps = BetaPsiCore(d, pa, pb, px, py);
  }

  // --- assembled per-lane direct side -------------------------------------
  BetaVal<D> val{zero, Dd<D>{zero, zero}};

  if (!op::AllFalse(d, m_r1) || !op::AllFalse(d, m_r2)) {
    // HUGE-parameter EXACT PRESCALE for the two DdMulD multiplicands
    // [2026-08-12, with BetaR1Series' huge-beta fix]: alpha or beta above
    // ops::ProdLow's 2^996 non-FMA Dekker ceiling breaks the split inside
    // DdMulD (a NaN, not a rounding -- and a NaN e defeats the saturation
    // clamp's comparison, so it REACHES the output). beta = 1e307 lanes
    // are live and unsaturated (l2 = beta*ln(y) ~ -7.6 on the corner
    // rows); alpha-huge lanes are always saturated but must still compute
    // a clean -inf-side e for the clamp to see. Scale-back by an exact
    // power of two; a true overflow scale-back lands on +-inf, which the
    // clamp handles (both logs are <= 0, so no inf - inf exists here).
    // Non-big lanes multiply by 1.0: bit-identical.
    const auto pba = op::Gt(alpha, op::Set(d, kBetaScaleAbove));
    const auto pbb = op::Gt(beta, op::Set(d, kBetaScaleAbove));
    const auto sa = op::IfThenElse(pba, op::Set(d, kBetaScaleDown), one);
    const auto ua = op::IfThenElse(pba, op::Set(d, kBetaScaleUp), one);
    const auto sb = op::IfThenElse(pbb, op::Set(d, kBetaScaleDown), one);
    const auto ub = op::IfThenElse(pbb, op::Set(d, kBetaScaleUp), one);
    const auto l1s = DdMulD(d, BetaLog(d, xi), op::Mul(alpha, sa));
    const Dd<D> l1{op::Mul(l1s.hi, ua), op::Mul(l1s.lo, ua)};
    const auto l2s = DdMulD(d, BetaLog(d, yv), op::Mul(beta, sb));
    const Dd<D> l2{op::Mul(l2s.hi, ub), op::Mul(l2s.lo, ub)};
    const auto la = BetaLog(d, alpha);

    // PA: -ln B = LgammaDiffDd(max, min) - lgamma(min). Scrubbed to (3, 2).
    //
    // lgamma(min) [TENTH correction]: for min <= kBetaPrTauMax it is the
    // analytic identity lgamma(min) = lgamma(1+min) - ln(min), with
    // lgamma(1+min) via the NINTH correction's LgammaDiffDd forms -- every
    // piece dd-relative, so the assembled lgamma(min) is dd-ABSOLUTE
    // (~2^-100; the cancellation at lgamma's zeros min = 1, 2 only ever
    // cancels to an absolute level dd carries). WHY: LgammaPosDd's
    // zone/Stirling budgets (3e-17 / 1e-18) are relative to lgamma
    // ITSELF with a double-class Horner floor, and a near-one R1 lane
    // hands its CMP side out at ulp(complement) ~ 2^-62..2^-64 --
    // measured 13 ULP at (a, b, x) = (0.5, 100, 0.05), implied
    // prefactor error 5e-18 = 5.1e-18 * |lgamma(0.5)|, squarely the
    // polynomial class (the R4-postroute/NINTH mechanism at its second
    // site). Near-one values require tau <= kBetaPrTauMax (the eighth-
    // correction pocket: above it the box-corner complement clears the
    // doctrine bound with margin), so min > 2.5 keeps LgammaPosDd; if
    // the G4 seam sweep exposes the tau > 2.5 corner, extend via
    // lgamma(min) = LgammaDiffDd(floor(min), frac) + ln((floor-1)!)
    // table (PLAN.md).
    Dd<D> epa{zero, zero};
    if (!op::AllFalse(d, m_pa)) {
      const auto mn = op::IfThenElse(m_pa, op::Min(alpha, beta), safe_a);
      const auto mx = op::IfThenElse(m_pa, op::Max(alpha, beta), safe_b);
      const auto lo = op::Ge(op::Set(d, detail::kBetaPrTauMax), mn);
      const auto mns = op::IfThenElse(lo, mn, one);  // scrub dead lanes
      const auto gt1 = op::Gt(mns, one);
      const auto lg1m = LgammaDiffDd(
          d, op::IfThenElse(gt1, op::Add(one, one), one),
          op::IfThenElse(gt1, op::Sub(mns, one), mns));
      const auto lga = DdSub(d, lg1m, BetaLog(d, mns));
      const auto lgp = LgammaPosDd(d, mn);
      const Dd<D> lgmn{op::IfThenElse(lo, lga.hi, lgp.hi),
                       op::IfThenElse(lo, lga.lo, lgp.lo)};
      epa = DdAdd(d, DdAdd(d, l1, l2),
                  DdSub(d, LgammaDiffDd(d, mx, mn), lgmn));
    }
    // PB: the analytic Stirling difference. Binet's derivative is <= 1/(12z^2)
    // here, so c.lo is negligible in Delta and c may even be +inf (Binet -> 0).
    Dd<D> epb{zero, zero};
    if (!op::AllFalse(d, m_pb)) {
      const auto delta =
          op::Sub(op::Add(BinetVal(d, alpha), BinetVal(d, beta)),
                  BinetVal(d, c_raw));
      auto e = DdSub(d, BetaLog(d, ps.nu),
                     Dd<D>{op::Set(d, kBetaLnTwoPiHi),
                           op::Set(d, kBetaLnTwoPiLo)});
      e = Dd<D>{op::Mul(e.hi, half), op::Mul(e.lo, half)};  // exact
      e = DdSub(d, e, ps.cpsi);
      epb = DdAddD(d, e, op::Neg(delta));
    }
    const Dd<D> efull{op::IfThenElse(m_pb, epb.hi, epa.hi),
                      op::IfThenElse(m_pb, epb.lo, epa.lo)};

    // R2 folds (-) log alpha; R1 additionally drops beta*ln y, which its own
    // series form does not carry. |beta*ln y| <= ~11 on R1's box (beta*xi <= 8
    // and xi <= 0.45), so this subtraction costs nothing.
    const auto e2 = DdSub(d, efull, la);
    const auto e1 = DdSub(d, e2, l2);

    if (!op::AllFalse(d, m_r1)) {
      const auto ea = BetaClampE(d, e1.hi, e1.lo);
      const auto as = op::IfThenElse(ea.sat, safe_a, alpha);
      const auto bs = op::IfThenElse(ea.sat, safe_b, beta);
      const auto xs = op::IfThenElse(ea.sat, safe_x, xi.hi);
      const auto ex = ExpDdFrac(d, ea.hi, ea.lo);
      auto r1 = BetaScale(d, DdMul(d, ex.m, BetaR1Series(d, as, bs, xs)), ex.e);
      r1.v = op::IfThenElse(ea.sat, zero, r1.v);
      r1.dd.hi = op::IfThenElse(ea.sat, zero, r1.dd.hi);
      r1.dd.lo = op::IfThenElse(ea.sat, zero, r1.dd.lo);
      val.v = op::IfThenElse(m_r1, r1.v, val.v);
      val.dd.hi = op::IfThenElse(m_r1, r1.dd.hi, val.dd.hi);
      val.dd.lo = op::IfThenElse(m_r1, r1.dd.lo, val.dd.lo);
    }
    if (!op::AllFalse(d, m_r2)) {
      const auto ea = BetaClampE(d, e2.hi, e2.lo);
      const auto as = op::IfThenElse(ea.sat, safe_a, alpha);
      const auto bs = op::IfThenElse(ea.sat, safe_b, beta);
      // c as the EXACT dd pair, per the design's exact-c rule: the CF's
      // (c+m) term is the one place c enters a product, and fl(a+b) alone
      // would put c*2^-53 into every odd d.
      const auto cd = TwoSum(d, alpha, beta);
      const Dd<D> cs{op::IfThenElse(ea.sat, op::Add(safe_a, safe_b), cd.hi),
                     op::IfThenElse(ea.sat, zero, cd.lo)};
      const Dd<D> xs{op::IfThenElse(ea.sat, safe_x, xi.hi),
                     op::IfThenElse(ea.sat, zero, xi.lo)};
      const auto ex = ExpDdFrac(d, ea.hi, ea.lo);
      auto r2 =
          BetaScale(d, DdMul(d, ex.m, BetaR2Cf(d, as, bs, cs, xs)), ex.e);
      r2.v = op::IfThenElse(ea.sat, zero, r2.v);
      r2.dd.hi = op::IfThenElse(ea.sat, zero, r2.dd.hi);
      r2.dd.lo = op::IfThenElse(ea.sat, zero, r2.dd.lo);
      val.v = op::IfThenElse(m_r2, r2.v, val.v);
      val.dd.hi = op::IfThenElse(m_r2, r2.dd.hi, val.dd.hi);
      val.dd.lo = op::IfThenElse(m_r2, r2.dd.lo, val.dd.lo);
    }
  }

  // --- SEVENTH ROUTING CORRECTION: near-one post-route ---------------------
  // An R1 lane whose evaluated dd value exceeds kBetaNearOne (1 - 2^-11)
  // would hand back a complement made of dd rounding noise (the
  // complement-slack doctrine's 1 - 2^-12 bound, with one bit of margin
  // for this compare being on the dd value). Such lanes are R4-SHAPED by
  // construction -- R1's box supplies R4's convergence caps (xi <= xi1,
  // beta*xi <= B1) and near-one puts the Expm1 argument in its ideal
  // zone -- so they fold into the R4 core's lane set below, SAME
  // orientation (alpha = min there: near-one requires the mean below
  // xi <= xi1, i.e. beta > alpha*(1-xi1)/xi1 > alpha -- holds at any
  // gate). The tau ceiling kBetaPrTauMax = 2.5 [EIGHTH correction; was
  // 1.5] is one lgamma recurrence step past the centre-2 edge (see
  // BetaR4Tiny's three-zone lgamma(1+tau)): check (e)'s extended pocket
  // found stay-R1 lanes at tau = 1.6 (b = 20, x = 0.4) whose complement
  // dips BELOW the 2^-12 doctrine bound -- the earlier "safe band"
  // claim was gamma-limit reasoning, wrong at moderate beta. With the
  // gate at 2.5, any doctrine-violating lane also exceeds the
  // kBetaNearOne bar (1 - 2^-11 < 1 - 2^-12) and therefore post-routes;
  // above 2.5 the pocket shows the box-corner complement clears 2^-12
  // with ~3x margin. This replaces the sixth correction's
  // opposite-orientation CF destination, which stalled at 2^-55.5 on
  // the CF's small-second-parameter weakness (generator check
  // (b)(viii), witness (0.0234, 1e6, 4e-6)).
  const auto i_pr = op::Mul(
      i_r1, op::Mul(BetaInd(d, op::Gt(val.dd.hi,
                                      op::Set(d, detail::kBetaNearOne))),
                    BetaInd(d, op::Ge(op::Set(d, detail::kBetaPrTauMax),
                                      alpha))));
  const auto m_pr = BetaIndMask(d, i_pr);
  const auto i_r4x = op::Max(i_r4, i_pr);
  const auto m_r4x = BetaIndMask(d, i_r4x);

  if (!op::AllFalse(d, m_r3)) {
    // The depth-extension indicator is the bmax gate alone (not i_glr's
    // nu part): every R3 slice lane gets the extension rows, including
    // nu >= T_ridge ones -- safe at any nu, see BetaR3Temme's comment.
    auto t3 = BetaR3Temme(
        d, ps, BetaInd(d, op::Ge(bmax, op::Set(d, detail::kBetaGammaLim))));
    t3.val.v = op::IfThenElse(t3.sat, zero, t3.val.v);
    t3.val.dd.hi = op::IfThenElse(t3.sat, zero, t3.val.dd.hi);
    t3.val.dd.lo = op::IfThenElse(t3.sat, zero, t3.val.dd.lo);
    val.v = op::IfThenElse(m_r3, t3.val.v, val.v);
    val.dd.hi = op::IfThenElse(m_r3, t3.val.dd.hi, val.dd.hi);
    val.dd.lo = op::IfThenElse(m_r3, t3.val.dd.lo, val.dd.lo);
    // is_p = XNOR(native, lambda >= 0): the routed triple's I is P when the
    // native orientation was evaluated, and R3 returns that I when lambda >= 0
    // and its complement otherwise.
    const auto lp = t3.lam_pos;
    const auto xn = op::MulAdd(i_nat, lp,
                               op::Mul(BetaIndNot(d, i_nat), BetaIndNot(d, lp)));
    is_p = op::IfThenElse(m_r3, xn, is_p);
  }

  if (!op::AllFalse(d, m_r4x)) {
    // R4's scrub point is its OWN interior (tau <= eps_R4), not the shared
    // (2, 3, 1/4): tau = 2 is outside this core's zones and would take
    // lgamma(1+tau) past the centre-2 edge. The mask is the EXTENDED set
    // m_r4x [SEVENTH correction]: routed R4 lanes plus near-one post-routed
    // R1 lanes -- for the latter, alpha/beta/xi ARE the fired orientation
    // and coincide with the tiny-first frame (alpha = min on every lane the
    // near-one bar can fire on), so lxt is the right log for both kinds.
    const auto t4 =
        op::IfThenElse(m_r4x, alpha, op::Set(d, detail::kBetaEpsR4));
    const auto b4 = op::IfThenElse(m_r4x, beta, safe_b);
    const auto x4 = op::IfThenElse(m_r4x, xi.hi, safe_x);
    // ln(1/4) = -2 ln 2, exact halving/doubling of log_dd's own dd pair, so
    // the scrubbed lane's log stays consistent with its scrubbed xi.
    const Dd<D> l4{
        op::IfThenElse(m_r4x, lxt.hi, op::Set(d, -2.0 * detail::kLogLn2Hi)),
        op::IfThenElse(m_r4x, lxt.lo, op::Set(d, -2.0 * detail::kLogLn2Lo))};
    const auto q4 = BetaR4Tiny(d, t4, b4, x4, l4);
    val.v = op::IfThenElse(m_r4x, DdToDouble(q4), val.v);
    val.dd.hi = op::IfThenElse(m_r4x, q4.hi, val.dd.hi);
    val.dd.lo = op::IfThenElse(m_r4x, q4.lo, val.dd.lo);
    // Post-routed lanes now hold the COMPLEMENT of the orientation they
    // evaluated in R1, so their is_p flips to the R4 convention.
    is_p = op::IfThenElse(m_pr, i_sw, is_p);
  }

  // --- (C) gamma-limit slice: R2's off-band remainder above B_GL -----------
  // (In-band traffic above B_GL went to R3 via the lowered ridge floor
  // kBetaGlRidgeMin.) One routed parameter is >= kBetaGammaLim = 2^59; the
  // beta CF is structurally degenerate there (d1 -> -(1 - tiny); mpmath's
  // own CF divides by zero at working precision -- G3 escalation (C)), so
  // the lane goes through the gamma limit: with s = the small routed
  // parameter and t = -(huge)*(ln of the huge side's own x-argument), dd:
  //   huge SECOND: I_xi(s, huge) ~ P_gamma(s, t),      t = -(huge)ln(1-xi)
  //   huge FIRST : I_xi(huge, s) ~ 1 - P_gamma(s, t),  t = -(huge)ln(xi)
  // relative correction O(1/huge), 2^-49-class at the B_GL pin (generator
  // overlap probe; gated at G4 as its own ULP row -- gammalim).
  // The sub-map mirrors gamma-inl.h's own routing verbatim (series for
  // s < kGammaAT and t <= s+1, or s >= kGammaAT and s >= 2t; CF otherwise;
  // gamma's in-band Temme case cannot occur here, those lanes are R3's).
  // val takes the NATURALLY COMPUTED side (series -> P_gamma, CF ->
  // Q_gamma) -- never a dd complement round-trip, which would hand a small
  // side back as absolute noise -- and is_p records which side of the
  // ORIGINAL pair that is.
  // BOTH-HUGE EXCLUSION (found by the first fresh-reference ULP run): the
  // slice's small/huge mapping is meaningless when BOTH parameters are
  // >= kBetaGammaLim -- s would itself be huge and the "shape" identity
  // breaks (witness (1.02e100, 5e101, 0.01), a band-edge float-fuzz
  // escapee of R3, which came back P = 1 for a true P ~ 0). Such lanes
  // stay in ordinary R2: off-band with nu >= 2^58 their cpsi exceeds the
  // saturation floor astronomically, so the PB prefactor's E-clamp
  // saturates them to the correct side by orientation alone.
  const auto i_gl = op::Mul(
      op::Mul(i_r2,
              BetaInd(d, op::Ge(bmax, op::Set(d, detail::kBetaGammaLim)))),
      BetaInd(d, op::Gt(op::Set(d, detail::kBetaGammaLim),
                        op::Min(alpha, beta))));
  const auto m_gl = BetaIndMask(d, i_gl);
  if (!op::AllFalse(d, m_gl)) {
    const auto hf = op::Ge(alpha, op::Set(d, detail::kBetaGammaLim));
    const auto i_hf = BetaInd(d, hf);
    const auto ss = op::IfThenElse(m_gl, op::IfThenElse(hf, beta, alpha), one);
    const auto huge = op::IfThenElse(hf, alpha, beta);
    const Dd<D> lx_gl =
        BetaLog(d, Dd<D>{op::IfThenElse(hf, xi.hi, yv.hi),
                         op::IfThenElse(hf, xi.lo, yv.lo)});
    // t > 0, dd. Same huge-parameter exact prescale as l1/l2 above: huge
    // exceeds ops::ProdLow's 2^996 non-FMA Dekker ceiling on the corner
    // rows' b = 1e307 gammalim lanes (bit-identical below 2^900).
    const auto gbig = op::Gt(huge, op::Set(d, kBetaScaleAbove));
    const auto gdn = op::IfThenElse(gbig, op::Set(d, kBetaScaleDown), one);
    const auto gup = op::IfThenElse(gbig, op::Set(d, kBetaScaleUp), one);
    const auto tds = DdMulD(d, lx_gl, op::Neg(op::Mul(huge, gdn)));
    const Dd<D> t_dd{op::Mul(tds.hi, gup), op::Mul(tds.lo, gup)};
    // E_g = s*ln t - t - lgamma(s); all the e^-t argument sensitivity is
    // absorbed HERE, in dd -- the cores below only see t.hi, whose 2^-53
    // relative slack enters their series/CF factors with O(1) sensitivity.
    auto e_g = DdMulD(d, BetaLog(d, t_dd), ss);
    e_g = DdSub(d, e_g, t_dd);
    e_g = DdSub(d, e_g, LgammaPosDd(d, ss));
    const auto th = t_dd.hi;
    const auto i_ser = op::Max(
        op::Mul(BetaInd(d, op::Lt(ss, op::Set(d, detail::kGammaAT))),
                BetaInd(d, op::Ge(op::Add(ss, one), th))),
        op::Mul(BetaInd(d, op::Ge(ss, op::Set(d, detail::kGammaAT))),
                BetaInd(d, op::Ge(ss, op::Add(th, th)))));
    const auto m_ser = BetaIndMask(d, i_ser);
    const auto e_s = DdSub(d, e_g, BetaLog(d, ss));  // gamma's R1 fold
    const Dd<D> e_pick{op::IfThenElse(m_ser, e_s.hi, e_g.hi),
                       op::IfThenElse(m_ser, e_s.lo, e_g.lo)};
    const auto ea = BetaClampE(d, e_pick.hi, e_pick.lo);
    const auto s_c = op::IfThenElse(ea.sat, one, ss);
    const auto t_c = op::IfThenElse(ea.sat, op::Set(d, 3.0),
                                    op::Min(th, op::Set(d, 4e306)));
    const auto ex = ExpDdFrac(d, ea.hi, ea.lo);
    const auto ser = GammaSeriesSum(d, s_c, t_c);
    const auto cfr = GammaCfRecip(d, s_c, t_c);
    const Dd<D> fac{op::IfThenElse(m_ser, ser.hi, cfr.hi),
                    op::IfThenElse(m_ser, ser.lo, cfr.lo)};
    auto gl = BetaScale(d, DdMul(d, ex.m, fac), ex.e);
    gl.v = op::IfThenElse(ea.sat, zero, gl.v);
    gl.dd.hi = op::IfThenElse(ea.sat, zero, gl.dd.hi);
    gl.dd.lo = op::IfThenElse(ea.sat, zero, gl.dd.lo);
    val.v = op::IfThenElse(m_gl, gl.v, val.v);
    val.dd.hi = op::IfThenElse(m_gl, gl.dd.hi, val.dd.hi);
    val.dd.lo = op::IfThenElse(m_gl, gl.dd.lo, val.dd.lo);
    // val holds P_gamma on series lanes, Q_gamma on CF lanes. P_gamma is
    // the ROUTED value iff the huge parameter is SECOND; the routed value
    // is P-of-original iff the orientation was native. XNOR twice:
    const auto agree = op::MulAdd(i_ser, BetaIndNot(d, i_hf),
                                  op::Mul(BetaIndNot(d, i_ser), i_hf));
    const auto isp_gl = op::MulAdd(agree, i_nat,
                                   op::Mul(BetaIndNot(d, agree), i_sw));
    is_p = op::IfThenElse(m_gl, isp_gl, is_p);
  }

  // The complement is formed from the dd BEFORE any rounding, so the single
  // rounding of 1 - direct is the only one it carries.
  const auto comp = DdToDouble(
      DdAdd(d, Dd<D>{one, zero}, Dd<D>{op::Neg(val.dd.hi), op::Neg(val.dd.lo)}));
  const auto m_p = BetaIndMask(d, is_p);
  auto res = kP ? op::IfThenElse(m_p, val.v, comp)
                : op::IfThenElse(m_p, comp, val.v);

  // --- specials (PLAN.md's pinned table) ----------------------------------
  // Applied last, in increasing priority. The doctrine is gamma's: one
  // degenerate parameter gets its limit; two degeneracies, or a degenerate
  // parameter meeting the x-boundary its mass sits on, give NaN.
  const auto p_hi = kP ? one : zero;   // "all the mass is at or below x"
  const auto p_lo = kP ? zero : one;   // "no mass at or below x"
  const auto a0 = op::Eq(a_in, zero);
  const auto b0 = op::Eq(b_in, zero);
  const auto ai = op::Eq(a_in, inf);
  const auto bi = op::Eq(b_in, inf);
  const auto x0 = op::Eq(x_in, zero);
  const auto x1 = op::Eq(x_in, one);

  res = op::IfThenElse(x0, p_lo, res);   // x = 0 -> P = 0
  res = op::IfThenElse(x1, p_hi, res);   // x = 1 -> P = 1
  // Mass at 0 (a = 0 or b = +inf): P = 1 on (0, 1]. Mass at 1 (b = 0 or
  // a = +inf): P = 0 on [0, 1).
  res = op::IfThenElse(a0, p_hi, res);
  res = op::IfThenElse(bi, p_hi, res);
  res = op::IfThenElse(b0, p_lo, res);
  res = op::IfThenElse(ai, p_lo, res);
  // A degenerate parameter meeting the boundary its own mass sits on.
  res = op::IfThenElse(op::Gt(op::Mul(BetaInd(d, a0), BetaInd(d, x0)), half),
                       qnan, res);
  res = op::IfThenElse(op::Gt(op::Mul(BetaInd(d, bi), BetaInd(d, x0)), half),
                       qnan, res);
  res = op::IfThenElse(op::Gt(op::Mul(BetaInd(d, b0), BetaInd(d, x1)), half),
                       qnan, res);
  res = op::IfThenElse(op::Gt(op::Mul(BetaInd(d, ai), BetaInd(d, x1)), half),
                       qnan, res);
  // Two degeneracies among {a in {0, inf}, b in {0, inf}}.
  {
    const auto da = op::Max(BetaInd(d, a0), BetaInd(d, ai));
    const auto db = op::Max(BetaInd(d, b0), BetaInd(d, bi));
    res = op::IfThenElse(op::Gt(op::Mul(da, db), half), qnan, res);
  }
  res = op::IfThenElse(op::Lt(a_in, zero), qnan, res);
  res = op::IfThenElse(op::Lt(b_in, zero), qnan, res);
  res = op::IfThenElse(op::Lt(x_in, zero), qnan, res);
  res = op::IfThenElse(op::Gt(x_in, one), qnan, res);
  res = op::IfThenElse(op::IsNaN(a_in), a_in, res);  // payload preserved
  res = op::IfThenElse(op::IsNaN(b_in), b_in, res);
  res = op::IfThenElse(op::IsNaN(x_in), x_in, res);
  return res;
}

// ------------------------------------------------------------------------
// lbeta(a, b) = ln B(a, b) = lgamma(a) + lgamma(b) - lgamma(a+b), on the
// positive-parameter domain (a, b > 0 finite; anything else is NaN --
// SciPy's betaln accepts negatives through |Gamma|, a documented deviation
// with no consumer need). This is the beta TU's third export because it is
// nothing but this TU's own prefactor machinery re-handed: with m = min,
// M = max,
//     lbeta = LgammaPosDd(m) - LgammaDiffDd(M, m),
// the exact pair PA assembles for -ln B, negated, rounded ONCE.
//
// TWO BANDS on m (one masked path, both computed, selected):
//  * m <= 2^990: the assembly above. This extends LgammaDiffDd's discharged
//    call-site bound (m <= 256) to 2^990 on the strength of the mechanism
//    audit, re-run for this caller [lbeta design, PLAN.md]: every
//    m-carrying TwoProd operand (m*ln z, (m-1/2)*log1p(w), the DdDivDd
//    numerator) stays below ops::ProdLow's 2^996 non-FMA Dekker ceiling up
//    to m = 2^990 with 64x headroom; the up-walk never fires (M >= m > Z0
//    whenever m > 256); w = m/z <= 1 sits inside Log1pmxDd's audited
//    domain; LgammaPosDd(m) and the internal m*ln z stay below DBL_MAX
//    (lgamma(2^990) ~ 7e300). Worst cancellation in the final subtraction
//    is ln(M)/ln(2) <= ~2^10 against the dd's ~2^-104 -- budgeted.
//  * m > 2^990: all three lgammas are deep-Stirling and the assembly's
//    terms would overflow before the RESULT does (ln B stays finite to
//    ~ -1.3e308 on the a = b ray). Direct grouped Stirling difference on
//    EXACTLY 2^-64-prescaled operands (the log RATIOS are scale-invariant;
//    the AGENTS.md 2^996 rule again):
//        ln B = -m*ln(c/m) - M*ln(c/M) + (1/2)*ln(c/(m*M)),
//    with c = m + M carried as an exact scaled TwoSum. No cancellation: both
//    leading terms are negative and |ln B| >= m*ln 2 here. The Binet terms
//    Binet(m)+Binet(M)-Binet(c) are DROPPED: each is < 1/(12*2^990), i.e.
//    < 2^-993 absolute against a result of magnitude >= 2^990*ln 2 --
//    relative < 2^-1980, unrepresentable. The (1/2)*ln term (magnitude
//    <= ~710) is likewise ~2^-990 relative but is kept -- it is one cheap
//    log. A result below -DBL_MAX saturates to -inf through the exact
//    final power-of-two scale-back's IEEE rounding, which is the correct
//    value-driven boundary (verified by bit-stepped reference rows).
inline constexpr double kLbetaBigMin = 0x1.0p+990;
inline constexpr double kLbetaScaleDown = 0x1.0p-64;
inline constexpr double kLbetaScaleUp = 0x1.0p+64;

template <class D>
HWY_NOINLINE op::V<D> LbetaVec(D d, op::V<D> a, op::V<D> b) {
  const auto zero = op::Zero(d);
  const auto one = op::Set(d, 1.0);
  const auto nan = op::Set(d, std::numeric_limits<double>::quiet_NaN());

  // Domain: a, b > 0 and finite. Everything else -> NaN (payload of a NaN
  // input is NOT preserved across the min/max scrub; the smoke test pins
  // quiet-NaN-out, matching beta_p's convention).
  const auto maxf = op::Set(d, (std::numeric_limits<double>::max)());
  const auto bad = op::Or(
      op::Or(op::Or(op::Ge(zero, a), op::Ge(zero, b)),
             op::Or(op::IsNaN(a), op::IsNaN(b))),
      op::Or(op::Gt(a, maxf), op::Gt(b, maxf)));

  // Scrub invalid lanes to a benign interior point before any arithmetic.
  auto as = op::IfThenElse(bad, one, a);
  auto bs = op::IfThenElse(bad, one, b);
  const auto m_raw = op::Min(as, bs);
  const auto mm_raw = op::Max(as, bs);

  const auto big = op::Gt(m_raw, op::Set(d, kLbetaBigMin));

  // --- main band (computed on every lane; big lanes scrubbed down) -------
  // ONLY the min is clamped: the scrub exists so big-band lanes stay
  // benign through this dead computation, and m <= 2^990 alone already
  // guarantees that (every m-carrying operand stays in range; see the
  // audit above). M is legitimate at ANY magnitude in the main band --
  // clamping it too corrupted every live lane with max > 2^990 (caught by
  // the reference gate's grid rows, 2026-08-11).
  const auto m1 = op::Min(m_raw, op::Set(d, kLbetaBigMin));
  const auto lgm = LgammaPosDd(d, m1);
  const auto dlg = LgammaDiffDd(d, mm_raw, m1);
  const auto main_dd = DdSub(d, lgm, dlg);

  // --- big band (m > 2^990; small lanes scrubbed up) ---------------------
  const auto m2 = op::Max(m_raw, op::Set(d, kLbetaBigMin));
  const auto mm2 = op::Max(mm_raw, op::Set(d, kLbetaBigMin));
  const auto dn = op::Set(d, kLbetaScaleDown);
  const auto ms = op::Mul(m2, dn);    // exact: operands normal, 2^-64 scale
  const auto mms = op::Mul(mm2, dn);  // exact
  const auto cs = TwoSum(d, mms, ms);  // c' = (m + M)*2^-64, exact dd
  const auto lr1 = BetaLog(d, DdDivDd(d, cs, Dd<D>{ms, zero}));   // ln(c/m)
  const auto lr2 = BetaLog(d, DdDivDd(d, cs, Dd<D>{mms, zero}));  // ln(c/M)
  auto acc = DdAdd(d, DdMulD(d, lr1, ms), DdMulD(d, lr2, mms));
  // Scale back exactly; overflow here IS the -inf boundary (see header).
  const auto up = op::Set(d, kLbetaScaleUp);
  acc = Dd<D>{op::Mul(acc.hi, up), op::Mul(acc.lo, up)};
  // + (1/2) ln(c/(m*M)) = (1/2) ln((c'/m'/M') * 2^-64), formed stepwise so
  // no intermediate overflows; ~2^-990 relative to the result, kept cheap.
  const auto rr = op::Mul(op::Div(op::Div(cs.hi, ms), mms), dn);
  const auto lhalf = BetaLog(d, rr);
  // If the scale-back overflowed, acc.hi is +inf and feeding it onward
  // would put inf-inf NaNs through TwoSum's error algebra -- detect it
  // HERE, hand the subtraction a benign clamped value, and select -inf.
  // No finite acc can overflow inside the subtraction instead: ulp(DBL_MAX)
  // is ~2^971, so subtracting the <=~710-magnitude half-log term cannot
  // push a finite sum across the rounding boundary to -inf.
  const auto inf = op::Set(d, std::numeric_limits<double>::infinity());
  const auto acc_ovf = op::Ge(acc.hi, inf);
  const Dd<D> acc_safe{op::IfThenElse(acc_ovf, one, acc.hi),
                       op::IfThenElse(acc_ovf, zero, acc.lo)};
  const auto big_dd =
      DdSub(d, Dd<D>{op::Mul(lhalf.hi, op::Set(d, 0.5)),
                     op::Mul(lhalf.lo, op::Set(d, 0.5))},
            acc_safe);
  const auto big_val = op::IfThenElse(acc_ovf, op::Neg(inf),
                                      op::Add(big_dd.hi, big_dd.lo));

  // --- combine, one rounding, specials last ------------------------------
  auto res = op::IfThenElse(big, big_val,
                            op::Add(main_dd.hi, main_dd.lo));
  res = op::IfThenElse(bad, nan, res);
  return res;
}

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#endif  // include guard
