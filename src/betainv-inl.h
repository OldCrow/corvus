// Inverse regularized incomplete beta: beta_p_inv(a, b, p) and
// beta_q_inv(a, b, q). Per-target include guard (Highway -inl.h idiom).
//
// Both public functions are ONE pipeline with one bit of orientation, per
// PLAN.md "P1 inverse incomplete beta -- detail design" (BINDING) and the
// parameters replay-pinned in src/betainv_data.h. The forward region cores
// (src/beta-inl.h) are consumed, never re-derived, and so is the gamma
// inverse's seed machinery (src/gammainv-inl.h) for the beta -> gamma limit;
// this file adds the frame, the seed stage, the log-space forward and the
// safeguarded step loop.
//
// THE TRANSFER-BUG RULE (PLAN.md's stage theme for this family). Every
// formula that LOOKS like gammainv's has been re-derived for beta at its
// definition site, because beta's second parameter changes the math and the
// generator stage caught five separate instances of a gamma formula carried
// over unchanged. The three re-derivations that live in THIS file are marked
// [TRANSFER SITE] and each carries its derivation; the load-bearing one is
// the Newton slope, which acquires a factor (1 - y) that gamma has no
// analogue for.
//
// ============================ THE INTERNAL FRAME ============================
// TWO EXACT RELABELLINGS COLLAPSE FOUR PROBLEMS INTO ONE, and everything
// downstream is simpler because of it.
//
//  (1) INPUT-SIDE FLIP. Whatever the caller asked for, the kernel solves
//      against sigma = min(s, 1 - s) <= 1/2. The subtraction is EXACT by
//      Sterbenz there, and the orientation bit records which member of the
//      pair sigma is.
//  (2) OUTPUT ORIENTATION SWAP. The swap identity I_x(a,b) = 1 - I_{1-x}(b,a)
//      relabels (a, b, side) without recomputing any complement, so the kernel
//      may solve for EITHER x or 1 - x and hand the other back with a single
//      subtraction. It solves for whichever of the two is the SMALL one, and
//      returns 1 (-) y when that was 1 - x.
//
// [SELF-CAUGHT, G3] THE SWAP IS NOT DECIDED BY WHICH SIDE sigma IS. The
// obvious frame -- always solve the variable whose probability is sigma <= 1/2
// -- is WRONG, and it is wrong in a way that looks right on paper: sigma = P
// puts y below the MEDIAN, and the median is near 1 whenever beta << alpha.
// The first reference run found it immediately: at (a, b, p) =
// (5.6e-4, 0.813, 0.658) the true x is SUBNORMAL while that frame solves for
// 1 - x = 1.0 and hands back an exact 0; at (1.4e39, 2.3e45, 1/2 + 1 ulp) it
// solves the near-one variable and the final subtraction costs 1.6e-10
// relative, 4.6e5 ULP. The frame must ask which END x is near, and that is a
// question about the DISTRIBUTION, not about the probability side.
//
// THE PROBE IS DEFINITIVE, NOT HEURISTIC: x > 1/2 exactly when
// I_{1/2}(a,b) < P, i.e. when m(1/2) < m_t in the kernel's own logit. So the
// driver evaluates its own forward once at y = 1/2 and reads the sign. It
// costs one forward evaluation out of twelve, needs no new machinery, and is
// exact wherever it matters (a lane the probe misclassifies has its root
// within an ulp of 1/2, where both orientations are equally good).
//
// After both relabellings the problem is
//         solve  I_y(alpha, beta) = tau,   y = min(x, 1-x),
// where tau is sigma or 1 - sigma. BOTH are available exactly -- sigma from
// the Sterbenz flip and 1 - sigma as the low word of an exact TwoSum -- so
// nothing is lost by tau landing on the large side, and the whole downstream
// stage generalizes by ONE sign carried on z = erfcinv(2 sigma):
//   * The logit objective m = ln P - ln Q is antisymmetric under each
//     relabelling, so one code path serves both exports and both orientations
//     and the target logit is just +-(ln sigma - ln(1 - sigma)).
//   * S1's ridge variable is zeta = erfcinv(2 tau)/sqrt(nu), and
//     erfcinv(2(1-sigma)) = -erfcinv(2 sigma) exactly, so the sign is carried
//     rather than a second inverse being evaluated. Likewise S5's normal
//     quantile.
//   * The deep-small closed form has ONE branch. The generator's q-side twin
//     y = 1 - exp((ln sigma + ln b + lnB(b,a))/b) IS this file's form after
//     the relabelling (alpha = b), including its cut -- the generator's
//     |1-a|*(1-y)/(1+b)*corr(1-y) is literally |1-beta|*y/(1+alpha)*corr(y)
//     here. The THIRD-correction disease (a q-side branch that no self-check
//     ever swept) cannot recur in a kernel that has no q-side branch.
//   * No complement of a rounded near-one value is ever formed. That is the
//     whole advantage of an inverse over the forward pair, whose complement
//     rounding is on the OUTPUT and cannot be undone.
// The caller's own escape hatch for a near-one answer is the same identity,
// documented in the public header: 1 - x at full relative precision is
// beta_p_inv(b, a, q).
//
// ========================= EVERYTHING IS A LOGIT =========================
// Newton runs on
//     m(y) = ln P(y) - ln Q(y)   against   m_t = ln sigma - ln(1 - sigma),
// never on F(y) - sigma and never on ln F alone -- gammainv's four reasons
// carry over verbatim (ln sigma reaches -745 so the difference must be dd;
// F itself underflows over a huge part of the domain while ln F stays an
// ordinary number; the solved side tends to 1 over the wrong half of the
// domain so its log saturates; and |ln min(P,Q)| jumps by 2 ln 2 exactly at
// the median, which is where a target of 1/2 puts its root).
//
// [TRANSFER SITE] THE SLOPE. dm/dy = f(y)/(P Q) with f the beta density, so
// the reciprocal the multiplicative step wants is
//     w = d(ln y)/dm = P Q / (y f(y)).
// With E = alpha ln y + beta ln(1-y) - ln B(alpha,beta) the density satisfies
//     y f(y) = y^alpha (1-y)^(beta-1) / B = e^E / (1 - y),
// hence
//     w = exp(ln u - E) * (1 - u) * (1 - y),   u = min(P, Q).
// GAMMA'S w HAS NO (1 - y) FACTOR: there x g(x) = x^a e^-x / Gamma(a) = e^E
// exactly, with no leftover. Carrying gamma's spelling over would under-step
// by (1 - y) at every iteration -- the same class of error the generator
// caught in its own m_and_w_mp (a missing 1/y factor there). Checked against
// the uniform case alpha = beta = 1, where m = ln y - ln(1-y) gives
// w = 1 - y exactly and the formula reproduces it.
//
// THE FORWARD IS beta's OWN ROUTER, IN LOG SPACE. BetaInvForward reproduces
// route_final and the four region cores of src/beta-inl.h, but returns
// ln(one side) + which side it is, never a probability:
//   * R1 and R2 never exponentiate at all -- ln I = E (-) ln alpha (-)
//     beta ln(1-xi) (+) ln(series) and E (-) ln alpha (+) ln(CF) are already
//     logs, so the whole underflow range stays live and beta's saturation
//     clamp (kBetaExpFloor) is simply absent.
//   * R3's tail branch is likewise exact in log space: its value is
//     e^-cpsi * bracket, so ln = -cpsi + ln(bracket).
//   * R4 and the gamma-limit slice return the genuinely small side, so their
//     log is taken directly.
// E is ORIENTATION-INVARIANT (swapping (alpha,beta,y) -> (beta,alpha,1-y)
// maps E to itself), which is why the shared cpsi machinery may be called in
// the kernel's own frame while R1/R2/R4 use the router's -- and why R3's
// lam_pos is directly this frame's "is the value P", with none of the
// forward's XNOR bookkeeping.
//
// SEEDS: FIVE FAMILIES, ALL GLOBAL, ONE RESIDUAL COMPARISON. Each candidate
// is evaluated wherever its own availability/stability gate passes and the
// one with the smallest |m(y0) - m_t| wins -- the same objective the steps
// then reduce, so seed choice and refinement cannot disagree about what
// "closer" means. The winner's forward evaluation is REUSED as the first
// step's. The five are S1 (beta-Temme ridge), S2 (small-y series inversion
// with Picard corrections), S3 (gamma-limit transfer, seeded by gammainv's
// OWN S1/S2/S3), S4 (exact-B leading-order closed form, either branch by
// sigma vs beta/(alpha+beta)) and S5 (logit-normal from exact digamma /
// trigamma moments). PLAN's FIRST correction is why all five are global:
// gating S4 on t_jt was the deviation that opened the moderate-tiny gap band.
//
// DEEP-SMALL CLOSED FORM. Where the dropped series factor S' is below the
// resolution of a double the equation collapses to
//     y = exp((ln sigma + ln alpha + lnB)/alpha),
// evaluated through ExpDdFrac's mantissa/exponent split so the power-of-two
// scaling is the LAST operation and a subnormal (or zero) answer carries
// exactly one rounding. The cut is the THIRD correction's, NOT any
// parameter*y form: what is dropped is |ln S'|/alpha whose leading term is
// |1-beta| y/(1+alpha) -- the OTHER side's parameter is the coefficient --
// times the exact closed-form multiplier corr(y) = -ln(1-y)/y that makes the
// bound sound at the widened gamma-limit corner. See src/betainv_data.h.
//
// ln sigma AND ln(1 - sigma), NEVER log of 1 (-) tiny. 1 - sigma is formed as
// an EXACT dd (TwoSum) and LogDdAny takes the pair: below sigma ~ 2^-53 the
// low word is the entire signal, and a native-float 1 - sigma would collapse
// to exactly 1.0 (the G2 completion round's own witness).
#if defined(CORVUS_BETAINV_INL_H_) == defined(HWY_TARGET_TOGGLE)
#ifdef CORVUS_BETAINV_INL_H_
#undef CORVUS_BETAINV_INL_H_
#else
#define CORVUS_BETAINV_INL_H_
#endif

#include <limits>

#include "src/beta-inl.h"
#include "src/beta_data.h"
#include "src/betainv_data.h"
#include "src/dd-inl.h"
#include "src/dd_special-inl.h"
#include "src/digamma-inl.h"
#include "src/erfinv-inl.h"
#include "src/exp_dd-inl.h"
#include "src/gamma-inl.h"
#include "src/gamma_data.h"
// gammainv-inl.h: ONLY for the three GammaInvSeed* template functions the S3
// gamma-limit transfer needs (cross-family include, the ErfcinvVec
// precedent). Templates instantiate what is called, so gammainv's own
// forward, its four region cores and its driver do NOT enter this TU.
#include "src/gammainv-inl.h"
#include "src/gammainv_data.h"
#include "src/lgamma-inl.h"
#include "src/log_dd-inl.h"
#include "src/ops-inl.h"
#include "src/trigamma-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

namespace op = ops;

// --- kernel-internal constants -------------------------------------------
// Structural (loop shapes, scrub points, finiteness sentinels), not fitted
// table data, so they live here rather than in the generated
// src/betainv_data.h -- the split erfinv, gamma and beta all make.

// Benign interior iterate every core is scrubbed to on a lane whose own
// candidate was rejected. Any point in (0, 1) does; 1/4 is beta's own
// kBetaSafeX and keeps the two files' scrub conventions identical.
inline constexpr double kBetaInvSafeY = 0.25;

// Newton iterations for the zeta -> lambda inversion inside S1, and the
// relative inset that keeps lambda strictly inside its bracket (-beta, alpha).
//
// [G3 LATITUDE TAKEN, PLAN's own flag] The generator needed niter = 100 and
// recorded that as a wart for this stage. The cause is not the arithmetic, it
// is the FORMULATION: it ran Newton on cpsi(lambda) - target, and cpsi has a
// QUADRATIC TANGENCY at lambda = 0 (cpsi ~ lambda^2/(2 nu), d cpsi/d lambda
// -> 0), so Newton there is not quadratic at all -- it halves the error per
// step, which is exactly the "99.99 -> 66.40 -> 88.53 -> converged" crawl the
// generator observed. Running the SAME Newton on
//     F(lambda) = sign(lambda) sqrt(cpsi(lambda)/nu)
// removes the tangency on paper (F ~ lambda/(nu sqrt 2) near zero, slope
// bounded away from zero) and restores quadratic convergence; eight steps from
// the leading-order start lambda_0 = sqrt(2) nu zeta cover the fitted zeta
// domain with room to spare, and a bisection safeguard makes the count an
// upper bound rather than an assumption. See BetaInvLamOfZeta.
inline constexpr int kBetaInvS1LamIters = 8;
inline constexpr double kBetaInvLamInset = 0x1.0p-40;

// Cap and freeze threshold for S2's plain-double replay of R1's series. The
// generator's own s2_series_S_double uses 80 terms and a 1e-18 relative
// freeze; this is a SEED, and its own rounding sits far below the seed floors
// the generator measured.
inline constexpr int kBetaInvS2SeriesN = 80;
inline constexpr double kBetaInvSeedEps = 0x1.0p-60;

// S3 availability gate: "neither parameter is usably huge" is the generator's
// own stated domain guard for the gamma-limit transfer, made explicit here so
// the three gamma-seed candidates (and the gammainv machinery they pull in)
// are skipped whole-vector when no lane can use them. The beta -> gamma limit
// carries a relative correction O(1/max(alpha,beta)), so at 2^20 the transfer
// is already a 10^-6 approximation and every other candidate outscores it;
// the gate is two orders below any point the gamma-limit stratum contains.
inline constexpr double kBetaInvS3Min = 0x1.0p+20;

// Sentinel replacing a non-finite log anywhere in the forward. Only lanes
// astronomically far from any root can reach one (alpha*ln y overflows to
// -inf for alpha past ~2.5e305 with |ln y| > 2, and lgamma(min) overflows
// past 2.556e305): such a lane must not poison its neighbours, must LOSE the
// seed comparison, and must freeze rather than step. A huge finite log does
// all three -- the residual it produces is enormous, and the w clamp below
// bounds the step it could ask for.
inline constexpr double kBetaInvLogSentinel = -0x1.0p+60;

// Clamp on the exponent of w = exp(ln u (-) E) * (1-u) * (1-y). exp saturates
// at +-745; clamping the argument at 700 keeps w finite for every input,
// including the ones whose E or ln u hit the sentinel above. |w| is the
// inverse's own condition number (|d ln y / d ln sigma| ~ 1/alpha in the
// power-law regime), and every input whose true y is a normal double has it
// far inside this bound -- the deep-small closed form owns the rest.
inline constexpr double kBetaInvWExpMax = 700.0;

// Residual reported for a lane whose forward could not be evaluated. Large
// enough to lose every comparison, small enough to keep the dd arithmetic in
// range.
inline constexpr double kBetaInvBigResid = 0x1.0p+100;

// ln 2, for the "is the computed side the smaller one" test.
inline constexpr double kBetaInvNegLn2 = -0.6931471805599453;

// ln(kBetaNearOne) = ln(1 - 2^-11): beta's near-one post-route bar, expressed
// in the log space this forward works in. Written out rather than derived at
// run time; the value is -(x + x^2/2 + x^3/3 + ...) at x = 2^-11, and the bar
// has three orders of margin either side of it, so the last digits are not
// load bearing.
inline constexpr double kBetaInvNearOneLog = -4.8840049810887e-4;

// sqrt(2), for S1's correction prefactor and S5's normal quantile.
inline constexpr double kBetaInvSqrt2 = 0x1.6a09e667f3bcdp+0;

// The two gates that keep E on its exact form (PA) wherever that form still
// has digits to spare -- see the prefactor block in BetaInvForward, which
// derives both. kBetaInvPaScaleMax is where PA's dd assembly runs out
// (2^40 * 2^-105 = 2^-65 absolute, below anything else here); kBetaInvPsiTMin
// is how close y may come to either end, as a fraction of the mean, before
// PB's cpsi loses its own accuracy to Log1pmxDd's unnormalized pair (t^3/3
// with t = lo/hi; at 2^-30 that term is 2^-92 before the alpha multiplier).
inline constexpr double kBetaInvPaScaleMax = 0x1.0p+40;
inline constexpr double kBetaInvPsiTMin = 0x1.0p-30;

// Dead band on the orientation probe: how far from 1/2 the root must be
// before the swap is worth taking. Inside it BOTH orientations give a value
// within 2^-10 of 1/2, where 1 (-) y is exact by Sterbenz -- so the tie goes
// to the cheaper orientation rather than to the forward's last bits. See the
// probe in BetaInvVec for the witness that made this necessary.
inline constexpr double kBetaInvSwapBand = 0x1.0p-10;

// ------------------------------------------------------------------------
// Small shared helpers.

// OUTLINED log, exp and sigmoid [MSVC BUILD-TIME GATE, AGENTS.md]. These are
// thin wrappers whose only purpose is the HWY_NOINLINE: log_dd, exp_dd and
// Highway's own Exp are each large, and this file reaches them from about
// twenty call sites -- one per region, per prefactor path, per seed and per
// step. Inlined, each of those becomes its own copy of a table gather plus a
// polynomial IN EVERY ONE OF THE COMPILED TARGETS, and cl.exe's optimizer is
// superlinear in function size: the un-outlined form of this TU ran past 45
// minutes and 7 GB on one MSVC invocation (see the final report), against a
// gate of ~4. Bit-identity is guaranteed by contraction-off and verified by
// byte-comparing the ULP tables across the change.
template <class D>
HWY_NOINLINE Dd<D> BetaInvLog(D d, Dd<D> x) {
  return LogDdAny(d, x);
}
template <class D>
HWY_NOINLINE Dd<D> BetaInvLog(D d, op::V<D> x) {
  return BetaInvLog(d, Dd<D>{x, op::Zero(d)});
}
template <class D>
HWY_NOINLINE Dd<D> BetaInvExpDd(D d, op::V<D> xh, op::V<D> xl) {
  return ExpDd(d, xh, xl);
}
template <class D>
HWY_NOINLINE op::V<D> BetaInvExp(D d, op::V<D> x) {
  return op::Exp(d, x);
}

// Replace a non-finite dd log by a finite sentinel so every downstream
// operation stays total. NaN fails |v| < inf as well, so one test covers both.
template <class D>
HWY_INLINE Dd<D> BetaInvFiniteLog(D d, Dd<D> v) {
  const auto ok = op::Lt(op::Abs(v.hi),
                         op::Set(d, std::numeric_limits<double>::infinity()));
  return Dd<D>{op::IfThenElse(ok, v.hi, op::Set(d, kBetaInvLogSentinel)),
               op::IfThenElse(ok, v.lo, op::Zero(d))};
}

// 1/(1 + e^-v), evaluated on the side that does not cancel: for v >= 0 the
// result is 1 (-) sigmoid(-v) with sigmoid(-v) <= 1/2, so the subtraction is
// a single rounding of a value >= 1/2 and never loses a bit that matters.
// Overflow-free by construction: only e^-|v| is ever formed.
template <class D>
HWY_NOINLINE op::V<D> BetaInvSigmoid(D d, op::V<D> v) {
  const auto one = op::Set(d, 1.0);
  const auto e = BetaInvExp(d, op::Neg(op::Abs(v)));
  const auto lo = op::Div(e, op::Add(one, e));  // sigmoid(-|v|) in (0, 1/2]
  return op::IfThenElse(op::Lt(v, op::Zero(d)), lo, op::Sub(one, lo));
}

// phi(w)/w^2 with phi(w) = w - log1p(w), in plain double, relatively accurate
// through w = 0 (where it is 1/2). Log1pmxDd is the primitive that is accurate
// RELATIVE to phi, which is exactly what makes the quotient safe; forming
// w - log1p(w) directly would cancel to nothing for small w. Seed-grade: only
// the high word is kept.
template <class D>
HWY_INLINE op::V<D> BetaInvPhiOverSq(D d, op::V<D> w) {
  const auto half = op::Set(d, 0.5);
  const auto ww = op::Mul(w, w);
  const auto phi = Log1pmxDd(d, Dd<D>{w, op::Zero(d)}).hi;
  // w == 0 is the limit 1/2; the select also covers a scrubbed lane's zero.
  return op::IfThenElse(op::Eq(ww, op::Zero(d)), half, op::Div(phi, ww));
}

// lgamma(1 + m) for m > 0, as a dd. Below kBetaPrTauMax this is beta's own
// NINTH-correction identity (LgammaDiffDd(1, m) or LgammaDiffDd(2, m-1),
// exact because lgamma(1) = lgamma(2) = 0), which stays relative-accurate as
// m -> 0 where fl(1 + m) rounds to exactly 1 and LgammaPosDd would return a
// flat zero. Above it lgamma(m) and ln m are both positive and neither is
// small, so lgamma(m) (+) ln m has nothing to cancel.
template <class D>
HWY_NOINLINE Dd<D> BetaInvLgamma1p(D d, op::V<D> m) {
  const auto one = op::Set(d, 1.0);
  const auto lo = op::Ge(op::Set(d, detail::kBetaPrTauMax), m);
  const auto ms = op::IfThenElse(lo, m, one);  // scrub the dead branch
  const auto g1 = op::Gt(ms, one);
  const auto diff = LgammaDiffDd(d, op::IfThenElse(g1, op::Add(one, one), one),
                                 op::IfThenElse(g1, op::Sub(ms, one), ms));
  const auto mh = op::IfThenElse(lo, op::Set(d, 3.0), m);  // scrub the other
  const auto gen = DdAdd(d, LgammaPosDd(d, mh), BetaInvLog(d, mh));
  return Dd<D>{op::IfThenElse(lo, diff.hi, gen.hi),
               op::IfThenElse(lo, diff.lo, gen.lo)};
}

// ------------------------------------------------------------------------
// Per-call context: everything the forward and the seeds need that does not
// depend on the iterate y. Computed ONCE per vector, then reused by up to
// eleven forward evaluations (seven seed candidates plus four steps).
//
// AGGREGATE-INITIALIZED AT EVERY CONSTRUCTION SITE, never default-constructed
// (AGENTS.md): an implicitly generated default constructor for a struct with
// vector members is instantiated OUTSIDE the per-target attribute region, and
// the NEON_BF16 slice then refuses to inline Vec128's always_inline
// constructors into it -- the macOS CI break of 2026-08-06, twice.
template <class D>
struct BetaInvCtx {
  op::V<D> alpha;
  op::V<D> beta;
  op::V<D> c_raw;   // alpha + beta, may be +inf; only ever compared or Binet'd
  Dd<D> mlnb;       // -ln B(alpha, beta)
  Dd<D> lb_a;       // ln alpha (+) ln B  = lgamma(1+alpha)+lgamma(beta)-lgamma(c)
  Dd<D> lb_b;       // ln beta  (+) ln B
  Dd<D> ln_alpha;   // ln alpha, ln beta: the forward's (-) ln alpha fold picks
  Dd<D> ln_beta;    //   whichever the router's orientation made first
  Dd<D> pbk;        // 1/2 ln(nu/2pi) (-) Delta: the PB prefactor's y-free part
  op::V<D> i_pb;    // 1.0 where E takes the Stirling-difference path
  op::V<D> nu;      // alpha*beta/c, unscaled (finite: nu <= min(alpha,beta))
  op::V<D> p_mean;  // alpha/c, on the exact prescale (finite for any pair)
  op::V<D> q_mean;  // beta/c
};

// The exact power-of-two prescale beta's own cpsi machinery uses, applied
// here for the same two reasons: c = alpha + beta is not representable for
// both parameters near DBL_MAX, and a parameter past 2^996 breaks
// ops::ProdLow's non-FMA Dekker split. p, q and nu/scale are 0- and
// 1-homogeneous respectively, so one multiply takes nu back.
template <class D>
HWY_NOINLINE BetaInvCtx<D> BetaInvPrepare(D d, op::V<D> alpha, op::V<D> beta) {
  const auto one = op::Set(d, 1.0);
  const auto half = op::Set(d, 0.5);

  const auto big = op::Gt(op::Max(alpha, beta), op::Set(d, kBetaScaleAbove));
  const auto dn = op::IfThenElse(big, op::Set(d, kBetaScaleDown), one);
  const auto up = op::IfThenElse(big, op::Set(d, kBetaScaleUp), one);
  const auto as = op::Mul(alpha, dn);
  const auto bs = op::Mul(beta, dn);
  const auto cs = op::Add(as, bs);  // finite by construction
  // nu = alpha*(beta/c): beta/c <= 1, so this ordering cannot overflow where
  // alpha*beta would, and nu <= min(alpha, beta) is always representable.
  const auto nu = op::Mul(op::Mul(as, op::Div(bs, cs)), up);

  // -ln B, via the analytic lgamma difference in BOTH terms, so the rounded
  // c = alpha + beta is never an lgamma argument (beta-inl.h's own PA rule,
  // and the reason lnB stays finite up to min(alpha,beta) ~ 2.5e305):
  //     ln B = lgamma(mn) (-) [lgamma(mx+mn) (-) lgamma(mx)]
  //          = lgamma(1+mn) (-) ln mn (-) LgammaDiffDd(mx, mn).
  // The two "ln alpha + ln B" forms the seeds and the closed form need follow
  // from it without ever forming ln alpha where alpha is tiny:
  //     LB(mn) = lgamma(1+mn) (-) D          [the +ln mn cancels exactly]
  //     LB(mx) = LB(mn) (+) ln mx (-) ln mn.
  const auto mn = op::Min(alpha, beta);
  const auto mx = op::Max(alpha, beta);
  const auto dlg = LgammaDiffDd(d, mx, mn);
  const auto lg1mn = BetaInvLgamma1p(d, mn);
  const auto lb_mn = DdSub(d, lg1mn, dlg);
  const auto lnmn = BetaInvLog(d, mn);
  const auto lnmx = BetaInvLog(d, mx);
  const auto lnb = DdSub(d, lb_mn, lnmn);
  const auto lb_mx = DdAdd(d, lb_mn, DdSub(d, lnmx, lnmn));
  const auto amin = op::Ge(beta, alpha);  // alpha is the min
  const Dd<D> lb_a{op::IfThenElse(amin, lb_mn.hi, lb_mx.hi),
                   op::IfThenElse(amin, lb_mn.lo, lb_mx.lo)};
  const Dd<D> lb_b{op::IfThenElse(amin, lb_mx.hi, lb_mn.hi),
                   op::IfThenElse(amin, lb_mx.lo, lb_mn.lo)};
  const Dd<D> ln_alpha{op::IfThenElse(amin, lnmn.hi, lnmx.hi),
                       op::IfThenElse(amin, lnmn.lo, lnmx.lo)};
  const Dd<D> ln_beta{op::IfThenElse(amin, lnmx.hi, lnmn.hi),
                      op::IfThenElse(amin, lnmx.lo, lnmn.lo)};

  // PB, the analytic Stirling difference, on beta's own gate. Binet's
  // arguments are clamped up to Z0 because this is computed unconditionally:
  // a subnormal parameter would send 1/z to infinity and NaN the whole pair,
  // and those lanes are selected away rather than masked off.
  const auto c_raw = op::Add(alpha, beta);
  const auto z0 = op::Set(d, detail::kBetaZ0);
  const auto delta = op::Sub(op::Add(BinetVal(d, op::Max(alpha, z0)),
                                     BinetVal(d, op::Max(beta, z0))),
                             BinetVal(d, op::Max(c_raw, z0)));
  auto pbk = DdSub(d, BetaInvLog(d, nu),
                   Dd<D>{op::Set(d, kBetaLnTwoPiHi),
                         op::Set(d, kBetaLnTwoPiLo)});
  pbk = Dd<D>{op::Mul(pbk.hi, half), op::Mul(pbk.lo, half)};  // exact
  pbk = DdAddD(d, pbk, op::Neg(delta));
  const auto i_pb =
      op::Mul(BetaInd(d, op::Gt(c_raw, op::Set(d, detail::kBetaClg))),
              BetaInd(d, op::Ge(mn, z0)));

  const BetaInvCtx<D> out{alpha, beta,     c_raw,          DdNeg(d, lnb),
                          lb_a,  lb_b,     ln_alpha,       ln_beta,
                          pbk,   i_pb,     nu,
                          op::Div(as, cs), op::Div(bs, cs)};
  return out;
}

// The same context with alpha and beta exchanged. Only five fields move: c,
// -ln B, nu, the PB constant and the PB gate are all SYMMETRIC in the two
// parameters (Delta and 1/2 ln(nu/2pi) manifestly so), which is why the
// orientation probe can run against one context and the winner be selected
// field-wise instead of the whole preparation being repeated.
template <class D>
HWY_INLINE BetaInvCtx<D> BetaInvSwapCtx(D, const BetaInvCtx<D>& cx,
                                        op::M<D> m) {
  const BetaInvCtx<D> out{
      op::IfThenElse(m, cx.beta, cx.alpha),
      op::IfThenElse(m, cx.alpha, cx.beta),
      cx.c_raw,
      cx.mlnb,
      Dd<D>{op::IfThenElse(m, cx.lb_b.hi, cx.lb_a.hi),
            op::IfThenElse(m, cx.lb_b.lo, cx.lb_a.lo)},
      Dd<D>{op::IfThenElse(m, cx.lb_a.hi, cx.lb_b.hi),
            op::IfThenElse(m, cx.lb_a.lo, cx.lb_b.lo)},
      Dd<D>{op::IfThenElse(m, cx.ln_beta.hi, cx.ln_alpha.hi),
            op::IfThenElse(m, cx.ln_beta.lo, cx.ln_alpha.lo)},
      Dd<D>{op::IfThenElse(m, cx.ln_alpha.hi, cx.ln_beta.hi),
            op::IfThenElse(m, cx.ln_alpha.lo, cx.ln_beta.lo)},
      cx.pbk,
      cx.i_pb,
      cx.nu,
      op::IfThenElse(m, cx.q_mean, cx.p_mean),
      op::IfThenElse(m, cx.p_mean, cx.q_mean)};
  return out;
}

// ------------------------------------------------------------------------
// The forward, in log space. Returns the kernel's OBJECTIVE and its slope --
// never a probability, which is what keeps the whole underflow range live.
template <class D>
struct BetaInvFwdOut {
  Dd<D> m;       // ln P(y) - ln Q(y) for the kernel's own (alpha, beta)
  op::V<D> w;    // d(ln y)/dm at this point
  op::V<D> unc;  // the residual's own noise floor (see the assembly below)
};

// HWY_NOINLINE from day one, like the region cores it calls: this function is
// reached from up to eleven call sites per export and inlines beta's four
// region cores, the shared cpsi machinery, LgammaDiffDd, Log1pmxDd, Expm1Dd,
// exp_dd, log_dd and both gamma cores. MSVC's optimizer is superlinear in
// function size (AGENTS.md).
template <class D>
HWY_NOINLINE BetaInvFwdOut<D> BetaInvForward(D d, const BetaInvCtx<D>& cx,
                                             op::V<D> y) {
  const auto one = op::Set(d, 1.0);
  const auto zero = op::Zero(d);
  const auto half = op::Set(d, 0.5);
  const auto safe_a = op::Set(d, kBetaSafeA);
  const auto safe_b = op::Set(d, kBetaSafeB);
  const auto safe_x = op::Set(d, kBetaSafeX);

  const auto alpha = cx.alpha;
  const auto beta = cx.beta;
  const auto yc = TwoSum(d, one, op::Neg(y));  // 1 - y, EXACT as a dd pair
  const auto lny = BetaInvLog(d, Dd<D>{y, zero});
  const auto lnyc = BetaInvLog(d, yc);

  // --- routing: beta's route_final with (a, b, x) = (alpha, beta, y) --------
  // Every predicate is written so an overflowing intermediate still decides
  // correctly (beta-inl.h's own note): b*x and B*xi_tau go to +inf and fail
  // their <= B1 test, nu is the harmonic form, the ratio band is xi*(1+b/a).
  const auto tau = op::Min(alpha, beta);
  const auto bmax = op::Max(alpha, beta);
  const auto sw4 = op::Lt(beta, alpha);  // tiny-first
  const Dd<D> xt{op::IfThenElse(sw4, yc.hi, y),
                 op::IfThenElse(sw4, yc.lo, zero)};
  // beta-inl.h takes a fresh LogDdAny here; both logs are already in hand.
  const Dd<D> lxt{op::IfThenElse(sw4, lnyc.hi, lny.hi),
                  op::IfThenElse(sw4, lnyc.lo, lny.lo)};

  const auto xi1 = op::Set(d, detail::kBetaXi1);
  const auto b1v = op::Set(d, detail::kBetaB1);
  const auto thr_t = op::Div(
      one, op::Add(one, op::Div(op::Add(bmax, one), op::Add(tau, one))));
  auto i_r4 = op::Mul(
      BetaInd(d, op::Ge(op::Set(d, detail::kBetaEpsR4), tau)),
      op::Mul(BetaInd(d, op::Ge(op::Set(d, detail::kBetaLn2),
                                op::Mul(tau, op::Abs(lxt.hi)))),
              op::Mul(op::Max(BetaInd(d, op::Ge(xi1, xt.hi)),
                              BetaInd(d, op::Lt(xt.hi, thr_t))),
                      BetaInd(d, op::Ge(b1v, op::Mul(bmax, xt.hi))))));

  const auto r1n = op::Mul(BetaInd(d, op::Ge(xi1, y)),
                           BetaInd(d, op::Ge(b1v, op::Mul(beta, y))));
  const auto r1s = op::Mul(BetaInd(d, op::Ge(xi1, yc.hi)),
                           BetaInd(d, op::Ge(b1v, op::Mul(alpha, yc.hi))));
  const auto i_r1 = op::Mul(op::Max(r1n, r1s), BetaIndNot(d, i_r4));

  const auto nu_r =
      op::Div(one, op::Add(op::Div(one, alpha), op::Div(one, beta)));
  const auto rat1 = op::Mul(y, op::Add(one, op::Div(beta, alpha)));      // xi/p
  const auto rat2 = op::Mul(yc.hi, op::Add(one, op::Div(alpha, beta)));  // y/q
  const auto lo_r = op::Set(d, detail::kBetaXiRatioLo);
  const auto hi_r = op::Set(d, detail::kBetaXiRatioHi);
  const auto band = op::Mul(
      op::Mul(BetaInd(d, op::Ge(rat1, lo_r)), BetaInd(d, op::Ge(hi_r, rat1))),
      op::Mul(BetaInd(d, op::Ge(rat2, lo_r)), BetaInd(d, op::Ge(hi_r, rat2))));
  const auto i_glr = op::Mul(
      BetaInd(d, op::Ge(bmax, op::Set(d, detail::kBetaGammaLim))),
      BetaInd(d, op::Ge(nu_r, op::Set(d, detail::kBetaGlRidgeMin))));
  const auto i_r3 = op::Mul(
      op::Mul(op::Max(BetaInd(d, op::Ge(nu_r, op::Set(d, detail::kBetaTRidge))),
                      i_glr),
              band),
      op::Mul(BetaIndNot(d, i_r4), BetaIndNot(d, i_r1)));

  const auto thr = op::Div(
      one, op::Add(one, op::Div(op::Add(beta, one), op::Add(alpha, one))));
  const auto i_r2 = op::Mul(BetaIndNot(d, i_r4),
                            op::Mul(BetaIndNot(d, i_r1), BetaIndNot(d, i_r3)));

  const auto i_sw = op::MulAdd(
      i_r4, BetaInd(d, sw4),
      op::MulAdd(i_r1, op::Mul(BetaIndNot(d, r1n), r1s),
                 op::MulAdd(i_r3, BetaInd(d, op::Gt(rat1, one)),
                            op::Mul(i_r2, BetaIndNot(
                                d, BetaInd(d, op::Lt(y, thr)))))));
  const auto m_sw = BetaIndMask(d, i_sw);

  const auto ra = op::IfThenElse(m_sw, beta, alpha);
  const auto rb = op::IfThenElse(m_sw, alpha, beta);
  const Dd<D> rxi{op::IfThenElse(m_sw, yc.hi, y),
                  op::IfThenElse(m_sw, yc.lo, zero)};
  const Dd<D> lrxi{op::IfThenElse(m_sw, lnyc.hi, lny.hi),
                   op::IfThenElse(m_sw, lnyc.lo, lny.lo)};
  const Dd<D> lryv{op::IfThenElse(m_sw, lny.hi, lnyc.hi),
                   op::IfThenElse(m_sw, lny.lo, lnyc.lo)};

  const auto m_r1 = BetaIndMask(d, i_r1);
  const auto m_r2 = BetaIndMask(d, i_r2);
  const auto m_r3 = BetaIndMask(d, i_r3);
  const auto m_r4 = BetaIndMask(d, i_r4);
  const auto i_nat = BetaIndNot(d, i_sw);
  auto is_p = op::IfThenElse(m_r4, i_sw, i_nat);

  // --- shared cpsi machinery, in THIS frame ---------------------------------
  // E and cpsi are both invariant under the orientation swap (u and v trade
  // places and alpha u + beta v = 0 identically), so the core may be called on
  // (alpha, beta, y, 1-y) whatever the router chose -- and R3's lam_pos is
  // then directly this frame's "the value is P", with none of beta's XNOR.
  // Scrubbed outside its own lanes: the prescale's no-underflow argument rests
  // on min(alpha, beta) >= Z0, which only holds where this is live.
  const auto m_psi = BetaIndMask(d, op::Max(i_r3, cx.i_pb));
  BetaPsi<D> ps{Dd<D>{zero, zero}, Dd<D>{one, zero}, zero, half, half, zero};
  if (!op::AllFalse(d, m_psi)) {
    const auto pa = op::IfThenElse(m_psi, alpha, safe_a);
    const auto pb = op::IfThenElse(m_psi, beta, safe_b);
    const Dd<D> px{op::IfThenElse(m_psi, y, safe_x), zero};
    const Dd<D> py{op::IfThenElse(m_psi, yc.hi, op::Set(d, 1.0 - kBetaSafeX)),
                   op::IfThenElse(m_psi, yc.lo, zero)};
    ps = BetaPsiCore(d, pa, pb, px, py);
  }

  // --- E = alpha ln y + beta ln(1-y) - ln B ---------------------------------
  // PA is the direct assembly (relative-accurate whenever at most one
  // parameter is large: the other's lgamma contribution is then bounded);
  // PB is beta's analytic Stirling difference, which removes the
  // both-parameters-large cancellation on paper. No saturation clamp: this is
  // a log, and its whole range is meaningful.
  const auto inf = op::Set(d, std::numeric_limits<double>::infinity());
  const auto e_pa =
      DdAdd(d, DdAdd(d, DdMulD(d, lrxi, ra), DdMulD(d, lryv, rb)), cx.mlnb);
  const auto e_pb = DdSub(d, cx.pbk, ps.cpsi);
  // PB IS USED ONLY WHERE PA CANNOT COPE [SELF-CAUGHT, G3 -- and the same
  // measurement is a defect in the SHIPPED FORWARD; see the final report].
  //
  // PA is the exact form: E = alpha ln y + beta ln(1-y) (-) ln B, every term
  // dd, and its only weakness is that those terms grow with the parameters
  // while E itself stays O(1) on the ridge -- a dd carries 2^-105 RELATIVE, so
  // PA is good to 2^-105 times the largest term. Below 2^40 that is 2^-65
  // absolute and nothing else in this kernel is that accurate. PB exists for
  // the range above it, and beta's own driver reaches for it much sooner
  // (c > 256 with both parameters past Z0).
  //
  // THAT IS TOO SOON, because PB's cpsi goes through Log1pmxDd at u = -lambda
  // /alpha, and u -> -1 whenever the iterate is deep in a tail (1 + u is
  // exactly c*y/alpha). Two failures follow, both measured:
  //   * 1 + u < 2^-53: u.hi rounds to EXACTLY -1, Log1pmxDd's large branch
  //     hands LogDdAny a pair whose high word is zero, and cpsi comes back
  //     NaN. Witness (19, 1e5, 5.2e-21): true I = 3.4e-308, shipped beta_p
  //     returns exactly 0, and the boundary sits at y ~ 2.1e-20.
  //   * 1 + u merely small: Log1pmxDd adds u.lo into the low word of an
  //     EXACT TwoSum, producing a pair whose ratio lo/hi is nowhere near the
  //     2^-53 that LogDd(Dd) assumes -- its correction keeps only log1p's
  //     quadratic term, so the cubic t^3/3 survives. At (19, 1e5, 1.73e-19),
  //     t = 0.028 and t^3/3 = 7.5e-6, which alpha multiplies to 1.4e-4 in E.
  //     The reference row is right and the shipped beta_p is 1.4e-4 out.
  // So: PB only above the scale where PA runs out of digits, only where its
  // own cpsi is finite, and only where neither u nor v has collapsed onto -1.
  const auto pa_scale = op::Max(op::Max(op::Abs(e_pa.hi), op::Abs(cx.mlnb.hi)),
                                op::Abs(DdMulD(d, lrxi, ra).hi));
  const auto t_u = op::Div(y, cx.p_mean);
  const auto t_v = op::Div(yc.hi, cx.q_mean);
  const auto small = op::Set(d, kBetaInvPsiTMin);
  const auto pb_ok = op::Mul(
      op::Mul(cx.i_pb, BetaInd(d, op::Lt(op::Abs(e_pb.hi), inf))),
      op::Mul(BetaInd(d, op::Gt(pa_scale, op::Set(d, kBetaInvPaScaleMax))),
              op::Mul(BetaInd(d, op::Ge(t_u, small)),
                      BetaInd(d, op::Ge(t_v, small)))));
  const auto m_use_pb = BetaIndMask(d, pb_ok);
  auto efull = Dd<D>{op::IfThenElse(m_use_pb, e_pb.hi, e_pa.hi),
                     op::IfThenElse(m_use_pb, e_pb.lo, e_pa.lo)};
  efull = BetaInvFiniteLog(d, efull);

  // Whether a lane's E is usable decides whether the series and the continued
  // fraction may see its arguments at all -- beta scrubs on its saturation
  // mask, this file on the same finiteness test the sentinel above applies.
  const auto esat = op::Ge(op::Set(d, kBetaInvLogSentinel), efull.hi);

  // --- ln of the directly computed side ------------------------------------
  Dd<D> lnf{zero, zero};
  const Dd<D> lra{op::IfThenElse(m_sw, cx.ln_beta.hi, cx.ln_alpha.hi),
                  op::IfThenElse(m_sw, cx.ln_beta.lo, cx.ln_alpha.lo)};
  const auto e2 = DdSub(d, efull, lra);                 // R2 folds (-) ln alpha
  const auto e1 = DdSub(d, e2, DdMulD(d, lryv, rb));    // R1 also drops beta*ln y

  if (!op::AllFalse(d, m_r1)) {
    const auto as = op::IfThenElse(esat, safe_a, ra);
    const auto bs = op::IfThenElse(esat, safe_b, rb);
    const auto xs = op::IfThenElse(esat, safe_x, rxi.hi);
    const auto l = DdAdd(d, e1, BetaInvLog(d, BetaR1Series(d, as, bs, xs)));
    lnf = Dd<D>{op::IfThenElse(m_r1, l.hi, lnf.hi),
                op::IfThenElse(m_r1, l.lo, lnf.lo)};
  }
  if (!op::AllFalse(d, m_r2)) {
    const auto as = op::IfThenElse(esat, safe_a, ra);
    const auto bs = op::IfThenElse(esat, safe_b, rb);
    const auto cd = TwoSum(d, ra, rb);  // c as the EXACT pair (beta's rule)
    const Dd<D> cs{op::IfThenElse(esat, op::Add(safe_a, safe_b), cd.hi),
                   op::IfThenElse(esat, zero, cd.lo)};
    const Dd<D> xs{op::IfThenElse(esat, safe_x, rxi.hi),
                   op::IfThenElse(esat, zero, rxi.lo)};
    const auto l = DdAdd(d, e2, BetaInvLog(d, BetaR2Cf(d, as, bs, cs, xs)));
    lnf = Dd<D>{op::IfThenElse(m_r2, l.hi, lnf.hi),
                op::IfThenElse(m_r2, l.lo, lnf.lo)};
  }
  if (!op::AllFalse(d, m_r3)) {
    auto t3 = BetaR3Temme(
        d, ps, BetaInd(d, op::Ge(bmax, op::Set(d, detail::kBetaGammaLim))));
    // THE TAIL BRANCH IS TAKEN IN LOG SPACE, NOT AS log(value) [MEASURED, G3].
    // Above cpsi = 36 the R3 value is e^-cpsi times a bracket, so its log is
    // -cpsi (+) ln(bracket) exactly. Taking the log of the ASSEMBLED value
    // instead loses the whole point of a log-space forward: past cpsi ~ 708
    // the value is subnormal, its mantissa is down to a few bits, and its log
    // is QUANTIZED -- two iterates 1e-8 apart return bit-identical logits and
    // the residual has a flat spot the Newton step cannot see across.
    // Witness: (a, b, sigma) = (2.9e5, 4.0e6, 2.0e-323), where the seed is
    // already correct to 1e-8 and no step can improve it. gamma's inverse
    // takes log(value) with a -a*phi fallback and gets away with it because
    // its own targets bottom out sooner; beta's do not. The core branch
    // (cpsi <= 36) keeps log(value): the value is >= ~1e-17 there, normal by
    // a wide margin, and the bracket form has no advantage.
    const auto lg3 = BetaInvLog(d, t3.val.dd);
    const auto lgt = DdAdd(d, DdNeg(d, ps.cpsi), BetaInvLog(d, t3.brk));
    const auto tail3 = op::Gt(ps.cpsi.hi, op::Set(d, kBetaZ2Split));
    const auto pos = op::Gt(t3.val.dd.hi, zero);
    auto l = Dd<D>{op::IfThenElse(pos, lg3.hi, op::Neg(ps.cpsi.hi)),
                   op::IfThenElse(pos, lg3.lo, op::Neg(ps.cpsi.lo))};
    l = Dd<D>{op::IfThenElse(tail3, lgt.hi, l.hi),
              op::IfThenElse(tail3, lgt.lo, l.lo)};
    lnf = Dd<D>{op::IfThenElse(m_r3, l.hi, lnf.hi),
                op::IfThenElse(m_r3, l.lo, lnf.lo)};
    is_p = op::IfThenElse(m_r3, t3.lam_pos, is_p);
  }

  // --- SEVENTH ROUTING CORRECTION: near-one post-route ----------------------
  // An R1 lane whose value exceeds kBetaNearOne would hand its complement back
  // as dd rounding noise. Such lanes are R4-SHAPED by construction, and R4
  // returns the SMALL side directly -- which for a logit is strictly better
  // than any complement, since the small side IS the limb that carries the
  // signal. The test is on the value, so it needs one exponential of a log
  // that is already bounded near zero there.
  const auto i_pr = op::Mul(
      i_r1,
      op::Mul(BetaInd(d, op::Gt(lnf.hi, op::Set(d, kBetaInvNearOneLog))),
              BetaInd(d, op::Ge(op::Set(d, detail::kBetaPrTauMax), ra))));
  const auto m_pr = BetaIndMask(d, i_pr);
  const auto i_r4x = op::Max(i_r4, i_pr);
  const auto m_r4x = BetaIndMask(d, i_r4x);

  if (!op::AllFalse(d, m_r4x)) {
    // R4's scrub point is its OWN interior (tau <= eps_R4), not the shared
    // (2, 3, 1/4): tau = 2 is outside this core's zones.
    const auto t4 =
        op::IfThenElse(m_r4x, ra, op::Set(d, detail::kBetaEpsR4));
    const auto b4 = op::IfThenElse(m_r4x, rb, safe_b);
    const auto x4 = op::IfThenElse(m_r4x, rxi.hi, safe_x);
    const Dd<D> l4{
        op::IfThenElse(m_r4x, lrxi.hi, op::Set(d, -2.0 * detail::kLogLn2Hi)),
        op::IfThenElse(m_r4x, lrxi.lo, op::Set(d, -2.0 * detail::kLogLn2Lo))};
    const auto q4 = BetaR4Tiny(d, t4, b4, x4, l4);
    const auto l = BetaInvLog(d, q4);
    lnf = Dd<D>{op::IfThenElse(m_r4x, l.hi, lnf.hi),
                op::IfThenElse(m_r4x, l.lo, lnf.lo)};
    is_p = op::IfThenElse(m_pr, i_sw, is_p);
  }

  // --- (C) gamma-limit slice: R2's off-band remainder above B_GL ------------
  // Transcribed from beta-inl.h, including the BOTH-HUGE exclusion; the only
  // change is that the value is kept in log space (e_pick + ln(core)) instead
  // of being exponentiated and clamped.
  const auto i_gl = op::Mul(
      op::Mul(i_r2,
              BetaInd(d, op::Ge(bmax, op::Set(d, detail::kBetaGammaLim)))),
      BetaInd(d, op::Gt(op::Set(d, detail::kBetaGammaLim), tau)));
  const auto m_gl = BetaIndMask(d, i_gl);
  if (!op::AllFalse(d, m_gl)) {
    const auto hf = op::Ge(ra, op::Set(d, detail::kBetaGammaLim));
    const auto i_hf = BetaInd(d, hf);
    const auto ss = op::IfThenElse(m_gl, op::IfThenElse(hf, rb, ra), one);
    const auto huge = op::IfThenElse(hf, ra, rb);
    const Dd<D> lx_gl{op::IfThenElse(hf, lrxi.hi, lryv.hi),
                      op::IfThenElse(hf, lrxi.lo, lryv.lo)};
    const auto t_dd = DdMulD(d, lx_gl, op::Neg(huge));  // t > 0, dd
    auto e_g = DdMulD(d, BetaInvLog(d, t_dd), ss);
    e_g = DdSub(d, e_g, t_dd);
    e_g = DdSub(d, e_g, LgammaPosDd(d, ss));
    const auto th = t_dd.hi;
    const auto i_ser = op::Max(
        op::Mul(BetaInd(d, op::Lt(ss, op::Set(d, detail::kGammaAT))),
                BetaInd(d, op::Ge(op::Add(ss, one), th))),
        op::Mul(BetaInd(d, op::Ge(ss, op::Set(d, detail::kGammaAT))),
                BetaInd(d, op::Ge(ss, op::Add(th, th)))));
    const auto m_ser = BetaIndMask(d, i_ser);
    const auto e_s = DdSub(d, e_g, BetaInvLog(d, ss));  // gamma's R1 fold
    Dd<D> e_pick{op::IfThenElse(m_ser, e_s.hi, e_g.hi),
                 op::IfThenElse(m_ser, e_s.lo, e_g.lo)};
    e_pick = BetaInvFiniteLog(d, e_pick);
    const auto gsat = op::Ge(op::Set(d, kBetaInvLogSentinel), e_pick.hi);
    const auto s_c = op::IfThenElse(gsat, one, ss);
    const auto t_c = op::IfThenElse(gsat, op::Set(d, 3.0),
                                    op::Min(th, op::Set(d, 4e306)));
    const auto ser = GammaSeriesSum(d, s_c, t_c);
    const auto cfr = GammaCfRecip(d, s_c, t_c);
    const Dd<D> fac{op::IfThenElse(m_ser, ser.hi, cfr.hi),
                    op::IfThenElse(m_ser, ser.lo, cfr.lo)};
    const auto l = DdAdd(d, e_pick, BetaInvLog(d, fac));
    lnf = Dd<D>{op::IfThenElse(m_gl, l.hi, lnf.hi),
                op::IfThenElse(m_gl, l.lo, lnf.lo)};
    // val holds P_gamma on series lanes, Q_gamma on CF lanes; P_gamma is the
    // ROUTED value iff the huge parameter is SECOND, and the routed value is
    // P-of-this-frame iff the orientation was native. XNOR twice.
    const auto agree = op::MulAdd(i_ser, BetaIndNot(d, i_hf),
                                  op::Mul(BetaIndNot(d, i_ser), i_hf));
    const auto isp_gl = op::MulAdd(agree, i_nat,
                                   op::Mul(BetaIndNot(d, agree), i_sw));
    is_p = op::IfThenElse(m_gl, isp_gl, is_p);
  }

  lnf = BetaInvFiniteLog(d, lnf);

  // --- assemble the logit ---------------------------------------------------
  // Put the SMALLER side in hand first: v > 1/2 <=> ln v > -ln 2. The region
  // map picks the smaller side by construction, so this fires only near a
  // region's own switch and on iterates far from the root, where the value
  // being complemented is >= ~0.4 and the subtraction costs a fraction of a
  // bit. The value is recovered by exponentiating its log rather than being
  // assembled a second time.
  {
    const auto m_cmp = op::Gt(lnf.hi, op::Set(d, kBetaInvNegLn2));
    if (!op::AllFalse(d, m_cmp)) {
      const auto vb = BetaInvExpDd(d, lnf.hi, lnf.lo);
      const auto cb = DdSub(d, Dd<D>{one, zero}, vb);
      const auto l = BetaInvFiniteLog(d, BetaInvLog(d, cb));
      lnf = Dd<D>{op::IfThenElse(m_cmp, l.hi, lnf.hi),
                  op::IfThenElse(m_cmp, l.lo, lnf.lo)};
      is_p = op::IfThenElse(m_cmp, BetaIndNot(d, is_p), is_p);
    }
  }
  const auto u = BetaInvExpDd(d, lnf.hi, lnf.lo);  // u = min(P, Q) <= 1/2
  const auto cmpl = DdSub(d, Dd<D>{one, zero}, u);
  const auto lnc = BetaInvLog(d, cmpl);
  const auto lg = DdSub(d, lnf, lnc);
  const auto m_p = BetaIndMask(d, is_p);
  Dd<D> m{op::IfThenElse(m_p, lg.hi, op::Neg(lg.hi)),
          op::IfThenElse(m_p, lg.lo, op::Neg(lg.lo))};

  // [TRANSFER SITE] w = exp(ln u (-) E) * (1 - u) * (1 - y). The last factor
  // is beta's own; gamma's slope has no analogue for it (file header).
  const auto wl = op::Max(op::Min(op::Sub(lnf.hi, efull.hi),
                                  op::Set(d, kBetaInvWExpMax)),
                          op::Set(d, -kBetaInvWExpMax));
  auto w = op::Mul(op::Mul(BetaInvExp(d, wl), cmpl.hi), yc.hi);

  // THE RESIDUAL'S OWN UNCERTAINTY, and the kernel needs it because the
  // forward is NOT dd-accurate near the median [MEASURED, G3].
  //   * ln u is taken straight from the log-space assembly, so its error is
  //     the dd's own, ~2^-100 relative to |m|.
  //   * ln(1 - u) is NOT: it needs the VALUE u, which comes from exp_dd, whose
  //     own documented budget is ~2^-70 relative (polynomial truncation 2^-72
  //     and the dropped r.lo at 2^-70.5; src/exp_dd-inl.h states it). That
  //     error lands on m amplified by u/(1-u) <= 1, since u is the smaller
  //     side by construction.
  // So the noise on m is ~2^-68*u + 2^-98*|m|: negligible in every tail (u is
  // tiny there and ln u carries the signal), and ~2^-69 near the median. It
  // matters because the inverse multiplies it by w, the condition number:
  // once w exceeds ~2^17 the step the residual asks for is pure noise, and
  // without this the step loop random-walks the answer AWAY from a good seed
  // (measured: a joint-tiny lane went from 1.3e-9 to 1.4e-7 relative over its
  // four steps, the trust bypass accepting every one). Reported by the
  // forward, consumed by the step loop's freeze.
  const auto unc = op::MulAdd(op::Set(d, 0x1.0p-68), op::Abs(u.hi),
                              op::Mul(op::Set(d, 0x1.0p-98), op::Abs(m.hi)));

  // A lane whose forward could not be evaluated reports a huge residual (so it
  // loses every seed comparison) and a zero slope (so it freezes rather than
  // stepping). That is the correct behaviour, not a fallback: the same thing
  // happens for real in the beyond-resolution regime, where w is genuinely
  // 2^-500-class and the iterate must not move.
  const auto okv = op::Mul(op::Mul(BetaInd(d, op::Lt(op::Abs(m.hi), inf)),
                                   BetaInd(d, op::Lt(op::Abs(m.lo), inf))),
                           BetaInd(d, op::Lt(op::Abs(w), inf)));
  const auto m_ok = BetaIndMask(d, okv);
  m = Dd<D>{op::IfThenElse(m_ok, m.hi, op::Set(d, kBetaInvBigResid)),
            op::IfThenElse(m_ok, m.lo, zero)};
  w = op::IfThenElse(m_ok, w, zero);
  return {m, w, op::IfThenElse(m_ok, unc, zero)};
}

// m(y) - m_t, rounded once. Both are dd, and near the root they cancel to
// nothing -- which is the entire reason they are dd.
template <class D>
HWY_INLINE op::V<D> BetaInvResid(D d, const BetaInvFwdOut<D>& f, Dd<D> mt) {
  return DdToDouble(DdSub(d, f.m, mt));
}

// ------------------------------------------------------------------------
// S1: the beta-Temme ridge seed.
//
// zeta -> lambda: invert cpsi(lambda) = zeta^2 nu on lambda in (-beta, alpha),
// with sign(lambda) = sign(zeta). NEWTON ON THE SQUARE ROOT, not on cpsi (see
// kBetaInvS1LamIters for why), with the analytic derivative
//     d cpsi/d lambda = lambda c / ((alpha - lambda)(beta + lambda))
// and the 0/0 at lambda = 0 removed algebraically rather than special-cased:
//     H := cpsi/lambda^2 = phi(u)/(alpha u^2) + phi(v)/(beta v^2),
//          u = -lambda/alpha,  v = lambda/beta
// (each term is phi(w)/w^2 divided by the parameter, and both are finite and
// positive through zero), so
//     F = lambda sqrt(H/nu),
//     lambda_new = lambda - (F - zeta) * 2 sqrt(nu H) (alpha-lambda)(beta+lambda)/c.
// At lambda = 0 that multiplier is exactly sqrt(2) nu, matching the
// leading-order relation lambda ~ sqrt(2) nu zeta the iteration starts from.
// The bracket is maintained and a step leaving it falls back to the midpoint,
// so the iteration count is an upper bound and not an assumption.
//
// EVERYTHING IS IN THE PRESCALED FRAME (as, bs, cs, nus): lambda and cpsi are
// 1-homogeneous, zeta and y are 0-homogeneous, so the answer is unchanged
// while c = alpha + beta stays representable for a pair near DBL_MAX and no
// operand exceeds ops::ProdLow's non-FMA Dekker ceiling. The generator has no
// prescale -- it lets huge parameters raise, and its caller discards S1 there.
template <class D>
HWY_NOINLINE op::V<D> BetaInvLamOfZeta(D d, op::V<D> zeta, op::V<D> as,
                                       op::V<D> bs, op::V<D> cs,
                                       op::V<D> nus) {
  const auto one = op::Set(d, 1.0);
  const auto zero = op::Zero(d);
  const auto half = op::Set(d, 0.5);
  const auto inset = op::Set(d, 1.0 - kBetaInvLamInset);
  const auto pos = op::Ge(zeta, zero);

  auto lo = op::IfThenElse(pos, zero, op::Neg(op::Mul(bs, inset)));
  auto hi = op::IfThenElse(pos, op::Mul(as, inset), zero);
  // Leading order: cpsi ~ lambda^2/(2 nu), hence zeta ~ lambda/(nu sqrt 2).
  // OUT OF BRACKET -> MIDPOINT, never a clamp to the endpoint [SELF-CAUGHT,
  // G3]. The leading-order guess overshoots for any zeta the ridge expansion
  // is not small at (sqrt2*nu*zeta exceeds alpha once zeta > p*sqrt(2 c/beta)
  // -- at nu = 228, zeta = 0.78 it is 250 against an alpha of 232), and
  // clamping lands the iteration exactly where its own derivative vanishes:
  // dlambda/dF carries a factor (alpha - lambda), so a start at the edge
  // crawls and eight steps do not recover. Measured at
  // (231.7, 14752.5, 4.7e-62): the clamped start returned 4.7e-4 for a true
  // 4.1e-3 and lost the residual comparison to a much weaker seed.
  const auto lam0 = op::Mul(op::Mul(op::Set(d, kBetaInvSqrt2), nus), zeta);
  const auto in0 = op::Mul(BetaInd(d, op::Gt(lam0, lo)),
                           BetaInd(d, op::Lt(lam0, hi)));
  auto lam = op::IfThenElse(BetaIndMask(d, in0), lam0,
                            op::Mul(op::Add(lo, hi), half));

  for (int i = 0; i < kBetaInvS1LamIters; ++i) {
    const auto u = op::Neg(op::Div(lam, as));
    const auto v = op::Div(lam, bs);
    const auto h = op::Add(op::Div(BetaInvPhiOverSq(d, u), as),
                           op::Div(BetaInvPhiOverSq(d, v), bs));
    const auto f = op::Mul(lam, op::Sqrt(op::Div(h, nus)));
    // F increases with lambda, so F < zeta makes lambda a lower bound.
    const auto below = op::Lt(f, zeta);
    lo = op::IfThenElse(below, op::Max(lo, lam), lo);
    hi = op::IfThenElse(below, hi, op::Min(hi, lam));
    // (alpha-lam)/c first, then times (beta+lam): the raw product
    // (alpha-lam)(beta+lam) overflows for two parameters near 2^900 while the
    // quotient form is bounded by beta.
    const auto mult = op::Mul(
        op::Mul(op::Add(one, one), op::Sqrt(op::Mul(nus, h))),
        op::Mul(op::Div(op::Sub(as, lam), cs), op::Add(bs, lam)));
    auto cand = op::Sub(lam, op::Mul(op::Sub(f, zeta), mult));
    // INCLUSIVE bounds, and the reason is not fussiness [SELF-CAUGHT, G3]: the
    // bracket was just tightened ONTO lambda itself, so a CONVERGED Newton
    // step -- which returns lambda unchanged -- sits exactly on an endpoint. A
    // strict test rejects it and substitutes the midpoint of a bracket whose
    // other end is still the initial 0 or -beta, throwing a converged iterate
    // halfway across the domain. Measured before the fix: at (a, b) =
    // (1.9e39, 1.0e42) the seed came back 0.00575 for a true 0.001848, and at
    // (19, 1e5) it walked to the bracket edge and stayed there. Every other
    // safeguard in this file has the same shape and the same reason.
    const auto inb = op::Mul(BetaInd(d, op::Ge(cand, lo)),
                             BetaInd(d, op::Ge(hi, cand)));
    cand = op::IfThenElse(BetaIndMask(d, inb), cand,
                          op::Mul(op::Add(lo, hi), half));
    lam = cand;
  }
  return lam;
}

// c_k(zeta, p) from the pinned 15x9 Chebyshev rows: Clenshaw in p inside each
// zeta row, then Clenshaw across the rows. The stored half is p <= 1/2; for
// p > 1/2 the table is evaluated at (-zeta, q) and NEGATED, exactly as R3's
// own e_k symmetry (the generator's _s1_ck_eval).
template <class D>
HWY_INLINE op::V<D> BetaInvS1Ck(D d, int k, op::V<D> t, op::V<D> t2,
                                op::V<D> u, op::V<D> u2) {
  auto b1 = op::Zero(d);
  auto b2 = op::Zero(d);
  for (int n = detail::kBetaInvS1NZ - 1; n >= 1; --n) {
    const double* row = detail::kBetaInvS1Cheb[k][n];
    auto c1 = op::Set(d, row[detail::kBetaInvS1NP - 1]);
    auto c2 = op::Zero(d);
    for (int mI = detail::kBetaInvS1NP - 2; mI >= 1; --mI) {
      const auto nb = op::Sub(op::MulAdd(u2, c1, op::Set(d, row[mI])), c2);
      c2 = c1;
      c1 = nb;
    }
    const auto rn = op::Sub(op::MulAdd(u, c1, op::Set(d, row[0])), c2);
    const auto nb = op::Sub(op::MulAdd(t2, b1, rn), b2);
    b2 = b1;
    b1 = nb;
  }
  const double* row0 = detail::kBetaInvS1Cheb[k][0];
  auto c1 = op::Set(d, row0[detail::kBetaInvS1NP - 1]);
  auto c2 = op::Zero(d);
  for (int mI = detail::kBetaInvS1NP - 2; mI >= 1; --mI) {
    const auto nb = op::Sub(op::MulAdd(u2, c1, op::Set(d, row0[mI])), c2);
    c2 = c1;
    c1 = nb;
  }
  const auto r0 = op::Sub(op::MulAdd(u, c1, op::Set(d, row0[0])), c2);
  return op::Sub(op::MulAdd(t, b1, r0), b2);
}

// S1 proper. zeta0 = z/sqrt(nu) with z = erfcinv(2 tau) -- NO factor of two,
// because beta's cpsi IS zeta^2 nu directly while gamma's is a eta^2/2 (the
// generator's own first transfer bug). z arrives already signed from the
// driver, which forms it as +-erfcinv(2 sigma): erfcinv(2(1-sigma)) =
// -erfcinv(2 sigma) exactly, so a target on the large side costs a sign and
// not an evaluation. The correction is one perturbative Newton step in zeta,
//     zeta <- zeta - sign(zeta) S(zeta, p)/(nu sqrt 2),  S = sum c_k/nu^k,
// gated to the table's own fitted domain |zeta| <= kBetaInvS1ZetaMax AND to
// nu >= kBetaInvS1NuMin (the 1/nu series is asymptotic; below nu ~ 2 the
// applied correction is unbounded -- measured).
template <class D>
HWY_NOINLINE op::V<D> BetaInvSeedS1(D d, const BetaInvCtx<D>& cx, op::V<D> z,
                                    op::V<D> as, op::V<D> bs, op::V<D> cs,
                                    op::V<D> nus) {
  const auto one = op::Set(d, 1.0);
  const auto zero = op::Zero(d);
  auto zeta = op::Div(z, op::Sqrt(cx.nu));
  zeta = op::IfThenElse(op::IsNaN(zeta), zero, zeta);

  const auto zmax = op::Set(d, detail::kBetaInvS1ZetaMax);
  const auto gate = op::Mul(
      BetaInd(d, op::Ge(zmax, op::Abs(zeta))),
      BetaInd(d, op::Ge(cx.nu, op::Set(d, detail::kBetaInvS1NuMin))));
  if (!op::AllFalse(d, BetaIndMask(d, gate))) {
    // Table symmetry: only p <= 1/2 is stored.
    const auto half = op::Set(d, 0.5);
    const auto pswap = op::Gt(cx.p_mean, half);
    const auto pe = op::IfThenElse(pswap, cx.q_mean, cx.p_mean);
    const auto ze = op::IfThenElse(pswap, op::Neg(zeta), zeta);
    const auto tsign = op::IfThenElse(pswap, op::Set(d, -1.0), one);
    const auto t = op::Mul(ze, op::Set(d, 1.0 / detail::kBetaInvS1ZetaMax));
    const auto t2 = op::Add(t, t);
    const auto u = op::Mul(op::Sub(pe, op::Set(d, detail::kBetaInvS1PMid)),
                           op::Set(d, 1.0 / detail::kBetaInvS1PHalf));
    const auto u2 = op::Add(u, u);
    const auto rnu = op::Div(one, cx.nu);
    auto s = zero;
    for (int k = detail::kBetaInvS1NCorr - 1; k >= 0; --k) {
      s = op::MulAdd(s, rnu, BetaInvS1Ck(d, k, t, t2, u, u2));
    }
    s = op::Mul(s, tsign);
    const auto sgn = op::IfThenElse(op::Ge(zeta, zero), one, op::Set(d, -1.0));
    const auto corr =
        op::Mul(op::Mul(sgn, s),
                op::Div(rnu, op::Set(d, kBetaInvSqrt2)));
    zeta = op::IfThenElse(BetaIndMask(d, gate), op::Sub(zeta, corr), zeta);
  }

  const auto lam = BetaInvLamOfZeta(d, zeta, as, bs, cs, nus);
  return op::Div(op::Sub(as, lam), cs);  // y = (alpha - lambda)/c
}

// ------------------------------------------------------------------------
// S2: small-y series inversion.
//
// R1's own value is P = y^alpha/B * sum_{n>=0} t_n/(alpha+n), whose n = 0 term
// is 1/alpha -- so P ~ y^alpha/(alpha B) as y -> 0 and the zeroth iterate is
//     y0 = exp((ln sigma + ln alpha + lnB)/alpha) = exp((ln sigma + LB)/alpha),
// with LB = ln alpha (+) ln B taken from the context, which is where the
// generator's "+ln alpha" bug (catastrophic at tiny alpha, and the same bug
// again in its deep-small twin) cannot recur: the term is never dropped
// because it is never separable -- LB IS lgamma(1+alpha)+lgamma(beta)-lgamma(c).
// The Picard correction is the EXACT fixed point of the defining equation,
//     ln y = (ln sigma + ln B - ln S)/alpha,  S = sum t_n/(alpha+n),
// which is why this seed alone recovers the tiny-alpha corner. It is written
// on R1's own alpha*S (= 1 + alpha sum_{n>=1}) rather than on S, so nothing
// forms 1/alpha and a subnormal alpha stays finite:
//     ln y = (ln sigma + LB - ln(alpha S))/alpha.
template <class D>
HWY_NOINLINE op::V<D> BetaInvS2Series(D d, op::V<D> a, op::V<D> b, op::V<D> y) {
  const auto one = op::Set(d, 1.0);
  const auto zero = op::Zero(d);
  const auto half = op::Set(d, 0.5);
  const auto eps = op::Set(d, kBetaInvSeedEps);
  auto t = one;
  auto s = zero;
  auto live = one;
  for (int n = 1; n <= kBetaInvS2SeriesN; ++n) {
    const auto nv = op::Set(d, static_cast<double>(n));
    t = op::Mul(t, op::Div(op::Mul(op::Sub(nv, b), y), nv));
    const auto contrib = op::Div(t, op::Add(a, nv));
    const auto lm = op::Gt(live, half);
    s = op::IfThenElse(lm, op::Add(s, contrib), s);
    live = op::IfThenElse(
        op::Lt(op::Abs(contrib), op::Mul(op::Abs(s), eps)), zero, live);
    if (op::AllTrue(d, op::Eq(live, zero))) break;
  }
  return op::MulAdd(a, s, one);
}

template <class D>
HWY_NOINLINE op::V<D> BetaInvSeedS2(D d, const BetaInvCtx<D>& cx,
                                    op::V<D> lsum) {
  const auto zero = op::Zero(d);
  const auto one = op::Set(d, 1.0);
  const auto inf = op::Set(d, std::numeric_limits<double>::infinity());
  auto y = BetaInvExp(d, op::Div(lsum, cx.alpha));
  for (int i = 0; i < detail::kBetaInvS2NCorr; ++i) {
    const auto ser = BetaInvS2Series(d, cx.alpha, cx.beta, y);
    const auto ly = op::Div(op::Sub(lsum, BetaInvLog(d, ser).hi), cx.alpha);
    const auto yn = BetaInvExp(d, ly);
    // The generator breaks out of the loop on a non-finite or out-of-range
    // iterate; per lane that is holding the previous value. Written through
    // the indicator arithmetic because the ops facade exposes no mask AND.
    const auto good = BetaIndMask(
        d, op::Mul(op::Mul(BetaInd(d, op::Gt(yn, zero)),
                           BetaInd(d, op::Lt(yn, one))),
                   BetaInd(d, op::Lt(op::Abs(ser), inf))));
    y = op::IfThenElse(good, yn, y);
  }
  return y;
}

// ------------------------------------------------------------------------
// One candidate seed: evaluate the forward at it, score it by |m - m_t|, and
// keep it if it beats the incumbent -- carrying the forward result along, so
// the winner's evaluation becomes the first Newton step's and no candidate is
// evaluated twice.
//
// A FREE FUNCTION, NOT A LAMBDA. A closure type's operator() is an
// implicitly-declared member instantiated outside the per-target attribute
// region -- the same shape as the aggregate-initialization hazard AGENTS.md
// records (NEON_BF16 refusing to inline Vec128's always_inline constructors
// into an unattributed function, which broke macOS CI twice).
template <class D>
HWY_NOINLINE void BetaInvConsider(D d, const BetaInvCtx<D>& cx, op::V<D> cand,
                                op::V<D> gate, Dd<D> mt, op::V<D>* best_y,
                                op::V<D>* best_r, BetaInvFwdOut<D>* best_f) {
  const auto zero = op::Zero(d);
  const auto one = op::Set(d, 1.0);
  const auto inf = op::Set(d, std::numeric_limits<double>::infinity());
  const auto okc = op::Mul(gate, op::Mul(BetaInd(d, op::Gt(cand, zero)),
                                         BetaInd(d, op::Lt(cand, one))));
  const auto mc = BetaIndMask(d, okc);
  // A rejected candidate is replaced by an ordinary interior point before the
  // forward runs: masked-off lanes execute every op, and a NaN y would
  // otherwise reach the series lengths and the region cores' gathers.
  const auto ys = op::IfThenElse(mc, cand, op::Set(d, kBetaInvSafeY));
  const auto f = BetaInvForward(d, cx, ys);
  auto r = op::Abs(BetaInvResid(d, f, mt));
  r = op::IfThenElse(mc, r, inf);
  r = op::IfThenElse(op::IsNaN(r), inf, r);
  const auto take = op::Lt(r, *best_r);
  *best_y = op::IfThenElse(take, ys, *best_y);
  *best_r = op::IfThenElse(take, r, *best_r);
  best_f->m.hi = op::IfThenElse(take, f.m.hi, best_f->m.hi);
  best_f->m.lo = op::IfThenElse(take, f.m.lo, best_f->m.lo);
  best_f->w = op::IfThenElse(take, f.w, best_f->w);
  best_f->unc = op::IfThenElse(take, f.unc, best_f->unc);
}

// t (the gamma-limit variable) back to y. On the beta >= alpha branch
// y = -expm1(-t/beta) keeps its relative accuracy for a tiny quotient, which
// is that branch's whole regime; on the other branch y = exp(-t/alpha) is
// near 1 and a plain exponential is what the answer can carry.
template <class D>
HWY_INLINE op::V<D> BetaInvGammaToY(D d, op::V<D> t, op::V<D> hg,
                                    op::M<D> blo) {
  const auto r = op::Div(t, hg);
  const auto e = Expm1Dd(d, Dd<D>{op::Neg(r), op::Zero(d)});
  return op::IfThenElse(blo, op::Neg(e.hi), BetaInvExp(d, op::Neg(r)));
}

// ------------------------------------------------------------------------
// The driver. kQ selects beta_q_inv (true) or beta_p_inv (false); the two
// differ only in how the input-side flip sets the orientation bit.
//
// HWY_NOINLINE for the AGENTS.md reason: it is inlined twice per export
// (full-vector and masked-tail call sites) and is the heaviest driver in the
// library.
template <bool kQ, class D>
HWY_NOINLINE op::V<D> BetaInvVec(D d, op::V<D> a_in, op::V<D> b_in,
                                 op::V<D> s_in) {
  const auto one = op::Set(d, 1.0);
  const auto zero = op::Zero(d);
  const auto half = op::Set(d, 0.5);
  const auto inf = op::Set(d, std::numeric_limits<double>::infinity());
  const auto qnan = op::Set(d, std::numeric_limits<double>::quiet_NaN());

  // --- scrub every lane the specials table decides -------------------------
  // Masked-off lanes still execute every op (AGENTS.md), and this kernel takes
  // logs of y and 1-y, divides by alpha and beta, multiplies (c+m) and reaches
  // value-derived gathers inside erfc's core: a NaN, a zero, a one or an
  // infinity left in place would propagate rather than sit quietly.
  const auto i_ok =
      op::Mul(op::Mul(BetaInd(d, op::Gt(a_in, zero)),
                      BetaInd(d, op::Lt(a_in, inf))),
              op::Mul(op::Mul(BetaInd(d, op::Gt(b_in, zero)),
                              BetaInd(d, op::Lt(b_in, inf))),
                      op::Mul(BetaInd(d, op::Gt(s_in, zero)),
                              BetaInd(d, op::Lt(s_in, one)))));
  const auto m_ok = BetaIndMask(d, i_ok);
  const auto a = op::IfThenElse(m_ok, a_in, op::Set(d, kBetaSafeA));
  const auto b = op::IfThenElse(m_ok, b_in, op::Set(d, kBetaSafeB));
  const auto s = op::IfThenElse(m_ok, s_in, op::Set(d, kBetaSafeX));

  // --- the two exact relabellings (file header) ----------------------------
  // (1) sigma = min(s, 1-s), the subtraction EXACT by Sterbenz. i_sq records
  // whether sigma is the Q of the CALLER's (a, b, x); that is a statement
  // about the probability, and it does NOT decide the orientation.
  const auto m_flip = op::Gt(s, half);
  const auto sigma = op::IfThenElse(m_flip, op::Sub(one, s), s);  // exact
  const auto i_flip = BetaInd(d, m_flip);
  const auto i_sq = kQ ? BetaIndNot(d, i_flip) : i_flip;
  const auto m_sq = BetaIndMask(d, i_sq);

  // ln sigma and ln(1 - sigma). 1 - sigma is formed as an EXACT dd: below
  // sigma ~ 2^-53 the low word is the entire signal, and a native-float
  // 1 - sigma would collapse to exactly 1.0.
  const auto lnsig = BetaInvLog(d, Dd<D>{sigma, zero});
  const auto lnsigc = BetaInvLog(d, TwoSum(d, one, op::Neg(sigma)));
  const Dd<D> lnp0{op::IfThenElse(m_sq, lnsigc.hi, lnsig.hi),
                   op::IfThenElse(m_sq, lnsigc.lo, lnsig.lo)};
  const Dd<D> lnq0{op::IfThenElse(m_sq, lnsig.hi, lnsigc.hi),
                   op::IfThenElse(m_sq, lnsig.lo, lnsigc.lo)};
  const auto mt0 = DdSub(d, lnp0, lnq0);  // target logit in the caller's frame

  // (2) THE ORIENTATION PROBE (file header). x > 1/2 exactly when
  // I_{1/2}(a,b) < P, i.e. when the forward's own logit at the midpoint falls
  // below the target. One forward evaluation, no new machinery, and the only
  // lanes it can misclassify have their root within an ulp of 1/2 -- where
  // both orientations are equally good. A lane whose midpoint forward is
  // unusable reports the huge positive sentinel and stays unswapped, which is
  // the right answer for the beyond-resolution pairs that produce it (their
  // root is the mean, and 1/2 to within an ulp when they are balanced).
  // The comparison carries a DEAD BAND of the probe's own uncertainty. Where
  // the root really is 1/2 the two sides are equal in exact arithmetic and the
  // sign of their difference is decided by the forward's last bits -- a coin
  // flip between two orientations that are equally correct but NOT equally
  // easy: (a, b, p) = (20, 1, 2^-20) has x = 1/2 exactly, and the unswapped
  // frame answers it in closed form (beta = 1 makes I_x = x^alpha, so the
  // deep-small route is exact) while the swapped one has to iterate there from
  // a 0.85-residual seed and lands 4e-12 out. Ties therefore go to "no swap".
  const auto cx0 = BetaInvPrepare(d, a, b);
  const auto f_half = BetaInvForward(d, cx0, half);
  const auto d_half = DdToDouble(DdSub(d, f_half.m, mt0));
  // The dead band is measured in the PROBE'S OWN STEP, not in the residual:
  // w*(m - m_t) is the relative distance from 1/2 to the root, so requiring it
  // to exceed kBetaInvSwapBand says "swap only when the answer is genuinely
  // away from the midpoint". Inside the band both orientations return a value
  // within 2^-10 of 1/2, where 1 (-) y is Sterbenz-exact and neither
  // orientation can lose anything -- so the tie may safely go to the cheaper
  // one, and the criterion needs no assumption about the forward's accuracy.
  const auto step_half = op::Abs(op::Mul(f_half.w, d_half));
  const auto m_swap = BetaIndMask(
      d, op::Mul(BetaInd(d, op::Lt(d_half, zero)),
                 BetaInd(d, op::Gt(step_half,
                                   op::Set(d, kBetaInvSwapBand)))));
  const auto i_swap = BetaInd(d, m_swap);
  const auto cx = BetaInvSwapCtx(d, cx0, m_swap);
  const auto alpha = cx.alpha;
  const auto beta = cx.beta;

  // The internal target tau = I_y(alpha, beta) and its complement, both taken
  // from the pair of exact logs above -- never recomputed, never complemented.
  const Dd<D> lntau{op::IfThenElse(m_swap, lnq0.hi, lnp0.hi),
                    op::IfThenElse(m_swap, lnq0.lo, lnp0.lo)};
  const Dd<D> lntauc{op::IfThenElse(m_swap, lnp0.hi, lnq0.hi),
                     op::IfThenElse(m_swap, lnp0.lo, lnq0.lo)};
  const auto mt = DdSub(d, lntau, lntauc);
  // tau IS sigma (hence tau <= 1/2) exactly when the two bits agree.
  const auto i_tis = op::MulAdd(i_swap, i_sq,
                                op::Mul(BetaIndNot(d, i_swap),
                                        BetaIndNot(d, i_sq)));
  const auto m_tis = BetaIndMask(d, i_tis);
  // z = erfcinv(2 tau), from the ONE inverse the frame can evaluate exactly:
  // erfcinv(2(1 - sigma)) = -erfcinv(2 sigma), and 2 sigma is exact.
  const auto zsig = ErfcinvVec(d, op::Add(sigma, sigma));
  const auto zq = op::IfThenElse(m_tis, zsig, op::Neg(zsig));

  // --- deep-small closed form ----------------------------------------------
  // y = exp((ln sigma + ln alpha + lnB)/alpha), with the power-of-two scaling
  // applied LAST so a subnormal or zero answer carries exactly one rounding.
  //
  // THE CUT IS THE THIRD CORRECTION'S (src/betainv_data.h), and it is NOT any
  // parameter*y form. What the closed form drops is the series factor S',
  // whose log divided by alpha has leading term |1-beta| y/(1+alpha): the
  // OTHER side's parameter supplies the coefficient (beta is the n = 1 series
  // coefficient; alpha enters only through the 1+alpha denominator), which is
  // exactly the dependence gamma cannot have and the reason its x0*(1+a) form
  // could not simply be carried over. The leading term alone under-predicts by
  // up to 13.8x at the widened gamma-limit corner, so it carries the exact
  // closed-form multiplier corr(y) = -ln(1-y)/y, which is exact in the
  // huge-other-side limit and a verified sound upper predictor elsewhere.
  // corr -> 1 as y -> 0, so the test reduces to the bare leading term exactly
  // where that was already valid.
  //
  // ONE FORM, BOTH ORIENTATIONS: the generator's q-side twin is this
  // expression after the frame's relabelling (file header). There is no second
  // branch for a self-check to miss.
  const auto ds_num = DdAdd(d, lntau, cx.lb_a);
  const auto ds_q = GammaInvDivD(d, ds_num, alpha);
  const auto ds_e = ExpDdFrac(d, ds_q.hi, ds_q.lo);
  const auto y_ds = ScaleTwo(d, DdToDouble(ds_e.m), ds_e.e);
  const auto corr_num = op::Neg(BetaInvLog(d, TwoSum(d, one, op::Neg(y_ds))).hi);
  const auto corr = op::IfThenElse(op::Lt(y_ds, op::Set(d, 1e-8)), one,
                                   op::Div(corr_num, y_ds));
  const auto ds_bound =
      op::Mul(op::Mul(op::Abs(op::Sub(one, beta)),
                      op::Div(y_ds, op::Add(one, alpha))), corr);
  const auto m_ds =
      BetaIndMask(d, op::Mul(BetaInd(d, op::Lt(ds_bound,
                                    op::Set(d, detail::kBetaInvDeepSmallCut))),
                             BetaInd(d, op::Lt(y_ds, one))));

  auto y = y_ds;
  if (!op::AllTrue(d, m_ds)) {
    // --- seed candidates ---------------------------------------------------
    // The prescaled frame S1 works in; also the source of the fallback
    // incumbent (the distribution's mean, an ordinary interior point for any
    // parameter pair, including one whose alpha + beta overflows).
    const auto big = op::Gt(op::Max(alpha, beta), op::Set(d, kBetaScaleAbove));
    const auto dn = op::IfThenElse(big, op::Set(d, kBetaScaleDown), one);
    const auto as = op::Mul(alpha, dn);
    const auto bs = op::Mul(beta, dn);
    const auto cs = op::Add(as, bs);
    const auto nus = op::Mul(as, op::Div(bs, cs));

    auto best_y = cx.p_mean;
    best_y = op::IfThenElse(
        BetaIndMask(d, op::Mul(BetaInd(d, op::Gt(best_y, zero)),
                               BetaInd(d, op::Lt(best_y, one)))),
        best_y, half);
    auto best_r = inf;
    // AGGREGATE-init, never default-construction (AGENTS.md; the NEON_BF16
    // implicit-ctor break of 2026-08-06).
    BetaInvFwdOut<D> best_f{Dd<D>{zero, zero}, zero, zero};

    // S1: beta-Temme, offered GLOBALLY [measured deviation, G3]. The
    // generator gates S1's whole candidacy on nu >= kBetaInvS1NuMin, but its
    // own stated reason is about the CORRECTION -- S(zeta,p)/nu is unbounded
    // as nu -> 0 -- and that correction is separately gated inside
    // BetaInvSeedS1, where the reason applies. The UNCORRECTED leading Temme
    // term is an ordinary seed at any nu, and offering it can only help,
    // because a seed that loses the residual comparison costs nothing but the
    // evaluation. Measured: at (20.000015, 2.1492, 4.9e-5), nu = 1.941 sits
    // just under the pinned 2 and S1 scores 0.49 where the best of the other
    // four scores 7.6 -- a factor of fifteen thrown away by a gate that was
    // never about candidacy.
    BetaInvConsider(d, cx, BetaInvSeedS1(d, cx, zq, as, bs, cs, nus), one, mt,
                    &best_y, &best_r, &best_f);

    // S2: series inversion with Picard corrections. Global.
    BetaInvConsider(d, cx, BetaInvSeedS2(d, cx, DdToDouble(ds_num)), one, mt,
                    &best_y, &best_r, &best_f);

    // S4: the exact-B leading-order closed form, branch by sigma vs
    // s* = beta/(alpha+beta) (PLAN's FIRST correction: the exp form, exact in
    // the logit of any sign and magnitude, with c(alpha,beta) DROPPED -- it
    // actively hurt once B is exact). The sigma <= s* branch is S2's own
    // zeroth iterate, i.e. exactly the deep-small quotient already in hand;
    // the other branch is its mirror through 1 - sigma and beta.
    {
      const auto v2 = GammaInvDivD(d, DdAdd(d, lntauc, cx.lb_b), beta);
      const auto y_lo = BetaInvSigmoid(d, DdToDouble(ds_q));
      const auto y_hi = BetaInvSigmoid(d, op::Neg(DdToDouble(v2)));
      // tau vs s* = beta/(alpha+beta), the branch predicate; tau itself is
      // never assembled as a value, so the comparison uses its log.
      const auto lnqm = BetaInvLog(d, cx.q_mean).hi;
      BetaInvConsider(d, cx,
                      op::IfThenElse(op::Ge(lnqm, lntau.hi), y_lo, y_hi), one,
                      mt, &best_y, &best_r, &best_f);
    }

    // S5: logit-normal from EXACT moments. logit(Y) for Y ~ Beta(alpha,beta)
    // is ln X1 - ln X2 for independent Gamma variates, so its mean and
    // variance are psi(alpha)-psi(beta) and psi'(alpha)+psi'(beta) exactly --
    // no asymptotics in either parameter. The normal quantile of tau is
    // -sqrt(2) erfcinv(2 tau), i.e. -sqrt(2) times the same signed z the ridge
    // seed uses.
    {
      const auto mu = op::Sub(DigammaVec(d, alpha), DigammaVec(d, beta));
      const auto var = op::Add(TrigammaVec(d, alpha), TrigammaVec(d, beta));
      const auto zz = op::Mul(op::Set(d, -kBetaInvSqrt2), zq);
      const auto v = op::MulAdd(zz, op::Sqrt(var), mu);
      BetaInvConsider(d, cx, BetaInvSigmoid(d, v), one, mt, &best_y, &best_r,
                      &best_f);
    }

    // S3: the gamma-limit transfer, seeded by gammainv's OWN S1/S2/S3 (each
    // offered as its own candidate and scored by THIS family's residual, which
    // is the only comparison that can decide between them here).
    //
    //   beta >= alpha: t = -beta ln(1-y),  gamma side p, y = -expm1(-t/beta)
    //   beta <  alpha: t = -alpha ln y,    gamma side q, y = exp(-t/alpha)
    //
    // -- mirroring src/beta-inl.h's own gamma-limit slice, with the sides
    // relabelled for the inverse. The whole block is skipped whole-vector when
    // no lane is in the transfer's domain.
    {
      const auto gate3 = BetaInd(d, op::Ge(op::Max(alpha, beta),
                                           op::Set(d, kBetaInvS3Min)));
      if (!op::AllFalse(d, BetaIndMask(d, gate3))) {
        const auto blo = op::Ge(beta, alpha);
        const auto i_blo = BetaInd(d, blo);
        const auto sg = op::IfThenElse(blo, alpha, beta);   // the gamma shape
        const auto hg = op::IfThenElse(blo, beta, alpha);   // the huge side
        // gamma's own p-value is tau when the huge parameter is SECOND (the
        // native limit I_y(alpha,beta) ~ P_gamma(alpha,t)) and 1 - tau when it
        // is first (I_y = 1 - I_{1-y}(beta,alpha) ~ 1 - P_gamma(beta,t)).
        const Dd<D> glnp{op::IfThenElse(blo, lntau.hi, lntauc.hi),
                         op::IfThenElse(blo, lntau.lo, lntauc.lo)};
        const auto lg1s = BetaInvLgamma1p(d, sg);
        const auto lsum = DdToDouble(DdAdd(d, glnp, lg1s));
        // gammainv's seeds solve against the SMALLER of the gamma pair, which
        // is sigma in every combination; its sgn is +1 when that solved side
        // is Q. gamma's P is tau (native) or 1 - tau (swapped), so the solved
        // side is P exactly when the two bits agree.
        const auto i_gp = op::MulAdd(i_blo, i_tis,
                                     op::Mul(BetaIndNot(d, i_blo),
                                             BetaIndNot(d, i_tis)));
        const auto sgn =
            op::IfThenElse(BetaIndMask(d, i_gp), op::Set(d, -1.0), one);
        op::V<D> eta0;
        const auto t1 = GammaInvSeedS1(d, sg, sigma, sgn, &eta0);
        const auto g1 = op::Mul(
            gate3,
            op::Mul(BetaInd(d, op::Ge(sg, op::Set(d,
                                        detail::kGammaInvS1AMin))),
                    BetaInd(d, op::Ge(op::Set(d, detail::kGammaInvEtaMax),
                                      op::Abs(eta0)))));
        const auto t2 = GammaInvSeedS2(d, sg, lsum);
        op::V<D> g3ok;
        const auto lga = LgammaPosDd(d, sg);
        // gamma's S3 is a far-q-tail fixed point taking ln q, and q = sigma
        // wherever it applies -- which is exactly where the gamma-solved side
        // is Q, its own availability condition.
        const auto t3 = GammaInvSeedS3(d, sg, lnsig.hi, lga.hi, &g3ok);
        const auto g3 = op::Mul(op::Mul(gate3, g3ok), BetaIndNot(d, i_gp));

        BetaInvConsider(d, cx, BetaInvGammaToY(d, t1, hg, blo), g1, mt,
                        &best_y, &best_r, &best_f);
        BetaInvConsider(d, cx, BetaInvGammaToY(d, t2, hg, blo), gate3, mt,
                        &best_y, &best_r, &best_f);
        BetaInvConsider(d, cx, BetaInvGammaToY(d, t3, hg, blo), g3, mt,
                        &best_y, &best_r, &best_f);
      }
    }

    // --- safeguarded logit-Newton, kBetaInvStepsN steps ---------------------
    // The gammainv G3 package carried whole, and the step count is PLAN's
    // SECOND correction (StepsN = 4): a bounded interior sub-band
    // (min(alpha,beta) ~ 0.02-0.5, skew 3-10x, y in 0.1-0.3) is short of the
    // gate after all five seed families, whose best there is 2-5 bits; the
    // convergence from that band's worst seed is cleanly quadratic
    // (2.12 -> 6.66 -> 16.48 -> 36.16 -> 75.51 bits) and step 4 clears the
    // gate by 20+ bits band-wide. Nothing is shaved: the safeguard makes the
    // fourth step idempotent for lanes already converged after the third.
    //
    // The step is MULTIPLICATIVE in y, y1 = y(1 + ls), ls = -w (m - m_t),
    // which is the form the generator's replay pinned four steps against. The
    // lower clamp is a floor under a runaway step; the upper side needs none,
    // because a candidate at or past 1 is outside the domain and the
    // acceptance test rejects it -- which is also what makes the backtrack
    // the right response there.
    //
    // EVERY STEP IS SAFEGUARDED: accepted only if it does not increase
    // |m - m_t|, and a rejected step is retried an eighth as long. Within
    // kBetaInvTrustResid of the solution the step is taken on trust instead,
    // because there the two residuals differ by the forward's own noise
    // rather than by anything about the step. Freeze by SELECT throughout,
    // never by adding a zero step.
    //
    // If ANY lane failed to select a candidate, its carried-over forward is
    // the zero-initialized one, which would make the first step read a
    // residual it never computed; recompute for the whole vector in that case.
    auto yn = best_y;
    BetaInvFwdOut<D> cur = best_f;
    auto rbest = best_r;
    if (!op::AllFalse(d, op::Eq(best_r, inf))) {
      cur = BetaInvForward(d, cx, yn);
      rbest = op::Abs(BetaInvResid(d, cur, mt));
    }
    auto scale = one;
    for (int k = 0; k < detail::kBetaInvStepsN; ++k) {
      const auto resid = BetaInvResid(d, cur, mt);
      auto ls = op::Neg(op::Mul(op::Mul(scale, cur.w), resid));
      ls = op::Max(ls, op::Set(d, -0.9));
      // FREEZE ON NOISE. A residual at or below the forward's own uncertainty
      // carries no information about where the root is, and stepping on it
      // moves the iterate by w times pure noise -- which for a large condition
      // number is a large move in the WRONG direction as often as the right
      // one, and the trust bypass below would accept every one of them. This
      // is the one place the safeguard package needed a beta-specific addition
      // (gammainv's w is bounded by ~2^10, so its noise steps are invisible;
      // beta's reaches 2^50 in the joint-tiny band).
      ls = op::IfThenElse(op::Ge(cur.unc, op::Abs(resid)), zero, ls);
      const auto cand = op::MulAdd(yn, ls, yn);

      const auto f = BetaInvForward(d, cx, cand);
      const auto rnew = op::Abs(BetaInvResid(d, f, mt));
      const auto trust =
          BetaInd(d, op::Gt(op::Set(d, detail::kBetaInvTrustResid), rbest));
      const auto acc = BetaIndMask(
          d, op::Mul(op::Mul(BetaInd(d, op::Gt(cand, zero)),
                             BetaInd(d, op::Lt(cand, one))),
                     op::Max(trust, BetaInd(d, op::Ge(rbest, rnew)))));
      yn = op::IfThenElse(acc, cand, yn);
      rbest = op::IfThenElse(acc, rnew, rbest);
      cur.m = Dd<D>{op::IfThenElse(acc, f.m.hi, cur.m.hi),
                    op::IfThenElse(acc, f.m.lo, cur.m.lo)};
      cur.w = op::IfThenElse(acc, f.w, cur.w);
      cur.unc = op::IfThenElse(acc, f.unc, cur.unc);
      scale = op::IfThenElse(acc, one, op::Mul(scale, op::Set(d, 0.125)));
    }
    y = op::IfThenElse(m_ds, y_ds, yn);
  }

  // --- back to the caller's variable ---------------------------------------
  // x = 1 (-) y on the swapped orientation, ONE rounding -- which is the best
  // a double x near 1 can carry, and exactly why the public header documents
  // the swap identity as the lossless route for the other end.
  auto res = op::IfThenElse(m_swap, op::Sub(one, y), y);

  // --- specials (beta's own table, read as quantiles) -----------------------
  // Applied last, in increasing priority. One degenerate parameter puts all
  // the mass at an endpoint and the quantile is that endpoint; two
  // degeneracies give NaN, as does a negative parameter or an out-of-range
  // probability.
  const auto at0 = kQ ? one : zero;  // s = 0
  const auto at1 = kQ ? zero : one;  // s = 1
  res = op::IfThenElse(op::Eq(s_in, zero), at0, res);
  res = op::IfThenElse(op::Eq(s_in, one), at1, res);

  const auto a0 = op::Eq(a_in, zero);
  const auto b0 = op::Eq(b_in, zero);
  const auto ai = op::Eq(a_in, inf);
  const auto bi = op::Eq(b_in, inf);
  // Mass at 0 (a = 0 or b = +inf): every quantile is 0. Mass at 1 (b = 0 or
  // a = +inf): every quantile is 1.
  res = op::IfThenElse(a0, zero, res);
  res = op::IfThenElse(bi, zero, res);
  res = op::IfThenElse(b0, one, res);
  res = op::IfThenElse(ai, one, res);
  {
    const auto da = op::Max(BetaInd(d, a0), BetaInd(d, ai));
    const auto db = op::Max(BetaInd(d, b0), BetaInd(d, bi));
    res = op::IfThenElse(op::Gt(op::Mul(da, db), half), qnan, res);
  }
  res = op::IfThenElse(op::Lt(a_in, zero), qnan, res);
  res = op::IfThenElse(op::Lt(b_in, zero), qnan, res);
  res = op::IfThenElse(op::Lt(s_in, zero), qnan, res);
  res = op::IfThenElse(op::Gt(s_in, one), qnan, res);
  res = op::IfThenElse(op::IsNaN(a_in), a_in, res);  // payload preserved
  res = op::IfThenElse(op::IsNaN(b_in), b_in, res);
  res = op::IfThenElse(op::IsNaN(s_in), s_in, res);
  return res;
}

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#endif  // include guard
