// Inverse regularized incomplete gamma: gamma_p_inv(a, p) and
// gamma_q_inv(a, q). Per-target include guard (Highway -inl.h idiom).
//
// Both public functions are ONE pipeline with one bit of orientation, per
// PLAN.md "P1 inverse incomplete gamma -- detail design" (BINDING) and the
// parameters replay-pinned in src/gammainv_data.h. The forward region cores
// (src/gamma-inl.h) are consumed, never re-derived; this file adds the seed
// stage, the residual stage and the routing that mirrors the forward's own
// region map.
//
// THE INPUT-SIDE FLIP. Whatever the caller asked for, the kernel solves for
// the side whose probability is <= 1/2: an input s > 1/2 is replaced by
// 1 (-) s, which is EXACT by Sterbenz there, and the orientation bit is
// flipped. So the target t is always in (0, 1/2] and the answer never comes
// from a subtraction. This is the inverse's whole advantage over the forward
// pair, whose complement rounding is on the OUTPUT and cannot be undone.
//
// EVERYTHING IS DONE IN LOG SPACE, AND THE RESIDUAL IS A LOGIT. Newton is run
// on
//     m(x) = log P(x) - log Q(x)   against   m_t = log p - log q,
// never on F(x) - t, and never on log F alone. Four separate things force
// this, and any one of them would:
//   * log t reaches -745. In double, log F (-) log t would carry an absolute
//     error of 2^-53*745 ~ 2^-43, which lands directly on x through the
//     step; taken as a dd difference it is 2^-95. (This is the same
//     cancellation erfinv's FAR branch solves the same way.)
//   * F itself underflows to zero over a huge part of the domain -- every
//     far tail with |E| > 745 -- while log F stays an ordinary number. A
//     residual-space kernel is blind there; a log-space one is not.
//   * The solved side tends to 1 over the wrong half of the domain, so its
//     log tends to 0 and every point out there scores alike. The logit's two
//     limbs cover the two halves between them and it saturates nowhere.
//   * |log min(P,Q)|, signed, would do that much but JUMPS by 2*log 2 at the
//     median -- exactly where a target of 1/2 puts its root. The logit is
//     continuous there.
// The step is x1 = x*(1 + ls), ls = -w*(m - m_t), with
//     w = exp(log u (-) E) * (1 - u),   u = min(P, Q),
//     E = log(x^a e^-x / Gamma(a)),
// because dm/d(log x) = x*g/(P*Q) and P*Q = u*(1-u). Note w is the inverse's
// own condition number |d log x / d log s| up to the O(1) factor (1-u);
// PLAN's conditioning adjudication bounds it by ~2^10 for every input whose
// true x is a normal double, and it goes to ZERO in the huge-a
// beyond-resolution regime -- which is exactly why no huge-a branch exists.
//
// E IS COMPUTED TWO WAYS, SPLIT AT a_T. Below a_T, directly:
// E = a*log(x) (-) x (-) lgamma(a), whose terms stay small enough that the dd
// cancellation is harmless. At and above a_T the direct form is USELESS --
// at a = 1e300 its terms are ~7e302 while E itself is O(100) on the ridge, so
// the dd difference has an absolute error of 1e271 -- and the Stirling form
//     E = -a*phi(lambda) + 1/2*log(a/(2*pi)) - mu(a),   phi = lambda-1-log lambda,
// is used instead: a*phi from Log1pmxDd (exactly as the forward's R3 does),
// and mu the Stirling remainder, which is lgamma's own kLgammaStirCoef tail
// (valid from kLgammaX0 = 8, and a_T = 20). This is also what makes the
// beyond-resolution rows work at all: the direct form's a*log(x) overflows to
// inf for a above ~2.5e305 and its residual to NaN, while the Stirling form
// at x = a gives phi = 0 and E = 1/2*log(a/2pi) exactly.
//
// SEEDS: THREE CANDIDATES, ONE RESIDUAL COMPARISON (PLAN.md FIRST
// CORRECTION -- the partition is by (side, lambda-regime) at ALL a, not by a
// alone). S1 is Temme's normal-quantile seed (erfcinv -> eta -> lambda(eta),
// plus the two pinned eps_k/a corrections); S2 is the p-form small-x seed
// exp((ln p + lgamma(1+a))/a) with kGammaInvS2NCorr Picard corrections, whose
// iteration is the EXACT fixed point of the defining equation and so is the
// only seed that recovers the tiny-a corner; S3 is the far-q-tail fixed point
// under its own stability gate. Each is evaluated where its domain guard
// passes, and the one with the smallest |m(x0) - m_t| wins -- the same
// objective the steps then reduce, so seed choice and refinement cannot
// disagree about what "closer" means. The winner's forward evaluation is
// REUSED as the first step's, so three candidates cost three forward
// evaluations and the three steps cost three more.
//
// DEEP-SMALL CLOSED FORM. Where the correction factor e^-x * Sum(a,x) is
// below the resolution of a double, the equation collapses to
// x = exp((ln p + lgamma(1+a))/a), evaluated with exp_dd's mantissa/exponent
// split so the power-of-two scaling is the LAST operation and a subnormal (or
// zero) answer carries exactly one rounding. That branch owns the entire
// tiny-a collapse zone and the whole subnormal-x band.
//
// ln p, NEVER log of 1 (-) tiny. On the q orientation the p-form seeds and the
// closed form need ln p = ln(1 - q). q is not tiny there by construction
// (q <= 1/2), but p can be within 2^-1000 of 1, so 1 (-) q is formed as an
// EXACT dd (TwoSum) and LogDdAny takes the pair: the low word is what carries
// the entire signal when fl(1 - q) rounds to 1.
#if defined(CORVUS_GAMMAINV_INL_H_) == defined(HWY_TARGET_TOGGLE)
#ifdef CORVUS_GAMMAINV_INL_H_
#undef CORVUS_GAMMAINV_INL_H_
#else
#define CORVUS_GAMMAINV_INL_H_
#endif

#include <limits>

#include "src/dd-inl.h"
#include "src/dd_special-inl.h"
#include "src/erfinv-inl.h"
#include "src/exp_dd-inl.h"
#include "src/gamma-inl.h"
#include "src/gamma_data.h"
#include "src/gammainv_data.h"
#include "src/lgamma-inl.h"
#include "src/lgamma_data.h"
#include "src/log_dd-inl.h"
#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

namespace op = ops;

// ------------------------------------------------------------------------
// OUTLINED log and exp [MSVC BUILD-TIME GATE, AGENTS.md]. These are thin
// wrappers whose only purpose is the HWY_NOINLINE: log_dd (via LogDdAny),
// exp_dd (via ExpDd) and Highway's own Exp are each large, and this file
// reaches them from about seventeen call sites -- the lambda-of-eta Newton,
// both S2 Picard loops, S3's fixed point, the forward's E/region assembly
// and the logit swap. Inlined, each of those becomes its own copy of a
// table gather plus a polynomial IN EVERY ONE OF THE COMPILED TARGETS, and
// cl.exe's optimizer is superlinear in function size: this is the same fix
// that took betainv.cpp (src/betainv-inl.h) from past 45 minutes and 7 GB
// on one MSVC invocation to 127 s, retro-applied here per PLAN.md's "MSVC
// build-time headroom" item. Bit-identity is guaranteed by contraction-off
// and verified by byte-comparing the ULP tables across the change.
template <class D>
HWY_NOINLINE Dd<D> GammaInvLog(D d, Dd<D> x) {
  return LogDdAny(d, x);
}
template <class D>
HWY_NOINLINE Dd<D> GammaInvLog(D d, op::V<D> x) {
  return GammaInvLog(d, Dd<D>{x, op::Zero(d)});
}
template <class D>
HWY_NOINLINE Dd<D> GammaInvExpDd(D d, op::V<D> xh, op::V<D> xl) {
  return ExpDd(d, xh, xl);
}
template <class D>
HWY_NOINLINE op::V<D> GammaInvExp(D d, op::V<D> x) {
  return op::Exp(d, x);
}

// The pin emits the two depth buckets separately so a future re-measure can
// split them without a header-format change; this kernel is the single-variant
// reading of "both true". If that ever stops holding, the assert is what makes
// it a build failure instead of a silent accuracy regression.
static_assert(detail::kGammaInvStepsLogResidualDeep &&
                  detail::kGammaInvStepsLogResidualShallow,
              "gammainv-inl.h implements the log-residual step for both depth "
              "buckets; a pin that splits the variants needs a kernel change");

// --- kernel-internal constants ------------------------------------------
// Closed-form, not fitted, so they live here rather than in the generated
// header -- the same split erfinv makes for kErfinvTwoOverSqrtPi.

// log(2*pi), as a dd. Halved at the use site (exact).
inline constexpr double kGammaInvLog2PiHi = 0x1.d67f1c864beb5p+0;
inline constexpr double kGammaInvLog2PiLo = -0x1.65b5a1b7ff5dfp-54;

// Freeze threshold for the seed's plain-double series (a seed, not a
// precision-bearing sum -- the forward's own kGammaFreezeEps is 2^-60).
inline constexpr double kGammaInvSeedEps = 0x1.0p-60;

// Ceiling on |S/a| in the deep-small closed form before the quotient is
// clamped. exp underflows at -745 and ExpDdFrac clamps its own argument at
// 1100, so any |S/a| past this produces the same exact 0 (or, for the sign
// that cannot occur with a tiny a, the same inf); clamping keeps the dd
// refinement from ever seeing an infinite first quotient.
inline constexpr double kGammaInvQuotMax = 0x1.0p+11;

// a below this has 1/a overflowing or losing its low word, so the deep-small
// quotient rescales both operands by an exact power of two first. 2^-900
// leaves room for the smallest subnormal a (2^-1074) to land at 2^-874.
inline constexpr double kGammaInvTinyA = 0x1.0p-900;
inline constexpr double kGammaInvTinyAScale = 0x1.0p+200;

// Residual below which a Newton step is taken on trust rather than tested
// against the previous one. See the step loop.
inline constexpr double kGammaInvTrustResid = 0.5;

// Scrub value for the CF's a: j*(j-a) overflows to -inf for a above ~4e306.
// No R2 lane with a that large is ever live (a*phi is hundreds of thousands
// of e-foldings past saturation), so this only keeps discarded lanes finite.
inline constexpr double kGammaInvCfAMax = 0x1.0p+996;

// (num.hi + num.lo) / den as a dd, for a den that may be subnormal and a
// quotient that may overflow.
//
// The classical two-step division: q0 = fl(num.hi/den), then the exact
// residual num (-) q0*den (TwoProd supplies q0*den exactly) divided again.
// Two guards make it total. The power-of-two RESCALE of both operands is
// exact and leaves the quotient unchanged, so a subnormal den -- routine
// here, since a is allowed to be subnormal -- never reaches the division as
// an operand whose reciprocal is infinite. The quotient CLAMP catches the
// case the rescale cannot: |num/den| genuinely exceeding the double range,
// which happens whenever the answer underflows to zero, and which would
// otherwise put an infinity into TwoProd and a NaN into the result. A clamped
// lane drops its refinement term (it is meaningless) and exponentiates to the
// same exact zero the unclamped value would have.
template <class D>
HWY_INLINE Dd<D> GammaInvDivD(D d, Dd<D> num, op::V<D> den) {
  const auto zero = op::Zero(d);
  const auto tiny = op::Lt(den, op::Set(d, kGammaInvTinyA));
  const auto sc =
      op::IfThenElse(tiny, op::Set(d, kGammaInvTinyAScale), op::Set(d, 1.0));
  const auto b = op::Mul(den, sc);
  const Dd<D> n{op::Mul(num.hi, sc), op::Mul(num.lo, sc)};

  const auto lim = op::Set(d, kGammaInvQuotMax);
  auto q0 = op::Div(n.hi, b);
  const auto over = op::Gt(op::Abs(q0), lim);
  q0 = op::IfThenElse(over, op::CopySign(lim, q0), q0);
  const auto pr = TwoProd(d, q0, b);
  const auto rem = DdAdd(d, n, Dd<D>{op::Neg(pr.hi), op::Neg(pr.lo)});
  const auto q1 = op::IfThenElse(over, zero, op::Div(DdToDouble(rem), b));
  return Fast2Sum(d, q0, q1);  // |q0| >> |q1| by construction
}

// ------------------------------------------------------------------------
// lambda(eta): the inverse of 1/2 eta^2 = lambda - 1 - log(lambda), with the
// branch chosen by sign(eta). Scheme and both pinned parameters come from
// tools/gen_gammainv_data.py self-checks (a)/(b): a Taylor series (derived
// there by exact-rational Lagrange reversion of the defining equation) for
// |eta| < kGammaInvEtaSeriesCut, and Newton in u = log(lambda) elsewhere,
// from the sign-dependent asymptotic guess that is exact to leading order in
// each of phi's two regimes. Both branches are plain double: this is a seed,
// and the generator measured its floor at 32.6 bits.
template <class D>
HWY_NOINLINE op::V<D> GammaInvLamOfEta(D d, op::V<D> eta) {
  const auto one = op::Set(d, 1.0);
  const auto zero = op::Zero(d);

  const auto* c = detail::kGammaInvLamSeriesCoef;
  auto s = op::Set(d, c[detail::kGammaInvLamSeriesOrder - 1]);
  for (int k = detail::kGammaInvLamSeriesOrder - 2; k >= 0; --k) {
    s = op::MulAdd(s, eta, op::Set(d, c[k]));
  }
  const auto lam_ser = op::MulAdd(eta, s, one);

  const auto m_ser = op::Lt(op::Abs(eta), op::Set(d, detail::kGammaInvEtaSeriesCut));
  if (op::AllTrue(d, m_ser)) return lam_ser;

  const auto target = op::Mul(op::Set(d, 0.5), op::Mul(eta, eta));
  // log1p(target) for the eta >= 0 branch: target >= 1/8 there, so the plain
  // log of 1 + target loses nothing a Newton start could notice.
  auto u = op::IfThenElse(op::Ge(eta, zero), LogDd(d, op::Add(one, target)).hi,
                          op::Neg(target));
  for (int i = 0; i < detail::kGammaInvLamNewtonIters; ++i) {
    const auto eu = GammaInvExp(d, u);
    const auto em1 = op::Sub(eu, one);
    const auto v = op::Sub(op::Sub(em1, u), target);
    // em1 == 0 only at u == 0, which is the series branch's own domain; the
    // select keeps a discarded lane from turning it into an infinity.
    u = op::IfThenElse(op::Eq(em1, zero), u, op::Sub(u, op::Div(v, em1)));
  }
  return op::IfThenElse(m_ser, lam_ser, GammaInvExp(d, u));
}

// c_k(eta0) from the pinned Chebyshev rows, Clenshaw in the mapped variable.
// The map is clamped to [-1, 1]: a no-op on every lane the S1 guard admits
// (|eta0| <= kGammaInvEtaMax), and it stops a discarded lane's huge eta0 from
// running the recurrence up to inf - inf = NaN.
template <class D>
HWY_INLINE op::V<D> GammaInvCk(D d, op::V<D> t, int k) {
  const double* row = detail::kGammaInvCkCheb[k];
  const auto t2 = op::Add(t, t);
  auto b1 = op::Zero(d);
  auto b2 = op::Zero(d);
  for (int j = detail::kGammaInvCkNCoef - 1; j >= 1; --j) {
    const auto nb = op::Sub(op::MulAdd(t2, b1, op::Set(d, row[j])), b2);
    b2 = b1;
    b1 = nb;
  }
  return op::Sub(op::MulAdd(t, b1, op::Set(d, row[0])), b2);
}

// S1: Temme's normal-quantile seed. z = erfcinv(2t) (2t is exact and lands in
// (0, 1], so z >= 0), eta0 = sgn*z*sqrt(2/a) with sgn = +1 when the solved
// side is Q, then the pinned eps_k(eta)/a^k correction to eta itself --
// side-SYMMETRIC, per the derivation in the generator -- and x0 = a*lambda.
template <class D>
HWY_NOINLINE op::V<D> GammaInvSeedS1(D d, op::V<D> a, op::V<D> t,
                                     op::V<D> sgn, op::V<D>* eta0_out) {
  const auto z = ErfcinvVec(d, op::Add(t, t));
  const auto eta0 =
      op::Mul(op::Mul(sgn, z), op::Sqrt(op::Div(op::Set(d, 2.0), a)));
  *eta0_out = eta0;

  const auto one = op::Set(d, 1.0);
  auto tc = op::Mul(eta0, op::Set(d, 1.0 / detail::kGammaInvEtaMax));
  tc = op::Max(op::Min(tc, one), op::Neg(one));
  const auto ra = op::Div(one, a);
  auto sc = op::Zero(d);
  for (int k = detail::kGammaInvSeedNCorr - 1; k >= 0; --k) {
    sc = op::MulAdd(sc, ra, GammaInvCk(d, tc, k));
  }
  const auto eta = op::MulAdd(sc, ra, eta0);
  return op::Mul(a, GammaInvLamOfEta(d, eta));
}

// Sum(a, x) = sum_{n>=0} x^n / ((a+1)...(a+n)), plain double, for the S2
// Picard correction only. Same shape as the forward's GammaSeriesSum (same
// cap, same clamp of x to the region's own a+1 bound so a discarded lane's
// enormous x cannot run the terms up to inf and never freeze), but the
// accumulator is a double: this feeds a seed, and its own rounding is far
// below the seed's measured floor.
template <class D>
HWY_INLINE op::V<D> GammaInvSeedSum(D d, op::V<D> a, op::V<D> x) {
  const auto one = op::Set(d, 1.0);
  const auto zero = op::Zero(d);
  const auto half = op::Set(d, 0.5);
  const auto eps = op::Set(d, kGammaInvSeedEps);
  const auto xc = op::Min(x, op::Add(a, one));

  auto t = one;
  auto s = one;
  auto live = one;
  for (int n = 1; n <= detail::kGammaSeriesN; ++n) {
    t = op::Mul(t, op::Div(xc, op::Add(a, op::Set(d, static_cast<double>(n)))));
    const auto lm = op::Gt(live, half);
    s = op::IfThenElse(lm, op::Add(s, t), s);
    live = op::IfThenElse(op::Lt(t, op::Mul(s, eps)), zero, live);
    if (op::AllTrue(d, op::Eq(live, zero))) break;
  }
  return s;
}

// S2: the p-form seed x0 = exp((ln p + lgamma(1+a))/a), refined by
// kGammaInvS2NCorr Picard steps of
//     log x = (ln p + lgamma(1+a) + x - log Sum(a,x)) / a,
// which is the EXACT fixed point of p = x^a e^-x Sum / Gamma(1+a) -- not an
// approximation being iterated, which is why this seed alone closes the
// tiny-a corner where the bare closed form is off by its whole dropped
// factor. lsum is (ln p + lgamma(1+a)) rounded once; the loop is double
// throughout, as measured.
template <class D>
HWY_NOINLINE op::V<D> GammaInvSeedS2(D d, op::V<D> a, op::V<D> lsum) {
  const auto zero = op::Zero(d);
  // The finiteness test is against infinity, never against DBL_MAX: a is
  // allowed to BE DBL_MAX, and so therefore is the answer.
  const auto dmax = op::Set(d, std::numeric_limits<double>::infinity());
  auto x = GammaInvExp(d, op::Div(lsum, a));
  for (int i = 0; i < detail::kGammaInvS2NCorr; ++i) {
    const auto sum = GammaInvSeedSum(d, a, x);
    const auto lx = op::Div(
        op::Sub(op::Add(lsum, x), LogDd(d, sum).hi), a);
    const auto xn = GammaInvExp(d, lx);
    // The generator breaks out of the loop on a non-finite or non-positive
    // iterate; per lane, that is holding the previous value. Written through
    // the indicator arithmetic because the ops facade exposes no mask AND.
    const auto good = IndMask(d, op::Mul(Ind(d, op::Gt(xn, zero)),
                                         Ind(d, op::Lt(xn, dmax))));
    x = op::IfThenElse(good, xn, x);
  }
  return x;
}

// S3: the far-q-tail fixed point x <- L + (a-1)*log(x), L = -ln(q*Gamma(a)),
// under the pinned stability gate L > kGammaInvS3StabilityMargin*|a-1| (the
// map's local contraction factor is (a-1)/x, so L > 0 alone is not
// sufficient -- the generator's own oscillation finding). Returns the
// candidate and, in *ok_out, the indicator saying it is inside its domain.
template <class D>
HWY_NOINLINE op::V<D> GammaInvSeedS3(D d, op::V<D> a, op::V<D> lnq,
                                     op::V<D> lga, op::V<D>* ok_out) {
  const auto one = op::Set(d, 1.0);
  const auto zero = op::Zero(d);
  const auto dmax = op::Set(d, std::numeric_limits<double>::infinity());
  const auto L = op::Sub(op::Neg(lnq), lga);
  const auto gate =
      op::Mul(Ind(d, op::Gt(L, op::Mul(op::Set(d,
                                               detail::kGammaInvS3StabilityMargin),
                                       op::Abs(op::Sub(a, one))))),
              Ind(d, op::Lt(L, dmax)));

  auto x = L;
  for (int i = 0; i < detail::kGammaInvS3NIter; ++i) {
    const auto xn = op::MulAdd(op::Sub(a, one), GammaInvLog(d, x).hi, L);
    const auto good = IndMask(d, op::Mul(Ind(d, op::Gt(xn, zero)),
                                         Ind(d, op::Lt(xn, dmax))));
    x = op::IfThenElse(good, xn, x);
  }
  *ok_out = gate;
  return x;
}

// ------------------------------------------------------------------------
// The forward, in log space. Returns the kernel's OBJECTIVE and its slope --
// never a probability, which is what keeps the whole underflow range live.
//
// THE OBJECTIVE IS THE LOGIT, m(x) = log P(x) - log Q(x). Three properties,
// and the kernel needs all three:
//   * MONOTONE and unbounded in both directions: m -> -inf as x -> 0 and
//     +inf as x -> inf, with no saturation anywhere. Scoring a candidate by
//     |log F_solve - log t| instead LOOKS equivalent and is not -- the solved
//     side tends to 1 over the wrong half of the domain, so its log tends to
//     0 and every point out there scores the same |log t|. That is not a
//     corner case: it is what ranked a candidate at 0.79*a ahead of a correct
//     one at a + 1 ulp (a = 1.9e34, q = 7.7e-53; measured 2e18 ULP), because
//     the solved Q at 0.79*a really is 1 to the last bit while log P there is
//     -4.8e32 and says exactly how wrong the point is.
//   * CONTINUOUS at the median, which |log min(P,Q)| with a sign is not: that
//     one jumps by 2*log 2 exactly where P = Q, i.e. exactly where a target
//     of 1/2 puts its root, and Newton then oscillates across the jump
//     (measured 5e14 ULP at s = 1/2 before the logit).
//   * Its two limbs are the two sides' own logs, so it degrades gracefully to
//     log t at either extreme and costs no accuracy at the median, where
//     m -> 0 while both limbs are -log 2 and the dd subtraction is exact.
// Its slope is dm/d(log x) = x*g*(1/P + 1/Q) = x*g/(P*Q), so the reciprocal
// the Newton step wants is w = P*Q/(x*g) = exp(log u - E) * (1 - u) with u
// the smaller side -- a product of two factors each of which stays in range
// however far into a tail the iterate sits.
//
// Whichever of the pair the region computed directly, the SMALLER one is
// whichever of {v, 1 - v} is <= 1/2, so the test is one comparison against
// log(1/2), and the complement is only ever formed FROM the larger one --
// never a subtraction that loses anything.
template <class D>
struct GammaInvFwdOut {
  Dd<D> m;     // log P - log Q
  op::V<D> w;  // d(log x)/dm at this point
};

// HWY_NOINLINE from day one, like the region cores it calls: this function is
// reached from six call sites per export (three candidate seeds, three steps)
// and, even with log_dd/exp_dd routed through GammaInvLog/GammaInvExpDd/
// GammaInvExp above, still inlines Log1pmxDd, the erfc core and all four
// region cores.
template <class D>
HWY_NOINLINE GammaInvFwdOut<D> GammaInvForward(D d, op::V<D> a, op::V<D> x,
                                               op::V<D> i_wp, Dd<D> lga,
                                               Dd<D> lna) {
  const auto one = op::Set(d, 1.0);
  const auto zero = op::Zero(d);
  const auto at = op::Set(d, detail::kGammaAT);
  const auto m_big = op::Ge(a, at);

  // --- E, split at a_T (see the file header) ------------------------------
  Dd<D> e{zero, zero};
  Dd<D> aphi{zero, zero};
  const bool need_small = !op::AllTrue(d, m_big);
  const bool need_big = !op::AllFalse(d, m_big);
  Dd<D> e_small{zero, zero};
  Dd<D> e_big{zero, zero};
  if (need_small) {
    const auto lx = GammaInvLog(d, x);
    e_small = DdAddD(d, DdMulD(d, lx, a), op::Neg(x));
    e_small = DdAdd(d, e_small, Dd<D>{op::Neg(lga.hi), op::Neg(lga.lo)});
  }
  if (need_big) {
    // a is clamped up to a_T for the lanes this branch is evaluated on but
    // does not own: 1/a and the Stirling remainder both need a bounded away
    // from zero, and a subnormal a among the discarded lanes would otherwise
    // reach them.
    const auto ac = op::Max(a, at);
    const auto u = DdMul(d, TwoSum(d, x, op::Neg(ac)), DdRecip(d, ac));
    aphi = DdMulD(d, Log1pmxDd(d, u), ac);
    // mu(a), the Stirling remainder, is lgamma's own tail polynomial in 1/a^2
    // divided by a -- the same coefficients LgammaStirling evaluates, valid
    // from kLgammaX0 = 8 and so over all of a >= a_T = 20.
    const auto ra2 = op::Div(one, op::Mul(ac, ac));
    auto pmu = op::Set(d, detail::kLgammaStirCoef[detail::kLgammaStirNCoef - 1]);
    for (int k = detail::kLgammaStirNCoef - 2; k >= 0; --k) {
      pmu = op::MulAdd(pmu, ra2, op::Set(d, detail::kLgammaStirCoef[k]));
    }
    const auto half = op::Set(d, 0.5);
    // 1/2*log(a) - 1/2*log(2*pi): both halvings exact.
    Dd<D> hl{op::Mul(lna.hi, half), op::Mul(lna.lo, half)};
    hl = DdAdd(d, hl,
               Dd<D>{op::Set(d, -0.5 * kGammaInvLog2PiHi),
                     op::Set(d, -0.5 * kGammaInvLog2PiLo)});
    e_big = DdAdd(d, Dd<D>{op::Neg(aphi.hi), op::Neg(aphi.lo)}, hl);
    e_big = DdAddD(d, e_big, op::Neg(op::Div(pmu, ac)));
  }
  e = Dd<D>{op::IfThenElse(m_big, e_big.hi, e_small.hi),
            op::IfThenElse(m_big, e_big.lo, e_small.lo)};

  // --- region map: GammaVec's, with kP carried per lane --------------------
  const auto i_small = Ind(d, op::Lt(a, at));
  const auto i_xle = Ind(d, op::Ge(op::Add(a, one), x));
  const auto i_lo = Ind(d, op::Ge(a, op::Add(x, x)));
  const auto i_hi = Ind(d, op::Ge(x, op::Add(a, a)));
  const auto i_box = op::Mul(Ind(d, op::Ge(op::Set(d, kGammaR4AMax), a)),
                             Ind(d, op::Ge(op::Set(d, kGammaR4XMax), x)));
  const auto i_wq = IndNot(d, i_wp);
  const auto i_bigr = IndNot(d, i_small);
  auto i_r1 = op::MulAdd(i_small, i_xle, op::Mul(i_bigr, i_lo));
  auto i_r2 = op::MulAdd(i_small, IndNot(d, i_xle), op::Mul(i_bigr, i_hi));
  const auto i_r3 = op::Mul(i_bigr, op::Mul(IndNot(d, i_lo), IndNot(d, i_hi)));
  // R4 takes the whole box when the solved side is Q, and only the x > a+1
  // part when it is P -- exactly gamma_q's and gamma_p's own rule.
  const auto i_r4 = op::Mul(i_box, IndNot(d, op::Mul(i_wp, i_xle)));
  i_r1 = op::Mul(i_r1, IndNot(d, op::Mul(i_wq, i_r4)));
  i_r2 = op::Mul(i_r2, IndNot(d, i_r4));

  const auto m_r1 = IndMask(d, i_r1);
  const auto m_r2 = IndMask(d, i_r2);
  const auto m_r3 = IndMask(d, i_r3);
  const auto m_r4 = IndMask(d, i_r4);

  // --- ln of the directly computed side ------------------------------------
  Dd<D> lnf{zero, zero};
  auto i_pdir = zero;
  if (!op::AllFalse(d, m_r1)) {
    // P = e^{E - log a} * Sum, so ln P needs no exponential at all.
    const auto l = DdAdd(d, DdAdd(d, e, Dd<D>{op::Neg(lna.hi), op::Neg(lna.lo)}),
                         LogDd(d, GammaSeriesSum(d, a, x)));
    lnf = Dd<D>{op::IfThenElse(m_r1, l.hi, lnf.hi),
                op::IfThenElse(m_r1, l.lo, lnf.lo)};
    i_pdir = op::IfThenElse(m_r1, one, i_pdir);
  }
  if (!op::AllFalse(d, m_r2)) {
    const auto as = op::Min(a, op::Set(d, kGammaInvCfAMax));
    const auto l = DdAdd(d, e, GammaInvLog(d, GammaCfRecip(d, as, x)));
    lnf = Dd<D>{op::IfThenElse(m_r2, l.hi, lnf.hi),
                op::IfThenElse(m_r2, l.lo, lnf.lo)};
  }
  if (!op::AllFalse(d, m_r3)) {
    const auto t3 = GammaTemme(d, a, x);
    // Past a*phi ~ 745 the assembled Temme value underflows to zero, and its
    // log would then be the log of the smallest normal rather than a number
    // with any bearing on the residual. -a*phi is the leading term of ln F
    // there and is exact to the dd; using it keeps the step's SIGN and
    // magnitude right, which is what pulls such a lane back toward the root
    // instead of stranding it. (No live lane is ever there: a target below
    // e^-745 is not a representable double.)
    const auto lg3 = GammaInvLog(d, t3.val.dd);
    const auto pos = op::Gt(t3.val.dd.hi, zero);
    const auto l = Dd<D>{op::IfThenElse(pos, lg3.hi, op::Neg(aphi.hi)),
                         op::IfThenElse(pos, lg3.lo, op::Neg(aphi.lo))};
    lnf = Dd<D>{op::IfThenElse(m_r3, l.hi, lnf.hi),
                op::IfThenElse(m_r3, l.lo, lnf.lo)};
    i_pdir = op::IfThenElse(m_r3, t3.is_p, i_pdir);
  }
  if (!op::AllFalse(d, m_r4)) {
    const auto l = GammaInvLog(d, GammaSmallQ(d, a, x));
    lnf = Dd<D>{op::IfThenElse(m_r4, l.hi, lnf.hi),
                op::IfThenElse(m_r4, l.lo, lnf.lo)};
    i_pdir = op::IfThenElse(m_r4, zero, i_pdir);
  }

  // --- assemble the logit --------------------------------------------------
  // First put the SMALLER side in hand: v > 1/2 <=> log v > -log 2. The
  // region map picks the smaller side by construction, so the swap below
  // fires only near a region's own switch (and on iterates far from the
  // root), and there the value being complemented is >= ~0.2 -- a subtraction
  // costing a fraction of a bit. The value is recovered by exponentiating its
  // log rather than assembled a second time: on that branch |log v| <= 0.7,
  // so the exponential's relative error is the dd's own.
  //
  // Then the larger side is 1 (-) u, whose log completes m = +-(log u -
  // log(1-u)). For a u that has underflowed to zero this is log 1 = 0
  // exactly, which is the correct limit and not a fudge: the whole logit is
  // then log u, and log u is the limb that was never formed as a value.
  {
    const auto ln_half = op::Set(d, -0.6931471805599453);
    const auto m_cmp = op::Gt(lnf.hi, ln_half);
    if (!op::AllFalse(d, m_cmp)) {
      const auto vb = GammaInvExpDd(d, lnf.hi, lnf.lo);
      const auto cb =
          DdAdd(d, Dd<D>{one, zero}, Dd<D>{op::Neg(vb.hi), op::Neg(vb.lo)});
      const auto l = GammaInvLog(d, cb);
      lnf = Dd<D>{op::IfThenElse(m_cmp, l.hi, lnf.hi),
                  op::IfThenElse(m_cmp, l.lo, lnf.lo)};
      i_pdir = op::IfThenElse(m_cmp, IndNot(d, i_pdir), i_pdir);
    }
  }
  const auto u = GammaInvExpDd(d, lnf.hi, lnf.lo);  // u = min(P, Q) <= 1/2
  const auto c = DdAdd(d, Dd<D>{one, zero}, Dd<D>{op::Neg(u.hi), op::Neg(u.lo)});
  const auto lnc = GammaInvLog(d, c);
  const auto lg = DdAdd(d, lnf, Dd<D>{op::Neg(lnc.hi), op::Neg(lnc.lo)});
  const auto m_p = IndMask(d, i_pdir);  // is the smaller side P?
  const Dd<D> m{op::IfThenElse(m_p, lg.hi, op::Neg(lg.hi)),
                op::IfThenElse(m_p, lg.lo, op::Neg(lg.lo))};
  const auto w = op::Mul(GammaInvExp(d, op::Sub(lnf.hi, e.hi)), c.hi);
  return {m, w};
}

// m(x) - m_t, rounded once. Both are dd, and near the root they cancel to
// nothing -- which is the entire reason they are dd.
template <class D>
HWY_INLINE op::V<D> GammaInvResid(D d, const GammaInvFwdOut<D>& f, Dd<D> mt) {
  return DdToDouble(DdAdd(d, f.m, Dd<D>{op::Neg(mt.hi), op::Neg(mt.lo)}));
}

// One candidate seed: evaluate the forward at it, score it by |m - m_t|, and
// keep it if it beats the incumbent -- carrying the forward result along, so
// the winner's evaluation becomes the first Newton step's and no candidate is
// evaluated twice.
template <class D>
HWY_INLINE void GammaInvConsider(D d, op::V<D> a, op::V<D> xc, op::V<D> i_ok,
                                 op::V<D> i_wp, Dd<D> lga, Dd<D> lna,
                                 Dd<D> mt, op::V<D>* best_x,
                                 op::V<D>* best_r, GammaInvFwdOut<D>* best_f) {
  const auto zero = op::Zero(d);
  const auto inf = op::Set(d, std::numeric_limits<double>::infinity());

  // Lt(., inf), not Lt(., DBL_MAX): a = DBL_MAX is a legal argument whose
  // answer is DBL_MAX, and rejecting that candidate strands the lane on the
  // fallback. (Found by the smoke test's beyond-resolution row.)
  const auto okv = op::Mul(i_ok, op::Mul(Ind(d, op::Gt(xc, zero)),
                                         Ind(d, op::Lt(xc, inf))));
  const auto m_ok = IndMask(d, okv);
  // A rejected candidate is replaced by an ordinary interior point before the
  // forward runs: masked-off lanes execute every op, and an x of NaN would
  // otherwise send the series to its full length and the region cores through
  // their gathers on a value nobody scrubbed.
  const auto xs = op::IfThenElse(m_ok, xc, op::Set(d, 3.0));
  const auto f = GammaInvForward(d, a, xs, i_wp, lga, lna);

  auto r = op::Abs(GammaInvResid(d, f, mt));
  r = op::IfThenElse(m_ok, r, inf);
  r = op::IfThenElse(op::IsNaN(r), inf, r);

  const auto take = op::Lt(r, *best_r);
  *best_x = op::IfThenElse(take, xs, *best_x);
  *best_r = op::IfThenElse(take, r, *best_r);
  best_f->m.hi = op::IfThenElse(take, f.m.hi, best_f->m.hi);
  best_f->m.lo = op::IfThenElse(take, f.m.lo, best_f->m.lo);
  best_f->w = op::IfThenElse(take, f.w, best_f->w);
}

// ------------------------------------------------------------------------
// The driver. kQ selects gamma_q_inv (true) or gamma_p_inv (false); the two
// differ only in how the input-side flip sets the orientation bit.
// HWY_NOINLINE for the reason in AGENTS.md: it is inlined twice per export
// (full-vector and masked-tail call sites) and is by far the heaviest driver
// in the library.
template <bool kQ, class D>
HWY_NOINLINE op::V<D> GammaInvVec(D d, op::V<D> a_in, op::V<D> s_in) {
  const auto one = op::Set(d, 1.0);
  const auto zero = op::Zero(d);
  const auto half = op::Set(d, 0.5);
  const auto inf = op::Set(d, std::numeric_limits<double>::infinity());
  const auto qnan = op::Set(d, std::numeric_limits<double>::quiet_NaN());

  // --- scrub every lane the specials table decides -------------------------
  // NaN compares false everywhere, so it falls out of i_ok with the rest.
  const auto i_ok =
      op::Mul(op::Mul(Ind(d, op::Gt(a_in, zero)), Ind(d, op::Lt(a_in, inf))),
              op::Mul(Ind(d, op::Gt(s_in, zero)), Ind(d, op::Lt(s_in, one))));
  const auto m_ok = IndMask(d, i_ok);
  const auto a = op::IfThenElse(m_ok, a_in, one);
  const auto s = op::IfThenElse(m_ok, s_in, op::Set(d, 0.25));

  // --- input-side flip: solve against t <= 1/2 -----------------------------
  const auto m_flip = op::Gt(s, half);
  const auto t = op::IfThenElse(m_flip, op::Sub(one, s), s);  // exact
  const auto i_flip = Ind(d, m_flip);
  const auto i_wq = kQ ? IndNot(d, i_flip) : i_flip;
  const auto i_wp = IndNot(d, i_wq);
  const auto m_wq = IndMask(d, i_wq);

  // --- quantities that do not depend on x ----------------------------------
  // Both sides' logs. 1 (-) t is formed EXACTLY, as a pair: when fl(1 - t)
  // rounds to 1 the low word is the entire signal, and it is what the tiny-a
  // q-side answers rest on.
  const auto lnt = GammaInvLog(d, t);  // log of the solved side; may be subnormal
  const auto lnct = GammaInvLog(d, TwoSum(d, one, op::Neg(t)));  // log(1 - t)
  const Dd<D> lnp{op::IfThenElse(m_wq, lnct.hi, lnt.hi),
                  op::IfThenElse(m_wq, lnct.lo, lnt.lo)};

  const auto lga = LgammaPosDd(d, a);
  const auto lna = GammaInvLog(d, a);
  // lgamma(1+a). For a <= 3/2 this MUST come from lgamma's own zone
  // polynomial at the exact shifted argument, exactly as GammaSmallQ does:
  // lga (+) lna is a cancellation of two ~|log a| terms whose dd residue is
  // 2^-105*|log a| ABSOLUTE, and the deep-small branch divides that by a. The
  // zone form is relative-accurate in a, which is what the tiny-a answers
  // need. Above 3/2 there is nothing to cancel and the identity is exact
  // enough (and covers every a up to the overflow boundary).
  const auto ac = op::Min(a, op::Set(d, kGammaR4AMax));
  const auto c1 = op::Ge(op::Set(d, detail::kLgammaZoneLo), a);
  const auto tz = op::IfThenElse(c1, ac, op::Sub(ac, one));  // exact both ways
  const auto lg1_zone = DdMulD(d, ZoneBracket(d, tz, c1), tz);
  const auto lg1_gen = DdAdd(d, lga, lna);
  const auto m_zone = op::Ge(op::Set(d, kGammaR4AMax), a);
  const Dd<D> lg1a{op::IfThenElse(m_zone, lg1_zone.hi, lg1_gen.hi),
                   op::IfThenElse(m_zone, lg1_zone.lo, lg1_gen.lo)};

  // --- deep-small closed form ---------------------------------------------
  // x = exp((ln p + lgamma(1+a))/a), with the power-of-two scaling applied
  // LAST so a subnormal or zero answer carries exactly one rounding.
  //
  // THE CUT IS ON x0*(1+a), NOT ON a*x0. What the closed form drops is the
  // factor e^-x * Sum(a,x), whose log is -a*h(x) with h(x) = x + O(x^2), so
  // the relative error it costs in x is h(x)/(1+a) ~ x/(1+a) -- INDEPENDENT
  // of a to leading order. The pinned constant is the right bound on that
  // quantity; a*x0 < cut is the same test only for a ~ 1, and for a << 1 it
  // is looser by 1/a. At a = 1e-4 the difference is real and reachable
  // (x0 ~ 9e-15 there satisfies a*x0 < 2^-60 while costing 2^-47 of relative
  // error), and it is reachable specifically from the q orientation, which
  // the generator's own self-check (g) sweep -- p from 1e-320 to 1e-4 -- does
  // not cover. Below x0*(1+a) < cut the dropped factor is under cut/(1+a)^2
  // for every a; above it the S2 seed's Picard iteration, which IS the exact
  // fixed point of the dropped factor, supplies the same answer to full
  // precision. See the final-report deviation note.
  const auto s0 = DdAdd(d, lnp, lg1a);
  const auto uds = GammaInvDivD(d, s0, a);
  const auto eds = ExpDdFrac(d, uds.hi, uds.lo);
  const auto x_ds = ScaleTwo(d, DdToDouble(eds.m), eds.e);
  const auto m_ds = op::Lt(op::Mul(x_ds, op::Add(one, a)),
                           op::Set(d, detail::kGammaInvDeepSmallCut));

  auto x = x_ds;
  if (!op::AllTrue(d, m_ds)) {
    // --- seed candidates -------------------------------------------------
    const auto sgn = op::IfThenElse(m_wq, one, op::Neg(one));
    op::V<D> eta0;
    const auto x1 = GammaInvSeedS1(d, a, t, sgn, &eta0);
    const auto i_ok1 =
        op::Mul(Ind(d, op::Ge(a, op::Set(d, detail::kGammaInvS1AMin))),
                Ind(d, op::Ge(op::Set(d, detail::kGammaInvEtaMax),
                              op::Abs(eta0))));
    const auto x2 = GammaInvSeedS2(d, a, DdToDouble(s0));
    op::V<D> i_ok3;
    const auto x3 = GammaInvSeedS3(d, a, lnt.hi, lga.hi, &i_ok3);
    i_ok3 = op::Mul(i_ok3, i_wq);

    // --- pick one by forward residual ------------------------------------
    // The fallback incumbent is x = a: finite, positive and inside every
    // region's domain, so a lane whose candidates all decline their guards
    // (which the generator's replay never observed) reports a defensible
    // number rather than a NaN.
    // The target's own logit: log t - log(1-t) when the solved side is P, and
    // its negation when it is Q, so that in both cases mt = log p - log q at
    // the root.
    const auto lgt = DdAdd(d, lnt, Dd<D>{op::Neg(lnct.hi), op::Neg(lnct.lo)});
    const Dd<D> mt{op::IfThenElse(m_wq, op::Neg(lgt.hi), lgt.hi),
                   op::IfThenElse(m_wq, op::Neg(lgt.lo), lgt.lo)};

    auto best_x = a;
    auto best_r = inf;
    // AGGREGATE-init, never default-construction: `GammaInvFwdOut<D> f;`
    // instantiates the struct's implicit ctor as a separate, UNATTRIBUTED
    // function, and the NEON_BF16 target then inlines Vec128's always_inline
    // ctor into it -- the exact macOS-CI break of 2026-08-06 (PLAN.md, beta
    // bf16 lesson). Brace-init has no such function; it is the fix, not the
    // hazard.
    GammaInvFwdOut<D> best_f{Dd<D>{zero, zero}, zero};
    GammaInvConsider(d, a, x1, i_ok1, i_wp, lga, lna, mt, &best_x, &best_r,
                     &best_f);
    GammaInvConsider(d, a, x2, one, i_wp, lga, lna, mt, &best_x, &best_r,
                     &best_f);
    GammaInvConsider(d, a, x3, i_ok3, i_wp, lga, lna, mt, &best_x, &best_r,
                     &best_f);

    // --- log-residual Newton, kGammaInvStepsN steps -----------------------
    // The step is Delta log x = -+ w*(ln F (-) ln t), applied MULTIPLICATIVELY
    // through expm1: x1 = x + x*expm1(step). Two reasons, and the first is
    // exactness, not robustness. For |step| below the expm1 series cut,
    // expm1(step) IS step to the last bit, so the converged step is the same
    // single rounding the additive form would have made -- nothing is given
    // up where it matters. Away from convergence the multiplicative form is
    // the correct one: ln F is close to linear in ln x in every tail (a*log x
    // in R1, -x in R2), and x*exp(step) cannot walk through zero, which
    // x*(1 - step) does the moment w*resid exceeds 1 -- a p-side deep-tail
    // seed with a residual of 40 does exactly that, and the positivity guard
    // would then strand the lane on its seed forever.
    //
    // EVERY STEP IS SAFEGUARDED: it is accepted only if it does not increase
    // |ln F (-) ln t|, and a rejected step is retried from the same point at
    // a shrunken length. This is not defensive dressing, it is what makes the
    // collapse zone right, and the ULP gate is what found it. For a above
    // ~2^105/|ln t| the whole P = 0 -> P = 1 transition happens inside one ulp
    // of x, so ln F is locally QUADRATIC in x - a (ln F ~ -(x-a)^2/(2a)) and
    // the iterate sits at its stationary point: Newton's linear model then
    // predicts a step twenty times too long, lands where the forward has
    // underflowed, and the next step -- reading a residual of 1e5 instead of
    // 1e3 -- throws the answer thousands of ulps to the far side. Measured
    // before the safeguard: 6626 ULP at a = 4.6e35, and 2e18 on the q side.
    // With it, the overshoot is simply not taken and the answer is the seed,
    // which is what PLAN's conditioning adjudication says it should be. The
    // test is inert wherever Newton behaves, i.e. everywhere the replay
    // measured, and it costs no extra work: the evaluation that scores step k
    // is the one step k+1 would have made anyway.
    //
    // Freeze by SELECT throughout, never by adding a zero step: a lane whose
    // forward is unusable (a saturated far tail, or the a > 2.5e305 band where
    // even a*log(x) overflows) keeps the value and the residual it already
    // has.
    //
    // If ANY lane failed to select a candidate, its carried-over forward is
    // the zero-initialized one, which would make the first step read a
    // residual it never computed. Recompute for the whole vector in that case
    // -- the fallback x = a is an ordinary point, so the recomputation is
    // meaningful -- and keep the reuse on the path the replay actually
    // exercises, where every lane selects something.
    auto xn_state = best_x;
    GammaInvFwdOut<D> cur = best_f;
    auto rbest = best_r;
    if (!op::AllFalse(d, op::Eq(best_r, inf))) {
      cur = GammaInvForward(d, a, xn_state, i_wp, lga, lna);
      rbest = op::Abs(GammaInvResid(d, cur, mt));
    }
    auto scale = one;
    for (int k = 0; k < detail::kGammaInvStepsN; ++k) {
      const auto resid = GammaInvResid(d, cur, mt);
      // m increases with x on BOTH orientations -- that is the point of using
      // it -- so the step needs no per-side sign, only the slope's reciprocal
      // w. The step is taken ADDITIVELY in x, x1 = x*(1 + ls), which is the
      // form the replay pinned three steps against, and it is the right one:
      // in the exponential far tail (Q ~ e^-x, so Q/g -> 1 and w -> 1/x) the
      // additive step lands on the root in ONE iteration, while the same
      // Newton written in log x converges merely quadratically -- measured
      // 49 bits after three steps against the additive form's 54 at
      // a = 1/2, q = 1/100, which is the difference between 12 ULP and
      // correctly rounded.
      //
      // The lower clamp is only a floor under a runaway step: any ls the
      // safeguard would accept is far inside it, and without it a residual
      // times w exceeding 1 would put the iterate at a negative x, which the
      // positivity guard then has to throw away entirely.
      auto ls = op::Neg(op::Mul(op::Mul(scale, cur.w), resid));
      ls = op::Max(ls, op::Set(d, -0.9));
      const auto cand = op::MulAdd(xn_state, ls, xn_state);

      const auto f = GammaInvForward(d, a, cand, i_wp, lga, lna);
      const auto rnew = op::Abs(GammaInvResid(d, f, mt));
      // The residual test only arbitrates the FAR phase. Within
      // kGammaInvTrustResid of the solution Newton is taken on trust: there
      // the two residuals differ by the forward's own noise rather than by
      // anything about the step, and letting noise veto the final refinement
      // measured 99 ULP in the small-a mid band -- a bucket that is otherwise
      // correctly rounded.
      const auto trust = Ind(d, op::Gt(op::Set(d, kGammaInvTrustResid), rbest));
      const auto acc = IndMask(
          d, op::Mul(op::Mul(Ind(d, op::Gt(cand, zero)),
                             Ind(d, op::Lt(cand, inf))),
                     op::Max(trust, Ind(d, op::Ge(rbest, rnew)))));
      xn_state = op::IfThenElse(acc, cand, xn_state);
      rbest = op::IfThenElse(acc, rnew, rbest);
      cur.m = Dd<D>{op::IfThenElse(acc, f.m.hi, cur.m.hi),
                    op::IfThenElse(acc, f.m.lo, cur.m.lo)};
      cur.w = op::IfThenElse(acc, f.w, cur.w);
      // Backtracking: a rejected step is retried an eighth as long.
      scale = op::IfThenElse(acc, one, op::Mul(scale, op::Set(d, 0.125)));
    }
    x = op::IfThenElse(m_ds, x_ds, xn_state);
  }

  // --- specials (SciPy limits; see the public header) ----------------------
  const auto at0 = kQ ? inf : zero;   // s = 0
  const auto at1 = kQ ? zero : inf;   // s = 1
  auto res = op::IfThenElse(op::Eq(s_in, zero), at0, x);
  res = op::IfThenElse(op::Eq(s_in, one), at1, res);
  res = op::IfThenElse(op::Lt(s_in, zero), qnan, res);
  res = op::IfThenElse(op::Gt(s_in, one), qnan, res);
  res = op::IfThenElse(op::Ge(zero, a_in), qnan, res);  // a <= 0, -0 included
  res = op::IfThenElse(op::Eq(a_in, inf), qnan, res);
  res = op::IfThenElse(op::IsNaN(a_in), a_in, res);  // payload preserved
  res = op::IfThenElse(op::IsNaN(s_in), s_in, res);
  return res;
}

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#endif  // include guard
