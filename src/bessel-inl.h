// Modified Bessel functions of the first kind, order 0 and 1, and their
// exponentially-scaled variants: i0(x), i1(x), i0e(x) = e^-|x|*I0(x),
// i1e(x) = sign(x)*e^-|x|*I1(|x|). Per-target include guard (Highway -inl.h
// idiom).
//
// Per PLAN.md's "P2 Bessel I0/I1 -- BINDING DESIGN" and the parameters
// pinned in src/bessel_data.h (whose header comment specifies the exact
// evaluation scheme this file follows): two regimes, split at
// kBesselSplit = 8, on ax = |x| -- both I0 and I1 are functions of ax alone
// (i0/i0e even, i1/i1e odd), so ONE table serves both signs and the odd
// forms reapply sign with CopySign at the very end, never implicitly.
//
// SERIES REGIME, ax <= kBesselSplit. I0(x) = S0(q), I1(x) = (ax/2)*S1(q),
// q = x^2/4, both series all-positive and perfectly conditioned. The scaled
// forms multiply by exp_dd(-ax) AFTER the unscaled series is known -- the
// unscaled value is the series itself, not a re-derived quantity, so no
// exponential is needed at all to answer i0/i1 in this regime.
//
// q IS CAPTURED EXACTLY. A rounded q = fl(x^2/4) costs ~x ULP through the
// series' own log-sensitivity q*S'/S = (x/2)*I1/I0 (grows with x -- this is
// the FIRST correction ratified at G1, not a probe artifact), so q is split
// into an exact hi/lo pair via ops::SquareLow's residual (the erfc/erfinv
// ssq/sl idiom) and the series is evaluated as
//     S(q_hi + q_lo) ~= S(q_hi) + S'(q_hi)*q_lo,
// a first-order derivative correction in PLAIN DOUBLE (the correction
// itself is already ~2^-53 relative of the result, so it needs no more).
// S(q_hi) itself keeps the lowest kBesselI{0,1}SeriesLead coefficients as
// dd (erfinv-inl.h's ErfinvCentralDd nested dd/double pattern); S'(q_hi)
// is one plain-double Horner over the *DCoef arrays, which are literally
// the termwise derivative of S (dcoef[k] = (k+1)*coef[k+1]).
//
// TAIL REGIME, ax > kBesselSplit. i_nu_e(x) = f_nu(t)/sqrt(2*pi*x), t = 1/x,
// f_nu a Chebyshev refit in normalized s = t*Scale + Shift (own nodes and
// budget, clean-room -- A&S coefficients were NOT ported, see PLAN.md).
// Plain-double Horner: t is flat/well-conditioned there and klead 0/1/2
// were measured to land at the same floor (erfc-tail precedent). The
// divide-by-sqrt(2*pi*x) is dd-assisted (DdSqrt then DdRecipDd), so the
// scaled result rounds exactly once; 2*pi itself is carried as a dd pair
// (kBesselTwoPiHi/Lo below) so its own double-rounding cannot spend any of
// that budget -- this is a kernel-internal MATH constant, not fitted table
// data, so it lives here rather than in the generated src/bessel_data.h.
// Unscaled forms in this regime multiply the scaled value by exp_dd(+ax)
// via ExpDdFrac's mantissa+exponent split, SCALING LAST so a result near
// the overflow boundary rounds once -- then saturate to +inf (i1: signed)
// past the EXACT pinned boundary (src/bessel_data.h), never relying on
// incidental IEEE overflow from the assembly, which could disagree with
// the boundary by a diagnosable-but-avoidable ULP.
//
// SPECIALS. x = 0: i0 = i0e = 1 (falls out of the series at q = 0, no
// override needed: S(0) = lead coefficient = 1 exactly for both nu); i1 =
// i1e = +-0 (series magnitude is exactly 0 at ax = 0, and CopySign carries
// the sign of x, including -0, onto it). x = +-inf: i0 = +inf, i0e = +0,
// i1 = +-inf, i1e = +-0 -- all obtained by scrubbing ax to a benign finite
// point for the arithmetic, then overriding the MAGNITUDE with the correct
// special before CopySign runs (i1/i1e) or directly (i0/i0e, which carry
// no sign). NaN propagates with its payload, applied last (after
// CopySign), matching every other corvus kernel's masked-lane-scrub
// convention (AGENTS.md).
#if defined(CORVUS_BESSEL_INL_H_) == defined(HWY_TARGET_TOGGLE)
#ifdef CORVUS_BESSEL_INL_H_
#undef CORVUS_BESSEL_INL_H_
#else
#define CORVUS_BESSEL_INL_H_
#endif

#include <limits>

#include "src/bessel_data.h"
#include "src/dd-inl.h"
#include "src/exp_dd-inl.h"
#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

namespace op = ops;

// 2*pi as a dd pair, ~2^-104 relative -- independently derived (not part of
// the generated fit data, which carries only the fitted series/tail
// parameters). Used solely by the tail's dd-assisted sqrt(2*pi*x); keeping
// it at dd precision rather than a single rounded double removes 2*pi's own
// ~2^-54 rounding error from the tail's budget entirely, for the cost of
// one extra DdMulD.
inline constexpr double kBesselTwoPiHi = 0x1.921fb54442d18p+2;
inline constexpr double kBesselTwoPiLo = 0x1.1a62633145c07p-52;

// EXACT power-of-two prescale/postscale for the tail's sqrt(2*pi*ax), per
// the AGENTS.md/lgamma-Stirling precedent ("a kernel whose operands can
// reach [a hazardous] range must scale by a power of two first and scale
// back after"). ax reaches DBL_MAX (the tail's own domain, per the
// reference set's DBL_MAX-neighborhood coverage), where 2*pi*ax OVERFLOWS
// (2*pi*DBL_MAX ~= 1.13e309). The BINDING constraint is NOT DBL_MAX but
// ops-inl.h ProdLow's Dekker-split operand bound: on the tiers WITHOUT
// native FMA (SSE4/SSSE3/SSE2) every dd TwoProd splits its operands via
// multiplication by 2^27+1, exact only for operands below ~2^996
// (~6.7e299). The original 2^-8 prescale kept the product under DBL_MAX --
// sufficient for the FMA tiers, where this was validated -- but left
// DdMulD's ax operand at up to ~7e305, far past the split bound, producing
// ~2^63-ULP garbage in i0e/i1e at the SSE4 tier for x ~> 4e302 (G4 sweep
// catch, 2026-08-11). 2^-32 (even exponent, so the sqrt's scale stays an
// exact power of two) puts the largest split operand at ~4.2e299 < 2^996
// with 1.6x headroom; the sqrt of a 2^-32 scale is 2^-16, so the postscale
// on the RESULT is 2^+16. Small side is safe: live tail lanes have
// ax >= 8, so the prescaled operand never approaches the subnormal range.
inline constexpr double kBesselSqrtPrescale = 0x1.0p-32;
inline constexpr double kBesselSqrtPostscale = 0x1.0p+16;

// Element counts of the *DCoef arrays, derived from their declared extent
// rather than duplicated as a literal: src/bessel_data.h fixes the array
// size but does not separately name it (it is the termwise derivative of
// *SeriesCoef plus the lead terms, one shorter than *SeriesCoef's own count
// plus lead -- deriving it here keeps this file from silently drifting out
// of sync with a future regeneration of the data file).
inline constexpr int kBesselI0SeriesNDCoef =
    static_cast<int>(sizeof(detail::kBesselI0SeriesDCoef) /
                     sizeof(detail::kBesselI0SeriesDCoef[0]));
inline constexpr int kBesselI1SeriesNDCoef =
    static_cast<int>(sizeof(detail::kBesselI1SeriesDCoef) /
                     sizeof(detail::kBesselI1SeriesDCoef[0]));

// --- outlined exp wrappers [MSVC BUILD-TIME GATE, AGENTS.md] --------------
// Thin wrappers whose only purpose is the HWY_NOINLINE: exp_dd's table
// gather and polynomial are reached from four call sites in this file (one
// per unscaled assembly, one per scaled-region exp fold), and cl.exe's
// optimizer is superlinear in function size -- outlining every heavy
// callee from day one is the betainv/gammainv/trigamma precedent (AGENTS.md,
// PLAN.md). Bit-identity is guaranteed by contraction-off.
template <class D>
HWY_NOINLINE Dd<D> BesselExpDd(D d, op::V<D> xh, op::V<D> xl) {
  return ExpDd(d, xh, xl);
}
template <class D>
HWY_NOINLINE ExpDdParts<D> BesselExpDdFrac(D d, op::V<D> xh, op::V<D> xl) {
  return ExpDdFrac(d, xh, xl);
}

// --- small shared polynomial evaluators (HWY_INLINE: erfinv's central-poly
// precedent -- genuinely small hot helpers stay inline, AGENTS.md) --------
// Own local copies rather than sharing digamma/trigamma's near-identical
// helpers: hoisting a helper out of a SHIPPED family invokes the
// byte-identity protocol on that family's whole gate set for no benefit to
// a family with no shared math (trigamma-inl.h's own stated rule -- "if a
// third consumer appears, hoist all of it at once").

// P(t) = lead[0] + lead[1]*t + ... + lead[n_lead-1]*t^(n_lead-1)
//        + t^n_lead * (coef[0] + coef[1]*t + ...), t a PLAIN DOUBLE
// argument, lead coefficients dd, tail coefficients plain double.
template <class D>
HWY_INLINE Dd<D> BesselLeadTailScalar(D d, op::V<D> t, const double* lead_hi,
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

// Plain-double Horner, ascending-degree coefficients, PLAIN DOUBLE argument
// and result: used both for the series' derivative correction (dcoef[k] =
// (k+1)*coef[k+1], the termwise derivative of the series in q) and for the
// tail's Chebyshev fit itself.
template <class D>
HWY_INLINE op::V<D> BesselPoly(D d, op::V<D> t, const double* coef, int n) {
  auto s = op::Set(d, coef[n - 1]);
  for (int k = n - 2; k >= 0; --k) {
    s = op::MulAdd(s, t, op::Set(d, coef[k]));
  }
  return s;
}

// --- series regime core ----------------------------------------------------
// S(q_hi + q_lo) via dd-lead/tail Horner at q_hi plus the first-order
// derivative correction S'(q_hi)*q_lo, both in the file header. HWY_NOINLINE
// like every other region core (AGENTS.md).
template <class D>
HWY_NOINLINE Dd<D> BesselSeriesS(D d, op::V<D> q_hi, op::V<D> q_lo,
                                 const double* lead_hi, const double* lead_lo,
                                 int n_lead, const double* coef, int n_coef,
                                 const double* dcoef, int n_dcoef) {
  const auto s = BesselLeadTailScalar(d, q_hi, lead_hi, lead_lo, n_lead, coef,
                                      n_coef);
  const auto sp = BesselPoly(d, q_hi, dcoef, n_dcoef);
  return DdAddD(d, s, op::Mul(sp, q_lo));
}

// --- tail regime core --------------------------------------------------
// i_nu_e(ax) = f_nu(t)/sqrt(2*pi*ax), as a dd (unrounded -- the caller
// decides whether this is the final scaled answer or an input to the
// unscaled exp_dd assembly, and either way there should be exactly one
// rounding). HWY_NOINLINE like every other region core.
template <class D>
HWY_NOINLINE Dd<D> BesselTailIveDd(D d, op::V<D> ax, const double* coef,
                                   int n_coef, double scale, double shift) {
  const auto one = op::Set(d, 1.0);
  const auto t = op::Div(one, ax);  // plain double: flat/well-conditioned
  const auto s = op::MulAdd(t, op::Set(d, scale), op::Set(d, shift));
  const auto poly = BesselPoly(d, s, coef, n_coef);

  const Dd<D> twopi{op::Set(d, kBesselTwoPiHi), op::Set(d, kBesselTwoPiLo)};
  const auto ax_ds = op::Mul(ax, op::Set(d, kBesselSqrtPrescale));  // exact
  const auto twopix = DdMulD(d, twopi, ax_ds);  // dd, safely bounded
  const auto denom_ds = DdSqrt(d, twopix);      // dd sqrt(2*pi*ax) * 2^-4
  const auto denom =
      DdMulD(d, denom_ds, op::Set(d, kBesselSqrtPostscale));  // exact rescale
  const auto recip = DdRecipDd(d, denom);    // dd 1/sqrt(2*pi*ax)
  return DdMulD(d, recip, poly);
}

// --- per-nu drivers ----------------------------------------------------
// One driver per order (nu = 0, 1) returns {unscaled, scaled} together: the
// two share every region core and the split/scrub logic, so exporting them
// as a pair keeps the per-target instantiation count at two (not four)
// while costing only one extra exp_dd fold and rounding for whichever half
// a given public export does not need. HWY_NOINLINE: the driver is inlined
// TWICE per export (full-vector and masked-tail call sites), so leaving it
// inline doubles the largest function the optimizer sees (AGENTS.md).
template <class D>
struct BesselPair {
  op::V<D> unscaled;
  op::V<D> scaled;
};

// i0/i0e: EVEN, no sign handling needed anywhere in this driver.
template <class D>
HWY_NOINLINE BesselPair<D> BesselNu0(D d, op::V<D> x) {
  const auto zero = op::Zero(d);
  const auto one = op::Set(d, 1.0);
  const auto quarter = op::Set(d, 0.25);
  const auto inf = op::Set(d, std::numeric_limits<double>::infinity());
  const auto maxf = op::Set(d, (std::numeric_limits<double>::max)());

  const auto ax_raw = op::Abs(x);
  const auto is_nan = op::IsNaN(x);
  const auto is_inf = op::Gt(ax_raw, maxf);

  // SCRUB. Masked-off lanes still execute every op (AGENTS.md): a NaN or
  // infinite ax must not reach the dd/exp arithmetic, so both are replaced
  // by a benign interior point before anything else runs; the true answer
  // is restored by the specials overrides below.
  auto ax = op::IfThenElse(is_inf, one, ax_raw);
  ax = op::IfThenElse(is_nan, one, ax);

  const auto split = op::Set(d, detail::kBesselSplit);
  const auto in_series = op::Ge(split, ax);  // ax <= split

  // --- series branch (computed on every lane; selected in on in_series) ---
  const auto ax_s = op::Min(ax, split);
  const auto ssq = op::Mul(ax_s, ax_s);
  const auto sl = op::SquareLow(d, ax_s, ssq);
  const auto q_hi = op::Mul(ssq, quarter);  // exact: *0.25 is a power of 2
  const auto q_lo = op::Mul(sl, quarter);

  const auto s0 = BesselSeriesS(
      d, q_hi, q_lo, detail::kBesselI0SeriesLeadHi,
      detail::kBesselI0SeriesLeadLo, detail::kBesselI0SeriesLead,
      detail::kBesselI0SeriesCoef, detail::kBesselI0SeriesNCoef,
      detail::kBesselI0SeriesDCoef, kBesselI0SeriesNDCoef);
  const auto exp_neg = BesselExpDd(d, op::Neg(ax_s), zero);
  const auto i0e_series = DdMul(d, s0, exp_neg);
  const auto i0_series_val = op::Add(s0.hi, s0.lo);

  // --- tail branch (computed on every lane; selected in on !in_series) ----
  const auto ax_t = op::Max(ax, split);
  const auto ive_tail = BesselTailIveDd(
      d, ax_t, detail::kBesselI0TailCoef, detail::kBesselI0TailNCoef,
      detail::kBesselI0TailScale, detail::kBesselI0TailShift);

  const auto parts = BesselExpDdFrac(d, ax_t, zero);
  const auto prod = DdMul(d, ive_tail, parts.m);
  const auto sc_hi = ScaleTwo(d, prod.hi, parts.e);
  const auto sc_lo = ScaleTwo(d, prod.lo, parts.e);
  auto i0_tail_val = op::Add(sc_hi, sc_lo);
  const auto boundary = op::Set(d, detail::kBesselI0OverflowX);
  i0_tail_val = op::IfThenElse(op::Gt(ax, boundary), inf, i0_tail_val);

  // --- combine, then specials, NaN last (payload preserved) --------------
  auto i0_val = op::IfThenElse(in_series, i0_series_val, i0_tail_val);
  auto i0e_val = op::IfThenElse(
      in_series, op::Add(i0e_series.hi, i0e_series.lo),
      op::Add(ive_tail.hi, ive_tail.lo));

  i0_val = op::IfThenElse(is_inf, inf, i0_val);
  i0e_val = op::IfThenElse(is_inf, zero, i0e_val);
  i0_val = op::IfThenElse(is_nan, x, i0_val);
  i0e_val = op::IfThenElse(is_nan, x, i0e_val);

  return {i0_val, i0e_val};
}

// i1/i1e: ODD. Every step above computes a MAGNITUDE on ax; CopySign(mag, x)
// is applied exactly once, at the very end, covering the regular case,
// +-0 and +-inf uniformly -- never an implicit sign carried through a dd
// assembly's internal Fast2Sum, which can turn -0 into +0 under
// round-to-nearest (the exact erfinv(-0) hazard src/erfinv-inl.h's
// ErfinvCentral documents and guards the same way; see src/bessel_data.h's
// header note).
template <class D>
HWY_NOINLINE BesselPair<D> BesselNu1(D d, op::V<D> x) {
  const auto zero = op::Zero(d);
  const auto one = op::Set(d, 1.0);
  const auto half = op::Set(d, 0.5);
  const auto quarter = op::Set(d, 0.25);
  const auto inf = op::Set(d, std::numeric_limits<double>::infinity());
  const auto maxf = op::Set(d, (std::numeric_limits<double>::max)());

  const auto ax_raw = op::Abs(x);
  const auto is_nan = op::IsNaN(x);
  const auto is_inf = op::Gt(ax_raw, maxf);

  auto ax = op::IfThenElse(is_inf, one, ax_raw);
  ax = op::IfThenElse(is_nan, one, ax);

  const auto split = op::Set(d, detail::kBesselSplit);
  const auto in_series = op::Ge(split, ax);

  // --- series branch -------------------------------------------------------
  const auto ax_s = op::Min(ax, split);
  const auto ssq = op::Mul(ax_s, ax_s);
  const auto sl = op::SquareLow(d, ax_s, ssq);
  const auto q_hi = op::Mul(ssq, quarter);
  const auto q_lo = op::Mul(sl, quarter);

  const auto s1 = BesselSeriesS(
      d, q_hi, q_lo, detail::kBesselI1SeriesLeadHi,
      detail::kBesselI1SeriesLeadLo, detail::kBesselI1SeriesLead,
      detail::kBesselI1SeriesCoef, detail::kBesselI1SeriesNCoef,
      detail::kBesselI1SeriesDCoef, kBesselI1SeriesNDCoef);
  const auto half_ax = op::Mul(ax_s, half);  // exact: *0.5 is a power of 2
  const auto mag_series = DdMulD(d, s1, half_ax);  // unscaled |I1|, dd
  const auto exp_neg = BesselExpDd(d, op::Neg(ax_s), zero);
  const auto i1e_series = DdMul(d, mag_series, exp_neg);
  const auto i1_series_mag = op::Add(mag_series.hi, mag_series.lo);

  // --- tail branch -----------------------------------------------------
  const auto ax_t = op::Max(ax, split);
  const auto ive_tail = BesselTailIveDd(
      d, ax_t, detail::kBesselI1TailCoef, detail::kBesselI1TailNCoef,
      detail::kBesselI1TailScale, detail::kBesselI1TailShift);

  const auto parts = BesselExpDdFrac(d, ax_t, zero);
  const auto prod = DdMul(d, ive_tail, parts.m);
  const auto sc_hi = ScaleTwo(d, prod.hi, parts.e);
  const auto sc_lo = ScaleTwo(d, prod.lo, parts.e);
  auto i1_tail_mag = op::Add(sc_hi, sc_lo);
  const auto boundary = op::Set(d, detail::kBesselI1OverflowX);
  i1_tail_mag = op::IfThenElse(op::Gt(ax, boundary), inf, i1_tail_mag);

  // --- combine magnitudes, specials on the magnitude, THEN sign ----------
  auto i1_mag = op::IfThenElse(in_series, i1_series_mag, i1_tail_mag);
  auto i1e_mag = op::IfThenElse(
      in_series, op::Add(i1e_series.hi, i1e_series.lo),
      op::Add(ive_tail.hi, ive_tail.lo));

  i1_mag = op::IfThenElse(is_inf, inf, i1_mag);
  i1e_mag = op::IfThenElse(is_inf, zero, i1e_mag);

  auto i1_val = op::CopySign(i1_mag, x);
  auto i1e_val = op::CopySign(i1e_mag, x);

  // NaN LAST, replacing the CopySign'd result entirely so the payload bits
  // survive exactly (CopySign only ever touches the sign bit of a NaN's
  // magnitude operand, which is not the same as preserving x's own bits).
  i1_val = op::IfThenElse(is_nan, x, i1_val);
  i1e_val = op::IfThenElse(is_nan, x, i1e_val);

  return {i1_val, i1e_val};
}

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#endif  // include guard
