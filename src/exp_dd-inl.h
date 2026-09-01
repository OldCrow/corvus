// Double-double exponential: corvus-owned, so no accuracy-critical kernel
// depends on the backend's exp. Internal (no public API): consumers are the
// erfc tail today, the incomplete gamma/beta prefactor later.
// Per-target include guard (Highway -inl.h idiom).
//
// METHOD
//   k = round(x * N/ln2),  j = k mod N,  e = (k - j)/N,  N = 128
//   r = x - k*(ln2/N),     |r| <= ln2/2N = 0.00271
//   exp(x) = 2^e * 2^(j/N) * e^r
// with 2^(j/N) tabulated as a dd pair and e^r - 1 from its Taylor series.
//
// ERROR BUDGET (relative, before the caller's final rounding). Every line is
// checked numerically by tools/gen_exp_table.py's self-check, which fails the
// generator rather than emitting a table that violates it:
//   argument reduction   2^-78.7  |k*L3|, the part of ln2/N no dd can hold
//   polynomial truncation 2^-72   r^7/5040 over |r| <= ln2/2N
//   r.lo dropped from r^2 2^-70.5 r*ulp(r)/2 through the quadratic term
//   table representation  2^-107  dd split of 2^(j/N)
//   dd assembly           2^-104  DdMul-class rounding in T*(1+p)
// Total ~2^-70, i.e. under 2^-17 ulp -- the result rounds to the correctly
// rounded double except within 2^-17 ulp of a tie.
//
// WHY THE REDUCTION IS EXACT
//   k*L1 is exact by construction: L1 keeps only 34 significant bits and
//   |k| <= 203132 needs 18, so the product fits in 53. The two subtractions
//   keep their residuals via TwoSum, and k*L2's residual via ops::ProdLow, so
//   the ONLY unrepresented part of k*(ln2/N) is k*L3 above. This matters more
//   than it looks: a naive 1/2-ulp error in r at |x| ~ 700 is a 1/2-ulp error
//   in the RESULT's exponent-scaled value, i.e. ~2^-53 relative -- the whole
//   point of the dd core is to not be that.
//
// DOMAIN
//   |x_hi| is clamped to kExpXMax = 1100; exp overflows above 709.79 and
//   flushes to zero below -745.2, so the clamp cannot change a representable
//   result. x_lo must be normalized against x_hi (|x_lo| <= ulp(x_hi)/2) --
//   callers assembling an exact product (erfc's a^2 = ssq + sl) satisfy this
//   by construction. NaN and Inf lanes are the CALLER's responsibility: the
//   table index is masked into range so nothing reads out of bounds, but the
//   value in such a lane is unspecified.
#if defined(CORVUS_EXP_DD_INL_H_) == defined(HWY_TARGET_TOGGLE)
#ifdef CORVUS_EXP_DD_INL_H_
#undef CORVUS_EXP_DD_INL_H_
#else
#define CORVUS_EXP_DD_INL_H_
#endif

#include "src/dd-inl.h"
#include "src/exp_dd_data.h"
#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

// exp(x) = (m.hi + m.lo) * 2^e, with m in ~[1, 2]. Kept unassembled so a
// consumer can fold its own factors in before the scaling rounds anything --
// that is what keeps erfc's subnormal tail at one rounding total.
template <class D>
struct ExpDdParts {
  Dd<D> m;
  op::V<op::SignedTag<D>> e;
};

// v * 2^e in two stages, so a result that lands in the subnormal range rounds
// exactly once. A single 2^e factor cannot express the range, and scaling all
// the way down in one step then multiplying would round twice.
// Precondition: |e| <= 1600 (guaranteed by ExpDdFrac's clamp on x), so both
// half-exponents stay inside the normal range and the first multiply is exact.
template <class D>
HWY_INLINE op::V<D> ScaleTwo(D d, op::V<D> v, op::V<op::SignedTag<D>> e) {
  const op::SignedTag<D> di;
  const auto e1 = op::ShiftRight<1>(e);  // floor(e/2), negatives included
  const auto e2 = op::Sub(e, e1);
  const auto bias = op::Set(di, int64_t{1023});
  const auto s1 = op::BitCast(d, op::ShiftLeft<52>(op::Add(e1, bias)));
  const auto s2 = op::BitCast(d, op::ShiftLeft<52>(op::Add(e2, bias)));
  return op::Mul(op::Mul(v, s1), s2);
}

// exp(x_hi + x_lo) in mantissa/exponent form. See the file header.
template <class D>
HWY_INLINE ExpDdParts<D> ExpDdFrac(D d, op::V<D> xh, op::V<D> xl) {
  const op::SignedTag<D> di;

  const auto lim = op::Set(d, detail::kExpXMax);
  xh = op::Min(op::Max(xh, op::Neg(lim)), lim);

  const auto kf = op::Round(op::Mul(xh, op::Set(d, detail::kExpInvL)));
  const auto ki = op::ConvertToInt(di, kf);  // exact: kf is integral, |kf| small

  // --- argument reduction, exact except for k*L3 ---
  const auto l2 = op::Set(d, detail::kExpL2);
  const auto kl1 = op::Mul(kf, op::Set(d, detail::kExpL1));  // exact
  const auto kl2 = op::Mul(kf, l2);
  const auto kl2_lo = op::ProdLow(d, kf, l2, kl2);  // exact residual
  const auto s0 = TwoSum(d, xh, op::Neg(kl1));
  const auto s1 = TwoSum(d, s0.hi, op::Neg(kl2));
  // Fast2Sum's ordering precondition USUALLY holds via |s1.hi| ~ 2^-8.5;
  // within ~2^-60 of a reduction grid point k*ln2/128 it can FAIL (s1.hi
  // shrinks below the correction operand). The result still stands, by a
  // different argument (#35 L5): when the precondition fails, both
  // operands are ~2^-60-scale, the pair error is bounded by an ulp of
  // their sum (~2^-112 ABSOLUTE in r), and that propagates through exp as
  // ~2^-112 relative -- orders below the 2^-70 budget. Verified by the
  // reference set's reduction-stress stratum and an adversarial
  // half-slot sweep at review.
  const auto r = Fast2Sum(
      d, s1.hi,
      op::Add(op::Sub(op::Add(s0.lo, s1.lo), kl2_lo), xl));

  // --- e^r - 1 = r + r^2*(1/2 + r/6 + r^2/24 + r^3/120 + r^4/720) ---
  const auto rh = r.hi;
  auto q = op::Set(d, 1.0 / 720.0);
  q = op::MulAdd(q, rh, op::Set(d, 1.0 / 120.0));
  q = op::MulAdd(q, rh, op::Set(d, 1.0 / 24.0));
  q = op::MulAdd(q, rh, op::Set(d, 1.0 / 6.0));
  q = op::MulAdd(q, rh, op::Set(d, 0.5));
  const auto p2 = op::Mul(op::Mul(rh, rh), q);
  // |rh| >= |p2| ~ rh^2/2 for |rh| < 1, and the two vanish together.
  const auto ps = Fast2Sum(d, rh, p2);
  const Dd<D> p{ps.hi, op::Add(ps.lo, r.lo)};

  // --- table lookup: j = k mod N, e = floor(k/N) ---
  // The mask also makes the gather index in-bounds for ANY input, including
  // the garbage ConvertToInt yields for a NaN lane.
  const auto j = op::And(ki, op::Set(di, detail::kExpN - 1));
  const auto e = op::ShiftRight<detail::kExpLog2N>(ki);
  const auto th = op::GatherIndex(d, detail::kExpTableHi, j);
  const auto tl = op::GatherIndex(d, detail::kExpTableLo, j);

  // --- m = T * (1 + p), T in [1, 2), |p| <= 0.0028 ---
  const auto u = TwoProd(d, th, p.hi);  // exact
  const auto a = Fast2Sum(d, th, u.hi);  // |th| >= 1 > |u.hi|
  auto lo = op::Add(op::Add(a.lo, u.lo), tl);
  lo = op::MulAdd(th, p.lo, lo);
  lo = op::MulAdd(tl, p.hi, lo);  // T.lo*p.lo ~ 2^-114 is dropped
  return {Fast2Sum(d, a.hi, lo), e};
}

// exp(x_hi + x_lo) assembled as a dd. Note the dd invariant necessarily
// degrades once the result is subnormal (lo has nowhere left to go); a
// consumer that cares about the subnormal range should fold its own factors
// into ExpDdFrac's mantissa and call ScaleTwo itself, as the erfc tail does.
template <class D>
HWY_INLINE Dd<D> ExpDd(D d, op::V<D> xh, op::V<D> xl) {
  const auto parts = ExpDdFrac(d, xh, xl);
  return {ScaleTwo(d, parts.m.hi, parts.e), ScaleTwo(d, parts.m.lo, parts.e)};
}

// ------------------------------------------------------------------------
// OUTLINED exp [MSVC BUILD-TIME GATE, AGENTS.md]. Thin wrappers whose only
// purpose is the HWY_NOINLINE: exp_dd (via ExpDd) and Highway's own Exp are
// each large, and gammainv, beta, betainv and bessel each reach ExpDd from
// several call sites (region assembly, prefactor folds, the unscaled
// Bessel exp fold), while gammainv and betainv separately reach op::Exp
// from their Newton/Picard steps and logit swap. Inlined, each call site
// becomes its own copy of a table gather plus a polynomial IN EVERY ONE OF
// THE COMPILED TARGETS, and cl.exe's optimizer is superlinear in function
// size -- at gammainv/betainv's TU scale the difference is minutes versus
// the better part of an hour on one MSVC invocation. Shared here (rather
// than once per family, as it was before A7) because the per-family
// wrappers were byte-identical modulo name; hosted beside ExpDd's other
// consumers since gammainv-inl.h and betainv-inl.h already include this
// file. Bit-identity across the outline is guaranteed by contraction-off.
template <class D>
HWY_NOINLINE Dd<D> OutlinedExpDd(D d, op::V<D> xh, op::V<D> xl) {
  return ExpDd(d, xh, xl);
}
template <class D>
HWY_NOINLINE op::V<D> OutlinedExp(D d, op::V<D> x) {
  return op::Exp(d, x);
}

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#endif  // include guard
