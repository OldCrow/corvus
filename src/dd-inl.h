// Double-double (dd) arithmetic primitives, shared by corvus's compensated
// transcendental kernels (exp_dd, log_dd) and by any kernel that needs more
// than working precision in an intermediate. Per-target include guard
// (Highway -inl.h idiom).
//
// A dd value is an unevaluated sum hi + lo of two doubles with
// |lo| <= ulp(hi)/2, so it carries ~106 significant bits. These are the
// classical algorithms (Dekker, Knuth, Kahan); they are exact-arithmetic
// identities, not approximations, and their correctness rests on:
//   * round-to-nearest, ties-to-even, with no double rounding, and
//   * no overflow/underflow of the intermediates.
// The second point is why the exp/log kernels keep their working values near
// 1 and apply the power-of-two scaling only at the very end.
//
// Nothing here uses hn:: directly: like every corvus kernel these are written
// against the ops:: facade, so the std::simd migration touches ops-inl.h only.
//
// FMA note: MulAdd is used freely below for ACCURACY (one rounding instead of
// two is a bonus, not a requirement), but never for exactness. Every exact
// residual goes through ops::ProdLow, which is capability-guarded -- Highway
// emulates MulSub as mul-then-sub on non-FMA targets, where the residual
// would silently come back zero. See the ops::SquareLow comment.
#if defined(CORVUS_DD_INL_H_) == defined(HWY_TARGET_TOGGLE)
#ifdef CORVUS_DD_INL_H_
#undef CORVUS_DD_INL_H_
#else
#define CORVUS_DD_INL_H_
#endif

#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

namespace op = ops;

// Unevaluated sum hi + lo. Lane-wise: each lane is its own dd value.
template <class D>
struct Dd {
  op::V<D> hi;
  op::V<D> lo;
};

// Exact sum of two doubles, REQUIRING |a| >= |b| (or b == 0). Cheaper than
// TwoSum by half; use it only where the ordering is guaranteed by
// construction, and say why at the call site.
template <class D>
HWY_INLINE Dd<D> Fast2Sum(D, op::V<D> a, op::V<D> b) {
  const auto s = op::Add(a, b);
  const auto z = op::Sub(s, a);  // the part of b that survived the rounding
  return {s, op::Sub(b, z)};
}

// Exact sum of two doubles, no ordering requirement (Knuth's TwoSum).
template <class D>
HWY_INLINE Dd<D> TwoSum(D, op::V<D> a, op::V<D> b) {
  const auto s = op::Add(a, b);
  const auto bv = op::Sub(s, a);
  const auto err = op::Add(op::Sub(a, op::Sub(s, bv)), op::Sub(b, bv));
  return {s, err};
}

// Exact product of two doubles.
template <class D>
HWY_INLINE Dd<D> TwoProd(D d, op::V<D> a, op::V<D> b) {
  const auto p = op::Mul(a, b);
  return {p, op::ProdLow(d, a, b, p)};
}

// dd + dd, accurate to ~2^-104 relative (Knuth's two-length addition).
template <class D>
HWY_INLINE Dd<D> DdAdd(D d, Dd<D> a, Dd<D> b) {
  const auto s = TwoSum(d, a.hi, b.hi);
  const auto t = TwoSum(d, a.lo, b.lo);
  const auto v = Fast2Sum(d, s.hi, op::Add(s.lo, t.hi));
  return Fast2Sum(d, v.hi, op::Add(v.lo, t.lo));
}

// dd + double.
template <class D>
HWY_INLINE Dd<D> DdAddD(D d, Dd<D> a, op::V<D> b) {
  const auto s = TwoSum(d, a.hi, b);
  return Fast2Sum(d, s.hi, op::Add(s.lo, a.lo));
}

// dd * dd, accurate to ~2^-104 relative. The dropped a.lo*b.lo term is
// O(2^-106) relative, below the representation's own resolution.
template <class D>
HWY_INLINE Dd<D> DdMul(D d, Dd<D> a, Dd<D> b) {
  const auto p = TwoProd(d, a.hi, b.hi);
  const auto lo = op::MulAdd(a.hi, b.lo, op::MulAdd(a.lo, b.hi, p.lo));
  // |p.hi| >= |lo|: lo is O(2^-53) relative to p.hi (both are zero together).
  return Fast2Sum(d, p.hi, lo);
}

// dd * double.
template <class D>
HWY_INLINE Dd<D> DdMulD(D d, Dd<D> a, op::V<D> b) {
  const auto p = TwoProd(d, a.hi, b);
  return Fast2Sum(d, p.hi, op::MulAdd(a.lo, b, p.lo));
}

// 1/a as a dd, from one Newton step on the rounded reciprocal.
//
// q = fl(1/a) is correct to 1/2 ulp, so p = fl(a*q) lies in [1-2^-52, 1+2^-52]
// and 1 - p is EXACT by Sterbenz. Adding the exact product residual gives the
// true 1 - a*q, hence the correction (1 - a*q)*q with relative error O(2^-105).
// Note this is precisely the shape that would be written fma(-a, q, 1) in
// scalar code -- routed through ops::ProdLow instead, because that spelling is
// silently zero on non-FMA targets.
template <class D>
HWY_INLINE Dd<D> DdRecip(D d, op::V<D> a) {
  const auto one = op::Set(d, 1.0);
  const auto q = op::Div(one, a);
  const auto p = op::Mul(a, q);
  const auto rem = op::Sub(op::Sub(one, p), op::ProdLow(d, a, q, p));
  return Fast2Sum(d, q, op::Mul(rem, q));
}

// 1/b for a dd b, as a dd. The double-argument DdRecip above with the one
// term it cannot have: q is built from b.hi alone, so the true residual is
// 1 - b*q = (1 - b.hi*q) - b.lo*q, and only the first bracket is what
// DdRecip forms. Same exactness argument as DdRecip for that bracket --
// p = fl(b.hi*q) lies in [1-2^-52, 1+2^-52] so 1 - p is EXACT by Sterbenz,
// and ops::ProdLow (never a bare MulSub: silently zero on non-FMA targets)
// supplies the product residual. The second term is already O(2^-53)
// relative, so its own rounding lands at O(2^-106).
//
// One Newton refinement then gives ~2^-105 relative, which is what the
// incomplete-gamma kernel needs from 1/(a+n), 1/sqrt(2*pi*a) and 1/Gamma(1+a):
// each of those is a factor of the RESULT, so a rounded reciprocal would put
// its own half ulp straight through.
template <class D>
HWY_INLINE Dd<D> DdRecipDd(D d, Dd<D> b) {
  const auto one = op::Set(d, 1.0);
  const auto q = op::Div(one, b.hi);
  const auto p = op::Mul(b.hi, q);
  const auto rem =
      op::Sub(op::Sub(op::Sub(one, p), op::ProdLow(d, b.hi, q, p)),
              op::Mul(b.lo, q));
  return Fast2Sum(d, q, op::Mul(rem, q));
}

// sqrt of a NON-NEGATIVE dd, ~2^-105 relative.
//
// s = fl(sqrt(a.hi)) is correct to half an ulp, so p = fl(s*s) is within a
// factor of two of a.hi and p - a.hi is EXACT by Sterbenz. Adding the exact
// product residual s^2 - p gives the true e = s^2 - a.hi with no error at
// all, and one Newton step on y^2 = a lands the correction
//     lo = (a.lo - e) / (2s),
// whose neglected term is O((a.lo - e)^2 / a^1.5) = O(2^-106) relative.
//
// The s^2 residual MUST go through ops::SquareLow, never a bare MulSub:
// Highway emulates MulSub as mul-then-sub on non-FMA targets, where the
// residual comes back identically zero and this quietly degrades to a plain
// Sqrt -- half an ulp of error where the caller asked for a hundred bits.
// (AGENTS.md, "any op whose CORRECTNESS depends on FMA fusion".)
//
// a.hi == 0 is the one input the Newton step cannot take (division by 2s);
// sqrt(0) = 0 exactly, so it is selected in. Negative a.hi is the caller's
// problem -- it produces NaN rather than trapping, which discarded lanes are
// allowed to do.
template <class D>
HWY_INLINE Dd<D> DdSqrt(D d, Dd<D> a) {
  const auto zero = op::Zero(d);
  const auto s = op::Sqrt(a.hi);
  const auto p = op::Mul(s, s);
  const auto e = op::Add(op::Sub(p, a.hi), op::SquareLow(d, s, p));
  const auto lo = op::Div(op::Sub(a.lo, e), op::Add(s, s));
  const auto r = Fast2Sum(d, s, lo);  // |s| >> |lo| by construction
  const auto z = op::Eq(a.hi, zero);
  return Dd<D>{op::IfThenElse(z, zero, r.hi), op::IfThenElse(z, zero, r.lo)};
}

// Round a dd to the nearest double (the single rounding a kernel should do
// exactly once, at the end).
template <class D>
HWY_INLINE op::V<D> DdToDouble(Dd<D> a) {
  return op::Add(a.hi, a.lo);
}

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#endif  // include guard
