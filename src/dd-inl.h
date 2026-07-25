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

// Renormalize a hi/lo pair whose lo may have grown past ulp(hi)/2.
template <class D>
HWY_INLINE Dd<D> DdNorm(D d, op::V<D> hi, op::V<D> lo) {
  return Fast2Sum(d, hi, lo);
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
