// SIMD op facade. Kernels use only corvus::HWY_NAMESPACE::ops — never hn::
// directly. To migrate to std::simd later, reimplement this file; kernels
// are untouched. Per-target include guard (Highway -inl.h idiom).
#if defined(CORVUS_OPS_INL_H_) == defined(HWY_TARGET_TOGGLE)
#ifdef CORVUS_OPS_INL_H_
#undef CORVUS_OPS_INL_H_
#else
#define CORVUS_OPS_INL_H_
#endif

#include "hwy/highway.h"
#include "hwy/contrib/math/math-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {
namespace ops {

namespace hn = hwy::HWY_NAMESPACE;

template <typename T> using ScalableTag = hn::ScalableTag<T>;
template <class D> using SignedTag = hn::RebindToSigned<D>;
template <class D> using V = hn::Vec<D>;
template <class D> using M = hn::Mask<D>;

using hn::Lanes;

template <class D> HWY_INLINE V<D> Load(D d, const double* p) { return hn::LoadU(d, p); }
template <class D> HWY_INLINE V<D> LoadN(D d, const double* p, size_t n) { return hn::LoadN(d, p, n); }
template <class D> HWY_INLINE void Store(V<D> v, D d, double* p) { hn::StoreU(v, d, p); }
template <class D> HWY_INLINE void StoreN(V<D> v, D d, double* p, size_t n) { hn::StoreN(v, d, p, n); }

// Lane type of a tag: Set/BitCast are used with both the double tag and the
// SignedTag (exponent assembly in exp_dd), so they are not double-only.
template <class D> using T = hn::TFromD<D>;

template <class D> HWY_INLINE V<D> Set(D d, T<D> x) { return hn::Set(d, x); }
template <class D> HWY_INLINE V<D> Zero(D d) { return hn::Zero(d); }

// Reinterpret lanes without converting: used to build 2^e from an integer
// exponent field. No-op at runtime.
template <class D, class VFrom> HWY_INLINE V<D> BitCast(D d, VFrom v) {
  return hn::BitCast(d, v);
}

template <class V> HWY_INLINE V Add(V a, V b) { return hn::Add(a, b); }
template <class V> HWY_INLINE V Sub(V a, V b) { return hn::Sub(a, b); }
template <class V> HWY_INLINE V Mul(V a, V b) { return hn::Mul(a, b); }
template <class V> HWY_INLINE V Div(V a, V b) { return hn::Div(a, b); }
template <class V> HWY_INLINE V MulAdd(V a, V b, V c) { return hn::MulAdd(a, b, c); }
// MulSub is NOT for exact residuals: on non-FMA targets Highway emulates it
// as mul-then-sub, which silently zeroes fma(a, b, -fl(a*b)). Exact residuals
// go through ProdLow/SquareLow below, which are capability-guarded. (No
// kernel currently uses MulSub; kept for the 1:1 hn:: mirror.)
template <class V> HWY_INLINE V MulSub(V a, V b, V c) { return hn::MulSub(a, b, c); }
template <class V> HWY_INLINE V Neg(V a) { return hn::Neg(a); }

template <class V> HWY_INLINE V Sqrt(V a) { return hn::Sqrt(a); }
template <class V> HWY_INLINE V Abs(V a) { return hn::Abs(a); }
template <class V> HWY_INLINE V Min(V a, V b) { return hn::Min(a, b); }
template <class V> HWY_INLINE V Max(V a, V b) { return hn::Max(a, b); }
template <class V> HWY_INLINE V CopySign(V magn, V sign) { return hn::CopySign(magn, sign); }

template <class V> HWY_INLINE auto Lt(V a, V b) { return hn::Lt(a, b); }
template <class V> HWY_INLINE auto Ge(V a, V b) { return hn::Ge(a, b); }
template <class V> HWY_INLINE auto Gt(V a, V b) { return hn::Gt(a, b); }
template <class V> HWY_INLINE auto Eq(V a, V b) { return hn::Eq(a, b); }
template <class V> HWY_INLINE auto IsNaN(V a) { return hn::IsNaN(a); }

// Round to nearest integral value, ties to even (matches FRINTN/VROUNDPD).
template <class V> HWY_INLINE V Round(V a) { return hn::Round(a); }
// Truncating float->int conversion (exact when the input is integral).
template <class DI, class V> HWY_INLINE hn::Vec<DI> ConvertToInt(DI di, V a) {
  return hn::ConvertTo(di, a);
}
template <int kBits, class V> HWY_INLINE V ShiftLeft(V a) {
  return hn::ShiftLeft<kBits>(a);
}
// Arithmetic shift on signed lanes: floor division by 2^kBits, negatives
// included (relied on to split an exponent into table index and power of 2).
template <int kBits, class V> HWY_INLINE V ShiftRight(V a) {
  return hn::ShiftRight<kBits>(a);
}
template <class V> HWY_INLINE V And(V a, V b) { return hn::And(a, b); }
template <class V> HWY_INLINE V Or(V a, V b) { return hn::Or(a, b); }
// Integer->double conversion (exact for the exponent-sized values corvus
// converts). The float->int direction is ConvertToInt above.
template <class D, class VI> HWY_INLINE V<D> ConvertToDouble(D d, VI a) {
  return hn::ConvertTo(d, a);
}
// out[i] = base[index[i]]; index in units of lanes, not bytes.
template <class D, class VI> HWY_INLINE V<D> GatherIndex(D d, const double* base, VI index) {
  return hn::GatherIndex(d, base, index);
}
template <class M, class V> HWY_INLINE V IfThenElse(M m, V yes, V no) { return hn::IfThenElse(m, yes, no); }
template <class D, class M> HWY_INLINE bool AllFalse(D d, M m) { return hn::AllFalse(d, m); }
template <class D, class M> HWY_INLINE bool AllTrue(D d, M m) { return hn::AllTrue(d, m); }

template <class D> HWY_INLINE V<D> Exp(D d, V<D> x) { return hn::Exp(d, x); }

// Low part of a*a given p = fl(a*a), i.e. the exact residual a^2 - p.
// With native FMA this is one MulSub; without it (SSE2/SSSE3/SSE4), the
// emulated mul-then-sub rounds a*a to exactly p and silently returns 0,
// so use Dekker's split (exact: a_hi has <= 26 significant bits, all
// partial products fit in a double).
template <class D> HWY_INLINE V<D> SquareLow(D d, V<D> a, V<D> p) {
#if HWY_NATIVE_FMA
  (void)d;
  return hn::MulSub(a, a, p);
#else
  const auto split = hn::Set(d, 134217729.0);  // 2^27 + 1
  const auto t = hn::Mul(a, split);
  const auto a_hi = hn::Sub(t, hn::Sub(t, a));
  const auto a_lo = hn::Sub(a, a_hi);
  const auto e = hn::Sub(hn::Mul(a_hi, a_hi), p);
  const auto cross = hn::Mul(hn::Set(d, 2.0), hn::Mul(a_hi, a_lo));
  return hn::Add(hn::Add(e, cross), hn::Mul(a_lo, a_lo));
#endif
}

// Low part of a*b given p = fl(a*b), i.e. the exact residual a*b - p; the
// two-operand form of SquareLow, and subject to the same FMA hazard, so it
// carries the same capability guard. Dekker's product needs both operands
// split. Exact for |a*b| within the normal range with room for the splits
// (|a|, |b| < 2^996); corvus's dd kernels stay far inside that.
template <class D> HWY_INLINE V<D> ProdLow(D d, V<D> a, V<D> b, V<D> p) {
#if HWY_NATIVE_FMA
  (void)d;
  return hn::MulSub(a, b, p);
#else
  const auto split = hn::Set(d, 134217729.0);  // 2^27 + 1
  const auto ta = hn::Mul(a, split);
  const auto a_hi = hn::Sub(ta, hn::Sub(ta, a));
  const auto a_lo = hn::Sub(a, a_hi);
  const auto tb = hn::Mul(b, split);
  const auto b_hi = hn::Sub(tb, hn::Sub(tb, b));
  const auto b_lo = hn::Sub(b, b_hi);
  auto e = hn::Sub(hn::Mul(a_hi, b_hi), p);
  e = hn::Add(e, hn::Mul(a_hi, b_lo));
  e = hn::Add(e, hn::Mul(a_lo, b_hi));
  return hn::Add(e, hn::Mul(a_lo, b_lo));
#endif
}

template <class D> HWY_INLINE double ReduceSum(D d, V<D> v) { return hn::ReduceSum(d, v); }
template <class D> HWY_INLINE double ReduceMax(D d, V<D> v) { return hn::ReduceMax(d, v); }

}  // namespace ops
}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#endif  // include guard
