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
template <class D> using V = hn::Vec<D>;
template <class D> using M = hn::Mask<D>;

using hn::Lanes;

template <class D> HWY_INLINE V<D> Load(D d, const double* p) { return hn::LoadU(d, p); }
template <class D> HWY_INLINE V<D> LoadN(D d, const double* p, size_t n) { return hn::LoadN(d, p, n); }
template <class D> HWY_INLINE void Store(V<D> v, D d, double* p) { hn::StoreU(v, d, p); }
template <class D> HWY_INLINE void StoreN(V<D> v, D d, double* p, size_t n) { hn::StoreN(v, d, p, n); }

template <class D> HWY_INLINE V<D> Set(D d, double x) { return hn::Set(d, x); }
template <class D> HWY_INLINE V<D> Zero(D d) { return hn::Zero(d); }

template <class V> HWY_INLINE V Add(V a, V b) { return hn::Add(a, b); }
template <class V> HWY_INLINE V Sub(V a, V b) { return hn::Sub(a, b); }
template <class V> HWY_INLINE V Mul(V a, V b) { return hn::Mul(a, b); }
template <class V> HWY_INLINE V Div(V a, V b) { return hn::Div(a, b); }
template <class V> HWY_INLINE V MulAdd(V a, V b, V c) { return hn::MulAdd(a, b, c); }
template <class V> HWY_INLINE V Neg(V a) { return hn::Neg(a); }

template <class V> HWY_INLINE V Sqrt(V a) { return hn::Sqrt(a); }
template <class V> HWY_INLINE V Abs(V a) { return hn::Abs(a); }
template <class V> HWY_INLINE V Min(V a, V b) { return hn::Min(a, b); }
template <class V> HWY_INLINE V Max(V a, V b) { return hn::Max(a, b); }
template <class V> HWY_INLINE V CopySign(V magn, V sign) { return hn::CopySign(magn, sign); }

template <class V> HWY_INLINE auto Lt(V a, V b) { return hn::Lt(a, b); }
template <class V> HWY_INLINE auto Eq(V a, V b) { return hn::Eq(a, b); }
template <class M, class V> HWY_INLINE V IfThenElse(M m, V yes, V no) { return hn::IfThenElse(m, yes, no); }

template <class D> HWY_INLINE V<D> Exp(D d, V<D> x) { return hn::Exp(d, x); }

template <class D> HWY_INLINE double ReduceSum(D d, V<D> v) { return hn::ReduceSum(d, v); }
template <class D> HWY_INLINE double ReduceMax(D d, V<D> v) { return hn::ReduceMax(d, v); }

}  // namespace ops
}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#endif  // include guard
