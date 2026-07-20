#include "corvus/corvus.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "src/erf.cpp"
#include "hwy/foreach_target.h"

#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

namespace op = ops;

// Abramowitz & Stegun 7.1.26 — PROVISIONAL bootstrap accuracy (~1.5e-7 max
// abs error). To be replaced by the max-1-ULP clean-room implementation.
template <class D>
HWY_INLINE op::V<D> ErfVec(D d, op::V<D> x) {
  const auto ax = op::Abs(x);
  const auto one = op::Set(d, 1.0);
  const auto t = op::Div(one, op::MulAdd(op::Set(d, 0.3275911), ax, one));
  auto poly = op::Set(d, 1.061405429);
  poly = op::MulAdd(poly, t, op::Set(d, -1.453152027));
  poly = op::MulAdd(poly, t, op::Set(d, 1.421413741));
  poly = op::MulAdd(poly, t, op::Set(d, -0.284496736));
  poly = op::MulAdd(poly, t, op::Set(d, 0.254829592));
  poly = op::Mul(poly, t);
  const auto e = op::Exp(d, op::Neg(op::Mul(ax, ax)));
  auto y = op::MulAdd(op::Neg(poly), e, one);
  // A&S coefficients don't sum to exactly 1, so force erf(+/-0) == +/-0.
  y = op::IfThenElse(op::Eq(ax, op::Zero(d)), op::Zero(d), y);
  return op::CopySign(y, x);
}

void ErfImpl(const double* in, double* out, size_t n) {
  const op::ScalableTag<double> d;
  const size_t N = op::Lanes(d);
  size_t i = 0;
  for (; i + N <= n; i += N) {
    op::Store(ErfVec(d, op::Load(d, in + i)), d, out + i);
  }
  if (i < n) {
    // Masked tail: same code path as the full lanes, no scalar fallback.
    op::StoreN(ErfVec(d, op::LoadN(d, in + i, n - i)), d, out + i, n - i);
  }
}

const char* TargetNameImpl() { return hwy::TargetName(HWY_TARGET); }

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace corvus {

HWY_EXPORT(ErfImpl);
HWY_EXPORT(TargetNameImpl);

void Erf(std::span<const double> in, std::span<double> out) {
  HWY_DYNAMIC_DISPATCH(ErfImpl)(in.data(), out.data(), in.size());
}

const char* ActiveTarget() {
  return HWY_DYNAMIC_DISPATCH(TargetNameImpl)();
}

}  // namespace corvus
#endif  // HWY_ONCE
