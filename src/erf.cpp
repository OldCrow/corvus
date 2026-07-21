#include "corvus/corvus.h"
#include "src/erf_data.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "src/erf.cpp"
#include "hwy/foreach_target.h"

#include "src/erf_core-inl.h"
#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

// See ErfTableCore (src/erf_core-inl.h) for the method. Max 1 ULP on
// validated tiers. |x| >= kErfAMax saturates to +/-1 (the bound is the
// smallest double whose erf rounds to 1, so this is correctly rounded);
// NaN is masked explicitly rather than relying on per-target Min/convert
// NaN semantics.
template <class D>
HWY_INLINE op::V<D> ErfVec(D d, op::V<D> x) {
  const auto one = op::Set(d, 1.0);
  const auto amax = op::Set(d, detail::kErfAMax);

  const auto ax = op::Abs(x);
  const auto nan = op::IsNaN(x);
  auto ac = op::Min(ax, amax);
  ac = op::IfThenElse(nan, op::Zero(d), ac);  // safe table index for NaN lanes

  const auto core = ErfTableCore(d, ac);
  auto res = op::Add(core.e_hi, core.small);

  res = op::IfThenElse(op::Ge(ax, amax), one, res);
  res = op::IfThenElse(nan, x, res);      // propagate NaN (payload preserved)
  return op::CopySign(res, x);            // erf is odd; also -0 -> -0
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

void erf(std::span<const double> in, std::span<double> out) {
  HWY_DYNAMIC_DISPATCH(ErfImpl)(in.data(), out.data(), in.size());
}

const char* active_target() {
  return HWY_DYNAMIC_DISPATCH(TargetNameImpl)();
}

}  // namespace corvus
#endif  // HWY_ONCE
