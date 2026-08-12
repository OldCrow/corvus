#include "corvus/corvus.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "src/trigamma.cpp"
#include "hwy/foreach_target.h"

#include "src/ops-inl.h"
#include "src/trigamma-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

static void TrigammaImpl(const double* in, double* out, size_t n) {
  const op::ScalableTag<double> d;
  const size_t N = op::Lanes(d);
  size_t i = 0;
  for (; i + N <= n; i += N) {
    op::Store(TrigammaVec(d, op::Load(d, in + i)), d, out + i);
  }
  if (i < n) {
    // Masked tail: same code path as the full lanes, no scalar fallback.
    op::StoreN(TrigammaVec(d, op::LoadN(d, in + i, n - i)), d, out + i, n - i);
  }
}

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace corvus {

HWY_EXPORT(TrigammaImpl);

void trigamma(std::span<const double> in, std::span<double> out) {
  HWY_DYNAMIC_DISPATCH(TrigammaImpl)(in.data(), out.data(), in.size());
}

}  // namespace corvus
#endif  // HWY_ONCE
