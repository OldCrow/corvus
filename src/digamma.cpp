#include "corvus/corvus.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "src/digamma.cpp"
#include "hwy/foreach_target.h"

#include "src/digamma-inl.h"
#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

static void DigammaImpl(const double* in, double* out, size_t n) {
  const op::ScalableTag<double> d;
  const size_t N = op::Lanes(d);
  size_t i = 0;
  for (; i + N <= n; i += N) {
    op::Store(DigammaVec(d, op::Load(d, in + i)), d, out + i);
  }
  if (i < n) {
    // Masked tail: same code path as the full lanes, no scalar fallback.
    op::StoreN(DigammaVec(d, op::LoadN(d, in + i, n - i)), d, out + i, n - i);
  }
}

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace corvus {

HWY_EXPORT(DigammaImpl);

// The dispatch calls stay INSIDE namespace corvus: with a single compiled
// target (the SSE2 cap) Highway collapses HWY_DYNAMIC_DISPATCH to
// N_SSE2::FUNC, and a globally qualified call would then name a namespace
// that does not exist. It compiles at every other tier, so only the cap
// sweep catches it.
void digamma(std::span<const double> in, std::span<double> out) noexcept {
  HWY_DYNAMIC_DISPATCH(DigammaImpl)(in.data(), out.data(), in.size());
}

}  // namespace corvus
#endif  // HWY_ONCE
