#include "corvus/corvus.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "src/gamma.cpp"
#include "hwy/foreach_target.h"

#include "src/gamma-inl.h"
#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

// One TU for both functions, per the AGENTS.md rule that the translation
// unit boundary is the SHARING boundary: gamma_p and gamma_q route onto the
// same four region cores and differ only in the router's two decisions, so
// splitting them would instantiate every core (and the Temme table's eleven
// Clenshaw passes) twice per target for nothing.
static void GammaPImpl(const double* a, const double* x, double* out, size_t n) {
  const op::ScalableTag<double> d;
  const size_t N = op::Lanes(d);
  size_t i = 0;
  for (; i + N <= n; i += N) {
    op::Store(GammaVec<true>(d, op::Load(d, a + i), op::Load(d, x + i)), d,
              out + i);
  }
  if (i < n) {
    // Masked tail: same code path as the full lanes, no scalar fallback.
    const size_t m = n - i;
    op::StoreN(
        GammaVec<true>(d, op::LoadN(d, a + i, m), op::LoadN(d, x + i, m)), d,
        out + i, m);
  }
}

static void GammaQImpl(const double* a, const double* x, double* out, size_t n) {
  const op::ScalableTag<double> d;
  const size_t N = op::Lanes(d);
  size_t i = 0;
  for (; i + N <= n; i += N) {
    op::Store(GammaVec<false>(d, op::Load(d, a + i), op::Load(d, x + i)), d,
              out + i);
  }
  if (i < n) {
    const size_t m = n - i;
    op::StoreN(
        GammaVec<false>(d, op::LoadN(d, a + i, m), op::LoadN(d, x + i, m)), d,
        out + i, m);
  }
}

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace corvus {

HWY_EXPORT(GammaPImpl);
HWY_EXPORT(GammaQImpl);

// The dispatch calls stay INSIDE namespace corvus: with a single compiled
// target (the SSE2 cap) Highway collapses HWY_DYNAMIC_DISPATCH to
// N_SSE2::FUNC, and a globally qualified call would then name a namespace
// that does not exist. It compiles at every other tier, so only the cap
// sweep catches it.
void gamma_p(std::span<const double> a, std::span<const double> x,
             std::span<double> out) noexcept {
  HWY_DYNAMIC_DISPATCH(GammaPImpl)(a.data(), x.data(), out.data(), a.size());
}

void gamma_q(std::span<const double> a, std::span<const double> x,
             std::span<double> out) noexcept {
  HWY_DYNAMIC_DISPATCH(GammaQImpl)(a.data(), x.data(), out.data(), a.size());
}

}  // namespace corvus
#endif  // HWY_ONCE
