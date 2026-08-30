#include "corvus/corvus.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "src/gamma.cpp"
#include "hwy/foreach_target.h"

#include "src/driver-inl.h"
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
static void GammaPImpl(std::span<const double> a, std::span<const double> x,
                       std::span<double> out) {
  DriveBinary([](auto d, auto va, auto vx) { return GammaVec<true>(d, va, vx); },
              a, x, out);
}

static void GammaQImpl(std::span<const double> a, std::span<const double> x,
                       std::span<double> out) {
  DriveBinary([](auto d, auto va, auto vx) { return GammaVec<false>(d, va, vx); },
              a, x, out);
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
  HWY_DYNAMIC_DISPATCH(GammaPImpl)(a, x, out);
}

void gamma_q(std::span<const double> a, std::span<const double> x,
             std::span<double> out) noexcept {
  HWY_DYNAMIC_DISPATCH(GammaQImpl)(a, x, out);
}

}  // namespace corvus
#endif  // HWY_ONCE
