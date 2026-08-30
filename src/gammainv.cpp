#include "corvus/corvus.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "src/gammainv.cpp"
#include "hwy/foreach_target.h"

#include "src/driver-inl.h"
#include "src/gammainv-inl.h"
#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

// One TU for both functions, per the AGENTS.md rule that the translation unit
// boundary is the SHARING boundary: gamma_p_inv and gamma_q_inv run the same
// seed stage, the same forward evaluator and the same Newton loop, and differ
// only in how the input-side flip sets the orientation bit. Splitting them
// would instantiate the whole pipeline -- which itself instantiates all four
// forward region cores, erfcinv and the dd transcendentals -- twice per
// target for nothing.
static void GammaPInvImpl(std::span<const double> a, std::span<const double> p,
                          std::span<double> out) {
  DriveBinary([](auto d, auto va, auto vp) { return GammaInvVec<false>(d, va, vp); },
              a, p, out);
}

static void GammaQInvImpl(std::span<const double> a, std::span<const double> q,
                          std::span<double> out) {
  DriveBinary([](auto d, auto va, auto vq) { return GammaInvVec<true>(d, va, vq); },
              a, q, out);
}

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace corvus {

HWY_EXPORT(GammaPInvImpl);
HWY_EXPORT(GammaQInvImpl);

// The dispatch calls stay INSIDE namespace corvus: with a single compiled
// target (the SSE2 cap) Highway collapses HWY_DYNAMIC_DISPATCH to
// N_SSE2::FUNC, and a globally qualified call would then name a namespace
// that does not exist. It compiles at every other tier, so only the cap
// sweep catches it.
void gamma_p_inv(std::span<const double> a, std::span<const double> p,
                 std::span<double> out) noexcept {
  HWY_DYNAMIC_DISPATCH(GammaPInvImpl)(a, p, out);
}

void gamma_q_inv(std::span<const double> a, std::span<const double> q,
                 std::span<double> out) noexcept {
  HWY_DYNAMIC_DISPATCH(GammaQInvImpl)(a, q, out);
}

}  // namespace corvus
#endif  // HWY_ONCE
