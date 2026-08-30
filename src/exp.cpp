#include "corvus/corvus.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "src/exp.cpp"
#include "hwy/foreach_target.h"

#include "src/driver-inl.h"
#include "src/exp-inl.h"
#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

static void ExpImpl(std::span<const double> in, std::span<double> out) {
  DriveUnary([](auto d, auto x) { return ExpVec(d, x); }, in, out);
}

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace corvus {

HWY_EXPORT(ExpImpl);

// The dispatch call stays INSIDE namespace corvus: with a single compiled
// target (the SSE2 cap) Highway collapses HWY_DYNAMIC_DISPATCH to
// N_SSE2::FUNC, and a globally qualified call would then name a namespace
// that does not exist. It compiles at every other tier, so only the cap
// sweep catches it.
void exp(std::span<const double> in, std::span<double> out) noexcept {
  HWY_DYNAMIC_DISPATCH(ExpImpl)(in, out);
}

}  // namespace corvus
#endif  // HWY_ONCE
