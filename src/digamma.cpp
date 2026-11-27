#include "corvus/corvus.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "src/digamma.cpp"
#include "hwy/foreach_target.h"

#include "src/digamma-inl.h"
#include "src/driver-inl.h"
#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

static void DigammaImpl(const double* pin, double* pout, size_t n) {
  const std::span<const double> in(pin, n);
  const std::span<double> out(pout, n);
  DriveUnary([](auto d, auto x) { return DigammaVec(d, x); }, in, out);
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
  HWY_DASSERT(in.size() == out.size());
  HWY_DYNAMIC_DISPATCH(DigammaImpl)(in.data(), out.data(), out.size());
}

}  // namespace corvus
#endif  // HWY_ONCE
