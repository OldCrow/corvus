// active_target()'s own TU (#8): it shares nothing with any function
// family, so hosting it here means a consumer calling only
// active_target() links this object alone rather than a family TU and
// its coefficient tables (it previously rode in erf.cpp for no reason
// beyond history).

#include "corvus/corvus.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "src/target.cpp"
#include "hwy/foreach_target.h"

#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {
namespace op = ops;

static const char* TargetNameImpl() { return op::TargetName(); }

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace corvus {

HWY_EXPORT(TargetNameImpl);

// The dispatch call stays INSIDE namespace corvus: with a single compiled
// target (the SSE2 cap) Highway collapses HWY_DYNAMIC_DISPATCH to
// N_SSE2::FUNC, and a globally qualified call would then name a namespace
// that does not exist. It compiles at every other tier, so only the cap
// sweep catches it.
const char* active_target() noexcept {
  return HWY_DYNAMIC_DISPATCH(TargetNameImpl)();
}

}  // namespace corvus
#endif  // HWY_ONCE
