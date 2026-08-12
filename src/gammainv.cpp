#include "corvus/corvus.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "src/gammainv.cpp"
#include "hwy/foreach_target.h"

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
static void GammaPInvImpl(const double* a, const double* p, double* out, size_t n) {
  const op::ScalableTag<double> d;
  const size_t N = op::Lanes(d);
  size_t i = 0;
  for (; i + N <= n; i += N) {
    op::Store(GammaInvVec<false>(d, op::Load(d, a + i), op::Load(d, p + i)), d,
              out + i);
  }
  if (i < n) {
    // Masked tail: same code path as the full lanes, no scalar fallback.
    const size_t m = n - i;
    op::StoreN(
        GammaInvVec<false>(d, op::LoadN(d, a + i, m), op::LoadN(d, p + i, m)),
        d, out + i, m);
  }
}

static void GammaQInvImpl(const double* a, const double* q, double* out, size_t n) {
  const op::ScalableTag<double> d;
  const size_t N = op::Lanes(d);
  size_t i = 0;
  for (; i + N <= n; i += N) {
    op::Store(GammaInvVec<true>(d, op::Load(d, a + i), op::Load(d, q + i)), d,
              out + i);
  }
  if (i < n) {
    const size_t m = n - i;
    op::StoreN(
        GammaInvVec<true>(d, op::LoadN(d, a + i, m), op::LoadN(d, q + i, m)), d,
        out + i, m);
  }
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
                 std::span<double> out) {
  HWY_DYNAMIC_DISPATCH(GammaPInvImpl)(a.data(), p.data(), out.data(), a.size());
}

void gamma_q_inv(std::span<const double> a, std::span<const double> q,
                 std::span<double> out) {
  HWY_DYNAMIC_DISPATCH(GammaQInvImpl)(a.data(), q.data(), out.data(), a.size());
}

}  // namespace corvus
#endif  // HWY_ONCE
