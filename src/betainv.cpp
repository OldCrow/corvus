#include "corvus/corvus.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "src/betainv.cpp"
#include "hwy/foreach_target.h"

#include "src/betainv-inl.h"
#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

// One TU for both functions, per the AGENTS.md rule that the translation unit
// boundary is the SHARING boundary: beta_p_inv and beta_q_inv run the same
// preparation, the same five seed families, the same log-space forward and the
// same Newton loop, and differ only in how the input-side flip sets the
// orientation bit. Splitting them would instantiate the whole pipeline --
// which itself instantiates all four of beta's region cores, the shared cpsi
// machinery, gammainv's three seed functions, erfcinv, digamma, trigamma and
// the dd transcendentals -- twice per target for nothing. This is the heaviest
// TU in the library and gets /d2ReducedOptimizeHugeFunctions on real MSVC
// (CMakeLists.txt) from day one.
static void BetaPInvImpl(const double* a, const double* b, const double* p,
                  double* out, size_t n) {
  const op::ScalableTag<double> d;
  const size_t N = op::Lanes(d);
  size_t i = 0;
  for (; i + N <= n; i += N) {
    op::Store(BetaInvVec<false>(d, op::Load(d, a + i), op::Load(d, b + i),
                                op::Load(d, p + i)),
              d, out + i);
  }
  if (i < n) {
    // Masked tail: same code path as the full lanes, no scalar fallback.
    const size_t m = n - i;
    op::StoreN(BetaInvVec<false>(d, op::LoadN(d, a + i, m),
                                 op::LoadN(d, b + i, m),
                                 op::LoadN(d, p + i, m)),
               d, out + i, m);
  }
}

static void BetaQInvImpl(const double* a, const double* b, const double* q,
                  double* out, size_t n) {
  const op::ScalableTag<double> d;
  const size_t N = op::Lanes(d);
  size_t i = 0;
  for (; i + N <= n; i += N) {
    op::Store(BetaInvVec<true>(d, op::Load(d, a + i), op::Load(d, b + i),
                               op::Load(d, q + i)),
              d, out + i);
  }
  if (i < n) {
    const size_t m = n - i;
    op::StoreN(BetaInvVec<true>(d, op::LoadN(d, a + i, m),
                                op::LoadN(d, b + i, m),
                                op::LoadN(d, q + i, m)),
               d, out + i, m);
  }
}

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace corvus {

HWY_EXPORT(BetaPInvImpl);
HWY_EXPORT(BetaQInvImpl);

// The dispatch calls stay INSIDE namespace corvus: with a single compiled
// target (the SSE2 cap) Highway collapses HWY_DYNAMIC_DISPATCH to
// N_SSE2::FUNC, and a globally qualified call would then name a namespace
// that does not exist. It compiles at every other tier, so only the cap
// sweep catches it.
void beta_p_inv(std::span<const double> a, std::span<const double> b,
                std::span<const double> p, std::span<double> out) {
  HWY_DYNAMIC_DISPATCH(BetaPInvImpl)
  (a.data(), b.data(), p.data(), out.data(), a.size());
}

void beta_q_inv(std::span<const double> a, std::span<const double> b,
                std::span<const double> q, std::span<double> out) {
  HWY_DYNAMIC_DISPATCH(BetaQInvImpl)
  (a.data(), b.data(), q.data(), out.data(), a.size());
}

}  // namespace corvus
#endif  // HWY_ONCE
