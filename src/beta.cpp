#include "corvus/corvus.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "src/beta.cpp"
#include "hwy/foreach_target.h"

#include "src/beta-inl.h"
#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

// One TU for both functions, per the AGENTS.md rule that the translation unit
// boundary is the SHARING boundary: beta_p and beta_q route onto the same four
// region cores, the same two prefactor paths and the same router, and differ
// only in the final handout -- splitting them would instantiate every core
// (and the R3 table's ten nested Clenshaw passes) twice per target for
// nothing. The TU does NOT include src/gamma-inl.h: R2's sweep covers the
// gamma-limit corner, so there is no gamma-core reuse to share.
void BetaPImpl(const double* a, const double* b, const double* x, double* out,
               size_t n) {
  const op::ScalableTag<double> d;
  const size_t N = op::Lanes(d);
  size_t i = 0;
  for (; i + N <= n; i += N) {
    op::Store(BetaVec<true>(d, op::Load(d, a + i), op::Load(d, b + i),
                            op::Load(d, x + i)),
              d, out + i);
  }
  if (i < n) {
    // Masked tail: same code path as the full lanes, no scalar fallback.
    const size_t m = n - i;
    op::StoreN(BetaVec<true>(d, op::LoadN(d, a + i, m), op::LoadN(d, b + i, m),
                             op::LoadN(d, x + i, m)),
               d, out + i, m);
  }
}

void BetaQImpl(const double* a, const double* b, const double* x, double* out,
               size_t n) {
  const op::ScalableTag<double> d;
  const size_t N = op::Lanes(d);
  size_t i = 0;
  for (; i + N <= n; i += N) {
    op::Store(BetaVec<false>(d, op::Load(d, a + i), op::Load(d, b + i),
                             op::Load(d, x + i)),
              d, out + i);
  }
  if (i < n) {
    const size_t m = n - i;
    op::StoreN(BetaVec<false>(d, op::LoadN(d, a + i, m), op::LoadN(d, b + i, m),
                              op::LoadN(d, x + i, m)),
               d, out + i, m);
  }
}

// Third export in this TU per the same sharing rule: lbeta is the PA
// prefactor's own LgammaPosDd/LgammaDiffDd assembly re-handed (see
// LbetaVec's header in beta-inl.h) -- a separate TU would re-instantiate
// that machinery per target for nothing.
void LbetaImpl(const double* a, const double* b, double* out, size_t n) {
  const op::ScalableTag<double> d;
  const size_t N = op::Lanes(d);
  size_t i = 0;
  for (; i + N <= n; i += N) {
    op::Store(LbetaVec(d, op::Load(d, a + i), op::Load(d, b + i)), d, out + i);
  }
  if (i < n) {
    const size_t m = n - i;
    op::StoreN(LbetaVec(d, op::LoadN(d, a + i, m), op::LoadN(d, b + i, m)), d,
               out + i, m);
  }
}

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace corvus {

HWY_EXPORT(BetaPImpl);
HWY_EXPORT(BetaQImpl);
HWY_EXPORT(LbetaImpl);

// The dispatch calls stay INSIDE namespace corvus: with a single compiled
// target (the SSE2 cap) Highway collapses HWY_DYNAMIC_DISPATCH to
// N_SSE2::FUNC, and a globally qualified call would then name a namespace that
// does not exist. It compiles at every other tier, so only the cap sweep
// catches it.
void beta_p(std::span<const double> a, std::span<const double> b,
            std::span<const double> x, std::span<double> out) {
  HWY_DYNAMIC_DISPATCH(BetaPImpl)
  (a.data(), b.data(), x.data(), out.data(), a.size());
}

void beta_q(std::span<const double> a, std::span<const double> b,
            std::span<const double> x, std::span<double> out) {
  HWY_DYNAMIC_DISPATCH(BetaQImpl)
  (a.data(), b.data(), x.data(), out.data(), a.size());
}

void lbeta(std::span<const double> a, std::span<const double> b,
           std::span<double> out) {
  HWY_DYNAMIC_DISPATCH(LbetaImpl)
  (a.data(), b.data(), out.data(), a.size());
}

}  // namespace corvus
#endif  // HWY_ONCE
