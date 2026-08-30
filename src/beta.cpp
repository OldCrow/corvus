#include "corvus/corvus.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "src/beta.cpp"
#include "hwy/foreach_target.h"

#include "src/beta-inl.h"
#include "src/driver-inl.h"
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
static void BetaPImpl(std::span<const double> a, std::span<const double> b,
                       std::span<const double> x, std::span<double> out) {
  DriveTernary(
      [](auto d, auto va, auto vb, auto vx) {
        return BetaVec<true>(d, va, vb, vx);
      },
      a, b, x, out);
}

static void BetaQImpl(std::span<const double> a, std::span<const double> b,
                       std::span<const double> x, std::span<double> out) {
  DriveTernary(
      [](auto d, auto va, auto vb, auto vx) {
        return BetaVec<false>(d, va, vb, vx);
      },
      a, b, x, out);
}

// Third export in this TU per the same sharing rule: lbeta is the PA
// prefactor's own LgammaPosDd/LgammaDiffDd assembly re-handed (see
// LbetaVec's header in beta-inl.h) -- a separate TU would re-instantiate
// that machinery per target for nothing.
static void LbetaImpl(std::span<const double> a, std::span<const double> b,
                       std::span<double> out) {
  DriveBinary([](auto d, auto va, auto vb) { return LbetaVec(d, va, vb); }, a, b,
              out);
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
            std::span<const double> x, std::span<double> out) noexcept {
  HWY_DYNAMIC_DISPATCH(BetaPImpl)(a, b, x, out);
}

void beta_q(std::span<const double> a, std::span<const double> b,
            std::span<const double> x, std::span<double> out) noexcept {
  HWY_DYNAMIC_DISPATCH(BetaQImpl)(a, b, x, out);
}

void lbeta(std::span<const double> a, std::span<const double> b,
           std::span<double> out) noexcept {
  HWY_DYNAMIC_DISPATCH(LbetaImpl)(a, b, out);
}

}  // namespace corvus
#endif  // HWY_ONCE
