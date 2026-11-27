#include "corvus/corvus.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "src/betainv.cpp"
#include "hwy/foreach_target.h"

#include "src/betainv-inl.h"
#include "src/driver-inl.h"
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
static void BetaPInvImpl(const double* pa, const double* pb, const double* pp,
                         double* pout, size_t n) {
  const std::span<const double> a(pa, n);
  const std::span<const double> b(pb, n);
  const std::span<const double> p(pp, n);
  const std::span<double> out(pout, n);
  DriveTernary([](auto d, auto va, auto vb,
                  auto vp) { return BetaInvVec<false>(d, va, vb, vp); },
               a, b, p, out);
}

static void BetaQInvImpl(const double* pa, const double* pb, const double* pq,
                         double* pout, size_t n) {
  const std::span<const double> a(pa, n);
  const std::span<const double> b(pb, n);
  const std::span<const double> q(pq, n);
  const std::span<double> out(pout, n);
  DriveTernary([](auto d, auto va, auto vb,
                  auto vq) { return BetaInvVec<true>(d, va, vb, vq); },
               a, b, q, out);
}

// Scalar entry points (#42): the same kernel on a broadcast point.
static double BetaPInvScalarImpl(double a, double b, double p) {
  return DriveScalarTernary(
      [](auto d, auto va, auto vb, auto vp) {
        return BetaInvVec<false>(d, va, vb, vp);
      },
      a, b, p);
}

static double BetaQInvScalarImpl(double a, double b, double q) {
  return DriveScalarTernary(
      [](auto d, auto va, auto vb, auto vq) {
        return BetaInvVec<true>(d, va, vb, vq);
      },
      a, b, q);
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
                std::span<const double> p, std::span<double> out) noexcept {
  HWY_DASSERT(a.size() == out.size() && b.size() == out.size() &&
              p.size() == out.size());
  HWY_DYNAMIC_DISPATCH(BetaPInvImpl)(a.data(), b.data(), p.data(), out.data(),
                                     out.size());
}

void beta_q_inv(std::span<const double> a, std::span<const double> b,
                std::span<const double> q, std::span<double> out) noexcept {
  HWY_DASSERT(a.size() == out.size() && b.size() == out.size() &&
              q.size() == out.size());
  HWY_DYNAMIC_DISPATCH(BetaQInvImpl)(a.data(), b.data(), q.data(), out.data(),
                                     out.size());
}

HWY_EXPORT(BetaPInvScalarImpl);
HWY_EXPORT(BetaQInvScalarImpl);

double beta_p_inv(double a, double b, double p) noexcept {
  return HWY_DYNAMIC_DISPATCH(BetaPInvScalarImpl)(a, b, p);
}

double beta_q_inv(double a, double b, double q) noexcept {
  return HWY_DYNAMIC_DISPATCH(BetaQInvScalarImpl)(a, b, q);
}

}  // namespace corvus
#endif  // HWY_ONCE
