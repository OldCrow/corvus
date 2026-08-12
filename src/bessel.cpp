#include "corvus/corvus.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "src/bessel.cpp"
#include "hwy/foreach_target.h"

#include "src/bessel-inl.h"
#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

// ONE TU, FOUR EXPORTS, per the AGENTS.md rule that the translation unit
// boundary is the SHARING boundary: i0/i0e share BesselNu0's series and
// tail cores, i1/i1e share BesselNu1's -- splitting by export would
// instantiate each shared core twice per target for nothing. i0/i1 and
// i0e/i1e likewise share the exp_dd wrappers and the dd primitives. See
// src/bessel-inl.h's file header for the evaluation scheme and
// PLAN.md's "P2 Bessel I0/I1" binding design.

void I0Impl(const double* in, double* out, size_t n) {
  const op::ScalableTag<double> d;
  const size_t N = op::Lanes(d);
  size_t i = 0;
  for (; i + N <= n; i += N) {
    op::Store(BesselNu0(d, op::Load(d, in + i)).unscaled, d, out + i);
  }
  if (i < n) {
    // Masked tail: same code path as the full lanes, no scalar fallback.
    op::StoreN(BesselNu0(d, op::LoadN(d, in + i, n - i)).unscaled, d, out + i,
              n - i);
  }
}

void I0eImpl(const double* in, double* out, size_t n) {
  const op::ScalableTag<double> d;
  const size_t N = op::Lanes(d);
  size_t i = 0;
  for (; i + N <= n; i += N) {
    op::Store(BesselNu0(d, op::Load(d, in + i)).scaled, d, out + i);
  }
  if (i < n) {
    op::StoreN(BesselNu0(d, op::LoadN(d, in + i, n - i)).scaled, d, out + i,
              n - i);
  }
}

void I1Impl(const double* in, double* out, size_t n) {
  const op::ScalableTag<double> d;
  const size_t N = op::Lanes(d);
  size_t i = 0;
  for (; i + N <= n; i += N) {
    op::Store(BesselNu1(d, op::Load(d, in + i)).unscaled, d, out + i);
  }
  if (i < n) {
    op::StoreN(BesselNu1(d, op::LoadN(d, in + i, n - i)).unscaled, d, out + i,
              n - i);
  }
}

void I1eImpl(const double* in, double* out, size_t n) {
  const op::ScalableTag<double> d;
  const size_t N = op::Lanes(d);
  size_t i = 0;
  for (; i + N <= n; i += N) {
    op::Store(BesselNu1(d, op::Load(d, in + i)).scaled, d, out + i);
  }
  if (i < n) {
    op::StoreN(BesselNu1(d, op::LoadN(d, in + i, n - i)).scaled, d, out + i,
              n - i);
  }
}

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace corvus {

HWY_EXPORT(I0Impl);
HWY_EXPORT(I0eImpl);
HWY_EXPORT(I1Impl);
HWY_EXPORT(I1eImpl);

// The dispatch calls stay INSIDE namespace corvus: with a single compiled
// target (the SSE2 cap) Highway collapses HWY_DYNAMIC_DISPATCH to
// N_SSE2::FUNC, and a globally qualified call would then name a namespace
// that does not exist. It compiles at every other tier, so only the cap
// sweep catches it.
void i0(std::span<const double> in, std::span<double> out) {
  HWY_DYNAMIC_DISPATCH(I0Impl)(in.data(), out.data(), in.size());
}

void i0e(std::span<const double> in, std::span<double> out) {
  HWY_DYNAMIC_DISPATCH(I0eImpl)(in.data(), out.data(), in.size());
}

void i1(std::span<const double> in, std::span<double> out) {
  HWY_DYNAMIC_DISPATCH(I1Impl)(in.data(), out.data(), in.size());
}

void i1e(std::span<const double> in, std::span<double> out) {
  HWY_DYNAMIC_DISPATCH(I1eImpl)(in.data(), out.data(), in.size());
}

}  // namespace corvus
#endif  // HWY_ONCE
