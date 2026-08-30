#include "corvus/corvus.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "src/bessel.cpp"
#include "hwy/foreach_target.h"

#include "src/bessel-inl.h"
#include "src/driver-inl.h"
#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

// ONE TU, FOUR EXPORTS, per the AGENTS.md rule that the translation unit
// boundary is the SHARING boundary: i0/i0e share BesselNu0's series and
// tail cores, i1/i1e share BesselNu1's -- splitting by export would
// instantiate each shared core twice per target for nothing. i0/i1 and
// i0e/i1e likewise share the exp_dd wrappers and the dd primitives.
// src/bessel-inl.h's file header and src/bessel_data.h document the
// evaluation scheme.

static void I0Impl(std::span<const double> in, std::span<double> out) {
  DriveUnary([](auto d, auto x) { return BesselNu0(d, x).unscaled; }, in, out);
}

static void I0eImpl(std::span<const double> in, std::span<double> out) {
  DriveUnary([](auto d, auto x) { return BesselNu0(d, x).scaled; }, in, out);
}

static void I1Impl(std::span<const double> in, std::span<double> out) {
  DriveUnary([](auto d, auto x) { return BesselNu1(d, x).unscaled; }, in, out);
}

static void I1eImpl(std::span<const double> in, std::span<double> out) {
  DriveUnary([](auto d, auto x) { return BesselNu1(d, x).scaled; }, in, out);
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
void i0(std::span<const double> in, std::span<double> out) noexcept {
  HWY_DYNAMIC_DISPATCH(I0Impl)(in, out);
}

void i0e(std::span<const double> in, std::span<double> out) noexcept {
  HWY_DYNAMIC_DISPATCH(I0eImpl)(in, out);
}

void i1(std::span<const double> in, std::span<double> out) noexcept {
  HWY_DYNAMIC_DISPATCH(I1Impl)(in, out);
}

void i1e(std::span<const double> in, std::span<double> out) noexcept {
  HWY_DYNAMIC_DISPATCH(I1eImpl)(in, out);
}

}  // namespace corvus
#endif  // HWY_ONCE
