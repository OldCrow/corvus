#include "corvus/corvus.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "src/erfinv.cpp"
#include "hwy/foreach_target.h"

#include "src/erfinv-inl.h"
#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

static void ErfinvImpl(const double* in, double* out, size_t n) {
  const op::ScalableTag<double> d;
  const size_t N = op::Lanes(d);
  size_t i = 0;
  for (; i + N <= n; i += N) {
    op::Store(ErfinvVec(d, op::Load(d, in + i)), d, out + i);
  }
  if (i < n) {
    // Masked tail: same code path as the full lanes, no scalar fallback.
    op::StoreN(ErfinvVec(d, op::LoadN(d, in + i, n - i)), d, out + i, n - i);
  }
}

static void ErfcinvImpl(const double* in, double* out, size_t n) {
  const op::ScalableTag<double> d;
  const size_t N = op::Lanes(d);
  size_t i = 0;
  for (; i + N <= n; i += N) {
    op::Store(ErfcinvVec(d, op::Load(d, in + i)), d, out + i);
  }
  if (i < n) {
    op::StoreN(ErfcinvVec(d, op::LoadN(d, in + i, n - i)), d, out + i, n - i);
  }
}

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace corvus {

HWY_EXPORT(ErfinvImpl);
HWY_EXPORT(ErfcinvImpl);

void erfinv(std::span<const double> in, std::span<double> out) {
  HWY_DYNAMIC_DISPATCH(ErfinvImpl)(in.data(), out.data(), in.size());
}

void erfcinv(std::span<const double> in, std::span<double> out) {
  HWY_DYNAMIC_DISPATCH(ErfcinvImpl)(in.data(), out.data(), in.size());
}

}  // namespace corvus
#endif  // HWY_ONCE
