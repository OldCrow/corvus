// The arity-generic SIMD driver (#6): the full-vector loop plus masked
// LoadN/StoreN tail that every family's dispatch Impl previously
// hand-rolled (20 copies across 11 TUs). One loop shape, written once,
// owning two pieces of doctrine:
//   * Vector and tail are ONE masked code path -- no scalar libm
//     fallback for the tail. The tail call is the same kernel invocation
//     as the full-vector body, under LoadN/StoreN masking.
//   * The debug-only span-length contract check (#5 S1): corvus.h
//     declares mismatched span lengths undefined behaviour; HWY_DASSERT
//     makes any mismatch fail loudly in debug builds. Zero release
//     cost. In release the loop bound is out.size(), so an out shorter
//     than its inputs truncates in-bounds; the dangerous mistake is an
//     INPUT shorter than out, which reads past the input's end (#34
//     S2-L5 -- an earlier version of this comment had the directions
//     swapped). HWY_DASSERT comes from hwy/base.h via ops-inl.h -- not
//     an hn:: symbol, so the facade rule holds.
//
// The kernel parameter is a callable (in practice a stateless lambda
// naming one per-target Vec function, e.g. GammaVec<true>) invoked as
// kernel(d, v...). Each family's exported Impl forwards here; the Impl
// is a dispatch-table root that is never inlined, so every instantiation
// is effectively one outlined driver per family per target -- the same
// codegen structure as the hand-rolled loops, and the MSVC
// outlining rule applies per instantiation exactly as before.
//
// Nothing here uses hn:: directly (facade rule; std::simd migration
// touches ops-inl.h only).
#if defined(CORVUS_DRIVER_INL_H_) == defined(HWY_TARGET_TOGGLE)
#ifdef CORVUS_DRIVER_INL_H_
#undef CORVUS_DRIVER_INL_H_
#else
#define CORVUS_DRIVER_INL_H_
#endif

#include <span>

#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {
namespace op = ops;

template <class Kernel>
static void DriveUnary(Kernel kernel, std::span<const double> in,
                       std::span<double> out) {
  HWY_DASSERT(in.size() == out.size());
  const op::ScalableTag<double> d;
  const size_t N = op::Lanes(d);
  const size_t n = out.size();
  const double* pi = in.data();
  double* po = out.data();
  size_t i = 0;
  for (; i + N <= n; i += N) {
    op::Store(kernel(d, op::Load(d, pi + i)), d, po + i);
  }
  if (i < n) {
    const size_t m = n - i;
    op::StoreN(kernel(d, op::LoadN(d, pi + i, m)), d, po + i, m);
  }
}

template <class Kernel>
static void DriveBinary(Kernel kernel, std::span<const double> a,
                        std::span<const double> b, std::span<double> out) {
  HWY_DASSERT(a.size() == out.size() && b.size() == out.size());
  const op::ScalableTag<double> d;
  const size_t N = op::Lanes(d);
  const size_t n = out.size();
  const double* pa = a.data();
  const double* pb = b.data();
  double* po = out.data();
  size_t i = 0;
  for (; i + N <= n; i += N) {
    op::Store(kernel(d, op::Load(d, pa + i), op::Load(d, pb + i)), d, po + i);
  }
  if (i < n) {
    const size_t m = n - i;
    op::StoreN(kernel(d, op::LoadN(d, pa + i, m), op::LoadN(d, pb + i, m)), d,
               po + i, m);
  }
}

template <class Kernel>
static void DriveTernary(Kernel kernel, std::span<const double> a,
                         std::span<const double> b, std::span<const double> c,
                         std::span<double> out) {
  HWY_DASSERT(a.size() == out.size() && b.size() == out.size() &&
              c.size() == out.size());
  const op::ScalableTag<double> d;
  const size_t N = op::Lanes(d);
  const size_t n = out.size();
  const double* pa = a.data();
  const double* pb = b.data();
  const double* pc = c.data();
  double* po = out.data();
  size_t i = 0;
  for (; i + N <= n; i += N) {
    op::Store(
        kernel(d, op::Load(d, pa + i), op::Load(d, pb + i), op::Load(d, pc + i)),
        d, po + i);
  }
  if (i < n) {
    const size_t m = n - i;
    op::StoreN(kernel(d, op::LoadN(d, pa + i, m), op::LoadN(d, pb + i, m),
                      op::LoadN(d, pc + i, m)),
               d, po + i, m);
  }
}

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#endif  // include guard
