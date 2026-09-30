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
// The masked tail pads its dead lanes with the tail's FIRST element rather
// than zero (#42). A zero lane is not a no-op for the region kernels: their
// specials scrub rewrites it to an interior safe point, whose region core
// then runs for the padding alone -- a span-of-1 gamma_p call cost 2.2-2.5x
// a full vector on AVX2 for that reason. A live element's copy shares its
// region and its convergence, so the tail costs what its live lanes cost.
// Results are unchanged: every lane is computed by its own region core and
// frozen on its own mask, so a neighbour's value never reaches it (the
// lane-mix determinism checks in the smoke tests are the guard).
//
// DriveScalar* are the scalar entry points (#42): one point broadcast to
// every lane, so a single call costs one vector in one region and nothing
// more, with no loop and no masked load/store. Bit-identical to the span
// form for the same reason the padding is.
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
    op::StoreN(kernel(d, op::LoadNOr(op::Set(d, pi[i]), d, pi + i, m)), d, po + i, m);
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
    op::StoreN(kernel(d, op::LoadNOr(op::Set(d, pa[i]), d, pa + i, m),
                      op::LoadNOr(op::Set(d, pb[i]), d, pb + i, m)),
               d, po + i, m);
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
    op::StoreN(kernel(d, op::LoadNOr(op::Set(d, pa[i]), d, pa + i, m),
                      op::LoadNOr(op::Set(d, pb[i]), d, pb + i, m),
                      op::LoadNOr(op::Set(d, pc[i]), d, pc + i, m)),
               d, po + i, m);
  }
}

template <class Kernel>
static double DriveScalarBinary(Kernel kernel, double a, double b) {
  const op::ScalableTag<double> d;
  return op::GetLane(kernel(d, op::Set(d, a), op::Set(d, b)));
}

template <class Kernel>
static double DriveScalarTernary(Kernel kernel, double a, double b, double c) {
  const op::ScalableTag<double> d;
  return op::GetLane(kernel(d, op::Set(d, a), op::Set(d, b), op::Set(d, c)));
}

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#endif  // include guard
