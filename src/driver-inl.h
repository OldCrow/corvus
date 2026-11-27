// The arity-generic SIMD driver (#6): the full-vector loop plus the tail
// that every family's dispatch Impl previously hand-rolled (20 copies
// across 11 TUs). One loop shape, written once, owning two pieces of
// doctrine:
//   * Vector and tail are ONE code path -- no scalar libm fallback for the
//     tail. The tail call is the same kernel invocation as the full-vector
//     body, on a vector whose spare lanes hold a live element
//     (LoadTail/StoreTail below).
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
// The tail pads its dead lanes with a live element (its last) rather than
// zero (#42). A zero lane is not a no-op for the region kernels: their
// specials scrub rewrites it to an interior safe point, whose region core
// then runs for the padding alone -- a span-of-1 gamma_p call cost 2.2-2.5x
// a full vector on AVX2 for that reason. A live element's copy shares its
// region and its convergence, so the tail costs what its live lanes cost.
// Results are unchanged: every lane is computed by its own region core and
// frozen on its own mask, so a neighbour's value never reaches it (the
// lane-mix determinism checks in the smoke tests are the guard).
//
// On 2-lane targets (NEON, the SSE tiers) DriveUnary runs two vectors per
// iteration (#43): the short elementary chains are latency-bound at two
// lanes, and a second independent chain fills the FP pipes (~11% on NEON
// cos/sin). Each lane's arithmetic is unchanged, so results are
// bit-identical. Wider targets keep one vector per iteration -- the extra
// inlined kernel copy costs code size and MSVC build time for no
// measured gain.
//
// SHORT SPANS (v1.2.1, #110). Two things made a call on fewer elements than
// a vector holds wait for the call before it, so a loop of them ran at the
// kernel's latency and not its throughput (64.7 against 29.3 ns for a
// central erfinv on AVX3_ZEN4, quiet run of 2026-11-26,
// docs/bench-evidence/2026-11-26-zen4-v1.2.1). Both are a load that
// cannot be forwarded from the stores it overlaps, and so waits for them to
// retire:
//   * The tail used to be a masked LoadN/StoreN pair: a full vector's worth
//     of masked bytes read from memory the caller had just written with
//     8-byte stores. The tail is now assembled in registers (LoadTail and
//     StoreTail below), and a span of ONE element skips even that: the value
//     is broadcast, the kernel runs once, lane 0 is stored. A tail of one
//     element takes the same path (StoreTail). The kernel sees the same
//     live lanes either way and returns the same bits.
//   * The public wrappers passed std::span objects to the dispatched Impl.
//     Under the Windows x64 convention a span travels as a pointer to a
//     copy, and the copy was one 16-byte load of what the caller wrote as
//     two 8-byte stores. The Impls now take pointers and one length (see
//     any src/<fn>.cpp); the wrappers keep the debug length check.
// Either one alone leaves the call serialised. tests/test_span1.cpp guards
// the result for spans of 1 to 9 elements; tests/bench_span1.cpp measures
// the caller shapes.
//
// DriveScalar* are the scalar entry points (#42): one point broadcast to
// every lane, so a single call costs one vector in one region and nothing
// more, with no loop and no tail. Bit-identical to the span
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

// The tail of 1 .. N-1 elements, loaded and stored without a wide access to
// the caller's memory (v1.2.1). A masked load of a whole vector from a span
// the caller has just written cannot be forwarded from those narrower
// stores, and a masked store cannot be forwarded to the caller's reads of
// the results; either one makes the call wait for the one before it (a
// span of one or two erfinv values cost 65 ns a call on AVX3_ZEN4 before
// v1.2.1, 30.5 and 34.2 after; quiet run of 2026-11-26,
// tests/bench_span1.cpp). So each element is loaded by itself and blended
// into place in registers, and the results leave through a local buffer one
// element at a time. Lanes past the tail hold its LAST element: a live
// value, so the padding shares a real lane's region (#42).
template <class D>
HWY_INLINE op::V<D> LoadTail(D d, const double* p, size_t m) {
  auto v = op::Set(d, p[0]);
  for (size_t k = 1; k < m; ++k) {
    v = op::IfThenElse(op::FirstN(d, k), v, op::Set(d, p[k]));
  }
  return v;
}

template <class D>
HWY_INLINE void StoreTail(op::V<D> v, D d, double* p, size_t m) {
  // A tail of one element leaves through GetLane, as a span of one does
  // (LoadTail has already made it a broadcast). The buffer round trip below
  // cost erf and exp 1.5-2.3 ns a call on NEON and 1.1-1.7 on SSE4, where
  // every odd length ends in such a tail; this path reads level with the
  // v1.2.0 masked tail on both (docs/bench-evidence/2026-11-27-m1-v1.2.1-
  // tail1 and 2026-11-27-ivybridge-v1.2.1, pass3). A 2-3 element tail that
  // follows a full vector still pays 2-8 ns over the masked form on AVX2
  // (2026-11-26-kaby-v1.2.1/n7): the full vector's load has already waited
  // for the caller's stores, so the register work buys nothing there. Left
  // for v1.3.0 with the rest of that record.
  if (m == 1) {
    p[0] = op::GetLane(v);
    return;
  }
  double buf[op::MaxLanes(d)];
  op::Store(v, d, buf);
  for (size_t k = 0; k < m; ++k) p[k] = buf[k];
}

template <class Kernel>
static void DriveUnary(Kernel kernel, std::span<const double> in,
                       std::span<double> out) {
  HWY_DASSERT(in.size() == out.size());
  const op::ScalableTag<double> d;
  const size_t N = op::Lanes(d);
  const size_t n = out.size();
  const double* pi = in.data();
  double* po = out.data();
  if (n == 1) {
    po[0] = op::GetLane(kernel(d, op::Set(d, pi[0])));
    return;
  }
  size_t i = 0;
  if constexpr (op::MaxLanes(d) == 2) {
    // Two independent vectors per iteration on 2-lane targets (#43).
    for (; i + 2 * N <= n; i += 2 * N) {
      const auto v0 = kernel(d, op::Load(d, pi + i));
      const auto v1 = kernel(d, op::Load(d, pi + i + N));
      op::Store(v0, d, po + i);
      op::Store(v1, d, po + i + N);
    }
  }
  for (; i + N <= n; i += N) {
    op::Store(kernel(d, op::Load(d, pi + i)), d, po + i);
  }
  if (i < n) {
    const size_t m = n - i;
    StoreTail(kernel(d, LoadTail(d, pi + i, m)), d, po + i, m);
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
  if (n == 1) {
    po[0] = op::GetLane(kernel(d, op::Set(d, pa[0]), op::Set(d, pb[0])));
    return;
  }
  size_t i = 0;
  for (; i + N <= n; i += N) {
    op::Store(kernel(d, op::Load(d, pa + i), op::Load(d, pb + i)), d, po + i);
  }
  if (i < n) {
    const size_t m = n - i;
    StoreTail(kernel(d, LoadTail(d, pa + i, m), LoadTail(d, pb + i, m)), d,
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
  if (n == 1) {
    po[0] = op::GetLane(
        kernel(d, op::Set(d, pa[0]), op::Set(d, pb[0]), op::Set(d, pc[0])));
    return;
  }
  size_t i = 0;
  for (; i + N <= n; i += N) {
    op::Store(kernel(d, op::Load(d, pa + i), op::Load(d, pb + i),
                     op::Load(d, pc + i)),
              d, po + i);
  }
  if (i < n) {
    const size_t m = n - i;
    StoreTail(kernel(d, LoadTail(d, pa + i, m), LoadTail(d, pb + i, m),
                     LoadTail(d, pc + i, m)),
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
