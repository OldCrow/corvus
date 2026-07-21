# corvus — Agent Guide

## Overview
corvus is a C++20 library of SIMD-vectorized statistical special functions
(erf/erfc, lgamma, incomplete gamma/beta, digamma, Bessel I0/I1, and their
inverses) with runtime multi-target dispatch via Google Highway. It fills the
gap between basic-transcendental SIMD libraries (SLEEF, Highway contrib math)
and SciPy-level special-function coverage. Design goals: audited accuracy
(documented ULP bounds per kernel per target), clean license provenance
(clean-room implementations only), and a swappable SIMD backend.

Highway is an implementation detail: public headers (`include/corvus/`)
expose only `std::span`/pointer APIs. All kernels are written against the
op facade in `src/ops-inl.h` — never against `hn::` directly — so the
backend can later be replaced (e.g. by `std::simd`) by reimplementing that
one file. Runtime dispatch (HWY_EXPORT/HWY_DYNAMIC_DISPATCH) has no
std::simd equivalent and will outlive the op-layer migration.

## Build/Test/Run Commands
```sh
cmake -B build                       # Release by default
cmake --build build
ctest --test-dir build --output-on-failure
```
- Highway: uses system install if `find_package(hwy)` succeeds, else
  FetchContent of pinned 1.2.0 (network needed on first configure).
- Tier capping for native per-tier validation:
  `cmake -B build-avx2 -DCORVUS_DISABLED_TARGETS="HWY_AVX3|HWY_AVX3_DL|HWY_AVX3_ZEN4|HWY_AVX3_SPR"`
  (pipe-separated HWY_* macros; same idea as libstats' LIBSTATS_MAX_SIMD_TIER).
- Naming: Highway's "AVX3" (HWY_AVX3 and its _DL/_ZEN4/_SPR variants) means
  AVX-512 — Highway-internal terminology, not an Intel ISA name. Expect it
  in build output, ActiveTarget() strings, and CORVUS_DISABLED_TARGETS.
- `install` target only exists when Highway came from find_package (see
  CMakeLists comment and PLAN.md).

## Architecture
- `include/corvus/corvus.h` — public API. Plain spans, no Highway types.
- `src/ops-inl.h` — the SIMD op facade (per-target include guard). The only
  file allowed to touch `hn::`. ~20 ops: load/store (incl. masked N
  variants), arithmetic, fma, abs/min/max/copysign, compare/select, exp,
  reductions.
- `src/<fn>.cpp` — one translation unit per function family. Pattern:
  `HWY_TARGET_INCLUDE` + `foreach_target.h`, kernel in
  `corvus::HWY_NAMESPACE` written against `ops::`, then `HWY_ONCE` section
  with `HWY_EXPORT` + public dispatch wrapper.
- `tests/` — ctest executables comparing against libm/reference values;
  test lengths deliberately non-multiples of lane counts to exercise the
  masked-tail path.

## Conventions
- C++20, zero public dependencies beyond std.
- Kernels never allocate; caller owns all buffers.
- Every kernel documents its approximation source and accuracy bound at the
  definition site. Provisional implementations are marked PROVISIONAL.
- Clean-room only: no ports of GPL/LGPL code (lesson from libstats issue
  #67, vector_erf_neon LGPL provenance).
- Accuracy claims are made per SIMD tier only after native validation on
  real silicon (lesson from libstats issue #74, SSE2 subnormal bug invisible
  under Rosetta). Fleet: M1 (NEON), Ryzen 7445 (AVX-512 down to SSE2 via
  tier capping), i7-7820HQ Kaby Lake (AVX2).
- Vector and tail paths must be the same code path (masked LoadN/StoreN),
  never a scalar libm fallback for the tail.
- Any op whose CORRECTNESS depends on FMA fusion (exact residuals like
  fma(a,b,-fl(a*b))) must be capability-guarded in the facade -- Highway
  emulates MulAdd/MulSub as mul-then-add on non-FMA targets (SSE2/SSSE3/
  SSE4), which silently zeroes exact residuals. See ops::SquareLow.

## Open Items
See PLAN.md.
