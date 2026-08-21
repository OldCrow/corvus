# corvus — Agent Guide

C++20 library of SIMD-vectorized statistical special functions (erf/erfc
and inverses, lgamma, lbeta, digamma/trigamma, incomplete gamma/beta and
their inverses, Bessel I0/I1 with scaled variants) with runtime
multi-target dispatch via Google Highway.
Design goals: audited accuracy (documented ULP bounds per kernel per
target), clean-room provenance (MIT), swappable SIMD backend. Fills the
gap between SLEEF/Highway-contrib transcendentals and SciPy-level
special-function coverage.

## Reading map — load on demand, not preemptively
- Kernel, generator, reference/oracle, or accuracy work →
  `docs/NUMERICAL-DOCTRINE.md` (hazard rules, test doctrine, generator
  recipes, oracle doctrine, effort routing). BINDING for that work.
- Building, validating, benchmarking, CMake/CI edits, anything
  machine- or compiler-specific → `docs/ENVIRONMENT.md` (fleet table,
  toolchain caveats, tier capping/sweep recipes, CMake standard, CI
  design).
- Session state, decisions, open items, shipped-family records →
  `PLAN.md` (full design texts live in its git history).
- Audited accuracy claims → `docs/ACCURACY.md`; update it in the same
  change set as any kernel or gate change.
- Performance numbers → `docs/PERFORMANCE.md`. PROVISIONAL, and not the
  same standing as ACCURACY.md: one machine, one libm, two families not
  yet reproducible. Never quote it as a claim, and never compare its two
  tables against each other — only erf/erfc/lgamma are timed against a
  libm; the rest are corvus against corvus called per element.
- Anything user-facing — examples, API advice, explaining a bound to a
  consumer → `docs/USER-GUIDE.md` (what corvus does/does not provide,
  the returned-vs-composed accuracy rule, pair-selection rules of thumb).
  Keep it in sync when a bound or the public surface changes.
- Orienting on where a change belongs, or on what the facade contains →
  `docs/ARCHITECTURE.md` (band diagram + layer-by-layer notes). VISUAL
  REFERENCE ONLY: it restates the Architecture section below in picture
  form and adds no constraint the non-negotiables do not already carry.

## Architecture
- `include/corvus/corvus.h` — public API: `std::span` in/out, std-only,
  Doxygen. Highway never appears in public headers.
- `src/ops-inl.h` — the ~40-op SIMD facade; the ONLY file allowed to touch
  `hn::`. All kernels are written against `ops::` (aliased `op::` inside
  every kernel header), which is what keeps the backend swappable
  (std::simd later = reimplement this file, plus the one `hwy::TargetName`
  call behind `active_target()` in `src/erf.cpp`).
- `src/<fn>.cpp` — one TU per function family; the TU boundary is the
  sharing/dependency boundary (families consuming the same cores share a
  TU with multiple HWY_EXPORTs). Per-target pattern: `HWY_TARGET_INCLUDE`
  + `foreach_target.h`, kernel in `corvus::HWY_NAMESPACE`, `HWY_ONCE`
  section with `HWY_EXPORT` + public dispatch wrapper.
- `src/dd-inl.h` / `src/dd_special-inl.h` / `src/<fn>_dd-inl.h` —
  double-double primitives, shared dd specials (Log1pmxDd/Expm1Dd), and
  corvus-owned transcendental cores (exp_dd/log_dd, mantissa+exponent
  form so scaling rounds last).
- `tests/` — ctest gates against checked-in reference sets, registered in
  dependency order; test lengths non-multiples of lane counts so the
  masked-tail path is always exercised.

## Build & test
```sh
cmake --preset release -G Ninja   # only windows-clang-cl pins a generator; pass -G Ninja otherwise
cmake --build build
ctest --test-dir build --output-on-failure
```
Release for every perf/accuracy number. Options are `CORVUS_`-prefixed
(`CORVUS_DISABLED_TARGETS` tier capping, `CORVUS_SANITIZE`,
`CORVUS_DEV_WARNINGS`, `CORVUS_BUILD_TESTS`/`CORVUS_BUILD_EXAMPLES`, …). Highway via `find_package` or pinned
FetchContent — bump the pin only with a revalidation pass. Presets,
capping recipes, sweep scripts: `docs/ENVIRONMENT.md`.

## Non-negotiables (traps that bite silently; detail in the reference docs)
- Never touch `hn::` outside `src/ops-inl.h`. Kernels never allocate.
  Vector and tail are ONE masked code path — no scalar libm fallback.
- Call `HWY_DYNAMIC_DISPATCH` from inside `namespace corvus` — a
  single-target cap build breaks otherwise, and only the SSE2 sweep
  catches it.
- Outline region cores AND per-lane drivers (`HWY_NOINLINE`) from day one;
  MSVC codegen time is superlinear in function size.
- A new function family goes, in dependency position, into all FOUR
  gating lists: tests/CMakeLists.txt, the three ULP-report steps in ci.yml,
  and the `$gates` array in tools/sweep_tiers.ps1. ctest auto-discovers;
  the other three silently omit. A fifth, non-gating list — the default
  `-Targets` array in tools/quiet_bench.ps1 — enumerates the benches the
  same way and needs the family's `bench_*` too.
- Assert the tier, never assume it: validate under
  `CORVUS_EXPECT_TARGET=<tier>` and confirm the active target before
  trusting any tier result. Windows validation numbers come from clang-cl
  ONLY — MSVC silently caps at AVX2, and mingw GCC miscompiles by-value
  vector calls at AVX2 and above (GCC PR 126741; safe only for the
  128-bit tiers).
- Accuracy claims are made per SIMD tier only after native-silicon
  validation.
- Clean-room only — no ports of GPL/LGPL code. FP contraction is OFF
  project-wide; fusion is requested in source (`ops::MulAdd`), never
  inferred.
- Generated tables and reference files are checked in; regenerate only
  when method or point selection changes, re-run ULP gates after.
  Generators self-check and exit non-zero on a missed budget — trust the
  printed budget over any comment.

## Conventions
- Public/non-SIMD code: snake_case functions, CamelCase types (house
  style). Kernel code and the facade follow Highway idiom — CamelCase
  ops mirroring `hn::` names 1:1 (that mapping is what makes the future
  std::simd swap mechanical). Constants kCamelCase. `.h`/`.cpp`;
  `-inl.h` for per-target headers with the toggle guard.
- `include/corvus/` is the installed public surface; everything under
  `src/` is implementation.
- Public headers carry Doxygen. Kernel internals carry prose derivation
  blocks (method, error bounds, exactness arguments) at the definition
  site — that math commentary IS the maintainer documentation, and short
  lane-variable names (`d`, `ax`, `ssq`) are the numerical-kernel idiom.
- Every kernel documents its approximation source and accuracy bound at
  the definition site; provisional work is marked PROVISIONAL.
- Layout: 2-space indent in `src/` and `tests/` (Highway idiom), 4-space in
  `examples/` (consumer style); no `.clang-format` is enforced. Test-harness
  helpers in `tests/` use CamelCase (Highway idiom, grandfathered) — the
  snake_case rule above applies to the public API and examples.

## Effort routing (full table: docs/NUMERICAL-DOCTRINE.md)
Design, error budgets, and oracle construction are frontier work;
settled-design implementation is mid-tier; documented recipes and
bookkeeping are small-model work. ESCALATE the moment a recipe task
surfaces a decision — a ULP gate trips, tiers disagree, a fit misses, or
Highway behaves contrary to assumptions. De-escalate when high-effort
work reaches pure execution of a settled plan.

## Open items
See PLAN.md.
