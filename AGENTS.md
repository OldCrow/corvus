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

## Development Fleet

| Machine | OS | CPU | SIMD | Validation role |
|---|---|---|---|---|
| MacBook Pro 2017 (Kaby Lake) | macOS Ventura | i7-7820HQ | AVX2+FMA | AVX2 native; SSE4/SSSE3/SSE2 via tier capping (no FMA on those) |
| Mac Mini M1 | macOS Tahoe | Apple M1 | NEON (native FMA) | NEON validation |
| Asus TUF A16 | Windows 11 | Ryzen 7 7445 (Zen 4) | AVX-512 | AVX3* native; every lower x86 tier via capping |

At session start, verify which machine you are on (`uname -m`,
`sysctl -n machdep.cpu.brand_string`) before interpreting SIMD dispatch,
accuracy, or benchmark results — the active tier and FMA availability
change per machine.

## Build/Test/Run Commands
```sh
cmake -B build -G Ninja              # Release by default; Ninja preferred
cmake --build build
ctest --test-dir build --output-on-failure
```
- Generator: Ninja preferred (faster, identical behavior across macOS/
  Linux/Windows-with-vcvars); Unix Makefiles works, nothing depends on it.
- Build types (single-config default Release — house rule: perf and
  accuracy numbers from optimized builds only):
  - `Release` — benchmarks, accuracy validation, distribution
  - `RelWithDebInfo` — profiling (symbols for Instruments/perf)
  - `Debug` — debugger sessions; pair with `-DCORVUS_SANITIZE=address;undefined`
- Options (all `CORVUS_`-prefixed): `CORVUS_BUILD_TESTS` (top-level only
  by default), `CORVUS_DEV_WARNINGS` (-Wall -Wextra -Wpedantic; top-level
  only, never exported), `CORVUS_WERROR` (CI), `CORVUS_DISABLED_TARGETS`
  (tier capping), `CORVUS_SANITIZE`.
- Highway: uses system install if `find_package(hwy)` succeeds, else
  FetchContent of a pinned version (network on first configure). The pin
  tracks the version the accuracy audit ran against — bump only with a
  revalidation pass.

### CI
`.github/workflows/ci.yml`. Runner minutes are a budgeted resource even on
a public repo (lesson from libstats' temporary private phase) — keep the
surface lean and justify every runner:
- Linux x86-64 (cheap): sequential tier sweep AVX2->SSE2 in ONE job
  (Highway builds once; a matrix would multiply setup and Highway builds),
  plus ASan+UBSan Debug build in the same job.
- macOS arm64 (expensive class): single config — the only runner producing
  new information (native NEON silicon; counts as real-silicon validation).
- Docs-only changes skip CI; concurrency cancellation; per-job timeouts.
- Deliberately absent: Windows/MSVC (never built there yet — add when MSVC
  support lands via the Ryzen box); required status checks (incompatible
  with direct-push-to-main workflow); caching (build is ~minutes; add only
  if minutes grow). AVX-512 cannot run on hosted runners — Ryzen stays a
  manual validation stop.
- Workflow security: GITHUB_TOKEN read-only (repo setting + workflow
  `permissions:`), no event-payload interpolation in `run:` blocks,
  Dependabot keeps action versions fresh.

### CMake standard
- Target-first: no directory-scope `include_directories`/`link_libraries`/
  global flags; interface vs build separation via generator expressions
  (`$<BUILD_INTERFACE:>`/`$<INSTALL_INTERFACE:>`).
- Requirements that consumers inherit are PUBLIC on the target and must
  survive export: C++20 travels as `target_compile_features(corvus PUBLIC
  cxx_std_20)` (the header uses std::span), never as directory variables.
- Warnings are PRIVATE and top-level-gated — a consumer building corvus
  via FetchContent must never inherit our -Werror.
- The static lib builds with POSITION_INDEPENDENT_CODE for future pybind11
  bindings (pycorvus, per the pylibhmm/pylibstats pattern).
- `compile_commands.json` exported when top-level (clangd; note clangd
  still can't model foreach_target self-inclusion — spurious N_SSE4/N_AVX3
  diagnostics in kernel TUs are expected and harmless).
- Minimum CMake 3.25 (PROJECT_IS_TOP_LEVEL, FetchContent SYSTEM keyword).
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

## Workflows

Generators need mpmath in a throwaway venv (network for pip only):
```sh
python3 -m venv /tmp/mpv && /tmp/mpv/bin/pip install mpmath
/tmp/mpv/bin/python tools/gen_erf_table.py        > src/erf_data.inc
/tmp/mpv/bin/python tools/gen_erfc_tail_poly.py   > src/erfc_tail_data.h
/tmp/mpv/bin/python tools/gen_erf_reference.py    > tests/data/erf_reference.txt
/tmp/mpv/bin/python tools/gen_erfc_reference.py   > tests/data/erfc_reference.txt
```
Reference files and generated tables are checked in; regenerate only when
the method or point selection changes, and re-run the ULP tests after.

Per-tier validation recipe (run on each machine; caps only remove tiers):
```sh
BASE="HWY_AVX3_SPR|HWY_AVX3_ZEN4|HWY_AVX3_DL|HWY_AVX3"
for CAP in "$BASE" "$BASE|HWY_AVX2" "$BASE|HWY_AVX2|HWY_SSE4" \
           "$BASE|HWY_AVX2|HWY_SSE4|HWY_SSSE3"; do
  cmake -B build-cap -DCORVUS_DISABLED_TARGETS="$CAP" && \
  cmake --build build-cap && ctest --test-dir build-cap --output-on-failure
done
```

Benchmarks (`bench_*`, not ctest-registered): Release build, quiet machine
only; numbers taken on a loaded machine are indicative and must be labeled
as such in PLAN.md / docs.

## Model & Effort Routing

Hints for agent harnesses (sub-agent partitioning) and for flagging when
to switch model/effort — derived from how this project was actually built
(scaffold, facade, erf/erfc kernels, and governance were Fable 5 at high
effort; much of the surrounding work did not need that).

**High effort, frontier model** — wrong judgment is expensive to discover:
- New function/kernel design: approximation method choice (table vs
  polynomial vs rational), region splits, interval and degree selection,
  series length — and the error-budget analysis justifying them
  (cancellation, relative-vs-absolute error metrics, exactness arguments:
  Sterbenz, Fast2Sum, Dekker, FMA dependence).
- Diagnosing accuracy regressions from symptom to root cause (precedents:
  the erfc series-truncation metric mismatch; the non-FMA MulSub
  zero-residual hazard — both multi-step numerical reasoning).
- Public API changes, facade/seam design, governance and policy decisions.

**Default effort, mid-tier model** — pattern-following with judgment:
- Implementing a function whose design is settled, through the
  established pipeline (generator → kernel via ops:: → smoke test →
  reference set → ULP gate → bench).
- New tests/generators/benches copying existing patterns; CMake/CI
  adjustments within documented policy.

**Low effort / small model** — documented recipes and mechanical edits:
- Running the per-tier validation recipe on another machine and recording
  results in ACCURACY.md/PLAN.md.
- Regenerating tables/references unchanged; renames; badge/link fixes;
  PLAN.md bookkeeping after decisions are made.

**Escalate** (or flag a model/effort switch to the user) the moment a
recipe task surfaces a decision: a ULP gate trips, measured values differ
across tiers, a fit misses its target, or Highway behaves contrary to
assumptions. **De-escalate** (flag it) when a high-effort session reaches
pure execution of a settled plan — don't spend frontier effort on rote
work.

## Conventions
- C++20, zero public dependencies beyond std.
- Naming (decided 2026-07-21): public API and non-SIMD code use snake_case
  functions and CamelCase types (house style, matches libstats/libhmm) —
  `corvus::erf`, `corvus::active_target`. Per-target kernel code and the
  ops facade follow Highway idiom (CamelCase functions), and facade op
  names mirror `hn::` names 1:1 deliberately — that mapping is what makes
  the future std::simd swap mechanical. Constants are kCamelCase.
- File extensions: `.h`/`.cpp` (matches libstats/libhmm and Highway;
  `-inl.h` for per-target headers with the toggle guard).
- Header policy: `include/corvus/` is the installed public surface —
  std-only, no Highway types, Doxygen-documented. Everything in `src/` is
  implementation. When function families grow, split public headers per
  family under `include/corvus/` with `corvus.h` as the umbrella.
- Comment style: public headers carry Doxygen; kernel internals carry prose
  derivation blocks (method, error bounds, exactness arguments) at the
  definition site — that math commentary is the maintainer documentation,
  and short lane-variable names (`d`, `ax`, `ssq`) are the numerical-kernel
  idiom, documented by those blocks rather than per-variable naming.
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

## Documentation Map
- `README.md` — user-facing overview, build, design bullets.
- `docs/ACCURACY.md` — the audit record: per-function method, measured ULP
  bounds per validated tier, oracle methodology. Update in the same change
  set as any kernel or gate change.
- `PLAN.md` — agent-facing session state, decisions, open items.
- This file — conventions and workflows.
Keep documentation minimal beyond these four; resist growth.

## Open Items
See PLAN.md.
