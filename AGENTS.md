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

| Machine | OS | CPU | SIMD | Compiler | Validation role |
|---|---|---|---|---|---|
| MacBook Pro 2017 (Kaby Lake) | macOS Ventura | i7-7820HQ | AVX2+FMA | Apple Clang | AVX2 native; SSE4/SSSE3/SSE2 via tier capping (no FMA on those) |
| Mac Mini M1 | macOS Tahoe | Apple M1 | NEON (native FMA) | Apple Clang | NEON validation |
| Asus TUF A16 | Windows 11 | Ryzen 7 7445 (Zen 4) | AVX-512 | **clang-cl or mingw GCC, not MSVC** | AVX3* native; every lower x86 tier via capping |

**corvus deviates from the house Windows default, and this is deliberate.**
The other projects in this fleet (libhmm, libstats, ewcalc and the Python
bindings) use MSVC on Windows, matching Apple Clang on macOS and Clang on
Linux. corvus does not, for accuracy work: MSVC cannot dispatch AVX-512 at
all (see below), so an MSVC build cannot validate or benchmark this
project's widest tier — the one the Ryzen box exists to cover. Use
`clang-cl` by preference, since it keeps the MSVC ABI and so stays
link-compatible with the MSVC-built rest of the fleet (relevant if
libstats/libhmm ever adopt corvus); mingw-w64 GCC is the second option and
is not ABI-compatible with MSVC. MSVC remains fully supported *as a
consumer toolchain* — it compiles clean, passes every gate, and is what the
CI Windows job uses precisely because it is the strictest diagnostic gate
and the likeliest consumer default. The exception is about which compiler
produces *validation numbers*, not about which compilers must work.

At session start, verify which machine you are on (`uname -m`,
`sysctl -n machdep.cpu.brand_string`) before interpreting SIMD dispatch,
accuracy, or benchmark results — the active tier and FMA availability
change per machine.

**On the Ryzen box, the compiler decides whether AVX-512 exists at all.**
Highway puts every AVX3\* target in `HWY_BROKEN_TARGETS` under MSVC, so an
MSVC build silently tops out at AVX2 while still looking like a clean pass
— never make an AVX-512 claim from one. Use GCC (mingw-w64/UCRT) or Clang
(clang-cl); both dispatch `AVX3_ZEN4` natively and agree point-for-point.
Confirm the active set before trusting any tier result:
`build/_deps/highway-build/hwy_list_targets`. Also note `AVX3_SPR`
(Intel Sapphire Rapids) and `AVX10_2` are not available on Zen 4, so
"AVX3\* native" means `AVX3`, `AVX3_DL`, and `AVX3_ZEN4` only.

The MSVC blocklist is one unconditional entry in `hwy/detect_targets.h`
(`HWY_BROKEN_MSVC`), citing a 2016 codegen bug with no compiler-version
floor — unlike every sibling entry in that file, which all gate on a
version. `-DCORVUS_MSVC_UNBLOCK_AVX512=ON` overrides it via Highway's own
sanctioned `#ifndef` escape hatch. It is OFF by default and should stay
that way for anything whose numbers get published: upstream declares the
path untested, so results from it are ours to defend, not Highway's. The
override is scoped `PRIVATE` to the `corvus` target and needs no Highway
rebuild — `ChosenTarget::GetIndex` masks with the *calling* TU's
`HWY_CHOSEN_TARGET_MASK_TARGETS` by design, so corvus may legitimately
compile a wider target set than the linked `libhwy`. Note that CMake's
`MSVC` variable is also true for clang-cl, which Highway does *not*
blocklist (`HWY_COMPILER_CLANGCL` is a separate macro); the option
therefore checks `CMAKE_CXX_COMPILER_ID` and warns if it is a no-op.

When running mingw-built test binaries, invoke them from PowerShell or put
the WinLibs `mingw64/bin` first on `PATH` — a Git Bash shell puts Git for
Windows' own `libstdc++-6.dll` ahead of it, and the ABI mismatch segfaults
the test before it prints anything.

## Build/Test/Run Commands
```sh
cmake --preset release -G Ninja      # Ninja preferred; presets set no generator
cmake --build build
ctest --test-dir build --output-on-failure
```
Manual alternative (no preset): `cmake -B build -DCMAKE_BUILD_TYPE=Release -G Ninja`.
- Generator: Ninja preferred (faster, identical behavior across macOS/
  Linux/Windows-with-vcvars); Unix Makefiles works, nothing depends on it.
  Presets never pin a generator — pass `-G Ninja` alongside `--preset` (CI
  does the same).
- Build types (single-config default Release — house rule: perf and
  accuracy numbers from optimized builds only):
  - `Release` — benchmarks, accuracy validation, distribution
  - `RelWithDebInfo` — profiling (symbols for Instruments/perf)
  - `Debug` — debugger sessions; pair with `-DCORVUS_SANITIZE=address;undefined`
- Options (all `CORVUS_`-prefixed): `CORVUS_BUILD_TESTS` (top-level only
  by default), `CORVUS_DEV_WARNINGS` (-Wall -Wextra -Wpedantic; top-level
  only, never exported), `CORVUS_WERROR` (CI), `CORVUS_DISABLED_TARGETS`
  (tier capping), `CORVUS_SANITIZE`, `CORVUS_MSVC_UNBLOCK_AVX512` (OFF;
  see the MSVC/AVX-512 caveat above — configuring with it ON emits a
  deliberate `message(WARNING)`).
- Highway: uses system install if `find_package(hwy)` succeeds, else
  FetchContent of a pinned version (network on first configure). The pin
  tracks the version the accuracy audit ran against — bump only with a
  revalidation pass.

### CI
Fleet-wide workflow rules: [CI House Style](https://github.com/OldCrow/standards/blob/main/CI-HOUSE-STYLE.md)
— several of them were derived from this repo's workflow.

`.github/workflows/ci.yml`. Runner minutes are a budgeted resource even on
a public repo (lesson from libstats' temporary private phase) — keep the
surface lean and justify every runner:
- Linux x86-64 (cheap): sequential tier sweep AVX2->SSE2 in ONE job
  (Highway builds once; a matrix would multiply setup and Highway builds),
  plus ASan+UBSan Debug build in the same job.
- macOS arm64 (expensive class): single config — the only runner producing
  new information (native NEON silicon; counts as real-silicon validation).
- Windows x86-64 / MSVC (added 2026-07-24, when the project first ran on
  Windows at all): **toolchain coverage, explicitly not tier coverage.** It
  is the only job that sees MSVC-only diagnostics (`/W4 /WX`) and the only
  one using a multi-config generator — a `set_property(CACHE
  CMAKE_BUILD_TYPE)` bug broke the VS generator while Ninja stayed green,
  invisible to both other jobs. No `-G` is passed on purpose: CMake's
  Windows default is the newest installed Visual Studio, so the job follows
  runner images without pinning a version, needs no vcvars shell, and needs
  no third-party setup action. `CORVUS_EXPECT_TARGET=AVX2` is stable there
  without capping because `HWY_BROKEN_MSVC` pins the ceiling. Never add
  `CORVUS_MSVC_UNBLOCK_AVX512=ON` to it — that manufactures a green
  "Windows passing" that reads as AVX-512 validation on an untested code
  path and unguaranteed runner silicon.
- Docs-only changes skip CI; concurrency cancellation; per-job timeouts.
- Deliberately absent: required status checks (incompatible with
  direct-push-to-main workflow); caching (build is ~minutes; add only
  if minutes grow). AVX-512 cannot run on hosted runners — Ryzen stays a
  manual validation stop.
- Workflow security: GITHUB_TOKEN read-only (repo setting + workflow
  `permissions:`), no event-payload interpolation in `run:` blocks,
  Dependabot keeps action versions fresh.

### CMake standard

Full rules: [CMake House Style](https://github.com/OldCrow/standards/blob/main/CMAKE-HOUSE-STYLE.md)
in the fleet standards repo — corvus is its reference implementation, so this
section restates rather than deviates. It is self-sufficient for this repo.

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
  `cmake -B build-avx2 -DCORVUS_DISABLED_TARGETS="HWY_AVX10_2|HWY_AVX3_SPR|HWY_AVX3_ZEN4|HWY_AVX3_DL|HWY_AVX3"`
  (pipe-separated HWY_* macros; same idea as libstats' LIBSTATS_MAX_SIMD_TIER).
  Cap the *whole* AVX-512 family including `HWY_AVX10_2` — leaving it out
  works only for as long as `HWY_BROKEN_AVX10_2`'s compiler-version gate
  holds, and that one expires.
- **Assert the tier, never assume it.** `CORVUS_EXPECT_TARGET=<name>` in the
  environment makes every test and bench fail (exit 2, before doing any work)
  if runtime dispatch did not land on that target; unset means report-only.
  Use it with every cap — a cap that fails to bite otherwise leaves a green
  suite measuring a tier nobody asked for, which is exactly what
  `HWY_BROKEN_MSVC` does silently. CI sets it on every sweep iteration.
- Naming: Highway's "AVX3" (HWY_AVX3 and its _DL/_ZEN4/_SPR variants) means
  AVX-512 — Highway-internal terminology, not an Intel ISA name. Expect it
  in build output, ActiveTarget() strings, and CORVUS_DISABLED_TARGETS.
- `install` target only exists when Highway came from find_package (see
  CMakeLists comment and PLAN.md).
- Presets (`CMakePresets.json`, schema 6, min CMake 3.25 — matches this
  repo's existing minimum): `release` → `build/`, `debug` → `build-debug/`,
  `rel-with-debug` → `build-relwithdebinfo/`, plus the `sanitize` extra →
  `build-san/` (Debug + `CORVUS_SANITIZE=address;undefined`, own binaryDir
  so it never leaves a sticky cache variable in `build/`). No `generator`
  field in any preset; pass `-G Ninja` alongside `--preset`.

## Architecture
- `include/corvus/corvus.h` — public API. Plain spans, no Highway types.
- `src/ops-inl.h` — the SIMD op facade (per-target include guard). The only
  file allowed to touch `hn::`. ~20 ops: load/store (incl. masked N
  variants), arithmetic, fma, abs/min/max/copysign, compare/select, exp,
  reductions.
- `src/<fn>.cpp` — one translation unit per function family. The TU
  boundary is the SHARING/DEPENDENCY boundary, not one-symbol-per-file:
  functions that consume the same kernel cores stay in one TU with
  multiple HWY_EXPORTs (erfinv + erfcinv both route onto both shared
  cores — splitting would instantiate every core twice per target for
  nothing); split within a family only when dependency sets differ
  materially (erf.cpp stays free of erfc's dd/exp_dd/tail-data
  dependencies, and a consumer linking only erf pulls only erf.o).
  Pattern:
  `HWY_TARGET_INCLUDE` + `foreach_target.h`, kernel in
  `corvus::HWY_NAMESPACE` written against `ops::`, then `HWY_ONCE` section
  with `HWY_EXPORT` + public dispatch wrapper. Call `HWY_DYNAMIC_DISPATCH`
  from *inside* `namespace corvus` — with a single compiled target (the
  SSE2 cap) Highway collapses it to `N_SSE2::FUNC`, and a globally
  qualified call then names a namespace that does not exist. It compiles
  at every other tier, so the cap sweep is what catches it.
- `src/dd-inl.h` — double-double primitives (Fast2Sum, TwoSum, TwoProd,
  DdAdd/DdMul/DdRecip) shared by the compensated kernels. Written against
  `ops::` like everything else. Exact residuals go through `ops::ProdLow`,
  never a bare `MulSub` — same FMA-capability hazard as `ops::SquareLow`.
- `src/<fn>_dd-inl.h` — corvus-owned transcendental cores (exp_dd, later
  log_dd), internal only, no public API. They return
  mantissa + exponent so a consumer folds its own factors in before the
  power-of-two scaling rounds anything; that is what keeps a subnormal
  result at one rounding.
- `tests/` — ctest executables comparing against libm/reference values;
  test lengths deliberately non-multiples of lane counts to exercise the
  masked-tail path. A test for an *internal* kernel compiles the kernel
  header itself through foreach_target and so uses
  `corvus_kernel_test_target()`, which links `hwy::hwy`, adds the source
  root, and — the part that matters — applies `CORVUS_HWY_TARGET_DEFS` so
  the test sees the same target set as the library. Such a test also
  asserts its own dispatched target equals `corvus::active_target()`.
  Where the kernel carries more than working precision, the reference file
  carries a double-double pair and the test measures relative error below
  the last bit of a double; rounding first would hide what is being
  tested.

## Workflows

Generators need mpmath in a throwaway venv (network for pip only):
```sh
python3 -m venv /tmp/mpv && /tmp/mpv/bin/pip install mpmath
/tmp/mpv/bin/python tools/gen_erf_table.py        > src/erf_data.inc
/tmp/mpv/bin/python tools/gen_erfc_tail_poly.py   > src/erfc_tail_data.h
/tmp/mpv/bin/python tools/gen_erf_reference.py    > tests/data/erf_reference.txt
/tmp/mpv/bin/python tools/gen_erfc_reference.py   > tests/data/erfc_reference.txt
/tmp/mpv/bin/python tools/gen_exp_table.py        > src/exp_dd_data.inc
/tmp/mpv/bin/python tools/gen_exp_dd_reference.py > tests/data/exp_dd_reference.txt
/tmp/mpv/bin/python tools/gen_log_table.py        > src/log_dd_data.inc
/tmp/mpv/bin/python tools/gen_log_dd_reference.py > tests/data/log_dd_reference.txt
/tmp/mpv/bin/python tools/gen_lgamma_data.py      > src/lgamma_data.h
/tmp/mpv/bin/python tools/gen_lgamma_reference.py > tests/data/lgamma_reference.txt
/tmp/mpv/bin/python tools/gen_erfinv_data.py      > src/erfinv_data.h
/tmp/mpv/bin/python tools/gen_erfinv_reference.py
```
`gen_erfinv_reference.py` writes both `tests/data/erfinv_reference.txt` and
`tests/data/erfcinv_reference.txt` directly (two output files, so no `>`
redirection) rather than printing one file to stdout.
Reference files and generated tables are checked in; regenerate only when
the method or point selection changes, and re-run the ULP tests after.
Table generators self-check on every run and exit non-zero rather than emit
a table that misses its error budget — `gen_exp_table.py` re-derives the
whole budget (reduction exactness, polynomial truncation, table dd error)
onto stderr. Trust that line over any claim in a comment.

Per-tier validation recipe (run on each machine; caps only remove tiers).
On Windows use `tools/sweep_tiers.ps1`, which does the same thing, runs the
gates individually, and aborts on the first configure/build/gate failure —
a build failure otherwise leaves the *previous* tier's binaries in place
and the next iteration re-measures them under the new tier's name:
```sh
BASE="HWY_AVX10_2|HWY_AVX3_SPR|HWY_AVX3_ZEN4|HWY_AVX3_DL|HWY_AVX3"
for TIER in AVX2 SSE4 SSSE3 SSE2; do
  case $TIER in
    AVX2)  CAP="$BASE" ;;
    SSE4)  CAP="$BASE|HWY_AVX2" ;;
    SSSE3) CAP="$BASE|HWY_AVX2|HWY_SSE4" ;;
    SSE2)  CAP="$BASE|HWY_AVX2|HWY_SSE4|HWY_SSSE3" ;;
  esac
  cmake -B build-cap -DCORVUS_DISABLED_TARGETS="$CAP" && \
  cmake --build build-cap && \
  CORVUS_EXPECT_TARGET="$TIER" ctest --test-dir build-cap --output-on-failure
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
- Masked-off lanes still EXECUTE every op, including gathers. Any gather
  whose index derives from lane VALUES (e.g. erf's round(ac*256); unlike
  log_dd's bit-masked slot index, which is bounded by construction) must
  have its input NaN/domain-scrubbed first -- a discarded lane's NaN
  otherwise reaches the index or not by platform accident (x86 minpd
  drops NaN, ARM fmin propagates it and fcvtzs(NaN) = 0), and Highway's
  debug bounds assert can trip. erfc.cpp's nan mask and erfinv's
  HalleyMid scrub are the pattern.
- The non-FMA fallback in ops::SquareLow/ProdLow is Dekker's split, and its
  intermediate a*(2^27+1) OVERFLOWS for |a| > 2^996 (~6.7e299). A kernel
  whose operands can reach that range must scale by a power of two first
  and scale back after -- exact, and linear in the operand, so it stays one
  code path for every target rather than a non-FMA special case. lgamma's
  Stirling product is the first kernel to need it; it read as correct on
  every FMA target and only the capped SSE sweep exposed it.

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
