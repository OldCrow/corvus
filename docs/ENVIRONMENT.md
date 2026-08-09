# corvus — Environment Reference

Loaded on demand via AGENTS.md's reading map. Covers the development
fleet, toolchains, build options, tier capping/validation recipes, the
CMake standard, and CI design. For kernel/generator/oracle rules see
`docs/NUMERICAL-DOCTRINE.md`; for session state see `PLAN.md`.

## Development Fleet

| Machine | OS | CPU | SIMD | Compiler | Validation role |
|---|---|---|---|---|---|
| MacBook Pro 2017 (Kaby Lake) | macOS Ventura | i7-7820HQ | AVX2+FMA | Apple Clang | AVX2 native; SSE4/SSSE3/SSE2 via tier capping (no FMA on those) |
| Mac Mini M1 | macOS Tahoe | Apple M1 | NEON (native FMA) | Apple Clang | NEON validation |
| Asus TUF A16 | Windows 11 | Ryzen 7 7445 (Zen 4) | AVX-512 | **clang-cl (GCC unsafe at AVX-512, MSVC can't dispatch it)** | AVX3* native; every lower x86 tier via capping |

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
— never make an AVX-512 claim from one. Use Clang (clang-cl), which
dispatches `AVX3_ZEN4` natively. mingw-w64 GCC also dispatches it but is
**disqualified for AVX-512 work as of 2026-07-28**: GCC 16.1 accesses the
ms_abi invisible-reference temporaries for 512-bit BY-VALUE arguments and
returns (any `__m512d` or wrapper struct passed to a non-inlined function
— exactly what HWY_NOINLINE outlining creates) with the aligned
`vmovapd`, while allocating them at plain `rsp`-relative offsets with no
realignment. The Windows ABI guarantees only 16-byte stack alignment, so
whether a binary faults is call-chain luck — `test_gamma_ulp` segfaulted
while the smoke test on the same kernel ran clean. Genuine register
spills are correctly `vmovupd`, and named over-aligned locals/return
slots in isolation get an aligned scratch pointer — only the argument
temporaries are broken. Reproduces at every -O level including -O0; no
flag rescues it (tested 2026-07-29: `-mstackrealign` and
`-mpreferred-stack-boundary=6` change nothing). clang-cl (and, per its
own escape-hatch caveat below, MSVC) is unaffected. Filed upstream
2026-08-08 as GCC PR 126741; minimal repro:
`C:\Users\gdwol\Development\gcc-zmm-mingw-repro\`. GCC remains fine for
capped tiers up to AVX2 (no zmm there), which is all
`tools/sweep_tiers.ps1` compiles — its `g++` default is safe for the
sweep, but the uncapped native build must be clang-cl (from a VS dev
shell so link.exe resolves; the sweep itself also runs clean under
`-CxxCompiler clang-cl -CCompiler clang-cl`).
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

## Build detail

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

## CMake standard

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

## Per-tier validation recipe

Run on each machine; caps only remove tiers. On Windows use
`tools/sweep_tiers.ps1`, which does the same thing, runs the gates
individually, and aborts on the first configure/build/gate failure — a
build failure otherwise leaves the *previous* tier's binaries in place
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

## CI

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
- House-style §5/§6 adopted at the first tagged release (v0.1.0,
  2026-08-06, closing issue #2): `lint-workflows.yml` carries the fleet
  actionlint + zizmor pair, and every action in both workflows is
  SHA-pinned with a `# vX.Y.Z` comment plus `persist-credentials: false`
  on every checkout. The pre-tag deferral (avoid the weekly Dependabot
  bump stream while the workflow surface was one rarely-changing file)
  is recorded in issue #2 for the rationale trail.
- Workflow security: GITHUB_TOKEN read-only (repo setting + workflow
  `permissions:`), no event-payload interpolation in `run:` blocks,
  Dependabot keeps action versions fresh across all SHA-pinned actions.
