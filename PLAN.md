# corvus — Plan / Session State

## Status [DERIVED] — end of session 2026-07-24 (Ryzen)
Scaffolded 2026-07-20 on the Kaby Lake Mac (AVX2); public at
github.com/OldCrow/corvus. Two functions shipped, both production-quality:

- **erf**: clean-room table + local-Taylor kernel (ported from libstats
  vector_erf_neon through the ops facade). Max 1 ULP over the full domain.
- **erfc**: two-region kernel (core reuses the erf table via compensated
  assembly; tail is a fitted e^{-a^2}*G(1/a)/a). Max 1 ULP for |x| <= 6
  and subnormal results; max 5 ULP normal-tail (bounded by backend Exp).

**Validated tiers, both functions, all within the same ULP bounds**: AVX2,
SSE4, SSSE3, SSE2 (Kaby Lake native + CORVUS_DISABLED_TARGETS capping) and
**NEON** (Apple Silicon GitHub Actions runner, native FMA — validated in CI
2026-07-21). NEON and AVX2 produce bit-identical ULP results point-for-
point on both reference sets (see docs/ACCURACY.md) — first evidence the
kernels are deterministic across ISA/compiler/OS on FMA-capable targets.
Bit-exactness holds among the FMA tiers specifically; the no-FMA SSE tiers
meet the same bounds but differ marginally in not-correctly-rounded counts
(quantified on the Ryzen 2026-07-24 — see the open item below).

**AVX-512 closed 2026-07-24 on the Ryzen 7445** — see the resolved open
item below and docs/ACCURACY.md. `AVX3_ZEN4` (native), `AVX3_DL`, and
`AVX3` all pass every gate with values identical to AVX2/NEON, reproduced
under both GCC and clang-cl; the capped sweep down through AVX2/SSE4/
SSSE3/SSE2 ran on that box as well. No accuracy tier gaps remain for erf
and erfc on available hardware. (`AVX3_SPR` is Intel-only and `AVX10_2`
absent on Zen 4 — neither is validatable on this fleet.)

**Repo infrastructure is fully stood up**: CI (Linux tier-sweep+sanitizers,
macOS arm64/NEON), branch protection, security scanning, topics, CMake
standard (exported C++20, private warnings, PIC), naming/docs conventions,
and model/effort routing hints (AGENTS.md) — all decided this session per
the explicit front-load-don't-evolve approach (see
[[frontload-project-conventions]] in the user's memory). Dependabot's
first PR (#1, actions/checkout v4->v7) merged clean after one rebase.

Facade ops added this session: SignedTag, Round, ConvertToInt, ShiftLeft,
GatherIndex, Ge, Gt, Eq, IsNaN, MulSub, SquareLow (FMA-capability-guarded),
AllFalse/AllTrue.

## Decisions
- Name: corvus (OldCrow tie-in). Namespace `corvus::`.
- Scope: statistical special functions only. Basic transcendentals
  (exp/log/trig/pow) are out of scope — Highway contrib math owns those.
  Planned families (P0 first): erf/erfc, erfinv/erfcinv, lgamma,
  regularized incomplete gamma P/Q, regularized incomplete beta; (P1)
  digamma, inverse incomplete gamma/beta, Bessel I0/I1.
- Backend: Highway, hidden behind `src/ops-inl.h` facade; public API is
  std-only. std::simd migration is a facade-reimplementation, deferred
  until implementations mature (as of mid-2026: GCC 16 partial — no
  simd.loadstore, partial simd.math; no libc++ implementation).
- Dependency model: find_package(hwy) preferred, FetchContent 1.2.0
  fallback.
- Ship model: static lib, Highway not exposed to consumers.
- License: MIT. Clean-room implementations only.
- Conventions audit (2026-07-21, user-driven): public API renamed to
  snake_case (corvus::erf/erfc/active_target — house style; the earlier
  CamelCase was accidentally inherited from Highway). Kernel/facade code
  stays Highway-idiom CamelCase, facade names deliberately mirror hn:: 1:1.
  .h/.cpp extensions confirmed. Single-public-header policy confirmed and
  documented (per-family split later). Doxygen on public header only;
  prose derivation blocks are the internal documentation. All recorded in
  AGENTS.md Conventions.
- Documentation set fixed at four files (README, docs/ACCURACY.md, PLAN.md,
  AGENTS.md) — resist growth. ACCURACY.md is the public audit record and
  must move with kernel/gate changes.
- Build-system audit (2026-07-21, user-driven): C++20 now propagates to
  consumers via target_compile_features PUBLIC (was a directory variable —
  invisible to the export, a real bug); dev warnings -Wall -Wextra
  -Wpedantic private + top-level-gated, CORVUS_WERROR for CI (build is
  warning-clean); PIC on for future pybind11 bindings; Highway FetchContent
  pin bumped 1.2.0 -> 1.4.0 to match the audited version; Ninja preferred
  generator; CORVUS_SANITIZE option added; build-type roles documented.
  CMake standard recorded in AGENTS.md. Deferred deliberately: LTO/IPO
  (profile first), shared-lib + symbol visibility (no demand yet),
  install-when-fetched (existing open item).
- Platform tiers: Tier 1 (accuracy-audited on real silicon) = NEON (M1),
  AVX-512/AVX2/SSE2 (Ryzen 7445 native + CORVUS_DISABLED_TARGETS capping),
  AVX2 (Kaby Lake). Tier 2 (compiles, unaudited) = SVE and anything else
  Highway emits.
- dd transcendental core (2026-07-21, resolved with user): build a
  corvus-owned double-double exp_dd + log_dd as its own phase BEFORE
  lgamma, with the erfc tail rewire as the acceptance test (existing
  reference sets and ULP gates; expected tail bound 5 -> ~1-2 ULP).
  Rationale: lgamma needs a compensated LOG, not exp (the earlier open
  item conflated these); the incomplete gamma/beta prefactor
  exp(a·log x − x − lgamma(a)) needs both; core-first avoids shipping
  lgamma bounded by hn::Log and revalidating later, and shrinks the
  audit's coupling to the Highway pin (contrib Exp/Log drop out of the
  accuracy-critical path). Placement: shared kernel headers
  (src/exp_dd-inl.h, src/log_dd-inl.h) against ops:: — NOT facade ops
  (facade stays a 1:1 hn:: mirror); one facade addition: ops::BitCast.
- lgamma v1 scope (2026-07-21, resolved with user): full real axis —
  negative x via reflection + sinpi, poles return +inf. Sign output
  (C's signgam) deferred; SciPy's gammaln offers none either.

## Erfc [DERIVED, 2026-07-21]
Two-region kernel: core |x| <= 6 reuses the erf table via compensated
1 -/+ erf assembly; tail 6 < a <= 28 is e^{-a^2}*G(1/a)/a with per-interval
polynomial fits of G (tools/gen_erfc_tail_poly.py; intervals [6,10],[10,17],
[17,28], degrees 11/10/8, coefficient-select + single Horner). Exact a^2
split via ops::SquareLow. Per-vector AllFalse/AllTrue branch skips the
unused region path.

Accuracy (Kaby Lake, all four tiers AVX2/SSE4/SSSE3/SSE2 identical):
core max 1 ULP; tail normal-result max 5 ULP (~59% not correctly rounded --
entirely the backend hn::Exp contribution); tail subnormal-result max 1 ULP.
Gates in test_erfc_ulp set to measured values, no margin.

Bench (Kaby Lake, session-loaded, indicative): core-dominated 3.5-3.8x vs
libm erfc; tail-only 2.0-2.4x (AVX2), 1.2-1.9x (SSE2); mixed 1.2-2.1x.

Two design findings worth remembering:
1. The 5-term erf series was tuned for erf's ABSOLUTE error (vs ulp(1));
   erfc needs the same series to RELATIVE erfc precision near a = 6, where
   the c6 truncation term is ~2e-13 relative. Fix: extended the shared
   series to d^8 (c6-c8 Hermite closed forms verified against the
   generator's recurrence; erf results unchanged, still max 1 ULP).
2. ops::MulSub is NOT an exact-residual primitive on non-FMA targets
   (SSE2/SSSE3/SSE4): Highway emulates it as mul-then-sub, which silently
   returns 0 for fma(a,a,-fl(a*a)) and reintroduced the amplified-argument
   error (501 ULP at a~25). Fix: ops::SquareLow is capability-guarded --
   HWY_NATIVE_FMA ? MulSub : Dekker split. Rule: any op whose CORRECTNESS
   (not just accuracy) depends on fusion must be guarded in the facade.

## dd transcendental core + lgamma [design resolved 2026-07-21]
Phase A — exp_dd / log_dd (internal only, dd in / dd out):
- exp_dd: Cody-Waite reduction k = round(x·N/ln2) with split constant
  L1+L2 so k·L1 is exact (|k| headroom far beyond erfc's a² <= ~750);
  N = 64 or 128 dd table of 2^(j/N) (two GatherIndex from separate
  hi[]/lo[] arrays — same pattern as the erf table); degree 5-6 poly for
  e^r − 1 on |r| <= ln2/2N; assemble T·(1+p) with one Fast2Sum;
  two-stage 2^k scaling through gradual underflow (the existing erfc
  subnormal reference band tests this for free). Budget ~2^-60 relative
  before final rounding — faithfully rounded.
- log_dd: exponent extraction (needs ops::BitCast) + table
  {R_j, L_hi, L_lo}; r = fma(R_j, m, −1); log x = k·ln2_dd + L_j_dd +
  log1p-poly(r), accumulated in dd.
- Generators tools/gen_exp_table.py / gen_log_table.py with mpmath
  self-check each run (house pipeline).
- Acceptance: erfc tail on exp_dd, retighten test_erfc_ulp tail gate to
  measured (expect <= 2 ULP); ACCURACY.md moves in the same change set.
  Bench may give back some of the 1.55-1.6x tail speedup (two gathers +
  dd math vs contrib poly) — measure and label.

Phase B — lgamma (public corvus::lgamma, span API):
- Regions: zeros zone [~0.75, ~2.5] — polys centered exactly at x = 1
  and x = 2 in exact t (Sterbenz), form t·(−γ + t·q(t)) with dd leading
  term, for RELATIVE accuracy at the zeros; middle (zone end, X0) —
  masked fixed-step product recurrence (<= 6 Select-multiply steps,
  P <= ~Γ(X0), then one log_dd(P)); Stirling x >= X0 —
  (x−½)·log_dd(x) − x + ½log(2π)_dd + φ, φ a Chebyshev fit in 1/x²
  (generator sweeps X0 in {8, 10, 13} × degree for 2^-60); (0,1) — one
  shift lgamma(x+1) − log_dd(x); negative axis — reflection
  log(π/|sin πx|) − lgamma(1−x), sinpi via exact u = x − round(x) plus
  poly, u == 0 mask -> +inf at poles.
- Specials: lgamma(1) = lgamma(2) = +0 exactly (falls out of exact-t
  form); +inf at poles, x = +inf, and overflow (x >~ 2.55e305); NaN
  propagates.
- Targets: <= 2 ULP relative on the positive axis including the zeros;
  negative axis measured per region, documented degradation near the
  |Γ| = 1 crossings (SciPy-level behavior, but documented).
- Considered and rejected (record in derivation blocks): Lanczos
  (plateaus ~1e-16 relative, can't reach the ULP target); per-lane
  iterated recurrence (lane divergence + per-step rounding).
- Left to generator/bench experiments: table sizes N, X0 and degrees,
  exact zone boundaries.

## GitHub Repo Settings [DERIVED, applied 2026-07-21 via gh api]
- Merge: all three styles allowed (matches libstats); auto-delete head
  branches on merge (matches libstats).
- Features: wiki and projects DISABLED (docs policy caps documentation at
  four files; a wiki would route around that). Issues on, discussions off.
- Topics: cpp, cpp20, simd, math, special-functions, avx2, avx512, neon,
  vectorization.
- Security: Dependabot alerts + automated security fixes, secret scanning
  + push protection, private vulnerability reporting — all enabled.
- Ruleset "protect-main": blocks force-push and deletion of the default
  branch; direct pushes remain allowed (solo workflow).
- Actions: default GITHUB_TOKEN read-only, cannot approve PRs (hardened
  before any CI exists).
- Deferred/manual: (a) add required-status-checks to the ruleset when CI
  lands; (b) signed-commits rule — all corvus commits so far verify (G),
  but confirm the M1 and the Ryzen/Windows box sign before enabling or
  their pushes will be rejected; (c) tag-protection ruleset for v* at
  first release.

## Open Items
- [OPEN, design resolved 2026-07-21] erfc tail 5-ULP bound is entirely
  hn::Exp's contribution. Resolution: dd transcendental core (exp_dd +
  log_dd) as Phase A before lgamma, erfc tail rewire as acceptance test —
  see the design section above. (Note: lgamma consumes log_dd, not
  exp_dd; exp_dd's next consumer after erfc is the incomplete gamma/beta
  prefactor.)
- [RESOLVED 2026-07-21] NEON validated in CI (Apple Silicon runner) for
  both erf and erfc, bit-identical to AVX2 results.
- [RESOLVED 2026-07-24] Native AVX-512 validation on the Ryzen 7445HS
  (Zen 4), Windows 11, Highway 1.4.0 — the last tier gap for erf/erfc is
  closed. Full sweep per the AGENTS.md recipe; ACCURACY.md matrix and
  cross-arch note updated. Results:
  - FMA tiers `AVX3_ZEN4` (native dispatch), `AVX3_DL`, `AVX3`, `AVX2`:
    erf max 1 ULP / 217 not-CR; erfc 1/5/1 ULP with 46/4917/12 not-CR.
    **Bit-identical to the NEON and Kaby-Lake/Linux AVX2 numbers**, so the
    cross-arch determinism claim now extends to AVX-512.
  - No-FMA tiers `SSE4`, `SSSE3`, `SSE2`: erf 218 not-CR, erfc 50/4916/12.
    Max ULP unchanged on every tier, so all gates hold. This is the
    Dekker/`ops::SquareLow` no-FMA path rounding differently — the
    previously unverified case, now measured rather than assumed.
  - Reproduced under two compilers, GCC 16.1 (mingw-w64/UCRT) and Clang
    22.1.8 (clang-cl, MSVC ABI), agreeing point-for-point.
  - **Toolchain constraint worth remembering: a default MSVC build cannot
    perform this validation.** Highway lists every AVX3\* target in
    `HWY_BROKEN_TARGETS` under MSVC, so an MSVC build silently caps at
    AVX2 and looks like a pass. Use GCC or Clang on this box for any
    AVX-512 claim. (Verify with
    `build/_deps/highway-build/hwy_list_targets`.) See the blocklist item
    below for why that cap is overridable but stays on by default.
  - Gotcha that cost time: running mingw-built test exes from a Git Bash
    shell segfaults in `libstdc++-6.dll` — Git for Windows' own mingw64
    runtime shadows the WinLibs one on PATH. Not a corvus bug; put the
    WinLibs `mingw64/bin` first, or run from PowerShell.
- [RESOLVED 2026-07-24] Windows/MSVC CI job added (`windows-msvc` in
  ci.yml) — scoped as toolchain coverage, not tier coverage, since
  `HWY_BROKEN_MSVC` pins MSVC dispatch to AVX2 regardless of runner silicon
  (which is what makes its `CORVUS_EXPECT_TARGET=AVX2` assertion stable
  without capping). Uses CMake's default Windows generator rather than a
  pinned `-G`: that yields the VS multi-config generator, tracks runner
  images without a version pin, and needs neither a vcvars shell nor a
  third-party setup action. Rationale for the job at all: this session
  produced three Windows-only defects nothing else could catch — the two
  C4566 Unicode-docstring warnings in the sibling Python bindings, and the
  multi-config configure error that Ninja-only CI was blind to.
  - **It justified itself on the first local replay**: `/WX` promoted MSVC's
    C4996 `std::getenv` deprecation in the then-new `tests/expect_target.h`
    into a hard error, breaking all six test/bench targets. Fixed with a
    locally scoped `#pragma warning(disable : 4996)` rather than
    `_CRT_SECURE_NO_WARNINGS`, so the deprecation stays live elsewhere. A
    GCC/Ninja `-Werror` build and the full VS-generator build+ctest+ULP
    sequence were both re-run green afterwards.
- [RESOLVED 2026-07-24] CI asserted its own tiers. Three linked defects, all
  the same shape as the MSVC blocklist — a cap that doesn't bite leaves a
  green suite measuring the wrong thing:
  - The sweep's first iteration used `CAP=""` (uncapped), so AVX2 coverage
    depended on the runner's CPU. On an AVX-512-capable hosted runner that
    iteration would have dispatched AVX3\* and **AVX2 would never have been
    exercised at all** — the next cap jumps to SSE4. Latent, not yet hit:
    logs from run 30052628221 show `AVX2`, so the draws so far have been
    AVX2-class. Now every tier is capped explicitly by name.
  - `BASE` omitted `HWY_AVX10_2`, so every capped iteration could be topped
    by AVX10.2. Harmless only while `HWY_BROKEN_AVX10_2`'s compiler-version
    gate holds — and unlike `HWY_BROKEN_MSVC`, that gate expires. Added.
  - Nothing verified the reached tier. New `CORVUS_EXPECT_TARGET` env var
    (`tests/expect_target.h`, used by all four tests and both benches) fails
    with exit 2 *before* doing any work if dispatch missed the expectation;
    unset = report-only, so local runs need no setup. Wired into both CI
    jobs, including macOS/NEON — that job is the sole source of the NEON row
    in ACCURACY.md, so a runner-image change moving dispatch to NEON_BF16
    must fail loudly rather than silently relabel the audit record.
    Verified: match passes, mismatch fails 4/4, unset reports only.
- [RESOLVED 2026-07-24] Highway's MSVC AVX-512 blocklist — investigated,
  and `CORVUS_MSVC_UNBLOCK_AVX512` added (default OFF) [DERIVED].
  - Mechanism: `HWY_BROKEN_MSVC` in `hwy/detect_targets.h` expands to
    `(HWY_AVX3 | (HWY_AVX3 - 1))`. Highway's target bits run *descending*
    in capability, so that mask removes AVX3 **and everything above it**
    (`AVX3_DL`, `AVX3_ZEN4`, `AVX3_SPR`, `AVX10_2`), leaving AVX2 as the
    ceiling. It is the only entry in that file with no compiler-version
    floor, and it cites a 2016-era bug report.
  - Measured with the block lifted, MSVC 19.51: builds clean (0 warnings),
    dispatches `AVX3_ZEN4`, and reproduces every documented ULP value
    exactly — so no correctness problem in *these* kernels. Throughput
    erf n=1e6: MSVC AVX2 9.02 → MSVC AVX3_ZEN4 5.49 ns/el (1.64x), erfc
    core 3.13 → 2.41. Separately note MSVC trails GCC ~2.8x at equal tier
    (GCC AVX3_ZEN4 1.95 ns/el), which is a codegen-quality gap, not an
    AVX-512 one.
  - Decision: keep Highway's default. Upstream declares the path untested,
    so numbers from it would be ours to defend; `clang-cl` gets AVX-512 on
    Windows with MSVC ABI and no override. The option exists so the
    override is a recorded, warned, reviewable switch rather than an
    ad-hoc `CMAKE_CXX_FLAGS` hack — the earlier form of this experiment
    clobbered CMake's MSVC defaults and lost `/EHsc`, which is exactly the
    failure mode a named option prevents.
  - Scoping detail: the define is `PRIVATE` to `corvus` and needs no
    Highway rebuild. `ChosenTarget::GetIndex` masks with the calling TU's
    `HWY_CHOSEN_TARGET_MASK_TARGETS` by design, so per-module target sets
    may differ — confirmed empirically against an unmodified `libhwy`.
  - [OPEN] Upstream path if ever worth it: add a version floor to
    `HWY_BROKEN_MSVC` like its siblings have. Needs Highway's *own* test
    suite passing under MSVC with AVX-512; two error-function kernels are
    not sufficient evidence for a PR.
- [RESOLVED 2026-07-24] `CMakeLists.txt` hard-errored under every
  multi-config generator: `set_property(CACHE CMAKE_BUILD_TYPE ...)` ran
  unconditionally, but the guard above it correctly skips creating that
  cache entry when `CMAKE_CONFIGURATION_TYPES` is set, so the property
  call referenced a nonexistent variable. `-G "Visual Studio 18 2026"`
  failed at configure; Ninja and the presets were unaffected, which is why
  it went unnoticed. Both statements now sit inside one
  `if(NOT CMAKE_CONFIGURATION_TYPES)`. Verified: VS generator configures
  clean, Ninja unchanged.
- [RESOLVED 2026-07-20] Gather performance on x86 (Kaby Lake, Release,
  session-loaded machine — treat as indicative) [DERIVED]: erf table-gather
  kernel beats scalar libm 3.5-4.8x on ALL tiers (AVX2/SSE4/SSSE3/SSE2,
  ~6-8.5 ns/el vs ~26-30), so NOT a repeat of libstats #33's null — scalar
  erf is expensive enough that SIMD wins regardless. BUT zero width
  scaling: AVX2 4-lane == SSE2 2-lane ns/el. Suspects: native AVX2 gather
  throughput (emulated SSE gathers are cheap scalar loads) and emulated
  f64->i64 ConvertTo below AVX-512DQ. Ship as-is; a non-gather x86 variant
  is a known ~2x AVX2 upside if ever needed.
- [RESOLVED 2026-07-24] Ryzen re-benchmark — **the "flat regardless of
  width" pattern does not hold on Zen 4; width scaling appears.** Ryzen
  7445HS, GCC 16.1, Release, machine had been running builds — indicative,
  not quiet-machine numbers [DERIVED]. erf ns/el at n=1e6: SSE2 (2-lane)
  5.49 → AVX2 (4-lane) 2.66 → AVX3_ZEN4 (8-lane) 1.95, i.e. ~2.1x then a
  further ~1.4x, vs libm 1.95–2.80 ns/el for 3.4–6.2x speedup. erfc at
  n=1e6 on AVX3_ZEN4: 9.05x core-dominated, 5.12x mixed, 8.38x tail-only
  — well above the Kaby Lake AVX2 figures recorded above. This supports
  the stated suspects: the Kaby Lake flatness was that CPU's gather
  throughput plus emulated f64->i64 below AVX-512DQ, not an inherent
  kernel limit. The non-gather x86 variant therefore remains a real but
  lower-priority upside, and mainly for pre-AVX-512 hardware. Worth
  re-running on a genuinely quiet machine before any published number.
- [OPEN] `CORVUS_SANITIZE` is not MSVC-aware: it emits `-fsanitize=<list>`
  unconditionally, which cl.exe does not accept (MSVC wants
  `/fsanitize=address`, and has no UBSan). Harmless today because sanitizer
  builds only run on Linux, but the option silently produces a broken
  command line if anyone tries it on Windows. Either branch on MSVC or
  reject the combination with a clear `message(FATAL_ERROR)`.
- [OPEN] bench_erf harness (tests/bench_erf.cpp, not ctest-registered) is
  the per-kernel benchmark pattern — reuse for erfc/lgamma.
- [OPEN] Install/export when Highway is FetchContent-built: currently
  disabled (exported target would dangle). Options: require system hwy for
  install (status quo), bundle hwy objects into libcorvus.a, or install a
  nested hwy. Decide before first tagged release.
- [OPEN] Generalize the ULP harness (gen_erf_reference.py + test_erf_ulp)
  into a per-kernel pattern as new functions land; reference files are
  checked in, generators need mpmath in a throwaway venv.
- [OPEN] Pre-release legal: binary artifacts that link Highway must carry
  its Apache-2.0 NOTICE; source-only distribution needs nothing. Handle
  when packaging/releases start.
- [RESOLVED 2026-07-21] CI: .github/workflows/ci.yml, designed around
  runner-minute economy (user lesson from libstats' private phase): one
  Linux job sweeps AVX2..SSE2 sequentially + ASan/UBSan in-job; one macOS
  arm64 job for native NEON. Windows/MSVC was absent only because the
  project had never run on Windows — a `windows-msvc` toolchain job was
  added 2026-07-24 (see the resolved item above); AVX-512 remains
  impossible on hosted runners (Ryzen stays manual);
  required status checks dropped (blocks direct-push workflow); no caching
  until minutes justify it. First green NEON run should update
  docs/ACCURACY.md's NEON column (gates may trip if hn::Exp accuracy
  differs on NEON — that would be a finding, not a nuisance).
- [OPEN] Decide whether libstats/libhmm adopt corvus as a dependency or
  keep their internal SIMD (migration is a separate project-level decision).
- [ILLUSTRATIVE] Possible future consumer: C++ port of multi-agent_sim
  (batch distance/trig), zeekhmm training pipelines.

## Build-Stack Standardization (2026-07-23) [DERIVED]
Cross-repo effort tracked in `~/Development/BUILD-STANDARDIZATION-PLAN.md`.
Commits: `3bbecf1` (pkg-config file, consumer example, installed-path CI
check), `c158765` (`find_package(hwy 1.4)` version floor + CI builds pinned
Highway 1.4.0 from source instead of apt, fixing a distro-libhwy-dev/
hn::ReduceMax mismatch found by the first Phase-1 CI run), `1c220b3`
(CMakePresets.json: release/debug/rel-with-debug/sanitize). AGENTS.md's
CMake-standard section checked post-Phase-3 and is still accurate (already
documents presets, the Highway find_package/FetchContent split, and the
install-when-system-Highway gate). The fetched-Highway install gate itself
stays tracked under "Open Items" above, not duplicated here.

## Next Steps
1. [RESOLVED 2026-07-24] Native AVX-512 validation on the Ryzen, per-tier
   recipe un-capped then capped down through AVX2/SSE4/SSSE3/SSE2; docs
   updated. See the resolved open item for results. (The reminder that
   held up: HWY_NATIVE_FMA follows the compiled HWY_TARGET, not the
   physical CPU, so Ryzen's capped SSE runs exercised ops::SquareLow's
   Dekker fallback exactly as Kaby Lake's do — the new information from
   Ryzen was AVX3\* itself, plus the first measurement of how far the
   no-FMA path diverges.)
2. [RESOLVED 2026-07-24] bench_erf / bench_erfc on the Ryzen — width
   scaling does appear on Zen 4; see the resolved re-benchmark item.
   Still worth one re-run on a genuinely quiet machine before publishing.
3. [RESOLVED 2026-07-21, M1 session] Sequencing decided with the user:
   dd transcendental core first (Phase A), then lgamma (Phase B), full
   real axis in v1 — design section above.
4. **Start here**: Phase A — exp_dd + log_dd + ops::BitCast + generators;
   erfc tail rewire and gate retightening as acceptance (ACCURACY.md in
   the same change set). New-kernel AVX-512 validation joins the Ryzen
   queue. No accuracy or tier work is outstanding for erf/erfc.
5. Phase B: lgamma per the design section (generator experiments pick
   X0, degrees, zone boundaries; new gen_lgamma_reference.py with
   zero-neighborhood, pole-neighborhood, and boundary-crossing points).
