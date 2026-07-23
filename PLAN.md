# corvus — Plan / Session State

## Status [DERIVED] — end of session 2026-07-21
Scaffolded 2026-07-20 on the Kaby Lake Mac (AVX2); public at
github.com/OldCrow/corvus. Two functions shipped, both production-quality:

- **erf**: clean-room table + local-Taylor kernel (ported from libstats
  vector_erf_neon through the ops facade). Max 1 ULP over the full domain.
- **erfc**: two-region kernel (core reuses the erf table via compensated
  assembly; tail is a fitted e^{-a^2}*G(1/a)/a). Max 1 ULP for |x| <= 6
  and subnormal results; max 5 ULP normal-tail (bounded by backend Exp).

**Validated tiers, both functions, all identical**: AVX2, SSE4, SSSE3,
SSE2 (Kaby Lake native + CORVUS_DISABLED_TARGETS capping) and **NEON**
(Apple Silicon GitHub Actions runner, native FMA — validated in CI
2026-07-21). NEON and AVX2 produce bit-identical ULP results point-for-
point on both reference sets (see docs/ACCURACY.md) — first evidence the
kernels are deterministic across ISA/compiler/OS on FMA-capable targets.

**Only remaining accuracy gap: AVX-512 (the four AVX3* variants) on the
Ryzen 7445** — structurally impossible to close via hosted CI runners;
this is the one gap that needs the physical machine. Tomorrow's
Ryzen session is expected to close it via the standard tier-sweep recipe
(AGENTS.md Workflows) run natively rather than capped, plus a
CORVUS_DISABLED_TARGETS sweep down through AVX2/SSE tiers on that box too
(second independent-hardware data point for those, on top of Kaby Lake +
Linux CI).

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
  both erf and erfc, bit-identical to AVX2 results. **[OPEN, Ryzen-bound]**
  Native AVX3*/AVX-512 validation is the only remaining tier gap for both
  functions — do this first on the Ryzen tomorrow (AGENTS.md per-tier
  recipe, run un-capped for AVX-512 itself, then capped down through
  AVX2/SSE tiers as a second hardware data point). Update docs/ACCURACY.md
  validation matrix and the cross-arch-reproducibility note with the
  result either way (a divergence from the NEON/AVX2 match would itself
  be a finding, not just a checkbox).
- [RESOLVED 2026-07-20] Gather performance on x86 (Kaby Lake, Release,
  session-loaded machine — treat as indicative) [DERIVED]: erf table-gather
  kernel beats scalar libm 3.5-4.8x on ALL tiers (AVX2/SSE4/SSSE3/SSE2,
  ~6-8.5 ns/el vs ~26-30), so NOT a repeat of libstats #33's null — scalar
  erf is expensive enough that SIMD wins regardless. BUT zero width
  scaling: AVX2 4-lane == SSE2 2-lane ns/el. Suspects: native AVX2 gather
  throughput (emulated SSE gathers are cheap scalar loads) and emulated
  f64->i64 ConvertTo below AVX-512DQ. Ship as-is; a non-gather x86 variant
  is a known ~2x AVX2 upside if ever needed. [OPEN, Ryzen-bound]
  Re-benchmark bench_erf/bench_erfc on the Ryzen: AVX-512 has native
  f64->i64 (AVX-512DQ) and different gather hardware, so this is the tier
  most likely to break the "flat regardless of width" pattern — worth
  running even though correctness validation is the primary reason to be
  on that machine.
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
  arm64 job for native NEON. Windows/MSVC deferred until MSVC support
  exists; AVX-512 impossible on hosted runners (Ryzen stays manual);
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

## Next Steps (tomorrow, expected on the Ryzen)
1. **Start here**: native AVX-512 validation via the AGENTS.md per-tier
   recipe, un-capped first (native AVX3* dispatch), then capped down
   through AVX2/SSE4/SSSE3/SSE2 for a second independent-hardware
   confirmation of the NEON/Kaby-Lake match. Update docs/ACCURACY.md
   matrix and the reproducibility note. (Reminder: HWY_NATIVE_FMA follows
   the compiled HWY_TARGET, not the physical CPU — Ryzen's capped SSE runs
   exercise ops::SquareLow's Dekker fallback exactly like Kaby Lake's do,
   not a new code path; the new information from Ryzen is AVX3* itself.)
2. bench_erf / bench_erfc on the Ryzen (see gather-performance open item
   above) — AVX-512 is the tier most likely to change the "flat regardless
   of width" finding from Kaby Lake.
3. [RESOLVED 2026-07-21, M1 session] Sequencing decided with the user:
   dd transcendental core first (Phase A), then lgamma (Phase B), full
   real axis in v1 — design section above.
4. Phase A: exp_dd + log_dd + ops::BitCast + generators; erfc tail
   rewire and gate retightening as acceptance (ACCURACY.md in the same
   change set). New-kernel AVX-512 validation joins the Ryzen queue.
5. Phase B: lgamma per the design section (generator experiments pick
   X0, degrees, zone boundaries; new gen_lgamma_reference.py with
   zero-neighborhood, pole-neighborhood, and boundary-crossing points).
