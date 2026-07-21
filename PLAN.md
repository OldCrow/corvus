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
- [OPEN] erfc tail 5-ULP bound is entirely hn::Exp's contribution. A
  corvus-owned compensated exp (double-double reduction) would tighten the
  tail toward 1-2 ULP and also serves future kernels (lgamma, incomplete
  gamma). Consider before or alongside lgamma.
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
3. Decide: corvus-owned compensated Exp before lgamma, or lgamma first
   (see erfc tail 5-ULP open item) — this is a [[frontload-project-
   conventions]]-flavored call worth making deliberately rather than by
   default-to-next-item.
4. lgamma — first P0 function without libstats prior art to port from.
