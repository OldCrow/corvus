# corvus — Plan / Session State

## Status [DERIVED]
Scaffolded 2026-07-20 on the Kaby Lake Mac (AVX2); pushed to
github.com/OldCrow/corvus (public) the same day. Erf is production-quality:
clean-room table + local-Taylor kernel ported from libstats
vector_erf_neon through the ops facade (gather-based, portable NaN
handling), **max 1 ULP validated vs mpmath oracle on AVX2 and tier-capped
SSE4/SSSE3/SSE2** (14,594-point reference incl. grid-residual and
saturation-boundary clusters; ~98.5% correctly rounded) — every x86 tier
Kaby Lake can express natively, via CORVUS_DISABLED_TARGETS capping.
Note Highway has no AVX1 target (AVX1-class hardware dispatches SSE4).
NEON (M1) and the AVX3* variants (Ryzen 7445) remain hardware-bound
pending — accuracy is claimed only for the four validated tiers until then.
Facade gained table-kernel ops: SignedTag, Round, ConvertToInt, ShiftLeft,
GatherIndex, Ge, IsNaN.

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

## Open Items
- [OPEN] erfc tail 5-ULP bound is entirely hn::Exp's contribution. A
  corvus-owned compensated exp (double-double reduction) would tighten the
  tail toward 1-2 ULP and also serves future kernels (lgamma, incomplete
  gamma). Consider before or alongside lgamma.
- [OPEN] Validate erf on NEON (M1) and native AVX3* (Ryzen 7445) before
  claiming those tiers — the only remaining gaps; all four x86 tiers
  expressible on Kaby Lake (AVX2/SSE4/SSSE3/SSE2) passed max 1 ULP
  2026-07-20.
- [RESOLVED 2026-07-20] Gather performance on x86 (Kaby Lake, Release,
  session-loaded machine — treat as indicative) [DERIVED]: erf table-gather
  kernel beats scalar libm 3.5-4.8x on ALL tiers (AVX2/SSE4/SSSE3/SSE2,
  ~6-8.5 ns/el vs ~26-30), so NOT a repeat of libstats #33's null — scalar
  erf is expensive enough that SIMD wins regardless. BUT zero width
  scaling: AVX2 4-lane == SSE2 2-lane ns/el. Suspects: native AVX2 gather
  throughput (emulated SSE gathers are cheap scalar loads) and emulated
  f64->i64 ConvertTo below AVX-512DQ. Ship as-is; a non-gather x86 variant
  is a known ~2x AVX2 upside if ever needed. Re-benchmark on Ryzen
  (AVX-512 has native f64->i64 and different gather hardware) and M1.
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
- [OPEN] CI: GitHub Actions matrix (macOS arm64/x86_64, Linux x86_64) with
  CORVUS_DISABLED_TARGETS jobs to natively exercise SSE2..AVX-512 tiers on
  the AVX-512 runner-or-self-hosted question TBD.
- [OPEN] Decide whether libstats/libhmm adopt corvus as a dependency or
  keep their internal SIMD (migration is a separate project-level decision).
- [ILLUSTRATIVE] Possible future consumer: C++ port of multi-agent_sim
  (batch distance/trig), zeekhmm training pipelines.

## Next Steps
1. Validate erf AND erfc on M1 (NEON) and Ryzen (native AVX3* tiers) —
   note NEON has native FMA so SquareLow takes the MulSub path there.
2. Decide: corvus-owned compensated Exp before lgamma, or lgamma first
   (see erfc tail open item).
3. lgamma — first P0 function without libstats prior art.
