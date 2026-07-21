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
- Platform tiers: Tier 1 (accuracy-audited on real silicon) = NEON (M1),
  AVX-512/AVX2/SSE2 (Ryzen 7445 native + CORVUS_DISABLED_TARGETS capping),
  AVX2 (Kaby Lake). Tier 2 (compiles, unaudited) = SVE and anything else
  Highway emits.

## Open Items
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
- [OPEN] CI: GitHub Actions matrix (macOS arm64/x86_64, Linux x86_64) with
  CORVUS_DISABLED_TARGETS jobs to natively exercise SSE2..AVX-512 tiers on
  the AVX-512 runner-or-self-hosted question TBD.
- [OPEN] Decide whether libstats/libhmm adopt corvus as a dependency or
  keep their internal SIMD (migration is a separate project-level decision).
- [ILLUSTRATIVE] Possible future consumer: C++ port of multi-agent_sim
  (batch distance/trig), zeekhmm training pipelines.

## Next Steps
1. Add erfc (shares the table infrastructure; needs its own compensated
   tail handling — erfc(x) for large x underflows differently).
2. Validate on M1 (NEON) and Ryzen (native AVX-512 + SSE2 cap).
3. Benchmark erf vs scalar std::erf per tier (see gather open item).
4. Then lgamma — first P0 function without libstats prior art to port.
