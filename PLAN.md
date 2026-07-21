# corvus — Plan / Session State

## Status [DERIVED]
Scaffolded 2026-07-20 on the Kaby Lake Mac (AVX2); pushed to
github.com/OldCrow/corvus (public) the same day. Erf is production-quality:
clean-room table + local-Taylor kernel ported from libstats
vector_erf_neon through the ops facade (gather-based, portable NaN
handling), **max 1 ULP validated vs mpmath oracle on AVX2 and tier-capped
SSE4** (14,594-point reference incl. grid-residual and saturation-boundary
clusters; ~98.5% correctly rounded). NEON (M1) and AVX-512 (Ryzen)
validation pending — accuracy is claimed only for AVX2/SSE4 until then.
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
- [OPEN] Validate erf on NEON (M1) and native AVX-512 (Ryzen 7445) before
  claiming those tiers; also run the SSE2-only cap (current lowest tested
  cap dispatches SSE4).
- [OPEN] Gather performance on x86 is unprofiled — libstats issue #33 found
  gather-based transcendentals a null win on AVX2/AVX-512 for exp/log.
  Correctness is validated; benchmark erf vs scalar before performance
  claims (may want a non-gather x86 variant if it disappoints).
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
