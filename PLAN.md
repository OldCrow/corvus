# corvus — Plan / Session State

## Status [DERIVED]
Scaffolded 2026-07-20 on the Kaby Lake Mac (AVX2). Repo proves the full
pipeline: op facade -> foreach_target kernel -> runtime dispatch -> masked
tail -> ctest. One function (Erf) implemented with a PROVISIONAL
approximation (A&S 7.1.26, ~1.5e-7 abs error) purely to validate plumbing.

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
- [OPEN] Erf is PROVISIONAL: replace A&S 7.1.26 with a clean-room max-ULP
  implementation (libstats' NEON clean-room erf is prior art — port the
  algorithm through the ops facade, verify provenance notes carry over).
- [OPEN] Install/export when Highway is FetchContent-built: currently
  disabled (exported target would dangle). Options: require system hwy for
  install (status quo), bundle hwy objects into libcorvus.a, or install a
  nested hwy. Decide before first tagged release.
- [OPEN] Accuracy harness: per-kernel, per-tier ULP measurement (mpmath
  reference like libstats' table generators) rather than ad-hoc abs-error
  checks in tests.
- [OPEN] CI: GitHub Actions matrix (macOS arm64/x86_64, Linux x86_64) with
  CORVUS_DISABLED_TARGETS jobs to natively exercise SSE2..AVX-512 tiers on
  the AVX-512 runner-or-self-hosted question TBD.
- [OPEN] Decide whether libstats/libhmm adopt corvus as a dependency or
  keep their internal SIMD (migration is a separate project-level decision).
- [ILLUSTRATIVE] Possible future consumer: C++ port of multi-agent_sim
  (batch distance/trig), zeekhmm training pipelines.

## Next Steps
1. Verify first build + ctest pass on this machine (AVX2) — in progress at
   scaffold time.
2. git init, initial commit, create OldCrow/corvus GitHub repo (user
   approval before push, SSH remote per global rules).
3. Replace provisional erf; add erfc alongside (shared kernel structure).
4. Build the accuracy harness before adding more functions — it gates all
   accuracy claims.
