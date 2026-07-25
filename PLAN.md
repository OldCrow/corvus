# corvus — Plan / Session State

## Status [DERIVED] — end of session 2026-07-25 (Ryzen)
**Phases A and B are complete.** Scaffolded 2026-07-20 on the Kaby Lake Mac
(AVX2); public at github.com/OldCrow/corvus. Three public functions plus
both internal transcendental cores, all production-quality:

- **lgamma**: full real axis, five regions, everything in dd with a single
  final rounding. Max 1 ULP across the positive axis including the exact
  zeros, correctly rounded throughout Stirling, 1 ULP / 2^-53 absolute on
  the negative axis. Identical on all five x86 tiers and under GCC and
  MSVC. See the Phase B section.

- **erf**: clean-room table + local-Taylor kernel (ported from libstats
  vector_erf_neon through the ops facade). Max 1 ULP over the full domain.
- **erfc**: two-region kernel (core reuses the erf table via compensated
  assembly; tail is a fitted e^{-a^2}*G(1/a)/a on exp_dd). Max 1 ULP for
  |x| <= 6 and subnormal results; max 2 ULP normal-tail (was 5 before the
  dd rewire; the residual is the G fit, decomposed in the open items).
- **exp_dd** (internal): 2^-68.45 relative. **log_dd** (internal):
  2^-67.88, correctly rounded on every reference point. No corvus kernel
  depends on Highway contrib math any more.

FP contraction is disabled project-wide (`CORVUS_FP_FLAGS`) — the dd layer's
exactness identities assume IEEE ops as written, and GCC's default was
silently shifting log_dd by 0.6 bits relative to MSVC. See the decision
section below.

**Next session starts at Phase C** — see Next Steps. The candidates are
erfinv/erfcinv and the regularized incomplete gamma/beta pair; no design
work has been done on either, so that is the first task.

**Validated tiers, all four kernels**: AVX3_ZEN4/AVX3_DL/AVX3, AVX2, SSE4,
SSSE3, SSE2 (Ryzen native + CORVUS_DISABLED_TARGETS capping, GCC 16.1;
erf/erfc also on Kaby Lake and Linux CI) and **NEON** (Apple Silicon CI
runner). erf/erfc are bit-identical across the FMA-capable targets and
differ marginally on the no-FMA SSE tiers in not-CR counts only; **exp_dd
and log_dd are identical on every tier and under all three compilers**
(GCC, MSVC, AppleClang) — same bound, same worst-case input. Counts per
tier are in docs/ACCURACY.md.

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
- Windows compiler (2026-07-25, resolved with user): **corvus takes a
  documented exception to the fleet standard.** The other projects default
  to MSVC on Windows (paired with Apple Clang on macOS, Clang on Linux);
  corvus prefers `clang-cl` there. Reason: Highway's `HWY_BROKEN_MSVC`
  removes every AVX-512 target, so MSVC cannot validate or benchmark the
  widest tier — which is the entire reason the Ryzen box is in the fleet.
  `clang-cl` is the preferred alternative rather than mingw GCC because it
  keeps the MSVC ABI and therefore stays link-compatible with the rest of
  the fleet, which matters for the open question of libstats/libhmm adopting
  corvus. Scope of the exception is narrow: it governs which compiler
  produces validation and benchmark numbers, not which compilers must work.
  MSVC stays a fully supported consumer toolchain and remains the CI
  Windows job's compiler, because it is both the strictest diagnostic gate
  and the likeliest consumer default. [OPEN] If a sibling project does adopt
  corvus, verify clang-cl-built corvus links cleanly into an MSVC-built
  consumer before relying on it — same-ABI is the design intent but has not
  been tested here.
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

## Phase A part 1 — exp_dd shipped [DERIVED, 2026-07-25 Ryzen]
`src/dd-inl.h` (dd primitives), `src/exp_dd-inl.h` + generated tables, and
the erfc tail rewired onto them. Acceptance met: **erfc tail-normal 5 -> 2
ULP**, gate retightened, same reference set, no other bound moved.

- exp_dd measures **2^-68.45 relative** (~2^-15 ULP) — 8 bits better than
  the 2^-60 design target — and is **bit-identical on all five validated
  tiers including the no-FMA ones**, which erf/erfc are not. The Dekker
  fallback in the new `ops::ProdLow` reproduces the fused path exactly
  here; same worst-case input (x = -696.994) everywhere.
- Facade grew `BitCast`, `And`, `ShiftRight`, `ProdLow` (the two-operand
  `SquareLow`, same FMA guard), and `Set` generalized to the tag's lane
  type so it works on the signed tag. `ops::Exp` is now unused by any
  shipped kernel — the contrib exp is out of the accuracy-critical path,
  as intended.
- Design point worth keeping: exp_dd returns **mantissa + exponent, not a
  scaled value** (`ExpDdFrac` + `ScaleTwo`). The consumer folds its own
  factors in first and the power-of-two scaling lands last, which is what
  keeps erfc's subnormal band at one rounding. `ExpDd` (fully scaled) is
  the convenience form; its dd invariant necessarily degrades once the
  result is subnormal.
- Cost, Ryzen/GCC/AVX3_ZEN4, n=1e6, session-loaded so indicative: tail-only
  8.38x -> **4.83x** vs libm (~1.7x slower tail path), mixed 5.12x ->
  3.60x. Core untouched (9.05x -> 9.80x, noise). This is the trade the
  design section predicted: two gathers plus dd math against one contrib
  Exp call. Re-measure on a quiet machine before publishing.
- Testing pattern established for internal kernels: `tests/test_exp_dd.cpp`
  compiles the kernel header itself through foreach_target instead of
  linking the public API, so it needs `hwy::hwy`, the source root, and
  **the same target defines as the library** — hoisted into
  `CORVUS_HWY_TARGET_DEFS` for exactly that reason. It also asserts its own
  dispatched target equals `corvus::active_target()`, so the two dispatch
  tables cannot silently disagree.
- Oracle pattern for dd kernels: the reference file carries the value as a
  **dd pair**, and the test measures relative error below the last bit of a
  double. Rounding first would hide what the kernel is for.

## Phase A part 2 — log_dd shipped [DERIVED, 2026-07-25 Ryzen]
`src/log_dd-inl.h` + generated tables + `tests/test_log_dd.cpp`, same
pipeline as exp_dd. No in-tree consumer until lgamma, so its own gate is the
acceptance test. **2^-68.48 relative on FMA tiers, 2^-67.88 on no-FMA;
correctly rounded on every one of the 7273 reference points, every tier.**

- Two design choices carry it, both worth not re-deriving: the mantissa is
  **centred on 1** (halving at the slot boundary 1+53/128, so the condition
  is exactly `j >= 53` and no slot straddles it), and the **two slots
  adjacent to m = 1 carry R_j = 1, L_j = 0 exactly**. Together they mean
  there is no special case for x near 1 and no cancellation on either side
  of it. Uncentred, x just below 1 would be ln2 + log(m) — ~45 bits of a
  dd's 106 burned exactly where relative accuracy matters most.
- `r = R_j*m - 1` is **exact**: p = fl(R_j*m) is in [1-2^-7, 1+2^-7] so
  p - 1 is exact by Sterbenz, and ops::ProdLow gives the product residual.
  A plain fma(R_j, m, -1) would carry 2^-62 (2^-53 without FMA).
- **The bug this cost, worth remembering: an exact dd pair is not
  necessarily a NORMALIZED one.** (p-1, p_lo) has |p_lo| up to 2^-53 against
  a hi of ~2^-8. Every term computed from `r.hi` alone then drops its share
  of `r.lo`; the r^3 term's share is r^2*r_lo ~ 2^-69, which showed up as
  2^-63.3 relative — 5 bits below target, and only near x = 1. Fix is one
  TwoSum (not Fast2Sum: r passes through zero inside every slot). The gate
  caught it; reading the code did not. Same shape as the MulSub/no-FMA
  hazard: a primitive that is *correct* used in a way that is *unsound*.
- Facade grew `Or` and `ConvertToDouble` (int->double). ops:: is now ~30 ops.

### FP contraction turned off project-wide [DECISION, 2026-07-25]
Found while cross-checking log_dd under MSVC: **GCC 2^-68.48 vs MSVC
2^-67.88 on identical source at the same tier.** Cause is GCC's default
`-ffp-contract=fast` fusing a Mul into an adjacent Add in the log1p path.
Confirmed by rebuilding GCC with `-ffp-contract=off` and landing exactly on
MSVC's number. Contraction was making the result *better* — which is why
nothing failed and why it would never have been noticed from one compiler.

Now `CORVUS_FP_FLAGS` sets `-ffp-contract=off` (clang-cl:
`/clang:-ffp-contract=off`; MSVC needs nothing, `/fp:precise` already
doesn't contract), PRIVATE to corvus and to the kernel test targets — a
test that compiles kernel code must measure the program that ships.

Rationale, and why this is not merely tuning: the dd layer is a set of
*exact* identities whose proofs assume each IEEE op is rounded as written.
A contraction landing inside a TwoSum makes `s` something other than
`fl(a+b)`, so the "exact" residual is the error of an operation that never
happened. corvus already claims cross-compiler reproducibility in
ACCURACY.md; leaving this to the optimizer contradicts that claim.
**erf and erfc are bit-for-bit unaffected** — every fusion they want is
already an explicit `ops::MulAdd`, which is the principle: fusion is
requested in the source, never inferred.

Throughput cost, same loaded machine, so indicative: erfc core 2.01 ->
2.11 ns/el, tail 4.50 -> 4.86, erf 1.95 -> 2.16. Up to ~8%, which is at
the edge of this machine's run-to-run spread — treat as "no larger than
8%" and re-measure on a quiet machine before publishing either number.

[OPEN] Worth a look at whether the sibling projects (libstats/libhmm) have
the same latent exposure in their compensated-summation paths — same
compilers, same class of algorithm, and the symptom is silence.

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

## Phase B — lgamma shipped [DERIVED, 2026-07-25 Ryzen]
`src/lgamma-inl.h` + `src/lgamma.cpp` + generated `src/lgamma_data.h`,
`tools/gen_lgamma_data.py`, `tools/gen_lgamma_reference.py`,
`tests/test_lgamma{,_ulp}.cpp`, `tests/bench_lgamma.cpp`. Bounds and method
are in docs/ACCURACY.md; what follows is only what would be expensive to
re-derive.

**The design section below was architecturally complete and was followed.**
Regions, the exact-t form at the zeros, the fixed-step masked recurrence,
Stirling with a swept X0, and the sinpi reflection all landed as written.
Four things it did not settle, resolved here:

1. **The recurrence direction is forced, and it is DOWN.** x − k is exactly
   representable for every integer k < x < 2^52 (x is a multiple of
   ulp(x) ≤ 1, and |x − k| < |x| needs no finer grid); x + k is a multiple
   of the same grid but larger, so it can need a bit the format does not
   have — x + 8 for x ∈ (2.5, 8) loses up to two. Walking up would put
   ~77 ULP into lgamma(2.5) unless every shifted argument were carried as a
   dd pair. Walking down also lands every lane in (1.5, 2.5], so one
   polynomial serves the whole region. This is why the design's
   "P ≤ ~Γ(X0)" bound is exactly right: 7·6·5·4·3·2 = Γ(8).
2. **The reflection uses Γ(1−x) = −x·Γ(−x), not Γ(x)Γ(1−x) directly.**
   1 − x is inexact and would have forced a dd argument through every
   region; −x is exact. Cost is one extra log. The payoff beyond exactness
   is structural: the positive pipeline then runs ONCE, on |x|, for lanes
   of either sign, and the reflection is a post-correction — no second
   evaluation, no lane divergence.
3. **π/|sin πx| is never formed.** It overflows at every pole. Writing
   sin(πu) = πu·(sin(πu)/(πu)) cancels π analytically and leaves
   −log|u| — the term that *should* diverge — plus a bounded correction
   fitted in u².
4. **lgamma(1) needs an explicit +0.0.** t·B(t) at t = 0 carries the sign
   of B, which is −γ at x = 1, so the natural result is −0 where C99 wants
   +0. Adding +0.0 is the identity on every finite nonzero value and costs
   one op; doing it that way rather than a select also removes any
   dependence on how a given target's Fast2Sum signs a zero.

**The bug worth remembering: `ops::ProdLow`'s non-FMA Dekker split
overflows for |a| > 2^996.** Its intermediate is a·(2^27+1). lgamma's
Stirling product x·(log x − 1) is the first corvus operand ever to reach
that range, so for x > ~6.7e299 the residual came back infinite on SSE2
through SSE4 while every FMA target was correct — 500 of 4242 points wrong,
max ULP ~8.6e16. Fixed by scaling x down 2^200 for the product and back up
after: exact, linear, one code path for all targets. **Only the capped
sweep could find this** — native AVX-512, AVX2, GCC and MSVC were all
green. Same family as the MulSub zero-residual hazard: a primitive that is
correct, used outside its stated precondition. The precondition is now in
AGENTS.md Conventions rather than only in the facade comment.

Two smaller notes:
- The Stirling assembly is grouped x·(log x − 1) − log(x)/2 + …, not the
  textbook (x − ½)·log x − x. Algebraically identical, but the textbook
  product exceeds the result by x, which overflows to inf — and its dd
  residual to NaN — over a ~0.1%-wide band of x below the true overflow
  threshold. Grouped this way the excess is ~log(x)/2 ≈ 351 against an ulp
  of 2e292, so product and result overflow at the same double.
- `mpmath` mis-sums a generator passed to `fsum`, and an ODD Chebyshev node
  count puts a node at t = 0 where `1 + t` rounds to 1 and B(t) collapses to
  0/0. Either one alone makes the coefficients plateau at ~0.007 instead of
  decaying, which reads as "this function is not smooth" rather than as a
  tooling bug. Both are commented in the generator.

Parameters chosen by the generator sweep, as the design intended: X0 = 8
(ψ degree 8; X0 = 13 would cut it to 6 but adds five recurrence steps),
zone split at 3/2, outer boundary at 1/2, degrees 34 and 21 with three dd
leading coefficients each.

Cost, Ryzen/GCC/AVX3_ZEN4, n=1e6, session-loaded so indicative: zone 2.60x,
recurrence 3.23x, Stirling 2.97x, mixed positive 1.89x, reflection 3.02x vs
libm. The zone is the slowest region — a degree-34 Horner with a per-lane
coefficient select is most of it. Re-measure on a quiet machine before
publishing.

## Phase C — broad design and build order [DECISION, 2026-07-25]
Order: **erfinv/erfcinv → regularized incomplete gamma P/Q → regularized
incomplete beta.** Why:
1. The inverse pair is the smallest candidate and unlocks the single
   most-demanded statistical inverse — the normal quantile
   (probit(p) = −√2·erfcinv(2p)) — immediately usable by libstats.
2. It settles the erfc-tail open item with numbers instead of a fear.
   The relative condition of x = erfcinv(s) in s is
   κ = s/(x·|erfc′(x)|) ≈ 1/(2x²) — the inverse of erfc's tail
   ill-conditioning — so a Newton step against erfc ATTENUATES the
   tail's 2-ULP error by ≥ 72× (x ≥ 6). The compounding the open item
   anticipated does not occur; expect the re-fit question to close as
   "not needed" once measured. The central region is the opposite:
   κ ≈ 1.0–1.2, a Newton against 1-ULP erf passes the error straight
   through, so the centre goes direct-polynomial. Same analysis, both
   directions, and it dictates the region structure below.
3. It establishes the seed + dd-Newton inverse pattern that the P1
   inverse incomplete gamma/beta will reuse, on its cheapest instance.
4. Gamma before beta: beta's hardest machinery (uniform asymptotics,
   Stirling-difference prefactor) is a superset of gamma's; build it
   once with two arguments before three.

### erfinv/erfcinv — broad parameters
Both public functions are thin routers over two shared cores, every
routed argument exact by Sterbenz:
- erfinv: |y| ≤ 1/2 → C(y); 1/2 < |y| < 1 → sign(y)·T(1 − |y|)
  (1 − |y| exact for |y| ∈ [1/2, 1]).
- erfcinv: z ∈ [1/2, 3/2] → C(1 − z) (1 − z exact on [1/2, 3/2]);
  z < 1/2 → T(z); z > 3/2 → −T(2 − z) (2 − z exact on [1, 2]).
  The zero at z = 1 gets relative accuracy for free, by the same
  mechanism as lgamma's zeros: an exact argument through an odd form.
- Core C (central, |x| ≤ 0.4769): x = y·Pc(y²), dd leading
  coefficient(s), relative target ~2^-55 — the lgamma-zone pattern.
  Direct fit, no Newton (κ note above). Modest degree expected (~16 by
  Bernstein estimate: nearest singularity y² = 1 vs interval [0, 1/4]).
  erfinv(±0) = ±0 falls out of the odd form.
- Core T (tail, s ∈ (0, 1/2) → x ∈ (0.4769, 27.217)): primary variable
  w = −log s via log_dd (LogDdAny covers subnormal s), t = √w.
  Far tail x ≥ ~6: seed fit in 1/t + one dd Newton/Halley step in log
  space — F(x) = log erfc(x) + w with
  log erfc(x) = −x² − log x + log G(1/x), reusing the erfc tail
  structure; x² exact via ProdLow (x ≤ 28, no 2^996 hazard). Halley is
  nearly free here (F″ = −2x·F′ − F′², no new transcendentals), which
  relaxes the seed to ~2^-19 for one cubic step.
  Mid region (0.4769, 6): either direct per-interval dd fits in t, or
  the same step against an internal dd erfc (the erfc core's
  compensated 1 ∓ erf assembly already carries dd internally — needs
  an ErfcDd exposure, no public change). Decide by generator
  experiment + bench; this is the main open design choice.
- Specials: erfinv(±1) = ±inf, |y| > 1 → NaN; erfcinv(0) = +inf,
  erfcinv(2) = −inf, z outside [0, 2] → NaN; NaN propagates.
- Target: ≤ 2 ULP everywhere, 1 ULP goal in C and the far tail.

### Phase C part 1 — erfinv/erfcinv design [resolved 2026-07-25]
Probes (mpmath, scratchpad `probe_erfinv.py`): central needs degree
~13–14 in y² for 2^-55; 2^-19 seeds need degree 6/8/4; direct
full-accuracy tail fits would need degree 41+ single-interval (20–29
split) AND are floored near 1 ULP anyway because t = √(−log s) arrives
rounded — unlike lgamma's exact t. That kills direct fits for the tail:
**seed + one dd Halley step everywhere in T.**

- **C** (|x| ≤ 0.4769): x = y·Pc(y²), degree ~14, dd leading
  coefficient(s) (count by generator replay). No log_dd on this path.
  erfinv(±0) = ±0 free; subnormal y fine (y² underflows, Pc → c0).
  erfcinv(1) = +0 falls out via C(1 − 1).
- **T** (s ∈ (0, ½)): w = −log s via LogDdAny (subnormal s reachable
  from erfcinv only), t = √(w.hi); seed x0 from three interval polys
  (t on [t_lo, ~2] and [~2, 6.195]; far in u = 1/t — the erfc-tail
  coefficient-select pattern), then one Halley step:
  - **far** (routed by t ≥ kTfar ⇔ x ≥ ~6): log space.
    F = (w ⊖ x0²)_dd − log(x0/G(1/x0)); x0² exact via ProdLow (x ≤ 28,
    no 2^996 hazard); G from erfc_tail_data; the log term needs only
    double-accuracy ABSOLUTE (2^-51 vs budget 2^-47.8 worst) — LogDd hi
    suffices. F′ = −2x0/(√π·G) — the e^{−x²} cancels in erfc′/erfc, no
    exp at all. F″ = −2x0·F′ − F′². After the dd residual, Halley runs
    in plain double: δ ~ 2^-19·x needs only 2^-40 relative.
    Residual space is IMPOSSIBLE here: for subnormal s the residual
    erfc(x0) − s sits below the 2^-1074 absolute floor; log space is
    immune and w is already paid for.
  - **mid** (x < ~6): residual space. F = ErfcDd(x0) ⊖ s (s exact; dd
    subtraction absorbs the ~2^-13-relative cancellation with ~2^-77 to
    spare). Needs ErfcDd: expose the erfc kernel's compensated
    1 − (e_hi + small) assembly pre-rounding — mechanical, erf_core
    already returns the pair; erfc public results must stay
    bit-identical (existing erfc gates are the regression guard).
    F′ = −2e^{−x0²}/√π with ops::Exp — backend exp is NOT
    accuracy-critical: its few-ULP error is attenuated by the 2^-19
    step to ~2^-69. Same Halley factor via erfc″ = −2x·erfc′.
  - x1 = fl(x0 + δ): one rounding; pre-rounding budget ≤ 2^-56,
    verified empirically by the generator self-check (house rule:
    trust the replay, not the derivation).
- Mid/far routing by t against a constant cut (deterministic in s);
  G extrapolated < 1e-5 interval-widths below x = 6 — generator
  confirms it stays inside slack.
- **erfinv never reaches the far tail**: max |erfinv| =
  erfcinv(2^-53) ≈ 5.86 < 6. Only erfcinv exercises x ∈ [6, 27.214];
  reference sets must cover the far tail through erfcinv.
- Specials: erfinv(±1) = ±inf, |y| > 1 → NaN; erfcinv(0) = +inf,
  erfcinv(2) = −inf, z ∉ [0, 2] → NaN; NaN propagates. Zero crossing
  at z = 1 gets bit-neighbourhood ULP testing like lgamma's zeros —
  the exact 1 − z argument is what makes relative accuracy hold there.
- Targets: C ≤ 1 ULP; T ≤ 1 ULP with low not-CR (near-CR expected from
  the 2^-56 budget). Gates set to measured values, no margin.
- Bench baseline: libm HAS no erfinv/erfcinv — pick and label a scalar
  baseline in implementation (scalar walk of our own kernel is
  acceptable; say so in the bench header).
- Left to the generator sweep: seed split points and degrees; Newton
  (deg ~10 seed) vs Halley (deg ~5 seed) — whichever meets the 2^-56
  replay check, cheaper per bench wins; Pc dd-lead count; kTfar.
- Rejected: direct dd tail fits (degree cost + rounded-argument
  floor); residual-space Newton in the far tail (underflow); any step
  against the public double-rounded erf/erfc (central κ ≈ 1 passes
  their error straight through, and double rounding floors the step).
- Oracle: mpmath erfinv for C; root-find on log erfc (probe pattern)
  for T — document in ACCURACY.md.

### Regularized incomplete gamma P/Q — broad parameters
Elementwise span API gamma_p(a, x, out) / gamma_q(a, x, out) in v1
(scalar-a broadcast overload later if profiling justifies). Always
compute the smaller of P/Q directly, complement the other.
- Region map: power series for P (x ≲ a+1, a below a_T) at FIXED
  length; continued fraction for Q (x above, a below a_T) at FIXED
  depth — both lengths proven at the region boundaries by the
  generator, not sampled; Temme uniform asymptotic for a ≥ a_T (~30):
  P/Q = ½erfc(∓η√(a/2)) + e^{−aη²/2}·S(η, 1/a)/√(2πa), fixed-length —
  the SIMD-friendly branch, and the payoff of Phases A/B: it consumes
  corvus erfc, exp_dd, log_dd, and lgamma's Stirling pieces.
- The ridge x ≈ a: η²/2 = φ(λ) = λ − 1 − log λ, λ = x/a, from
  u = (x − a)/a — x − a is exact by Sterbenz exactly where it matters
  (a/2 ≤ x ≤ 2a) — with a dd log1p. No naive a·log x − x cancellation
  anywhere.
- Prefactor x^a e^{−x}/Γ(a): for a ≥ 8,
  a·log x − x − lgamma(a) = −a·φ(λ) + ½log(a/(2π)) − φst(a), reusing
  kLgammaStirCoef verbatim; for a < 8 use lgamma's internal
  pre-rounding dd (LgammaPosDd — needs internal exposure, no public
  change). exp via ExpDdFrac with late scaling so subnormal tails stay
  one rounding.
- Oracle: mpmath regularized gammainc. The 2D reference strategy is
  genuinely new: per-region grids, ridge stress lines, boundary
  brackets, corners.
- Target: set from a ridge survey during detail design; expectation
  ≤ 4 ULP on the directly computed side, documented per region. Domain
  caps (huge/tiny a) decided and documented in detail design.

### Regularized incomplete beta — broad parameters (thin by intent)
I_x(a,b): symmetry I_x(a,b) = 1 − I_{1−x}(b,a) with "compute the side
≤ ½" routing (1 − x exact for x ∈ [1/2, 1]); CF at fixed depth as the
core; power series for small b·x; Temme uniform regimes for large a,b;
prefactor via Stirling-difference — the lgamma(a) + lgamma(b) −
lgamma(a+b) cancellation removed analytically with the same φ/φst
machinery, never by subtracting three rounded lgammas. Detail design
deferred until gamma ships and its η/φ dd utilities exist to reuse.

### Routing for Phase C [per AGENTS.md]
- Detail design + error budgets per family: frontier model, high
  effort.
- Implementation of a settled design (generator → kernel → tests →
  sweep): mid-tier model, default effort, with the AGENTS.md
  escalation triggers restated in the brief.
- Tier-sweep bookkeeping folds into the implementation run rather than
  a separate low-effort spawn.

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
- [RESOLVED 2026-07-25] erfc tail 5-ULP bound was entirely hn::Exp's
  contribution. exp_dd shipped and the tail rewired onto it: 5 -> 2 ULP.
  See the Phase A part 1 section above. log_dd (lgamma's dependency) is
  still outstanding; exp_dd's next consumer after erfc is the incomplete
  gamma/beta prefactor.
- [OPEN, decomposed 2026-07-25] **The erfc tail's 2 ULP and 48% not-CR are
  two different problems, and only one of them is cheap.** Measured by
  replaying the tail formula in mpmath against 3875 points in [6, 27.2],
  adding one error source at a time (exp exact throughout, since exp_dd
  contributes 2^-68):

  | model | max ULP | not-CR |
  |---|---|---|
  | A: stored coefficients, exact evaluation, exact args | 1 | 52.9% |
  | B: A + Horner in double | 2 | 51.8% |
  | C: B + u and s rounded | 3 | 54.6% |
  | shipped kernel | 2 | 48.5% |

  Readings:
  - **The not-CR rate is set by the fit and its double coefficients, and
    nothing else.** Model A is already 52.9% with everything downstream
    exact. No amount of compensated evaluation touches it; only a re-fit
    with dd coefficients would, and then the whole path has to stay dd.
  - **The max ULP is set by the evaluation.** Model A reaches 1 ULP, so a
    dd Horner would take the shipped kernel 2 -> 1 while leaving not-CR at
    ~53%. Cost: ~11 dd steps on a tail path already 1.7x slower than the
    hn::Exp version. Poor trade — it moves a bound that is already
    documented and does not move the quality signal that matters.
  - The kernel sitting at 48.5% (below model C's 54.6%, above B) confirms
    its dd argument and assembly work is buying something real.

  Decision: leave it. 2 ULP documented is honest and is better than most
  production libms in this region. **Revisit when erfcinv is scheduled** —
  Newton refinement there would be bounded by exactly this, which is the
  first time the tail's accuracy compounds into another function rather
  than sitting at a leaf. (2026-07-25: now scheduled, and the Phase C
  condition analysis says the bound goes the other way — erfcinv
  attenuates the tail error by 1/(2x²), so no compounding. Measure during
  the erfinv/erfcinv work, then close this as "no re-fit needed.")
- [OPEN] `HWY_DYNAMIC_DISPATCH` must be invoked from **inside** the
  namespace holding the per-target functions. With a single compiled
  target (the SSE2 cap) Highway collapses it to `N_SSE2::FUNC`, so a
  call written `HWY_DYNAMIC_DISPATCH(corvus::Foo)` from global scope
  becomes `N_SSE2::corvus::Foo` and fails to compile — at that tier only.
  Found in test_exp_dd during the 2026-07-25 sweep and fixed there; the
  shipped kernels already do it correctly. Worth a look if a future
  test/bench is written from global scope. Related: the SSE2 build failure
  left the SSSE3 binaries in place and the sweep script happily re-ran
  them — `CORVUS_EXPECT_TARGET` caught it, and `tools/sweep_tiers.ps1` now
  aborts on the first build failure rather than relying on that.
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
  docs/ACCURACY.md's NEON column. (That column's original caveat — "gates
  may trip if hn::Exp accuracy differs on NEON" — expired on 2026-07-25:
  no corvus kernel calls hn::Exp any more. The live NEON question is now
  exp_dd/log_dd, which CI picks up automatically since both are
  ctest-registered.)
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
4. [RESOLVED 2026-07-25] Phase A part 1 — exp_dd, the dd primitive layer,
   the facade additions, and the erfc tail rewire; validated on all five
   x86 tiers on the Ryzen and documented. See the Phase A section above.
5. [RESOLVED 2026-07-25] Phase A part 2 — log_dd, validated on all five
   x86 tiers and documented. **Phase A is complete**: corvus owns both
   transcendental cores it needs, and no accuracy-critical kernel depends
   on Highway contrib any more.
6. [RESOLVED 2026-07-25] Phase B — lgamma, full real axis, validated on all
   five x86 tiers under GCC and MSVC. See the Phase B section.
7. **Start here**: Phase C detail design, erfinv/erfcinv first. Broad
   design and build order are resolved and recorded in the Phase C
   section above (2026-07-25): erfinv/erfcinv → incomplete gamma P/Q →
   incomplete beta. Detail design is frontier-model work; the main open
   choice is direct dd fits vs seed + dd-Newton for the tail core's mid
   region. Implementation after that is mid-tier per the routing
   section.
8. [RESOLVED 2026-07-25] NEON — CI run 30163646136 (all three jobs green)
   reproduced every number point-for-point: exp_dd 2^-68.45 at the same
   worst-case input, log_dd 2^-67.88, erfc tail at the new 2-ULP gate,
   erf/erfc not-CR counts unchanged. ACCURACY.md's NEON row is filled.
   Worth noting AppleClang's default contraction (`on`) differs from GCC's
   (`fast`) and MSVC's (none), so this agreement doubles as a check that
   CORVUS_FP_FLAGS reaches every compiler.
