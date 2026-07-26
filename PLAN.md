# corvus — Plan / Session State

## Status [DERIVED] — 2026-07-25 (Ryzen)
**Phases A and B are complete; Phase C is in progress.** Scaffolded
2026-07-20; public at github.com/OldCrow/corvus. Shipped and
production-quality (per-tier audit record: docs/ACCURACY.md):

- **erf** max 1 ULP over the full domain. **erfc** max 1 ULP for
  |x| <= 6 and subnormal results, 2 ULP in the normal tail
  (fit-limited; decomposition in Open Items). **lgamma** over the full
  real axis: max 1 ULP on the positive axis including the exact zeros,
  correctly rounded throughout Stirling, 1 ULP / 2^-53 absolute on the
  negative axis. Internal dd cores **exp_dd** (2^-68.45 relative) and
  **log_dd** (2^-67.88, correctly rounded on every reference point).
  No accuracy-critical path depends on Highway contrib math.
- Validated on real silicon: AVX3_ZEN4/AVX3_DL/AVX3, AVX2, SSE4, SSSE3,
  SSE2 (Ryzen native + capping, under GCC, MSVC and clang-cl), NEON
  (Apple Silicon CI), Kaby Lake AVX2. The dd cores are bit-identical on
  every tier and compiler; erf/erfc differ on the no-FMA tiers in
  not-CR counts only.
- **Phase C part 1 (erfinv/erfcinv) shipped 2026-07-25.** Max 1 ULP
  everywhere on all five validated x86 tiers, including erfcinv's far
  tail down to the smallest subnormal result. Full design, three bugs
  found by the reference sweep, and measured numbers are in the Phase C
  part 1 section under Shipped phases below. Incomplete-gamma detail
  design is the next frontier task.

## Next Steps
1. **Start here**: detail design for regularized incomplete gamma P/Q (frontier
   effort). The broad section below lists what it must settle: a_T and
   the series/CF boundaries with PROVEN fixed lengths, the ridge
   accuracy survey that sets the public target, domain caps.
2. Regularized incomplete beta after gamma (reuses its η/φ dd
   machinery).
3. Quiet-machine bench pass before publishing any performance number —
   everything so far is session-loaded and labeled indicative.

## Open Items
- [RESOLVED 2026-07-25] **The erfc tail's 2 ULP and 48% not-CR are two
  different problems, and only one is cheap.** mpmath replay of the
  tail formula over 3875 points in [6, 27.2], adding one error source
  at a time (exp exact throughout — exp_dd contributes 2^-68):

  | model | max ULP | not-CR |
  |---|---|---|
  | A: stored coefficients, exact evaluation, exact args | 1 | 52.9% |
  | B: A + Horner in double | 2 | 51.8% |
  | C: B + u and s rounded | 3 | 54.6% |
  | shipped kernel | 2 | 48.5% |

  The not-CR rate is set by the fit's double coefficients alone (model
  A); the max ULP is set by the evaluation (a dd Horner would buy
  2 -> 1 at ~11 dd steps on an already-1.7x-slower path — poor trade).
  Decision: leave at 2 ULP, documented. **Closed by measurement**:
  erfcinv's far-tail Halley step (Phase C part 1) reuses the same G fit
  and lands at max 1 ULP end to end, confirming the condition-number
  analysis — erfcinv attenuates the tail's error rather than compounding
  it. No re-fit needed.
- [OPEN] CORVUS_SANITIZE is not MSVC-aware: emits `-fsanitize=<list>`
  unconditionally, which cl.exe rejects. Harmless while sanitizer
  builds are Linux-only; branch on MSVC or reject with FATAL_ERROR.
- [OPEN] Install/export when Highway is FetchContent-built is disabled
  (the exported target would dangle). Decide before the first tagged
  release: require system hwy (status quo), bundle hwy objects into
  libcorvus.a, or install a nested hwy.
- [OPEN] Pre-release legal: binary artifacts linking Highway must carry
  its Apache-2.0 NOTICE; source-only distribution needs nothing. Handle
  when packaging starts.
- [OPEN] Decide whether libstats/libhmm adopt corvus as a dependency or
  keep their internal SIMD (separate project-level decision).
- [OPEN] If a sibling project adopts corvus, verify clang-cl-built
  corvus links cleanly into an MSVC-built consumer before relying on
  it — same-ABI is the design intent but is untested.
- [OPEN] Upstream path for HWY_BROKEN_MSVC: add a compiler-version
  floor like its siblings have. Needs Highway's own suite passing under
  MSVC with AVX-512; two kernels are not sufficient evidence for a PR.
- [OPEN, low priority] Non-gather x86 kernel variant: ~2x upside on
  gather-weak pre-AVX-512 CPUs (Kaby Lake class); Zen 4 scales fine
  without it (see Resolved log).
- [OPEN] **lgamma is below scalar-libm parity at 2 lanes — by design
  economics, not defect.** Measured 2026-07-25 on Ryzen/SSE2 (loaded,
  indicative): 0.46–0.82x vs mingw libm (zone worst — the degree-34
  Horner with per-lane selects; reflection best at 0.82x). M1/NEON
  reported 0.2–0.4x, consistent with the same 2-lane cost divided by
  Apple's faster vendor lgamma. The dd-heavy 1-ULP design pays in
  flops and wins by lane count: 8 lanes 1.9–3.2x, 2 lanes < 1x.
  Options if narrow-width parity ever matters: cheaper zone path
  (split intervals to cut degree), or simply document that 2-lane
  targets favor scalar libm for lgamma. erf/erfc expected unaffected
  (light kernels vs expensive scalar) — verify with bench_erf/erfc on
  the M1 before concluding anything is wrong there.
  **M1/NEON cross-check run 2026-07-25 (loaded, indicative): confirms
  lane economics, no NEON-specific defect.** erf 5.5–6.5x and erfc
  core 3.6–3.7x vs Apple libm — healthy, and erf's table gathers cost
  nothing measurable (2.35 ns/el including emulated gathers), ruling
  out NEON's lack of hardware gather as a suspect. erfc tail-only
  0.95–0.97x decomposes exactly: corvus per-element NEON/AVX3_ZEN4 =
  10.96/4.86 = 2.25x, Apple-vs-mingw erfc baseline = 10.6/~23.3 =
  2.2x, product ≈ the full 4.8x→0.96x ratio gap. Per-VECTOR the dd
  tail is ~70 cycles on M1 vs ~180 on Zen 4 — NEON does the same dd
  work in fewer cycles and loses only by doing 4x fewer elements per
  vector. lgamma NEON 0.22–0.43x vs SSE2's 0.46–0.82x is the same
  kernel cost divided by Apple's 2–3x faster vendor lgamma
  (8.5–12 ns/el positive axis, 31–34 negative). No tuning action on
  NEON; the option list above stands unchanged.
- [ILLUSTRATIVE] Possible future consumers: C++ port of multi-agent_sim
  (batch distance/trig), zeekhmm training pipelines.

## Decisions
- Name: corvus (OldCrow tie-in). Namespace `corvus::`.
- Scope: statistical special functions only; basic transcendentals
  (exp/log/trig/pow) belong to Highway contrib. Families (P0 first):
  erf/erfc, erfinv/erfcinv, lgamma, regularized incomplete gamma P/Q,
  regularized incomplete beta; (P1) digamma, inverse incomplete
  gamma/beta, Bessel I0/I1.
- Backend: Highway behind the `src/ops-inl.h` facade; public API is
  std-only. std::simd migration = facade reimplementation, deferred
  until implementations mature (mid-2026: GCC 16 partial, no libc++).
- Dependency model: find_package(hwy) preferred, FetchContent fallback
  pinned to the audited version (1.4.0 — bump only with revalidation).
  Ship model: static lib, PIC on, Highway not exposed. MIT, clean-room
  implementations only.
- Naming, file-extension, header and documentation conventions, and the
  CMake standard: AGENTS.md is canonical (user-driven audits
  2026-07-21). Documentation capped at four files. Deferred
  deliberately: LTO/IPO (profile first), shared lib + symbol visibility
  (no demand).
- **FP contraction off project-wide** (2026-07-25): GCC's default
  `-ffp-contract=fast` fused inside log_dd's dd identities and shifted
  it 0.6 bits vs MSVC on identical source — and made it BETTER, which
  is why nothing failed and one compiler alone could never notice. The
  dd layer's exactness proofs assume IEEE ops as written, and
  ACCURACY.md claims cross-compiler reproducibility, so contraction
  cannot be left to the optimizer. `CORVUS_FP_FLAGS` sets
  `-ffp-contract=off` (clang-cl: `/clang:-ffp-contract=off`; MSVC
  `/fp:precise` already doesn't contract), PRIVATE to corvus and to the
  kernel test targets — a test measuring kernel code must measure the
  program that ships. Principle: fusion is requested in source
  (`ops::MulAdd`), never inferred. Cost <= ~8%, indicative.
- Windows compiler (2026-07-25, resolved with user): documented
  exception to the fleet MSVC default — Windows validation and bench
  numbers come from clang-cl (preferred: MSVC ABI) or mingw GCC,
  because `HWY_BROKEN_MSVC` caps MSVC dispatch at AVX2. MSVC remains a
  fully supported consumer toolchain and runs the CI Windows job. Full
  rationale and scope in AGENTS.md.
- Platform tiers: Tier 1 (accuracy-audited on real silicon) = NEON
  (M1), AVX-512/AVX2/SSE2 family (Ryzen native + capping), AVX2 (Kaby
  Lake). Tier 2 (compiles, unaudited) = SVE and anything else Highway
  emits.
- Sequencing (2026-07-21, resolved with user): dd cores before lgamma
  (Phase A -> B), with the erfc tail rewire as Phase A's acceptance
  test. Both landed as planned.
- lgamma v1 scope (2026-07-21, resolved with user): full real axis,
  poles return +inf, sign output (signgam) deferred — SciPy's gammaln
  offers none either.

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

## Shipped phases — what would be expensive to re-derive
Full method and measured bounds live in docs/ACCURACY.md; kernel
derivation blocks carry the math. This section keeps only the design
points and bugs worth not re-learning.

### Phase A — exp_dd + log_dd [shipped 2026-07-25]
Acceptance met: erfc tail 5 -> 2 ULP on the unchanged reference set,
gate retightened. exp_dd 2^-68.45 relative, bit-identical on ALL tiers
including no-FMA (the Dekker ProdLow fallback reproduces the fused path
exactly here); log_dd 2^-67.88, correctly rounded on all 7273 reference
points, every tier, all three compilers.
- exp_dd returns **mantissa + exponent** (ExpDdFrac + ScaleTwo), not a
  scaled value: consumers fold their own factors in first and the
  power-of-two scaling lands last — that is what keeps subnormal
  results at one rounding.
- log_dd's mantissa is **centred on 1** (halving at slot boundary
  1 + 53/128) and the two slots adjacent to m = 1 carry R = 1, L = 0
  exactly: no special case and no cancellation on either side of x = 1.
  r = R_j·m − 1 is exact (Sterbenz + ProdLow); a plain fma would carry
  2^-62.
- **An exact dd pair is not necessarily a NORMALIZED one.** (p−1, p_lo)
  had |lo| up to 2^-53 against hi ~2^-8; terms computed from hi alone
  dropped their share of lo (surfaced as 2^-63.3, only near x = 1).
  Fix: one TwoSum — not Fast2Sum, r crosses zero inside every slot.
  The gate caught it; reading the code did not.
- The internal-kernel test pattern and the dd-pair oracle pattern
  established here are documented in AGENTS.md (Architecture).
- Cost of the erfc rewire, indicative: tail-only 8.4x -> 4.8x vs libm —
  the predicted two-gathers-plus-dd-math trade.

### Phase B — lgamma [shipped 2026-07-25]
The 2026-07-21 design (full text in this file's git history) was
architecturally complete and was followed: five regions, exact-t zeros
form, masked fixed-step recurrence, swept-X0 Stirling, sinpi
reflection. Four things it did not settle:
1. **The recurrence direction is forced, and it is DOWN**: x − k is
   exact for every integer k < x < 2^52; x + k is not (walking up would
   put ~77 ULP into lgamma(2.5) without dd shifts). Down also lands
   every lane in (1.5, 2.5] — one polynomial — and makes the design's
   P ≤ Γ(X0) bound exact: 7·6·5·4·3·2 = Γ(8).
2. **Reflection via Γ(1−x) = −x·Γ(−x)**: −x is exact where 1 − x is
   not, and the positive pipeline runs ONCE on |x| with reflection as a
   post-correction — no second evaluation, no lane divergence.
3. **π/|sin πx| is never formed** (it overflows at every pole): factor
   sin(πu) = πu·sinc to cancel π analytically, leaving −log|u| — the
   term that should diverge — plus a bounded u² fit.
4. **lgamma(1) needs an explicit +0.0**: t·B(t) at t = 0 carries B's
   sign (−γ). Adding +0.0 is the identity on every finite nonzero value
   and also removes dependence on how a target's Fast2Sum signs zero.
The bug worth remembering — **ops::ProdLow's non-FMA Dekker split
overflows for |a| > 2^996** — is now an AGENTS.md convention. Stirling's
x·(log x − 1) was the first corvus operand to reach that range; only
the capped SSE sweep caught it (500/4242 points wrong while every FMA
target, GCC and MSVC stayed green). Fix: scale 2^-200/2^200 around the
product — exact, one code path for all targets.
Also: Stirling is grouped x·(log x − 1) − log(x)/2, not the textbook
(x − ½)·log x, whose product exceeds the result by x and overflows —
dd residual to NaN — over a ~0.1%-wide band below the true threshold.
And two generator traps, both commented in gen_lgamma_data.py: mpmath's
fsum mis-sums a generator argument, and an odd Chebyshev node count
puts a node where 1 + t rounds to 1 (coefficients plateau, reading as
"not smooth" instead of "tooling bug").
Parameters from the sweep: X0 = 8 (ψ degree 8), zone split 3/2, outer
boundary 1/2, degrees 34/21 with three dd leading coefficients each.
Cost, indicative (loaded Ryzen): 1.9–3.2x vs libm; the degree-34 zone
Horner is the slow region.

### Phase C part 1 — erfinv/erfcinv [shipped 2026-07-25]
Design (this file's git history, "Phase C part 1 — erfinv/erfcinv design")
was followed as written: C/T routing by Sterbenz-exact arguments, central
direct-polynomial fit, tail seed + one dd Halley step split at the same
x = 6 threshold erfc's own core/tail split uses. **Halley over Newton**,
per the design's own resolution (seed degrees 7/9/5 against a 2^-19 target,
cheaper than a Newton-sized ~10-degree seed for the same 2^-56 end-to-end
budget) — confirmed by `gen_erfinv_data.py`'s replay self-check, not
re-swept independently. Measured: **max 1 ULP on every region, every
validated x86 tier** (AVX3_ZEN4/AVX3_DL/AVX3, AVX2, SSE4, SSSE3, SSE2);
not-CR counts and method detail in docs/ACCURACY.md.

Three bugs, all found by the reference set's boundary-neighbourhood points
rather than by reasoning about the code, and all the same shape as the
lgamma/log_dd ones: **an identity that is exact in general does the wrong
thing at one distinguished point, and only testing that exact point catches
it.**
1. **ErfcCoreDd's table clamp needed widening.** It clamped its argument to
   exactly 6.0 -- safe for erfc.cpp (which only ever keeps results for
   |x| <= 6, discarding anything computed above it), but erfinv's mid-region
   seed can legitimately land a few 2^-17-ish PAST 6 when the true root sits
   right at the mid/far seam (seed error ~2^-19 relative; 6*2^-19 is many
   ulps at that magnitude). Clamping silently evaluated erfc at 6.0 instead
   of the seed's actual x0, decoupling f from f' and producing a ~1e-6 error
   (not 1 ULP) exactly at the boundary reference points built to probe it.
   Fixed by widening to `kErfcCoreSafeMax = 6 + 1/1024` -- still 500x inside
   the erf table's real safe-extrapolation limit (6 + 1/512, where the grid
   index would go out of bounds) and providably a no-op for erfc's own
   public results (every lane the wider clamp changes is one erfc.cpp always
   discards via its own tail/core select).
2. **Sign of zero, two places.** `ErfinvCentral`'s dd assembly adds a +0 and
   a -0 partial sum together internally; IEEE 754 defines (-0)+(+0) as +0 in
   round-to-nearest, so `erfinv(-0)` silently came back +0 without an
   explicit trailing `CopySign`. Separately, erfcinv's C-branch argument had
   been written `Neg(z - 1)` rather than `Sub(1, z)` -- identical for every
   z except z = 1, where `x - x` is always +0 but negating that +0 gives -0,
   so `erfcinv(1)` came back -0 instead of the design's stated +0.
3. **The reference oracle's own initial guess blew up at s = 1.** Root-
   finding `log(erfc(x)) - log(s)` for `erfcinv`'s reference points used an
   initial guess built from `log(-log(s))`; at s = 1, `-log(1) = 0` and
   `log(0) = -inf` fed into a `max(..., inf)`, sending the guess to infinity
   and making mpmath's solver return near-zero garbage (~1e-227) rather than
   exactly 0. Fixed by using `mpmath.erfinv(1 - s)` directly whenever
   s >= 1/2 (safe there: 1 - s is never subnormal-relative to 1 in that
   branch) and reserving the log-space solve for s < 1/2, where it is
   needed for precision (forming 1 - s for a subnormal s would need ~1075
   bits just to distinguish it from 1).

A fourth defect was caught in review, not by any test, because it is
invisible in results: **HalleyMid passed a possibly-NaN x0 into
ErfcCoreDd's table gather.** Discarded lanes (erfinv(NaN), erfcinv(z < 0),
±inf inputs) route s ≤ 0 into ErfcInvCore, where sqrt(−log s) goes NaN —
and ErfcCoreDd's gather index is round(ac·256), NOT masked. It never
misbehaved only by unrelated platform accidents (x86 minpd returns the
non-NaN second operand; ARM fcvtzs(NaN) = 0), neither of which is a
guarantee, and Highway's debug-mode gather bounds assert could trip.
Fixed with the same one-op NaN scrub erfc.cpp itself uses; measured
values unchanged on every tier (the scrub only affects discarded lanes).
The general rule is now in AGENTS.md Conventions: masked-off lanes still
EXECUTE gathers, so any value-derived gather index must be scrubbed.

Bench, Ryzen/GCC/AVX3_ZEN4, scalar-walk-of-own-kernel baseline (libm has no
erfinv/erfcinv), session-loaded so indicative: erfinv central 22x, erfinv
mixed [-0.999,0.999] 7.3-7.7x, erfcinv central 53-54x, erfcinv T-mid 11-11.5x,
erfcinv T-far (subnormal z) 12x. The "scalar" baseline calls the public API
with length-1 spans per element, so it is not a clean scalar-vs-SIMD
comparison — it includes per-call dispatch/span overhead the real kernel
loop doesn't pay, which inflates every ratio above; treat these as upper
bounds, not the kernel's true SIMD speedup.

### erf + erfc [shipped 2026-07-20/21]
erf: table + local-Taylor kernel (clean-room port of libstats
vector_erf_neon through the facade), max 1 ULP. erfc: core reuses the
erf table via compensated 1 ∓ erf assembly; tail is
e^{−a²}·G(1/a)/a on exp_dd, three fitted intervals with coefficient
select. Lessons that outlived the phase: the shared erf series must
meet erfc's RELATIVE precision near a = 6, not erf's absolute one
(series extended to d^8, erf unchanged); and the MulSub/no-FMA
exact-residual hazard that created ops::SquareLow's capability guard —
now an AGENTS.md convention.

## GitHub repo settings [applied 2026-07-21 via gh api]
Merge: all three styles, auto-delete head branches (matches libstats).
Wiki and projects DISABLED (four-file docs policy), issues on,
discussions off. Topics set. Security: Dependabot alerts + auto fixes,
secret scanning + push protection, private vulnerability reporting.
Ruleset "protect-main": blocks force-push/deletion, direct pushes
allowed (solo workflow). Actions GITHUB_TOKEN read-only, cannot approve
PRs. Deferred: signed-commits rule (confirm the M1 and Ryzen boxes sign
before enabling), tag-protection ruleset for v* at first release.
Required status checks deliberately absent (incompatible with
direct-push workflow).

## Build-stack standardization (2026-07-23) [DERIVED]
Cross-repo effort tracked in ~/Development/BUILD-STANDARDIZATION-PLAN.md.
corvus commits: pkg-config file + consumer example + installed-path CI
check; find_package(hwy 1.4) version floor with CI building pinned
Highway 1.4.0 from source; CMakePresets.json. AGENTS.md's CMake section
verified still accurate afterwards.

## Resolved log
One line per closed item; detail lives in this file's git history,
AGENTS.md, and docs/ACCURACY.md.
- 2026-07-25 Phase C part 1 (erfinv/erfcinv) shipped (0ed13ab), max
  1 ULP on all five validated x86 tiers AND NEON — CI run 30180799151,
  all three jobs green, NEON point-identical to the x86 FMA tiers;
  ACCURACY.md NEON row filled. erfc-tail open item closed by
  measurement (attenuated, not compounded, as the condition analysis
  predicted). One watch item: the Windows/MSVC CI job took 16 min vs
  its usual ~5.5 (the new erfinv TU is the heaviest yet, plus hosted
  2-core runner variance) — worth a look only if it repeats.
- 2026-07-25 Sibling FP-contraction audit delegated to its own trackers:
  OldCrow/libstats#84 and OldCrow/libhmm#70 carry the full context and
  own the question from here.
- 2026-07-25 Phase B (lgamma) shipped; NEON row filled from CI run
  30170907111.
- 2026-07-25 Phase A (exp_dd, log_dd, dd primitives, erfc tail rewire
  5 -> 2 ULP) shipped.
- 2026-07-25 HWY_DYNAMIC_DISPATCH must be called from inside
  `namespace corvus` (single-target collapse breaks global-scope calls
  at the SSE2 cap only) — in AGENTS.md; sweep_tiers.ps1 now aborts on
  first build failure so stale binaries can't be re-measured.
- 2026-07-24 AVX-512 validated natively on the Ryzen (GCC + clang-cl,
  point-identical to AVX2/NEON); no-FMA divergence measured: not-CR
  counts only, no bound moved. Git-Bash/libstdc++ gotcha in AGENTS.md.
- 2026-07-24 Windows/MSVC CI job added as toolchain coverage; /WX
  caught C4996 in expect_target.h, and the VS generator exposed the
  set_property(CACHE CMAKE_BUILD_TYPE) multi-config bug — both fixed.
- 2026-07-24 CI asserts its tier: CORVUS_EXPECT_TARGET everywhere after
  three silent-cap defects; every sweep iteration capped by name
  including HWY_AVX10_2.
- 2026-07-24 HWY_BROKEN_MSVC investigated end-to-end;
  CORVUS_MSVC_UNBLOCK_AVX512 exists, default OFF — mechanism, scoping
  and policy in AGENTS.md.
- 2026-07-21 CI designed around runner-minute economy (Linux tier sweep
  + sanitizers, macOS NEON); NEON validated for erf/erfc; repo settings
  applied; conventions + build-system audits recorded in AGENTS.md.
- 2026-07-20 x86 gather performance: Kaby Lake is flat across widths
  (its gather + emulated f64->i64); Zen 4 scales properly (SSE2 5.49 ->
  AVX2 2.66 -> AVX3_ZEN4 1.95 ns/el on erf). Non-gather variant stays a
  noted upside for pre-AVX-512 hardware only.
- ULP-harness generalization: the per-kernel generator/reference/gate
  pattern is established across all five shipped kernels.
