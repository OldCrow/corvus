# corvus — Plan / Session State

Detail policy: this file holds STATE — what's decided, what's open, the
next concrete step, and lessons too expensive to re-learn. Full designs,
session narratives, and measurement play-by-play live in this file's git
history (compacted 2026-08-06), docs/ACCURACY.md, and the kernel/generator
source, which are the official record for finished work.

## Status [DERIVED] — 2026-08-06

**v0.1.0 RELEASED** (tag at b4eaeea, immutable under protect-tags;
https://github.com/OldCrow/corvus/releases/tag/v0.1.0). P0 is complete:
erf/erfc, erfinv/erfcinv, lgamma, gamma_p/q, beta_p/q — all shipped,
audited per tier on real silicon (docs/ACCURACY.md), reference oracles
independently certified (oracle-trust directive discharged for both
gamma and beta). CI green on Linux (tier sweep + sanitizers + install
contract), macOS arm64 (NEON), Windows (MSVC); lint-workflows adopted
(issue #2 closed). One branch (main).

## Next Steps
1. **NEXT SESSION (fresh fork): Bessel I0/I1 + lbeta (P2)** — staged
   in full below (probe questions, API decision points, effort
   routing); open at frontier with the staging block. BETA_P_INV /
   BETA_Q_INV SHIPPED [2026-08-10] (ledger: probe 0 / G1 3 ratified
   corrections incl. one orchestrator-review catch / G2 0 design +
   2 scope continuations + 3 frontier rulings / G3 1 adjudicated
   escalation + nine accepted deviations / G4-G5 1 CI warning fix).
   GAMMA_P_INV / GAMMA_Q_INV SHIPPED [2026-08-09] (ledger: probe 0 /
   G1 2 ratified corrections + frontier takeover at chain depth 3 /
   G2 0 / G3 0 with nine accepted deviations). P1 COMPLETE — corvus
   now covers every PDF/CDF/quantile/MLE need in the libhmm/libstats
   inventories except von Mises (Bessel, P2).
   Pipeline template: probe→design→G1/G2→G3→G4/G5 (digamma and
   trigamma both shipped through it, one session each).
   Escalation-density rule [2026-08-06, user]: ~3 chained
   escalations in a delegated stage → default back to frontier.
   Pipeline ledgers: digamma probe 0 / G1 1 / G2 0 / G3 0; trigamma
   probe 0 / G1+G2 1 (+1 process fault) / G3 0 — both families' sole
   design escalations were probe-grid artifacts caught by the next
   stage's replay, and the edge-refined-sampling rule from
   trigamma's FIRST correction is binding for all future families.
   DIGAMMA + TRIGAMMA SHIPPED [2026-08-08]: with them, every MLE
   update in libhmm/libstats has its ψ/ψ′ pair; after the inverse
   pair lands, corvus covers every PDF/CDF/quantile/MLE need in both
   consumers' current inventories except von Mises (Bessel, last).
2. **ROADMAP GAP ANALYSIS [2026-08-08, vs the actual libhmm +
   libstats distribution inventories]**: one REQUIRED addition —
   **trigamma** (public ψ₁): every second-order MLE in both consumers
   needs it alongside digamma (Gamma shape Newton, Beta bivariate
   Newton on ψ′(a)/ψ′(b)/ψ′(a+b), NegBin r-update). Added to P1;
   cheapest family left (digamma architecture minus the root
   product-form — ψ₁ has no positive zero; reflection is the SUM
   identity ψ₁(x) + ψ₁(1−x) = π²/sin²(πx) on existing sinc-pair
   machinery; the internal rough-trigamma is the structural sketch).
   Slot at user's discretion: warm-up before the inverse pair or
   immediately after. REFINEMENT to the planned Bessel item: must
   include exponentially-scaled variants (i0e/i1e or log-I0) — I0
   overflows past κ ≈ 713 and von Mises log-density at large κ is a
   primary consumer use. lbeta PROMOTED to committed P2 for v1.0.0
   [2026-08-10, user decision after the milestone/issue sweep]: the
   BetaBinomial PMF (libstats v2.4.0, #62) needs two ln B per point
   in its hot path on top of the F/StudentT/Binomial delegations;
   consumers currently form ln B as three lgammas — the a+b
   cancellation hazard corvus solved internally via LgammaDiffDd.
   Nearly free: expose the internal machinery as one thin public
   kernel; slot after Bessel, possibly same session. [OPEN, P2
   candidate, not required]: erfcx (nearly free from the erfc tail
   machinery; speculative consumer found 2026-08-10 — TruncatedNormal
   far-truncation moments/MLE run on the Mills ratio = erfcx — but no
   filed need yet). Sweep also confirmed: libstats #47 (A&S Bessel
   fallback capping VonMises at ~1e-7) is exactly what P2 Bessel
   retires, and libstats #52 (slow Binomial CDF) is beta_p — an
   integration note for libstats, not a corvus gap.
2. **Quiet-machine bench_beta re-run** [bench SHIPPED 2026-08-06 —
   numbers below are loaded/indicative]: re-run on an idle Ryzen for
   publishable numbers, and fold into the Kaby bench pass when that
   machine leg happens.
3. **Kaby Lake legs when the machine is available** [OPEN, machine
   access; ruled NON-GATING 2026-08-06, user decision]: beta AVX2-native
   + capped sweep (additional cross-machine check, not a claim gap —
   ACCURACY.md dagger note), and the quiet-machine Kaby bench pass (its
   gamma bench numbers are loaded/indicative only).

## Open Items
- [OPEN, PRIORITY — 2026-08-10, found by betainv G3] TWO defects in
  the SHIPPED beta forward (beta_p/beta_q), PB prefactor's cpsi →
  Log1pmxDd at u → −1 (u = −λ/α): (1) 1+u < 2⁻⁵³ → u.hi rounds to
  exactly −1 → LogDdAny of a zero-high pair → NaN path → beta_p(19,
  1e5, 5.204222470155122e-21) returns EXACTLY 0, truth 3.36e-308
  (boundary y ≈ 2.1e-20); (2) 1+u merely small → Log1pmxDd adds u.lo
  into an exact TwoSum's low word and LogDd(Dd) keeps only the
  quadratic, t³/3 survives → 1.4e-4 error in E at (19, 1e5,
  1.73e-19), verified against a correct reference row. The forward
  reference set never sampled this corner (deep tail at moderate-a/
  huge-b). betainv is IMMUNE (routes E to PA except where gated).
  Needs its own fix arc: kernel correction + reference rows covering
  the corner + full revalidation. Scheduling at user's discretion.
- [OPEN, enhancement, low priority — 2026-08-10] exp_dd accuracy
  bump (one more polynomial term + keep r.lo through the quadratic)
  would raise betainv's y-ULP κ-horizon from 2¹⁸ toward the design's
  2⁵²; shared-core change, full-fleet revalidation; only worth it if
  a consumer needs y-ULP in the plateau band κ ∈ (2¹⁸, 2⁵²) — those
  rows already meet the backward contract at 0.000 ulp(σ).
- [OPEN, no rush — user, 2026-08-09] docs/ARCHITECTURE.md (layering
  diagram, added by user): decide whether README (users/maintainers)
  and/or AGENTS.md's reading map (agents) should reference it; if the
  reading map takes it, as a load-on-demand visual-reference entry,
  not always-read content.
- [WATCH] Ryzen box stability: two GPU-stack bugchecks 2026-07-29 (0x9F
  power-IRP, 0x10E video memory), root-caused [DERIVED] to a
  half-committed NVIDIA driver install; DDU clean reinstall the same
  night. Crash-free since, including the full beta validation ladder
  2026-08-05. Residual watch: whether the NEXT NVIDIA App driver update
  completes cleanly (chronic installer freezes implicate accumulated App
  state; manual driver-only installs are the fallback). Recurrence of
  either bugcheck on the clean stack flips suspicion to VRAM/hardware.
- [WATCH] mingw GCC 16.1 AVX-512 by-value-argument misalignment bug:
  filed upstream 2026-08-08 as GCC PR 126741
  (https://gcc.gnu.org/bugzilla/show_bug.cgi?id=126741). Local repro
  kept at `C:\Users\gdwol\Development\gcc-zmm-mingw-repro\`. Re-qualify
  mingw GCC for AVX-512 work only after the fix lands. Technical
  detail lives in AGENTS.md (Development Fleet).
- [OPEN] CORVUS_SANITIZE is not MSVC-aware (emits `-fsanitize=<list>`
  unconditionally). Harmless while sanitizer builds are Linux-only;
  branch on MSVC or reject with FATAL_ERROR.
- [OPEN] Pre-release legal: BINARY artifacts linking Highway must carry
  its Apache-2.0 NOTICE; source-only distribution needs nothing (v0.1.0
  is source-only). Handle when packaging starts.
- [OPEN] Decide whether libstats/libhmm adopt corvus as a dependency or
  keep their internal SIMD (separate project-level decision).
- [OPEN] If a sibling adopts corvus, verify clang-cl-built corvus links
  cleanly into an MSVC-built consumer (same-ABI is the design intent,
  untested).
- [OPEN] Upstream path for HWY_BROKEN_MSVC (add a compiler-version
  floor): needs Highway's own suite passing under MSVC with AVX-512;
  two kernels are not sufficient evidence for a PR.
- [OPEN, low] Non-gather x86 kernel variant: ~2x upside on gather-weak
  pre-AVX-512 CPUs (Kaby class); Zen 4 scales fine without it.
- [OPEN] lgamma performance-positioning wording: settle against a
  same-libm comparison before publishing any claim. Measured picture
  (loaded, indicative; full data in git history): the dd-heavy 1-ULP
  design wins by lane count — after the 2026-07-25 all-zone fast path
  (bit-identical, verified every tier), zone is 5.6–6.5x at 8 lanes,
  ~1.0x at 2 lanes; recurrence is the floor everywhere (0.2–0.6x vs
  fast vendor libms) and its cost is genuine work (X0 = 8 is
  accuracy-forced). "Wins from 4 lanes up" did NOT survive the Kaby
  4-lane measurement — the honest form is per-region and per-libm.
  Interval splitting (halve zone again) deferred; trigger = profiling
  gamma/beta end-to-end shows the zone Horner as a real bottleneck.
- [OPEN, gamma] Add a ≥ 2^998 Temme witness rows at gamma's next
  reference touch — the kGammaTwoPiAClamp Dekker-ceiling defect was
  latent for lack of them (clamp already fixed to 2^900, 2026-08-05).
- [OPEN, low] mingw-g++-built test binaries crash at process EXIT
  (0xC0000005 AFTER full PASS output, from PowerShell — not the known
  Git-Bash DLL-shadowing signature). clang-cl sweeps unaffected.
  Diagnose on a quiet day.
- [OPEN, pre-v1.0.0] Source-comment trim [2026-08-08, user]: many
  implementation comments embed operational history (stage tags,
  correction ordinals, session dates) whose authoritative record is
  PLAN.md/ACCURACY.md/git history — a drift hazard as code evolves.
  Sweep src/ + tools/ keeping math derivations, error bounds, and
  constraints at definition sites (those ARE the maintainer docs) while
  pruning narrative history to pointers. First instances fixed
  2026-08-08: three stale PROVISIONAL markers in beta-inl.h/beta_data.h
  that predated the G4 gate pinning. SAME PASS [added 2026-08-09,
  user]: review/compress docs/ENVIRONMENT.md and
  docs/NUMERICAL-DOCTRINE.md — the 2026-08-09 carve moved AGENTS.md
  content near-verbatim for diff auditability; the references still
  carry dated incident narrative that can shrink to rules + pointers
  once the split has proven itself in agent use.
- [OPEN] Signed-commits ruleset: confirm the M1 and Ryzen boxes sign
  before enabling.
- [ILLUSTRATIVE] Possible future consumers: C++ port of multi-agent_sim
  (batch distance/trig), zeekhmm training pipelines.

## Decisions
- Name: corvus (OldCrow tie-in). Namespace `corvus::`.
- Scope: statistical special functions only; basic transcendentals
  belong to Highway contrib. P0 (done): erf/erfc, erfinv/erfcinv,
  lgamma, incomplete gamma P/Q, incomplete beta. P1: digamma, inverse
  incomplete gamma/beta, Bessel I0/I1.
- Backend: Highway behind the `src/ops-inl.h` facade; public API
  std-only. std::simd migration = facade reimplementation, deferred
  until implementations mature (mid-2026: GCC 16 partial, no libc++).
- Dependency model: find_package(hwy) preferred, FetchContent fallback
  pinned to the audited version (1.4.0 — bump only with revalidation).
  Static lib, PIC on, Highway not exposed. MIT, clean-room only.
- Naming/extension/header/doc conventions and the CMake standard:
  AGENTS.md is canonical. Doc policy REVISED [2026-08-09, user]: the
  former four-file cap had concentrated growth into a 28 KB AGENTS.md
  — a recurring per-session context cost. Now: compact always-read
  AGENTS.md core (~4 KB: architecture map, non-negotiable traps,
  conventions, reading map) + on-demand references
  docs/ENVIRONMENT.md (fleet/build/CMake/CI) and
  docs/NUMERICAL-DOCTRINE.md (kernel hazards, test doctrine,
  generators, oracle doctrine, effort routing), carved by
  WHEN-needed, not topic. Same intent as the old cap (bounded context
  cost); resist growth in the core. Deferred deliberately: LTO/IPO
  (profile first), shared lib (no demand).
- **FP contraction off project-wide** (2026-07-25): GCC's default
  contraction fused inside log_dd's dd identities and shifted it 0.6
  bits vs MSVC — BETTER, which is why one compiler alone could never
  notice. The dd exactness proofs assume IEEE ops as written;
  CORVUS_FP_FLAGS sets -ffp-contract=off (clang-cl:
  /clang:-ffp-contract=off; MSVC /fp:precise already doesn't), PRIVATE
  to corvus and kernel-test targets. Fusion is requested in source
  (ops::MulAdd), never inferred. Cost ≤ ~8%, indicative.
- Windows compiler (2026-07-25, user): documented exception to the
  fleet MSVC default — validation/bench numbers from clang-cl
  (preferred, MSVC ABI) because HWY_BROKEN_MSVC caps MSVC at AVX2. MSVC
  stays a fully supported consumer toolchain and runs the CI Windows
  job. Full rationale in AGENTS.md.
- Platform tiers: Tier 1 (accuracy-audited on real silicon) = NEON
  (M1/CI), AVX-512+AVX2+SSE family (Ryzen native + capping), AVX2
  (Kaby). Tier 2 (compiles, unaudited) = SVE and anything else Highway
  emits.
- Phase order (2026-07-25): erfinv/erfcinv → incomplete gamma →
  incomplete beta. Rationale: smallest first unlocking probit; the
  inverse settles the erfc-tail question by measurement (it attenuated,
  κ ≈ 1/(2x²)); establishes the seed + dd-Newton pattern P1 inverses
  reuse; gamma's Temme machinery is beta's subset — build it with two
  arguments before three. All four predictions held.
- lgamma v1 scope (2026-07-21): full real axis, poles +inf, signgam
  deferred (SciPy's gammaln offers none either).
- **Oracle-trust doctrine** (2026-08-05, user): reference-oracle
  construction for a function WITHOUT a trusted library baseline is
  FRONTIER work (now in AGENTS.md routing — the beta oracle cost more
  sessions than its kernel). References are trusted only after an
  INDEPENDENTLY constructed verification harness passes clean with
  baked-in negative controls; tools/verify_beta_reference.py is the
  pattern. Exposure audit of the other oracles 2026-08-05: gamma was
  the only other one with custom logic — its a > 1e4 exact-asymptotic
  branch spot-checked 42/42 on 2026-08-06 via three Temme-independent
  evaluators (series / peak-normalized log-quad / high-depth Legendre
  CF, layered dps, overlap self-certification); the rest are thin
  wrappers over gold-standard single-argument mpmath calls.
- **First release** (2026-08-06, user): v0.1.0; Kaby non-gating;
  install/export status quo ratified — `cmake --install` requires a
  system Highway (find_package(hwy 1.4)), FetchContent builds are
  build-tree-only; no bundling, no nested install (keeps source tags
  free of the NOTICE obligation); revisit only if packaging starts.
  protect-tags ruleset makes v* tags immutable. Every future v-tag gets
  a Release object (page-coherence cadence).

## Shipped families — what would be expensive to re-derive
Full method and measured bounds: docs/ACCURACY.md. Math: kernel
derivation blocks at the definition sites. Full design texts: this
file's git history (pre-2026-08-06 compaction). Below: design points,
pinned parameters, and bugs worth not re-learning.

### erf + erfc [2026-07-20/21]
erf: table + local-Taylor (clean-room port of libstats vector_erf_neon
through the facade), max 1 ULP. erfc: core reuses the erf table via
compensated 1 ∓ erf assembly; tail e^{−a²}·G(1/a)/a on exp_dd, three
fitted intervals. Lessons that outlived the phase: the shared erf series
must meet erfc's RELATIVE precision near a = 6, not erf's absolute one;
the MulSub/no-FMA exact-residual hazard that created ops::SquareLow's
capability guard (AGENTS.md convention). The erfc normal tail's 2 ULP is
fit-limited (closed by decomposition 2026-07-25: not-CR rate is set by
the fit's double coefficients alone; a dd Horner buys 2→1 at a poor
speed trade — documented, accepted).

### Phase A — exp_dd + log_dd [2026-07-25]
exp_dd 2^-68.45 relative; log_dd 2^-67.88, correctly rounded on every
reference point; both bit-identical on all tiers and compilers.
- exp_dd returns mantissa + exponent (ExpDdFrac + ScaleTwo): consumers
  fold factors in first, power-of-two scaling lands LAST — that is what
  keeps subnormal results at one rounding. Every later family uses it.
- log_dd's mantissa is centred on 1; the slots adjacent to m = 1 carry
  R = 1, L = 0 exactly — no special case at x = 1.
- **An exact dd pair is not necessarily NORMALIZED**: (p−1, p_lo) had
  |lo| up to 2^-53 against hi ~ 2^-8; fix is one TwoSum (not Fast2Sum —
  r crosses zero). The gate caught it; reading the code did not.

### Phase B — lgamma [2026-07-25]
Five regions, exact-t zeros form, masked fixed-step DOWNWARD recurrence
(x − k exact; down also lands every lane in one polynomial), swept-X0
Stirling (X0 = 8), sinpi reflection via Γ(1−x) = −x·Γ(−x) (−x exact
where 1 − x is not; positive pipeline runs once). π/|sinπx| never
formed (overflows at poles) — sinc factoring cancels π analytically.
lgamma(1) needs an explicit +0.0 (t·B(t) carries −γ's sign at t = 0).
THE bug worth remembering — **ops::ProdLow's non-FMA Dekker split
overflows for |a| > 2^996** — is an AGENTS.md convention; Stirling's
x·(log x − 1) was the first operand there, and only the capped SSE
sweep caught it. Stirling grouped x·(log x − 1) − log(x)/2, never
(x − ½)·log x (product overflows a band below threshold). Generator
traps commented in gen_lgamma_data.py (mpmath fsum; odd node count
putting a node where 1 + t rounds to 1).

### Phase C part 1 — erfinv/erfcinv [2026-07-25]
C/T routing by Sterbenz-exact arguments; central direct dd-lead
polynomial (κ ≈ 1: a Newton would pass erf's error straight through);
tail seed + ONE dd Halley step (Halley beat Newton on seed degree
economics), log-space in the far tail (residual space underflows for
subnormal s), residual-space mid (needs ErfcDd exposure). Max 1 ULP
everywhere, every validated tier. Four bugs, all the shape "an identity
exact in general is wrong at one distinguished point, and only testing
that exact point catches it": ErfcCoreDd clamp at exactly 6.0 (seed can
land past 6; widened to 6 + 1/1024, provably no-op for erfc's public
results); (−0)+(+0) = +0 needed a trailing CopySign; `Neg(z − 1)` vs
`Sub(1, z)` differ only at z = 1 (−0); the reference oracle's own guess
blew up at s = 1 (log(−log 1)). Fourth defect caught in review, not
test: HalleyMid passed possibly-NaN x0 into a value-derived gather —
now the AGENTS.md masked-lane scrub rule. Bench (Ryzen, loaded,
indicative, scalar-walk baseline — an UPPER BOUND: the baseline pays
per-call dispatch/span overhead): central 22x/53x, mixed 7–12x.

### Phase C part 2 — incomplete gamma P/Q [2026-07-28]
Region map (λ = x/a, a_T = 20): R1 series-P fixed N=64; R2 CF-Q
BACKWARD fixed depth N=44 (backward contracts rounding — forward Lentz
was 24 ULP worst and false-converges); R3 Temme uniform asymptotic on
the ridge ratio band λ ∈ [½, 2] (clean-room c_k(η) extracted from the
oracle by Vandermonde solve in 1/a, disjoint-sample validated); R4
small-a Q-direct via Expm1Dd assembly (closes the small-a complement
corner with full relative accuracy). Small-side-direct routing with
1 ⊖ dd complements; prefactor E = a·LogDd(x) ⊖ x ⊖ LgammaPosDd(a);
shared primitives Log1pmxDd/Expm1Dd (hoisted to dd_special-inl.h
2026-07-29, byte-identity protocol), DdSqrt, DdRecipDd. Measured: 2 ULP
direct everywhere (one 4-ULP complement corner), identical cell-for-cell
on every fleet tier. What the stage gates caught (process record):
generator sup-proof self-checks found probe blind spots (N_cf 40→44,
R4 cap 30→36); the reference set caught a φ-coefficient double-rounding
worth 12 ULP through e^{−aφ} at aφ ≈ 740 (six leading φ coefficients now
dd); freeze masks must SELECT the accumulator, never add zero
(DdAddD renormalizes — lane-mix determinism test polices). Oracle traps
now institutional: compute the SMALL side directly in the oracle too;
mpmath's lower gammainc hangs/diverges for large a near the ridge —
exact-asymptotic oracle above a = 1e4 (spot-checked independently
2026-08-06, 42/42). Bench: Ryzen QUIET-machine at 7b52ed1: 49–74 ns/el
per region, 8.9–21× scalar walk (upper bound); Kaby numbers indicative
only.

### Incomplete beta P/Q [G1 2026-07-31 → shipped 2026-08-05/06]
The hardest family; eleven routing/assembly corrections and six oracle
defect classes. Kernel: src/beta-inl.h (own TU; consumes dd,
dd_special, exp_dd, log_dd, lgamma internals, erfc core, plus exactly
two gamma template cores for the gamma-limit slice). Everything below
is the surviving state; the correction-by-correction narrative is in
git history.

**Architecture.** Small-side-direct with 1 ⊖ dd complements; the
≤ 1 − 2⁻¹² doctrine applies to the EVALUATED side, orientation is
region-driven and decoupled from the direct/complement handout.
Routing order: (0) min ≤ ε_R4 tiny-first → R4 if inside its caps, else
fall through; (1) R1 power series if either orientation has ξ ≤ ξ₁ ∧
βξ ≤ B₁; (2) R3 Temme if ν ≥ T_ridge ∧ ratio band ξ/p, (1−ξ)/q ∈
[½, 2] (gamma's band shape; saturation a separate overlay); (3) R2 CF,
orientation by ξ < (α+1)/(c+2). Near-one R1 lanes (value > kBetaNearOne
= 1 − 2⁻¹¹, provably τ ≤ kBetaPrTauMax = 2.5) POST-ROUTE into R4's
analytic assembly in the fired orientation. Gamma-limit slice: max
param ≥ kBetaGammaLim = 2⁵⁹ routes R2-family lanes through gamma's
series/CF template cores with beta-side dd prefactor; R3's ridge floor
drops to kBetaGlRidgeMin = 20 in-band via the p→0-edge depth extension
kBetaR3GlExt (k = 10..12; extraction at p = 2⁻²⁰, NOT 2⁻⁵⁰ — the
ladder's c = ν/p must stay below the CF ground-truth ceiling ~2⁶¹).
Prefactor lnB never forms rounded c as an lgamma argument:
LgammaDiffDd(max, min) ⊖ lgamma(min), the latter via the
NINTH/TENTH-correction identities (lgamma(min) = LgammaDiffDd(1|2,·) −
ln(min)) for min ≤ 2.5 — a component-relative zone-poly budget is VOID
once the assembly cancels below it. Pinned constants (src/beta_data.h,
generator self-checked): N₁ = 64, N₂ = 64, T_ridge = 32, B₁ = 8,
ξ₁ = 0.45, ε_R4 = 2⁻⁶, R4 N = 48, Z₀ = 10, K_B = 16, C_lg = 256,
E_floor = −800, ζ_max = √(3·ln2/2), R3 K = 10 @ 25×15 + gl-ext
(29.88 KiB), η_γ = −ζ√2 gamma-limit mapping, e_k(ζ,p) = −e_k(−ζ,q),
kBetaTwoPiNuClamp = 2⁹⁰⁰ (Dekker ceiling; gamma's twin fixed too),
subnormal-τ rescale constants (2⁻⁹⁵⁰/2⁷⁰⁰ reframing + both-tiny
closed form Q̃ = r/(1+r)).

**Correction ledger** (each has a witness in git history): (1) G1a —
orientation decoupled from handout, R4 tiny-first caps added; (2) G1b —
tiny-min guard hoisted ABOVE R1 (rule order); (3) G1c — R3 membership
is the ratio band, not a cψ cap (mis-remembered gamma; no 32 KB fit
spans ζ ∈ [−5,5]); (4) G3 — R4 window cap max(ξ₁, thr_τ); (5) R1 λ ≥ 0
rule — REVERTED, too blunt; (6) near-one post-route → opposite CF —
destination failed its own check; (7) post-route destination = R4
assembly, fired orientation; (8) kBetaPrTauMax 1.5 → 2.5 + third
lgamma zone (the "safe band" claim was gamma-limit reasoning, wrong at
moderate β); (9) postroute lg1 via LgammaDiffDd identities (55/209 ULP
→ 1/1); (10) same disease at PA's lnB (R1 cmp 13 → 1); (11) non-FMA
subnormal-τ reframings + the 2¹⁰⁰⁰ → 2⁹⁰⁰ clamp (both only visible on
capped no-FMA tiers).

**Oracle/harness record.** Six oracle defect classes vs ZERO kernel
defects in adjudication — the origin of the oracle-trust doctrine.
Disease classes (all now guarded at enforcement sites in
gen_beta_reference.py): (a) truncation-at-ambient-dps — mpf ops
truncate higher-precision operands, mp.mpf() on an mpf RE-ROUNDS,
1 + τ = 1, a + b collapses in lnB, log1p(−near-1) = −inf; exact
complements carried, never recomputed; (b) component-relative error
budgets voided by downstream cancellation (Taylor branch deleted:
required accuracy is set by the RESULT's cancellation depth, unknowable
at the component site); (c) mp.quad returns √eps-scale noise when the
integral sits below working epsilon — normalize the integrand by peak
log-magnitude; (d) false saturation certificates from collapsed
complements; (e) mpmath betainc returns internally-CONSISTENT garbage
in the gammalim corner (layered-dps agreement cannot catch it — hard
excluded); (f) small-side-direct applies to the oracle too.
tools/verify_beta_reference.py is the INDEPENDENT harness (no oracle
import; layered series + half-split log-quad evaluator; exact analytic
lines; saturation log-bounds): 154 negative-control failures were
adjudicated two-of-three — 139 the harness's own quad, 4 the already-
fixed oracle rows (independent confirmation), 11 resolved by hand in
stored's favor. The 4 adjudicated rows are BAKED IN as a negative
control: every harness run must reject them and accept their
corrections before judging anything (exit 2 otherwise). Clean pass
2026-08-05 = the trust gate for the shipped references (37,099
rows/side + 251 specials, zero drops).

**Shipped numbers** (pinned gates, no margin, max over p/q; identical
gate cells on AVX3_ZEN4 native, AVX2/SSE4/SSSE3/SSE2 capped, Linux CI,
NEON CI run 31066128952): R1 1/1, R2 0/0, R3 3/1, R4 2/0, postroute
1/0, gammalim 0/1, specials exact. Monotonicity post-pass (3,278
(a,b)-groups, kernel dip slack 4 ULP) and ten 4001-point seam sweeps
(each must CROSS its boundary or fail): zero violations. Both live in
test_beta_ulp — no new binary, so tier coverage and the four-list rule
are automatic. BETA_ULP_DUMP env prints every not-CR row at/above a
threshold (permanent gate-pinning tool). bf16 lesson (2026-08-06 CI):
never default-construct vector-member structs — implicit ctors
instantiate outside the per-target attribute region; aggregate-init
(now everywhere).

**Bench** [2026-08-06, tests/bench_beta.cpp — Sonnet agent, reviewed]:
per-region point sets with a router-replica membership diagnostic
(printed, never gated; all six sets 100% in-region), scalar-walk
upper-bound baseline (bench_gamma pattern). Ryzen AVX3_ZEN4,
LOADED/INDICATIVE: R1 292, R2 448, R3 320, R4-tiny 278, postroute 543,
gammalim 598 ns/el (7.4–15.1× the walk). Quiet-machine re-run owed
before publishing (Next Steps).

**Process lessons** (tooling, not math): mp.dps must be set INSIDE
every computation layer — worker, subprocess-send, subprocess-parse
(bitten three times). Probe with exact hex from dumps, never
display-rounded values (phantom 2e-3 "errors" at d(lnI)/dx ~ 6800).
Windows multiprocessing: current_process().name, not parent_process();
probes need `if __name__ == '__main__'` guards or spawn bootstrapping
masquerades as fast failure. Sub-agent briefs: name the MECHANISM, not
the concept ("never set run_in_background on any tool call");
foreground only, chunk sweeps ≤ ~5-min re-runnable commands. Verify
exponent arithmetic by exponent SUM (8e-100·1e100 = 8, not 0.8 — a
false kernel-bug hypothesis cost an hour).

## P1 digamma — detail design [2026-08-06, frontier; BINDING]
Probe-validated (Sonnet probe agent, scratchpad digamma/p1–p6 scripts,
layered dps 60/100 clean; orchestrator-reviewed). [DERIVED, empirical]
unless noted.

**API**: digamma(x, out), full real axis. Specials (scipy 1.17.1
parity, probed): ψ(±0) = ∓inf (signed-zero pole convention); negative
integers → NaN (every double ≤ −2^53 is an integer → NaN); +inf →
+inf; −inf → NaN; NaN propagates; subnormal x → ∓inf (the −1/x term).

**Positive pipeline** (every shifted argument via exact TwoSum, never
a bare subtraction; x₀ dd: hi 0x1.762d86356be3fp+0, lo
0x1.b86a722197829p-54; trigamma(x₀) ≈ 0.9677):
- **Zone [1, 2)** — the full width-1 recurrence landing interval
  (probe: a narrower window cannot be a landing target under integer
  steps): product form ψ = t ⊗ P(t), t = dd shift vs x₀. Probe:
  all-double coefficients PLATEAU at 2^-53.85..2^-54.7 (degree 19–21)
  → dd LEADING coefficients required, the lgamma-zone pattern; degree
  and dd-lead count pinned by generator replay against 2^-55 relative,
  including 1e-14 neighborhoods of x₀.
- **(0, 1)** — one up-step WITHOUT forming 1+x (fl(1+x) would cost
  ~2^-52.5 relative near x → 1⁻): the same P evaluated at t₁ = dd
  shift vs (x₀−1) (x₀.hi − 1 exact; lo unchanged), then ψ =
  t₁⊗P(t₁) ⊖ DdRecip(x). Probe cancellation ratio < 1.74 (~1 bit) —
  dd absorbs trivially. −1/x dominates as x → 0 (no cancellation at
  the pole; ratio → 1).
- **[2, X0)** — masked fixed-step down-walk to [1,2): ψ = zone ⊕
  Σ_{j=1..m} 1/(x−j) in dd, x−j exact (Sterbenz, x < X0), weights
  DdRecipDd of exact pairs (gamma-R4 pattern), freeze-by-select.
- **[X0, ∞), X0 = 8** [probe table: K = 9 Bernoulli terms → 2^-56.5
  relative, the best margin of the sweep; X0 ≤ 5 cannot reach 2^-55
  at all, 6 is marginal]: ψ = LogDd(x) ⊖ DdRecip(2x) ⊖ x⁻²·S(x⁻²),
  dd head + double Horner tail, replay pins the dd-head count.
**Negative axis**: ψ(x) = ψ(1−x) − π·cot(πx), assembled in dd.
y_dd = TwoSum(1, −x) EXACT feeds the positive pipeline at dd argument
(lo correction via a rough-trigamma poly, ~2^-40 suffices — beta's
c.lo·ψ̃ pattern). cot from exact reduction u = x − round(x) (exact for
every double; π cancelled analytically via sinc-pair fits sin(πu)/πu
and cos(πu), lgamma-reflection pattern), ratio in dd. Probe HEADLINE:
at the nearest double to each of the first 20 negative zeros a plain-
double assembly loses 47.8–49.0 bits (3–4 correct bits — these are
legitimate ULP-sweep inputs, not edge cases). FIRST DESIGN CORRECTION
[2026-08-06, G1 escalation, chain depth 1]: the probe's "dd assembly
retains ~55+ bits" assumed IDEAL-dd components; with the design's own
2^-55-class fits the difference caps at fit-precision − cancellation
(replay: 5–20 relative bits at the adversarial doubles, components'
identity verified to ~249 bits in exact arithmetic — architecture
sound, metric wrong). Near-relative at the zeros would need
2^-104-class fits (zone degree ~40 all-dd + raised X0) — REJECTED on
cost for a contract lgamma's negative axis doesn't offer either.
Check (c) is the doctrine's own dual metric: (c1) ABSOLUTE ≤ 2^-56 at
the 20 adversarial doubles (replay measured 2.4e-18 ≈ 2^-58.5 —
margin held); (c2) RELATIVE ≤ 2^-52 at dense negative-axis samples
with |ψ| ≥ 1 (n = 1..20 intervals + log-spaced far intervals).
Accuracy doctrine (lgamma analog): relative where |ψ| ≥ 1, else
2^-53-class absolute; per-zero band width W ≈ target/|trigamma(z₀)|
(|trigamma| grows 8.9 → 18.9 over n = 1..20); pin to measured at G4.
Ratified G1 deviations: cos(πu) fit target is ABSOLUTE 2^-58 (cos has
a zero at u = ½; the cot assembly inherits the zero-band doctrine
there); rough-trigamma uses a floor-6 recurrence walk (~5 cheap
double steps — the KERNEL must mirror it); Chebyshev fits use
per-degree matched-node fitting (truncating a high-degree fit cancels
catastrophically in the monomial conversion). NO reflection domain ceiling: u is exact everywhere and
|x| ≥ 2^53 negatives are all integers → NaN (scipy's NaN at
−1e300+0.5 is input rounding to an integer, not a formula limit).
**Seeds**: DigammaRough/kBetaDigammaCoef NOT reusable — plain value
fit whose measured floor (2^-43.8) sits exactly at x₀; prior art only;
stays in beta_data.h documented unused.
**Oracle**: mpmath.digamma — gold-standard single-argument thin
wrapper; standard generator/reference pattern, no independent harness
(inside the oracle-trust doctrine's trusted-baseline scope).
Reference set: per-region grids; x₀ bit-neighborhoods (offsets to
1e-14 both sides); the 20 adversarial nearest-double negative-zero
points ± offsets; pole neighborhoods −n ± ulp-scale; zone/X0 boundary
brackets; subnormals; huge x; the specials table.
**Targets** [ILLUSTRATIVE until measured; pin to measured, no
margin]: positive axis ≤ 1 ULP relative; negative axis ≤ 1–2 ULP
where |ψ| ≥ 1 and 2^-53-class absolute inside the zero bands.
**Kernel/TU**: src/digamma-inl.h + digamma.cpp, own TU (consumes dd,
log_dd, ops; no gamma/beta cores); HWY_NOINLINE day-one on cores AND
driver; test_digamma_ulp + smoke registered in DEPENDENCY position in
all FOUR lists.
**Stage record**: G1 SHIPPED [2026-08-06, 9203b1f] — pinned by
replay: zone degree 21 / 2 dd leads (2.44e-17), asymptotic K = 9 /
1 dd head (9.82e-18), sinc 8/3 rel + cos 9/3 abs at 2^-58,
rough-trigamma K = 8 floor-6 walk (3.47e-13 vs 2^-40); reflection
(c1) 4.83e-18 abs worst at the 20 adversarial doubles, (c2) 7.19e-18
rel worst over 855 |ψ| ≥ 1 samples. One escalation (depth 1 → FIRST
correction); three self-caught tooling bugs (spurious π² in the cot
ratio, Chebyshev truncate-from-high-degree cancellation → per-degree
matched-node fitting, array-emission double-subscript). CI green.
G2 SHIPPED [2026-08-06, 8f23e99] — 15,709 rows, lgamma two-hex-double
format; every row layered-dps 60/100; root + 20 negative zeros
independently recomputed (6.8e-67 worst layer disagreement);
independent hand-derived 25-row spot rederivation worst 6.0e-17; 139
oracle-overflow points excluded (subnormal-x → ∓inf is smoke
doctrine); walk-step brackets added (325 rows). Zero escalations.
G3 SHIPPED [2026-08-07, a931228, Opus agent, zero escalations] —
kernel + smoke + ULP + bench + four-list registration. Measured
(AVX3_ZEN4 native AND SSE2 capped, IDENTICAL gate cells): 1 ULP max
in all five relative buckets ((0,1)/zone/walk/asym/neg-|ψ|≥1),
1.00×2⁻⁵³ absolute in the negative zero bands; all 15 prior gates
byte-identical. Three deviations reviewed and ACCEPTED: (i) large-x
cut kDigammaAsymCut = 2⁸⁵ (brief's 2⁵⁵-class sketch would flip
rounding at 2⁻⁶¹-relative dropped terms; derivation at the site),
(ii) 2⁻⁹⁶⁰ direct-reciprocal shortcuts on (0,1) and the cot ratio
(quotient-side Dekker ceiling; delivers the subnormal ∓inf doctrine —
infinities cannot flow through TwoSum), (iii) w = (1/x)² never 1/x²
(x·x overflows past 1.3e154). Bench loaded/indicative: 7.0–27.4
ns/el, 7.2–25.2× scalar walk.
G4 COMPLETE [2026-08-08]: gates PINNED to measured, no margin (1 ULP
in all five relative buckets, 1.0 × 2⁻⁵³ absolute band). Full ladder:
AVX3_ZEN4 native; AVX2/SSE4/SSSE3/SSE2 capped (clang-cl sweep, all
tiers passed); Linux CI 4-tier sweep + sanitizers; NEON (CI run
31231085106 — table identical to native INCLUDING not-CR counts and
worst-x points); Windows MSVC. The near-pole shortcut-coverage watch
item CLOSED by analysis: |u| < 2⁻⁹⁶⁰ is reachable only through
x ∈ (−2⁻⁹⁶⁰, 0) — near a pole −n the smallest representable |u| is
ulp-scale ≥ 2⁻⁵³ — and that band's doctrine answer (±inf) is
smoke-gated; no reference rows needed.
G5 COMPLETE [2026-08-08]: ACCURACY.md matrix row (dagger extended to
digamma) + full family section; README status/bullet/example — same
change set as the gate pinning. **DIGAMMA SHIPPED** — the first P1
family, via the agent pipeline: probe (Sonnet) → G1 (Sonnet, one
escalation → FIRST correction) → G2 (Sonnet, zero) → G3 (Opus, zero,
three reviewed-and-accepted deviations) → G4/G5 (orchestrator).

## P1 trigamma — detail design [2026-08-08, frontier; BINDING]
Probe-validated (scratchpad trigamma/p1–p5, layered dps 60/100 clean).
Digamma-shaped MINUS the hard parts; ALL-RELATIVE everywhere.

**API**: trigamma(x, out), full real axis. ψ₁(x) = Σ 1/(x+n)² is a sum
of squares — positive wherever finite (the no-zeros PROOF), so no
absolute band exists anywhere; single relative metric. Negative-axis
global min 8.933 at x ≈ −0.4957; per-interval minima rise
monotonically to π²; worst reflection cancellation ratio 1.107
(≈ 0.15 bit) at x ≈ −0.455.
**Specials (scipy parity, probed — NOT digamma's convention)**: every
pole is a DOUBLE pole, sign-unambiguous → **+inf** at ±0, at every
negative integer, and at every negative double |x| ≥ 2^53 (all
integers); +inf → +0; −inf → +inf (scipy); NaN propagates;
subnormals of both signs → +inf.
**Positive pipeline** (shifts via exact TwoSum, digamma mechanisms):
- Zone [1,2): PLAIN value fit (no product form — no zero). FIRST
  CORRECTION [2026-08-08, G1 escalation, depth 1]: the probe's
  "degree 24 / 1 dd-lead / 2^-55.16" was a GRID ARTIFACT — the true
  worst points sit within ~1e-10 of the interval edges (coherent
  Chebyshev coefficient-rounding), and edge-refined bit-stepped
  sampling shows 1 lead plateaus at 2^-53.7..54.0, 2 leads at
  2^-54.9. PINNED: **degree 27, 3 dd-leads** (2^-56.5, >1 bit
  margin). BINDING RULE from the root cause: every replay self-check
  in this family uses edge-refined bit-stepped boundary sampling —
  the gamma probe1 / beta R3-lens disease at its third occurrence.
- (0,1): up-step without forming 1+x (zone at shifted centre) ⊕
  dd(1/x²). Covers tiny x naturally — the zone term → ψ₁(1) = π²/6,
  the Laurent constant (probe REFUTED the 1/x term: coefficient is
  exactly 0; π²/6 stops mattering below ~2^-28). Deep-tiny guard:
  below ~2^-480 the dd 1/x² alone is the result (zone < 2^-950
  relative), overflowing to +inf at x ≤ 2^-512 = 1/√DBL_MAX. PROBE
  WARNING for G3: naive double (1/x)² or 1/(x·x) is NOT reliably CR
  (24–46% 1-ULP misses) — the reciprocal-square stays dd end-to-end
  down to the overflow boundary; if Dekker-split limbs land subnormal
  in the deep-tiny lane, use exact power-of-two rescaling (beta
  ELEVENTH pattern). ESCALATE if dd cannot reach the gate there.
- [2,8): down-walk ≤ 6 exact steps subtracting dd 1/(x−j)²
  (x−j exact), freeze-by-select.
- [8, 2^89): Bernoulli asymptotic in the DIRECT (unfactored) sum form
  (probe: conditions better than 1/x-factored), K = 11, dd-head count
  by generator replay. NOTE: log-free — trigamma does not consume
  log_dd at all.
- x ≥ kTrigammaAsymCut = 2^89 (conservative analytic cut per the
  digamma doctrine; empirical crossover ~2^55): fl(1/x) alone
  (dropped part < 2^-90 relative); retires every large-operand dd op
  below the non-FMA Dekker ceiling.
**Negative axis**: ψ₁(x) = π²/sin²(πx) − ψ₁(1−x) in dd: exact
u = x − round(x); π²/sin² = 1/(u·sinc(u))² from the sinc fit
(re-emitted into trigamma's own header by the same fit procedure —
~10 duplicated constants beat invoking the hoist/byte-identity
protocol on shipped digamma); y_dd = TwoSum(1, −x) exact, lo
correction y.lo·ψ₂(y.hi) via a CRUDE tetragamma (bound analysis:
whole correction ≤ ~2^-55.9 relative because ψ₁ ≥ 8.93 — a ~2^-30
floor-walk asymptotic fit is ample; pattern-identical to digamma's
rough-trigamma, far looser target). Probe (c): simulated dd assembly
worst 2^-54.7 with 2^-55 components — all-relative gating safe.
**Oracle**: mpmath trigamma/polygamma(1,·) — trusted single-argument
baseline; layered dps on every row; NO adversarial-zero stratum
(no zeros exist). Independent spot rederivation: direct Σ 1/(x+n)²
+ Euler–Maclaurin tail (trivially mpmath-digamma-independent).
**Reference strata**: region grids ((0,1) log incl. the 2^-512
overflow boundary, the ~2^-480 guard, and π²/6-crossover ~2^-28
brackets; zone dense; walk + step brackets; asym log to the 2^89 cut
both sides and on to 1e308); negative: dense (−50, 0), global-min
neighborhood, near-pole ulp-offset brackets (n = 1..20, 100, 1e3,
1e6-class), far log-spaced to 2^52. Specials excluded (smoke).
**Targets** [ILLUSTRATIVE until measured; pin at G4, no margin]:
≤ 1 ULP relative everywhere, single metric.
**Kernel/TU**: src/trigamma-inl.h + trigamma.cpp, own TU (consumes dd
+ ops ONLY); HWY_NOINLINE day-one; tests at the END of all four
lists; bench_trigamma per-region.
**Process**: combined G1+G2 in ONE Sonnet tooling agent (generator +
reference set, both self-check families, single review gate —
justified by family simplicity); G3 Opus; G4/G5 orchestrator.
Escalation-density rule applies.
**Stage record — TRIGAMMA SHIPPED [2026-08-08]**: G1+G2 (a3bcbf6,
combined Sonnet agent): zone 27/3 dd-leads (9.6e-18, edge-refined),
asym K=11 head B₂ dd (1.4e-17), sinc 8/3 bit-identical to digamma's,
crude tetragamma floor-6 K=6, deep-tiny guard 2⁻⁴⁸⁰ derived,
reflection replay 2⁻⁵⁹·², 14,928 rows all layered-dps, oracle
fast-path via reflection for |x| > 50 (mpmath polygamma(1,·) is
O(|x|) on the negative axis — verified 1e-99), 25-row independent
spot check. Ledger: one design escalation (grid artifact → FIRST
correction + the edge-refined-sampling binding rule), one SECOND
correction self-diagnosed (derived recurrence replay target,
ψ₁(1)/ψ₁(8) ≈ 12.4× amplification), one process fault (parked on a
Monitor — recovered by resume; memory updated to name every
background door). G3 (eb0a556, Opus, zero escalations, seven
reviewed-accepted deviations — the standout: deep-tiny via exact
2⁵¹² rescale whose lower clamp itself delivers +inf, one rounding,
one code path). G4/G5: gate PINNED at 1 ULP single relative metric
(a hardcoded "gate 8" display string fixed — enforcement was always
kMaxUlp); ladder identical cells everywhere — AVX3_ZEN4 native,
AVX2/SSE4/SSSE3/SSE2 capped, Linux CI sweep, NEON (identical incl.
not-CR counts), Windows MSVC (12m11s — trigamma.cpp needs NO /d2
flag); (0,1) bucket correctly rounded; walk amplification cost no
bit (shows as not-CR 2.97%/4.44% FMA/non-FMA). ACCURACY.md + README
in the change set. Bench indicative: 5.7–27 ns/el.

## P1 inverse incomplete gamma — detail design [2026-08-08, frontier; BINDING]
Probe-validated (Sonnet probe agent, scratchpad gammainv/p1–p6 + common.py,
layered dps 60/100; orchestrator-reviewed). [DERIVED, empirical] unless noted.

**API**: gamma_p_inv(a, p, out) → x with P(a,x) = p; gamma_q_inv(a, q, out)
→ x with Q(a,x) = q. Full [0,1] contract on BOTH sides: input s > 1/2 flips
to the complement side via 1 − s, EXACT by Sterbenz for s ≥ 1/2 — the
inverse's complement transform is on the INPUT and costs nothing (cleaner
than the forward pair, whose complement rounding is on the output). One
shared core pipeline, two exports (erfinv/erfcinv TU pattern).
**Specials (scipy parity, probed P1c)**: p=0 → 0, p=1 → +inf (q mirrored:
0 → +inf, 1 → 0); s outside [0,1] → NaN; a ≤ 0, a = +inf, NaN → NaN.
**Conditioning adjudication (frontier review of probe P1 — the probe's two
"severe collapse" findings DISSOLVE, neither weakens the contract):**
- Tiny-a "zero-bit collapse" (κ_p = 1/a unbounded): lies ENTIRELY inside
  the output-underflow region — P1b measured that for a ≤ ~9.3e-4 the whole
  small-p side maps below DBL_MIN_NORMAL, and beyond-round-to-zero for most
  of it (CR answer 0, exact). The probe's high-κ tiny-a interior points all
  have their SMALL side on q with κ ≈ O(1). With the exact input-side flip,
  every input whose true x is a normal double has κ ≤ ~2^10.1 (κ = 2^10
  contour at a = 2^-10, self-limited: |ln p|/a ≤ 745 wherever x is
  representable).
- Huge-a "non-injectivity" (a ≳ 3e34: whole transition < 1 ULP of x): κ → 0
  there (measured 5e-91 at a = 1e90) — the Temme seed alone is CR-class and
  Newton steps SELF-FREEZE (Δ/x ~ κ). Test stratum, not a branch.
**Architecture** (erfinv seed+dd-step precedent; forward-core reuse per
beta's gamma-limit precedent — the inverse TU assembles prefactor ⊗ region
core in dd itself, mirroring GammaVec's internal assembly unrounded):
1. Side selection: solve against small side s ≤ 1/2 (exact flip above).
2. SEED (double precision, per-region; parameters replay-pinned at G1 with
   edge-refined bit-stepped sampling):
   - S1 (a ≥ a_T ≈ 20): Temme normal-quantile — z = erfcinv(2s) (sign by
     side), η₀ = -z·√(2/a), invert ½η² = λ − 1 − ln λ (series near η = 0,
     Newton elsewhere; scheme pinned at G1), x₀ = a·λ(η₀), ε_k(η)/a^k
     corrections (count by replay; probe's S1 one-Newton anomalies at
     a ~ 100 were its own correction formula — G1 re-derives from the
     published Temme 1992 expansion, clean-room paper math).
   - S2 (p-side, a < a_T): x₀ = exp((ln p + lnΓ(1+a))/a) + Picard
     corrections (count by replay; probe: wins everywhere on p-side at
     a < 1, 60 bits at small p).
   - S3 (q-side, a < a_T): L = -ln(q·Γ(a)) fixed-point iterations (count
     by replay).
   - Weak-seed middle band (a ∈ [~0.1, a_T), s near 1/2 — probe S4: best
     seed ~4–6 bits): curvature is benign there — covered by step count 3;
     if G1 replay shows 3 steps insufficient anywhere, ESCALATE (a
     dedicated 2D fit is a design change).
3. STEPS (dd residual against forward template cores GammaSeriesSum /
   GammaCfRecip / GammaTemme / GammaSmallQ + prefactor e^E machinery;
   routing by (a, x_seed) mirrors the forward region map;
   freeze-by-select; per-region count pinned by replay, max 3):
   - Plain dd Newton Δ = (P_dd(x) ⊖ p)/g; g from the forward prefactor
     (dP/dx = e^{E}-class, dQ/dx = −g).
   - Log-residual Newton in the far q-tail — MANDATORY, not an
     optimization (P3: 65 bits from a 14-bit seed vs 30 plain): Δ =
     (ln Q_dd ⊖ ln q)·Q/g in dd; also the recovery tool on the ridge
     curvature band (a ≫ 1, λ within O(1/√a) of 1) if replay wants it.
   - Halley (analytic g'/g = (a−1)/x − 1) available; G3 may trade 1
     Halley vs 2 Newton on bench, gates must hold either way.
   Probe P3/P4: 2 steps from a ≥ 20-bit seed reach the internal-dd floor
   (57–69 bits) in every interior; internal budgets 2^-56..2^-58 (fit-
   limited, not dd-representation-limited); R3 ridge does NOT inherit the
   forward's external 2 ULP.
4. DEEP-SMALL closed form (p-side, when a·x₀ < 2^-60 — the correction
   series is dead): x = exp_dd((LogDd(p) ⊕ dd lnΓ(1+a)) ⊘ a), mantissa +
   exponent scaling LAST = one rounding into subnormals/zero. Owns the
   entire tiny-a collapse zone and the subnormal-x band. Amplification
   argument: rel-x error = |ΔS|/a ≤ 745·rel_S (self-limited by the exp
   underflow range |S/a| ≤ 745), so dd S ⇒ ≤ 2^-90-class — CR throughout.
**Targets** [ILLUSTRATIVE until measured; pin at G4, no margin]: 1–2 ULP
relative, both sides, full domain including subnormal outputs. scipy
baseline (probe P1c vs prototype oracle): median 3 ULP, p99 341, max 793.
**Oracle (G2; frontier-specified, oracle-trust doctrine — no library
inverse exists anywhere)**: per-row BRACKET CERTIFICATION at layered dps
60/100 (dps 30 measured under-certifying 2/20 — never lower): root-find
x* seeded by S1/S2/S3, round to xd, then certify sign(P − p) flips across
the two half-ulp midpoints of xd as exact mpf, forward evals via mpmath
gammainc with the a > 1e4 exact-asymptotic branch reused from
gen_gamma_reference.py. Deep-small rows: NO root-find — closed-form
log-space oracle (ln x = S/a in mpf), certified in log space against the
subnormal/zero boundary midpoints. Huge-a rows: exact-asymptotic oracle
as the independent second route (probe's elementary-series cross-check
cannot converge there — expected, not a defect). Negative controls baked
in (beta doctrine): known-bad rows the certifier must reject, exit 2
otherwise. Measured cost 29 ms/row median → 15–40k rows ≈ 10–90 min.
**Reference strata (G2)**: per-region (a, s, side) grids; underflow
p-threshold brackets (P1b table); subnormal-x band; deep p-tail to
subnormal-min p; far q-tail to subnormal-min q; ridge band λ ∈ 1 ±
O(1/√a); huge-a {1e16 … 1.7e308} incl. a·φ saturation edges; a_T
brackets; weak-seed middle band dense; specials excluded (smoke).
**Kernel/TU**: src/gammainv-inl.h + gammainv.cpp, both exports one TU
(shared cores); consumes gamma-inl.h template cores + dd/dd_special +
exp_dd/log_dd + lgamma internals. HWY_NOINLINE day-one on cores AND
driver; /d2ReducedOptimizeHugeFunctions on gammainv.cpp from day one
(real MSVC only — the TU instantiates the heavy gamma cores twice per
export; gamma.cpp precedent). Tests at the END of all FOUR lists.
**Process**: G1 (Sonnet, generator gen_gammainv_data.py →
src/gammainv_data.h) and G2 (Sonnet, oracle + references, SEPARATE agent
— the oracle is the risk item) → G3 (Opus kernel) → G4/G5 orchestrator.
Escalation-density rule applies. Beta inverse follows as its own
pipeline after this one ships (P6 scoping: swap identity
I_x(a,b) = 1 − I_{1−x}(b,a) gives lossless near-1 output via argument
swap; tiny-a AND tiny-b are independent collapse triggers — needs its
own probe).
**FIRST CORRECTION [2026-08-08, G1 escalation, chain depth 1]**: the
seed partition is by (side, λ-regime) at ALL a, not by a alone — the
probe's "S2 wins p-side / S3 wins q-side" was a λ-regime truth tested
only at a < 1. G1's replay caught the two corners the a-gated partition
leaves uncovered: (i) deep p-tail at a ≥ a_T (a=20, λ=0.02: S1's
η ≈ −2.4 weak tail, 47.98 bits at 3 steps — while S2's Picard
contraction x/(a+1) ≈ 0.02 there seeds ~17 bits and converges easily);
(ii) small-a mid band (best seed 4–6 bits; Halley topped at 54.5,
half a bit under margin — fix the seed, don't shave margin). RATIFIED:
tri-candidate seed {S1 if |η| in domain, S2 p-form (usable from either
input side via the exact complement), S3 under its stability gate}
selected per lane by cheap forward-residual comparison (G1's own
mechanism, now global); S2 Picard count re-pinned by replay; a_T
governs only the central/ridge band. Ratified G1 deviations: S1
corrections K=2 (c₂ Vandermonde extraction unstable at sane node/dps
budgets, marginal seed-bit gain); S3 stability gate L > 3·|a−1| (the
design's L > 0 was necessary, not sufficient — contraction factor is
(a−1)/x); step variant is per-depth-bucket (log-residual deep, plain
shallow), not global.
**SECOND CORRECTION [2026-08-08, G1 re-escalation, chain depth 2]**:
the replay's UNIFORM internal-dd noise model (2^-56/2^-58) is wrong at
the shallow small-x class (a ≈ 0.1–0.3, moderate s, x tiny) — those
budgets are region-worst fit/series-LENGTH bounds that bind near the
forward region's far boundary (x ≈ a+1, slowest convergence), while at
tiny x the R1/R4 series is super-converged and the true per-point
forward error is component-limited via the prefactor's ABSOLUTE-E error
(LogDd 2^-67.9, lgamma dd core 2^-68-class → ~2^-66 at |E| ≈ 3.5).
Diagnostic confirmation: the measured 54.68-bit floor is EXACTLY
58 − log2(1/a) at a = 0.1 — pure eps·κ, seeds and steps already
optimal. RATIFIED: per-point analytic eps in the replay at that class —
max(series-tail bound at N, dd-accumulation, prefactor-component bound,
safety floor 2^-64) — uniform 2^-56/2^-58 retained everywhere fits and
lengths genuinely bind. Gate stays ≥ 55 bits, unshaved; predicted floor
at the class ≈ 60+ bits. G3's silicon gates remain the arbiter — the
replay is design-sanity, not proof. A further escalation hits chain
depth 3 — takeover assessed under the escalation-density rule, which
is a JUDGMENT call, not a hard ceiling [2026-08-08, user refinement]:
take over on churn, or when sequential local fixes signal a global
flaw a larger-context pass should assess whole; clean,
precisely-diagnosed chains (as both of this stage's have been) may
continue delegated.
**Stage record**: PROBE COMPLETE [2026-08-08, Sonnet]: 6 self-caught
tooling bugs (worst: native-float bisection bounds silently capping the
oracle at double precision — caught by bracket-certification failures),
0 design escalations; both "severe" findings adjudicated at frontier
review as metric-framing artifacts (see Conditioning above); P2's
largest-a S1 rows predate the bisection fix — indicative only, G1 replay
re-measures. G1 SHIPPED [2026-08-09, Sonnet + frontier takeover at
chain depth 3]: gen_gammainv_data.py + src/gammainv_data.h. Sonnet
carried both ratified corrections and 14 self-caught bugs (worst:
fixed-dps erfcinv rounding 1−y to 1 for y < 1e-30 corrupting every
deep-tail S1 seed; the probe's S1 correction-formula side-sign bug —
the a ~ 100 anomaly's root cause; the direct-region-vs-small-
probability side confusion in the replay harness). TAKEOVER (judgment
per the refined rule — not the count): the last two fixes changed
measurement infrastructure and the numbers moved in ways the agent
could not characterize. Frontier root cause, one line: the main
replay loop solved against the UNROUNDED mpf forward value while
being measured against oracle_x's root of the ROUNDED double s — a
κ·2⁻⁵⁴ basis mismatch (50.95 bits at κ ≈ 2^3.3, exactly as
predicted). Second frontier fix: the per-point eps model extended to
the q-side twins of the SECOND-correction class (same components,
P→Q ratio conversion max(1, (1−s)/s), floor stays absolute 2⁻⁶⁴) —
those q-side points (κ ≈ 9–10) were the residual 54.68 floor. FINAL
PINS: a_T = 20 (= kGammaAT, referenced not duplicated); λ(η) series
order 12 |η| < 1/2 + 6-iter log-space Newton to |η| ≤ 9.24; S1 K=2
Chebyshev 2×25 wide-η; S2 Picard 6; S3 fixed-point 3 gated
L > 3|a−1|; S1_A_MIN = 0.3; tri-candidate residual-compared seed at
all a < a_T; steps: nsteps = 3 shared, LOG-RESIDUAL NEWTON BOTH
BUCKETS (the plain-wins-shallow read was a basis-mismatch artifact;
plain leaves shallow worst cases at 20–30 bits — G3 gets ONE step
variant, a kernel simplification); deep-small cut 2⁻⁶⁰; replay
floors deep 59.68 / shallow 57.54 vs 55-bit gate, per-point eps
audited (198 uses, floor-term wins 171, max 2⁻⁵⁹·¹).
G2 SHIPPED [2026-08-09, Sonnet, zero escalations]:
gen_gammainv_reference.py + 14,926 rows (p 8,406 / q 6,520,
three-hex-double a s xd), every row bracket-certified at layered dps
60→100; deep-small rows via the log-space closed form with the
dropped-correction bound folded into certification slack; a ≥ 1e16
rows dual-route (two independently-anchored Temme fits, kext 15/13,
distinct nodes+anchors, overlap-validated 1.7e-35 / 5.8e-30); 4
negative controls REJECTED on every run (they caught the certifier's
own first-draft deep-small hint bypass — the doctrine paying for
itself); 40/14,966 boundary-ladder points correctly declined rather
than guessed. KEY FINDING (test design input for G3/G4): the
beyond-resolution collapse is domain-wide for a ≳ 3e34 — 1 ULP off
x = a already saturates a·φ past 800 — so ~half the random-grid rows
certify xd = a exactly; correct but trivial. The ULP test MUST
bucket beyond-resolution rows separately so their zeros don't dilute
real-region statistics; the resolvable domain carries ~7.0k rows
across all named strata. Orchestrator review: 12/12 independent
mpmath bracket spot-check (direct-side evaluation — the reviewer's
own first check hit complement collapse, disease class (a), fixed).
Ratified deviations: no in-file header (matches all 13 existing
reference files + the raw-tokenizer reader pattern); route 2 as an
independently-anchored Temme fit (no unrelated trusted method exists
at that scale).
G3 SHIPPED [2026-08-09, Opus, zero escalations, NINE reviewed-and-
accepted deviations — the standout family of the arc]: kernel + smoke
+ ULP + bench + four-list registration. Measured 1 ULP max in every
bucket, deep-small/subnormal/x=0 correctly rounded, IDENTICAL cells
(incl. not-CR and worst points) on clang-cl AVX3_ZEN4 native, g++
SSE2-capped, and MSVC AVX2 — three toolchains, FMA and non-FMA.
Accepted deviations, each measured against the pinned alternative:
(1) deep-small cut is x₀(1+a) < 2⁻⁶⁰, NOT a·x₀ — dropped-factor
error ~x/(1+a) is a-independent; the pinned form is 1/a looser below
a=1 (90 ULP reachable at a=1e-4 via the q orientation; G1
self-check (g) swept only the p orientation — [OPEN] annotate the
generator); (2) Newton objective is the LOGIT m = lnP − lnQ (solved-
side log saturates: 2e18 ULP at a=1.9e34; signed ln min(P,Q) jumps
2ln2 at the median: 5e14 ULP; logit is continuous, saturates
nowhere, first-order-identical step, nsteps=3 holds); (3)
safeguarded Newton (reject residual-increasing steps, 1/8 backtrack,
bypass when |resid| < 1/2) — the design's "steps self-freeze" fails
in the collapse zone where lnF is locally quadratic (6626 ULP
measured), and the bypass is needed the other way (99 ULP);
(4) additive-x step (log-x converges slower: 12 ULP vs CR), −0.9
relative-step floor; (5) E dual-form split at a_T (direct form's
terms ~7e302 at a=1e300; a·ln x overflows past ~2.5e305), Stirling
μ from kLgammaStirCoef; (6) forward returns logit+slope, NO
saturation clamp — keeps the whole underflow range live; (7)
tri-candidate seeds at ALL a per FIRST correction (the G3 brief's
"S1 alone ≥ a_T" was the brief's error); (8) lnΓ(1+a) via lgamma
zone poly at exact shifted arg for a ≤ 3/2; (9) non-finite
candidates score +inf. Six self-caught bugs (worst: a=DBL_MAX
candidate rejected by Lt(x, DBL_MAX) → uninitialized-forward step to
5.5e307). Bench indicative: 43 ns/el deep-small fast path, 520–1056
elsewhere, 6.3–12.3× scalar walk (cost is structural: 3 seeds + 6
forward evals). MSVC WATCH: gammainv.cpp 2.0 min with /d2, library
9.3 min on the Ryzen box — largest single Windows-build jump; CI
timeout 25 min.
G4/G5 COMPLETE — **GAMMA_P_INV / GAMMA_Q_INV SHIPPED [2026-08-09]**:
gate PINNED to measured, no margin (1 ULP every bucket;
deep-small/subnormal/x=0 bands CR). Full ladder asserting: AVX3_ZEN4
native; AVX2/SSE4/SSSE3/SSE2 capped clang-cl sweep (all tiers
passed, cells identical); Linux CI sweep + sanitizers; NEON (CI run
31324636938); Windows MSVC (13.7 min — watch item stands, timeout
25). ACCURACY.md matrix rows + family section, README, in the gate-
pinning change set. POSTSCRIPT (fifth family to teach a lesson at
the last leg): the first pinned-gate CI run failed on NEON — G3 had
default-constructed GammaInvFwdOut, the 2026-08-06 bf16 lesson with
its site comment INVERTED (the agent cited the precedent but
memorized the fix backwards; only the macOS job materializes the
attribute mismatch). Fixed 88a980b, aggregate-init. Also: a
background CI watch attached to the previous still-running run and
reported a false green — verify the watched run's SHA, always.
Next session opens the BETA inverse (own probe → design pipeline;
P6 scoping notes in the probe record).

## P1 inverse incomplete beta — detail design [2026-08-09, frontier; BINDING]
Probe COMPLETE [2026-08-09, Sonnet]: 5 self-caught tooling bugs, 0 design
escalations; 4 flagged open questions, ALL adjudicated here. The full
record was scratchpad-only; everything binding survives in this section.
**API**: beta_p_inv(a, b, p, x) solves I_x(a,b) = p; beta_q_inv(a, b, q,
x) solves 1 − I_x(a,b) = q. Spans, one TU two exports. The swap identity
I_x(a,b) = 1 − I_{1−x}(b,a) is the documented lossless-near-1 mechanism:
1 − x at full relative precision = beta_p_inv(b, a, q). Doxygen states it.
**Conditioning adjudications (probe B-P1; one dissolves, one is REAL,
one is a stratum)**:
- Single-tiny parameter: DISSOLVES per gamma's precedent with the
  boundary GENERALIZED — κ = 1/a exactly in the power-law regime,
  self-limiting boundary a*(b) ≈ 1/(1074 − log2(b)); the gamma-limit
  corner (b → 1e300) WIDENS the collapse zone 14× (a* ≈ 1.3e-2 there).
  Input-side flip + deep-small closed form own it, as in gamma.
- JOINT-tiny plateau (both a, b tiny) — REAL, does not dissolve
  [measured]: interior density f(1/2) ≈ 4·min(a,b) ⇒ κ ~ 1/min(a,b) at
  interior REPRESENTABLE x where NEITHER probability side is small
  (plateau value s* = b/(a+b) is interior). dd (2^-105-class) resolves
  y to 1 ULP only for κ ≤ 2^52 (min(a,b) ≳ 2^-52; measured threshold
  1.1e-16). ADJUDICATION — dedicated joint-tiny route (S4 below) plus a
  CONTRACT SPLIT: rows with κ ≤ 2^52 stay under the y-ULP gate; rows
  above carry a BACKWARD-ERROR contract (forward value of the returned
  y within ~2 ulp of s — the statistically meaningful guarantee: the
  returned quantile inverts a probability indistinguishable at double
  precision). Metric-band precedent: lgamma/digamma absolute bands.
  G2 computes κ per row and buckets. The route achieves the
  information limit of dd precision in ONE evaluation — no iteration
  scheme at dd precision can beat κ·2^-105, so the split is honest.
- Huge-ν beyond-resolution: whole transition < 1 ulp of x once the
  SHAPE-side parameter reaches ~1e33–7e34 at every skew tested — the
  beta→gamma limit reproducing gamma's own ~3e34 threshold. Test
  stratum, not a branch (gamma precedent); bucketed separately.
**Architecture** (gammainv is the pattern; its G3-proven mechanisms are
carried as house doctrine, not re-derived):
1. Input side: solve against s = min(p, 1−p), exact Sterbenz flip.
2. Output orientation: per-lane swap so the solved variable y is the
   end x is near — I_y(α,β) = target with (α,β,side) relabeled by the
   swap identity, σ carried EXACT through both flips (the swap
   re-labels p↔q; no complement is ever recomputed). The logit
   objective is antisymmetric under both flips — one code path.
3. SEEDS — quad-candidate, cheap forward-residual global selection at
   ALL (a,b,s) (the FIRST-correction mechanism, now standard; probe
   measured every candidate everywhere):
   - S1 beta-Temme: z = erfcinv(2σ), invert beta's ridge mapping
     (clean-room from the published Temme beta expansion; forward R3's
     e_k(ζ,p) machinery direction-reversed). Probe floor (plain CLT)
     already wins the balanced ridge.
   - S2 small-y series inversion: y₀ = exp((ln σ + ln α + lnB)/α) +
     Picard via the R1 series (count replay-pinned); swapped twin
     covers the other end free. Wins R1-tiny (11/12) and the moderate
     plurality.
   - S3 gamma-limit transfer: for huge β, map t = −β·log1p(−y) and
     seed via the EXISTING GammaInvSeedS1/S2/S3 template functions
     (src/gammainv-inl.h; cross-family include, ErfcinvVec precedent),
     inverting y = −expm1(−t/β). Probe: wins gamma-limit by 15–52
     bits with a sharp measured seam at α ≈ 20 = kGammaAT (reference
     the constant, never duplicate it).
   - S4 joint-tiny logit closed form: logit(y) = (s − s*)/w + c(α,β),
     w = αβ/(α+β), s* = β/(α+β) — from the u = logit(t) substitution
     B_y = ∫^{logit y} exp(−α·ln(1+e^-u) − β·ln(1+e^u)) du (integrand
     → min(e^{αu}, e^{−βu}); c = 0 at α = β by symmetry). G1 derives
     c(α,β) and the correction order on paper, pins the route gate
     (max(α,β) < t_jt, which MUST own min(a,b) ≤ 2^-52 with margin)
     and the large-|logit| seam onto the power-law/deep-small form.
4. STEPS: safeguarded logit-Newton m = lnP − lnQ, 3 shared steps —
   the gammainv G3 package carried whole (reject residual-increasing
   steps, 1/8 backtrack, bypass |resid| < 1/2, additive-y step with
   relative floor; no saturation clamp, forward returns logit+slope).
   Forward: dd assembly of beta's region cores (R1 series / R2 CF /
   R3 Temme / gamma-limit via gamma cores), lnB via the LgammaDiffDd
   identities, prefactor in log space (params to 1e308 — the E
   dual-form lesson applies). Probe: 3 steps from a 6-bit seed reach
   the noise floor in EVERY region; the ridge inherits NO external
   penalty from the forward's 3 ULP. G1 replay: per-point analytic
   eps wherever series super-converge (SECOND-correction lesson,
   q-side twins included from day one), always solving against the
   root of the ROUNDED double s.
5. DEEP-SMALL closed form at both ends: y = exp_dd((LogDd(σ) ⊕ ln α ⊕
   lnB_dd)/α), mantissa + exponent, scaling last; cut on the
   DROPPED-FACTOR error < 2^-60 measured in BOTH orientations from
   the start (gammainv G3 deviation-1: the single-orientation
   self-check left 90 ULP reachable).
**Targets** [ILLUSTRATIVE until G4; pin to measured, no margin]: 1–2 ULP
relative, both sides, full domain, EXCLUDING the two named buckets:
plateau κ > 2^52 (backward-error ≤ 2-ulp-class contract) and huge-ν
beyond-resolution (xd at least as close as either neighbor). scipy
betaincinv baseline [probe B-P1d]: median 2.7 ULP, p99 4.7e11, max
5.4e14 — it collapses near s = 1; the lossless near-1 story via the
swap identity is exactly the gap.
**Oracle (G2; frontier-specified — THREE binding constructions beyond
the gammainv pattern)**:
1. FAST-PATH forward evaluator for R1-tiny/joint-tiny certification:
   plain mpf series at target dps, bypassing small_side_direct's
   escalation ladder — measured 100× per-call cost there (400–524 ms
   vs 3–6 ms; 43 s/row unseeded ⇒ infeasible at any stratum size).
   Validate fast-vs-full on a stratum sample, then certify with the
   fast path plus layered spot-checks. Seed every root-find from the
   winning seed candidate (iteration savings are real but secondary).
2. GUARD the reused gamma-corner route AT THE ENFORCEMENT SITE:
   small_side_direct HANGS for both params ≳ 1e17 balanced
   (gamma_corner_value feeds min(a,b) to mpmath.gammainc as a shape
   argument unconditionally — untested at huge shape). Bound the
   shape argument; route both-huge-balanced traffic through an
   R3-Temme extraction (gen_beta_data.py gamma_ck machinery),
   dual-anchored per gammainv G2's route-2 — this same route is the
   huge-ν stratum's independent certification.
3. Plateau rows: κ per row; κ ≤ 2^52 → normal bracket certification;
   above → BACKWARD-ERROR certification (forward of the stored y at
   dps 100 within the contract) — no y-bracket exists to certify
   there. Deep-small rows: log-space certification (gammainv pattern).
Everything else per gammainv G2: half-ulp midpoint sign-flip bracket
certification, layered dps 60→100, negative controls baked in with
exit 2 (probe prototype already rejected 3/3).
**Reference strata (G2)** [probe B-P5; ~14–21k logical rows]: R1-tiny
both orientations (4–6k); ridge balanced + skewed sub-bands + the
S1/S3 skew seam (3–4k); gamma-limit dense at the α ≈ 20 seam (2–3k);
joint-tiny plateau band, SEPARATE BUCKET (1.5–2.5k); underflow
thresholds both ends across the widened a*(b) boundary (1–1.5k);
subnormal-y both ends (0.8–1.2k); huge-ν beyond-resolution, SEPARATE
BUCKET (1–1.5k); a_T-seam bit-stepped bracket (0.5–0.8k); specials
smoke (~250). The swap identity HALVES orientation coverage (one of
(a,b)/(b,a) per logical point) EXCEPT near-diagonal and plateau rows,
where the swap maps s ↔ 1−s — both orientations needed there.
**Kernel/TU**: src/betainv-inl.h + betainv.cpp, both exports one TU;
consumes beta-inl.h region cores + gammainv-inl.h seed machinery +
dd/dd_special + exp_dd/log_dd + lgamma internals. HWY_NOINLINE day one
on cores AND driver; /d2ReducedOptimizeHugeFunctions day one (real
MSVC only). MSVC BUILD-TIME GATE: heaviest TU yet (beta cores AND
gammainv seeds, instantiated twice per export) — if the Windows CI
build pushes past ~18 min, ESCALATE before G4 (mitigations: audit
which cores the TU actually instantiates, TU split). Tests at the END
of all FOUR lists.
**FIRST CORRECTION [2026-08-09, G1 escalation, chain depth 1 —
RATIFIED, resolution in flight]**: G1's replay found a seed-coverage
gap — near-symmetric moderate-tiny (a,b) between t_jt = 2⁻⁸ and S1's
ν ≥ 2 boundary with a moderate target (witness (a,b,y) =
(0.02, 0.02, 1e-4), 0.00 bits post-steps; bucket tracked un-gated,
correctly escalated). Frontier ruling: t_jt gates the CLOSED-FORM
ROUTE, not S4's candidacy — the contract's candidates are global, and
excluding S4 above t_jt is the suspected primary cause (at C = 0.04
its leading form still seeds several bits). Fix hierarchy ratified:
(a) offer S4 globally; (b) pin next-order S4 correction in C
(measured residual coeff ≈ −0.4817·C; L-dependence to be determined)
holding ≥6 seed bits as far up in ν as it reaches; (c) lower S1's
ν-gate via a log-space λ(ζ) Newton (also retires the generator's
niter=100 wart rather than bequeathing it to G3); (d) fifth-seed fit
only as last resort. Mechanism must be named before the fix is
pinned (S4-not-offered vs residual-selection failure in the log-flat
band vs step traversal — the last is a NEW escalation, not a silent
step-count bump). Acceptance: the gap bucket JOINS the hard 55-bit
gate; only plateau-contract and beyond-resolution stay outside the
y-ULP gate.
**FIRST-correction resolution [2026-08-09]**: mechanism (i) confirmed
at the witness — S4 was never offered (candidate-gated at t_jt, a
contract deviation; fallback S3 seeded 3.8e-14 bits). Fixes landed:
S4 global candidacy (t_jt now gates only the closed-form route); the
S4 linear form's error measured O(C·|L|) — replaced by the EXACT
leading-order relation (exp-form with exact lnB, error O(C) uniform
in L; c(α,β) DROPPED — it hurts once B is exact, 3.95b → 0.17b at
C=2), which is S2's own zeroth iterate and thereby exposed two
independent S2 bugs (missing +ln α term, catastrophic at tiny α; a
q-side orientation bug complementing σ with unswapped a,b). Fix
level (c) tested and RULED OUT by measurement (S1's floor at ν=0.18
is asymptotic O(1/√ν), not Newton convergence). Fifth seed ratified
and added: logit-normal via exact ψ/ψ₁ moments of logit(Y)
(Gamma-ratio identity, clean-room); Cornish-Fisher skewness tested
and rejected. Also caught: eps_for's region-name gate missed
plateau-adjacent points labeled "gap" (74.38b → 25.38b regression),
generalized to min(a,b). Floors after cycle: S1 59.27 / S2 59.49 /
S4 74.38 / plateau 67.73.
**SECOND CORRECTION [2026-08-09, G1 re-escalation, chain depth 2 —
RATIFIED]**: residual gap sub-band (min(a,b) ≈ 0.02–0.5, skew
3–10×, y interior 0.1–0.3) is STEP TRAVERSAL, not seed quality —
the named tripwire, correctly re-escalated rather than silently
fixed. Seeds top out at 2–5 bits after five families honestly
exhausted; measured convergence is cleanly quadratic
(2.12→6.66→16.48→36.16→75.51 bits), step 4 clearing the gate by
20+ bits band-wide. RULING: StepsN = 4 shared. The gammainv "fix
the seed, don't shave margin" precedent does not apply: the seed
side is exhausted, nothing is shaved (full 75+-bit margin
restored), and the safeguard package makes a fourth step idempotent
for converged lanes — cost bounded at one forward eval, strictly
cheaper than a sixth fitted seed. G3 latitude: MAY add a
whole-vector all-converged skip after step 3 (bench call; gates
must hold either way). Closure requires full replay at 4 steps
everywhere (no floor may regress), band-wide basin verification
from the worst seed, gap bucket joining the hard 55-bit gate.
SECOND-correction closure verified: 4-step floors S1 59.27 / S2
95.83 / S4 74.38 / plateau 67.74 / former-gap 74.75 (joins hard
gate); 2,592-point band sweep, zero below 55b, worst seed 2.06b →
84.73b; byte-reproducible, rc=0.
**THIRD CORRECTION [2026-08-09, orchestrator review of closed
deliverables — RATIFIED, resolution in flight]**: the deep-small
piece carries the gammainv G3 deviation-1 disease in beta form,
three related defects confined to that subsystem. (1) deep_small_y
missing +ln a (the exact bug the agent fixed in seed_S2 this cycle,
surviving in the twin function; masked as a seed candidate by
residual selection). (2) The pinned cut a·y < 2⁻⁶⁰ has NO
b-dependence, but the dropped-factor error is |1−b|·y/(1+a) — the
OTHER side's parameter is the leading coefficient (gamma had no
second parameter; that's why its fixed form was x₀(1+a)). Reachable
witness: a=0.9, b=1e5, y=2⁻⁶⁰ (σ ≈ 2e-12, all normal doubles) —
route fires, ships ~410 ULP. (3) Self-check (f) validated none of
it: grid never samples the cut boundary (violates the binding
edge-refined rule), b fixed at 5, and the q-side loop is DEAD CODE
(pre-swapped args + side="p" make ax0 ≈ b, always skipped) — "both
orientations" swept zero q rows, the single-orientation hole the
contract named verbatim. Fixes: correct formula both branches;
re-derive + re-pin the cut from the dropped-factor bound (class
|1−b|·y/(1+a) < 2⁻⁶⁰ + exact q twin), measured boundary tightness
both orientations across b to 1e300; rebuild (f) with bit-stepped
boundary sampling; directed audit of checks (a)–(h) for the two
disease classes (boundary never sampled; dead orientation branch).
**THIRD-correction resolution [2026-08-09]**: all three defects
confirmed and fixed. The agent's measurement then found the
orchestrator's leading-order bound ITSELF insufficient at the
widened gamma-limit corner (true/bound ratio to 13.8 — the ln(S′)
linearization needs the leading term small, which fails at huge
other-side parameter); resolved with an exact closed-form
multiplier corr(y′) = −ln(1−y′)/y′, exact in the huge-other-side
limit (S′ → (1−y′)^(other−1)), verified sound (ratio ≤ 1, worst
1.0000000004 = boundary float noise) across 1,141 bit-stepped
boundary points, both orientations, b < 1 through b = 1e300. Final
route: P |1−β|·y/(1+α)·corr(y) < 2⁻⁶⁰, Q twin mirrored. Witness
re-tested: correctly rejected (true error 205-ULP-class). Directed
audit found ONE more real instance: check (b)'s "deep-small both
orientations" block appended 5-tuples into a list dispatched on
len==4 — 20 points silently dropped every run (removed, superseded
by rebuilt (f)); remaining checks clean or N/A, one low-risk note
((a) approaches but does not bit-step its domain edge; continuous
check, not a route decision). Floors unchanged, rc=0,
byte-reproducible.
**Stage record**: G1 SHIPPED [2026-08-09, Sonnet, three ratified
corrections — FIRST chain depth 1→2 (S4 candidacy + exact-B form +
two S2 bugs + fifth seed S5 logit-normal via exact ψ/ψ₁ moments;
Cornish-Fisher tested-and-rejected), SECOND (StepsN=4, the named
step-traversal tripwire, correctly re-escalated), THIRD
(orchestrator review: deep-small transfer-bug cluster)]:
gen_betainv_data.py + src/betainv_data.h (162 lines). Final floors
S1 59.27 / S2 95.83 / S4 74.38 / plateau 67.74 / former-gap 74.75
vs 55-bit gate; 2,592-point gap-band sweep zero below gate. Pinned:
S1 K=2 15×9 Chebyshev (ζ ≤ 3.5, ν ≥ 2 — the 1/ν series diverges
below, measured), S2 Picard 6, S4 exact-B global candidate (t_jt =
2⁻⁸ PROVISIONAL as closed-form-route gate only), S5 tableless, 
StepsN=4 TrustResid=1/2, deep-small cut per THIRD correction.
STAGE THEME for the G3 brief: every G1 defect except the SECOND
correction was a TRANSFER BUG — gamma formulas or gamma-shaped
assumptions carried where beta's second parameter changes the math
(√2 scaling, missing ν factor, missing +ln α twice, b-independent
deep-small cut, single-orientation checks). G3 must treat every
gammainv-inherited formula as UNVERIFIED until re-derived for beta.
**G2 IN FLIGHT — paused mid-stage [2026-08-10, session pause; two
agent cycles done, deliverables UNCOMMITTED→committed as WIP this
change set]**: gen_betainv_reference.py + 9,352 certified rows
(6,084 p / 3,268 q; five-hex rows a b sigma yd marker, marker
N/P/B — P = plateau backward-error contract, B =
beyond-resolution). ACCEPTED as implemented and validated: the
three binding constructions — fast path (worst fast-vs-full
disagreement 4.1e-59 on 40-point sample), gamma-corner hang guard
(hang confirmed at shape as low as 1.6e17 when near gammainc's OWN
ridge — proximity-dependent, not magnitude alone; guard threshold
1e10 + exception rescue, external wrapper only) with dual-anchored
R3-Temme route 2, per-row κ-split (9 plateau-contract rows);
negative controls 5/5 rejected first on every run; route-2
dual-check caught 1 real disagreement → row correctly DECLINED.
OPEN ITEMS for the completion round (next session, same agent or
fresh with this record): (1) huge-ν calibration harness DEFECT —
its "balanced resolvable to ν=1e60 under z=3 probe" claim is
impossible (frontier direct check 2026-08-10: balanced
central-band collapse starts between ν=1e31 and 1e32; at 1e60 even
σ=0.3 rounds to yd=0.5 exactly); root-cause the probe (suspect
mpf-space comparison or residual ζ-vs-z confusion). The SHIPPED
bucketing survives review anyway: B-floor ν=1e35 sits inside
central-band collapse and near the max-z full-collapse boundary
(~5e35 by direct estimate; probe's 3e33 was the ±6σ convention) —
B rows are certified under neighbor semantics and safe. (2)
N-marked huge-ν rows between collapse onset (~1e32) and 1e35 are
trivially-satisfiable dilution rows — G4's ULP test must bucket
huge-ν separately regardless (gammainv dilution lesson); relabel
or bucket-by-formula decision owed. (3) r1-tiny p/q imbalance
17:1 (1872:110) despite nominally-symmetric construction —
root-cause owed. (4) gamma-limit-seam declines: 365/409 declines
concentrate there (root-finder converges to a wrong point at
b ~ 1e111–1e250; bracket certification correctly rejects — no bad
rows shipped, a YIELD problem). (5) scale 9.4k → design's 14–21k;
subnormal-y/underflow strata still single-sided by construction.
G2 self-caught bug ledger so far: 11 across both cycles (worst:
both-huge bisection at 500 ms/call inside a 150-iteration loop;
q-side round-to-one threshold off by ~700 decades; blind-σ
deep-small construction 77% drop → inversion-first fix → ~12%).
FRONTIER RULINGS [2026-08-10, open items 1–3 CLOSED; completion
agent launched with them]:
(1) Calibration probe root cause: its "resolvable" criterion tested
forward-VALUE saturation (P/Q == exactly 0.0/1.0 in double at the
z-probe point) — not y-space collapse; the balanced case
degenerates because the probed y rounds to the mean double 0.5
where P = 0.5 by symmetry, never saturating at ANY ν (the
"resolvable to 1e60" artifact); the skewed 2e36–4e36 numbers
measure a third quantity (rounding perturbation exceeding the
whole z-range). CORRECT criterion, per-(a,b): y(z = ±Z_MAX ≈
±38.5, the subnormal-σ limit) double-rounded ≤ 1 ulp apart ⇒
beyond-resolution. Anchors: balanced central-band onset between
ν = 1e31 and 1e32 (measured directly); balanced full collapse
~5e35-class. Shipped B rows (ν ≥ 1e35) accepted as-is.
(2) Dilution: NO relabeling — G4's ULP test buckets huge-ν
statistics BY FORMULA from (a,b); markers carry certification
semantics only. Rule recorded in the generator's format docstring.
(3) r1-tiny 17:1: STRUCTURALLY EXPECTED, not a bug — log-uniform-y
sampling with tiny a puts the P ≤ 1/2 crossover at y_med ≈
exp(−ln2/a) (a = 0.01 → ~1e-30), so nearly all sampled decades are
p-side. Coverage intent still binds: scale-up constructs q-side
rows DIRECTLY (σ-targeted, inversion-first). Same for the
single-sided subnormal-y/underflow strata.
G2 SHIPPED [2026-08-10, two agents: original (3 cycles) + fresh
completion round carrying the frontier rulings]: 16,883 certified
rows (8,603 p / 8,280 q, 99.32% of constructed; 115 declines, all
boundary-ladder/root-find classes), five-hex format a b sigma yd
marker. Rulings implemented: Z_MAX = 27.2005633… derived by mpf
bisection in the file's own erfc z-convention (≡ the frontier's
38.5 normal-quantile figure), both sanity anchors reproduced
(central-band onset 1e31–1e32 exact; full collapse 6.0e34–1.3e35);
bucket-by-formula note in the format docstring; q-side coverage
r1-tiny 1.68:1 / subnormal-y 1.26:1 / underflow 1.90:1 (was 17:1 /
q=0 / q=0). Completion round's own catches: (i) q-branch
construction built near-1 arguments as native-float 1.0−y —
collapses to exactly 1.0 for y < ~2⁻⁵³ (96.5% of the r1-tiny draw
range; witness float(1.0−1e-150).hex() == 0x1.0p+0) — the
escalation was raised rather than folded in silently; (ii) seam
declines root-caused DEEPER than the ruling: gb.route_final
silently misroutes a neighborhood just past the true root into
small_val_via_cf, which is INVALID at extreme skew and returns
catastrophically wrong values with no exception (witness: true
P=0.5207 vs computed 1e-4541) — fixed by bisecting against the
audited evaluator directly (oracle_y_audited, S3-seeded); seam
declines 66% → 2.5%. SHARED-MACHINERY CAVEAT recorded: the
route_final/small_val_via_cf silent-garbage combination is
reachable when shipped generator machinery is driven at
out-of-domain points; the shipped beta FORWARD reference set is
believed unaffected (its constructions stay in-domain) but any
future reuse must guard routes the way this generator now does.
ACCEPTED SHORTS: underflow 806 vs 1–1.5k target; huge-ν B-bucket
646 vs 1–1.5k (doctrine-correct: the Z_MAX criterion bounds the
outer envelope, specific σ draws inside it legitimately certify N
— those rows still enter the huge-ν formula bucket at G4).
Negative controls 5/5 on all ~21 invocations; ~2.5 h compute,
checkpointed. Orchestrator review: 12/12 independent mpmath
bracket spot-check (moderate-param rows, dps 60), zero format
defects.
G3 SHIPPED [2026-08-10, Opus, ONE escalation adjudicated + NINE
measured deviations, 8 self-caught bugs]: kernel + smoke + ULP +
bench + four-list registration (confirmed, dependency position).
Measured, IDENTICAL tables clang-cl AVX3_ZEN4 native and g++
SSE2-capped: deep-small 0/0; huge-ν formula bucket (ν ≥ 1e31,
from (a,b) alone per ruling 2) 1/1; ridge 1/1; gamma-limit 1/1;
small-param remainder 1/1; B rows 1/2 (neighbor semantics); P rows
+ κ-bucket backward contract 0.000 ulp(σ); subnormal-x and x=1
cross-cuts 0. ESCALATION ADJUDICATED (frontier): the κ contract
boundary is 2¹⁸, not the design's 2⁵² — the design assumed a
dd-accurate forward, but near the median the logit's ln(1−u) chain
rides exp_dd's ~2⁻⁷⁰ budget; κ·2⁻⁷⁰ crosses 2⁻⁵³ at κ ≈ 2¹⁸.
ACCEPTED: the 731 affected rows all meet the backward contract at
0.000 ulp(σ); semantics unchanged, band boundary moves; exp_dd
upgrade recorded as optional enhancement (Open Items). Key
deviations (full list in the G3 report/git history): orientation
frame by definitive median probe (the "σ ≤ 1/2 ⇒ y below median"
brief wording was wrong when β ≪ α); S1 λ-inversion in √-space (8
iters vs the generator's niter=100 raw-λ oscillation); S1 offered
globally (the ν ≥ 2 gate was the correction table's, not the
seed's); one ADDITIVE field on BetaR3Out (tail bracket in log
space — no arithmetic changed) so the inverse's residual has no
subnormal flat spot; residual-uncertainty freeze (beta's w reaches
2⁵⁰ vs gamma's 2¹⁰ — trust bypass was accepting noise steps);
StepsN=4 kept, all-converged skip DECLINED (seed stage dominates).
MSVC BUILD GATE RESOLVED, not deferred: betainv.cpp >45 min/7 GB
before outlining ~20 log/exp call sites through HWY_NOINLINE
wrappers → 127 s, lighter than beta.cpp; ULP tables byte-identical
across the change. Also FOUND (not fixed): two defects in the
shipped beta forward at u → −1 (Open Items, PRIORITY). Bench
indicative: 0.59–11.4 µs/el; cost dominated by up to 12
region-routed forward evals (probe + 7 candidates + 4 steps);
candidate count is the throughput lever if wanted. Orchestrator
review: 23/23 ctest re-run verified on the agent's tree; BetaR3Out
single construction site, no default-construction anywhere in the
new TU (bf16 pattern grep clean).
G4/G5 COMPLETE — **BETA_P_INV / BETA_Q_INV SHIPPED [2026-08-10]**:
gates PINNED to measured, no margin (1 ULP every y-bucket both
sides; B rows 2 vs certified answer; backward contract 1.0 ulp(σ),
measured 0.000; deep-small/subnormal/x=1 CR). Full ladder
asserting under CORVUS_EXPECT_TARGET: AVX3_ZEN4 native;
AVX2/SSE4/SSSE3/SSE2 capped clang-cl sweep (all tiers, run under
pwsh in the VS dev env); Linux CI sweep + sanitizers 10.3 min;
NEON; Windows MSVC 17.4 min (watch item stands — heaviest yet,
timeout 25; betainv.cpp itself 127 s post-outlining). CI run
31448781077 verified by SHA and per-job conclusions. ACCURACY.md
matrix rows + family section (κ contract band, swap-identity
lossless-near-1, oracle record), README, in the gate-pinning
change set. POSTSCRIPT lessons: (i) the G3 commit failed
Linux/macOS CI on three -Werror unused-variable hits — the G3
agent's g++ leg had CORVUS_DEV_WARNINGS=OFF (build-cap trees carry
it OFF; only CI's dev-warnings build catches this class) — fixed
after verifying all three were genuinely dead, swept to zero under
a dev-warnings g++ build; (ii) sweep_tiers.ps1 invoked under
Windows PowerShell 5.1 silently fails to apply the pipe-delimited
cap (cap didn't bite, AVX3_ZEN4 ran under the AVX2 name — the
script's own expect-target assertion caught it); use pwsh; (iii)
the g++-default sweep hits the standing mingw exit-crash item
(gamma prints PASS then segfaults at teardown → false gate
failure) — the validated sweep compiler on this box is clang-cl,
matching the ENVIRONMENT.md rule.
Next: Bessel I0/I1 + lbeta (P2, staged below) from a fresh fork.
G3 BRIEF INGREDIENTS (compose at launch): transfer-bug stage theme
(every gammainv-inherited formula UNVERIFIED until re-derived —
G1's ledger is the witness list); marker-column semantics (N/P/B;
P = backward-error contract rows — the test verifies |forward(yd)
− σ| ≤ contract, NOT y-ULP; B = neighbor semantics; huge-ν
statistics bucketed by formula per ruling 2); StepsN = 4 with the
whole-vector all-converged skip latitude; deep-small cut is the
THIRD-correction form (other-side coefficient × corr(y′), NOT
a·y); S5 needs the kernel's own digamma/trigamma cores (tableless);
S3 calls GammaInvSeed* cross-family (include gammainv-inl.h);
/d2 flag + MSVC build-time gate ~18 min; four-list registration at
END; betainv-inl.h + betainv.cpp one TU two exports.
**Process**: G1 (Sonnet, gen_betainv_data.py → src/betainv_data.h:
replay with per-point analytic eps, edge-refined bit-stepped sampling,
both-orientation deep-small validation, c(α,β) derivation, t_jt +
seam pins) → G2 (Sonnet, SEPARATE agent — the oracle is the risk item;
the three binding constructions above are its brief) → G3 (Opus
kernel) → G4/G5 orchestrator. Escalation-density judgment rule
applies. Probe stage record: 5 self-caught bugs (linear-space
bisection unusable at y ~ 1e-300, fixed to ln-space; a positional-arg
swap feeding y_true into the target slot — plausible-looking wrong
numbers, caught only by hand-deriving one point; the
small_side_direct hang, two orphaned PIDs killed per process rule;
CF non-convergence boundary mapped at ν ~ 1e18 — real R3 territory,
not a bug; one grid point silently in the wrong regime, left
documented as a non-fit rather than dropped).

## P2 Bessel I0/I1 — STAGING [2026-08-10, pre-probe; design next session]
Last planned family; smallest since erf/erfc. Estimate: ONE session
for the full pipeline, ~half trigamma's effort.
**Consumer requirements** (the reason the family exists): von Mises
log-density needs log I0(κ) at large κ — plain I0 overflows past
x ≈ 713.5 (e^x/√(2πx) vs DBL_MAX) — and the von Mises MLE κ-update
needs the ratio A(κ) = I1/I0. Plain I0/I1 for densities at small κ.
**API candidates [design decides]**: i0, i1 (unscaled, overflow past
~713); i0e = e⁻ˣI0, i1e (scaled — the kernel-natural pair, full
axis); log-I0 direct vs composed log(i0e(x)) + x — composition is
accuracy-fine (the log term is O(ln x) against x, its error vanishes
relatively; verify at probe) so the question is API convenience only.
A(κ) = i1e/i0e composes exactly (scaling cancels) — likely document,
not ship. Negative axis free: I0 even, I1 odd.
**Structure** (two regimes, both precedented): small-x power series
Σ(x²/4)ᵏ/(k!)² — ALL-POSITIVE terms, no cancellation, perfectly
conditioned; large-x scaled asymptotic → clean-room Chebyshev refit
of e⁻ˣI0(x)·√x-class in 1/x via the standard generator approach
(erfc-tail machinery is the direct precedent; A&S-form coefficients
are NOT to be ported — refit from scratch, own nodes and budgets).
Unscaled assembly: i0e · e^x via exp_dd mantissa+exponent (scaling
last), inheriting exp's conditioning honestly — document like the
erfc tail. i1e ~ (x/2)e⁻ˣ near 0, relative-clean, no hazard.
**Conditioning**: i0e/i1e bounded mild condition numbers, no zeros
on the axis, no reflection, no ill-conditioned bands. The first
family since erf with NO adversarial conditioning story.
**Oracle**: mpmath besseli has no known pathology classes —
erf/lgamma-difficulty verification (layered dps + ONE independent
cross-check route, e.g. own-series at high dps vs besseli), no
bracket certification. Doctrine minimums still apply (negative
controls in the generator self-check, oracle-trust posture).
**Probe questions (light — probe + design in one frontier pass)**:
(a) series/fit split point and fit degree/budget per function
(1-ULP target expected; the fit region may carry 2 like erfc's
tail); (b) composed-log-I0 accuracy measurement (close the API
question with numbers); (c) i1 sign/odd-extension and specials
policy (x=0: I0=1 exact, I1=0 exact; NaN propagation); (d) overflow
boundary exactness for unscaled forms (saturate to +inf past the
measured cut, erfc-underflow precedent in reverse); (e) [from the
2026-08-10 milestone sweep] von Mises CDF (libstats #51) needs
Σ I_j(κ)·sin(jθ)/j — higher-order I_j via standard backward
recurrence from I0/I1: DESIGN DECIDES whether i0e/i1e alone
suffice (libstats recurs on its own; document the recipe) or a
recurrence helper is worth shipping — measure the recurrence's
stability/accuracy from corvus seeds before deciding.
**lbeta (committed P2, after Bessel — possibly same session)**:
public ln B(a,b) exposing the internal LgammaDiffDd assembly (the
a+b cancellation hazard is already solved in-house); consumer
drivers: BetaBinomial PMF hot path (2 evals/point), F/StudentT/
Binomial delegations. Trivial oracle (mpmath directly), thin TU or
co-located with an existing lgamma-family TU per the dependency-
boundary rule — design decides placement. erfcx stays optional:
Mills-ratio consumer (TruncatedNormal far truncation) is
speculative, no filed need.
**Kernel/TU**: src/bessel-inl.h + bessel.cpp, one TU, exports for
the chosen API set (shared series/fit cores); consumes exp_dd only.
HWY_NOINLINE day one; four-list registration at END. MSVC: expected
LIGHT (no heavy core instantiation) — no /d2 unless measured.
**Effort routing**: probe+design frontier (one pass); G1 fits
generator + G2 references Sonnet (G2 is light here); G3 mid-tier
candidate per the effort table (settled two-regime design — Sonnet
with escalation rights; Opus only if the design pass flags risk);
G4/G5 orchestrator. Transfer-bug theme does NOT carry (no gamma/beta
formula inheritance) — but the erfc-tail fit-budget discipline does.

## GitHub repo settings [applied 2026-07-21 via gh api]
Merge: all three styles, auto-delete head branches (PR merges only —
a direct fast-forward push bypasses it; prune manually). Wiki and
projects DISABLED (four-file docs policy), issues on, discussions off.
Security: Dependabot alerts + auto fixes, secret scanning + push
protection, private vulnerability reporting. Ruleset "protect-main":
blocks force-push/deletion, direct pushes allowed (solo workflow).
Ruleset "protect-tags" [2026-08-06, id 20491885]: v* tags immutable
(deletion + update blocked, no bypass actors; creation open). Actions
GITHUB_TOKEN read-only, cannot approve PRs. Deferred: signed-commits
rule (see Open Items). Required status checks deliberately absent
(incompatible with direct-push workflow).

## Build-stack standardization (2026-07-23) [DERIVED]
Cross-repo effort tracked in the fleet standards repo
([record](https://github.com/OldCrow/standards/blob/main/records/BUILD-STANDARDIZATION-PLAN.md)).
corvus commits: pkg-config file + consumer example + installed-path CI
check; find_package(hwy 1.4) version floor with CI building pinned
Highway 1.4.0 from source; CMakePresets.json.

## Resolved log
One line per closed item; detail in this file's git history, AGENTS.md,
and docs/ACCURACY.md.
- 2026-08-09 gen_beta_data.py kBetaGammaLim standing "frontier review
  owed" flag (2^-49 target deviation, from the beta G1's escalation (C))
  reviewed and RATIFIED: routing threshold not truncation depth (end-to-
  end ULP bar applies), pin midpoint target-invariant (log-linear
  frontiers), empirically confirmed by beta's shipped gammalim gates +
  boundary-crossing seam sweeps; flag text retired at the print site.
- 2026-08-06 v0.1.0 first release: tag at b4eaeea gated on full CI
  green; issue #2 closed (fleet lint-workflows.yml, SHA-pinned actions,
  persist-credentials: false); protect-tags ruleset; install/export
  decision ratified; Release object published; beta/g3-kernel pruned;
  README beta example + release badge.
- 2026-08-06 gamma a > 1e4 oracle branch independently spot-checked
  42/42 — oracle-trust directive fully discharged (both families).
- 2026-08-05/06 incomplete beta shipped end-to-end: eleven corrections,
  gates pinned to measured, monotonicity + ten seam sweeps clean, full
  tier ladder + NEON gate-cell-identical, harness-certified references,
  ACCURACY.md/README/AGENTS.md updated in the change set.
- 2026-08-06 macOS CI bf16 failure fixed (aggregate-init of
  vector-member structs) + clang-cl -Wall/=/Weverything flag mapping
  corrected to /W4 -Wpedantic (frontend-variant gate); local clang-cl
  build warning-free.
- 2026-07-29 dd_special hoist (Log1pmxDd/Expm1Dd → dd_special-inl.h)
  under the byte-identity protocol; reference file renamed verbatim.
- 2026-07-29 MSVC codegen blowup closed: outline region cores AND
  per-lane drivers from day one (AGENTS.md rule); gamma + erfinv
  outlined, CI Windows job 25-min timeout → 6m55s; erfinv's central
  poly deliberately left inline (measured 0.5–0.7 ns/el cost for no
  relief).
- 2026-07-28 gamma_p/q shipped: 2 ULP direct on every region, four x86
  tiers cell-identical; Sonnet-tooling + Opus-kernel split with
  orchestrator review gates.
- 2026-07-25 erfinv/erfcinv shipped (max 1 ULP everywhere); erfc-tail
  open item closed by measurement (attenuated, as predicted).
- 2026-07-25 Phase B lgamma and Phase A (exp_dd/log_dd, erfc tail
  rewire 5 → 2 ULP) shipped; NEON rows from CI.
- 2026-07-25 HWY_DYNAMIC_DISPATCH must be called inside namespace
  corvus (SSE2-cap collapse); sweep_tiers.ps1 aborts on first failure.
- 2026-07-24 AVX-512 validated natively on Ryzen; Windows/MSVC CI job
  added (toolchain coverage); CORVUS_EXPECT_TARGET everywhere after
  three silent-cap defects; HWY_BROKEN_MSVC investigated,
  CORVUS_MSVC_UNBLOCK_AVX512 exists default-OFF.
- 2026-07-21 CI designed around runner-minute economy; repo settings
  applied; conventions audits recorded in AGENTS.md.
- 2026-07-20 x86 gather perf: Kaby flat across widths, Zen 4 scales
  (erf 5.49 → 1.95 ns/el SSE2 → AVX3_ZEN4); non-gather variant noted
  as pre-AVX-512 upside only.
