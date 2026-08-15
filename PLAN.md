# corvus — Plan / Session State

Detail policy: this file holds STATE — what's decided, what's open, the
next concrete step, and lessons too expensive to re-learn. Full designs,
session narratives, stage ledgers, and measurement play-by-play live in
this file's git history (compacted 2026-08-06 and 2026-08-13),
docs/ACCURACY.md, and the kernel/generator source, which are the official
record for finished work. Binding cross-family engineering rules live in
docs/NUMERICAL-DOCTRINE.md, not here.

## Status [DERIVED] — 2026-08-14

**v0.5.0 RELEASED — CORE/GENERATOR/TEST FREEZE IN EFFECT** (signed
tag at 73eaee0, CI-gated per-job on runs 31849654346 + 31850384519;
https://github.com/OldCrow/corvus/releases/tag/v0.5.0): the
pre-v1.0.0 documentation trim + the huge-parameter Dekker audit
(three defect classes fixed, corner-row gated — Resolved log).
Anything after v0.5.0 (user-focused examples and other additions)
builds on the frozen base; core/generator/test changes now require
explicit unfreeze.

**v0.4.0 RELEASED** (tag at 96d181d, CI-gated per-job on run
31652276608;
https://github.com/OldCrow/corvus/releases/tag/v0.4.0): the
beta-forward u → −1 fix arc + the zero-warning static-analysis pass +
the CORVUS_SANITIZE MSVC guard.

**Pre-v0.5.0 documentation trim COMPLETE [2026-08-13]**: all four
phases (frontier hot kernels + data-header emits; Sonnet small
families/tests/ENVIRONMENT/README; ACCURACY compression + PLAN
compaction; generators under the lighter standard) shipped in
commits 6f050c5/c9d264f/ed7a759. Validation: all 8 data headers
regenerated from trimmed emit code, numeric sections byte-identical,
generator self-checks green; full ctest 27/27 on the clang-cl tree;
CI green per-job on run 31758023543 (SHA-verified). THE ARCHAEOLOGY
TEST [policy, ratified, applied]: a comment survives only if it
states a constraint, bound, derivation, or hazard the CODE cannot
show; who found it, when, in which pass, and what it measured before
the fix all go; rewrite, don't delete — keep each site's math half.
The lighter generator standard: every enforcement-site/self-check
rationale stays, timeless.

Remaining before v1.0.0: consumer-integration notes (libstats
#47/#51/#52), ARCHITECTURE.md referencing decision, Kaby machine
legs (non-gating), and the user-focused examples phase on the frozen
base.

Release history: v0.1.0 2026-08-06 (P0: erf/erfc, erfinv/erfcinv,
lgamma, gamma P/Q, beta P/Q); v0.2.0 2026-08-10 (P1: digamma,
trigamma, gamma_p_inv/gamma_q_inv, beta_p_inv/beta_q_inv); v0.3.0
2026-08-12 (P2: i0/i1/i0e/i1e, lbeta — with it corvus covers every
PDF/CDF/quantile/MLE special-function need in the libhmm/libstats
inventories, von Mises included); v0.4.0 2026-08-12. All tags
CI-gated, immutable under protect-tags, each with a Release object.

## Next Steps
1. **Consumer-integration notes** for libstats: #47 (A&S Bessel
   fallback capping VonMises at ~1e-7 — retired by i0/i1/i0e/i1e),
   #51 (von Mises CDF — Miller recurrence recipe documented in the
   bessel design record below), #52 (slow Binomial CDF — beta_p is
   the integration answer, a libstats note, not a corvus gap).
   NOW SOURCED FROM THE ADOPTION SPIKE (Open Items) rather than
   written cold: #47's note is stage S4's deliverable, #52's rides
   along as a one-liner, #51's waits on the adoption verdict since
   the recurrence recipe only pays off inside a real consumer.
2. **Kaby Lake legs when the machine is available** [OPEN, machine
   access; ruled NON-GATING 2026-08-06, user decision]: beta
   AVX2-native + capped sweep (additional cross-machine check, not a
   claim gap — ACCURACY.md dagger note), and the quiet-machine Kaby
   bench pass (its gamma bench numbers are loaded/indicative only).

## Open Items
- [OPEN, low] A few stage-tag strings remain inside generator stderr
  banners (runtime diagnostics only, never emitted into headers);
  clean opportunistically at the next generator touch.
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
- [OPEN, P2 candidate, not required] erfcx (nearly free from the erfc
  tail machinery; speculative consumer found 2026-08-10 —
  TruncatedNormal far-truncation moments/MLE run on the Mills ratio =
  erfcx — but no filed need yet).
- [OPEN, gamma] Add a ≥ 2^998 Temme witness rows at gamma's next
  reference touch — the kGammaTwoPiAClamp Dekker-ceiling defect was
  latent for lack of them (clamp already fixed to 2^900, 2026-08-05).
- [WATCH] Ryzen box stability: two GPU-stack bugchecks 2026-07-29 (0x9F
  power-IRP, 0x10E video memory), root-caused [DERIVED] to a
  half-committed NVIDIA driver install; DDU clean reinstall the same
  night. Crash-free since. Residual watch: whether the NEXT NVIDIA App
  driver update completes cleanly (chronic installer freezes implicate
  accumulated App state; manual driver-only installs are the fallback).
  Recurrence of either bugcheck on the clean stack flips suspicion to
  VRAM/hardware.
- [WATCH] mingw GCC 16.1 by-value-vector-argument misalignment bug,
  AVX2 and above: filed upstream 2026-08-08 as GCC PR 126741
  (https://gcc.gnu.org/bugzilla/show_bug.cgi?id=126741) against
  512-bit; an independently contributed repro, verified locally, shows
  the identical defect at 256-bit. Repros:
  `C:\Users\gdwol\Development\gcc-zmm-mingw-repro\` (512-bit),
  `C:\Users\gdwol\Development\gcc-ymm-mingw-repro\` (256-bit). Docs
  re-scoped 2026-08-10: mingw GCC qualified for 128-bit tiers only;
  past GCC AVX2 numbers stay valid (fault, never corruption).
  Re-qualify for AVX2+ only after the fix lands. The mingw test-binary
  "exit crash" is this same bug at 256-bit (initial-stack-residue
  trigger — an after-PASS fault is structurally impossible in these
  tests); technical detail in docs/ENVIRONMENT.md.
- [OPEN] Pre-release legal: BINARY artifacts linking Highway must carry
  its Apache-2.0 NOTICE; source-only distribution needs nothing (all
  releases so far are source-only). Handle when packaging starts.
- [OPEN] Decide whether libstats/libhmm adopt corvus as a dependency or
  keep their internal SIMD (separate project-level decision). SPIKE
  AUTHORIZED 2026-08-15 — libstats branch `spike/corvus-bessel` wires
  i0/i1/i0e behind stats::detail::bessel_* (its #47 path: three
  functions, eight call sites, all scalar, all in von_mises.cpp, none
  hot — the payoff is retiring its 1.6e-7 A&S tier, NOT throughput).
  corvus is consumed AT THE v0.5.0 TAG AND NOT EDITED; holding the
  freeze against a real consumer is part of what the spike tests, so a
  corvus-side need it surfaces is a FINDING, not a fix. Execution state
  (stages, next step) lives in libstats/PLAN.md In Progress and is
  deliberately not restated here — same rule as the pylibhmm pin: a
  copy here cannot notice when it goes wrong. This entry records the
  VERDICT, which has two INDEPENDENT halves: (a) does clang-cl-built
  corvus link into an MSVC consumer — corvus's question, ANSWERED YES
  2026-08-15 by stage S1, see the Resolved log; (b) does libstats adopt
  — libstats's question, still open. A "no" on (b) does not invalidate
  (a), and the #47 consumer-integration note is produced either way.
  SPIKE CLOSED 2026-08-15, recommendation ADOPT-BUT-NOT-FOR-BESSEL-ALONE:
  the mechanism works and corvus beats both of libstats's existing tiers
  measurably, but eight scalar call sites do not justify a dependency on
  their own; the justification is the wider surface (#47 retired, #51's
  Miller recurrence, #52's beta_p, plus erfinv and incomplete gamma/beta).
  Report: https://claude.ai/code/artifact/ab912f14-d920-4eb2-b60f-5d92222bb33f
  Rationale and costs in libstats/PLAN.md In Progress; the real decision
  still wants the macOS leg, since Tier 2 is where #47's users are.
  Three defects the spike found in libstats are filed there as #92/#93/#94
  and are independent of adoption — corvus owes nothing for them.
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
- [OPEN] Signed-commits ruleset: confirm the M1 and Ryzen boxes sign
  before enabling.
- [ILLUSTRATIVE] Possible future consumers: C++ port of multi-agent_sim
  (batch distance/trig), zeekhmm training pipelines.

## Decisions
- Name: corvus (OldCrow tie-in). Namespace `corvus::`.
- Scope: statistical special functions only; basic transcendentals
  belong to Highway contrib. P0 (done): erf/erfc, erfinv/erfcinv,
  lgamma, incomplete gamma P/Q, incomplete beta. P1 (done): digamma,
  trigamma, inverse incomplete gamma/beta. P2 (done): Bessel I0/I1 +
  scaled variants, lbeta.
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
  FRONTIER work (routing in AGENTS.md — the beta oracle cost more
  sessions than its kernel). References are trusted only after an
  INDEPENDENTLY constructed verification harness passes clean with
  baked-in negative controls; tools/verify_beta_reference.py is the
  pattern, gammainv/betainv per-row bracket certification the other.
  Exposure audit of the other oracles 2026-08-05/06: gamma's a > 1e4
  exact-asymptotic branch spot-checked 42/42 via three
  Temme-independent evaluators; the rest are thin wrappers over
  gold-standard single-argument mpmath calls.
- **Pipeline template** (P1/P2 families): probe → design →
  G1 (generator/data) → G2 (oracle/references) → G3 (kernel + tests)
  → G4 (gate pinning + tier ladder) → G5 (docs). Escalation-density
  rule [2026-08-06, user; refined 2026-08-08]: takeover is a JUDGMENT
  call, not a counter — take over on churn or local-fix myopia; clean,
  precisely-diagnosed escalation chains may continue delegated.
- **First release** (2026-08-06, user): v0.1.0; Kaby non-gating;
  install/export status quo ratified — `cmake --install` requires a
  system Highway (find_package(hwy 1.4)), FetchContent builds are
  build-tree-only; no bundling, no nested install (keeps source tags
  free of the NOTICE obligation); revisit only if packaging starts.
  protect-tags ruleset makes v* tags immutable. Every future v-tag gets
  a Release object (page-coherence cadence).

## Shipped families — what would be expensive to re-derive
Full method and measured bounds: docs/ACCURACY.md. Math: kernel
derivation blocks at the definition sites. Full design texts and stage
ledgers: this file's git history. Below: design points, pinned
parameters, and bugs worth not re-learning.

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
exact-asymptotic oracle above a = 1e4. Bench: Ryzen quiet-machine
49–74 ns/el per region, 8.9–21× scalar walk (upper bound).

### Incomplete beta P/Q [shipped 2026-08-05/06; corner arc 2026-08-12]
The hardest family; the routing/assembly corrections and six oracle
defect classes are recorded correction-by-correction in git history.
Kernel: src/beta-inl.h (own TU; consumes dd, dd_special, exp_dd,
log_dd, lgamma internals, erfc core, plus exactly two gamma template
cores for the gamma-limit slice). The surviving architecture, routing
order, prefactor identities, saturation/scrub rules, and pinned
constants are documented at their definition sites in beta-inl.h and
src/beta_data.h (generator self-checked); the derivations there are
the maintainer record.

**Oracle/harness record.** Six oracle defect classes vs ZERO kernel
defects in adjudication — the origin of the oracle-trust doctrine.
Disease classes (all now guarded at enforcement sites in
gen_beta_reference.py): (a) truncation-at-ambient-dps — mpf ops
truncate higher-precision operands, mp.mpf() on an mpf RE-ROUNDS,
1 + τ = 1, a + b collapses in lnB, log1p(−near-1) = −inf; exact
complements carried, never recomputed; (b) component-relative error
budgets voided by downstream cancellation; (c) mp.quad returns
√eps-scale noise when the integral sits below working epsilon —
normalize the integrand by peak log-magnitude; (d) false saturation
certificates from collapsed complements; (e) mpmath betainc returns
internally-CONSISTENT garbage in the gammalim corner (layered-dps
agreement cannot catch it — hard excluded); (f) small-side-direct
applies to the oracle too. tools/verify_beta_reference.py is the
INDEPENDENT harness (no oracle import; layered series + half-split
log-quad evaluator; exact analytic lines; saturation log-bounds); its
4 adjudicated rows are BAKED IN as a negative control (exit 2 unless
rejected). Clean pass = the trust gate for the shipped references
(37,099 rows/side + 251 specials, zero drops).

**Shipped numbers** (pinned gates, no margin, max over p/q; identical
gate cells on AVX3_ZEN4 native, AVX2/SSE4/SSSE3/SSE2 capped, Linux CI,
NEON): R1 1/1, R2 0/0, R3 3/1, R4 2/0, postroute 1/0, gammalim 0/1
(re-pinned dir 0 → 1 at the 2026-08-12 corner arc), specials exact.
Monotonicity post-pass (3,278 (a,b)-groups, kernel dip slack 4 ULP)
and ten 4001-point seam sweeps (each must CROSS its boundary or fail):
zero violations. Both live in test_beta_ulp. BETA_ULP_DUMP env prints
every not-CR row at/above a threshold (permanent gate-pinning tool).
bf16 lesson (2026-08-06 CI): never default-construct vector-member
structs — implicit ctors instantiate outside the per-target attribute
region; aggregate-init everywhere.

**Process lessons** (tooling, not math): mp.dps must be set INSIDE
every computation layer — worker, subprocess-send, subprocess-parse
(bitten three times). Probe with exact hex from dumps, never
display-rounded values. Windows multiprocessing:
current_process().name, not parent_process(); probes need
`if __name__ == '__main__'` guards or spawn bootstrapping masquerades
as fast failure. Sub-agent briefs: name the MECHANISM, not the concept
("never set run_in_background on any tool call"); foreground only,
chunk sweeps ≤ ~5-min re-runnable commands. Verify exponent arithmetic
by exponent SUM (8e-100·1e100 = 8, not 0.8).

### digamma [SHIPPED 2026-08-08 — first family through the agent pipeline]
Full real axis. Positive pipeline: product-form zone fit on [1,2)
around the root x₀ (dd shift, dd leading coefficients), one up-step for
(0,1) WITHOUT forming 1+x, masked down-walk [2, X0), Bernoulli
asymptotic at X0 = 8; negative axis by reflection ψ(1−x) − π·cot(πx)
with exact u = x − round(x), sinc-pair fits, y_dd = TwoSum(1, −x) with
a rough-trigamma lo correction. All mechanisms documented at definition
sites in src/digamma-inl.h and src/digamma_data.h.
**Accuracy doctrine** (the family's lasting contribution): relative
metric where |ψ| ≥ 1, 2^-53-class ABSOLUTE inside the per-zero bands
(near-relative at the 20 adversarial negative-zero doubles would need
2^-104-class fits — REJECTED on cost for a contract lgamma's negative
axis doesn't offer either). Gates pinned 1 ULP in all five relative
buckets, 1.0 × 2⁻⁵³ absolute band; full ladder identical cells
including NEON not-CR counts. Specials: ψ(±0) = ∓inf, every double
≤ −2^53 is an integer → NaN, subnormal x → ∓inf.
Bench loaded/indicative: 7.0–27.4 ns/el, 7.2–25.2× scalar walk.

### trigamma [SHIPPED 2026-08-08]
Digamma-shaped MINUS the hard parts; ALL-RELATIVE everywhere (ψ₁ is a
sum of squares — positive wherever finite, no zeros ⇒ single metric;
negative-axis global min 8.933, worst reflection cancellation ~0.15
bit). Every pole is a DOUBLE pole → +inf (both signed zeros, negative
integers, subnormals, −inf; +inf → +0; scipy parity). Zone [1,2) plain
value fit (degree 27, 3 dd-leads — pinned by edge-refined bit-stepped
replay, the origin of that BINDING rule, now in
docs/NUMERICAL-DOCTRINE.md); (0,1) up-step ⊕ dd(1/x²) with deep-tiny
guard (dd 1/x² alone below ~2^-480, +inf at x ≤ 2^-512); down-walk
[2,8); direct-form Bernoulli asymptotic K = 11 (log-free — trigamma
consumes no log_dd); fl(1/x) alone above kTrigammaAsymCut = 2^89.
Reflection π²/sin²(πx) − ψ₁(1−x) on a sinc fit (squaring removes
parity — no cos table); crude-tetragamma lo correction. Gate pinned
1 ULP relative, single metric, full ladder. Bench indicative: 5.7–27
ns/el.

### gamma_p_inv / gamma_q_inv [SHIPPED 2026-08-09]
One pipeline, one bit of orientation: solve against the small side
s = min(p, 1−p) (exact Sterbenz flip — the inverse's complement
transform is on the INPUT). Logit-Newton objective m = lnP − lnQ
(monotone, unbounded, continuous at the median — the solved-side log
saturates and signed ln min(P,Q) jumps at the median; both alternatives
fail measurably). Tri-candidate seeds at ALL a — S1 Temme
normal-quantile, S2 p-form Picard fixed point (owns the tiny-a corner),
S3 far-q-tail fixed point under stability gate L > 3|a−1| — selected
per lane by cheap forward-residual comparison; partition is by (side,
λ-regime), never by a alone. Safeguarded log-residual Newton, 3 steps
(reject residual-increasing, 1/8 backtrack, trust bypass near the
root); deep-small closed form x = exp_dd((LogDd(p) ⊕ lnΓ(1+a))/a) cut
on the dropped-factor bound x₀(1+a) < 2⁻⁶⁰ measured in BOTH
orientations. E dual-form split at a_T (direct below, Stirling above —
direct terms overflow at huge a). Conditioning adjudication: tiny-a
"collapse" lies inside output underflow (κ ≤ ~2^10 wherever x is a
normal double); huge-a beyond-resolution (a ≳ 3e34) has κ → 0 — a test
STRATUM bucketed separately (its trivially-exact rows must not dilute
real-region statistics), not a branch. Oracle: per-row BRACKET
CERTIFICATION at layered dps 60/100 (sign(P−p) flips across the
half-ulp midpoints of xd), deep-small rows certified in log space,
huge-a dual-route via independently-anchored Temme fits; negative
controls exit-2. Gates pinned 1 ULP every bucket,
deep-small/subnormal/x=0 CR, three toolchains identical cells.
Bench indicative: 43 ns/el deep-small path, 520–1056 elsewhere
(structural: 3 seeds + 6 forward evals).

### beta_p_inv / beta_q_inv [SHIPPED 2026-08-10]
The gammainv pipeline generalized; every re-derivation forced by
beta's second parameter is marked [TRANSFER SITE] at its definition
site in src/betainv-inl.h (the load-bearing one: the Newton slope's
(1 − y) factor). Internal frame: input-side Sterbenz flip ⊗ output
orientation swap via I_x(a,b) = 1 − I_{1−x}(b,a); the swap is decided
by a DEFINITIVE median probe (forward at y = 1/2), not by which side
sigma is — the median is near 1 whenever β ≪ α. Log-space forward
reproducing beta's router (R1/R2 never exponentiate; R3 tail in
−cpsi ⊕ ln bracket form — log of an assembled subnormal is QUANTIZED
and flat-spots the residual); PA/PB prefactor split gated by where
PA's dd assembly genuinely runs out (2^40 scale) AND u,v away from −1
(PB's Log1pmxDd degenerates there). FIVE global seed families
(S1 beta-Temme, S2 series inversion, S3 gamma-limit transfer via
gammainv's own seed machinery, S4 exact-B closed form, S5 logit-normal
from exact ψ/ψ₁ moments); candidacy gates are NEVER domain heuristics
— the residual comparison judges. StepsN = 4 (sized by the interior
band where seeds top out at 2–5 bits; quadratic convergence clears the
gate by 20+ bits); residual-uncertainty freeze (beta's condition
number w reaches 2^50 vs gamma's 2^10 — stepping on noise
random-walks). Deep-small cut |1−β|·y/(1+α)·corr(y) < 2⁻⁶⁰, corr(y) =
−ln(1−y)/y (the OTHER side's parameter is the coefficient; both
orientations swept bit-stepped).
**Contract split** (κ = condition number ~1/min(a,b) on the
joint-tiny plateau): y-ULP gate for κ ≤ 2¹⁸ (boundary set by exp_dd's
~2⁻⁷⁰ budget through the forward's ln(1−u), not the design's dd
assumption — measured, adjudicated); above it a BACKWARD-ERROR
contract (forward of returned y within ulp-class of σ; measured 0.000
ulp(σ)). Huge-ν beyond-resolution bucketed by formula from (a,b)
(ν ≥ 1e31; markers carry certification semantics only). Oracle: three
binding constructions beyond the gammainv pattern — fast-path forward
evaluator (100× per-call), gamma-corner hang guard at the enforcement
site, per-row κ-split with backward-error certification where no
y-bracket exists. SHARED-MACHINERY CAVEAT: driving shipped generator
machinery at out-of-domain points can silently misroute
(route_final/small_val_via_cf returns garbage with no exception) —
any future reuse must guard routes. Gates pinned 1 ULP every y-bucket
both sides, B rows 2 (neighbor semantics), backward contract 1.0
ulp(σ), deep-small/subnormal/x=1 CR; full ladder. Bench indicative:
0.59–11.4 µs/el (up to 12 region-routed forward evals; candidate
count is the throughput lever if wanted).

### i0 / i1 / i0e / i1e [SHIPPED 2026-08-11]
Four exports, one TU; even/odd in x so one table serves both signs
(sign reapplied by CopySign at the end — dd Fast2Sum can turn −0 into
+0). Two regimes split at x_s = 8: series in q = x²/4 (exact q via
TwoProd + first-order S′·q_lo correction — q's single rounding under
the series' log-sensitivity (x/2)·I1/I0 was the family's one ratified
design amendment) and tail Chebyshev in 1/x of e^{−x}I_ν(x)√(2πx),
÷√(2πx) in dd with an exact 2^-32/2^16 prescale/postscale (Dekker
ceiling). Unscaled forms via exp_dd mantissa+exponent, saturating at
the EXACT bisected last-finite boundaries. NO log_i0 export (log(i0e)
+ x composes at < 1 ulp relative for x ≳ 2, absolute 3.3e-16
everywhere — documented in ACCURACY with the small-x caveat); NO
recurrence helper (von Mises CDF: Miller backward recurrence from
j_max + ~15 on scaled values normalized by i0e — recipe documented for
libstats #51; forward recurrence is unusable). Non-FMA lesson (now
BINDING in docs/NUMERICAL-DOCTRINE.md): the G1 replay sim originally
modeled FMA-only semantics and only the scaled assemblies — the G4
capped sweep caught 2–3 ULP unscaled non-FMA rows and a
Dekker-ceiling break; klead re-derived by the honest sweep (4/4).
Gates pinned 1 ULP every bucket, all four functions. Bench
indicative: 11.5–19.0 ns/el.

### lbeta [SHIPPED 2026-08-11]
ln B(a,b) as the beta TU's third export: PA's own assembly re-handed,
lbeta = LgammaPosDd(min) − LgammaDiffDd(max, min), one rounding; big
band (min > 2^990) by grouped Stirling difference on 2^-64-prescaled
operands (Binet terms provably unrepresentable, dropped). Domain
a, b > 0 finite, else NaN (SciPy betaln's |Γ| negatives: documented
deviation, no consumer need). Metric: relative where |ln B| ≥ 1,
2^-53-class ABSOLUTE near the zero manifold through (1,1) (lgamma
negative-axis precedent). MEASURED CORRECTLY ROUNDED on every row of
every band (gates 0 ULP / 0.5·2^-53); −inf boundary 44/44 exact (only
the a = b ray reaches −inf — entropy factor H(a/c)). Two kernel
defects caught by tests worth remembering: big-band overflow NaN
through TwoSum's inf−inf error algebra (clamp-and-select with the
empty-sliver proof via ulp(DBL_MAX)); main-band scrub must clamp ONLY
min — clamping max too corrupts every live lane with max > 2^990, and
the smoke symmetry check is blind to symmetric corruption. Oracle
lesson: an ungated hi==0 ∧ lo==0 "agreement" shortcut let ~89-digit
cancellation zero BOTH dps tiers identically — fixed with a
cancellation-free term-magnitude probe scaling dps to the floor.

## GitHub repo settings [applied 2026-07-21 via gh api]
Merge: all three styles, auto-delete head branches (PR merges only —
a direct fast-forward push bypasses it; prune manually). Wiki and
projects DISABLED (docs policy), issues on, discussions off.
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
- 2026-08-15 clang-cl → MSVC ABI item CLOSED, and the install/export path
  exercised for the FIRST TIME. Adoption-spike stage S1 (libstats branch
  `spike/corvus-bessel`): an MSVC 19.51 consumer links the clang-cl-built,
  INSTALLED corvus and runs it correctly — every i0/i1/i0e/i1e row agrees
  with the consumer's OWN std::cyl_bessel_i to ≤ 2.3e-16, both regime paths
  (series below x_s = 8, Chebyshev tail above) exercised. Values were
  checked, not just the link: an ABI fault links cleanly and returns
  garbage, so a clean link alone would have been no evidence. Three
  findings beyond the pass. (1) Install/export had never run before —
  every build tree to date shows hwy_DIR-NOTFOUND, i.e. fetched Highway,
  which disables install BY DESIGN; corvus-config/targets/.pc emit and
  resolve correctly against a clang-cl-built Highway 1.4.0 prefix.
  (2) The dispatch tier is a function of the COMPILER THAT BUILDS CORVUS's
  TUs, NOT of the delivery mechanism. Configs A and B had varied compiler
  and delivery together, so a third de-confounding run was needed:
  FetchContent + clang-cl → AVX3_ZEN4 in 177 s. Full matrix —
  FetchContent+MSVC AVX2, FetchContent+clang-cl AVX3_ZEN4,
  installed-clang-cl + MSVC-consumer AVX3_ZEN4. So the installed path buys
  DECOUPLING (corvus on clang-cl while the consumer stays pinned to
  cl.exe), not tier: a consumer willing to build everything with clang-cl
  gets AVX3 from plain FetchContent with no prefix, and an MSVC-pinned one
  can alternatively flip CORVUS_MSVC_UNBLOCK_AVX512 (measured working,
  gates pass, deliberately unsupported upstream). Secondary datum: clang-cl
  builds corvus 3.9x faster than cl.exe on the same FetchContent config
  (177 s vs 682 s). (3) The exported
  cxx_std_20 requirement travels: the consumer sets no standard and
  std::span still compiles. Config A (all-MSVC FetchContent at the v0.5.0
  tag) also passes, at AVX2, in 682 s wall. The two configs' probe output
  is byte-identical EXCEPT the active_target line, so the A-vs-B adoption
  choice is performance-only, not accuracy.
- 2026-08-14 betainv huge-parameter Dekker audit CLOSED, findings wider
  than the item: (1) the three named DdMulD sites fixed by exact
  prescales (reachable unsaturated, e.g. the Beta(2, 1e307) median on
  non-FMA tiers); (2) NEW — shared BetaR4Tiny carried R1's two huge-B
  hazards unfixed (intermediate overflow near the B·ξ = 8 cap on every
  tier + non-FMA Dekker break; affected the shipped beta forward),
  fixed by R1's twin prescale; (3) NEW — betainv's u = exp(lnf) site
  fed exp_dd an unclamped log: wildly skewed pairs (first parameter
  ≳ 1e69-class, moderate second) NaN'd the orientation probe's
  midpoint forward on every tier and returned the 0.5 fallback; fixed
  by a value-identical argument clamp at the underflow floor (gammainv
  got the same defensive clamp at its latent twin site). Permanent
  gating: 215 R4-huge-B beta rows (harness-verified, 0 declines) +
  203 betainv huge-corner rows (bracket-certified, 86 unprovable rows
  declined); two oracle bugs found and fixed in gen_beta_reference's
  cross-check dps escalation en route. All existing gate cells
  byte-identical pre/post-fix; kappa-bucket comment/code inconsistency
  resolved in the code's favor (35, the shipped replay-validated
  threshold).
- 2026-08-12 quiet-machine Ryzen bench pass DONE (supersedes the
  loaded/indicative v0.3.0-session table as the publishable Ryzen
  set). Protocol: detached self-gating runner (build-clangcl/
  quiet_bench.ps1) refuses to start until ambient < 5% total over two
  consecutive 10 s samples and logs noise between every target; this
  pass gated at 3.6% and held 2.6–6.2% avg throughout (floor ~3%).
  Getting there required killing APSDaemon (~1 core, constant) and
  temporarily stopping WSearch (0.83 core) + LightingService (0.13);
  Performance power plan, clocks 105–108% base. Per-family speedups
  (AVX3_ZEN4; libm-baseline for erf/erfc/lgamma, scalar-walk
  otherwise): erf 3.7–6.4, erfc 2.8–9.5, lgamma 1.4–5.5, erfinv
  5.2–14.9, gamma 8.9–24.0, digamma 7.1–25.2, trigamma 7.6–30.9, beta
  7.3–15.1, gammainv 6.2–142.0, betainv 4.0–55.4, bessel 5.3–9.2,
  lbeta 8.0/12.7. Raw per-band logs: quiet_bench_bench_*.txt (local
  build tree, not checked in; the noise-annotated quiet_bench.log
  alongside is the evidence chain).
- 2026-08-12 v0.4.0 released: tag at 96d181d gated on full CI green
  (run 31652276608, per-job verified); ships the beta-forward fix arc
  and the QC sweep below. Version sweep clean — only the two enforced
  sites carry the number (README badge is dynamic), configure-time
  check re-verified on the bumped tree.
- 2026-08-12 static-analysis sweep COMPLETE (cppcheck 2.21 + clang-tidy
  22.1.3 + lizard 1.22.1, check set pinned in .clang-tidy): 73 raw
  findings triaged to a small hygiene batch — `static` on all
  per-target Impl/dispatch functions and in-TU template kernels (25),
  dead DdNorm removed, two confusable-identifier renames, two
  const-correctness sites, three test const-pointer sites,
  NOLINTNEXTLINE on the CF's intentional floor k/2. Zero
  kernel-arithmetic changes; post-batch: 0 build warnings, 0 tidy
  warnings, 27/27 ctest. ACCEPTED (won't-fix, recorded): bench raw
  fill loops (useStlAlgorithm ×19), aggregate-table
  uninitMemberVarNoCtor ×17 (always brace-initialized), the
  betainv-ulp gate ternary (distinct gates currently equal),
  Case::why documentation member, and the Windows CI D9025 pair —
  Highway's own deliberate /EHs-c- retraction of CMake's default
  exception flags, upstream-intentional and confined to the
  FetchContent MSVC path. CCN measured: max 14 (BetaVec), avg 1.8;
  lizard tripwires CCN 15 / length 400 (tripwires against drift, not
  targets). CORVUS_SANITIZE Open Item also closed: MSVC-ABI configure
  now FATAL_ERRORs instead of letting cl D9002-ignore -fsanitize
  (negative-controlled both directions).
- 2026-08-12 beta-forward u → −1 defect pair FIXED (the PRIORITY item
  open since 2026-08-10, disclosed in every release since v0.1.0):
  closed-form 1+u = c·ξ/α / 1+v = c·y/β on corner lanes in BetaPsiCore
  (be55662; non-corner lanes bit-identical; Log1pmxDd gained the 2-arg
  w-overload + a definition-site u → −1 hazard rule). The fix's own
  786-row pb-corner reference family (d275bcd, --corner-append mode
  with point-bits-digest checkpoint sig) found a THIRD latent defect —
  R1 series recurrence overflow at β ~ 1e307 plus non-FMA
  Dekker-ceiling breaks at l1/l2 and the gammalim t·, all fixed by
  exact power-of-two prescales (68deb66, bit-identical below 2^900).
  Generator cross-check repaired en route (compared raw CF where the
  oracle's routing replaces it; betainc's internal 1−x truncation →
  dps escalated by −log10 x). Both witnesses now CORRECTLY ROUNDED;
  gates re-pinned at the one moved cell (gammalim dir 0 → 1 ULP,
  6/430 rows); full-suite + tier sweep + CI revalidation in this
  arc's commits.
- 2026-08-12 v0.3.0 released: tag at 0bebf95 gated on full CI green
  (run 31561265568, per-job verified); P2 complete — i0/i1/i0e/i1e
  (1 ULP everywhere) + lbeta (correctly rounded everywhere); version
  sweep found corvus.h kVersion* still at 0.1.0 (missed at v0.2.0) —
  fixed + configure-time consistency check added (negative-controlled).
  Major code-writing closed; next phase is analysis/clean-up.
- 2026-08-11 MSVC build-time headroom retro-outlining shipped: log/exp
  HWY_NOINLINE wrappers in beta/lgamma/gammainv -inl.h (11+3+17 call
  sites). Isolated Ninja+cl per-TU times: beta 339.0→96.5 s (3.51×),
  lgamma 193.0→72.1 s, gammainv 75.6→52.1 s, betainv (untouched)
  130.3→95.9 s downstream via smaller lgamma instantiation. All 7 ULP
  gates byte-identical at native AVX3_ZEN4; 4-tier capped sweep green
  with expect-target held; zero new warnings with DEV_WARNINGS=ON.
  NOTE: the historic 547–904/681/486 s priors were MSBuild-batched
  measurements — the Ninja isolated-edge numbers are the comparable
  baseline going forward. CI confirmed on d3e09aa: Windows
  17.4 → 10.6 min, Linux 10.3 → 7.0 min, all jobs green.
- 2026-08-09 gen_beta_data.py kBetaGammaLim standing "frontier review
  owed" flag (2^-49 target deviation) reviewed and RATIFIED: routing
  threshold not truncation depth (end-to-end ULP bar applies), pin
  midpoint target-invariant (log-linear frontiers), empirically
  confirmed by beta's shipped gammalim gates + boundary-crossing seam
  sweeps; flag text retired at the print site.
- 2026-08-10 v0.2.0 released: tag at b0221d1 gated on full CI green
  (run 31449939039, per-job verified); P1 complete — digamma,
  trigamma, gamma_p_inv/gamma_q_inv, beta_p_inv/beta_q_inv; Release
  notes disclose the beta-forward u→−1 known issue.
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
