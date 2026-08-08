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
1. **DIGAMMA SHIPPED [2026-08-08] — next P1 family: inverse
   incomplete gamma/beta** (erfinv's seed + dd-Newton pattern; then
   Bessel I0/I1). Detail design + error budgets = frontier work,
   probe stage first (the digamma probe→G1→G2→G3 pipeline is the
   template). Escalation-density rule [2026-08-06, user]: if
   resolving an escalation spawns a new one more than ~3 deep in a
   chain, the stage defaults back to frontier hands-on work. Digamma
   pipeline ledger (final): probe 0, G1 one (depth 1 → FIRST
   correction), G2 zero, G3 zero (three reviewed-accepted
   deviations).
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
- [WATCH] Ryzen box stability: two GPU-stack bugchecks 2026-07-29 (0x9F
  power-IRP, 0x10E video memory), root-caused [DERIVED] to a
  half-committed NVIDIA driver install; DDU clean reinstall the same
  night. Crash-free since, including the full beta validation ladder
  2026-08-05. Residual watch: whether the NEXT NVIDIA App driver update
  completes cleanly (chronic installer freezes implicate accumulated App
  state; manual driver-only installs are the fallback). Recurrence of
  either bugcheck on the clean stack flips suspicion to VRAM/hardware.
- [OPEN — user action] File the mingw GCC 16.1 AVX-512
  by-value-argument misalignment bug upstream. Repro + bugzilla-ready
  draft: `C:\Users\gdwol\Development\gcc-zmm-mingw-repro\` (repro.cpp,
  repro-struct.cpp, report.md; cite PR 110273 and PR 49001). Blocked on
  GNU bugzilla account approval (pending since 2026-07-29). Technical
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
  AGENTS.md is canonical. Documentation capped at four files. Deferred
  deliberately: LTO/IPO (profile first), shared lib (no demand).
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
