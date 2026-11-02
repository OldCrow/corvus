# corvus — Plan / Session State

Detail policy: this file holds STATE — what's decided, what's open, the
next concrete step, and lessons too expensive to re-learn. Full designs,
session narratives, stage ledgers, and measurement play-by-play live in
this file's git history (compacted 2026-08-06, 2026-08-13 and 2026-11-02;
the last trim dropped the v1.1.0 order, the per-row narratives of both
orders and the Resolved log — read them at `adef4bd`), the GitHub issues,
docs/ACCURACY.md, and the kernel/generator source, which are the official
record for finished work. Binding cross-family engineering rules live in
docs/NUMERICAL-DOCTRINE.md, not here.

## Status [DERIVED] — 2026-11-02

**v1.2.0: row D (release tail) at freeze head `2f863be`** (kernels; later
commits are docs only). Milestone v1.2.0 has 0 open issues. All three
machine legs are DONE; bench numbers are in docs/PERFORMANCE.md and the
commit messages of `0a3ee74`, `ac472c3` and `d2a9d20`.
- Zen 4 (2026-11-01): PERFORMANCE.md §4, §5 and §8.2's Zen 4 column
  rewritten from two gated passes plus four repeats; `i1`/`i1e` read with
  `i0`/`i0e`. Evidence `docs/bench-evidence/2026-11-01-zen4-v1.2.0/`.
- Kaby Lake (2026-11-01, `0a3ee74`): `CORVUS_WERROR=ON`, 38/38
  tier-asserted on AVX2 native and SSE4/SSSE3/SSE2 capped, 0 warnings;
  every CI 41348222134 AVX2 report line reproduced; the three no-FMA
  tiers bit-identical. Evidence `docs/bench-evidence/2026-11-01-kaby-v1.2.0/`.
- M1 (2026-11-02, `ac472c3`): 38/38 NEON-asserted, 0 warnings; the ULP
  report equals CI 41348222134's NEON job line for line (639 lines).
  §8.2 has its M1 column. Evidence `docs/bench-evidence/2026-11-01-m1-v1.2.0/`.
- Beta R4 tiny direct reads max 1 ULP (real 0.5241/0.6138) natively on
  all four x86 tiers and NEON, so the 1-ULP bound decided at C3b holds.
- §5 speed modes: Zen 4 shows three bands 1.50× slower in 1–2 of 6
  processes; K and M1 show no 1.50× mode, but the M1 shows seven bands
  stepping 1.07–1.30× between processes. Cause unknown; #97 (v1.3.0)
  owns it with all three machines' data. [DECIDED 2026-11-01, user]
  v1.2.0 does not wait for it.

**Docs pass DONE 2026-11-02:** PERFORMANCE.md and the README performance
section restated at v1.2.0 on all three machines; ACCURACY.md, the
user-facing docs, examples and corvus.h Doxygen audited and fixed; this
file trimmed; RELEASING.md gained the PLAN-trim, docs-pass and
merge-to-main preconditions (5–7). Release-notes draft is in the session
scratchpad; its content is the inputs block below plus the v1.1.0-notes
correction (scalar overloads shipped in v1.1.0, unlisted there).
**Merged 2026-11-02:** version pair 1.2.0 (`2c3e686`; fresh `CORVUS_WERROR`
build, 38/38 tier-asserted on Zen 4); PR #98 green on every job and merged
to `main` as `148bb64`.
**Next:** CI green on `main`, then RELEASING.md from the signed tag (step
11) with the user; then the libstats and libhmm bumps (Cross-machine
handshakes below).
- [OPEN] Housekeeping: remove the `../corvus-c1` and `../corvus-c3`
  worktrees on Kaby Lake and the M1 (A trees of the C1/C3 A/Bs).

**Release-note inputs for v1.2.0 [DERIVED]** (row D):
- New exports: span `expm1` and `log1pmx`, correctly rounded (#45).
- Bound tightenings (MINOR): beta_p/beta_q R4 tiny direct 2 → 1 ULP
  (#94) [DECIDED 2026-11-01, user: a floor for all of 1.x]; lgamma
  table band (5/2, 8) claims 0 ULP outside a 2^-12-ulp tie window
  (`lgamma.table`, #31); lbeta(1, x) = −ln x at 0 ULP (`lbeta.unit`,
  #80 part 2).
- Throughput: #42 lever 1 (beta CF precision split), #31 (lgamma table
  band, single reflection log); per-machine ratios in PERFORMANCE.md
  §12–§14.
- Consumer note: the scalar entry points replace libstats' deducing
  `corvus_scalar` wrapper (the family's names are overload sets).
- Gather-weak x86: v1.2.0 leans harder on gathers (11 per lgamma table
  vector); cost is unmeasured on parts running the Downfall (GDS,
  CVE-2022-40982) mitigation — Skylake through Ice Lake/Tiger Lake — and
  on AMD before Zen 3. #22 (non-gather variant) is on v1.3.0. Secondary
  sources only, enough for one sentence; the README says the same.
- Errata: none recorded for v1.2.0.

## v1.2.0 issue order [DECIDED 2026-10-26, user]

Milestone 7, "Throughput & new exports", on `dev/v1.2.0` (cut from `main`
at `c668080`, the v1.1.0 tag plus the R3 plan commit).

Rules: correctness before bit-changing performance in the same files;
guards written fail-first and landed with their fix; the forward before
its inverse; generators hardened before reference rows are regenerated,
and regeneration done once. Two more:
- **The libstats/libhmm pin bump is the critical path.** What libstats
  adopts in v3.0.0 — #42 lever 1, #45, #31 — gates the release. Anything
  else comes off the milestone if it would hold the release.
- **Guards that protect the performance rewrite land before it.** Lever 1
  is measured against the tightened guards (A4), the stated floors (B1)
  and the marked twins (B2).

| Row | Issue | Done | Commit | CI |
|---|---|---|---|---|
| A1 | #81 tooling chore | 2026-10-26 | `d54ab11` | 40437030819 |
| A2 | #78 fast-math guard | 2026-10-26 | `b01b162` | 40441318445 |
| A3 | #87, #89 generators | 2026-10-27 | `1780304` | (generators only; closed at A4) |
| A4 | #79 guards + one regeneration (#87, #89 closed; #94 filed) | 2026-10-27 | `78c1bbf`, `f9835f0` | 40611128749 |
| B1 | #86 `// FLOOR:` everywhere, lint `--strict` in CI | 2026-10-28 | `537ba3b` | 40782653042 |
| B2 | #77 TWIN markers | 2026-10-28 | `11ab298` | 40788554389 |
| B3 | #80 part 1 (erf j = 0, beta R3 Reinsch form) | 2026-10-28 | `9429cf3` | 40804958773 |
| B4 | #68 Log1pmxDd cut: measured negligible, documented | 2026-10-29 | comment-only | 40804958773 (identical code) |
| C1 | #42 lever 1; K 2026-10-30, M1 2026-10-31; closed 2026-11-01 | 2026-10-29 | `6027e2d`, `0967994` | 40985511199 |
| C2 | #45 span expm1, log1pmx | 2026-10-29 | `19bcff8`, `747dce9` | 40993748302 |
| C3 | #31 + #80 part 2; K 2026-10-30, M1 2026-10-31; closed 2026-11-01 | 2026-10-30 | `5fef537`, `2001fc8`, `040da41`, `5612bdb` | 41149000287 |
| — | Outlining rule applied (lgamma regions; 28 definitions in 17 headers) | 2026-10-30 | `d5c731d`, `ca066a2` | 41157925670 |
| C3b | #94 beta R4 cap corner (BinetDiffDd k = 2 term in dd) | 2026-11-01 | `8a6533f` | 41348222134 |
| C4 | #22 non-gather x86 | MOVED to v1.3.0 2026-11-01 [user] | — | — |
| D | Release tail: PERFORMANCE.md re-run and quiet bench on three machines (DONE); native M1 + Kaby validation at the freeze head (DONE); release notes (inputs above); README line on gather-weak x86 and #22; RELEASING.md | open | — | — |

After D: libstats and libhmm bump to v1.2.0 (standards workplan tasks 2
and 7).

Gated, worked only if the gate opens: #18 (no filed need), #19 (filed
need, P3, frontier oracle; libstats mitigates locally), #20 (opens only
if libstats takes Zipf's closed-form CDF at its v3.1.0 planning), #21 (no
filed need). #44 stays unmilestoned research, answered at C1. #28 and #29
stay off every milestone.

**Decided 2026-10-26 [user]:**
1. Phases A–B run before C1, as ordered above.
2. #89: the kernels keep their tie-to-even result at the subnormal
   exact-midpoint ties of i1/i1e. 1 ULP there is accepted and stated in
   ACCURACY.md. A3 fixes the reference rows and the generators' tie rule
   only.
3. #80 is split: part 1 in B3, part 2 after #31 (C3).
4. #22 is conditional: it comes off the milestone if it would hold the
   release.
5. [2026-10-30] #94 is row C3b, diagnosis on the milestone and the fix
   conditional as above. #95 (R3 k = 0 compensated sum) is off every
   milestone.

**Milestones re-cut 2026-11-01 [DECIDED, user].** v1.2.0 is rows C3b and
D. v1.3.0 holds #19, #22 and #97. Gated work no longer rolls between
versions: #18, #20, #21, #28, #29, #44, #95 and #96 sit off every
milestone with their trigger stated on the issue.

## Errata by release (RELEASING.md step 14)

**v1.1.0 errata (release notes, RELEASING step 14) [DERIVED 2026-10-25].**
Each was never met by a release and never enforced by a gate; each has
an "Erratum" paragraph in docs/ACCURACY.md naming its witness.
- gamma_p/gamma_q R4 complement: 0 → 1 ULP (row 12).
- beta_p_inv/beta_q_inv deep-small closed form and subnormal/x = 1 results:
  0 → 1 ULP (row 9).
- lbeta: 0 → 1 ULP relative, 2^-54 → 2^-53 absolute (#76, row 12j).
- beta R2 complement: 0 → 1 ULP (#76 M-9).
- trigamma on (0, 1): CR → 1 ULP; the gate was already 1 (#76 M-8).
- gamma_p_inv/gamma_q_inv deep-small closed form and subnormal/zero results:
  CR → 1 ULP; the gate was already 1 (#76 M-6, L-1).
Disclosed as v1.0.x defects but fixed, bound kept (PATCH-class, list
under notable fixes): the i0/i1/i0e/i1e tail (#52, row 10); erfc's top
subnormal binade (#74, dd G chosen over an erratum). [DERIVED]
Not errata (checked against v1.0.0/v1.0.1, which published 1 ULP; only
dev/v1.1.0 had tightened them to 0): gamma R2 complement (#56 tightened
it; row 12 mislabelled the restatement) and beta's gamma-limit slice
(#54 tightened it; #90 restates it). Nothing to disclose. [DERIVED]

**v1.2.0:** none recorded (see the release-note inputs under Status).

## Cross-machine handshakes [OPEN]

- libstats and libhmm bump their corvus pin to v1.2.0 (standards
  workplan tasks 2 and 7; fleet coherence rule: libhmm adopts and spikes,
  libhmm #106, against the same release as libstats).
- libstats #156 (open): its v2.5.0 throughput gate waits on v1.2.0.
  #42's NEON targets (gamma_p ≤ 40, beta_p ≤ 30 ns/el) are out of reach
  under the accuracy contract (gamma's prefactor alone is 29 ns/el on 8
  lanes); told to libstats on #156 (2026-10-29). Re-baselining its gate is
  libstats' decision. libstats also exempts or re-baselines its log
  speedup gate (corvus log stays ~0.7× scalar libm, Decisions).
- libstats #107/#108 and libhmm #99 and the kTrigDMax half of libhmm
  #101 (all open) close against each consumer's corvus swap, citing the
  v0.8.0 consumer cross-validation record on corvus #32.
- libstats #62 (open) decides #20: Hurwitz zeta opens only if libstats
  takes Zipf's closed-form CDF.

## Open items

**Open issues** (each carries its trigger and detail):
- #18 erfcx (P2; no filed need).
- #19 von Mises ratio complement 1 − I₁/I₀ in dd (v1.3.0; filed need,
  frontier oracle).
- #20 Hurwitz zeta (conditional on libstats #62).
- #21 exp_dd accuracy bump, betainv κ-horizon 2¹⁸ → 2⁵² (no filed need).
- #22 non-gather x86 variant (v1.3.0; sized by the Kaby C3 leg: the
  table band runs 5.2× slower per element than Zen 4 against 2.9× (zone)
  and 3.6× (Stirling), so a non-gather variant recovers at most 1.5–1.8×
  there; needs a selection mechanism and a new accuracy audit).
- #28 upstream Highway HWY_BROKEN_MSVC version floor (needs Highway's own
  suite passing under MSVC/AVX-512 first).
- #29 watch GCC PR 126741; re-qualify mingw GCC for AVX2+ when it lands.
- #44 promotion candidates from libstats (the quadrature question is
  answered no, Decisions).
- #95 beta R3 compensated k = 0 S-sum (halves the S error, replay max
  0.69 → 0.33 bracket ulp; costs R3 throughput; the bound holds).
- #96 untraced 1–3% cost on the lgamma zone and lbeta main band after C3,
  same sign on three machines; a record, not worth a session now.
- #97 cause of the process-to-process speed modes (v1.3.0).

**Unfiled, still open** (recorded only in closed #82 or here):
- [OPEN] Four kernel routes neither reached nor argued unreachable, each
  with a `// ROUTE-OPEN:` comment: gammainv/newton/floor,
  gammainv/bracket/cap, betainv/seed/none, betainv/special/nonfinite_log.
  Deciding one means proving and deleting the branch (a kernel change,
  so a fix row) or keeping it with a written guard rationale.
- [OPEN] Run the route census in CI; it needs a census build leg (about
  one extra configure).
- [OPEN] Oracle gaps: `beta_inv_root` fails at extreme tiny pairs (a, b ~
  1e-94/1e-189, root within e^-1e131 of 1); beta inverses off the median
  with shapes ≥ 1e100 (gigabytes per row; the cross stratum stops at
  1e15); both shapes tiny (I_x flat; the two forms agreed on a wrong
  root, so the cross stratum excludes max(a, b) ≤ 1e-20).
- [OPEN] exp's ExpDdVec still rounds twice into subnormals: its
  documented 0.75-ulp design bound. A candidate only.
- [OPEN] LogDd(Dd) pays a Div even for x.lo = 0; a double-only LogDdAny
  would speed every consumer's fall-back. Touches lgamma/gamma/beta.
- [OPEN] Rare-path regions are not audited for plain HWY_NOINLINE (the
  2026-10-30 outlining audit covered hot regions); MSVC runtime cost of
  the outlining split is unmeasured.
- [OPEN] Size drift to watch, not refactor (2026-10-10 sweep): BetaVec
  (511 lines) and BetaInvForward (413) exceed the 400-line tripwire;
  BesselSeriesS (10) and GammaInvConsider (11) exceed 8 parameters.
- [OPEN] AGENTS.md is 10.4 KB against the ~4 KB core budget (Decisions);
  trim or re-budget.
- [OPEN] `gen_betainv_reference.status()` appends to an old scratchpad
  file outside the repo.
- [DEFERRED until the CMake minimum reaches 3.28, user 2026-10-11] Add
  `EXCLUDE_FROM_ALL` to the Highway `FetchContent_Declare` so MSVC stops
  building `hwy_test`/`hwy_list_targets` (CI minutes, half the D9025
  warnings; benign). Do not raise the minimum (3.25) for this alone; then
  check `hwy`/`hwy_contrib` still build and install gating is unchanged.
- [WATCH, 2026-07-29] Ryzen box: two GPU-stack bugchecks traced to a
  half-committed NVIDIA driver install, clean since the DDU reinstall.
  Recurrence on the clean stack points at VRAM/hardware.
- [ILLUSTRATIVE] Possible future consumers: a C++ port of
  multi-agent_sim (batch distance/trig), zeekhmm training pipelines.

## Decisions

**Project shape**
- Name corvus, namespace `corvus::`. Scope: statistical special functions
  plus the elementary family (exp/expm1/log/log1p/log1pmx/cos/sin).
- Backend: Highway behind `src/ops-inl.h`; public API std-only.
  std::simd migration = facade reimplementation, deferred until
  implementations mature (mid-2026: GCC 16 partial, no libc++).
- Dependencies: find_package(hwy) preferred, FetchContent fallback pinned
  to the audited 1.4.0 (bump only with revalidation). Static lib, PIC,
  MIT, clean-room only. Deferred deliberately: LTO/IPO (profile first),
  shared lib (no demand).
- Install: `cmake --install` requires a system Highway
  (find_package(hwy 1.4)); FetchContent builds are build-tree only; no
  bundling, no nested install (keeps source tags free of the NOTICE
  obligation). Distribution is source-only (RELEASING.md). Every v-tag
  gets a GitHub Release object; protect-tags makes v* tags immutable.
- Platform tiers: Tier 1 (accuracy-audited on silicon) = NEON (M1, CI;
  covers NEON_BF16, verified on CI 41149000287), AVX-512/AVX2/SSE family
  (Zen 4 native + capping), AVX2 (Kaby Lake). Tier 2 (compiles,
  unaudited) = SVE and anything else Highway emits.
- Doc policy [2026-08-09, user]: compact always-read AGENTS.md core
  (~4 KB) plus on-demand references carved by when they are needed
  (ENVIRONMENT.md, NUMERICAL-DOCTRINE.md). Resist growth in the core.
- Comment archaeology test [ratified 2026-08-13]: a comment survives only
  if it states a constraint, bound, derivation or hazard the code cannot
  show; who found it, when, and what it measured before the fix go.
  Rewrite, don't delete — keep each site's math half.
- Conventions, the FP-contraction rule, the Windows toolchain rule, the
  outlining split, public naming (#34) and the four gating lists:
  AGENTS.md. Bounds are contractual floors, errata clause: VERSIONING.md.
  Twin-drift, FLOOR and oracle-trust rules: NUMERICAL-DOCTRINE.md.
- Reassess the outlining split entirely if MSVC becomes properly
  supportable under Highway, AVX-512 dispatch included [user] (AGENTS.md
  points here).

**Process**
- Milestones: gated work stays off every milestone with its trigger on
  the issue; it does not roll between versions [2026-11-01, user].
- Acceptance rule for performance work: bounds, not bits [2026-10-29,
  user]. Every ULP gate green on every tier, no published bound moves,
  AND per region the real maximum on every runnable tier rises by no more
  than tie noise (0.01 ulp); a routed hunt at the conditioning corners;
  the error mechanism stated at the definition site. A rise is
  stop-and-ask, not a re-pin. Lane-mix determinism is unchanged.
- A same-sign sub-spread bench reading gets an issue only when a further
  machine reproduces it [2026-10-30 rule; applied as #96].
- Oracle memory guard [2026-10-25, user]: `refgen_common.cap_process_memory`
  caps each oracle worker (default 1024 MB, `CORVUS_ORACLE_WORKER_MB`) and
  `core.plan_workers` sizes pools to free memory minus MIN_FREE_MB. Ad-hoc
  sweeps use the same initializer, never a bare Pool.
- Pipeline template (new families): probe → design → G1 generator/data →
  G2 oracle/references → G3 kernel + tests → G4 gate pinning + tier
  ladder → G5 docs. Takeover from delegated work is a judgment call
  (churn or local-fix myopia), not a counter [2026-08-08, user].
- CI [2026-10-11, user]: CI workflow changes land on `main` first
  (`workflow_dispatch` works only once the trigger is on `main`), then
  cherry-pick to the dev branch; dev branches without a PR run CI via
  `gh workflow run ci.yml --ref <branch>`. The Windows job builds
  clang-cl capped at AVX2 and tier-asserted besides MSVC. The sanitizer
  step builds Debug at -O1 [2026-10-25, user]. The macOS job builds native
  (gates unasserted, logs `FEAT_BF16`) and `(HWY_NEON-1)`-capped
  (asserted NEON gates and the ULP report).
- Release split precedent [2026-10-20]: a correctness release is measured
  and tagged before bit-changing performance work, so each speedup's
  accuracy check is one diff against a tagged, guarded baseline.

**Performance decisions that close a question**
- log stays at ~0.7× scalar libm [2026-10-11, user]. A correctly rounded
  log on this structure cannot reach 1.5× libm (exact r, the error-free
  assembly and the rounding test are ~45 of ~60 ops; the fall-back's cost
  forbids a looser window). The 129-slot log table is closed (~6% at best).
  A non-CR opt-in log entry point was considered and not taken.
- #43 [2026-10-10, user]: no libstats exp/log import (it would loosen a
  validated bound, MAJOR); corvus got its own CR fast paths instead. What
  did NOT help, so nobody re-tries it: paired-load gather for erf
  (Highway's NEON GatherIndex is fine); Estrin for log's series; a chunked
  driver keeping the fall-back call out of the hot loop, or inlining the
  fall-back (the cost was the rounding test's own ~10 ops); HWY_LIKELY on
  the skip branch. DriveUnary runs two vectors per iteration on 2-lane
  targets only.
- #42 [2026-10-29, user]: lever 1 is gated per lane (deep CF steps in
  double, last 16 in dd, only where F < 16; F reaches ~2^43 in R2, so a
  fixed split errs by ~F·2^-53). Lever 2 (per-lane late start) dropped:
  gamma R2 needs all 44 steps below x ~ 5. #44 answered no: a dd
  quadrature node costs ~10 ns/el on AVX3_ZEN4, 24 nodes ~251 ns/el
  against gamma R2's whole CF at 26.
- #31: the Estrin zone leg was retired [2026-08-31] — the Horner is not
  the zone's dominant cost on any machine.
- Broadcast (constant-argument) overloads: libstats fills spans on its
  side [2026-09-29, user]. File one only if libstats measures the fill.
- libhmm consumes transcendentals in-register; corvus offers no
  register-level entry point by design (it would contradict the std-only
  public surface). No corvus-side work proposed.

## Lessons too expensive to re-learn

- Reference sets prove the grid, not the domain: every 2026-10-10 hunt
  defect sat in a band the generators never sampled, often because the
  generator shared the kernel's formula (#54) or capped its own
  self-check (#61).
- One clean oracle row never refutes a band claim; re-check across the
  band. A derivation comment is a claim to test, not evidence (#52, #56
  bounds were never achievable by their method).
- Any huge-N oracle must budget the loggamma cancellation or check at
  doubled dps (#65 was a dps-60 oracle artefact).
- Store probe oracle roots with `mp.nstr(v, 40)`, never `str(mpf)`
  (15 digits read as 10–30 ulp "failures" in #88's sweep).
- A/B builds restored with Copy-Item keep the old mtime, so ninja skips
  the rebuild; hash-check both binaries (C1's first A/B measured a stale
  build).
- An integer-ULP gate fails only when real error passes B + 0.5; real-error
  gates per region and tier class (`tests/data/real_pins.txt`) are what
  see a margin erode (#82 T1).

## Shipped families — what would be expensive to re-derive

Measured bounds, per-tier cells and method: docs/ACCURACY.md. Math:
kernel derivation blocks at the definition sites. Below: design choices
and bugs worth not re-learning.

### erf + erfc
erf: table + local Taylor (clean-room port of libstats vector_erf_neon
through the facade). erfc: core reuses the erf table via compensated
1 ∓ erf assembly; tail e^{−a²}·G(1/a)/a on exp_dd. The shared erf series
must meet erfc's RELATIVE precision near a = 6, not erf's absolute one.
The MulSub/no-FMA exact-residual hazard created ops::SquareLow's
capability guard.

### exp_dd + log_dd
- exp_dd returns mantissa + exponent (ExpDdFrac + ScaleTwo): consumers
  fold factors in first and power-of-two scaling lands LAST, which keeps
  subnormal results at one rounding. Every later family uses it.
- log_dd's mantissa is centred on 1; the slots adjacent to m = 1 carry
  R = 1, L = 0 exactly — no special case at x = 1.
- An exact dd pair is not necessarily NORMALIZED: (p−1, p_lo) had |lo|
  up to 2^-53 against hi ~ 2^-8; fix is one TwoSum (not Fast2Sum — r
  crosses zero). The gate caught it; reading the code did not.

### lgamma
Masked fixed-step DOWNWARD recurrence (x − k exact; down lands every lane
in one polynomial), Stirling at X0 = 8, sinpi reflection via
Γ(1−x) = −x·Γ(−x) (−x exact where 1 − x is not). π/|sinπx| is never
formed (overflows at poles); sinc factoring cancels π analytically.
lgamma(1) needs an explicit +0.0. ops::ProdLow's non-FMA Dekker split
overflows for |a| > 2^996; Stirling groups x·(log x − 1) − log(x)/2,
never (x − ½)·log x. Since #31: a table band on (5/2, 8) (h = 1/32,
degree 8, 15.6 KB, 11 gathers per vector) and a single reflection log.

### erfinv / erfcinv
Central direct dd-lead polynomial (κ ≈ 1: a Newton would pass erf's
error straight through); tail seed + ONE dd Halley step, log-space in
the far tail. Four bugs, all "an identity exact in general is wrong at
one distinguished point": ErfcCoreDd clamp at exactly 6.0; (−0)+(+0) = +0
needing a trailing CopySign; `Neg(z − 1)` vs `Sub(1, z)` at z = 1; the
oracle's own guess at s = 1. HalleyMid's possibly-NaN x0 into a gather is
the origin of the masked-lane scrub rule.

### Incomplete gamma P/Q
Region map (λ = x/a, a_T = 20): R1 series-P; R2 CF-Q BACKWARD fixed depth
44 (forward Lentz was 24 ULP worst and false-converges); R3 Temme on
λ ∈ [½, 2] (clean-room c_k(η) by Vandermonde solve in 1/a); R4 small-a
Q-direct via Expm1Dd. Small side computed directly, 1 ⊖ dd complements.
Freeze masks must SELECT the accumulator, never add zero (DdAddD
renormalizes). Oracle traps: compute the small side directly in the
oracle too; mpmath's lower gammainc hangs near the ridge for large a —
exact-asymptotic oracle above a = 1e4.

### Incomplete beta P/Q
Architecture, routing order, prefactor identities and pinned constants
are documented in src/beta-inl.h and src/beta_data.h. Six oracle defect
classes against zero kernel defects in adjudication — the origin of the
oracle-trust doctrine: (a) truncation at ambient dps (mp.mpf() on an mpf
re-rounds; carry exact complements, never recompute); (b) component
budgets voided by downstream cancellation; (c) mp.quad returns √eps noise
when the integral sits below working epsilon (normalize by peak
log-magnitude); (d) false saturation certificates from collapsed
complements; (e) mpmath betainc returns internally consistent garbage in
the gamma-limit corner (layered dps cannot catch it); (f) small-side-direct
applies to the oracle too. tools/verify_beta_reference.py is the
independent harness, with its adjudicated rows baked in as negative
controls. Never default-construct vector-member structs (implicit ctors
instantiate outside the per-target attribute region; bf16 CI).
Tooling: set mp.dps inside every computation layer (worker,
subprocess-send, subprocess-parse); probe with exact hex, never
display-rounded values; Windows multiprocessing needs `__main__` guards.

### digamma / trigamma
digamma: relative metric where |ψ| ≥ 1, 2^-53-class ABSOLUTE inside the
per-zero bands (near-relative at the negative zeros would need
2^-104-class fits — rejected on cost). trigamma: all-relative (no zeros);
every pole is a double pole → +inf; direct-form Bernoulli asymptotic
(log-free); fl(1/x) alone above 2^89. The edge-refined bit-stepped replay
rule (NUMERICAL-DOCTRINE.md) comes from trigamma's zone fit.

### gamma_p_inv / gamma_q_inv
Solve against the small side s = min(p, 1−p) (exact Sterbenz flip on the
INPUT). Logit-Newton objective m = lnP − lnQ (monotone, continuous at the
median; both alternatives fail measurably). Tri-candidate seeds at all a,
selected per lane by forward residual; partition by (side, λ-regime),
never by a alone. Since #48: a bracketed Newton phase in log x (rtsafe
rule) between seed selection and the three safeguarded steps. Huge-a
beyond-resolution (a ≳ 3e34, κ → 0) is a test stratum, not a branch.
Oracle: per-row bracket certification at layered dps.

### beta_p_inv / beta_q_inv
The gammainv pipeline generalized; each re-derivation is marked
[TRANSFER SITE] in src/betainv-inl.h. The output orientation swap is
decided by a definitive median probe (forward at y = 1/2), with a logit
tie band (#55). Log-space forward reproducing beta's router (log of an
assembled subnormal is quantized and flat-spots the residual). Five seed
families; candidacy is judged by residual, never domain heuristics.
Contract split at κ = 2¹⁸ (exp_dd's budget through ln(1−u)): y-ULP below,
backward error above, gated non-circularly by
`tests/data/betainv_backward_reference.txt` (#60). Driving shipped
generator machinery at out-of-domain points can silently misroute; any
reuse must guard routes.

### i0 / i1 / i0e / i1e
Even/odd in x, one table both signs (CopySign last — dd Fast2Sum can turn
−0 into +0). Series in q = x²/4 below 8 (exact q via TwoProd + first-order
q_lo correction), tail Chebyshev in 1/x above, ÷√(2πx) in dd with exact
prescale. No log_i0 export (log(i0e) + x composes at < 1 ulp for x ≳ 2);
no recurrence helper (von Mises CDF: Miller backward recurrence on
scaled values, recipe given to libstats #51). The non-FMA replay rule
(NUMERICAL-DOCTRINE.md) comes from this family's G4 sweep.

### lbeta
lbeta = LgammaPosDd(min) − LgammaDiffDd(max, min), one rounding; big band
(min > 2^990) by grouped Stirling difference on prescaled operands.
Relative where |ln B| ≥ 1, absolute near the zero manifold through
(1, 1). Since #80 part 2: lbeta(1, x) = −ln x from one log behind an
AllFalse skip. Bugs worth remembering: big-band overflow NaN through
TwoSum's inf−inf algebra; the main-band scrub must clamp ONLY min (the
smoke symmetry check is blind to symmetric corruption); an ungated
hi == 0 ∧ lo == 0 oracle "agreement" shortcut let cancellation zero both
dps tiers identically.
