# corvus — Plan / Session State

Detail policy: this file holds STATE — what's decided, what's open, the
next concrete step, and lessons too expensive to re-learn. Full designs,
session narratives, stage ledgers, and measurement play-by-play live in
this file's git history (compacted 2026-08-06, 2026-08-13, 2026-11-02 and
2026-11-27; the 2026-11-02 trim dropped the v1.1.0 order, read it at
`adef4bd`; the 2026-11-27 trim dropped the v1.2.0 order table and the
v1.2.1 per-machine leg narratives, read them at `6d4ddb6`), the GitHub
issues, docs/ACCURACY.md, docs/bench-evidence/, and the kernel/generator
source, which are the official record for finished work. Binding
cross-family engineering rules live in docs/NUMERICAL-DOCTRINE.md, not
here.

## Status [DERIVED] — 2026-11-27

**v1.2.0: released 2026-11-02.** Signed tag `v1.2.0` on `691db43`,
GitHub release published, milestone closed. Rows A1 to D of its issue
order are done; #22 moved to v1.3.0. Nothing is owed for it.

**v1.2.1: released 2026-11-27.** "One value through the span API": a
one-value or short-span call no longer waits on the call before it. No
result moves and the API is unchanged, so it is a patch (VERSIONING.md).
PR #117 merged to `main` as `77bd92b`; signed tag `v1.2.1` on that
commit, signature verified, `check_release_version.py` OK; CI green
there on every job (run 46049888250); GitHub release published;
milestone v1.2.1 closed. Nothing is owed for it but the handshake below.
- Validation, at `6940de9`, the last code commit: Z (clang-cl 22.1.3)
  40/40 on AVX3_ZEN4 and the AVX2, SSE4, SSSE3, SSE2 caps,
  `masked_const` PASS; K (AppleClang 15) 40/40 on AVX2 and the SSE4,
  SSSE3, SSE2 caps; M 40/40 on NEON; I (clang 18.1.3) 40/40 on SSE4
  native with SSSE3 and SSE2 from CI's sweep. `span1` passes on every
  tier. The release head passes 40/40 on K from a fresh `-Werror` build;
  the docs pass's comment edits left all 19 library objects
  byte-identical.
- **Next:** the libstats handshake (Cross-machine handshakes below),
  then milestone v1.3.0.

**What v1.2.1 shipped [DERIVED]** (the published notes say the same)
- What changed, all in how a call reaches a kernel (no kernel
  arithmetic):
  1. A span of one element skips the tail: broadcast, kernel, store one
     lane (`src/driver-inl.h`).
  2. The 27 dispatched `Impl` functions take pointers and one length,
     not `std::span` by value.
  3. `ErfcinvVec`'s and `Log1pVec`'s special-value selects sit behind a
     branch.
  4. A tail of 2 to N-1 elements is assembled in registers, with no
     masked load or store; a tail of one element leaves through
     `GetLane` (`0dfe170`).
- Why: two loads that could not be forwarded from the stores they
  overlapped made each short call wait for the one before it (the masked
  tail's load, and under the Windows x64 convention the wrapper's copy
  of a `std::span`); and a merge-masked load of a constant with an empty
  mask took a slow path while its page was unmapped.
- Figures, `main` against the branch, quiet runs, two build trees:
  - Zen 4, Windows (`docs/bench-evidence/2026-11-26-zen4-v1.2.1/`): one
    value through a span, `erfinv` central 64.7 -> 29.3 ns, `erf` 32.9 ->
    15.7, `exp` 35.8 -> 14.4; `erfcinv` 132.0 -> 29.4 one value and
    14.7 -> 3.4 per element in a batch; `log1p` 3.2x.
  - Off Windows the gain is small and at one value: Kaby Lake AVX2
    0.84x to 1.01x (layout-matched, see #120), Ivy Bridge SSE4 0.89x to
    0.98x, M1 NEON 4 to 10% less per call. On NEON the serialisation
    never existed.
- Guards: `tests/test_span1.cpp` (27 span functions; every length 1 to
  9 at every offset equals the same values in a long span, bit for bit)
  and `tests/test_masked_const.cpp` (Windows: a band reads the same with
  the image's pages trimmed and mapped). The suite is 40 tests.
- No bound moved; no API change; errata: none.
- Known and left for v1.3.0: #120 (Kaby Lake: a load after a gather is
  delayed at one page offset) and #121 (the 2-to-3 element tail after a
  full vector, the full-vector load, `bench_span1`'s caller).
- Consumer note: libstats v3.0.0 pins v1.2.1 [DECIDED, user]. Its scalar
  Gaussian and LogNormal quantiles are one `erfcinv` call each (72 ->
  106 ns on Z when it adopted corvus, OldCrow/libstats#156); the effect
  there is not yet measured.

**v1.2.1: what is not known [OPEN]**
- Which change produces the gain off Windows. All three Unix legs gain a
  little at one value and none separates the changes.
- The masked-load cost is shown only on Zen 4 under Windows, its
  mechanism inferred from on/off experiments. About 1,400 loads of the
  same shape remain in the AVX3_ZEN4 functions; `ErfcinvVec` and
  `Log1pVec` are changed. `tests/bench_cold.cpp` flags nothing on the
  branch, which says only that this link exposes no such constant. A
  static count per kernel is owed (v1.3.0).
- Whether Linux or macOS has an unmapped-page cost: `bench_cold` runs on
  Windows only, and nobody has written a way to unmap image pages
  elsewhere.
- #97 has three leads, none tested against its evidence: the two above
  and #120. I's `erf` and `exp` short-span cells read at two levels
  between processes on both trees (`main` `erf` n = 7: 76.2, and 95.8
  and 84.8 in two of eight).

**Milestones [DECIDED, user, re-cut 2026-11-26].** v1.3.0 is "Scalar
paths and ISA tuning", performance only: #110, #111, #22, #97, #120,
#121. #110 stays as the scalar entry points' API: after v1.2.1 a
length-1 span costs what a broadcast costs, so what is left is a uniform
surface, not speed. The branch `exp/narrow-lane-scalar` (`23398fa`, not
for merge) stays on origin until #110 is done [DECIDED, user,
2026-11-27]: it holds the only copy of the harness that ran `erfcinv` at
one and two lanes (bit-identical, and little faster; figures on #110).
v1.4.0 is "von Mises complement" (#19), a new export
with its own unfreeze. Gated work does not roll between versions: #18,
#20, #21, #28, #29, #44, #95 and #96 sit off every milestone with their
trigger on the issue.

**Fourth machine, Ivy Bridge [DERIVED 2026-11-11]** (MacBook Pro 9,1,
i7-3820QM, Linux Mint 22.3, glibc 2.39). It dispatches `SSE4`, so it is
not a new tier: it adds the non-FMA Dekker path at its real cost, a
glibc baseline, and `SSE4` as a native top tier. An uncapped,
tier-asserted ctest run is one of the release validation legs; the quiet
bench re-runs only when the dd path or the table kernels change.
PERFORMANCE.md §15 and ACCURACY.md's validation-leg line are written.
- [OPEN] No capped-`SSE4` ULP report is published (CI prints AVX2 only),
  so the native leg is not diffed against a capped one. A candidate:
  have CI's sweep print the report per tier, or save one from Kaby Lake.
- [OPEN] The Dekker pre-split is unfiled; file it only if a consumer
  reports pre-Haswell hardware.
- [OPEN] Housekeeping on Kaby Lake: remove the `../corvus-c1` and
  `../corvus-c3` worktrees, and v1.2.1's `../corvus-main`,
  `../corvus-exp` and `../corvus-plan`.

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

**v1.2.0:** none recorded.

**v1.2.1:** none recorded.

## Cross-machine handshakes [OPEN]

- libstats v3.0.0 pins v1.2.1 [DECIDED, user]; libhmm adopts and spikes
  (libhmm #106) against the same release (standards workplan tasks 2 and
  7; fleet coherence rule). v1.2.1 is tagged (2026-11-27), so the pin
  can ride along with libstats' correctness leg that follows its pool
  fix; a sweep with no row moved is the per-tier proof that no result
  changed. libstats has not been told yet [OPEN].
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

## Reach and distribution [discussed 2026-11-03]

Strategy for building a user base (opensource.guide/finding-users). Nothing
is built or filed yet.

**Decided [user]**
- Paper: arXiv first (cs.MS, cross-list math.NA and stat.CO; the user
  already publishes there, no endorsement needed), linked from any site.
  JOSS is nice to have, pursued only once its 2026 scope is met (evidence
  of impact, 6 months of public history, AI disclosure). ACM TOMS is the
  better peer-reviewed venue. Add a Zenodo DOI and CITATION.cff for
  citability.
- C++ first. No per-layer Python bindings (corvus, linalg, rng, libpar
  each bound separately would be a maintenance chore). [DERIVED,
  recommendation] The Python surface follows user concepts: corvus as
  `pylibstats.special`, real NumPy ufuncs (PyUFunc_FromFuncAndData),
  float64 only, GIL released.
- [DERIVED, 2026-11-03] No object-level interop between pylibstats and
  pylibhmm: the peer rule (standards `records/FLEET-DESIGN.md` §2) keeps
  libhmm off libstats, and the hmm+stats bridge (§6) is the only link.
  A Python bridge converts parameters. So separate wheels, each linking
  its own private corvus; the only published corvus Python surface is
  `pylibstats.special`.

**Facts that constrain the plan [DERIVED]**
- GitHub Pages cannot be private on a personal account; refine a site
  locally or in a private repo, publish when ready.
- Tag tarballs come from `git archive`, so `.gitattributes export-ignore`
  decides whether site files ship in releases.
- vcpkg, Conan, conda-forge, Spack and Homebrew build from source, so
  corvus signs no binaries. Homebrew core needs 30 forks, 30 watchers or
  75 stars (225 to self-submit).
- Once libstats links corvus, the pylibstats and pylibhmm wheels contain
  corvus and Highway binaries: NOTICE and Apache-2.0 obligations apply to
  those wheels.

**[OPEN]**
- Whether to build a site at all, and where (Pages from this repo or a
  separate one).
- Accuracy-contract scope for third-party builds (which tiers and
  compilers a port's bound claims cover), Highway version policy for
  ports, MSVC/Windows handling in ports, and the maintenance budget per
  channel.
- Upstream SciPy reports (outward-facing, user's call): erfcinv of the
  smallest subnormal = inf; i0/i1 overflow just below the true
  threshold; gammaln(subnormal) = inf; gammainc with subnormal a;
  polygamma(1, x) costs O(|x|) for x < 0.

**scipy.special comparison [DERIVED, 2026-11-03]** — SciPy 1.18.1 scored
on corvus's reference sets at `b2848d4` with the gate metric, against
corvus's per-function published maximum; about 150 rows cross-checked in
mpmath (no reference off by > 0.5 ULP). Write-up, scripts and outputs:
`../corvus-scratch/scipy-comparison-2026-11-03/` (outside the repo;
REPORT.md). Rerun takes about a minute.
- Parity: exp, log, cos, sin, erf, erfinv, log1p, expm1 (max 1–3 ULP;
  corvus's log1p/expm1 are CR, SciPy's are not).
- SciPy a few ULP worse: i0e/i1e (max 7/11), trigamma (7), erfcinv on
  normal inputs (5).
- SciPy far worse or wrong: erfc tail (moderate-input median 12 ULP, 0
  where the result is subnormal); lgamma and lbeta inf/NaN on subnormal
  arguments; digamma near the negative poles (to 3.5e15 ULP); lbeta
  cancellation (to 1.4e7·2⁻⁵³); gamma_q median 24 ULP; incomplete
  gamma/beta 0 or NaN at extreme parameters; inverses on moderate inputs
  53–59% within corvus's bound (max 1e4–3e5 ULP).
- Limits: the reference sets are adversarial, so full-set percentages
  overstate typical-use error. No speed ratio is supported; SciPy's
  betaincinv averaged ~350 µs/row on these sets, everything else 60 ns–
  5 µs. A speed claim needs a typical-input bench on all three machines.

## Open items

**Open issues** (each carries its trigger and detail):
- #18 erfcx (P2; no filed need).
- #19 von Mises ratio complement 1 − I₁/I₀ in dd (v1.4.0; filed need,
  frontier oracle).
- #20 Hurwitz zeta (conditional on libstats #62).
- #21 exp_dd accuracy bump, betainv κ-horizon 2¹⁸ → 2⁵² (no filed need).
- #22 non-gather x86 variant (v1.3.0; sized by the Kaby C3 leg: the
  table band runs 5.2× slower per element than Zen 4 against 2.9× (zone)
  and 3.6× (Stirling), so a non-gather variant recovers at most 1.5–1.8×
  there; needs a selection mechanism and a new accuracy audit).
  [DERIVED 2026-11-11] Sizing step before any code: cap Kaby Lake to
  `SSE4` and time the table bands against its AVX2 build. The Ivy
  Bridge rows (PERFORMANCE.md §15.3, clang both sides: lgamma table
  0.66×, exp 0.84×, erf 1.13× of Kaby AVX2 per element) support the
  premise but are cross-machine, so silicon is still confounded.
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
- #110 scalar code paths for the span-only functions (v1.3.0; libstats
  calls twelve of them one value at a time).
- #111 `lbeta` and `LgammaDiffDd` cost on the non-FMA tier (v1.3.0;
  three questions, no design).
- #120 on Kaby Lake a load after an AVX2 gather is delayed at page
  offset 0x7c0 to 0x7cf, so a gather kernel's one-value cost depends on
  where the linker puts its 8-byte literals: 10 to 30 ns per call, at
  about 3 of 64 positions (v1.3.0; three questions: other machines,
  which PERFORMANCE.md Kaby figures are affected, a remedy).
- #121 short spans after v1.2.1 (v1.3.0): the 2-to-3 element tail after
  a full vector is 2 to 8 ns slower than the masked one on AVX2; the
  full-vector load still waits on the caller's stores (about 30 ns on
  `erfinv`); `bench_span1`'s write loop is a `memcpy` call on AppleClang.

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

**CI and tooling [2026-11-11]** (rules and incidents: CI House Style §1, §6, §8)
- Runner labels all float on `*-latest`; nothing is pinned.
- The canary runs weekly (Monday 03:45 UTC). A `canary-report` job opens
  a "Canary failed" issue, or comments on the open one (#103).
- Pre-commit is the formatting and hygiene gate (#104): the `precommit` CI
  job runs every hook on all files; the local hook is opt-in per clone
  (`python -m pre_commit install`). Dependabot updates the hook revs.
- `.clang-format` is Google base, 2-space, 80 columns, chosen because the
  code already matched it (about 1,800 changed lines against about 15,900
  for either sibling style). Includes are not sorted (order is
  load-bearing around `hwy/foreach_target.h`); comments are not reflowed.
  `examples/` overrides the indent to 4.
- Coefficient and test-vector tables sit inside `// clang-format off/on`
  (15 pairs, 10 files): clang-format rewrites a braced table to one
  element per line. New tables need the same guard.
- Generated data (`src/*_data.*`), `tests/data/` and
  `docs/bench-evidence/` are excluded from every hook.
- The reformat is `b839820`, listed in `.git-blame-ignore-revs`.

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
- A speed claim about the library comes from a program that links only
  the library, built in its own tree. v1.2.1 began with three wrong
  conclusions from a harness that compiled its own kernel copy, and from
  a restored file whose old timestamp skipped the rebuild
  (OldCrow/standards, PERFORMANCE-MEASUREMENT.md).
- Two trees of the same kernel can differ by 10 to 30 ns a call on Kaby
  Lake because of where each link put one literal (#120). Before
  believing a before-and-after on a gather kernel there, shift one
  tree's constants 64 bytes and see whether the difference survives
  (`docs/bench-evidence/2026-11-26-kaby-v1.2.1/n7/`).
- Read the bench's own disassembly. `bench_span1`'s "write n values"
  loop was a `memcpy` call, so two machines' short-span tables measured
  a caller nobody meant (#121).
- The first cells a process measures read high: 5 to 12% on Kaby Lake,
  and on Linux under `schedutil` a whole first band until the core is
  warmed for 2 s. Warm first, or do not compare a first cell with a
  later one.

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
