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
**Tagged 2026-11-02:** CI green on `main` at `148bb64`; signed tag
`v1.2.0` on `691db43` (docs only after `148bb64`), signature verified,
`check_release_version.py` OK.
**Released 2026-11-02:** GitHub release v1.2.0 published; milestone
v1.2.0 closed. corvus has nothing owed for v1.2.0.
**Next:** libstats' and libhmm's PLAN.md carry the bump note
(Cross-machine handshakes below); each library picks up its bump in its
own session. corvus's next work is v1.2.1 (one value through the span
API; draft PR #117), then milestone v1.3.0 (#110, #111, #22, #97); both
under "Milestones re-cut 2026-11-26" below.
**Fourth machine, Ivy Bridge [DERIVED 2026-11-11]** (MacBook Pro 9,1,
i7-3820QM, Linux Mint 22.3, glibc 2.39; #101 holds the indicative GCC
pass, #100 the `quiet_bench.sh` Linux bugs it found):
- It dispatches `SSE4` (Highway has no AVX1 target), so it is not a new
  tier. It adds three things: the non-FMA Dekker path at its real cost,
  the first glibc baseline, and `SSE4` as a native top tier, so an
  uncapped build runs the real dispatch decision on an AVX-without-AVX2
  CPU.
- Workflow: a fleet-table row in ENVIRONMENT.md; an uncapped,
  tier-asserted ctest run joins the release validation legs; the quiet
  bench re-runs only when the dd path or the table kernels change, not
  at every release.
- Accuracy: no new table or column. The SSE4 target compiles without FMA
  on any host, so the capped Kaby, Zen 4 and CI Linux/GCC legs already
  run this code path. One validation-leg line in ACCURACY.md once a
  tier-asserted gate run and ULP report from this machine are recorded.
- Performance [DERIVED 2026-11-11, DONE]: PERFORMANCE.md §15, §10.1 and
  the Ivy Bridge column of §11, from two clang 18.1.3 passes at `2565e85`
  (gates 1.00/3.00% and 3.00/1.00%; pass 2 quoted). Evidence
  `docs/bench-evidence/2026-11-11-ivybridge-v1.2.0/`. The README gained
  the non-FMA sentence beside the gather-weak line: the dd bands cost
  3.55× (lgamma zone) and 2.65× (erfc tail) Kaby Lake per element where
  lane count predicts 2× (clang against clang; #101's GCC pass read
  3.75×/2.75×). The three highest ratios, 3.7–4.1×, are §5's three bands
  (i0 series, digamma and trigamma up-step); cause not known.
- No unfreeze. The pass shows no defect; both code candidates (#22, and
  pre-split constant operands for Dekker on non-FMA tiers) are
  speed-only. The libstats/libhmm bump to v1.2.0 proceeds unchanged.
- [DERIVED 2026-11-12, DONE] Order of work, all on the Mint machine: fix
  #100 (Linux side done 2026-11-11), then the recordable #101 run (done
  2026-11-11, row above), then the ctest + ULP report (done 2026-11-11:
  38/38 and 18/18 report programs, `SSE4`-asserted, clang 18.1.3 and GCC
  13.3 byte-identical; ACCURACY.md has the validation-leg line). #100's
  Darwin branch confirmed on the M1 2026-11-12 (below).
- [OPEN] No capped-`SSE4` ULP report is published (CI prints AVX2 only),
  so the native leg is not diffed against a capped one. A candidate:
  have CI's sweep print the report per tier, or save one from Kaby Lake.
- [OPEN] The Dekker pre-split is unfiled; file it only if a consumer
  reports pre-Haswell hardware.
- [OPEN] Housekeeping: remove the `../corvus-c1` and `../corvus-c3`
  worktrees on Kaby Lake (A trees of the C1/C3 A/Bs). M1: already gone
  (checked 2026-11-27); its `../corvus-main` v1.2.1 worktree goes once
  the M line here is committed.

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

**Milestones re-cut 2026-11-26 [DECIDED, user].** v1.3.0 is "Scalar
paths and ISA tuning", performance only: #110, #111, #22 and #97. #19
moves to v1.4.0 ("von Mises complement"), a new export with its own
unfreeze. The off-milestone list above is unchanged.
- **v1.2.1 is in review [DECIDED, user, 2026-11-26]: one value through
  the span API.** Draft PR #117, branch `dev/v1.2.1`. No result moves
  and the API is unchanged, so it is a patch (VERSIONING.md). Five
  changes, each found by measuring through the shipped library; CI is
  green on `653758c` (Linux, Windows, macOS):
  1. A span of one element skips the tail (`src/driver-inl.h`).
  2. The dispatched `Impl` functions take pointers and one length, not
     `std::span` by value (all 27).
  3. `ErfcinvVec`'s special-value selects sit behind a branch.
  4. `Log1pVec`'s likewise: its band read 3.2x slower with the image's
     pages unmapped (15.8 against 4.9 ns per value).
  5. The tail of 2 to N-1 elements is assembled in registers, with no
     masked load or store: a 2-element `erfinv` span cost 62 ns against
     31 for one element, and now costs 34. corvus is a public library,
     so a short span is an ordinary call whether or not libstats makes
     it [user, 2026-11-26].
- **What was wrong [DERIVED, Z, 2026-11-26, indicative: machine in
  use].** Two loads that could not be forwarded from the stores they
  overlap made each one-value call wait for the one before it: the
  masked tail's load from a slot the caller had just written, and the
  wrappers' 16-byte copy of a `std::span` (the Windows x64 convention).
  Either alone leaves the call at the kernel's latency. Separately,
  `erfcinv`'s unconditional selects compile to merge-masked loads of
  their constants; with an all-false mask the load never maps its page
  and costs about 90 ns a vector while it is unmapped, which every
  library-only link hit. Library-only bench, `main` against the branch,
  ns per one-value call: `erfcinv` central 131 -> 30, probit shape
  164 -> 69, `erfinv` central 65 -> 29; `erfcinv` central batch 14.6 ->
  3.4 per element.
- **This replaces two earlier notes here, both wrong.** The first said
  a length-1 fast path alone saved 33 ns: an experiment
  (`exp/narrow-lane-scalar`) measured that in a test program with its
  own kernel copy and raw pointers, and the library did not reproduce
  it. Before that, a note proposed a branchy scalar path on the premise
  that a one-value call ran every region; the driver has padded a short
  span with its live element since #42. The rules this produced are in
  `docs/NUMERICAL-DOCTRINE.md` (Kernel construction),
  `docs/PERFORMANCE.md` §6 and, as a proposed fleet standard,
  OldCrow/standards `PERFORMANCE-MEASUREMENT.md`.
- **Validated so far.** Z, clang-cl 22.1.3: ctest 40/40 on AVX3_ZEN4; the
  capped sweep on AVX2, SSE4, SSSE3 and SSE2. Guards:
  `tests/test_span1.cpp` (27 span functions; every span length from 1
  to 9, at every offset, equals the same values in a long span, bit for
  bit; also on the AVX2 cap) and
  `tests/test_masked_const.cpp` (Windows: a band reads the same with
  the image's pages trimmed and mapped; 4.0x on `main`, 1.00x on the
  branch; central and three tail bands). `tests/bench_cold.cpp` runs
  the same method over 41 bands of every span function; it found
  `log1p` and flags nothing on the branch.
- **Before the tag [OPEN].**
  - The fleet legs below, on K, M and I.
  - Z's quiet run is done [DERIVED, 2026-11-26]: `bench_span1` and
    `bench_cold`, `main` against the branch, on PR #117 as `64629ed`
    (`docs/bench-evidence/2026-11-26-zen4-v1.2.1/`). It confirms the
    indicative figures within about 1 ns; PERFORMANCE.md §6 cites it.
  - Before merging #117: the comments in `src/driver-inl.h` and
    `src/erfinv-inl.h` marked their figures "indicative" (`log-inl.h`
    never did). Corrected in the one-element-tail commit (M line
    below), which is the code change they were waiting for.
  - **Re-gate `0dfe170`** (the one-element tail; the legs validated
    `653758c`, and the driver is compiled per tier). No re-bench
    anywhere: the only new path is a one-element tail, measured with
    this exact diff on both 2-lane machines, and no K or Z cell meets
    one (n = 5 or 9 on four lanes, 9 on eight); every other path is
    byte-identical to what they timed.
    - **K owes:** AVX2 native and the SSE4/SSSE3/SSE2 caps, each
      tier-asserted, `CORVUS_WERROR=ON`; `span1` PASS on all four.
    - **Z: DONE 2026-11-27 at `6940de9`** (which contains `0dfe170`),
      clang-cl 22.1.3. AVX3_ZEN4 native with `CORVUS_WERROR=ON`: a full
      rebuild with no warning, 40/40 tier-asserted, `span1` PASS on all
      27, `masked_const` PASS on its five bands (ratios 0.97 to 1.01).
      `tools/sweep_tiers.ps1`: AVX2, SSE4, SSSE3 and SSE2 all pass,
      each tier-asserted on every gate run, `span1` PASS and
      `masked_const` at 1.00 on each. No bench, as directed.
    - I: done (pass 3 gated the identical diff on SSE4 native; its
      SSSE3/SSE2 caps are CI's Linux sweep). M: done (40/40 NEON).
    - CI on the push: green at `6940de9` (Linux sweep, Windows MSVC +
      clang-cl, macOS NEON, pre-commit), read 2026-11-27.
    Report each in the fleet-legs lines below.
  - Version, release notes, RELEASING.md.
- **v1.2.1 fleet legs [OPEN; one per machine: K (AVX2, macOS), M (NEON,
  macOS), I (SSE4, Linux)].** Session-start check first (AGENTS.md).
  Report in the status line below; commit and push only when the user
  asks.
  1. **Gate, native.** `git fetch && git switch dev/v1.2.1` (head
     `653758c` or later). Release build, then
     `CORVUS_EXPECT_TARGET=<tier> ctest --test-dir build
     --output-on-failure`: 40 of 40, with the tier printed. Two of the
     40 are new: `span1` must print PASS for all 27 functions;
     `masked_const` prints "not run: Windows only" and passes.
  2. **Gate, capped.** The capped sweep for every tier below the native
     one, as for v1.2.0 (docs/ENVIRONMENT.md).
  3. **Bench, quiet, before and after.** This answers an open question:
     whether Linux and macOS gain anything. Build `main` in a second
     checkout (`git worktree add ../corvus-main main`), copy
     `tests/bench_span1.cpp` into it and add
     `add_executable(bench_span1 bench_span1.cpp)` and
     `corvus_test_target(bench_span1)` to its `tests/CMakeLists.txt`.
     Then on a quiet machine, alternating, three times each:
     `tools/quiet_bench.sh -b <build> -t <tier> bench_span1`. Report
     the `local` column and the short-span table (n = 1, 2, 3, 4, 7)
     for both trees. Never compare by swapping files in one tree
     (OldCrow/standards, PERFORMANCE-MEASUREMENT.md, rule 2).
  4. **Not asked:** `bench_cold`. It runs on Windows only. Whether the
     unmapped-page cost exists on Linux or macOS needs a method for
     unmapping image pages there, which nobody has written.
  - K: DONE 2026-11-26 at `653758c`, on PR #117 as `ba23fac` and
    `02944a2` (`docs/bench-evidence/2026-11-26-kaby-v1.2.1/`). Gates:
    40/40 on AVX2 native and on SSE4, SSSE3 and SSE2 capped, each
    tier-asserted, `CORVUS_WERROR=ON`, 0 warnings; `span1` PASS on all
    four. Bench, six runs gated at 2.7–4.5% with the user's background
    processes stopped around each, median of three, `main` -> branch:
    `local` 0.57x (`erf`) to 0.94x (`erfinv` central); short spans 0.64x
    to 0.98x. So macOS on AVX2 gains; which change does it is not
    separated. `main`'s runs spread up to 1.21x on short-span cells, so
    a ratio within about 10% of 1 is not established.
    - One cell is slower: `erfinv` central at n = 7, 111.7 -> 118.2 ns
      per call. An on/off build shows the register-assembled tail
      (change 5) costs it about 6.5 ns when the tail follows a full
      vector, and `exp` 2.5 ns; the same tail gains 5 ns on `erf` and 42
      to 47 ns on `erfcinv` there, and 7 to 36 ns on all four at n = 2
      and 3. [DECIDED, user, 2026-11-27] Accepted as recorded; the
      driver is unchanged for v1.2.1.
    - [DERIVED, K, 2026-11-27; `.../2026-11-26-kaby-v1.2.1/n7/`, on PR
      #117 as `a954181`]
      The sign does not differ by kernel: with the layout matched the
      register tail is 2 to 8 ns slower at n = 7 on all four kernels and
      6 to 12 ns faster at n = 2 and 3. `bench_span1` read `erf` and
      `erfcinv` the other way because of a CPU effect, shown by a
      standalone program: on this Kaby Lake (microcode 0xF8) a load that
      follows an AVX2 gather and touches page offset 0x7c0..0x7cf is
      delayed. Where the linker puts the kernels' 8-byte literals decides
      whether a gather kernel pays 10 to 30 ns per call; a 64-byte shift
      turns it on or off, and the rule predicts 64 of 64 positions.
      `main`'s bench binary had a shared literal there and the branch's
      did not, so K's `local` ratios for the gather kernels are
      overstated: layout-matched they are 0.84x (`erf`) to 1.01x
      (`erfinv` central), not 0.57x to 0.94x. Still a gain.
    - [OPEN, from the above] (1) Run `n7/gather_alias2.cpp` on Z and any
      other AVX2 machine: is there a window, and where. (2) Which Kaby
      Lake figures in docs/PERFORMANCE.md were taken with a literal in
      the window. (3) A remedy is not designed: constants loaded from a
      table whose alignment corvus controls, or #22's non-gather
      variant. (4) `bench_span1`'s "write n values" loop compiles to a
      `memcpy` call on AppleClang, so its short-span table does not
      measure 8-byte stores. (5) The full-vector load still waits for the
      caller's stores at n >= N (about 30 ns on `erfinv`); v1.2.1 fixes
      only the tail. None of these blocks v1.2.1; all are v1.3.0
      candidates for the user to file or drop.
    - For the hold on M's one-element-tail diff: K's result is in. The
      n = 7 cost is the register tail's own work (on AVX2 a blend loop
      and a `memcpy` call) and K found no reason to change the 2-to-3
      element tail. K did not build M's diff or measure a one-element
      tail (n = 5 or 9 on AVX2). M's guess that the blend chain and the
      copy loop explain the sign by kernel is not what K found: there is
      no sign difference once the layout is matched.
    - The plan's `git worktree add ../corvus-main main` fails where
      `main` is already checked out; a detached worktree of
      `origin/main` works.
  - I: DONE 2026-11-27 at `653758c`, on PR #117 as `954f852`
    (`docs/bench-evidence/2026-11-27-ivybridge-v1.2.1/`). Gates: 40/40
    on SSE4 native and on SSSE3 and SSE2 capped, each tier-asserted,
    clang 18.1.3, `CORVUS_WERROR=ON`, 0 warnings; `span1` PASS on all
    three. Bench, two passes, every run gated at 0 to 4%; pass 2 quoted
    (each bench pinned to one core warmed for 2 s, eight alternations),
    median, `main` -> branch: `local` 0.89x (`erf`) to 0.98x (`erfinv`);
    n = 1 calls 0.89x to 0.95x; n = 2 to 7 0.96x to 1.03x. So Linux on
    SSE4 gains a little at one value and about nothing on a longer
    span; which change does it is not separated.
    - Three cells are slower, by 1.1 to 1.7 ns per call: `erfinv`
      central at n = 3 (172.4 -> 173.9) and `erf` at n = 3 (52.0 ->
      53.7) in all eight pairs, `erf` at n = 7 (76.2 -> 77.3) in six.
      K's slower cell, `erfinv` central at n = 7, is level here (312.7
      -> 312.8).
    - They are the one-element tail, and M's saved diff removes the
      cost [DERIVED, I, 2026-11-27; `pass3/`, on PR #117 as `80da7e3`].
      Four trees, eight alternating processes each, pinned and warmed,
      gated at 0 to 3%, `main` / branch / branch with K's masked tail /
      branch with M's diff: `erfinv` central n = 3 172.4 / 173.9 /
      172.4 / 172.2; `erf` n = 3 52.0 / 53.8 / 51.8 / 52.2; `erf` n = 7
      76.2 / 77.3 / 74.7 / 75.5; `exp` n = 7 76.0 / 78.2 / 74.4 / 74.9.
      n = 1, 2 and 4 and the long table read the same on the branch and
      both experiment trees (within 0.6 ns). The patched tree gates
      40/40 on SSE4, `span1` PASS. So the two 2-lane legs agree on the
      defect and the fix. The masked one-element tail is no slower on
      SSE4: 0.5 to 0.8 ns under M's diff on `erf` and `exp` at n = 7.
      Nothing here bears on K's n = 7 cell. Whether M's diff goes on
      the branch stays with the user [OPEN], as on the M line.
    - Pass 1 (the plan's protocol, three alternations) read its first
      band slow on both trees (batch 51 to 57 ns against 37): under
      `schedutil` the bench starts on a core still coming up from idle.
      Warming the core first removes it. A quiet Linux run of a bench
      this short needs that warm-up.
    - `erf` and `exp` short-span cells read at two levels from process
      to process on both trees (`main` `erf` n = 7: 76.2, and 95.8 and
      84.8 in two of eight). Cause not known; a lead for #97.
  - M: DONE 2026-11-26/27 at `653758c` (evidence
    `docs/bench-evidence/2026-11-26-m1-v1.2.1/` and
    `2026-11-27-m1-v1.2.1-tail1/`, on the M1, not yet committed).
    Gate: 40/40 NEON-asserted, `span1` PASS on all 27, `masked_const`
    not run (Windows only); no tier below NEON to cap on arm64. Bench:
    `bench_span1`, `main` `c1d6640` in a worktree against the branch,
    three quiet passes each in 35 min with the user's indexing, media
    and sync processes suspended (libstats' SIGSTOP method, approved by
    the user; gates 3.2–3.8%). Medians of three: (a) the serialisation
    is absent here — on `main` the `local` column is within 2–8% of
    `array` on every band, where Zen 4 read 2x; (b) a span of one costs
    4–10% less per call (erf/exp −1.1 ns, erfinv −1.8, probit −2.3),
    `local` erf/exp 9.7/9.5 → 8.4/8.3, erfcinv central 30.0 → 28.1,
    probit 59.5 → 56.0, lgamma 41.5 → 40.3; (c) erfcinv central batch
    12.6 → 11.1 (change 3, the selects behind a branch); (d) the
    register-built tail costs +1.5–2.3 ns per call on erf/exp at n = 3
    and 7 (+6–9%, repeatable to 0.1 ns), noise on the inverses. On
    NEON that tail is always ONE element (every odd length), a case
    K's n = 2, 3 and 7 cells never meet (their tails are 2–3 elements;
    a 1-element tail there is n = 5 or 9, unmeasured). Second sitting,
    00:49–01:10: `StoreTail` with m = 1 writes `GetLane(v)`, so a
    trailing single element takes the length-1 path; gate 40/40,
    `span1` PASS; the cost is gone (erf n = 3: `main` 18.8, branch
    20.3, patched 19.1; exp n = 7: 29.8 / 31.5 / 29.5; the inverses
    0.4–3.9 ns under `main`), n = 1, 2, 4 and the long table unchanged
    within noise. Eight lines in `src/driver-inl.h`, uncommitted on the
    M1 (diff beside the evidence). It changes no cell K or Z measured
    and does not touch the 2–3-element tail the 2026-11-27 decision
    kept, but it IS a driver change after "the driver is unchanged for
    v1.2.1". **[DECIDED, user, 2026-11-27] Adopted for v1.2.1.** I
    reproduced the three slower cells as this case and measured the
    identical diff level with `main` on SSE4 (I line above); K found no
    sign by kernel once the layout was matched and no reason to change
    the 2–3-element tail (K line above), so M's earlier guess that the
    blend chain and copy loop explained K's sign is withdrawn — they
    are the register tail's cost, the sign was layout. The commit on
    `dev/v1.2.1` (`StoreTail` m = 1 through `GetLane`, eight lines, plus
    the "indicative" comment corrections in `src/driver-inl.h` and
    `src/erfinv-inl.h` quoting Z's quiet run) is `0dfe170` on PR #117;
    gated 40/40 NEON, `span1` PASS. What it owes is the "Re-gate
    `0dfe170`" item under "Before the tag" above (K and Z; no
    re-bench). The 2–3-element tail after a full vector (K's
    n = 7, 2–8 ns) is the v1.3.0 tail design, with K's items 1–5.
    Caveat on M's short-span table (K's item 4): AppleClang compiled
    the bench's write loop to `memcpy`, so n >= 2 cells on M and K did
    not write separate 8-byte stores; I's pass 3 is the clean reading
    and agrees. erfinv's two batch cells read in modes (§5); not
    settled by three passes.
- **Not known [OPEN].** Which change produces the gain off Windows.
  All three Unix legs gain a little at one value (K, I, M above) and
  none separates the five changes; on NEON the `local`/`array` ratio on
  `main` shows the serialisation never existed there. The
  System V convention passes a span in registers, and the masked-load
  cost is shown only on Zen 4 under Windows, its mechanism inferred
  from on/off experiments. About 1,400 loads of the same shape remain
  in the AVX3_ZEN4 functions; `ErfcinvVec` and `Log1pVec` are changed.
  A clean survey says only that this link exposes no such constant, so
  a static count of these loads per kernel is still owed (v1.3.0). #97 has
  both leads, neither tested against its evidence.
- **libstats v3.0.0 pins v1.2.1** [DECIDED, user]: it stays on the v1.2
  line. Tag before libstats' correctness leg that follows its pool fix,
  so the new pin rides along; a sweep with no row moved is the per-tier
  proof that no result changed. libstats' scalar Gaussian and LogNormal
  quantiles are one `erfcinv` call each (72 -> 106 ns on Z when it
  adopted corvus, OldCrow/libstats#156); the effect there is not yet
  measured.
- **#110 stays on v1.3.0** as the scalar entry points' API. After
  v1.2.1 a length-1 span costs what a broadcast costs, so what is left
  is a uniform surface, not speed.
- #112 closed 2026-11-26 (`b8c401f`): the Windows runner paired MSVC
  19.52, whose standard library requires Clang 22, with a standalone
  LLVM 20.1.8. The clang-cl job now uses Visual Studio's bundled
  clang-cl (22.1.3 on that image).

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
- #100 `quiet_bench.sh` on Linux (ps listing, shebang, machine line);
  `tools/` only, no release needed. [DERIVED 2026-11-11] Fixed and tested
  on the Mint machine (failing gate lists five consumers, target lines
  carry `top:`, machine line names CPU and distribution, runs as
  `./tools/quiet_bench.sh`; now POSIX `#!/bin/sh`, tested under dash and
  `bash --posix`). Darwin branch confirmed on the M1 2026-11-12 (macOS
  27.2.1, `/bin/sh` = bash 3.2 POSIX mode, script at `38026b6`): runs
  as `./tools/quiet_bench.sh`, machine line unchanged, forced gate
  failure lists five consumers, target lines carry `top:`, `-t NEON`
  passes and `-t AVX2` fails rc=2. Closes.
- #101 Ivy Bridge SSE4/glibc machine in PERFORMANCE.md: written
  2026-11-11 (Status above); closes with its PR.

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
