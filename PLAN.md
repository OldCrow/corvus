# corvus — Plan / Session State

## Status [DERIVED] — 2026-07-30 (Ryzen)
**Incomplete beta detail design COMPLETE (section below); next step is
G1 — spawn the Sonnet tooling sub-agent on gen_beta_data.py.** CI went
green on the hoist commit fbefba0 with test_dd_special confirmed running
on all three runners; zero crashes/WHEA on this box since the driver
remediation (stability watch continues).

Previous session (2026-07-29 evening): dd_special hoist landed.
Log1pmxDd/Expm1Dd hoisted out of gamma-inl.h into `src/dd_special-inl.h`
(+ `src/dd_special_data.h` via the new `tools/gen_dd_special_data.py`,
self-checks carried over as (a)/(b)); the gate renamed
test_gamma_util → test_dd_special and moved to its dependency slot right
after log_dd in all four lists; the reference file renamed VERBATIM
(generation stays in gen_gamma_reference.py — shared seeded rng stream, see
AGENTS.md). Validated on this box under clang-cl: native AVX3_ZEN4 gamma and
lgamma ULP outputs byte-identical to the pre-hoist baseline, dd_special
output byte-identical to the old gamma_util output, 13/13 native ctest, and
the full AVX2→SSE2 capped sweep green (`build-cap-cc`). AVX3_DL/AVX3 capped
re-runs were skipped this session (machine stability, open item below) —
optional due diligence for the next Ryzen session given the byte-identity.
Beta API resolved with user 2026-07-29: **beta_p / beta_q (a, b, x, out)**,
mirroring the gamma pair; hoist-first sequencing also user-confirmed.

**Phase C part 2 (incomplete gamma P/Q) SHIPPED and fully validated.**
gamma_p/gamma_q public; max 2 ULP direct side everywhere (4 ULP on one
complement corner), gates pinned to measured. Validated identical cell for
cell on every fleet tier: AVX2/SSE4/SSSE3/SSE2 (Kaby Lake native + capped),
NEON (CI run 30414669451; reproduced on the fleet M1 Mini 2026-07-29 at
efa98e1 — clean build, 13/13 gates, measured ULP values identical cell for
cell), and AVX-512 (Ryzen 2026-07-28: AVX3_ZEN4 native
+ AVX3_DL/AVX3 capped + the four lower tiers re-swept, all under clang-cl —
mingw GCC 16.1 is disqualified at AVX-512 by a misaligned-zmm-spill codegen
bug found this session; see ACCURACY.md and AGENTS.md). Built by the
two-agent split decided 2026-07-26 (Sonnet tooling, Opus kernel) with
orchestrator review gates between stages — the workflow notes at the end of
the Phase C part 2 section record what the review gates caught.

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
1. **Start here — resume point 2026-08-02. Everything lives on branch
   `beta/g3-kernel`** (pushed; CI is main/PR-only so the branch is
   silent; merge to main only after G4 goes green). G1+G2 shipped on
   main; G3 kernel + the G1/G2 revision are branch commits. Read the
   dated subsections of the beta design section in order — the SIXTH
   routing correction and kBetaGammaLim = 2⁵⁹ subsection is the
   latest binding state. Resume sequence:
   a. **Tooling revision cycle 2** (Sonnet agent): in
      tools/gen_beta_data.py implement the sixth correction in
      route_final (REVERT the fifth/λ≥0 rule; model the near-one
      post-route: R1-evaluated value > 1 − 2⁻¹¹ → R2-CF opposite
      orientation), swap check (b)'s lattice (vi) for the re-route
      band (viii), re-prove check (e) with the post-route, emit
      kBetaGammaLim = 2⁵⁹ (existing tables MUST stay bit-identical —
      diff the header); in tools/gen_beta_reference.py add the
      betainc-with-timeout rescue for the 5,978 small-τ-guard drops,
      regenerate references, re-audit histogram/coverage under the
      new routing. Values in the current checked-in reference files
      are routing-independent and valid as a checkpoint.
   b. **Kernel revision** (Opus agent — FRESH spawn; brief it from
      the PLAN subsections G3 + G1/G2-revision + sixth-correction,
      NOT from memory): the committed kernel is at FOURTH-correction
      state. Changes: implement the sixth-correction post-route
      (mask update between the R1 and CF passes); add the
      gamma-limit slice (max CF-param ≥ kBetaGammaLim → t_dd =
      −β·log1p(−ξ) dd, E = α·LogDdAny(t) ⊖ t ⊖ LgammaPosDd(α), core
      via gamma-inl.h's GammaCfRecip/GammaSeriesSum templates —
      instantiate only those two); DELETE the kRefDefectCutoff hatch
      in tests/test_beta_ulp.cpp; rebuild (build-clangcl, vcvars
      import pattern in AGENTS/old notes), re-run smoke + ULP vs the
      regenerated references; expectations: direct ≤ 3–4 ULP/region,
      R2-gammalim reported as its own row.
   c. **G4**: pin gates to measured (R2-gammalim row included; decide
      R4 depth 48 vs 56 from the widened-window data; C_lg 256→128
      only if measured hot), Kaby Lake native+caps → NEON CI → Ryzen
      AVX3 native + capped sweep (sweep_tiers.ps1 -CxxCompiler
      clang-cl), then merge branch → main. **G5**: four-list audit,
      ACCURACY.md + README in the same change set, PLAN close-out.
   Sub-agent brief rules that earned their place: state the MECHANISM
   "never set run_in_background on any tool call" (three agents
   parked on the concept phrasing); mp.dps inside every layer;
   multiprocessing child-detection via current_process().name.
   Old-session scratchpad (agent NOTES.md, oracle checkpoints, probe
   scripts — still on disk, path-stable): C:\Users\gdwol\AppData\
   Local\Temp\claude\C--Users-gdwol-Development-corvus\
   4bada365-0af6-4c65-8090-e31c82fb549a\scratchpad\beta\
   The design already bakes in the lessons list (day-one HWY_NOINLINE
   on cores AND driver, masked-lane scrubs incl. the CF, freeze masks
   select-not-add, small-side rule in the ORACLE too, mpmath ceiling →
   exact-asymptotic fallback, fixed lengths proven at boundaries by
   generator self-checks, four-list registration, sub-agent briefs
   finish-and-report).
2. Quiet-machine bench pass on Kaby Lake before publishing its performance
   numbers (its gamma bench was session-loaded, labeled indicative). The
   Ryzen gamma bench IS quiet-machine and at 7b52ed1 (2026-07-29,
   AVX3_ZEN4, clang-cl): simd ns/el R1 64–68, R2 56–58, R3 49–56,
   R4 67–74 (8.9–21× the scalar-walk baseline, an upper bound per its
   caveat). Identical within noise to the pre-outlining build — the
   HWY_NOINLINE driver costs nothing measurable at 8 lanes.

## Open Items
- [OPEN — stability watch only] **Ryzen box: two kernel crashes within
  ~45 min on 2026-07-29 evening; remediated the same night.** #1 ~20:55:
  bugcheck 0x9F DRIVER_POWER_STATE_FAILURE (param1=0x3, a driver blocked
  a power IRP) after Modern Standby churn while a corvus build ran
  unattended — the build was collateral, not cause (user-mode code cannot
  block power IRPs). #2 21:25: bugcheck 0x10E
  VIDEO_MEMORY_MANAGEMENT_INTERNAL (param1=0x2D) during ordinary
  interactive use. Both GPU-stack domain; NO WHEA events either time.
  Cause [DERIVED]: the user installed NVIDIA driver 32.0.16.1088 that
  evening and the installer hung per a months-long pattern (NVIDIA App
  error dialog, frozen until reboot), leaving a half-committed driver
  stack; the driver store held SIX coexisting nvpcf.inf generations
  (Optimus power management), four nvhda, three nvppc.
  **Remediation DONE 2026-07-29 late evening**: NVIDIA App state wiped,
  DDU safe-mode clean (NVIDIA only; AMD iGPU untouched), clean reinstall
  of 32.0.16.1088 — and the installer COMPLETED WITHOUT FREEZING for the
  first time in months, implicating the accumulated App/NvContainer state
  in the chronic install failures. NVIDIA App reinstalled by choice
  (game-settings use). Verified after: driver active and healthy, zero
  crashes/WHEA since, store down to one current generation per component
  plus four INERT 04/2025 leftovers (nvpcf/nvhda/nvppc/nvam — not
  referenced by any device; optional prune via elevated
  `pnputil /delete-driver oemNN.inf`, no /force). Remaining watch: a few
  days crash-free confirms the driver theory; recurrence of EITHER
  bugcheck on the clean stack flips suspicion to VRAM/HW (stress test,
  ASUS support). Watch also whether the NEXT NVIDIA App driver update
  completes cleanly — if the freeze returns, the App is implicated and
  driver-only manual installs are the fallback. Until a few clean days
  accumulate, treat long unattended sweeps here as slightly at-risk.
- [OPEN — repro ready, user action to file] **File the mingw GCC AVX-512
  by-value-argument misalignment bug upstream.** Root cause pinned
  2026-07-29 with a 60-line freestanding repro: GCC 16.1
  (x86_64-w64-mingw32) stores the ms_abi invisible-reference temporaries
  for 512-bit by-value arguments/returns with aligned vmovapd but never
  64-byte-aligns them (prologue is a bare sub; ABI guarantees 16) — one
  legal rsp residue in four is safe. NOT a regalloc-spill bug: real
  spills are correctly vmovupd (the gamma TU's 620), the 120 aligned
  accesses are argument/return temps, and HWY_NOINLINE outlining is what
  created them. Crashes at every -O level incl. -O0; clang-cl and MSVC
  handle the identical source; -mstackrealign/-mpreferred-stack-boundary=6
  change nothing. Deliverables in
  `C:\Users\gdwol\Development\gcc-zmm-mingw-repro\` (repro.cpp,
  repro-struct.cpp, report.md — bugzilla-ready draft; reference
  PR 110273, the i686 sibling, and PR 49001). Blocked only on the user's
  GCC bugzilla account, which is pending GNU's manual account approval
  (anti-spam) as of 2026-07-29 — no further agent work until it clears.
- [RESOLVED 2026-07-29] **Windows CI MSVC-codegen blowup — fixed for
  gamma (7b52ed1) and erfinv (1202273); CI Windows job 22m48s → 6m55s,
  green.** Durable rule now in AGENTS.md (Architecture): outline region
  cores AND the driver of every heavy family TU from day one; keep small
  hot helpers inline. Full history below. History: 5.5 min typical → 16 min with erfinv
  (watch item) → >25 min (killed) on all four gamma pushes 2026-07-28/29,
  including 03e80c9 which outlined the four region cores but left GammaVec
  inlining LgammaPosDd, Log1pmxDd, ExpDdFrac and the erfc core into both
  export loops (MSVC's optimizer is superlinear in function size; runner
  image and MSVC version verified unchanged — drift ruled out). Fix
  7b52ed1 (Kaby Lake): GammaVec HWY_NOINLINE + MSVC-only
  /d2ReducedOptimizeHugeFunctions on gamma.cpp. Confirmed: its CI run is
  green on all three jobs, Windows at 22m48s of the 25-min ceiling.
  Bit-identity verified on AVX2 (Kaby Lake) and on this box's full
  re-validation at 7b52ed1 (native AVX3_ZEN4 table byte-identical, all
  six capped tiers green). [DERIVED] Local MSVC timing, Ryzen Ninja+cl
  19.51, pre-fix source: gamma.cpp ~15 min of codegen, erfinv.cpp the
  long pole (finished between 15:18 and 22:24 total) — on the 2-core CI
  runner either TU alone exceeds the budget, so the pre-fix timeouts were
  deterministic, and erfinv.cpp is now what keeps the job near 23 min.
  Post-fix re-time on the same box, same recipe: gamma.cpp ~15 min →
  ~2:16 (7×), library total 22:24 → 18:13 with erfinv.cpp the last
  ~13 min of it. **erfinv treated the same way on the Ryzen box
  2026-07-29** (ErfcInvCore + both drivers HWY_NOINLINE; ErfinvCentral
  deliberately left inline — outlining it cost a measured 0.5–0.7 ns/el
  on the all-central fast path for negligible codegen weight, and the
  final config's only cost is ~0.3 ns/el on that path from the driver
  call, everything else flat): library total 18:13 → **3:59**, erfinv.cpp
  from ~17 min to under 4 (lgamma.cpp is now the last TU to finish).
  Validated before push: native AVX3_ZEN4 gates 13/13 with the erfinv
  ULP table byte-identical, full clang-cl sweep AVX2→SSE2 + AVX3_DL/AVX3
  green, MSVC-built AVX2 gates green. No /d2 escape hatch needed for
  erfinv.cpp. CLOSED by CI run 30422030773 (2026-07-29): all three jobs
  green, Windows **6m55s** (25+ timeout → 22:48 → 6:55), and the Windows
  ULP report step confirmed emitting the measured values (identical to
  every other validated platform). Incomplete beta outlines its driver
  from day one per the AGENTS.md rule.
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
  **Options evaluated 2026-07-25 (M1), two bit-identical skip guards
  shipped in lgamma-inl.h:** an all-zone fast path in LgammaLow (when
  every lane sits in [1/2, 5/2] the recurrence multiplies P by exact
  ones and the log runs on P = 1, whose slots carry R = 1, L = 0
  exactly — both contribute exactly zero) and a monotone early break
  out of the recurrence walk. ULP printout byte-identical before and
  after. Zone on NEON: 41.2 → 19.6 ns/el (0.29x → 0.61x); other
  regions unchanged. Width-independent — expect Ryzen zone to lift
  too (SSE2 0.46x toward parity, AVX3 2.60x toward ~4x); re-measure
  next Ryzen session. Remaining zone cost is the degree-34
  select-Horner chain itself. **Interval splitting (option A) is
  DEFERRED, trigger recorded:** revisit only if profiling incomplete
  gamma P/Q end-to-end (its prefactor consumes lgamma) shows the zone
  Horner as a real bottleneck; it would roughly halve zone again but
  costs a generator rework, per-coefficient select-vs-gather redesign
  whose trade flips per tier, and fleet revalidation. Recurrence
  region (0.22x NEON, worst) has no lever: its cost is genuine work,
  and X0 = 8 is accuracy-forced. Docs positioning note (2-lane
  targets favor scalar libm for lgamma; corvus wins from 4 lanes up)
  rides the next docs pass.
  **Ryzen re-measure done 2026-07-25 (loaded, indicative),
  post-fast-path:** ULP stats byte-identical on AVX3_ZEN4 native AND
  all four capped tiers (FMA row 30/5 and no-FMA row 29/3 both
  unchanged), so the bit-identity claim now holds on every validated
  x86 tier, not just NEON. Zone AVX3_ZEN4: 2.60x -> 5.6–6.5x
  (5.8 ns/el — every 8-lane vector in the zone range qualifies for
  the skip). Zone SSE2: 0.46–0.50x -> 0.98–1.07x — parity at 2 lanes,
  as predicted. Recurrence/Stirling/mixed/reflection unchanged within
  noise on both. Updated 2-lane picture: zone ~1.0x, recurrence
  ~0.6x, Stirling ~0.7x, reflection ~0.8x vs mingw libm.
  **Kaby Lake 4-lane measurement 2026-07-26 (loaded, indicative) — the
  pending docs claim "corvus wins from 4 lanes up" does not survive it.**
  AVX2 native on the i7-7820HQ vs Apple libm: zone 1.02x, Stirling
  0.77–0.90x, mixed 0.57x, reflection 0.49–0.51x, recurrence 0.28–0.31x.
  Only the zone reaches parity at 4 lanes, so the positioning note needs
  rewording before the docs pass — the honest form is per-region and
  per-libm, not a lane-count threshold. Not a regression and not a new
  defect: it is the lane economics already recorded, divided by a fast
  vendor libm (Apple lgamma 12.9–25.0 ns/el here, the same 2–3x-over-
  mingw advantage seen on the M1), on 2017 silicon. The zone fast path
  is confirmed working at a third width (2 lanes 0.61x NEON / ~1.0x
  SSE2, 4 lanes 1.02x, 8 lanes 5.6–6.5x). [OPEN] Settle the wording
  against a same-libm comparison before publishing any lgamma
  positioning claim.
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

### Phase C part 2 — regularized incomplete gamma P/Q detail design
### [resolved 2026-07-26; SHIPPED 2026-07-28 — see the shipped-record
### subsection at the end for deltas and findings]

Probe-validated on the Kaby Lake box (mpmath oracle + python-float double
replay, the erfc-tail methodology; scripts were session-scratchpad
`probe1/3c/4/4c/5*.py` [ILLUSTRATIVE — generators must re-derive and
self-check every number below anyway]).

**Public API**: `gamma_p(a, x, out)` / `gamma_q(a, x, out)`, elementwise
spans, same length. Specials (SciPy limits): x<0 or a<0 → NaN; NaN
propagates; (x=0, a>0) → P=0, Q=1; (a=0, x>0) → P=1, Q=0; (0,0) → NaN;
x=+inf → P=1 (a finite); a=+inf → Q=1 (x finite); (+inf,+inf) → NaN. No
domain caps; saturation handled by an E<−800 underflow mask.

**Region map** (λ = x/a, a_T = 20):
- R1 series-P: {a<20, 0<x≤a+1} ∪ {a≥20, λ≤½}. Fixed cap N=64 (worst need
  52 at a→20⁻, x=a+1; ≤58 at λ=½ for any a), per-lane freeze mask
  t < s.hi·2⁻⁶⁰ (monotone), vector break on AllFalse. Terms double, sum dd.
  Replay: ≤3 ULP direct-P over the whole region incl. a=1e-300.
- R4 small-a Q-direct: {0<a≤3/2, 0<x≤4}:
  Γ(a,x) = [(Γ(1+a)−1) − (x^a−1)]/a − x^a·Σ_{n≥1}(−x)ⁿ/(n!(a+n));
  Q = a·Γ(a,x)/Γ(1+a). Γ(1+a)−1 = Expm1Dd of lgamma's OWN zone poly at
  EXACT argument (t=a centre-1 for a≤½; t=a−1 centre-2 for a≤3/2 — never
  form 1+a); x^a−1 = Expm1Dd(a·LogDd x). Σ: dd terms (static dd 1/n table
  n≤36, weights DdRecipDd(TwoSum(a,n))), alternating, cap 36 + freeze
  [cap was 30 until 2026-07-27: the 4^n/n! tail at x=4 needs 34 terms
  for 2^-58 — generator self-check catch].
  Replay: correctly rounded on the full probe grid incl. a=1e-300 and the
  x=e^−γ cancellation line. This region is what closes the small-a
  complement corner: Q ~ a·E1(x) keeps full relative accuracy.
- R2 CF-Q: {a<20, x>a+1} ∪ {a≥20, λ≥2}, minus R4. **Backward** fixed-depth
  N=44, no convergence test: K = x+2N+1−a; for j=N..1:
  K = (x+2j−1−a) − j(j−a)/K; Q = e^E·recip(K). Backward contracts rounding
  (replay ≤4 ULP worst vs 24 for forward Lentz; forward also
  false-converges — (30,133) stopped at N=8 with 70 ULP).
  [N was 40 until 2026-07-27: the generator's sup sweep found the true
  worst point at a→1.5⁺, x=a+1 — a blind spot in probe1's a-grid, which
  skipped (1,2). N=42 meets 2^-56; 44 carries margin.]
- R3 Temme: {a≥20, ½<λ<2}, direct side by sign(x−a):
  smaller = ½erfc(z) ± R, z = √(aφ), R = e^{−aφ}/√(2πa)·S(η,1/a),
  S = Σ_{k<11} c_k(η)/a^k. c_k: Chebyshev fits over η ∈ [−0.6215, 0.7834],
  degrees 16–20 (pad to uniform table [11][deg+1]); truncation+fit
  2^−57.1 exact-arith worst at (a_T=20, K=11). Replay end-to-end (kernel
  arithmetic, real dd ops): ≤3 ULP direct, ≤1 ULP complement, a up to
  1e250 incl. λ=1±2⁻⁵⁰ and x==a.
- Routing differs per function where R1/R4 overlap: gamma_p uses R1 (P-direct)
  for x≤a+1 — R4's complement would destroy tiny-P relative accuracy;
  gamma_q uses R4 (Q-direct) for all {a≤3/2, x≤4}. Everything else:
  complement = 1 ⊖ direct_dd before the single rounding. By construction
  every complement is ≥ ~0.4, so ALL gates are relative — no lgamma-style
  absolute band needed.

**Prefactor** E = ln(x^a e^{−x}/Γ(a)) [AS SHIPPED, simplified 2026-07-28]:
E = a·LogDd(x) ⊖ x ⊖ LgammaPosDd(a) for ALL a in R1/R2 — LgammaPosDd
covers every positive a internally (its own 0<a<½ log shift, its own
2^200-scaled Stirling), so no separate large-a φ-form is needed; worst
extra cost ~1 ULP-class at the extreme non-saturated corner a ≈ 4e3. R3
needs no lgamma at all: the extracted c_k absorb the Stirling remainder
and 1/√(2πa) is formed directly (the design's original "a ≥ 8 Stirling
form" paragraph is superseded). e^E via ExpDdFrac; every factor folds into
mantissa space; power-of-two scale last (erfc pattern) so subnormal
results take one rounding. P-prefactor folds ⊖ LogDdAny(a) — Γ(a) →
Γ(a+1) without forming 1+a and without a 1/a that overflows for
subnormal a.

**New shared primitives** (beta will reuse all four):
- Log1pmxDd(u_dd): φ(u)=u−log1p(u). |u|≤1/16: φ = u²·T(u), T = dd Horner,
  18 double coeffs (−1)^k/(k+2) in u.hi, u² = DdMul(u,u). |u|>1/16:
  u ⊖ LogDdAny(TwoSum(1,u.hi), lo+=u.lo). Amplification 2⁻⁶⁸·2/u ≤ 2⁻⁶⁴ at
  the cut; err(aφ) ≤ 800·2⁻⁶⁴ in budget. NEVER the naive u ⊖ LogDd(fl(1+u)):
  small-u cancellation amplifies log_dd's 2⁻⁶⁸ by 2/u — fatal at large a.
- DdSqrt: s=Sqrt(hi); residual = (s·s−hi Sterbenz) + ops::SquareLow(s,·) —
  the s² residual MUST be capability-guarded, never bare MulSub (no-FMA
  hazard, AGENTS.md); lo = (al−e)/(2s).
- DdRecipDd: Newton from 1/hi, residual through ops::ProdLow (like DdRecip)
  plus the −r0·b.lo term.
- Expm1Dd (R4 only): |w|<2⁻¹⁰ dd series to k=6, else ExpDd(w) ⊖ 1.

**Temme kernel details**: η = CopySign(DdSqrt(2φ).hi, x−a); S in plain
double (11 Clenshaw passes + Horner in r=1/a — replay says that suffices);
rv = DdRecipDd(DdSqrt(2π_dd·a)); z≤6 (select on aφ vs exact 36.0):
½ErfcCoreDd(z.hi) ⊕ [sgn·S·rv ⊖ z.lo/√π] ⊗ e^{−aφ} — the z.lo term is the
first-order erfc correction and e^{−z²} comes from the dd aφ, NEVER by
re-squaring rounded z (2⁻⁴⁸ error near z=6 otherwise); z>6: reuse
erfc_tail_data's G: [G(1/z)/(2z) ⊕ sgn·S·rv] ⊗ e^{−aφ}, scale last
(z ≤ √800 ≈ 28.3 is inside G's fitted range; deep-tail gate inherits G's
~2 ULP class, same source as erfc's own tail bound).

**Masked-lane hygiene** (AGENTS.md rules apply): specials and E<−800
saturation lanes are masked AND their (a,x) scrubbed to (1,3) before
series/CF — j(j−a) overflows at a ≳ 4.5e306 otherwise; R3's z NaN-scrubbed
before ErfcCoreDd's value-derived gather (erfinv HalleyMid pattern).

**Clean-room Temme coefficients** (generator method, validated): extract
c_k(η) from the ORACLE by Vandermonde solve in 1/a at fixed η
(a_j = 512·2^j, j=0..14, dps 100; disjoint sample sets agree to 1e-19
through c_10). No recursion ported from anything. THE ORACLE TRAP: for
η<0 extract via the P-side identity R = ½erfc(−η√(a/2)) − P with P =
regularized LOWER computed directly — the Q-side is a 1-minus-tiny
cancellation needing ~aφ·log₁₀e digits (first attempt died exactly there).
Same small-side-direct rule for every reference point. Second trap:
mpmath's lower-gammainc (hyp1f1) fails to CONVERGE for large a near the
ridge — use mpmath both sides only for a ≤ 1e4 and exact-arithmetic
Temme (full-degree fits, mpf) as the oracle above; the fits themselves are
validated against true gammainc on overlapping range a ≤ ~2e6.

**Expected gates** (replay, ideal dd cores; pin to real measured, no
margin): R1 ≤3 / R2 ≤4 / R4 ~CR / R3 ≤3 ULP direct; complements ≤~6 ULP
relative. 2⁻⁶⁸ dd-core injection moved nothing. P+Q=1 within 1 ULP —
cheap smoke invariant. Data header src/gamma_data.h (Temme cheb table,
dd 1/n table, φ coeffs, constants incl. kGammaAT=20, N=64/44/36,
PHI_CUT=1/16, −800). References: gamma_p/q_reference.txt (`a x P Q`,
~20k pts: per-region grids, ridge lines λ=1±2⁻ᵏ, x==a at several binades,
boundary bit-brackets x=a+1/λ=½/λ=2/a=20/a=3/2/x=4/aφ=36, subnormal band
aφ ∈ [700,760], huge a incl. 2⁵³ neighborhood, tiny a, specials) plus
dd_special_reference.txt (renamed verbatim from gamma_util_reference.txt in
the 2026-07-29 hoist; dd triples for the Log1pmxDd micro-gate,
corvus_kernel_test_target pattern). Tests: smoke (specials, P+Q=1,
lane-mix determinism probing the freeze masks) + per-region ULP gates +
bench with scalar-walk-of-own-kernel baseline (libm has no gammainc;
label it, erfinv precedent).

#### Phase C part 2 shipped record [2026-07-28]

Measured (Kaby Lake AVX2 native + SSE4/SSSE3/SSE2 caps, max ULP identical
cell for cell on all four; two not-CR counts moved by ±1): direct side
R1 = 2, R2 = 2, R3 = 2, R4 = 1 ULP; complements ≤ 1 except gamma_q's R1
complement corner, 4 ULP at a = 1.5+ulp. Gates pinned, no margin.
ACCURACY.md has the full table; NEON/AVX-512 pending (CI / Ryzen).
Bench, loaded Kaby Lake, indicative, scalar-walk baseline: R1 213, R2 182,
R3 152, R4 199 ns/el (8–11× the walk).

What the stage gates caught — kept for the next family's process:
- Generator self-checks vs the design's own probes: probe1's CF a-grid
  skipped (1,2), hiding the true worst point (a→1.5⁺, x=a+1): N_cf 40→44.
  The R4 cap was under-sized from a casual estimate: 30→36. Both caught by
  the sup-proof self-checks before anything was emitted.
- Kernel review vs the reference set: the φ-series coefficient rounding
  (1/3 to double = 2^-55.9 abs) is invisible in φ and worth 12 ULP through
  e^{−aφ} at a·φ ≈ 740 (a = 3.8e5, λ = 1.062) — a band the reference set
  under-sampled and now covers; fixed by carrying the six leading φ
  coefficients as dd (kGammaPhiCoefLo, generator-emitted + self-checked).
- Kernel-stage deviations, all sound and documented at their sites: R1
  folds ⊖ LogDdAny(a) instead of ×(1/a) (subnormal-a safe); R4 multiplies
  through by a (no 1/a at all, two fewer roundings); a·log x → NaN
  overflow above a ≈ 2.5e305 is caught by the E-floor clamp; 2πa clamped
  at 2^1000 before DdSqrt (x == a at a ~ 1e308 is a legitimate unsaturated
  point); freeze masks select the accumulator rather than adding zero
  (DdAddD(s, 0) renormalizes — lane-mix determinism test is what polices
  this).
- Oracle traps now institutional: mpmath's regularized LOWER gammainc
  hangs at a ~ 1e250 and raises NoConvergence for a ≳ 1e7 near the ridge —
  reference generation switches to the exact-asymptotic oracle for every
  a > 1e4 (truncation < 1e-40 there, provable from |c_11|).
- Sub-agent workflow lesson (harness, not math): a STOPPED sub-agent is
  not woken by its own background shells or monitors — twice the tooling
  agent parked itself waiting on a watcher that fired into a stopped
  transcript. Briefs for long-running sub-agent work must say: wait
  synchronously or finish-and-report in the same turn.

### Regularized incomplete beta — detail design [2026-07-30, frontier]
Supersedes the earlier thin broad-parameters section; its decisions
(symmetry routing, fixed-depth CF core, small series, Temme regimes,
analytic Stirling-difference prefactor) all survive and are made precise
here. API settled with user: **beta_p / beta_q (a, b, x, out)**, four
equal-length spans. Everything below is [DERIVED] unless tagged.

#### Complement doctrine and routing
P = I_x(a,b), Q = I_{1−x}(b,a), P + Q = 1. One internal evaluation
computes the DIRECT side to dd and the other as 1 ⊖ dd before the single
rounding (GammaVec / GammaScale pattern, mantissa+exponent so subnormal
directs take one rounding). Complement budget: with the direct side at dd
relative ~2⁻⁶⁶, the complement's relative error is ~2⁻⁶⁶·P̃/(1−P̃), so
routing only has to keep the direct side ≤ 1 − 2⁻¹² (complement then
≤ ~¼-ulp class) — NOT ≤ ½. That slack is what makes a simple predicate
provable.

Baseline predicate: sign of λ = a − c·x with c = a + b, i.e. x vs the
mean p = a/c — computed as λ_dd = a ⊖ x ⊗ c_dd with c_dd = TwoSum(a,b);
no rounded p is ever formed. For min(a,b) ≥ 1 [claim]: I_mean ∈
[~1/e, ~1−1/e] (limits I_x(a,1) = x^a and I_x(1,b) = 1−(1−x)^b both give
exactly 1/e-class extremes) — ample slack. The predicate FAILS as
min → 0 (Q at the mean ~ α·ln(β/α) → 0, direct side → 1): below ε_R4 the
R4 logic routes instead, with the flip to the ξ^α side once
α·|ln ξ| ≳ ln 2. Generator self-check (e) proves the whole tree:
max direct-side value over a dense (α,β,ξ) boundary lattice incl.
endpoints ≤ 1 − 2⁻¹² (the gamma probe1 lesson: never skip edge grid
values). [G1a CORRECTION applies: the mean predicate alone strands
β ≪ α, ξ → 1 points — orientation is region-driven and decoupled from
the direct/complement handout; see the G1a subsection below.]

#### Region map (canonical direct side (α, β, ξ) after routing;
numbering mirrors gamma)
- **R0 specials**: table below. Masked lanes scrubbed to (2, 3, ¼)
  before ANY core — CF products, series weights, and R3's value-derived
  ErfcCoreDd gather index (AGENTS masked-lane rules).
- **R1 power series** (BPSER analog):
  I = [ξ^α / B(α,β)] · Σ_{n≥0} (1−β)_n ξⁿ / (n! (α+n)).
  Membership β·ξ ≤ B₁ and ξ ≤ ξ₁ [targets B₁ = 8, ξ₁ = 0.45 —
  ILLUSTRATIVE until sup-probe]. Fixed depth N₁; terms rise until
  n ~ βξ ≤ B₁ then decay — probe (a) proves the sup on the full
  boundary. Freeze below kBetaFreezeEps by SELECT, never add-zero
  (DdAddD(s,0) renormalizes; lane-mix test polices). 1/(α+n) via
  DdRecipDd of the EXACT TwoSum(α,n) (gamma R4 precedent). Prefactor
  folds ⊖ LogDdAny(α), never 1/α (subnormal-α safe); the series is then
  1 + α·Σ_{n≥1}, benign as α → 0.
- **R2 continued fraction** (BFRAC analog, DLMF 8.17.22):
  d_{2m} = m(β−m)ξ / ((α+2m−1)(α+2m)),
  d_{2m+1} = −(α+m)(c+m)ξ / ((α+2m)(α+2m+1)); fixed depth N₂ backward
  (GammaCfRecip pattern). Order-of-operations rule: each d formed as
  (ratio ≤ ~1) ⊗ (bounded factor), never the raw (α+m)(c+m) product —
  overflows at c ~ 1e308. Covers the middle band at moderate min; ALL
  off-ridge large-c (prefactor carries the smallness; CF is fast
  off-ridge); and the gamma-limit corner (α small, βξ > B₁) on the
  swapped side. Probe (b) must sweep the α→0, β→∞, βξ ∈ [B₁, ∞)
  gamma-limit line AND the near-ridge boundary at min ≈ T_ridge — the
  depth risk. Expected N₂ ~ 44–72 [OPEN — pinned by probe; if the
  near-ridge sup demands N₂ > ~80, lower T_ridge instead and let R3's
  fit hold to smaller ν].
- **R3 Temme erf-form ridge**: membership ν ≥ T_ridge and cψ below
  saturation [G1c CORRECTION: membership is the RATIO BAND
  ξ/p ∈ [½, 2] ∧ (1−ξ)/q ∈ [½, 2] (linked caps), saturation a separate
  overlay — gamma's shipped table has the same shape; see the G1c
  subsection], where ν = αβ/c (≤ min(α,β); the true expansion
  parameter) and
  cψ = α·φ(u) ⊕ β·φ(v),  φ = Log1pmxDd,  u = −λ/α,  v = +λ/β —
  ONE dd λ serves both, u,v in dd via DdMul with DdRecip. Then
  z = DdSqrt(cψ) with λ's sign: z ≤ 6 → ½ErfcCoreDd(z.hi) plus the z.lo
  first-order correction, e^{−z²} = e^{−cψ} from the DD cψ (never
  re-square rounded z); z > 6 → erfc_tail_data G, scale last. Gamma R3
  verbatim with cψ replacing aφ; z² = cψ exactly parallels z² = aφ.
  rv = DdRecipDd(DdSqrt(2π·ν_dd)), 2πν clamped (kGammaTwoPiAClamp
  analog). Correction series S = Σ_k e_k(ζ, p)/ν^k with ζ² = cψ/ν
  (λ's sign) and p = α/c. In the gamma limit p→0: ν→α, ζ→η_γ, so
  e_k(ζ, 0⁺) MUST reproduce gamma's c_k(η) — mandatory cross-check (d),
  and the reason the ansatz is trustworthy. Extraction is the gamma
  clean-room protocol plus one dimension: Vandermonde solve in 1/ν at
  fixed (ζ, p) against the oracle, then per-order tensor Chebyshev in
  (ζ, 2p−1); symmetry e_k(ζ,p) = ±e_k(−ζ,q) halves the table (check
  (h)). Δδ (Binet difference) is ABSORBED by the extracted e_k — its
  1/α^{2k−1} = (q/ν)^{2k−1}-type terms are p-smooth in the ν-power
  frame, same as gamma's c_k absorbing the Stirling remainder. Budget:
  table ≤ 32 KB [K = 11 orders, degrees ~(20, 12) ILLUSTRATIVE];
  runtime ~2–3× gamma R3's S, ridge lanes only. T_ridge target 32
  [OPEN — co-pinned with N₂].
- **R4 tiny-min** (APSER/FPSER analog): min(α,β) ≤ ε_R4 [target 2⁻⁶,
  OPEN], βξ ≤ B₁, α|ln ξ| ≲ ln 2 — the α-scaled side ~ α·J is direct.
  Assembly is gamma-R4 verbatim in beta clothing: nothing is ever
  "1 − rounded near-1" — ξ^α − 1 = Expm1Dd(α ⊗ LogDdAny(ξ)),
  Γ-ratio − 1 = Expm1Dd(exact-argument lgamma combination), products of
  (1+Eᵢ) expanded so the 1s cancel analytically. Exact-argument
  lgammas near 1 come from lgamma's ZoneBracket at the exact shifted t
  (the GammaSmallQ mechanism — already exposed by lgamma-inl.h). The
  box's complement cases route away by construction: α|ln ξ| large → R1
  (ξ^α side direct), βξ > B₁ → R2 swapped CF. Probe (f) proves the
  α-expansion truncation over the closed box. [G1a CORRECTION: the box
  is stated in the TINY-FIRST orientation and carries explicit caps
  ξτ ≤ ξ₁ and B·ξτ ≤ B₁ — see the G1a subsection.]

#### Prefactor E = α ln ξ + β ln(1−ξ) − ln B(α,β)
The central hazard. y_dd = TwoSum(1, −ξ) is EXACT and feeds LogDdAny's
dd-argument path (the Log1pmxDd log-side mechanism) — a plain fl(1−ξ)
for ξ < ½ costs β·2⁻⁵³ absolute in E, fatal at large β. ln B is NEVER
three independently rounded lgammas at large c. A new intrinsic hazard
vs gamma: c = a + b is the first ROUNDED lgamma argument in the family
(gamma's arguments were all exact inputs); fl(c) alone costs
ψ(c)·c·2⁻⁵³ absolute — e.g. 2⁻⁴²·class at c = 256. Every path below
carries c as the exact TwoSum dd pair. Three assembly paths, all ending
in the GammaClampE-analog → ExpDdFrac → factors folded in mantissa
space → ScaleTwo LAST (one rounding into subnormals):
- **P1** (c ≤ C_lg [target 256]):
  E = α·LogDdAny(ξ) ⊕ β·LogDdAny(y_dd) ⊖ LgammaPosDd(α) ⊖
  LgammaPosDd(β) ⊕ [LgammaPosDd(c.hi) ⊕ c.lo·ψ̃(c.hi)] — the c.lo
  correction uses DigammaRough (new, ~2⁻⁴⁰ suffices: residual
  c.lo·ψ·2⁻⁴⁰ ≪ budget). Cancellation budget: |lgamma| ≤ ~1400 at
  c ≤ 512-class → 3·1400·2⁻⁶⁸ ≈ 2⁻⁵⁶ absolute in E ✓.
- **P2** (c > C_lg, α ≥ Z₀ [target 10]): the analytic
  Stirling-difference — algebra collapses to
  E = −cψ ⊕ ½·LogDdAny(ν_dd/2π) ⊖ Δδ,  Δδ = BinetDd(α) ⊕ BinetDd(β)
  ⊖ BinetDd(c.hi) (Binet's derivative ≤ 1/(12z²) makes c.lo
  negligible here). cψ is the SAME dd machinery as R3 — one
  implementation. Derivation: Stirling on all three lgammas leaves
  a ln(x/p) + b ln(y/q) = au + bv − (aφ(u)+bφ(v)) and au + bv = 0
  identically (u = −λ/a, v = λ/b) — the cancellation is removed on
  paper, not in floats. Budget: cψ ≤ ~750 pre-saturation, error
  ≤ 2·750·2⁻⁶⁴ ≈ 2⁻⁵³·½ — the same accepted ½-ulp-class extreme as
  gamma's aφ budget line.
- **P3** (c > C_lg, α < Z₀): E = α·LogDdAny(ξ) ⊕ β·LogDdAny(y_dd) ⊖
  LgammaPosDd(α) ⊕ LgammaDiffDd(β, α), where LgammaDiffDd(β,α) =
  lgamma(β+α) − lgamma(β) is computed analytically (α·ln β, a
  Log1pmx-decomposed (β+α−½)·log1p(α/β) term, −α, Binet difference) —
  β ≥ c/2 > 128 ≥ Z₀ always holds here, and c is never formed as an
  lgamma argument. Cancellation scale is capped by α < Z₀:
  pieces ≤ ~10·|ln| ≤ ~7100 → error ≤ 2⁻⁵⁵ absolute ✓.
Overflow guard: c > ~2¹⁰²² lanes get α, β, λ halved (exact powers of
two; ψ is 0-homogeneous in the pair through u, v, so (c/2)ψ ⊗ 2 is
exact) [precise guard site at implementation]. Saturation: E ≤ E_floor
(kGammaExpFloor reuse, −800) → direct saturates to 0, complement to 1,
lane scrubbed.

#### Specials [pinned now; gamma-consistent doctrine: one degenerate
parameter gets its limit, two degeneracies (or a degenerate parameter
meeting the x-boundary its mass sits on) → NaN]
- NaN anywhere → NaN (payload preserved). x ∉ [0,1] → NaN. a < 0 or
  b < 0 → NaN. Two of {a ∈ {0,∞}, b ∈ {0,∞}} degenerate → NaN.
- x = 0 → P = +0, Q = 1;  x = 1 → P = 1, Q = +0  (a, b > 0 finite).
- a = 0 (mass at 0): P = 1 for x ∈ (0,1]; NaN at x = 0.
- b = 0 (mass at 1): P = +0 for x ∈ [0,1); NaN at x = 1.
- b = +inf (mass at 0): P = 1 for x ∈ (0,1]; NaN at x = 0.
- a = +inf (mass at 1): P = +0 for x ∈ [0,1); NaN at x = 1
  (maps to gamma's P(∞,∞) = NaN under x=1 ↔ x=∞).
Exact/brutal invariants for smoke + reference: I_½(a,a) = ½ EXACTLY;
P + Q = 1 within 1 ulp; analytic lines I_x(a,1) = x^a,
I_x(1,b) = 1−(1−x)^b, I_x(½,½) = (2/π)·asin(√x); lane-mix determinism.

#### New shared primitives (and homes)
- **BinetDd(z)**, z ≥ Z₀: Stirling tail Σ B_{2k}/(2k(2k−1)z^{2k−1}).
  Expose from lgamma-inl.h if its Stirling zone factors cleanly, else
  fresh with generator-emitted coefficients [OPEN at implementation —
  ANY touch of lgamma-inl.h re-runs the full lgamma+gamma byte-identity
  guard, the dd_special-hoist protocol].
- **LgammaDiffDd(β, α)** (P3 above): lives in beta-inl.h until a second
  consumer (inverse incomplete beta) hoists it — the Log1pmxDd history
  pattern, noted at the definition site.
- **DigammaRough(z)** on (0, 2Z₀]: small poly + recurrence, ~2⁻⁴⁰; also
  the natural seed for the future public digamma (P1 roadmap).
- **cψ/λ/ν machinery** shared by P2 and R3 — single implementation.

#### Oracle and reference set
gen_beta_reference.py: mpmath betainc regularized, SMALL SIDE DIRECT
always (swap the argument triple; 1−x is exact in mpf) — the gamma
oracle rule, enforced at the same kind of site comment. Determine
mpmath's convergence/latency ceiling in (a,b) near the ridge
empirically [OPEN, G1]; above it, the exact-asymptotic oracle is our
own erf-form at full degree in mpf, validated on an overlap band
(gamma protocol verbatim). Secondary disjoint-sample cross-check via
high-dps power series. NEW generator with its OWN seed — shares nothing
with gamma's rng stream (that stream stays frozen; AGENTS note).
Point economics [target ~40k pts, ~4 MB total]: per-region
log-lattices in (α,β) × ξ grids; ridge lines λ = 0± at binades of ν
crossed with p ∈ {tiny, ¼, ½, ¾, 1−tiny}; bit-brackets BOTH sides of
EVERY boundary (routing predicate, B₁, ξ₁, T_ridge, ε_R4, C_lg, Z₀,
z = 6, E_floor); the x = ½, a = b diagonal; the three analytic lines;
subnormal band E ∈ [−745, −700]; huge/tiny parameters incl. subnormal
α and c near overflow; the full specials table.

#### Data, generator self-checks (exit nonzero, budgets to stderr)
src/beta_data.h ← tools/gen_beta_data.py: routing/region constants,
N₁/N₂, R3 tensor tables + K, Binet coefficients + Z₀, DigammaRough
poly, E_floor, clamps. Self-checks: (a) R1 truncation sup on the full
membership boundary incl. α→0 and the (B₁, ξ₁) edges; (b) N₂ sup incl.
gamma-limit line and near-ridge min ≈ T_ridge; (c) R3 fit residual +
total S truncation over a (ζ, p, ν) lattice vs the mpf oracle;
(d) e_k(ζ, p→0⁺) == gamma c_k(η); (e) routing-safety max direct side
≤ 1 − 2⁻¹²; (f) R4 α-expansion truncation over its closed box;
(g) Binet truncation at Z₀ vs dd target; (h) e_k symmetry identity.

#### Kernel, TU, build
src/beta-inl.h + src/beta.cpp, OWN TU — dependency set differs
materially from gamma's (consumes dd, dd_special, exp_dd, log_dd,
lgamma-inl (LgammaPosDd, ZoneBracket, Binet exposure), erfc_core, ops;
does NOT include gamma-inl.h — the R2 sweep covers the gamma-limit
corner, so no gamma-core reuse; decided here). HWY_NOINLINE from day
one on every region core, every prefactor path, AND the per-lane
driver (it inlines twice per export). /d2ReducedOptimizeHugeFunctions
for beta.cpp only if MSVC CI times degrade (gamma's threshold story).
HWY_DYNAMIC_DISPATCH called from inside namespace corvus (SSE2-cap
collapse rule).

#### Gates and tests [ILLUSTRATIVE targets; pin to measured, no margin]
Direct side ≤ 3–4 ULP per region (R4 ~1–2); complements relative
(≥ 2⁻¹²-ish by routing) ≤ ~6 ULP; P + Q = 1 ≤ 1 ulp; I_½(a,a) exact.
test_beta_smoke (specials table, invariants, lane-mix) +
test_beta_ulp (public API, `a b x P Q` reference files) — registered
after test_gamma_ulp in ALL FOUR lists (tests/CMakeLists.txt, three
ci.yml report steps, sweep_tiers.ps1 $gates). bench_beta with the
scalar-walk-of-own-kernel baseline, labeled (no libm betainc; erfinv
precedent). ACCURACY.md + README in the SAME change set as the kernel.

#### Process: sub-agent split and review gates
Sonnet tooling + Opus kernel, briefs composed from this section at
spawn time, finish-and-report (never park on a watcher — the stopped-
subagent lesson).
- **G1** (Sonnet, generator): gen_beta_data.py per this design; the
  [OPEN] constants (N₁, N₂, T_ridge, B₁, ξ₁, ε_R4, C_lg, Z₀, K, fit
  degrees, mpmath ceiling) are pinned by its self-checks. Orchestrator
  reviews the stderr budget lines against this section BEFORE any
  table is committed. Escalate to frontier if: any self-check cannot
  meet budget at target constants; cross-check (d) fails; table
  exceeds 32 KB.
- **G2** (Sonnet, references): gen_beta_reference.py + reference set;
  orchestrator spot-audits vs independent oracle samples incl. all
  analytic lines.
- **G3** (Opus, kernel): reviewed against this design BEFORE the first
  ULP run — checklist: freeze-by-select, masked-lane scrubs (CF
  products, R3 gather), d-term order of operations, ⊖LogDdAny(α)
  folds, exact c_dd handling in every path, scale-last.
- **G4**: measured ULP tables → gates pinned. Tier order: Kaby Lake
  native + caps first, NEON via CI, Ryzen native AVX3 last (stability
  watch).
- **G5**: four-list registration audit, CI green, ACCURACY.md/README/
  PLAN in the change set.

#### G1a probe results and design corrections [2026-07-30]
Sonnet probe agent, orchestrator-reviewed at frontier effort. Pinned
[DERIVED, empirical]: **N₁ = 64** (B₁ = 8, ξ₁ = 0.45 unchanged; worst
case the β→0 geometric corner, 2⁻⁷⁴·⁹ margin); **N₂ = 64, T_ridge =
32** (binding case the symmetric middle (32,32,ξ≈0.70), sharp elbow to
2⁻¹⁸⁷ by N=64; gamma-limit line needs only 48); **Z₀ = 10, K_B = 16**;
**ε_R4 = 2⁻⁶ confirmed** (routing predicate safe above it — max direct
0.9120 — and fails ONLY below it, checked to min ~1e-42 incl. the
(ε_R4, 1) band); **C_lg = 256 PROVISIONAL** (budget passes by 0.23
bits; drop to 128 if G4 measures hot); **mpmath betainc ceiling**:
reliable but latency-bound — safe to min(a,b) ~2000 near the ridge,
seconds by ~4000, hangs beyond; failure mode is latency ONLY, never
wrong digits. Oracle hygiene lesson from the probe harness itself: a
stale mp.dps (set outside the worker loop) produced fake convergence
plateaus — generators must set dps INSIDE every computation context.

Two escalations, resolved here:
- **R4/routing coverage gap (real design flaw).** As first written, R4's
  box had no ξ cap and the mean predicate alone strands points: the
  probe's witness (8, 2⁻⁶, 1−9.5e-7) — a legitimate direct value 0.16 —
  was outside every region, and the β ≪ α, ξ → 1 family generalizes it.
  Root cause: the design conflated ORIENTATION (which argument triple
  the kernel evaluates) with the direct/complement HANDOUT. Correction:
  they are decoupled. The ≤ 1 − 2⁻¹² doctrine applies to the EVALUATED
  side; the handout is a final select. Orientation is REGION-DRIVEN,
  first match:
  1. R1 if either orientation satisfies ξ ≤ ξ₁ ∧ βξ ≤ B₁ (at most one
     orientation has ξ ≤ 0.45; R1 handles tiny first-param via its
     1 + α·Σ form down to α ~ 2⁻¹¹, complement-slack covers the rest).
  2. R4 if min ≤ ε_R4, evaluated TINY-FIRST (τ, B, ξτ) with
     **ξτ ≤ ξ₁ AND B·ξτ ≤ B₁** (the missing caps; the probe confirms
     N = 48–56 inside them). ξτ = 1−x is exact for x ≥ ½ (Sterbenz),
     and for x < ½ the 2⁻⁵⁴ abs error is damped by ∂I/∂ξ ~ τ.
  3. R3 if ν ≥ T_ridge ∧ cψ ≤ saturation (mean-predicate orientation;
     ridge symmetric under swap).
  4. Else R2 CF, orientation = the side where the CF converges at
     N₂ — candidate crisp rule ξ < (α+1)/(c+2) (the DLMF-fast side,
     reachable in exactly one orientation up to an overlap band), to be
     PINNED by a G1b probe measuring the depth surface in BOTH
     orientations over the gap zones (second-param ≤ 1 with
     ξ ∈ (0.45, mean]; tiny-first with Bξτ ∈ (B₁, 800]).
  Saturation (E ≤ E_floor) remains the net under every region in every
  orientation. Self-check (e) must be RE-PROVEN under this final
  orientation rule (evaluated side ≤ 1 − 2⁻¹²), not just the mean
  predicate.
- **R3 extraction ill-conditioning (protocol fix, ansatz intact).** The
  spot-probe stabilized only k ≤ 1: the mpmath ceiling capped the
  1/ν ladder at ν ≤ 2048 (7 rungs), far short of gamma's 15-rung span —
  a conditioning problem, not evidence against the (ζ, p, 1/ν) form.
  The e₀ sign flip vs gamma's c₀ is a convention artifact suspect (the
  design's ζ carries λ = a − cx's sign; gamma's η carries x − a's) and
  the ~7% magnitude gap smells of rv normalization — both to be nailed,
  not fitted around. G1b protocol, in order: (1) build a tanh–sinh
  quad oracle for ridge points (split at the mode, exact log-prefactor
  at high dps; validate vs betainc below the ceiling), sidestepping the
  hyp2f1 latency wall and restoring the tall ladder; (2) anchor at the
  GAMMA LIMIT first — p = 2⁻⁴⁰, extract e_k along ν = 512·2^j,
  j = 0..12+, and match gamma's c_k (conventions read from
  gen_gamma_data.py's own extraction) to ≤ 1e-15 through k = 5 BEFORE
  any 2D work — converting the inconclusive cross-check (d) into the
  anchor; (3) only then the 2D spots, solved as regularized least
  squares in a Chebyshev-in-1/ν basis on disjoint oversampled ladders.

#### G1b probe results and second routing correction [2026-07-31]
- **R2 orientation rule PINNED**: evaluate on the side with
  ξ < (α+1)/(c+2), else swap — confirmed as literally stated over all
  256 gap-zone points, zero failures, N₂ = 64 reconfirmed (worst
  predicted-side error 2.8e-26 ≪ 2⁻⁶⁸). One lattice extension owed:
  the second-param band was probed only down to β = 2⁻⁶ — generator
  self-check (b) must extend it to β → 1e-300.
- **R3 extraction VIABLE with wide margin.** Quad oracle validated
  (63/63, ≤ 9.2e-47; 0.04–0.09 s/pt to min ~1e8; needs the u = t^a
  substitution for a < 1 — endpoint singularity gave WRONG answers
  before it); the gamma-limit anchor resolved both G1a artifacts: the
  sign flip was a branch-order bug, and the magnitude gap is the exact
  variable mapping **η_γ = −ζ·√2 as p → 0** — which is what the
  definitions predict (ζ² = cψ/ν → φ = η²/2 at the limit; λ = a − cx
  opposes x − a), so cross-check (d) is now stated as this mapping.
  Anchor match ≤ 1e-15 through k = 5 at p = 2⁻⁵⁰ (residual at 2⁻⁴⁰ is
  O(p), slope-confirmed). 2D spots: all 18 (p, ζ) combos stable
  through k = 8 (disjoint-set agreement 2.7e-30..7.8e-27 vs a 1e-13
  bar); extraction = regularized least squares in Chebyshev-in-1/ν,
  extract K = 15 keep ~9 (truncation-bias lesson). The anchor-ladder
  oracle is the CF (quad needs infeasible dps for that shape).
- **Task C ESCALATE, resolved here (second routing flaw, rule-ORDER).**
  Witness (1e-20, 1, 0.4): R1-native fires (its box has no α floor)
  and evaluates 1.0 — yet G1a's own analysis had derived R1's tiny-α
  validity floor (complement-slack needs α·J ≥ 2⁻¹²); the floor never
  made it into the decision list. Correction: hoist the tiny-min guard
  ABOVE R1 — final orientation order:
  0. min(a,b) ≤ ε_R4 → tiny-first (τ, B, ξτ): if τ|ln ξτ| ≤ ln 2 ∧
     ξτ ≤ ξ₁ ∧ B·ξτ ≤ B₁ → **R4**; else fall through (its ξτ^τ-side
     and CF cases land correctly in R1/R2 below — spot-verified at the
     witness's neighbors).
  1. R1 if either orientation has ξ ≤ ξ₁ ∧ βξ ≤ B₁.
  2. R3 if ν ≥ T_ridge ∧ cψ ≤ saturation.
  3. R2, orientation by the pinned ξ < (α+1)/(c+2) rule.
  Self-check (e) runs against THIS order in the generator and must
  pass there — it is the regression guard for both routing flaws.
- Probe-harness hygiene, now thrice-learned: mp.dps must be set inside
  EVERY computation layer (worker, subprocess-send, subprocess-parse —
  G1b found the same trap at three layers). And both G1a and G1b agents
  parked on background jobs despite briefs forbidding it — future
  briefs say: foreground only, sweeps chunked into ≤ ~5-minute
  re-runnable commands, no Monitors.

#### G1c generator results and third correction [2026-07-31]
First cut of tools/gen_beta_data.py landed (uncommitted): checks (a),
(b), (d), (f), (g), (h), (i) all PASS with margin (R1 2⁻⁷⁴·⁹, CF
2⁻⁷⁶·⁷, anchor 6.6e-16, Binet 2⁻⁷⁴·⁶, DigammaRough 2⁻⁵⁰),
bit-identical across runs, ~310 s. Two escalations, both resolved:
- **R3 fit stalled at 2⁻¹⁶ → root cause was a DESIGN error in R3's
  membership, found by reading gen_gamma_data.py.** Gamma's shipped
  Temme table spans only η ∈ [−√(2φ(½)), +√(2φ(2))] ≈ [−0.62, +0.78] —
  the ridge RATIO band λ ∈ [½, 2] — with NNODES = 33 and a 2⁻⁵⁶ replay
  target; outside the band, series/CF cover even at huge a. This
  design's "R3 = ν ≥ T ∧ cψ ≤ 800" mis-remembered that, and the ζ ∈
  [−5, 5] fit domain it implied is what no 32 KB fit can span at dd
  level (a 2⁻¹⁶ residual would cost ~2⁻¹⁹ relative on ridge values —
  unshippable). CORRECTED R3 membership, mirroring gamma exactly:
  **ν ≥ T_ridge ∧ ξ/p ∈ [½, 2] ∧ (1−ξ)/q ∈ [½, 2]** (the caps are
  linked by pu + qv = 0; at p = ½ the joint cap binds at u = ±½),
  giving a derived ζ band ⊂ ~[−0.76, +0.76] (generator derives and
  states the exact sup). Saturation stays a separate overlay; the
  z ∈ (6, √800] G-tail is unchanged (live wide-z lanes only occur at
  moderate ν where cψ ≤ 800). Check (c) target becomes gamma-class
  2⁻⁵⁶, not pinned-to-measured. NEW risk transferred to R2, to be
  probed in the revision's check (b): the CF now owns the band edge at
  ALL ν — sweep u = ±caps at ν ∈ {32, 1e4, 1e8, 1e12, 1e16} to confirm
  N₂ = 64 is ν-independent at fixed ratio distance; growth with ν is
  an ESCALATE.
- **Check (e) "failures" in R4 are a category error, accepted as the
  agent gated it.** R4's contract is analytic small-side assembly (the
  Expm1Dd product expansion) — it never complements a rounded near-1
  evaluation, so the ≤ 1 − 2⁻¹² doctrine is satisfied by construction
  there. (e) formally gates R1/R2/R3; R4's guarantee is a G3 kernel-
  review checklist item (confirm the assembly produces the small
  quantity directly) plus the smoke invariants.

#### G1c revision shipped — G1 COMPLETE [2026-07-31]
tools/gen_beta_data.py + src/beta_data.h landed (orchestrator-reviewed,
independently re-run and hash-verified). All self-checks pass at their
fixed targets: (a) 2⁻⁷⁴·⁹, (b) 2⁻⁷⁶·⁷ plus the transferred-risk
ratio-cap × ν sweep — CF margin IMPROVES with ν (2⁻¹⁵⁹ at ν = 32,
2⁻³¹⁴ at 1e4, exact beyond), so the band handoff is safe at all ν;
(c) 2⁻⁵⁷·³⁵ vs the fixed 2⁻⁵⁶ target over the corrected band;
(d) anchor 6.5e-16; (e) R1/R2/R3 max 0.9973, zero failures, both
escalation witnesses in-lattice; (f) R4 N = 48 → 2⁻⁶⁰·⁵; (g) 2⁻⁷⁴·⁶;
(h) 2.6e-43 (sign resolved: e_k(ζ,p) = −e_k(−ζ,1−p)); (i) 2⁻⁵⁰.
R3 table: K = 10 rows × 25 × 15 tensor-Chebyshev, 29.3 KiB. Exact band
sup DERIVED: ζ_max = √(3·ln 2/2) ≈ 1.0197 — both ratio caps bind
jointly at p = ⅓ (closed form verified; the G1c subsection's ~0.76 was
the p = ½ slice, superseded). Two findings worth keeping: the binding
fit axis was N_p (p-direction), not N_ζ; and R3 test points must be
scaled to the p-local reachable lens (zeta_max_at_p), not the global
sup — sampling outside the lens fakes a fit plateau. Oracle note: the
CF is the R3 extraction oracle exclusively (quad needs infeasible dps
for the scaled-remainder ladder shape; deviation documented in the
generator). Generator runs ~7 min, bit-identical across runs and
machines-of-invocation; mpmath sits in system python (the throwaway
venv creation half-failed — rebuild the venv per AGENTS.md before any
regeneration on a clean machine).

#### G2 shipped — reference set [2026-08-01]
tools/gen_beta_reference.py → tests/data/beta_p/q_reference.txt
(35,478 rows each, 3.7 MB, byte-identical pair per the gamma
convention; `a b x P Q`, fresh seed, nothing shared with gamma's rng
stream). Oracle: CF (G1-validated) primary, mpmath betainc
cross-check (worst 7.5e-24), plus an APSER-style small-τ branch
Q̃ = −expm1(w + ln S) for min(a,b) ≤ 2⁻⁴ — added after the first cut
dropped 1,614 points whose single root cause was CF failure in the
small-τ/small-ξ corner (drop histogram audit; the lost sets included
all nine ε_R4/ln2-wall brackets). The branch is validated at 2⁻¹⁸⁸
worst diff on the CF overlap band and doubles as independent
validation of R4's own kernel expansion. Final drops: 2 (one huge
symmetric ridge bracket a=b=1e10, x=½±ulp, CF non-convergent at every
dps — the exact-½ point itself survives). Region histogram
R1/R2/R3/R4 = 13,484/5,937/2,906/12,900, all ≥ 2000; P+Q=1 bit-exact;
small-side-direct enforced; analytic lines at 1.5e-31. Four generator
bugs found by its own checks, incl. the checker itself computing a
forbidden 1−(near-1) — fixed with expm1/acos identities.

#### G3 kernel results, fourth routing correction, escalation
resolutions [2026-08-01]
Opus kernel agent delivered src/beta-inl.h + beta.cpp + both tests +
four-list registration (on branch beta/g3-kernel, NOT main — three
routing/oracle issues below must close before merge). clang-cl
AVX3_ZEN4: beta.cpp compiles in 19.5 s (no codegen blowup), ctest
15/15 green (existing gates untouched — lgamma-inl.h NOT modified;
BinetDd implemented fresh from kBetaBinetCoef), smoke fully green
incl. the whole specials table and I_½(a,a) exact. Measured
PROVISIONAL direct sides: R1 = 1, R2 = 0, R3 = 3, R4 = 2–3 ULP — all
inside targets. Six kernel bugs found by the reference sweep, each
with a witness (R3 orientation×λ sign XNOR; freeze eps 2⁻⁶⁰→2⁻¹⁰⁵;
R2 d-terms in dd; R2 numerator overflow — two-bounded-factors rule
restored; LgammaDiffDd compensation term removed in favor of the
φ-form; and the FOURTH ROUTING CORRECTION, ratified: R4's ξτ cap is
max(ξ₁, thr_τ) — for B < 0.24 there is a ξτ ∈ (ξ₁, thr_τ) window
where neither R1 orientation fires and R2 would evaluate the near-one
side; measured N=48 truncation in the widened window is 2⁻⁵⁷·³ (vs
2⁻⁷¹·⁹ in-box) — G4 decides from measurement whether R4 depth bumps
to 56). Design deviations ratified: P1 (three-lgamma) path replaced
by P3's LgammaDiffDd(max,min) ⊖ LgammaPosDd(min) everywhere — never
forms rounded c as an lgamma argument, and kBetaDigammaCoef is
therefore UNUSED (keep in header for the future public digamma; G5
cleanup decision); P2/P3 split on min(α,β), not α (α=255, β=2 broke
the literal predicate); c-overflow guard is an exact 2⁻²⁰⁰ prescale
inside the cψ core; λ = α·y ⊖ β·ξ with TwoSum (algebraically a − cx,
symmetric); R2 folds ⊖LogDdAny(α) (the CF prefactor's own form).

Escalations, resolved at frontier effort:
- **(A) G2 reference defect, kernel adjudicated RIGHT (237 vs 3).**
  small_tau_oracle forms lnΓ(B+τ) − lnΓ(B) at fixed dps: for τ ≪ B
  the difference is identically 0 in mpf and the emitted value goes
  INDEPENDENT of the max parameter (witness family (a, 1.4e-300,
  1−2⁻⁵²): constant P across a ∈ [2⁻²⁰, 2⁻⁶] where ψ(a) ~ −1/a should
  swing it by 2¹⁴). FIX (G2 revision): τ ≤ B·1e-8 → analytic Taylor
  τ·ψ(B) + τ²/2·ψ₁(B) (with a τ³ bound check); else direct difference
  at dps ≥ base + log10(B/τ) + 20. ~2394 rows regenerate; the
  kernel's kRefDefectCutoff hatch is DELETED after regeneration.
- **(B) FIFTH routing correction: R1 needs the λ ≥ 0 side.** R1's box
  admits an evaluated side of 1 − 1.2e-6 (witness (0.158, 20, 0.396);
  429 ULP complement — E1 error already at the 2⁻⁷⁰ floor, so this is
  routing, not arithmetic). Correction: R1 fires only in the
  orientation whose exact λ_dd ≥ 0 (evaluated ξ ≤ its mean; exactly
  one orientation qualifies). Displaced traffic (ξ ∈ (mean, ξ₁],
  βξ ≤ B₁) lands in R2 via the orientation rule — spot-verified sound
  incl. both witnesses; generator check (b) gains the moved-traffic
  lattice and check (e) gains the pocket its lattice MISSED (dense
  sampling near the βξ = B₁ edge with α ∈ (ε_R4, 1) — 0.9973 was a
  lattice artifact, not the true sup).
- **(C) Gamma-limit slice: the no-gamma-core decision is PARTIALLY
  REVERSED.** At (0.05, 1e100, 2e-99) the beta CF is structurally
  degenerate (d₁ → −(1 − 2e-99); mpmath's own CF divides by zero at
  dps 40) — no depth or precision rescues it. For max-param ≥
  kBetaGammaLim [target 2⁸⁰, ILLUSTRATIVE — generator pins by
  BOTH-SIDES overlap: beta-CF validity probed upward from the 2⁴⁰
  G1a-validated line, gamma-form correction error ≤ 2⁻⁶⁰ probed
  downward] R2 routes to a gamma-limit path: t_dd = −β·log1p(−ξ) in
  dd, E = α·LogDdAny(t) ⊖ t ⊖ LgammaPosDd(α), core via gamma-inl.h's
  own GammaCfRecip / GammaSeriesSum templates (instantiating exactly
  two gamma cores in beta.o — templates cost only what is called; the
  TU-boundary objection was about wholesale instantiation). The G2
  oracle for that corner likewise switches to mpmath gammainc.
- **(D)** covered by the ratified deviations above.

Revision sequence: one Sonnet tooling agent updates gen_beta_data.py
(route_final: λ≥0 R1 rule, R4 window cap, gamma-limit slice + B_gl
pin; checks (b)/(e)/(f) lattice extensions; header regen) and
gen_beta_reference.py (small-τ Taylor/adaptive-dps fix, gamma-corner
gammainc oracle, full regeneration + drop audit); then the Opus
kernel agent applies the two routing changes + deletes the defect
hatch; then G4.

#### G1/G2 revision results — SIXTH routing correction, B_GL pinned
[2026-08-02]
The fifth correction (R1 requires λ ≥ 0) is REVERTED as too blunt:
its own displaced traffic broke check (b) — witness (0.158, 1000,
0.00251) needs CF depth ~512, yet R1-native served it perfectly
(evaluated 0.994, complement at 2⁻⁵⁸). **SIXTH correction (the
precise form): R1 fires in either orientation as before; lanes whose
R1-evaluated dd value exceeds 1 − 2⁻¹¹ are POST-ROUTED into the R2-CF
pass in the opposite orientation** (a mask update between passes the
kernel already runs; the small side is then computed directly in the
G1b-validated CF band). Both G3 witnesses re-route correctly; the
complement-slack doctrine is enforced by construction. Generator:
check (b) drops the moved-traffic lattice (vi), gains the re-route
band (viii) (swapped triples (β, α, 1−ξ) with α ∈ (ε_R4, ~1.2],
near-one condition α|ln ξ| ≲ 1.5, βξ ≤ B₁); check (e) re-proved with
the post-route applied.
**kBetaGammaLim = 2⁵⁹ RATIFIED, provisional-to-measurement.** The
literal 2⁻⁶⁰ overlap bar is unachievable (CF-dd conditioning ceiling
2⁵⁰·³ vs gamma-form floor 2⁶⁸·⁶ — an 18-binade inversion); the agent's
2⁻⁴⁹ bar (= the design's own ~6-ULP complement gate class) gives a
3.75-binade overlap, pin 2⁵⁹. G4 reports R2-gammalim as its own
region row; measured > 6-ULP class ⇒ follow-up item: first-order
Temme correction term on the gamma form (not a threshold change).
Reference set: oracle fix (A) verified (defect family now varies with
the max parameter, matches independent Taylor derivation); blast
radius 10,268 changed rows (the fix touches ALL min ≤ 2⁻⁴ points, not
just the witness family — expected on reflection). OPEN in the next
cycle: 5,978 drops (the new cancellation guard correctly rejects what
the old oracle silently corrupted, but the (τ ≤ 2⁻⁴, B ~ 10³) gap
needs a betainc-with-timeout rescue — reliable exactly there,
off-ridge below the latency ceiling); reference VALUES are
routing-independent so the current files are a valid checkpoint, but
the histogram/coverage audit re-runs under the sixth correction.
Windows multiprocessing lesson (generator): parent_process() is
unreliable under spawn on this box — use current_process().name.

#### SEVENTH correction + slice IMPLEMENTED (frontier, hands-on)
[2026-08-03 — the orchestrator took over implementation per user
decision; the escalation density showed every stage needs
design-boundary judgment]
The sixth correction's CF destination FAILED its own check (b)(viii)
at 2⁻⁵⁵·⁵ (witness (0.0234, 1e6, 4e-6)) — the opposite-orientation CF
inherits the small-second-parameter weakness. **SEVENTH correction
(supersedes the sixth's destination): near-one R1 lanes fold into
R4's analytic assembly in the FIRED orientation** — they are R4-shaped
by construction (R1's box supplies R4's caps; near-one puts the Expm1
argument below ~2⁻¹⁰; ε_R4 was a routing threshold, never a validity
bound; post-routed lanes provably have α = min and τ ≤ ~1.35, gated
at kBetaPrTauMax = 1.5 = lgamma's centre-2 edge). Measured: post-route
domain truncation 2⁻⁶⁷·⁵ at N = 48 (vs the CF's 2⁻⁵⁵·⁵ — structurally
better), check (e) PASSES (worst evaluated 0.99948), both witnesses
post-route. IMPLEMENTED in: gen_beta_data.py (route_final, checks
(b)/(e)/(f), kBetaNearOne/kBetaPrTauMax/kBetaGlRidgeMin emission),
src/beta-inl.h (BetaR4Tiny two-zone lgamma; driver post-route mask
into m_r4x; is_p flip), tests/test_beta_ulp.cpp (router sync incl.
R4-postroute + R2-gammalim buckets; kRefDefectCutoff hatch DELETED).
The (C) slice is ALSO implemented: R3 ridge floor drops to
kBetaGlRidgeMin = 20 for max ≥ B_GL in-band (check (c) extension
lattice proves the 1/ν extrapolation at the anchor p; gamma's own
table has the same shape); off-band slice lanes go through the gamma
limit via gamma-inl.h's GammaSeriesSum/GammaCfRecip templates with
beta-side dd prefactor (E_g = s·ln t ⊖ t ⊖ lgamma(s), t dd; val =
naturally-computed side, is_p double-XNOR — never a dd complement
round-trip of a small side).

#### R3 depth extension + EIGHTH correction (frontier, hands-on)
[2026-08-03, post-compaction session — both found by the verification
run the seventh-correction commit left in flight]
1. **Check (c) extension lattice FAILED at 2⁻⁵⁰·²⁷ (ν = 20, p-edge)**
   — the K = 10 table's 1/ν truncation alone (probe-matched exactly),
   NOT a fit or extrapolation defect. Fix: **p→0-edge DEPTH EXTENSION**
   kBetaR3GlExt[3][25] — e_k(ζ, 2⁻²⁰) for k = 10..12, 1D Chebyshev in
   the same t = ζ/ζ_max, applied by BetaR4Temme→BetaR3Temme on
   gamma-limit lanes only (masked add ·r¹⁰; safe at every slice ν
   since non-negligible ⇒ p tiny). K = 13 measured 2⁻⁶⁰·¹ worst
   all-interpolated (4-bit margin); table now 29.88 KiB of 32.
   TRAPS learned: (a) extraction at p = 2⁻⁵⁰ is INVALID — the ladder's
   c = ν/p crosses the CF-ground-truth ceiling (~2⁶¹, B_GL derivation's
   own finding) and LSQ fits the noise into high orders; extract at
   2⁻²⁰ (ladder c ≤ 2³⁵; e_k p-slope ~1e-7, ÷ν¹⁰ ⇒ 2⁻⁶⁶-class).
   (b) POINTWISE extraction at small nonzero ζ (|ζ| ~ 0.03) is
   ill-conditioned at ANY p (e_9: 5.9e-4 → 19; ~1/ζ growth signature)
   — node extraction + interpolation is clean and probe-validated at
   exactly those ζ (2⁻⁶⁰·⁶..2⁻⁶²). Probes: scratchpad\beta\probe9*.
2. **Check (e) extended pocket FAILED (first actual run — the prior
   pass died at (c) before reaching it): (1.6, 20, 0.4) value 0.99985,
   complement 1.52e-4 < 2⁻¹² doctrine bound, stay-R1 under the 1.5
   τ-gate.** The "τ > 1.5 safe band" claim was gamma-limit reasoning,
   wrong at moderate β (second time this class of error bit — see the
   τ-gate lattice bug). **EIGHTH correction: kBetaPrTauMax 1.5 → 2.5**
   + BetaR4Tiny's lgamma(1+τ) grows a third zone (one recurrence step:
   lgamma(τ) + ln τ, τ−2 Sterbenz-exact, LogDdAny; sum cannot cancel —
   lgamma ≥ −0.1215 vs ln τ ≥ 0.405). Bar (1−2⁻¹¹) < doctrine bound
   (1−2⁻¹²) ⇒ every violating τ ≤ 2.5 lane post-routes by
   construction; pocket alphas 2.6/2.8/3/4 police the remainder
   (box-corner Q ~ 7e-4 at 2.6, ~3× margin). Rejected alternatives:
   R1 swap-orientation (violates ξ₁ cap: 0.6⁶⁴ = 2⁻⁴⁷); routing-time
   R2-swap diversion (reintroduces the sixth-correction CF stall shape
   at ξ→1). LgammaDiffDd audit: m ≤ 256 covers τ ≤ 2.5; α = min holds
   at any gate (near-one ⇒ β > α(1−ξ₁)/ξ₁ > 1.22α).

RESUME STATE [2026-08-03, second breakpoint]:
- DONE this session: generator (P_EXT_GL/K_GL_EXT/BETA_PR_TAU_MAX,
  build_r3_gl_ext, eval_r3_S gl_ext, check (c) gl lattice w/ ext +
  ν 45/128 rows, check (f) α lattice to 2.5 + (1.6,20,0.4) witness,
  check (e) pocket + 2.6/2.8, kBetaR3GlExt/kBetaPrTauMax emission);
  kernel (BetaR3Temme i_gl arg + ext Clenshaw block, BetaR4Tiny
  three-zone lgamma, driver comment + call site); pool_rescue3.py
  (windowed-parallel, works: ~340/chunk, format verified).
- IN FLIGHT: generator run 3 (scratchpad\beta_data_final3.h +
  beta_gen_final3_stderr.txt) — expect ALL checks green now; then
  diff vs committed src/beta_data.h: pre-existing tables bit-identical,
  NEW kBetaR3GlExt[3][25] + kBetaPrTauMax now 2.5 + kBetaGammaLim/
  kBetaNearOne/kBetaGlRidgeMin as before. Install + commit with the
  kernel edits (they compile only together — kBetaR3GlExt).
  ALSO in flight: rescue loop (pool_rescue3.py × ≤16 chunks; was
  4,405 FAILED after two chunks from 4,872).
- NEXT: (1) finish rescue loop (leave genuinely-dead rows FAILED).
  (2) Invalidate checkpoint rows whose oracle dispatch changed:
  min ∈ (2⁻⁴, 2.5] & route_final tag R4-postroute — REMOVE the rows
  (compact the file, keep latest-per-idx) so the MAIN loop recomputes
  via the small-τ path, NOT by marking FAILED (prewarm would betainc
  them instead). NOTE the gate is now 2.5, not 1.5. (3) Reference
  regen chunks to completion; self-checks + coverage audit (histogram
  rows for R4-postroute and R2-gammalim). (4) Rebuild build-clangcl
  (vcvars64 import), smoke + ULP vs regenerated references — expect
  R1-cmp fixed (the 429-ULP cell was the routing hole), R4-postroute/
  R2-gammalim/R3-floor as own rows. (5) G4: pin gates to measured,
  tier order Kaby native+caps → NEON CI → Ryzen last, merge to main;
  G5 four-list audit + ACCURACY.md/README same change set.

#### Decisions made here / still open
Decided: no gamma-core dependency EXCEPT the (C) gamma-limit slice
(two template cores + the R3 ridge-floor extension); R2 covers the
gamma limit below B_GL AND the off-band far ridge at all ν (margin
improving with ν, measured); R1 near-one lanes post-route to R4's
assembly in the fired orientation (SEVENTH correction, supersedes the
sixth);
erf-form in (ζ, p) with 1/ν powers and symmetry e_k(ζ,p) = −e_k(−ζ,q)
over the RATIO-BAND domain; complement-slack doctrine applied to the
EVALUATED side with region-driven orientation in the G1b final order;
exact-c_dd rule for every lgamma-argument path; specials table above;
own TU. Pinned: N₁ = 64, N₂ = 64, T_ridge = 32, B₁ = 8, ξ₁ = 0.45,
ε_R4 = 2⁻⁶, R4 N = 48, Z₀ = 10, K_B = 16, C_lg = 256 (provisional),
E_floor = −800, R2 orientation rule ξ < (α+1)/(c+2), ratio caps
[½, 2], ζ_max = √(3 ln 2/2), R3 K = 10 @ 25×15 + gl-ext K 10..12 @ 25
(29.88 KiB total), η_γ = −ζ√2 mapping, kBetaPrTauMax = 2.5 (EIGHTH),
kBetaGammaLim = 2⁵⁹, kBetaNearOne = 1−2⁻¹¹, kBetaGlRidgeMin = 20,
P_EXT_GL = 2⁻²⁰ (extraction-side only). [OPEN — later gates]: BinetDd's kernel home (G3;
generator emits fresh coefficients either way); the c-overflow guard
site (G3); C_lg drop to 128 if G4 measures hot; G2 next — reference
generator + point sets.

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
Cross-repo effort tracked in the fleet standards repo
([record](https://github.com/OldCrow/standards/blob/main/records/BUILD-STANDARDIZATION-PLAN.md)).
corvus commits: pkg-config file + consumer example + installed-path CI
check; find_package(hwy 1.4) version floor with CI building pinned
Highway 1.4.0 from source; CMakePresets.json. AGENTS.md's CMake section
verified still accurate afterwards.

## Resolved log
One line per closed item; detail lives in this file's git history,
AGENTS.md, and docs/ACCURACY.md.
- 2026-07-28 Phase C part 2 (gamma_p/gamma_q) shipped: 2 ULP direct side
  on every region, all-relative gates, four x86 tiers cell-identical;
  Sonnet-tooling + Opus-kernel split with orchestrator review gates —
  what each gate caught is recorded in the Phase C part 2 shipped record.
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
