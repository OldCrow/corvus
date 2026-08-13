# corvus — Numerical Doctrine

Loaded on demand via AGENTS.md's reading map. BINDING for kernel,
generator, reference/oracle, and accuracy work. For machines, builds,
tier capping, and CI see `docs/ENVIRONMENT.md`; for session state and
family design records see `PLAN.md`.

## Kernel construction

- **TU boundary = sharing/dependency boundary**, not one-symbol-per-file:
  functions that consume the same kernel cores stay in one TU with
  multiple HWY_EXPORTs (erfinv + erfcinv both route onto both shared
  cores — splitting would instantiate every core twice per target for
  nothing); split within a family only when dependency sets differ
  materially (erf.cpp stays free of erfc's dd/exp_dd/tail-data
  dependencies, and a consumer linking only erf pulls only erf.o).
- Per-target pattern: `HWY_TARGET_INCLUDE` + `foreach_target.h`, kernel in
  `corvus::HWY_NAMESPACE` written against `ops::`, then `HWY_ONCE` section
  with `HWY_EXPORT` + public dispatch wrapper. Call `HWY_DYNAMIC_DISPATCH`
  from *inside* `namespace corvus` — with a single compiled target (the
  SSE2 cap) Highway collapses it to `N_SSE2::FUNC`, and a globally
  qualified call then names a namespace that does not exist. It compiles
  at every other tier, so the cap sweep is what catches it.
- **Outline heavy kernels from day one**: region cores AND the per-lane
  driver of any non-trivial family are `HWY_NOINLINE` (gamma and erfinv
  are the pattern). MSVC's optimizer is superlinear in function size —
  fully inlined, each export becomes one enormous function per target
  (the driver inlines TWICE per export: full-vector + masked-tail call
  sites) and cl.exe codegen for one such TU reached ~15–17 minutes,
  killing the CI Windows job at its 25-minute timeout twice (2026-07-29)
  before outlining brought the whole library to ~4 minutes. Bit-identity
  is guaranteed by contraction-off (verify anyway: ULP tables
  byte-compare across the boundary change). Keep genuinely small hot
  helpers inline — outlining erfinv's central polynomial cost a measured
  0.5–0.7 ns/el on its 3 ns/el fast path for no meaningful codegen
  relief, and was reverted. MSVC additionally gets
  `/d2ReducedOptimizeHugeFunctions` on the heaviest TUs (gamma.cpp,
  beta.cpp; CMakeLists; real MSVC only, clang-cl rejects /d2 flags) —
  erfinv needed no such flag once outlined.
- dd layer: `src/dd-inl.h` holds the double-double primitives (Fast2Sum,
  TwoSum, TwoProd, DdAdd/DdMul/DdRecip, DdSqrt, DdRecipDd), written
  against `ops::` like everything else. Exact residuals go through
  `ops::ProdLow`, never a bare `MulSub` — same FMA-capability hazard as
  `ops::SquareLow`. `src/dd_special-inl.h` holds shared dd specials one
  layer up (Log1pmxDd, Expm1Dd), consumed by gamma and beta, gated by
  `test_dd_special`. `src/<fn>_dd-inl.h` are corvus-owned transcendental
  cores (exp_dd, log_dd), internal only: they return mantissa + exponent
  so a consumer folds its own factors in before the power-of-two scaling
  rounds anything — that is what keeps a subnormal result at one rounding.

## Numerical hazard rules

- Vector and tail paths must be the same code path (masked LoadN/StoreN),
  never a scalar libm fallback for the tail.
- Any op whose CORRECTNESS depends on FMA fusion (exact residuals like
  fma(a,b,-fl(a*b))) must be capability-guarded in the facade — Highway
  emulates MulAdd/MulSub as mul-then-add on non-FMA targets (SSE2/SSSE3/
  SSE4), which silently zeroes exact residuals. See ops::SquareLow.
- Masked-off lanes still EXECUTE every op, including gathers. Any gather
  whose index derives from lane VALUES (e.g. erf's round(ac*256); unlike
  log_dd's bit-masked slot index, which is bounded by construction) must
  have its input NaN/domain-scrubbed first — a discarded lane's NaN
  otherwise reaches the index or not by platform accident (x86 minpd
  drops NaN, ARM fmin propagates it and fcvtzs(NaN) = 0), and Highway's
  debug bounds assert can trip. erfc.cpp's nan mask and erfinv's
  HalleyMid scrub are the pattern.
- The non-FMA fallback in ops::SquareLow/ProdLow is Dekker's split, and its
  intermediate a*(2^27+1) OVERFLOWS for |a| > 2^996 (~6.7e299). A kernel
  whose operands can reach that range must scale by a power of two first
  and scale back after — exact, and linear in the operand, so it stays one
  code path for every target rather than a non-FMA special case. lgamma's
  Stirling product is the first kernel to need it; it read as correct on
  every FMA target and only the capped SSE sweep exposed it.
- FP contraction is OFF project-wide (`CORVUS_FP_FLAGS`): the dd exactness
  proofs assume IEEE ops as written; fusion is requested in source
  (`ops::MulAdd`), never inferred by the compiler.

## Test doctrine

- ctest executables compare against libm/reference values; test lengths
  are deliberately non-multiples of lane counts to exercise the
  masked-tail path.
- **Tests are registered in DEPENDENCY order, not development order** (dd
  cores first, then families in the order they consume each other), so a
  shared-core regression fails under its own name rather than a
  consumer's.
- **A new function family must be added in its dependency position to all
  FOUR explicit lists**: tests/CMakeLists.txt, the three `ULP report`
  steps in .github/workflows/ci.yml, and the `$gates` array in
  tools/sweep_tiers.ps1 — ctest picks new tests up automatically, but the
  report steps and the Ryzen sweep enumerate binaries explicitly and
  silently omit anything not added (the gamma pair shipped one commit
  without its CI reports before this rule existed).
- A test for an *internal* kernel compiles the kernel header itself
  through foreach_target and so uses `corvus_kernel_test_target()`, which
  links `hwy::hwy`, adds the source root, and — the part that matters —
  applies `CORVUS_HWY_TARGET_DEFS` so the test sees the same target set as
  the library. Such a test also asserts its own dispatched target equals
  `corvus::active_target()`.
- Where the kernel carries more than working precision, the reference file
  carries a double-double pair and the test measures relative error below
  the last bit of a double; rounding first would hide what is being
  tested.

## Generators and reference sets

Generators need mpmath in a throwaway venv (network for pip only):
```sh
python3 -m venv /tmp/mpv && /tmp/mpv/bin/pip install mpmath
/tmp/mpv/bin/python tools/gen_erf_table.py        > src/erf_data.inc
/tmp/mpv/bin/python tools/gen_erfc_tail_poly.py   > src/erfc_tail_data.h
/tmp/mpv/bin/python tools/gen_erf_reference.py    > tests/data/erf_reference.txt
/tmp/mpv/bin/python tools/gen_erfc_reference.py   > tests/data/erfc_reference.txt
/tmp/mpv/bin/python tools/gen_exp_table.py        > src/exp_dd_data.inc
/tmp/mpv/bin/python tools/gen_exp_dd_reference.py > tests/data/exp_dd_reference.txt
/tmp/mpv/bin/python tools/gen_log_table.py        > src/log_dd_data.inc
/tmp/mpv/bin/python tools/gen_log_dd_reference.py > tests/data/log_dd_reference.txt
/tmp/mpv/bin/python tools/gen_lgamma_data.py      > src/lgamma_data.h
/tmp/mpv/bin/python tools/gen_lgamma_reference.py > tests/data/lgamma_reference.txt
/tmp/mpv/bin/python tools/gen_erfinv_data.py      > src/erfinv_data.h
/tmp/mpv/bin/python tools/gen_erfinv_reference.py
/tmp/mpv/bin/python tools/gen_dd_special_data.py  > src/dd_special_data.h
/tmp/mpv/bin/python tools/gen_gamma_data.py       > src/gamma_data.h
/tmp/mpv/bin/python tools/gen_gamma_reference.py
/tmp/mpv/bin/python tools/gen_beta_data.py        > src/beta_data.h
/tmp/mpv/bin/python tools/gen_beta_reference.py
/tmp/mpv/bin/python tools/gen_digamma_data.py     > src/digamma_data.h
/tmp/mpv/bin/python tools/gen_digamma_reference.py > tests/data/digamma_reference.txt
/tmp/mpv/bin/python tools/gen_trigamma_data.py    > src/trigamma_data.h
/tmp/mpv/bin/python tools/gen_trigamma_reference.py > tests/data/trigamma_reference.txt
/tmp/mpv/bin/python tools/gen_gammainv_data.py    > src/gammainv_data.h
/tmp/mpv/bin/python tools/gen_gammainv_reference.py
/tmp/mpv/bin/python tools/gen_betainv_data.py     > src/betainv_data.h
/tmp/mpv/bin/python tools/gen_betainv_reference.py
/tmp/mpv/bin/python tools/gen_bessel_data.py      > src/bessel_data.h
/tmp/mpv/bin/python tools/gen_bessel_reference.py
/tmp/mpv/bin/python tools/gen_lbeta_reference.py
```
`gen_erfinv_reference.py` writes both `tests/data/erfinv_reference.txt` and
`tests/data/erfcinv_reference.txt` directly (two output files, so no `>`
redirection) rather than printing one file to stdout; `gen_gamma_reference.py`
likewise writes its three files directly (`gamma_p`/`gamma_q`/`dd_special`
reference — the last gates `Log1pmxDd`, which lives in `src/dd_special-inl.h`
but keeps its reference generation in the gamma generator because it shares
the seeded rng stream with the P/Q point sets; carving it out would silently
change every point). The gamma reference oracle has two non-obvious rules baked in —
compute the SMALL side of P/Q directly (never 1 − a near-1 value), and use
the exact-asymptotic oracle rather than mpmath's gammainc for a > 1e4
(mpmath's regularized lower hangs or fails to converge at large a) — both
documented at the enforcement sites in the generator.
`gen_beta_reference.py` writes both beta reference files directly
(checkpointed and resumable — re-run until it writes the files), and any
beta reference regeneration MUST end with a clean
`tools/verify_beta_reference.py` pass before the files are trusted: the
harness is oracle-independent and carries a baked-in negative control
(known-bad rows that must be rejected, exit 2 otherwise) — the
oracle-trust doctrine's enforcement point (see PLAN.md Decisions).
`gen_gammainv_reference.py` and `gen_betainv_reference.py` write their
reference files directly (checkpointed/resumable) with per-row bracket
certification and baked-in negative controls — same exit-2 doctrine.
`gen_bessel_reference.py` and `gen_lbeta_reference.py` likewise write
their files directly (layered-dps oracles, no bracket certification —
both functions have trusted single-argument mpmath baselines).

**Replay self-checks (binding, both rules paid for in shipped defects):**
- Every generator replay/self-check samples DOMAIN BOUNDARIES with
  edge-refined, bit-stepped points — uniform or random grids miss worst
  points that sit within ~1e-10 of an interval edge (coherent Chebyshev
  coefficient rounding puts them there), and under-pin dd-lead counts.
- A replay sim must model UNFUSED MulAdd semantics as well as FMA, and
  must replay EVERY shipped assembly (scaled and unscaled) — non-FMA
  tiers (SSE4/SSSE3/SSE2) are first-class, and FMA-only validation is
  the "assert the tier, never assume it" trap in numerical form.
- A truncation/rounding budget proven relative to a COMPONENT's own
  scale is void once the consuming assembly cancels below it — required
  accuracy is set by the result's cancellation depth, not the
  component's.

Reference files and generated tables are checked in; regenerate only when
the method or point selection changes, and re-run the ULP tests after.
A resumable generator's checkpoint signature must bind to the POINT
IDENTITIES (a bits-digest), not just the count: an edit that preserves N
replays stale oracle values under new point identities — caught live
2026-08-12 when a b-column swap in the pb-corner family kept N = 786 and
the first re-run served b = 1e20 values as b = 1e30 rows.
Table generators self-check on every run and exit non-zero rather than emit
a table that misses its error budget — `gen_exp_table.py` re-derives the
whole budget (reduction exactness, polynomial truncation, table dd error)
onto stderr. Trust that line over any claim in a comment.

**Oracle-trust doctrine** (canonical statement in PLAN.md Decisions):
reference-oracle construction for a function WITHOUT a trusted library
baseline is FRONTIER work. References are trusted only after independent
verification — a separate harness with baked-in negative controls
(tools/verify_beta_reference.py is the pattern) or per-row certification
designed at the frontier (the gammainv bracket-certification pattern).

## Model & effort routing

Hints for agent harnesses (sub-agent partitioning) and for flagging when
to switch model/effort — derived from how this project was actually built
(scaffold, facade, erf/erfc kernels, and governance were Fable 5 at high
effort; much of the surrounding work did not need that).

**High effort, frontier model** — wrong judgment is expensive to discover:
- New function/kernel design: approximation method choice (table vs
  polynomial vs rational), region splits, interval and degree selection,
  series length — and the error-budget analysis justifying them
  (cancellation, relative-vs-absolute error metrics, exactness arguments:
  Sterbenz, Fast2Sum, Dekker, FMA dependence).
- Diagnosing accuracy regressions from symptom to root cause (precedents:
  the erfc series-truncation metric mismatch; the non-FMA MulSub
  zero-residual hazard — both multi-step numerical reasoning).
- Reference-oracle construction wherever no trusted library baseline
  exists (precedent: the beta oracle — mpmath betainc returns
  internally-consistent garbage in whole parameter regions, and six
  oracle defect classes cost more sessions than the kernel itself;
  "it's just a generator script" was the misjudgment that caused it).
- Public API changes, facade/seam design, governance and policy decisions.

**Default effort, mid-tier model** — pattern-following with judgment:
- Implementing a function whose design is settled, through the
  established pipeline (generator → kernel via ops:: → smoke test →
  reference set → ULP gate → bench).
- New tests/generators/benches copying existing patterns; CMake/CI
  adjustments within documented policy.

**Low effort / small model** — documented recipes and mechanical edits:
- Running the per-tier validation recipe on another machine and recording
  results in ACCURACY.md/PLAN.md.
- Regenerating tables/references unchanged; renames; badge/link fixes;
  PLAN.md bookkeeping after decisions are made.

**Escalate** (or flag a model/effort switch to the user) the moment a
recipe task surfaces a decision: a ULP gate trips, measured values differ
across tiers, a fit misses its target, or Highway behaves contrary to
assumptions. **De-escalate** (flag it) when a high-effort session reaches
pure execution of a settled plan — don't spend frontier effort on rote
work.
