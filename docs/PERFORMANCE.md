# corvus — Performance

**Status: MEASURED, three microarchitectures (v0.9.0, 2026-08-31).**
Zen 4 / UCRT and Kaby Lake / Apple x86 libm under the 5% ambient gate;
Apple M1 / Apple arm64 libm under the §10 gate (ratified 2026-08-31),
stable rows only. **Start at §11** — the cross-machine synthesis is the
quotable part; §§1–10 are the dated measurement ledger behind it, kept
because a number without its provenance cannot be audited later.

What promotion does NOT change: `docs/ACCURACY.md` is the audited
document and nothing here has that standing; every ratio names the libm
and machine it was measured against, and ratios against different libms
are not comparable (§6); §5's two families still do not reproduce
run-to-run and stay provisional inside an otherwise settled document.

---

## 1. What was measured

| | |
|---|---|
| machine | Ryzen 7 (Zen 4), Windows 11 |
| SIMD tier | `AVX3_ZEN4`, asserted per target via `CORVUS_EXPECT_TARGET` |
| compiler | clang-cl, Release |
| libm | UCRT (`std::lgamma` etc. as clang-cl links them) |
| harness | `tools/quiet_bench.ps1` |
| date | 2026-08-15 |

The harness refuses to start until ambient CPU is under 5% across two
consecutive 10-second windows, and samples noise between every target. The
sweep below gated at 4.45% and held 3.06–9.10% (average 5.10%); the peak was
the first sample, with prior activity still decaying.

Each figure is the median of repeated timings after a warm-up call, on seeded
inputs — so the same binary on the same data gives the same answer, run to run,
except where §5 says otherwise.

**The v0.8.0 elementary family** (`exp`, `log`, `log1p`, `cos`, `sin`)
was measured in the 2026-08-30 quiet pass on this machine — see §9.

---

## 2. Two baselines, and why they must not be compared

This is the most important thing on the page.

**`erf`, `erfc` and `lgamma` are timed against the system libm.** Those ratios
mean what you would expect: corvus versus another implementation of the same
function.

**Every other family is timed against corvus called one element at a time.**
That is not a comparison with any other library. It measures how much you gain
by batching versus paying per-call overhead — corvus against itself. The
benchmarks label that column **"upper bnd"** rather than "speedup", and the
distinction is deliberate.

So a 141× on `gammainv` and a 1.96× on `lgamma` are not commensurable, and
ranking families across the two tables below would be meaningless. There is no
vendor implementation of the inverse incomplete gamma in the C runtime to
compare against, which is why the baseline differs in the first place.

---

## 3. Against the system libm

Real comparisons. n = 10⁶ unless noted.

### erf

No region split — this varies with **array size**, not with the domain:

| n | corvus ns/el | libm ns/el | ratio |
|---:|---:|---:|---:|
| 1 000 | 1.70 | 6.30 | 3.71× |
| 10 000 | 1.65 | 7.05 | 4.27× |
| 100 000 | 1.66 | 10.55 | 6.36× |
| 1 000 000 | 1.70 | 10.77 | 6.33× |

corvus is flat across sizes; the libm slows down as the array leaves cache.
Much of the spread here is that, not vectorisation.

### erfc

| band | corvus ns/el | libm ns/el | ratio |
|---|---:|---:|---:|
| [−6, 6] core-dominated | 1.94 | 18.16 | 9.34× |
| [6, 28] tail-only | 4.10 | 19.92 | 4.86× |
| [−6.5, 28] mixed | 6.27 | 19.57 | 3.12× |

The mixed band is slowest because lanes in a vector can land in different
regions, and the kernel pays for both.

### lgamma

| band | corvus ns/el | libm ns/el | ratio |
|---|---:|---:|---:|
| [0.5, 2.5] zone (the zeros) | 5.25 | 28.75 | 5.47× |
| [0.01, 100] mixed positive | 13.34 | 56.06 | 4.20× |
| [2.5, 8] recurrence | 13.24 | 37.55 | 2.84× |
| [8, 1000] Stirling | 6.81 | 13.36 | 1.96× |
| [−30, −0.01] reflection | 33.25 | 48.83 | 1.47× |

Best to worst spans nearly four times, against one baseline on one machine.
That spread is the point of publishing per region at all: a single number for
"how much faster is lgamma" would be hiding it.

Note the reflection band's absolute cost (33 ns/el) — the negative axis runs
the positive pipeline and then a reflection on top, so it is genuinely more
work, not a missed optimisation.

---

## 4. Batching gain (corvus against corvus, per element)

Upper bounds, not comparisons with other libraries. Range across each family's
bands at the largest n.

| family | upper bound |
|---|---|
| `gamma_p` / `gamma_q` | 8.79 – 23.89× |
| `beta_p` / `beta_q` | 7.28 – 15.14× |
| `digamma` | 7.13 – 25.12× |
| `trigamma` | 7.29 – 30.93× |
| `gamma_p_inv` / `gamma_q_inv` | 6.21 – 141.57× |
| `beta_p_inv` / `beta_q_inv` | 3.89 – 54.18× |
| `lbeta` | 7.99× main band, 12.73× big band |
| `erfinv` / `erfcinv` | 5.17 – 14.94× — see §5 |
| `i0` / `i1` / `i0e` / `i1e` | 8.31 – 8.45× — see §5 |

The very large figures at the top of the `gammainv` range are a property of the
baseline, not of the kernel: an iterative solve has a lot of per-call setup to
amortise, so calling it one element at a time is especially wasteful.

---

## 5. Two families are not reproducible yet

`erfinv` and `bessel` do not give the same answer run to run, and the cause is
not known. Three runs of each:

| measurement | run A (in sweep) | run B | run C |
|---|---:|---:|---:|
| `erfcinv` central [0.5, 1.5] | **36.78×** | 14.94× | 14.71× |
| `bessel` i0, series band | 8.36× | **5.62×** | 8.34× |

Two of three agree in each case and the third is well outside. The outliers
land in different runs for the two families, so whatever causes them is
sporadic rather than systematic.

What it is **not**: inputs are seeded, so every run consumes identical data;
the harness warms up and takes a median, so neither a cold start nor a single
stray sample survives; and in the same binary seconds apart, `erfinv`'s own
central band reproduced to three digits while `erfcinv`'s moved four-fold.
Ambient load, clock and thermal effects would move every band together. These
did not.

Until that is understood, **treat both rows in §4 as provisional even by the
standards of this provisional document.** The figures quoted are the modal
values, not the extremes.

---

## 6. What would change these numbers

Everything here is one point in a space with at least four axes, and three of
them are unexplored.

- **Another microarchitecture.** An earlier "wins from N lanes up" claim died
  on a Kaby Lake measurement. §8 is the second machine: a 4-lane AVX2 part
  from 2017, where the §4 batching gains run at roughly half these figures.
- **Another libm.** UCRT's `lgamma` is slow — 27–56 ns/el in two of its bands.
  Against glibc or Apple's the margins in §3 would narrow, and could plausibly
  invert somewhere. **It does: §8.** Apple's `lgamma` is 12–21 ns/el in the
  same bands and corvus loses to it everywhere but the zone. This is the
  largest single source of error on this page, as predicted.
- **Another compiler.** clang-cl only. MSVC cannot reach AVX-512 at all and
  would be measuring a different tier under the same name.
- **More runs.** §5 exists because two runs looked like agreement until a third
  disagreed.

A published performance claim needs at least the first two settled. **As of
v0.9.0 they are** — three microarchitectures, three vendor libms (§11) —
which is what promoted this document. The compiler axis and §5 remain open;
the README quotes ranges with their libm named, never a single number.

---

## 7. Reproducing this

```powershell
tools\quiet_bench.ps1 -Targets bench_lgamma -ExpectTarget AVX3_ZEN4
```

The harness will refuse to run on a busy machine and tell you what to quiet
first. Raw per-target output and the noise-annotated log land in the build
directory; keep them together, since a number without its noise context cannot
be audited later.

Benchmarks must come from a Release build. `corvus::active_target()` — which
every bench prints — is the only reliable way to know which tier produced a
number.

## 8. Second machine: Kaby Lake AVX2, Apple libm (2026-08-23)

| | |
|---|---|
| machine | i7-7820HQ (Kaby Lake, 4 cores / 8 threads), macOS 13.7 |
| SIMD tier | `AVX2`, asserted per target via `CORVUS_EXPECT_TARGET` |
| compiler | Apple Clang 15, Release (`release` preset, `-G Ninja`) |
| libm | Apple libSystem (`std::lgamma` etc. as AppleClang links them) |
| harness | `tools/quiet_bench.sh` (the Unix counterpart of `quiet_bench.ps1`, same protocol) |
| noise | gate 4.46% / 2.59%; ten targets measured at 1.5–3.7% ambient; `erf`/`erfc` re-run separately under a 2.13% / 2.57% gate at 2.65% / 2.38% (their first samples caught a `diagnostics_agent` burst at 6–10%). A prior full pass at 10–20% ambient (`mdworker` indexing) was discarded; its rows differed from the quiet pass by a median 8.9% per row and up to 44% on the first targets, which is the reason the gate exists. |

### 8.1 Against Apple's libm

| function | band | corvus ns/el | libm ns/el | ratio | Zen 4 / UCRT ratio (§3) |
|---|---|---:|---:|---:|---:|
| `erf` | n = 10⁶ | 6.71 | 23.88 | 3.56× | 6.33× |
| `erfc` | [−6, 6] core | 7.10 | 24.56 | 3.46× | 9.34× |
| `erfc` | [6, 28] tail | 16.26 | 25.11 | 1.54× | 4.86× |
| `erfc` | [−6.5, 28] mixed | 20.16 | 25.63 | 1.27× | 3.12× |
| `lgamma` | [0.5, 2.5] zone | 16.58 | 15.46 | 0.93× | 5.47× |
| `lgamma` | [0.01, 100] mixed | 34.99 | 19.95 | **0.57×** | 4.20× |
| `lgamma` | [2.5, 8] recurrence | 43.48 | 12.20 | **0.28×** | 2.84× |
| `lgamma` | [8, 1000] Stirling | 23.83 | 20.05 | **0.84×** | 1.96× |
| `lgamma` | [−30, −0.01] reflection | 102.48 | 45.78 | **0.45×** | 1.47× |

(n = 10⁶ rows; the n = 10⁴ rows agree within 10%.)

Two effects stack, and they pull in the same direction for `lgamma`:

- **Apple's `lgamma` is fast.** 12–21 ns/el in the bands where UCRT's is
  27–56. The vendor baseline moved by 2–3×; corvus's own per-element cost
  moved by the lane count.
- **Four lanes, not eight, on a 2017 core.** corvus `lgamma` costs 17–43
  ns/el here against 5–13 on Zen 4, a 3.2–3.3× ratio that is what one
  expects from half the lanes on a slower clock. `erf`/`erfc` show the
  same 3.5–4× per-element ratio against their Zen 4 cells.

So the §3 statement "lgamma is faster than the vendor libm in every band"
is a **Zen 4 / UCRT statement**. On Kaby Lake against Apple's libm, corvus
`lgamma` is slower in every band but the zone — by 3.6× in the recurrence
band. `erf`/`erfc` stay ahead of Apple's libm on both machines, by smaller
margins here. Any published positioning must name the libm it was measured
against; the README's decision to quote no figures stands.

### 8.2 Batching gain, corvus against corvus

| family | upper bound (Kaby Lake, AVX2) | Zen 4 (§4) |
|---|---|---|
| `gamma_p` / `gamma_q` | 4.83 – 11.07× | 8.79 – 23.89× |
| `beta_p` / `beta_q` | 3.85 – 7.95× | 7.28 – 15.14× |
| `digamma` | 4.27 – 11.92× | 7.13 – 25.12× |
| `trigamma` | 4.11 – 15.95× | 7.29 – 30.93× |
| `gamma_p_inv` / `gamma_q_inv` | 3.69 – 60.73× | 6.21 – 141.57× |
| `beta_p_inv` / `beta_q_inv` | 2.35 – 27.67× | 3.89 – 54.18× |
| `lbeta` | 4.27× main band, 6.18× big band | 7.99×, 12.73× |
| `erfinv` / `erfcinv` | 3.15 – 20.28× | 5.17 – 14.94× — see §5 |
| `i0` / `i1` / `i0e` / `i1e` | 4.75 – 6.49× | 8.31 – 8.45× — see §5 |

Roughly half of Zen 4 across the board, as the lane count predicts. The
ordering of families and of bands within a family is preserved, so the
structural conclusions of §4 hold on the second microarchitecture; only
the magnitudes are lane-bound.

### 8.3 §5 on this machine

`erfcinv` central [0.5, 1.5]: 21.27× (discarded noisy pass) / 20.28×
(quiet pass). `i0` series band: 4.90× / 4.79×. Neither family produced an
outlier here in two runs, under very different ambient load. That neither
confirms nor clears §5 — two runs looked like agreement on Zen 4 too — but
it does say the sporadic outlier is not a property of the kernels alone.

---

## 9. Same machine, second quiet pass (2026-08-30) — v0.9.0 S1

Same box, tier, compiler and libm as §1; harness gate passed at 4.02%,
noise held 3.46–14.91% (average 5.74%). Raw per-target outputs and the
runner log: `docs/bench-evidence/2026-08-30-zen4-quiet/`. Full-surface
pass (all fifteen family benches plus the #30 component instrument),
per the one-pass-per-machine rule on #24.

### 9.1 Elementary family, first numbers (against UCRT — §3 baseline)

| function | corvus ns/el | speedup |
|---|---:|---:|
| `cos` / `sin` | 1.6 | 5.5× streaming; 2.2–2.4× at n ≤ 10⁴ |
| `exp` | 2.1 | 1.19–1.22× |
| `log` | 4.9 | **0.61–0.64× — slower than UCRT** |
| `log1p` | 4.9 | 0.86–0.92× |

`log` below 1.0 is the accuracy trade stated plainly: corvus `log` is
correctly rounded on every reference row (see ACCURACY.md), UCRT's is a
fast table method. There is no configuration in which corvus wins this
comparison without giving that up, and it does not intend to.

### 9.2 lgamma per-region rerun

Zone 5.2×, recurrence 2.9×, Stirling 1.9×, mixed 4.2×, reflection 1.5× —
no band below 1.0×, confirming §3's quiet finding against the retracted
loaded-set recurrence figure. The #30 per-band component profile
(`bench_lgamma_components`, same evidence directory) attributes the
costs: the zone's dd lead ladder outweighs its 32-coefficient Horner
(2.2 vs 1.3 of 4.65 ns/el); the recurrence band splits between the dd
log and the zone floor (4.4 + 4.6 of 12.6); reflection's two extra dd
logs are 10.0 of 33.1 ns/el against an 18.3 ns positive pipeline. The
gather-weak Kaby half of that profile is what decides #31.

### 9.3 Reproducibility against §3/§4 (fifteen days apart)

Where the two passes overlap, figures agree to a few percent: gamma
49–74 → 49–77 ns/el, digamma 7.0–27.4 → 7.0–27.4, trigamma 5.7–27 →
5.6–26.9, lbeta main 106 → 106. On §5's two families: `erfcinv`
central reproduced the modal value (15.2× / 15.4× vs modal ~14.8×,
no outlier); `bessel` i0 series came in at 9.9× vs modal 8.35× —
outside the few-percent class, consistent with §5's caution. §5
stands unchanged.

---

## 10. Apple Silicon ambient gate — RATIFIED

**Status: RATIFIED 2026-08-31 (v0.9.0 S4, user).** The 5% gate stays the
universal standard; this section is the documented deviation for
heterogeneous-core Apple Silicon only. The M1 rows in §11 are published
under it, subject to its own two-run agreement requirement — which is
why the M1 elementary rows and the lgamma zone band are flagged rather
than published there.

**Proposal.** On Apple Silicon, gate at 10% ambient (same two-window
protocol), with two additions: per-target noise annotation stays
mandatory (the runner already records it), and a row is publishable
only if it reproduces within a few percent across two gated runs.
Rows that fail the agreement check are flagged, not published.

**Why the 5% gate is unreachable there.** macOS Tahoe's on-demand ML
services (`mediaanalysisd`, `photoanalysisd`) are SIP-protected,
relaunch on demand despite launchd disables, and set a practical
ambient floor. ~30 h of gate sampling on the Mac Mini M1 across
2026-08-28..31 (evidence: `docs/bench-evidence/2026-08-31-m1-indicative/`,
gate-abort logs included) produced a minimum of 4.30% and never two
consecutive sub-5% windows. Kaby Lake on macOS 13 passes the 5% gate
(§8), so this is a modern-macOS-services floor, not a macOS one.

**Why 10% is equivalent rigor, not lesser.** Global CPU% on M-class
parts averages the performance and efficiency clusters, and macOS
schedules background/utility-QoS work onto the E-cores: one
background task saturating an E-core reads as ~12.5% ambient while
the single-threaded bench, at default QoS, runs on a P core
essentially uncontended. On SMT Intel, the same 10% means contention
for the very cores, shared L3 and ring the bench occupies — which is
what §8's discarded 10–20% pass measured (median 8.9% per-row shift).
Empirically on the M1, three independent runs (2026-08-23 twice,
2026-08-31 once) agree to 2–3% on every stable row under ambient
anywhere from 6% to 50%.

**Caveats the proposal carries.** Unified memory is the leak in the
isolation: E-core video decode still consumes shared bandwidth, so
streaming rows (n = 10⁶) are more exposed than cache-resident ones —
the likely mechanism behind the one persistently unstable row (the
lgamma zone band: 0.35–0.95× across runs). The two-run agreement
requirement exists to catch exactly this class. A more principled
gate — P-cluster idle via `powermetrics` — needs root and new
tooling; noted as a possible refinement, not part of this proposal.

---

## 11. Fleet synthesis — v0.9.0 (2026-08-31)

The quotable section. Sources: §9 (Zen 4, 5% gate,
`docs/bench-evidence/2026-08-30-zen4-quiet/`), the Kaby v0.9.0 pass
(5% gate, `docs/bench-evidence/2026-08-30-kaby-quiet/`), the M1 pass
(§10 gate, `docs/bench-evidence/2026-08-31-m1-indicative/`). All rows
n = 10⁶, tier-asserted (`AVX3_ZEN4` / `AVX2` / `NEON`). Ratios in one
COLUMN share a libm; ratios across columns do not, and §2's rule
against comparing them still applies row by row.

### 11.1 Elementary family, three libms

| function | Zen 4 vs UCRT | Kaby vs Apple x86 | M1 vs Apple arm64 ‡ |
|---|---:|---:|---:|
| `cos` | 5.50× (1.6 ns/el) | 4.87× (4.4) | 3.27× (4.0) |
| `sin` | 5.58× (1.6) | 4.98× (4.4) | 3.21× (4.1) |
| `exp` | 1.20× (2.1) | 0.65× (8.4) | 1.28× (4.7) |
| `log` | 0.64× (4.9) | 0.27× (19.4) | 0.21× (13.5) |
| `log1p` | 0.92× (4.9) | 0.40× (19.4) | 0.32× (13.8) |

‡ Single gated run — the §10 two-run agreement requirement is not yet
met for the elementary family (it did not exist at the 2026-08-23 M1
runs), so this column is REPORTED, NOT PUBLISHED; it firms up at the
next M1 pass. Its shape is consistent across sizes 10⁵–10⁶.

Three positioning facts, stable across all three libms:

- **`cos`/`sin` win everywhere** — 3.2–5.6× against every vendor libm
  measured, full double range, 1 ULP. The strongest performance row
  corvus has.
- **`log` loses everywhere (0.2–0.64×), and that is the product
  working as designed**: corvus `log` is correctly rounded on every
  reference row; every vendor libm measured trades exactly that away
  for speed. `log1p` sits below parity for the same reason. There is
  no configuration in which corvus wins this without giving up the
  bound.
- **`exp` is libm-dependent in both directions** (1.2× vs UCRT and
  Apple arm64, 0.65× vs Apple's x86 vectorized exp) — the cleanest
  demonstration that the baseline, not the kernel, moves these ratios.

### 11.2 lgamma across the fleet, and the Apple-libm question

§8 left open whether the lgamma inversion was an Apple-libm fact or a
Kaby-lane fact. **Answered: it is an Apple-libm fact.** M1 stable rows
(two-run agreement met): recurrence 0.20×, Stirling 0.39×, mixed
0.34×, reflection 0.39× — the same inversion as Kaby's 0.28×/0.77–
0.84×/0.55×/0.44× on a different microarchitecture with eight-lane
NEON-class arithmetic replaced by Apple's arm64 libm at 8.4–32 ns/el.
(M1 zone: flagged unstable under §10, not published.) Meanwhile Zen 4
vs UCRT stays above 1.0× in every band (§9.2). One function, one
kernel, and the vendor baseline alone decides which side of 1.0 it
lands on — this is why no README figure ever omits the libm name.

### 11.3 #30 component attribution, three machines (ns/el, n = 10⁶)

| component | Zen 4 | Kaby | M1 ‡ |
|---|---:|---:|---:|
| zone: full band | 4.7 | 14.9 | 17.7 |
| zone: 32-coef Horner (with selects) | 1.3 | 5.6 | 5.7 |
| zone: dd lead ladder | 2.2 | 6.7 | 6.2 |
| recurrence: full band | 12.6 | 44.7 | 42.1 |
| recurrence: masked walk-down | 1.4 | 3.8 | 5.0 |
| recurrence: outlined dd log | 4.4 | 18.5 | 12.4 |
| Stirling: full band | 6.8 | 24.6 | 19.9 |
| Stirling: dd log | 3.0 | 11.6 | 6.7 |
| Stirling: 1/x² remainder | 0.5 | 1.3 | 0.6 |
| reflection: full band | 33.1 | 107.6 | 81.7 |
| reflection: positive pipeline | 18.3 | 53.8 | 38.3 |
| reflection: two extra dd logs | 10.0 | 40.4 | 27.3 |
| reflection: LogSinc | 2.0 | 6.2 | 6.1 |

‡ M1 column single-run (same caveat as §11.1); its ratios agree with
the other machines' shape.

The #30 exit questions, answered on every machine measured:

- **(a) Does the zone Horner dominate the zone? No, anywhere.** The dd
  lead-term ladder outweighs the Horner on all three machines; the
  per-lane coefficient selects are a minor share of the Horner itself.
  This retires the Estrin-zone leg of #31.
- **(b) Does the dd log dominate the recurrence band? On the
  gather-weak machine, nearly** — 41% on Kaby (18.5 of 44.7, one dd
  log costing 18–20 ns/el against Zen 4's ~5), 35% on Zen 4, 29% on
  M1. The rest is the zone floor; the walk itself is cheap everywhere.
- **(c) Do the two extra logs dominate reflection?** They are the
  largest addressable share — 30% / 38% / 33% — atop an irreducible
  positive pipeline at ~55% / 50% / 47%. A single-log rewrite (#31)
  saves ~15% on Zen 4 and ~18% on Kaby in this band.

These answers scope #31: the table-driven (5/2, X0) band and the
single reflection log proceed (both attack the dd-log share, largest
exactly where corvus is weakest — Apple-libm machines); the Estrin
zone leg is retired by (a).
