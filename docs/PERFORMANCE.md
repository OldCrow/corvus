# corvus — Performance (provisional)

**Status: FIRST PASS. These numbers will change.** One machine, one compiler,
one libm, and two families whose figures are not yet reproducible. This is the
working record, not a claim — see §6 before quoting anything from it.

`docs/ACCURACY.md` is the audited document. Nothing here has that standing.

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

A published performance claim needs at least the first two settled. Until then
this page is a record of what one machine did on one afternoon, and the README
deliberately quotes no figures at all.

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
