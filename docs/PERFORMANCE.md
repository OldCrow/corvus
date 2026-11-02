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
are not comparable (§6); §5's bands read in two speed modes from process
to process, cause unknown, so a single reading of them is provisional.

**v1.1.0 (2026-10-25): partly stale.** The figures below were measured on
v0.9.0 code. v1.1.0's accuracy fixes made several families slower, and
its elementary-family work (#43) made that family faster. The only v1.1.0
measurement is one ad-hoc Zen 4 run of every family bench on the freeze
code (`0508ef1`, `AVX3_ZEN4`, clang-cl Release) at 20–69% ambient load.
It is not a quiet-bench pass and is INDICATIVE only: differences below
about 10% are within its noise. Against §3/§9 (ns/el, corvus only):

| family | v0.9.0 quiet → v1.1.0 ad-hoc | cause |
|---|---|---|
| `erf`, `erfinv`/`erfcinv`, `lgamma`, `digamma`, `trigamma`, `lbeta` | unchanged within noise | — |
| `cos` / `sin` | 1.6 → 1.1 | #43 |
| `exp` | 2.1 → 1.5 | #43 |
| `log` | 4.9 → 2.7 (0.64× → 1.34–1.51× against UCRT) | #43 |
| `log1p` | 4.9 → 4.9–5.4 (0.84–0.93× against UCRT) | — |
| `erfc` tail / mixed | 4.1 → 9.9 / 6.3 → 14.2 | #74 (§3) |
| `gamma_p`/`gamma_q` R1 / R2 / R3 / R4 | 66 / 54 / 49 / 77 → 109 / 68 / 57 / 85 | #56, #57, #64, #67, #70 |
| `beta_p`/`beta_q` | R1, R2 and gamma-limit +32–37%; R3, R4, post-route +9–12% | #47, #54, #69, #84, #85 |
| `gamma_p_inv`/`gamma_q_inv` | +15–55% by band | the forward, plus #48, #72, #83 |
| `beta_p_inv`/`beta_q_inv` | +10–25% by band; joint-tiny plateau 529 → 829 | the forward, plus #50, #55, #60 |
| `i0`/`i1`/`i0e`/`i1e` | series 13.9–14.2 → 14.0–25.4; tail 13.9–14.2 → 6.4; mixed 13.9–14.2 → 19.3–19.8 | #75 (§4) |

The beta row ran at the 69% load peak, so its size is the least certain.
Kaby Lake and M1 were not re-measured; the same code changes apply there.
The §4 and §8.2 upper bounds are stale on two counts, these changes and
#42's per-element baseline (§2). Row D's quiet full-surface pass replaced
§4 on Zen 4 (`2f863be`, 2026-11-01) and §8.2's Kaby Lake column
(`0a3ee74`, same kernels, 2026-11-01); M1 is pending, and this note goes
when it is in.

**v1.2.0 row C1 (2026-10-29):** #42 lever 1 made `beta_p`/`beta_q` and
their inverses faster on Zen 4; quiet A/B in §12. Gamma is unchanged.
Reproduced on Kaby Lake 2026-10-30 (§12.1) and on the M1 2026-10-31
(§12.2). **Row C3 (2026-10-30):** #31's lgamma table band and single-log
reflection, Zen 4 in §13, Kaby Lake in §13.1, M1 in §13.2; the §8.1 and
§11.2 lgamma-against-Apple-libm ratios are superseded there.

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

**The elementary family** (`exp`, `log`, `log1p`, `cos`, `sin`, shipped
in v0.8.0) was measured in the 2026-08-30 quiet pass on this machine —
see §9.

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

**The per-element baseline changed on 2026-09-30 (#42).** Every "upper
bnd" on this page was measured while a length-1 call zero-padded its dead
lanes; the specials scrub turned those zeros into an interior safe point
whose region core then ran alongside the live lane's, so a single call cost
up to ~2.5× one vector (Kaby Lake AVX2: `gamma_p` 1,258 → 506 ns, `lgamma`
274 → 122 ns per call after the driver pads with a live element). The
batching-gain tables in §4 and §8.2 therefore overstate what a caller gains
today by roughly that factor for the families whose safe point lands in a
different region (`gamma_p`, `lgamma`, `digamma`/`trigamma` to be
re-measured; `beta_p` and the inverses were already at one vector per call).
They are left as measured; a v1.1.0 re-run replaces them.

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

**Stale since row 12h (2026-10-25).** The tail now evaluates the dd G
(#74, docs/ACCURACY.md erfc). An ad-hoc Zen 4 run, not a quiet-bench run,
put the tail-only band at ~9.8 ns/el against ~5.5 for the plain-double G,
and the mixed band at ~12.7 against ~8.8. The core band is unchanged. The
tail and mixed rows above, and in §8.1, predate the change; re-measure at
the next quiet bench.

### lgamma

| band | corvus ns/el | libm ns/el | ratio |
|---|---:|---:|---:|
| [0.5, 2.5] zone (the zeros) | 5.25 | 28.75 | 5.47× |
| [0.01, 100] mixed positive | 13.34 | 56.06 | 4.20× |
| [2.5, 8] recurrence | 13.24 | 37.55 | 2.84× |
| [8, 1000] Stirling | 6.81 | 13.36 | 1.96× |
| [−30, −0.01] reflection | 33.25 | 48.83 | 1.47× |

Best to worst spans nearly a factor of four, against one baseline on one
machine. That spread is the point of publishing per region at all: a single
number for "how much faster is lgamma" would hide it.

Note the reflection band's absolute cost (33 ns/el) — the negative axis runs
the positive pipeline and then a reflection on top, so it is genuinely more
work, not a missed optimisation.

---

## 4. Batching gain (corvus against corvus, per element)

Upper bounds, not comparisons with other libraries. Range across each family's
bands at the largest n. Zen 4 at the v1.2.0 freeze head (`2f863be`, `AVX3_ZEN4`,
clang-cl Release): two full `tools/quiet_bench.ps1` passes, both gated (evidence
in `docs/bench-evidence/2026-11-01-zen4-v1.2.0/`). Bands in §5's fast mode.

| family | upper bound |
|---|---|
| `gamma_p` / `gamma_q` | 6.02 – 8.45× |
| `beta_p` / `beta_q` | 6.23 – 8.05× |
| `digamma` | 6.09 – 12.66× |
| `trigamma` | 5.42 – 13.23× |
| `gamma_p_inv` / `gamma_q_inv` | 3.83 – 8.07× |
| `beta_p_inv` / `beta_q_inv` | 3.12 – 7.77× |
| `lbeta` | 6.0× main band, 8.4× big band |
| `erfinv` / `erfcinv` | 6.36 – 17.77× |
| `i0` / `i1` / `i0e` / `i1e` | 6.23 – 12.17× |

The per-element baseline changed with #42 (§2), so these do not compare with
the v0.9.0 figures, which ran to 141× against a slower baseline.

The `i1`/`i1e` series now reads with `i0`/`i0e`: 14.2 against 13.9 ns/el. The
v1.1.0 ad-hoc 25.4/20.8 was not a property of the kernel.

The very large figures at the top of the `gammainv` range are a property of the
baseline, not of the kernel: an iterative solve has a lot of per-call setup to
amortise, so calling it one element at a time is especially wasteful.

---

## 5. Bands read in two speed modes

Some bands read in one of two speeds from process to process, and the slow
one is 1.50× the fast one. The cause is not known. Six fresh processes of the
same binary at `2f863be`, Zen 4 (SIMD ns/el, evidence as §4):

| band | fast | slow | runs slow |
|---|---:|---:|---:|
| `digamma` (0,1) up-step | 8.77 | 13.17 | 2 of 6 |
| `trigamma` (0,1) up-step | 10.40 | 15.47 | 1 of 6 |
| `i0` series | 13.92 | 20.89 | 1 of 6 |

The scalar walk over the same data holds within 0.5% in those runs. Inputs are
seeded, so every run consumes identical data, and the harness warms up and
takes a median, so neither a cold start nor a stray sample survives. Different
bands land slow in different runs, in the same process, so ambient load, clock
and thermal effects are ruled out: those would move every band together.
What varies between processes is where the buffers land. Address-dependent
behaviour is the leading suspect, not a finding. The v0.9.0 `i0` outlier here
(5.62× against 8.36×) is the same ratio; `erfcinv`'s old outlier (36.78×
against 14.94×) has not recurred in these six runs (14.04–15.12 ns/el).

Read a band's fast mode as its speed, and repeat a surprising reading before
trusting it. #97 tracks the cause. Kaby Lake showed no slow mode in six
processes; the M1 showed none at 1.50× in thirteen, but five bands step
1.07–1.30× between processes with the scalar walk stepping on its own
(§8.3).

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

All three columns are v1.2.0 at the freeze kernels. Kaby Lake: `0a3ee74`
(kernels as `2f863be`), `AVX2`, AppleClang 15 Release, two full
`tools/quiet_bench.sh` passes gated at 4.04/4.00% and 2.03/2.29% (evidence in
`docs/bench-evidence/2026-11-01-kaby-v1.2.0/`). The column is pass 2, the
quieter pass; pass 1 agrees within 6% on every band but four, named below.
M1: `ac472c3` (same kernels), `NEON`, AppleClang 21 Release, eight full
passes under the §10 gate (evidence in
`docs/bench-evidence/2026-11-01-m1-v1.2.0/`). The column is pass 2 (gate
4.03/4.64%); pass 11 (2.40/3.14%) agrees within 3.5% on every band, five more
gated passes within 8.2%, and pass 1 within 6% on every band but eight, all
named below.

| family | upper bound (Kaby Lake, AVX2) | Zen 4 (§4) | M1 (NEON) |
|---|---|---|---|
| `gamma_p` / `gamma_q` | 3.13 – 3.92× | 6.02 – 8.45× | 1.78 – 2.00× |
| `beta_p` / `beta_q` | 3.21 – 4.02× | 6.23 – 8.05× | 1.73 – 2.00× |
| `digamma` | 2.88 – 4.77× | 6.09 – 12.66× | 1.70 – 2.10× |
| `trigamma` | 2.58 – 4.88× | 5.42 – 13.23× | 1.60 – 2.11× |
| `gamma_p_inv` / `gamma_q_inv` | 2.29 – 4.03× | 3.83 – 8.07× | 1.46 – 2.00× |
| `beta_p_inv` / `beta_q_inv` | 2.06 – 3.88× | 3.12 – 7.77× | 1.35 – 1.91× |
| `lbeta` | 3.16× main band, 3.97× big band | 6.0×, 8.4× | 1.66 – 1.70× main band, 2.02× big band |
| `erfinv` / `erfcinv` | 2.36 – 5.02× | 6.36 – 17.77× | 1.50 – 2.44× |
| `i0` / `i1` / `i0e` / `i1e` | 2.91 – 4.53× | 6.23 – 12.17× | 1.45 – 2.25× |

Pass 1 ran slower than pass 2 on `digamma` (0,1) up-step (SIMD 27.0 against
24.3 ns/el, scalar unchanged), `trigamma` [2,8) down-walk and [8,2^89)
asymptotic (54.1/20.6 against 47.9/16.7, scalar moved with it) and
`beta_p_inv` weak-seed (5,008 against 4,528). The four §8.3 repeats read
`digamma` and `trigamma` with pass 2 to 1%, so pass 1 is the slow reading.
`trigamma`'s [2^89,inf) 1/x-only band is left out of the range: pass 2 read
6.50× from a scalar walk at 108 ns/el against pass 1's 86 (4.87×), with the
SIMD side at 16.6 in both.

M1: pass 1 gated at 7% but from `bench_digamma` on sampled 50% with
`corespotlightd` at three to four cores, and every one of its eight bands
outside 6% is that load: `trigamma` [8,2^89) asymptotic 2.00/2.15×;
`gamma_p_inv` a ≥ a_T deep tail 1.65/1.42× and beyond resolution 1.76/1.99×
(pass 1's scalar walk 18% and 6% slower); `beta_p_inv` deep-small 1.83/1.95×
(`p_inv` 1,285 against 1,205 ns/el, the band §12.2 flags); `i0` series
2.23/1.91× and `i0e` series 2.07/1.92×; `i1` tail 2.42/2.24×; lbeta main
band 1.43/1.66× (356 against 308 ns/el). Each of those bands but the two
Bessel series reads within 3% of pass 2 in the six other gated passes (the
deep-tail `gamma_p_inv` bands within 0.3%, `beta_p_inv` deep-small
1,205–1,238 ns/el, lbeta main 300–308, `i1` tail 2.21–2.25×, `trigamma`
asymptotic 2.11–2.25× on a 31.2/32.9 step of the walk), so the §13.2 lbeta
flag did not recur in seven passes and the column's range ends stand; the
`i0` and `i0e` series are §8.3's modes. The one range end that is mode-dependent is
`digamma`'s top: 2.10× in the four passes whose scalar walk read 59.9–60.2
ns/el on the [1,2) zone and 2.25× in the three at 64.6, the scalar walk's
two modes (§8.3) with the SIMD side at 28.77 in all seven. Passes 8 and 12 ran
with `mediaanalysisd` at one to two cores before every bench (13–22%, the
E-core class §10 describes) and still agree with pass 2 within 7.8%, the
excursions being the `i0e` series mode and, in pass 8, the `i1e` mixed band
at 65.78 ns/el against 61.3–61.5 in eleven processes with the walk steady
at 90.4: the same stepping (§8.3).

Half of Zen 4 or less across the board, as the lane count predicts, and
further below it on `digamma`, `trigamma` and the Bessel family, where the
Zen 4 gains run to 12–13×. The ordering of families and of bands within a
family is preserved, so the structural conclusions of §4 hold on the second
microarchitecture; only the magnitudes are lane-bound. Against this machine's
v0.9.0 column (4.83–11.07× down to 2.35–27.67×, old baseline), every range
is lower for the reason §2 gives: #42 made the per-element baseline faster.
The M1 column is a quarter of Zen 4's, give or take, as two lanes against
eight predict, and the family and band ordering is preserved there too:
`erfinv` central and the Bessel tail at the top, the inverses' deep and
skewed bands and the mixed-sign Bessel band at the bottom.

### 8.3 §5 on this machine

v1.2.0 (`0a3ee74`, 2026-11-01): the two full passes plus four repeats of
`bench_digamma bench_trigamma bench_erfinv bench_bessel`, six fresh processes,
all gated (evidence as §8.2, `rep1`–`rep4`). §5's three bands (SIMD ns/el):

| band | six readings | slow/fast |
|---|---|---:|
| `digamma` (0,1) up-step | 27.03, 24.27, 24.26, 24.27, 24.34, 24.27 | 1.11× |
| `trigamma` (0,1) up-step | 30.63, 29.07, 29.06, 29.06, 29.06, 29.06 | 1.05× |
| `i0` series | 42.87, 40.76, 52.80, 40.67, 40.75, 40.67 | 1.30× |

No 1.50× mode. The `digamma` 1.11× is pass 1, which ran slow across
`digamma` and `trigamma` (§8.2). The `i0` 52.80 is repeat 1, whose gate
sampled 10.7% (`contactsd`) just before the Bessel bench and whose scalar
walk moved with it (235 against 185 ns/el, and the tail bands up to 1.45×):
ambient, not §5, which holds the scalar walk. The other three repeats agree
with pass 2 to 0.3% on all three bands. On Zen 4 the slow mode appeared in
1–2 of 6 processes; six here showed none. That does not clear §5 on Kaby
Lake, and it says nothing about the cause.

M1, v1.2.0 (`ac472c3`, 2026-11-01/02): eight full passes and six short
processes of `bench_digamma bench_trigamma bench_erfinv bench_bessel`, all
gated (evidence as §8.2; the short ones are `rep1`, `rep3`, `rep4`, `rep9`,
`rep11`, `rep12`; the three day-run repeats sampled under 4.2% before every
bench). Pass 1 ran its second half under `corespotlightd` at three to four
cores (§8.2) and is left out. Thirteen fresh processes, SIMD ns/el:

| band | thirteen readings | slow/fast |
|---|---|---:|
| `digamma` (0,1) up-step | 30.76–30.78 ×10, 25.75, 31.43, 33.37 | 1.30× |
| `digamma` [1,2) zone | 28.76–28.83 ×12, 25.51 | 1.13× |
| `trigamma` (0,1) up-step | 38.09–38.14 ×12, 41.36 | 1.08× |
| `i0` series | 59.16–59.31 ×9, 64.19–64.23 ×4 | 1.09× |
| `i0e` series | 51.32–51.51 ×9, 55.32–55.34 ×4 | 1.08× |
| `i1` series | 57.69–57.94 ×12, 61.69 | 1.07× |
| `i1e` series | 53.52–53.72 ×13 | 1.00× |
| `i1e` mixed sign | 61.34–61.46 ×11, 63.05, 65.78 | 1.07× |
| `erfinv` central, n = 10⁶ | 11.57–11.62 ×11, 11.82, 12.05 | 1.04× |

No 1.50× mode, and a mode nonetheless. Each band above but the `i1e`
series and `erfinv` takes two or three values across processes, the three repeats that
sampled under 4.2% before every bench included, and the values repeat to
the hundredth:
`i0` reads 59.2 or 64.2, never between. The scalar walk has modes of its own and is not in step
with the SIMD side: the `digamma` up-step walk reads 64.6 in five
processes, 69.1 in seven and 57.27 in one (the process with the 25.75), the
[1,2) zone walk 59.9–60.2 or 64.6–64.9, the `trigamma` up-step walk
75.7 ×9, 76.3, 77.0 and 80.95 ×2, and the `i0` walk 114.5–114.8 ×8 or 122.5
×5. The four processes with `i0` at 64.2 are the four with `i0e` at 55.3
and the `i0` walk at 122.5, so some modes travel together and others
(`digamma` SIMD against its walk) do not. Zen 4's three bands read 1.50×
with the walk steady; here the steps are 1.07–1.30× and the walk steps
too, which fits the buffer-address suspicion of §5 (the walk and the SIMD
loop read the same buffers) better than anything per-kernel. `erfinv` and the
`i1e` series never move. Pass 1's `i0` at 51.41, the one reading below 59, came
from the loaded pass and is not counted. Recorded for #97, not published as
speeds; the column uses the ratios, which move by at most the walk's step.


v0.9.0 (2026-08-23): `erfcinv` central [0.5, 1.5]: 21.27× (discarded noisy
pass) / 20.28× (quiet pass). `i0` series band: 4.90× / 4.79×. Neither family
produced an outlier here in two runs, under very different ambient load.

---

## 9. Same machine, second quiet pass (2026-08-30, v0.9.0)

Same box, tier, compiler and libm as §1; harness gate passed at 4.02%,
noise held 3.46–14.91% (average 5.74%). Raw per-target outputs and the
runner log: `docs/bench-evidence/2026-08-30-zen4-quiet/`. Full-surface
pass (all fifteen family benches plus the lgamma component instrument,
`bench_lgamma_components`), per the one-pass-per-machine rule: a
machine's published numbers come from a single full-surface run.

### 9.1 Elementary family, first numbers (against UCRT — §3 baseline)

| function | corvus ns/el | speedup |
|---|---:|---:|
| `cos` / `sin` | 1.6 | 5.5× streaming; 2.2–2.4× at n ≤ 10⁴ |
| `exp` | 2.1 | 1.19–1.22× |
| `log` | 4.9 | **0.61–0.64× — slower than UCRT** |
| `log1p` | 4.9 | 0.86–0.92× |

**Stale for v1.1.0 (2026-10-25).** #43 moved the elementary kernels. The
ad-hoc v1.1.0 run (Status note) reads `cos`/`sin` 1.1 ns/el, `exp` 1.5
and `log` 2.7, which is 1.34–1.51× against UCRT: faster than UCRT, where
this table has it slower. `log1p` is unchanged. The `log` paragraph
below describes v0.9.0.

`log` below 1.0 is the accuracy trade stated plainly: corvus `log` is
correctly rounded on every reference row (see ACCURACY.md), UCRT's is a
fast table method. There is no configuration in which corvus wins this
comparison without giving that up, and it does not intend to.

### 9.2 lgamma per-region rerun

Zone 5.2×, recurrence 2.9×, Stirling 1.9×, mixed 4.2×, reflection 1.5× —
no band below 1.0×, confirming §3's quiet-machine finding against an
earlier recurrence figure that was measured on a loaded machine and
retracted. The per-band component profile
(`bench_lgamma_components`, same evidence directory) attributes the
costs: the zone's dd lead ladder outweighs its 32-coefficient Horner
(2.2 vs 1.3 of 4.65 ns/el); the recurrence band splits between the dd
log and the zone floor (4.4 + 4.6 of 12.6); reflection's two extra dd
logs are 10.0 of 33.1 ns/el against an 18.3 ns positive pipeline. The
gather-weak Kaby half of that profile is what decides the follow-up
lgamma optimization work (§11.3).

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

**Status: RATIFIED 2026-08-31 by the maintainer.** The 5% gate stays the
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

**v0.9.0 figures (2026-10-25 note).** On Zen 4, v1.1.0 `log` beats UCRT
(1.34–1.51×, ad-hoc; Status note), so §11.1's "`log` loses everywhere"
no longer holds for that column. The Kaby and M1 columns were not
re-measured.

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

### 11.3 lgamma component attribution, three machines (ns/el, n = 10⁶)

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

The component study's exit questions, answered on every machine measured:

- **(a) Does the zone Horner dominate the zone? No, anywhere.** The dd
  lead-term ladder outweighs the Horner on all three machines; the
  per-lane coefficient selects are a minor share of the Horner itself.
  This retires the proposed Estrin rewrite of the zone Horner.
- **(b) Does the dd log dominate the recurrence band? On the
  gather-weak machine, nearly** — 41% on Kaby (18.5 of 44.7, one dd
  log costing 18–20 ns/el against Zen 4's ~5), 35% on Zen 4, 29% on
  M1. The rest is the zone floor; the walk itself is cheap everywhere.
- **(c) Do the two extra logs dominate reflection?** They are the
  largest addressable share — 30% / 38% / 33% — atop an irreducible
  positive pipeline at ~55% / 50% / 47%. A single-log rewrite
  saves ~15% on Zen 4 and ~18% on Kaby in this band.

These answers scope the follow-up lgamma optimization work: the
table-driven (5/2, X0) band and the
single reflection log proceed (both attack the dd-log share, largest
exactly where corvus is weakest — Apple-libm machines); the Estrin
zone leg is retired by (a).

---

## 12. Row C1, #42 lever 1: Zen 4 quiet A/B (2026-10-29, v1.2.0 dev)

Two commits on `dev/v1.2.0`: `6027e2d` (beta R2's backward CF in double
for its deep steps behind a per-lane F gate) and `0967994` (one
LgammaDiffDd for both analytic lgamma(min) forms). Both keep every ULP
gate cell's real maximum unchanged on five tiers (docs/ACCURACY.md, beta).

**Protocol.** `AVX3_ZEN4`, clang-cl Release. A = the four family benches
built from `3712cc3` (before both commits), B = the same benches at
`64a085b`; binaries hash-checked distinct. Six `tools/quiet_bench.ps1`
passes alternating A, B, A, B, A, B: gated at 3.5–4.6%, per-pass noise
averages 2.9–4.9% (one sample each in passes 4 and 6 reached 11%).
Median of three per variant; per-band spread under 1% except where noted.

| band (ns/el, n = 4001) | A | B | A/B |
|---|---:|---:|---:|
| `beta_p` R1 series | 354.4 | 268.8 | 1.32× |
| `beta_q` R2 backward CF | 517.1 | 279.9 | 1.85× |
| `beta_q` R2 gamma-limit | 687.0 | 450.6 | 1.52× |
| `beta_p` R3, `beta_q` R4, `beta_p` post-route | 476.5 / 264.0 / 522.7 | 476.7 / 264.2 / 518.4 | 1.00× |
| `beta_p_inv` by band | 532–9,740 | 380–8,370 | 1.05–1.40× |
| `beta_q_inv` by band | 810–10,680 | 658–10,790 | 0.99–1.23× |

Controls: every `gamma_p`/`gamma_q` band and every `gamma_p_inv`/`gamma_q_inv`
band reads 1.00× (0.99× for one q_inv band), and so do the beta regions
neither commit touches. Two `beta_q_inv` bands (skewed CF,
gamma-limit transfer) read 0.99–1.00×; their q-side forward work sits
outside both changes.

**#44 answered here.** A dd Gauss–Legendre node (log_dd + exp_dd + dd
multiply-accumulate) costs ~10 ns/el on this machine, so a 24-node rule
(251 ns/el, indicative) costs ten times gamma R2's whole 44-step CF.
corvus's region cores already run at fixed cost; the depth, not the
method, was the open question.

**Owed:** ~~the same A/B on M1 (NEON, the tier #42's targets name) and on
Kaby Lake (AVX2)~~ (§12.1, §12.2). #42's NEON targets (gamma_p <= 40,
beta_p <= 30 ns/el) are out of reach under the accuracy contract; gamma's
prefactor alone costs 29 ns/el here across 8 lanes, and on NEON B reads
`gamma_p` R1 at 272 and `beta_p` R1 at 823 ns/el (§12.2).

### 12.1 Kaby Lake AVX2 (2026-10-30)

Evidence: `docs/bench-evidence/2026-10-30-kaby-c1-ab/`.

**Protocol.** `AVX2`, MacBook Pro 2017 (i7-7820HQ, macOS 13), AppleClang
15.0.0 Release. A = the four family benches built from `3712cc3`, B = the
same benches at `0708fed` (the `dev/v1.2.0` head, so B also carries rows C2
and C3); binaries hash-checked distinct. Six `tools/quiet_bench.sh` passes
alternating A, B, A, B, A, B: gated at 3.2–3.6%, per-bench noise 3.1–5.5%
with one sample at 8.7%. Median of three per variant; spread under 5% on
the rows below except `beta_q` R2 gamma-limit (B 10.6%) and `beta_p` R1
(A 6.3%).

| band (ns/el, n = 4001) | A | B | A/B | Zen 4 (§12) |
|---|---:|---:|---:|---:|
| `beta_p` R1 series | 1,025.8 | 785.3 | 1.31× | 1.32× |
| `beta_q` R2 backward CF | 1,462.2 | 787.5 | 1.86× | 1.85× |
| `beta_q` R2 gamma-limit | 2,010.9 | 1,284.7 | 1.57× | 1.52× |
| `beta_p` R3, `beta_q` R4, `beta_p` post-route | 1,285.5 / 757.0 / 1,545.2 | 1,237.1 / 733.9 / 1,514.3 | 1.02–1.04× | 1.00× |
| `beta_p_inv` by band | 1,534–21,509 | 1,114–17,840 | 1.05–1.38× | 1.05–1.40× |
| `beta_q_inv` by band | 2,194–23,768 | 1,805–23,113 | 1.03–1.22× | 0.99–1.23× |

The three lever-1 bands reproduce Zen 4's ratios within 0.05×, on four
lanes instead of eight. Per element, Kaby Lake runs these bands 2.8–2.9×
slower than Zen 4.

Controls: `gamma_p`/`gamma_q` read 1.00–1.06× and `gamma_p_inv`/`gamma_q_inv`
1.03–1.12× by band. Those are not lever 1's: the C3 A/B (§13.1), whose B
is the same binary, reads the same 1.04–1.05× on gamma R1/R2 and
1.03–1.11× on the gammainv bands against an A that already has lever 1,
so they belong to #31's lgamma, which gamma's prefactor consumes. The beta
regions neither commit touches read 1.02–1.04×, within their spread.

### 12.2 Mac Mini M1 NEON (2026-10-31)

Evidence: `docs/bench-evidence/2026-10-31-m1-c1-ab/`. The tier #42's
targets name.

**Protocol.** `NEON`, Mac Mini M1 (macOS 27.2), AppleClang 21.0.0, Highway
1.4.0, Release. A = the four family benches built from `3712cc3`, B = the
same benches at `4f78527` (the `dev/v1.2.0` head; kernels and benches
identical to Kaby's `0708fed`, so B also carries rows C2 and C3); binaries
hash-checked distinct. Six `tools/quiet_bench.sh` passes attempted under
the §10 gate (10%, 30 min per pass), alternating A, B, A, B, A, B. Three
ran (1 A, 4 B, 5 A), gated at 3.1–6.7%, per-bench noise 2.7–9.2% in passes
1 and 4 and 11.7–15.7% in pass 5 (`mediaanalysisd` at one E-core); the
other three never met the gate (`mediaanalysisd` with VTDecoder/VTEncoder,
`corespotlightd`). Median of two for A. B is one run here; its §10 second
reading is the C3 A/B's B (§13.2: the same four binaries by hash, gated
twice more), which agrees within 0.5% on the beta rows and within 1% on
the betainv rows except the two cheapest `beta_p_inv` bands (deep-small,
R1-tiny) at 5–7%. A spread at most 1.5% on the rows below except
`beta_q_inv` gamma-limit transfer and large-side (7%).

| band (ns/el, n = 4001) | A | B | A/B | Zen 4 (§12) | Kaby (§12.1) |
|---|---:|---:|---:|---:|---:|
| `beta_p` R1 series | 1,098.8 | 823.5 | 1.33× | 1.32× | 1.31× |
| `beta_q` R2 backward CF | 1,622.4 | 850.4 | 1.91× | 1.85× | 1.86× |
| `beta_q` R2 gamma-limit | 2,066.8 | 1,308.5 | 1.58× | 1.52× | 1.57× |
| `beta_p` R3, `beta_q` R4, `beta_p` post-route | 1,433.2 / 844.7 / 1,764.5 | 1,431.2 / 845.0 / 1,748.2 | 1.00–1.01× | 1.00× | 1.02–1.04× |
| `beta_p_inv` by band | 1,725–20,741 | 1,271–18,335 | 1.05–1.36× | 1.05–1.40× | 1.05–1.38× |
| `beta_q_inv` by band | 2,250–22,809 | 1,854–20,790 | 1.05–1.21× | 0.99–1.23× | 1.03–1.22× |

The three lever-1 bands reproduce Zen 4's ratios within 0.06×, on two
lanes. Per element, B on the M1 runs these bands 2–8% slower than Kaby
Lake's four lanes and 2.9–3.1× slower than Zen 4's eight.

Controls: `gamma_p` R3 / `gamma_q` R4 read 1.00× / 0.99×, and every
`gamma_p_inv`/`gamma_q_inv` band 0.99–1.01× except the deep-small `p_inv`
band, flagged under §10 (184 / 278 ns/el over the two A runs, 266 for B;
the bench's shortest timing). `gamma_p` R1 and `gamma_q` R2 read 1.09× and
1.14×: not lever 1's, as on Kaby (§12.1). The C3 A/B (§13.2) reads the
same 1.09× / 1.14× against an A that already has lever 1, so it is #31's
lgamma through gamma's prefactor, a larger share on this machine than
Kaby's 1.04–1.05×. The beta regions neither commit touches read
1.00–1.01×.

## 13. Row C3, #31: Zen 4 quiet A/B (2026-10-30, v1.2.0 dev)

`5fef537` replaces lgamma's walk-down and dd log on 5/2 < x < 8 with a
Taylor table (11 gathers per vector) and takes the reflection's two logs as
one. Evidence: `docs/bench-evidence/2026-10-30-zen4-c3-ab/`.

**Protocol.** `AVX3_ZEN4`, clang-cl Release. A = the six family benches
built from `1fe08b1` (before #31), B = the same benches at `2001fc8`, whose
kernels match `5fef537`; binaries hash-checked distinct. Six
`tools/quiet_bench.ps1` passes alternating A, B, A, B, A, B: gated at
3.5–4.5%, per-pass noise averages 3.9–4.9%, except pass 2 (B) at 7.9%, with
one sample reaching 13.6%. Median of three per variant; spread at most 2.2%
on the rows below, except `gamma_q` R2 (3.8%).

| band (ns/el, n = 10^6) | A | B | A/B | UCRT |
|---|---:|---:|---:|---:|
| `lgamma` [1/2, 5/2], zone | 5.50 | 5.67 | 0.97× | 29.0 |
| `lgamma` [5/2, 8], table | 13.05 | 5.35 | 2.44× | 37.6 |
| `lgamma` [8, 1000], Stirling | 6.90 | 6.88 | 1.00× | 13.3 |
| `lgamma` [0.01, 100], mixed | 13.31 | 11.05 | 1.20× | 55.7 |
| `lgamma` [−30, −0.01], reflection | 33.11 | 24.16 | 1.37× | 48.5 |
| `gamma_p` R1 series | 102.4 | 97.7 | 1.05× | |
| `gamma_q` R2 backward CF | 63.3 | 58.4 | 1.09× | |
| `lbeta` main band (n as the bench sets it) | 106.8 | 109.8 | 0.97× | |

Controls: `gamma_p` R3, `gamma_q` R4, every `beta_p`/`beta_q` band, and every
inverse band read 0.99–1.02×. The UCRT column reads 0.99–1.00× between
variants, which bounds the machine's drift across the run.

**Costs.** The zone and lbeta's main band each read 3% slower. The cause is
not traced.

**Owed:** ~~the same A/B on Kaby Lake (11 gathers on a gather-weak part)~~
(§13.1) ~~and on M1~~ (§13.2).

### 13.1 Kaby Lake AVX2 (2026-10-30)

Evidence: `docs/bench-evidence/2026-10-30-kaby-c3-ab/`. The gather check
the Zen 4 row owed: Kaby Lake is the fleet's gather-weak part.

**Protocol.** `AVX2`, MacBook Pro 2017 (i7-7820HQ, macOS 13), AppleClang
15.0.0 Release, Apple libm. A = the six family benches built from
`1fe08b1` (before #31), B = the same benches at `0708fed`; binaries
hash-checked distinct. Six `tools/quiet_bench.sh` passes alternating A, B,
A, B, A, B: gated at 3.2–4.8%; pass 1 (A) saw 12–15% noise before its
beta, lbeta and gammainv benches, the rest 3.0–6.6%. Median of three per
variant; spread at most 6% on the lgamma rows below (B zone 5.3%, A mixed
6.1%); the beta rows carry 8–14% A spread from pass 1 and are listed as
controls only.

| band (ns/el, n = 10^6) | A | B | A/B | Apple libm | B / libm | §8.1 (v0.9.0) |
|---|---:|---:|---:|---:|---:|---:|
| `lgamma` [1/2, 5/2], zone | 16.47 | 16.59 | 0.99× | 15.8 | 0.95× | 0.93× |
| `lgamma` [5/2, 8], table | 44.69 | 27.78 | 1.61× | 12.2 | 0.44× | 0.28× |
| `lgamma` [8, 1000], Stirling | 25.20 | 24.95 | 1.01× | 20.9 | 0.84× | 0.84× |
| `lgamma` [0.01, 100], mixed | 36.94 | 34.98 | 1.06× | 20.2 | 0.58× | 0.57× |
| `lgamma` [−30, −0.01], reflection | 106.34 | 80.62 | 1.32× | 47.0 | 0.58× | 0.45× |
| `gamma_p` R1 series | 298.1 | 286.0 | 1.04× | | | |
| `gamma_q` R2 backward CF | 195.9 | 186.7 | 1.05× | | | |
| `lbeta` main band (n as the bench sets it) | 315.5 | 323.1 | 0.98× | | | |

The libm column reads 0.97–1.02× between variants, bounding the drift.

**The gathers did not lose.** The table band gains 1.61× here against
2.44× on Zen 4: the 11 gathers per vector cost on this part, but the band
still leaves the walk-down and dd log well behind. The accuracy gain
stands either way. Relative to Zen 4, Kaby Lake runs the table band 5.2×
slower per element, the zone 2.9× and Stirling 3.6×; the ~1.5–1.8× excess
on the table band over the gather-free bands is the gather penalty, and
it sizes #22 (row C4): a non-gather variant can recover at most that much
on this band on this class of part.

**Against Apple's libm.** The §8.1 statement that corvus `lgamma` is 3.6×
slower than Apple's in the recurrence band is superseded: 2.3× (0.44×)
after #31, and the reflection moves from 0.45× to 0.58×. Apple's lgamma
is still ahead in every band but the zone on this machine; the README's
no-figures positioning stands.

**The Zen 4 0.97× readings.** The zone reads 0.99× (spread 2.5% / 5.3%)
and lbeta's main band 0.98× (3.1% / 5.1%): the same sign as Zen 4, both
within spread on both machines. Recorded as within spread, no issue
(PLAN.md, 2026-10-30 decisions); a third machine with the same sign
reopens it.

Controls: `gamma_p` R3, `gamma_q` R4, lbeta's big band and every
`beta_p_inv`/`beta_q_inv` band read 0.98–1.02×; `gamma_p_inv`/`gamma_q_inv`
read 1.03–1.11× by band, the same uplift as §12.1 sees, through lgamma.


### 13.2 Mac Mini M1 NEON (2026-10-31)

Evidence: `docs/bench-evidence/2026-10-31-m1-c3-ab/`. Apple's arm64 libm,
two lanes: the machine on which corvus lgamma was furthest behind (§11.2).

**Protocol.** `NEON`, Mac Mini M1 (macOS 27.2), AppleClang 21.0.0, Highway
1.4.0, Release, Apple arm64 libm. A = the six family benches built from
`1fe08b1` (before #31), B = the same benches at `4f78527`; binaries
hash-checked distinct. Six `tools/quiet_bench.sh` passes attempted under
the §10 gate, alternating A, B, A, B, A, B. Four ran (2 B, 3 A, 4 B, 5 A),
gated at 2.5–9.5%, per-bench noise 2.9–5.8% before bench_lgamma and
7–45% before the later benches (pass 2 7–25%, passes 3–5 15–45%:
`mediaanalysisd` at one core, `corespotlightd` bursts, §10's E-core class
of load); passes 1 and 6 never met the gate. A seventh gated A pass from
an abandoned earlier run (`extra-A-1828/`, outside the medians) is the
third A reading: within 2% on the zone, mixed and reflection bands,
within 3% on every beta, gammainv and betainv band (the flagged gammainv
band apart), and 6–7% above on the table and Stirling bands and gamma
R1. Median of two per variant; spread at most 2.9% on the
lgamma rows below except reflection (B 6.9%: 54.2 / 58.1 ns/el).

| band (ns/el, n = 10^6) | A | B | A/B | Apple libm | B / libm | §11.2 (v0.9.0) | Zen 4 A/B | Kaby A/B |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `lgamma` [1/2, 5/2], zone | 24.16 | 24.50 | 0.99× | 11.58 | 0.47× | flagged | 0.97× | 0.99× |
| `lgamma` [5/2, 8], table | 42.61 | 6.16 | 6.92× | 8.59 | 1.39× | 0.20× | 2.44× | 1.61× |
| `lgamma` [8, 1000], Stirling | 21.46 | 21.30 | 1.01× | 8.62 | 0.40× | 0.39× | 1.00× | 1.01× |
| `lgamma` [0.01, 100], mixed | 27.95 | 23.42 | 1.19× | 9.61 | 0.41× | 0.34× | 1.20× | 1.06× |
| `lgamma` [−30, −0.01], reflection | 80.38 | 56.12 | 1.43× | 32.08 | 0.57× | 0.39× | 1.37× | 1.32× |
| `gamma_p` R1 series | 296.1 | 272.0 | 1.09× | | | | 1.05× | 1.04× |
| `gamma_q` R2 backward CF | 229.8 | 201.3 | 1.14× | | | | 1.09× | 1.05× |
| `lbeta` main band (n as the bench sets it) | flagged | flagged | — | | | | 0.97× | 0.98× |

The libm column reads 0.97–1.01× between variants, bounding the drift.

**The table band is the fleet's largest gain, and it passes Apple's
libm.** 42.61 → 6.16 ns/el. On two lanes the 11 gathers per vector are
22 scalar loads, cheap, while what the table replaces costs most on this
machine (§11.3: the dd log 12.4 ns/el and the walk 5.0 atop the zone's
floor). Per element, B runs the band within 15% of Zen 4's 5.35 on eight
lanes, where the gather-free bands run 2.3–4.3× slower than Zen 4 (zone
4.3×, Stirling 3.1×, reflection 2.3×). corvus lgamma is 1.39× faster
than Apple's arm64 lgamma on this band; it was 0.20× in v0.9.0 (§11.2).
Reflection moves from 0.39× to 0.57× of libm, mixed 0.34× → 0.41×.
Apple's is still ahead in every other band; the README's no-figures
positioning stands.

**The zone, third machine.** The zone reads 0.99× here (A 24.16 / 24.16,
B 24.50 / 24.50, extra A 24.62): the same sign as Zen 4 (0.97×) and Kaby
Lake (0.99×), with the B median inside A's three-reading range. §13.1's
rule (a third machine with the same sign reopens it) is met by sign and
not by magnitude; filed as #96, off milestone, with the bisect recipe. The zone row
itself reproduced to 2% across three A and two B readings, the first
time it has under §10, which names it the unstable row. lbeta's main band
is flagged on this machine: bimodal regardless of variant (A 293.8 /
302.9 / 347.9, B 292.4 / 347.2), so the M1 says nothing about its Zen 4
0.97×; the big band reads 1.00×.

Controls: `gamma_p` R3, `gamma_q` R4, lbeta's big band and every
`beta_p`/`beta_q` band read 1.00–1.02×; every inverse band 0.99–1.02×
except `beta_p_inv` deep-small at 1.05× (spreads 1.4 / 1.2%) and the
flagged `gamma_p_inv` deep-small band (§12.2). `gamma_p` R1 / `gamma_q` R2
read 1.09× / 1.14×, #31's lgamma through gamma's prefactor; the C1 A/B
(§12.2) reads the same from an A without lever 1, so lever 1 adds nothing
to gamma here either.

## 14. Row C3b, #94: Zen 4 quiet A/B (2026-11-01, v1.2.0 dev)

Evidence: `docs/bench-evidence/2026-11-01-zen4-c3b-ab/`; read it with
`docs/bench-evidence/2026-10-30-zen4-c3-ab/summarise.py <dir>`.

**Protocol.** `AVX3_ZEN4`, clang-cl Release, the six family benches. A is
`6684208`, B the same tree with BinetDiffDd's k = 2 term in dd; binaries
hash-checked distinct. Six `tools/quiet_bench.ps1` passes alternating A, B;
all six met the gate. Median of three per variant, ns/element.

| Band | A | B | A/B |
|---|---:|---:|---:|
| `beta_p` R1 | 268.91 | 273.53 | 0.98× |
| `beta_q` R2 | 281.00 | 284.58 | 0.99× |
| `beta_p` R3 | 476.58 | 476.56 | 1.00× |
| `beta_q` R4 tiny-min | 263.91 | 268.63 | 0.98× |
| `beta_p` R4 near-one | 518.50 | 528.32 | 0.98× |
| `beta_q` gamma-limit | 446.39 | 450.16 | 0.99× |
| `lbeta` main band | 110.37 | 112.32 | 0.98× |
| `lbeta` big band | 71.46 | 74.71 | 0.96× |

**Reading.** The dd term costs 1–2% wherever BinetDiffDd runs. The
`lbeta` big band's 0.96× has a 5.6% spread on A, so it says little beyond
that. `betainv` bands read 0.99–1.01× and `gammainv`, `gamma_p`/`gamma_q`
and `lgamma` read 1.00–1.01×. Not measured on Kaby Lake or M1: the
change adds no gather, and the cost is a few dd operations per lane.
