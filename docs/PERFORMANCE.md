# corvus — Performance

**Status: MEASURED at v1.2.0, three microarchitectures (2026-11-02).**
Zen 4 / UCRT and Kaby Lake / Apple x86 libm under the 5% ambient gate;
Apple M1 / Apple arm64 libm under the §10 gate. §§3, 4, 8 and 11 come from
quiet full-surface passes at the v1.2.0 freeze kernels (`2f863be`), two or
more per machine. A fourth machine, Ivy Bridge / glibc at `SSE4`, joined
on 2026-11-11 at the same kernels (§15). **Start at §11**, the
cross-machine synthesis. The rest is
the dated measurement ledger behind it, kept because a number without its
provenance cannot be audited later; §§12–14 are A/B measurements of single
v1.2.0 changes.

What this page is not: `docs/ACCURACY.md` is the audited document and
nothing here has that standing; every ratio names the libm and machine it
was measured against, and ratios against different libms are not
comparable (§6); §5's bands read in two speed modes from process to
process, cause unknown (#97), so a single reading of them is provisional.

---

## 1. What was measured

| | |
|---|---|
| machine | Ryzen 7 (Zen 4), Windows 11 |
| SIMD tier | `AVX3_ZEN4`, asserted per target via `CORVUS_EXPECT_TARGET` |
| compiler | clang-cl, Release |
| libm | UCRT (`std::lgamma` etc. as clang-cl links them) |
| harness | `tools/quiet_bench.ps1` |
| code | v1.2.0 freeze kernels, `2f863be` |
| date | 2026-11-01 |

The harness refuses to start until ambient CPU is under 5% across two
consecutive 10-second windows, and samples noise between every target. The
two full passes gated at 4.77% and 3.92% and held 2.3–7.3% (averages 3.7%
and 3.9%); they agree within 5% on every §3 row. Evidence:
`docs/bench-evidence/2026-11-01-zen4-v1.2.0/`; §3 and §4 quote pass 2.

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

So an 8× on `gammainv` and a 1.9× on `lgamma` are not commensurable, and
ranking families across the two tables below would be meaningless. There is no
vendor implementation of the inverse incomplete gamma in the C runtime to
compare against, which is why the baseline differs in the first place.

**The per-element baseline changed on 2026-09-30 (#42).** Before it, a
length-1 call zero-padded its dead lanes; the specials scrub turned those
zeros into an interior safe point whose region core then ran alongside the
live lane's, so a single call cost up to ~2.5× one vector (Kaby Lake AVX2:
`gamma_p` 1,258 → 506 ns, `lgamma` 274 → 122 ns per call after the driver
pads with a live element). The v1.2.0 tables in §4 and §8.2 use the new
baseline; the v0.9.0 figures, which ran to 141×, do not compare with them.

---

## 3. Against the system libm

Real comparisons. n = 10⁶ unless noted.

### erf

No region split — this varies with **array size**, not with the domain:

| n | corvus ns/el | libm ns/el | ratio |
|---:|---:|---:|---:|
| 1 000 | 1.70 | 6.20 | 3.65× |
| 10 000 | 1.67 | 6.92 | 4.14× |
| 100 000 | 1.67 | 10.44 | 6.24× |
| 1 000 000 | 1.67 | 10.51 | 6.28× |

corvus is flat across sizes; the libm slows down as the array leaves cache.
Much of the spread here is that, not vectorisation.

### erfc

| band | corvus ns/el | libm ns/el | ratio |
|---|---:|---:|---:|
| [−6, 6] core-dominated | 2.41 | 18.01 | 7.47× |
| [6, 28] tail-only | 9.74 | 19.84 | 2.04× |
| [−6.5, 28] mixed | 12.68 | 19.48 | 1.54× |

The tail evaluates a double-double G (#74, docs/ACCURACY.md erfc), which
costs about 1.8× the plain-double G it replaced in v1.1.0. The mixed band
is slowest because lanes in a vector can land in different regions, and the
kernel pays for both.

### lgamma

| band | corvus ns/el | libm ns/el | ratio |
|---|---:|---:|---:|
| [0.5, 2.5] zone (the zeros) | 5.67 | 28.66 | 5.06× |
| [2.5, 8] table | 5.33 | 37.69 | 7.07× |
| [8, 1000] Stirling | 6.93 | 13.31 | 1.92× |
| [0.01, 100] mixed positive | 11.07 | 55.64 | 5.03× |
| [−30, −0.01] reflection | 24.54 | 48.38 | 1.97× |

Best to worst spans a factor of 3.7, against one baseline on one machine.
That spread is the point of publishing per region at all: a single number
for "how much faster is lgamma" would hide it.

The reflection band costs most because the negative axis runs the positive
pipeline and then a reflection on top. Since #31 the reflection takes one dd
log, not two, and the [2.5, 8] band reads a Taylor table instead of walking
down to the zone (§13).

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
processes; the M1 showed none at 1.50× in thirteen, but seven bands step
1.07–1.30× between processes with the scalar walk stepping on its own
(§8.3).

---

## 6. What would change these numbers

Everything here is one point in a space with at least four axes, and three of
them are unexplored.

- **Another microarchitecture.** An earlier "wins from N lanes up" claim died
  on a Kaby Lake measurement. §8 is the second machine: a 4-lane AVX2 part
  from 2017, where the §4 batching gains run at roughly half these figures.
- **Another libm.** UCRT's `lgamma` is slow — 29–56 ns/el in three of its
  bands. Against Apple's the margins in §3 narrow and invert: Apple's x86
  `lgamma` is 12–20 ns/el in the same bands and corvus loses to it in every
  band on Kaby Lake (§8.1); on the M1 corvus wins only the table band (§11.2).
  This is the largest single source of error on this page.
- **Another compiler.** clang-cl only. MSVC cannot reach AVX-512 at all and
  would be measuring a different tier under the same name. Every fleet
  column is a clang-family build; §15.4 has the one GCC reading there is.
- **More runs.** §5 exists because two runs looked like agreement until a third
  disagreed.

A published performance claim needs at least the first two settled. **They
are** — four microarchitectures, four system libms (§11), three of each
first at v0.9.0 and again at v1.2.0, the fourth at v1.2.0 (§15). The compiler axis and §5 remain open;
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

## 8. Second machine: Kaby Lake AVX2, Apple libm

| | |
|---|---|
| machine | i7-7820HQ (Kaby Lake, 4 cores / 8 threads), macOS 13.7 |
| SIMD tier | `AVX2`, asserted per target via `CORVUS_EXPECT_TARGET` |
| compiler | Apple Clang 15, Release (`release` preset, `-G Ninja`) |
| libm | Apple libSystem (`std::lgamma` etc. as AppleClang links them) |
| harness | `tools/quiet_bench.sh` (the Unix counterpart of `quiet_bench.ps1`, same protocol) |
| code | v1.2.0 freeze kernels (`0a3ee74`, kernels as `2f863be`) |
| noise | two full passes gated at 4.04/4.00% and 2.03/2.29%; pass 1 ran 5–30% slower on some bands, so the tables quote pass 2. Evidence: `docs/bench-evidence/2026-11-01-kaby-v1.2.0/`. |

### 8.1 Against Apple's libm

| function | band | corvus ns/el | libm ns/el | ratio | Zen 4 / UCRT ratio (§3) |
|---|---|---:|---:|---:|---:|
| `erf` | n = 10⁶ | 6.69 | 24.34 | 3.64× | 6.28× |
| `erfc` | [−6, 6] core | 7.78 | 24.96 | 3.21× | 7.47× |
| `erfc` | [6, 28] tail | 29.95 | 25.54 | **0.85×** | 2.04× |
| `erfc` | [−6.5, 28] mixed | 41.00 | 26.11 | **0.64×** | 1.54× |
| `lgamma` | [0.5, 2.5] zone | 16.55 | 15.09 | **0.91×** | 5.06× |
| `lgamma` | [2.5, 8] table | 27.20 | 11.93 | **0.44×** | 7.07× |
| `lgamma` | [8, 1000] Stirling | 24.74 | 20.52 | **0.83×** | 1.92× |
| `lgamma` | [0.01, 100] mixed | 34.02 | 19.92 | **0.59×** | 5.03× |
| `lgamma` | [−30, −0.01] reflection | 79.49 | 46.86 | **0.59×** | 1.97× |

(n = 10⁶ rows.)

Two effects stack, and they pull in the same direction:

- **Apple's libm is fast.** Its `lgamma` is 12–21 ns/el outside the
  reflection band, 2–3× faster than UCRT's in the zone, table and mixed
  bands (Stirling is the exception: UCRT is faster there). corvus's own
  per-element cost moved by the lane count.
- **Four lanes, not eight, on a 2017 core.** corvus costs 2.9–4.0× more per
  element here than on Zen 4; the table band, at 5.1×, also pays for its
  gathers (§13.1).

So "lgamma is faster than the vendor libm in every band" is a **Zen 4 /
UCRT statement**. On Kaby Lake against Apple's libm, corvus `lgamma` is
slower in every band, and `erfc`'s double-double tail (#74) now loses there
too. `erf` and `erfc`'s core stay ahead of Apple's libm. Any published
positioning must name the libm it was measured against.

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

## 9. Retired

The v0.9.0 second Zen 4 pass (2026-08-30) lived here; §§1–4 and §11 now
carry v1.2.0 figures. Its text is in git history (tag `v1.1.0`), its
evidence in `docs/bench-evidence/2026-08-30-zen4-quiet/`.

---

## 10. Apple Silicon ambient gate — RATIFIED

**Status: RATIFIED 2026-08-31 by the maintainer.** The 5% gate stays the
universal standard; this section is the documented deviation for
heterogeneous-core Apple Silicon only. The M1 rows in §§8.2 and 11 are
published under it, subject to its own two-run agreement requirement. At
v1.2.0 that requirement is met on every band: passes 2 and 11 agree within
3.5% (§8.2). Getting there took 32 attempts over two nights and a day,
14 of them gated; `mediaanalysisd` and `corespotlightd` held the rest
above 10%.

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

### 10.1 Linux desktop: quit Electron apps, locking is not enough

Not a deviation: the Ivy Bridge machine (§15) meets the 5% gate. It meets
it only with every Electron app quit. With the Claude desktop app open, a
five-minute test read 9–14% busy unlocked and 6–7% with the screen locked
(#101): the app keeps rendering behind the Cinnamon lock screen and the
compositor keeps working. On 2026-11-11 the same app held the machine at
11–28% while in use; with it quit, both §15 passes gated on their first two
windows at 1–3% and sampled 0.00% before every bench. The Kaby Lake recipe
(screen locked) is therefore not enough on Linux/X11. Since #100 a failing
gate names the consumers on Linux as it does on macOS.

---

## 11. Fleet synthesis — v1.2.0 (2026-11-02; Ivy Bridge column 2026-11-11)

The quotable section. Sources: Zen 4 pass 2 (§1), Kaby Lake pass 2 (§8),
M1 pass 2 (§8.2; pass 11 agrees within 3.5%), Ivy Bridge pass 2 (§15). All
rows n = 10⁶, tier-asserted (`AVX3_ZEN4` / `AVX2` / `NEON` / `SSE4`), ratio
first, corvus ns/el in brackets. Ratios in one COLUMN share a libm; ratios across columns do
not, and §2's rule against comparing them still applies row by row.

### 11.1 Elementary family, four libms

| function | Zen 4 vs UCRT | Kaby vs Apple x86 | M1 vs Apple arm64 | Ivy Bridge vs glibc |
|---|---:|---:|---:|---:|
| `cos` | 8.44× (1.06) | 7.20× (3.05) | 5.21× (2.50) | 2.81× (8.46) |
| `sin` | 8.52× (1.06) | 6.68× (3.26) | 5.02× (2.61) | 2.79× (8.51) |
| `exp` | 1.85× (1.36) | 0.70× (8.32) | 2.59× (2.27) | 1.69× (7.01) |
| `expm1` | 1.07× (3.07) | 0.97× (11.07) | 0.45× (6.78) | 0.44× (24.75) |
| `log` | 1.50× (2.12) | 0.44× (10.28) | 0.66× (4.11) | 0.56× (13.86) |
| `log1p` | 0.91× (4.93) | 0.39× (18.81) | 0.32× (13.76) | 0.43× (33.92) |
| `log1pmx` † | 0.50× (11.36) | 0.26× (33.66) | 0.23× (27.11) | 0.19× (56.67) |

† Against `log1p(x) - x`, a cost reference only: it cancels near 0 and is
not an accurate log1pmx (#45).

Three positioning facts:

- **`cos`/`sin` win everywhere** — 5.0–8.5× against every vendor libm
  measured and 2.8× against glibc on two non-FMA lanes, full double range,
  1 ULP. The strongest performance row corvus has.
- **The log family trails the Apple libms and glibc (0.2–0.7×), and that is the
  product working as designed**: corvus `log` is correctly rounded on
  every reference row and `log1p`/`log1pmx` carry dd cores; the vendor
  libms trade exactly that away for speed. Against UCRT, `log` wins
  (1.50×) since #43.
- **`exp` is libm-dependent in both directions** (1.9× vs UCRT, 2.6× vs
  Apple arm64, 1.7× vs glibc, 0.70× vs Apple's x86 exp) — the cleanest demonstration that
  the baseline, not the kernel, moves these ratios.

### 11.2 lgamma and erf/erfc across the fleet

| function | band | Zen 4 vs UCRT | Kaby vs Apple x86 | M1 vs Apple arm64 | Ivy Bridge vs glibc |
|---|---|---:|---:|---:|---:|
| `erf` | n = 10⁶ | 6.28× (1.67) | 3.64× (6.69) | 6.52× (1.99) | 3.61× (7.54) |
| `erfc` | [−6, 6] core | 7.47× (2.41) | 3.21× (7.78) | 3.28× (3.91) | 2.49× (11.11) |
| `erfc` | [6, 28] tail | 2.04× (9.74) | 0.85× (29.95) | 0.38× (26.68) | 0.48× (79.25) |
| `erfc` | [−6.5, 28] mixed | 1.54× (12.68) | 0.64× (41.00) | 0.42× (27.80) | 0.39× (87.74) |
| `lgamma` | [0.5, 2.5] zone | 5.06× (5.67) | 0.91× (16.55) | 0.47× (25.11) | 0.32× (58.82) |
| `lgamma` | [2.5, 8] table | 7.07× (5.33) | 0.44× (27.20) | 1.39× (6.31) | 1.46× (17.97) |
| `lgamma` | [8, 1000] Stirling | 1.92× (6.93) | 0.83× (24.74) | 0.40× (21.84) | 0.41× (47.09) |
| `lgamma` | [0.01, 100] mixed | 5.03× (11.07) | 0.59× (34.02) | 0.40× (23.32) | 0.38× (53.10) |
| `lgamma` | [−30, −0.01] reflection | 1.97× (24.54) | 0.59× (79.49) | 0.59× (54.15) | 0.92× (138.70) |

`erf` and `erfc`'s core beat every libm measured. `lgamma` beats UCRT in
every band and loses to both Apple libms and to glibc in every band but
the table band on the M1 and on Ivy Bridge, where #31's Taylor table
(§13.2) passes Apple's arm64 `lgamma` and glibc's. One
function, one kernel, and the vendor baseline alone decides which side of
1.0 it lands on — this is why no README figure omits the libm name.
`erfc`'s tail and mixed bands lose to both Apple libms and to glibc since
#74 put a double-double G in the tail (§3).

The v0.9.0 synthesis, including the lgamma component study that scoped #31
(PERFORMANCE.md §11.3 at tag `v1.1.0`), is in git history.

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

| band (ns/el, n = 10^6) | A | B | A/B | Apple libm | B / libm | v0.9.0 |
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

**Against Apple's libm.** v0.9.0's statement that corvus `lgamma` is 3.6×
slower than Apple's in the recurrence band is superseded: 2.3× (0.44×)
after #31, and the reflection moves from 0.45× to 0.58×. Apple's lgamma
is still ahead in every band on this machine (§8.1).

**The Zen 4 0.97× readings.** The zone reads 0.99× (spread 2.5% / 5.3%)
and lbeta's main band 0.98× (3.1% / 5.1%): the same sign as Zen 4, both
within spread on both machines. The M1 read the same sign (§13.2), and
#96 carries the bisect recipe.

Controls: `gamma_p` R3, `gamma_q` R4, lbeta's big band and every
`beta_p_inv`/`beta_q_inv` band read 0.98–1.02×; `gamma_p_inv`/`gamma_q_inv`
read 1.03–1.11× by band, the same uplift as §12.1 sees, through lgamma.


### 13.2 Mac Mini M1 NEON (2026-10-31)

Evidence: `docs/bench-evidence/2026-10-31-m1-c3-ab/`. Apple's arm64 libm,
two lanes: the machine on which corvus lgamma was furthest behind at v0.9.0.

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

| band (ns/el, n = 10^6) | A | B | A/B | Apple libm | B / libm | v0.9.0 | Zen 4 A/B | Kaby A/B |
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
machine (v0.9.0's component study: the dd log 12.4 ns/el and the walk 5.0
atop the zone's floor). Per element, B runs the band within 15% of Zen 4's 5.35 on eight
lanes, where the gather-free bands run 2.3–4.3× slower than Zen 4 (zone
4.3×, Stirling 3.1×, reflection 2.3×). corvus lgamma is 1.39× faster
than Apple's arm64 lgamma on this band; it was 0.20× in v0.9.0.
Reflection moves from 0.39× to 0.57× of libm, mixed 0.34× → 0.41×.
Apple's is still ahead in every other band (§11.2).

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

---

## 15. Fourth machine: Ivy Bridge SSE4, glibc (2026-11-11)

| | |
|---|---|
| machine | MacBook Pro 9,1 (2012), i7-3820QM (Ivy Bridge, 4 cores / 8 threads): AVX, **no AVX2, no FMA** |
| OS | Linux Mint 22.3 (Ubuntu 24.04 base), kernel 7.0.0 |
| SIMD tier | `SSE4`, asserted per target via `CORVUS_EXPECT_TARGET`; the CPU's native top tier, uncapped (Highway has no AVX1 target) |
| compiler | Ubuntu clang 18.1.3, Release (`release` preset, `-G Ninja`, `CC=clang-18 CXX=clang++-18`) |
| libm | glibc 2.39 |
| harness | `tools/quiet_bench.sh` after #100 |
| code | `2565e85` (kernels as `2f863be`, the v1.2.0 freeze) |
| noise | two full passes gated at 1.00/3.00% and 3.00/1.00%, 0.00% before every bench; the gates tie and the tables quote pass 2 (below). Evidence: `docs/bench-evidence/2026-11-11-ivybridge-v1.2.0/`. |

What this machine adds: two lanes with no FMA on real silicon, so the
Dekker-split `TwoProd` (`src/ops-inl.h`) is timed at its real cost; the
first glibc baseline; and `SSE4` reached by the real dispatch decision on an
AVX-without-AVX2 CPU. Every other `SSE4` leg in the fleet is a capped build
on a CPU that has FMA. It is not a new tier and adds no accuracy claim.

**Which pass.** The n = 10⁶ rows and the fixed-n bands agree within 1%
between passes but for the rows named in the evidence README; the median
difference over all 242 timed values is 0.17%. Pass 1's disturbed rows are
`digamma`'s first two bands, which have no second n to read; pass 2's are
n = 10 000 rows of `exp`, `expm1` and the `lgamma` zone, each with an
n = 10⁶ row beside it. So the tables quote pass 2.

### 15.1 Against glibc

n = 10⁶ rows. The elementary family and the `erf`/`erfc`/`lgamma` rows are
in §11's Ivy Bridge column; this is the same data with the libm cost shown.

| function | band | corvus ns/el | glibc ns/el | ratio |
|---|---|---:|---:|---:|
| `erf` | n = 10⁶ | 7.54 | 27.20 | 3.61× |
| `erfc` | [−6, 6] core | 11.11 | 27.61 | 2.49× |
| `erfc` | [6, 28] tail | 79.25 | 38.13 | **0.48×** |
| `erfc` | [−6.5, 28] mixed | 87.74 | 34.15 | **0.39×** |
| `lgamma` | [0.5, 2.5] zone | 58.82 | 18.73 | **0.32×** |
| `lgamma` | [2.5, 8] table | 17.97 | 26.29 | 1.46× |
| `lgamma` | [8, 1000] Stirling | 47.09 | 19.15 | **0.41×** |
| `lgamma` | [0.01, 100] mixed | 53.10 | 19.94 | **0.38×** |
| `lgamma` | [−30, −0.01] reflection | 138.70 | 127.62 | **0.92×** |
| `cos` / `sin` | | 8.46 / 8.51 | 23.80 / 23.76 | 2.81× / 2.79× |
| `exp` / `expm1` | | 7.01 / 24.75 | 11.88 / 10.87 | 1.69× / **0.44×** |
| `log` / `log1p` / `log1pmx` † | | 13.86 / 33.92 / 56.67 | 7.74 / 14.69 / 11.02 | **0.56×** / **0.43×** / **0.19×** |

† As §11.1: against `log1p(x) - x`, a cost reference only.

glibc is fast where Apple's libms are: `lgamma` outside the table band at
19–20 ns/el, `log` at 7.7, `expm1` at 10.9. On two non-FMA lanes corvus
loses those rows and keeps `erf`, `erfc`'s core, `cos`/`sin`, `exp` and the
`lgamma` table band. §8.1's rule extends to glibc unchanged: any published
positioning names the libm.

### 15.2 Batching gain, corvus against corvus

Pass 2, beside the other 2-lane tier and the 4-lane x86 one (§8.2).

| family | upper bound (Ivy Bridge, SSE4) | M1 (NEON) | Kaby Lake (AVX2) |
|---|---|---|---|
| `gamma_p` / `gamma_q` | 1.76 – 1.98× | 1.78 – 2.00× | 3.13 – 3.92× |
| `beta_p` / `beta_q` | 1.75 – 2.00× | 1.73 – 2.00× | 3.21 – 4.02× |
| `digamma` | 1.63 – 2.13× | 1.70 – 2.10× | 2.88 – 4.77× |
| `trigamma` | 1.57 – 2.24× | 1.60 – 2.11× | 2.58 – 4.88× |
| `gamma_p_inv` / `gamma_q_inv` | 1.21 – 1.99× | 1.46 – 2.00× | 2.29 – 4.03× |
| `beta_p_inv` / `beta_q_inv` | 1.24 – 1.97× | 1.35 – 1.91× | 2.06 – 3.88× |
| `lbeta` | 1.76× main band, 2.01× big band | 1.66 – 1.70× main band, 2.02× big band | 3.16×, 3.97× |
| `erfinv` / `erfcinv` | 1.48 – 2.27× | 1.50 – 2.44× | 2.36 – 5.02× |
| `i0` / `i1` / `i0e` / `i1e` | 1.52 – 2.19× | 1.45 – 2.25× | 2.91 – 4.53× |

The column tracks the lane count and its range ends read within 7% of the
M1's on every family but the inverses, whose bottom ends are lower here:
`gamma_p_inv` a ≥ a_T deep tail 1.21×, `beta_p_inv` skewed (CF) 1.24×, both
branch-divergent bands. `beta_p_inv`'s top is pass 1's deep-small reading:
pass 2 read 2.16× there from a scalar walk 10% slower than pass 1's, with
the SIMD side at 2,937 and 2,940 ns/el.

### 15.3 Per-element cost against Kaby Lake

Both columns are pass 2 of a clang-family build at the same kernels, SIMD
ns/el. Half the lanes predicts 2×; this is also a 2012 core against a 2017
one, without FMA, so 2× is a floor for arithmetic-bound kernels, not an
expectation.

| row | Ivy Bridge SSE4 | Kaby Lake AVX2 | ratio |
|---|---:|---:|---:|
| `i0` series, n = 4001 | 166.41 | 40.76 | 4.08× |
| `trigamma` (0,1) up-step, n = 4001 | 108.88 | 29.06 | 3.75× |
| `digamma` (0,1) up-step, n = 4001 | 89.48 | 24.27 | 3.69× |
| `lgamma` zone [0.5, 2.5] | 58.82 | 16.55 | 3.55× |
| `erfinv` central [−0.5, 0.5] | 36.83 | 10.47 | 3.52× |
| `cos` / `sin` | 8.46 / 8.51 | 3.05 / 3.26 | 2.77× / 2.61× |
| `erfc` tail [6, 28] | 79.25 | 29.95 | 2.65× |
| `expm1` | 24.75 | 11.07 | 2.24× |
| `gamma_p_inv` / `gamma_q_inv` deep-small, n = 40 001 | 316.45 / 5,904.08 | 155.62 / 2,812.04 | 2.03× / 2.10× |
| `lgamma` Stirling [8, 1000] | 47.09 | 24.74 | 1.90× |
| `log1p` / `log1pmx` | 33.92 / 56.67 | 18.81 / 33.66 | 1.80× / 1.68× |
| `log` | 13.86 | 10.28 | 1.35× |
| `erf` | 7.54 | 6.69 | 1.13× |
| `exp` | 7.01 | 8.32 | **0.84×** |
| `lgamma` table [2.5, 8] | 17.97 | 27.20 | **0.66×** |

(n = 10⁶ unless the row says otherwise.)

Two readings, neither of them a decomposition:

- **The double-double bands sit well above 2×.** The `lgamma` zone runs dd
  near the zeros and `erfc`'s tail runs a dd G (#74); they cost 3.55× and
  2.65× Kaby Lake. Without FMA every `TwoProd` is a Dekker split, which is
  the cost this machine was added to expose, but no A/B isolates it here:
  core generation is in the same ratio. The three rows above them are §5's
  three bands; why they in particular cost most is not known.
- **The table kernels run at or below Kaby Lake's cost on half the lanes.**
  `SSE4` has no hardware gather, so Highway emulates it with per-lane
  loads: in effect the non-gather variant #22 proposes. `lgamma`'s table
  band (11 gathers per vector) is 1.5× faster per element here than on
  Kaby Lake's AVX2, `exp` 1.2× faster, and `erf` within 13%. That supports
  #22's premise that Kaby Lake's hardware gathers are its bottleneck. It is
  cross-machine evidence only; the sizing step is still the same-silicon
  A/B, Kaby Lake capped to `SSE4` against its AVX2 build.

§5 on this machine: two processes are too few to say. `digamma` (0,1)
up-step read 101.32 and 89.48 ns/el (1.13×) with its scalar walk moving
further (211.83 against 160.43), which is not §5's signature of a steady
walk; `trigamma` (0,1) up-step read 109.98 and 108.88, `i0` series 166.91
and 166.41. Recorded for #97.

### 15.4 Compiler side note: GCC 13.3 on the same silicon

Not a fleet column. One indicative GCC 13.3 pass on this machine (gate
3.00/0.00%, 2026-11-07, `c83e6f8`, same kernels) preceded this run; its raw
output is in #101. Against it, the clang build is faster on almost every
row, by 2–8% on most and by 15–18% on `gamma_p` R3 Temme, the
`gamma_p_inv` ridge and beyond-resolution bands and `trigamma`'s asymptotic
bands. `expm1` goes the other way: GCC reads 21.96 ns/el against clang's
24.75. One pass against two, so these are sizes, not rankings; they are the
only reading of §6's compiler axis so far.
