# v1.2.1, M leg: bench_span1 on the M1, main against dev/v1.2.1 (2026-11-26)

Machine: Mac Mini M1 (Macmini9,1; 4 performance + 4 efficiency cores), macOS
27.2.1, AppleClang 21.0.0, system Highway 1.4.0, `NEON` asserted
(`CORVUS_EXPECT_TARGET`). Why: PLAN.md, "v1.2.1 fleet legs", step 3: whether
macOS gains anything from the branch. Both trees `cmake --preset release -G
Ninja -DCORVUS_WERROR=ON`; `environment.txt` has the hashes and sizes (they
differ, OldCrow/standards PERFORMANCE-MEASUREMENT.md rule 2).

- `main`: `c1d6640` in a worktree (`../corvus-main`) with the branch's
  `tests/bench_span1.cpp` copied in and registered; not committed.
- dev: `dev/v1.2.1` at `653758c`, clean.
- Gate first (step 1): `CORVUS_EXPECT_TARGET=NEON ctest` 40/40 on the dev
  tree, `span1` PASS on all 27 span functions, `masked_const` "not run:
  Windows only". Step 2 (capped) has no tier below NEON on arm64; not run,
  as for v1.2.0.

## Protocol

`quiet_run.sh.txt`, launched under `caffeinate`, 22:56–23:31 local
(03:56–04:31 UTC 2026-11-27). Three passes; each pass runs `main` then dev
(dev first on pass 2). Each run is `tools/quiet_bench.sh -b <tree>/build -t
NEON -m 8 -w 5 bench_span1`: ambient CPU (100 − `top` idle) below 8% over
two consecutive 10 s windows, then the bench in a fresh process. Before each
run fourteen of the user's processes (indexing, media analysis, sync; the
list in the script and `environment.txt`, eight of them present) were
suspended with SIGSTOP and resumed with SIGCONT after each pass — the
libstats 2026-11-24 method, approved for corvus by the user 2026-11-26.
The user was at the machine for the first half hour.

- `run.txt`: the sitting log. The first gate took 29 minutes (four 5-minute
  waits, `aborted/`; `WindowServer`, the Vallum endpoint-security
  extension, `mediaanalysisd` between suspensions, and Spotlight's `mds` /
  `mdworker_shared`, which run as `_spotlight` and are out of a user
  signal's reach, indexing the two fresh build trees). The other five gates
  opened in 20–90 s. Every run started at 3.2–3.8% ambient.
- `<tree>_passN/quiet_bench_bench_span1.txt`: that run's output.
  `<tree>_passN/quiet_bench.log` is the tree's cumulative log as
  `quiet_bench.sh` appends it, so each copy repeats the earlier runs, and
  the dev tree's copy begins with the v1.2.0 row-D attempts (15:51 lines,
  10% gate) made in the same build directory on 2026-11-01. Lines from
  22:56 on 2026-11-26 are this sitting.

## Results (median of three; min–max in `run_medians.txt`)

ns per value, 4096 values. `batch`: one span call. `array`: 4096 calls of
length 1 over arrays. `local`: value in, result out, as a scalar wrapper
does it. `main` → dev (dev/`main`).

| band | batch | array | local |
|---|---|---|---|
| erfinv [−0.5, 0.5] | 15.5 → 12.4 (0.80)† | 30.2 → 28.6 (0.95) | 33.6 → 29.4 (0.87)‡ |
| erfinv [−0.999, 0.999] | 37.3 → 42.7 (1.14)† | 60.6 → 55.9 (0.92)‡ | 60.1 → 58.7 (0.98)‡ |
| erfcinv [0.5, 1.5] | 12.6 → 11.1 (0.88) | 29.0 → 27.5 (0.95) | 30.0 → 28.1 (0.94) |
| erfcinv (0, 2) probit | 38.0 → 37.0 (0.97) | 60.7 → 54.6 (0.90) | 59.5 → 56.0 (0.94) |
| erf [−3, 3] | 2.0 → 2.0 (1.00) | 8.9 → 7.6 (0.85) | 9.7 → 8.4 (0.87) |
| lgamma [1, 21] | 25.2 → 25.2 (1.00) | 39.7 → 38.9 (0.98) | 41.5 → 40.3 (0.97) |
| exp [−20, 20] | 2.3 → 2.3 (1.00) | 9.1 → 7.8 (0.86) | 9.5 → 8.3 (0.87) |

ns per CALL, a span of n values written into a local array just before
the call and read back after it (`LocalN`). NEON has 2 lanes, so n = 2 and
4 are whole vectors, n = 3 and 7 end in a 1-element tail.

| band | n=1 | n=2 | n=3 | n=4 | n=7 |
|---|---|---|---|---|---|
| erfinv [−0.5, 0.5] | 33.9 → 32.1 (0.95) | 36.3 → 36.1 | 57.8 → 58.9 (1.02) | 58.4 → 58.6 | 107.5 → 109.2 (1.02) |
| erfcinv (0, 2) probit | 62.4 → 60.1 (0.96) | 88.6 → 87.3 | 142.2 → 141.0 | 168.2 → 163.9 (0.97) | 299.1 → 293.8 (0.98) |
| erf [−3, 3] | 11.6 → 10.5 (0.91) | 13.2 → 13.3 | 18.8 → 20.3 (1.08) | 17.4 → 17.3 | 28.8 → 31.1 (1.08) |
| exp [−20, 20] | 11.6 → 10.4 (0.90) | 13.3 → 13.3 | 18.6 → 20.3 (1.09) | 17.8 → 18.1 | 29.8 → 31.5 (1.06) |

Spread. Cells without a mark repeat to within 0.1–0.6 ns over the three
passes in both trees. † erfinv central batch reads 11.5 or 15.5 in both
trees (`main` 11.5/15.5/15.5, dev 12.4/11.6/15.5): a two-mode band
(PERFORMANCE.md §5), not a difference between the trees. erfinv wide batch
read 37.0/47.3/42.7 on dev against 37.3/37.6/37.3 on `main`; the erfinv
kernel is unchanged on the branch, so this is read as a mode too, but three
passes do not settle it. ‡ `main` pass 1, the first run of the sitting,
read 72.5 and 95.0 in the two erfinv `local` cells and 83.5 in the wide
`array` cell, against 29.7–33.6, 58.0–60.1 and 58.1–60.6 in passes 2 and
3; the medians stand on the other two passes.

## What it says

1. **The serialisation defect is not on this platform.** On `main`, the
   `local` column equals the `array` column within 2–8% on every band
   (erf 9.7 against 8.9, lgamma 41.5 against 39.7, probit 59.5 against
   60.7). On Zen 4 under Windows `local` read about twice `array`. The
   AAPCS64 convention passes a span in two registers, and NEON's masked
   tail did not serialise the calls either.
2. **macOS gains a little, from the length-1 path and the erfcinv branch.**
   A span of one costs 1.1–1.2 ns less per call on erf and exp (10%), 1.8
   ns on erfinv central (5%) and 2.3 ns on the probit band (4%). In the
   `local` shape the cheap kernels read 13% less, the inverses 3–13%. The
   erfcinv central batch, where the driver does nothing different, reads
   12% less (12.6 → 11.1 ns per value): that is change 3, the special-value
   selects behind a branch, and the only kernel change these bands touch.
3. **The register-built 1-element tail costs about 1.5–2.3 ns more per
   call on the cheap kernels** (erf and exp at n = 3 and n = 7: +6–9%,
   repeatable to 0.1 ns), and is within noise on the inverses. A 2-lane
   target meets this tail on every odd length. The old tail was one masked
   `StoreN`; the new one stores the vector to a stack buffer and copies one
   element out. A tail of one element could take the length-1 path
   (broadcast, `GetLane`) instead. Tried the same night:
   `../2026-11-27-m1-v1.2.1-tail1/` — the cost is gone, nothing else
   moves; whether it goes on the branch is the user's decision [OPEN].

Indicative rather than published: six fresh processes, three per tree, on
one evening; the §5 mode bands above did not settle. No number here goes
into PERFORMANCE.md's tables without the user's say.

Caveat added 2026-11-27 (K's finding, `2026-11-26-kaby-v1.2.1/n7/`):
AppleClang compiles `bench_span1`'s write-n-values loop to a `memcpy`
call, so the short-span table for n >= 2 here, as on Kaby Lake, did not
write the inputs as separate 8-byte stores. The n = 1 column and the long
table are unaffected. Ivy Bridge (clang 18) is the clean reading of the
short-span shape, and its pass 3 confirmed the one-element-tail finding
and the fix with the same diff.
