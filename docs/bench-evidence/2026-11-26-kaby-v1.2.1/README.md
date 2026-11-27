# v1.2.1 fleet leg, Kaby Lake (2026-11-26)

Machine: MacBookPro14,3, i7-7820HQ (Kaby Lake), macOS 13.7.8, AppleClang
15.0.0, system Highway 1.4.0. Release, `CORVUS_WERROR=ON`. Branch
`dev/v1.2.1` at `653758c`; `main` at `c1d6640`. `environment.txt` has the
rest.

## Legs 1 and 2: gates (`gates/`)

| Tier | Build | ctest | `span1` | Warnings |
|---|---|---|---|---|
| AVX2 | native (`cmake --preset release -G Ninja -DCORVUS_WERROR=ON`) | 40/40 | PASS, 27 functions | 0 |
| SSE4 | capped | 40/40 | PASS, 27 functions | 0 |
| SSSE3 | capped | 40/40 | PASS, 27 functions | 0 |
| SSE2 | capped | 40/40 | PASS, 27 functions | 0 |

Every run was asserted with `CORVUS_EXPECT_TARGET=<tier>`; each
`span1_*.txt` opens with the active target. `masked_const` prints "not
run: Windows only" and passes. Caps are the ENVIRONMENT.md recipe
(`build-cap`).

## Leg 3: `bench_span1`, `main` against the branch

Two build trees: this checkout and a detached worktree of `origin/main`
at `/Users/wolfman/Development/corvus-main`, with `tests/bench_span1.cpp`
copied in and its two CMake lines appended. `libcorvus.a` is 5,065,216
bytes on `main` and 5,166,472 on the branch.

Protocol (`run_span1.sh`): `tools/quiet_bench.sh -b build -t AVX2 -m 5 -s
10 -w 5 bench_span1`, `main` then branch, three times, each a fresh
process. The user's background processes named in `environment.txt` were
stopped (SIGSTOP) before each gate and resumed (SIGCONT) after each run
[user, 2026-11-26]. Screen locked, AC power. All six runs gated on the
first try, between 2.66% and 4.46%. `run.log` has the order and times.

Noise: `main2` sampled 8.10% just before the bench (`signpost_reporter`);
`main3`'s process line shows `biomesyncd` at 96% with ambient at 3.82%.
The three `main` processes spread more than the branch's: 1.51x on batch
`erfinv` central (16.3 in the first process, 11.2 and 10.8 after), up to
1.21x on the short-span cells (`exp` n = 7: 78.8, 88.4, 95.7) and under
1.12x elsewhere. The three branch processes agree within 6% on every
cell. So a ratio within about 10% of 1 is not established by these runs.

Median of three, ns, `main` -> branch (`summary.txt`, from `sum.py`):

| Band | batch | array | local |
|---|---|---|---|
| erfinv [-0.5, 0.5] | 11.2 -> 11.2 | 56.4 -> 52.5 | 54.3 -> 51.2 (0.94x) |
| erfinv [-0.999, 0.999] | 53.0 -> 48.3 | 124.8 -> 103.6 | 131.4 -> 105.0 (0.80x) |
| erfcinv [0.5, 1.5] | 11.2 -> 10.2 | 55.9 -> 46.7 | 53.6 -> 48.1 (0.90x) |
| erfcinv (0, 2) (probit) | 50.9 -> 44.6 | 121.3 -> 102.2 | 127.7 -> 103.6 (0.81x) |
| erf [-3, 3] | 6.7 -> 6.4 | 44.5 -> 30.8 | 55.5 -> 31.9 (0.57x) |
| lgamma [1, 21] | 49.2 -> 49.2 | 110.6 -> 103.1 | 114.5 -> 104.7 (0.91x) |
| exp [-20, 20] | 7.7 -> 7.7 | 37.6 -> 34.1 | 48.7 -> 37.1 (0.76x) |

ns per call, span in place:

| Band | n = 1 | n = 2 | n = 3 | n = 4 | n = 7 |
|---|---|---|---|---|---|
| erfinv [-0.5, 0.5] | 70.2 -> 63.4 | 83.0 -> 73.1 | 89.9 -> 75.7 | 83.6 -> 81.6 | 111.7 -> 118.4 |
| erfcinv (0, 2) (probit) | 132.0 -> 108.9 | 203.2 -> 167.8 | 237.9 -> 194.0 | 243.4 -> 209.4 | 422.2 -> 374.1 |
| erf [-3, 3] | 50.1 -> 35.7 | 61.1 -> 40.7 | 69.0 -> 44.3 | 51.6 -> 49.4 | 87.9 -> 75.7 |
| exp [-20, 20] | 43.6 -> 39.3 | 53.2 -> 43.4 | 59.8 -> 47.9 | 54.1 -> 49.1 | 88.4 -> 81.3 |

What it says:
- macOS on AVX2 gains. The `local` shape is 0.57x to 0.94x of `main`; the
  short-span calls are 0.64x to 0.98x, with one exception.
- The exception: `erfinv` central at n = 7 reads 118.4 on all three
  branch processes against 111.7, 111.7 and 121.1 on `main`, 6% slower at
  the median. It is real; see "The n = 7 cell" below.
- Batch throughput is unchanged or better on every band.
- The mechanism is not shown here. The System V convention passes a span
  in registers, so the Windows 16-byte copy does not apply; which of the
  five changes produces the gain on this machine is not separated.

Not run: `bench_cold` (Windows only).

## The n = 7 cell (`n7/`, 2026-11-27)

`erfinv` central at n = 7 is 6% slower on the branch, about 6.5 ns per
call, and the register-assembled tail (change 5) is what costs it. Shown
by turning it off: a third tree, the branch at `653758c` with only
`DriveUnary`'s tail put back to `main`'s `LoadNOr`/`StoreN`
(`n7/exp_tail.diff`; `test_span1` passes on it).

`n7/three_way.log`: `main`, branch and that tree, alternating, 12 fresh
processes each, background processes stopped, ambient 3.30% before and
3.62% after. Every binary reads one of two levels 5.5% apart from process
to process (PERFORMANCE.md's speed modes; the masked-tail tree sat on the
higher one in 7 of 12 processes), so the table gives each binary's lower
level, and the n = 7 / n = 4 ratio inside one process, which does not move
with the level.

| ns per call | n = 2 | n = 3 | n = 4 | n = 7 | n = 7 / n = 4 |
|---|---|---|---|---|---|
| erfinv, `main` | 82.7 | 84.0 | 81.6 | 111.7 | 1.37 |
| erfinv, branch | 73.1 | 75.4 | 81.6 | 118.2 | 1.45 |
| erfinv, branch with the masked tail | 82.6 | 83.9 | 81.6 | 112.4 | 1.38 |
| exp, `main` | 52.9 | 54.6 | 49.2 | 78.8 | 1.60 |
| exp, branch | 43.4 | 47.7 | 49.1 | 81.3 | 1.66 |
| exp, branch with the masked tail | 53.0 | 54.5 | 48.8 | 78.6 | 1.61 |
| erf, `main` | 60.7 | 62.8 | 49.2 | 80.9 | 1.64 |
| erf, branch | 40.6 | 44.2 | 49.4 | 75.7 | 1.53 |
| erf, branch with the masked tail | 61.2 | 62.6 | 49.4 | 80.9 | 1.64 |
| erfcinv probit, `main` | 202.6 | 231.9 | 240.7 | 421.0 | 1.75 |
| erfcinv probit, branch | 167.8 | 194.0 | 209.4 | 373.8 | 1.79 |
| erfcinv probit, branch with the masked tail | 201.2 | 230.4 | 239.6 | 415.7 | 1.73 |

- n = 4 has no tail and reads the same in all three trees.
- A tail by itself (n = 2, 3) is 7 to 36 ns faster assembled in
  registers, on all four bands.
- A tail after a full vector (n = 7) is 6.5 ns slower for `erfinv`
  central and 2.5 ns slower for `exp`, 5 ns faster for `erf`, and 42 to
  47 ns faster for `erfcinv` probit.
- Why the sign differs by kernel is not known. No change is proposed
  here: sending a tail that follows a full vector back to the masked form
  would give up the `erf` and `erfcinv` gains at n = 7.

`n7/alternate.log` is an earlier two-way run on a machine in use
(ambient 11.4% and 7.5%), indicative only; it shows the same sign in 14
of 15 pairs.
