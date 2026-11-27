# v1.2.1 fleet leg, Ivy Bridge (2026-11-26/27)

Machine: MacBook Pro 9,1 (2012), i7-3820QM (Ivy Bridge, 4 cores / 8
threads; AVX, no AVX2, no FMA), Linux Mint 22.3, kernel 7.0.0, glibc 2.39,
Ubuntu clang 18.1.3, Highway 1.4.0 (fetched). Release, `CORVUS_WERROR=ON`
(`CC=clang-18 CXX=clang++-18 cmake --preset release -G Ninja -B
build-clang -DCORVUS_WERROR=ON`). Branch `dev/v1.2.1` at `653758c`; `main`
at `c1d6640`. `SSE4` is this CPU's native top tier, so the native build is
uncapped. `environment.txt` has the rest.

## Legs 1 and 2: gates (`gates/`)

| Tier | Build | ctest | `span1` | Warnings |
|---|---|---|---|---|
| SSE4 | native | 40/40 | PASS, 27 functions | 0 |
| SSSE3 | capped | 40/40 | PASS, 27 functions | 0 |
| SSE2 | capped | 40/40 | PASS, 27 functions | 0 |

Every run was asserted with `CORVUS_EXPECT_TARGET=<tier>`. `masked_const`
prints "not run: Windows only" and passes. Caps are the ENVIRONMENT.md
recipe, each from an empty `build-cap`. The `span1_*.txt` files open with
the active target; there is none for SSSE3 (its build directory was
reused for SSE2 before the output was saved), so SSSE3's `span1` is the
passing line in `ctest_SSSE3_capped.txt`.

## Leg 3: `bench_span1`, `main` against the branch

Two build trees: this checkout and a worktree of `main` at
`/home/wolfman/Development/corvus-main`, with `tests/bench_span1.cpp`
copied in and its two CMake lines appended, built with the same command.
`libcorvus.a` is 8,558,858 bytes on `main` and 8,706,130 on the branch;
each pass's `head.txt` has both trees' commits and hashes.

Both passes: `tools/quiet_bench.sh -b build-clang -t SSE4 -w 30`, `main`
then branch, each a fresh process, from a plain terminal, AC power,
governor `schedutil`. Redaction: 13 `top:` entries that named user
applications read `[user-app]` (each at 1.0% but one at 3.9%); nothing
else in the logs is edited. One of those applications was still running
during both passes; ambient read 0.00% before every bench.

### Pass 1 (`pass1/`): the plan's protocol, three alternations

`bench_span1` as the plan gives it. All six runs gated on their first two
windows (1.00/0.00%, 2.00/0.00%, the rest 0.00/0.00%).

Two of its readings are not usable, and pass 2 exists because of them:
- The first band the program times, `erfinv` central, read slow on both
  trees: batch 50.8 to 57.2 ns per value against 37 in pass 2, and the
  branch's `local` 88.5, 104.9 and 103.9. The kernel behind that batch
  figure is unchanged by the branch.
- `erf` on the branch read at two levels across processes: `local` 24.4
  in the first and 21.2 and 21.3 after; n = 7 at 99.1, then 77.2 and
  78.1.

Outside those two bands and the batch column, every cell agrees with
pass 2 within 3%, median against median.

### Pass 2 (`pass2/`): pinned and warmed, eight alternations

Two changes from pass 1 and no others: each bench starts behind a wrapper
(`bench_span1_warm.branch.txt`) that pins to CPU 2, spins that CPU for
2 s, then execs the bench on it; and eight alternations, not three. Each
output opens with the core's frequency before the exec (3.69 GHz in all
16). Fourteen runs gated at 0.00/0.00%, one at 4.00/0.00%, one at
0.00/1.00%.

With the core at speed the first band reads the same on both trees (batch
37.1 and 37.0), so pass 1's slow first row was the bench starting on a
core still coming up from idle. Demonstrated by turning it off, on this
machine and governor only.

The two `erf` levels are not a property of the branch. In pass 2 the
branch never showed the higher one in eight processes, and `main` showed
it in two (`local` 28.0 once; n = 7 at 95.8 and 84.8 against 76.2). `exp`
at n = 3 and n = 7 does the same on both trees (`main` 61.6 to 61.8 in
three of eight processes against 50.5 to 50.8; the branch 61.3 in one).
Cause not known; PERFORMANCE.md §5 and #97 describe the same behaviour.

Median of eight, ns per value, `main` -> branch (`summary.txt`, from
`sum.py`):

| Band | batch | array | local |
|---|---|---|---|
| erfinv [-0.5, 0.5] | 37.1 -> 37.0 | 80.6 -> 80.3 | 82.8 -> 81.4 (0.98x) |
| erfinv [-0.999, 0.999] | 80.7 -> 80.6 | 118.3 -> 117.2 | 120.6 -> 118.6 (0.98x) |
| erfcinv [0.5, 1.5] | 38.6 -> 37.2 | 83.6 -> 79.5 | 85.3 -> 81.0 (0.95x) |
| erfcinv (0, 2) (probit) | 82.2 -> 80.9 | 121.3 -> 117.1 | 123.8 -> 118.8 (0.96x) |
| erf [-3, 3] | 7.5 -> 7.6 | 21.2 -> 19.8 | 24.0 -> 21.2 (0.89x) |
| lgamma [1, 21] | 60.5 -> 60.5 | 88.7 -> 85.9 | 90.5 -> 88.0 (0.97x) |
| exp [-20, 20] | 7.4 -> 7.5 | 20.5 -> 19.2 | 22.2 -> 21.1 (0.95x) |

ns per call, span in place:

| Band | n = 1 | n = 2 | n = 3 | n = 4 | n = 7 |
|---|---|---|---|---|---|
| erfinv [-0.5, 0.5] | 90.5 -> 86.0 | 91.3 -> 90.4 | 172.4 -> 173.9 | 165.2 -> 164.7 | 312.7 -> 312.8 |
| erfcinv (0, 2) (probit) | 132.5 -> 123.1 | 189.0 -> 184.8 | 305.2 -> 302.8 | 352.9 -> 348.1 | 630.5 -> 621.5 |
| erf [-3, 3] | 29.8 -> 28.2 | 30.4 -> 30.1 | 52.0 -> 53.7 | 43.6 -> 42.8 | 76.2 -> 77.3 |
| exp [-20, 20] | 26.4 -> 23.5 | 28.1 -> 27.0 | 50.8 -> 51.6 | 42.1 -> 41.2 | 78.7 -> 77.6 |

Outside the excursions named above, most cells repeat within 2% from
process to process; single readings sit up to 6% out, and one 12%
(`exp` at n = 2 on the branch, 30.2 against 27.0) (`summary.txt`, the
max/min table). The n = 1 and `local` gains have the same sign in all
eight pairs.

What it says:
- Linux on SSE4 gains a little on a one-value call and close to nothing
  on a longer span. The `local` shape is 0.89x (`erf`) to 0.98x of
  `main`; n = 1 calls are 0.89x to 0.95x, 1.6 to 9.4 ns per call.
- From n = 2 the ratios are 0.96x to 1.03x. On two lanes n = 2 and n = 4
  have no tail. The cells with one (n = 3, n = 7) read 0.99x to 1.03x.
- Three cells are slower, by 1.1 to 1.7 ns per call. Two are slower in
  all eight pairs: `erfinv` central at n = 3 (172.4 -> 173.9) and `erf`
  at n = 3 (52.0 -> 53.7). `erf` at n = 7 (76.2 -> 77.3) is slower in the
  six pairs where `main` sat on its lower level. Kaby Lake's slower
  cell, `erfinv` central at n = 7, is level here (312.7 -> 312.8).
- Batch throughput is unchanged, but for `erfcinv` central (38.6 ->
  37.2) and probit (82.2 -> 80.9).
- The mechanism is not shown here. The System V convention passes a span
  in registers, so the Windows 16-byte copy does not apply; which of the
  five changes produces the gain on this machine is not separated.

Not run: `bench_cold` (Windows only).

## The slower cells (`pass3/`, 2026-11-27)

The three slower cells are the one-element tail, and M's patch removes
the cost. On a 2-lane tier every odd length ends in a tail of one
element, which the branch's `StoreTail` sends through a stack buffer.
Shown by two on/off builds beside `main` and the branch, all four clang
18.1.3 Release from their own trees:

- **masked**: `653758c` with `DriveUnary`'s tail put back to `main`'s
  `LoadNOr`/`StoreN` (K's `n7/exp_tail.diff`; `pass3/exp_masked.diff`).
- **lane0**: `653758c` with M's saved diff, `StoreTail` with m = 1 writes
  `GetLane(v)` (`2026-11-27-m1-v1.2.1-tail1/driver-inl.h.diff`;
  `pass3/exp_lane0.diff`). Gated here first: 40/40 under
  `CORVUS_EXPECT_TARGET=SSE4`, `span1` PASS on all 27, 0 warnings
  (`pass3/ctest_lane0_SSE4.txt`).

Protocol as pass 2 (pinned to CPU 2, warmed 2 s, `quiet_bench.sh -t SSE4`),
the four trees alternating, eight fresh processes each. 31 runs gated at
0.00/0.00% and one at 3.00/0.00%; the core read 3.67 to 3.69 GHz before
every bench. `main`'s library is the binary of passes 1 and 2 (its
worktree's HEAD moved to `c826564`, PLAN.md only; `pass3/head.txt` has
the hashes). 20 `top:` entries are redacted as above.

Median of eight, ns per call, `main` / branch / masked / lane0
(`summary_pass3.txt`, from `sum4.py`, with min and max):

| Band | n = 1 | n = 2 | n = 3 | n = 4 | n = 7 |
|---|---|---|---|---|---|
| erfinv [-0.5, 0.5] | 90.4 / 85.9 / 86.0 / 86.0 | 91.2 / 90.3 / 90.4 / 90.4 | 172.4 / 173.9 / 172.4 / **172.2** | 165.1 / 164.6 / 164.9 / 164.7 | 312.7 / 312.8 / 311.2 / **312.0** |
| erfcinv (0, 2) (probit) | 132.6 / 122.7 / 122.9 / 122.8 | 188.9 / 184.2 / 184.2 / 184.3 | 305.3 / 302.8 / 300.4 / **301.1** | 352.8 / 347.6 / 347.6 / 347.5 | 630.7 / 621.4 / 620.0 / **621.6** |
| erf [-3, 3] | 29.7 / 28.2 / 28.1 / 28.1 | 30.4 / 30.0 / 30.1 / 30.1 | 52.0 / 53.8 / 51.8 / **52.2** | 43.5 / 42.8 / 42.8 / 42.9 | 76.2 / 77.3 / 74.7 / **75.5** |
| exp [-20, 20] | 26.1 / 23.8 / 23.3 / 23.6 | 27.6 / 27.8 / 26.8 / 27.0 | 50.6 / 51.7 / 50.3 / **50.0** | 42.0 / 41.6 / 41.0 / 41.2 | 76.0 / 78.2 / 74.4 / **74.9** |

- **The cost is gone with M's patch.** At n = 3 and n = 7 the patched
  tree reads at or under `main` on every band but `erf` at n = 3, which
  is 0.2 ns over (52.2 against 52.0), and 1.6 to 3.3 ns under the branch
  on `erfinv` n = 3, `erf` and `exp`. The branch reproduces pass 2's
  slower cells to 0.1 ns (173.9, 53.8, 77.3).
- **Nothing else moved.** n = 1, 2 and 4 read the same on the branch and
  both experiment trees within 0.6 ns, bar `exp` at n = 2 (27.8 against
  26.8 and 27.0, inside that cell's spread). The long table's three
  columns agree across the three within 0.6 ns on every band.
- **The masked tail is as fast here, or faster.** The masked tree reads
  0.5 to 0.8 ns under the patched one on `erf` and `exp` at n = 7, and
  within 0.4 ns of it at n = 3. So on SSE4 the masked one-element tail
  is no slower than either register form; the patched form gives up at
  most 0.8 ns to it.
- This is the 2-lane case only. It says nothing on K's n = 7 cell, a
  three-element tail after a full vector on four lanes.
- Excursions as in pass 2: single processes read high on single cells in
  every tree (the branch's first band once at 98.7 and 99.1; `main`
  `erf` n = 7 once at 91.3; `exp` n = 7 at 86.3 and 89.0), so the
  medians are quoted and the ranges are in `summary_pass3.txt`.
