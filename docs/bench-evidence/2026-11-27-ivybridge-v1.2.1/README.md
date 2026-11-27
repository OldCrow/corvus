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
