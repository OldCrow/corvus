# Kaby Lake re-gate of the one-element tail (2026-11-27)

`dev/v1.2.1` at `6940de9`, which contains `0dfe170` (a tail of one
element leaves through `GetLane`). Machine and toolchain as the parent
folder: i7-7820HQ, macOS 13.7.8, AppleClang 15.0.0, system Highway 1.4.0.
Release, `CORVUS_WERROR=ON`, each build from clean
(`cmake --build <dir> --clean-first`), each run asserted with
`CORVUS_EXPECT_TARGET=<tier>`.

| Tier | Build | ctest | `span1` | Warnings |
|---|---|---|---|---|
| AVX2 | native | 40/40 | PASS, 27 functions, lengths 1 to 9 | 0 |
| SSE4 | capped | 40/40 | PASS, 27 functions, lengths 1 to 9 | 0 |
| SSSE3 | capped | 40/40 | PASS, 27 functions, lengths 1 to 9 | 0 |
| SSE2 | capped | 40/40 | PASS, 27 functions, lengths 1 to 9 | 0 |

`masked_const` prints "not run: Windows only" and passes. Caps are the
ENVIRONMENT.md recipe. No bench, as PLAN.md directs: no cell K timed has
a one-element tail (n = 5 or 9 on four lanes).

Files: `ctest_*.txt`, `span1_*.txt` (each opens with the active target),
`build_*.log`, `head.txt` (commit, and the count of modified tracked
files at build time: 0).
