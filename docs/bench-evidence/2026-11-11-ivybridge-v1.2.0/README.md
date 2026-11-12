# Fourth machine, Ivy Bridge SSE4 / glibc, at v1.2.0 (2026-11-11)

Machine: MacBook Pro 9,1 (2012), i7-3820QM (Ivy Bridge, 4 cores / 8 threads;
AVX, no AVX2, no FMA), Linux Mint 22.3 (Ubuntu 24.04 base, kernel 7.0.0),
glibc 2.39, Ubuntu clang 18.1.3, Highway 1.4.0 (fetched), Release
(`CC=clang-18 CXX=clang++-18 cmake --preset release -G Ninja -B build-clang`),
0 warnings, 38/38 ctest under `CORVUS_EXPECT_TARGET=SSE4`. `SSE4` is this
CPU's native top tier (Highway has no AVX1 target), so the build is uncapped.
Code: `head.txt` (`2565e85`; kernels as `2f863be`).

Protocol: `./tools/quiet_bench.sh -b build-clang -t SSE4 -w 30`, run from a
plain terminal with the Claude desktop app quit (PERFORMANCE.md §10.1: an
Electron app holds this machine above the gate even behind the lock screen).
`pass1`, `pass2`: the default target list (15 benches). `status.txt`: exit
codes. Each `quiet_bench.log` carries the gate samples and the per-bench
noise: both passes gated on their first two windows (1.00/3.00% and
3.00/1.00%) and sampled 0.00% before every bench. These are the first logs
written by the script after #100, so the machine line and the `top:` fields
are populated. Redaction: four `top:` entries that named user applications
read `[user-app]` (each at 1.0%); kernel threads, system daemons and the
desktop's own processes are as logged. Nothing else in the logs is edited.

Noise notes: the gates tie. The n = 10⁶ rows and the fixed-n = 4001 bands
agree within 1% between passes but for those named here; the median
difference over all 242 timed values is 0.17%. Pass 1 ran slow on `digamma`
(0,1) up-step and [1,2) zone (SIMD 101.32/84.15 against 89.48/71.06 ns/el,
the scalar walk moving with it) and on `log`/`log1p`/`log1pmx` at n = 1000.
Pass 2 ran slow on `exp`/`expm1` at n = 10 000 (1.17×/1.30×), on the `lgamma`
zone at n = 10 000 (1.28×) and on `erf` libm at n = 10 000 (1.25×); its
`beta_p_inv` deep-small and gamma-limit scalar walks read 10% and 7.5% above
pass 1 with the SIMD side within 0.1%. `erf` at n = 10⁶ read 8.05 (pass 1)
and 7.54 (pass 2). `erfinv` central at n = 10 000 read 48.3 ns/el in both
passes against 36.8 at n = 10⁶; it is the first band that bench times.

Validation leg (ACCURACY.md, the v1.2.0 validation list): `ulp_report.sh
<build-dir>` runs the program list of CI's "ULP report" step, in its order,
printing a `### <program> <reference files>` line before each.
`ulp-report-clang.txt` is the clang 18.1.3 build above and
`ulp-report-gcc.txt` a GCC 13.3 Release build of the same tree (0 warnings,
38/38 ctest), both under `CORVUS_EXPECT_TARGET=SSE4`: 18 of 18 programs
pass in each and the two files are byte-identical. Against the AVX2 report
of CI run 43129858616 (same commit) the 639 result lines keep their order
and row counts; 149 of 153 max-ULP readings are equal, none passes its
gate, and 109 of 121 real-error pin values are equal, the rest moving as
ACCURACY.md's no-FMA note describes (`erf`: 37 not correctly rounded and
1.10 real, its published SSE4 figures). No capped-`SSE4` report is
published to diff against, so identity with a capped build on an FMA host
is expected, not shown.

`sum.py` (copied from the Kaby Lake folder) compares pass 1 and pass 2 per
band. Published: PERFORMANCE.md §15 and the Ivy Bridge column of §11, from
pass 2. The indicative GCC 13.3 pass that preceded this run is in #101, not
here: it is a compiler side note, not a fleet column.
