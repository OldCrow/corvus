# Row D, Kaby Lake leg at the v1.2.0 freeze kernels (2026-11-01)

Machine: i7-7820HQ (Kaby Lake), macOS 13, AppleClang 15.0.0, Release
(`cmake --preset release -G Ninja -DCORVUS_WERROR=ON`), `AVX2` asserted via
`CORVUS_EXPECT_TARGET`. Code: `head.txt` (`0a3ee74`; kernels as `2f863be`).

Protocol: `tools/quiet_bench.sh -b build -t AVX2 -m 5 -w 20`, screen locked,
non-security services stopped. `pass1`, `pass2`: the default target list
(15 benches). `rep1`–`rep4`: `bench_digamma bench_trigamma bench_erfinv
bench_bessel`, four more fresh processes for PERFORMANCE.md §5/§8.3.
`status.txt`, `status-rep.txt`: exit codes. Every `quiet_bench.log` carries
the gate samples and the per-bench noise; all six gated below 5%.

Noise notes: pass 2 sampled 17% (`biomesyncd`, `mds`) before `bench_exp` and
`bench_log`; rep 1 sampled 10.7% (`contactsd`) before `bench_bessel`, and its
Bessel scalar walk moved with it. Pass 1 (gate 4.0/4.0%) ran 5–30% slower
than pass 2 on `digamma` and `trigamma`; the repeats agree with pass 2.

`sum.py` (copied from the Zen 4 folder) compares pass 1 and pass 2 per band.
Published: PERFORMANCE.md §8.2 (Kaby column) and §8.3.
