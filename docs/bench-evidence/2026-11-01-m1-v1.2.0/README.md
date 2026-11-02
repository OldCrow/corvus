# Row D, M1 leg at the v1.2.0 freeze kernels (2026-11-01/02)

Machine: Mac Mini M1, macOS 27.2, AppleClang 21.0.0, Highway 1.4.0 (system),
Release (`cmake --preset release -G Ninja -DCORVUS_WERROR=ON`), `NEON`
asserted via `CORVUS_EXPECT_TARGET`. Code: `head.txt` (`ac472c3`; kernels
as `2f863be`).

Protocol: `tools/quiet_bench.sh -b build -t NEON -m 10` (PERFORMANCE.md
§10), screen locked, Time Machine stopped, Spotlight indexing disabled on
both volumes, Backblaze and the other LaunchControl-managed agents
unloaded. `passN`: the default target list (15 benches). `repN`:
`bench_digamma bench_trigamma bench_erfinv bench_bessel`, fresh processes
for PERFORMANCE.md §5/§8.3. `status.txt`, `status-rep.txt`: every
attempt's exit code, 32 attempts in all. A pass that did not see two
consecutive sub-10% windows within its wait (30, 60 or 45 min by attempt)
aborted with exit 3; its gate log is kept under `aborted/<pass>-<HHMM
started>/`. Night of 2026-11-01 (three drivers, 30 then 60 min waits):
passes 1 and 2, reps 1, 3 and 4 gated. Day run 2026-11-02 07:18–16:19
(`day-run.log`; one driver alternating full and short attempts, 45 min
waits, stopping on one full pass and three short processes with noise
under 10% before every bench): passes 5, 7, 8, 9, 11 and 12 and reps 9, 11
and 12 gated; the three reps met the per-bench rule, no full pass did
(`mediaanalysisd` or `corespotlightd` above 10% before 5 to 15 of the
benches); pass 13 was stopped by hand when the user returned. Gated: 8
full passes, 6 short processes.

Gates: pass 1 7.17/7.02%, pass 2 4.03/4.64%, pass 5 5.70/2.49%, pass 7
2.32/2.70%, pass 8 3.85/5.25%, pass 9 2.48/2.81%, pass 11 2.40/3.14%,
pass 12 3.37/2.26%; rep 1 3.74/3.83%, rep 3 8.71/4.86%, rep 4 6.60/6.50%,
rep 9 5.76/3.64%, rep 11 4.13/1.80%, rep 12 2.24/2.07%. Per-bench noise
before each bench is in every `quiet_bench.log`; `status*.txt` carries
each day-run pass's maximum.

Noise notes. Pass 1 gated at 20:38 and from `bench_digamma` on sampled
50–53% with `corespotlightd` at 280–340% (three to four cores, so P cores
too) and `mediaanalysisd` at 40–180%; its `digamma` and `trigamma` bands
read 4–33% slower than passes 2, rep 1 and rep 4, SIMD and scalar walk
together (up-step 40.99/84.71 against 30.76/64.61 ns/el), and its
`gamma_p_inv`/`gamma_q_inv` deep-tail bands moved the same way. The
upper-bound ratios survive that (both sides move), which is why pass 1
still agrees with pass 2 within 6% on every band but the eight named in
PERFORMANCE.md §8.2. Pass 2 sampled 15–19% before five benches
(`mediaanalysisd` at one core, the E-core class §10 describes) and 3–8%
before the rest. Rep 1 sampled 11–24% (`mediaanalysisd`, VideoToolbox);
rep 4 17–31% (`photoanalysisd` and `mediaanalysisd` at a core each).
Rep 3, gated on its third attempt, sampled 3–9% before every bench. Rep 2
aborted three times (minimum windows 17.4%, 14.6% and 8.38%, the last a
single sub-10% window with no second). Day run: pass 11 sampled under 5.5%
before nine benches and 13–15% (`mediaanalysisd` at a core) before the
other six, and agrees with pass 2 within 3.5% on every band; passes 5, 7
and 9 sampled under 10% before ten benches each; passes 8 and 12 ran with
`mediaanalysisd` at one to two cores throughout (13–22%, pass 12 32–50%
before `digamma`, `trigamma` and `gammainv` with `corespotlightd`). Reps 9,
11 and 12 sampled under 4.2% before every bench.

The set holds fourteen gated processes with the four short benches, which
is what PERFORMANCE.md §8.3's mode table counts (thirteen, pass 1 apart).

`sum.py` (copied from the Zen 4 folder) compares pass 1 and pass 2 per
band; `python3 sum.py` with the pass numbers edited compares any two.
Published: PERFORMANCE.md §8.2 (M1 column) and §8.3.
