# 2026-08-31 M1 full-surface run — INDICATIVE (10% gate), not a quiet pass

Mac Mini M1 (Apple M1, macOS Tahoe 26.6.2, AppleClang 21, Homebrew
Highway 1.4.0), `build-m1` Release, all 16 targets tier-asserted NEON
(`CORVUS_EXPECT_TARGET=NEON`), 16/16 rc=0. Gate armed at 9.02/9.50%
under `-m 10` per the ratified S2 fallback; per-target noise ran
13.8–50.5% (`mediaanalysisd` + `VTDecoderXPCService` resumed ~10 s
after gating, the same pattern as the 2026-08-23 runs). The
publish/withhold call is deferred to v0.9.0 S4.

Why the 10% gate: the 5% gate is structurally unreachable on this box
while the Photos first-run video-analysis backlog persists. This
session added two more failed windows on top of the ~20 h recorded
2026-08-28/29 (PLAN.md #24 entry):

- `gate_abort_2026-08-30_60min.log` — 341 samples, min 13.43%,
  `mediaanalysisd` in 327/341 failure snapshots.
- `gate_abort_2026-08-30_overnight_8h.log` — 2,724 samples, min 4.62%
  (a single sub-5% sample at 05:03, never two consecutive), median
  25.2%; `mediaanalysisd` in 94% of failure snapshots; `bztransmit`
  self-resumed overnight despite a manual pause (third night running),
  and `bzfilelist` scans even while transfers are paused.

`runner_log.txt` is the passing run (gate, per-target noise
annotations, rc lines). Per-target raw output:
`quiet_bench_bench_*.txt`.

Cross-check against the 2026-08-23 indicative runs (PLAN.md Resolved
log): stable rows agree — lgamma recurrence 0.20x (was 0.22x),
Stirling 0.39x (same), mixed 0.33–0.34x (same), erf 5.90x at n=1e6
(was 5.88–5.90x). The lgamma zone rows remain the noise-sensitive
outlier (0.35x/0.48x here vs 0.95x/0.61x then). First M1 elementary
rows: cos/sin ~3.2–3.3x streaming, exp 1.21–1.28x (inverts Kaby's
0.61x — Apple's fast x86 exp does not carry to M1), log 0.19–0.24x,
log1p ~0.31x (publish-honestly rows, same shape as Kaby).
