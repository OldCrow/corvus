# Mac Mini M1 quiet A/B for row C3 (#31), 2026-10-31

Six `tools/quiet_bench.sh` passes attempted on `NEON` (Mac Mini M1, macOS
27.2, AppleClang 21.0.0, Highway 1.4.0, Release), alternating A (benches
built from `1fe08b1`, before #31) and B (built from `4f78527`, the
`dev/v1.2.0` head, kernels and benches identical to Kaby's `0708fed`).
Gated at 10% ambient (PERFORMANCE.md §10), 30 min per pass. Four passes
ran: 2 (B), 3 (A), 4 (B) and 5 (A), gated at 2.5–9.5%; per-bench noise
2.9–5.8% before bench_lgamma, then 7–45% before the later benches (pass
2 7–25%, passes 3–5 15–45%: `mediaanalysisd` at one core, `corespotlightd`
and `spotlightknowledged.updater` bursts, §10's E-core class of load).
Passes 1 and 6 aborted (no two consecutive sub-10% windows in 30 min,
`mediaanalysisd`); their folders hold the gate logs only. `status.txt`
records every exit code; `hashes.txt` the A and B binary hashes.

`extra-A-1828/` is a gated A pass (3.5 / 5.2%, exit 0) from an earlier,
abandoned run of the same driver on the same binaries; `summarise.py`
skips it (not `pass*`), and it serves as a third A reading for the §10
agreement check: within 2% of the A medians on the zone, mixed and
reflection bands, within 3% on every beta, gammainv and betainv band
(the flagged gammainv band apart); 6–7% above on the table and Stirling
bands and gamma R1.

Medians of two per variant. `python3 ../2026-10-30-zen4-c3-ab/summarise.py .`
prints the table that docs/PERFORMANCE.md §13.2 quotes. Spread at most
2.9% on the lgamma rows except reflection (B 6.9%: 54.2 / 58.1 ns/el).
Flagged, not published: the lbeta main band, bimodal regardless of
variant (A 293.8 / 302.9 and 347.9 in the extra pass; B 292.4 / 347.2);
the gammainv deep-small `p_inv` band (B 5.4% spread; 46% in the C1 A/B).
