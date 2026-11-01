# Mac Mini M1 quiet A/B for row C1 (#42 lever 1), 2026-10-31

Six `tools/quiet_bench.sh` passes attempted on `NEON` (Mac Mini M1, macOS
27.2, AppleClang 21.0.0, Highway 1.4.0, Release), alternating A (benches
built from `3712cc3`, before `6027e2d` and `0967994`) and B (built from
`4f78527`, the `dev/v1.2.0` head, kernels and benches identical to Kaby's
`0708fed`; B also carries rows C2 and C3). Gated at 10% ambient
(PERFORMANCE.md §10), 30 min per pass. Three passes ran: 1 (A), 4 (B) and
5 (A), gated at 3.1–6.7%, per-bench noise 2.7–4.5% (pass 1), 3.2–9.2%
(pass 4) and 11.7–15.7% (pass 5, `mediaanalysisd` at one core
throughout). Passes 2, 3 and 6 aborted: the machine never held two
consecutive sub-10% windows in 30 min (`mediaanalysisd` 50–180% with
VTDecoder/VTEncoder, `corespotlightd` bursts); their folders hold the
gate logs only and carry no bench output. `status.txt` records every
exit code; `hashes.txt` the A and B binary hashes.

Medians: A over two runs, B over one here. B's four binaries are the
same files as the C3 A/B's B (`hashes.txt` in both folders), which ran
gated twice more in `2026-10-31-m1-c3-ab` (passes 2 and 4), so the §10
two-run check for B uses those: the beta rows agree within 0.5%, the
betainv rows within 1% except the two cheapest `beta_p_inv` bands
(deep-small closed form, R1-tiny deep target), which read 5–7% apart.
`python3 ../2026-10-30-zen4-c3-ab/summarise.py .` prints the table that
docs/PERFORMANCE.md §12.2 quotes. Flagged, not published: the gammainv
deep-small `p_inv` band (184 / 278 ns/el over the two A runs, 266 for B;
the bench's shortest timing).
