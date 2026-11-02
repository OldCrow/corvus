# Kaby Lake quiet A/B for row C1 (#42 lever 1), 2026-10-30

Six `tools/quiet_bench.sh` passes on `AVX2` (MacBook Pro 2017, i7-7820HQ,
macOS 13, AppleClang 15.0.0, Release), alternating A (benches built from
`3712cc3`, before `6027e2d` and `0967994`) and B (built from `0708fed`, the
`dev/v1.2.0` head, which also carries rows C2 and C3). Gated at 3.2–3.6%
ambient. Each pass directory holds the gate log and the raw per-bench
output; `hashes.txt` has the A and B binary hashes; `summarise.py` prints
the median-of-three table that docs/PERFORMANCE.md §12.1 quotes.
