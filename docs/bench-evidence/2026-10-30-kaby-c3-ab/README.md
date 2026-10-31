# Kaby Lake quiet A/B for row C3 (#31), 2026-10-30

Six `tools/quiet_bench.sh` passes on `AVX2` (MacBook Pro 2017, i7-7820HQ,
macOS 13, AppleClang 15.0.0, Release), alternating A (benches built from
`1fe08b1`, before #31) and B (built from `0708fed`, the `dev/v1.2.0` head).
Gated at 3.2–4.8% ambient; pass 1 (A) saw 12–15% noise before its beta,
lbeta and gammainv benches, which widens those rows' A spread. Each pass
directory holds the gate log and the raw per-bench output; `hashes.txt`
has the A and B binary hashes; `summarise.py` prints the median-of-three
table that docs/PERFORMANCE.md §13.1 quotes.
