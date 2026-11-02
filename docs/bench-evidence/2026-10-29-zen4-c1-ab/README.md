# Zen 4 quiet A/B for row C1 (#42 lever 1), 2026-10-29

Six `tools/quiet_bench.ps1` passes on `AVX3_ZEN4` (clang-cl Release),
alternating A (benches built from `3712cc3`, before `6027e2d` and
`0967994`) and B (built from `64a085b`). Each pass directory holds the
gate log and the raw per-bench output. `summarise.py` prints the median-of-three table that
docs/PERFORMANCE.md §12 quotes.
