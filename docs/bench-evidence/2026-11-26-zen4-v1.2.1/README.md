# Zen 4, v1.2.1: one value through the span API, quiet before and after

2026-11-26, Asus TUF A16 (Ryzen 7 7445HS, Windows 11), active target
`AVX3_ZEN4`, clang-cl 22.1.3, Release. Power plan Turbo, screen locked.

- **Before:** `main` at `ad5c9d5` (v1.2.0 code), built in its own checkout.
- **After:** `dev/v1.2.1` at `653758c`.

Both benches link the library only. Three rounds, the two trees
alternating. Every measured process was launched directly and opted out of
Windows power throttling, which a locked screen otherwise applies (about
1.5x slower). Before each launch a gate waited for total CPU under 8% in
two consecutive 5-second samples and a fixed integer loop timed the clock.
`run.log` has the gate and probe of every launch: gates 2.1 to 7.6%,
probes 664 to 710 ms. `quiet.ps1.txt` is the script.

## One value, value in and result out (`r*_span1.txt`, column `local`)

ns per call, median of three rounds.

| | before | after |
|---|---:|---:|
| `erfcinv` [0.5, 1.5] | 132.0 | 29.4 |
| `erfcinv` uniform on (0, 2) | 165.8 | 67.3 |
| `erfinv` [-0.5, 0.5] | 64.7 | 29.3 |
| `erfinv` [-0.999, 0.999] | 97.9 | 67.8 |
| `erf` [-3, 3] | 32.9 | 15.7 |
| `exp` [-20, 20] | 35.8 | 14.4 |
| `lgamma` [1, 21] | 83.5 | 56.5 |

## A short span built in place (`r*_span1.txt`, second table)

`erfinv` [-0.5, 0.5], ns per call.

| | n=1 | n=2 | n=3 | n=4 | n=7 |
|---|---:|---:|---:|---:|---:|
| before | 65.4 | 65.8 | 66.4 | 67.8 | 69.9 |
| after | 30.5 | 34.2 | 35.0 | 35.9 | 39.3 |

## Batch of 4096 (`r*_span1.txt`, column `batch`)

ns per element. `erfinv` [-0.5, 0.5] reads 3.4 to 3.8 after; before it
read 3.8 to 6.4, the widest spread in the run. `erfcinv` [0.5, 1.5]:
14.7 before, 3.4 after. `erfcinv` uniform on (0, 2): 26.8 before, 15.6
after.

## Image pages unmapped against mapped (`r*_cold.txt`)

41 bands of the 27 span functions; a band is flagged when it is slower
unmapped in both of two passes.

| | before | after |
|---|---|---|
| bands flagged, each of three rounds | 3 | 0 |
| `erfcinv` [0.5, 1.5] | 3.97x to 4.09x | 0.99 |
| `erfcinv` 2^[-26, -1] | 1.86x to 1.90x | none |
| `log1p` [-0.5, 2] | 3.16x to 3.31x | 1.00 |

After, no band in any round reads more than 10% from 1.00 in either pass.

## Not covered

One machine and one compiler. The warm-up outputs are not kept. Each bench
takes a median of 41 passes (`bench_span1`) or a minimum of 30
(`bench_cold`); the rounds repeat within about 1 ns.
