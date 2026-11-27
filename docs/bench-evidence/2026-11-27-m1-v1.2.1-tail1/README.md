# v1.2.1, M leg, second sitting: the one-element tail through GetLane (2026-11-27)

Follows `../2026-11-26-m1-v1.2.1/` (same machine, trees, protocol and
script, `quiet_run.sh.txt`; read that README first). That sitting found
the register-built tail costing erf and exp 1.5–2.3 ns a call at n = 3 and
n = 7 on NEON, where every odd length ends in a 1-element tail (its item
3). This one tests the fix proposed there: `StoreTail` with m = 1 writes
`GetLane(v)` and returns, so a trailing single element takes the same ops
as a span of one (`LoadTail` with m = 1 was already a broadcast). The
change is `driver-inl.h.diff`, eight lines, UNCOMMITTED on `dev/v1.2.1`
at `653758c` when this ran; `environment.txt` reads `dev: dirty=1` and the
library hash differs from the first sitting's (`1e8afa01…` against
`22b0a838…`). `main` is the same binary as the first sitting
(`3967d9cc…`), rebuilt by nobody. The branch moved to `02944a2` the same
night with K's and Z's evidence and a PERFORMANCE.md §6 edit, no code.

Gate first: the patched tree passes `CORVUS_EXPECT_TARGET=NEON ctest`
40/40, `span1` PASS on all 27 span functions (tails of 1 to 9 bit-identical
to a long span).

Protocol: 00:49–01:10 local (05:49–06:10 UTC), screen locked, user away.
The first gate took 15 min (two 5-min waits, `aborted/`: `mediaanalysisd`
between suspensions, the Vallum extension); the other five opened in
40–110 s. Gates 3.2–7.2%. **dev pass 2 ran at 12.3% ambient**
(`Bartender`, `Siri`, `mediaanalysisd` woke between the gate and the
bench); its long-table cells are the wide ones below. The other five runs
started at 3.2–6.8%.

## Results

`run_medians.txt` has every cell four ways: `main` on the 26th, `main` on
the 27th, dev1 (`653758c`), dev2 (`653758c` + this change); median of
three with min–max. Here, ns per CALL for a span of n values in place
(the table the change targets), `main` (27th) / dev1 / dev2:

| band | n=1 | n=2 | n=3 | n=4 | n=7 |
|---|---|---|---|---|---|
| erfinv [−0.5, 0.5] | 33.9 / 32.1 / 32.1 | 36.3 / 36.1 / 36.1 | 57.8 / 58.9 / **57.4** | 58.4 / 58.6 / 58.6 | 107.5 / 109.2 / **107.7** |
| erfcinv (0, 2) probit | 62.6 / 60.1 / 60.0 | 88.6 / 87.3 / 87.3 | 142.2 / 141.0 / **138.3** | 168.3 / 163.9 / 163.9 | 299.1 / 293.8 / **292.2** |
| erf [−3, 3] | 11.7 / 10.5 / 10.6 | 13.2 / 13.3 / 13.3 | 18.8 / 20.3 / **19.1** | 17.4 / 17.3 / 17.3 | 28.8 / 31.1 / **29.0** |
| exp [−20, 20] | 11.6 / 10.4 / 10.4 | 13.3 / 13.3 / 13.3 | 18.6 / 20.3 / **18.6** | 17.8 / 18.1 / 18.1 | 29.8 / 31.5 / **29.5** |

- **The tail cost is gone.** At n = 3 and n = 7, dev2 reads what `main`
  reads within 0.3 ns on erf and exp (18.6 against 18.6; 29.5 against
  29.8; 19.1 against 18.8; 29.0 against 28.8), and 0.4–3.9 ns under
  `main` on the inverses. Against dev1 that is −6 to −8% on erf/exp
  (0.92–0.94). Repeatable: the min of three equals the median in 19 of
  20 dev2 cells.
- **Nothing else moved.** n = 1 is unchanged from dev1 (the path is the
  same code), n = 2 and 4 are whole vectors and read the same in all
  three. `main` on the 27th reproduces `main` on the 26th to 0.1 ns in
  every short-table cell: same binary, same reading, one night apart.
- Long table (ns per value): dev2 equals dev1 within noise on every
  band's batch and array columns; the `local` column is 0.0–0.3 ns off on
  erf, exp and lgamma. The cells that read 1.04–1.13 (erfinv wide
  `local`, erfcinv central `local`, probit all three) have dev pass 2's
  12.3%-ambient reading as their max and one of the other two passes
  equal to dev1; not a difference between the trees. The `main` rows
  carry the same kind of outlier (erfcinv central `local` 78.0 once).

The first sitting's three findings stand; its item 3 is now answered:
the fix costs nothing anywhere measured and removes the 6–9%. Whether it
goes on the branch is the user's decision [OPEN]; if it does, K and I
validate the new head, and `test_span1` already guards it.
