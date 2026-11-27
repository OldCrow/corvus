# Why the n = 7 sign differed by kernel (Kaby Lake, 2026-11-27)

Machine and toolchain as the parent folder. i7-7820HQ, family 6 model 158
stepping 9, microcode revision 248 (0xF8). All figures are ns, the lowest
reading over the fresh processes of a run, with the user's background
processes stopped (SIGSTOP) for the run and ambient CPU between 2% and 4%
unless a line says otherwise. Each claim below says whether an on/off
experiment shows it or it is inferred.

## Answer

1. **The sign does not differ by kernel.** With the caller and the layout
   controlled, the register-assembled tail is 2 to 8 ns slower than the
   masked tail at n = 7 on all four kernels, and 6 to 12 ns faster at
   n = 2 and 3. [shown]
2. **`bench_span1` showed `erf` and `erfcinv` the other way because of
   where the linker put an 8-byte constant.** On this CPU a load that
   touches the 16 bytes at offset 0x7c0 to 0x7cf of its 4 KB page, and
   follows an AVX2 gather, is delayed. `main`'s and the masked-tail
   variant's executables had such a constant on the `erf` and inverse-erf
   paths; the branch's did not. Moving the constants 64 bytes turns the
   10 to 30 ns difference off, and on in the branch's library. [shown]
3. **So the committed leg-3 table overstates the gain on the gather
   kernels.** With the constant out of the window, `main` -> branch on the
   one-value `local` shape is 0.84x to 1.01x, not 0.57x to 0.94x.
   `matched_layout.txt`; table below. [shown]

## The tail by itself

What the two tails are in the AVX2 code (disassembly of `ErfinvImpl`):
- `main`'s `LoadNOr` is not a masked load on AVX2. Highway emits a branch
  on the count: for three elements a 16-byte load plus an 8-byte load,
  for two a 16-byte load. `StoreN` is a real `vmaskmovpd`.
- The branch's `LoadTail` is a loop of compare, broadcast and blend per
  element. `StoreTail`'s copy loop is compiled to a `memcpy` call.
- Every outlined kernel call passes its 256-bit argument and result
  through stack memory (the translation unit's baseline is SSE2, so the
  System V classification of a 256-bit vector is MEMORY).

Four libraries built from one source tree, differing only in
`DriveUnary`'s tail: branch (`v0`), masked load with `StoreTail` (`v1`),
`LoadTail` with masked store (`v2`), `main`'s tail (`v3`), `LoadTail` with
a store loop that cannot become `memcpy` (`v4`). `bench_tail2.cpp` writes
the inputs as n separate 8-byte stores, calls, and reads the results;
`tail_store_width.txt` has every cell at eight alignments.

| write, call, read (8-byte stores) | n = 2 | n = 3 | n = 4 | n = 7 | n = 8 |
|---|---|---|---|---|---|
| erfinv central, masked tail (`v3`) | 81.2 | 82.9 | 81.2 | 124.3 | 121.2 |
| erfinv central, branch (`v0`) | 69.2 | 73.9 | 81.2 | 128.2 | 121.4 |
| erfinv central, `LoadTail` + masked store (`v2`) | 66.3 | 70.4 | 80.9 | 127.9 | 121.1 |
| erf, masked tail | 49.3 | 51.5 | 49.2 | 77.7 | 74.2 |
| erf, branch | 39.3 | 42.8 | 49.1 | 80.0 | 74.3 |
| exp, masked tail | 51.5 | 52.9 | 48.1 | 83.1 | 77.9 |
| exp, branch | 42.1 | 46.6 | 48.5 | 85.6 | 78.6 |

- The gain at n = 2 and 3 is on the load side: `v2` has it, `v1` does
  not. The 16-byte load of two values the caller wrote separately cannot
  be forwarded; the per-element loads can. [shown]
- At n = 7 the full vector's 32-byte load has the same problem in every
  variant, so the tail's forwarding buys nothing and the extra work of
  the blend loop and the `memcpy` call shows. With nothing written or read
  around the call (`tail_shapes.txt`, shape `none`) the register tail
  costs `erfinv` 8 to 16 ns at n = 2, 3 and 7. [shown]
- n = 4 and n = 8 are about 30 ns above the same call with nothing written
  beforehand (81 against 49, 121 against 91): the full-vector load still
  waits for the caller's stores. v1.2.1 does not address that. [shown]
- Array alignment, the result array's position (`bench_tail3.cpp`) and
  16 to 48 bytes of code padding change none of this. [shown]

## The gather window

`gather_alias.cpp` and `gather_alias2.cpp` contain no corvus code: one
`vgatherqpd` (or `vgatherdpd`) from a small table, then one 8-byte load,
both feeding the next iteration.

| Probe | Median | Load touching page offset 0x7c0..0x7cf |
|---|---|---|
| 64-bit-index gather, then load | 13.67 | 18.2 |
| 32-bit-index gather, then load | 13.17 | 18.7 |
| plain 32-byte load, then load | 1.03 | 1.03 |

- The slow start offsets are 0x7b9 to 0x7cf at one-byte resolution: every
  8-byte load that overlaps those 16 bytes. The window is the same in
  both pages swept and for five positions of the gathered table,
  including one that crosses a page. [shown; `gather_alias2.txt`]
- It behaves as a load held behind a store 4096*k bytes away (4K
  aliasing), with the gather standing in for the store. What the gather
  does internally is not known, nor whether this is the Downfall (GDS)
  microcode or the instruction as designed. [inferred]

In corvus the loads are the 8-byte literals of the polynomial that follows
the table lookup. `window_check.py` lists, for an executable, the constant
loads after the first gather on the one-value path that fall in the
window.
- `literal_sweep.txt`: the constants moved through all 64 positions of a
  page in 64-byte steps, `erf` called with one fixed input. Three
  positions cost 44 to 47 ns against 33.6. The rule predicts all 64.
  [shown]
- At a slow position the cost is the same at 256 stack positions, with
  20 sizes of environment, for every input from 0 to 100, with the table
  moved 64 bytes, and with the code at the same addresses as a fast
  build. Only moving the literals changes it. [shown; these fixed-input
  probes ran without the ambient gate, and their levels repeat to
  0.1 ns]
- Calls chained through their results slow by the same 10.5 ns
  (41.5 -> 52.5) as independent ones (33.6 -> 44.2): the delay is on the
  kernel's own path. [shown]
- In `bench_span1` itself (`const_shift_span1.txt`): the masked-tail
  library at its built position reads `erf` local 41.9 and `erfcinv`
  n = 4 239.7; 64 bytes away, 31.9 and 209.7. The branch's library 64
  bytes the other way reads 41.9 and 239.5. [shown]

One 8-byte literal is shared by many kernels, because the linker merges
equal literals. In `main`'s `bench_span1` it sits at 0x7c8 and 19 loads
in 10 AVX2 functions read it.

## Leg 3 with the layout matched

`bench_span1` against `main`'s library as built, against the same library
with its constants moved 1024 bytes (no load in the window), and against
the branch's as built (none on these paths).

| `local`, one value | `main` as built | `main`, clean | branch | clean ratio |
|---|---|---|---|---|
| erfinv [-0.5, 0.5] | 54.9 | 50.7 | 51.2 | 1.01x |
| erfinv [-0.999, 0.999] | 121.2 | 108.1 | 104.7 | 0.97x |
| erfcinv [0.5, 1.5] | 53.1 | 53.1 | 48.1 | 0.91x |
| erfcinv (0, 2) (probit) | 123.5 | 110.6 | 103.7 | 0.94x |
| erf [-3, 3] | 47.8 | 38.0 | 31.9 | 0.84x |
| lgamma [1, 21] | 109.7 | 110.0 | 104.6 | 0.95x |
| exp [-20, 20] | 40.9 | 41.2 | 35.2 | 0.85x |

| per call, clean `main` -> branch | n = 1 | n = 2 | n = 3 | n = 4 | n = 7 |
|---|---|---|---|---|---|
| erfinv [-0.5, 0.5] | 66.7 -> 59.9 | 82.9 -> 73.1 | 84.5 -> 75.6 | 81.9 -> 81.6 | 111.7 -> 118.4 |
| erfcinv (0, 2) (probit) | 119.3 -> 108.8 | 179.2 -> 167.8 | 204.1 -> 194.1 | 210.6 -> 209.5 | 375.9 -> 374.1 |
| erf [-3, 3] | 40.1 -> 35.7 | 51.5 -> 40.7 | 53.3 -> 44.3 | 49.4 -> 49.3 | 71.3 -> 75.8 |
| exp [-20, 20] | 43.6 -> 39.3 | 53.2 -> 43.4 | 55.1 -> 47.9 | 49.4 -> 49.1 | 78.9 -> 81.3 |

macOS on AVX2 still gains on every one-value band but `erfinv` central,
which is level. Batch throughput moves with the layout too (`erfinv`
wide: 54.5 as built, 45.8 clean).

## Ruled out, each by an experiment

- Jitter: 15 and 12 alternating pairs (`alternate.log`, `three_way.log`).
- The build tree: `bench_span1` built in the tree and outside it agree to
  0.2 ns on every cell.
- Code position: 16 to 2048 bytes of padding after the library's code.
- Stack position: `span1_slide.cpp` at 14 positions, `bench_fixed.cpp` at
  256, and the environment at 20 sizes.
- The table's position, and the input.
- The branch predictor learning the 512 repeated calls: 64 to 65,536
  distinct calls per repetition move `erfcinv` by 3 to 4%
  (`bench_tail5.cpp`; machine not gated).
- A constant that straddles a page or cache line: none does.

## Also found

- **`bench_span1`'s caller does not write 8-byte stores.** Its "write n
  values" loop compiles to a `memcpy` call on AppleClang 15, so the
  driver's loads meet whatever stores `memmove` chooses for n*8 bytes.
  With real 8-byte stores n = 7 and n = 8 cost `erfinv` 10 ns more
  (`tail_store_width.txt`, `s8` against `mc`).
- **A gather can also alias a real stack store.** With one fixed input,
  the table row whose address is 4096*k bytes from a stack slot the call
  has just written costs 4 to 7 ns more; the slow rows are exactly 4096
  bytes apart and move with the table (`bench_rowmap.cpp`). Over random
  inputs it averages out, but it gives bands of stack positions where
  `erf` at n = 4 costs up to 30% more (`stack_slide.txt`).
- **The first cells a process measures read 5 to 12% high.** A cell
  measured first is not comparable with one measured later.

## Not known

- Whether other Intel parts, other microcode revisions, AMD or any
  AVX-512 gather have a window, and where. `gather_alias2.cpp` answers it
  in a minute on any AVX2 machine.
- Whether this is behind #97's speed modes. On macOS the page offset of a
  constant is fixed per executable, so it varies between builds, not
  between processes; #97's modes vary between processes.
- Which figures in docs/PERFORMANCE.md's Kaby Lake columns were taken with
  a constant in the window.
- A remedy. The literal pool's position is the linker's choice. Loading
  constants from a table whose alignment corvus controls, or a
  non-gather variant (#22), would remove the exposure; neither is
  designed.

## Files

- Benches (each links `libcorvus.a` only; build line in
  `bench_tail.cpp`): `bench_tail.cpp` (caller shapes), `bench_tail2.cpp`
  (store width, alignment), `bench_tail3.cpp` (result position),
  `bench_tail5.cpp` (distinct calls), `bench_slide.cpp`,
  `span1_slide.cpp`, `bench_fixed.cpp` (stack position),
  `bench_rowmap.cpp` (cost per input), `span1_local.cpp`.
- Standalone: `gather_alias.cpp`, `gather_alias2.cpp` and their outputs.
- `exp_tail.diff`: the first masked-tail variant. The four-variant macro
  is the same edit with `CORVUS_TAILV` selecting the load and the store.
- `run_tail.sh` (the stop/resume runner), `agg.py`, `agg2.py`,
  `aggslide.py`, `window_check.py`.
- Tables: `tail_shapes.txt`, `tail_store_width.txt`, `stack_slide.txt`,
  `literal_sweep.txt`, `const_shift_span1.txt`, `matched_layout.txt`.
- Raw per-process outputs are not kept here.
