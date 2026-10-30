#!/usr/bin/env python3
"""Generate tests/data/expm1_reference.txt -- correctly rounded expm1 oracle.

Each line: <input-hex-double> <expm1(input)-hex-double>, output rounded once
to nearest double by refgen_common.round_resolving_ties (a value within
1e-12 ulp of a midpoint is re-evaluated at 400 digits, or the row is
dropped). Specials (+-0 -> +-0 signed, +-inf, NaN) are covered by the smoke
test, not here; every input is a finite nonzero double.

Point selection (fixed seed, reproducible), by the kernel's own seams
(src/exp-inl.h Expm1Vec):
  - series band |x| <= 2^-6, both signs: log-spaced down to the subnormal
    inputs (expm1(x) = x there), plus full-mantissa draws per binade
    2^-60..2^-7 (power-of-two anchors have near-empty mantissas, which hide
    rounding in a series -- #35 H1)
  - the 2^-6 series/exponential cut: +-128-ulp bit ladders at +-2^-6 and
    dense draws in [2^-7, 2^-4] both signs, where the exponential branch's
    amplified error peaks
  - the exponential branch: uniform on (-38, 0) and log-uniform up to
    709.78
  - the -38 shortcut: a ladder around -38 and draws on [-45, -30], where
    expm1 meets -1 (the last value above -1 sits near x = -37.43)
  - the overflow boundary: +-160-ulp ladder around the largest x with a
    finite result (exp's; expm1 differs from exp by 1, far below an ulp
    there), so rows above it are exact +inf

Usage:
    python tools/gen_expm1_reference.py > tests/data/expm1_reference.txt
"""

import math
import random
import struct
import sys

import mpmath as mp

from refgen_common import round_resolving_ties

DPS = 60
mp.mp.dps = DPS

SEED = 20261029


def as_bits(x: float) -> int:
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def from_bits(u: int) -> float:
    return struct.unpack("<d", struct.pack("<Q", u))[0]


def ref_value(x):
    return round_resolving_ties(
        mp.expm1(mp.mpf(x)),
        lambda _dps, x=x: mp.expm1(mp.mpf(x)))


def emit(points):
    seen = set()
    rows = []
    for x in points:
        if not (math.isfinite(x) and x != 0.0):
            continue
        b = as_bits(x)
        if b in seen:
            continue
        seen.add(b)
        y = ref_value(x)
        if y is None:  # unresolved near-tie (#89): decline the row
            continue
        rows.append((x, y))
    return rows


def self_check(rows):
    """Cut rows (|x| within a factor 2 of 2^-6), overflow-ladder rows and
    every 97th row overall, recomputed at 2x dps bit-identical. Per-arm
    caps, never a truncating global cap (#35 M2/L3)."""
    cut = [i for i, (x, _) in enumerate(rows)
           if 2.0 ** -7 <= abs(x) <= 2.0 ** -5][::7][:600]
    over = [i for i, (x, _) in enumerate(rows) if x > 709.0]
    every_97th = list(range(0, len(rows), 97))
    sample = sorted(set(cut + over + every_97th))
    bad = []
    with mp.workdps(2 * DPS):
        for i in sample:
            x, y = rows[i]
            y2 = ref_value(x)
            if y2 is None:
                y2 = float("nan")
            if as_bits(y2) != as_bits(y):
                bad.append((x, y, y2))
    if bad:
        for x, y, y2 in bad:
            print(f"self-check MISMATCH: x={x.hex()} stored={y.hex()} "
                  f"recomputed={y2.hex()}", file=sys.stderr)
        return False
    n_inf = sum(1 for _, y in rows if math.isinf(y))
    if n_inf == 0:
        print("self-check: no +inf overflow rows", file=sys.stderr)
        return False
    print(f"self-check: {len(sample)} rows ({len(cut)} cut, {len(over)} "
          f"overflow ladder, {n_inf} +inf) re-verified at dps={2 * DPS}: OK",
          file=sys.stderr)
    return True


def ladder(anchor, k):
    b = as_bits(anchor)
    return [from_bits(b + j) for j in range(-k, k + 1)]


def main():
    rng = random.Random(SEED)
    pts = []

    # Series band, both signs, down into the subnormal inputs.
    for _ in range(3072):
        x = 10.0 ** rng.uniform(-320, math.log10(2.0 ** -6))
        pts.append(rng.choice([x, -x]))
    pts += [from_bits(rng.randrange(1, 1 << 52)) for _ in range(128)]
    pts += [-from_bits(rng.randrange(1, 1 << 52)) for _ in range(128)]
    for e in range(-60, -6):
        exp_bits = (1023 + e) << 52
        for _ in range(48):
            x = from_bits(exp_bits | rng.randrange(1 << 52))
            pts += [x, -x]

    # The 2^-6 cut.
    for a in (2.0 ** -6, -(2.0 ** -6)):
        pts += ladder(a, 128)
    for _ in range(2048):
        x = rng.uniform(2.0 ** -7, 2.0 ** -4)
        pts.append(rng.choice([x, -x]))

    # Exponential branch.
    pts += [rng.uniform(-38.0, 0.0) for _ in range(2048)]
    for _ in range(3072):
        pts.append(10.0 ** rng.uniform(math.log10(2.0 ** -4), math.log10(709.7)))

    # The -38 shortcut and where expm1 meets -1.
    pts += ladder(-38.0, 64)
    pts += [rng.uniform(-45.0, -30.0) for _ in range(512)]

    # Overflow boundary (exp's): largest x with a finite exp is
    # 0x1.62e42fefa39efp9.
    pts += ladder(float.fromhex("0x1.62e42fefa39efp9"), 160)

    rows = emit(pts)
    if not self_check(rows):
        return 1
    for x, y in rows:
        print(f"{x.hex()} {y.hex() if math.isfinite(y) else 'inf'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
