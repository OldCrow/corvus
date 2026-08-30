#!/usr/bin/env python3
"""Generate tests/data/log1p_reference.txt -- correctly rounded log1p oracle.

Each line: <input-hex-double> <log1p(input)-hex-double>, output rounded once
to nearest double by refgen_common.round_to_double. Specials (x = -1 ->
-inf, x < -1 -> NaN, +/-0 -> +/-0 signed, +/-inf, NaN) are covered by the
smoke/edge test, not here; every input is a finite double > -1, excluding
zero.

Point selection (fixed seed, reproducible):
  - near -1 from above: -1 + 2^-k for k <= 53 and the -1 + j*2^-53 ulp
    ladder. NOTE: no representable double x > -1 has a subnormal 1 + x --
    the closest approach is nextafter(-1, 0) with 1 + x = 2^-53 EXACTLY
    (log1p ~ -36.7), so the "deep corner" bottoms out there by IEEE
    grid construction, not by point-selection choice
  - the Sterbenz-exact zone [-1, -0.5] dense (1 + x is exact there; the
    kernel claim being exercised is TwoSum(1, x) feeding LogDdAny)
  - near 0, both signs: |x| log-spaced 1e-320..1 incl. subnormal inputs
    (log1p(x) ~ x - x^2/2)
  - the |x| ~ 2^-53 seam, bit-neighborhoods (TwoSum lo-only regime border)
  - moderate/huge x log-uniform to the top of the double range
    (log1p(x) ~ log(x))

Usage:
    python tools/gen_log1p_reference.py > tests/data/log1p_reference.txt
"""

import math
import random
import struct
import sys

import mpmath as mp

from refgen_common import round_to_double

DPS = 60
mp.mp.dps = DPS

SEED = 20260830


def as_bits(x: float) -> int:
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def from_bits(u: int) -> float:
    return struct.unpack("<d", struct.pack("<Q", u))[0]


def emit(points):
    seen = set()
    rows = []
    for x in points:
        if not (math.isfinite(x) and x > -1.0 and x != 0.0):
            continue
        b = as_bits(x)
        if b in seen:
            continue
        seen.add(b)
        rows.append((x, round_to_double(mp.log1p(mp.mpf(x)))))
    return rows


def self_check(rows):
    """Every corner row (x within 2^-45 of -1: the -1 + j*2^-53 ladder,
    where 1 + x is smallest and cancellation deepest) plus every 97th row,
    capped at 600, re-verified at 2x dps bit-identical (issue #13 N14.2
    pattern)."""
    corner = [i for i, (x, _) in enumerate(rows)
              if x < -1.0 + 2.0 ** -45]
    corner_set = set(corner)
    every_97th = [i for i in range(0, len(rows), 97) if i not in corner_set]
    sample = (corner + every_97th)[:600]
    bad = []
    with mp.workdps(2 * DPS):
        for i in sample:
            x, y = rows[i]
            y2 = round_to_double(mp.log1p(mp.mpf(x)))
            if as_bits(y2) != as_bits(y):
                bad.append((x, y, y2))
    if bad:
        for x, y, y2 in bad:
            print(f"self-check MISMATCH: x={x.hex()} stored={y.hex()} "
                  f"recomputed={y2.hex()}", file=sys.stderr)
        return False
    print(f"self-check: {len(sample)} rows ({len(corner)} deep-corner) "
          f"re-verified at dps={2 * DPS}: OK", file=sys.stderr)
    return True


def main():
    rng = random.Random(SEED)
    pts = []

    # Cancellation corner: -1 + 2^-k (representable only for k <= 53; the
    # double grid bottoms out at 1 + x = 2^-53) and random points riding it.
    for k in range(1, 54):
        pts.append(-1.0 + math.ldexp(1.0, -k))
    for _ in range(1536):
        k = rng.uniform(1.0, 52.0)
        pts.append(-1.0 + math.ldexp(rng.uniform(1.0, 2.0), -int(k)))
    # bit-neighborhood of -1 from above.
    m1_bits = as_bits(-1.0)
    pts += [from_bits(m1_bits - k) for k in range(1, 257)]  # toward zero

    # Sterbenz-exact zone [-1, -0.5], dense.
    pts += [rng.uniform(-1.0, -0.5) for _ in range(2048)]

    # Near zero, both signs, down into the subnormals.
    for _ in range(3072):
        e = rng.uniform(-320, 0)
        x = 10.0 ** e
        pts.append(rng.choice([x, -x]))
    pts += [from_bits(rng.randrange(1, 1 << 52)) for _ in range(128)]
    pts += [-from_bits(rng.randrange(1, 1 << 52)) for _ in range(128)]

    # The |x| ~ 2^-53 seam.
    for anchor in (2.0 ** -53, 2.0 ** -52, 2.0 ** -54):
        b = as_bits(anchor)
        for k in range(-64, 65):
            x = from_bits(b + k)
            pts += [x, -x]

    # Moderate and huge arguments.
    for _ in range(3072):
        e = rng.uniform(0, 308)
        pts.append(10.0 ** e * rng.uniform(1.0, 9.99))
    pts += [rng.uniform(-0.5, 4.0) for _ in range(1024)]

    rows = emit(pts)
    if not self_check(rows):
        return 1
    for x, y in rows:
        print(f"{x.hex()} {y.hex()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
