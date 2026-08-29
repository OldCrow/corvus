#!/usr/bin/env python3
"""Generate tests/data/erfc_reference.txt -- correctly rounded erfc oracle.

Each line: <input-hex-double> <erfc-hex-double>, output rounded to nearest
by mpmath at 50 digits. Specials (0, inf, NaN) are covered by test_erfc.

Point selection (fixed seed, reproducible):
  - uniform over [-6.5, 29] (core, tail, and both saturated ends)
  - clustered at erf-table grid points (core region worst case)
  - log-spaced small magnitudes (erfc ~ 1 - 2x/sqrt(pi))
  - dense tail coverage in a in [6, 28], plus bit-neighborhoods of the
    region boundaries 6, 10, 17 and the erf saturation bound
  - the subnormal-result band a in [26.2, 27.9]

Usage:
    python3 tools/gen_erfc_reference.py > tests/data/erfc_reference.txt
"""

import random
import struct
import sys

import mpmath as mp

from refgen_common import round_to_double

DPS = 50
mp.mp.dps = DPS

SEED = 20260721
_MIN_NORMAL = 2.0 ** -1022


def as_bits(x: float) -> int:
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def from_bits(u: int) -> float:
    return struct.unpack("<d", struct.pack("<Q", u))[0]


def emit(points):
    seen = set()
    rows = []
    for x in points:
        b = as_bits(x)
        if b in seen:
            continue
        seen.add(b)
        y = round_to_double(mp.erfc(mp.mpf(x)))
        rows.append((x, y))
    return rows


def self_check(rows):
    """Re-verify a deterministic sample at double dps (issue #13 N14.2).

    Sample = every subnormal-output row plus every 97th row, capped at 500
    rows total (no rng -- deterministic by construction). Recomputes the
    oracle at 2*DPS and requires a bit-identical double, compared by packed
    bit pattern so -0.0 vs +0.0 is caught.
    """
    subnormal = [i for i, (_, y) in enumerate(rows) if 0 < abs(y) < _MIN_NORMAL]
    subnormal_set = set(subnormal)
    every_97th = [i for i in range(0, len(rows), 97) if i not in subnormal_set]
    sample = (subnormal + every_97th)[:500]
    bad = []
    with mp.workdps(2 * DPS):
        for i in sample:
            x, y = rows[i]
            y2 = round_to_double(mp.erfc(mp.mpf(x)))
            if as_bits(y2) != as_bits(y):
                bad.append((x, y, y2))
    if bad:
        for x, y, y2 in bad:
            print(f"self-check MISMATCH: x={x.hex()} stored={y.hex()} "
                  f"recomputed={y2.hex()}", file=sys.stderr)
        return False
    print(f"self-check: {len(sample)} rows re-verified at dps={2 * DPS}: OK",
          file=sys.stderr)
    return True


def main():
    rng = random.Random(SEED)
    pts = []

    # Uniform sweep.
    pts += [rng.uniform(-6.5, 29.0) for _ in range(6144)]

    # Core-region grid clusters (both signs).
    h = 1.0 / 256.0
    for _ in range(3072):
        j = rng.randint(0, 1536)
        kind = rng.random()
        if kind < 0.5:
            d = rng.uniform(-0.5, 0.5) * h
        else:
            d = rng.choice([-1.0, 1.0]) * h * 0.5 * rng.random() ** 8
        x = j * h + d
        pts.append(rng.choice([x, -x]))

    # Small magnitudes.
    for _ in range(1536):
        e = rng.uniform(-300, -1)
        x = 10.0 ** e
        pts.append(rng.choice([x, -x]))

    # Tail: uniform in a and uniform in log(a).
    pts += [rng.uniform(6.0, 28.0) for _ in range(3072)]
    import math
    pts += [math.exp(rng.uniform(math.log(6.0), math.log(28.0))) for _ in range(1024)]

    # Subnormal-result band.
    pts += [rng.uniform(26.2, 27.9) for _ in range(1024)]

    # Bit-neighborhoods of region boundaries and the erf saturation bound.
    for b0 in (6.0, 10.0, 17.0, 5.921587195794507):
        b = as_bits(b0)
        for k in range(-64, 65):
            x = from_bits(b + k)
            pts += [x, -x]

    rows = emit(pts)
    if not self_check(rows):
        return 1
    for x, y in rows:
        print(f"{x.hex()} {y.hex()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
