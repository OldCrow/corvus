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

import mpmath as mp

mp.mp.dps = 50

SEED = 20260721


def as_bits(x: float) -> int:
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def from_bits(u: int) -> float:
    return struct.unpack("<d", struct.pack("<Q", u))[0]


def emit(points):
    seen = set()
    for x in points:
        b = as_bits(x)
        if b in seen:
            continue
        seen.add(b)
        y = float(mp.erfc(mp.mpf(x)))  # float() rounds to nearest
        print(f"{x.hex()} {y.hex()}")


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

    emit(pts)
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
