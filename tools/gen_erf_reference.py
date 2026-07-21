#!/usr/bin/env python3
"""Generate tests/data/erf_reference.txt -- correctly rounded erf oracle.

Each line: <input-hex-double> <erf-hex-double>, where the output is
round-to-nearest(erf(input)) computed by mpmath at 50 digits. Used by
test_erf_ulp to measure the kernel's max ULP deviation. Specials
(0, inf, NaN) are covered by test_erf, not here.

Point selection (fixed seed, reproducible):
  - uniform over [-6.5, 6.5] (covers the saturated region on both sides)
  - clustered at grid points r_j = j/256 +/- small residuals (worst case for
    the table's local-Taylor correction)
  - log-spaced magnitudes down to 1e-300 (small-x linear regime)
  - bit-neighborhoods of the saturation bound kErfAMax

Usage:
    python3 tools/gen_erf_reference.py > tests/data/erf_reference.txt
"""

import random
import struct

import mpmath as mp

mp.mp.dps = 50

SEED = 20260720


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
        y = float(mp.erf(mp.mpf(x)))  # float() rounds to nearest
        print(f"{x.hex()} {y.hex()}")


def main():
    rng = random.Random(SEED)
    pts = []

    # Uniform sweep, both signs, into the saturated region.
    pts += [rng.uniform(-6.5, 6.5) for _ in range(8192)]

    # Grid-point neighborhoods: r_j plus residuals across the full |d| <= h/2
    # range and tiny offsets that stress the compensated E_hi + E_lo add.
    h = 1.0 / 256.0
    for _ in range(4096):
        j = rng.randint(0, 1536)
        kind = rng.random()
        if kind < 0.5:
            d = rng.uniform(-0.5, 0.5) * h
        else:
            d = rng.choice([-1.0, 1.0]) * h * 0.5 * rng.random() ** 8
        x = j * h + d
        pts.append(rng.choice([x, -x]))

    # Small magnitudes: erf(x) ~ (2/sqrt(pi)) x.
    for _ in range(2048):
        e = rng.uniform(-300, -1)
        x = 10.0 ** e
        pts.append(rng.choice([x, -x]))

    # Saturation bound bit-neighborhood.
    amax = 5.921587195794507
    b = as_bits(amax)
    for k in range(-64, 65):
        x = from_bits(b + k)
        pts += [x, -x]

    emit(pts)
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
