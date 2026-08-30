#!/usr/bin/env python3
"""Generate tests/data/log_reference.txt -- correctly rounded log oracle.

Each line: <input-hex-double> <log(input)-hex-double>, output rounded once
to nearest double by refgen_common.round_to_double. Specials (+/-0 -> -inf,
x<0 -> NaN, +inf, NaN) are covered by the smoke/edge test, not here
(erf-reference convention); every input in this file is a positive finite
double.

Point selection (fixed seed, reproducible):
  - log-uniform over the whole positive normal range (random exponent in
    [-1022, 1023], random mantissa)
  - subnormal inputs: random subnormal bit patterns plus denorm_min and
    min-normal bit-neighborhoods (the prescale seam)
  - near 1 from both sides: +/-512-ulp bit-neighborhood of 1.0, plus
    1 +/- 10^-e sprays (relative accuracy is hardest here; the log_dd
    table's centred-mantissa design is what is being exercised)
  - log_dd slot-boundary stress: x = (1 + j/128) * 2^k +/- small residuals
    (mirrors the erf grid-point stratum; 128 = kLogN slots)
  - sqrt(2) mantissa-centering boundary bit-neighborhoods across exponents

Usage:
    python tools/gen_log_reference.py > tests/data/log_reference.txt
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
        if not (math.isfinite(x) and x > 0.0):
            continue
        b = as_bits(x)
        if b in seen:
            continue
        seen.add(b)
        rows.append((x, round_to_double(mp.log(mp.mpf(x)))))
    return rows


def self_check(rows):
    """Every subnormal-INPUT row plus every 97th row, capped at 600,
    re-verified at 2x dps bit-identical (issue #13 N14.2 pattern)."""
    subnormal = [i for i, (x, _) in enumerate(rows) if x < 2.0 ** -1022]
    subnormal_set = set(subnormal)
    every_97th = [i for i in range(0, len(rows), 97) if i not in subnormal_set]
    sample = (subnormal + every_97th)[:600]
    bad = []
    with mp.workdps(2 * DPS):
        for i in sample:
            x, y = rows[i]
            y2 = round_to_double(mp.log(mp.mpf(x)))
            if as_bits(y2) != as_bits(y):
                bad.append((x, y, y2))
    if bad:
        for x, y, y2 in bad:
            print(f"self-check MISMATCH: x={x.hex()} stored={y.hex()} "
                  f"recomputed={y2.hex()}", file=sys.stderr)
        return False
    print(f"self-check: {len(sample)} rows ({len(subnormal)} subnormal-input) "
          f"re-verified at dps={2 * DPS}: OK", file=sys.stderr)
    return True


def main():
    rng = random.Random(SEED)
    pts = []

    # Whole positive normal range, log-uniform.
    for _ in range(6144):
        e = rng.randint(-1022, 1023)
        m = 1.0 + rng.random()
        pts.append(math.ldexp(m, e))

    # Subnormals: random patterns + the two boundary neighborhoods.
    pts += [from_bits(rng.randrange(1, 1 << 52)) for _ in range(2048)]
    for k in range(1, 129):
        pts.append(from_bits(k))                       # denorm_min ladder
    min_normal_bits = as_bits(2.0 ** -1022)
    pts += [from_bits(min_normal_bits + k) for k in range(-128, 129)]

    # Near 1, both sides: dense bit-neighborhood + relative sprays.
    one_bits = as_bits(1.0)
    pts += [from_bits(one_bits + k) for k in range(-512, 513)]
    for _ in range(2048):
        e = rng.uniform(-16, -1)
        d = 10.0 ** e
        pts.append(rng.choice([1.0 + d, 1.0 - d]))

    # log_dd slot boundaries: (1 + j/128) * 2^k with small residuals.
    for _ in range(2048):
        j = rng.randint(0, 128)
        k = rng.randint(-64, 64)
        kind = rng.random()
        if kind < 0.5:
            d = rng.uniform(-0.5, 0.5) / 128.0
        else:
            d = rng.choice([-1.0, 1.0]) * rng.random() ** 8 / 256.0
        pts.append(math.ldexp(1.0 + j / 128.0 + d, k))

    # sqrt(2) centering boundary across exponents.
    for k in range(-16, 17):
        b = as_bits(math.ldexp(math.sqrt(2.0), k))
        pts += [from_bits(b + j) for j in range(-32, 33)]

    rows = emit(pts)
    if not self_check(rows):
        return 1
    for x, y in rows:
        print(f"{x.hex()} {y.hex()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
