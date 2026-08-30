#!/usr/bin/env python3
"""Generate tests/data/exp_reference.txt -- correctly rounded exp oracle.

Each line: <input-hex-double> <exp(input)-hex-double>, output rounded once
to nearest double by refgen_common.round_to_double (single rounding matters
throughout the subnormal-output band). Specials (+/-0 -> 1, +/-inf, NaN)
are covered by the smoke/edge test, not here (erf-reference convention).

Point selection (fixed seed, reproducible):
  - uniform over [-750, 712]: past both saturation points
  - log-spaced tiny |x| to 1e-300 incl. subnormal inputs (exp(x) ~ 1 + x)
  - overflow boundary bit-neighborhood: the largest double with a finite
    correctly rounded exp (~709.7827)
  - underflow-to-zero boundary bit-neighborhood (~-745.1332): where
    RN(exp(x)) crosses from the smallest subnormal to 0
  - SUBNORMAL-OUTPUT band, dense: x in [-745.14, -708.39] -- the one-
    effective-rounding claim of the ScaleTwo assembly is tested here
  - normal/subnormal output boundary bit-neighborhood (~-708.3964)
  - reduction stress: x near k*(ln2/128) half-slot points (table-boundary
    neighborhoods of the exp_dd reduction)

Usage:
    python tools/gen_exp_reference.py > tests/data/exp_reference.txt
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
_MIN_NORMAL = 2.0 ** -1022


def as_bits(x: float) -> int:
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def from_bits(u: int) -> float:
    return struct.unpack("<d", struct.pack("<Q", u))[0]


def emit(points):
    seen = set()
    rows = []
    for x in points:
        if not math.isfinite(x):
            continue
        b = as_bits(x)
        if b in seen:
            continue
        seen.add(b)
        rows.append((x, round_to_double(mp.exp(mp.mpf(x)))))
    return rows


def self_check(rows):
    """Every subnormal-output row plus every 97th row, capped at 600,
    re-verified at 2x dps bit-identical (issue #13 N14.2 pattern)."""
    subnormal = [i for i, (_, y) in enumerate(rows) if 0 < abs(y) < _MIN_NORMAL]
    subnormal_set = set(subnormal)
    every_97th = [i for i in range(0, len(rows), 97) if i not in subnormal_set]
    sample = (subnormal + every_97th)[:600]
    bad = []
    with mp.workdps(2 * DPS):
        for i in sample:
            x, y = rows[i]
            y2 = round_to_double(mp.exp(mp.mpf(x)))
            if as_bits(y2) != as_bits(y):
                bad.append((x, y, y2))
    if bad:
        for x, y, y2 in bad:
            print(f"self-check MISMATCH: x={x.hex()} stored={y.hex()} "
                  f"recomputed={y2.hex()}", file=sys.stderr)
        return False
    print(f"self-check: {len(sample)} rows ({len(subnormal)} subnormal-output) "
          f"re-verified at dps={2 * DPS}: OK", file=sys.stderr)
    return True


def main():
    rng = random.Random(SEED)
    pts = []

    pts += [rng.uniform(-750.0, 712.0) for _ in range(6144)]

    for _ in range(2048):
        e = rng.uniform(-300, 2)
        x = 10.0 ** e
        pts.append(rng.choice([x, -x]))
    pts += [from_bits(rng.randrange(1, 1 << 52)) for _ in range(128)]
    pts += [-from_bits(rng.randrange(1, 1 << 52)) for _ in range(128)]

    # Boundary bit-neighborhoods. The anchors are the nearest doubles to
    # the saturation points; +/-160 ulps blankets the exact crossovers so
    # no hand-derived threshold is trusted.
    for anchor in (709.782712893384, -745.1332191019412, -708.3964185322641):
        b = as_bits(anchor)
        pts += [from_bits(b + k) for k in range(-160, 161)]

    # Subnormal-output band, dense.
    pts += [rng.uniform(-745.14, -708.39) for _ in range(3072)]

    # Reduction/table stress: x near k*(ln2/128) and near half-slot points.
    l = math.log(2.0) / 128.0
    for _ in range(2048):
        k = rng.randint(-95500, 90800)
        kind = rng.random()
        if kind < 0.5:
            d = rng.uniform(-0.5, 0.5) * l
        else:
            d = rng.choice([-1.0, 1.0]) * l * 0.5 * rng.random() ** 8
        pts.append(k * l + d)

    rows = emit(pts)
    if not self_check(rows):
        return 1
    for x, y in rows:
        print(f"{x.hex()} {y.hex()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
