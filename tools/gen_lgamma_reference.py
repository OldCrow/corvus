#!/usr/bin/env python3
"""Generate tests/data/lgamma_reference.txt -- correctly rounded lgamma oracle.

Each line: <input-hex-double> <lgamma-hex-double>, rounded to nearest by mpmath
at 50 digits. Specials (poles, infinities, NaN) are covered by test_lgamma.

Point selection (fixed seed, reproducible) targets what the kernel's structure
makes fragile rather than what is convenient to sample:
  - bit-neighbourhoods of BOTH zeros (x = 1, x = 2), where the claim is
    relative accuracy and the exact-t form is what delivers it
  - bit-neighbourhoods of every region boundary (1/2, 3/2, 5/2, X0) and of
    each recurrence step threshold 5/2 + k, where a lane can take a different
    number of product steps than its neighbour
  - the negative axis, split between generic points and the neighbourhoods of
    the |Gamma| = 1 crossings, where lgamma has zeros with no closed form and
    only absolute accuracy is claimed
  - pole approaches from both sides, and the subnormal-argument band where
    lgamma(x) = -log x and the log has to prescale
  - log-spaced large x up to the overflow threshold

Usage:
    python3 tools/gen_lgamma_reference.py > tests/data/lgamma_reference.txt
"""

import math
import random
import struct
import sys

import mpmath as mp

from refgen_common import round_to_double

DPS = 50
mp.mp.dps = DPS

SEED = 20260725
X0 = 8.0
_MIN_NORMAL = 2.0 ** -1022


def as_bits(x: float) -> int:
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def from_bits(u: int) -> float:
    return struct.unpack("<d", struct.pack("<Q", u))[0]


def neighbourhood(x0: float, k: int = 48):
    b = as_bits(x0)
    return [from_bits(b + j) for j in range(-k, k + 1)]


def gamma_one_crossings():
    """Points on the negative axis where |Gamma| = 1, i.e. lgamma = 0.

    Two per interval (-n-1, -n) for n >= 2, bracketing the local |Gamma|
    minimum. These are the only places relative accuracy is unattainable, so
    the reference has to contain them rather than sample around them by luck.
    """
    out = []
    for n in range(2, 12):
        lo, hi = -(n + 1), -n
        # The minimum of |Gamma| sits near the midpoint; bisect each side.
        mid = mp.findroot(lambda z: mp.digamma(z), mp.mpf(lo) + mp.mpf("0.5"))
        for a, b in ((lo + 1e-6, float(mid)), (float(mid), hi - 1e-6)):
            fa, fb = mp.loggamma(mp.mpf(a)), mp.loggamma(mp.mpf(b))
            if mp.sign(fa) == mp.sign(fb):
                continue
            for _ in range(200):
                m = (mp.mpf(a) + mp.mpf(b)) / 2
                if mp.sign(mp.loggamma(m)) == mp.sign(fa):
                    a = float(m)
                else:
                    b = float(m)
            out.append(float(a))
    return out


def oracle(xm):
    """The lgamma oracle value at mpf xm (xm != 0, not a pole)."""
    return mp.loggamma(xm) if xm > 0 else mp.log(abs(mp.gamma(xm)))


def emit(points):
    seen = set()
    rows = []
    for x in points:
        if not math.isfinite(x) or x == 0.0:
            continue
        b = as_bits(x)
        if b in seen:
            continue
        seen.add(b)
        xm = mp.mpf(x)
        if xm < 0 and xm == mp.floor(xm):
            continue  # pole; test_lgamma covers it
        y = oracle(xm)
        if not mp.isfinite(y):
            continue
        yf = round_to_double(y)
        if not math.isfinite(yf):
            continue
        rows.append((x, yf))
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
            x, yf = rows[i]
            y2 = round_to_double(oracle(mp.mpf(x)))
            if as_bits(y2) != as_bits(yf):
                bad.append((x, yf, y2))
    if bad:
        for x, yf, y2 in bad:
            print(f"self-check MISMATCH: x={x.hex()} stored={yf.hex()} "
                  f"recomputed={y2.hex()}", file=sys.stderr)
        return False
    print(f"self-check: {len(sample)} rows re-verified at dps={2 * DPS}: OK",
          file=sys.stderr)
    return True


def main():
    rng = random.Random(SEED)
    pts = []

    # Both zeros, to the bit. Relative accuracy here is the whole point of the
    # zero-centred polynomials.
    pts += neighbourhood(1.0, 96)
    pts += neighbourhood(2.0, 96)

    # Region boundaries and every recurrence step threshold.
    for b in (0.5, 1.5, 2.5, X0) + tuple(2.5 + k for k in range(1, 7)):
        pts += neighbourhood(b)

    # Positive axis, broad.
    pts += [rng.uniform(0.0, 2.5) for _ in range(4096)]
    pts += [rng.uniform(2.5, X0) for _ in range(2048)]
    pts += [rng.uniform(X0, 200.0) for _ in range(2048)]
    pts += [math.exp(rng.uniform(math.log(X0), math.log(1e300)))
            for _ in range(1536)]
    # Up to the overflow threshold: the Stirling product is closest to
    # overflowing here, and the grouping that avoids it is only exercised here.
    pts += [math.exp(rng.uniform(math.log(1e300), math.log(2.55e305)))
            for _ in range(512)]

    # Small and subnormal arguments: lgamma(x) -> -log x, and the log has to
    # prescale to reach a subnormal at all.
    pts += [10.0 ** rng.uniform(-323.0, -1.0) for _ in range(1536)]

    # Negative axis, generic.
    pts += [-rng.uniform(0.0, 20.0) for _ in range(3072)]
    pts += [-math.exp(rng.uniform(math.log(20.0), math.log(1e15)))
            for _ in range(1024)]

    # Negative axis, near the poles from both sides (|u| small).
    for n in range(1, 25):
        for _ in range(24):
            du = rng.choice([-1.0, 1.0]) * 0.5 * rng.random() ** 6
            pts.append(-n + du)

    # Negative axis, near the |Gamma| = 1 crossings.
    for z in gamma_one_crossings():
        pts += neighbourhood(z, 24)
        pts += [z + rng.choice([-1.0, 1.0]) * 10.0 ** rng.uniform(-12, -2)
                for _ in range(24)]

    rows = emit(pts)
    if not self_check(rows):
        return 1
    for x, yf in rows:
        print(f"{x.hex()} {yf.hex()}")
    print(f"emitted {len(rows)} points", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
