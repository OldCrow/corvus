#!/usr/bin/env python3
"""Generate tests/data/erfinv_reference.txt and erfcinv_reference.txt --
correctly rounded oracles for corvus::erfinv / corvus::erfcinv.

Each line: <input-hex-double> <output-hex-double>, rounded to nearest by
mpmath at 60 digits. Specials (0, +-1/+-2, out-of-domain, NaN) are covered by
test_erfinv.cpp, not here.

Oracle: mpmath's erfinv directly for erfinv points (accurate across its
whole domain); for erfcinv points, root-finding on log(erfc(x)) - log(s) --
mpmath has no erfcinv, and this is also what the design's own tail model is
checked against (see tools/gen_erfinv_data.py).

Point selection (fixed seed, reproducible), per the design's own list of
what is fragile:
  - erfinv dense over (-1, 1), including bit-neighbourhoods of the C/T
    routing boundary y = +-0.5, the extremes +-(1 - k*2^-52), subnormal y,
    and +-0
  - erfcinv log-spaced z down to the smallest subnormal double, so the far
    tail (only reachable through erfcinv) is densely covered
  - bit-neighbourhoods of z = 1 (the zero crossing -- relative accuracy
    there is a design CLAIM, not a convenience), z = 0.5 and z = 1.5 (the
    C/T routing boundary), z near 2, and the mid/far Halley split (kTfar,
    i.e. x = 6)
  - erfinv NEVER reaches the far tail (max |erfinv| = erfcinv(2^-53) < 6),
    so far-tail coverage comes only through erfcinv, by construction here

Usage:
    python3 tools/gen_erfinv_reference.py
"""

import random
import struct
import sys

import mpmath as mp

mp.mp.dps = 60

SEED = 20260725


def as_bits(x: float) -> int:
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def from_bits(u: int) -> float:
    return struct.unpack("<d", struct.pack("<Q", u))[0]


def neighbourhood(x0: float, k: int = 48):
    b = as_bits(abs(x0))
    sign = -1.0 if x0 < 0 or (x0 == 0.0 and mp.sign(x0) < 0) else 1.0
    return [sign * from_bits(b + j) for j in range(-k, k + 1) if b + j >= 0]


def erfcinv_mp(s):
    """x with erfc(x) = s, s in (0, 2).

    Two regimes, matching where each is numerically safe rather than where
    the kernel happens to route (though they end up close):
      s >= 1/2: x = erfinv(1 - s) directly. 1 - s is a NORMAL-magnitude
        mpmath value here (no precision loss forming it), and mpmath's own
        erfinv is robust across it -- including exactly 0 at s = 1, which a
        log-space Newton iteration gets badly wrong (w = -log(1) = 0 makes
        the initial-guess formula's log(pi*w) blow up to -inf/inf and the
        solver returns garbage near zero rather than exactly zero).
      s < 1/2 (down to the smallest subnormal double): root-find in log
        space on w = -log(s). This is NOT merely a style choice: forming
        1 - s directly for subnormal s would need ~1075 bits (mp.dps ~ 324)
        to distinguish it from 1 at all, since s can be as small as 2^-1074;
        working in w keeps everything at ordinary magnitude.
    """
    s = mp.mpf(s)
    if s > 1:
        return -erfcinv_mp(2 - s)  # erfc(x) = 2 - erfc(-x)
    if s >= mp.mpf("0.5"):
        return mp.erfinv(1 - s)
    w = -mp.log(s)
    x0 = mp.sqrt(max(w - mp.log(mp.pi * w) / 2, mp.mpf(1) / 4))
    return mp.findroot(lambda x: mp.log(mp.erfc(x)) - mp.log(s), x0)


def emit(path, points, oracle):
    seen = set()
    n = 0
    with open(path, "w") as f:
        for x in points:
            if not isinstance(x, float) or x != x or abs(x) == float("inf"):
                continue
            b = as_bits(x)
            if b in seen:
                continue
            seen.add(b)
            try:
                y = oracle(x)
            except Exception:
                continue
            if y is None or not mp.isfinite(y):
                continue
            yf = float(y)
            import math
            if not math.isfinite(yf):
                continue
            f.write(f"{x.hex()} {yf.hex()}\n")
            n += 1
    print(f"{path}: emitted {n} points", file=sys.stderr)
    return n


def gen_erfinv():
    rng = random.Random(SEED)
    pts = []

    # Dense over the whole open domain.
    pts += [rng.uniform(-1.0, 1.0) for _ in range(8192)]

    # C/T routing boundary, to the bit, both signs.
    pts += neighbourhood(0.5, 96)
    pts += neighbourhood(-0.5, 96)

    # Extremes: 1 - k*2^-52 and its neighbours, both signs.
    for k in range(0, 64):
        x = 1.0 - k * 2.0 ** -52
        pts += [x, -x]
    pts += neighbourhood(1.0, 64)
    pts += neighbourhood(-1.0, 64)

    # Subnormal y and +-0's immediate neighbourhood.
    pts += [rng.uniform(0, 1) * 2.0 ** rng.randint(-1074, -1000) for _ in range(512)]
    for _ in range(512):
        v = rng.uniform(0, 1) * 2.0 ** rng.randint(-1074, -1000)
        pts.append(rng.choice([v, -v]))
    pts += neighbourhood(0.0, 32)

    # A few bit-neighbourhoods scattered across (0, 0.5) and (0.5, 1) so the
    # central polynomial's degree and the seed's own accuracy both get
    # boundary-level (not just interior) coverage.
    for b in (0.1, 0.25, 0.4, 0.6, 0.75, 0.9, 0.95, 0.99):
        pts += neighbourhood(b, 16)
        pts += neighbourhood(-b, 16)

    return emit("tests/data/erfinv_reference.txt", pts, lambda x: mp.erfinv(mp.mpf(x)))


def gen_erfcinv():
    rng = random.Random(SEED + 1)
    pts = []

    # Log-spaced z down to the smallest subnormal double: this is the ONLY
    # way the far tail (x up to ~27.2) is exercised at all.
    pts += [2.0 ** rng.uniform(-1074.0, -1.0) for _ in range(4096)]
    pts += [rng.uniform(1e-8, 0.5) for _ in range(2048)]

    # z = 1: the zero crossing. Relative accuracy there is a design CLAIM
    # (the exact 1 - z argument feeding the central polynomial), so it gets
    # the same bit-level treatment as lgamma's zeros.
    pts += neighbourhood(1.0, 96)

    # C/T routing boundaries (0.5, 1.5) and the mid/far Halley split.
    pts += neighbourhood(0.5, 96)
    pts += neighbourhood(1.5, 96)
    t_far_x = float(mp.sqrt(-mp.log(mp.erfc(mp.mpf(6)))))
    s_at_far = float(mp.e ** (-mp.mpf(t_far_x) ** 2))
    pts += neighbourhood(s_at_far, 64)

    # z near 2 (erfc saturation on the negative-x side) and z in (1.5, 2).
    pts += neighbourhood(2.0, 96)
    pts += [2.0 - 2.0 ** rng.uniform(-52.0, -1.0) for _ in range(512)]
    pts += [rng.uniform(1.5, 2.0) for _ in range(1024)]

    # Dense over the rest of (0, 2).
    pts += [rng.uniform(0.5, 1.5) for _ in range(4096)]

    return emit("tests/data/erfcinv_reference.txt", pts, erfcinv_mp)


def main():
    n1 = gen_erfinv()
    n2 = gen_erfcinv()
    return 0 if (n1 > 5000 and n2 > 5000) else 1


if __name__ == "__main__":
    sys.exit(main())
