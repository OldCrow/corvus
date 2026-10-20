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

Self-check (erfinv rows, issue #61): every row is re-solved independently at
dps 120 (Newton on erf, log-space Newton on erfc -- not mp.erfinv) and the
stored double must be its correct rounding; a negative control nudges stored
values by one ulp and requires each nudge to be flagged. Failure exits
non-zero before any file is written.

Usage:
    python3 tools/gen_erfinv_reference.py
"""

import math
import random
import struct
import sys

import mpmath as mp

from refgen_common import round_to_double

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


HP_DPS = 120
_NUDGE_ROWS = 6


def erfinv_hp(y):
    """x with erf(x) = y, solved at HP_DPS by Newton from a double guess.

    Independent of the generator's oracle (mp.erfinv): plain Newton on
    erf(x) - y for |y| <= 1/2, and log-space Newton on log(erfc(x)) -
    log(1 - y) for |y| > 1/2 (1 - |y| is exact at this precision, and the
    log form keeps the iteration well scaled in the far tail). Returns
    None if it does not converge.
    """
    with mp.workdps(HP_DPS):
        yy = mp.mpf(y)
        if yy == 0:
            return mp.mpf(0)
        neg = yy < 0
        a = abs(yy)
        c = 2 / mp.sqrt(mp.pi)
        # Seed: sqrt(pi)/2 * a is exact enough for small a; for a > 1/2
        # start from the asymptotic w - log(pi*w)/2 with w = -log(1 - a).
        if a <= mp.mpf("0.5"):
            x = a * mp.sqrt(mp.pi) / 2
            for _ in range(200):
                dx = (mp.erf(x) - a) / (c * mp.exp(-x * x))
                x -= dx
                if abs(dx) <= mp.mpf(10) ** (-(HP_DPS - 20)) * abs(x):
                    return -x if neg else x
            return None
        s = 1 - a
        w = -mp.log(s)
        x = mp.sqrt(max(w - mp.log(mp.pi * w) / 2, mp.mpf(1) / 4))
        ls = mp.log(s)
        for _ in range(200):
            ec = mp.erfc(x)
            dx = (mp.log(ec) - ls) / (-c * mp.exp(-x * x) / ec)
            x -= dx
            if abs(dx) <= mp.mpf(10) ** (-(HP_DPS - 20)) * abs(x):
                return -x if neg else x
        return None


def erfinv_row_ok(y, xs):
    """True iff double xs is the correct rounding of erfinv(y)."""
    if y == 0.0:
        return as_bits(xs) == as_bits(y)
    r = erfinv_hp(y)
    if r is None:
        return False
    with mp.workdps(HP_DPS):
        return as_bits(round_to_double(r)) == as_bits(xs)


def erfinv_self_check(rows):
    """Re-solve every erfinv row at HP_DPS and require the stored double to
    be the correct rounding (issue #61), then prove the check can fail: nudge
    sampled stored values by one ulp (math.nextafter) and require each
    nudged row to be flagged. Returns True on success."""
    bad = [(x, yf) for x, yf in rows if not erfinv_row_ok(x, yf)]
    if bad:
        for x, yf in bad[:20]:
            print(f"erfinv self-check MISMATCH: y={x.hex()} stored={yf.hex()}",
                  file=sys.stderr)
        print(f"erfinv self-check: {len(bad)} of {len(rows)} rows FAILED",
              file=sys.stderr)
        return False
    # Negative control: evenly spread nonzero rows, nudged up and down.
    nz = [(x, yf) for x, yf in rows if yf != 0.0]
    step = max(1, len(nz) // _NUDGE_ROWS)
    picks = nz[::step][:_NUDGE_ROWS]
    # Always include the smallest-|y| (subnormal) and the largest-|y| rows.
    picks += [min(nz, key=lambda r: abs(r[0])), max(nz, key=lambda r: abs(r[0]))]
    missed = 0
    for x, yf in picks:
        for direction in (math.inf, -math.inf):
            if erfinv_row_ok(x, math.nextafter(yf, direction)):
                missed += 1
                print(f"erfinv negative control NOT flagged: y={x.hex()} "
                      f"nudged from {yf.hex()} toward {direction}",
                      file=sys.stderr)
    if missed:
        return False
    print(f"erfinv self-check: {len(rows)} rows re-solved at dps={HP_DPS}: OK; "
          f"negative control: {2 * len(picks)} one-ulp nudges all flagged",
          file=sys.stderr)
    return True


def emit(path, points, oracle, in_domain, check=None):
    seen = set()
    rows = []
    for x in points:
        if not isinstance(x, float) or x != x or abs(x) == float("inf"):
            continue
        if not in_domain(x):
            # Bit-neighbourhoods of the domain endpoints deliberately
            # step OUTSIDE the domain (x just past +-1, z just past 2);
            # the old blanket try/except dropped them as a side effect.
            # An explicit filter keeps the drop intentional while a
            # non-converged findroot inside the domain stays FATAL.
            continue
        b = as_bits(x)
        if b in seen:
            continue
        seen.add(b)
        # No try/except here (issue #13 N14.2): a non-converged
        # mp.findroot must abort the run, not silently drop the point.
        y = oracle(x)
        if y is None or not mp.isfinite(y):
            # Legitimate infinite oracle values at domain endpoints
            # (erfinv(+-1) = +-inf), not a convergence failure.
            continue
        yf = round_to_double(y)
        if x == 0.0:
            # mpmath's erfinv(mpf(-0.0)) loses the sign (-> mpf(0)); the
            # double input carries it, and corvus.h promises
            # erfinv(+-0) = +-0.
            yf = x
        if not math.isfinite(yf):
            continue
        rows.append((x, yf))
    # Verify before writing, so a failed check leaves the checked-in file
    # untouched.
    if check is not None and not check(rows):
        return -1
    with open(path, "w") as f:
        for x, yf in rows:
            f.write(f"{x.hex()} {yf.hex()}\n")
    print(f"{path}: emitted {len(rows)} points", file=sys.stderr)
    return len(rows)


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

    # Both zeros as deliberate point-set members (issue #13 N6): dedup is by
    # bit pattern, so +0.0 and -0.0 both survive; neighbourhood(0.0, ...)
    # alone never reaches -0.0 (its b+j >= 0 filter excludes j < 0 at b = 0).
    pts += [0.0, -0.0]

    # A few bit-neighbourhoods scattered across (0, 0.5) and (0.5, 1) so the
    # central polynomial's degree and the seed's own accuracy both get
    # boundary-level (not just interior) coverage.
    for b in (0.1, 0.25, 0.4, 0.6, 0.75, 0.9, 0.95, 0.99):
        pts += neighbourhood(b, 16)
        pts += neighbourhood(-b, 16)

    # Closed domain [-1, 1]: mp.erfinv RAISES outside it; +-1 exactly give
    # +-inf and fall to the non-finite filter as before.
    return emit("tests/data/erfinv_reference.txt", pts,
                lambda x: mp.erfinv(mp.mpf(x)),
                in_domain=lambda x: -1.0 <= x <= 1.0,
                check=erfinv_self_check)


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

    # Open domain (0, 2): at z = 2 exactly erfcinv_mp reduces to -log(0)
    # and can only fail; z <= 0 is out of domain outright.
    return emit("tests/data/erfcinv_reference.txt", pts, erfcinv_mp,
                in_domain=lambda z: 0.0 < z < 2.0)


def main():
    n1 = gen_erfinv()
    if n1 < 0:
        return 1
    n2 = gen_erfcinv()
    return 0 if (n1 > 5000 and n2 > 5000) else 1


if __name__ == "__main__":
    sys.exit(main())
