#!/usr/bin/env python3
"""Generate tests/data/cos_reference.txt / sin_reference.txt -- correctly
rounded cos/sin oracle over the FULL double range.

Each line: <input-hex-double> <fn(input)-hex-double>, output rounded once to
nearest double by refgen_common.round_to_double. The point set is identical
for both functions (same seed, selection independent of --fn), so the two
files share their inputs row-for-row. Specials (+/-0, +/-inf, NaN) are
covered by the smoke/edge test, not here (erf-reference convention).

Point selection (fixed seed, reproducible):
  small region (|x| <= 2^23, the exact-split reduction):
  - uniform over +/-2^23
  - log-spaced small magnitudes to 1e-300 incl. subnormals, BOTH signs
    (sin(x) ~ x; negative subnormals added by #35 L8)
  - bit-neighborhoods of k*pi/2 (random k <= n_max = 5,340,354): the
    reduction-boundary stress rows
  - per-exponent CF worst cases for e in [2, 22] (deepest cancellation
    reachable inside the small region; tools/trig_common.py)
  large region (|x| > 2^23, the Payne-Hanek reduction):
  - per-exponent CF worst cases, e = 23..1023, top 2 each -- includes
    binary64's literature worst case (m=0x16ac5b262ca1ff, e=849,
    |r| = 2^-60.89)
  - random mantissas on an exponent stride (breadth between worst cases)
  - the consumer specials' finite points verbatim: 2^23, nextafter(2^23),
    2^24, 1e9, 1e300 (libstats/libhmm trig specials gates)
  Both signs throughout.

Oracle precision scales with the exponent: cos/sin of x ~ 2^e needs
~0.302*e digits just to survive the reduction, so each point is evaluated
at workdps(ceil(0.302*e) + 80); the self-check re-verifies a deterministic
sample at twice that and requires bit-identical doubles.

Usage:
    python tools/gen_trig_reference.py --fn cos > tests/data/cos_reference.txt
    python tools/gen_trig_reference.py --fn sin > tests/data/sin_reference.txt
"""

import argparse
import math
import random
import struct
import sys

import mpmath as mp

from refgen_common import round_to_double
from trig_common import PH_E_MAX, PH_E_MIN, worst_m_for_exponent

SEED = 20260830
D_MAX = float(2 ** 23)
N_MAX = 5340354


def as_bits(x: float) -> int:
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def from_bits(u: int) -> float:
    return struct.unpack("<d", struct.pack("<Q", u))[0]


def work_dps(x: float) -> int:
    _, ex = math.frexp(abs(x))
    return max(60, int(0.302 * max(ex, 0)) + 80)


def oracle(fn, x: float, dps_scale: int = 1) -> float:
    with mp.workdps(dps_scale * work_dps(x)):
        return round_to_double(fn(mp.mpf(x)))


def emit(fn, points):
    seen = set()
    rows = []
    for x in points:
        if not math.isfinite(x):
            continue
        b = as_bits(x)
        if b in seen:
            continue
        seen.add(b)
        rows.append((x, oracle(fn, x)))
    return rows


def self_check(fn, rows):
    """Deterministic sample re-verified at 2x dps, bit-identical (issue #13
    N14.2 pattern). EVERY huge-region row with |x| >= 1e250 is in the
    sample (they lean hardest on the oracle's internal reduction; these
    rows also sit above both consumers' 2^23 domain, so this replay is
    their only oracle cross-check), plus every 97th row. No truncating
    global cap: #35 M2 found the old (huge + every_97th)[:600] silently
    dropped the 334 HIGHEST-exponent huge rows and the entire every-97th
    arm."""
    huge = [i for i, (x, _) in enumerate(rows) if abs(x) >= 1e250]
    huge_set = set(huge)
    every_97th = [i for i in range(0, len(rows), 97) if i not in huge_set]
    sample = huge + every_97th
    bad = []
    for i in sample:
        x, y = rows[i]
        y2 = oracle(fn, x, dps_scale=2)
        if as_bits(y2) != as_bits(y):
            bad.append((x, y, y2))
    if bad:
        for x, y, y2 in bad:
            print(f"self-check MISMATCH: x={x.hex()} stored={y.hex()} "
                  f"recomputed={y2.hex()}", file=sys.stderr)
        return False
    print(f"self-check: {len(sample)} rows ({len(huge)} huge, all of them) "
          f"re-verified at 2x dps: OK",
          file=sys.stderr)
    return True


def build_points():
    rng = random.Random(SEED)
    pts = []

    # --- small region ---
    pts += [rng.uniform(-D_MAX, D_MAX) for _ in range(4096)]

    for _ in range(1536):
        e = rng.uniform(-300, 6)
        x = 10.0 ** e
        pts.append(rng.choice([x, -x]))
    pts += [from_bits(rng.randrange(1, 1 << 52)) for _ in range(256)]  # subnormals
    # Negative subnormals (#35 L8): the stratum above never negated, so the
    # kernel's |x| symmetry met no negative-subnormal reference row. A
    # dedicated rng stream keeps every row above byte-identical (pure
    # append, same pattern as gen_log1p's H1 stratum).
    rng_l8 = random.Random(SEED ^ 0x1B8)
    pts += [-from_bits(rng_l8.randrange(1, 1 << 52)) for _ in range(256)]

    with mp.workdps(60):
        for _ in range(3584):
            k = rng.randint(1, N_MAX)
            x0 = float(mp.pi / 2 * k)
            j = rng.choice([0, 1, -1, 2, -2, 4, -4, 8, -8])
            x = from_bits(as_bits(x0) + j)
            pts.append(rng.choice([x, -x]))

    for e in range(2, PH_E_MIN):
        for m, _depth in worst_m_for_exponent(e, 2):
            x = math.ldexp(m, e - 52)
            if abs(x) <= D_MAX:
                pts += [x, -x]

    # --- large region ---
    for e in range(PH_E_MIN, PH_E_MAX + 1):
        for m, _depth in worst_m_for_exponent(e, 2):
            x = math.ldexp(m, e - 52)
            pts += [x, -x]

    for e in range(PH_E_MIN, PH_E_MAX + 1, 5):
        for _ in range(2):
            x = math.ldexp(rng.randrange(1 << 52, 1 << 53), e - 52)
            pts += [x, -x]

    for x in (D_MAX, math.nextafter(D_MAX, math.inf), 2.0 * D_MAX, 1e9, 1e300):
        pts += [x, -x]

    return pts


def main():
    ap = argparse.ArgumentParser(allow_abbrev=False)
    ap.add_argument("--fn", required=True, choices=["cos", "sin"])
    args = ap.parse_args()
    fn = mp.cos if args.fn == "cos" else mp.sin

    rows = emit(fn, build_points())
    if not self_check(fn, rows):
        return 1
    for x, y in rows:
        print(f"{x.hex()} {y.hex()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
