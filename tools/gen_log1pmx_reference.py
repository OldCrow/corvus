#!/usr/bin/env python3
"""Generate tests/data/log1pmx_reference.txt -- correctly rounded
log1pmx(x) = log1p(x) - x oracle.

Each line: <input-hex-double> <log1pmx(input)-hex-double>, output rounded
once to nearest double by refgen_common.round_resolving_ties. Specials
(+-0 -> +0, x = -1 -> -inf, x < -1 -> NaN, +inf -> -inf, NaN) are covered
by the smoke test, not here; every input is a finite nonzero double > -1.
Rows whose value underflows to -0 are kept: they certify the tiny branch's
rounding onto zero, and the gate holds them to bit identity.

The oracle cancels near 0 (log1p(x) ~ x), so it runs at 60 digits plus
twice the decimal exponent of 1/|x|: the result ~ x^2/2 keeps 60 digits.

Point selection (fixed seed, reproducible), by the kernel's own seams
(src/log-inl.h Log1pmxVec, src/dd_special-inl.h Log1pmxDd):
  - near 0 both signs: |x| log-spaced 1e-320..1/16, plus full-mantissa draws
    per binade 2^-60..2^-5
  - the 2^-480 tiny/series cut: +-64-ulp ladders and full-mantissa draws in
    [2^-490, 2^-470], both signs
  - subnormal results: full-mantissa |x| in [2^-540, 2^-505], both signs,
    which spans the normal/subnormal result boundary (|x| = 2^-511) and the
    rounding to zero (|x| ~ 2^-537.5)
  - the 1/16 series/direct cut (kPhiCut): +-128-ulp ladders at +-1/16 and
    dense draws in [1/32, 1/8] both signs
  - the approach to -1: -1 + 2^-k, the ulp ladder above -1, random points
    riding it, and the Sterbenz zone [-1, -1/2] dense
  - moderate and huge x up to DBL_MAX (log1pmx ~ -x)

Usage:
    python tools/gen_log1pmx_reference.py > tests/data/log1pmx_reference.txt
"""

import math
import random
import struct
import sys

import mpmath as mp

from refgen_common import round_resolving_ties

DPS = 60
mp.mp.dps = DPS

SEED = 20261030


def as_bits(x: float) -> int:
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def from_bits(u: int) -> float:
    return struct.unpack("<d", struct.pack("<Q", u))[0]


def _extra(x):
    """Digits lost to cancellation: log1p(x) - x ~ x^2/2 against x."""
    if abs(x) >= 0.5:
        return 10
    return 10 + 2 * int(math.ceil(-math.log10(abs(x))))


def exact(x, dps):
    with mp.workdps(dps + _extra(x)):
        m = mp.mpf(x)
        return +(mp.log1p(m) - m)


def ref_value(x):
    return round_resolving_ties(exact(x, mp.mp.dps),
                                lambda dps, x=x: exact(x, dps))


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
        y = ref_value(x)
        if y is None:  # unresolved near-tie (#89): decline the row
            continue
        rows.append((x, y))
    return rows


def self_check(rows):
    """Tiny-cut and subnormal-result rows, 1/16-cut rows, deep-corner rows
    near -1 and every 97th row overall, recomputed at 2x dps bit-identical.
    Per-arm caps (#35 M2/L3)."""
    tiny = [i for i, (x, _) in enumerate(rows) if abs(x) < 2.0 ** -460][::5][:600]
    cut = [i for i, (x, _) in enumerate(rows)
           if 2.0 ** -5 <= abs(x) <= 2.0 ** -3][::7][:600]
    corner = [i for i, (x, _) in enumerate(rows) if x < -1.0 + 2.0 ** -45][:300]
    every_97th = list(range(0, len(rows), 97))
    sample = sorted(set(tiny + cut + corner + every_97th))
    bad = []
    with mp.workdps(2 * DPS):
        for i in sample:
            x, y = rows[i]
            y2 = ref_value(x)
            if y2 is None:
                y2 = float("nan")
            if as_bits(y2) != as_bits(y):
                bad.append((x, y, y2))
    if bad:
        for x, y, y2 in bad:
            print(f"self-check MISMATCH: x={x.hex()} stored={y.hex()} "
                  f"recomputed={y2.hex()}", file=sys.stderr)
        return False
    n_zero = sum(1 for _, y in rows if y == 0.0)
    n_sub = sum(1 for _, y in rows if 0.0 < abs(y) < 2.0 ** -1022)
    if n_zero == 0 or n_sub == 0:
        print(f"self-check: missing tiny-branch rows (zero {n_zero}, "
              f"subnormal {n_sub})", file=sys.stderr)
        return False
    print(f"self-check: {len(sample)} rows ({len(tiny)} tiny, {len(cut)} cut, "
          f"{len(corner)} corner; {n_sub} subnormal and {n_zero} zero results) "
          f"re-verified at dps={2 * DPS}: OK", file=sys.stderr)
    return True


def ladder(anchor, k):
    b = as_bits(anchor)
    return [from_bits(b + j) for j in range(-k, k + 1)]


def full_mantissa_binades(rng, elo, ehi, per):
    out = []
    for e in range(elo, ehi + 1):
        exp_bits = (1023 + e) << 52
        for _ in range(per):
            x = from_bits(exp_bits | rng.randrange(1 << 52))
            out += [x, -x]
    return out


def main():
    rng = random.Random(SEED)
    pts = []

    # Near 0, both signs, down into the subnormal inputs.
    for _ in range(3072):
        x = 10.0 ** rng.uniform(-320, math.log10(1.0 / 16))
        pts.append(rng.choice([x, -x]))
    pts += full_mantissa_binades(rng, -60, -5, 32)

    # The 2^-480 tiny/series cut.
    for a in (2.0 ** -480, -(2.0 ** -480)):
        pts += ladder(a, 64)
    pts += full_mantissa_binades(rng, -490, -470, 16)

    # Subnormal results and the rounding to zero.
    pts += full_mantissa_binades(rng, -540, -505, 40)

    # The 1/16 cut.
    for a in (1.0 / 16, -1.0 / 16):
        pts += ladder(a, 128)
    for _ in range(2048):
        x = rng.uniform(1.0 / 32, 1.0 / 8)
        pts.append(rng.choice([x, -x]))

    # The approach to -1 (1 + x is exact on [-1, -1/2]; the grid bottoms
    # out at 1 + x = 2^-53).
    for k in range(1, 54):
        pts.append(-1.0 + math.ldexp(1.0, -k))
    pts += [from_bits(as_bits(-1.0) - k) for k in range(1, 257)]
    for _ in range(1024):
        pts.append(-1.0 + math.ldexp(rng.uniform(1.0, 2.0),
                                     -int(rng.uniform(1.0, 52.0))))
    pts += [rng.uniform(-1.0, -0.5) for _ in range(2048)]

    # Moderate and huge arguments.
    pts += [rng.uniform(-0.5, 4.0) for _ in range(1024)]
    for _ in range(3072):
        pts.append(10.0 ** rng.uniform(0, 308) * rng.uniform(1.0, 1.79))
    pts.append(sys.float_info.max)

    rows = emit(pts)
    if not self_check(rows):
        return 1
    for x, y in rows:
        print(f"{x.hex()} {y.hex()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
