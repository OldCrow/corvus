#!/usr/bin/env python3
"""Generate the two dd_special gate sets that tests/data/dd_special_reference.txt
does not cover (#58):

  tests/data/expm1_dd_reference.txt   w.hi w.lo  want.hi want.lo
      Expm1Dd(w) = e^w - 1 for a dd w, to dd precision. Expm1Dd had no gate
      at all; #58 found its exponential branch at ~2^-61.6 relative, not the
      ~2^-95 its comment claimed.
  tests/data/log1pmx_ulo_reference.txt  u.hi u.lo  want.hi want.lo
      Log1pmxDd(u) = u - log1p(u) with a NONZERO u.lo. The existing set
      (written by gen_gamma_reference.py) drives u.lo = 0 only, so the u.lo
      correction and the 2-arg overload never ran under a gate.

The want pair is the exact value rounded to a dd: hi = RN(v), lo = RN(v - hi).
The input pair (hi, lo) is an exact dd by construction (|lo| <= ulp(hi)/2),
and the oracle evaluates at hi + lo exactly in mpf.

Point sets (own deterministic RNG, seed 58; no other generator's stream):
  expm1: w.hi of both signs, log-uniform over [2^-60, 2^4], denser across
         [2^-12, 2^-4] (the series/exponential seams the kernel has used:
         2^-10 before #58, 2^-6 after), plus w in [-40, -2^4] where e^w - 1
         approaches -1. Half the rows carry a nonzero w.lo.
  log1pmx: u.hi of both signs over [2^-40, 1/16] (the series branch, where
         the u.lo correction lives) and (1/16, 1] / [-1/2, -1/16) (the
         LogDdAny branch, which takes w = 1 + u with u.lo folded in). Every
         row carries a nonzero u.lo.

Self-checks (exit 1 on any failure; files untouched):
  (a) dps ladder: every value computed at dps 60 and dps 100 must agree to
      2^-120 relative.
  (b) every want pair is normalized: |lo| <= ulp(hi)/2 and hi = RN(hi + lo).
  (c) every input pair is normalized the same way.
  (d) negative control: a want pair perturbed by 2^-100 relative must be
      rejected by the same agreement test (c) uses.
  (e) row counts at least the planned sizes (no silent shrinkage).

Usage (from the repo root):
    python tools/gen_dd_special_gates.py
Writes both files directly.
"""
import math
import random
import sys

import mpmath as mp

SEED = 58
EXPM1_PATH = "tests/data/expm1_dd_reference.txt"
LOG1PMX_PATH = "tests/data/log1pmx_ulo_reference.txt"
AGREE = mp.mpf(2) ** -120


def hexd(x):
    return float(x).hex()


def ulp(x):
    return math.ulp(x)


def random_lo(rng, hi):
    """A nonzero lo with |lo| <= ulp(hi)/2 that leaves (hi, lo) normalized."""
    while True:
        lo = (rng.random() - 0.5) * ulp(hi)
        if lo != 0.0 and hi + lo == hi:
            return lo


def dd_of(v):
    hi = float(v)
    lo = float(v - mp.mpf(hi))
    return hi, lo


def normalized(hi, lo):
    return abs(lo) <= ulp(hi) / 2 and hi + lo == hi


def expm1_mp(whi, wlo, dps):
    with mp.workdps(dps):
        return mp.expm1(mp.mpf(whi) + mp.mpf(wlo))


def log1pmx_mp(uhi, ulo, dps):
    with mp.workdps(dps):
        u = mp.mpf(uhi) + mp.mpf(ulo)
        return u - mp.log1p(u)


def loguni(rng, lo, hi):
    return math.exp(rng.uniform(math.log(lo), math.log(hi)))


def expm1_points(rng):
    pts = []
    for _ in range(1600):  # whole range, both signs
        w = loguni(rng, 2.0 ** -60, 16.0) * rng.choice((-1.0, 1.0))
        pts.append(w)
    for _ in range(1200):  # the seam band
        w = loguni(rng, 2.0 ** -12, 2.0 ** -4) * rng.choice((-1.0, 1.0))
        pts.append(w)
    for _ in range(200):  # e^w - 1 -> -1
        pts.append(-rng.uniform(16.0, 40.0))
    for k in (-10, -6):  # both seams, bit-stepped
        c = 2.0 ** k
        for s in (-1.0, 1.0):
            for x in (c, math.nextafter(c, 0.0), math.nextafter(c, 1.0)):
                pts.append(s * x)
    out = []
    for i, whi in enumerate(pts):
        wlo = random_lo(rng, whi) if i % 2 else 0.0
        out.append((whi, wlo))
    return out


def log1pmx_points(rng):
    out = []
    for _ in range(1400):  # series branch
        u = loguni(rng, 2.0 ** -40, 1.0 / 16.0) * rng.choice((-1.0, 1.0))
        out.append(u)
    for _ in range(300):  # LogDdAny branch, positive
        out.append(rng.uniform(1.0 / 16.0, 1.0))
    for _ in range(300):  # LogDdAny branch, negative
        out.append(-rng.uniform(1.0 / 16.0, 0.5))
    for s in (-1.0, 1.0):  # the 1/16 cut, bit-stepped
        c = 1.0 / 16.0
        for x in (c, math.nextafter(c, 0.0), math.nextafter(c, 1.0)):
            out.append(s * x)
    return [(u, random_lo(rng, u)) for u in out]


def build(points, fn, label):
    rows, bad = [], 0
    for hi_in, lo_in in points:
        if not normalized(hi_in, lo_in):
            print(f"  FAILED (c): {label} input ({hi_in!r}, {lo_in!r}) "
                  f"not normalized", file=sys.stderr)
            bad += 1
            continue
        v60 = fn(hi_in, lo_in, 60)
        v100 = fn(hi_in, lo_in, 100)
        with mp.workdps(100):
            if v100 == 0 or abs((v60 - v100) / v100) > AGREE:
                print(f"  FAILED (a): {label} dps ladder at "
                      f"({hi_in!r}, {lo_in!r})", file=sys.stderr)
                bad += 1
                continue
            hi, lo = dd_of(v100)
        if not normalized(hi, lo):
            print(f"  FAILED (b): {label} want pair at ({hi_in!r}, {lo_in!r})",
                  file=sys.stderr)
            bad += 1
            continue
        rows.append((hi_in, lo_in, hi, lo, v100))
    return rows, bad


def negative_control(rows, label):
    """A want pair 2^-100 relative off must fail the 2^-120 agreement test."""
    hi_in, lo_in, hi, lo, v = rows[len(rows) // 2]
    with mp.workdps(100):
        perturbed = (mp.mpf(hi) + mp.mpf(lo)) * (1 + mp.mpf(2) ** -100)
        caught = abs((perturbed - v) / v) > AGREE
    if not caught:
        print(f"  FAILED (d): {label} negative control not rejected",
              file=sys.stderr)
    return caught


def write(path, rows, header):
    with open(path, "w", newline="\n") as f:
        for line in header:
            f.write(f"# {line}\n")
        for hi_in, lo_in, hi, lo, _ in rows:
            f.write(f"{hexd(hi_in)} {hexd(lo_in)} {hexd(hi)} {hexd(lo)}\n")


def main():
    rng = random.Random(SEED)
    e_pts = expm1_points(rng)
    l_pts = log1pmx_points(rng)
    e_rows, e_bad = build(e_pts, expm1_mp, "expm1")
    l_rows, l_bad = build(l_pts, log1pmx_mp, "log1pmx")
    ok = e_bad == 0 and l_bad == 0
    ok &= negative_control(e_rows, "expm1")
    ok &= negative_control(l_rows, "log1pmx")
    for label, rows, want in (("expm1", e_rows, len(e_pts)),
                              ("log1pmx", l_rows, len(l_pts))):
        if len(rows) < want:
            print(f"  FAILED (e): {label} {len(rows)} rows < {want}",
                  file=sys.stderr)
            ok = False
    if not ok:
        print("Self-checks FAILED -- not touching reference files.",
              file=sys.stderr)
        return 1
    write(EXPM1_PATH, e_rows,
          ["Expm1Dd gate (#58): w.hi w.lo want.hi want.lo",
           "tools/gen_dd_special_gates.py, seed 58, mpmath dps 60/100"])
    write(LOG1PMX_PATH, l_rows,
          ["Log1pmxDd nonzero-u.lo gate (#58): u.hi u.lo want.hi want.lo",
           "tools/gen_dd_special_gates.py, seed 58, mpmath dps 60/100"])
    print(f"expm1: {len(e_rows)} rows; log1pmx: {len(l_rows)} rows; "
          f"self-checks passed", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
