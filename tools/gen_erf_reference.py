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
  - subnormal x (APPEND mode, issue #61): 3001 rows, both signs, the 64
    smallest subnormals, the 2^-1022 boundary, log-uniform fill; added to
    the existing file without regenerating it

Usage:
    python3 tools/gen_erf_reference.py > tests/data/erf_reference.txt
    python3 tools/gen_erf_reference.py --subnormal-append [path]
        (appends in place; refuses if the file already has subnormal x)
"""

import math
import random
import struct
import sys

import mpmath as mp

from refgen_common import round_to_double

DPS = 50
mp.mp.dps = DPS

SEED = 20260720
_MIN_NORMAL = 2.0 ** -1022


def as_bits(x: float) -> int:
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def from_bits(u: int) -> float:
    return struct.unpack("<d", struct.pack("<Q", u))[0]


def emit(points):
    seen = set()
    rows = []
    for x in points:
        b = as_bits(x)
        if b in seen:
            continue
        seen.add(b)
        y = round_to_double(mp.erf(mp.mpf(x)))
        rows.append((x, y))
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
            x, y = rows[i]
            y2 = round_to_double(mp.erf(mp.mpf(x)))
            if as_bits(y2) != as_bits(y):
                bad.append((x, y, y2))
    if bad:
        for x, y, y2 in bad:
            print(f"self-check MISMATCH: x={x.hex()} stored={y.hex()} "
                  f"recomputed={y2.hex()}", file=sys.stderr)
        return False
    print(f"self-check: {len(sample)} rows re-verified at dps={2 * DPS}: OK",
          file=sys.stderr)
    return True


SUBNORMAL_SEED = SEED + 1000
SUBNORMAL_ROWS = 3001  # appended rows; ODD (house rule), see append_subnormal


def subnormal_points(rng, n, taken):
    """n distinct subnormal-x doubles (both signs) not already in `taken`.

    Anchors first: the 64 smallest subnormals and the 2^-1022 boundary
    (64 subnormals below it plus the smallest normal itself), both signs.
    Then log-uniform fill: m = floor(2^U), U ~ uniform(0, 52), so every
    binade of the subnormal range is hit equally often.
    """
    pts = []
    for k in range(1, 65):
        pts += [from_bits(k), -from_bits(k)]
    base = as_bits(_MIN_NORMAL)
    for j in range(-64, 1):
        x = from_bits(base + j)
        pts += [x, -x]
    out, seen = [], set(taken)
    for x in pts:
        if as_bits(x) not in seen:
            seen.add(as_bits(x))
            out.append(x)
    while len(out) < n:
        m = int(2.0 ** rng.uniform(0.0, 52.0))
        if m < 1:
            m = 1
        x = from_bits(m)
        x = x if rng.random() < 0.5 else -x
        if as_bits(x) not in seen:
            seen.add(as_bits(x))
            out.append(x)
    return out[:n]


def verify_appended(rows):
    """Re-verify EVERY appended row at 2*DPS (bit-identical double), then a
    negative control: a one-ulp nudge (math.nextafter) of sampled stored
    values must be flagged by the same comparison."""
    with mp.workdps(2 * DPS):
        want = [round_to_double(mp.erf(mp.mpf(x))) for x, _ in rows]
    bad = [(x, y, w) for (x, y), w in zip(rows, want) if as_bits(y) != as_bits(w)]
    if bad:
        for x, y, w in bad[:20]:
            print(f"subnormal self-check MISMATCH: x={x.hex()} stored={y.hex()} "
                  f"recomputed={w.hex()}", file=sys.stderr)
        return False
    flagged = 0
    picks = list(range(0, len(rows), max(1, len(rows) // 16)))[:16]
    for i in picks:
        for direction in (math.inf, -math.inf):
            if as_bits(math.nextafter(rows[i][1], direction)) == as_bits(want[i]):
                print(f"subnormal negative control NOT flagged: x={rows[i][0].hex()}",
                      file=sys.stderr)
                return False
            flagged += 1
    print(f"subnormal self-check: {len(rows)} rows re-verified at dps={2 * DPS}: OK; "
          f"negative control: {flagged} one-ulp nudges flagged", file=sys.stderr)
    return True


def append_subnormal(path):
    """Append a subnormal-x stratum to an existing erf reference file.

    Existing rows are NOT regenerated or rewritten. Refuses to run if the
    file already holds a subnormal x (idempotency guard) or if the new total
    row count would be even (house rule: never a multiple of 2/4/8).
    """
    with open(path, "rb") as f:
        raw = f.read()
    lines = raw.decode("ascii").splitlines()
    taken = set()
    for ln in lines:
        x = float.fromhex(ln.split()[0])
        if x != 0.0 and abs(x) < _MIN_NORMAL:
            print(f"{path} already holds subnormal x; refusing to append",
                  file=sys.stderr)
            return 1
        taken.add(as_bits(x))
    if (len(lines) + SUBNORMAL_ROWS) % 2 == 0:
        print("total row count would be even", file=sys.stderr)
        return 1
    rng = random.Random(SUBNORMAL_SEED)
    pts = subnormal_points(rng, SUBNORMAL_ROWS, taken)
    rows = emit(pts)
    if len(rows) != SUBNORMAL_ROWS or not verify_appended(rows):
        return 1
    eol = b"\r\n" if b"\r\n" in raw else b"\n"
    out = b"" if raw.endswith(b"\n") else eol
    out += b"".join(f"{x.hex()} {y.hex()}".encode("ascii") + eol for x, y in rows)
    with open(path, "ab") as f:
        f.write(out)
    print(f"{path}: appended {len(rows)} subnormal-x rows "
          f"({len(lines)} -> {len(lines) + len(rows)})", file=sys.stderr)
    return 0


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--subnormal-append":
        return append_subnormal(sys.argv[2] if len(sys.argv) > 2
                                else "tests/data/erf_reference.txt")
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

    rows = emit(pts)
    if not self_check(rows):
        return 1
    for x, y in rows:
        print(f"{x.hex()} {y.hex()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
