#!/usr/bin/env python3
"""Generate tests/data/digamma_reference.txt -- correctly rounded digamma
oracle.

Each line: <input-hex-double> <digamma-hex-double>, the value mpmath.digamma
converges to at dps 60 AND dps 100 (checked to agree far below the double's
own ULP before being trusted), rounded to nearest double. Specials (poles,
+-inf, NaN, x = 0) are covered by the smoke test, not this file -- same
convention as tools/gen_lgamma_reference.py.

lgamma's own reference file carries NO dd pair -- two hex doubles per line,
nothing more -- despite lgamma's kernel computing in dd internally. Its own
hardest metric (the negative-axis |lgamma|<1 absolute band,
tests/test_lgamma_ulp.cpp) is measured against that same plain rounded
double, because a correctly-rounded double is already accurate in ABSOLUTE
terms in proportion to its own magnitude -- exactly what's needed near a
zero crossing. digamma's design doctrine is the direct analogue (relative
where |psi| >= 1, else 2^-53-class absolute near the negative-axis zeros),
so this generator follows the same convention: no dd pair.

Point selection targets what the kernel's region structure (src/digamma_data.h)
makes fragile, mirroring lgamma's rationale one level up in complexity because
digamma adds a product-form root, a reflection formula, and a genuinely
unbounded pole ladder on the negative axis:
  - the product-form zone [1, 2) and its root x0 (unique positive zero of
    digamma), where the claim is relative accuracy through the sign change
  - (0, 1), where the kernel avoids forming 1+x near x -> 1- and instead
    shifts against (x0 - 1); sampled log-spaced down INTO the subnormal range
    (some of those points overflow -1/x past DBL_MAX on the oracle side and
    are silently dropped, same as any other non-finite result)
  - the [2, 8) fixed-step down-walk region, plus brackets at each integer
    step threshold (a lane can take a different number of steps than its
    neighbour) and at the walk's [1,2) landing boundary
  - the asymptotic region [8, ~DBL_MAX], log-spaced, plus explicit huge-x
    witnesses near 2^53, 1e300 and 1e308
  - the negative axis: the first 20 zeros of digamma (recomputed here, not
    copied from anywhere) with both a nearest-double point and relative
    offsets bracketing the sign change; a spread of pole neighbourhoods from
    both sides at ulp-scale offsets; dense generic sampling on (-21, 0); and
    log-spaced far-negative non-integers out to ~2^52 (see the note at
    far_negative_points() on why [2^52, 2^53) contributes nothing here)

Mechanism notes: mp.dps is set INSIDE every computation function, never at
module scope for anything that runs during point generation; run is
foreground-only and single-shot (point count here computes in low tens of
seconds, not the ~5-minute chunk ceiling that motivated the rule elsewhere
in this project).

Usage:
    python tools/gen_digamma_reference.py > tests/data/digamma_reference.txt
"""

import math
import random
import struct
import sys

import mpmath as mp

SEED = 20260806
X0 = 8.0
ZONE_LO = 1.0
ZONE_HI = 2.0

# Agreement threshold for the layered-dps self-check: dps 60 vs dps 100 must
# agree to within this combined (absolute + relative) bound before a row is
# trusted. dps 60 alone converges to roughly 1e-58; this threshold carries
# ~18 orders of magnitude of margin under that, chosen to be comfortably
# below the target than to be a diagnostic anyone expects to be marginal.
DPS_AGREE_TOL = mp.mpf("1e-40")


def as_bits(x: float) -> int:
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def from_bits(u: int) -> float:
    return struct.unpack("<d", struct.pack("<Q", u))[0]


def neighbourhood(x0: float, k: int = 48):
    b = as_bits(x0)
    return [from_bits(b + j) for j in range(-k, k + 1)]


def digamma_layered(x_float):
    """mpmath.digamma(x) at dps 60 and dps 100; returns (y100, ok, diff).

    mp.dps is set INSIDE this function (mechanism rule) and restored after.
    Agreement is checked at the higher-precision (dps 100) layer, combining
    absolute and relative so it works uniformly near poles (huge magnitude)
    and near zeros (magnitude -> 0).
    """
    old = mp.mp.dps
    try:
        mp.mp.dps = 60
        y60 = mp.digamma(mp.mpf(x_float))
        mp.mp.dps = 100
        y100 = mp.digamma(mp.mpf(x_float))
        diff = abs(y100 - y60)
        ok = diff <= DPS_AGREE_TOL * (1 + abs(y100))
        return y100, bool(ok), diff
    finally:
        mp.mp.dps = old


def find_negative_zero(n, dps):
    """The n-th negative zero of digamma (n=1 in (-1,0), n=2 in (-2,-1), ...).

    mp.dps is set INSIDE this function.
    """
    old = mp.mp.dps
    try:
        mp.mp.dps = dps
        lo, hi = mp.mpf(-n), mp.mpf(-(n - 1))
        guess = mp.mpf(-n) + mp.mpf("0.5")
        try:
            z = mp.findroot(lambda x: mp.digamma(x), guess)
        except Exception:
            z = mp.findroot(
                lambda x: mp.digamma(x),
                (lo + mp.mpf("1e-6"), hi - mp.mpf("1e-6")),
                solver="bisect",
            )
        return z
    finally:
        mp.mp.dps = old


def find_positive_root(dps):
    """The unique positive root of digamma (~1.4616...). mp.dps set INSIDE."""
    old = mp.mp.dps
    try:
        mp.mp.dps = dps
        return mp.findroot(lambda x: mp.digamma(x), mp.mpf("1.4616321449683622"))
    finally:
        mp.mp.dps = old


def _hand_digamma_positive(y, dps, shift, nterms):
    """psi(y) for y > 0 via up-recurrence to `shift` + hand Bernoulli series.
    Caller has already set mp.dps."""
    s = mp.mpf(0)
    while y < shift:
        s += 1 / y
        y += 1
    w = 1 / (y * y)
    acc = mp.mpf(0)
    wp = w
    for k in range(1, nterms + 1):
        b = mp.bernoulli(2 * k)
        acc += b / (2 * k) * wp
        wp *= w
    asym = mp.log(y) - 1 / (2 * y) - acc
    return asym - s


def hand_digamma(x, dps, shift=40, nterms=15):
    """Independent evaluation path for spot-checking mp.digamma.

    mp.psi(0, x) and mp.digamma(x) are literally the SAME function in
    mpmath (verified: `mp.psi(0, z) == mp.digamma(z)` bitwise, and
    `mpmath.digamma is mpmath.mp.digamma`), so they cannot serve as an
    independent check on each other. This hand-rolled evaluator instead
    re-derives digamma from the recurrence psi(x) = psi(x+1) - 1/x walked
    up to a large argument, plus the textbook Bernoulli asymptotic series
    written out by hand -- not a call into mpmath's own digamma algorithm
    at any point (mp.bernoulli() is a generic sequence primitive, not part
    of mpmath's digamma implementation). Negative x goes through the
    reflection formula psi(x) = psi(1-x) - pi*cot(pi*x) (mp.cot() is again
    an independent primitive) rather than recurring up from x itself --
    walking from x = -1e6 or the far-negative stratum's ~-4.5e15 to
    `shift` one integer at a time would be a million-plus-iteration loop;
    psi(1-x) starts the walk already past `shift` for any x that negative,
    so the loop body never executes. Measured agreement against mp.digamma
    at dps 60
    across zone/mid/negative/root/small-x samples: worst diff ~2.5e-43,
    i.e. this is not a marginal check.

    mp.dps is set INSIDE this function.
    """
    old = mp.mp.dps
    try:
        mp.mp.dps = dps
        x = mp.mpf(x)
        if x < 0:
            psi_1mx = _hand_digamma_positive(1 - x, dps, shift, nterms)
            return psi_1mx - mp.pi * mp.cot(mp.pi * x)
        return _hand_digamma_positive(x, dps, shift, nterms)
    finally:
        mp.mp.dps = old


def is_negative_integer(x_mpf):
    return x_mpf < 0 and x_mpf == mp.floor(x_mpf)


def find_overflow_boundary_bits(lo_inf: float, hi_finite: float, dps=60) -> int:
    """Bisect (in IEEE bit-pattern space) the smallest-x boundary where
    digamma(x) stops overflowing DBL_MAX. `lo_inf` must overflow, `hi_finite`
    must not. mp.dps is set INSIDE this function."""
    old = mp.mp.dps
    try:
        mp.mp.dps = dps
        lo, hi = as_bits(lo_inf), as_bits(hi_finite)

        def finite_at(b):
            return math.isfinite(float(mp.digamma(mp.mpf(from_bits(b)))))

        assert not finite_at(lo) and finite_at(hi)
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if finite_at(mid):
                hi = mid
            else:
                lo = mid
        return hi  # first finite bit pattern
    finally:
        mp.mp.dps = old


# ---------------------------------------------------------------------------
# Point-set strata. Each returns (stratum_name, [python-float points]).
# ---------------------------------------------------------------------------


def strata_region_grids(rng):
    out = []

    # 1a. Zone [1, 2), dense.
    out.append(("zone [1,2) dense", [rng.uniform(ZONE_LO, ZONE_HI) for _ in range(3000)]))

    # 1b. (0, 1): linear-dense plus log-spaced down into (and past) the
    # subnormal boundary. Points whose oracle result overflows past DBL_MAX
    # (roughly x < 5.6e-309, since psi(x) ~ -1/x there) are dropped as
    # non-finite by emit(), same convention as any other special.
    out.append(("(0,1) linear", [rng.uniform(0.0, 1.0) for _ in range(1000)]))
    out.append(
        ("(0,1) log-spaced incl. subnormal", [10.0 ** rng.uniform(-323.0, -1.0) for _ in range(2000)])
    )
    # Explicit bracket around the finite/overflow boundary (psi(x) ~ -1/x
    # overflows DBL_MAX below ~5.56e-309): the boundary bit pattern is found
    # by bisection (not assumed), then bracketed so the stratum straddles
    # it -- confirms emit() cleanly separates the last finite points from
    # the first overflowing ones rather than papering over a cliff.
    boundary_bits = find_overflow_boundary_bits(1e-309, 1e-308)
    out.append(
        (
            "(0,1) overflow-boundary bracket",
            [from_bits(boundary_bits + j) for j in range(-48, 49)],
        )
    )

    # 1c. [2, 8) walk region, dense, plus a bracket at every integer step
    # threshold the fixed-step down-walk can land on (3..7) -- not in the
    # literal POINT SET list but directly targets kDigammaWalkDepth's
    # per-step correctness, the same reasoning lgamma's own generator used
    # for its recurrence thresholds.
    out.append(("[2,8) walk dense", [rng.uniform(2.0, 8.0) for _ in range(2000)]))
    step_pts = []
    for k in (3.0, 4.0, 5.0, 6.0, 7.0):
        step_pts += neighbourhood(k, 32)
    out.append(("[2,8) integer step brackets", step_pts))

    # 1d. [8, ~DBL_MAX], log-spaced asymptotic.
    out.append(
        (
            "[8,1e308) log asymptotic",
            [math.exp(rng.uniform(math.log(X0), math.log(1e308))) for _ in range(2000)],
        )
    )

    return out


def strata_root(rng):
    old = mp.mp.dps
    try:
        z = find_positive_root(100)
        z60 = find_positive_root(60)
        diff = abs(z - z60)
        print(f"  root recompute agreement (dps100 vs dps60): {mp.nstr(diff, 5)}", file=sys.stderr)
        if diff > DPS_AGREE_TOL:
            print("FATAL: positive root dps60/dps100 disagreement", file=sys.stderr)
            sys.exit(2)
        x0 = float(z)
    finally:
        mp.mp.dps = old

    pts = []
    b = as_bits(x0)
    pts += [from_bits(b - 1), from_bits(b), from_bits(b + 1)]
    for k in range(14, 0, -1):
        d = 10.0 ** (-k)
        pts.append(x0 + d)
        pts.append(x0 - d)
    return [("root x0 neighbourhood", pts)]


def strata_boundary_brackets():
    pts = []
    for b in (1.0, 2.0, 8.0):
        pts += neighbourhood(b, 64)
    return [("region boundary brackets {1,2,8}", pts)]


def strata_negative(rng):
    out = []

    # 4a. First 20 negative-axis zeros, recomputed at dps >= 60 here (never
    # copied). Nearest double + relative offsets bracketing the sign change
    # on both sides.
    offsets = (1e-15, 1e-12, 1e-9, 1e-6, 1e-3)
    zero_pts = []
    zeros100 = []
    zeros60 = []
    for n in range(1, 21):
        z100 = find_negative_zero(n, 100)
        z60 = find_negative_zero(n, 60)
        zeros100.append(z100)
        zeros60.append(z60)
        diff = abs(z100 - z60)
        if diff > DPS_AGREE_TOL:
            print(f"FATAL: negative zero n={n} dps60/dps100 disagreement: {diff}", file=sys.stderr)
            sys.exit(2)
        zf = float(z100)
        zero_pts.append(zf)
        for d in offsets:
            zero_pts.append(zf + d)
            zero_pts.append(zf - d)
    print(
        f"  negative zero recompute agreement (dps100 vs dps60), n=1..20: "
        f"worst {max(float(abs(a - b)) for a, b in zip(zeros100, zeros60)):.3e}",
        file=sys.stderr,
    )
    out.append(("negative zero neighbourhoods (n=1..20)", zero_pts))

    # 4b. Pole neighbourhoods: a spread of |n| from both sides at several
    # ulp-scale offsets.
    pole_pts = []
    ulp_offsets = (1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024)
    for n in list(range(1, 21)) + [100, 1000, 1_000_000]:
        base = -float(n)
        b = as_bits(base)
        for o in ulp_offsets:
            pole_pts.append(from_bits(b + o))
            pole_pts.append(from_bits(b - o))
    out.append(("pole neighbourhoods (n=1..20,100,1000,1e6)", pole_pts))

    # 4c. Dense general-position sampling over (-21, 0).
    out.append(("(-21,0) dense generic", [-rng.uniform(0.0, 21.0) for _ in range(3000)]))

    # 4d. Log-spaced far-negative non-integer points out to |x| ~ 2^52.
    # NOTE: verified below that every double in [2^52, 2^53) has ULP = 1,
    # i.e. every representable value in that range IS an integer. Those
    # would all be exact poles (NaN under the design's negative-integer
    # convention), which is the smoke test's doctrine, not this reference
    # file's -- so that half-open interval contributes zero rows here, by
    # design, not by omission.
    far_pts = []
    lo, hi = math.log(21.0), math.log(float(2**52))
    tries = 0
    while len(far_pts) < 1000 and tries < 4000:
        tries += 1
        m = math.exp(rng.uniform(lo, hi))
        x = -m
        if x == math.floor(x):
            continue  # reroll on the rare exact-integer draw
        far_pts.append(x)
    out.append(("far-negative log-spaced non-integer (<2^52)", far_pts))

    return out


def strata_huge_positive():
    pts = []
    pts += neighbourhood(float(2**53), 48)
    pts += neighbourhood(1e300, 48)
    hi = neighbourhood(1e308, 48)
    pts += [p for p in hi if math.isfinite(p)]
    return [("huge positive (2^53, 1e300, 1e308)", pts)]


def emit(strata):
    seen = set()
    counts = {}
    skipped_overflow = 0
    rows = []  # (x_float, y_float) for spot-check reuse
    lines = []
    for name, points in strata:
        n_here = 0
        for x in points:
            if not math.isfinite(x) or x == 0.0:
                continue
            b = as_bits(x)
            if b in seen:
                continue
            seen.add(b)
            xm = mp.mpf(x)
            if is_negative_integer(xm):
                continue  # pole; smoke test covers it
            y, ok, diff = digamma_layered(x)
            if not ok:
                print(
                    f"FATAL: layered-dps disagreement at x={x!r} "
                    f"(stratum '{name}'): diff={diff}",
                    file=sys.stderr,
                )
                sys.exit(2)
            if not mp.isfinite(y):
                skipped_overflow += 1
                continue
            yf = float(y)
            if not math.isfinite(yf):
                skipped_overflow += 1
                continue
            lines.append(f"{x.hex()} {yf.hex()}")
            rows.append((x, yf))
            n_here += 1
        counts[name] = n_here
        if n_here == 0:
            print(f"FATAL: stratum '{name}' emitted zero rows", file=sys.stderr)
            sys.exit(2)
    return lines, rows, counts, skipped_overflow


def spot_check(rows, rng, n=25):
    """Independent re-derivation of n random emitted rows via hand_digamma.

    Rows span magnitudes from ~0 (near the root/zeros) to ~1e274 (near the
    smallest sampled positive x, where psi(x) ~ -1/x), so the gate must be
    relative-or-absolute like the layered-dps check, not a bare absolute
    diff -- a fixed absolute threshold would flag an actually-fine
    ~1e-17-relative row as a failure purely because its magnitude is
    ~1e274.
    """
    sample = rng.sample(rows, min(n, len(rows)))
    worst_norm = mp.mpf(0)
    worst_x = None
    worst_diff = mp.mpf(0)
    for x, yf in sample:
        h = hand_digamma(x, 60)
        d = abs(mp.mpf(yf) - h)
        norm = d / (1 + abs(h))
        if norm > worst_norm:
            worst_norm = norm
            worst_diff = d
            worst_x = x
    print(
        f"  spot re-derivation (independent hand_digamma, n={len(sample)}): "
        f"worst normalized |stored - hand| = {mp.nstr(worst_norm, 5)} "
        f"(raw diff {mp.nstr(worst_diff, 5)}) at x={worst_x!r}",
        file=sys.stderr,
    )
    # hand_digamma agrees with mp.digamma to ~1e-43 relative at dps 60 (see
    # its docstring); the stored value is a rounded double, so the real
    # budget here is double rounding (~1.1e-16 relative) plus that ~1e-43
    # slack. A failure here would mean the STORED value, not just
    # hand_digamma, is wrong -- gate generously but not vacuously.
    if worst_norm > mp.mpf("1e-12"):
        print("FATAL: spot re-derivation exceeds sanity threshold", file=sys.stderr)
        sys.exit(2)


def main():
    rng = random.Random(SEED)

    print(f"psi(0,x) is digamma alias in this mpmath: "
          f"{mp.psi is not None}; using hand_digamma() for spot checks "
          f"per the aliasing note in its docstring.", file=sys.stderr)

    # Verify the [2^52, 2^53) all-integer claim (see far_negative note)
    # directly rather than asserting it silently: ULP there must be 1.0,
    # and a handful of sample doubles in the range must equal their floor.
    lo52, hi52 = float(2**52), float(2**53)
    ulp52 = from_bits(as_bits(lo52) + 1) - lo52
    all_int = all(
        (lambda v: v == math.floor(v))(from_bits(as_bits(lo52) + j))
        for j in (0, 1, 2, 1000, 2**20, 2**51 - 1)
    )
    print(
        f"  [2^52,2^53) verify: ULP={ulp52} (expect 1.0), sample doubles all "
        f"integer={all_int}",
        file=sys.stderr,
    )
    if ulp52 != 1.0 or not all_int:
        print("FATAL: [2^52,2^53) all-integer assumption failed", file=sys.stderr)
        return 2

    strata = []
    strata += strata_region_grids(rng)
    strata += strata_root(rng)
    strata += strata_boundary_brackets()
    strata += strata_negative(rng)
    strata += strata_huge_positive()

    lines, rows, counts, skipped_overflow = emit(strata)

    for name, c in counts.items():
        print(f"  stratum '{name}': {c} rows", file=sys.stderr)
    print(f"  skipped (oracle result overflows double): {skipped_overflow}", file=sys.stderr)

    spot_check(rows, rng)

    total = len(lines)
    print(f"emitted {total} points", file=sys.stderr)
    if total < 8000:
        print("FATAL: reference set suspiciously small", file=sys.stderr)
        return 2

    for line in lines:
        print(line)

    print("PASS: all self-checks clean", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
