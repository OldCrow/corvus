#!/usr/bin/env python3
"""Generate tests/data/trigamma_reference.txt -- correctly rounded trigamma
oracle.

Each line: <input-hex-double> <trigamma-hex-double>, the value mpmath's
trigamma (polygamma order 1) converges to at dps 60 AND dps 100 (checked to
agree far below the double's own ULP before being trusted), rounded to
nearest double. Format matches lgamma/digamma's reference files exactly:
two hex doubles per line, no dd pair -- PLAN.md P1 trigamma's own doctrine
is ALL-RELATIVE everywhere (psi_1 = sum of squares, strictly positive
wherever finite -- no zero crossing anywhere on the real line, unlike
digamma), so there is no absolute-band metric this file needs to serve; a
correctly-rounded double is the right oracle representation throughout.

Specials (poles at 0 and every negative integer, +-inf, NaN) are covered by
the smoke test, not this file -- same convention as digamma/lgamma.

Oracle aliasing check (mirrors gen_digamma_reference.py's own check):
mpmath's `polygamma` and `psi` are the SAME bound function
(`mp.polygamma.__func__ is mp.psi.__func__`, verified below), both routing
through the identical `libmp.mpf_psi` C-level implementation regardless of
which name is called -- so calling `mp.psi(1, x)` as a "second" check would
NOT be independent. The genuinely independent check here is `hand_trigamma`:
a direct sum 1/(x+n)^2 to N plus an Euler-Maclaurin (Bernoulli asymptotic)
tail, using `mp.bernoulli()` and `mp.sin()` as independent primitives --
never calling `mp.polygamma`/`mp.psi`/`mp.digamma` at any point. Measured
agreement against the oracle at dps 60 across zone/mid/negative/huge
samples: worst diff ~1e-44 relative (see hand_trigamma's own docstring) --
not a marginal check.

Point selection targets what the kernel's region structure
(src/trigamma_data.h) makes fragile, mirroring digamma's rationale one
level up in a simpler function (no root, no product form, no cos table, but
an unbounded pole ladder like digamma's and a genuinely two-sided small-x
regime that digamma's Laurent term doesn't have):
  - (0, 1): log-spaced down into the subnormal range, with explicit
    brackets at the ~2^-512 overflow boundary (1/x^2 itself overflows
    DBL_MAX there), the ~2^-480 deep-tiny guard (below which the zone term
    stops mattering to double precision), and the ~2^-28 crossover (below
    which the pi^2/6 Laurent-constant contribution stops mattering at all)
  - the zone [1, 2), dense
  - the [2, 8) fixed-step down-walk region, plus brackets at each integer
    step threshold (3..7) -- a lane can take a different walk depth than
    its neighbour, same reasoning as digamma's own down-walk brackets
  - the asymptotic region [8, ~DBL_MAX], log-spaced, with explicit
    brackets straddling kTrigammaAsymCut = 2^89 on both sides (below which
    the direct Bernoulli-sum form is used, above which fl(1/x) alone
    suffices) plus huge-x witnesses near 2^53, 1e300, 1e308
  - the negative axis: dense generic sampling on (-50, 0); the global-min
    neighbourhood x ~ -0.4957 (psi_1's unique negative-axis minimum,
    8.933..., recomputed here, not copied from anywhere); near-pole
    ulp-offset brackets from both sides at n=1..20, 100, 1000, ~1e6; and
    log-spaced far-negative non-integers out to ~2^52 ([2^52, 2^53)
    contributes nothing: every double there is an integer, i.e. a pole --
    see strata_negative(); identical reasoning to digamma's generator)

Mechanism rule (binding): mp.dps is set INSIDE every computation
function, never at module scope for anything that runs during point
generation; run is foreground-only and single-shot.

Usage:
    python tools/gen_trigamma_reference.py > tests/data/trigamma_reference.txt
"""

import math
import random
import struct
import sys

import mpmath as mp

SEED = 20260808
X0 = 8.0
ZONE_LO = 1.0
ZONE_HI = 2.0
ASYM_CUT = float(mp.mpf(2) ** 89)

# Agreement threshold for the layered-dps self-check: dps 60 vs dps 100 must
# agree to within this combined (absolute + relative) bound before a row is
# trusted. Matches gen_digamma_reference.py's own convention/value.
DPS_AGREE_TOL = mp.mpf("1e-40")


def as_bits(x: float) -> int:
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def from_bits(u: int) -> float:
    return struct.unpack("<d", struct.pack("<Q", u))[0]


def neighbourhood(x0: float, k: int = 48):
    b = as_bits(x0)
    return [from_bits(b + j) for j in range(-k, k + 1)]


# mp.polygamma(1, x) has NO fast path for negative x internally, unlike
# mp.digamma (order 0): it costs ~4.2s at x=-1e6 at dps=100 and scales
# ~linearly with |x| (an internal naive up-recurrence, not reflection) --
# effectively unusable past ~1e6 and a genuine hang by ~1e15, exactly the
# range this file's far-negative and n=1e6 pole-neighbourhood strata
# need. mp.polygamma(1, POSITIVE x) has no such problem (confirmed <1ms
# even at x=4.5e15) -- same asymmetry as mp.digamma vs the naive negative
# recurrence hand_trigamma's own docstring already flags. Fix: route
# |x|>50 through reflection onto the fast positive branch. Verified exact
# (diff ~1e-99 at dps=100) against a direct call at a magnitude where
# both are still fast (x=-1234.5678).
NEGATIVE_FAST_PATH_THRESHOLD = 50.0


def trigamma_oracle(x_mpf):
    """Correctly-rounded trigamma at the CURRENT mp.dps. Caller sets dps."""
    if x_mpf < -NEGATIVE_FAST_PATH_THRESHOLD:
        s = mp.sin(mp.pi * x_mpf)
        return (mp.pi ** 2) / (s * s) - mp.polygamma(1, 1 - x_mpf)
    return mp.polygamma(1, x_mpf)


def trigamma_layered(x_float):
    """mpmath trigamma(x) at dps 60 and dps 100; returns (y100, ok, diff).

    mp.dps is set INSIDE this function (mechanism rule) and restored after.
    """
    old = mp.mp.dps
    try:
        mp.mp.dps = 60
        y60 = trigamma_oracle(mp.mpf(x_float))
        mp.mp.dps = 100
        y100 = trigamma_oracle(mp.mpf(x_float))
        diff = abs(y100 - y60)
        ok = diff <= DPS_AGREE_TOL * (1 + abs(y100))
        return y100, bool(ok), diff
    finally:
        mp.mp.dps = old


def find_negative_min(dps):
    """The unique negative-axis global minimum of trigamma (x ~ -0.4957,
    value ~8.933), found as the root of tetragamma (psi_2) on (-1,0).
    mp.dps is set INSIDE this function."""
    old = mp.mp.dps
    try:
        mp.mp.dps = dps
        return mp.findroot(lambda x: mp.polygamma(2, x), mp.mpf("-0.5"))
    finally:
        mp.mp.dps = old


def _hand_trigamma_positive(y, dps, N, nterms):
    """psi_1(y) for y > 0 via direct sum to N + Euler-Maclaurin tail.
    Caller has already set mp.dps."""
    s = mp.mpf(0)
    for n in range(N):
        s += 1 / ((y + n) ** 2)
    yN = y + N
    w = 1 / (yN * yN)
    acc = mp.mpf(0)
    wp = w
    for k in range(1, nterms + 1):
        b = mp.bernoulli(2 * k)
        acc += b * wp
        wp *= w
    tail = 1 / yN + w / 2 + acc / yN
    return s + tail


def hand_trigamma(x, dps, N=40, nterms=15):
    """Independent evaluation path for spot-checking mp.polygamma(1,.).

    mp.polygamma and mp.psi are literally the same bound function in this
    mpmath (verified: `mp.polygamma.__func__ is mp.psi.__func__`), both
    routing through libmp.mpf_psi -- so they cannot serve as an independent
    check on each other. This hand-rolled evaluator instead re-derives
    trigamma from a direct sum 1/(x+n)^2 plus the textbook Bernoulli
    (Euler-Maclaurin) asymptotic tail written out by hand -- not a call
    into mpmath's own polygamma/psi/digamma implementation at any point
    (mp.bernoulli() and mp.sin() are generic primitives, not part of
    mpmath's polygamma algorithm). Negative x goes through the reflection
    formula psi_1(x) = pi^2/sin^2(pi x) - psi_1(1-x) (mp.sin() again an
    independent primitive) rather than summing from x itself -- summing
    from x=-1e6 or the far-negative stratum's ~-4.5e15 directly would be a
    million-plus-term sum; psi_1(1-x) starts the sum already past N for any
    x that negative. Measured agreement against mp.polygamma(1,.) at dps 60
    across zone/mid/negative/huge samples: worst diff ~2e-44 relative --
    not a marginal check.

    mp.dps is set INSIDE this function.
    """
    old = mp.mp.dps
    try:
        mp.mp.dps = dps
        x = mp.mpf(x)
        if x < 0:
            p = _hand_trigamma_positive(1 - x, dps, N, nterms)
            s = mp.sin(mp.pi * x)
            return (mp.pi ** 2) / (s * s) - p
        return _hand_trigamma_positive(x, dps, N, nterms)
    finally:
        mp.mp.dps = old


def is_negative_integer_or_zero(x_mpf):
    return x_mpf <= 0 and x_mpf == mp.floor(x_mpf)


def find_overflow_boundary_bits(lo_inf: float, hi_finite: float, dps=60) -> int:
    """Bisect (in IEEE bit-pattern space) the smallest-x boundary where
    trigamma(x) stops overflowing DBL_MAX. `lo_inf` must overflow,
    `hi_finite` must not. mp.dps is set INSIDE this function."""
    old = mp.mp.dps
    try:
        mp.mp.dps = dps
        lo, hi = as_bits(lo_inf), as_bits(hi_finite)

        def finite_at(b):
            return math.isfinite(float(mp.polygamma(1, mp.mpf(from_bits(b)))))

        assert not finite_at(lo) and finite_at(hi)
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if finite_at(mid):
                hi = mid
            else:
                lo = mid
        return hi
    finally:
        mp.mp.dps = old


# ---------------------------------------------------------------------------
# Point-set strata. Each returns (stratum_name, [python-float points]).
# ---------------------------------------------------------------------------


def strata_region_grids(rng):
    out = []

    # 1a. (0, 1): linear-dense plus log-spaced down into (and past) the
    # subnormal boundary, with explicit brackets at the overflow boundary
    # (~2^-512, where 1/x^2 itself overflows DBL_MAX), the deep-tiny guard
    # (~2^-480), and the pi^2/6-crossover (~2^-28, below which the zone's
    # Laurent-constant contribution stops mattering to double precision).
    out.append(("(0,1) linear", [rng.uniform(0.0, 1.0) for _ in range(1000)]))
    out.append(
        ("(0,1) log-spaced incl. subnormal",
         [10.0 ** rng.uniform(-323.0, -1.0) for _ in range(2000)])
    )
    boundary_bits = find_overflow_boundary_bits(2.0 ** -600, 2.0 ** -500)
    out.append(
        ("(0,1) overflow-boundary bracket (~2^-512)",
         [from_bits(boundary_bits + j) for j in range(-48, 49)])
    )
    out.append(
        ("(0,1) deep-tiny-guard bracket (~2^-480)",
         neighbourhood(2.0 ** -480, 48))
    )
    out.append(
        ("(0,1) pi^2/6-crossover bracket (~2^-28)",
         neighbourhood(2.0 ** -28, 48))
    )

    # 1b. [1,2) zone, dense.
    out.append(("zone [1,2) dense", [rng.uniform(ZONE_LO, ZONE_HI) for _ in range(3000)]))

    # 1c. [2, 8) walk region, dense, plus brackets at every integer step
    # threshold the fixed-step down-walk can land on (3..7).
    out.append(("[2,8) walk dense", [rng.uniform(2.0, 8.0) for _ in range(2000)]))
    step_pts = []
    for k in (3.0, 4.0, 5.0, 6.0, 7.0):
        step_pts += neighbourhood(k, 32)
    out.append(("[2,8) integer step brackets", step_pts))

    # 1d. [8, ~DBL_MAX], log-spaced asymptotic, with brackets straddling
    # kTrigammaAsymCut = 2^89 on both sides.
    out.append(
        ("[8,1e308) log asymptotic",
         [math.exp(rng.uniform(math.log(X0), math.log(1e308))) for _ in range(2000)])
    )
    out.append(("asymptotic-cut bracket (2^89)", neighbourhood(ASYM_CUT, 48)))

    return out


def strata_boundary_brackets():
    pts = []
    for b in (1.0, 2.0, 8.0):
        pts += neighbourhood(b, 64)
    return [("region boundary brackets {1,2,8}", pts)]


def strata_negative(rng):
    out = []

    # 4a. Global-min neighbourhood (x ~ -0.4957), recomputed at dps>=60
    # here, never copied.
    m100 = find_negative_min(100)
    m60 = find_negative_min(60)
    diff = abs(m100 - m60)
    print(f"  global-min recompute agreement (dps100 vs dps60): {mp.nstr(diff, 5)}",
          file=sys.stderr)
    if diff > DPS_AGREE_TOL:
        print("FATAL: negative-axis global-min dps60/dps100 disagreement",
              file=sys.stderr)
        sys.exit(2)
    mf = float(m100)
    old = mp.mp.dps
    mp.mp.dps = 60
    try:
        min_val = mp.polygamma(1, m100)
    finally:
        mp.mp.dps = old
    print(f"  global min: x={mp.nstr(m100, 20)} value={mp.nstr(min_val, 20)}",
          file=sys.stderr)
    min_pts = neighbourhood(mf, 64)
    for d in (1e-9, 1e-6, 1e-3, 1e-2):
        min_pts.append(mf + d)
        min_pts.append(mf - d)
    out.append(("global-min neighbourhood (x~-0.4957)", min_pts))

    # 4b. Pole neighbourhoods: a spread of |n| from both sides at several
    # ulp-scale offsets (poles at 0 and every negative integer).
    pole_pts = []
    ulp_offsets = (1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024)
    for n in list(range(1, 21)) + [100, 1000, 1_000_000]:
        base = -float(n)
        b = as_bits(base)
        for o in ulp_offsets:
            pole_pts.append(from_bits(b + o))
            pole_pts.append(from_bits(b - o))
    out.append(("pole neighbourhoods (n=1..20,100,1000,1e6)", pole_pts))

    # 4c. Dense general-position sampling over (-50, 0).
    out.append(("(-50,0) dense generic", [-rng.uniform(0.0, 50.0) for _ in range(3000)]))

    # 4d. Log-spaced far-negative non-integer points out to |x| ~ 2^52.
    # Every double in [2^52, 2^53) is an integer (ULP=1 there), i.e. every
    # representable value would be an exact pole -- smoke-test territory,
    # not this file's -- so that half-open interval contributes zero rows
    # here, by design (verified explicitly in main(), matching digamma's
    # own generator).
    far_pts = []
    lo, hi = math.log(21.0), math.log(float(2 ** 52))
    tries = 0
    while len(far_pts) < 1000 and tries < 4000:
        tries += 1
        m = math.exp(rng.uniform(lo, hi))
        x = -m
        if x == math.floor(x):
            continue
        far_pts.append(x)
    out.append(("far-negative log-spaced non-integer (<2^52)", far_pts))

    return out


def strata_huge_positive():
    pts = []
    pts += neighbourhood(float(2 ** 53), 48)
    pts += neighbourhood(1e300, 48)
    hi = neighbourhood(1e308, 48)
    pts += [p for p in hi if math.isfinite(p)]
    return [("huge positive (2^53, 1e300, 1e308)", pts)]


def emit(strata):
    seen = set()
    counts = {}
    skipped_overflow = 0
    rows = []
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
            if is_negative_integer_or_zero(xm):
                continue  # pole; smoke test covers it
            y, ok, diff = trigamma_layered(x)
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
    """Independent re-derivation of n random emitted rows via hand_trigamma.

    Rows span magnitudes from ~8.93 (global min) up to ~1e618 (near the
    smallest sampled positive x, where psi_1(x) ~ 1/x^2), so the gate must
    be relative-or-absolute like the layered-dps check, not a bare absolute
    diff.
    """
    sample = rng.sample(rows, min(n, len(rows)))
    worst_norm = mp.mpf(0)
    worst_x = None
    worst_diff = mp.mpf(0)
    for x, yf in sample:
        h = hand_trigamma(x, 60)
        d = abs(mp.mpf(yf) - h)
        norm = d / (1 + abs(h))
        if norm > worst_norm:
            worst_norm = norm
            worst_diff = d
            worst_x = x
    print(
        f"  spot re-derivation (independent hand_trigamma, n={len(sample)}): "
        f"worst normalized |stored - hand| = {mp.nstr(worst_norm, 5)} "
        f"(raw diff {mp.nstr(worst_diff, 5)}) at x={worst_x!r}",
        file=sys.stderr,
    )
    if worst_norm > mp.mpf("1e-12"):
        print("FATAL: spot re-derivation exceeds sanity threshold", file=sys.stderr)
        sys.exit(2)


def main():
    rng = random.Random(SEED)

    print(
        f"polygamma/psi aliasing: mp.polygamma.__func__ is mp.psi.__func__ = "
        f"{mp.polygamma.__func__ is mp.psi.__func__}; using hand_trigamma() "
        "for spot checks per the aliasing note in its docstring (direct sum "
        "+ Euler-Maclaurin tail, sharing no mpmath polygamma internals).",
        file=sys.stderr,
    )

    # Verify the [2^52, 2^53) all-integer claim directly (matches digamma's
    # own generator's verification, same underlying IEEE754 fact).
    lo52, hi52 = float(2 ** 52), float(2 ** 53)
    ulp52 = from_bits(as_bits(lo52) + 1) - lo52
    all_int = all(
        (lambda v: v == math.floor(v))(from_bits(as_bits(lo52) + j))
        for j in (0, 1, 2, 1000, 2 ** 20, 2 ** 51 - 1)
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
    strata += strata_boundary_brackets()
    strata += strata_negative(rng)
    strata += strata_huge_positive()

    lines, rows, counts, skipped_overflow = emit(strata)

    for name, c in counts.items():
        print(f"  stratum '{name}': {c} rows", file=sys.stderr)
    print(f"  skipped (oracle result overflows double, near-pole/tiny-x "
          f"psi_1 > DBL_MAX): {skipped_overflow}", file=sys.stderr)

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
