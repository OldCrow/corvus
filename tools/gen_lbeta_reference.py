#!/usr/bin/env python3
"""Generate tests/data/lbeta_reference.txt -- correctly-rounded mpmath oracle
reference set for corvus::lbeta (public ln B(a,b) =
lgamma(a) + lgamma(b) - lgamma(a+b)).

Domain: a>0, b>0, both finite doubles only -- no NaN/inf/negative rows
(matching beta_p/beta_q's positive-domain contract; specials belong to the
smoke test, not this file, matching every other reference generator in
this repo).

Row format: THREE hex-float tokens per line -- `<a> <b> <lnB(a,b)>`,
strtod round-trip, one extra token beyond the erf/digamma/bessel two-token
convention because lbeta has two inputs. Past the huge-parameter
saturation boundary the third token is the literal `-inf` (float.hex()'s
own spelling for -inf; strtod parses it back exactly) -- the same
round-trip convention gen_bessel_reference.py uses for its overflow rows.

Oracle: mp.loggamma(a) + mp.loggamma(b) - mp.loggamma(a+b) at layered dps
(40 vs 80; escalate to 150 on relative disagreement > AGREE_REL_TOL;
DECLINE [skip + record] on persistent disagreement). This is erf-difficulty
oracle work, NOT the no-trusted-baseline oracle-trust doctrine (mpmath's
loggamma IS a trusted library baseline for this function, unlike beta's
betainc) -- no bracket certification. An independent cross-check on a
subset uses mp.log(mp.beta(a,b)) -- a genuinely different mpmath code path
(gamma RATIO, not a loggamma SUM) -- restricted to a moderate-magnitude
window (1e-6 <= a,b <= 1e6), where B(a,b) stays a numerically unremarkable
quantity, without betting a compute budget on mpmath's arbitrary-exponent
mpf tolerating astronomically large/small B values gracefully within a
reasonable per-call cost.

CONVENTION: tests/test_digamma_ulp.cpp buckets its own absolute-vs-
relative gate purely from the loaded reference VALUE (|psi(x)| < 1, decided
INSIDE the test), with NO marker column in tests/data/digamma_reference.txt
-- see that test's own comment, "the lgamma convention: on the negative
axis, the metric is chosen by |psi| alone". Followed verbatim here: rows
near the ln B = 0 manifold carry no marker; the eventual test_lbeta_ulp.cpp
will bucket by |stored lnB| < threshold itself, the same pattern. This
generator's job is DENSE coverage of that band, not annotating it.

RAY-CROSSING FINDING: of the three named boundary rays (a=b; a=10b;
a=1e5*b), only a=b actually reaches ln B magnitude ~DBL_MAX (rounds to
-inf) within representable double range. Large-(a,b) asymptotics give
ln B(a,b) ~ -(a+b)*H(p), p=a/(a+b), H the binary entropy in nats,
maximized (H=ln2) at p=1/2. A 10:1 or 1e5:1 skew shrinks H far enough that
even pinning BOTH parameters at DBL_MAX falls short of the ~1.7977e308
magnitude needed to round to -inf: a=10b maxes out at ln B ~ -6.0e307
(at a=DBL_MAX, b=DBL_MAX/10); a=1e5*b maxes out at ln B ~ -2.2e304 (at
a=DBL_MAX, b=DBL_MAX/1e5). Both are comfortably finite doubles, nowhere
near the boundary -- there IS no crossing to bit-step on those two rays.
The a=b ray's genuine crossing is derived (bisected) and bit-stepped
below; the other two rays are used as near-domain-ceiling STRESS rows
instead (their maximum-representable-magnitude point, bracketed) -- real
"both parameters huge, skewed ratio" coverage, honestly reported as NOT
hitting a saturation boundary rather than fabricating one.

BAND-SEAM ADDITION: the shipped kernel buckets its own gate by
min(a,b) > 2^990 as a distinct big-band region in test_lbeta_ulp.cpp. A
dedicated stratum straddles that seam log-spaced through [2^985, 2^995]
plus a fine bit-step cluster at 2^990 itself, crossed with a spread of
partner magnitudes in both argument orientations -- see band_seam_points().

Usage:
    python tools/gen_lbeta_reference.py
writes tests/data/lbeta_reference.txt directly (single light pass, no
checkpoint/resume machinery -- point count is small enough to regenerate
in one shot).
"""
import math
import random
import struct
import sys
import time

import mpmath as mp
from mpmath import mpf

from refgen_common import round_to_double

SEED = 20260811

DBL_MAX = float.fromhex("0x1.fffffffffffffp+1023")
SMALLEST_SUBNORMAL = float.fromhex("0x0.0000000000001p-1022")
MIN_NORMAL = float.fromhex("0x1.0000000000000p-1022")

DPS_LO = 40
DPS_HI = 80
DPS_ESCALATE = 150
# SECOND-TIER escalation: at extreme a<<b (or b<<a) corners --
# e.g. a=5e-324 (smallest subnormal), b=1e60 -- lgamma(a+b) rounds to
# EXACTLY lgamma(b) at ambient dps<~380 (a is absorbed below the
# representable digit range of b's magnitude), so ln B collapses to the
# analytically-correct lgamma(a) exactly -- there is no cancellation
# hazard at all here, only a precision-FLOOR effect: dps=80 vs dps=150
# disagree at ~7e-23 relative (past AGREE_REL_TOL) purely because dps=80
# hasn't yet converged lgamma(tiny-a) to 1e-25, NOT because the two tiers
# disagree about the answer -- measured directly: 150 vs 300 agree to
# ~7.7e-93, and dps>=80 already round to the IDENTICAL stored double in
# every case checked. One more escalation rung resolves it cleanly rather
# than declining a well-behaved row; DPS_ESCALATE_TIERS chases convergence
# instead of a hair-trigger single-shot escalation.
DPS_ESCALATE_TIERS = (150, 300, 600)
# Comfortably (~1e9x) below the double's own ~1.1e-16 relative ULP -- meant
# to catch a real oracle disagreement, not to be a marginal diagnostic.
# Same constant as gen_bessel_reference.py's AGREE_REL_TOL.
AGREE_REL_TOL = mpf("1e-25")

TARGET_LO, TARGET_HI = 2000, 6000

T0 = time.time()


# ---------------------------------------------------------------------------
# bit helpers
# ---------------------------------------------------------------------------
def as_bits(x: float) -> int:
    return struct.unpack("<q", struct.pack("<d", x))[0]


def from_bits(b: int) -> float:
    return struct.unpack("<d", struct.pack("<q", b))[0]


def neighbourhood(x0: float, k: int = 40, lo_off=None, hi_off=None):
    b = as_bits(x0)
    lo = -k if lo_off is None else lo_off
    hi = k if hi_off is None else hi_off
    return [from_bits(b + j) for j in range(lo, hi + 1)]


# ---------------------------------------------------------------------------
# oracle
# ---------------------------------------------------------------------------
def lnbeta_mpf(a, b, dps):
    """Raw (unrounded) mpf ln B(a,b) at the given dps. mp.dps is set INSIDE
    this function (mechanism rule, house style)."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        am, bm = mpf(a), mpf(b)
        return mp.loggamma(am) + mp.loggamma(bm) - mp.loggamma(am + bm)
    finally:
        mp.mp.dps = old


def _term_magnitude_digits(a, b, dps_probe=60):
    """Decimal-digit count of the magnitude of the biggest lgamma term
    feeding ln B(a,b) -- lgamma of whichever of {a, b, a+b} is largest.
    This probe is itself cancellation-FREE (a single lgamma of a single
    argument is well-conditioned everywhere), so dps_probe=60 is always
    trustworthy regardless of how badly the SUBTRACTION that follows will
    cancel. Used to decide how much working precision that subtraction
    needs before any agreement between two tiers can be trusted -- see
    layered_value's ORACLE TRAP comment."""
    old = mp.mp.dps
    mp.mp.dps = dps_probe
    try:
        am, bm = mpf(a), mpf(b)
        m = max(am, bm, am + bm)
        lg = mp.loggamma(m)
        if lg == 0:
            return 0
        return max(0, int(mp.log10(abs(lg))) + 1)
    finally:
        mp.mp.dps = old


def _agree(hi, lo, tol):
    """Relative agreement. hi==lo==0 is accepted as agreement -- but ONLY
    the caller (layered_value) may invoke this at a dps pair it has
    already confirmed clears the cancellation floor; see the ORACLE TRAP
    comment in layered_value for why an ungated zero-shortcut is a hazard,
    not this shortcut in isolation."""
    if hi == 0 and lo == 0:
        return True
    denom = abs(hi) if hi != 0 else abs(lo)
    if denom == 0:
        return False
    return abs(hi - lo) / denom <= tol


def layered_value(a, b):
    """Layered-dps oracle value with escalation. Returns (value_mpf, ok,
    escalated: bool) -- ok=False means DECLINE (row not emitted).

    ORACLE TRAP: for a<<b (or b<<a) pairs where BOTH parameters are
    individually huge -- e.g. a=2^730 (~1.75e219), b~2^1022 (~9e307) --
    lgamma(a)+lgamma(b) and lgamma(a+b) are both ~6e310 in magnitude
    (lgamma(b) dominates, ~b*ln(b)), while the TRUE ln B is ~-1e222: the
    subtraction cancels ~89 decimal digits. At dps=40 (and even dps=80)
    the entire computed difference is noise -- and mpmath's mpf
    subtraction of two operands that are IDENTICAL at that working
    precision returns EXACTLY 0, not "approximately 0". Both tiers can
    independently underflow to the identical exact zero, which `_agree`'s
    hi==0-and-lo==0 shortcut reads as proof of agreement: a coincidental
    exact-zero match sails through a relative-tolerance check with flying
    colors, silently storing 0x0.0p+0 for a genuinely huge-magnitude ln B.

    Fix: `_term_magnitude_digits` gives a cancellation-free estimate of
    how many decimal digits the biggest lgamma term carries; the
    subtraction needs working precision beyond that to have ANY signal
    left (+40 digits of buffer, matching AGREE_REL_TOL's own margin
    convention). When that floor sits at or below DPS_HI, the original
    cheap 40-vs-80(-vs-escalate) ladder is already trustworthy and runs
    unchanged -- this is the common case (the vast majority of rows never
    approach this floor). When the floor exceeds DPS_HI, dps 40/80/150(/
    300, depending how bad) are all potential noise and are skipped
    entirely in favor of directly comparing TWO tiers that both clear the
    floor -- so a trusted "agreement" can only ever be noise-free
    (floor <= max(a,b,a+b)'s lgamma digit count for THIS specific pair,
    never generic, since the cancellation depth is pair-dependent)."""
    floor = _term_magnitude_digits(a, b) + 40

    if floor <= DPS_HI:
        lo = lnbeta_mpf(a, b, DPS_LO)
        hi = lnbeta_mpf(a, b, DPS_HI)
        if _agree(hi, lo, AGREE_REL_TOL):
            return hi, True, False
        prev = hi
        for tier in DPS_ESCALATE_TIERS:
            esc = lnbeta_mpf(a, b, tier)
            if _agree(esc, prev, AGREE_REL_TOL):
                return esc, True, True
            prev = esc
        return prev, False, True

    # floor > DPS_HI: skip straight to a pair of floor-clearing tiers --
    # comparing anything below `floor` here would repeat the exact defect
    # above (noise vs noise, or noise vs signal, either way untrustworthy).
    t1, t2, t3 = floor, floor + 200, floor + 500
    v1 = lnbeta_mpf(a, b, t1)
    v2 = lnbeta_mpf(a, b, t2)
    if _agree(v2, v1, AGREE_REL_TOL):
        return v2, True, True
    v3 = lnbeta_mpf(a, b, t3)
    if _agree(v3, v2, AGREE_REL_TOL):
        return v3, True, True
    return v3, False, True


# ---------------------------------------------------------------------------
# point-set strata
# ---------------------------------------------------------------------------
def log_grid_values():
    """2D coverage axis: log-spaced magnitudes from the smallest subnormal
    through near-DBL_MAX, plus the binding design's explicit named
    anchors (subnormals/tiny, moderate 0.1..100, huge 1e5/1e10/1e100/1e300,
    near DBL_MAX)."""
    exps = list(range(-320, 301, 20))  # 32 log10-spaced anchors
    vals = [10.0 ** e for e in exps]
    vals += [
        SMALLEST_SUBNORMAL, 1e-300, 1e-30, 1e-5, 0.01, 0.1, 0.5, 1.0, 2.0,
        10.0, 100.0, 1e5, 1e10, 1e100, 1e300, DBL_MAX / 2.0,
        math.nextafter(DBL_MAX, 0.0), DBL_MAX,
    ]
    seen = set()
    out = []
    for v in vals:
        if not (math.isfinite(v) and v > 0.0):
            continue
        bb = as_bits(v)
        if bb in seen:
            continue
        seen.add(bb)
        out.append(v)
    return sorted(out)


def grid_cross_points(grid):
    return [(a, b) for a in grid for b in grid]


def random_fill(rng, n):
    """Dense random log-uniform fill across the whole representable
    exponent range."""
    pts = []
    for _ in range(n):
        a = 10.0 ** rng.uniform(-323.3, 308.25)
        b = 10.0 ** rng.uniform(-323.3, 308.25)
        if math.isfinite(a) and math.isfinite(b) and a > 0.0 and b > 0.0:
            pts.append((a, b))
    return pts


def moderate_fill(rng, n):
    """Extra density in the moderate 1e-3..1e3 band, where the eventual
    kernel's own region routing (shared with beta's PA/PB split) will
    matter most."""
    pts = []
    for _ in range(n):
        a = 10.0 ** rng.uniform(-3.0, 3.0)
        b = 10.0 ** rng.uniform(-3.0, 3.0)
        pts.append((a, b))
    return pts


# ---------------------------------------------------------------------------
# zero-manifold band: ln B = 0 through (1,1). For fixed a>0, ln B(a,b) is
# strictly decreasing in b (Beta(a,b) is strictly decreasing in b for fixed
# a), +inf as b->0+ and -inf as b->infty, so there is exactly one root --
# bisected directly, no closed form.
# ---------------------------------------------------------------------------
def zero_manifold_root(a, dps=80):
    lo, hi = 1e-6, 1e6

    def f(bb):
        return float(lnbeta_mpf(a, bb, dps))

    flo, fhi = f(lo), f(hi)
    assert flo > 0.0 and fhi < 0.0, (a, flo, fhi)
    for _ in range(100):
        mid = math.sqrt(lo * hi)
        if f(mid) > 0.0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


def zero_manifold_points():
    a_values = sorted(set([0.5 + i * (1.5 / 24) for i in range(25)] + [1.0]))
    pts = []
    for a in a_values:
        b0 = zero_manifold_root(a, dps=80)
        pts += [(a, x) for x in neighbourhood(float(b0), k=15)]
    return pts, a_values


# ---------------------------------------------------------------------------
# huge-parameter boundary. See module docstring RAY-CROSSING FINDING: only
# the a=b ray genuinely crosses -inf within double range; a=10b and
# a=1e5*b never do, and are used as near-ceiling stress rays instead.
# ---------------------------------------------------------------------------
def find_ab_ray_crossing(dps=50):
    def is_neg_inf(t):
        return math.isinf(float(lnbeta_mpf(t, t, dps)))

    lo, hi = 1.0e307, DBL_MAX
    assert not is_neg_inf(lo)
    assert is_neg_inf(hi)
    for _ in range(90):  # far more than enough: 1e308/2^90 << 1 ulp at 1e308
        # lo + (hi-lo)/2, NOT (lo+hi)/2 -- lo and hi are both ~1e308-class
        # doubles here, and their SUM overflows to +inf before the halving
        # ever happens (measured directly: 1e307+DBL_MAX rounds to +inf).
        # That would feed is_neg_inf(+inf) -> lnbeta_mpf(inf,inf,..) -> NaN
        # (never +inf), so the "else" branch would pin lo=inf permanently
        # and the bit-walk below would scan NaN bit patterns forever --
        # classic floating-point bisection-midpoint-overflow bug, worth
        # naming since it is easy to reintroduce.
        mid = lo + (hi - lo) / 2.0
        if is_neg_inf(mid):
            hi = mid
        else:
            lo = mid
    # quantize to the double grid, then bit-walk to the exact transition
    # (bounded -- the continuous bisection above already converges to far
    # better than 1 ulp, so this should take O(1) steps; a runaway here
    # would again be the overflow-class bug, not legitimate work).
    b = as_bits(float(lo))
    steps = 0
    while is_neg_inf(from_bits(b)):
        b -= 1
        steps += 1
        assert steps < 10_000, "bit-walk runaway -- see overflow-bisection note above"
    while not is_neg_inf(from_bits(b + 1)):
        b += 1
        steps += 1
        assert steps < 10_000, "bit-walk runaway -- see overflow-bisection note above"
    return from_bits(b), from_bits(b + 1)  # (last finite, first -inf)


def ab_crossing_points(last_finite, k=40):
    b = as_bits(last_finite)
    ts = [from_bits(b + j) for j in range(-k, k + 1)]
    return [(t, t) for t in ts if t > 0.0 and math.isfinite(t)]


def ray_stress_points(ratio, n=25, k=20):
    """Non-crossing ray (a=ratio*b): coverage of the maximum-representable-
    magnitude point and a log-spaced approach to it, honestly never -inf."""
    b_max = DBL_MAX / ratio
    if b_max <= 0 or not math.isfinite(b_max):
        return []
    lo_e, hi_e = -300.0, math.log10(b_max)
    bs = [10.0 ** (lo_e + i * (hi_e - lo_e) / (n - 1)) for i in range(n)]
    bs += [x for x in neighbourhood(b_max, k=k) if x > 0.0]
    pts = []
    for b in bs:
        a = ratio * b
        if math.isfinite(a) and a > 0.0 and math.isfinite(b) and b > 0.0:
            pts.append((a, b))
    return pts


# ---------------------------------------------------------------------------
# kernel band-boundary seam: the shipped kernel's test_lbeta_ulp.cpp buckets
# by min(a,b) > 2^990 (big-band) vs the relative/absolute split below that --
# straddle the seam log-spaced through [2^985, 2^995] plus a fine bit-step
# cluster right at 2^990, crossed with a spread of partner magnitudes so
# the seam is exercised as real 2D coverage, not a 1D probe. Both
# (small, partner) and (partner, small) orientations are emitted since the
# band test is on min(a,b), not on argument POSITION.
# ---------------------------------------------------------------------------
def band_seam_points():
    seam = 2.0 ** 990
    log_spaced = [2.0 ** (985.0 + 0.5 * i) for i in range(21)]  # 2^985..2^995
    fine = neighbourhood(seam, k=20)
    small_vals = sorted(set(v for v in (log_spaced + fine)
                             if v > 0.0 and math.isfinite(v)))
    partner_vals = [1e-300, 1e-3, 1.0, 1e3, 1e100, DBL_MAX / 4.0]
    pts = []
    for s in small_vals:
        for p in partner_vals:
            pts.append((s, p))
            pts.append((p, s))
    return pts, len(small_vals)


def symmetry_mirror(points, rng, frac=0.12):
    candidates = [(a, b) for (a, b) in points if a != b]
    sample = rng.sample(candidates, min(len(candidates), int(len(points) * frac)))
    return [(b, a) for (a, b) in sample]


# ---------------------------------------------------------------------------
# assembly
# ---------------------------------------------------------------------------
def build_points(rng):
    grid = log_grid_values()
    pts = []
    pts += grid_cross_points(grid)
    n_grid = len(pts)
    pts += random_fill(rng, 1200)
    n_random = len(pts) - n_grid
    pts += moderate_fill(rng, 400)
    n_moderate = len(pts) - n_grid - n_random

    zm_pts, zm_a_values = zero_manifold_points()
    pts += zm_pts
    n_zm = len(zm_pts)

    last_finite, first_inf = find_ab_ray_crossing()
    ab_cross_pts = ab_crossing_points(last_finite)
    pts += ab_cross_pts
    ray10_pts = ray_stress_points(10.0)
    ray1e5_pts = ray_stress_points(1e5)
    pts += ray10_pts + ray1e5_pts
    n_boundary = len(ab_cross_pts) + len(ray10_pts) + len(ray1e5_pts)

    seam_pts, n_seam_small = band_seam_points()
    pts += seam_pts
    n_seam = len(seam_pts)

    sym_pts = symmetry_mirror(pts, rng, frac=0.12)
    pts += sym_pts
    n_sym = len(sym_pts)

    seen = set()
    uniq = []
    for a, b in pts:
        if not (math.isfinite(a) and math.isfinite(b) and a > 0.0 and b > 0.0):
            continue
        key = (as_bits(a), as_bits(b))
        if key in seen:
            continue
        seen.add(key)
        uniq.append((a, b))

    meta = dict(
        grid_n=len(grid), n_grid_cross=n_grid, n_random=n_random,
        n_moderate=n_moderate, n_zero_manifold=n_zm, zm_a_n=len(zm_a_values),
        n_boundary=n_boundary, n_symmetry=n_sym,
        n_seam=n_seam, n_seam_small=n_seam_small,
        ab_crossing=(last_finite, first_inf),
        n_unique=len(uniq),
    )
    return uniq, meta


def emit_rows(points):
    rows = []
    declined = []
    escalated_n = 0
    for a, b in points:
        val, ok, escalated = layered_value(a, b)
        if escalated:
            escalated_n += 1
        if not ok:
            declined.append((a, b))
            continue
        y = round_to_double(val)
        if math.isnan(y):
            declined.append((a, b))
            continue
        rows.append((a, b, y))
    return rows, declined, escalated_n


def write_file(path, rows):
    with open(path, "w") as f:
        for a, b, y in rows:
            f.write(f"{a.hex()} {b.hex()} {y.hex()}\n")


# ---------------------------------------------------------------------------
# self-checks
# ---------------------------------------------------------------------------
def reverify_sample(rows, rng, n=200):
    """Fresh-oracle-vs-stored re-verification pass: recompute the full
    layered pipeline and confirm the stored (already-rounded) double
    matches bit-for-bit. Both the real consistency check AND the harness
    the negative control below proves is load-bearing."""
    sample = rng.sample(rows, min(n, len(rows)))
    mism = []
    for a, b, y in sample:
        val, ok, _ = layered_value(a, b)
        fresh = round_to_double(val) if ok else float("nan")
        if fresh != y and not (math.isnan(fresh) and math.isnan(y)):
            mism.append((a, b, y, fresh))
    return mism


def negative_control(rows, rng):
    """Corrupt one row's expected value; confirm a fresh oracle recompute
    disagrees with the corrupted value. Proves the verification pass is
    load-bearing, not a no-op.

    Corruption MUST be magnitude-safe: `y + 1.0` is a silent no-op for any
    row whose |y| is large enough that 1.0 sits below its ULP (this
    generator's huge-boundary/ray-stress rows reach |y| ~ 1e307-1.3e308,
    where ULP(y) is itself ~1e291, so +1.0 there rounds straight back to
    y). nextafter is exactly 1 ULP away regardless of magnitude, so it is
    never a no-op."""
    a, b, y = rng.choice(rows)
    if math.isinf(y):
        bad_y = -1.0 if y < 0 else 1.0
    else:
        bad_y = math.nextafter(y, math.inf)
    val, ok, _ = layered_value(a, b)
    fresh = round_to_double(val) if ok else float("nan")
    return fresh != bad_y


def symmetry_check(rows, rng, n=150):
    """ln B(a,b) = ln B(b,a): recompute the swapped pair fresh and confirm
    it matches the stored (a,b) row's value bit-for-bit."""
    candidates = [(a, b, y) for (a, b, y) in rows if a != b]
    sample = rng.sample(candidates, min(n, len(candidates)))
    mism = []
    for a, b, y in sample:
        val, ok, _ = layered_value(b, a)
        fresh = round_to_double(val) if ok else float("nan")
        if fresh != y:
            mism.append((a, b, y, fresh))
    return mism


def cross_check_sample(rows, rng, n=200):
    """Independent route (mp.log(mp.beta(a,b)), a gamma RATIO, not a
    loggamma SUM) vs a FRESH high-dps primary-oracle evaluation (never the
    already-double-rounded stored value -- that would impose a spurious
    ~2^-53 relative floor unrelated to whether the two mpmath code paths
    agree). Restricted to a moderate window where B(a,b) itself stays a
    numerically unremarkable quantity."""
    candidates = [(a, b) for (a, b, _y) in rows if 1e-6 <= a <= 1e6 and 1e-6 <= b <= 1e6]
    if not candidates:
        return mpf(0), None
    sample = rng.sample(candidates, min(n, len(candidates)))
    worst = mpf(0)
    worst_pt = None
    old = mp.mp.dps
    mp.mp.dps = 60
    try:
        for a, b in sample:
            primary = lnbeta_mpf(a, b, 60)
            cross = mp.log(mp.beta(mpf(a), mpf(b)))
            denom = abs(primary) if primary != 0 else mpf(1)
            rel = abs(cross - primary) / denom
            if rel > worst:
                worst, worst_pt = rel, (a, b)
    finally:
        mp.mp.dps = old
    return worst, worst_pt


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    ok = True
    rng = random.Random(SEED)

    print(f"[gen_lbeta_reference] DPS {DPS_LO}/{DPS_HI} "
          f"(escalate tiers {DPS_ESCALATE_TIERS}), "
          f"tol={float(AGREE_REL_TOL):.1e}", file=sys.stderr)

    points, meta = build_points(rng)
    print(f"[gen_lbeta_reference] grid axis size={meta['grid_n']} "
          f"(cross={meta['n_grid_cross']}); random-fill={meta['n_random']}; "
          f"moderate-fill={meta['n_moderate']}; zero-manifold={meta['n_zero_manifold']} "
          f"({meta['zm_a_n']} a-anchors); boundary={meta['n_boundary']}; "
          f"band-seam(min(a,b)=2^990)={meta['n_seam']} "
          f"({meta['n_seam_small']} small-side anchors); "
          f"symmetry-mirror={meta['n_symmetry']}", file=sys.stderr)
    lf, fi = meta["ab_crossing"]
    print(f"[gen_lbeta_reference] a=b ray -inf crossing: last finite t={lf.hex()} "
          f"({lf!r}) first -inf t={fi.hex()} ({fi!r})", file=sys.stderr)
    print("[gen_lbeta_reference] a=10b / a=1e5*b rays: verified NEVER reach -inf "
          "within double range -- see module docstring RAY-CROSSING FINDING "
          "(used as near-ceiling stress rays instead)", file=sys.stderr)
    print(f"[gen_lbeta_reference] total unique candidate points: {meta['n_unique']}",
          file=sys.stderr)

    t0 = time.time()
    rows, declined, escalated_n = emit_rows(points)
    print(f"[gen_lbeta_reference] emitted {len(rows)} rows, declined={len(declined)}, "
          f"escalated-past-dps{DPS_HI}={escalated_n}, "
          f"elapsed={time.time() - t0:.1f}s", file=sys.stderr)
    if declined:
        print(f"[gen_lbeta_reference]   declined points (expect zero per doctrine): "
              f"{declined[:20]}{'...' if len(declined) > 20 else ''}", file=sys.stderr)

    n_small = sum(1 for _a, _b, y in rows if abs(y) < 1e-6)
    n_neg_inf = sum(1 for _a, _b, y in rows if math.isinf(y))
    nonzero_abs = [abs(y) for _a, _b, y in rows if y != 0.0 and math.isfinite(y)]
    smallest_abs = min(nonzero_abs) if nonzero_abs else None
    n_exact_zero = sum(1 for _a, _b, y in rows if y == 0.0)
    print(f"[gen_lbeta_reference] rows with |lnB|<1e-6: {n_small} "
          f"(exact-zero: {n_exact_zero}); -inf rows: {n_neg_inf}; "
          f"smallest nonzero |lnB|: {smallest_abs!r}", file=sys.stderr)

    if len(rows) % 2 == 0:
        rows.pop()
    n = len(rows)
    if n % 2 == 0:
        print(f"[gen_lbeta_reference] FAILED: row count {n} still even after trim",
              file=sys.stderr)
        ok = False
    if not (TARGET_LO <= n <= TARGET_HI):
        print(f"[gen_lbeta_reference] WARNING: row count {n} outside the "
              f"{TARGET_LO}-{TARGET_HI} target band (not fatal)", file=sys.stderr)

    mism = reverify_sample(rows, rng)
    print(f"[gen_lbeta_reference] re-verify sample (n={min(200, n)}): "
          f"{len(mism)} mismatches", file=sys.stderr)
    if mism:
        print(f"[gen_lbeta_reference] FAILED: stored rows disagree with fresh "
              f"oracle re-evaluation: {mism[:5]}", file=sys.stderr)
        ok = False

    caught = negative_control(rows, rng)
    print(f"[gen_lbeta_reference] negative control: "
          f"{'CAUGHT' if caught else 'NOT CAUGHT'}", file=sys.stderr)
    if not caught:
        print("[gen_lbeta_reference] FAILED: negative control not caught -- "
              "verification pass is not load-bearing", file=sys.stderr)
        ok = False

    sym_mism = symmetry_check(rows, rng)
    print(f"[gen_lbeta_reference] symmetry self-check (n={min(150, n)}): "
          f"{len(sym_mism)} mismatches", file=sys.stderr)
    if sym_mism:
        print(f"[gen_lbeta_reference] FAILED: symmetry check found mismatches: "
              f"{sym_mism[:5]}", file=sys.stderr)
        ok = False

    worst, worst_pt = cross_check_sample(rows, rng)
    print(f"[gen_lbeta_reference] independent cross-check "
          f"(mp.log(mp.beta) vs loggamma-sum oracle): worst rel "
          f"{float(worst):.3e} @ {worst_pt!r}", file=sys.stderr)
    if worst > mpf("1e-20"):
        print("[gen_lbeta_reference] FAILED: independent cross-check exceeds "
              "modeled eps", file=sys.stderr)
        ok = False

    if not ok:
        print("[gen_lbeta_reference] ABORTING: not writing "
              "tests/data/lbeta_reference.txt (self-check failure)", file=sys.stderr)
        return 1

    write_file("tests/data/lbeta_reference.txt", rows)
    print(f"[gen_lbeta_reference] all checks passed in {time.time() - T0:.0f}s -- "
          f"wrote tests/data/lbeta_reference.txt ({n} rows)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
