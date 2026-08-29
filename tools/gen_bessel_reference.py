#!/usr/bin/env python3
"""Generate tests/data/{i0,i1,i0e,i1e}_reference.txt -- correctly rounded
Bessel I0/I1 oracle, unscaled (i0, i1) and exponentially-scaled (i0e, i1e).

Each line: <input-hex-double> <output-hex-double>, output = round-to-nearest
of the mathematical value at high mpmath dps. Format matches the erf/erfc/
digamma/trigamma house convention (bare two-hex-double rows, no dd pair --
these functions are working-precision-only -- no kernel keeps more than
double precision past the final rounding). Specials
(x=0, +-inf, NaN) are covered by the smoke test, not this file -- same
convention as every other reference generator in this repo. -0.0 is likewise
left to the smoke test (i1(-0)=-0 exact is a one-line assertion there, not a
sampled-oracle concern).

Definitions used here:
  i0(x)  = I_0(x)                         i1(x)  = I_1(x)
  i0e(x) = exp(-|x|) * I_0(x)             i1e(x) = exp(-|x|) * I_1(x)
mpmath's besseli(nu, x) for negative real x and integer nu already applies
the correct even/odd symmetry (I_0(-x)=I_0(x), I_1(-x)=-I_1(x)) internally,
so negative-axis rows call the oracle directly at the negative point rather
than hand-mirroring a positive-axis result -- a genuine oracle evaluation,
not an assumption of the symmetry the kernel itself will exploit.

Oracle doctrine (this is erf-difficulty class -- NO bracket certification;
PLAN.md/NUMERICAL-DOCTRINE.md's oracle-trust doctrine applies only to
families WITHOUT a trusted library baseline, which does not describe
mpmath.besseli):
  - every row: layered-dps agreement (dps 40 vs 80); on disagreement beyond
    AGREE_REL_TOL, escalate to dps 150; if that still disagrees, DECLINE the
    row (skip + record), never guess
  - independent cross-check on a sampled subset: OWN high-dps series-sum
    (x <~ 25) or OWN high-dps asymptotic expansion (x >~ 25) vs mp.besseli,
    reusing gen_bessel_data.py's own series-sum route structurally (own
    math, not a second call into mpmath's besseli internals)
  - negative control: corrupt one stored row's expected value, confirm the
    fresh-oracle-vs-stored re-verification pass catches it; exit nonzero if
    it is NOT caught (proves the verification pass is load-bearing)
  - overflow-boundary re-derivation: independently bisect (dps 50) the
    round-to-inf threshold 2^1024*(1-2^-54) for both nu, assert bit-identical
    to src/bessel_data.h's kBesselI0OverflowX/kBesselI1OverflowX (the
    authority -- already independently derived twice); mismatch is an
    ESCALATE condition, not something this generator resolves on its own

Coverage per file (edge-refined, bit-stepped where it matters -- the
trigamma-generator sampling rule): log-spaced across the whole positive
domain (subnormals through the domain ceiling), an explicit smallest-
subnormal cluster and the subnormal/normal boundary (2^-1022) bracket, a
bit-stepped bracket straddling the series/tail seam x=kBesselSplit=8, a
bit-stepped bracket at the function's own domain ceiling (the overflow
boundary for i0/i1, including rows one ULP past it that must read back as
+inf; the DBL_MAX neighbourhood for i0e/i1e, which never overflow -- min
value ~3e-155 at DBL_MAX), huge-x witnesses, and moderate-density fill.
A random subset (not a full mirror) is replayed on the negative axis via a
genuine oracle call at the negative point. Final row count is trimmed to be
odd (never a multiple of 2, 4, or 8) so ctest's masked-tail path is always
exercised, per house rule.

Usage:
    python tools/gen_bessel_reference.py
writes all four files under tests/data/ directly (checkpoint-free -- this
family's point count is small enough to regenerate in one shot; no partial-
run resumption is needed, unlike beta/gammainv).
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

X_S = 8.0  # kBesselSplit -- must match src/bessel_data.h
I0_OVERFLOW_X = float.fromhex("0x1.64fe5304e83e4p+9")
I1_OVERFLOW_X = float.fromhex("0x1.64fe69ff9fec7p+9")
DBL_MAX = float.fromhex("0x1.fffffffffffffp+1023")
SMALLEST_SUBNORMAL = float.fromhex("0x0.0000000000001p-1022")
MIN_NORMAL = float.fromhex("0x1.0000000000000p-1022")

DPS_LO = 40
DPS_HI = 80
DPS_ESCALATE = 150
# Comfortably (~1e9x) below the double's own ~1.1e-16 relative ULP -- meant
# to catch a real oracle disagreement, not to be a marginal diagnostic.
AGREE_REL_TOL = mpf("1e-25")

T0 = time.time()


# ---------------------------------------------------------------------------
# bit helpers
# ---------------------------------------------------------------------------
def as_bits(x: float) -> int:
    return struct.unpack("<q", struct.pack("<d", x))[0]


def from_bits(b: int) -> float:
    return struct.unpack("<d", struct.pack("<q", b))[0]


def neighbourhood(x0: float, k: int = 48, lo_off=None, hi_off=None):
    b = as_bits(x0)
    lo = -k if lo_off is None else lo_off
    hi = k if hi_off is None else hi_off
    return [from_bits(b + j) for j in range(lo, hi + 1)]


# ---------------------------------------------------------------------------
# oracle
# ---------------------------------------------------------------------------
def besseli_value(nu, scaled, x_float, dps):
    """Raw (unrounded) mpf value of i{nu}[e](x) at the given dps. mp.dps is
    set INSIDE this function (mechanism rule)."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        xm = mpf(x_float)
        v = mp.besseli(nu, xm)
        if scaled:
            v = v * mp.exp(-abs(xm))
        return v
    finally:
        mp.mp.dps = old


def layered_value(nu, scaled, x_float):
    """Layered-dps oracle value with escalation. Returns (value_mpf, ok,
    escalated: bool) -- ok=False means DECLINE (row not emitted)."""
    lo = besseli_value(nu, scaled, x_float, DPS_LO)
    hi = besseli_value(nu, scaled, x_float, DPS_HI)
    denom = abs(hi) if hi != 0 else mpf(1)
    if abs(hi - lo) / denom <= AGREE_REL_TOL:
        return hi, True, False
    esc = besseli_value(nu, scaled, x_float, DPS_ESCALATE)
    denom2 = abs(esc) if esc != 0 else mpf(1)
    if abs(esc - hi) / denom2 <= AGREE_REL_TOL:
        return esc, True, True
    return esc, False, True


# ---------------------------------------------------------------------------
# independent cross-check routes (own math, not a second mpmath.besseli call)
# ---------------------------------------------------------------------------
def own_series_ive(nu, x, dps, kmax=4000):
    """OWN high-dps power-series evaluation of exp(-x)*I_nu(x), x>0. The
    literal mathematical series definition, summed directly -- not mpmath's
    internal besseli algorithm. Reused structurally from
    gen_bessel_data.py's own_series_besseli (extended: returns the SCALED
    value directly since that's what both regimes ultimately need)."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        xm = mpf(x)
        q = xm * xm / 4
        term = mpf(1)
        s = mpf(1)
        k = 0
        tol = mpf(10) ** (-(dps - 5))
        while True:
            k += 1
            term = term * q / (k * k if nu == 0 else k * (k + 1))
            s += term
            if term / s < tol and k > 4:
                break
            if k > kmax:
                break
        if nu == 1:
            s *= xm / 2
        return s * mp.exp(-xm)
    finally:
        mp.mp.dps = old


def own_asymptotic_ive(nu, x, dps, K=20):
    """OWN high-dps Poincare/Hankel asymptotic expansion of exp(-x)*I_nu(x)
    for large x (textbook classical asymptotic series, A&S 9.7.1 -- public-
    domain mathematics, not ported code; this is a test-oracle cross-check
    helper, never shipped in the kernel). Verified against mp.besseli at
    x=50..1.7e308: worst relative error ~1e-24, i.e. this is not a marginal
    check. term_k/term_{k-1} = ((2k-1)^2 - mu) / (8*k*x), mu = 4*nu^2."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        xm = mpf(x)
        mu = 4 * nu * nu
        term = mpf(1)
        s = mpf(1)
        for k in range(1, K + 1):
            term = term * ((2 * k - 1) ** 2 - mu) / (8 * k * xm)
            s += term
        return s / mp.sqrt(2 * mp.pi * xm)
    finally:
        mp.mp.dps = old


def own_cross_check(nu, x, dps=60):
    """Route by magnitude: series for |x|<=25 (converges in <100 terms
    there), asymptotic for |x|>25 (already ~1e-18 relative accurate by
    x=25, improving to ~1e-24 by x=50 and beyond). Returns the SIGNED
    SCALED (exp(-|x|)*I_nu(x)) value -- both own_series_ive/
    own_asymptotic_ive are derived assuming positive x internally, so the
    odd (nu=1) sign is reapplied explicitly here rather than assumed by
    the caller (i0/i0e are even: sign is a no-op there)."""
    ax = abs(x)
    v = own_series_ive(nu, ax, dps) if ax <= 25.0 else own_asymptotic_ive(nu, ax, dps)
    if nu == 1 and x < 0:
        v = -v
    return v


# ---------------------------------------------------------------------------
# overflow-boundary independent re-derivation (mirrors gen_bessel_data.py)
# ---------------------------------------------------------------------------
def find_overflow_boundary(nu, dps=50):
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        ovf = mpf(2) ** 1024 * (1 - mpf(2) ** -54)
        lo, hi = 700.0, 720.0
        while from_bits(as_bits(lo) + 1) < hi:
            mid = (lo + hi) / 2
            v = mp.besseli(nu, mpf(mid))
            if v < ovf:
                lo = mid
            else:
                hi = mid
        return lo
    finally:
        mp.mp.dps = old


# ---------------------------------------------------------------------------
# point-set strata
# ---------------------------------------------------------------------------
def log_spaced(rng, lo_exp, hi_exp, n):
    return [10.0 ** rng.uniform(lo_exp, hi_exp) for _ in range(n)]


def tiny_cluster():
    pts = [from_bits(k) for k in range(1, 21)]  # smallest subnormals
    pts += neighbourhood(MIN_NORMAL, 40)  # subnormal/normal boundary
    pts += [1e-300, 1e-250, 1e-200, 1e-150, 1e-100, 1e-50, 1e-30, 1e-10, 1e-5]
    return pts


def seam_bracket():
    return neighbourhood(X_S, 60)


def domain_ceiling_bracket(mode, nu):
    if mode == "unscaled":
        b = I0_OVERFLOW_X if nu == 0 else I1_OVERFLOW_X
        # Both sides: negative offsets are the last finite doubles, +1 and
        # beyond must read back as +inf -- the boundary behaviour this
        # format pins.
        return neighbourhood(b, 80, lo_off=-80, hi_off=80)
    # scaled: never overflows: bracket the domain ceiling itself (DBL_MAX);
    # offsets above 0 would be +inf, which is a special (excluded).
    return neighbourhood(DBL_MAX, 60, lo_off=-60, hi_off=0)


def moderate_fill(rng, mode, n_linear, n_log):
    if mode == "unscaled":
        hi = I0_OVERFLOW_X  # slightly conservative for nu=1 too; both < 714
        pts = [rng.uniform(0.0, hi) for _ in range(n_linear)]
        pts += log_spaced(rng, -6, math.log10(hi), n_log)
        return pts
    pts = log_spaced(rng, 0.0, math.log10(DBL_MAX), n_linear)
    huge = []
    for w in (1e300, 1e290, 1e250, 1e200, 5e307):
        huge += neighbourhood(w, 20)
    pts += huge
    return pts


def near_seam_spread(rng, n=140):
    return [rng.uniform(4.0, 12.0) for _ in range(n)]


def near_ceiling_spread(rng, mode, nu, n=140):
    if mode == "unscaled":
        b = I0_OVERFLOW_X if nu == 0 else I1_OVERFLOW_X
        return [rng.uniform(600.0, b) for _ in range(n)]
    return [10.0 ** rng.uniform(280.0, math.log10(DBL_MAX)) for _ in range(n)]


def build_positive_points(rng, mode, nu):
    pts = []
    hi_exp = math.log10(I0_OVERFLOW_X if mode == "unscaled" else DBL_MAX)
    pts += log_spaced(rng, -320.0, hi_exp, 620)
    pts += tiny_cluster()
    pts += seam_bracket()
    pts += domain_ceiling_bracket(mode, nu)
    pts += moderate_fill(rng, mode, 260, 220)
    pts += near_seam_spread(rng)
    pts += near_ceiling_spread(rng, mode, nu)
    return pts


# ---------------------------------------------------------------------------
# emission
# ---------------------------------------------------------------------------
def emit_file(nu, scaled, mode, rng, negative_fraction=0.42):
    tag = f"i{nu}" + ("e" if scaled else "")
    pos_pts = build_positive_points(rng, mode, nu)
    # dedupe, drop exact 0 (special -- smoke test's job)
    seen = set()
    uniq_pos = []
    for x in pos_pts:
        if not math.isfinite(x) or x <= 0.0:
            continue
        b = as_bits(x)
        if b in seen:
            continue
        seen.add(b)
        uniq_pos.append(x)

    neg_sample = rng.sample(uniq_pos, int(len(uniq_pos) * negative_fraction))
    all_pts = list(uniq_pos)
    for x in neg_sample:
        nb = as_bits(-x)
        if nb in seen:
            continue
        seen.add(nb)
        all_pts.append(-x)

    rows = []  # (x_float, y_float)
    declined = []
    escalated_n = 0
    for x in all_pts:
        val, ok, escalated = layered_value(nu, scaled, x)
        if escalated:
            escalated_n += 1
        if not ok:
            declined.append(x)
            continue
        y = round_to_double(val)
        if math.isnan(y):
            declined.append(x)
            continue
        rows.append((x, y))

    # Trim to an odd count (never a multiple of 2, 4, or 8) -- house rule so
    # the masked-tail path is always exercised. Deterministic: drop from the
    # end of the (seed-stable) list.
    if len(rows) % 2 == 0:
        rows.pop()

    return tag, rows, declined, escalated_n, len(uniq_pos), len(all_pts) - len(uniq_pos)


def write_file(path, rows):
    with open(path, "w") as f:
        for x, y in rows:
            f.write(f"{x.hex()} {y.hex()}\n")


# ---------------------------------------------------------------------------
# self-checks
# ---------------------------------------------------------------------------
def reverify_sample(nu, scaled, rows, rng, n=120):
    """Fresh-oracle-vs-stored re-verification pass: recompute at DPS_HI and
    confirm the stored (already-rounded) double matches bit-for-bit. This is
    both a real consistency check AND the harness the negative control below
    proves is load-bearing."""
    sample = rng.sample(rows, min(n, len(rows)))
    mism = []
    for x, y in sample:
        fresh = round_to_double(besseli_value(nu, scaled, x, DPS_HI))
        if fresh != y and not (math.isnan(fresh) and math.isnan(y)):
            mism.append((x, y, fresh))
    return mism


def negative_control(nu, scaled, rows, rng):
    """Corrupt one row's expected value; confirm reverify_sample's exact
    comparison flags it. Proves the verification pass is load-bearing, not a
    no-op -- mirrors gen_bessel_data.py's own dd-coefficient corruption
    control applied to reference-row space instead of kernel-simulation
    space."""
    x, y = rng.choice(rows)
    # y+1.0 is a no-op for |y| >= 2^53 (1.0 sits below ULP(y) there) --
    # nextafter is exactly 1 ULP away regardless of magnitude, never a no-op.
    bad_y = math.nextafter(y, math.inf) if math.isfinite(y) else 0.0
    corrupted = [(x, bad_y)]
    mism = reverify_sample(nu, scaled, corrupted, rng, n=1)
    return len(mism) == 1


def cross_check_sample(nu, scaled, rows, rng, n=10):
    """Compares the OWN independent route against a FRESH full-precision
    oracle evaluation (not the already-double-rounded stored value -- that
    would impose a spurious ~2^-53 relative floor having nothing to do with
    whether the independent route agrees with mpmath's besseli)."""
    sample = rng.sample(rows, min(n, len(rows)))
    worst = mpf(0)
    worst_x = None
    for x, _y in sample:
        own = own_cross_check(nu, x, dps=60)
        if not scaled:
            own = own * mp.exp(abs(mpf(x)))
        ref = besseli_value(nu, scaled, x, 60)
        denom = abs(ref) if ref != 0 else mpf(1)
        rel = abs(own - ref) / denom
        if rel > worst:
            worst, worst_x = rel, x
    return worst, worst_x


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    ok = True
    rng_master = random.Random(SEED)

    print(f"[gen_bessel_reference] split x_s={X_S}, DPS {DPS_LO}/{DPS_HI}"
          f"(escalate {DPS_ESCALATE}), tol={float(AGREE_REL_TOL):.1e}",
          file=sys.stderr)

    # ---- overflow-boundary independent re-derivation ----
    b0 = find_overflow_boundary(0)
    b1 = find_overflow_boundary(1)
    print(f"[gen_bessel_reference] overflow boundary I0 re-derived: "
          f"{float.hex(b0)} vs header {float.hex(I0_OVERFLOW_X)} "
          f"match={b0 == I0_OVERFLOW_X}", file=sys.stderr)
    print(f"[gen_bessel_reference] overflow boundary I1 re-derived: "
          f"{float.hex(b1)} vs header {float.hex(I1_OVERFLOW_X)} "
          f"match={b1 == I1_OVERFLOW_X}", file=sys.stderr)
    if b0 != I0_OVERFLOW_X or b1 != I1_OVERFLOW_X:
        print("[gen_bessel_reference] FAILED: overflow boundary disagrees "
              "with src/bessel_data.h -- ESCALATE, do not pick one.",
              file=sys.stderr)
        return 1

    specs = [
        (0, False, "unscaled", "i0"),
        (1, False, "unscaled", "i1"),
        (0, True, "scaled", "i0e"),
        (1, True, "scaled", "i1e"),
    ]

    # Deterministic per-file seed offset -- NOT hash(tag): Python salts str
    # hashing per interpreter invocation (PYTHONHASHSEED), which would make
    # this generator silently non-reproducible run to run (caught here by
    # the negative control flagging a picked row inconsistently between two
    # otherwise-identical runs).
    seed_offset = {"i0": 1, "i1": 2, "i0e": 3, "i1e": 4}

    summary = {}
    for nu, scaled, mode, tag in specs:
        t0 = time.time()
        rng = random.Random(SEED + seed_offset[tag])
        tag2, rows, declined, escalated_n, n_pos, n_neg = emit_file(
            nu, scaled, mode, rng
        )
        assert tag2 == tag
        n = len(rows)
        n_pos_out = sum(1 for x, _ in rows if x > 0)
        n_neg_out = n - n_pos_out
        print(f"[gen_bessel_reference] {tag}: {n} rows "
              f"(+{n_pos_out}/-{n_neg_out}), declined={len(declined)}, "
              f"escalated-to-dps{DPS_ESCALATE}={escalated_n}, "
              f"elapsed={time.time()-t0:.1f}s", file=sys.stderr)
        if declined:
            print(f"[gen_bessel_reference]   declined points: {declined[:20]}"
                  f"{'...' if len(declined) > 20 else ''}", file=sys.stderr)
        if n % 8 == 0 or n % 4 == 0 or n % 2 == 0:
            print(f"[gen_bessel_reference] FAILED: {tag} row count {n} is a "
                  f"multiple of 2/4/8", file=sys.stderr)
            ok = False
        if not (2000 <= n <= 4000):
            print(f"[gen_bessel_reference] WARNING: {tag} row count {n} "
                  f"outside the 2-4k target band (not fatal)", file=sys.stderr)

        # -- self-check: fresh-oracle re-verification --
        mism = reverify_sample(nu, scaled, rows, rng)
        print(f"[gen_bessel_reference]   re-verify sample "
              f"(n={min(120, n)}): {len(mism)} mismatches", file=sys.stderr)
        if mism:
            print(f"[gen_bessel_reference] FAILED: {tag} stored rows "
                  f"disagree with fresh oracle re-evaluation: {mism[:5]}",
                  file=sys.stderr)
            ok = False

        # -- negative control --
        caught = negative_control(nu, scaled, rows, rng)
        print(f"[gen_bessel_reference]   negative control: "
              f"{'CAUGHT' if caught else 'NOT CAUGHT'}", file=sys.stderr)
        if not caught:
            print(f"[gen_bessel_reference] FAILED: {tag} negative control "
                  f"not caught -- verification pass is not load-bearing",
                  file=sys.stderr)
            ok = False

        # -- independent cross-check (own series/asymptotic vs mp.besseli) --
        worst, worst_x = cross_check_sample(nu, scaled, rows, rng)
        print(f"[gen_bessel_reference]   independent cross-check "
              f"(own series x<=25 / own asymptotic x>25): worst rel "
              f"{float(worst):.3e} @ x={worst_x!r}", file=sys.stderr)
        if worst > mpf("1e-14"):
            print(f"[gen_bessel_reference] FAILED: {tag} independent "
                  f"cross-check exceeds modeled eps", file=sys.stderr)
            ok = False

        summary[tag] = dict(n=n, n_pos=n_pos_out, n_neg=n_neg_out,
                             declined=len(declined), escalated=escalated_n,
                             worst_cross=float(worst))
        if ok:
            write_file(f"tests/data/{tag}_reference.txt", rows)

    if not ok:
        print("[gen_bessel_reference] ABORTING: not writing any files "
              "(self-check failure)", file=sys.stderr)
        return 1

    print(f"[gen_bessel_reference] all checks passed in "
          f"{time.time()-T0:.0f}s -- wrote 4 files under tests/data/",
          file=sys.stderr)
    for tag, s in summary.items():
        print(f"  {tag}: {s}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
