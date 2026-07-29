#!/usr/bin/env python3
"""Generate tests/data/gamma_{p,q,util}_reference.txt -- correctly rounded
oracles for corvus::gamma_p / corvus::gamma_q / the Log1pmxDd micro-gate.

Per PLAN.md "Phase C part 2", region map (lambda = x/a, a_T = 20):
  R1 series-P:    {a<20, 0<x<=a+1} u {a>=20, lambda<=1/2}
  R4 small-a Q:   {0<a<=3/2, 0<x<=4}
  R2 backward-CF: {a<20, x>a+1} u {a>=20, lambda>=2}, minus R4
  R3 Temme:       {a>=20, 1/2<lambda<2}

Oracle rules (both hard requirements):
  SMALL-SIDE-DIRECT: x>=a -> compute Q via regularized upper gammainc,
  P=1-Q; x<a -> compute P via regularized lower, Q=1-P. Never subtract a
  near-1 value.

  ORACLE TRAP: mpmath's lower-gammainc path (hyp1f1) genuinely hangs/fails
  to converge for large a (measured: NoConvergence for a>=1e7 near the
  ridge; multi-minute hangs for a~1e250 even AWAY from the ridge). Use
  mpmath only for a<=1e4. For a>1e4: if a*phi(lambda) > 800 the small side
  has already underflowed past any double's representable range (e^-800 <<
  the smallest subnormal, e^-745ish) -- return the exact saturated pair
  (0,1)/(1,0) with no further computation. Otherwise (a>1e4, a*phi<=800)
  lambda is necessarily within roughly [0.5, 1.5] of the ridge (phi<=
  800/1e4=0.08 bounds it), safely inside the Temme fit's validated band
  [1/2, 2] -- use the exact-arithmetic Temme oracle: mpf evaluation of
  1/2*erfc(+-eta*sqrt(a/2)) +- e^{-a*phi}/sqrt(2*pi*a)*sum_{k<11} c_k(eta)/a^k
  with FULL-DEGREE (untruncated) Chebyshev fits, dps>=100. Truncation
  error at a>1e4 is bounded by |c_11|*a^-11 <~ 2e-3 * 1e-44 << 1e-40 (c_11
  was the first coefficient this generator's own extraction drops at the
  2^-60 tail cut, and every extracted |c_k| for k<=11 stays under 2e-3 --
  see gen_gamma_data.py's stderr degree/coefficient report). Self-checked
  below on the a in [5e3,1e4] overlap band against mpmath directly.

Usage:
    python3 tools/gen_gamma_reference.py
"""

import math
import random
import struct
import sys

import mpmath as mp

mp.mp.dps = 100

SEED = 20260727
A_SWITCH = mp.mpf("1e4")
APHI_SAT = mp.mpf(800)

# --- Temme extraction (duplicated from gen_gamma_data.py -- generators are
# independent scripts by project convention; this one only needs the
# FULL-DEGREE, untruncated fits as an oracle, not the emitted table). -------
KEXT = 15
K = 11
NNODES = 33
A0 = 512


def phi_lam(lam):
    return lam - 1 - mp.log(lam)


ETA_LO = -mp.sqrt(2 * phi_lam(mp.mpf("0.5")))
ETA_HI = mp.sqrt(2 * phi_lam(mp.mpf(2)))
ETA_MID = (ETA_HI + ETA_LO) / 2
ETA_HALF = (ETA_HI - ETA_LO) / 2


def lam_of_eta(eta):
    if eta == 0:
        return mp.mpf(1)
    target = eta * eta / 2
    if eta > 0:
        lo, hi = mp.mpf(1), mp.mpf(2)
        while hi - 1 - mp.log(hi) < target:
            hi *= 2
    else:
        lo, hi = mp.mpf("1e-30"), mp.mpf(1)
    for _ in range(400):
        mid = (lo + hi) / 2
        v = mid - 1 - mp.log(mid) - target
        if eta > 0:
            lo, hi = (mid, hi) if v < 0 else (lo, mid)
        else:
            lo, hi = (lo, mid) if v < 0 else (mid, hi)
    return (lo + hi) / 2


def r_exact(a, eta, lam):
    """Same P-side/Q-side oracle-trap split as gen_gamma_data.py."""
    if eta >= 0:
        Q = mp.gammainc(a, lam * a, regularized=True)
        base = Q - mp.erfc(eta * mp.sqrt(a / 2)) / 2
    else:
        P = mp.gammainc(a, 0, lam * a, regularized=True)
        base = mp.erfc(-eta * mp.sqrt(a / 2)) / 2 - P
    return base * mp.sqrt(2 * mp.pi * a) * mp.exp(a * eta * eta / 2)


def extract_c(eta, lam, a0=A0, kext=KEXT):
    A = mp.matrix(kext, kext)
    b = mp.matrix(kext, 1)
    for j in range(kext):
        a = mp.mpf(a0) * 2 ** j
        v = 1 / a
        for k in range(kext):
            A[j, k] = v ** k
        b[j] = r_exact(a, eta, lam)
    return mp.lu_solve(A, b)


def cheb_coeffs_from_vals(vals):
    n = len(vals)
    out = []
    for j in range(n):
        s = mp.fsum([vals[i] * mp.cos(j * mp.pi * (2 * i + 1) / (2 * n))
                      for i in range(n)])
        out.append(s * 2 / n if j else s / n)
    return out


def clenshaw(coefs, t):
    b1 = b2 = mp.mpf(0)
    for j in range(len(coefs) - 1, 0, -1):
        b1, b2 = 2 * t * b1 - b2 + coefs[j], b1
    return t * b1 - b2 + coefs[0]


def extract_temme_full():
    """Full-degree (untruncated), high-precision Chebyshev fits c_0..c_10."""
    cvals = []
    for i in range(NNODES):
        t = mp.cos(mp.pi * (2 * i + 1) / (2 * NNODES))
        eta = ETA_MID + ETA_HALF * t
        lam = lam_of_eta(eta)
        cvals.append(extract_c(eta, lam))
    fits = []
    for k in range(K):
        vals = [cvals[i][k] for i in range(NNODES)]
        fits.append(cheb_coeffs_from_vals(vals))
    return fits


def temme_pq(a, lam, phi, fits):
    """Full-precision Temme evaluation; valid where a*phi <= ~800 (forces
    lambda near the ridge, well inside the fitted band for a > A_SWITCH)."""
    eta = mp.sqrt(2 * phi)
    if lam < 1:
        eta = -eta
    t = (eta - ETA_MID) / ETA_HALF
    ck = [clenshaw(row, t) for row in fits]
    S = mp.mpf(0)
    for k in range(K - 1, -1, -1):
        S = S / a + ck[k]
    R = mp.exp(-a * phi) / mp.sqrt(2 * mp.pi * a) * S
    z = eta * mp.sqrt(a / 2)
    if lam >= 1:
        Q = mp.erfc(z) / 2 + R
        P = 1 - Q
    else:
        P = mp.erfc(-z) / 2 - R
        Q = 1 - P
    return P, Q


def oracle_pq(a, x, fits):
    """SMALL-SIDE-DIRECT with the a<=1e4 mpmath / a>1e4 Temme-or-saturated
    switch. a, x: python float or mpf, both > 0 and finite."""
    a_m, x_m = mp.mpf(a), mp.mpf(x)
    if a_m <= A_SWITCH:
        if x_m >= a_m:
            Q = mp.gammainc(a_m, x_m, regularized=True)
            P = 1 - Q
        else:
            P = mp.gammainc(a_m, 0, x_m, regularized=True)
            Q = 1 - P
        return P, Q
    lam = x_m / a_m
    phi = phi_lam(lam)
    if a_m * phi > APHI_SAT:
        return (mp.mpf(1), mp.mpf(0)) if lam >= 1 else (mp.mpf(0), mp.mpf(1))
    return temme_pq(a_m, lam, phi, fits)


def check_oracle_overlap(fits):
    """Self-check: on a in [5e3,1e4] with lambda in the Temme band, the
    mpmath path (actually used there) and the Temme path (used just above
    the switch) must agree to 2^-60 -- validates the switch has no seam."""
    worst = mp.mpf(0)
    worst_at = None
    a_vals = [mp.mpf(v) for v in ("5000", "6000", "7500", "9000", "9999", "10000")]
    lam_vals = [mp.mpf(v) for v in
                ("0.55", "0.7", "0.85", "0.97", "1.0", "1.05", "1.2", "1.5", "1.9")]
    for a in a_vals:
        for lam in lam_vals:
            phi = phi_lam(lam)
            if a * phi > APHI_SAT:
                continue
            x = lam * a
            if x >= a:
                Qd = mp.gammainc(a, x, regularized=True)
                Pd = 1 - Qd
            else:
                Pd = mp.gammainc(a, 0, x, regularized=True)
                Qd = 1 - Pd
            Pt, Qt = temme_pq(a, lam, phi, fits)
            small_d, small_t = (Qd, Qt) if lam >= 1 else (Pd, Pt)
            if small_d == 0:
                continue
            rel = abs((small_d - small_t) / small_d)
            if rel > worst:
                worst, worst_at = rel, (float(a), float(lam))
    return worst, worst_at


# --- helpers -----------------------------------------------------------------
def as_bits(x):
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def hexd(x):
    """Hex-float of a python float, mpf, or int -- round-trip exact."""
    return float(x).hex()


def log_grid(lo, hi, n):
    llo, lhi = math.log10(lo), math.log10(hi)
    return [10 ** (llo + (lhi - llo) * i / (n - 1)) for i in range(n)]


NEXT_UP = lambda v: math.nextafter(v, math.inf)
NEXT_DN = lambda v: math.nextafter(v, -math.inf)


class PointSet:
    """Collects (a, x, keep_if_saturated, avoid_min_subnormal) tuples,
    deduped by (a,x) bits."""

    def __init__(self):
        self.seen = set()
        self.pts = []

    def add(self, a, x, keep_if_saturated=True, avoid_min_subnormal=False):
        a, x = float(a), float(x)
        if not (math.isfinite(a) and math.isfinite(x) and a > 0 and x > 0):
            return
        key = (as_bits(a), as_bits(x))
        if key in self.seen:
            return
        self.seen.add(key)
        self.pts.append((a, x, keep_if_saturated, avoid_min_subnormal))


# --- region point generation --------------------------------------------------
def gen_r1(ps, rng):
    n0 = len(ps.pts)
    a_list = sorted(set(log_grid(1e-300, 20.0, 70) + [
        1e-300, 1e-100, 1e-30, 1e-8, 0.01, 0.3, 0.9,
        NEXT_DN(1.5), 1.5, NEXT_UP(1.5), 2.0, 3.7,
        NEXT_DN(8.0), 8.0, NEXT_UP(8.0), 15.0, 19.99,
    ]))
    fracs = [0.02, 0.1, 0.3, 0.6, 0.9, 0.99, 1.0 - 2.0 ** -52]
    for a in a_list:
        for f in fracs:
            ps.add(a, f * (a + 1.0))
    # Random fill, dense, to broaden coverage beyond the structured grid.
    # keep_if_saturated=False: incidental full-saturation from pure random
    # noise adds no information beyond the structured/designed points above.
    for _ in range(5000):
        a = 10.0 ** rng.uniform(-300.0, math.log10(20.0))
        f = rng.uniform(0.0, 1.0)
        ps.add(a, f * (a + 1.0), keep_if_saturated=False)
    # a>=20 cross set, filtered by a*phi<=800.
    for lam in (0.01, 0.05, 0.2, 0.35, 0.499, 0.5):
        for a in (20.0, 100.0, 4e3, 1e8):
            phi = float(phi_lam(mp.mpf(lam)))
            if a * phi > 800.0:
                continue
            ps.add(a, lam * a)
    print(f"  R1: {len(ps.pts) - n0} points", file=sys.stderr)


def gen_r4(ps, rng):
    n0 = len(ps.pts)
    gamma_const = float(mp.euler)
    e_neg_gamma = math.exp(-gamma_const)
    a_list = sorted(set(log_grid(1e-300, 1.5, 40) + [
        1e-300, 1e-100, 1e-16, 1e-8, 1e-4, 0.01, 0.1, 0.3, 0.561, 0.9,
        1.2, NEXT_DN(1.5), 1.5, NEXT_UP(1.5),
    ]))
    x_list = sorted(set([
        1e-30, 0.01, 0.1, 0.3, 0.561, 0.9, 1.5, 2.5, 3.3,
        NEXT_DN(e_neg_gamma), e_neg_gamma, NEXT_UP(e_neg_gamma),
        NEXT_DN(4.0), 4.0, NEXT_UP(4.0),
    ]))
    for a in a_list:
        for x in x_list:
            if x <= 4.0:
                ps.add(a, x)
        # x spanning (a+1, 4) and x <= a+1, per spec, for a where a+1 < 4.
        if a + 1.0 < 4.0:
            for f in (0.3, 0.7, 0.99, 1.0, 1.01, 1.5, 2.0):
                x = f * (a + 1.0)
                if 0 < x <= 4.0:
                    ps.add(a, x)
    for _ in range(4200):
        a = 10.0 ** rng.uniform(-300.0, math.log10(1.5))
        x = rng.uniform(1e-30, 4.0)
        ps.add(a, x, keep_if_saturated=False)
    print(f"  R4: {len(ps.pts) - n0} points", file=sys.stderr)


def gen_r2(ps, rng):
    n0 = len(ps.pts)
    a_list = sorted(set(log_grid(1e-300, 1e6, 60) + [
        1e-300, 1e-100, 1e-8, 0.01, 0.3, 0.9, 1.5, 2.0, 3.7, 8.0, 15.0, 19.99,
        20.0, 100.0, 1000.0,
    ]))
    for a in a_list:
        base = a + 1.0
        x_candidates = [base, NEXT_UP(base), a + 1.5, 2.0 * base, 5.0 * base]
        # deep tail up to Q-underflow: push aphi well past 800.
        x_candidates.append(a + 40.0 * max(a, 1.0) + 800.0)
        for x in x_candidates:
            if a <= 1.5 and x <= 4.0:
                continue  # a<=3/2 uses x>4 only here (R4 owns x<=4)
            ps.add(a, x)
    # a<=3/2 bracket, x=4+1ulp only (the R2 sliver just past R4's boundary).
    for a in (1e-300, 1e-8, 0.1, 0.5, 1.0, 1.5):
        for x in (NEXT_UP(4.0), NEXT_UP(NEXT_UP(4.0))):
            ps.add(a, x)
    # Ridge-adjacent lambda cross set for moderate a.
    for a in (20.0, 100.0, 1000.0, 2600.0):
        for lam in (2.0, NEXT_UP(2.0), 2.6, 4.0, 10.0):
            ps.add(a, lam * a)
    for _ in range(4900):
        a = 10.0 ** rng.uniform(-300.0, 6.0)
        f = rng.uniform(1.0001, 20.0)
        ps.add(a, f * (a + 1.0), keep_if_saturated=False)
    print(f"  R2: {len(ps.pts) - n0} points", file=sys.stderr)


def gen_r3(ps, rng):
    n0 = len(ps.pts)
    a_list = [
        NEXT_DN(20.0), 20.0, NEXT_UP(20.0), 25.0, 100.0, 1000.0, 3.7e4, 1e8,
        1e10, 9e15, math.nextafter(2.0 ** 53, -math.inf), 2.0 ** 53,
        math.nextafter(2.0 ** 53, math.inf), 1e250,
    ]
    lam_band = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.93, 0.96,
                0.98, 0.99, 1.0, 1.01, 1.02, 1.05, 1.1, 1.2, 1.3, 1.4, 1.5,
                1.6, 1.8, 2.0]
    ridge_ks = list(range(1, 21)) + [26, 30, 40, 50, 52]

    def maybe_add(a, lam):
        phi = float(phi_lam(mp.mpf(lam)))
        if a * phi > 800.0:
            return
        ps.add(a, lam * a)

    for a in a_list:
        for lam in lam_band:
            maybe_add(a, lam)
        for k in ridge_ks:
            maybe_add(a, 1.0 - 2.0 ** -k)
            maybe_add(a, 1.0 + 2.0 ** -k)
        ps.add(a, a)  # x == a exactly
        # z=6 seam: solve a*phi(lam)=36 for lam on both sides, bracket +-2ulp.
        for target_lam_guess, sign in ((1.3, 1), (0.7, -1)):
            lo, hi = (1.0, 4.0) if sign > 0 else (1e-6, 1.0)
            for _ in range(60):
                mid = (lo + hi) / 2
                phi_mid = float(phi_lam(mp.mpf(mid)))
                if (a * phi_mid - 36.0) * sign < 0:
                    lo = mid
                else:
                    hi = mid
            x0 = mid * a
            for dx in (-2, -1, 0, 1, 2):
                xv = x0
                for _ in range(abs(dx)):
                    xv = NEXT_UP(xv) if dx > 0 else NEXT_DN(xv)
                ps.add(a, xv)
        # Subnormal band: aphi in {700,720,740,745,750} on both lambda sides.
        for target in (700.0, 720.0, 740.0, 745.0, 750.0):
            for sign in (1, -1):
                lo, hi = (1.0, 4.0) if sign > 0 else (1e-6, 1.0)
                mid = 1.0
                for _ in range(60):
                    mid = (lo + hi) / 2
                    phi_mid = float(phi_lam(mp.mpf(mid)))
                    if (a * phi_mid - target) * sign < 0:
                        lo = mid
                    else:
                        hi = mid
                ps.add(a, mid * a)
    # A few explicitly saturated points per side (a*phi ~ 850): kept even
    # though both P and Q round to exactly 1/0 -- this is the point.
    for a in (100.0, 1000.0, 1e8):
        for target in (820.0, 835.0, 850.0):
            for sign in (1, -1):
                lo, hi = (1.0, 4.0) if sign > 0 else (1e-6, 1.0)
                mid = 1.0
                for _ in range(60):
                    mid = (lo + hi) / 2
                    phi_mid = float(phi_lam(mp.mpf(mid)))
                    if (a * phi_mid - target) * sign < 0:
                        lo = mid
                    else:
                        hi = mid
                ps.add(a, mid * a, keep_if_saturated=True)
    print(f"  R3: {len(ps.pts) - n0} points", file=sys.stderr)


def gen_r3_deep_tail(ps):
    """Deep-Temme-tail band [added after kernel review]: a*phi ~ 700-800 at
    moderately large a is exactly where the phi-series coefficient rounding
    (1/3, 1/5, 1/6, 1/7 not exact in double) gets amplified through
    e^{-a*phi} -- measured ~12 ULP at a=3.79e5, lambda=1.062 before
    kGammaPhiCoefLo existed. All these a exceed A_SWITCH, so oracle_pq
    already routes them through the exact-Temme path unconditionally.
    """
    n0 = len(ps.pts)
    a_list = (1e5, 3.79e5, 1e6, 3e6, 1e7)
    deltas = (0.03, 0.045, 0.055, 0.0625, 0.08)  # 0.0625 = the |u|=1/16 cut
    for a in a_list:
        for d in deltas:
            for lam in (1.0 + d, 1.0 - d):
                phi = float(phi_lam(mp.mpf(lam)))
                if a * phi > 800.0:
                    continue
                # avoid_min_subnormal=True: keep subnormal-result points,
                # but not ones landing on the single smallest subnormal,
                # where ULP distance downward is meaningless.
                ps.add(a, lam * a, keep_if_saturated=True, avoid_min_subnormal=True)
    print(f"  R3 deep-tail (phi-series coefficient band): "
          f"{len(ps.pts) - n0} points", file=sys.stderr)


# --- oracle evaluation + emission ---------------------------------------------
MIN_SUBNORMAL = math.ldexp(1.0, -1074)  # smallest positive double


def compute_and_write(ps, fits):
    rows = []
    n_pruned = 0
    n_min_subnormal = 0
    for a, x, keep_sat, avoid_min_sub in ps.pts:
        try:
            P, Q = oracle_pq(a, x, fits)
        except Exception:
            continue
        if not (mp.isfinite(P) and mp.isfinite(Q)):
            continue
        Pf, Qf = float(P), float(Q)
        if not (math.isfinite(Pf) and math.isfinite(Qf)):
            continue
        if not keep_sat:
            if Pf in (0.0, 1.0) and Qf in (0.0, 1.0):
                n_pruned += 1
                continue
        if avoid_min_sub and (Pf == MIN_SUBNORMAL or Qf == MIN_SUBNORMAL):
            # At exactly the smallest subnormal, ULP distance below is
            # undefined -- not useful for the deep-tail gate this point
            # exists to feed.
            n_min_subnormal += 1
            continue
        rows.append((a, x, Pf, Qf))
    print(f"  total rows after oracle eval: {len(rows)} "
          f"(pruned {n_pruned} incidental saturations, "
          f"{n_min_subnormal} at the minimum subnormal)", file=sys.stderr)
    for path in ("tests/data/gamma_p_reference.txt", "tests/data/gamma_q_reference.txt"):
        with open(path, "w") as f:
            for a, x, P, Q in rows:
                f.write(f"{hexd(a)} {hexd(x)} {hexd(P)} {hexd(Q)}\n")
        print(f"  wrote {path}: {len(rows)} points", file=sys.stderr)
    return len(rows)


def gen_util_reference(rng):
    """gamma_util_reference.txt: u phi_hi phi_lo, dd pairs for Log1pmxDd.

    phi(u) = u - log1p(u), computed at dps=60 per PLAN.md's own spec for
    this table (independent of the dps=100 used for the Temme oracle above).
    """
    with mp.workdps(60):
        pts = set()
        # Series side: |u| in [2^-60, 1/16), log-spaced, both signs.
        for _ in range(1000):
            mag = 2.0 ** rng.uniform(-60.0, math.log2(1.0 / 16.0))
            pts.add(mag)
            pts.add(-mag)
        # Boundary brackets around the 1/16 cut.
        cut = 1.0 / 16.0
        for v in (NEXT_DN(cut), cut, NEXT_UP(cut), NEXT_UP(NEXT_UP(cut))):
            pts.add(v)
            pts.add(-v)
        # Log side, negative: (1/16, 0.99...], incl named brackets, u>-1.
        for _ in range(700):
            pts.add(-rng.uniform(1.0 / 16.0, 0.999999999))
        for v in (0.9, 0.99, 0.999999, 1.0 - 2.0 ** -40):
            pts.add(-v)
        # Log side, positive: (1/16, 100], incl named points.
        for _ in range(700):
            pts.add(rng.uniform(1.0 / 16.0, 100.0))
        for v in (1.0, 3.5, 10.0, 100.0):
            pts.add(v)

        rows = []
        seen_bits = set()
        for u in pts:
            if not (math.isfinite(u) and u > -1.0):
                continue
            b = as_bits(u)
            if b in seen_bits:
                continue
            seen_bits.add(b)
            um = mp.mpf(u)
            phi = um - mp.log1p(um)
            phi_hi = float(phi)
            phi_lo = float(phi - mp.mpf(phi_hi))
            rows.append((u, phi_hi, phi_lo))

    with open("tests/data/gamma_util_reference.txt", "w") as f:
        for u, phi_hi, phi_lo in rows:
            f.write(f"{hexd(u)} {hexd(phi_hi)} {hexd(phi_lo)}\n")
    print(f"  wrote tests/data/gamma_util_reference.txt: {len(rows)} points",
          file=sys.stderr)
    return len(rows)


def main():
    rng = random.Random(SEED)

    print("extracting full-degree Temme fits for the a>1e4 oracle ...",
          file=sys.stderr)
    fits = extract_temme_full()

    print("self-check: oracle overlap on a in [5e3,1e4] ...", file=sys.stderr)
    worst, worst_at = check_oracle_overlap(fits)
    print(f"  worst rel diff {float(worst):.3e} "
          f"(2^{float(mp.log(worst, 2)):.2f}) at (a,lam)={worst_at}, "
          f"target 2^-60", file=sys.stderr)
    if worst > mp.mpf(2) ** -60:
        print("  FAILED: oracle overlap exceeds 2^-60 -- aborting.",
              file=sys.stderr)
        return 1

    print("generating region point sets ...", file=sys.stderr)
    ps = PointSet()
    gen_r1(ps, rng)
    gen_r4(ps, rng)
    gen_r2(ps, rng)
    gen_r3(ps, rng)
    gen_r3_deep_tail(ps)
    print(f"  total distinct (a,x) points: {len(ps.pts)}", file=sys.stderr)

    print("evaluating oracle + writing gamma_p/q_reference.txt ...",
          file=sys.stderr)
    n_pq = compute_and_write(ps, fits)

    print("generating gamma_util_reference.txt ...", file=sys.stderr)
    n_util = gen_util_reference(rng)

    if n_pq < 10000 or n_util < 1500:
        print(f"FAILED: point counts too low (pq={n_pq}, util={n_util})",
              file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
