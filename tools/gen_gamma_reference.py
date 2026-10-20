#!/usr/bin/env python3
"""Generate tests/data/gamma_{p,q,util}_reference.txt -- correctly rounded
oracles for corvus::gamma_p / corvus::gamma_q / the Log1pmxDd micro-gate.

Region map (lambda = x/a, a_T = 20):
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
  ridge; multi-minute hangs for a~1e250 even AWAY from the ridge). mpmath
  gammainc is called in exactly two places: (1) oracle_pq evaluates rows
  with a<=1e4 directly; (2) the Temme coefficient extraction (r_exact /
  extract_c) evaluates it at the 15 nodes a = 512*2^j, j=0..14 (512 ..
  8388608 ~ 8.4e6), on the 33 Chebyshev eta nodes of the fitted band
  [eta(1/2), eta(2)] only -- every such call converges (a < 1e7). No row
  with a>1e4 ever calls gammainc. For a>1e4: if a*phi(lambda) > 800 the small side
  has already underflowed past any double's representable range (e^-800 <<
  the smallest subnormal, e^-745ish) -- return the exact saturated pair
  (0,1)/(1,0) with no further computation. Otherwise (a>1e4, a*phi<=800)
  lambda is necessarily within roughly [0.5, 1.5] of the ridge (phi<=
  800/1e4=0.08 bounds it), safely inside the Temme fit's validated band
  [1/2, 2] -- use the exact-arithmetic Temme oracle: mpf evaluation of
  1/2*erfc(+-eta*sqrt(a/2)) +- e^{-a*phi}/sqrt(2*pi*a)*sum_{k<11} c_k(eta)/a^k
  with FULL-DEGREE (untruncated) Chebyshev fits, dps>=100. The k<11
  truncation alone is bounded at a>1e4 by |c_11|*a^-11 <~ 2e-3 * 1e-44
  (every extracted |c_k| for k<=11 stays under 2e-3 -- see
  gen_gamma_data.py's stderr degree/coefficient report), but it is not the
  dominant term: the measured total oracle error, which includes the
  extraction and 33-node Chebyshev fit error, is ~1e-35 relative on the
  small side (still ~19 decimal orders below double resolution).

  Self-checks (every mode, exit non-zero on a miss; run_oracle_self_checks):
  - overlap: Temme vs mpmath on a in [5e3,1e4] -- shows the switch has no
    seam, but those a lie between the extraction's own nodes, so it does
    not by itself certify the a>1e4 oracle;
  - certification: at off-node a in (1e4, 1e8] (past the last extraction
    node, 8.4e6), the small side from oracle_pq is compared with an
    independent route that never calls gammainc -- own power series for P
    (x<a), own Legendre continued fraction for Q (x>=a), at two working
    precisions that must agree with each other -- to tolerance 2^-100.
    a>1e8 is not certified directly (series/CF cost grows as sqrt(a)); it
    rests on the same fits with every c_k/a^k (k>=1) term smaller still;
  - negative control: the certification must reject an oracle whose c_0
    fit is perturbed by a relative 2^-40, a shift that moves the small side
    by more than 2^-53 at some certified point.

Usage:
    python3 tools/gen_gamma_reference.py               # full regeneration
    python3 tools/gen_gamma_reference.py --check-only  # self-checks only, writes nothing
    python3 tools/gen_gamma_reference.py --dd-special-append  # row-12 Log1pmxDd rows
"""

import argparse
import hashlib
import math
import os
import random
import struct
import sys
import tempfile
import time

import mpmath as mp

from refgen_common import round_to_double

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


# --- independent certification of the a > 1e4 Temme oracle --------------------
# Off-node a (the extraction nodes are 512*2^j; 8388608 is the last one and is
# kept as the endpoint) spanning (1e4, 1e8]. Per a, x is placed at a*phi in
# CERT_APHI on both sides of the ridge (0 = the ridge itself), so every point
# routes through temme_pq, never through saturation.
CERT_A = (12345.678, 3.7e4, 1e5, 5e5, 2e6, 8388608.0, 8.4e6, 3e7, 1e8)
CERT_APHI = (0, 1, 36, 300, 700)
CERT_DPS = (130, 160)
CERT_CONV_TOL = mp.mpf(10) ** -110  # the two precisions must agree this well
CERT_TOL = mp.mpf(2) ** -100        # oracle vs independent route
NEGCTL_K = 0                        # perturb c_0's fit ...
NEGCTL_REL = mp.mpf(2) ** -40       # ... by this relative amount
DOUBLE_RES = mp.mpf(2) ** -53


def indep_small_side(a, x, dps):
    """Small side of (P, Q) without gammainc: P by the power series
    P = x^a e^-x / Gamma(a+1) * sum_n x^n / ((a+1)...(a+n))   (x < a),
    Q by the Legendre continued fraction (modified Lentz)
    Q = x^a e^-x / Gamma(a) * 1/(x+1-a- 1(1-a)/(x+3-a- 2(2-a)/(x+5-a- ...)))
    (x >= a). Stops when the next term / CF factor changes the result by
    less than 10^-(dps-5). Returns (value, iterations)."""
    with mp.workdps(dps):
        a, x = mp.mpf(a), mp.mpf(x)
        eps = mp.mpf(10) ** -(dps - 5)
        cap = 10 ** 7
        if x >= a:
            tiny = mp.mpf(10) ** -(2 * dps)
            b = x + 1 - a
            c = 1 / tiny
            d = 1 / b
            h = d
            n = 0
            while True:
                n += 1
                an = -n * (n - a)
                b += 2
                d = an * d + b
                if d == 0:
                    d = tiny
                c = b + an / c
                if c == 0:
                    c = tiny
                d = 1 / d
                delta = d * c
                h *= delta
                if abs(delta - 1) < eps:
                    break
                if n > cap:
                    raise RuntimeError(f"CF did not converge at a={a} x={x}")
            pref = mp.exp(a * mp.log(x) - x - mp.loggamma(a))
            return +(pref * h), n
        term = mp.mpf(1)
        s = mp.mpf(1)
        n = 0
        while True:
            n += 1
            term *= x / (a + n)
            s += term
            # x < a + n: terms are strictly decreasing from here on.
            if term < s * eps and x < a + n:
                break
            if n > cap:
                raise RuntimeError(f"series did not converge at a={a} x={x}")
        pref = mp.exp(a * mp.log(x) - x - mp.loggamma(a + 1))
        return +(pref * s), n


def build_temme_cert_set():
    """(a, x, is_q_side, indep_value) per certification point, plus the
    worst two-precision disagreement. Raises if a point fails to converge
    or the two precisions disagree beyond CERT_CONV_TOL."""
    pts = []
    worst_conv = mp.mpf(0)
    for a in CERT_A:
        for tgt in CERT_APHI:
            for sgn in ((1,) if tgt == 0 else (1, -1)):
                lam = lam_of_eta(sgn * mp.sqrt(2 * mp.mpf(tgt) / a))
                x = float(lam * a)
                v_lo, _ = indep_small_side(a, x, CERT_DPS[0])
                v_hi, _ = indep_small_side(a, x, CERT_DPS[1])
                conv = abs((v_hi - v_lo) / v_hi)
                if conv > CERT_CONV_TOL:
                    raise RuntimeError(
                        f"independent route unconverged at a={a!r} x={x!r}: "
                        f"dps {CERT_DPS} disagree by {float(conv):.3e}")
                worst_conv = max(worst_conv, conv)
                pts.append((a, x, mp.mpf(x) >= mp.mpf(a), v_hi))
    return pts, worst_conv


def temme_cert_worst(pts, fits):
    """Worst relative small-side gap between oracle_pq(fits) and the
    independent values. Asserts each point actually takes the Temme path."""
    worst = mp.mpf(0)
    worst_at = None
    for a, x, q_side, ref in pts:
        a_m, x_m = mp.mpf(a), mp.mpf(x)
        assert a_m > A_SWITCH and a_m * phi_lam(x_m / a_m) <= APHI_SAT
        P, Q = oracle_pq(a, x, fits)
        got = Q if q_side else P
        rel = abs((got - ref) / ref)
        if rel > worst:
            worst, worst_at = rel, (a, x)
    return worst, worst_at


def run_oracle_self_checks(fits):
    """Overlap check, independent a>1e4 certification, and its negative
    control. Returns 0 on success, else a non-zero exit code."""
    t0 = time.time()
    print("self-check: oracle overlap on a in [5e3,1e4] ...", file=sys.stderr)
    worst, worst_at = check_oracle_overlap(fits)
    print(f"  worst rel diff {float(worst):.3e} "
          f"(2^{float(mp.log(worst, 2)):.2f}) at (a,lam)={worst_at}, "
          f"target 2^-60", file=sys.stderr)
    if worst > mp.mpf(2) ** -60:
        print("  FAILED: oracle overlap exceeds 2^-60 -- aborting.",
              file=sys.stderr)
        return 1

    print(f"self-check: independent certification of the a>1e4 Temme oracle "
          f"(own series/CF, dps {CERT_DPS}) ...", file=sys.stderr)
    pts, worst_conv = build_temme_cert_set()
    worst, worst_at = temme_cert_worst(pts, fits)
    print(f"  {len(pts)} points, a in {CERT_A}; independent route "
          f"self-agreement {float(worst_conv):.3e}", file=sys.stderr)
    print(f"  worst rel diff {float(worst):.3e} "
          f"(2^{float(mp.log(worst, 2)):.2f}) at (a,x)={worst_at}, "
          f"target 2^-100", file=sys.stderr)
    if worst > CERT_TOL:
        print("  FAILED: Temme oracle misses the independent route -- aborting.",
              file=sys.stderr)
        return 3

    bad = [list(row) for row in fits]
    bad[NEGCTL_K][0] *= 1 + NEGCTL_REL
    shift = mp.mpf(0)
    for a, x, q_side, _ref in pts:
        P0, Q0 = oracle_pq(a, x, fits)
        P1, Q1 = oracle_pq(a, x, bad)
        s0, s1 = (Q0, Q1) if q_side else (P0, P1)
        shift = max(shift, abs((s1 - s0) / s0))
    bad_worst, bad_at = temme_cert_worst(pts, bad)
    print(f"  negative control: c_{NEGCTL_K} fit perturbed by rel "
          f"2^{float(mp.log(NEGCTL_REL, 2)):.0f} -> max small-side shift "
          f"{float(shift):.3e} (double resolution 2^-53 = "
          f"{float(DOUBLE_RES):.3e}); certification gap {float(bad_worst):.3e} "
          f"at (a,x)={bad_at}", file=sys.stderr)
    if shift <= DOUBLE_RES or bad_worst <= CERT_TOL:
        print("  FAILED: negative control not flagged (or perturbation "
              "below double resolution) -- aborting.", file=sys.stderr)
        return 4
    print(f"  negative control rejected as required; self-checks "
          f"{time.time() - t0:.1f}s", file=sys.stderr)
    return 0


# --- helpers -----------------------------------------------------------------
def as_bits(x):
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def from_bits(b):
    return struct.unpack("<d", struct.pack("<Q", b))[0]


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
    """Deep-Temme-tail band: a*phi ~ 700-800 at moderately large a is
    exactly where the phi-series coefficient rounding (1/3, 1/5, 1/6, 1/7
    not exact in double) gets amplified through e^{-a*phi} -- this band
    exists to guard against that amplification. All these a exceed
    A_SWITCH, so oracle_pq already routes them through the exact-Temme path
    unconditionally.
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


def gen_diag_huge(ps):
    """Deterministic huge-a diagonal stratum (#12 regression rows).

    Draws NO rng -- called after the seeded-rng builders and before
    gen_dd_special_reference(rng) so the shared rng call sequence (and
    therefore that file's points) stays unchanged; see main().

    At a in {2^996, 2^997, 2^1000, DBL_MAX}, x = a +/- a few ULPs: k=0 is
    the exact diagonal (lambda=1 exactly, P=Q~0.5 via the Temme oracle at
    eta=0); k!=0 puts lambda within ~k*2^-52 of 1, so
    phi(lambda) ~ (lambda-1)^2/2 ~ k^2*2^-105 and a*phi ~ k^2*2^891 --
    astronomically past APHI_SAT=800 -- so oracle_pq's existing saturation
    branch returns the exact (1,0)/(0,1) pair with no Temme evaluation at
    all. No special-casing needed here beyond point selection.
    """
    n0 = len(ps.pts)
    for a0 in (2.0 ** 996, 2.0 ** 997, 2.0 ** 1000, sys.float_info.max):
        for k in range(-3, 4):
            x = from_bits(as_bits(a0) + k)
            if not math.isfinite(x):
                continue  # DBL_MAX, k>0 -> inf: skip
            ps.add(a0, x, keep_if_saturated=True)
    print(f"  diag-huge (#12 regression): {len(ps.pts) - n0} points",
          file=sys.stderr)


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
        Pf, Qf = round_to_double(P), round_to_double(Q)
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


def gen_dd_special_reference(rng):
    """dd_special_reference.txt: u phi_hi phi_lo, dd pairs for Log1pmxDd.

    phi(u) = u - log1p(u), computed at dps=60 (independent of the dps=100
    used for the Temme oracle above).

    Log1pmxDd itself lives in src/dd_special-inl.h, but its reference
    generation stays HERE deliberately: it consumes the same seeded rng
    stream as the P/Q point sets above, so carving it into its own
    generator would silently change every point in the checked-in file.
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
            phi_hi = round_to_double(phi)
            phi_lo = round_to_double(phi - mp.mpf(phi_hi))
            rows.append((u, phi_hi, phi_lo))

    with open("tests/data/dd_special_reference.txt", "w") as f:
        for u, phi_hi, phi_lo in rows:
            f.write(f"{hexd(u)} {hexd(phi_hi)} {hexd(phi_lo)}\n")
    print(f"  wrote tests/data/dd_special_reference.txt: {len(rows)} points",
          file=sys.stderr)
    return len(rows)


# --- incremental append: Log1pmxDd's LogDdAny branch (row 12) ------------------
# The base set above samples |u| > 1/16 at 700 points per side over wide
# ranges and never landed near u = 0.1166, where #58's nonzero-u.lo set found
# the primitive's worst (2^-71.41; the same u.hi with u.lo = 0 measures the
# same). This mode appends dense LogDdAny-branch rows from its OWN rng, so the
# main SEED stream -- and every gamma_{p,q} and existing dd_special row -- is
# untouched; existing lines stay byte-identical and rows already present are
# skipped (a re-run appends nothing). Same oracle and format as
# gen_dd_special_reference: phi = u - log1p(u) at dps 60, rounded to a dd.
#
# Strata (u.lo = 0, as the whole file):
#   near-cut   u in (1/16, 1/4], uniform -- the 2/u amplification is largest
#   pos        u in (1/4, 1], uniform
#   pos-wide   u in (1, 100], log-uniform
#   neg        u in [-1/2, -1/16), uniform
#   neg-deep   t = 1 + u log-uniform in [2^-52, 1/2], u = fl(t - 1); with
#              u in [-1, -1/2], 1 + u is exact (Sterbenz), so the 1-arg
#              form's w is exact -- its hazard near -1 is u.lo, which is 0
#   worst-nbhd the #58 worst u.hi = 0.11658359198367813 +- 256 ulps
#              bit-stepped, plus 400 uniform in +-2^-12 relative around it
#   cut-dense  u in (1/16, 0.064] (4000) and [-0.064, -1/16) (1000), uniform,
#              plus 0.0625809 +- 128 ulps: the true worst sits here (row 12)
#
# Self-checks (exit non-zero, file untouched):
#   (a) re-verify: every row recomputed at dps 120 must round to the SAME dd
#       pair, bit for bit;
#   (b) every want pair is normalized: |lo| <= ulp(hi)/2 and hi + lo == hi;
#   (c) negative control: a want pair with lo moved by one ulp(lo) -- the
#       smallest change a dd can express, ~2^-106 relative -- must fail (a);
#   (d) every appended u is in the LogDdAny branch (|u| > 1/16), and the
#       stratum counts meet the plan (no silent shrinkage).
DDSPEC12_SEED = 1212
DDSPEC12_WORST = 0.11658359198367813
DDSPEC12_CUT_PEAK = 0.0625809
DDSPEC12_PLAN = (("near-cut", 1500), ("pos", 800), ("pos-wide", 200),
                 ("neg", 800), ("neg-deep", 300), ("worst-nbhd", 913),
                 ("cut-dense", 5257))
DD_SPECIAL_PATH = os.path.join("tests", "data", "dd_special_reference.txt")


def build_ddspec12_points():
    """Returns [(u, stratum)] in file order, duplicates within the list removed."""
    rng = random.Random(DDSPEC12_SEED)
    cut = 1.0 / 16.0
    pts = []
    for _ in range(1500):
        pts.append((rng.uniform(cut, 0.25), "near-cut"))
    for _ in range(800):
        pts.append((rng.uniform(0.25, 1.0), "pos"))
    for _ in range(200):
        pts.append((10.0 ** rng.uniform(0.0, 2.0), "pos-wide"))
    for _ in range(800):
        pts.append((-rng.uniform(cut, 0.5), "neg"))
    # Near t = 2^-52 the 2^-53 grid of u collapses draws; draw to 300 distinct.
    deep = set()
    while len(deep) < 300:
        u = 2.0 ** rng.uniform(-52.0, -1.0) - 1.0
        if u not in deep:
            deep.add(u)
            pts.append((u, "neg-deep"))
    up = dn = DDSPEC12_WORST
    pts.append((DDSPEC12_WORST, "worst-nbhd"))
    for _ in range(256):
        up = math.nextafter(up, math.inf)
        dn = math.nextafter(dn, -math.inf)
        pts.append((up, "worst-nbhd"))
        pts.append((dn, "worst-nbhd"))
    for _ in range(400):
        pts.append((DDSPEC12_WORST * (1.0 + rng.uniform(-2.0 ** -12, 2.0 ** -12)),
                    "worst-nbhd"))
    # cut-dense (drawn last, so the strata above keep their draws): the
    # error peaks just above the cut, where log_dd's error is amplified by
    # ~2/u; 100k-point probes found 2^-67.63 near u = 0.0625809.
    for _ in range(4000):
        pts.append((rng.uniform(cut, 0.064), "cut-dense"))
    for _ in range(1000):
        pts.append((-rng.uniform(cut, 0.064), "cut-dense"))
    up = dn = DDSPEC12_CUT_PEAK
    pts.append((DDSPEC12_CUT_PEAK, "cut-dense"))
    for _ in range(128):
        up = math.nextafter(up, math.inf)
        dn = math.nextafter(dn, -math.inf)
        pts.append((up, "cut-dense"))
        pts.append((dn, "cut-dense"))
    out, seen = [], set()
    for u, tag in pts:
        b = as_bits(u)
        if b not in seen:
            seen.add(b)
            out.append((u, tag))
    return out


def phi_dd(u, dps):
    with mp.workdps(dps):
        um = mp.mpf(u)
        phi = um - mp.log1p(um)
        hi = round_to_double(phi)
        lo = round_to_double(phi - mp.mpf(hi))
    return hi, lo


def dd_normalized(hi, lo):
    return abs(lo) <= math.ulp(hi) / 2 and hi + lo == hi


def ddspec12_append():
    pts = build_ddspec12_points()
    with open(DD_SPECIAL_PATH, "rb") as f:
        old = f.read()
    if old and not old.endswith(b"\n"):
        print("FAILED: dd_special_reference.txt lacks a trailing newline -- "
              "nothing written.", file=sys.stderr)
        return 6
    eol = b"\r\n" if b"\r\n" in old else b"\n"
    present = {ln.split()[0] for ln in old.decode().splitlines() if ln.split()}

    bad = 0
    rows, counts, skipped = [], {}, 0
    for u, tag in pts:
        if not (abs(u) > 1.0 / 16.0 and u > -1.0 and math.isfinite(u)):
            print(f"  FAILED (d): u={u!r} ({tag}) outside the LogDdAny branch",
                  file=sys.stderr)
            bad += 1
            continue
        if hexd(u) in present:
            skipped += 1
            continue
        hi, lo = phi_dd(u, 60)
        if (hi, lo) != phi_dd(u, 120):
            print(f"  FAILED (a): u={u!r} dps 60/120 dd pairs differ",
                  file=sys.stderr)
            bad += 1
            continue
        if not dd_normalized(hi, lo):
            print(f"  FAILED (b): u={u!r} want pair not normalized",
                  file=sys.stderr)
            bad += 1
            continue
        rows.append((u, hi, lo))
        counts[tag] = counts.get(tag, 0) + 1

    # (c) negative control on the middle row: one ulp of lo must be caught.
    if rows:
        u, hi, lo = rows[len(rows) // 2]
        if (hi, math.nextafter(lo, math.inf)) == phi_dd(u, 120):
            print("  FAILED (c): 1-ulp(lo) perturbation not rejected",
                  file=sys.stderr)
            bad += 1
    # (d) plan sizes, unless this is a re-run that found everything present.
    if skipped == 0:
        for tag, want in DDSPEC12_PLAN:
            if counts.get(tag, 0) < want:
                print(f"  FAILED (d): stratum {tag} has {counts.get(tag, 0)} "
                      f"rows < {want}", file=sys.stderr)
                bad += 1
    if bad:
        print(f"Self-checks FAILED ({bad}) -- nothing written.", file=sys.stderr)
        return 1

    new = [f"{hexd(u)} {hexd(hi)} {hexd(lo)}".encode() + eol
           for u, hi, lo in rows]
    with open(DD_SPECIAL_PATH, "wb") as f:
        f.write(old + b"".join(new))
    print(f"  appended {len(new)} rows to dd_special_reference.txt "
          f"({skipped} already present); by stratum: {counts}; "
          f"self-checks (a)-(d) passed", file=sys.stderr)
    return 0


# --- incremental append of the #56 ridge / #57 Temme-tail families ------------
# Own RNG (Random(5657)); the main SEED stream is never touched. Only the new
# family is certified (same oracle_pq, same finite / saturation-prune /
# min-subnormal row filters as compute_and_write, preceded by the same
# extraction + oracle self-checks main() runs), then the rows are
# appended to both data files, leaving every existing line byte-identical.
# Ordering convention: file order == point-set insertion order, so new rows
# go at the end. Rows whose (a, x) is already present are skipped.
RIDGE56_SEED = 5657
RIDGE56_N = 600
RIDGE56_STEP_A = 60
TAIL57_N = 300
TAIL57_PHI_MIN = 36.5
RIDGE56_WITNESSES = (
    (15.171008041303773, 14.337519187856447),
    (1.5964951531210116, 2.355985895630986),
    (14.580163310119733, 15.71328391163746),
    (14.509912867374792, 15.54067538287069),
)
TAIL57_WITNESSES = (
    (221.53038557579066, 436.1874278050619),
    (2293.9196984916325, 4212.8856201041035),
    (4493.060743078933, 2958.983641334607),
)
RIDGE56_CKPT_PATH = os.path.join(tempfile.gettempdir(),
                                 f"corvus_gamma_ref_ckpt_ridge56_{RIDGE56_SEED}.tsv")


def build_ridge56_point_set():
    """Returns (PointSet, tags) with tags[(abits, xbits)] -> family label."""
    rng = random.Random(RIDGE56_SEED)
    ps = PointSet()
    tags = {}

    def add(a, x, tag, keep_sat):
        n = len(ps.pts)
        ps.add(a, x, keep_if_saturated=keep_sat)
        if len(ps.pts) > n:
            tags[(as_bits(float(a)), as_bits(float(x)))] = tag
            return True
        return False

    for a, x in RIDGE56_WITNESSES:
        add(a, x, "ridge-witness", True)
    for a, x in TAIL57_WITNESSES:
        add(a, x, "tail-witness", True)
    lo, hi = math.log10(1.5), math.log10(20.0)
    a_ridge = []
    while len(a_ridge) < RIDGE56_N:
        a = 10.0 ** rng.uniform(lo, hi)
        if not (1.5 <= a < 20.0):
            continue
        x = a * rng.uniform(0.8, 1.25)
        if add(a, x, "ridge", True):
            a_ridge.append(a)
    for a in rng.sample(a_ridge, RIDGE56_STEP_A):
        x = a + 1.0
        xs = {0: x}
        up = dn = x
        for k in (1, 2, 3):
            up = math.nextafter(up, math.inf)
            dn = math.nextafter(dn, -math.inf)
            xs[k], xs[-k] = up, dn
        for k in range(-3, 4):
            add(a, xs[k], "ridge", True)
    n_tail = 0
    lo, hi = math.log10(20.0), math.log10(3000.0)
    while n_tail < TAIL57_N:
        a = 10.0 ** rng.uniform(lo, hi)
        lam = rng.uniform(0.5, 2.0)
        if not (20.0 <= a <= 3000.0 and 0.5 < lam < 2.0):
            continue
        if a * (lam - 1.0 - math.log(lam)) <= TAIL57_PHI_MIN:
            continue
        if add(a, a * lam, "tail", False):
            n_tail += 1
    return ps, tags


def ridge56_append():
    t0 = time.time()
    print("extracting full-degree Temme fits for the a>1e4 oracle ...",
          file=sys.stderr)
    fits = extract_temme_full()
    rc = run_oracle_self_checks(fits)
    if rc:
        return rc

    ps, tags = build_ridge56_point_set()
    total = len(ps.pts)
    dig = hashlib.sha256()
    for a, x, _k, _m in ps.pts:
        dig.update(struct.pack("<QQ", as_bits(a), as_bits(x)))
    sig = f"ridge56-v1-{dig.hexdigest()[:16]} SEED={RIDGE56_SEED} N={total}"
    print(f"  ridge56/tail57 point set: {total} points", file=sys.stderr)

    done = {}
    mode = "w"
    if os.path.exists(RIDGE56_CKPT_PATH):
        with open(RIDGE56_CKPT_PATH, "r") as f:
            if f.readline().strip() == sig:
                mode = "a"
                for line in f:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) >= 2:
                        done[int(parts[0])] = parts[1:]
            else:
                print("  checkpoint signature mismatch -- starting fresh.",
                      file=sys.stderr)
    print(f"  checkpoint: {len(done)}/{total} points already computed",
          file=sys.stderr)
    with open(RIDGE56_CKPT_PATH, mode) as fh:
        if mode == "w":
            fh.write(sig + "\n")
            fh.flush()
        for idx, (a, x, _k, _m) in enumerate(ps.pts):
            if idx in done:
                continue
            try:
                P, Q = oracle_pq(a, x, fits)
                if not (mp.isfinite(P) and mp.isfinite(Q)):
                    raise ValueError("non-finite oracle value")
                if abs(P + Q - 1) > mp.mpf(10) ** -80:
                    print(f"\nFAILED: P+Q != 1 at a={a!r} x={x!r}: "
                          f"P={float(P)!r} Q={float(Q)!r} -- nothing written.",
                          file=sys.stderr)
                    return 7
                Pf, Qf = round_to_double(P), round_to_double(Q)
                if not (math.isfinite(Pf) and math.isfinite(Qf)):
                    raise ValueError("non-finite rounded value")
                rec = [hexd(Pf), hexd(Qf)]
            except Exception as e:
                rec = ["FAILED", f"{type(e).__name__}:{e}"]
            fh.write(str(idx) + "\t" + "\t".join(rec) + "\n")
            fh.flush()
            if (idx + 1) % 100 == 0:
                print(f"    ... {idx + 1}/{total} ({time.time() - t0:.0f}s)",
                      file=sys.stderr)
        # reload
    done = {}
    with open(RIDGE56_CKPT_PATH, "r") as f:
        f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            done[int(parts[0])] = parts[1:]

    rows = []
    n_pruned = 0
    n_min_sub = 0
    declined = []
    for idx, (a, x, keep_sat, avoid_min_sub) in enumerate(ps.pts):
        rec = done[idx]
        if rec[0] == "FAILED":
            declined.append((a, x, rec[1]))
            continue
        Pf, Qf = float.fromhex(rec[0]), float.fromhex(rec[1])
        if not keep_sat and Pf in (0.0, 1.0) and Qf in (0.0, 1.0):
            n_pruned += 1
            continue
        if avoid_min_sub and (Pf == MIN_SUBNORMAL or Qf == MIN_SUBNORMAL):
            n_min_sub += 1
            continue
        rows.append((a, x, Pf, Qf))
    for a, x, why in declined:
        print(f"  DECLINED (not emitted): a={a!r} x={x!r} ({why})", file=sys.stderr)
    got = {(a, x): (Pf, Qf) for a, x, Pf, Qf in rows}
    for w in RIDGE56_WITNESSES + TAIL57_WITNESSES:
        if w not in got:
            print(f"\nFAILED: witness {w!r} was not certified -- nothing written.",
                  file=sys.stderr)
            return 5
    print(f"  certified {len(rows)} rows; {len(declined)} declined, "
          f"{n_pruned} incidental saturations pruned, {n_min_sub} min-subnormal",
          file=sys.stderr)

    plan = []
    for name, col in (("gamma_p_reference.txt", 2), ("gamma_q_reference.txt", 3)):
        path = os.path.join("tests", "data", name)
        with open(path, "rb") as f:
            old = f.read()
        if old and not old.endswith(b"\n"):
            print(f"\nFAILED: {name} lacks a trailing newline -- nothing written.",
                  file=sys.stderr)
            return 6
        eol = b"\r\n" if b"\r\n" in old else b"\n"
        present = set()
        for ln in old.decode().splitlines():
            t = ln.split()
            if len(t) == 4:
                present.add((t[0], t[1]))
        new = []
        counts = {}
        for a, x, Pf, Qf in rows:
            if (hexd(a), hexd(x)) in present:
                print(f"  already present, skipped: {name} a={hexd(a)} x={hexd(x)}",
                      file=sys.stderr)
                continue
            new.append(f"{hexd(a)} {hexd(x)} {hexd(Pf)} {hexd(Qf)}".encode() + eol)
            fam = tags[(as_bits(a), as_bits(x))].split("-")[0]
            counts[fam] = counts.get(fam, 0) + 1
        plan.append((name, path, old, new, counts))
    for name, path, old, new, counts in plan:
        with open(path, "wb") as f:
            f.write(old + b"".join(new))
        print(f"  appended {len(new)} rows to {name} (by family: {counts})",
              file=sys.stderr)
    for a, x in RIDGE56_WITNESSES + TAIL57_WITNESSES:
        Pf, Qf = got[(a, x)]
        print(f"  witness a={a!r} x={x!r} -> P={Pf!r} ({hexd(Pf)}) "
              f"Q={Qf!r} ({hexd(Qf)})", file=sys.stderr)
    print(f"ridge56 runtime (this invocation): {time.time()-t0:.1f}s",
          file=sys.stderr)
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                     allow_abbrev=False)
    parser.add_argument("--ridge56-append", action="store_true",
                        help="certify only the #56 ridge / #57 Temme-tail "
                             "families and append them to "
                             "tests/data/gamma_{p,q}_reference.txt")
    parser.add_argument("--check-only", action="store_true",
                        help="extract the Temme fits and run the oracle "
                             "self-checks only; write no files")
    parser.add_argument("--dd-special-append", action="store_true",
                        help="append the row-12 LogDdAny-branch rows to "
                             "tests/data/dd_special_reference.txt (own rng; "
                             "no Temme work, gamma files untouched)")
    args = parser.parse_args()
    if args.ridge56_append:
        return ridge56_append()
    if args.dd_special_append:
        return ddspec12_append()
    rng = random.Random(SEED)

    print("extracting full-degree Temme fits for the a>1e4 oracle ...",
          file=sys.stderr)
    fits = extract_temme_full()

    rc = run_oracle_self_checks(fits)
    if rc or args.check_only:
        return rc

    print("generating region point sets ...", file=sys.stderr)
    ps = PointSet()
    gen_r1(ps, rng)
    gen_r4(ps, rng)
    gen_r2(ps, rng)
    gen_r3(ps, rng)
    gen_r3_deep_tail(ps)
    gen_diag_huge(ps)  # no rng draw -- rng stream into gen_dd_special_reference below is unaffected
    print(f"  total distinct (a,x) points: {len(ps.pts)}", file=sys.stderr)

    print("evaluating oracle + writing gamma_p/q_reference.txt ...",
          file=sys.stderr)
    n_pq = compute_and_write(ps, fits)

    print("generating dd_special_reference.txt ...", file=sys.stderr)
    n_util = gen_dd_special_reference(rng)

    if n_pq < 10000 or n_util < 1500:
        print(f"FAILED: point counts too low (pq={n_pq}, util={n_util})",
              file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
