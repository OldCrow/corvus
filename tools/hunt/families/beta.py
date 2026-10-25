"""beta_p / beta_q hunt strata (#82).

Bounds (docs/ACCURACY.md, beta_p and beta_q), per region and side, the
direct side being the one the kernel computes (orientation replicated
from src/beta-inl.h's router, below):

  region                 direct   complement
  R1 power series          1         1
  R2 continued fraction    1         0  (CR claim; #76 disputes it)
  R3 Temme ridge           3         1
  R4 tiny-min corner       2         0  (CR claim)
  R4 near-one post-route   1         0  (CR claim)
  R2 gamma-limit slice     0         0  (CR claim)

The 0 cells are scored at bound 0 with tie_exempt = 0: #76 is correcting
several false CR claims (R2 complement among them), so violations there
are expected and still counted. Rows whose routing the float replica
cannot decide (a predicate within ~1e-12 of its threshold, the near-one
bar within 2^-40) take the larger bound of the candidate regions.
"""
from __future__ import annotations

import math
import random

import mpmath as mp
from mpmath import mpf

from core import Stratum, full_mantissa, log_uniform, neighbours, subnormal
from ._beta_oracles import (approx_I, betainc_pq, betainc_pq_cf)

# ---- kernel constants (src/beta_data.h, src/beta-inl.h, src/gamma_data.h) ----
XI1 = float.fromhex("0x1.ccccccccccccdp-2")     # kBetaXi1 = 0.45
B1 = 8.0                                        # kBetaB1
EPS_R4 = 2.0 ** -6                              # kBetaEpsR4
T_RIDGE = 32.0                                  # kBetaTRidge
Z0 = 10.0                                       # kBetaZ0
C_LG = 256.0                                    # kBetaClg
LN2 = float.fromhex("0x1.62e42fefa39efp-1")     # kBetaLn2
GAMMA_LIM = 2.0 ** 47                           # kBetaGammaLim
NEAR_ONE = float.fromhex("0x1.ffc0000000000p-1")  # kBetaNearOne = 1 - 2^-11
PR_TAU_MAX = 2.5                                # kBetaPrTauMax
GL_RIDGE_MIN = 20.0                             # kBetaGlRidgeMin
GAMMA_AT = 20.0                                 # kGammaAT (gamma-limit sub-map)
Z2_SPLIT = 25.0                                 # kBetaZ2Split (cpsi; 36 until #70)
EXP_FLOOR = -800.0                              # kBetaExpFloor
TAU_SUBN = 2.0 ** -950                          # kBetaTauSubn
B_TINY_CUT = 2.0 ** -140                        # kBetaBTinyCut
SCALE_ABOVE = 2.0 ** 900                        # kBetaScaleAbove
DIV_FLOOR = 2.0 ** -900                         # kBetaDivFloor
TWO_PI_NU_CLAMP = 2.0 ** 900                    # kBetaTwoPiNuClamp
BINET_Z_CLAMP = 2.0 ** 500                      # kBetaBinetZClamp

SEAMS = [
    {"fn": "beta", "name": "xi1_x", "coord": "x", "value": XI1, "source": "kBetaXi1 (R1 native box)"},
    {"fn": "beta", "name": "xi1_y", "coord": "1-x", "value": XI1, "source": "kBetaXi1 (R1 swap box)"},
    {"fn": "beta", "name": "b1_bx", "coord": "b*x", "value": B1, "source": "kBetaB1 (R1 box)"},
    {"fn": "beta", "name": "epsr4_min", "coord": "min(a,b)", "value": EPS_R4, "source": "kBetaEpsR4"},
    {"fn": "beta", "name": "r4_ln2", "coord": "tau*|ln xi|", "value": LN2, "source": "kBetaLn2 (R4 window)"},
    {"fn": "beta", "name": "r2_thr", "coord": "x", "value": "(a+1)/(a+b+2)", "source": "R2 orientation / #47 R4 window"},
    {"fn": "beta", "name": "near_one", "coord": "I (R1 value)", "value": NEAR_ONE, "source": "kBetaNearOne post-route"},
    {"fn": "beta", "name": "pr_tau_max", "coord": "alpha", "value": PR_TAU_MAX, "source": "kBetaPrTauMax"},
    {"fn": "beta", "name": "lgamma_table_max", "coord": "min(a,b)", "value": 16.0, "source": "kBetaLgammaAnalyticMax: PA lgamma(min) table form above kBetaPrTauMax (#85)"},
    {"fn": "beta", "name": "r4_tau1", "coord": "alpha", "value": 1.0, "source": "BetaR4Tiny gt1 (lgamma(1+tau) zone)"},
    {"fn": "beta", "name": "tridge_nu", "coord": "nu", "value": T_RIDGE, "source": "kBetaTRidge"},
    {"fn": "beta", "name": "band_rat1_lo", "coord": "x/p", "value": 0.5, "source": "kBetaXiRatioLo"},
    {"fn": "beta", "name": "band_rat1_hi", "coord": "x/p", "value": 2.0, "source": "kBetaXiRatioHi"},
    {"fn": "beta", "name": "band_rat2_lo", "coord": "y/q", "value": 0.5, "source": "kBetaXiRatioLo"},
    {"fn": "beta", "name": "band_rat2_hi", "coord": "y/q", "value": 2.0, "source": "kBetaXiRatioHi"},
    {"fn": "beta", "name": "z2split", "coord": "cpsi", "value": Z2_SPLIT, "source": "kBetaZ2Split (R3 core/tail)"},
    {"fn": "beta", "name": "r3_sat", "coord": "cpsi", "value": -EXP_FLOOR, "source": "kBetaExpFloor (R3 saturation)"},
    {"fn": "beta", "name": "r3_pswap", "coord": "p=a/(a+b)", "value": 0.5, "source": "BetaR3Temme pswap"},
    {"fn": "beta", "name": "clg_c", "coord": "a+b", "value": C_LG, "source": "kBetaClg (PA/PB prefactor)"},
    {"fn": "beta", "name": "z0_min", "coord": "min(a,b)", "value": Z0, "source": "kBetaZ0 (PA/PB prefactor)"},
    {"fn": "beta", "name": "gammalim", "coord": "max(a,b)", "value": GAMMA_LIM, "source": "kBetaGammaLim"},
    {"fn": "beta", "name": "gl_ridge", "coord": "min(a,b) (max>=2^47)", "value": GL_RIDGE_MIN, "source": "kBetaGlRidgeMin, gated on s since #84"},
    {"fn": "beta", "name": "gl_ser_cf", "coord": "t (s<20)", "value": "s+1", "source": "gamma-limit sub-map (kGammaAT)"},
    {"fn": "beta", "name": "gl_at", "coord": "s", "value": GAMMA_AT, "source": "kGammaAT in the gamma-limit slice"},
    {"fn": "beta", "name": "gl_sat", "coord": "E_g", "value": EXP_FLOOR, "source": "BetaClampE / kBetaGlSatT (slice saturation)"},
    {"fn": "beta", "name": "scale_above", "coord": "max(a,b)", "value": SCALE_ABOVE, "source": "kBetaScaleAbove"},
    {"fn": "beta", "name": "div_floor", "coord": "max(a,b)", "value": DIV_FLOOR, "source": "kBetaDivFloor (DdDivDd, #59)"},
    {"fn": "beta", "name": "tau_subn", "coord": "min(a,b)", "value": TAU_SUBN, "source": "kBetaTauSubn (R4)"},
    {"fn": "beta", "name": "b_tiny_cut", "coord": "max(a,b) in R4", "value": B_TINY_CUT, "source": "kBetaBTinyCut (R4)"},
    {"fn": "beta", "name": "two_pi_nu_clamp", "coord": "nu", "value": TWO_PI_NU_CLAMP, "source": "kBetaTwoPiNuClamp"},
    {"fn": "beta", "name": "binet_z_clamp", "coord": "max(a,b)", "value": BINET_Z_CLAMP, "source": "kBetaBinetZClamp (LgammaDiffDd)"},
]

# Regions -> (direct bound, complement bound)
BOUNDS = {"R1": (1, 1), "R2": (1, 1), "R3": (3, 1), "R4": (2, 0),  # R2 cmp 1 since #76
          "R4pr": (1, 0), "GL": (1, 1)}  # GL 1 since #90


# ---- router replica (src/beta-inl.h BetaVec routing, float-exact) -------------
def _near(v, t, rel=1e-12):
    return abs(v - t) <= rel * abs(t)


def route(a, b, x):
    """[(region, swapped), ...] candidate routings: one entry normally, more
    when a float predicate sits within ~1e-12 of its threshold."""
    y0 = 1.0 - x
    tau, bmax = min(a, b), max(a, b)
    sw4 = b < a
    xt = y0 if sw4 else x
    if sw4:
        lxt = float(mp.log(mp.fsub(1, mpf(x), exact=True))) if x < 1 else -math.inf
    else:
        lxt = math.log(x) if x > 0 else -math.inf
    thr = 1.0 / (1.0 + (b + 1.0) / (a + 1.0))
    r2tf = (x < thr) != sw4
    amb = False
    p4 = tau * abs(lxt)
    if tau <= EPS_R4 and _near(p4, LN2):
        amb = True
    i_r4 = tau <= EPS_R4 and LN2 >= p4 and (xt <= XI1 or r2tf) and B1 >= bmax * xt
    r1n = x <= XI1 and B1 >= b * x
    r1s = y0 <= XI1 and B1 >= a * y0
    i_r1 = (r1n or r1s) and not i_r4
    with_inf = 1.0 / a + 1.0 / b
    nu = 1.0 / with_inf if with_inf > 0 else math.inf
    rat1 = x * (1.0 + b / a)
    rat2 = y0 * (1.0 + a / b)
    band = 0.5 <= rat1 <= 2.0 and 0.5 <= rat2 <= 2.0
    glr = bmax >= GAMMA_LIM and min(a, b) >= GL_RIDGE_MIN  # on s, not nu (#84)
    i_r3 = (nu >= T_RIDGE or glr) and band and not i_r4 and not i_r1
    if i_r4:
        out = [("R4", sw4)]
    elif i_r1:
        out = [("R1", (not r1n) and r1s)]
    elif i_r3:
        out = [("R3", rat1 > 1.0)]
    else:
        sw = not (x < thr)
        gl = bmax >= GAMMA_LIM and min(a, b) < GAMMA_LIM
        out = [("GL" if gl else "R2", sw)]
    if amb:
        out.append(("R1", False))
        out.append(("R2", not (x < thr)))
    return out


def bound_for(fn, a, b, x, v):
    """Integer ULP bound for beta_p/beta_q at (a, b, x) with exact value v."""
    with mp.workdps(60):
        P = v if fn == "beta_p" else 1 - v
        lam = mpf(a) * mp.fsub(1, mpf(x), exact=True) - mpf(b) * mpf(x)
    best = 0
    for reg, sw in route(a, b, x):
        cands = []
        if reg == "R1":
            alpha = b if sw else a
            routed = P if not sw else 1 - P
            if alpha <= PR_TAU_MAX and routed > NEAR_ONE * (1 - 2.0 ** -40):
                cands.append(("R4pr", sw))      # post-route: complement of routed
            if not (alpha <= PR_TAU_MAX and routed > NEAR_ONE * (1 + 2.0 ** -40)):
                cands.append(("R1", not sw))
        elif reg == "R4":
            cands.append(("R4", sw))           # R4 computes the complement
        elif reg == "R3":
            cands.append(("R3", lam >= 0))
        else:
            cands.append((reg, not sw))
        for r, direct_is_p in cands:
            d, c = BOUNDS[r]
            is_direct = direct_is_p == (fn == "beta_p")
            best = max(best, d if is_direct else c)
    return best


def region_tag(a, b, x):
    return route(a, b, x)[0][0]


# ---- oracles -------------------------------------------------------------------
def _oracle(fn):
    k = 0 if fn == "beta_p" else 1

    def f(args, dps):
        return betainc_pq(args[0], args[1], args[2], dps)[k]
    return f


def _oracle2(fn):
    k = 0 if fn == "beta_p" else 1

    def f(args, dps):
        return betainc_pq_cf(args[0], args[1], args[2], dps)[k]
    return f


def _bound(fn):
    return lambda args, v: bound_for(fn, args[0], args[1], args[2], v)


# ---- float helpers for point placement ----------------------------------------
def _phi(t):
    return t * t * (0.5 - t / 3 + t * t / 4) if abs(t) < 1e-4 else t - math.log1p(t)


def cpsi_f(a, b, x):
    lam = a * (1.0 - x) - b * x
    return a * _phi(-lam / a) + b * _phi(lam / b)


def x_for_cpsi(a, b, target, upper):
    """x on the lower (upper=False) or upper side of the mean with
    cpsi(a, b, x) ~ target (float bisection)."""
    m = a / (a + b)
    lo, hi = (m, 1.0) if upper else (0.0, m)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mid <= lo or mid >= hi:
            break
        c = cpsi_f(a, b, mid)
        if (c < target) == upper:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def shapes_nu_ratio(nu, r):
    """(a, b) with a/b = r and a*b/(a+b) = nu."""
    b = nu * (1.0 + r) / r
    return nu * (1.0 + r), b


def ln_gamma_p_small(s, t):
    """ln P(s,t) (t < s) or ln Q(s,t) (t > s), float leading forms."""
    if t < s:
        return s * math.log(t) - t - math.lgamma(s + 1) - math.log1p(-t / (s + 1))
    return (s - 1) * math.log(t) - t - math.lgamma(s) - math.log1p(-max(s - 1, 0) / t)


def t_for_lnP(s, target, upper):
    lo, hi = (s * 1.0001 + 1, 1e6) if upper else (1e-300, s * 0.9999)
    for _ in range(300):
        mid = math.exp(0.5 * (math.log(lo) + math.log(hi))) if not upper else 0.5 * (lo + hi)
        v = ln_gamma_p_small(s, mid)
        if (v > target) == upper:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def ln_I_left(a, b, x):
    """ln I_x(a,b) leading form for x below the mean (float; placement only)."""
    lnb = math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)
    r = min((a + b) * x / (a + 1), 0.999)
    return a * math.log(x) + b * math.log1p(-x) - math.log(a) - lnb - math.log1p(-r) * 0.5


def x_for_lnI(a, b, target):
    lo, hi = 1e-320, a / (a + b)
    for _ in range(300):
        mid = math.exp(0.5 * (math.log(lo) + math.log(hi)))
        if ln_I_left(a, b, mid) < target:
            lo = mid
        else:
            hi = mid
    return math.exp(0.5 * (math.log(lo) + math.log(hi)))


def _ok(a, b, x):
    return (math.isfinite(a) and math.isfinite(b) and a > 0 and b > 0
            and 0 < x < 1)


def _take(rng, pts, n):
    pts = [p for p in pts if _ok(*p)]
    if n < len(pts):
        pts = rng.sample(pts, n)
    return pts


def _fill(rng, n, draw, tries=40):
    out = []
    for _ in range(n * tries):
        if len(out) >= n:
            break
        p = draw(rng)
        if p is not None and _ok(*p):
            out.append(p)
    return out


def _sw(rng, a, b, x):
    """Randomly mirror (a, b, x) -> (b, a, 1-x) (exact in 1-x only when
    x >= 1/4 or tiny; used where the seam is orientation-symmetric)."""
    if rng.random() < 0.5:
        return a, b, x
    return b, a, 1.0 - x


# ---- generators: one per seam (bits, band) -------------------------------------
def g_xi1_x_bits(rng, n):
    combos = [(0.5, 3.0), (2.0, 11.0), (9.0, 17.5), (25.0, 1.5)]
    return _take(rng, [(a, b, x) for a, b in combos for x in neighbours(XI1, 64)], n)


def g_xi1_x_band(rng, n):
    return _fill(rng, n, lambda r: (log_uniform(r, 0.03, 60), log_uniform(r, 0.05, 17.6),
                                    full_mantissa(r, XI1 * 0.93, XI1 * 1.07)))


def g_xi1_y_bits(rng, n):
    combos = [(3.0, 0.5), (11.0, 2.0), (17.5, 9.0), (1.5, 25.0)]
    x0 = 1.0 - XI1
    return _take(rng, [(a, b, x) for a, b in combos for x in neighbours(x0, 64)], n)


def g_xi1_y_band(rng, n):
    def d(r):
        y = full_mantissa(r, XI1 * 0.93, XI1 * 1.07)
        return log_uniform(r, 0.05, 17.6), log_uniform(r, 0.03, 60), 1.0 - y
    return _fill(rng, n, d)


def g_b1_bits(rng, n):
    pts = []
    for x, a in ((0.05, 0.5), (0.2, 3.0), (0.44, 1.2), (0.01, 7.0)):
        for b in neighbours(B1 / x, 64):
            pts.append((a, b, x))
    return _take(rng, pts, n)


def g_b1_band(rng, n):
    def d(r):
        x = log_uniform(r, 1e-6, 0.45)
        b = full_mantissa(r, B1 / x * 0.93, B1 / x * 1.07)
        return _sw(r, log_uniform(r, 0.02, 40), b, x)
    return _fill(rng, n, d)


def g_epsr4_bits(rng, n):
    pts = []
    for b, x in ((0.3, 0.05), (2.0, 0.3), (7.0, 0.9), (0.9, 0.6)):
        for t in neighbours(EPS_R4, 64):
            pts.append((t, b, x))
    return _take(rng, [_sw(rng, *p) for p in pts], n)


def g_epsr4_band(rng, n):
    return _fill(rng, n, lambda r: _sw(r, full_mantissa(r, EPS_R4 * 0.9, EPS_R4 * 1.1),
                                       log_uniform(r, 0.05, 60), r.uniform(1e-4, 0.9999)))


def g_r4ln2_bits(rng, n):
    pts = []
    for tau, b in ((EPS_R4, 1.0), (0.005, 5.0), (0.0011, 0.4)):
        x0 = math.exp(-LN2 / tau)
        for x in neighbours(x0, 64):
            pts.append((tau, b, x))
    return _take(rng, pts, n)


def g_r4ln2_band(rng, n):
    def d(r):
        tau = log_uniform(r, 0.00095, EPS_R4)
        x = math.exp(-LN2 / tau * r.uniform(0.9, 1.1))
        return tau, log_uniform(r, 0.05, 30), x
    return _fill(rng, n, d)


def g_r2thr_bits(rng, n):
    pts = []
    for a, b in ((3.5, 7.25), (12.0, 5.0), (0.7, 30.0), (1e-3, 0.2), (0.15, 1e-30)):
        t = 1.0 / (1.0 + (b + 1.0) / (a + 1.0))
        for x in neighbours(t, 64):
            pts.append((a, b, x))
    return _take(rng, pts, n)


def g_r2thr_band(rng, n):
    def d(r):
        a, b = log_uniform(r, 0.3, 60), log_uniform(r, 0.3, 60)
        t = 1.0 / (1.0 + (b + 1.0) / (a + 1.0))
        return a, b, full_mantissa(r, t * 0.95, min(t * 1.05, 0.999))
    return _fill(rng, n, d)


def _near_one_point(r, a_lo, a_hi, rel):
    """R1-native lane whose value sits at 1 - 2^-11 * (1 +- rel)."""
    a = r.uniform(a_lo, a_hi)
    x = r.uniform(0.08, 0.449)
    target = 1.0 - (1.0 - NEAR_ONE) * (1.0 + r.uniform(-rel, rel))
    blo, bhi = a, B1 / x
    if approx_I(a, bhi, x) < target:
        return None
    for _ in range(80):
        bm = 0.5 * (blo + bhi)
        if approx_I(a, bm, x) < target:
            blo = bm
        else:
            bhi = bm
    return (a, 0.5 * (blo + bhi), x)


def g_nearone_band(rng, n):
    return _fill(rng, n, lambda r: _near_one_point(r, 0.02, 2.5, 0.05))


def g_nearone_bits(rng, n):
    pts = []
    r = random.Random(rng.getrandbits(32))
    while len(pts) < 4:
        p = _near_one_point(r, 0.05, 2.4, 0.0)
        if p:
            pts.append(p)
    return _take(rng, [(a, b, xx) for a, b, x in pts for xx in neighbours(x, 32)], n)


def _box_near_one(r, a):
    x = r.uniform(0.3, 0.449)
    b = B1 / x * r.uniform(0.75, 1.0)
    return a, b, x


def g_prtau_bits(rng, n):
    r = random.Random(rng.getrandbits(32))
    pts = []
    for _ in range(3):
        _, b, x = _box_near_one(r, 2.5)
        pts += [(a, b, x) for a in neighbours(PR_TAU_MAX, 64)]
    return _take(rng, pts, n)


def g_prtau_band(rng, n):
    return _fill(rng, n, lambda r: _box_near_one(r, full_mantissa(r, 2.3, 2.7)))


def g_lgtab_bits(rng, n):
    """#85: PA's lgamma(min) switches from the analytic table to LgammaPosDd
    at min = 16. a + b <= 216 keeps the lane on PA (C_lg = 256)."""
    r = random.Random(rng.getrandbits(32))
    pts = [(a, b, x) for a in neighbours(16.0, 24) for b, x in
           ((r.uniform(16.0, 200.0), r.uniform(0.05, 0.9)) for _ in range(2))]
    return _take(rng, pts, n)


def g_lgtab_band(rng, n):
    def d(r):
        return (full_mantissa(r, 15.2, 16.8), r.uniform(16.0, 200.0), r.uniform(0.02, 0.95))
    return _fill(rng, n, d)


def g_tau1_bits(rng, n):
    r = random.Random(rng.getrandbits(32))
    pts = []
    for _ in range(3):
        _, b, x = _box_near_one(r, 1.0)
        pts += [(a, b, x) for a in neighbours(1.0, 64)]
    return _take(rng, pts, n)


def g_tau1_band(rng, n):
    return _fill(rng, n, lambda r: _box_near_one(r, full_mantissa(r, 0.93, 1.07)))


def g_tridge_bits(rng, n):
    pts = []
    for r_ in (0.2, 1.0, 7.0):
        a, b = shapes_nu_ratio(T_RIDGE, r_)
        m = a / (a + b)
        for bb in neighbours(b, 43):
            for xm in (0.9, 1.05):
                pts.append((a, bb, m * xm))
    return _take(rng, pts, n)


def g_tridge_band(rng, n):
    def d(r):
        a, b = shapes_nu_ratio(T_RIDGE * r.uniform(0.95, 1.05), log_uniform(r, 0.05, 20))
        m = a / (a + b)
        return a, b, m * r.uniform(0.6, 1.6)
    return _fill(rng, n, d)


def _band_edge_point(r, which, rel):
    a, b = shapes_nu_ratio(log_uniform(r, 33, 1e5), log_uniform(r, 0.05, 20))
    p, q = a / (a + b), b / (a + b)
    f = 1.0 + r.uniform(-rel, rel)
    if which == "rat1_lo":
        x = 0.5 * p * f
    elif which == "rat1_hi":
        x = 2.0 * p * f
    elif which == "rat2_lo":
        x = 1.0 - 0.5 * q * f
    else:
        x = 1.0 - 2.0 * q * f
    return a, b, x


def g_band_bits(which):
    def g(rng, n):
        r = random.Random(rng.getrandbits(32))
        pts = []
        for _ in range(3):
            a, b, x = _band_edge_point(r, which, 0.0)
            if _ok(a, b, x):
                pts += [(a, b, xx) for xx in neighbours(x, 64)]
        return _take(rng, pts, n)
    return g


def g_band_band(which):
    return lambda rng, n: _fill(rng, n, lambda r: _band_edge_point(r, which, 0.05))


def _cpsi_point(r, lo, hi, nu_lo=33, nu_hi=1e6):
    a, b = shapes_nu_ratio(log_uniform(r, nu_lo, nu_hi), log_uniform(r, 0.03, 30))
    x = x_for_cpsi(a, b, r.uniform(lo, hi), r.random() < 0.5)
    return a, b, x


def g_z2_bits(rng, n):
    r = random.Random(rng.getrandbits(32))
    pts = []
    for _ in range(4):
        a, b, x = _cpsi_point(r, Z2_SPLIT, Z2_SPLIT)
        pts += [(a, b, xx) for xx in neighbours(x, 64)]
    return _take(rng, pts, n)


def g_cpsi_band(lo, hi, nu_lo=33, nu_hi=1e6):
    return lambda rng, n: _fill(rng, n, lambda r: _cpsi_point(r, lo, hi, nu_lo, nu_hi))


def g_r3sat_bits(rng, n):
    r = random.Random(rng.getrandbits(32))
    pts = []
    for c in (-EXP_FLOOR, 745.13, 708.4):
        a, b, x = _cpsi_point(r, c, c, 40, 1e4)
        pts += [(a, b, xx) for xx in neighbours(x, 40)]
    return _take(rng, pts, n)


def g_pswap_bits(rng, n):
    pts = []
    for bb in (100.0, 1e5, 1e12):
        for a in neighbours(bb, 64):
            for z in (0.0, 3.0):
                sd = math.sqrt(0.25 / (2 * bb))
                pts.append((a, bb, 0.5 + z * sd))
    return _take(rng, pts, n)


def g_pswap_band(rng, n):
    def d(r):
        b = log_uniform(r, 40, 1e9)
        a = b * r.uniform(0.95, 1.05)
        sd = math.sqrt(a * b / ((a + b) ** 3))
        return a, b, a / (a + b) + r.uniform(-6, 6) * sd
    return _fill(rng, n, d)


def _offband_x(r, a, b):
    m = a / (a + b)
    if r.random() < 0.5:
        return m * r.uniform(0.05, 0.45)
    return 1.0 - (1.0 - m) * r.uniform(0.05, 0.45)


def g_clg_bits(rng, n):
    r = random.Random(rng.getrandbits(32))
    pts = []
    for a in (10.0, 40.0, 128.0):
        b0 = C_LG - a
        x = _offband_x(r, a, b0)
        pts += [(a, b, x) for b in neighbours(b0, 64)]
    return _take(rng, pts, n)


def g_clg_band(rng, n):
    def d(r):
        c = full_mantissa(r, C_LG * 0.93, C_LG * 1.07)
        a = r.uniform(Z0, c - Z0)
        b = c - a
        return a, b, _offband_x(r, a, b)
    return _fill(rng, n, d)


def g_z0_bits(rng, n):
    r = random.Random(rng.getrandbits(32))
    pts = []
    for b in (300.0, 5e3, 1e7):
        x = _offband_x(r, Z0, b)
        pts += [_sw(r, a, b, x) for a in neighbours(Z0, 64)]
    return _take(rng, pts, n)


def g_z0_band(rng, n):
    def d(r):
        a = full_mantissa(r, Z0 * 0.93, Z0 * 1.07)
        b = log_uniform(r, 250, 1e9)
        return _sw(r, a, b, _offband_x(r, a, b))
    return _fill(rng, n, d)


def _gl_point(r, N, s, t):
    """(s, N, x) with t = -(N + (s-1)/2) log1p(-x) ~ t, either orientation."""
    x = t / (N + (s - 1) / 2)
    if x <= 0 or x >= 1:
        return None
    if r.random() < 0.5 or x < 2.0 ** -40:
        return s, N, x          # huge-first mirror needs 1-x to carry x's bits
    return N, s, 1.0 - x


def g_gl_bits(rng, n):
    r = random.Random(rng.getrandbits(32))
    pts = []
    for s, tf in ((0.5, 1.0), (3.0, 0.3), (25.0, 3.0), (400.0, 0.4), (3000.0, 1.3)):
        for N in neighbours(GAMMA_LIM, 40):
            x = s * tf / N
            pts.append((s, N, x) if (r.random() < 0.5 or x < 2.0 ** -40) else (N, s, 1.0 - x))
    return _take(rng, pts, n)


def g_gl_band(rng, n):
    def d(r):
        N = full_mantissa(r, GAMMA_LIM * 0.9, GAMMA_LIM * 1.1)
        s = log_uniform(r, 0.05, 3000)
        return _gl_point(r, N, s, s * log_uniform(r, 0.05, 6))
    return _fill(rng, n, lambda r: d(r))


def g_glridge_bits(rng, n):
    pts = []
    for N in (2.0 ** 47, 2.0 ** 60, 2.0 ** 200):
        s0 = GL_RIDGE_MIN * N / (N - GL_RIDGE_MIN)
        for s in neighbours(s0, 40):
            for tf in (0.7, 1.4):
                pts.append((s, N, s * tf / N))
    return _take(rng, pts, n)


def g_glridge_s20_bits(rng, n):
    """#84: the ridge is gated on s = min(a, b) at 20. Bit-stepped across
    s = 20 for N from 2^47 to 2^1000, at the t/s ratios the issue's rows
    used (strata_append keeps 48 rows of a :bits stratum)."""
    r = random.Random(rng.getrandbits(32))
    pts = []
    for N in (2.0 ** 47, 2.0 ** 48, 2.0 ** 60, 2.0 ** 200, 2.0 ** 1000):
        for s in neighbours(GL_RIDGE_MIN, 24):
            for tf in (0.7, 1.4):
                p = _gl_point(r, N, s, s * tf)
                if p is not None:
                    pts.append(p)
    return _take(rng, pts, n)


def g_issue84_sliver(rng, n):
    """#84: inside the sliver nu < 20 <= s, s in [20, 20N/(N-20)), where the
    s-gated router sends lanes to R3 and the nu-gated one sent them to the
    slice's CF. N where the sliver is several ulps of 20 wide."""
    r = random.Random(rng.getrandbits(32))
    pts = []
    for N in (2.0 ** 47, 2.0 ** 48, 2.0 ** 50):
        s1 = GL_RIDGE_MIN * N / (N - GL_RIDGE_MIN)
        for _ in range(16):
            s = r.uniform(GL_RIDGE_MIN, s1)
            for tf in (0.55, 0.7, 0.85, 1.2, 1.4, 1.9):
                p = _gl_point(r, N, s, s * tf)
                if p is not None:
                    pts.append(p)
    return _take(rng, pts, n)


def g_glridge_band(rng, n):
    def d(r):
        N = log_uniform(r, GAMMA_LIM, 1e100)
        s = full_mantissa(r, 18.5, 21.5)
        return _gl_point(r, N, s, s * r.uniform(0.55, 1.9))
    return _fill(rng, n, d)


def g_glentry_band(rng, n):
    """R3 entered through the gamma-limit ridge (#80): b ~ 1e13-1e16 either
    side of kBetaGammaLim, a ~ 20-5000, x in the R3 band. The worst beta_q
    R3 real error (1.356 FMA) sat here with no stratum aimed at it."""
    def d(r):
        N = log_uniform(r, 1e13, 1e16)
        s = log_uniform(r, GL_RIDGE_MIN, 5000.0)
        return _gl_point(r, N, s, s * r.uniform(0.55, 1.9))
    return _fill(rng, n, d)


def g_glsercf_bits(rng, n):
    pts = []
    for s in (0.5, 2.0, 7.0, 15.0):
        for N in (2.0 ** 48, 2.0 ** 100):
            x0 = (s + 1) / (N + (s - 1) / 2)
            pts += [(s, N, x) for x in neighbours(x0, 24)]
    return _take(rng, pts, n)


def g_glsercf_band(rng, n):
    def d(r):
        N = log_uniform(r, GAMMA_LIM, 1e200)
        s = log_uniform(r, 0.05, 19.9)
        return _gl_point(r, N, s, (s + 1) * r.uniform(0.9, 1.1))
    return _fill(rng, n, d)


def g_glat_bits(rng, n):
    pts = []
    for tf in (0.2, 0.45, 2.5, 4.0):
        for s in neighbours(GAMMA_AT, 32):
            N = 2.0 ** 52
            pts.append((s, N, s * tf / N))
    return _take(rng, pts, n)


def g_glat_band(rng, n):
    def d(r):
        N = log_uniform(r, GAMMA_LIM, 1e150)
        s = full_mantissa(r, 18.6, 21.4)
        tf = r.choice((r.uniform(0.05, 0.49), r.uniform(2.05, 6.0)))
        return _gl_point(r, N, s, s * tf)
    return _fill(rng, n, d)


def _glsat_point(r, lnlo, lnhi):
    s = log_uniform(r, 1.0, 3000)
    upper = r.random() < 0.5
    t = t_for_lnP(s, r.uniform(lnlo, lnhi), upper)
    N = log_uniform(r, GAMMA_LIM, 1e250)
    if t / N < 1e-300:
        return None
    return _gl_point(r, N, s, t)


def g_glsat_band(lo, hi):
    return lambda rng, n: _fill(rng, n, lambda r: _glsat_point(r, lo, hi))


def g_scale_bits(rng, n):
    pts = []
    for s, tf in ((2.0, 0.3), (50.0, 3.0), (5.0, 1.0)):
        for N in neighbours(SCALE_ABOVE, 32):
            pts.append((s, N, s * tf / N))
    for rr in (1.0, 3.0):
        for a in neighbours(SCALE_ABOVE, 16):
            b = a * rr
            pts.append((a, b, a / (a + b)))
    return _take(rng, pts, n)


def g_scale_band(rng, n):
    def d(r):
        N = full_mantissa(r, SCALE_ABOVE * 0.9, SCALE_ABOVE * 1.1)
        s = log_uniform(r, 0.1, 5000)
        return _gl_point(r, N, s, s * log_uniform(r, 0.1, 5))
    return _fill(rng, n, d)


def g_divfloor_bits(rng, n):
    pts = []
    for tau, x in ((2.0 ** -1000, 0.1), (1e-310, 0.5), (2.0 ** -905, 0.9)):
        for M in neighbours(DIV_FLOOR, 40):
            pts.append(_sw(rng, tau, M, x))
    return _take(rng, pts, n)


def g_divfloor_band(rng, n):
    def d(r):
        M = full_mantissa(r, DIV_FLOOR * 0.9, DIV_FLOOR * 1.1)
        tau = log_uniform(r, 5e-324, M)
        return _sw(r, tau, M, r.uniform(1e-3, 0.999))
    return _fill(rng, n, d)


def g_tausubn_bits(rng, n):
    pts = []
    for b, x in ((0.5, 0.05), (3.0, 0.6), (2.0 ** -920, 0.97)):
        for t in neighbours(TAU_SUBN, 40):
            pts.append(_sw(rng, t, b, x))
    return _take(rng, pts, n)


def g_tausubn_band(rng, n):
    def d(r):
        t = full_mantissa(r, TAU_SUBN * 0.9, TAU_SUBN * 1.1)
        b = log_uniform(r, 1e-300, 7.9)
        x = r.uniform(1e-3, min(0.999, 8.0 / b))
        return _sw(r, t, b, x)
    return _fill(rng, n, d)


def g_btiny_bits(rng, n):
    pts = []
    for tau, x in ((1e-200, 0.3), (2.0 ** -145, 0.7), (1e-160, 0.5)):
        for M in neighbours(B_TINY_CUT, 40):
            pts.append(_sw(rng, tau, M, x))
    return _take(rng, pts, n)


def g_btiny_band(rng, n):
    def d(r):
        M = full_mantissa(r, B_TINY_CUT * 0.9, B_TINY_CUT * 1.1)
        tau = log_uniform(r, 1e-300, M)
        return _sw(r, tau, M, r.uniform(1e-3, 0.999))
    return _fill(rng, n, d)


def g_twopinu_bits(rng, n):
    pts = []
    a0 = 2.0 * TWO_PI_NU_CLAMP
    for a in neighbours(a0, 32):
        pts.append((a, a, 0.5))
        pts.append((a, a0, 0.5))
    for x in neighbours(0.5, 8):
        pts.append((a0, a0, x))
    return _take(rng, pts, n)


def g_binet_bits(rng, n):
    pts = []
    for s, t in ((0.5, 0.3), (5.0, 3.0), (2.0, 7.9)):
        for N in neighbours(BINET_Z_CLAMP, 40):
            pts.append((s, N, t / N))
    return _take(rng, pts, n)


def g_binet_band(rng, n):
    def d(r):
        N = full_mantissa(r, BINET_Z_CLAMP * 0.9, BINET_Z_CLAMP * 1.1)
        s = log_uniform(r, 0.01, 9.9)
        return _gl_point(r, N, s, log_uniform(r, 0.01, 7.9))
    return _fill(rng, n, d)


# ---- subnormal strata ------------------------------------------------------------
def g_sub_a(rng, n):
    return _fill(rng, n, lambda r: _sw(r, subnormal(r), log_uniform(r, 1e-3, 1e4),
                                       r.choice((r.random(), log_uniform(r, 1e-300, 0.5)))))


def g_sub_x(rng, n):
    return _fill(rng, n, lambda r: _sw(r, log_uniform(r, 1e-3, 1.5), log_uniform(r, 0.05, 1e6),
                                       subnormal(r)))


def g_sub_pair(rng, n):
    return _fill(rng, n, lambda r: (subnormal(r), subnormal(r), r.random()))


def _subout_point(r, alo, ahi, blo, bhi):
    a, b = log_uniform(r, alo, ahi), log_uniform(r, blo, bhi)
    tgt = r.uniform(-1074, -1000) * math.log(2)
    x = x_for_lnI(a, b, tgt)
    return _sw(r, a, b, x) if x > 0 else None


def g_subout(alo, ahi, blo, bhi):
    return lambda rng, n: _fill(rng, n, lambda r: _subout_point(r, alo, ahi, blo, bhi))


def g_subout_r4(rng, n):
    def d(r):
        tau = log_uniform(r, 5e-324, 1e-300)
        b = log_uniform(r, 0.05, 5.0)
        return _sw(r, tau, b, r.uniform(0.01, 0.9))
    return _fill(rng, n, d)


# ---- dense populations ----------------------------------------------------------
def g_r4_dense_b8(rng, n):
    def d(r):
        tau = log_uniform(r, 1e-30, EPS_R4)
        B = r.uniform(4.0, 16.0)
        lo = math.exp(-LN2 / tau)
        hi = min(XI1, B1 / B)
        if lo >= hi:
            return None
        return _sw(r, tau, B, log_uniform(r, max(lo, 1e-300), hi))
    return _fill(rng, n, d)


def g_r2_comp_mid(rng, n):
    def d(r):
        a, b = log_uniform(r, 3, 40), log_uniform(r, 3, 40)
        if a * b / (a + b) >= T_RIDGE:
            return None
        sd = math.sqrt(a * b / ((a + b) ** 2 * (a + b + 1)))
        return a, b, a / (a + b) + r.uniform(-3, 3) * sd
    return _fill(rng, n, d)


def g_witness(rng, n):
    w = (float.fromhex("0x1.7ebfd5597637bp+4"), float.fromhex("0x1.390d72284b741p+4"),
         float.fromhex("0x1.1e9d7495b0b32p-1"))
    return [w] + [(w[0], w[1], x) for x in neighbours(w[2], 16) if x != w[2]]


# ---- assembly -------------------------------------------------------------------
def _pair(name, gen, n, seam=None, note=""):
    out = []
    for fn in ("beta_p", "beta_q"):
        out.append(Stratum(fn, name, gen, n, _oracle(fn), _bound(fn),
                           oracle2=_oracle2(fn), seam=seam, note=note))
    return out


def strata():
    S = []
    BITS, BAND = 258, 500
    S += _pair("seam:xi1_x:bits", g_xi1_x_bits, 4 * 129, "xi1_x")
    S += _pair("seam:xi1_x:band", g_xi1_x_band, BAND, "xi1_x")
    S += _pair("seam:xi1_y:bits", g_xi1_y_bits, 4 * 129, "xi1_y")
    S += _pair("seam:xi1_y:band", g_xi1_y_band, BAND, "xi1_y")
    S += _pair("seam:b1_bx:bits", g_b1_bits, 4 * 129, "b1_bx")
    S += _pair("seam:b1_bx:band", g_b1_band, BAND, "b1_bx")
    S += _pair("seam:epsr4_min:bits", g_epsr4_bits, 4 * 129, "epsr4_min")
    S += _pair("seam:epsr4_min:band", g_epsr4_band, BAND, "epsr4_min")
    S += _pair("seam:r4_ln2:bits", g_r4ln2_bits, 3 * 129, "r4_ln2")
    S += _pair("seam:r4_ln2:band", g_r4ln2_band, BAND, "r4_ln2")
    S += _pair("seam:r2_thr:bits", g_r2thr_bits, 5 * 129, "r2_thr", "complement side is a 0-ULP cell")
    S += _pair("seam:r2_thr:band", g_r2thr_band, 800, "r2_thr", "complement side is a 0-ULP cell")
    S += _pair("seam:near_one:bits", g_nearone_bits, 4 * 65, "near_one")
    S += _pair("seam:near_one:band", g_nearone_band, BAND, "near_one")
    S += _pair("seam:pr_tau_max:bits", g_prtau_bits, 3 * 129, "pr_tau_max")
    S += _pair("seam:pr_tau_max:band", g_prtau_band, BAND, "pr_tau_max")
    S += _pair("seam:lgamma_table_max:bits", g_lgtab_bits, 49 * 2, "lgamma_table_max")
    S += _pair("seam:lgamma_table_max:band", g_lgtab_band, BAND, "lgamma_table_max")
    S += _pair("seam:r4_tau1:bits", g_tau1_bits, 3 * 129, "r4_tau1")
    S += _pair("seam:r4_tau1:band", g_tau1_band, BAND, "r4_tau1")
    S += _pair("seam:tridge_nu:bits", g_tridge_bits, 3 * 87 * 2, "tridge_nu")
    S += _pair("seam:tridge_nu:band", g_tridge_band, BAND, "tridge_nu")
    for w in ("rat1_lo", "rat1_hi", "rat2_lo", "rat2_hi"):
        S += _pair(f"seam:band_{w}:bits", g_band_bits(w), 3 * 129, f"band_{w}")
        S += _pair(f"seam:band_{w}:band", g_band_band(w), 300, f"band_{w}")
    S += _pair("seam:z2split:bits", g_z2_bits, 4 * 129, "z2split")
    S += _pair("seam:z2split:band", g_cpsi_band(23.75, 26.25), BAND, "z2split")
    S += _pair("dense:r3_old_split", g_cpsi_band(33.0, 39.0), BAND, None,
               "around the pre-#70 split cpsi = 36, now tail lanes")
    S += _pair("seam:z2split:core_below", g_cpsi_band(20.0, 36.0), 1500, "z2split",
               "R3 Temme core just below cpsi=36 (gamma's twin 3 ULP, #70)")
    S += _pair("seam:z2split:core_big_nu", g_cpsi_band(25.0, 36.0, 1e6, 1e15), 400, "z2split")
    S += _pair("seam:r3_sat:bits", g_r3sat_bits, 3 * 81, "r3_sat")
    S += _pair("seam:r3_sat:subnormal_out", g_cpsi_band(700.0, 750.0, 40, 1e5), 800, "r3_sat",
               "R3 results in the subnormal band")
    S += _pair("seam:r3_pswap:bits", g_pswap_bits, 3 * 129 * 2, "r3_pswap")
    S += _pair("seam:r3_pswap:band", g_pswap_band, BAND, "r3_pswap")
    S += _pair("seam:clg_c:bits", g_clg_bits, 3 * 129, "clg_c")
    S += _pair("seam:clg_c:band", g_clg_band, BAND, "clg_c")
    S += _pair("seam:z0_min:bits", g_z0_bits, 3 * 129, "z0_min")
    S += _pair("seam:z0_min:band", g_z0_band, BAND, "z0_min")
    S += _pair("seam:gammalim:bits", g_gl_bits, 5 * 81, "gammalim", "1-ULP cells since #90")
    S += _pair("seam:gammalim:band", g_gl_band, 800, "gammalim", "1-ULP cells since #90")
    S += _pair("seam:gl_ridge:bits", g_glridge_bits, 3 * 81 * 2, "gl_ridge")
    S += _pair("seam:gl_ridge:band", g_glridge_band, BAND, "gl_ridge")
    S += _pair("seam:gl_ridge_s20:bits", g_glridge_s20_bits, 5 * 49 * 2, "gl_ridge",
               "s across 20, N 2^47..2^1000 (#84)")
    # #85 witnesses: min shape just past kBetaPrTauMax with Q ~ 2^-10.5, so no
    # post-route; the R1 complement amplified LgammaPosDd's 2^-62 absolute
    # error to 2 ULP on SSE4/SSE2. Bit neighbourhoods in x, both witnesses.
    S += _pair("issue85:bits", lambda r, n: _take(r, [
        (a, b, xx)
        for a, b, x0 in ((float.fromhex("0x1.4000000000040p+1"), float.fromhex("0x1.16d2f94b501a8p+4"),
                          float.fromhex("0x1.bbb1361c00740p-2")),
                         (float.fromhex("0x1.4497c15e49c2bp+1"), float.fromhex("0x1.1d9cd393a33dep+4"),
                          float.fromhex("0x1.c0256ad5870c5p-2")))
        for xx in neighbours(x0, 24)], n), 2 * 49, "pr_tau_max", "#85 witness neighbourhoods")
    S += _pair("issue84:sliver", g_issue84_sliver, 3 * 16 * 6, "gl_ridge",
               "nu < 20 <= s: the slice took the CF here (#84)")
    S += _pair("ridge:gl_entry", g_glentry_band, 1500, "gl_ridge",
               "R3 via the gamma-limit ridge, b ~ 1e13-1e16 (#80)")
    S += _pair("seam:gl_ser_cf:bits", g_glsercf_bits, 8 * 49, "gl_ser_cf")
    S += _pair("seam:gl_ser_cf:band", g_glsercf_band, BAND, "gl_ser_cf")
    S += _pair("seam:gl_at:bits", g_glat_bits, 4 * 65, "gl_at")
    S += _pair("seam:gl_at:band", g_glat_band, BAND, "gl_at")
    S += _pair("seam:gl_sat:subnormal_out", g_glsat_band(-745.1, -708.4), 800, "gl_sat",
               "gamma-limit slice, subnormal results")
    S += _pair("seam:gl_sat:clamp", g_glsat_band(-830.0, -770.0), 300, "gl_sat")
    S += _pair("seam:scale_above:bits", g_scale_bits, 3 * 65 + 2 * 33, "scale_above",
               "dps ~300 oracle rows")
    S += _pair("seam:scale_above:band", g_scale_band, 200, "scale_above")
    S += _pair("seam:div_floor:bits", g_divfloor_bits, 3 * 81, "div_floor")
    S += _pair("seam:div_floor:band", g_divfloor_band, 400, "div_floor")
    S += _pair("seam:tau_subn:bits", g_tausubn_bits, 3 * 81, "tau_subn")
    S += _pair("seam:tau_subn:band", g_tausubn_band, 400, "tau_subn")
    S += _pair("seam:b_tiny_cut:bits", g_btiny_bits, 3 * 81, "b_tiny_cut")
    S += _pair("seam:b_tiny_cut:band", g_btiny_band, 400, "b_tiny_cut")
    S += _pair("seam:two_pi_nu_clamp:bits", g_twopinu_bits, 2 * 65 + 17, "two_pi_nu_clamp")
    S += _pair("seam:binet_z_clamp:bits", g_binet_bits, 3 * 81, "binet_z_clamp")
    S += _pair("seam:binet_z_clamp:band", g_binet_band, 300, "binet_z_clamp")
    # subnormal inputs / outputs
    S += _pair("subnormal:a_or_b", g_sub_a, 1200)
    S += _pair("subnormal:x", g_sub_x, 1000)
    S += _pair("subnormal:a_and_b", g_sub_pair, 600)
    S += _pair("subnormal_out:R1", g_subout(0.3, 30, 0.05, 30), 1000)
    S += _pair("subnormal_out:PB", g_subout(30, 3000, 30, 3000), 800)
    S += _pair("subnormal_out:R4", g_subout_r4, 600)
    # dense populations
    S += _pair("dense:r4_tiny_min_B8", g_r4_dense_b8, 3000, None,
               "R4 tiny-min corner around B ~ 8 (known one-signed -1.60 ulp)")
    S += _pair("dense:r2_complement_mid", g_r2_comp_mid, 3000, "r2_thr",
               "R2 moderate shapes around the mean; complement 0-ULP cell (#76)")
    S += _pair("witness:r2_complement", g_witness, 33, None,
               "beta_p(0x1.7ebfd5597637bp+4, 0x1.390d72284b741p+4, 0x1.1e9d7495b0b32p-1) 1 ULP in a 0 cell")
    # T2 completion (#82): two_pi_nu_clamp had bits only. a = b = 2 nu at
    # x = 1/2 puts lambda at 0, so the band stays out of saturation.
    S += _pair("seam:two_pi_nu_clamp:band",
               lambda rng, n: [(lambda a: (a, a, 0.5))(2.0 * full_mantissa(rng, TWO_PI_NU_CLAMP * 0.95,
                                                                             TWO_PI_NU_CLAMP * 1.05))
                               for _ in range(n)], 200, "two_pi_nu_clamp")
    # #90 witnesses: gamma-limit near-tie misses (1 ULP, 0.501-0.504 real, all
    # three tiers) in cells that were gated 0 on dev/v1.1.0; v1.0.x published 1.
    S += _pair("issue90:witness", lambda r, n: [tuple(float.fromhex(v) for v in w) for w in (
        ("0x1.e000000000000p+3", "0x1.0000000000000p+48", "0x1.fffffffffff14p-45"),
        ("0x1.93d7ee129117cp+3", "0x1.8344b91246540p+546", "0x1.36fc9942edc90p-543"),
        ("0x1.47e1175a8ef1cp+4", "0x1.ab5ee2ef0bebap+230", "0x1.4d5ca4b82e4dep-225"),
        ("0x1.3614b84676e08p+4", "0x1.84103d08cbd29p+467", "0x1.4c55bfa7da863p-462"),
        ("0x1.3668a40581ecdp+4", "0x1.93c3031cebd06p+120", "0x1.1d3ff08c79bb7p-114"),
        ("0x1.fdb4fb2d48518p+3", "0x1.03d33000fec35p+650", "0x1.065c045774ee0p-646"))],
        6, "gl_at", "#90 gamma-limit near-tie witnesses")
    return S


# ---- seam-audit fields (#82 T2; tools/hunt/seam_audit.py docstring) ----------
def _regions(r):
    try:
        return dict(route(r["a"], r["b"], r["x"]))
    except (ZeroDivisionError, ValueError, OverflowError):
        return {}


def _nu(r):
    return 1.0 / (1.0 / r["a"] + 1.0 / r["b"])


def _r4_window(r):
    a, b, x = r["a"], r["b"], r["x"]
    tau = min(a, b)
    xt = 1.0 - x if b < a else x
    return tau * abs(math.log(xt)) if tau <= EPS_R4 and 0 < xt < 1 else None


def _r1_frame(r):
    """(alpha, R1 value) in R1's own frame, for R1 rows only."""
    reg = _regions(r)
    if "R1" not in reg:
        return None
    a, b, x = r["a"], r["b"], r["x"]
    return (b, approx_I(b, a, 1.0 - x)) if reg["R1"] else (a, approx_I(a, b, x))


def _gl(r):
    """(s, N, x_s) for gamma-limit-slice rows: s the small shape."""
    a, b, x = r["a"], r["b"], r["x"]
    if not (max(a, b) >= GAMMA_LIM > min(a, b)):
        return None
    return (a, b, x) if a < b else (b, a, 1.0 - x)


def _gl_t(r):
    g = _gl(r)
    if g is None or g[0] >= GAMMA_AT:
        return None
    s, N, xs = g
    return xs * (N + (s - 1.0) / 2.0) / (s + 1.0)


_AUDIT = {
    "r4_ln2": {"audit_f": _r4_window},
    "r2_thr": {"audit_f": lambda r: r["x"] / ((r["a"] + 1.0) / (r["a"] + r["b"] + 2.0)),
               "audit_value": 1.0},
    "near_one": {"audit_f": lambda r: (_r1_frame(r) or (None, None))[1]},
    "pr_tau_max": {"audit_f": lambda r: (_r1_frame(r) or (None, None))[0]},
    "lgamma_table_max": {"audit_f": lambda r: min(r["a"], r["b"])},
    "r4_tau1": {"audit_f": lambda r: (_r1_frame(r) or (None, None))[0]},
    "tridge_nu": {"audit_f": _nu},
    "two_pi_nu_clamp": {"audit_f": _nu},
    "band_rat1_lo": {"audit_f": lambda r: r["x"] * (r["a"] + r["b"]) / r["a"]},
    "band_rat1_hi": {"audit_f": lambda r: r["x"] * (r["a"] + r["b"]) / r["a"]},
    "band_rat2_lo": {"audit_f": lambda r: (1.0 - r["x"]) * (r["a"] + r["b"]) / r["b"]},
    "band_rat2_hi": {"audit_f": lambda r: (1.0 - r["x"]) * (r["a"] + r["b"]) / r["b"]},
    "z2split": {"audit_f": lambda r: cpsi_f(r["a"], r["b"], r["x"])},
    "r3_sat": {"audit_f": lambda r: cpsi_f(r["a"], r["b"], r["x"])},
    "gl_ridge": {"audit_f": lambda r: min(r["a"], r["b"]) if max(r["a"], r["b"]) >= GAMMA_LIM else None},
    "gl_ser_cf": {"audit_f": _gl_t, "audit_value": 1.0},
    "gl_at": {"audit_f": lambda r: (_gl(r) or (None,))[0]},
    "gl_sat": {"route": "beta/GL/saturated"},
    "b_tiny_cut": {"audit_f": lambda r: max(r["a"], r["b"]) if "R4" in _regions(r) else None},
}
for _s in SEAMS:
    _s.update(_AUDIT.get(_s["name"], {}))
