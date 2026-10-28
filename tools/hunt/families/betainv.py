"""beta_p_inv / beta_q_inv hunt strata (#82).

Bounds (docs/ACCURACY.md, beta_p_inv and beta_q_inv): 1 ULP in every regime
(deep-small, seeds, ridge, gamma-limit transfer, subnormal and x = 1
results); 2 ULP from the huge-nu onset nu = ab/(a+b) >= 1e31 (beyond-
resolution rows, neighbour semantics). Rows whose condition number
kappa = s/(xi f(x)) exceeds 2^18 (s = min(p, 1-p), xi = min(x, 1-x)) carry a
BACKWARD-error contract, which this harness cannot score: they get bound
2^40 (never flagged). Generators keep such rows rare (no both-tiny pairs).

Oracle: the exact root by safeguarded Newton in u = ln(xi) against the
bounded forward (_beta_oracles.betainc_pq), solved in the frame whose root
is <= 1/2 and on the side whose target is <= 1/2, so roots near 1 are
judged in x's own ulps (x = 1 - xi formed at dps + 30). oracle2: one
Newton correction of that root through the CF forward (betainc_pq_cf), so
a violation is certified against a second forward form.

Placement uses float approximations of P (approx_I, normal form for huge
shapes): rows land NEAR a seam, not on it bit-exactly, except where the
seam is a p-coordinate.
"""
from __future__ import annotations

import math
import random

import mpmath as mp
from mpmath import mpf

from core import Stratum, full_mantissa, log_uniform, neighbours, subnormal
from ._beta_oracles import (approx_I, beta_inv_root, betainc_pq_cf, ln_beta)
from .beta import (B1, GAMMA_LIM, GL_RIDGE_MIN, NEAR_ONE, XI1, ln_I_left,
                   cpsi_f, shapes_nu_ratio, x_for_cpsi)

DEEP_SMALL_CUT = 2.0 ** -60          # kBetaInvDeepSmallCut (betainv_data.h)
S3_MIN = 2.0 ** 20                   # kBetaInvS3Min
S1_NU_MIN = 2.0                      # kBetaInvS1NuMin
TIE_BAND = 2.0 ** -40                # kBetaInvTieBand (orientation probe)
PA_SCALE_MAX = 2.0 ** 40             # kBetaInvPaScaleMax
HUGE_NU = 1e31                       # huge-nu band onset (ACCURACY.md)
KAPPA_MAX = 2.0 ** 18

SEAMS = [
    {"fn": "betainv", "name": "fwd_xi1", "coord": "root x", "value": XI1, "source": "kBetaXi1 (forward router twin)"},
    {"fn": "betainv", "name": "fwd_b1", "coord": "b*x", "value": B1, "source": "kBetaB1 (forward router twin)"},
    {"fn": "betainv", "name": "fwd_r2_thr", "coord": "root x", "value": "(a+1)/(a+b+2)", "source": "R2 orientation (forward twin)"},
    {"fn": "betainv", "name": "fwd_z2split", "coord": "cpsi(root)", "value": 25.0, "source": "kBetaZ2Split (forward core/tail; 36 until #70)"},
    {"fn": "betainv", "name": "fwd_gammalim", "coord": "max(a,b)", "value": GAMMA_LIM, "source": "kBetaGammaLim (forward twin)"},
    {"fn": "betainv", "name": "fwd_gl_ridge", "coord": "min(a,b) (max>=2^47)", "value": GL_RIDGE_MIN, "source": "kBetaGlRidgeMin, gated on s since #84"},
    {"fn": "betainv", "name": "deep_small", "coord": "|1-b| y/(1+a) corr(y)", "value": DEEP_SMALL_CUT, "source": "kBetaInvDeepSmallCut"},
    {"fn": "betainv", "name": "s3_min", "coord": "max(a,b)", "value": S3_MIN, "source": "kBetaInvS3Min (gamma-limit seed gate)"},
    {"fn": "betainv", "name": "s1_nu_min", "coord": "nu", "value": S1_NU_MIN, "source": "kBetaInvS1NuMin (ridge seed gate)"},
    {"fn": "betainv", "name": "tie_band", "coord": "p - P(1/2)", "value": TIE_BAND, "source": "kBetaInvTieBand (orientation probe)"},
    {"fn": "betainv", "name": "near_one_log", "coord": "ln f", "value": -4.8840049810887e-4, "source": "kBetaInvNearOneLog"},
    {"fn": "betainv", "name": "pa_scale", "coord": "max(|E|,|lnB|)", "value": PA_SCALE_MAX, "source": "kBetaInvPaScaleMax"},
    {"fn": "betainv", "name": "huge_nu", "coord": "nu", "value": HUGE_NU, "source": "huge-nu band onset (ACCURACY.md)"},
    {"fn": "betainv", "name": "lnsigc_log1p", "coord": "sigma = min(s, 1-s)", "value": 2.0 ** -7, "source": "ln(1 - sigma) as Log1pDd below 2^-7 (#88)"},
    {"fn": "betainv", "name": "ds_scaled_a", "coord": "alpha (deep-small frame)", "value": 2.0 ** -900, "source": "kGammaInvTinyA: deep-small scaled numerator (#88)"},
    {"fn": "betainv", "name": "lgamma1p_table_max", "coord": "min(a,b)", "value": 15.0, "source": "BetaInvLgamma1p table form on (2.5, 15] (#85)"},
    {"fn": "betainv", "name": "ds_scaled_b", "coord": "beta (deep-small frame)", "value": 2.0 ** -500, "source": "deep-small scaled numerator needs beta >= 2^-500 (#88)"},
]

FNS = ("beta_p_inv", "beta_q_inv")


# ---- oracle / bound ----------------------------------------------------------
def _side(fn):
    return "p" if fn == "beta_p_inv" else "q"


def _oracle(fn):
    side = _side(fn)
    return lambda args, dps: beta_inv_root(args[0], args[1], args[2], side, dps)


def _oracle2(fn):
    side = _side(fn)

    def f(args, dps):
        a, b, p = args
        r = beta_inv_root(a, b, p, side, dps)
        with mp.workdps(dps + 30):
            if r <= 0 or r >= 1:
                return r
            P, Q = betainc_pq_cf(a, b, r, dps)
            lnB = ln_beta(a, b, dps + 10)
            dens = mp.exp((a - 1) * mp.log(r) + (b - 1) * mp.log1p(-r) - lnB)
            resid = (P - p) if side == "p" else (p - Q)
            return r - resid / dens
    return f


def _bound(fn):
    side = _side(fn)

    def bnd(args, v):
        a, b, p = args
        nu = a * b / (a + b) if a + b > 0 else 0.0
        base = 2 if nu >= HUGE_NU else 1
        if v <= 0 or v >= 1:
            return base
        wd = 25 + int(math.log10(a + b + 1))   # lnf cancels ~ (a+b) digits
        with mp.workdps(wd):
            s = min(mpf(p), 1 - mpf(p))
            xi = min(v, 1 - v)
            lnB = ln_beta(a, b, wd)
            lnf = (a - 1) * mp.log(v) + (b - 1) * mp.log1p(-v) - lnB
            lk = mp.log(s) - mp.log(xi) - lnf
        if lk > math.log(KAPPA_MAX):
            return 1 << 40
        return base
    return bnd


# ---- placement helpers ----------------------------------------------------------
def _target(r, a, b, x, fn):
    """p (beta_p_inv) or q (beta_q_inv) for root ~x, float; None if out of range."""
    if not (0 < x < 1):
        return None
    P = approx_I(a, b, x)
    if fn == "beta_p_inv":
        t = P if P < 0.5 else 1.0 - approx_I(b, a, 1.0 - x)
    else:
        t = 1.0 - P if P < 0.5 else approx_I(b, a, 1.0 - x)
    if not (0 < t < 1) or not math.isfinite(t):
        return None
    return t


def _ok(a, b, p):
    return (math.isfinite(a) and math.isfinite(b) and a > 0 and b > 0
            and 0 < p < 1)


def _fill(fn, rng, n, draw, tries=40):
    out = []
    for _ in range(n * tries):
        if len(out) >= n:
            break
        q = draw(rng)
        if q is None:
            continue
        a, b, x = q
        t = _target(rng, a, b, x, fn)
        if t is not None and _ok(a, b, t):
            out.append((a, b, t))
    return out


def _fill_p(rng, n, draw, tries=40):
    """draw returns (a, b, p) directly."""
    out = []
    for _ in range(n * tries):
        if len(out) >= n:
            break
        q = draw(rng)
        if q is not None and _ok(*q):
            out.append(q)
    return out


def _take(rng, pts, n):
    pts = [p for p in pts if _ok(*p)]
    if n < len(pts):
        pts = rng.sample(pts, n)
    return pts


def _flip(r, a, b, x):
    return (a, b, x) if r.random() < 0.5 else (b, a, 1.0 - x)


# ---- generators -------------------------------------------------------------------
def d_xi1(r):
    return (log_uniform(r, 0.3, 60), log_uniform(r, 0.3, 17.6),
            full_mantissa(r, XI1 * 0.95, XI1 * 1.05))


def d_b1(r):
    x = log_uniform(r, 1e-4, 0.45)
    return _flip(r, log_uniform(r, 0.3, 40), B1 / x * r.uniform(0.95, 1.05), x)


def d_thr(r):
    a, b = log_uniform(r, 0.5, 60), log_uniform(r, 0.5, 60)
    t = (a + 1) / (a + b + 2)
    return a, b, t * r.uniform(0.95, 1.05)


def d_cpsi(lo, hi, nlo=33, nhi=1e6):
    def d(r):
        a, b = shapes_nu_ratio(log_uniform(r, nlo, nhi), log_uniform(r, 0.03, 30))
        return a, b, x_for_cpsi(a, b, r.uniform(lo, hi), r.random() < 0.5)
    return d


def d_gl(r):
    N = full_mantissa(r, GAMMA_LIM * 0.9, GAMMA_LIM * 1.1)
    s = log_uniform(r, 0.05, 3000)
    x = s * log_uniform(r, 0.05, 6) / N
    return (s, N, x) if r.random() < 0.5 else (N, s, 1.0 - x)


def d_glridge(r):
    N = log_uniform(r, GAMMA_LIM, 1e60)
    s = full_mantissa(r, 18.5, 21.5)
    return s, N, s * r.uniform(0.55, 1.9) / N


def d_deep_small(r):
    al = log_uniform(r, 0.02, 12)
    be = log_uniform(r, 0.02, 1e5)
    if abs(1 - be) < 1e-3:
        return None
    y = DEEP_SMALL_CUT * (1 + al) / abs(1 - be) * r.uniform(0.9, 1.1)
    if not 0 < y < 0.5:
        return None
    return _flip(r, al, be, y)


def d_s3(r):
    M = full_mantissa(r, S3_MIN * 0.9, S3_MIN * 1.1)
    s = log_uniform(r, 0.05, 200)
    x = s * log_uniform(r, 0.05, 8) / M
    return (s, M, x) if r.random() < 0.5 else (M, s, 1.0 - x)


def d_s1(r):
    a, b = shapes_nu_ratio(S1_NU_MIN * r.uniform(0.95, 1.05), log_uniform(r, 0.05, 20))
    m = a / (a + b)
    return a, b, min(0.999, m * log_uniform(r, 0.01, 3))


def d_near_one(fn):
    def d(r):
        a, b = log_uniform(r, 0.3, 300), log_uniform(r, 0.3, 300)
        t = (1.0 - NEAR_ONE) * r.uniform(0.95, 1.05)
        return (a, b, t) if fn == "beta_q_inv" else (a, b, 1.0 - t)
    return d


def d_pa_scale(fn):
    def d(r):
        a, b = shapes_nu_ratio(log_uniform(r, 1e10, 1e13), log_uniform(r, 0.1, 10))
        z = r.uniform(0, 37)
        t = 0.5 * math.erfc(z / math.sqrt(2))
        return a, b, t if r.random() < 0.5 else 1.0 - t
    return d


def d_huge_nu(lo, hi, tail):
    def d(r):
        a, b = shapes_nu_ratio(log_uniform(r, lo, hi), log_uniform(r, 0.01, 100))
        if tail:
            t = log_uniform(r, 1e-300, 1e-3)
        else:
            t = r.uniform(0.001, 0.999)
        return a, b, t if r.random() < 0.5 else (1.0 - t if t > 1e-16 else t)
    return d


def d_sub_target(alo, ahi, blo, bhi):
    def d(r):
        return log_uniform(r, alo, ahi), log_uniform(r, blo, bhi), subnormal(r)
    return d


def d_sub_root(r):
    a = log_uniform(r, 0.03, 1.1)
    b = log_uniform(r, 0.1, 1e3)
    x = subnormal(r)
    lnp = ln_I_left(a, b, x)
    if lnp < -744 or lnp > -1e-3:
        return None
    return a, b, math.exp(lnp)


def d_moderate(r):
    a, b = log_uniform(r, 0.5, 1e3), log_uniform(r, 0.5, 1e3)
    t = log_uniform(r, 1e-300, 0.5) if r.random() < 0.5 else r.uniform(0.001, 0.999)
    return a, b, t


def g_tie_bits(rng, n):
    pts = []
    for a in (3.0, 1e3, 1e8):
        pts += [(a, a, p) for p in neighbours(0.5, 64)]
    return _take(rng, pts, n)


def g_tie_band(rng, n):
    def d(r):
        b = log_uniform(r, 0.5, 1e8)
        a = b * (1 + r.uniform(-1e-6, 1e-6))
        P = approx_I(a, b, 0.5)
        return a, b, P * (1 + r.uniform(-1e-9, 1e-9))
    return _fill_p(rng, n, d)


def g_deep_small_bits(fn):
    def g(rng, n):
        r = random.Random(rng.getrandbits(32))
        pts = []
        while len(pts) < 3 * 129:
            q = d_deep_small(r)
            if q is None:
                continue
            t = _target(r, *q, fn)
            if t is None or not 1e-300 < t < 0.5:
                continue
            pts += [(q[0], q[1], p) for p in neighbours(t, 64)]
        return _take(rng, pts, n)
    return g


# ---- assembly ----------------------------------------------------------------------
# Row budget: the inverse oracle costs ~0.1-0.6 s per row (Newton through the
# bounded forward); BUDGET scales the population strata so the full run stays
# near 20 min on 6 workers. Bit-stepped strata are not scaled.
BUDGET = 0.68


def _xstrat(name, draw, n, seam=None, note=""):
    n = max(1, int(n * BUDGET))
    out = []
    for fn in FNS:
        gen = (lambda fn_: (lambda rng, k: _fill(fn_, rng, k, draw)))(fn)
        out.append(Stratum(fn, name, gen, n, _oracle(fn), _bound(fn),
                           oracle2=_oracle2(fn), seam=seam, note=note))
    return out


def _pstrat(name, draw_or_gen, n, seam=None, note="", per_fn=False, raw=False, ref_cap=None):
    if not name.endswith(":bits"):
        n = max(1, int(n * BUDGET))
    out = []
    for fn in FNS:
        if raw:
            gen = draw_or_gen(fn) if per_fn else draw_or_gen
        else:
            dr = draw_or_gen(fn) if per_fn else draw_or_gen
            gen = (lambda dr_: (lambda rng, k: _fill_p(rng, k, dr_)))(dr)
        out.append(Stratum(fn, name, gen, n, _oracle(fn), _bound(fn),
                           oracle2=_oracle2(fn), seam=seam, note=note, ref_cap=ref_cap))
    return out


def strata():
    S = []
    S += _xstrat("seam:fwd_xi1:band", d_xi1, 600, "fwd_xi1")
    S += _xstrat("seam:fwd_b1:band", d_b1, 500, "fwd_b1")
    S += _xstrat("seam:fwd_r2_thr:band", d_thr, 800, "fwd_r2_thr")
    S += _xstrat("seam:fwd_z2split:band", d_cpsi(22.0, 28.0), 500, "fwd_z2split")
    S += _xstrat("dense:fwd_r3_old_split", d_cpsi(30.0, 42.0), 500, None,
                 "around the pre-#70 split cpsi = 36, now tail lanes")
    S += _xstrat("seam:fwd_z2split:core_below", d_cpsi(20.0, 36.0), 800, "fwd_z2split",
                 "forward R3 core just below cpsi=36 (#70 twin)")
    S += _xstrat("seam:fwd_r3_tail", d_cpsi(25.0, 700.0), 500, "fwd_z2split")
    S += _xstrat("seam:fwd_gammalim:band", d_gl, 500, "fwd_gammalim")
    S += _xstrat("seam:fwd_gl_ridge:band", d_glridge, 500, "fwd_gl_ridge")
    S += _xstrat("seam:deep_small:band", d_deep_small, 800, "deep_small")
    S += _pstrat("seam:deep_small:bits", g_deep_small_bits, 3 * 129, "deep_small",
                 per_fn=True, raw=True)
    S += _xstrat("seam:s3_min:band", d_s3, 500, "s3_min")
    S += _xstrat("seam:s1_nu_min:band", d_s1, 600, "s1_nu_min")
    S += _pstrat("seam:tie_band:bits", g_tie_bits, 3 * 129, "tie_band", raw=True)
    S += _pstrat("seam:tie_band:band", g_tie_band, 300, "tie_band", raw=True)
    S += _pstrat("seam:near_one_log:band", d_near_one, 500, "near_one_log", per_fn=True)
    S += _pstrat("seam:pa_scale:band", d_pa_scale, 300, "pa_scale", per_fn=True)
    S += _pstrat("seam:huge_nu:onset", d_huge_nu(3e30, 3e31, False), 150, "huge_nu",
                 "oracle ~5 s/row (bisection at ~45+log10(nu)/2 digits)")
    S += _pstrat("seam:huge_nu:onset_tail", d_huge_nu(3e30, 3e31, True), 150, "huge_nu")
    # subnormal targets / roots
    S += _pstrat("subnormal:target_large_shapes", d_sub_target(20, 1e6, 20, 1e6), 1000, None,
                 "subnormal p (q) with a, b in [20, 1e6]")
    S += _pstrat("subnormal:target_wide", d_sub_target(0.05, 1e6, 0.05, 1e6), 500)
    S += _pstrat("subnormal:target_huge", d_sub_target(1e6, 1e30, 0.5, 1e30), 200)
    S += _pstrat("subnormal:root", d_sub_root, 800, None, "subnormal x results")
    # dense
    # #79: 400 reference rows, not 120 -- the file's overall not-CR rate
    # (3.93%) hid this population's (~15%).
    S += _pstrat("dense:huge_nu_tail", d_huge_nu(1e20, 1e40, True), 600, "huge_nu",
                 "huge-nu tail-p rows", ref_cap=400)
    S += _pstrat("dense:huge_nu_mid", d_huge_nu(1e20, 1e40, False), 150, "huge_nu")
    S += _pstrat("dense:moderate", d_moderate, 1200)
    return S + _t2_strata()


# ---- T2 completion strata (#82): bit neighbourhoods the audit found empty ----
def _exact_target(a, b, x, fn):
    """The rounded p or q whose root is x to within an ulp or so: the target
    from the exact forward (the float approx_I cannot hold a root to 64 ulp)."""
    from ._beta_oracles import betainc_pq
    with mp.workdps(40):
        P, Q = betainc_pq(a, b, x, 30)
        t = float(P if fn == "beta_p_inv" else Q)
    return t if 0 < t < 1 else None


def _root_bits(points):
    """points: [(a, b, x0)] -> per-fn generator of (a, b, target) rows whose
    roots step through x0's bit neighbourhood."""
    def per_fn(fn):
        def g(rng, n):
            rows = []
            for a, b, x0 in points:
                for x in neighbours(x0, 12):
                    t = _exact_target(a, b, x, fn)
                    if t is not None:
                        rows.append((a, b, t))
            return _take(rng, rows, n)
        return g
    return per_fn


def _shape_bits(shapes, targets):
    """[(a, b)] stepped by the caller -> rows crossed with fixed targets."""
    return lambda rng, n: _take(rng, [(a, b, s) for a, b in shapes for s in targets], n)


def _t2_strata():
    S = []
    S += _pstrat("seam:fwd_xi1:bits", _root_bits([(2.0, 5.0, XI1), (0.7, 3.0, XI1), (5.0, 12.0, XI1)]),
                 3 * 25, "fwd_xi1", per_fn=True, raw=True)
    S += _pstrat("seam:fwd_b1:bits", _root_bits([(1.5, 16.0, 0.5), (3.0, 20.0, 0.4), (1.2, 40.0, 0.2)]),
                 3 * 25, "fwd_b1", per_fn=True, raw=True)
    z2 = [(a, b, x_for_cpsi(a, b, 25.0, up)) for a, b, up in ((40.0, 60.0, False), (200.0, 90.0, True),
                                                             (1e4, 3e4, False))]
    S += _pstrat("seam:fwd_z2split:bits", _root_bits(z2), 3 * 25, "fwd_z2split", per_fn=True, raw=True)
    tg = (1e-5, 0.3, 0.7)
    S += _pstrat("seam:s3_min:bits",
                 _shape_bits([(m, w) for w in (0.5, 30.0) for m in neighbours(2.0 ** 20, 16)]
                             + [(w, m) for w in (1e4,) for m in neighbours(2.0 ** 20, 16)], tg),
                 3 * 99, "s3_min", raw=True)
    S += _pstrat("seam:s1_nu_min:bits",
                 _shape_bits([(a, 4.0) for a in neighbours(4.0, 16)]
                             + [(a, 6.0) for a in neighbours(3.0, 16)], tg),
                 3 * 66, "s1_nu_min", raw=True)
    S += _pstrat("seam:huge_nu:bits",
                 _shape_bits([(a, 2e31) for a in neighbours(2e31, 16)], (0.3, 0.5, 1e-5)),
                 3 * 33, "huge_nu", raw=True)
    a20 = 20.0 / (1.0 - 20.0 / 2.0 ** 60)
    S += _pstrat("seam:fwd_gl_ridge:bits",
                 _shape_bits([(a, 2.0 ** 60) for a in neighbours(a20, 16)], (1e-5, 0.3, 0.7)),
                 3 * 33, "fwd_gl_ridge", raw=True)
    # ln(1 - sigma) switches to Log1pDd at sigma = 2^-7 (#88).
    t7 = 2.0 ** -7
    S += _pstrat("seam:lgamma1p_table_max:bits",
                 _shape_bits([(a, b) for a in neighbours(15.0, 16) for b in (20.0, 80.0)],
                             (1e-5, 0.3, 0.7)),
                 3 * 66, "lgamma1p_table_max", raw=True)
    S += _pstrat("seam:lgamma1p_table_max:band",
                 lambda r: (full_mantissa(r, 14.2, 15.8), log_uniform(r, 15.0, 200.0),
                            log_uniform(r, 1e-6, 0.5)),
                 300, "lgamma1p_table_max")
    S += _pstrat("seam:lnsigc_log1p:bits",
                 _shape_bits([(2.0, 5.0), (0.3, 40.0), (1e3, 1e3)], neighbours(t7, 64)),
                 3 * 129, "lnsigc_log1p", raw=True)
    S += _pstrat("seam:lnsigc_log1p:band",
                 lambda r: (log_uniform(r, 0.05, 1e6), log_uniform(r, 0.05, 1e6),
                            full_mantissa(r, t7 * 0.95, t7 * 1.05)),
                 400, "lnsigc_log1p")
    return S


# ---- seam-audit fields (#82 T2; tools/hunt/seam_audit.py docstring) ----------
def _nu_r(r):
    return 1.0 / (1.0 / r["a"] + 1.0 / r["b"])


def _p_of(r):
    """The lower-tail target of a reference row (betainv_q files store q)."""
    return 1.0 - r["s"] if r["file"].endswith("_q") else r["s"]


def _tie(r):
    """|target - side(1/2)|; float first, the exact forward near the seam
    (the float approx_I cannot resolve 64 ulp of 9e-13)."""
    a, b = r["a"], r["b"]
    if a + b > 1e7:
        return None
    side_q = r["file"].endswith("_q")
    c = abs(r["s"] - (1.0 - approx_I(a, b, 0.5) if side_q else approx_I(a, b, 0.5)))
    if abs(c - TIE_BAND) > 1e-11:
        return c
    from ._beta_oracles import betainc_pq
    with mp.workdps(40):
        P, Q = betainc_pq(a, b, 0.5, 30)
        return float(abs(mp.mpf(r["s"]) - (Q if side_q else P)))


_AUDIT = {
    "fwd_r2_thr": {"audit_f": lambda r: r["root x"] / ((r["a"] + 1.0) / (r["a"] + r["b"] + 2.0)),
                   "audit_value": 1.0},
    "fwd_z2split": {"audit_f": lambda r: cpsi_f(r["a"], r["b"], r["root x"])},
    "fwd_gl_ridge": {"audit_f": lambda r: min(r["a"], r["b"]) if max(r["a"], r["b"]) >= GAMMA_LIM else None},
    "deep_small": {"route": "betainv/deep-small"},
    "s1_nu_min": {"audit_f": _nu_r},
    "huge_nu": {"audit_f": _nu_r},
    # The coordinate moves in steps of ulp(p) ~ 2^-53 near 1/2, far coarser
    # than 64 ulp of 2^-40: the bits window is 64 of those steps.
    "tie_band": {"audit_f": _tie, "window": 64 * 2.0 ** -53},
    "near_one_log": {"route": "betainv/fwd/near-one"},
    # sigma = 1 - s on a row near 1 moves in steps of 2^-53: window as tie_band.
    "lnsigc_log1p": {"audit_f": lambda r: min(r["s"], 1.0 - r["s"]), "window": 64 * 2.0 ** -53},
    # Both select inside the deep-small closed form (itself route-covered):
    # the census marks the scaled branch; the #88 band sweep checked it
    # against an independent quadrature.
    "lgamma1p_table_max": {"audit_f": lambda r: min(r["a"], r["b"])},
    "ds_scaled_a": {"route": "betainv/deep-small/scaled"},
    "ds_scaled_b": {"route": "betainv/deep-small/scaled"},
    "pa_scale": {"internal": "max(|E|, |ln B|) is a prescale test inside the forward's "
                             "PA assembly; the census marks the PA form, not this threshold"},
}
for _s in SEAMS:
    _s.update(_AUDIT.get(_s["name"], {}))
