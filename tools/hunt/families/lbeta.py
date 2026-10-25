"""lbeta hunt strata (#82).

Bound (docs/ACCURACY.md, lbeta): 1 ULP wherever |ln B| >= 1 (relative
metric), and absolute 1 * 2^-53 in the band |ln B| < 1 around ln B's zero
curve through (1, 1). It was documented as correctly rounded until #76
(order row 12j) found ~2% not-CR on [0.25, 4]^2; tie_exempt = 0.

Oracle: ln B = lgamma(m) + lgamma(M) - lgamma(m + M) at cancellation-aware
precision (_beta_oracles.ln_beta: absolute accuracy 10^-dps, digits for the
a+b cancellation and for |lgamma(M)| ~ M ln M; m psi(M) Taylor form for
m < 10^-(dps/2+8); exact (m+M)/(mM) split for M <= 1). oracle2:
log(beta(a, b)) at the same raised precision (mpmath's gamma-ratio path).
True values below -DBL_MAX round to -inf (the oracle's finite value
converts to -inf and is scored as exact agreement).
"""
from __future__ import annotations

import math

import mpmath as mp
from mpmath import mpf

from core import Stratum, full_mantissa, log_uniform, neighbours, subnormal
from ._beta_oracles import ln_beta

BIG_MIN = 2.0 ** 990          # kLbetaBigMin (src/beta-inl.h)
SCALE_DOWN = 2.0 ** -64       # kLbetaScaleDown (exact prescale in the big band)
SCALE_UP = 2.0 ** 64          # kLbetaScaleUp
Z0 = 10.0                     # kBetaZ0 (LgammaDiffDd up-walk target)
DIV_FLOOR = 2.0 ** -900       # kBetaDivFloor (DdDivDd, #59)
SCALE_ABOVE = 2.0 ** 900      # kBetaScaleAbove (LgammaDiffDd z prescale)
BINET_Z_CLAMP = 2.0 ** 500    # kBetaBinetZClamp
LG_ZONES = (0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.0)  # lgamma_data.h zones, X0

SEAMS = [
    {"fn": "lbeta", "name": "big_min", "coord": "min(a,b)", "value": BIG_MIN, "source": "kLbetaBigMin (grouped Stirling band, 2^-64/2^64 prescale)"},
    {"fn": "lbeta", "name": "neg_inf", "coord": "a=b", "value": "~1.2968e308", "source": "ln B = -DBL_MAX crossing (big band scale-back overflow)"},
    {"fn": "lbeta", "name": "lg_zones", "coord": "min(a,b)", "value": list(LG_ZONES), "source": "kLgammaZoneLo/Mid/Hi + mid steps, kLgammaX0 (LgammaPosDd)"},
    {"fn": "lbeta", "name": "walk", "coord": "max(a,b)", "value": [float(k) for k in range(1, 11)], "source": "kBetaZ0 up-walk (LgammaDiffDd)"},
    {"fn": "lbeta", "name": "div_floor", "coord": "max(a,b)", "value": DIV_FLOOR, "source": "kBetaDivFloor (DdDivDd)"},
    {"fn": "lbeta", "name": "scale_above", "coord": "max(a,b)", "value": SCALE_ABOVE, "source": "kBetaScaleAbove (LgammaDiffDd)"},
    {"fn": "lbeta", "name": "binet_z_clamp", "coord": "max(a,b)", "value": BINET_Z_CLAMP, "source": "kBetaBinetZClamp"},
    {"fn": "lbeta", "name": "abs_band", "coord": "|ln B|", "value": 1.0, "source": "ACCURACY.md absolute band edge"},
]

ABS_TOL = mpf(2) ** -53  # #76: was 0.5 * 2^-53


def _oracle(args, dps):
    return ln_beta(args[0], args[1], dps)


def _oracle2(args, dps):
    a, b = mpf(args[0]), mpf(args[1])
    m, M = min(a, b), max(a, b)
    extra = max(0, int(mp.log10(M / m))) + 20
    if M > 10:
        extra += int(mp.log10(M * mp.log(M))) + 1
    with mp.workdps(dps + extra):
        return mp.log(mp.beta(a, b))


def _bound(args, v):
    return 1  # #76


def _abs_tol(args, v):
    return ABS_TOL if abs(v) < 1 else None


def _st(name, gen, n, seam=None, note=""):
    return Stratum("lbeta", name, gen, n, _oracle, _bound, oracle2=_oracle2,
                   abs_tol=_abs_tol, seam=seam, note=note)


def _ok(a, b):
    return math.isfinite(a) and math.isfinite(b) and a > 0 and b > 0


def _take(rng, pts, n):
    pts = [p for p in pts if _ok(*p)]
    if n < len(pts):
        pts = rng.sample(pts, n)
    return pts


def _fill(rng, n, draw):
    out = []
    for _ in range(n * 20):
        if len(out) >= n:
            break
        p = draw(rng)
        if p is not None and _ok(*p):
            out.append(p if rng.random() < 0.5 else (p[1], p[0]))
    return out


def _neg_inf_a0():
    """a with ln B(a, a) = -DBL_MAX (bisection on the Stirling form)."""
    with mp.workdps(40):
        target = -mpf(2) ** 1024 * (1 - mpf(2) ** -53)
        lo, hi = mpf("1.2e308"), mpf("1.4e308")
        for _ in range(200):
            mid = (lo + hi) / 2
            v = -2 * mid * mp.log(2) + mp.log(mp.pi / mid) / 2
            if v > target:
                lo = mid
            else:
                hi = mid
        return float(lo)


def g_neg_inf_bits(rng, n):
    a0 = _neg_inf_a0()
    pts = [(a, a) for a in neighbours(a0, 64)]
    pts += [(a, math.nextafter(a, math.inf)) for a in neighbours(a0, 32)]
    return _take(rng, pts, n)


def g_lg_zones_bits(rng, n):
    pts = []
    for z in LG_ZONES:
        for o in (z * 1.3, 50.0, 1e10):
            pts += [(m, o) for m in neighbours(z, 24)]
    return _take(rng, pts, n)


def g_lg_zones_band(rng, n):
    def d(r):
        z = r.choice(LG_ZONES)
        return full_mantissa(r, z * 0.95, z * 1.05), log_uniform(r, 0.3, 1e12)
    return _fill(rng, n, d)


def g_walk_bits(rng, n):
    pts = []
    for M in (10.0, 9.0, 5.0, 2.0, 1.0):
        for m in (0.3, 0.9 * M, 1e-5):
            pts += [(m, x) for x in neighbours(M, 24)]
    return _take(rng, pts, n)


def g_walk_band(rng, n):
    def d(r):
        M = r.randint(1, 10) * r.uniform(0.95, 1.05)
        return log_uniform(r, 1e-6, M), M
    return _fill(rng, n, d)


def g_zero_curve(rng, n):
    """Points with |ln B| < 1 (abs-judged band): around B(a, b) = 1."""
    def d(r):
        a = log_uniform(r, 0.2, 6)
        # B(a, b) = 1 -> solve b by bisection on lgamma (float)
        lo, hi = 1e-3, 50.0
        f = lambda b: math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)
        if (f(lo) > 0) == (f(hi) > 0):
            return None
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if (f(mid) > 0) == (f(lo) > 0):
                lo = mid
            else:
                hi = mid
        b = 0.5 * (lo + hi) * (1 + r.uniform(-0.3, 0.3))
        return a, b
    return _fill(rng, n, d)


def g_witness(rng, n):
    a, b = float.fromhex("0x1.5c788bfb641a0p+1"), float.fromhex("0x1.5f0d82e6d1d90p+0")
    return [(a, b)] + [(x, b) for x in neighbours(a, 8) if x != a]


def strata():
    S = []
    S.append(_st("seam:big_min:bits", lambda rng, n: _take(rng, [
        (m, o) for m in neighbours(BIG_MIN, 64)
        for o in (m, 2.0 * BIG_MIN, 1.7e308, m * 1.000001)], n),
        4 * 129, "big_min", "min crosses 2^990: main band -> grouped Stirling"))
    S.append(_st("seam:big_min:band", lambda rng, n: _fill(rng, n, lambda r: (
        full_mantissa(r, BIG_MIN * 0.9, BIG_MIN * 1.1),
        full_mantissa(r, BIG_MIN * 0.9, 1.7e308) if r.random() < 0.5 else log_uniform(r, BIG_MIN * 1.1, 1.7e308))),
        2000, "big_min"))
    S.append(_st("seam:neg_inf:bits", g_neg_inf_bits, 129 + 65, "neg_inf"))
    S.append(_st("seam:lg_zones:bits", g_lg_zones_bits, len(LG_ZONES) * 3 * 49, "lg_zones"))
    S.append(_st("seam:lg_zones:band", g_lg_zones_band, 4000, "lg_zones"))
    S.append(_st("seam:walk:bits", g_walk_bits, 5 * 3 * 49, "walk"))
    S.append(_st("seam:walk:band", g_walk_band, 3000, "walk"))
    S.append(_st("seam:div_floor:bits", lambda rng, n: _take(rng, [
        (M, m) for M in neighbours(DIV_FLOOR, 64) for m in (5e-324, 2.0 ** -1000, M * 0.5)], n),
        3 * 129, "div_floor"))
    S.append(_st("seam:div_floor:band", lambda rng, n: _fill(rng, n, lambda r: (
        full_mantissa(r, DIV_FLOOR * 0.9, DIV_FLOOR * 1.1), log_uniform(r, 5e-324, DIV_FLOOR))),
        1500, "div_floor"))
    S.append(_st("seam:scale_above:bits", lambda rng, n: _take(rng, [
        (M, m) for M in neighbours(SCALE_ABOVE, 64) for m in (0.5, 1e10, 2.0 ** 880)], n),
        3 * 129, "scale_above"))
    S.append(_st("seam:scale_above:band", lambda rng, n: _fill(rng, n, lambda r: (
        full_mantissa(r, SCALE_ABOVE * 0.9, SCALE_ABOVE * 1.1), log_uniform(r, 1e-300, 2.0 ** 899))),
        1500, "scale_above"))
    S.append(_st("seam:binet_z_clamp:bits", lambda rng, n: _take(rng, [
        (M, m) for M in neighbours(BINET_Z_CLAMP, 64) for m in (0.5, 30.0, 1e100)], n),
        3 * 129, "binet_z_clamp"))
    S.append(_st("seam:binet_z_clamp:band", lambda rng, n: _fill(rng, n, lambda r: (
        full_mantissa(r, BINET_Z_CLAMP * 0.9, BINET_Z_CLAMP * 1.1), log_uniform(r, 1e-300, 1e140))),
        1000, "binet_z_clamp"))
    S.append(_st("seam:abs_band:zero_curve", g_zero_curve, 4000, "abs_band",
                 "|ln B| < 1: judged absolutely (0.5 * 2^-53)"))
    # subnormal
    S.append(_st("subnormal:one", lambda rng, n: _fill(rng, n, lambda r: (
        subnormal(r), log_uniform(r, 1e-300, 1e300))), 3000))
    S.append(_st("subnormal:both", lambda rng, n: _fill(rng, n, lambda r: (
        subnormal(r), subnormal(r))), 1500))
    # crossed extremes
    S.append(_st("crossed:tiny_huge", lambda rng, n: _fill(rng, n, lambda r: (
        log_uniform(r, 1e-320, 1e-200), log_uniform(r, 1e200, 1.7e308))), 1500))
    S.append(_st("crossed:ratio_sweep", lambda rng, n: _fill(rng, n, lambda r: (
        log_uniform(r, 1e-300, 1e300), log_uniform(r, 1e-300, 1e300))), 4000))
    # dense
    S.append(_st("dense:quarter_to_four", lambda rng, n: _fill(rng, n, lambda r: (
        full_mantissa(r, 0.25, 4.0), full_mantissa(r, 0.25, 4.0))), 20000, None,
        "[0.25, 4]^2 (known ~2% not CR, #76)"))
    S.append(_st("dense:moderate", lambda rng, n: _fill(rng, n, lambda r: (
        log_uniform(r, 1e-3, 1e6), log_uniform(r, 1e-3, 1e6))), 6000))
    S.append(_st("witness", g_witness, 17, None,
                 "lbeta(0x1.5c788bfb641a0p+1, 0x1.5f0d82e6d1d90p+0) 1 ULP vs 0"))
    # T2 completion (#82): bucket gaps the audit found.
    for k in (3.0, 4.0, 6.0, 7.0, 8.0):  # g_walk_bits steps 10, 9, 5, 2, 1 only
        S.append(_st(f"seam:walk_{int(k)}:bits", lambda rng, n, k=k: _take(rng, [
            (m, x) for m in (0.3, 0.9 * k, 1e-5) for x in neighbours(k, 16)], n), 3 * 33, "walk"))
    a0 = _neg_inf_a0()
    S.append(_st("seam:neg_inf:band", lambda rng, n: [(lambda a: (a, a))(full_mantissa(rng, a0 * 0.95, a0 * 1.05))
                                                    for _ in range(n)], 300, "neg_inf"))

    def _b_for_lnB(a, target):  # ln B(a, b) decreases in b
        lo, hi = 1e-300, 1e300
        for _ in range(400):
            mid = math.sqrt(lo * hi)
            v = math.lgamma(a) + math.lgamma(mid) - math.lgamma(a + mid)
            if v > target:
                lo = mid
            else:
                hi = mid
        return math.sqrt(lo * hi)
    ab = [(a, _b_for_lnB(a, t)) for a in (0.3, 2.0, 7.0) for t in (1.0, -1.0)]
    S.append(_st("seam:abs_band:bits", lambda rng, n: _take(rng, [
        (a, b) for a, b0 in ab for b in neighbours(b0, 12)], n), 6 * 25, "abs_band",
        "|ln B| = 1, the absolute-band edge"))
    return S


# ---- seam-audit fields (#82 T2; tools/hunt/seam_audit.py docstring) ----------
_AUDIT = {
    "neg_inf": {"audit_f": lambda r: r["a"] if r["a"] == r["b"] else None,
                "audit_value": _neg_inf_a0()},
    "abs_band": {"audit_f": lambda r: abs(r["lnB"])},
}
for _s in SEAMS:
    _s.update(_AUDIT.get(_s["name"], {}))
