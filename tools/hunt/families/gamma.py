"""gamma_p / gamma_q hunt strata (#82).

Bounds (docs/ACCURACY.md, gamma_p and gamma_q): 1 ULP everywhere except the
DIRECT side of R3 (Temme ridge: a >= 20, 1/2 < x/a < 2), published 2 ULP.
Direct side in R3 is P for x < a and Q for x >= a (GammaTemme's sgn rule:
CopySign(1, x - a), so x == a is Q-direct). R1/R2/R4 directs and every
complement are 1 ULP (R2/R4 complements by the row-12 erratum). No absolute
band. Oracle: _gamma_oracles.pq (small side direct, capped series/CF/R4
identity); oracle2 = Gauss CF / backward Legendre CF / mpmath upper.
Not covered: a > 1e7 except saturated rows (a*phi > 790) and the exact
diagonal x == a for a >= 2^100.
"""
import math

import mpmath as mp

from core import Stratum, full_mantissa, log_uniform, neighbours, subnormal

from . import _gamma_oracles as G

FNS = ("gamma_p", "gamma_q")
A_T = 20.0
# Sampled (band) strata are scaled by this so the module totals ~80k rows at
# --scale 1 (bit-stepped strata are fixed-size and subsample below it).
BAND_SCALE = 0.5

SEAMS = [
    # region map (gamma-inl.h GammaVec, gamma_data.h)
    {"fn": "gamma_p/q", "name": "aT", "coord": "a", "value": 20.0, "source": "kGammaAT (gamma_data.h): R1/R2 vs R3 split"},
    {"fn": "gamma_p/q", "name": "x_eq_a_plus_1", "coord": "x-(a+1)", "value": 0.0, "source": "GammaVec i_xle: R1/R2 split for a < 20"},
    {"fn": "gamma_p/q", "name": "lam_half", "coord": "x/a", "value": 0.5, "source": "GammaVec i_lo (2x <= a): R1/R3; = kGammaEtaMid-kGammaEtaHalf"},
    {"fn": "gamma_p/q", "name": "lam_two", "coord": "x/a", "value": 2.0, "source": "GammaVec i_hi (x >= 2a): R3/R2; = kGammaEtaMid+kGammaEtaHalf"},
    {"fn": "gamma_p/q", "name": "lam_one", "coord": "x/a", "value": 1.0, "source": "GammaTemme sgn = sign(x-a): R3 direct-side switch"},
    {"fn": "gamma_p/q", "name": "r4_amax", "coord": "a", "value": 1.5, "source": "kGammaR4AMax = kLgammaZoneHi-1 (R4 box)"},
    {"fn": "gamma_p/q", "name": "r4_xmax", "coord": "x", "value": 4.0, "source": "kGammaR4XMax (R4 box)"},
    {"fn": "gamma_p/q", "name": "r4_zone_half", "coord": "a", "value": 0.5, "source": "kLgammaZoneLo: GammaSmallQ lgamma(1+a) centre switch; LgammaPosDd zone"},
    {"fn": "gamma_p/q", "name": "r4_log1p_cut", "coord": "a", "value": 2.0 ** -7, "source": "GammaSmallQ r_cut 0x1p-7 (Log1pDd domain)"},
    {"fn": "gamma_p/q", "name": "lgamma_zone_hi", "coord": "a", "value": 2.5, "source": "kLgammaZoneHi (LgammaPosDd, prefactor E for R1/R2)"},
    {"fn": "gamma_p/q", "name": "lgamma_x0", "coord": "a", "value": 8.0, "source": "kLgammaX0 (LgammaPosDd Stirling switch)"},
    # Temme (R3)
    {"fn": "gamma_p/q", "name": "temme_split", "coord": "a*phi", "value": 36.0, "source": "kGammaTemmeZ2Split (z = 6 = kGammaTailZLo, erfc core/tail)"},
    {"fn": "gamma_p/q", "name": "expfloor_r3", "coord": "a*phi", "value": 800.0, "source": "-kGammaExpFloor: R3 saturation and live_tail"},
    {"fn": "gamma_p/q", "name": "tail_zhi", "coord": "a*phi", "value": 812.25, "source": "kGammaTailZHi^2 (28.5^2): inside saturation, result 0 -- no stratum"},
    {"fn": "gamma_p/q", "name": "twopi_a_clamp", "coord": "a", "value": 2.0 ** 900, "source": "kGammaTwoPiAClamp = kGammaScaleAbove (prescale)"},
    # R1/R2 prefactor
    {"fn": "gamma_p/q", "name": "expfloor_r1", "coord": "E-log a", "value": -800.0, "source": "kGammaExpFloor on R1's exponent"},
    {"fn": "gamma_p/q", "name": "expfloor_r2", "coord": "E", "value": -800.0, "source": "kGammaExpFloor on R2's exponent"},
    # outputs
    {"fn": "gamma_p/q", "name": "out_min_normal", "coord": "result", "value": 2.0 ** -1022, "source": "derived: subnormal output band"},
    {"fn": "gamma_p/q", "name": "out_zero", "coord": "result", "value": 2.0 ** -1075, "source": "derived: zero output"},
    # not mappable to an input coordinate
    {"fn": "gamma_p/q", "name": "freeze_eps", "coord": "term/sum", "value": 2.0 ** -60, "source": "kGammaFreezeEps: per-term freeze, no input seam (unmapped)"},
    {"fn": "gamma_p/q", "name": "loop_shapes", "coord": "-", "value": 0, "source": "kGammaSeriesN/CfN/R4N/CfDdSteps/TemmeK: fixed depths, worst at the a->20-, x=a+1 corner (covered by aT x a+1 rows)"},
]


# ---- oracles and bound ------------------------------------------------------
def _oracle(fn):
    return lambda args, dps: G.gamma_value(fn, args[0], args[1], dps)


def _oracle2(fn):
    return lambda args, dps: G.gamma_value(fn, args[0], args[1], dps, two=True)


def region(a, x):
    """GammaVec's region map, replayed in double arithmetic."""
    small = a < A_T
    xle = a + 1.0 >= x
    lo = a >= x + x
    hi = x >= a + a
    if small:
        return "R1" if xle else "R2"
    if lo:
        return "R1"
    if hi:
        return "R2"
    return "R3"


def _bound(fn):
    def b(args, v):
        a, x = args
        if region(a, x) == "R3":
            p_direct = x < a
            if (fn == "gamma_p") == p_direct:
                return 2
        return 1
    return b


# ---- float helpers for the generators ---------------------------------------
def _phi(u):
    """phi(1+u) = u - log1p(u)."""
    return u - math.log1p(u)


def lam_for_aphi(a, L, sgn):
    """lambda with a*phi(lambda) = L on the side sgn (+1: lambda > 1);
    float accuracy (band placement only)."""
    target = L / a
    if sgn > 0:
        lo, hi = 0.0, 1.0
        while _phi(hi) < target:
            hi *= 2
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if _phi(mid) < target:
                lo = mid
            else:
                hi = mid
        return 1.0 + 0.5 * (lo + hi)
    # lambda < 1: bisect y = ln(lambda) in [-745, 0], phi = expm1(y) - y
    lo, hi = -745.0, 0.0
    if math.expm1(lo) - lo < target:
        return None
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if math.expm1(mid) - mid > target:
            lo = mid
        else:
            hi = mid
    return math.exp(0.5 * (lo + hi))


def x_for_aphi(a, L, sgn):
    lam = lam_for_aphi(a, L, sgn)
    return None if lam is None else lam * a


def x_for_aphi_exact(a, L, sgn):
    """Seam location to < 1 ulp (bit-stepped strata): mpmath bisection."""
    with mp.workdps(40):
        a_m, L = mp.mpf(a), mp.mpf(L)
        lo, hi = (mp.mpf(a), mp.mpf(a) * 64) if sgn > 0 else (mp.mpf(0), mp.mpf(a))
        for _ in range(300):
            mid = (lo + hi) / 2
            v = a_m * (mid / a_m - 1 - mp.log(mid / a_m))
            if (v < L) == (sgn > 0):
                lo = mid
            else:
                hi = mid
        return float(lo)


def x_for_E_exact(a, target, r1):
    with mp.workdps(40):
        a_m = mp.mpf(a)
        lg = mp.loggamma(a_m) + (mp.log(a_m) if r1 else 0)
        f = lambda x: a_m * mp.log(x) - x - lg
        if r1:
            lo, hi = mp.mpf(2) ** -1074, a_m
        else:
            lo, hi = a_m, a_m * 4 + 4000
        if (f(lo) - target) * (f(hi) - target) > 0:
            return None
        inc = f(hi) > f(lo)
        for _ in range(400):
            mid = (lo + hi) / 2
            if (f(mid) < target) == inc:
                lo = mid
            else:
                hi = mid
        return float(lo)


def _E(a, x, r1):
    e = a * math.log(x) - x - math.lgamma(a)
    return e - math.log(a) if r1 else e


def x_for_E(a, target, r1):
    """x with (R1) E - log a = target on x <= a, or (R2) E = target on x >= a."""
    if r1:
        lo, hi = math.log(5e-324), math.log(a)
        if _E(a, math.exp(lo), True) > target or _E(a, a, True) < target:
            return None
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if _E(a, math.exp(mid), True) < target:
                lo = mid
            else:
                hi = mid
        return math.exp(0.5 * (lo + hi))
    lo, hi = max(a, 1e-300), max(2 * a, 10.0)
    while _E(a, hi, False) > target:
        hi *= 2
        if hi > 1e300:
            return None
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if _E(a, mid, False) > target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _sub(r):
    """Subnormal, half uniform in mantissa bits, half uniform in exponent."""
    return subnormal(r) if r.random() < 0.5 else subnormal(r, 52)


def _take(rng, rows, n):
    rows = [r for r in rows if r is not None and all(v > 0 and math.isfinite(v) for v in r)]
    return rows if n >= len(rows) else rng.sample(rows, n)


def _many(rng, n, f, tries=20):
    out = []
    for _ in range(n * tries):
        if len(out) >= n:
            break
        r = f(rng)
        if r is not None and all(v > 0 and math.isfinite(v) for v in r):
            out.append(r)
    return out


# ---- generators --------------------------------------------------------------
def _bits_a(t, xfuncs):
    def gen(rng, n):
        rows = [(a, xf(a)) for a in neighbours(t, 64) for xf in xfuncs]
        return _take(rng, rows, n)
    return gen, 129 * len(xfuncs)


def _band_a(t, delta, side, xsampler):
    def gen(rng, n):
        lo, hi = (t * (1 - delta), t) if side < 0 else (t, t * (1 + delta))
        return _many(rng, n, lambda r: (lambda a: (a, xsampler(r, a)))(full_mantissa(r, lo, hi)))
    return gen


def _lam_x(lams):
    return [lambda a, l=l: l * a for l in lams] + [lambda a: a + 1.0]


def strata():
    out = []

    def add(name, gen, n, seam=None, note=""):
        if not (name.endswith(":bits") or name.startswith(("witness", "extreme:diag"))):
            n = max(1, int(n * BAND_SCALE))
        for fn in FNS:
            out.append(Stratum(fn, name, gen, n, _oracle(fn), _bound(fn),
                               oracle2=_oracle2(fn), seam=seam, note=note))

    # 1. a seams: bit-stepped a crossed with x, and full-mantissa bands.
    lam_cross = _lam_x([0.3, 0.5, 0.9, 1.0, 1.6, 2.0, 3.0])
    box_cross = [lambda a, v=v: v for v in (1e-3, 0.3, 1.0, 3.0, 4.0, 5.0)] + [lambda a: a + 1.0]
    a_seams = [
        ("aT", A_T, lam_cross, lambda r, a: a * log_uniform(r, 0.25, 4.0)),
        ("r4_amax", 1.5, box_cross, lambda r, a: full_mantissa(r, 1e-6, 5.0)),
        ("r4_zone_half", 0.5, box_cross, lambda r, a: full_mantissa(r, 1e-6, 5.0)),
        ("r4_log1p_cut", 2.0 ** -7, [lambda a, v=v: v for v in (1e-300, 1e-3, 0.3, 1.0, 2.0, 4.0)],
         lambda r, a: log_uniform(r, 1e-300, 4.0) if r.random() < 0.3 else full_mantissa(r, 1e-6, 4.0)),
        ("lgamma_zone_hi", 2.5, _lam_x([0.3, 0.8, 1.2, 2.0, 3.0]), lambda r, a: a * log_uniform(r, 0.1, 6.0)),
        ("lgamma_x0", 8.0, _lam_x([0.3, 0.8, 1.2, 2.0, 3.0]), lambda r, a: a * log_uniform(r, 0.1, 6.0)),
    ]
    for name, t, xf, samp in a_seams:
        g, nb = _bits_a(t, xf)
        add(f"seam:{name}:bits", g, nb, seam=name)
        for side, tag in ((-1, "lo"), (1, "hi")):
            add(f"seam:{name}:band_{tag}", _band_a(t, 0.05, side, samp), 900, seam=name)

    # 2. x seam R4XMax = 4.
    a_set = (1e-300, 1e-20, 1e-5, 0.1, 0.43, 0.5, 1.0, 1.5)
    add("seam:r4_xmax:bits",
        lambda r, n, a_set=a_set: _take(r, [(a, x) for x in neighbours(4.0, 64) for a in a_set], n),
        129 * len(a_set), seam="r4_xmax")
    for lo, hi, tag in ((3.8, 4.0, "lo"), (4.0, 4.2, "hi")):
        add(f"seam:r4_xmax:band_{tag}",
            lambda r, n, lo=lo, hi=hi: _many(r, n, lambda q: (
                log_uniform(q, 1e-8, 1.5) if q.random() < 0.8 else log_uniform(q, 1e-308, 1e-8),
                full_mantissa(q, lo, hi))),
            1000, seam="r4_xmax")

    # 3. x = a + 1 (R1/R2 for a < 20; the series cap / CF depth corner).
    a_set = (0.37, 1.5, 1.7, 3.3, 7.5, 12.1, 19.75, math.nextafter(20.0, 0))
    add("seam:x_eq_a_plus_1:bits",
        lambda r, n, a_set=a_set: _take(r, [(a, x) for a in a_set for x in neighbours(a + 1.0, 64)], n),
        129 * len(a_set), seam="x_eq_a_plus_1")
    for sgn, tag in ((-1, "lo"), (1, "hi")):
        add(f"seam:x_eq_a_plus_1:band_{tag}",
            lambda r, n, sgn=sgn: _many(r, n, lambda q: (lambda a: (a, (a + 1) * (1 + sgn * 0.04 * q.random())))(
                full_mantissa(q, 0.05, 20.0))),
            1200, seam="x_eq_a_plus_1")

    # 4. lambda seams 1/2, 2 (region) and 1 (R3 side switch).
    a_set = (20.0, 20.5, 33.3, 117.0, 1e3, 1e4, 1e5, 1e6)
    for name, lam in (("lam_half", 0.5), ("lam_two", 2.0), ("lam_one", 1.0)):
        add(f"seam:{name}:bits",
            lambda r, n, lam=lam, a_set=a_set: _take(r, [(a, x) for a in a_set for x in neighbours(lam * a, 64)], n),
            129 * len(a_set), seam=name)
        for sgn, tag in ((-1, "lo"), (1, "hi")):
            if lam == 1.0:
                f = (lambda q, sgn=sgn: (lambda a: (a, a * (1 + sgn * 3.0 * q.random() / math.sqrt(a))))(
                    log_uniform(q, 20.0, 1e6)))
            else:
                f = (lambda q, sgn=sgn, lam=lam: (lambda a: (a, a * lam * (1 + sgn * 0.05 * q.random())))(
                    log_uniform(q, 20.0, 1e6)))
            add(f"seam:{name}:band_{tag}", lambda r, n, f=f: _many(r, n, f), 1200, seam=name)

    # 5. Temme core/tail split a*phi = 36 (R3 had 1 reference row in [30, 36)).
    a_set = (150.0, 300.0, 1e3, 1e4, 1e5, 1e6)
    add("seam:temme_split:bits",
        lambda r, n, a_set=a_set: _take(r, [(a, x) for a in a_set for s in (-1, 1)
                               for x in neighbours(x_for_aphi_exact(a, 36.0, s), 64)], n),
        129 * 2 * len(a_set), seam="temme_split")
    for lo, hi, tag in ((30.0, 36.0, "lo"), (36.0, 40.0, "hi")):
        add(f"seam:temme_split:band_{tag}",
            lambda r, n, lo=lo, hi=hi: _many(r, n, lambda q: (lambda a: (a, x_for_aphi(a, q.uniform(lo, hi), q.choice((-1, 1)))))(
                log_uniform(q, 100.0, 1e6))),
            2000, seam="temme_split")
    add("corner:temme_core_lam2",
        lambda r, n: _many(r, n, lambda q: (lambda a: (a, x_for_aphi(a, q.uniform(28.0, 36.0), 1)))(
            full_mantissa(q, 117.0, 190.0))),
        1500, seam="temme_split", note="#70: R3 core just below a*phi = 36 near lambda = 2")

    # 6. Saturation a*phi = 800 (R3) and E = -800 (R1, R2).
    a_set = (5e3, 1e4, 1e5, 1e6, 1e7)
    add("seam:expfloor_r3:bits",
        lambda r, n, a_set=a_set: _take(r, [(a, x) for a in a_set for s in (-1, 1)
                               for x in neighbours(x_for_aphi_exact(a, 800.0, s), 64)], n),
        129 * 2 * len(a_set), seam="expfloor_r3")
    for lo, hi, tag in ((760.0, 800.0, "lo"), (800.0, 840.0, "hi")):
        add(f"seam:expfloor_r3:band_{tag}",
            lambda r, n, lo=lo, hi=hi: _many(r, n, lambda q: (lambda a: (a, x_for_aphi(a, q.uniform(lo, hi), q.choice((-1, 1)))))(
                log_uniform(q, 5e3, 1e7))),
            1000, seam="expfloor_r3")
    for name, r1, arange in (("expfloor_r1", True, (1.2, 1e3)), ("expfloor_r2", False, (1e-3, 1e3))):
        add(f"seam:{name}:band",
            lambda r, n, r1=r1, arange=arange: _many(r, n, lambda q: (lambda a: (a, x_for_E(a, q.uniform(-840.0, -760.0), r1)))(
                log_uniform(q, *arange))),
            1500, seam=name)
        a_bits = (2.0, 7.5, 19.5, 40.0) if r1 else (0.5, 3.0, 19.5, 40.0)
        add(f"seam:{name}:bits",
            lambda r, n, r1=r1, a_bits=a_bits: _take(r, [(a, x) for a in a_bits
                                                         for x in neighbours(x_for_E_exact(a, -800.0, r1) or 1.0, 64)], n),
            129 * 4, seam=name)

    # 7. Subnormal outputs: R3 a*phi in [700, 760]; R1/R2 exponent in [-750, -700].
    add("subnormal_out:r3",
        lambda r, n: _many(r, n, lambda q: (lambda a: (a, x_for_aphi(a, q.uniform(700.0, 760.0), q.choice((-1, 1)))))(
            log_uniform(q, 5e3, 1e7))),
        2500, seam="out_min_normal")
    for name, r1, arange in (("r1", True, (1.2, 1e4)), ("r2", False, (1e-3, 1e4))):
        add(f"subnormal_out:{name}",
            lambda r, n, r1=r1, arange=arange: _many(r, n, lambda q: (lambda a: (a, x_for_E(a, q.uniform(-750.0, -700.0), r1)))(
                log_uniform(q, *arange))),
            1500, seam="out_min_normal")

    # 8. Subnormal a (R4 with Q ~ a*E1(x) subnormal; #71) and tiny a at the box edges.
    add("subnormal_a:x_log",
        lambda r, n: [(_sub(r), log_uniform(r, 1e-6, 60.0)) for _ in range(n)], 2000)
    add("subnormal_a:box",
        lambda r, n: [(_sub(r), full_mantissa(r, 1e-3, 4.0)) for _ in range(n)], 2500,
        note="#71: subnormal a, x in the R4 box")
    add("tiny_a:x_near_4",
        lambda r, n: [(log_uniform(r, 1e-308, 1e-8), full_mantissa(r, 3.8, 4.2)) for _ in range(n)], 800,
        seam="r4_xmax")
    add("tiny_a:x_near_a_plus_1",
        lambda r, n: [(log_uniform(r, 1e-308, 1e-8), full_mantissa(r, 0.95, 1.05)) for _ in range(n)], 800,
        seam="x_eq_a_plus_1")

    # 9. Crossed extremes.
    add("extreme:big_a_ridge",
        lambda r, n: _many(r, n, lambda q: (lambda a: (a, x_for_aphi(a, q.uniform(0.0, 40.0), q.choice((-1, 1)))))(
            log_uniform(q, 1e6, 1e7))),
        600, note="a in [1e6, 1e7] on the ridge (costliest oracle rows)")
    add("extreme:huge_a_saturated",
        lambda r, n: _many(r, n, lambda q: (lambda a: (a, a * (log_uniform(q, 1e-300, 0.45) if q.random() < 0.5
                                                               else log_uniform(q, 2.2, 1e300))))(
            math.ldexp(1.0 + q.random(), q.randint(24, 1022)))),
        800)
    add("extreme:diag_huge",
        lambda r, n: _take(r, [(a, x) for k in range(100, 1024, 3) for a in (math.ldexp(1.0 + r.random() * 0.999, k),)
                               for x in ((a,) if k < 120 else (a, math.nextafter(a, math.inf), math.nextafter(a, 0)))], n),
        3 * len(range(100, 1024, 3)),
        note="x == a for a >= 2^100: Q = 1/2 - 1/(3 sqrt(2 pi a)); +-1 ulp (saturated) from 2^120")
    add("seam:twopi_a_clamp:bits",
        lambda r, n: _take(r, [(a, x) for a in neighbours(2.0 ** 900, 64)
                               for x in (a, math.nextafter(a, math.inf), math.nextafter(a, 0))], n),
        129 * 3, seam="twopi_a_clamp")
    add("extreme:subnormal_x",
        lambda r, n: [(log_uniform(r, 1e-3, 30.0), _sub(r)) for _ in range(n)], 800,
        seam="out_zero")
    add("extreme:huge_x",
        lambda r, n: [(log_uniform(r, 1e-300, 1e6), log_uniform(r, 1e10, 1.7e308)) for _ in range(n)], 400)

    # 10. Generic full-mantissa coverage per region.
    add("generic:r1_small",
        lambda r, n: _many(r, n, lambda q: (lambda a: (a, full_mantissa(q, 1e-3, a + 1.0)))(full_mantissa(q, 1.5, 20.0))), 1000)
    add("generic:r2_small",
        lambda r, n: _many(r, n, lambda q: (lambda a: (a, (a + 1.0) * log_uniform(q, 1.0, 50.0)))(full_mantissa(q, 0.01, 20.0))), 1000)
    add("generic:r3",
        lambda r, n: _many(r, n, lambda q: (lambda a: (a, a * full_mantissa(q, 0.5, 2.0)))(log_uniform(q, 20.0, 1e5))), 1500)
    add("generic:r4",
        lambda r, n: [(log_uniform(r, 1e-12, 1.5), full_mantissa(r, 1e-6, 4.0)) for _ in range(n)], 1000)

    # 11. Known witnesses (validation of the strata).
    wit = [(float.fromhex("0x1.eb4e156a22059p+6"), float.fromhex("0x1.e35ffc4f3e367p+7")),
           (float.fromhex("0x0.0000024de7ffep-1022"), float.fromhex("0x1.db1b6d2972468p-2"))]
    add("witness", lambda r, n: list(wit), len(wit), note="#70 (R3 core 3 ULP), #71 (subnormal a, SSE4/SSE2)")
    return out
