"""trigamma hunt strata. Bounds: ACCURACY.md trigamma -- 1 ULP everywhere
(relative on both axes, no absolute band). 0 < x < 1 is documented CR, but
#76 corrects that cell to 1 ULP; scored at bound 1 here (see notes).
+inf at +-0, at negative integers, and for every subnormal (1/x^2 overflows)."""
import math

import mpmath as mp

from core import Stratum, full_mantissa, log_uniform, neighbours, subnormal

_TRI_MIN = float.fromhex("-0x1.fb9ce3f1d2f24p-2")  # global minimum on x < 0

SEAMS = [
    {"fn": "trigamma", "name": "zone_lo", "coord": "x", "value": 1.0, "source": "kTrigammaZoneLo (up-step|zone)"},
    {"fn": "trigamma", "name": "zone_centre", "coord": "x", "value": 1.5, "source": "kTrigammaZoneCentre"},
] + [
    {"fn": "trigamma", "name": f"walk_{j}", "coord": "x", "value": float(j), "source": f"walk step fire x >= {j}"}
    for j in range(2, 8)
] + [
    {"fn": "trigamma", "name": "x0", "coord": "x", "value": 8.0, "source": "kTrigammaX0 (walk|asymptotic)"},
    {"fn": "trigamma", "name": "asym_cut", "coord": "x", "value": 2.0 ** 89, "source": "kTrigammaAsymCut (fl(1/x) alone)"},
    {"fn": "trigamma", "name": "sq_ovf", "coord": "x", "value": 2.0 ** 512, "source": "x*x overflow; 1/x^2 subnormal output"},
    {"fn": "trigamma", "name": "tiny_guard", "coord": "x", "value": 2.0 ** -480, "source": "kTrigammaDeepTinyGuard"},
    {"fn": "trigamma", "name": "neg_tiny_guard", "coord": "x", "value": -2.0 ** -480, "source": "kTrigammaDeepTinyGuard, |u| side"},
    {"fn": "trigamma", "name": "out_ovf", "coord": "x", "value": 2.0 ** -512, "source": "1/x^2 overflows (2^-512 = 1/sqrt(DBL_MAX))"},
    {"fn": "trigamma", "name": "neg_out_ovf", "coord": "x", "value": -2.0 ** -512, "source": "1/u^2 overflows, x < 0"},
    {"fn": "trigamma", "name": "min_normal_in", "coord": "x", "value": 2.0 ** -1022, "source": "subnormal input edge"},
    {"fn": "trigamma", "name": "neg_tie_0_5", "coord": "x", "value": -0.5, "source": "u = x - round(x) sign flip"},
    {"fn": "trigamma", "name": "neg_tie_1_5", "coord": "x", "value": -1.5, "source": "u = x - round(x) sign flip"},
    {"fn": "trigamma", "name": "neg_min", "coord": "x", "value": _TRI_MIN, "source": "global minimum ~8.933 on x < 0"},
    {"fn": "trigamma", "name": "pole_1", "coord": "x", "value": -1.0, "source": "double pole (+inf)"},
    {"fn": "trigamma", "name": "pole_2", "coord": "x", "value": -2.0, "source": "double pole (+inf)"},
    {"fn": "trigamma", "name": "rough_floor", "coord": "x", "value": -5.5, "source": "kTrigammaRoughTetraFloor: y = 1-x = 6.5"},
    {"fn": "trigamma", "name": "neg_half_grid", "coord": "x", "value": -2.0 ** 52, "source": "|x| >= 2^52: integers and halves"},
    {"fn": "trigamma", "name": "neg_int_grid", "coord": "x", "value": -2.0 ** 53, "source": "|x| >= 2^53: every double a pole"},
]

_BANDS = {"neg_half_grid": (-2.0 ** 52 * 1.05, -2.0 ** 52 * 0.95)}


def _special(a):
    if a == 0 or (a < 0 and a == math.floor(a)):
        return mp.inf
    return None


def _oracle(args, dps):
    sp = _special(args[0])
    if sp is not None:
        return sp
    x = mp.mpf(args[0])
    with mp.workdps(dps + 15):
        if x < 0:  # mpmath's negative-axis polygamma is O(|x|); reflect
            return (mp.pi / mp.sinpi(x)) ** 2 - mp.polygamma(1, 1 - x)
        return mp.polygamma(1, x)


def _oracle2(args, dps):
    sp = _special(args[0])
    if sp is not None:
        return sp
    x = mp.mpf(args[0])
    with mp.workdps(dps + 25):
        if x < 0:
            return (mp.pi / mp.sin(mp.pi * x)) ** 2 - mp.zeta(2, 1 - x)
        return mp.zeta(2, x)  # Hurwitz zeta


def _bound(args, v):
    return 1


def _st(name, gen, n, seam=None, note=""):
    return Stratum("trigamma", name, gen, n, _oracle, _bound, oracle2=_oracle2, seam=seam, note=note)


def _near_poles(r, n):
    out = []
    for _ in range(n):
        k = r.randint(1, 30)
        d = log_uniform(r, 2.0 ** -52 * k, 1e-2)
        out.append((-k + (d if r.random() < 0.5 else -d),))
    return out


def strata():
    out = []
    for s in SEAMS:
        t = s["value"]
        lo, hi = _BANDS.get(s["name"], tuple(sorted((t * 0.95, t * 1.05))))
        note = "#76: (0,1) documented CR, scored at 1 ULP" if 0 < t < 1 else ""
        out.append(_st(f"seam:{s['name']}:bits", lambda r, n, t=t: [(x,) for x in neighbours(t, 64)],
                       129, s["name"], note))
        out.append(_st(f"seam:{s['name']}:band",
                       lambda r, n, lo=lo, hi=hi: [(full_mantissa(r, lo, hi),) for _ in range(n)],
                       800, s["name"], note))
    out += [
        _st("pos_upstep", lambda r, n: [(full_mantissa(r, 0.0, 1.0),) for _ in range(n)], 4000,
            note="#76: documented CR, scored at 1 ULP"),
        _st("pos_tiny", lambda r, n: [(log_uniform(r, 2.0 ** -1022, 1e-3),) for _ in range(n)], 2500,
            note="#76 cell; deep-tiny 1/x^2 path below 2^-480"),
        _st("pos_zone", lambda r, n: [(full_mantissa(r, 1.0, 2.0),) for _ in range(n)], 3000),
        _st("zone_edges", lambda r, n: [((1.0 if r.random() < 0.5 else 2.0)
                                         + (1 if r.random() < 0.5 else -1) * log_uniform(r, 2.0 ** -52, 1e-6),)
                                        for _ in range(n)], 1500, note="fit worst points within ~1e-10 of edges"),
        _st("pos_walk", lambda r, n: [(full_mantissa(r, 2.0, 8.0),) for _ in range(n)], 3000),
        _st("pos_asym", lambda r, n: [(log_uniform(r, 8.0, 1.7e308),) for _ in range(n)], 2500,
            note="subnormal outputs above ~2^1022"),
        _st("subnormal", lambda r, n: [((subnormal(r) if r.random() < 0.5 else -subnormal(r)),) for _ in range(n)],
            1000, note="+inf"),
        _st("neg_unit_intervals", lambda r, n: [(-r.randint(0, 19) - full_mantissa(r, 0.0, 1.0),) for _ in range(n)],
            6000, note="between poles, -20 < x < 0"),
        _st("neg_near_poles", _near_poles, 2000),
        _st("neg_large", lambda r, n: [(-log_uniform(r, 20.0, 2.0 ** 52),) for _ in range(n)], 2000),
    ]
    return out
