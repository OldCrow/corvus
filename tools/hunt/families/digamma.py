"""digamma hunt strata. Bounds: ACCURACY.md digamma --
  x > 0                          1 ULP
  x < 0, |psi| >= 1              1 ULP
  x < 0, |psi| < 1               2^-53 absolute
psi(+-0) = -+inf, negative integers -> NaN. ACCURACY.md says positive
subnormals -> -inf; the real value -1/x is finite for x >= ~2^-1024, so the
top two subnormal binades are scored against the finite exact value."""
import math

import mpmath as mp

from core import Stratum, full_mantissa, log_uniform, neighbours, subnormal

_ROOT = float.fromhex("0x1.762d86356be3fp+0")  # kDigammaRootHi
# Negative zeros, one per (-k-1, -k), mpmath findroot at dps 40.
_NEG_ZEROS = [float.fromhex(h) for h in (
    "-0x1.02172b05ee260p-1", "-0x1.92d0cbc289d4ap+0", "-0x1.4e2c19f679e5ap+1",
    "-0x1.d1514b041b2a8p+1", "-0x1.29cea5c1ccbd0p+2", "-0x1.6ab2ca18e6ce3p+2")]

SEAMS = [
    {"fn": "digamma", "name": "zone_lo", "coord": "x", "value": 1.0, "source": "kDigammaZoneLo (up-step|zone)"},
    {"fn": "digamma", "name": "root", "coord": "x", "value": _ROOT, "source": "kDigammaRootHi (t = x - x0)"},
] + [
    {"fn": "digamma", "name": f"walk_{j}", "coord": "x", "value": float(j), "source": f"walk step fire x >= {j}"}
    for j in range(2, 8)
] + [
    {"fn": "digamma", "name": "x0", "coord": "x", "value": 8.0, "source": "kDigammaX0 (walk|asymptotic)"},
    {"fn": "digamma", "name": "sq_ovf", "coord": "x", "value": 2.0 ** 512, "source": "x*x overflow avoided (w = (1/x)^2)"},
    {"fn": "digamma", "name": "asym_cut", "coord": "x", "value": 2.0 ** 85, "source": "kDigammaAsymCut (log_dd alone)"},
    {"fn": "digamma", "name": "direct_recip", "coord": "x", "value": 2.0 ** -960, "source": "kDigammaDirectRecip (-1/x alone)"},
    {"fn": "digamma", "name": "neg_direct_recip", "coord": "x", "value": -2.0 ** -960, "source": "kDigammaDirectRecip, |u| side"},
    {"fn": "digamma", "name": "out_ovf", "coord": "x", "value": 2.0 ** -1024, "source": "-1/x overflows below ~2^-1024"},
    {"fn": "digamma", "name": "min_normal_in", "coord": "x", "value": 2.0 ** -1022, "source": "subnormal input edge"},
    {"fn": "digamma", "name": "neg_tie_0_5", "coord": "x", "value": -0.5, "source": "u = x - round(x) sign flip"},
    {"fn": "digamma", "name": "neg_tie_1_5", "coord": "x", "value": -1.5, "source": "u = x - round(x) sign flip"},
    {"fn": "digamma", "name": "pole_1", "coord": "x", "value": -1.0, "source": "pole (NaN)"},
    {"fn": "digamma", "name": "pole_2", "coord": "x", "value": -2.0, "source": "pole (NaN)"},
    {"fn": "digamma", "name": "rough_floor", "coord": "x", "value": -5.5, "source": "kDigammaRoughTrigammaFloor: y = 1-x = 6.5"},
    {"fn": "digamma", "name": "neg_half_grid", "coord": "x", "value": -2.0 ** 52, "source": "|x| >= 2^52: integers and halves"},
    {"fn": "digamma", "name": "neg_int_grid", "coord": "x", "value": -2.0 ** 53, "source": "|x| >= 2^53: every double a pole"},
] + [
    {"fn": "digamma", "name": f"neg_zero_{i}", "coord": "x", "value": z, "source": "psi(x) = 0, x < 0 (absolute band)"}
    for i, z in enumerate(_NEG_ZEROS)
]

_BANDS = {"neg_half_grid": (-2.0 ** 52 * 1.05, -2.0 ** 52 * 0.95),
          "root": (_ROOT - 0.02, _ROOT + 0.02)}
for _i, _z in enumerate(_NEG_ZEROS):
    _BANDS[f"neg_zero_{_i}"] = (_z - 2e-3, _z + 2e-3)


def _special(a):
    """Exact value for the specials, or None."""
    if a == 0:
        return -mp.inf if math.copysign(1.0, a) > 0 else mp.inf
    if a < 0 and a == math.floor(a):
        return mp.nan
    return None


def _extra(x):
    if abs(x - _ROOT) < 0.05:
        return 40
    if x < 0:
        return 25   # reflection cancels ~49 bits near the negative zeros
    return 15


def _oracle(args, dps):
    sp = _special(args[0])
    if sp is not None:
        return sp
    x = mp.mpf(args[0])
    with mp.workdps(dps + _extra(x)):
        if x < 0:
            return mp.digamma(1 - x) - mp.pi * mp.cospi(x) / mp.sinpi(x)
        return mp.digamma(x)


def _psi_series(x):
    """psi(x), x > 0: upward recurrence to y >= N, then the Bernoulli series.
    Shares no mpmath digamma internals."""
    n_min = int(mp.mp.dps * 1.2) + 10
    s = mp.mpf(0)
    y = x
    if y < n_min:
        k = int(mp.ceil(n_min - y))
        for j in range(k):
            s += 1 / (x + j)
        y = x + k
    r = mp.log(y) - 1 / (2 * y)
    y2 = y * y
    yp = y2
    eps = mp.mpf(2) ** -(mp.mp.prec + 10)
    for k in range(1, 200):
        term = mp.bernoulli(2 * k) / (2 * k * yp)
        r -= term
        if abs(term) < eps * abs(r):
            return r - s
        yp *= y2
    raise RuntimeError("psi series did not converge")


def _oracle2(args, dps):
    sp = _special(args[0])
    if sp is not None:
        return sp
    x = mp.mpf(args[0])
    with mp.workdps(dps + _extra(x) + 10):
        if x < 0:
            return _psi_series(1 - x) - mp.pi / mp.tan(mp.pi * x)
        return _psi_series(x)


def _bound(args, v):
    return 1


_ABS = mp.mpf(2) ** -53


def _abs_tol(args, v):
    if args[0] < 0 and mp.isfinite(v) and abs(v) < 1:
        return _ABS
    return None


def _st(name, gen, n, seam=None, note=""):
    return Stratum("digamma", name, gen, n, _oracle, _bound, oracle2=_oracle2,
                   abs_tol=_abs_tol, seam=seam, note=note)


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
        out.append(_st(f"seam:{s['name']}:bits", lambda r, n, t=t: [(x,) for x in neighbours(t, 64)],
                       129, s["name"]))
        out.append(_st(f"seam:{s['name']}:band",
                       lambda r, n, lo=lo, hi=hi: [(full_mantissa(r, lo, hi),) for _ in range(n)],
                       800, s["name"]))
    out += [
        _st("pos_upstep", lambda r, n: [(full_mantissa(r, 0.0, 1.0),) for _ in range(n)], 2500),
        _st("pos_tiny", lambda r, n: [(log_uniform(r, 2.0 ** -1022, 1e-3),) for _ in range(n)], 1500),
        _st("pos_zone", lambda r, n: [(full_mantissa(r, 1.0, 2.0),) for _ in range(n)], 3000),
        _st("near_root", lambda r, n: [(_ROOT + (1 if r.random() < 0.5 else -1) * log_uniform(r, 2.0 ** -52, 1e-2),)
                                       for _ in range(n)], 2000, note="relative bound at the root"),
        _st("pos_walk", lambda r, n: [(full_mantissa(r, 2.0, 8.0),) for _ in range(n)], 3000),
        _st("pos_asym", lambda r, n: [(log_uniform(r, 8.0, 1.7e308),) for _ in range(n)], 2500),
        _st("subnormal", lambda r, n: [((subnormal(r) if r.random() < 0.5 else -subnormal(r)),) for _ in range(n)],
            1500, note="documented +-inf except top binades (finite exact)"),
        _st("subnormal_top", lambda r, n: [(subnormal(r, 2),) for _ in range(n)], 600,
            note="x >= 2^-1024: -1/x finite"),
        _st("neg_unit_intervals", lambda r, n: [(-r.randint(0, 19) - full_mantissa(r, 0.0, 1.0),) for _ in range(n)],
            6000, note="between poles, -20 < x < 0"),
        _st("neg_near_poles", _near_poles, 2000),
        _st("neg_large", lambda r, n: [(-log_uniform(r, 20.0, 2.0 ** 52),) for _ in range(n)], 2000),
    ]
    return out
