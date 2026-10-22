"""erfinv / erfcinv hunt strata. Bound: ACCURACY.md erfinv and erfcinv -- max
1 ULP everywhere, both functions, far tail and subnormal arguments included.
Known: erfinv of a subnormal is 2 ULP on SSE4/SSE2 (#73,
erfinv(-0x0.2884b16852f8ep-1022)); the subnormal strata below reproduce it."""
import math

import mpmath as mp

from core import Stratum, full_mantissa, log_uniform, neighbours, subnormal

_YSPLIT = 0.5                                       # kErfinvYSplit
_TSPLIT_S = math.exp(-4.0)                          # t = kErfinvTSplit = 2 <=> s = e^-4
_TFAR = float.fromhex("0x1.8c7a44436a6f6p+2")       # kErfinvTFar
_TFAR_S = float.fromhex("0x1.8cf81557d209ap-56")    # s = exp(-TFar^2)
_ERFC6 = float.fromhex("0x1.8cf81557d20b6p-56")     # erfc(6): root x = 6 (kErfcCoreSafeMax clamp)
_ISSUE73 = -float.fromhex("0x0.2884b16852f8ep-1022")

SEAMS = [
    {"fn": "erfinv", "name": "ysplit", "coord": "y", "value": _YSPLIT, "source": "kErfinvYSplit (C|T)"},
    {"fn": "erfinv", "name": "tsplit", "coord": "y", "value": 1 - _TSPLIT_S, "source": "kErfinvTSplit t=2 <=> s=e^-4 (seed interval 0|1)"},
    {"fn": "erfinv", "name": "y_to_1", "coord": "y", "value": 1 - 2.0 ** -53, "source": "largest y < 1 (s = 2^-53, x ~ 5.86)"},
    {"fn": "erfinv", "name": "min_normal_in", "coord": "y", "value": 2.0 ** -1022, "source": "subnormal input edge"},
    {"fn": "erfcinv", "name": "c_lo", "coord": "z", "value": 0.5, "source": "erfcinv C|T split z = 1/2"},
    {"fn": "erfcinv", "name": "c_hi", "coord": "z", "value": 1.5, "source": "erfcinv C|T split z = 3/2"},
    {"fn": "erfcinv", "name": "zero_cross", "coord": "z", "value": 1.0, "source": "erfcinv(1) = +0, C(1 - z)"},
    {"fn": "erfcinv", "name": "tsplit_lo", "coord": "z", "value": _TSPLIT_S, "source": "kErfinvTSplit t=2 <=> s=e^-4"},
    {"fn": "erfcinv", "name": "tsplit_hi", "coord": "z", "value": 2 - _TSPLIT_S, "source": "kErfinvTSplit via 2 - z"},
    {"fn": "erfcinv", "name": "tfar", "coord": "z", "value": _TFAR_S, "source": "kErfinvTFar (Halley mid|far, seed interval 1|2)"},
    {"fn": "erfcinv", "name": "x6_root", "coord": "z", "value": _ERFC6, "source": "erfc(6): kErfcCoreSafeMax clamp / erfc core|tail"},
    {"fn": "erfcinv", "name": "z_to_2", "coord": "z", "value": 2 - 2.0 ** -52, "source": "largest z < 2 (s = 2^-52)"},
    {"fn": "erfcinv", "name": "min_normal_in", "coord": "z", "value": 2.0 ** -1022, "source": "LogDdAny subnormal s edge"},
    {"fn": "erfcinv", "name": "min_subnormal", "coord": "z", "value": 2.0 ** -1074, "source": "smallest s, x ~ 27.2"},
]

# Bands per seam (lo, hi); default is +-5% of the seam value.
_BANDS = {
    "tsplit": (1 - 1.1 * _TSPLIT_S, 1 - 0.9 * _TSPLIT_S),
    "y_to_1": (1 - 1e-13, 1 - 2.0 ** -53),
    "tsplit_hi": (2 - 1.1 * _TSPLIT_S, 2 - 0.9 * _TSPLIT_S),
    "tfar": (_TFAR_S * 0.5, _TFAR_S * 2),
    "x6_root": (_ERFC6 * 0.9, _ERFC6 * 1.1),
    "z_to_2": (2 - 1e-12, 2 - 2.0 ** -52),
    "min_subnormal": (2.0 ** -1074, 2.0 ** -1060),
}


# ---- oracles ----------------------------------------------------------------
def _tail_guess(s):
    w = -mp.log(s)
    if w > 2:
        return mp.sqrt(w - mp.log(mp.pi * w) / 2)
    with mp.workdps(20):
        return mp.erfinv(1 - s)


def _t_log(s, dps):
    """Positive root of erfc(x) = s, Newton on log(erfc(x)) - log(s)."""
    with mp.workdps(dps + 15):
        s = mp.mpf(s)
        x, ls = _tail_guess(s), mp.log(s)
        tol = mp.mpf(10) ** -(dps + 8)
        for _ in range(100):
            e = mp.erfc(x)
            dx = (mp.log(e) - ls) * mp.sqrt(mp.pi) * e / (2 * mp.exp(-x * x))
            x += dx
            if abs(dx) <= tol * x:
                return x
    raise RuntimeError("log-space Newton did not converge")


def _t_lin(s, dps):
    """Same root, Newton on erfc(x) - s (linear space; mpmath erfc is
    relative-accurate in the tail)."""
    with mp.workdps(dps + 15):
        s = mp.mpf(s)
        x = _tail_guess(s)
        tol = mp.mpf(10) ** -(dps + 8)
        for _ in range(200):
            dx = (mp.erfc(x) - s) * mp.sqrt(mp.pi) / (2 * mp.exp(-x * x))
            x += dx
            if abs(dx) <= tol * x:
                return x
    raise RuntimeError("linear Newton did not converge")


def _erf_newton(y, dps):
    """erfinv(y), |y| <= 1/2, by Newton on erf from a low-precision start."""
    with mp.workdps(dps + 15):
        y = mp.mpf(y)
        if y == 0:
            return y
        with mp.workdps(15):
            x = mp.erfinv(y)
        tol = mp.mpf(10) ** -(dps + 8)
        for _ in range(60):
            dx = (mp.erf(x) - y) * mp.sqrt(mp.pi) / (2 * mp.exp(-x * x))
            x -= dx
            if abs(dx) <= tol * abs(x):
                return x
    raise RuntimeError("erf Newton did not converge")


def _erfinv(args, dps):
    y = mp.mpf(args[0])
    if abs(y) >= 1:
        return mp.nan if abs(y) > 1 else mp.sign(y) * mp.inf
    with mp.workdps(dps + 10):
        return mp.erfinv(y)


def _erfinv2(args, dps):
    y = mp.mpf(args[0])
    if abs(y) >= 1:
        return _erfinv(args, dps)
    if abs(y) <= 0.5:
        return _erf_newton(y, dps)
    return mp.sign(y) * _t_lin(1 - abs(y), dps)  # 1 - |y| exact (Sterbenz)


def _erfcinv(args, dps):
    z = mp.mpf(args[0])
    if z < 0 or z > 2:
        return mp.nan
    if z == 0:
        return mp.inf
    if z == 2:
        return -mp.inf
    if z < 0.5:
        return _t_log(z, dps)
    if z > 1.5:
        return -_t_log(2 - z, dps)                 # 2 - z exact in mpf
    with mp.workdps(dps + 10):
        return mp.erfinv(1 - z)                    # 1 - z exact (Sterbenz)


def _erfcinv2(args, dps):
    z = mp.mpf(args[0])
    if z <= 0 or z >= 2:
        return _erfcinv(args, dps)
    if z < 0.5:
        return _t_lin(z, dps)
    if z > 1.5:
        return -_t_lin(2 - z, dps)
    return _erf_newton(1 - z, dps)


def _bound(args, v):
    return 1


# ---- strata -----------------------------------------------------------------
def _sgn(r, x):
    return -x if r.random() < 0.5 else x


def _st(fn, name, gen, n, seam=None, note=""):
    o, o2 = (_erfinv, _erfinv2) if fn == "erfinv" else (_erfcinv, _erfcinv2)
    return Stratum(fn, name, gen, n, o, _bound, oracle2=o2, seam=seam, note=note)


def _y_tail(r, lo_s, hi_s):
    """y = +-(1 - s), s log-uniform; y rounds to the 2^-53 grid below 1."""
    return _sgn(r, 1 - log_uniform(r, lo_s, hi_s))


def strata():
    out = []
    for s in SEAMS:
        fn, t = s["fn"], s["value"]
        lo, hi = _BANDS.get(s["name"], (t * 0.95, t * 1.05))
        sym = fn == "erfinv"
        ok = (lambda x: -1 < x < 1) if sym else (lambda x: 0 < x < 2)

        def bits_gen(r, n, t=t, sym=sym, ok=ok):
            pts = [x for x in neighbours(t, 64) if ok(x)]
            return [(x,) for x in pts] + ([(-x,) for x in pts] if sym else [])

        def band_gen(r, n, lo=lo, hi=hi, sym=sym, ok=ok):
            pts = [full_mantissa(r, lo, hi) for _ in range(n)]
            return [((_sgn(r, x) if sym else x),) for x in pts if ok(x)]

        out.append(_st(fn, f"seam:{s['name']}:bits", bits_gen, 129, s["name"]))
        out.append(_st(fn, f"seam:{s['name']}:band", band_gen, 1000, s["name"]))
    out += [
        # erfinv
        _st("erfinv", "central", lambda r, n: [(_sgn(r, full_mantissa(r, 0.0, 0.5)),) for _ in range(n)], 5000),
        _st("erfinv", "central_tiny", lambda r, n: [(_sgn(r, log_uniform(r, 2.0 ** -1022, 1e-3)),) for _ in range(n)], 2000),
        _st("erfinv", "tail_mid", lambda r, n: [(_sgn(r, full_mantissa(r, 0.5, 1 - _TSPLIT_S)),) for _ in range(n)], 4000),
        _st("erfinv", "tail_deep", lambda r, n: [(_y_tail(r, 2.0 ** -53, _TSPLIT_S),) for _ in range(n)], 4000),
        _st("erfinv", "subnormal", lambda r, n: [(_sgn(r, subnormal(r)),) for _ in range(n)], 6000,
            note="#73: 2 ULP on SSE4/SSE2"),
        _st("erfinv", "subnormal_top", lambda r, n: [(_sgn(r, subnormal(r, 4)),) for _ in range(n)], 3000,
            note="#73 neighbourhood: top subnormal binades"),
        _st("erfinv", "issue73:bits", lambda r, n: [(x,) for x in neighbours(_ISSUE73, 64)], 129,
            note="#73 witness erfinv(-0x0.2884b16852f8ep-1022)"),
        # erfcinv
        _st("erfcinv", "central", lambda r, n: [(full_mantissa(r, 0.5, 1.5),) for _ in range(n)], 5000),
        _st("erfcinv", "near_one", lambda r, n: [(1 + _sgn(r, log_uniform(r, 2.0 ** -53, 1e-3)),) for _ in range(n)], 2000,
            note="results down to ~1e-16"),
        _st("erfcinv", "tail_seed0", lambda r, n: [(log_uniform(r, _TSPLIT_S, 0.5),) for _ in range(n)], 3000),
        _st("erfcinv", "tail_seed1", lambda r, n: [(log_uniform(r, _TFAR_S, _TSPLIT_S),) for _ in range(n)], 4000),
        _st("erfcinv", "tail_far", lambda r, n: [(log_uniform(r, 2.0 ** -1022, _TFAR_S),) for _ in range(n)], 4000),
        _st("erfcinv", "upper", lambda r, n: [(2 - log_uniform(r, 2.0 ** -52, 0.5),) for _ in range(n)], 3000),
        _st("erfcinv", "subnormal", lambda r, n: [(subnormal(r),) for _ in range(n)], 3000),
    ]
    return out
