"""lgamma hunt strata. Bounds: ACCURACY.md lgamma --
  0 < x <= 5/2                   1 ULP
  5/2 < x < 8 (table)            0 ULP (correctly rounded, near-tie window TIE_TABLE)
  x >= 8 (Stirling)              0 ULP (correctly rounded, near-tie window TIE)
  x < 0, |lgamma| >= 1           1 ULP
  x < 0, |lgamma| < 1            2^-53 absolute
Poles (x <= 0 integer, and |x| > kLgammaMaxArg) give +inf."""
import mpmath as mp

from core import Stratum, full_mantissa, log_uniform, neighbours, subnormal

_MAXARG = float.fromhex("0x1.754d9278b51a7p+1014")   # kLgammaMaxArg
# Negative zeros of lgamma (|Gamma(x)| = 1), mpmath findroot at dps 40.
_NEG_ZEROS = [float.fromhex(h) for h in (
    "-0x1.3a7fc9600f86cp+1", "-0x1.5fb410a1bd901p+1", "-0x1.9260dbc9e59afp+1",
    "-0x1.fa471547c2fe5p+1", "-0x1.0284e78599581p+2", "-0x1.3f7577a6eeafdp+2",
    "-0x1.4086a57f0b6d9p+2")]

SEAMS = [
    {"fn": "lgamma", "name": "zone_lo", "coord": "x", "value": 0.5, "source": "kLgammaZoneLo (small|zone)"},
    {"fn": "lgamma", "name": "zero_1", "coord": "x", "value": 1.0, "source": "zero at x=1, t = x-1 exact"},
    {"fn": "lgamma", "name": "zone_mid", "coord": "x", "value": 1.5, "source": "kLgammaZoneMid (centre 1|2)"},
    {"fn": "lgamma", "name": "zero_2", "coord": "x", "value": 2.0, "source": "zero at x=2, t = x-2 exact"},
    {"fn": "lgamma", "name": "zone_hi", "coord": "x", "value": 2.5, "source": "kLgammaZoneHi (zone|table)"},
] + [
    {"fn": "lgamma", "name": f"cell_edge_{k}", "coord": "x", "value": 2.5 + k + 1 / 64, "source": "table cell edge, |r| = 1/64 (LgammaTable)"}
    for k in range(0, 5)
] + [
    {"fn": "lgamma", "name": "x0", "coord": "x", "value": 8.0, "source": "kLgammaX0 (table|Stirling)"},
    {"fn": "lgamma", "name": "stir_scale", "coord": "x", "value": 2.0 ** 200, "source": "LgammaStirling 2^-200 product scaling"},
    {"fn": "lgamma", "name": "dekker", "coord": "x", "value": 2.0 ** 996, "source": "non-FMA ProdLow Dekker split ceiling 2^996"},
    {"fn": "lgamma", "name": "max_arg", "coord": "x", "value": _MAXARG, "source": "kLgammaMaxArg (overflow)"},
    {"fn": "lgamma", "name": "min_normal_in", "coord": "x", "value": 2.0 ** -1022, "source": "subnormal input edge"},
    {"fn": "lgamma", "name": "neg_min_normal_in", "coord": "x", "value": -2.0 ** -1022, "source": "subnormal input edge, x<0"},
    {"fn": "lgamma", "name": "neg_tiny_scale", "coord": "x", "value": -2.0 ** -400, "source": "reflection product scales by 2^600 below |x| = 2^-400"},
    {"fn": "lgamma", "name": "neg_tie_0_5", "coord": "x", "value": -0.5, "source": "u = x - round(x) sign flip"},
    {"fn": "lgamma", "name": "neg_tie_2_5", "coord": "x", "value": -2.5, "source": "u = x - round(x) sign flip"},
    {"fn": "lgamma", "name": "neg_tie_3_5", "coord": "x", "value": -3.5, "source": "u = x - round(x) sign flip"},
    {"fn": "lgamma", "name": "pole_1", "coord": "x", "value": -1.0, "source": "pole, u = 0"},
    {"fn": "lgamma", "name": "pole_3", "coord": "x", "value": -3.0, "source": "pole, u = 0"},
    {"fn": "lgamma", "name": "neg_half_grid", "coord": "x", "value": -2.0 ** 52, "source": "|x| >= 2^52: only integers and halves"},
    {"fn": "lgamma", "name": "neg_int_grid", "coord": "x", "value": -2.0 ** 53, "source": "|x| >= 2^53: every double a pole"},
] + [
    {"fn": "lgamma", "name": f"neg_zero_{i}", "coord": "x", "value": z, "source": "lgamma(x) = 0, x < 0 (absolute band)"}
    for i, z in enumerate(_NEG_ZEROS)
]

_BANDS = {
    "max_arg": (_MAXARG * 0.95, _MAXARG),          # above it: +inf by construction
    "neg_half_grid": (-2.0 ** 52 * 1.05, -2.0 ** 52 * 0.95),
}
for _i, _z in enumerate(_NEG_ZEROS):
    _BANDS[f"neg_zero_{_i}"] = (_z - 2e-3, _z + 2e-3)


def _is_pole(x):
    return x <= 0 and x == mp.floor(x)


def _extra(x):
    # Relative accuracy near the positive zeros needs the cancellation covered.
    return 30 if (abs(x - 1) < 0.25 or abs(x - 2) < 0.25) else 10


def _oracle(args, dps):
    x = mp.mpf(args[0])
    if _is_pole(x):
        return mp.inf
    with mp.workdps(dps + _extra(x)):
        return mp.re(mp.loggamma(x))


def _oracle2(args, dps):
    x = mp.mpf(args[0])
    if _is_pole(x):
        return mp.inf
    with mp.workdps(dps + 2 * _extra(x)):
        if x > 0 and x < 1e5:
            return mp.log(mp.gamma(x))
        if x < 0 and x > -1e5:
            return mp.log(abs(mp.gamma(x)))
        # Large |x|: explicit Stirling series (+ reflection), no mpmath loggamma.
        if x < 0:
            return (mp.log(mp.pi) - mp.log(abs(mp.sinpi(x))) - mp.log(-x)
                    - _stirling(-x))
        return _stirling(x)


def _stirling(x):
    """lgamma(x) for x >= 1e5 by the Stirling series (terms shrink past 2^-400)."""
    r = (x - mp.mpf(0.5)) * mp.log(x) - x + mp.log(2 * mp.pi) / 2
    xp, x2 = x, x * x
    eps = mp.mpf(2) ** -(mp.mp.prec + 10)
    for k in range(1, 60):
        term = mp.bernoulli(2 * k) / (2 * k * (2 * k - 1) * xp)
        r += term
        if abs(term) < eps * abs(r):
            return r
        xp *= x2
    raise RuntimeError("Stirling series did not converge")


def _bound(args, v):
    return 0 if args[0] > 2.5 else 1


TIE = 2.0 ** -15  # the x >= 8 cell's near-tie exception (corvus docs, #31)
TIE_TABLE = 2.0 ** -12  # the table's: 2^-65 relative is <= 2^-12 ulp


def _tie(args):
    return TIE_TABLE if args[0] < 8 else TIE


_ABS = mp.mpf(2) ** -53


def _abs_tol(args, v):
    if args[0] < 0 and mp.isfinite(v) and abs(v) < 1:
        return _ABS
    return None


def _st(name, gen, n, seam=None, note=""):
    return Stratum("lgamma", name, gen, n, _oracle, _bound, oracle2=_oracle2,
                   abs_tol=_abs_tol, tie_exempt=_tie, seam=seam, note=note)


def _near_poles(r, n):
    out = []
    for _ in range(n):
        k = r.randint(1, 30)
        d = log_uniform(r, 2.0 ** -52 * k, 1e-2)
        out.append((-k + (d if r.random() < 0.5 else -d),))
    return out


def _half_grid(r, n):
    return [(-(r.randint(2 ** 52, 2 ** 53 - 1) + 0.5 * r.randint(0, 1)),) for _ in range(n)]


def _table_edges(r, n):
    out = []
    for i in range(176):
        edge = 2.5 + (i + 0.5) / 32
        out.extend((x,) for x in neighbours(edge, 2) if 2.5 < x < 8.0)
    return out


def _table_edge_mantissa(r, n):
    out = []
    for _ in range(n):
        edge = 2.5 + (r.randrange(176) + 0.5) / 32
        out.append((full_mantissa(r, edge - 2.0 ** -12, edge + 2.0 ** -12),))
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
        _st("pos_small", lambda r, n: [(log_uniform(r, 2.0 ** -1022, 0.5),) for _ in range(n)], 2000),
        _st("pos_zone", lambda r, n: [(full_mantissa(r, 0.5, 2.5),) for _ in range(n)], 3000),
        _st("near_zero_1", lambda r, n: [(1 + (1 if r.random() < 0.5 else -1) * log_uniform(r, 2.0 ** -52, 0.1),)
                                         for _ in range(n)], 1500, note="relative bound at the zero"),
        _st("near_zero_2", lambda r, n: [(2 + (1 if r.random() < 0.5 else -1) * log_uniform(r, 2.0 ** -51, 0.1),)
                                         for _ in range(n)], 1500, note="relative bound at the zero"),
        _st("pos_table", lambda r, n: [(full_mantissa(r, 2.5, 8.0),) for _ in range(n)], 3000),
        _st("stirling_tie_witness",
            lambda r, n: [(float.fromhex("0x1.9bdb8fed86e38p+3"),)], 1,
            note="1 ULP, 1.8e-5 ulp from a midpoint: the 0-ULP cell's near-tie"),
        _st("table_tie_witness",
            lambda r, n: [(float.fromhex("0x1.a1f3dc0403336p+1"),)], 1,
            note="1 ULP, 3.7e-5 ulp from a midpoint: the widest of 6M inputs"),
        _st("table_edges", _table_edges, 1408, note="all 176 cell edges, +-2 ulps"),
        _st("table_edge_mantissa", _table_edge_mantissa, 2000,
            note="full-mantissa rows within 2^-12 of a random cell edge"),
        _st("pos_stirling", lambda r, n: [(log_uniform(r, 8.0, _MAXARG),) for _ in range(n)], 3000,
            note="CR cell"),
        _st("pos_stirling_low", lambda r, n: [(full_mantissa(r, 8.0, 64.0),) for _ in range(n)], 2000,
            note="CR cell, tightest end"),
        _st("subnormal", lambda r, n: [((subnormal(r) if r.random() < 0.5 else -subnormal(r)),) for _ in range(n)],
            1500),
        _st("neg_unit_intervals", lambda r, n: [(-r.randint(0, 19) - full_mantissa(r, 0.0, 1.0),) for _ in range(n)],
            6000, note="between poles, -20 < x < 0"),
        _st("neg_near_poles", _near_poles, 2000),
        _st("neg_large", lambda r, n: [(-log_uniform(r, 20.0, 2.0 ** 52),) for _ in range(n)], 2000),
        _st("neg_half_grid_rand", _half_grid, 400, note="2^52 <= |x| < 2^53"),
        _st("overflow_over", lambda r, n: [(log_uniform(r, _MAXARG, 1.7e308),) for _ in range(n)], 200,
            note="+inf"),
    ]
    return out


# ---- seam-audit fields (#82 T2; tools/hunt/seam_audit.py docstring) ----------
for _s in SEAMS:
    if _s["name"] == "max_arg":
        # Above it the result is +inf by construction, an exact-match special
        # (route lgamma/special/overflow, exercised by the census grid).
        _s["sides"] = "below"
