"""Elementary-family hunt strata: exp, expm1, log, log1p, log1pmx, cos,
sin (#82; expm1 and log1pmx #45).

Bounds (docs/ACCURACY.md "exp, log, log1p" and "cos and sin"):
  exp    normal results 0 ULP (CR) except within ~2^-15 ulp of a rounding tie
         (include/corvus/corvus.h); subnormal results 1 ULP
  log    0 ULP (same near-tie statement)
  log1p  0 ULP (same; ~2^-23 ulp below |x| = 2^-30, we exempt 2^-15 only
         outside the tiny band)
  expm1  0 ULP, except within ~2^-11 ulp of a tie above |x| = 2^-6 (the
         exponential branch amplifies the core up to ~65x) and ~2^-15 below
  log1pmx 0 ULP (same near-tie statement as log)
  cos    1 ULP over the full double range
  sin    1 ULP over the full double range
Seam constants come from src/exp-inl.h (|x| <= 708 fast path, ScaleTwo
thresholds), src/exp_dd_data.h (N = 128, kExpXMax = 1100), src/log-inl.h
(fast path on normals, log1p cut 2^-30), src/log_dd_data.h (128 slots,
centring at slot 53), src/trig_data.h (kTrigDMax = 2^23), src/trig_ph_data.h
(Payne-Hanek rows e = 23..1023).
"""
import math
import random
from functools import lru_cache

import mpmath as mp
from mpmath import mpf

from core import (Stratum, bits, from_bits, full_mantissa, neighbours,
                  subnormal)

DBL_MAX = 1.7976931348623157e308
MIN_NORMAL = 2.2250738585072014e-308
TIE = 2.0 ** -15          # documented near-tie exception, in ulps
TIE_LOG1P_TINY = 2.0 ** -23
TIE_EXPM1_EXP = 2.0 ** -11  # expm1 above the 2^-6 cut (corvus.h)


def _f(v):
    with mp.workdps(60):
        return float(v)


with mp.workdps(60):
    _EXP_OVER = float(mp.log(mpf(2) ** 1024 * (1 - mpf(2) ** -54)))  # exp(x) -> inf above
    _EXP_MINN = float(-1022 * mp.log(2))                              # exp(x) = 2^-1022
    _EXP_ZERO = float(-1075 * mp.log(2))                              # exp(x) < 2^-1075 -> 0
    _LN2_256 = float(mp.log(2) / 256)                                 # k rounding 0 -> 1
    _PI4 = float(mp.pi / 4)
    _PI2 = float(mp.pi / 2)
    _PI = float(mp.pi)
    _3PI2 = float(3 * mp.pi / 2)
    _2PI = float(2 * mp.pi)
_WORST_E849 = math.ldexp(6381956970095103, 797)  # |r| = 2^-60.89 (ACCURACY.md)

# ---- seam registry ----------------------------------------------------------
SEAMS = [
    # exp
    {"fn": "exp", "name": "fast_hi", "coord": "x", "value": 708.0, "source": "ExpFast domain |x| <= 708 (exp-inl.h)"},
    {"fn": "exp", "name": "fast_lo", "coord": "x", "value": -708.0, "source": "ExpFast domain |x| <= 708 (exp-inl.h)"},
    {"fn": "exp", "name": "overflow", "coord": "x", "value": _EXP_OVER, "source": "exp(x) = 2^1024(1-2^-54): RN to +inf above"},
    {"fn": "exp", "name": "min_normal_out", "coord": "x", "value": _EXP_MINN, "source": "exp(x) = 2^-1022 (subnormal-output edge)"},
    {"fn": "exp", "name": "zero_out", "coord": "x", "value": _EXP_ZERO, "source": "exp(x) = 2^-1075 (underflow to 0)"},
    {"fn": "exp", "name": "one_up", "coord": "x", "value": 2.0 ** -53, "source": "exp(x) = 1 + 2^-53 + 2^-107: tie + tiny"},
    {"fn": "exp", "name": "one_down", "coord": "x", "value": -(2.0 ** -54), "source": "exp(x) rounds to 1 below 1-2^-54"},
    {"fn": "exp", "name": "reduce_k1", "coord": "x", "value": _LN2_256, "source": "kExpInvL rounding of k, ln2/(2N) (exp_dd_data.h)"},
    {"fn": "exp", "name": "clamp_hi", "coord": "x", "value": 1100.0, "source": "kExpXMax (exp_dd_data.h)"},
    {"fn": "exp", "name": "clamp_lo", "coord": "x", "value": -1100.0, "source": "kExpXMax (exp_dd_data.h)"},
    # log
    {"fn": "log", "name": "one", "coord": "x", "value": 1.0, "source": "slots 0/127 R=1, L=0 (log_dd-inl.h)"},
    {"fn": "log", "name": "two", "coord": "x", "value": 2.0, "source": "binade edge, k changes (log_dd-inl.h)"},
    {"fn": "log", "name": "split", "coord": "x", "value": 1.0 + 53.0 / 128.0, "source": "kLogSplit = 1 + 53/128 (log_dd_data.h)"},
    {"fn": "log", "name": "slot1", "coord": "x", "value": 1.0 + 1.0 / 128.0, "source": "first slot boundary (kLogN = 128)"},
    {"fn": "log", "name": "min_normal", "coord": "x", "value": MIN_NORMAL, "source": "LogFast domain x >= 2^-1022; 2^600 prescale below"},
    {"fn": "log", "name": "dbl_max", "coord": "x", "value": DBL_MAX, "source": "LogFast domain x <= DBL_MAX"},
    # log1p
    {"fn": "log1p", "name": "tiny_hi", "coord": "x", "value": 2.0 ** -30, "source": "tiny-series cut |x| <= 2^-30 (log-inl.h)"},
    {"fn": "log1p", "name": "tiny_lo", "coord": "x", "value": -(2.0 ** -30), "source": "tiny-series cut |x| <= 2^-30 (log-inl.h)"},
    {"fn": "log1p", "name": "minus_one", "coord": "x", "value": -1.0, "source": "domain edge, 1+x = 0"},
    {"fn": "log1p", "name": "minus_half", "coord": "x", "value": -0.5, "source": "1+x = 1/2, Sterbenz edge of TwoSum"},
    {"fn": "log1p", "name": "split_m1", "coord": "x", "value": 53.0 / 128.0, "source": "1+x = kLogSplit"},
    {"fn": "log1p", "name": "plus_one", "coord": "x", "value": 1.0, "source": "1+x = 2, binade edge"},
    {"fn": "log1p", "name": "ulp_half", "coord": "x", "value": 2.0 ** -53, "source": "1+x rounding edge; TwoSum lo only"},
    {"fn": "log1p", "name": "dbl_max", "coord": "x", "value": DBL_MAX, "source": "largest finite input"},
    # expm1 (src/exp-inl.h Expm1Vec)
    {"fn": "expm1", "name": "cut_hi", "coord": "x", "value": 2.0 ** -6, "source": "kExpm1Cut series/exponential (dd_special-inl.h)"},
    {"fn": "expm1", "name": "cut_lo", "coord": "x", "value": -(2.0 ** -6), "source": "kExpm1Cut series/exponential (dd_special-inl.h)"},
    {"fn": "expm1", "name": "minus38", "coord": "x", "value": -38.0, "source": "x <= -38 returns -1 exactly (exp-inl.h)"},
    {"fn": "expm1", "name": "overflow", "coord": "x", "value": _EXP_OVER, "source": "e^x - 1 = 2^1024(1-2^-54): RN to +inf above"},
    # log1pmx (src/log-inl.h Log1pmxVec)
    {"fn": "log1pmx", "name": "tiny_hi", "coord": "x", "value": 2.0 ** -480, "source": "tiny/series cut (log-inl.h)"},
    {"fn": "log1pmx", "name": "tiny_lo", "coord": "x", "value": -(2.0 ** -480), "source": "tiny/series cut (log-inl.h)"},
    {"fn": "log1pmx", "name": "subnormal_out", "coord": "x", "value": 2.0 ** -511, "source": "x^2/2 = 2^-1023: normal/subnormal result"},
    {"fn": "log1pmx", "name": "cut_hi", "coord": "x", "value": 1.0 / 16, "source": "kPhiCut series/direct (dd_special_data.h)"},
    {"fn": "log1pmx", "name": "cut_lo", "coord": "x", "value": -1.0 / 16, "source": "kPhiCut series/direct (dd_special_data.h)"},
    {"fn": "log1pmx", "name": "minus_half", "coord": "x", "value": -0.5, "source": "1+x = 1/2, Sterbenz edge of TwoSum"},
    {"fn": "log1pmx", "name": "minus_one", "coord": "x", "value": -1.0, "source": "domain edge, 1+x = 0"},
    {"fn": "log1pmx", "name": "dbl_max", "coord": "x", "value": DBL_MAX, "source": "largest finite input"},
    # trig (cos and sin share the seam values)
] + [
    {"fn": fn, "name": nm, "coord": "x", "value": val, "source": src}
    for fn in ("cos", "sin")
    for nm, val, src in [
        ("small_tiny", 2.0 ** -26.5, "1 - x^2/2 / x - x^3/6 rounding edge (cores)"),
        ("pi4", _PI4, "Round(x*2/pi) switches 0 -> 1 (kTrigTwoOverPi)"),
        ("pio2", _PI2, "cos zero / sin = 1"),
        ("pi", _PI, "sin zero / cos = -1"),
        ("3pio2", _3PI2, "cos zero"),
        ("2pi", _2PI, "sin zero"),
        ("ph_threshold", 2.0 ** 23, "kTrigDMax: exact split -> Payne-Hanek"),
        ("worst_e849", _WORST_E849, "global worst cancellation |r| = 2^-60.89"),
        ("dbl_max", DBL_MAX, "largest finite, last PH row e = 1023"),
    ]
]

_DOMAIN = {
    "exp": lambda x: math.isfinite(x),
    "log": lambda x: math.isfinite(x) and x > 0,
    "log1p": lambda x: math.isfinite(x) and x > -1,
    "expm1": lambda x: math.isfinite(x),
    "log1pmx": lambda x: math.isfinite(x) and x > -1,
    "cos": lambda x: math.isfinite(x),
    "sin": lambda x: math.isfinite(x),
}


# ---- oracles ----------------------------------------------------------------
def _exp_o(args, dps):
    x = args[0]
    if x > 800:
        return mp.ldexp(mpf(1), 5000)
    if x < -800:
        return mp.ldexp(mpf(1), -5000)
    with mp.workdps(dps + 10):
        return mp.exp(mpf(x))


def _exp_o2(args, dps):
    x = args[0]
    if x > 800 or x < -800:
        return _exp_o(args, dps)
    with mp.workdps(dps + 40):
        h = mp.exp(mpf(x) / 2)
        return h * h


def _log_o(args, dps):
    x = args[0]
    if x == 0:
        return mp.ninf
    if x < 0:
        return mp.nan
    with mp.workdps(dps + 10):
        m = mpf(x)
        if abs(m - 1) < 0.5:
            return mp.log1p(m - 1)  # m - 1 exact
        return mp.log(m)


def _log_o2(args, dps):
    x = args[0]
    if x <= 0:
        return _log_o(args, dps)
    with mp.workdps(dps + 60):
        m = mpf(x)
        return mp.log(m)


def _log1p_o(args, dps):
    x = args[0]
    if x < -1:
        return mp.nan
    if x == -1:
        return mp.ninf
    if x == 0:
        return mpf(x)
    with mp.workdps(dps + 10):
        return mp.log1p(mpf(x))


def _log1p_o2(args, dps):
    x = args[0]
    if x <= -1 or x == 0:
        return _log1p_o(args, dps)
    extra = 60 + (int(-2 * math.log10(abs(x))) if abs(x) < 1 else 0)
    with mp.workdps(dps + extra):
        return mp.log(1 + mpf(x))


def _expm1_o(args, dps):
    x = args[0]
    if x > 800:
        return mp.ldexp(mpf(1), 5000)
    if x == 0:
        return mpf(x)
    with mp.workdps(dps + 10):
        return mp.expm1(mpf(x))


def _expm1_o2(args, dps):
    """Independent form: exp(x) - 1 at enough extra digits to absorb the
    cancellation (log10(1/|x|) digits) for small x; 2*sinh(x/2)*e^(x/2)
    is not used, so the two oracles share no mpmath routine."""
    x = args[0]
    if x > 800 or x == 0:
        return _expm1_o(args, dps)
    extra = 60 + (int(-math.log10(abs(x))) if abs(x) < 1 else 0)
    with mp.workdps(dps + extra):
        return mp.exp(mpf(x)) - 1


def _log1pmx_extra(x):
    return 10 + (2 * int(math.ceil(-math.log10(abs(x)))) if abs(x) < 0.5 else 0)


def _log1pmx_o(args, dps):
    x = args[0]
    if x < -1:
        return mp.nan
    if x == -1:
        return mp.ninf
    if x == 0:
        return mpf(0)
    with mp.workdps(dps + _log1pmx_extra(x)):
        m = mpf(x)
        return mp.log1p(m) - m


def _log1pmx_o2(args, dps):
    """Independent form: the alternating series -x^2/2 + x^3/3 - ... for
    |x| < 1/4 (no cancellation), log(1 + x) - x at extra digits above."""
    x = args[0]
    if x <= -1 or x == 0:
        return _log1pmx_o(args, dps)
    if abs(x) < 0.25:
        with mp.workdps(dps + 20):
            m = mpf(x)
            tot, term, k = mpf(0), m * m, 2
            eps = mpf(10) ** -(dps + 15)
            while True:
                t = term / k * (1 if k % 2 else -1)
                tot += t
                if abs(t) < eps * abs(tot):
                    break
                term *= m
                k += 1
            return tot
    with mp.workdps(dps + 40 + _log1pmx_extra(x)):
        m = mpf(x)
        return mp.log(1 + m) - m


def _trig_extra(x):
    e = math.frexp(x)[1] if x else 0
    return 50 + (int(e * 0.30103) if e > 0 else 0)


def _cos_o(args, dps):
    with mp.workdps(dps + _trig_extra(args[0])):
        return mp.cos(mpf(args[0]))


def _sin_o(args, dps):
    with mp.workdps(dps + _trig_extra(args[0])):
        return mp.sin(mpf(args[0]))


def _cos_o2(args, dps):
    with mp.workdps(dps + _trig_extra(args[0]) + 40):
        return mp.cos(mpf(args[0]))


def _sin_o2(args, dps):
    with mp.workdps(dps + _trig_extra(args[0]) + 40):
        return mp.sin(mpf(args[0]))


ORACLES = {"exp": (_exp_o, _exp_o2), "log": (_log_o, _log_o2),
           "log1p": (_log1p_o, _log1p_o2), "cos": (_cos_o, _cos_o2),
           "sin": (_sin_o, _sin_o2), "expm1": (_expm1_o, _expm1_o2),
           "log1pmx": (_log1pmx_o, _log1pmx_o2)}

_MIN_N = mpf(2) ** -1022


def _bound_exp(args, v):
    return 1 if abs(v) < _MIN_N else 0


def _bound0(args, v):
    return 0


def _bound1(args, v):
    return 1


BOUNDS = {"exp": _bound_exp, "log": _bound0, "log1p": _bound0,
          "cos": _bound1, "sin": _bound1, "expm1": _bound0,
          "log1pmx": _bound0}


# ---- generators -------------------------------------------------------------
def _rand_exp(rng, elo, ehi, sign=False):
    """Double with unbiased exponent uniform in [elo, ehi], random mantissa."""
    e = rng.randint(elo, ehi)
    x = from_bits(((e + 1023) << 52) | rng.getrandbits(52))
    return -x if sign and rng.random() < 0.5 else x


def _band(fn, t, side=0):
    """Full-mantissa offsets around seam t at log-uniform distances
    |t| * 2^-u, u in [3, 52]; points outside the function's domain are
    mirrored back."""
    dom = _DOMAIN[fn]

    def gen(rng, n):
        out = []
        while len(out) < n:
            u = rng.uniform(3, 52)
            off = abs(t) * 2.0 ** -u * rng.uniform(1, 2)
            s = rng.choice((-1, 1)) if side == 0 else side
            for sg in (s, -s):
                x = t + sg * off
                if dom(x):
                    out.append((x,))
                    break
        return out
    return gen


def _seam_strata(fn):
    out = []
    o, o2 = ORACLES[fn]
    for s in SEAMS:
        if s["fn"] != fn:
            continue
        t, nm = s["value"], s["name"]
        out.append(Stratum(fn, f"seam:{nm}:bits",
                           lambda r, n, t=t: [(x,) for x in neighbours(t, 64)],
                           129, o, BOUNDS[fn], o2, tie_exempt=TIE, seam=nm))
        out.append(Stratum(fn, f"seam:{nm}:band", _band(fn, t), 150,
                           o, BOUNDS[fn], o2, tie_exempt=TIE, seam=nm))
    return out


def _tie_dist(v):
    """Distance (ulps) of a normal-range mpf to the nearest rounding midpoint,
    computed on the mpf's integer mantissa (no logs)."""
    sign, man, exp, bc = v._mpf_
    sh = bc - 53
    if sh <= 0:
        return 0.5
    rem = man & ((1 << sh) - 1)
    return abs(rem / (1 << sh) - 0.5)


def _near_tie(fn, sample, cap_factor=80):
    """Dense random candidates screened at 28 dps; keeps rows within 0.01 ulp
    of a midpoint (about 2% of candidates), bounded at cap_factor*n tries."""
    o = ORACLES[fn][0]

    def gen(rng, n):
        out, tries = [], 0
        with mp.workdps(28):
            while len(out) < n and tries < cap_factor * n:
                tries += 1
                x = sample(rng)
                if not _DOMAIN[fn](x):
                    continue
                v = +o((x,), 15)
                if not v or not mp.isfinite(v) or abs(v) < _MIN_N:
                    continue
                if _tie_dist(v) < 0.01:
                    out.append((x,))
        return out
    return gen


# --- trig: nearest doubles to multiples of pi/4, pi/2 -------------------------
def _multiples(rng, n_centres, k_hi_log2, step_num, step_den, offset, nb):
    """Doubles around (k*step + offset)*pi: k small and log-uniform."""
    out = []
    with mp.workdps(80):
        kmax = 2 ** k_hi_log2
        ks = list(range(0, min(32, n_centres // 4)))
        while len(ks) < n_centres:
            ks.append(int(2 ** rng.uniform(5, k_hi_log2 - 0.01)))
        for k in ks:
            v = (mpf(k) * step_num + offset) * mp.pi / step_den
            x = float(v)
            if x > 0:
                out.extend((y,) for y in neighbours(x, nb) if y > 0)
    return out


def _gen_near_pi4(rng, n):  # x ~ (2k+1) * pi/4
    return _multiples(rng, max(1, n // 7), 22, 2, 4, 1, 3)


def _gen_near_pio2(rng, n):  # x ~ k * pi/2
    return _multiples(rng, max(1, n // 5), 22, 1, 2, 0, 2)


@lru_cache(maxsize=1)
def _cf_worst():
    """Per exponent e = 23..1023, mantissas M in [2^52, 2^53) minimising the
    distance of M*2^(e-52)*(2/pi) to an integer (mixed parity) and to an even
    integer, from continued-fraction convergents of the fractional part and
    their smallest in-range multiples."""
    rows = []
    prec_bits = 1400
    with mp.workprec(prec_bits):
        twoopi = 2 / mp.pi
        for e in range(23, 1024):
            for divisor in (1, 2):
                beta = twoopi * mpf(2) ** (e - 52) / divisor
                P = e + 250
                frac = beta - mp.floor(beta)
                num = int(mp.floor(frac * mpf(2) ** P))
                den = 1 << P
                # CF of num/den
                p0, q0, p1, q1 = 0, 1, 1, 0
                a, b = num, den
                best = None
                while b:
                    t = a // b
                    a, b = b, a - t * b
                    p0, q0, p1, q1 = p1, q1, t * p1 + p0, t * q1 + q0
                    q = q1
                    if q == 0:
                        continue
                    if q >= (1 << 53):
                        break
                    err = abs(q * num - p1 * den)  # |q*alpha - p| * den
                    m = -(-(1 << 52) // q)
                    if m * q >= (1 << 53):
                        continue
                    cost = m * err
                    if best is None or cost < best[0]:
                        best = (cost, m * q)
                if best:
                    rows.append(math.ldexp(best[1], e - 52))
    return tuple(rows)


def _gen_huge_cf(rng, n):
    out = [(x,) for x in _cf_worst()]
    out += [(y,) for y in neighbours(_WORST_E849, 3)]
    return out


def _gen_trig_small(rng, n):
    return [(_rand_exp(rng, -27, 22, True),) for _ in range(n)]


def _gen_trig_huge(rng, n):
    return [(_rand_exp(rng, 23, 1023, True),) for _ in range(n)]


def _gen_trig_tiny(rng, n):  # below the polynomial's visible x^2/x^3 terms
    return [(_rand_exp(rng, -60, -27, True),) for _ in range(n)]


def _gen_subnormal_in(rng, n):
    out = [(subnormal(rng) * rng.choice((-1, 1)),) for _ in range(n)]
    out += [(from_bits(1),), (from_bits(0x000FFFFFFFFFFFFF),), (-from_bits(1),),
            (from_bits(0x0008000000000000),), (-from_bits(0x000FFFFFFFFFFFFF),)]
    return out


# ---- strata -----------------------------------------------------------------
def _stratum(fn, name, gen, n, **kw):
    o, o2 = ORACLES[fn]
    kw.setdefault("tie_exempt", TIE)
    return Stratum(fn, name, gen, n, o, BOUNDS[fn], o2, **kw)


def _exp_strata():
    S = [_stratum("exp", "normal_uniform",
                  lambda r, n: [(full_mantissa(r, _EXP_MINN, _EXP_OVER),) for _ in range(n)], 3000),
         _stratum("exp", "small_logu",
                  lambda r, n: [(_rand_exp(r, -60, -1, True),) for _ in range(n)], 2000),
         _stratum("exp", "near_zero_unit",
                  lambda r, n: [(full_mantissa(r, -1.0, 1.0),) for _ in range(n)], 1000),
         _stratum("exp", "subnormal_out_band",
                  lambda r, n: [(full_mantissa(r, -745.14, _EXP_MINN),) for _ in range(n)], 3000,
                  note="x in [-745.14, -708.4]"),
         _stratum("exp", "subnormal_out_top",
                  lambda r, n: [(full_mantissa(r, _EXP_MINN - 0.7, _EXP_MINN),) for _ in range(n)], 1000,
                  note="top subnormal binade, 0.75-ulp construction bound"),
         _stratum("exp", "near_tie", _near_tie("exp", lambda r: (
             full_mantissa(r, -700, 700) if r.random() < 0.5 else full_mantissa(r, -2, 2))), 800,
             note="screened within 0.01 ulp of a midpoint"),
         _stratum("exp", "pow2_inputs",
                  lambda r, n: [(s * 2.0 ** -k,) for k in r.sample(range(0, 80), min(80, max(1, n))) for s in (-1, 1)],
                  120),
         _stratum("exp", "subnormal_in", _gen_subnormal_in, 300),
         _stratum("exp", "extremes",
                  lambda r, n: [(v,) for v in (0.0, DBL_MAX, -DBL_MAX, 1e300, -1e300, 1e10, -1e10,
                                               1e5, -1e5, 800.5, -800.5, 746.0, -746.0, 709.9, 710.0)],
                  15)]
    # reduction-grid midpoints (k + 1/2) * ln2/128, where Round(x*N/ln2) flips
    def grid(r, n):
        out = []
        with mp.workdps(60):
            L = mp.log(2) / 128
            for _ in range(max(1, n // 5)):
                k = r.randint(-5400, 5400)
                out.extend((y,) for y in neighbours(float((k + mpf(1) / 2) * L), 2))
        return out
    S.append(_stratum("exp", "reduction_midpoints", grid, 1000))
    return S + _seam_strata("exp")


def _log_strata():
    def near_one(r, n):
        out = []
        for _ in range(n):
            u = r.uniform(1, 53)
            d = 2.0 ** -u * r.uniform(1, 2)
            out.append((1.0 + d if r.random() < 0.5 else 1.0 - d / 2,))
        return out

    def slots(r, n):
        out = []
        for _ in range(max(1, n // 9)):
            j = r.randint(0, 128)
            k = r.randint(-30, 30) if r.random() < 0.7 else r.randint(-1021, 1022)
            out.extend((math.ldexp(x, k),) for x in neighbours(1.0 + j / 128.0, 4))
        return out

    def pow2(r, n):
        out = []
        for k in r.sample(range(-1074, 1024), min(n // 3, 2098)):
            x = math.ldexp(1.0, k)
            out.extend((y,) for y in neighbours(x, 1) if y > 0)
        return out

    S = [_stratum("log", "general_exp_uniform",
                  lambda r, n: [(_rand_exp(r, -1022, 1023),) for _ in range(n)], 3000),
         _stratum("log", "near_one", near_one, 3000),
         _stratum("log", "subnormal_x",
                  lambda r, n: [(subnormal(r),) for _ in range(n // 2)]
                  + [(subnormal(r, 6),) for _ in range(n - n // 2)], 1500),
         _stratum("log", "slot_boundaries", slots, 1000,
                  note="1+j/128 times 2^k, 4 bit-neighbours each side"),
         _stratum("log", "pow2_neighbours", pow2, 900),
         _stratum("log", "unit_range",
                  lambda r, n: [(full_mantissa(r, 0.5, 4.0),) for _ in range(n)], 1000),
         _stratum("log", "near_tie", _near_tie("log", lambda r: (
             _rand_exp(r, -1022, 1023) if r.random() < 0.3 else full_mantissa(r, 0.1, 10.0))), 800,
             note="screened within 0.01 ulp of a midpoint")]
    return S + _seam_strata("log")


def _log1p_strata():
    def general(r, n):
        out = []
        for i in range(n):
            if i % 3 == 0:
                out.append((_rand_exp(r, -60, 1023),))
            elif i % 3 == 1:
                out.append((-full_mantissa(r, 0.0, 1.0),))
            else:
                out.append((-_rand_exp(r, -60, -1),))
        return [a for a in out if a[0] > -1]

    def near_minus_one(r, n):
        return [(-1.0 + 2.0 ** -r.uniform(1, 53) * r.uniform(1, 2),) for _ in range(n)]

    def tiny(r, n):  # [2^-60, 2^-26], both signs
        return [(_rand_exp(r, -60, -27, True),) for _ in range(n)]

    S = [_stratum("log1p", "general", general, 3000),
         _stratum("log1p", "tiny_band", tiny, 2500, tie_exempt=TIE_LOG1P_TINY,
                  note="|x| in [2^-60, 2^-26]: the fixed v0.8 defect band"),
         _stratum("log1p", "tiny_band_2pow53", lambda r, n: [(_rand_exp(r, -55, -45, True),) for _ in range(n)],
                  1500, tie_exempt=TIE_LOG1P_TINY, note="defect core [2^-53, 2^-47]"),
         _stratum("log1p", "near_minus_one", near_minus_one, 1500),
         _stratum("log1p", "unit_range",
                  lambda r, n: [(full_mantissa(r, -0.99, 3.0),) for _ in range(n)], 1000),
         _stratum("log1p", "near_tie", _near_tie("log1p", lambda r: (
             _rand_exp(r, -29, 1023)
             if r.random() < 0.5 else -full_mantissa(r, 0.0, 0.999))), 800,
             note="screened within 0.01 ulp of a midpoint"),
         _stratum("log1p", "subnormal_in", _gen_subnormal_in, 300,
                  tie_exempt=TIE_LOG1P_TINY)]
    return S + _seam_strata("log1p")


def _expm1_strata():
    def cut_band(r, n):  # [2^-7, 2^-4] both signs: the amplified-error zone
        return [(full_mantissa(r, 2.0 ** -7, 2.0 ** -4) * r.choice((-1, 1)),)
                for _ in range(n)]

    S = [_stratum("expm1", "series_logu",
                  lambda r, n: [(_rand_exp(r, -60, -7, True),) for _ in range(n)], 2000),
         _stratum("expm1", "cut_band", cut_band, 2000, tie_exempt=TIE_EXPM1_EXP,
                  note="|x| in [2^-7, 2^-4]: exponential branch, ~65x amplification"),
         _stratum("expm1", "exp_branch",
                  lambda r, n: [(full_mantissa(r, -38.0, _EXP_OVER),) for _ in range(n)], 2000,
                  tie_exempt=TIE_EXPM1_EXP),
         _stratum("expm1", "neg_to_minus_one",
                  lambda r, n: [(full_mantissa(r, -45.0, -30.0),) for _ in range(n)], 500,
                  tie_exempt=TIE_EXPM1_EXP, note="where expm1 meets -1"),
         _stratum("expm1", "near_tie", _near_tie("expm1", lambda r: (
             full_mantissa(r, -30, 700) if r.random() < 0.5
             else _rand_exp(r, -40, -7, True))), 800,
             tie_exempt=TIE_EXPM1_EXP, note="screened within 0.01 ulp of a midpoint"),
         _stratum("expm1", "subnormal_in", _gen_subnormal_in, 300)]
    out = S
    for st in _seam_strata("expm1"):
        if st.seam in ("cut_hi", "cut_lo", "minus38", "overflow"):
            st.tie_exempt = TIE_EXPM1_EXP
        out.append(st)
    return out


def _log1pmx_strata():
    def general(r, n):
        out = []
        for i in range(n):
            if i % 3 == 0:
                out.append((_rand_exp(r, -4, 1023),))
            elif i % 3 == 1:
                out.append((-full_mantissa(r, 1.0 / 16, 1.0),))
            else:
                out.append((-_rand_exp(r, -60, -5),))
        return [a for a in out if a[0] > -1]

    def near_minus_one(r, n):
        return [(-1.0 + 2.0 ** -r.uniform(1, 53) * r.uniform(1, 2),) for _ in range(n)]

    S = [_stratum("log1pmx", "general", general, 3000),
         _stratum("log1pmx", "series_band",
                  lambda r, n: [(_rand_exp(r, -480, -5, True),) for _ in range(n)], 2500,
                  note="|x| in [2^-480, 2^-4]: the series branch"),
         _stratum("log1pmx", "cut_band",
                  lambda r, n: [(full_mantissa(r, 1.0 / 32, 1.0 / 8) * r.choice((-1, 1)),)
                                for _ in range(n)], 2000,
                  note="|x| in [1/32, 1/8]: both sides of kPhiCut"),
         _stratum("log1pmx", "subnormal_out",
                  lambda r, n: [(_rand_exp(r, -540, -505, True),) for _ in range(n)], 1500,
                  note="results subnormal or rounding to zero"),
         _stratum("log1pmx", "near_minus_one", near_minus_one, 1500),
         _stratum("log1pmx", "near_tie", _near_tie("log1pmx", lambda r: (
             _rand_exp(r, -4, 1023) if r.random() < 0.5
             else -full_mantissa(r, 1.0 / 16, 0.999))), 800,
             note="screened within 0.01 ulp of a midpoint"),
         _stratum("log1pmx", "subnormal_in", _gen_subnormal_in, 300)]
    return S + _seam_strata("log1pmx")


def _trig_strata(fn):
    # #79: small_random is the generic population above the tiny band; 2000
    # reference rows, not 120 (83% of the base rows sit within 2^-53 of 0 or
    # +-1, where any kernel is correctly rounded).
    S = [_stratum(fn, "small_random", _gen_trig_small, 3000, tie_exempt=0.0,
                  note="full-mantissa, exponent uniform in [-27, 22]", ref_cap=2000),
         _stratum(fn, "huge_random", _gen_trig_huge, 3000, tie_exempt=0.0,
                  note="full-mantissa, exponent uniform in [23, 1023]"),
         _stratum(fn, "tiny_random", _gen_trig_tiny, 500, tie_exempt=0.0),
         _stratum(fn, "near_pi4_odd", _gen_near_pi4, 500, tie_exempt=0.0),
         _stratum(fn, "near_k_pio2", _gen_near_pio2, 500, tie_exempt=0.0),
         _stratum(fn, "huge_cf_worst", _gen_huge_cf, 2000, tie_exempt=0.0,
                  note="per-exponent worst cancellation (mixed + even n) and e849 global worst"),
         _stratum(fn, "subnormal_in", _gen_subnormal_in, 300, tie_exempt=0.0)]
    for s in _seam_strata(fn):
        s.tie_exempt = 0.0
        S.append(s)
    return S


def strata():
    out = _exp_strata() + _log_strata() + _log1p_strata()
    out += _expm1_strata() + _log1pmx_strata()
    out += _trig_strata("cos") + _trig_strata("sin")
    return out


# ---- seam-audit fields (#82 T2; tools/hunt/seam_audit.py docstring) ----------
_AUDIT = {
    (fn, "dbl_max"): {"sides": "below"}
    for fn in ("log", "log1p", "log1pmx", "cos", "sin")
}
_AUDIT[("log1p", "minus_one")] = {"sides": "above"}  # x < -1 is outside the domain
_AUDIT[("log1pmx", "minus_one")] = {"sides": "above"}
for _s in SEAMS:
    _s.update(_AUDIT.get((_s["fn"], _s["name"]), {}))
