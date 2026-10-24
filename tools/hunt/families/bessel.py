"""Bessel hunt strata: i0, i1, i0e, i1e. Bound: ACCURACY.md Bessel -- max 1 ULP
in every region on every tier. i0/i0e are even, i1/i1e odd: strata run on
both signs. Kernel: series for |x| <= 8 (kBesselSplit), Chebyshev tail
f(1/x)/sqrt(2*pi*x) above (src/bessel-inl.h, src/bessel_data.h)."""
import math

import mpmath as mp

from core import Stratum, bits, from_bits, full_mantissa, log_uniform, neighbours, subnormal

FNS = ("i0", "i1", "i0e", "i1e")
DBL_MAX = 1.7976931348623157e308
I0_OVF = float.fromhex("0x1.64fe5304e83e4p+9")
I1_OVF = float.fromhex("0x1.64fe69ff9fec7p+9")
OVF_MAG = mp.mpf(2) ** 1024 * (1 - mp.mpf(2) ** -54)  # round-to-inf threshold

# (name, value, source, applies-to)
_SEAM_DEFS = [
    ("split", 8.0, "kBesselSplit (bessel_data.h): series/tail switch", FNS),
    ("ovf_i0", I0_OVF, "kBesselI0OverflowX (bessel_data.h)", ("i0", "i0e")),
    ("ovf_i1", I1_OVF, "kBesselI1OverflowX (bessel_data.h)", ("i1", "i1e")),
    ("exp_ovf", 709.782712893384, "ln(DBL_MAX): exp(ax) range split in tail",
     FNS),
    ("tail_700", 700.0, "tail sweep upper edge of reference set", FNS),
    ("one_tie", 2.0 ** -25.5, "i0: x^2/4 = 2^-53, 1 + q rounding tie", ("i0", "i0e")),
    ("one_q54", 2.0 ** -26, "i0: x^2/4 = 2^-54, below half-ulp of 1", ("i0", "i0e")),
    ("q_underflow", 2.0 ** -510, "x^2/4 reaches 2^-1022 (q subnormal)", FNS),
    ("min_normal", 2.0 ** -1022, "smallest normal x; x/2 subnormal in i1", FNS),
    ("half_normal", 2.0 ** -1021, "x/2 = 2^-1022: i1 half_ax exactness edge", ("i1", "i1e")),
    ("dbl_max", DBL_MAX, "tail domain top: prescale 2^-32 / postscale 2^16",
     ("i0e", "i1e")),
    ("prescale_split", 2.0 ** 996, "ProdLow Dekker bound 2^996 (no-FMA tiers)", ("i0e", "i1e")),
    ("prescale_pre", 2.0 ** 1020, "prescaled operand x*2^-32 = 2^988 region", ("i0e", "i1e")),
]

SEAMS = [{"fn": fn, "name": f"{n}", "coord": "x", "value": v, "source": src}
         for (n, v, src, fns) in _SEAM_DEFS for fn in fns]


def _series(nu, x):
    """Own power series sum (x/2)^(2k+nu)/(k!(k+nu)!) at the ambient dps."""
    h = x / 2
    q = h * h
    term = h ** nu / mp.factorial(nu)
    s = term
    k = 0
    while True:
        k += 1
        term = term * q / (k * (k + nu))
        s += term
        if abs(term) < abs(s) * mp.mpf(2) ** -(mp.mp.prec + 8) or k > 400:
            break
    if k > 400:
        raise RuntimeError("series cap")
    return s


def _mag(fn, ax, dps, own):
    nu = 1 if fn.startswith("i1") else 0
    scaled = fn.endswith("e")
    ax = mp.mpf(ax)
    if own:
        with mp.workdps(dps + 30):
            v = _series(nu, ax)
            if scaled:
                v = v * mp.exp(-ax)
        return +v
    with mp.workdps(dps + 40):
        v = mp.besseli(nu, ax)
        if scaled:
            v = v * mp.exp(-ax)
    return +v


def _make(fn, own):
    odd = fn.startswith("i1")
    scaled = fn.endswith("e")

    def oracle(args, dps):
        x = mp.mpf(args[0])
        ax = abs(x)
        if own and ax > 30:
            raise RuntimeError("own series out of range")
        v = _mag(fn, ax, dps, own)
        if not scaled and v >= OVF_MAG:
            # core.quantum() cannot take inf: return the finite round-to-inf
            # threshold (round_double maps it to inf; float() gives inf).
            return OVF_MAG if (x > 0 or not odd) else -OVF_MAG
        if odd and x < 0:
            v = -v
        return v
    return oracle


def _bound(args, v):
    return 1


def _stratum(fn, name, gen, n, seam=None, note=""):
    # oracle2 independent only where the own series is usable.
    return Stratum(fn, name, gen, n, _make(fn, False), _bound,
                   oracle2=_make(fn, True), seam=seam, note=note)


def _signed(rng, x):
    return -x if rng.random() < 0.5 else x


def _stratum2(fn, name, gen, n, seam=None, note="", own=False):
    s = Stratum(fn, name, gen, n, _make(fn, False), _bound, seam=seam, note=note)
    if own:
        s.oracle2 = _make(fn, True)
    return s


def strata():
    out = []
    for (sname, t, src, fns) in _SEAM_DEFS:
        for fn in fns:
            own = t <= 30
            top = t >= DBL_MAX / 1.1

            def bits_gen(r, n, t=t):
                pts = neighbours(t, 64)
                return [(x if i % 2 == 0 else -x,) for i, x in enumerate(pts)]

            def lo_gen(r, n, t=t):
                return [(_signed(r, full_mantissa(r, t * 0.97, t)),) for _ in range(n)]

            def hi_gen(r, n, t=t, top=top):
                hi = DBL_MAX if top else t * 1.03
                return [(_signed(r, min(full_mantissa(r, t, hi), DBL_MAX)),) for _ in range(n)]

            out.append(_stratum2(fn, f"seam:{sname}:bits", bits_gen, 129, sname, own=own))
            out.append(_stratum2(fn, f"seam:{sname}:below", lo_gen, 150, sname, own=own))
            if not top:
                out.append(_stratum2(fn, f"seam:{sname}:above", hi_gen, 150, sname, own=own))
    for fn in FNS:
        for lo, hi in ((1, 4), (4, 6), (6, 8)):
            out.append(_stratum2(
                fn, f"series:[{lo},{hi})",
                lambda r, n, lo=lo, hi=hi: [(_signed(r, full_mantissa(r, lo, hi)),) for _ in range(n)],
                2200, note="#75 band" if lo == 6 else "", own=True))
        out.append(_stratum2(
            fn, "series:small",
            lambda r, n: [(_signed(r, log_uniform(r, 1e-8, 1.0)),) for _ in range(n)],
            600, own=True))
        out.append(_stratum2(
            fn, "tiny:log",
            lambda r, n: [(_signed(r, log_uniform(r, 1e-300, 1e-8)),) for _ in range(n)],
            400, own=True))
        out.append(_stratum2(
            fn, "tiny:subnormal",
            lambda r, n: [(_signed(r, subnormal(r)),) for _ in range(n)], 400, own=True,
            note="subnormal x; i1/i1e subnormal outputs"))
        out.append(_stratum2(
            fn, "tiny:subnormal_top",
            lambda r, n: [(_signed(r, subnormal(r, 3)),) for _ in range(n)], 200, own=True))
        out.append(_stratum2(
            fn, "tail:(8,700]",
            lambda r, n: [(_signed(r, log_uniform(r, 8.0000001, 700.0)),) for _ in range(n)],
            1500, own=False))
        out.append(_stratum2(
            fn, "tail:[8,16)",
            lambda r, n: [(_signed(r, full_mantissa(r, 8.0, 16.0)),) for _ in range(n)],
            800, own=True))
        top = I0_OVF if fn in ("i0", "i0e") else I1_OVF
        out.append(_stratum2(
            fn, "tail:(700,ovf]",
            lambda r, n, top=top: [(_signed(r, full_mantissa(r, 700.0, top)),) for _ in range(n)],
            600, own=False))
        if fn.endswith("e"):
            out.append(_stratum2(
                fn, "tail:(ovf,DBL_MAX]",
                lambda r, n: [(_signed(r, log_uniform(r, 720.0, DBL_MAX)),) for _ in range(n)],
                600, own=False))
    return out


# ---- seam-audit fields (#82 T2; tools/hunt/seam_audit.py docstring) ----------
for _s in SEAMS:
    if _s["name"] == "dbl_max":
        _s["sides"] = "below"  # no double above DBL_MAX
