"""erf hunt strata. Bound: ACCURACY.md erf -- max 1 ULP over the full domain,
subnormal inputs included (row 11, #61)."""
import mpmath as mp

from core import Stratum, bits, from_bits, full_mantissa, log_uniform, neighbours, subnormal

# Seams: src/erf_core-inl.h table index Round(256*a) (cell edges at
# (j + 1/2)/256), src/erf_data.inc kErfAMax saturation, the subnormal edge.
SEAMS = [
    {"fn": "erf", "name": "cell_0_1", "coord": "x", "value": 1 / 512, "source": "ErfTableCore Round(a*256): j=0|1 edge"},
    {"fn": "erf", "name": "cell_127_128", "coord": "x", "value": 127.5 / 256, "source": "ErfTableCore cell edge near 1/2"},
    {"fn": "erf", "name": "cell_255_256", "coord": "x", "value": 255.5 / 256, "source": "ErfTableCore cell edge near 1"},
    {"fn": "erf", "name": "cell_511_512", "coord": "x", "value": 511.5 / 256, "source": "ErfTableCore cell edge near 2"},
    {"fn": "erf", "name": "cell_1023_1024", "coord": "x", "value": 1023.5 / 256, "source": "ErfTableCore cell edge near 4"},
    {"fn": "erf", "name": "amax", "coord": "x", "value": float.fromhex("0x1.7afb48dc96627p+2"), "source": "kErfAMax (erf_data.inc)"},
    {"fn": "erf", "name": "min_normal_in", "coord": "x", "value": 2.0 ** -1022, "source": "subnormal input edge"},
]


def _oracle(args, dps):
    with mp.workdps(dps + 10):
        return mp.erf(mp.mpf(args[0]))


def _oracle2(args, dps):
    # Independent form: regularized lower incomplete gamma P(1/2, x^2).
    x = mp.mpf(args[0])
    with mp.workdps(dps + 10):
        if x == 0:
            return x
        return mp.sign(x) * mp.gammainc(mp.mpf(0.5), 0, x * x, regularized=True)


def _bound(args, v):
    return 1


def _sgn(r, x):
    return -x if r.random() < 0.5 else x


def _cell_edges(r, n):
    out = []
    for _ in range(n):
        e = (r.randint(0, 1535) + 0.5) / 256
        out.append((_sgn(r, from_bits(bits(e) + r.randint(-(1 << 20), 1 << 20))),))
    return out


def _st(name, gen, n, seam=None, note=""):
    return Stratum("erf", name, gen, n, _oracle, _bound, oracle2=_oracle2, seam=seam, note=note)


def strata():
    out = []
    for s in SEAMS:
        t = s["value"]
        out.append(_st(f"seam:{s['name']}:bits",
                       lambda r, n, t=t: [(x,) for x in neighbours(t, 64)] + [(-x,) for x in neighbours(t, 64)],
                       258, s["name"]))
        out.append(_st(f"seam:{s['name']}:band",
                       lambda r, n, t=t: [(_sgn(r, full_mantissa(r, t * 0.95, t * 1.05)),) for _ in range(n)],
                       1500, s["name"]))
    out += [
        _st("cell_j0_uniform", lambda r, n: [(_sgn(r, full_mantissa(r, 0.0, 1 / 512)),) for _ in range(n)],
            5000, note="table cell j=0: dd 2/sqrt(pi) path since #80"),
        _st("cell_j0_log", lambda r, n: [(_sgn(r, log_uniform(r, 2.0 ** -1022, 1 / 512)),) for _ in range(n)],
            3000, note="tiny normals, erf ~ 2x/sqrt(pi)"),
        _st("cell_edges_all", _cell_edges, 6000, note="+-2^20 ulps of every (j+1/2)/256 edge"),
        _st("interior_0_2", lambda r, n: [(_sgn(r, full_mantissa(r, 0.0, 2.0)),) for _ in range(n)], 8000),
        _st("interior_2_6", lambda r, n: [(_sgn(r, full_mantissa(r, 2.0, 6.0)),) for _ in range(n)], 6000),
        _st("saturated", lambda r, n: [(_sgn(r, log_uniform(r, 5.9216, 1e300)),) for _ in range(n)], 500),
        _st("subnormal", lambda r, n: [(_sgn(r, subnormal(r)),) for _ in range(n)], 4000,
            note="subnormal in and (mostly) subnormal out"),
        _st("subnormal_top", lambda r, n: [(_sgn(r, subnormal(r, 3)),) for _ in range(n)], 2000,
            note="top subnormal binades; outputs cross 2^-1022"),
    ]
    return out

