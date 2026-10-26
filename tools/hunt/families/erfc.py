"""erfc hunt strata. Bounds: ACCURACY.md erfc -- 1 ULP for |x| <= 6, 2 ULP for
normal results with x > 6, 1 ULP for subnormal results."""
import mpmath as mp

from core import Stratum, full_mantissa, neighbours

# Seams: src/erfc_tail_data.h tail breakpoints and clamp; derived output seams.
SEAMS = [
    {"fn": "erfc", "name": "core_tail", "coord": "x", "value": 6.0, "source": "kErfcTailLo (erfc_tail_data.h)"},
    {"fn": "erfc", "name": "tail_b1", "coord": "x", "value": 10.0, "source": "kErfcTailBoundN (erfc_tail_data.h)"},
    {"fn": "erfc", "name": "tail_b2", "coord": "x", "value": 17.0, "source": "kErfcTailBoundN (erfc_tail_data.h)"},
    {"fn": "erfc", "name": "min_normal_out", "coord": "x", "value": 26.543258234, "source": "erfc(x) = 2^-1022"},
    {"fn": "erfc", "name": "zero_out", "coord": "x", "value": 27.226017, "source": "erfc(x) < 2^-1075"},
    {"fn": "erfc", "name": "clamp", "coord": "x", "value": 28.0, "source": "kErfcTailHi (erfc_tail_data.h)"},
]


def _oracle(args, dps):
    with mp.workdps(dps):
        return mp.erfc(mp.mpf(args[0]))


def _bound(args, v):
    x = args[0]
    if abs(x) <= 6:
        return 1
    return 1 if abs(v) < mp.mpf(2) ** -1022 else 2


def strata():
    out = []
    for s in SEAMS:
        t = s["value"]
        out.append(Stratum("erfc", f"seam:{s['name']}:bits",
                           lambda r, n, t=t: [(x,) for x in neighbours(t, 64)],
                           129, _oracle, _bound, seam=s["name"]))
        out.append(Stratum("erfc", f"seam:{s['name']}:band",
                           lambda r, n, t=t: [(full_mantissa(r, t * 0.97, t * 1.03),) for _ in range(n)],
                           3000, _oracle, _bound, seam=s["name"]))
    out.append(Stratum("erfc", "top_subnormal_binade",
                       lambda r, n: [(full_mantissa(r, 26.5433, 26.556),) for _ in range(n)],
                       6000, _oracle, _bound, note="#74 band"))
    return out
