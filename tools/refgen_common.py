"""Shared reference-writer helpers (#13 N14).

float(mpf) double-rounds subnormal results: mpmath's to_float rounds the
mantissa to 53 bits, then ldexp rounds AGAIN onto the 2^-1074 grid — two
roundings, up to 1 ulp off correctly rounded in the subnormal band.
round_to_double rounds once.
"""

import math

import mpmath as mp

_MIN_NORMAL = 2.0 ** -1022


def round_to_double(v):
    """Nearest double to mpf/int/float v, single rounding, ties-to-even.

    Normal range: float(mpf) IS a single correct rounding (53-bit mantissa,
    exponent representable), so it is used directly. Subnormal range: scale
    by 2^1074 (exact — pure exponent shift in binary fp), round to an
    integer once (ties-to-even, matching IEEE), ldexp back (exact: the
    integer is <= 2^52).
    """
    v = mp.mpf(v)
    if not mp.isfinite(v):
        return float(v)
    if v == 0:
        # An mpf zero carries no usable sign; callers owning a signed-zero
        # contract (erfinv(-0) = -0) must emit that sign explicitly.
        return 0.0
    if abs(v) >= _MIN_NORMAL:
        return float(v)
    # t = |v| * 2^1074 < 2^52. frac = t - n is exact at this workprec:
    # t's mantissa spans <= prec bits below an MSB <= 2^51, so the
    # fractional part needs < prec bits — no bits are dropped, and the
    # tie comparison against exactly 0.5 is therefore decisive.
    with mp.workprec(mp.mp.prec + 64):
        t = abs(v) * mp.mpf(2) ** 1074
        n = int(mp.floor(t))
        frac = t - n
        if frac > 0.5 or (frac == 0.5 and (n & 1)):
            n += 1
    return math.copysign(math.ldexp(float(n), -1074), 1.0 if v > 0 else -1.0)


def _self_test():
    # Tie at the +0/min-subnormal boundary: 2^-1075 is exactly half the
    # smallest subnormal; ties-to-even rounds to 0 (even).
    assert round_to_double(mp.mpf(2) ** -1075) == 0.0
    # Just above the tie rounds up to the min subnormal.
    assert round_to_double(mp.mpf(2) ** -1075 * (1 + mp.mpf(2) ** -40)) == 2.0 ** -1074
    # Just below rounds down to 0.
    assert round_to_double(mp.mpf(2) ** -1075 * (1 - mp.mpf(2) ** -40)) == 0.0
    # Odd-multiple tie rounds to the EVEN neighbor (up): 3 * 2^-1075 -> 2 * 2^-1074.
    assert round_to_double(mp.mpf(3) * mp.mpf(2) ** -1075) == 2 * 2.0 ** -1074
    # Even-multiple neighborhood: 5 * 2^-1075 -> tie between 2 and 3 ulps -> 2 (even...
    # 5/2 = 2.5 -> even = 2).
    assert round_to_double(mp.mpf(5) * mp.mpf(2) ** -1075) == 2 * 2.0 ** -1074
    # Subnormal/normal boundary: rounding up across it is exact.
    below = mp.mpf(2) ** -1022 * (1 - mp.mpf(2) ** -70)
    assert round_to_double(below) == 2.0 ** -1022
    # Normal range unchanged versus float().
    assert round_to_double(mp.mpf(1) / 3) == float(mp.mpf(1) / 3)
    # Sign preserved through the subnormal path.
    assert round_to_double(-(mp.mpf(3) * mp.mpf(2) ** -1075)) == -(2 * 2.0 ** -1074)
    # The double-rounding defect this module exists to fix: a value whose
    # 53-bit rounding lands exactly on a half-ulp of the subnormal grid
    # (second rounding then goes the wrong way ~half the time). Construct
    # one: v = (2^-1070)*(1 + 2^-53 + 2^-90). 53-bit rounding of the
    # mantissa gives 1 + 2^-52 exactly (tie broken up by the 2^-90 bit is
    # already above the tie): float(v) via to_float rounds mantissa to
    # 1 + 2^-52, then ldexp onto the 2^-1074 grid rounds 16 + 2^-48-ish...
    # the ONE-rounding answer: v*2^1074 = 16*(1 + 2^-53 + 2^-90) =
    # 16 + 2^-49 + 2^-86, nearest integer = 16.
    v = mp.mpf(2) ** -1070 * (1 + mp.mpf(2) ** -53 + mp.mpf(2) ** -90)
    assert round_to_double(v) == math.ldexp(16.0, -1074)


_self_test()
