"""Shared trig-reduction math for corvus generators.

Used by gen_trig_ph_table.py (the Payne-Hanek window table and its
self-check) and gen_trig_reference.py (worst-case point selection), so the
two scripts cannot drift apart on what "the bits of 2/pi" or "the worst
M for exponent e" mean.

Clean-room: everything here is textbook -- Payne-Hanek windowing (Payne &
Hanek 1983; Ng 1992), continued-fraction best rational approximation
(three-distance / convergent theory). mpmath supplies high-precision VALUES
of pi only.

Frames used throughout:
  x = M * 2^(e-52),  M integer in [2^52, 2^53), e = unbiased exponent.
  W_e = (2^(e-52) * (2/pi)) mod 4   -- the per-exponent window, in [0, 4).
    Bits of 2/pi at weight 2^-j with j <= e-54 contribute exact multiples
    of 4 to x*(2/pi) (M * 2^(e-52-j), e-52-j >= 2), so they are dropped.
  f = (M * W_e) mod 4 == (x * 2/pi) mod 4;  n = round(f);  quadrant n mod 4;
  r = (f - n) * pi/2.
"""

import math
from fractions import Fraction

import mpmath as mp

# Bits of 2/pi kept per window: 4 chunks x 53 bits, spanning weights
# 2^1 .. 2^-210 of W_e. Truncation in f is bounded by M * 2^-210 < 2^-157.
PH_FRAC_BITS = 210
PH_CHUNKS = 4
PH_E_MIN = 23     # smallest exponent routed to the large-argument region
PH_E_MAX = 1023   # largest finite-double exponent

# Working precision for the master bit string: enough for e = 1023 windows
# (K >= e + 158) with slack.
_K = 1280


def two_over_pi_bits() -> int:
    """floor((2/pi) * 2^_K): the master bit string, as a Python int."""
    with mp.workprec(_K + 64):
        return int(mp.floor(mp.mpf(2) / mp.pi * mp.mpf(2) ** _K))


_B = two_over_pi_bits()


def window_int(e: int) -> int:
    """W_e * 2^PH_FRAC_BITS truncated to an integer (2 integer bits +
    PH_FRAC_BITS fractional bits of the window). The TABLE uses only
    e >= PH_E_MIN; smaller exponents are legal here because the worst-case
    search (gen_trig_reference.py) also probes the small region, where the
    same frame applies with no dropped bits."""
    assert 0 <= e <= PH_E_MAX
    if e <= 54:
        kept = _B  # no bits dropped
    else:
        kept = _B % (1 << (_K - (e - 54)))  # drop weights 2^-j, j <= e-54
    shift = _K - (e - 52) - PH_FRAC_BITS
    assert shift >= 0, e
    return kept >> shift


def window_chunks(e: int):
    """The four 53-bit chunk doubles c0..c3 of W_e (exactly representable:
    53-bit integer times a power of two)."""
    v = window_int(e)
    assert v < 1 << (PH_CHUNKS * 53)  # 2 integer + 210 fractional bits = 212
    chunks = []
    for i in range(PH_CHUNKS):
        sh = (PH_CHUNKS - 1 - i) * 53
        sl = (v >> sh) & ((1 << 53) - 1)
        c = math.ldexp(sl, sh - PH_FRAC_BITS)
        assert int(c * 2.0 ** (PH_FRAC_BITS - sh)) == sl, (e, i)
        chunks.append(c)
    return chunks


def window_fraction(e: int) -> Fraction:
    """W_e truncated at PH_FRAC_BITS, as an exact rational."""
    return Fraction(window_int(e), 1 << PH_FRAC_BITS)


def theta_fraction(e: int) -> Fraction:
    """theta_e = W_e mod 1 (exact rational of the truncated window). The
    distance of M*theta_e to the nearest integer equals the distance of
    M*W_e (mod 4) to the nearest integer, which is what cancellation
    depth in the reduction is about."""
    w = window_int(e)
    return Fraction(w % (1 << PH_FRAC_BITS), 1 << PH_FRAC_BITS)


def _convergent_denominators(theta: Fraction, qmax: int):
    """Denominators of the continued-fraction convergents of theta, up to
    qmax. Standard recurrence; theta is an exact rational so the expansion
    is exact and finite."""
    num, den = theta.numerator, theta.denominator
    qs = []
    q_prev2, q_prev1 = 1, 0   # q_{-2}, q_{-1}
    while den:
        a = num // den
        num, den = den, num - a * den
        q = a * q_prev1 + q_prev2
        q_prev2, q_prev1 = q_prev1, q
        if q > qmax:
            break
        if q > 1:
            qs.append(q)
    return qs


def dist_to_nearest_int(fr: Fraction) -> Fraction:
    n = round(fr)
    return abs(fr - n)


def worst_m_for_exponent(e: int, count: int = 2):
    """Mantissa integers M in [2^52, 2^53) with the deepest reduction
    cancellation for exponent e: minimize ||M * theta_e||. Candidates are
    the CF convergent denominators of theta_e (lifted into range by an
    integer multiplier when small). Returns [(M, depth_bits), ...] sorted
    deepest-first; depth_bits = -log2 ||M * theta_e||."""
    theta = theta_fraction(e)
    if theta == 0:
        return []
    lo, hi = 1 << 52, (1 << 53) - 1
    best = {}
    for q in _convergent_denominators(theta, hi):
        cands = []
        if q >= lo:
            cands.append(q)
        else:
            for k in (-(-lo // q), -(-lo // q) + 1):  # ceil(lo/q), +1
                if lo <= k * q <= hi:
                    cands.append(k * q)
        for m in cands:
            d = dist_to_nearest_int(m * theta)
            if d > 0:
                # depth in bits, via integer bit lengths (d can underflow
                # a float): -log2(p/q) = len(q) - len(p) +- 1.
                best[m] = d.denominator.bit_length() - d.numerator.bit_length()
    ranked = sorted(best.items(), key=lambda kv: -kv[1])
    return ranked[:count]


def reduce_exact(m: int, e: int):
    """Simulate the table-backed reduction EXACTLY (rational arithmetic):
    returns (n mod 4, r_table as Fraction * (pi/2 NOT applied), i.e. the
    exact f - n). The kernel's fp rounding is the gates' concern; this
    isolates the table truncation."""
    f = (m * window_fraction(e)) % 4
    n = round(f)
    return n % 4, f - n
