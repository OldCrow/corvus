#!/usr/bin/env python3
"""Generate src/trig_ph_data.inc -- Payne-Hanek windows for the cos/sin
kernels' large-argument region (|x| > 2^23).

Clean-room: the windowing scheme is textbook Payne-Hanek (Payne & Hanek
1983; Ng 1992, "Argument Reduction for Huge Arguments: Good to the Last
Bit"; Muller et al., Handbook of Floating-Point Arithmetic). mpmath supplies
high-precision VALUES of pi only; no third-party reduction implementation
was consulted. The shared frame lives in tools/trig_common.py.

Method (kernel-side derivation lives in src/trig-inl.h when S3 lands it):
    x = M * 2^(e-52), M in [2^52, 2^53), e in [23, 1023].
    Bits of 2/pi at weight 2^-j with j <= e-54 contribute exact multiples
    of 4 to x*(2/pi) and are dropped; the remainder, scaled to
    W_e = (2^(e-52) * 2/pi) mod 4 in [0, 4), is stored per exponent as
    FOUR 53-bit chunks at FIXED weights (2^1..2^-51, 2^-52..2^-104,
    2^-105..2^-157, 2^-158..2^-210) -- every chunk is a 53-bit integer
    times a power of two, exactly representable, and every product
    M * chunk is an EXACT double-double via TwoProd (no 26-bit splitting).
    The kernel accumulates the products as an exact expansion (integer
    stripping hi - 4*Round(hi/4) is exact for any double), so cancellation
    costs nothing; only the final dd compression rounds, RELATIVE to the
    already-cancelled f - n. Window truncation contributes
    |M| * 2^-210 < 2^-157 absolute in f -- at binary64's deepest
    cancellation (|f - n| ~ 2^-61.5, e = 849 row) that is ~2^-96 relative
    in r, verified below at that exact point.

Self-check (exits non-zero on any failure):
    1. every chunk reconstructs its integer slice exactly, and the four
       chunks sum to the truncated window exactly (all 1001 rows);
    2. window truncation is in [0, 2^-210) against 1400-bit 2/pi;
    3. end-to-end: for every row, random mantissas plus the row's
       continued-fraction worst cases -- the exact-rational table reduction
       must match the 1400-bit true reduction to <= 2^-150 absolute in f,
       and to <= 2^-88 relative in r at every probed point;
    4. the literature worst case (M = 0x16ac5b262ca1ff, e = 849,
       |r| = 2^-60.89) is probed explicitly.

Usage:
    python tools/gen_trig_ph_table.py > src/trig_ph_data.inc
"""

import random
import sys
from fractions import Fraction

import mpmath as mp

from trig_common import (
    PH_CHUNKS,
    PH_E_MAX,
    PH_E_MIN,
    PH_FRAC_BITS,
    window_chunks,
    window_fraction,
    window_int,
    worst_m_for_exponent,
)

SEED = 20260830
ROWS = PH_E_MAX - PH_E_MIN + 1
PROBES_PER_ROW = 2   # random mantissas per row, plus the CF worst cases
WORKPREC = 1400


def as_hex(x: float) -> str:
    return float(x).hex()


def true_f_and_r(m: int, e: int):
    """(x*2/pi) mod 4 and the reduced argument, at WORKPREC bits."""
    x = mp.mpf(m) * mp.mpf(2) ** (e - 52)
    t = x * 2 / mp.pi
    f = t - 4 * mp.floor(t / 4)
    n = mp.nint(f)
    return f, (f - n) * mp.pi / 2, int(n) % 4


def self_check():
    rng = random.Random(SEED)
    fails = []
    worst_abs = mp.mpf(0)
    worst_rel = mp.mpf(0)
    deepest = (0, 0, 0)  # (depth_bits, m, e)

    with mp.workprec(WORKPREC):
        two_over_pi = mp.mpf(2) / mp.pi
        for e in range(PH_E_MIN, PH_E_MAX + 1):
            # 1. chunk exactness / reconstruction (pure integer identity).
            chunks = window_chunks(e)
            if sum(Fraction(c) for c in chunks) != window_fraction(e):
                fails.append(f"e={e}: chunks do not reconstruct the window")
                continue

            # 2. truncation of the window itself.
            w_true = (mp.mpf(2) ** (e - 52) * two_over_pi) % 4
            tail = w_true - mp.mpf(window_int(e)) / mp.mpf(2) ** PH_FRAC_BITS
            if not (0 <= tail < mp.mpf(2) ** -PH_FRAC_BITS):
                fails.append(f"e={e}: window truncation out of range: {tail}")

            # 3. end-to-end probes: random + CF worst cases.
            probes = [rng.randrange(1 << 52, 1 << 53)
                      for _ in range(PROBES_PER_ROW)]
            worst = worst_m_for_exponent(e, 2)
            probes += [m for m, _ in worst]
            if worst and worst[0][1] > deepest[0]:
                deepest = (worst[0][1], worst[0][0], e)
            for m in probes:
                f_tab = m * window_fraction(e) % 4
                f_true, r_true, q_true = true_f_and_r(m, e)
                # wrap-aware absolute error in f (both live in [0, 4)).
                d = (mp.mpf(f_tab.numerator) / f_tab.denominator - f_true) % 4
                d = min(d, 4 - d)
                worst_abs = max(worst_abs, d)
                if d > mp.mpf(2) ** -150:
                    fails.append(f"e={e} m={m:#x}: |f_tab - f_true| = "
                                 f"2^{float(mp.log(d, 2)):.1f} > 2^-150")
                n_tab = round(f_tab)
                r_tab = (mp.mpf((f_tab - n_tab).numerator)
                         / (f_tab - n_tab).denominator) * mp.pi / 2
                if r_true != 0:
                    rel = abs(r_tab / r_true - 1)
                    worst_rel = max(worst_rel, rel)
                    if rel > mp.mpf(2) ** -88:
                        fails.append(
                            f"e={e} m={m:#x}: r rel err "
                            f"2^{float(mp.log(rel, 2)):.1f} > 2^-88 "
                            f"(|r| = 2^{float(mp.log(abs(r_true), 2)):.1f})")
                if n_tab % 4 != q_true:
                    # legitimate only if f sits essentially on a half-integer
                    half_dist = abs(f_true - mp.nint(f_true))
                    if abs(half_dist - mp.mpf('0.5')) > mp.mpf(2) ** -60:
                        fails.append(f"e={e} m={m:#x}: quadrant "
                                     f"{n_tab % 4} != {q_true}")

        # 4. the literature worst case, explicitly.
        m_lit, e_lit = 0x16AC5B262CA1FF, 849
        f_tab = m_lit * window_fraction(e_lit) % 4
        n_tab = round(f_tab)
        _, r_true, _ = true_f_and_r(m_lit, e_lit)
        r_tab = (mp.mpf((f_tab - n_tab).numerator)
                 / (f_tab - n_tab).denominator) * mp.pi / 2
        rel = abs(r_tab / r_true - 1)
        if not (mp.mpf(2) ** -61 < abs(r_true) < mp.mpf(2) ** -60.5):
            fails.append("literature worst case: |r| not at 2^-60.89 -- "
                         "wrong point or wrong math")
        if rel > mp.mpf(2) ** -88:
            fails.append(f"literature worst case: rel err "
                         f"2^{float(mp.log(rel, 2)):.1f} > 2^-88")

    return fails, worst_abs, worst_rel, deepest


def main():
    fails, worst_abs, worst_rel, deepest = self_check()
    for f in fails[:40]:
        print(f"self-check FAILED: {f}", file=sys.stderr)
    if fails:
        return 1
    print(
        f"self-check OK: {ROWS} rows x {PH_CHUNKS} chunks; worst |f| error "
        f"2^{float(mp.log(worst_abs, 2)):.1f}, worst r rel error "
        f"2^{float(mp.log(worst_rel, 2)):.1f}; deepest cancellation probed "
        f"{deepest[0]} bits (m={deepest[1]:#x}, e={deepest[2]}); literature "
        f"worst case (e=849) verified",
        file=sys.stderr,
    )

    with mp.workprec(WORKPREC):
        pio2 = mp.pi / 2
        pio2_hi = float(pio2)
        pio2_lo = float(pio2 - mp.mpf(pio2_hi))
        dd_err = abs((mp.mpf(pio2_hi) + mp.mpf(pio2_lo)) / pio2 - 1)
        if dd_err > mp.mpf(2) ** -104:
            print(f"self-check FAILED: pi/2 dd pair error 2^"
                  f"{float(mp.log(dd_err, 2)):.1f}", file=sys.stderr)
            return 1

    print("// Auto-generated by tools/gen_trig_ph_table.py. DO NOT EDIT.")
    print("// Payne-Hanek windows for the cos/sin large-argument region:")
    print("// W_e = (2^(e-52) * 2/pi) mod 4 as four 53-bit chunks at fixed")
    print("// weights (2^1..2^-51 / ..2^-104 / ..2^-157 / ..2^-210), one row")
    print("// per unbiased exponent e = 23..1023, indexed by e - kTrigPhEMin.")
    print("// Method, bounds, and the self-check: tools/gen_trig_ph_table.py.")
    print(f"static_assert(kTrigPhEMin == {PH_E_MIN}, "
          "\"regenerate src/trig_ph_data.inc\");")
    print(f"static_assert(kTrigPhRows == {ROWS}, "
          "\"regenerate src/trig_ph_data.inc\");")
    print(f"const double kTrigPio2DdHi = {as_hex(pio2_hi)};")
    print(f"const double kTrigPio2DdLo = {as_hex(pio2_lo)};")
    all_chunks = [window_chunks(e) for e in range(PH_E_MIN, PH_E_MAX + 1)]
    for c in range(PH_CHUNKS):
        print(f"alignas(64) const double kTrigPhTab{c}[{ROWS}] = {{")
        for row in all_chunks:
            print(f"    {as_hex(row[c])},")
        print("};")
    return 0


if __name__ == "__main__":
    sys.exit(main())
