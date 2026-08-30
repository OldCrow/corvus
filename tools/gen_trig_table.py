#!/usr/bin/env python3
"""Generate src/trig_data.h -- constants for the corvus cos/sin kernels'
small-argument region (|x| <= 2^23).

PORT of libstats scripts/gen_trig_cleanroom_table.py (same owner, MIT,
clean-room derived -- see libstats docs/NEON_TRIG_DERIVATION.md for the
mathematics and docs/NEON_TRIG_DIVERGENCE_AUDIT.md for the point-by-point
comparison against ARM optimized-routines advsimd sin.c/cos.c confirming no
shared expression; no third-party source). The math here is byte-for-byte
that generator's math -- same precision, same splits, same fits -- so the
emitted constants are REQUIRED to be bit-identical to the donor values both
consumers already ship; the EXPECTED table below freezes those values and
the self-check exits non-zero on any deviation. A port that "almost"
reproduces the donor constants would silently invalidate the donor kernels'
validation history, which is what S3 transcribes against.

Design (donor DERIVATION.md sections 1-4):
    Quadrant reduction x = n*(pi/2) + r, n = round(x * 2/pi), |r| <= pi/4;
    quadrant q = n mod 4 selects/negates the two parity cores
        sin(r) = r + r*(u*P(u)),  cos(r) = 1 + u*Q(u),  u = r^2.
    pi/2 split into 4 parts, 30 significant bits each (last part full
    precision): every product n*p_k is EXACT for |n| < 2^(53-30) = 2^23
    WITHOUT FMA (30+23 bits fit in 53) -- the non-FMA safety of the
    reduction is structural. Domain of this region: |x| <= 2^23
    (n_max = 5,340,354 < 2^23). Beyond it corvus routes to the
    Payne-Hanek region (src/trig_ph_data.h), NOT to a scalar fixup.
    P, Q are degree-6 near-minimax fits (mpmath chebyfit at 320-bit
    precision) on u in [0, (pi/4 * (1+1e-6))^2]; Q[0] == -1/2 exactly is
    load-bearing (the kernel's exact 1 - u/2 head/tail split).

Usage:
    python tools/gen_trig_table.py > src/trig_data.h
"""

import math
import sys

import mpmath as mp

mp.mp.prec = 320

SIG_BITS = 30   # per-part significant bits -> n*part exact for |n| < 2^23
N_PARTS = 4     # donor measured: 3 parts fail (1.6 ULP stress); 2 catastrophic
D_MAX = float(2 ** 23)
DEG = 6
FIT_PAD = 1.0 + 1e-6

PIO2 = mp.pi / 2
TWO_OVER_PI = float(mp.mpf(2) / mp.pi)
U_MAX = float((mp.pi / 4 * mp.mpf(FIT_PAD)) ** 2)

# The certified donor constants (libstats/libhmm trig_cleanroom_data.inc,
# payload-identical in both consumers, audit chain in PLAN.md S1). The port
# must reproduce these bit-for-bit.
EXPECTED = {
    "TwoOverPi": ["0x1.45f306dc9c883p-1"],
    "Pio2": [
        "0x1.921fb54800000p+0",
        "-0x1.de973dc800000p-31",
        "-0x1.9d9cceb800000p-62",
        "-0x1.1fc8f8cbb5bf7p-93",
    ],
    "SinC": [
        "-0x1.5555555555555p-3",
        "0x1.1111111111110p-7",
        "-0x1.a01a01a019938p-13",
        "0x1.71de3a5460952p-19",
        "-0x1.ae645412c46d3p-26",
        "0x1.61217f0abf087p-33",
        "-0x1.ab17d393166dcp-41",
    ],
    "CosC": [
        "-0x1.0000000000000p-1",
        "0x1.5555555555551p-5",
        "-0x1.6c16c16c15d79p-10",
        "0x1.a01a019de130cp-16",
        "-0x1.27e4f8e4a1e4cp-22",
        "0x1.1eea7f24ce4e2p-29",
        "-0x1.8ff9d3c0d9835p-37",
    ],
}


def as_hex(x: float) -> str:
    if x == 0.0 and math.copysign(1.0, x) < 0:
        return "-0x0p+0"
    return float(x).hex()


def split_pio2(nparts: int, sig: int):
    parts = []
    R = mp.mpf(PIO2)
    for _ in range(nparts - 1):
        _, e = mp.frexp(R)
        g = int(e) - sig
        q = int(mp.nint(R / mp.mpf(2) ** g))
        assert abs(q) < 2 ** (sig + 1)
        p = float(mp.mpf(q) * mp.mpf(2) ** g)
        assert mp.mpf(p) == mp.mpf(q) * mp.mpf(2) ** g
        assert abs(q) * (2 ** (53 - sig) - 1) < 2 ** 53 or q == 0
        parts.append(p)
        R = R - mp.mpf(p)
    parts.append(float(R))
    resid = abs(R - mp.mpf(float(R)))
    return parts, resid


def h_sin(u):
    u = mp.mpf(u)
    if u < mp.mpf(2) ** -80:
        return mp.mpf(-1) / 6 + u / 120 - u ** 2 / 5040
    r = mp.sqrt(u)
    return (mp.sin(r) - r) / (u * r)


def h_cos(u):
    u = mp.mpf(u)
    if u < mp.mpf(2) ** -80:
        return mp.mpf(-1) / 2 + u / 24 - u ** 2 / 720
    return (mp.cos(mp.sqrt(u)) - 1) / u


def fit_poly(h, deg):
    coeffs = mp.chebyfit(h, [0, U_MAX], deg + 1)  # highest-degree first
    return list(reversed([float(c) for c in coeffs]))


def poly_mp(cs, u):
    acc = mp.mpf(0)
    for c in reversed(cs):
        acc = acc * u + mp.mpf(c)
    return acc


def main():
    fails = []

    parts, resid = split_pio2(N_PARTS, SIG_BITS)
    if not resid < mp.mpf(2) ** -140:
        fails.append(f"split truncation 2^{float(mp.log(resid, 2)):.1f}, want < 2^-140")
    nmax = int(mp.nint(mp.mpf(D_MAX) * 2 / mp.pi))
    if not nmax < 2 ** (53 - SIG_BITS):
        fails.append(f"n_max = {nmax} breaks the exact-product lemma")

    sin_c = fit_poly(h_sin, DEG)
    cos_c = fit_poly(h_cos, DEG)

    worst_s = mp.mpf(0)
    worst_c = mp.mpf(0)
    for i in range(2001):
        u = mp.mpf(U_MAX) * i / 2000
        r = mp.sqrt(u)
        if r > 0:
            s_approx = r + r * u * poly_mp(sin_c, u)
            worst_s = max(worst_s, abs(s_approx / mp.sin(r) - 1))
        c_approx = 1 + u * poly_mp(cos_c, u)
        worst_c = max(worst_c, abs(c_approx / mp.cos(r) - 1))
    if not worst_s < mp.mpf(2) ** -56:
        fails.append(f"sin fit 2^{float(mp.log(worst_s, 2)):.1f}, want < 2^-56")
    if not worst_c < mp.mpf(2) ** -56:
        fails.append(f"cos fit 2^{float(mp.log(worst_c, 2)):.1f}, want < 2^-56")
    if cos_c[0] != -0.5:
        fails.append(f"CosC[0] = {cos_c[0].hex()} != -0.5 (head/tail split needs it exact)")

    generated = {
        "TwoOverPi": [TWO_OVER_PI],
        "Pio2": parts,
        "SinC": sin_c,
        "CosC": cos_c,
    }
    for key, vals in generated.items():
        exp = EXPECTED[key]
        if len(exp) != len(vals):
            fails.append(f"{key}: length {len(vals)} != expected {len(exp)}")
            continue
        for i, (g, ex) in enumerate(zip(vals, exp)):
            if g.hex() != float.fromhex(ex).hex():
                fails.append(f"{key}[{i}]: generated {g.hex()} != donor {ex}")

    for f in fails:
        print(f"self-check FAILED: {f}", file=sys.stderr)
    if fails:
        return 1
    print(
        f"self-check OK: {N_PARTS}-part split (resid 2^{float(mp.log(resid, 2)):.1f}), "
        f"deg {DEG}/{DEG}, sin fit 2^{float(mp.log(worst_s, 2)):.1f}, "
        f"cos fit 2^{float(mp.log(worst_c, 2)):.1f}; all constants bit-identical "
        f"to the certified donor values",
        file=sys.stderr,
    )

    print("// Auto-generated by tools/gen_trig_table.py. DO NOT EDIT.")
    print("// Small-argument-region constants for the cos/sin kernels")
    print("// (src/trig-inl.h). Clean-room provenance and the derivation live in")
    print("// the generator's docstring; the values are bit-identical to the")
    print("// certified donor constants (generator self-check).")
    print("#ifndef CORVUS_TRIG_DATA_H_")
    print("#define CORVUS_TRIG_DATA_H_")
    print()
    print("namespace corvus::detail {")
    print()
    print("// Small-region domain bound: n*part_k products are exact only for")
    print("// |n| < 2^23. NOT a public domain limit -- beyond it the kernel")
    print("// routes to the Payne-Hanek region (trig_ph_data.h).")
    print(f"inline constexpr double kTrigDMax = {as_hex(D_MAX)};  // 2^23")
    print(f"inline constexpr double kTrigTwoOverPi = {as_hex(TWO_OVER_PI)};")
    print()
    print("// pi/2 split: parts 0-2 carry 30 significant bits each, so n*part_k")
    print("// is exact for |n| < 2^23 WITHOUT FMA; part 3 is a full 53-bit tail")
    print("// whose product rounds (<= 2^-124 absolute, below the reduction's")
    print("// ~2^-60 error floor). Subtract in order, compensated.")
    print(f"inline constexpr double kTrigPio2[{N_PARTS}] = {{")
    for p in parts:
        print(f"    {as_hex(p)},")
    print("};")
    print()
    print("// Degree-6 near-minimax parity cores on u = r^2 in [0, (pi/4)^2")
    print("// (padded)]: sin(r) = r + r*(u*P(u)); cos(r) = 1 + u*Q(u).")
    print("// kTrigCosC[0] == -1/2 EXACTLY (generator-asserted): the kernel's")
    print("// exact 1 - u/2 head/tail split depends on it.")
    print(f"inline constexpr double kTrigSinC[{DEG + 1}] = {{")
    for c in sin_c:
        print(f"    {as_hex(c)},")
    print("};")
    print(f"inline constexpr double kTrigCosC[{DEG + 1}] = {{")
    for c in cos_c:
        print(f"    {as_hex(c)},")
    print("};")
    print()
    print("}  // namespace corvus::detail")
    print()
    print("#endif  // CORVUS_TRIG_DATA_H_")
    return 0


if __name__ == "__main__":
    sys.exit(main())
