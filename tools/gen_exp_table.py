#!/usr/bin/env python3
"""Generate src/exp_dd_data.inc -- tables for corvus's double-double exp kernel.

Clean-room: the argument reduction is textbook Cody-Waite and the table is
2^(j/N) evaluated by mpmath. No third-party exp implementation was consulted;
mpmath is used only as a high-precision oracle for VALUES.

Method (see src/exp_dd-inl.h for the kernel-side derivation):

    k = round(x * N/ln2),  j = k mod N,  e = (k - j)/N
    r = x - k*(ln2/N)      computed with ln2/N split as L1 + L2
    exp(x) = 2^e * 2^(j/N) * e^r

L1 carries only the top 34 significant bits, so k*L1 is EXACT for every |k|
this kernel can see (|x| <= 1100 => |k| < 2^18). L2 captures the next 53 bits;
the unrepresented remainder L3 = ln2/N - L1 - L2 is bounded below, and its
amplified contribution |k*L3| is what limits the reduction.

Emitted:
    kExpN            table size N (and the reduction's N/ln2 scale)
    kExpInvL         fl(N/ln2)
    kExpL1, kExpL2   the Cody-Waite split of ln2/N
    kExpTableHi/Lo   2^(j/N) as a dd pair, j = 0..N-1, two flat arrays so the
                     kernel does two stride-1 gathers (same shape as the erf
                     table's field gathers)

Usage:
    python3 tools/gen_exp_table.py > src/exp_dd_data.inc
"""

import struct
import sys

import mpmath as mp

mp.mp.dps = 60

N = 128            # table entries per octave; |r| <= ln2/(2N) = 0.00271
L1_ZERO_BITS = 19  # low mantissa bits cleared in L1 => |k| < 2^19 stays exact
X_MAX = 1100.0     # kernel clamps |x| here (exp saturates well inside it)
POLY_TERMS = 6     # e^r - 1 = r + r^2*(1/2 + ... + r^4/720); truncation r^7/5040


def as_hex(x: float) -> str:
    return float.hex(x)


def as_bits(x: float) -> int:
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def from_bits(u: int) -> float:
    return struct.unpack("<d", struct.pack("<Q", u))[0]


def split_ln2_over_n():
    """L1 (top bits, exactly multipliable by k) + L2 (next 53 bits)."""
    exact = mp.log(2) / N
    l1 = from_bits(as_bits(float(exact)) & ~((1 << L1_ZERO_BITS) - 1))
    l2 = float(exact - mp.mpf(l1))
    l3 = exact - mp.mpf(l1) - mp.mpf(l2)
    return l1, l2, l3, exact


def self_check(l1, l2, l3, exact):
    """Verify every claim the kernel's derivation block makes."""
    fails = []
    k_max = int(X_MAX * N / float(mp.log(2))) + 1

    # 1. k*L1 exact for all reachable k: L1's significand times k must still
    #    fit in 53 bits. Checked structurally (bit widths), not by sampling.
    l1_sig = as_bits(l1) & ((1 << 52) - 1) | (1 << 52)
    while l1_sig and l1_sig % 2 == 0:
        l1_sig >>= 1
    if l1_sig.bit_length() + k_max.bit_length() > 53:
        fails.append(
            f"k*L1 not exact: L1 has {l1_sig.bit_length()} significant bits, "
            f"k up to {k_max} needs {k_max.bit_length()}"
        )

    # 2. Amplified reduction remainder |k*L3|: the absolute error floor of r,
    #    and so the relative error floor of exp.
    kl3 = abs(mp.mpf(k_max) * l3)
    if kl3 > mp.mpf(2) ** -70:
        fails.append(f"|k*L3| = 2^{float(mp.log(kl3, 2)):.1f}, want <= 2^-70")

    # 3. Polynomial truncation over the reduced range.
    r_max = mp.log(2) / (2 * N)
    trunc = mp.mpf(0)
    for i in range(2001):
        r = -r_max + 2 * r_max * i / 2000
        series = r
        for n in range(2, POLY_TERMS + 1):
            series += r ** n / mp.factorial(n)
        trunc = max(trunc, abs(mp.e ** r - 1 - series))
    if trunc > mp.mpf(2) ** -70:
        fails.append(f"poly truncation 2^{float(mp.log(trunc, 2)):.1f}, want <= 2^-70")

    # 4. Table dd pairs represent 2^(j/N) to dd precision.
    worst_tbl = mp.mpf(0)
    for j in range(N):
        v = mp.mpf(2) ** (mp.mpf(j) / N)
        hi = float(v)
        lo = float(v - mp.mpf(hi))
        worst_tbl = max(worst_tbl, abs((mp.mpf(hi) + mp.mpf(lo) - v) / v))
    if worst_tbl > mp.mpf(2) ** -104:
        fails.append(f"table dd error 2^{float(mp.log(worst_tbl, 2)):.1f}")

    return fails, dict(k_max=k_max, kl3=kl3, trunc=trunc, worst_tbl=worst_tbl,
                       l1_bits=l1_sig.bit_length(), r_max=r_max, exact=exact)


def main():
    l1, l2, l3, exact = split_ln2_over_n()
    fails, info = self_check(l1, l2, l3, exact)
    for f in fails:
        print(f"self-check FAILED: {f}", file=sys.stderr)
    if fails:
        return 1
    print(
        f"self-check OK: L1 {info['l1_bits']} sig bits, |k| <= {info['k_max']}, "
        f"|k*L3| <= 2^{float(mp.log(info['kl3'], 2)):.1f}, "
        f"poly truncation <= 2^{float(mp.log(info['trunc'], 2)):.1f}, "
        f"table dd <= 2^{float(mp.log(info['worst_tbl'], 2)):.1f}",
        file=sys.stderr,
    )

    inv_l = float(N / mp.log(2))

    print("// Auto-generated exp tables for corvus. DO NOT EDIT --")
    print("// regenerate with tools/gen_exp_table.py.")
    print("// Cody-Waite reduction k = round(x*N/ln2) with ln2/N = L1 + L2 + O(2^-95);")
    print(f"// L1's low {L1_ZERO_BITS} mantissa bits are zero, so k*L1 is exact for")
    print(f"// |k| < 2^{L1_ZERO_BITS} (the kernel clamps |x| <= {X_MAX:g}, giving |k| <= {info['k_max']}).")
    print("// Table: 2^(j/N) as a double-double, split across two flat arrays.")
    # N drives compile-time shift counts in the kernel, so it lives in the
    # header as a constexpr; the .inc only pins it to what was generated.
    print(f"static_assert(kExpN == {N}, \"regenerate src/exp_dd_data.inc\");")
    print(f"static_assert(kExpXMax == {X_MAX:g}, \"regenerate src/exp_dd_data.inc\");")
    print(f"const double kExpInvL = {as_hex(inv_l)};   // fl(N/ln2)")
    print(f"const double kExpL1 = {as_hex(l1)};   // ln2/N, top {info['l1_bits']} bits")
    print(f"const double kExpL2 = {as_hex(l2)};")
    print(f"alignas(64) const double kExpTableHi[{N}] = {{")
    for j in range(N):
        print(f"    {as_hex(float(mp.mpf(2) ** (mp.mpf(j) / N)))},")
    print("};")
    print(f"alignas(64) const double kExpTableLo[{N}] = {{")
    for j in range(N):
        v = mp.mpf(2) ** (mp.mpf(j) / N)
        print(f"    {as_hex(float(v - mp.mpf(float(v))))},")
    print("};")
    return 0


if __name__ == "__main__":
    sys.exit(main())
