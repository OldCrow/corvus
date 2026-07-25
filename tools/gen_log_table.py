#!/usr/bin/env python3
"""Generate src/log_dd_data.inc -- tables for corvus's double-double log kernel.

Clean-room: table-plus-log1p is the textbook shape and every value comes from
mpmath. No third-party log implementation was consulted.

Method (see src/log_dd-inl.h for the kernel-side derivation):

    x = 2^k * m,  m in [0.70703125, 1.4140625)
    j = top 7 mantissa bits of x     (the slot index, before centring)
    r = R_j*m - 1                    exact as a dd: R_j*m via TwoProd, and
                                     p - 1 is exact by Sterbenz
    log(x) = k*ln2 + L_j + log1p(r),  L_j = -log(R_j)

Two choices here do the real work:

* The mantissa is CENTRED on 1 (m in [1/sqrt2-ish, sqrt2-ish], the split
  taken at the slot boundary 1+53/128 rather than at sqrt2 itself so the
  condition is exactly "j >= 53" and no slot straddles it). Without the
  centring, x slightly below 1 lands at k = -1 and log(x) comes out as
  ln2 + log(m) -- two terms near 0.693 cancelling to ~1e-16, which throws
  away most of a double-double's headroom precisely where relative accuracy
  matters most.

* Slots 0 and 127 -- the two adjacent to m = 1 -- get R_j = 1 and L_j = 0
  EXACTLY. Then r = m - 1 exactly (Sterbenz) and log(x) is log1p(r) alone,
  with no cancellation at all on either side of 1. This is why the kernel
  needs no special case for x near 1.

Emitted:
    kLogSplit               1 + 53/128, the m0 threshold for halving
    kLogLn2Hi/Lo            ln2 as a dd
    kLogTableR              R_j
    kLogTableLHi/LLo        -log(R_j) as a dd

Usage:
    python3 tools/gen_log_table.py > src/log_dd_data.inc
"""

import struct
import sys

import mpmath as mp

mp.mp.dps = 60

N = 128          # slots over m0 in [1, 2)
SPLIT_SLOT = 53  # j >= SPLIT_SLOT halves m0; 1+53/128 = 1.4140625 ~ sqrt(2)
POLY_LAST = 11   # log1p series carried through r^POLY_LAST


def as_hex(x: float) -> str:
    return float.hex(x)


def slot_range(j):
    """(lo, hi) of m for slot j, after the centring halve."""
    lo = mp.mpf(1) + mp.mpf(j) / N
    hi = mp.mpf(1) + mp.mpf(j + 1) / N
    if j >= SPLIT_SLOT:
        lo, hi = lo / 2, hi / 2
    return lo, hi


def recip(j):
    """R_j: 1 exactly next to m = 1, else the rounded reciprocal of the centre."""
    if j == 0 or j == N - 1:
        return 1.0
    lo, hi = slot_range(j)
    return float(1 / ((lo + hi) / 2))


def self_check(rs):
    fails = []

    # 1. |r| bound, and p = R_j*m staying inside [0.5, 2] so that p - 1 is
    #    exact by Sterbenz -- the claim that makes r exact as a dd.
    worst_r = mp.mpf(0)
    for j in range(N):
        lo, hi = slot_range(j)
        for i in range(129):
            m = lo + (hi - lo) * i / 128
            p = mp.mpf(rs[j]) * m
            if not (mp.mpf(0.5) <= p <= 2):
                fails.append(f"slot {j}: p = {float(p)} outside [0.5, 2]")
                break
            worst_r = max(worst_r, abs(p - 1))
    r_max = worst_r

    # 2. log1p truncation over |r| <= r_max, measured RELATIVE to |r| because
    #    that is the accuracy the near-1 slots must deliver.
    trunc = mp.mpf(0)
    for i in range(1, 2001):
        r = -r_max + 2 * r_max * i / 2000
        if r == 0:
            continue
        series = mp.mpf(0)
        for n in range(POLY_LAST, 0, -1):
            series = r * (series + mp.mpf((-1) ** (n + 1)) / n)
        trunc = max(trunc, abs(mp.log(1 + r) - series) / abs(r))
    if trunc > mp.mpf(2) ** -70:
        fails.append(f"log1p truncation 2^{float(mp.log(trunc, 2)):.1f} relative, want <= 2^-70")

    # 3. L_j = -log(R_j) representable as a dd to ~2^-105 relative.
    worst_l = mp.mpf(0)
    for j in range(N):
        v = -mp.log(mp.mpf(rs[j]))
        if v == 0:
            continue
        hi = float(v)
        lo = float(v - mp.mpf(hi))
        worst_l = max(worst_l, abs((mp.mpf(hi) + mp.mpf(lo) - v) / v))
    if worst_l > mp.mpf(2) ** -104:
        fails.append(f"table dd error 2^{float(mp.log(worst_l, 2)):.1f}")

    # 4. The two slots adjacent to m = 1 must be exactly (R, L) = (1, 0), or
    #    the no-cancellation argument for x near 1 does not hold.
    if rs[0] != 1.0 or rs[N - 1] != 1.0:
        fails.append("slots 0 and 127 must have R_j == 1 exactly")

    return fails, dict(r_max=r_max, trunc=trunc, worst_l=worst_l)


def main():
    rs = [recip(j) for j in range(N)]
    fails, info = self_check(rs)
    for f in fails:
        print(f"self-check FAILED: {f}", file=sys.stderr)
    if fails:
        return 1
    print(
        f"self-check OK: |r| <= 2^{float(mp.log(info['r_max'], 2)):.2f}, "
        f"log1p truncation <= 2^{float(mp.log(info['trunc'], 2)):.1f} relative, "
        f"table dd <= 2^{float(mp.log(info['worst_l'], 2)):.1f}",
        file=sys.stderr,
    )

    ln2 = mp.log(2)
    ln2_hi = float(ln2)
    ln2_lo = float(ln2 - mp.mpf(ln2_hi))
    split = 1.0 + SPLIT_SLOT / N

    print("// Auto-generated log tables for corvus. DO NOT EDIT --")
    print("// regenerate with tools/gen_log_table.py.")
    print(f"// {N} slots on the top {N.bit_length() - 1} mantissa bits; m centred on 1 by halving")
    print(f"// when m0 >= {split} (i.e. j >= {SPLIT_SLOT}), so x near 1 needs no cancelling")
    print("// k*ln2 term. Slots 0 and 127 carry R_j = 1, L_j = 0 exactly.")
    print(f"static_assert(kLogN == {N}, \"regenerate src/log_dd_data.inc\");")
    print(f"static_assert(kLogSplitSlot == {SPLIT_SLOT}, \"regenerate src/log_dd_data.inc\");")
    print(f"const double kLogSplit = {as_hex(split)};   // 1 + {SPLIT_SLOT}/{N}")
    print(f"const double kLogLn2Hi = {as_hex(ln2_hi)};")
    print(f"const double kLogLn2Lo = {as_hex(ln2_lo)};")
    print(f"alignas(64) const double kLogTableR[{N}] = {{")
    for j in range(N):
        print(f"    {as_hex(rs[j])},")
    print("};")
    for name, part in (("Hi", 0), ("Lo", 1)):
        print(f"alignas(64) const double kLogTableL{name}[{N}] = {{")
        for j in range(N):
            v = -mp.log(mp.mpf(rs[j]))
            hi = float(v)
            print(f"    {as_hex(hi if part == 0 else float(v - mp.mpf(hi)))},")
        print("};")
    return 0


if __name__ == "__main__":
    sys.exit(main())
