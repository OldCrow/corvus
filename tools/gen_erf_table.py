#!/usr/bin/env python3
"""Generate src/erf_data.inc -- lookup table for corvus's erf kernel.

Ported from libstats scripts/gen_neon_erf_table.py (same author). The method
was INDEPENDENTLY DERIVED from first principles -- local Taylor expansion
about uniform grid points, with Hermite-polynomial derivatives of exp(-x^2);
see libstats docs/NEON_ERF_DERIVATION.md. No third-party erf implementation
was consulted; mpmath is used only as a high-precision oracle for VALUES.

    erf(a) = erf(r) + S * ( d + c2 d^2 + ... + c5 d^5 ),   d = a - r,
    S = erf'(r) = (2/sqrt(pi)) e^{-r^2},
    c_k(r) = (-1)^{k-1} H_{k-1}(r) / k!.

Layout difference vs libstats: 3 doubles per entry {e_hi, s, e_lo} in a flat
array (r_j is reconstructed in the kernel from the grid index, exactly, since
h = 2^-8). Emits the saturation bound kErfAMax: the smallest double whose
erf rounds to 1.

Usage:
    python3 tools/gen_erf_table.py > src/erf_data.inc
"""

import random
import struct
import sys

import mpmath as mp

mp.mp.dps = 50

LOG2_INV_H = 8           # grid spacing h = 2^-8 = 1/256
GRID_MAX = 6.0           # table covers [0, 6]
N_TERMS = 8              # series length used by the kernel (c6-c8 for erfc relative accuracy)


def as_hex(x: float) -> str:
    return float.hex(x)


def as_bits(x: float) -> int:
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def from_bits(u: int) -> float:
    return struct.unpack("<d", struct.pack("<Q", u))[0]


def hermite(n_max):
    """Physicists' Hermite polynomials H_0..H_{n_max} as coefficient lists.

    H_0 = 1, H_1 = 2x, H_{n+1}(x) = 2x H_n(x) - 2n H_{n-1}(x).
    """
    H = [[mp.mpf(1)], [mp.mpf(0), mp.mpf(2)]]
    for n in range(1, n_max):
        two_x_Hn = [mp.mpf(0)] + [2 * c for c in H[n]]
        prev = H[n - 1] + [mp.mpf(0)] * (len(two_x_Hn) - len(H[n - 1]))
        H.append([a - 2 * n * b for a, b in zip(two_x_Hn, prev)])
    return H


def c_k(k, r, H):
    """c_k(r) = (-1)^{k-1} H_{k-1}(r) / k!."""
    coeffs = H[k - 1]
    val = sum(c * r ** i for i, c in enumerate(coeffs))
    return (-1) ** (k - 1) * val / mp.factorial(k)


def self_check():
    """Sum the N-term series against mpmath.erf at random (r, d) pairs.

    A wrong coefficient blows this up; a correct derivation lands at the
    truncation floor (<< 1 ULP for |d| <= h/2)."""
    H = hermite(N_TERMS)
    two_over_sqrtpi = 2 / mp.sqrt(mp.pi)
    half_h = mp.mpf(2) ** (-LOG2_INV_H) / 2
    rng = random.Random(20260720)
    worst = mp.mpf(0)
    for _ in range(400):
        r = mp.mpf(rng.uniform(0.0, GRID_MAX))
        d = mp.mpf(rng.uniform(-1.0, 1.0)) * half_h
        S = two_over_sqrtpi * mp.e ** (-r * r)
        series = d  # c_1 = 1
        for k in range(2, N_TERMS + 1):
            series += c_k(k, r, H) * d ** k
        approx = mp.erf(r) + S * series
        true = mp.erf(r + d)
        if abs(true) > 0:
            worst = max(worst, abs(approx - true) / abs(true))
    return float(worst)


def find_saturation():
    """Smallest double a with round-to-nearest(erf(a)) == 1.0, by bit bisection."""
    def rounded(x):
        return float(mp.erf(mp.mpf(x)))
    lo, hi = 5.0, 6.5
    assert rounded(lo) < 1.0 and rounded(hi) == 1.0
    blo, bhi = as_bits(lo), as_bits(hi)
    while bhi - blo > 1:
        mid = (blo + bhi) // 2
        if rounded(from_bits(mid)) == 1.0:
            bhi = mid
        else:
            blo = mid
    return from_bits(bhi)


def main():
    worst = self_check()
    if worst > 1e-15:
        print(f"self-check FAILED: worst relative series error {worst:g}", file=sys.stderr)
        return 1
    print(f"self-check OK: worst relative N={N_TERMS} series error {worst:.2e}", file=sys.stderr)

    a_max = find_saturation()
    two_over_sqrtpi = 2 / mp.sqrt(mp.pi)
    h = 2.0 ** (-LOG2_INV_H)
    n = int(GRID_MAX / h)  # last index (r_n = 6)

    print("// Auto-generated erf table for corvus. DO NOT EDIT --")
    print("// regenerate with tools/gen_erf_table.py.")
    print("// Method independently derived (local Taylor expansion about grid points,")
    print("// Hermite-polynomial derivatives of exp(-x^2)); see libstats")
    print("// docs/NEON_ERF_DERIVATION.md (same author). mpmath used only as a")
    print("// high-precision oracle for values.")
    print(f"// Grid r_j = j/{1 << LOG2_INV_H}, j = 0..{n} (covers [0, 6]). Flat layout, 3 doubles")
    print("// per entry: {E_hi, S, E_lo}. E_hi = RN(erf(r_j)); E_lo = RN(erf(r_j) - E_hi);")
    print("// S = RN((2/sqrt(pi)) exp(-r_j^2)). r_j is reconstructed in the kernel.")
    print(f"// kErfAMax = {a_max!r} (bits {as_bits(a_max):#018x}): the smallest double")
    print("// whose erf rounds to 1.0; |x| >= kErfAMax saturates to +/-1 (correctly rounded).")
    print(f"const double kErfAMax = {as_hex(a_max)};")
    print(f"const int64_t kErfTableLastIndex = {n};")
    print(f"alignas(64) const double kErfTable[{(n + 1) * 3}] = {{")
    for j in range(n + 1):
        r = float(j) * h  # exact: h = 2^-8, j <= 1536
        E = mp.erf(mp.mpf(r))
        e_hi = float(E)
        e_lo = float(E - mp.mpf(e_hi))
        s = float(two_over_sqrtpi * mp.e ** (-mp.mpf(r) ** 2))
        print(f"    {as_hex(e_hi)}, {as_hex(s)}, {as_hex(e_lo)},")
    print("};")
    return 0


if __name__ == "__main__":
    sys.exit(main())
