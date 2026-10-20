#!/usr/bin/env python3
"""Generate src/bessel_data.h -- every table corvus::i0/i1/i0e/i1e need.

Two regimes, split at kBesselSplit, for BOTH nu=0 and nu=1, per the
frontier probe (PROBE-RECORD.md):

  SERIES  x in [0, kBesselSplit]: truncated power series in q = x^2/4,
      I0: sum q^k/(k!)^2 ; I1: (x/2) * sum q^k/(k!(k+1)!). All-positive,
      perfectly conditioned IN q -- but q itself is NOT well-conditioned
      to compute from x: I_nu's own logarithmic sensitivity to x is
      O(x) (elasticity x*I1(x)/I0(x) -> x for large x), so a single
      rounding of q=x*x/4 alone costs ~x ULP in the series result
      (measured directly: 2-7 ULP depending on split, independent of how
      many leading coefficients are promoted to dd -- NOT a Horner-
      accumulation artifact, a q-INPUT artifact). Fix: capture x*x
      EXACTLY via TwoProd (the ops::SquareLow/ProdLow exact-residual
      idiom NUMERICAL-DOCTRINE.md documents, same technique erfc.cpp
      uses for its own ssq/sl before ExpDdFrac), then apply a first-
      order derivative correction S(q_hi+q_lo) ~= S(q_hi)+S'(q_hi)*q_lo
      (S' cheap and plain-double -- the correction itself is already
      ~2^-53-relative, doesn't need its own extra precision). On TOP of
      that fix, the kBesselI0SeriesLead/kBesselI1SeriesLead lowest-degree
      (dominant-magnitude) coefficients are kept as dd pairs and folded
      via the SAME nested dd/double Horner pattern src/erfinv-inl.h's
      ErfinvCentralDd uses (L0+v*(L1+v*(L2+v*S(v))), leading terms dd,
      remainder one plain-double Horner pass) -- this generator pins the
      minimal depth (smallest that reaches the measured floor; beyond it
      more dd depth was measured to buy nothing further).

  TAIL  x >= kBesselSplit: clean-room Chebyshev refit (own nodes, own
      budget -- A&S/Cephes coefficients NOT consulted) of
      f_nu(t) = e^{-x} I_nu(x) sqrt(2*pi*x), t = 1/x, on [0, 1/kBesselSplit],
      converted to monomial-in-normalized-s Horner coefficients
      (gen_erfc_tail_poly.py's own pattern, reused structurally). Evaluated
      as plain-double Horner except the LAST step c0 + s*acc, which is dd
      (TwoProd + DdAddD, #52), so the poly is never rounded to a double.
      (t is flat/well-conditioned, no squaring hazard; this generator does
      not sweep deeper dd promotion for the tail.) Assembly: i_nu_e =
      DdMul(1/sqrt(2*pi*x), poly) with the reciprocal sqrt in dd (DdSqrt,
      DdRecipDd), rounded once; unscaled i_nu = that dd times exp_dd(+x)
      in ExpDdFrac's mantissa+exponent form, scaled by 2^e last.

  OVERFLOW BOUNDARIES: last-finite doubles for unscaled i0/i1, bisected
      independently against the round-to-inf threshold 2^1024*(1-2^-54)
      at dps=50; cross-checked bit-identical against PROBE-RECORD.md's
      independent dps=40 derivation (0x1.64fe5304e83e4p+9 /
      0x1.64fe69ff9fec7p+9) -- any future disagreement here is an
      ESCALATE, not something to resolve by picking one side.

Self-checks (mandatory; stderr budget lines; ANY miss -> exit nonzero,
emit nothing):
  (a) series truncation: measured tail fraction at q_max < FIT_TOL, both
      nu.
  (b) tail Chebyshev truncation: measured dense-sample relative error
      of the fitted monomial polynomial < FIT_TOL, both nu.
  (c) end-to-end replay: simulates the kernel arithmetic of
      src/bessel-inl.h (exact-q TwoProd + derivative correction + dd-lead
      Horner on the series side; Horner with a dd last step + DdMul by a
      dd 1/sqrt(2*pi*x), then the ExpDdFrac fold for unscaled, on the tail
      side) under fused AND unfused MulAdd, for BOTH the scaled and the
      unscaled assembly. exp_dd is modeled with an injected +/-2^-70
      relative perturbation, never assumed exact. Series sample: log-
      spaced coverage, bit-stepped refinement below the split, tiny-x/
      subnormal probes. Tail sample: 20k full-mantissa points in
      (kBesselSplit, min overflow boundary] (both assemblies), 1k log-
      uniform points above it up to DBL_MAX (scaled only), plus anchors
      (bit-stepped above the split and around both overflow boundaries).
      The replay measures the worst-case ULP floor
      (magnitude ULP distance against the correctly-rounded mpmath
      truth; NOT relative-error bits, which artificially blow up near
      the subnormal floor even when the double result IS the correctly
      rounded nearest value, which relative error would wrongly flag near
      the subnormal floor). Series and tail targets are both 1 ULP
      (ULP_TARGET_SERIES, ULP_TARGET_TAIL).
  (d) negative controls: corrupt one dd-lead coefficient's lo word and
      confirm the series floor measurement explodes (proves the dd
      machinery is load-bearing, not a no-op); add 4 ulp to the tail's c0
      and confirm the tail replay exceeds ULP_TARGET_TAIL (proves it
      resolves a few-ulp defect). Exit nonzero if either is NOT caught.
  (e) seam check: both routes evaluated at x=kBesselSplit, each within
      its own ULP target of the truth (they need not agree bit-for-bit;
      the kernel sends x=kBesselSplit to the series).
  (f) independent oracle cross-check: OWN high-dps series-sum (the
      literal mathematical definition, dps=80, NOT mpmath's internal
      besseli algorithm) vs mp.besseli at a handful of points, layered
      dps 40/80 -- agreement recorded, no oracle defect class found
      (unlike the beta family; this is erf-difficulty class per the
      binding design).

Usage:
    python3 tools/gen_bessel_data.py > src/bessel_data.h
"""
import math
import random
import struct
import sys
import time

import mpmath as mp
from mpmath import mpf

mp.mp.dps = 60

X_S = 8.0                    # kBesselSplit
FIT_TOL = mpf(2) ** -60      # both regimes' truncation target (design's stated budget)
ULP_TARGET_SERIES = 1.0
ULP_TARGET_TAIL = 1.0        # the published bound (#52: 2.0 hid a ~1.5 ULP kernel)
TAIL_FIT_REL_TARGET = 2.5e-16  # dense-verification gate for the RAW fitted
# polynomial's double-precision evaluation (erfc_tail's own REL_TARGET,
# gen_erfc_tail_poly.py's exact choice, reused): this is a ~1.5 ULP
# evaluation-rounding floor, NOT the mpmath-side truncation cut (FIT_TOL,
# checked separately at full precision before any double rounding enters).
KLEAD = {}                   # DERIVED during the replay self-check:
                             # smallest dd depth whose worst replay ULP over
                             # BOTH assemblies (unscaled bare hi+lo AND
                             # scaled i_nu_e) and BOTH MulAdd semantics
                             # (fused / unfused -- SSE4/SSSE3/SSE2 ship
                             # without FMA) meets ULP_TARGET_SERIES. An
                             # idealized sweep restricted to FMA-only,
                             # scaled-only assembly can hide a real ULP
                             # miss on the non-FMA, unscaled tiers.
TAIL_NODES = 50

T0 = time.time()


# ============================================================================
# hex/dd emission helpers
# ============================================================================
def hexf(x):
    return float.hex(float(x))


def round_to_dd(mpf_val):
    hi = float(mpf_val)
    lo = float(mpf_val - mp.mpf(hi))
    return hi, lo


def emit_array_1d(name, vals, ncols=6):
    print(f"inline constexpr double {name}[{len(vals)}] = {{")
    for i in range(0, len(vals), ncols):
        print("    " + ", ".join(hexf(v) for v in vals[i:i + ncols]) + ",")
    print("};")


# ============================================================================
# Part 0: exact-double primitives (mpmath at dps=60 models correctly-rounded
# FMA/mul/add/div/sqrt the same way ops:: delivers it on an FMA-capable
# target -- far more headroom than a*b+c for three doubles ever needs, so
# round(mpf) here IS what op::MulAdd computes).
# ============================================================================
DPS_EXACT = 60


def fma_d(a, b, c):
    return float(mp.mpf(a) * mp.mpf(b) + mp.mpf(c))


def mul_d(a, b):
    return float(mp.mpf(a) * mp.mpf(b))


def two_prod(a, b):
    p = mul_d(a, b)
    err = float(mp.mpf(a) * mp.mpf(b) - mp.mpf(p))
    return p, err


def horner_plain(coeffs_desc, x):
    acc = coeffs_desc[0]
    for c in coeffs_desc[1:]:
        acc = fma_d(acc, x, c)
    return acc


# ---- faithful dd/float-semantics primitives ----
# The primitives above model correctly-rounded FMA arithmetic -- true on the
# FMA-capable tiers (AVX2 and up) ONLY. SSE4/SSSE3/SSE2 are SHIPPING tiers
# without FMA: there op::MulAdd is an unfused mul-then-add (two roundings)
# and the dd layer's TwoProd takes ops-inl.h's Dekker path (exact in range,
# so two_prod above validly models BOTH). The kernel assemblies must
# therefore be replayed under BOTH semantics, with the dd algorithms
# mirrored from dd-inl.h op-for-op (TwoSum/Fast2Sum chains) rather than
# idealized exact-dd folds. Python float arithmetic IS IEEE double
# round-to-nearest, so plain expressions give unfused semantics directly.
def _mafold(a, b, c, fused):
    return fma_d(a, b, c) if fused else (a * b) + c


def two_sum_f(a, b):
    s = a + b
    bv = s - a
    return s, (a - (s - bv)) + (b - bv)


def fast2sum_f(a, b):
    s = a + b
    return s, (a - s) + b


def dd_add_f(ah, al, bh, bl):
    sh, sl = two_sum_f(ah, bh)
    th, tl = two_sum_f(al, bl)
    vh, vl = fast2sum_f(sh, sl + th)
    return fast2sum_f(vh, vl + tl)


def dd_addd_f(ah, al, b):
    sh, sl = two_sum_f(ah, b)
    return fast2sum_f(sh, sl + al)


def dd_muld_f(ah, al, b, fused):
    ph, pl = two_prod(ah, b)
    return fast2sum_f(ph, _mafold(al, b, pl, fused))


def dd_mul_f(ah, al, bh, bl, fused):
    ph, pl = two_prod(ah, bh)
    lo = _mafold(ah, bl, _mafold(al, bh, pl, fused), fused)
    return fast2sum_f(ph, lo)


def horner_plain_sem(coeffs_desc, x, fused):
    acc = coeffs_desc[0]
    for c in coeffs_desc[1:]:
        acc = _mafold(acc, x, c, fused)
    return acc


def horner_dd_lead(coeffs_asc_hi, coeffs_asc_lo, klead, x):
    """Mirrors src/erfinv-inl.h's ErfinvCentralDd: the lowest `klead`
    coefficients folded in dd, the remainder collapsed to one plain-
    double Horner pass first. Returns (hi, lo)."""
    n = len(coeffs_asc_hi)
    if klead == 0:
        return horner_plain(list(reversed(coeffs_asc_hi)), x), 0.0
    acc_hi, acc_lo = coeffs_asc_hi[klead - 1], coeffs_asc_lo[klead - 1]
    if klead < n:
        tail_desc = list(reversed(coeffs_asc_hi[klead:n]))
        s = horner_plain(tail_desc, x)
        exact = mp.mpf(acc_hi) + mp.mpf(acc_lo) + mp.mpf(mul_d(s, x))
        acc_hi = float(exact)
        acc_lo = float(exact - mp.mpf(acc_hi))
    for k in range(klead - 2, -1, -1):
        prod = (mp.mpf(acc_hi) + mp.mpf(acc_lo)) * mp.mpf(x)
        exact = mp.mpf(coeffs_asc_hi[k]) + mp.mpf(coeffs_asc_lo[k]) + prod
        acc_hi = float(exact)
        acc_lo = float(exact - mp.mpf(acc_hi))
    return acc_hi, acc_lo


def exp_dd_sim(x_mpf, rel_noise_ulp, sign):
    """Model exp_dd(x) (mantissa+exponent form): inject the documented
    ~2^-70 relative budget as a worst-case perturbation rather than
    assume the primitive exact."""
    true = mp.exp(mp.mpf(x_mpf))
    pert = true * (1 + sign * mp.mpf(rel_noise_ulp))
    e = int(mp.floor(mp.log(pert, 2))) + 1
    m = pert / mp.mpf(2) ** e
    hi = float(m)
    lo = float(m - mp.mpf(hi))
    return hi, lo, e


def ulp_distance_pos(true_val, got):
    """Magnitude ULP distance (valid for got, true_val >= 0: IEEE-754
    double bit patterns are monotonic in value there, covering subnormals
    uniformly -- unlike relative error, which blows up near the subnormal
    floor even for a correctly-rounded result)."""
    if not math.isfinite(got) or got < 0:
        return float("inf")
    ref = float(true_val)
    if got == ref:
        return 0.0
    b_ref = struct.unpack('<q', struct.pack('<d', ref))[0]
    b_got = struct.unpack('<q', struct.pack('<d', got))[0]
    return float(abs(b_got - b_ref))


def next_after(x, up=True):
    b = struct.unpack('<q', struct.pack('<d', x))[0]
    b += 1 if up else -1
    return struct.unpack('<d', struct.pack('<q', b))[0]


# ============================================================================
# Part 1: overflow boundary (independent re-derivation)
# ============================================================================
def find_overflow_boundary(nu, dps=50):
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        ovf = mpf(2) ** 1024 * (1 - mpf(2) ** -54)
        lo, hi = 700.0, 720.0
        while next_after(lo) < hi:
            mid = (lo + hi) / 2
            v = mp.besseli(nu, mpf(mid))
            if v < ovf:
                lo = mid
            else:
                hi = mid
        return lo
    finally:
        mp.mp.dps = old


# ============================================================================
# Part 2: series fit (raw reciprocal-factorial, truncated). Kept RAW rather
# than economized/refit: all-positive, exactly the mathematical series (no
# fit-generation risk), and the term count (~22) is already compact -- an
# economized Chebyshev-in-q refit was not pursued given the raw series
# already meets budget with a clean, trivially-verified construction.
# ============================================================================
def series_coeffs_raw(nu, kmax=60, dps=60):
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        coeffs = [mpf(1)]
        term = mpf(1)
        for k in range(1, kmax + 1):
            term = term / (k * k if nu == 0 else k * (k + 1))
            coeffs.append(term)
        return coeffs
    finally:
        mp.mp.dps = old


def series_truncate(coeffs, q_max, tol):
    old = mp.mp.dps
    mp.mp.dps = 60
    try:
        q_max = mpf(q_max)
        terms = [c * q_max ** k for k, c in enumerate(coeffs)]
        total = sum(terms)
        tail_after = [mpf(0)] * (len(terms) + 1)
        for k in range(len(terms) - 1, -1, -1):
            tail_after[k] = tail_after[k + 1] + terms[k]
        for N in range(len(terms)):
            frac = tail_after[N + 1] / total
            if frac < tol:
                return coeffs[: N + 1], float(frac)
        return coeffs[:], float(tail_after[0] / total)
    finally:
        mp.mp.dps = old


def deriv_coeffs(coeffs_asc):
    """d/dq of sum c_k q^k, ascending, index0..N-1."""
    return [k * c for k, c in enumerate(coeffs_asc)][1:]


# ============================================================================
# Part 3: tail fit (clean-room Chebyshev-in-t refit, own nodes/budget)
# ============================================================================
def cheb_coeffs_mp(f, a, b, n_nodes):
    a, b = mpf(a), mpf(b)
    c = (a + b) / 2
    h = (b - a) / 2
    nodes = [mp.cos(mp.pi * (j + mpf(1) / 2) / n_nodes) for j in range(n_nodes)]
    vals = [f(c + h * s) for s in nodes]
    coeffs = []
    for k in range(n_nodes):
        acc = mpf(0)
        for j in range(n_nodes):
            acc += vals[j] * mp.cos(mp.pi * k * (j + mpf(1) / 2) / n_nodes)
        a_k = 2 * acc / n_nodes
        coeffs.append(a_k / 2 if k == 0 else a_k)
    return coeffs, c, h


def cheb_truncate(coeffs, tol):
    deg = len(coeffs) - 1
    while deg > 0 and abs(coeffs[deg]) < tol:
        deg -= 1
    return coeffs[: deg + 1]


def cheb_to_monomial(coeffs):
    t_prev = [mpf(1)]
    t_cur = [mpf(0), mpf(1)]
    mono = [mpf(0)] * len(coeffs)

    def add(poly, w):
        for i, p in enumerate(poly):
            mono[i] += w * p

    add(t_prev, coeffs[0])
    if len(coeffs) > 1:
        add(t_cur, coeffs[1])
    for k in range(2, len(coeffs)):
        t_next = [mpf(0)] + [2 * x for x in t_cur]
        for i, x in enumerate(t_prev):
            t_next[i] -= x
        add(t_next, coeffs[k])
        t_prev, t_cur = t_cur, t_next
    return mono


def f_scaled(nu):
    def f(t):
        if t == 0:
            return mpf(1)
        x = 1 / t
        return mp.besseli(nu, x) * mp.exp(-x) * mp.sqrt(2 * mp.pi * x)
    return f


def build_tail_fit(nu, x_s, tol, n_nodes=TAIL_NODES):
    lo, hi = mpf(0), 1 / mpf(x_s)
    coeffs, c, h = cheb_coeffs_mp(f_scaled(nu), lo, hi, n_nodes)
    coeffs_t = cheb_truncate(coeffs, tol)
    mono = cheb_to_monomial(coeffs_t)
    scale = 1 / h
    shift = -c / h
    return mono, scale, shift


# ============================================================================
# Part 4: simulated kernel assembly (mirrors the planned C++ exactly enough
# to measure the real ULP floor -- see module docstring).
# ============================================================================
def series_assemble_sim(nu, coeffs_hi, coeffs_lo, dcoef_desc, klead, x,
                        exp_sign, fused):
    """Faithful replay of src/bessel-inl.h's series branch: BOTH assemblies
    (unscaled i_nu = bare hi+lo / (x/2)-scaled hi+lo; scaled i_nu_e via a
    faithful DdMul against the exp_dd pair) under fused OR unfused MulAdd
    semantics, with the dd algorithms mirrored from dd-inl.h op-for-op.
    Replaying only the scaled assembly under FMA-only, exact-arithmetic
    folds can hide a real ULP miss on the unscaled path at non-FMA tiers.
    Returns (unscaled, scaled)."""
    ssq, sl = two_prod(x, x)             # SquareLow: exact on BOTH paths
    q_hi, q_lo = ssq * 0.25, sl * 0.25   # exact power-of-two scales
    n = len(coeffs_hi)
    # BesselLeadTailScalar, faithfully
    if klead == 0:
        s_hi, s_lo = horner_plain_sem(list(reversed(coeffs_hi)), q_hi,
                                      fused), 0.0
    else:
        acc = (coeffs_hi[klead - 1], coeffs_lo[klead - 1])
        if klead < n:
            s = horner_plain_sem(list(reversed(coeffs_hi[klead:])), q_hi,
                                 fused)
            acc = dd_addd_f(acc[0], acc[1], mul_d(s, q_hi))
        for k in range(klead - 2, -1, -1):
            acc = dd_muld_f(acc[0], acc[1], q_hi, fused)
            acc = dd_add_f(acc[0], acc[1], coeffs_hi[k], coeffs_lo[k])
        s_hi, s_lo = acc
    # BesselSeriesS derivative correction
    if dcoef_desc:
        dS = horner_plain_sem(dcoef_desc, q_hi, fused)
        s_hi, s_lo = dd_addd_f(s_hi, s_lo, mul_d(dS, q_lo))
    # exp_dd(-x) as a dd VALUE (the series regime uses ExpDd directly, not
    # the mantissa+exponent form), documented ~2^-70 budget injected
    # worst-case in the given sign
    ev = mp.exp(-mp.mpf(x)) * (1 + exp_sign * mp.mpf(2) ** -70)
    eh = float(ev)
    el = float(ev - mp.mpf(eh))
    if nu == 0:
        unscaled = s_hi + s_lo
        sch, scl = dd_mul_f(s_hi, s_lo, eh, el, fused)
    else:
        mag_h, mag_l = dd_muld_f(s_hi, s_lo, x * 0.5, fused)
        unscaled = mag_h + mag_l
        sch, scl = dd_mul_f(mag_h, mag_l, eh, el, fused)
    return unscaled, sch + scl


def tail_ive_dd_sim(coeffs_desc, scale, shift, x, fused=True):
    """Mirrors BesselTailIveDd and returns its UNROUNDED dd result: t =
    fl(1/x), s = MulAdd(t, scale, shift), plain-double Horner down to c1,
    then the last step c0 + s*acc in dd (TwoProd + DdAddD, #52), then
    DdMul(recip, poly).

    MODEL ASSUMPTION: recip = 1/sqrt(2*pi*x) enters as the correctly rounded
    dd of the exact value. The kernel builds it from a dd 2*pi, an exact
    power-of-two prescale, DdSqrt and DdRecipDd, each ~2^-104 relative; that
    chain is not replayed op-for-op, on the argument that its error sits
    ~2^-50 ulp below the result's rounding. The final DdMul IS replayed
    (fused or unfused), since its lo-word MulAdds differ between tiers."""
    u_hi = 1.0 / x
    s = _mafold(u_hi, scale, shift, fused)
    acc = coeffs_desc[0]
    for c in coeffs_desc[1:-1]:
        acc = _mafold(acc, s, c, fused)
    ph, pl = two_prod(acc, s)
    poly_h, poly_l = dd_addd_f(ph, pl, coeffs_desc[-1])
    rh, rl = round_to_dd(1 / mp.sqrt(2 * mp.pi * mp.mpf(x)))
    return dd_mul_f(rh, rl, poly_h, poly_l, fused)


def tail_ive_sim(coeffs_desc, scale, shift, x, fused=True):
    """Scaled tail assembly: the kernel's Add(ive_tail.hi, ive_tail.lo)."""
    h, l = tail_ive_dd_sim(coeffs_desc, scale, shift, x, fused)
    return h + l


EXP_INV_L = float(128 / mp.log(2))   # exp_dd's N/ln2, N = 2^kExpLog2N = 128


def exp_dd_frac_sim(x, sign):
    """Model ExpDdFrac(x, 0) as (m_hi, m_lo, e) with exp(x) = m * 2^e.

    The exponent split mirrors the kernel (k = round(x*N/ln2), e = k >> 7,
    m ~ [1, 2)); with no over/underflow inside ScaleTwo the split cannot
    change the result, since every scaling by 2^e is exact.

    MODEL ASSUMPTION: m is the correctly rounded dd of exp(x)*(1 +
    sign*2^-70)/2^e. That is, exp_dd's whole error is a FIXED-SIGN relative
    perturbation at its documented budget bound (exp_dd-inl.h header,
    ~2^-70), replayed at both signs. The real table/polynomial error varies
    in sign and magnitude below that bound; it is not reproduced."""
    k = int(round(x * EXP_INV_L))
    e = k >> 7
    pert = mp.exp(mp.mpf(x)) * (1 + sign * mp.mpf(2) ** -70)
    return round_to_dd(pert / mp.mpf(2) ** e) + (e,)


def scale_two(v, e):
    """ScaleTwo: v * 2^floor(e/2) * 2^(e - floor(e/2)), two IEEE multiplies
    (Python float * overflows to inf exactly as the kernel's Mul does)."""
    e1 = e >> 1
    return (v * 2.0 ** e1) * 2.0 ** (e - e1)


def tail_iv_sim(ive_dd, x, exp_sign, fused):
    """Unscaled tail assembly of BesselNu0/BesselNu1: prod = DdMul(ive_tail,
    parts.m), each word scaled by ScaleTwo(., parts.e), then ONE Add -- the
    scaling is last, so a result near the overflow boundary rounds once."""
    mh, ml, e = exp_dd_frac_sim(x, exp_sign)
    ph, pl = dd_mul_f(ive_dd[0], ive_dd[1], mh, ml, fused)
    return scale_two(ph, e) + scale_two(pl, e)


def truth_ive(nu, x, dps=70):
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        xm = mpf(x)
        return mp.besseli(nu, xm) * mp.exp(-xm)
    finally:
        mp.mp.dps = old


def truth_iv_ive(nu, x, dps=70):
    """(I_nu(x), e^-x I_nu(x)) from one besseli evaluation at dps."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        xm = mpf(x)
        iv = mp.besseli(nu, xm)
        return iv, iv * mp.exp(-xm)
    finally:
        mp.mp.dps = old


# ============================================================================
# Part 5: sample construction (edge-refined, bit-stepped, per doctrine)
# ============================================================================
def bits_to_float(b):
    return struct.unpack('<d', struct.pack('<q', b))[0]


def struct_bits(x):
    return struct.unpack('<q', struct.pack('<d', x))[0]


def series_sample_points(x_s, n=300):
    pts = []
    for i in range(n):
        t = (i + 0.5) / n
        pts.append(x_s * math.exp(-10 * (1 - t)))
    b = struct_bits(float(x_s))
    for k in range(-200, 1):   # wide bracket: the
        pts.append(bits_to_float(b + k))   # unfused-noise worst sits ~13
        # ulps below the split; 200 gives deep margin either side of it
    pts += [1e-300, 1e-150, 1e-30, 1e-10, 1e-5, 5e-324, 2.0 ** -1074,
            2.0 ** -1022, 2.0 ** -1021]
    return sorted(set(p for p in pts if 0 < p <= x_s))


def tail_sample_points(x_s, n=300, seed=20260811):
    rng = random.Random(seed)
    pts = []
    for _ in range(n):
        u = rng.random()
        pts.append(x_s + u * 400)
    b = struct_bits(float(x_s))
    for k in range(0, 200):
        pts.append(bits_to_float(b + k))
    for nu_bound in (I0_BOUND_CACHE, I1_BOUND_CACHE):
        b2 = struct_bits(nu_bound)
        for k in range(-30, 30):
            pts.append(bits_to_float(b2 + k))
    return sorted(set(p for p in pts if p >= x_s and math.isfinite(p)))


def tail_replay_points(x_s, x_hi, n_value=10000, n_bits=10000, n_large=1000,
                       seed=20261020):
    """Full-mantissa tail replay sample (#61).

    dense: points in (x_s, x_hi], half uniform in value (low 20 mantissa
        bits re-drawn so every bit is random), half uniform in bit pattern
        (equal weight per binade, so the 8..64 end is not thinned out).
        x_hi is the smaller unscaled overflow boundary, so every dense point
        is valid for BOTH assemblies.
    large: scaled-only points above x_hi, uniform in bit pattern (log-
        uniform) up to DBL_MAX, plus DBL_MAX and its neighbours."""
    rng = random.Random(seed)
    b_lo, b_hi = struct_bits(x_s) + 1, struct_bits(x_hi)
    dense = set()
    while len(dense) < n_value:
        x0 = x_s + rng.random() * (x_hi - x_s)
        x = bits_to_float((struct_bits(x0) & ~0xFFFFF) | rng.getrandbits(20))
        if x_s < x <= x_hi:
            dense.add(x)
    while len(dense) < n_value + n_bits:
        dense.add(bits_to_float(rng.randint(b_lo, b_hi)))
    b_max = struct_bits(sys.float_info.max)
    large = set(bits_to_float(b_max - k) for k in range(4))
    while len(large) < n_large:
        large.add(bits_to_float(rng.randint(b_hi + 1, b_max)))
    return sorted(dense), sorted(large)


# ============================================================================
# main
# ============================================================================
def main():
    ok = True

    print(f"[gen_bessel_data] split x_s={X_S}, fit tol={float(FIT_TOL):.3e}, "
          f"klead=derived-by-sweep (both semantics, both assemblies)",
          file=sys.stderr)

    # ---- (overflow boundaries) ----
    global I0_BOUND_CACHE, I1_BOUND_CACHE
    i0_bound = find_overflow_boundary(0)
    i1_bound = find_overflow_boundary(1)
    PROBE_I0 = float.fromhex('0x1.64fe5304e83e4p+9')
    PROBE_I1 = float.fromhex('0x1.64fe69ff9fec7p+9')
    print(f"[gen_bessel_data] overflow boundary I0: {i0_bound!r} ({float.hex(i0_bound)}) "
          f"probe={float.hex(PROBE_I0)} match={i0_bound == PROBE_I0}", file=sys.stderr)
    print(f"[gen_bessel_data] overflow boundary I1: {i1_bound!r} ({float.hex(i1_bound)}) "
          f"probe={float.hex(PROBE_I1)} match={i1_bound == PROBE_I1}", file=sys.stderr)
    if i0_bound != PROBE_I0 or i1_bound != PROBE_I1:
        print("[gen_bessel_data] FAILED: overflow boundary disagrees with the frontier "
              "probe -- ESCALATE, do not pick one.", file=sys.stderr)
        return 1
    I0_BOUND_CACHE, I1_BOUND_CACHE = i0_bound, i1_bound

    # ---- (series fits) ----
    series = {}
    for nu in (0, 1):
        raw = series_coeffs_raw(nu)
        q_max = mpf(X_S) ** 2 / 4
        trunc, tail_frac = series_truncate(raw, q_max, FIT_TOL)
        n = len(trunc) - 1
        print(f"[gen_bessel_data] series nu={nu}: N={n} (terms 0..{n}), "
              f"measured tail fraction={tail_frac:.3e} (target {float(FIT_TOL):.1e})",
              file=sys.stderr)
        if tail_frac >= float(FIT_TOL):
            print(f"[gen_bessel_data] FAILED: series nu={nu} truncation misses budget",
                  file=sys.stderr)
            ok = False
        dcoeffs = deriv_coeffs(trunc)
        series[nu] = dict(coeffs=trunc, n=n, dcoeffs=dcoeffs, tail_frac=tail_frac)

    # ---- (tail fits) ----
    tail = {}
    for nu in (0, 1):
        mono, scale, shift = build_tail_fit(nu, X_S, FIT_TOL)
        n = len(mono) - 1
        mono_d = [float(c) for c in mono]
        # dense verification of the fitted polynomial itself (own budget,
        # gen_erfc_tail_poly.py precedent) before it ever enters the
        # end-to-end replay below.
        rng = random.Random(20260811 + nu)
        worst = 0.0
        for _ in range(4000):
            t = rng.uniform(0.0, 1.0 / X_S)
            s = t * float(scale) + float(shift)
            got = horner_plain(list(reversed(mono_d)), s)
            want = f_scaled(nu)(mpf(t)) if t > 0 else mpf(1)
            rel = abs((mpf(got) - want) / want) if want != 0 else abs(mpf(got))
            worst = max(worst, float(rel))
        print(f"[gen_bessel_data] tail nu={nu}: degree={n}, dense max rel err={worst:.3e} "
              f"(evaluation-rounding gate {TAIL_FIT_REL_TARGET:.1e}, erfc_tail's own)",
              file=sys.stderr)
        if worst > TAIL_FIT_REL_TARGET:
            print(f"[gen_bessel_data] FAILED: tail nu={nu} fit exceeds budget",
                  file=sys.stderr)
            ok = False
        tail[nu] = dict(mono=mono_d, scale=float(scale), shift=float(shift), n=n)

    if not ok:
        return 1

    # ---- (end-to-end replay self-check + klead derivation) ----
    print("[gen_bessel_data] replay: series regime "
          "(fused+unfused x unscaled+scaled; klead swept)", file=sys.stderr)
    for nu in (0, 1):
        coeffs = series[nu]["coeffs"]
        coeffs_hi = [round_to_dd(c)[0] for c in coeffs]
        coeffs_lo = [round_to_dd(c)[1] for c in coeffs]
        dcoef_desc = [float(c) for c in reversed(series[nu]["dcoeffs"])]
        pts = series_sample_points(X_S)
        truth_u = {}
        truth_e = {}
        for p in pts:
            te = truth_ive(nu, p)
            truth_e[p] = te
            truth_u[p] = te * mp.exp(mp.mpf(p))   # unscaled I_nu (same dps)
        chosen = None
        for klead in range(1, 13):
            worst = -1.0
            worst_pt = None
            worst_tag = None
            for p in pts:
                for fused in (True, False):
                    for sign in (1.0, -1.0):
                        un, sc = series_assemble_sim(
                            nu, coeffs_hi, coeffs_lo, dcoef_desc, klead, p,
                            sign, fused)
                        du = ulp_distance_pos(truth_u[p], un)
                        ds = ulp_distance_pos(truth_e[p], sc)
                        if du > worst:
                            worst, worst_pt = du, p
                            worst_tag = ("unscaled", fused)
                        if ds > worst:
                            worst, worst_pt = ds, p
                            worst_tag = ("scaled", fused)
            print(f"[gen_bessel_data]   nu={nu} klead={klead}: worst "
                  f"{worst:.3f} ULP @ x={worst_pt!r} "
                  f"[{worst_tag[0]}, {'fused' if worst_tag[1] else 'unfused'}]",
                  file=sys.stderr)
            if worst <= ULP_TARGET_SERIES:
                chosen = klead
                series[nu]["worst_ulp"] = worst
                break
        if chosen is None:
            print(f"[gen_bessel_data] FAILED: series nu={nu} misses "
                  f"{ULP_TARGET_SERIES} ULP at every klead <= 12",
                  file=sys.stderr)
            ok = False
        else:
            KLEAD[nu] = chosen
            print(f"[gen_bessel_data]   nu={nu}: klead PINNED at {chosen}",
                  file=sys.stderr)
    if not ok:
        return 1

    # Tail replay (#61): fused/unfused x scaled/unscaled. Unscaled points
    # stop at the nu's own overflow boundary (the kernel saturates past it);
    # the unscaled exp_dd fold is replayed at both noise signs.
    print("[gen_bessel_data] replay: tail regime "
          "(fused+unfused x unscaled+scaled)", file=sys.stderr)
    t_tail = time.time()
    dense, large = tail_replay_points(X_S, min(I0_BOUND_CACHE, I1_BOUND_CACHE))
    anchors = tail_sample_points(X_S)
    tail_pts = sorted(set(dense) | set(anchors) | set(large))
    for nu in (0, 1):
        mono_d = tail[nu]["mono"]
        scale, shift = tail[nu]["scale"], tail[nu]["shift"]
        coeffs_desc = list(reversed(mono_d))
        bound = I0_BOUND_CACHE if nu == 0 else I1_BOUND_CACHE
        worst = {}
        count = {"unscaled": 0, "scaled": 0}
        for p in tail_pts:
            unscaled_ok = p <= bound
            if unscaled_ok:
                truth_u, truth_e = truth_iv_ive(nu, p)
            else:
                truth_e = truth_ive(nu, p)
            count["scaled"] += 1
            count["unscaled"] += unscaled_ok
            for fused in (True, False):  # FMA tiers and SSE4/SSSE3/SSE2
                ive = tail_ive_dd_sim(coeffs_desc, scale, shift, p, fused)
                got = [("scaled", ulp_distance_pos(truth_e, ive[0] + ive[1]))]
                if unscaled_ok:
                    for sign in (1.0, -1.0):
                        un = tail_iv_sim(ive, p, sign, fused)
                        got.append(("unscaled", ulp_distance_pos(truth_u, un)))
                for asm, u in got:
                    key = (asm, fused)
                    if key not in worst or u > worst[key][0]:
                        worst[key] = (u, p)
        nu_worst = -1.0
        for asm in ("unscaled", "scaled"):
            for fused in (True, False):
                u, p = worst[(asm, fused)]
                nu_worst = max(nu_worst, u)
                print(f"[gen_bessel_data]   nu={nu} {asm:8s} "
                      f"{'fused  ' if fused else 'unfused'}: worst {u:.3f} ULP "
                      f"@ x={p!r} over {count[asm]} pts "
                      f"(target <= {ULP_TARGET_TAIL})", file=sys.stderr)
        if nu_worst > ULP_TARGET_TAIL:
            print(f"[gen_bessel_data] FAILED: tail nu={nu} exceeds ULP target",
                  file=sys.stderr)
            ok = False
        tail[nu]["worst_ulp"] = nu_worst
    print(f"[gen_bessel_data]   tail replay time {time.time()-t_tail:.0f}s",
          file=sys.stderr)

    # ---- (seam continuity) ----
    print("[gen_bessel_data] seam continuity at x=kBesselSplit", file=sys.stderr)
    for nu in (0, 1):
        coeffs = series[nu]["coeffs"]
        coeffs_hi = [round_to_dd(c)[0] for c in coeffs]
        coeffs_lo = [round_to_dd(c)[1] for c in coeffs]
        dcoef_desc = [float(c) for c in reversed(series[nu]["dcoeffs"])]
        _, ser_val = series_assemble_sim(nu, coeffs_hi, coeffs_lo, dcoef_desc,
                                         KLEAD[nu], X_S, 1.0, True)
        coeffs_desc = list(reversed(tail[nu]["mono"]))
        tail_val = tail_ive_sim(coeffs_desc, tail[nu]["scale"], tail[nu]["shift"], X_S)
        truth = truth_ive(nu, X_S)
        u_s = ulp_distance_pos(truth, ser_val)
        u_t = ulp_distance_pos(truth, tail_val)
        print(f"[gen_bessel_data]   nu={nu}: series={ser_val!r} ({u_s:.2f} ULP) "
              f"tail={tail_val!r} ({u_t:.2f} ULP)", file=sys.stderr)
        if u_s > ULP_TARGET_SERIES or u_t > ULP_TARGET_TAIL:
            print("[gen_bessel_data] FAILED: seam mismatch", file=sys.stderr)
            ok = False

    # ---- (negative control) ----
    print("[gen_bessel_data] negative control: corrupt one dd-lead coefficient",
          file=sys.stderr)
    for nu in (0, 1):
        coeffs = series[nu]["coeffs"]
        coeffs_hi = [round_to_dd(c)[0] for c in coeffs]
        coeffs_lo = [round_to_dd(c)[1] for c in coeffs]
        dcoef_desc = [float(c) for c in reversed(series[nu]["dcoeffs"])]
        klead = KLEAD[nu]
        idx = klead - 1 if klead > 0 else 0
        bad_lo = list(coeffs_lo)
        bad_lo[idx] = bad_lo[idx] + 1e-3
        worst_bad = -1.0
        for p in series_sample_points(X_S)[:300]:
            _, got = series_assemble_sim(nu, coeffs_hi, bad_lo, dcoef_desc,
                                         klead, p, 1.0, True)
            truth = truth_ive(nu, p)
            worst_bad = max(worst_bad, ulp_distance_pos(truth, got))
        caught = worst_bad > 100.0
        print(f"[gen_bessel_data]   nu={nu}: corrupted floor {worst_bad:.3e} ULP "
              f"({'CAUGHT' if caught else 'NOT CAUGHT'})", file=sys.stderr)
        if not caught:
            print("[gen_bessel_data] FAILED: negative control not caught -- the dd "
                  "machinery is not load-bearing as implemented", file=sys.stderr)
            ok = False

    # Tail negative control: c0 + TAIL_NEG_ULPS*ulp(c0) must push the tail
    # replay past ULP_TARGET_TAIL, or the replay cannot see a few-ulp
    # coefficient defect. Every 20th dense point, both semantics and both
    # assemblies (one noise sign), to keep it cheap.
    TAIL_NEG_ULPS = 4
    print(f"[gen_bessel_data] negative control: tail c0 + {TAIL_NEG_ULPS} ulp",
          file=sys.stderr)
    for nu in (0, 1):
        bad = list(reversed(tail[nu]["mono"]))
        bad[-1] += TAIL_NEG_ULPS * math.ulp(bad[-1])
        worst_bad = -1.0
        for p in dense[::20]:
            truth_u, truth_e = truth_iv_ive(nu, p)
            for fused in (True, False):
                ive = tail_ive_dd_sim(bad, tail[nu]["scale"], tail[nu]["shift"],
                                      p, fused)
                worst_bad = max(worst_bad,
                                ulp_distance_pos(truth_e, ive[0] + ive[1]),
                                ulp_distance_pos(truth_u,
                                                 tail_iv_sim(ive, p, 1.0, fused)))
        caught = worst_bad > ULP_TARGET_TAIL
        print(f"[gen_bessel_data]   nu={nu}: corrupted tail floor {worst_bad:.3f} "
              f"ULP ({'CAUGHT' if caught else 'NOT CAUGHT'})", file=sys.stderr)
        if not caught:
            print("[gen_bessel_data] FAILED: tail negative control not caught -- "
                  "the tail replay cannot resolve a few-ulp defect",
                  file=sys.stderr)
            ok = False

    # ---- (independent oracle cross-check) ----
    print("[gen_bessel_data] independent oracle cross-check (own high-dps series "
          "vs mp.besseli)", file=sys.stderr)

    def own_series_besseli(nu, x, dps):
        old = mp.mp.dps
        mp.mp.dps = dps
        try:
            x = mpf(x)
            q = x * x / 4
            term = mpf(1)
            s = mpf(1)
            k = 0
            tol = mpf(10) ** (-(dps - 5))
            while True:
                k += 1
                term *= q / (k * k if nu == 0 else k * (k + 1))
                s += term
                if term / s < tol and k > 4:
                    break
                if k > 3000:
                    break
            if nu == 1:
                s *= x / 2
            return s
        finally:
            mp.mp.dps = old

    for nu in (0, 1):
        for xv in (0.5, 3.0, 7.9, 20.0, 100.0):
            for dps in (40, 80):
                own = own_series_besseli(nu, xv, dps + 20)
                old = mp.mp.dps
                mp.mp.dps = dps
                try:
                    ref = mp.besseli(nu, mpf(xv))
                finally:
                    mp.mp.dps = old
                rel = abs((own - ref) / ref) if ref != 0 else abs(own)
                if rel > mpf(10) ** (-(dps - 8)):
                    print(f"[gen_bessel_data] FAILED: independent cross-check nu={nu} "
                          f"x={xv} dps={dps} rel={rel}", file=sys.stderr)
                    ok = False
    print("[gen_bessel_data]   independent cross-check: all points agree within the "
          "modeled dps floor", file=sys.stderr)

    if not ok:
        print("[gen_bessel_data] ABORTING: emitting nothing (self-check failure)",
              file=sys.stderr)
        return 1

    # ---- emit header ----
    print(f"[gen_bessel_data] all checks passed in {time.time()-T0:.0f}s -- emitting "
          f"src/bessel_data.h", file=sys.stderr)

    print("// Auto-generated by tools/gen_bessel_data.py. DO NOT EDIT.")
    print("// Fits and their error budgets are documented in the generator;")
    print("// the kernel that consumes them is src/bessel-inl.h.")
    print("//")
    print("// All fits below are functions of q=x*x/4 or t=1/|x| -- both EVEN")
    print("// in x -- so a single table serves both signs: the kernel takes")
    print("// ax=|x| into every formula here and reapplies sign only at the very")
    print("// end for the odd i1/i1e (CopySign, not implicit propagation -- the")
    print("// dd assembly's internal Fast2Sum can turn -0 into +0 under IEEE")
    print("// round-to-nearest, the exact erfinv(-0) hazard src/erfinv-inl.h's")
    print("// ErfinvCentral already documents and guards the same way).")
    print("#ifndef CORVUS_BESSEL_DATA_H_")
    print("#define CORVUS_BESSEL_DATA_H_")
    print()
    print("namespace corvus::detail {")
    print()
    print("// Series/tail regime split, pinned by replay (probe range [8,12]):")
    print("// x_s=8 is the ONLY candidate in range that reaches the 1 ULP series")
    print("// target without deep dd promotion -- x_s=10/12 were measured to")
    print("// plateau at 2-7 ULP regardless of dd-lead depth (larger x_s widens")
    print("// the series' own hump-shaped term profile near the domain edge,")
    print("// worsening Horner conditioning there faster than it shortens the")
    print("// tail fit).")
    print(f"inline constexpr double kBesselSplit = {hexf(X_S)};")
    print()
    print("// Overflow boundaries (unscaled i0/i1): last-finite double x,")
    print("// independently re-derived (bisected against the round-to-inf")
    print("// threshold 2^1024*(1-2^-54) at dps=50).")
    print(f"inline constexpr double kBesselI0OverflowX = {hexf(i0_bound)};")
    print(f"inline constexpr double kBesselI1OverflowX = {hexf(i1_bound)};")
    print()
    print("// SERIES regime, x in [0, kBesselSplit]. i0e = S0(q)*exp_dd(-x);")
    print("// i1e = (x/2)*S1(q)*exp_dd(-x), q = x*x/4 captured EXACTLY (TwoProd/")
    print("// SquareLow-style residual -- q's own rounding costs ~x ULP here,")
    print("// NOT attenuated the way erfc's flat-in-u tail is) with a first-order")
    print("// derivative correction S(q_hi+q_lo)~=S(q_hi)+S'(q_hi)*q_lo using the")
    print("// *DCoef arrays below (plain double, cheap: the correction itself is")
    print("// already ~2^-53-relative). The lowest kBesselI{0,1}SeriesLead")
    print("// coefficients are dd (Hi/Lo); the remainder is one plain-double")
    print("// Horner pass (erfinv-inl.h ErfinvCentralDd's nested dd/double")
    print("// pattern) -- *SeriesCoef holds degrees [Lead..N] ascending.")
    for nu in (0, 1):
        klead = KLEAD[nu]
        coeffs = series[nu]["coeffs"]
        hi = [round_to_dd(c)[0] for c in coeffs]
        lo = [round_to_dd(c)[1] for c in coeffs]
        dcoef = [float(c) for c in series[nu]["dcoeffs"]]
        tag = f"I{nu}"
        ntail = len(coeffs) - klead
        print()
        print(f"inline constexpr int kBessel{tag}SeriesLead = {klead};")
        print(f"inline constexpr int kBessel{tag}SeriesNCoef = {ntail};")
        emit_array_1d(f"kBessel{tag}SeriesLeadHi", hi[:klead])
        emit_array_1d(f"kBessel{tag}SeriesLeadLo", lo[:klead])
        emit_array_1d(f"kBessel{tag}SeriesCoef", hi[klead:])
        emit_array_1d(f"kBessel{tag}SeriesDCoef", dcoef)
    print()
    print("// TAIL regime, x >= kBesselSplit. i_nu_e(x) = f_nu(t)/sqrt(2*pi*x),")
    print("// t = 1/x, f_nu Chebyshev-refit (own nodes/budget) in normalized")
    print("// s = t*Scale + Shift, s in [-1,1]. Plain-double Horner down to c1,")
    print("// then the last step c0 + s*acc in dd (#52), times the dd")
    print("// 1/sqrt(2*pi*x) (DdSqrt, DdRecipDd), so the result rounds once.")
    for nu in (0, 1):
        tag = f"I{nu}"
        mono_d = tail[nu]["mono"]
        n = len(mono_d)
        print()
        print(f"inline constexpr int kBessel{tag}TailNCoef = {n};")
        print(f"inline constexpr double kBessel{tag}TailScale = {hexf(tail[nu]['scale'])};")
        print(f"inline constexpr double kBessel{tag}TailShift = {hexf(tail[nu]['shift'])};")
        emit_array_1d(f"kBessel{tag}TailCoef", mono_d)
    print()
    print("}  // namespace corvus::detail")
    print()
    print("#endif  // CORVUS_BESSEL_DATA_H_")
    return 0


if __name__ == "__main__":
    sys.exit(main())
