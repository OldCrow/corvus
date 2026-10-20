#!/usr/bin/env python3
"""Generate src/gamma_data.h -- every table the gamma_p/gamma_q kernel needs.

Three independent pieces:

  Temme table   c_0..c_10 (K=11), each a Chebyshev fit in eta over the band
                eta in [-sqrt(2*phi(1/2)), +sqrt(2*phi(2))], phi(l)=l-1-log(l).
                Clean-room extraction (NO recursion ported from any book or
                library): at each of 33 Chebyshev nodes in eta, solve a
                15x15 Vandermonde in v=1/a (samples a_j=512*2^j, j=0..14,
                dps=100) for R(a,eta), the scaled difference between the
                true regularized tail and its leading erfc term. THE ORACLE
                TRAP (see R_exact docstring): for eta<0 this MUST be
                extracted via the P-side identity, never the Q-side -- the
                Q-side is a 1-minus-tiny cancellation that silently returns
                garbage.
                Kernel evaluates each row by Clenshaw in eta, then Horner in
                r=1/a across the 11 rows -- so coefficients stay in
                CHEBYSHEV form (unlike lgamma_data.h's monomial tables).

  1/n table     dd pairs for n=1..36, R4's alternating-series weights.

  dd constants  2*pi and 1/sqrt(pi), for Temme's sqrt(2*pi*a) and the z.lo
                erfc correction.

  Tail G        gamma-private dd fit of G(z) = z e^{z^2} erfc(z) over the
                Temme tail's z range (#64), separate from erfc's own
                erfc_tail_data.h so erfc/erfinv bits stay put.

(The phi series for Log1pmxDd moved with the shared dd primitives to
src/dd_special-inl.h; it is emitted by tools/gen_dd_special_data.py,
self-checks included.)

Self-checks (mandatory, budget lines to stderr; ANY miss -> exit nonzero,
emit nothing):
  (a) disjoint-sample re-extraction (a_j=768*2^j) agrees with the primary
      extraction (a_j=512*2^j) to the stated per-coefficient tolerance.
  (b) cross-check against a saved probe JSON, when present.
  (c) end-to-end replay: EMITTED (float64, truncated+padded) coefficients,
      otherwise exact mpf arithmetic, vs the small-side oracle.
  (d) series length N<=64 (worst-case exact-arithmetic tail bound).
  (e) backward-CF depth N=kGammaCfN (worst-case boundary region; includes a
      permanent fine sweep over (1.5, 2.5] at x=a+1 -- a coarse discrete
      a-grid can miss this boundary's true supremum entirely, since it
      sits at a->1.5+ rather than at any conveniently round sample point).
  (f) R4 alternating-series length N=kGammaR4N.
  (g) tail G: op-for-op replay of GammaTailG (and the gt = G/(2z) product
      GammaTemme forms from it) under FUSED and UNFUSED MulAdd, over
      uniform and edge-bit-stepped arms, vs mpmath; budget 2^-63 relative.
      Negative control: dropping c0's low word must fail the same budget.
  (Former checks (g)/(h) covered the phi series and moved with it to
  gen_dd_special_data.py as its checks (a)/(b).)

Usage:
    python3 tools/gen_gamma_data.py > src/gamma_data.h
"""

import math
import sys
from fractions import Fraction

import mpmath as mp

mp.mp.dps = 100

# --- Design constants -------------------------------------------------------
KEXT = 15          # Vandermonde order: extracts c_0..c_14 per node.
K = 11             # table rows kept/emitted: c_0..c_10.
NNODES = 33        # Chebyshev nodes spanning the eta band.
A0_PRIMARY = 512   # primary sample grid a_j = A0 * 2^j, j=0..KEXT-1.
A0_CHECK = 768     # disjoint sample grid for self-check (a).
TAIL_CUT = mp.mpf(2) ** -60

SERIES_NMAX = 64
SERIES_TARGET = mp.mpf(2) ** -58
# CF_N=40 measured 2^-54.71 at the a->1.5+, x=a+1 boundary (self-check (e)).
# 44 gives ~2^-57-class margin there; four extra backward-CF divisions are
# cheap.
CF_N = 44
CF_TARGET = mp.mpf(2) ** -56
# R4_NMAX=30 measured 2^-47.7 at x=4 (minimum needed is ~34); 36 = 34 + 2
# margin.
R4_NMAX = 36
R4_TARGET = mp.mpf(2) ** -58
REPLAY_TARGET = mp.mpf(2) ** -56

# Gamma-private erfc tail factor G(z) = z e^{z^2} erfc(z) (#64; see the
# "Tail G" block below). One interval over the Temme tail's live z range,
# TAIL_DEG + 1 monomial coefficients in s = 1/z - C, the first TAIL_NDD of
# them dd and evaluated by dd Horner steps. TAIL_ZHI sits above
# sqrt(-kGammaExpFloor) = sqrt(800) ~ 28.28, the largest z of a live tail lane.
TAIL_ZLO = 6.0
TAIL_ZHI = 28.5
TAIL_DEG = 15
TAIL_NDD = 3
TAIL_NNODES = 96
TAIL_TARGET = mp.mpf(2) ** -63
EXP_FLOOR = -800.0

# --- Temme band + lambda(eta) ------------------------------------------------
def phi_lam(lam):
    """phi(lambda) = lambda - 1 - log(lambda); phi(1) = 0, phi >= 0."""
    return lam - 1 - mp.log(lam)


ETA_LO = -mp.sqrt(2 * phi_lam(mp.mpf("0.5")))
ETA_HI = mp.sqrt(2 * phi_lam(mp.mpf(2)))
ETA_MID = (ETA_HI + ETA_LO) / 2
ETA_HALF = (ETA_HI - ETA_LO) / 2


def lam_of_eta(eta):
    """Invert phi(lambda) = eta^2/2 by bisection -- monotone on each side of 1."""
    if eta == 0:
        return mp.mpf(1)
    target = eta * eta / 2
    if eta > 0:
        lo, hi = mp.mpf(1), mp.mpf(2)
        while hi - 1 - mp.log(hi) < target:
            hi *= 2
    else:
        lo, hi = mp.mpf("1e-30"), mp.mpf(1)
    for _ in range(400):
        mid = (lo + hi) / 2
        v = mid - 1 - mp.log(mid) - target
        if eta > 0:
            lo, hi = (mid, hi) if v < 0 else (lo, mid)
        else:
            lo, hi = (lo, mid) if v < 0 else (mid, hi)
    return (lo + hi) / 2


def R_exact(a, eta, lam):
    """R(a,eta) = (tail - leading-erfc-term) * sqrt(2*pi*a) * exp(a*phi).

    THE ORACLE TRAP: for eta >= 0 (lambda >= 1) both Q and erfc(.)/2 are
    same-scale tinies, so the difference is a safe relative cancellation.
    For eta < 0 (lambda < 1), that same difference computed from the Q side
    is ONE MINUS TWO TINIES -- it needs ~a*phi*log10(e) extra decimal
    digits to survive and silently returns garbage/inf otherwise. The fix
    is the P-side identity:
    R = (1/2*erfc(-eta*sqrt(a/2)) - P) * sqrt(2*pi*a) * exp(a*phi), with P
    the regularized LOWER incomplete gamma computed directly -- never
    derived as 1 - Q.
    """
    if eta >= 0:
        Q = mp.gammainc(a, lam * a, regularized=True)
        base = Q - mp.erfc(eta * mp.sqrt(a / 2)) / 2
    else:
        P = mp.gammainc(a, 0, lam * a, regularized=True)
        base = mp.erfc(-eta * mp.sqrt(a / 2)) / 2 - P
    return base * mp.sqrt(2 * mp.pi * a) * mp.exp(a * eta * eta / 2)


def extract_c(eta, lam, a0, kext=KEXT):
    """Solve the kext x kext Vandermonde in v=1/a for R's power-series coeffs."""
    A = mp.matrix(kext, kext)
    b = mp.matrix(kext, 1)
    for j in range(kext):
        a = mp.mpf(a0) * 2 ** j
        v = 1 / a
        for k in range(kext):
            A[j, k] = v ** k
        b[j] = R_exact(a, eta, lam)
    return mp.lu_solve(A, b)


def cheb_coeffs_from_vals(vals):
    """Standard DCT-II Chebyshev coefficients from values at the NNODES nodes.

    A list comprehension, not a generator, feeds fsum: mpmath's fsum
    mis-sums a generator argument (gen_lgamma_data.py's own documented trap).
    """
    n = len(vals)
    out = []
    for j in range(n):
        s = mp.fsum([vals[i] * mp.cos(j * mp.pi * (2 * i + 1) / (2 * n))
                      for i in range(n)])
        out.append(s * 2 / n if j else s / n)
    return out


def clenshaw(coefs, t):
    """Clenshaw evaluation of a Chebyshev series (coefs[0] full weight)."""
    b1 = b2 = mp.mpf(0)
    for j in range(len(coefs) - 1, 0, -1):
        b1, b2 = 2 * t * b1 - b2 + coefs[j], b1
    return t * b1 - b2 + coefs[0]


# --- doubles ------------------------------------------------------------------
def rd(x):
    return float(x)


def dd_split(x):
    hi = rd(x)
    return hi, rd(mp.mpf(x) - mp.mpf(hi))


def hexf(x):
    return float.hex(float(x))


# --- extraction driver ---------------------------------------------------------
def extract_all_nodes(a0):
    """Return (nodes, cvals): nodes[i]=(t_i, eta_i, lam_i), cvals[i]=len-KEXT mpf vec."""
    nodes = []
    cvals = []
    for i in range(NNODES):
        t = mp.cos(mp.pi * (2 * i + 1) / (2 * NNODES))
        eta = ETA_MID + ETA_HALF * t
        lam = lam_of_eta(eta)
        nodes.append((t, eta, lam))
        cvals.append(extract_c(eta, lam, a0))
    return nodes, cvals


def cheb_fits_from_cvals(cvals, k_rows=K):
    """Per-row (k=0..k_rows-1) Chebyshev coefficients over the NNODES nodes."""
    fits = []
    for k in range(k_rows):
        vals = [cvals[i][k] for i in range(NNODES)]
        fits.append(cheb_coeffs_from_vals(vals))
    return fits


def truncate(coeffs):
    d = len(coeffs) - 1
    while d > 0 and abs(coeffs[d]) < TAIL_CUT:
        d -= 1
    return coeffs[: d + 1]


# --- replay side: R1 series / R2 backward CF / R4 sum / phi series ------------
def series_len_bound(a, x, nmax=SERIES_NMAX, target=SERIES_TARGET):
    a, x = mp.mpf(a), mp.mpf(x)
    t = mp.mpf(1)
    s = mp.mpf(1)
    for n in range(1, nmax + 1):
        t *= x / (a + n)
        s += t
        r = x / (a + n + 1)
        if r < 1:
            tail = t * r / (1 - r)
            if tail < target * s:
                return n
    return None


def cf_backward_val(a, x, n=CF_N):
    a, x = mp.mpf(a), mp.mpf(x)
    kk = x + 2 * n + 1 - a
    for j in range(n, 0, -1):
        kk = (x + 2 * j - 1 - a) - j * (j - a) / kk
    return 1 / kk


def r4_series_contribution(a, x, nterms):
    """Partial sum of sum_{n>=1} (-x)^n/(n! (a+n)), dd terms modelled exactly."""
    a, x = mp.mpf(a), mp.mpf(x)
    term = mp.mpf(1)
    s = mp.mpf(0)
    for n in range(1, nterms + 1):
        term *= -x / n
        s += term / (a + n)
    return s


# --- Tail G: gamma-private dd erfc tail factor (#64) --------------------------
# GammaTemme's tail (a*phi in (36, 800], z = sqrt(a*phi) in (6, 28.29]) used
# erfc's plain-double G from erfc_tail_data.h: ~1.7 x 2^-53 relative (fit
# 0.81 + Horner 0.89), entering gamma_q at weight up to (lambda-1)/eta = 1.28
# near lambda = 2 -- worth 2.1 of the 2.9 ulp at #64's witness. This is a
# separate fit so erfc's and erfinv's bits do not move.
#
# Variable: u = 1/z (erfc's choice -- G is a function of u^2 analytic off the
# imaginary u axis, so on [1/28.5, 1/6] the nearest trouble is u = 0), used
# UNSCALED as s = u - C with C = the interval midpoint rounded to double: the
# kernel forms s from the dd u = 2 * DdRecip(2z) by one DdAddD, so s is a dd
# with no scale-constant rounding at all. Chebyshev interpolation on
# TAIL_NNODES nodes, truncated at degree TAIL_DEG (|a_16|/a_0 ~ 2e-21), then
# re-expanded exactly into monomials in s; coefficients 0..TAIL_NDD-1 are dd,
# the rest doubles.
#
# Why 3 dd steps (measured with this replay, both MulAdd semantics, deg 15):
# 1 dd step (c0 dd only) 2^-58.6; 2 steps 2^-60.9 (floor set by the double
# step-2 rounding, |s q_2| / G ~ 2^-10); 3 steps 2^-65.4; 4 steps 2^-67.3.
# Degree at 3 steps: 13 -> 2^-59.3, 14 -> 2^-62.8, 15 -> 2^-65.4, 16 -> same.
# One interval at degree 15 beats erfc's three at degree 11: no per-lane
# coefficient selects at all.
def tail_g_exact(z):
    z = mp.mpf(z)
    return z * mp.exp(z * z) * mp.erfc(z)


def tail_g_fit():
    """(C, hi[TAIL_DEG+1], lo[TAIL_NDD], cheb) for the emitted table."""
    ulo, uhi = 1 / mp.mpf(TAIL_ZHI), 1 / mp.mpf(TAIL_ZLO)
    c, h = (ulo + uhi) / 2, (uhi - ulo) / 2
    n = TAIL_NNODES
    th = [mp.pi * (j + mp.mpf(1) / 2) / n for j in range(n)]
    vals = [tail_g_exact(1 / (c + h * mp.cos(t))) for t in th]
    cheb = []
    for k in range(TAIL_DEG + 1):
        acc = mp.fsum(vals[j] * mp.cos(k * th[j]) for j in range(n))
        cheb.append(acc / n if k == 0 else 2 * acc / n)
    # Chebyshev -> monomial in t = (u - c)/h.
    tp, tc = [mp.mpf(1)], [mp.mpf(0), mp.mpf(1)]
    mono_t = [mp.mpf(0)] * len(cheb)
    for i, p in enumerate(tp):
        mono_t[i] += cheb[0] * p
    for i, p in enumerate(tc):
        mono_t[i] += cheb[1] * p
    for k in range(2, len(cheb)):
        tn = [mp.mpf(0)] + [2 * x for x in tc]
        for i, x in enumerate(tp):
            tn[i] -= x
        for i, p in enumerate(tn):
            mono_t[i] += cheb[k] * p
        tp, tc = tc, tn
    # t = alpha*s + beta, s = u - C exactly: re-expand in s.
    cd = rd(c)
    alpha, beta = 1 / h, (mp.mpf(cd) - c) / h
    mono_s = [mp.mpf(0)] * len(mono_t)
    for k, m in enumerate(mono_t):
        for j in range(k + 1):
            mono_s[j] += m * mp.binomial(k, j) * alpha ** j * beta ** (k - j)
    hi = [rd(m) for m in mono_s]
    lo = [rd(mono_s[k] - mp.mpf(hi[k])) for k in range(TAIL_NDD)]
    return cd, hi, lo, cheb


# Op-for-op replay of GammaTailG and its caller's DdRecip / DdMul, in IEEE
# double (python floats) with MulAdd modelled FUSED (one rounding, exact
# rationals) or UNFUSED (mul then add -- Highway's emulation on SSE2/SSSE3/
# SSE4). ProdLow is exact on both (FMA residual or Dekker split; no operand
# here is near the 2^996 split ceiling or the 2^-968 residual floor).
def _fma(a, b, c, fused):
    if fused:
        return float(Fraction(a) * Fraction(b) + Fraction(c))
    return (a * b) + c


def _prodlow(a, b, p):
    return float(Fraction(a) * Fraction(b) - Fraction(p))


def _two_sum(a, b):
    s = a + b
    bv = s - a
    return s, (a - (s - bv)) + (b - bv)


def _fast2sum(a, b):
    s = a + b
    return s, b - (s - a)


def _two_prod(a, b):
    p = a * b
    return p, _prodlow(a, b, p)


def _dd_add(a, b):
    s = _two_sum(a[0], b[0])
    t = _two_sum(a[1], b[1])
    v = _fast2sum(s[0], s[1] + t[0])
    return _fast2sum(v[0], v[1] + t[1])


def _dd_add_d(a, b):
    s = _two_sum(a[0], b)
    return _fast2sum(s[0], s[1] + a[1])


def _dd_mul(a, b, fused):
    p = _two_prod(a[0], b[0])
    lo = _fma(a[0], b[1], _fma(a[1], b[0], p[1], fused), fused)
    return _fast2sum(p[0], lo)


def _dd_recip(a):
    q = 1.0 / a
    p = a * q
    rem = (1.0 - p) - _prodlow(a, q, p)
    return _fast2sum(q, rem * q)


def tail_g_replay(z, tab, fused):
    """(G, gt) as the kernel computes them at the double z: G dd, and the
    gt = G/(2z) dd GammaTemme consumes."""
    cd, hi, lo = tab
    zc = min(max(z, TAIL_ZLO), TAIL_ZHI)
    r = _dd_recip(zc + zc)
    u = (2.0 * r[0], 2.0 * r[1])
    s = _dd_add_d(u, -cd)
    q = hi[TAIL_DEG]
    for k in range(TAIL_DEG - 1, TAIL_NDD - 1, -1):
        q = _fma(q, s[0], hi[k], fused)
    k = TAIL_NDD - 1
    g = _dd_add((hi[k], lo[k]), _two_prod(s[0], q))
    for k in range(TAIL_NDD - 2, -1, -1):
        g = _dd_add((hi[k], lo[k]), _dd_mul(s, g, fused))
    return g, _dd_mul(g, r, fused)


def tail_g_points():
    """Per-arm sample: uniform in z and in u, plus bit-stepped neighbourhoods
    of the domain edges (6, sqrt(800), 28.5) and of u = C (s ~ 0, where
    DdAddD's Fast2Sum ordering is tightest)."""
    import random
    import struct

    def bits(x):
        return struct.unpack("<Q", struct.pack("<d", x))[0]

    def unbits(b):
        return struct.unpack("<d", struct.pack("<Q", b))[0]

    rng = random.Random(6464)
    zlive = math.sqrt(-EXP_FLOOR)
    arms = {
        "uniform-z": [rng.uniform(TAIL_ZLO, zlive) for _ in range(2000)],
        "uniform-u": [1.0 / rng.uniform(1.0 / zlive, 1.0 / TAIL_ZLO)
                      for _ in range(2000)],
    }
    for name, e in (("edge-6", TAIL_ZLO), ("edge-sqrt800", zlive),
                    ("edge-28.5", TAIL_ZHI)):
        arms[name] = [unbits(bits(e) + k) for k in range(-96, 97)]
    cd = rd((1 / mp.mpf(TAIL_ZHI) + 1 / mp.mpf(TAIL_ZLO)) / 2)
    arms["s~0"] = [1.0 / unbits(bits(cd) + k) for k in range(-96, 97)]
    for name in arms:
        arms[name] = [z for z in arms[name] if TAIL_ZLO <= z <= TAIL_ZHI]
    return arms


def tail_g_worst(tab, arms):
    """Worst relative error of G and of gt over every arm, both semantics."""
    worst = {}
    for fused in (True, False):
        wg = wt = mp.mpf(0)
        at = None
        for pts in arms.values():
            for z in pts:
                g, gt = tail_g_replay(z, tab, fused)
                want = tail_g_exact(z)
                eg = abs((mp.mpf(g[0]) + mp.mpf(g[1]) - want) / want)
                wgt = want / (2 * mp.mpf(z))
                et = abs((mp.mpf(gt[0]) + mp.mpf(gt[1]) - wgt) / wgt)
                if eg > wg:
                    wg, at = eg, z
                wt = max(wt, et)
        worst[fused] = (wg, wt, at)
    return worst


# --- self-checks ---------------------------------------------------------------
def check_a_disjoint_extraction(cvals_primary, cvals_check):
    print("(a) disjoint-sample re-extraction (a0=512 vs a0=768):", file=sys.stderr)
    rc = 0
    for k in range(K):
        worst = mp.mpf(0)
        for i in range(NNODES):
            av, bv = cvals_primary[i][k], cvals_check[i][k]
            if av != 0:
                worst = max(worst, abs((av - bv) / av))
        thresh = mp.mpf("1e-15") if k <= 5 else mp.mpf("1e-10")
        status = "OK" if worst <= thresh else "FAIL"
        print(f"    c_{k:2d}: worst rel diff {float(worst):.3e} "
              f"(threshold {float(thresh):.1e}) {status}", file=sys.stderr)
        if worst > thresh:
            rc = 1
    return rc


# The former "(b) cross-check against seed JSON" compared the fits against a
# temme_chebs.json from the original design session's machine-local
# scratchpad; the file never shipped, so the check had been a silent
# warn-and-skip no-op since commit (#17). Removed rather than repointed: the
# surviving checks re-derive the budget from scratch each run, which is the
# stronger guarantee.
def emit_hex_array_1d(name, vals):
    print(f"inline constexpr double {name}[{len(vals)}] = {{")
    print("    " + ", ".join(hexf(v) for v in vals) + ",")
    print("};")


def emit_hex_array_2d(name, rows):
    ncols = len(rows[0])
    print(f"inline constexpr double {name}[{len(rows)}][{ncols}] = {{")
    for row in rows:
        print("    {" + ", ".join(hexf(v) for v in row) + "},")
    print("};")


def main():
    rc = 0

    print(f"eta band: [{float(ETA_LO):.6f}, {float(ETA_HI):.6f}] "
          f"mid={float(ETA_MID):.6f} half={float(ETA_HALF):.6f}", file=sys.stderr)

    # --- Temme extraction, primary grid -------------------------------------
    print(f"extracting Temme coefficients: {NNODES} nodes x {KEXT} samples "
          f"(a0={A0_PRIMARY}), dps={mp.mp.dps} ...", file=sys.stderr)
    _, cvals_primary = extract_all_nodes(A0_PRIMARY)
    fits_primary = cheb_fits_from_cvals(cvals_primary, K)

    # --- self-check (a): disjoint re-extraction -----------------------------
    print(f"re-extracting with disjoint grid (a0={A0_CHECK}) for self-check (a) ...",
          file=sys.stderr)
    _, cvals_check = extract_all_nodes(A0_CHECK)
    rc |= check_a_disjoint_extraction(cvals_primary, cvals_check)

    if rc:
        print("Self-check (a) failed -- aborting before emission.",
              file=sys.stderr)
        return rc

    # --- truncate + pad rows -------------------------------------------------
    truncated = [truncate(row) for row in fits_primary]
    maxdeg = max(len(row) - 1 for row in truncated)
    emitted_rows = [[rd(c) for c in row] + [0.0] * (maxdeg - (len(row) - 1))
                    for row in truncated]
    for k, row in enumerate(truncated):
        print(f"  c_{k:2d}: truncated degree {len(row) - 1}", file=sys.stderr)
    print(f"padded table: [{K}][{maxdeg + 1}]", file=sys.stderr)

    eta_mid_d, eta_half_d = rd(ETA_MID), rd(ETA_HALF)

    # --- self-check (c): end-to-end replay -----------------------------------
    print("(c) end-to-end truncation+fit replay:", file=sys.stderr)
    worst_c = mp.mpf(0)
    worst_c_at = None
    lams = [mp.mpf(s) for s in ("0.5", "0.6", "0.75", "0.9", "0.97")]
    lams += [1 - mp.mpf(2) ** -20, 1 + mp.mpf(2) ** -20]
    lams += [mp.mpf(s) for s in ("1.03", "1.12", "1.3", "1.6", "1.85", "2")]
    a_vals = [mp.mpf(s) for s in ("20", "27.4", "40", "80", "320", "5120", "1.31e6")]
    rows_mpf = [[mp.mpf(c) for c in row] for row in emitted_rows]
    mid_m, half_m = mp.mpf(eta_mid_d), mp.mpf(eta_half_d)
    for lam in lams:
        phi = phi_lam(lam)
        sign = 1 if lam >= 1 else -1
        eta = sign * mp.sqrt(2 * phi)
        t = (eta - mid_m) / half_m
        ck = [clenshaw(row, t) for row in rows_mpf]
        for a in a_vals:
            aphi = a * phi
            if aphi > 800:
                continue
            S = mp.mpf(0)
            for k in range(K - 1, -1, -1):
                S = S / a + ck[k]
            R = mp.exp(-aphi) / mp.sqrt(2 * mp.pi * a) * S
            z = eta * mp.sqrt(a / 2)
            if lam >= 1:
                direct = mp.erfc(z) / 2 + R
            else:
                direct = mp.erfc(-z) / 2 - R
            if a <= mp.mpf("1e4"):
                oracle = (mp.gammainc(a, lam * a, regularized=True) if lam >= 1
                          else mp.gammainc(a, 0, lam * a, regularized=True))
            else:
                # Full-degree (untruncated) fit at dps=100 -- 1e-40-class
                # residual, far below the 2^-56 replay budget being measured.
                Sf = mp.mpf(0)
                fits_full = fits_primary  # untruncated Chebyshev coefficients
                ckf = [clenshaw(fits_full[k], t) for k in range(K)]
                for k in range(K - 1, -1, -1):
                    Sf = Sf / a + ckf[k]
                Rf = mp.exp(-aphi) / mp.sqrt(2 * mp.pi * a) * Sf
                oracle = (mp.erfc(z) / 2 + Rf) if lam >= 1 else (mp.erfc(-z) / 2 - Rf)
            if oracle == 0:
                continue
            rel = abs((direct - oracle) / oracle)
            if rel > worst_c:
                worst_c, worst_c_at = rel, (float(lam), float(a))
    print(f"    worst rel err {float(worst_c):.3e} (2^{float(mp.log(worst_c, 2)):.2f}) "
          f"at lam={worst_c_at[0]:.6g} a={worst_c_at[1]:.6g}, "
          f"target 2^-56", file=sys.stderr)
    if worst_c > REPLAY_TARGET:
        print("    FAILED: exceeds 2^-56", file=sys.stderr)
        rc = 1

    # --- self-check (d): series length ---------------------------------------
    print("(d) series fixed length N<=64:", file=sys.stderr)
    worst_d, worst_d_at = 0, None
    for i in range(200):
        a = mp.mpf(19) + mp.mpf(i) / 200  # dense [19,20)
        n = series_len_bound(a, a + 1)
        if n is None:
            n = SERIES_NMAX + 1
        if n > worst_d:
            worst_d, worst_d_at = n, ("a in [19,20)", float(a))
    for a in (mp.mpf("20"), mp.mpf("1e3"), mp.mpf("1e8"), mp.mpf("1e100"), mp.mpf("1e300")):
        n = series_len_bound(a, a / 2)
        if n is None:
            n = SERIES_NMAX + 1
        if n > worst_d:
            worst_d, worst_d_at = n, ("lam=1/2", float(a))
    print(f"    worst N={worst_d} at {worst_d_at}, cap {SERIES_NMAX}", file=sys.stderr)
    if worst_d > 62:
        print(f"    FAILED: worst N={worst_d} exceeds budget 62", file=sys.stderr)
        rc = 1

    # --- self-check (e): backward CF depth ------------------------------------
    print(f"(e) backward CF depth N={CF_N}:", file=sys.stderr)
    worst_e, worst_e_at = mp.mpf(0), None

    def cf_check(a, x):
        nonlocal worst_e, worst_e_at
        a, x = mp.mpf(a), mp.mpf(x)
        exact = mp.gammainc(a, x) * mp.exp(x) * x ** (-a)
        if exact == 0:
            return
        got = cf_backward_val(a, x)
        rel = abs((got - exact) / exact)
        if rel > worst_e:
            worst_e, worst_e_at = rel, (float(a), float(x))

    for a in (1.6, 2.0, 5.0, 10.0, 15.0, 19.0, 19.999, 20.0):
        x0 = a + 1.0
        for x in (x0, math.nextafter(x0, math.inf), math.nextafter(x0, -math.inf)):
            cf_check(a, x)
    # Permanent fine sweep over (1.5, 2.5] at x=a+1: the true supremum of
    # the {a in (3/2,20], x=a+1} case sits at a->1.5+, not at any of the
    # discrete a values in the list above (at N=40 it measures 2^-54.71
    # there, worse than the discrete list's a=1.6 sample suggests), so a
    # coarse discrete grid alone cannot be trusted to find it.
    for i in range(1, 201):
        a = 1.5 + i * (2.5 - 1.5) / 200.0
        x0 = a + 1.0
        for x in (x0, math.nextafter(x0, math.inf), math.nextafter(x0, -math.inf)):
            cf_check(a, x)
    for a in ("1e-300", "1e-8", "0.1", "0.5", "1.0", "1.5"):
        af = float(a)
        for x in (4.0, math.nextafter(4.0, math.inf), math.nextafter(4.0, -math.inf)):
            cf_check(af, x)
    for a in (20.0, 100.0, 1000.0, 2600.0):
        cf_check(a, 2 * a)
    print(f"    worst rel err {float(worst_e):.3e} (2^{float(mp.log(worst_e, 2)):.2f}) "
          f"at (a,x)={worst_e_at}, target 2^-56", file=sys.stderr)
    if worst_e > CF_TARGET:
        print("    FAILED: exceeds 2^-56", file=sys.stderr)
        rc = 1

    # --- self-check (f): R4 series length -------------------------------------
    print(f"(f) R4 alternating-series length N={R4_NMAX}:", file=sys.stderr)
    worst_f, worst_f_at = mp.mpf(0), None
    for a in ("1e-300", "0.5", "1.5"):
        af, x = mp.mpf(a), mp.mpf(4)
        exact_sum = r4_series_contribution(af, x, 400)
        partial_sum = r4_series_contribution(af, x, R4_NMAX)
        remainder = exact_sum - partial_sum
        gamma_ax = mp.gammainc(af, x)
        contribution = abs(x ** af * remainder)
        rel = contribution / abs(gamma_ax)
        if rel > worst_f:
            worst_f, worst_f_at = rel, float(a)
    print(f"    worst rel contribution {float(worst_f):.3e} "
          f"(2^{float(mp.log(worst_f, 2)):.2f}) at a={worst_f_at}, target 2^-58",
          file=sys.stderr)
    if worst_f > R4_TARGET:
        print("    FAILED: exceeds 2^-58", file=sys.stderr)
        rc = 1

    # --- self-check (g): tail G replay, fused AND unfused, + negative control -
    print(f"(g) tail G: z in [{TAIL_ZLO}, {TAIL_ZHI}], degree {TAIL_DEG}, "
          f"{TAIL_NDD} dd steps, {TAIL_NNODES} nodes:", file=sys.stderr)
    tail_c, tail_hi, tail_lo, tail_cheb = tail_g_fit()
    tail_tab = (tail_c, tail_hi, tail_lo)
    arms = tail_g_points()
    print("    arms: " + ", ".join(f"{k} {len(v)}" for k, v in arms.items())
          + f"; |a_{TAIL_DEG}|/a_0 = "
          f"{float(abs(tail_cheb[-1] / tail_cheb[0])):.2e}", file=sys.stderr)
    tail_w = tail_g_worst(tail_tab, arms)
    for fused in (True, False):
        wg, wt, at = tail_w[fused]
        print(f"    {'fused  ' if fused else 'unfused'}: G worst rel "
              f"2^{float(mp.log(wg, 2)):.2f} at z={at!r}; gt = G/(2z) worst "
              f"2^{float(mp.log(wt, 2)):.2f}; target 2^-63", file=sys.stderr)
        if max(wg, wt) > TAIL_TARGET:
            print("    FAILED: exceeds 2^-63", file=sys.stderr)
            rc = 1
    # Negative control: the same replay must reject a table whose leading
    # coefficient is a plain double (c0's low word dropped).
    bad_lo = [0.0] + tail_lo[1:]
    bad_w = tail_g_worst((tail_c, tail_hi, bad_lo), arms)
    bad = max(max(bad_w[f][0], bad_w[f][1]) for f in (True, False))
    caught = bad > TAIL_TARGET
    print(f"    negative control (c0 low word dropped): worst "
          f"2^{float(mp.log(bad, 2)):.2f} -> "
          f"{'REJECTED' if caught else 'NOT REJECTED'}", file=sys.stderr)
    if not caught:
        print("    FAILED: negative control not rejected", file=sys.stderr)
        rc = 1

    if rc:
        print("One or more self-checks failed -- emitting nothing.", file=sys.stderr)
        return rc

    # --- emit ------------------------------------------------------------------
    recip_n = [dd_split(mp.mpf(1) / n) for n in range(1, R4_NMAX + 1)]
    two_pi = dd_split(2 * mp.pi)
    inv_sqrt_pi = dd_split(1 / mp.sqrt(mp.pi))

    print("// Auto-generated by tools/gen_gamma_data.py. DO NOT EDIT.")
    print("// Fits and their error budgets are documented in the generator;")
    print("// the kernel that consumes them is src/gamma-inl.h.")
    print("#ifndef CORVUS_GAMMA_DATA_H_")
    print("#define CORVUS_GAMMA_DATA_H_")
    print()
    print("namespace corvus::detail {")
    print()
    print("// Temme uniform-asymptotic coefficient table: row k holds the")
    print("// Chebyshev series (NOT monomial -- the kernel Clenshaw-evaluates")
    print("// each row directly) for c_k(eta), eta in [kGammaEtaMid -")
    print("// kGammaEtaHalf, kGammaEtaMid + kGammaEtaHalf] i.e. lambda in")
    print("// [1/2, 2]. S(eta, 1/a) = sum_{k=0}^{10} c_k(eta) / a^k, evaluated")
    print("// by Horner in r=1/a across the 11 rows.")
    print(f"inline constexpr int kGammaTemmeK = {K};")
    print(f"inline constexpr int kGammaTemmeNCoef = {maxdeg + 1};")
    print(f"inline constexpr double kGammaEtaMid = {hexf(eta_mid_d)};")
    print(f"inline constexpr double kGammaEtaHalf = {hexf(eta_half_d)};")
    emit_hex_array_2d("kGammaTemmeCheb", emitted_rows)
    print()
    print(f"// dd pairs of 1/n, n = 1..{R4_NMAX} -- R4's alternating-series weights.")
    emit_hex_array_1d("kGammaRecipNHi", [p[0] for p in recip_n])
    emit_hex_array_1d("kGammaRecipNLo", [p[1] for p in recip_n])
    print()
    print("// dd constants shared by the Temme prefactor and z.lo correction.")
    print(f"inline constexpr double kGammaTwoPiHi = {hexf(two_pi[0])};")
    print(f"inline constexpr double kGammaTwoPiLo = {hexf(two_pi[1])};")
    print(f"inline constexpr double kGammaInvSqrtPiHi = {hexf(inv_sqrt_pi[0])};")
    print(f"inline constexpr double kGammaInvSqrtPiLo = {hexf(inv_sqrt_pi[1])};")
    print()
    print("// Region-map and fixed-length constants (all probe-validated).")
    print(f"inline constexpr double kGammaAT = {hexf(20.0)};")
    print(f"inline constexpr int kGammaSeriesN = {SERIES_NMAX};")
    print(f"inline constexpr int kGammaCfN = {CF_N};")
    print(f"inline constexpr int kGammaR4N = {R4_NMAX};")
    print(f"inline constexpr double kGammaExpFloor = {hexf(EXP_FLOOR)};")
    print()
    print("// Gamma-private erfc tail factor G(z) = z e^{z^2} erfc(z) (#64), for")
    print(f"// z in [{TAIL_ZLO}, {TAIL_ZHI}]: monomial coefficients (ascending) in")
    print("// s = 1/z - kGammaTailC; the first kGammaTailNDd are dd")
    print("// (kGammaTailCoef holds the high words, kGammaTailCoefLo the low).")
    print(f"inline constexpr double kGammaTailZLo = {hexf(TAIL_ZLO)};")
    print(f"inline constexpr double kGammaTailZHi = {hexf(TAIL_ZHI)};")
    print(f"inline constexpr double kGammaTailC = {hexf(tail_c)};")
    print(f"inline constexpr int kGammaTailNCoef = {TAIL_DEG + 1};")
    print(f"inline constexpr int kGammaTailNDd = {TAIL_NDD};")
    emit_hex_array_1d("kGammaTailCoef", tail_hi)
    emit_hex_array_1d("kGammaTailCoefLo", tail_lo)
    print()
    print("}  // namespace corvus::detail")
    print()
    print("#endif  // CORVUS_GAMMA_DATA_H_")
    return 0


if __name__ == "__main__":
    sys.exit(main())
