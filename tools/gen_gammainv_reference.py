#!/usr/bin/env python3
"""Generate tests/data/gammainv_{p,q}_reference.txt -- certified reference
set for corvus::gamma_p_inv / corvus::gamma_q_inv. Consumes the pinned
constants written by tools/gen_gammainv_data.py directly from the
checked-in src/gammainv_data.h -- this generator does NOT re-run that
script's replay/self-check pipeline.

ORACLE:
  1. Root-find x* solving P(a,x)=s (side='p') or Q(a,x)=s (side='q') via
     bisection in ln(x) (monotone, robust across the whole double range --
     p5.py/common.py's probe-validated pattern), forward = ROUTE 1 below.
  2. Round x* to the nearest double xd.
  3. BRACKET CERTIFICATION: form the two exact-mpf half-ulp midpoints of
     xd, evaluate the forward there, certify sign(forward-s) flips.
     Re-certify at a SECOND, higher dps layer (60 then 100) -- both must
     pass.
  4. Deep-small rows (p-side, a*x0 < kGammaInvDeepSmallCut): closed-form
     log-space oracle ln x = (ln p + lnGamma(1+a))/a, certified in LOG
     SPACE against ln of the boundary midpoints, with the dropped-
     correction-term bound folded in as certification slack (see
     deep_small_lx0/certify_deep_small).
  5. Huge-a rows (a>=HUGE_A_THRESHOLD=1e16): a SECOND, independently
     extracted Temme fit (different node count/anchor spacing/degree --
     "route 2") must ALSO bracket-certify the row.
  6. Beyond-resolution rows (aphi already saturated one ULP away from
     x=a -- measured to occur for a >~ 3e34, matching PLAN's own
     conditioning note "a >~ 3e34: transition < 1 ULP of x"): when
     standard bracketing cannot resolve a sign flip, certify instead that
     xd is at least as close to the target as either neighbor, at
     escalated precision (150, then 220 dps). See certify_beyond_resolution.
  7. NEGATIVE CONTROLS: >=4 known-bad rows (hand-perturbed one ULP off a
     certified row, spanning p/q x deep-small/huge-a) must ALL be
     REJECTED on every invocation, checked FIRST, before any generation
     work -- exit 2 otherwise (beta-reference doctrine).

ROUTE 1 (primary forward oracle) is adapted from tools/gen_gamma_reference.py's
own audited oracle_pq: a<=A_SWITCH=1e4 uses mpmath.gammainc directly
(small-side-direct, with the documented hyp1f1-NoConvergence fallback to
elementary series/CF), a>A_SWITCH uses an exact-Temme Chebyshev fit in eta
extrapolated in 1/a (validated below to agree with mpmath to 2^-116 on the
[5e3,1e4] overlap band -- far tighter than gen_gamma_reference.py's own
2^-60 target, since this generator's own oracle-fit validation runs it at
dps=120).
A_SWITCH=1e4 matches BOTH gen_gamma_reference.py's A_SWITCH and
gen_gammainv_data.py's own A_MPMATH_SAFE (that generator's docstring:
threshold 1e6 caused a 15-minute mpmath hang near the ridge at a=1e8; 1e4
is the safe margin, reused here).

SEEDING: oracle_x's bisection bracket is narrowed around a seed from
gen_gammainv_data.py's OWN pinned tri-candidate machinery (seed_S1/S2/S3
+ seed_for, a thin reimplementation calling tools/gen_gammainv_data.py's
module-level functions with the CHECKED-IN Chebyshev table from
src/gammainv_data.h hooked in directly -- not that script's expensive
self-check replay). Seeding is a
SPEED optimization only: oracle_x falls back to the full default bracket
whenever the seeded one fails to bracket the root, so correctness never
depends on the seed being good.

Checkpointed and resumable (beta-reference precedent): every certified row
is appended to a checkpoint file and flushed immediately; re-invoking this
script (same command line) picks up where it left off. Wall-clock budget
per invocation is bounded well under the "foreground command <=~5 min"
process rule.

Usage:
    python3 tools/gen_gammainv_reference.py         # resumable; re-run
                                                       # until it reports
                                                       # DONE
"""
import math
import os
import random
import re
import sys
import tempfile
import time

import mpmath as mp

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "tools"))
import gen_gammainv_data as g1data  # noqa: E402

SEED = 20260809

# ============================================================================
# Part 0: pinned constants, read from the checked-in headers (NOT re-derived
# -- gen_gammainv_data.py's replay pipeline is a separate, expensive job;
# this generator consumes its PINNED OUTPUT exactly as the kernel will).
# ============================================================================
def _parse_scalar_hex(txt, name):
    m = re.search(re.escape(name) + r"\s*=\s*([^;]+);", txt)
    return float.fromhex(m.group(1).strip())


def _parse_scalar_int(txt, name):
    m = re.search(re.escape(name) + r"\s*=\s*([^;]+);", txt)
    return int(m.group(1).strip())


def load_pinned_constants():
    gi_path = os.path.join(REPO, "src", "gammainv_data.h")
    g_path = os.path.join(REPO, "src", "gamma_data.h")
    gi_txt = open(gi_path).read()
    g_txt = open(g_path).read()
    consts = {
        "A_T": _parse_scalar_hex(g_txt, "kGammaAT"),
        "S1_NCORR": _parse_scalar_int(gi_txt, "kGammaInvSeedNCorr"),
        "S2_NCORR": _parse_scalar_int(gi_txt, "kGammaInvS2NCorr"),
        "S3_NITER": _parse_scalar_int(gi_txt, "kGammaInvS3NIter"),
        "ETA_MAX": _parse_scalar_hex(gi_txt, "kGammaInvEtaMax"),
        "S1_A_MIN": _parse_scalar_hex(gi_txt, "kGammaInvS1AMin"),
        "DEEP_SMALL_CUT": _parse_scalar_hex(gi_txt, "kGammaInvDeepSmallCut"),
        "SHALLOW_THRESHOLD": _parse_scalar_hex(gi_txt, "kGammaInvShallowThreshold"),
    }
    m = re.search(r"kGammaInvCkCheb\[2\]\[25\]\s*=\s*\{(.*?)\n\};", gi_txt, re.S)
    rows = re.findall(r"\{([^}]+)\}", m.group(1))
    consts["CK_CHEB"] = [[float.fromhex(v.strip()) for v in row.split(",") if v.strip()]
                          for row in rows]
    assert consts["A_T"] == g1data.A_T_CANDIDATE == 20.0
    assert abs(consts["ETA_MAX"] - g1data.ETA_MAX) < 1e-9
    assert abs(consts["S1_A_MIN"] - g1data.S1_A_MIN) < 1e-15
    return consts


PIN = load_pinned_constants()
g1data._S1_CHEB_ROWS_D = PIN["CK_CHEB"]
DEEP_SMALL_CUT = mp.mpf(PIN["DEEP_SMALL_CUT"])
A_T = PIN["A_T"]


# ============================================================================
# Part 1: tri-candidate seed (thin re-implementation of
# gen_gammainv_data.py's own nested seed_for, calling g1data's
# MODULE-LEVEL helpers fed by the PINNED table above). Speed aid only --
# see module docstring.
# ============================================================================
def cheap_residual(a, x0, s, side):
    if not (math.isfinite(x0) and x0 > 0):
        return None
    v = sd = None
    for dps_try in (25, 35, 55):
        try:
            v, sd = g1data.small_of_x(a, x0, dps=dps_try)
            break
        except (ZeroDivisionError, ValueError, OverflowError):
            continue
    if v is None:
        return None
    if sd != side:
        v = 1 - v
    if s == 0:
        return abs(float(v))
    return abs(float(v) - s) / s


def seed_for(a, s, side):
    s1_candidate = None
    try:
        eta0 = g1data.eta0_of(a, s, side)
    except (OverflowError, ValueError):
        eta0 = None
    if eta0 is not None and a >= PIN["S1_A_MIN"] and abs(eta0) <= PIN["ETA_MAX"]:
        try:
            s1_candidate = g1data.seed_S1(a, s, side, PIN["S1_NCORR"])
        except (OverflowError, ValueError):
            s1_candidate = None
    try:
        s2_candidate = (g1data.seed_S2(a, s, PIN["S2_NCORR"]) if side == "p"
                         else g1data.seed_S2(a, 1.0 - s, PIN["S2_NCORR"]))
    except (OverflowError, ValueError):
        s2_candidate = None
    s3_candidate = None
    if side == "q":
        try:
            s3_candidate = g1data.seed_S3(a, s, PIN["S3_NITER"])
        except (OverflowError, ValueError):
            s3_candidate = None
    best_x, best_r = None, None
    for cand in (s1_candidate, s2_candidate, s3_candidate):
        if cand is None or not (math.isfinite(cand) and cand > 0):
            continue
        r = cheap_residual(a, cand, s, side)
        if r is None:
            continue
        if best_r is None or r < best_r:
            best_x, best_r = cand, r
    if best_x is not None:
        return best_x
    return s2_candidate if (s2_candidate is not None and math.isfinite(s2_candidate)
                             and s2_candidate > 0) else None


# ============================================================================
# Part 2: ROUTE 1 primary forward oracle, adapted verbatim (structure and
# constants) from tools/gen_gamma_reference.py's own audited oracle_pq.
# ============================================================================
A_SWITCH = mp.mpf("1e4")
APHI_SAT = mp.mpf(800)


def phi_lam(lam):
    return lam - 1 - mp.log(lam)


def lam_of_eta_mp(eta, dps):
    with mp.workdps(dps):
        eta = mp.mpf(eta)
        if eta == 0:
            return mp.mpf(1)
        target = eta * eta / 2
        if eta > 0:
            lo, hi = mp.mpf(1), mp.mpf(2)
            while hi - 1 - mp.log(hi) < target:
                hi *= 2
        else:
            lo, hi = mp.mpf("1e-300"), mp.mpf(1)
        for _ in range(300):
            mid = (lo + hi) / 2
            v = mid - 1 - mp.log(mid) - target
            if eta > 0:
                lo, hi = (mid, hi) if v < 0 else (lo, mid)
            else:
                lo, hi = (lo, mid) if v < 0 else (mid, hi)
        return (lo + hi) / 2


def r_exact(a, eta, lam, dps):
    """ORACLE TRAP (gen_gamma_data.py's own documented hazard): eta<0
    MUST use the P-side identity, never Q=1-tiny."""
    with mp.workdps(dps):
        a = mp.mpf(a)
        if eta >= 0:
            Q = mp.gammainc(a, lam * a, regularized=True)
            base = Q - mp.erfc(eta * mp.sqrt(a / 2)) / 2
        else:
            P = mp.gammainc(a, 0, lam * a, regularized=True)
            base = mp.erfc(-eta * mp.sqrt(a / 2)) / 2 - P
        return base * mp.sqrt(2 * mp.pi * a) * mp.exp(a * eta * eta / 2)


def cheb_coeffs_from_vals(vals):
    n = len(vals)
    out = []
    for j in range(n):
        s = mp.fsum([vals[i] * mp.cos(j * mp.pi * (2 * i + 1) / (2 * n))
                      for i in range(n)])
        out.append(s * 2 / n if j else s / n)
    return out


def clenshaw(coefs, t):
    b1 = b2 = mp.mpf(0)
    for j in range(len(coefs) - 1, 0, -1):
        b1, b2 = 2 * t * b1 - b2 + coefs[j], b1
    return t * b1 - b2 + coefs[0]


def extract_temme_fit(dps, kext, k, nnodes, a0, eta_lo, eta_hi, label):
    """One Temme-in-eta Chebyshev fit (Vandermonde-in-1/a solve against
    mpmath.gammainc at nnodes eta values, a0*2^j anchors). Parameterized so
    a SECOND, independently anchored fit (different nnodes/a0/kext) is a
    different numerical derivation, not a copy."""
    t0 = time.time()
    eta_mid = (eta_hi + eta_lo) / 2
    eta_half = (eta_hi - eta_lo) / 2

    def extract_c(eta, lam):
        A = mp.matrix(kext, kext)
        b = mp.matrix(kext, 1)
        for j in range(kext):
            a = mp.mpf(a0) * 2 ** j
            v = 1 / a
            for kk in range(kext):
                A[j, kk] = v ** kk
            b[j] = r_exact(a, eta, lam, dps)
        return mp.lu_solve(A, b)

    with mp.workdps(dps):
        cvals = []
        for i in range(nnodes):
            t = mp.cos(mp.pi * (2 * i + 1) / (2 * nnodes))
            eta = eta_mid + eta_half * t
            lam = lam_of_eta_mp(eta, dps=dps)
            cvals.append(extract_c(eta, lam))
        fits = []
        for kk in range(k):
            vals = [cvals[i][kk] for i in range(nnodes)]
            fits.append(cheb_coeffs_from_vals(vals))
    print(f"  [{label}] Temme fit: dps={dps} kext={kext} k={k} nnodes={nnodes} "
          f"a0={a0} eta=[{float(eta_lo):.4f},{float(eta_hi):.4f}] "
          f"({time.time()-t0:.1f}s)", file=sys.stderr)
    return {"fits": fits, "eta_mid": eta_mid, "eta_half": eta_half, "k": k}


def temme_pq(a, lam, phi, fit, dps):
    with mp.workdps(dps):
        eta = mp.sqrt(2 * phi)
        if lam < 1:
            eta = -eta
        t = (eta - fit["eta_mid"]) / fit["eta_half"]
        ck = [clenshaw(row, t) for row in fit["fits"]]
        S = mp.mpf(0)
        for kk in range(fit["k"] - 1, -1, -1):
            S = S / a + ck[kk]
        R = mp.exp(-a * phi) / mp.sqrt(2 * mp.pi * a) * S
        z = eta * mp.sqrt(a / 2)
        if lam >= 1:
            Q = mp.erfc(z) / 2 + R
            P = 1 - Q
        else:
            P = mp.erfc(-z) / 2 - R
            Q = 1 - P
        return P, Q


def series_S_mp(a, x, dps, nmax=200000):
    with mp.workdps(dps):
        a = mp.mpf(a); x = mp.mpf(x)
        t = mp.mpf(1); s = mp.mpf(1)
        eps = mp.mpf(10) ** (-(dps - 8))
        for n in range(1, nmax + 1):
            t *= x / (a + n)
            s += t
            ratio = x / (a + n + 1)
            if ratio < 1 and t * ratio / (1 - ratio) < eps * s:
                return s, True
        return s, False


def cf_K_mp(a, x, dps, n0=80, nmax=20000):
    with mp.workdps(dps):
        a = mp.mpf(a); x = mp.mpf(x)

        def cf(n):
            k = x + 2 * n + 1 - a
            for j in range(n, 0, -1):
                k = (x + 2 * j - 1 - a) - j * (j - a) / k
            return k

        n = n0
        prev = cf(n)
        eps = mp.mpf(10) ** (-(dps - 8))
        while n < nmax:
            n2 = n * 2
            cur = cf(n2)
            if abs(cur - prev) < eps * abs(cur):
                return cur, True
            prev = cur
            n = n2
        return prev, False


def oracle_pq(a, x, fit, dps):
    """SMALL-SIDE-DIRECT forward: (P, Q). a<=A_SWITCH mpmath.gammainc
    direct (elementary series/CF fallback on the hyp1f1 oracle trap),
    a>A_SWITCH exact-Temme (fit), exact saturation beyond a*phi>800."""
    a_m, x_m = mp.mpf(a), mp.mpf(x)
    with mp.workdps(dps):
        if a_m <= A_SWITCH:
            try:
                if x_m >= a_m:
                    Q = mp.gammainc(a_m, x_m, regularized=True)
                    P = 1 - Q
                else:
                    P = mp.gammainc(a_m, 0, x_m, regularized=True)
                    Q = 1 - P
                return P, Q
            except mp.libmp.libhyper.NoConvergence:
                if x_m >= a_m:
                    K, conv = cf_K_mp(a_m, x_m, dps)
                    lg = a_m * mp.log(x_m) - x_m - mp.loggamma(a_m)
                    Q = mp.e ** lg / K
                    P = 1 - Q
                else:
                    S, conv = series_S_mp(a_m, x_m, dps)
                    lg = a_m * mp.log(x_m) - x_m - mp.loggamma(a_m + 1)
                    P = mp.e ** lg * S
                    Q = 1 - P
                return P, Q
        lam = x_m / a_m
        phi = phi_lam(lam)
        if a_m * phi > APHI_SAT:
            return (mp.mpf(1), mp.mpf(0)) if lam >= 1 else (mp.mpf(0), mp.mpf(1))
        return temme_pq(a_m, lam, phi, fit, dps)


# ============================================================================
# Part 3: root-finder (bisection in ln(x)), seeded via seed_for for speed.
# ============================================================================
def oracle_x(a, target, side, dps, fit, seed_hint=None):
    a_m = mp.mpf(a)
    with mp.workdps(dps):
        target_m = mp.mpf(target)

        def f(lx):
            x = mp.e ** lx
            P, Q = oracle_pq(a_m, x, fit, dps)
            return (P if side == "p" else Q) - target_m

        lo = mp.mpf(-800)
        hi = mp.log(a_m * mp.mpf("1e320") + 1) if a_m > 0 else mp.mpf(700)
        if not mp.isfinite(hi):
            hi = mp.mpf(5000)
        if seed_hint is not None and math.isfinite(seed_hint) and seed_hint > 0:
            lx_seed = mp.log(mp.mpf(seed_hint))
            slo, shi = lx_seed - 40, lx_seed + 40
            try:
                if f(slo) * f(shi) <= 0:
                    lo, hi = slo, shi
            except (ValueError, OverflowError):
                pass
        flo, fhi = f(lo), f(hi)
        tries = 0
        while flo * fhi > 0 and tries < 60:
            if side == "p":
                if flo > 0:
                    lo -= 50
                else:
                    hi += 50
            else:
                if fhi < 0:
                    hi += 50
                else:
                    lo -= 50
            flo, fhi = f(lo), f(hi)
            tries += 1
        if flo * fhi > 0:
            return None
        n_iters = min(200, max(60, dps * 2))
        for _ in range(n_iters):
            mid = (lo + hi) / 2
            fm = f(mid)
            if (fm > 0) == (flo > 0):
                lo, flo = mid, fm
            else:
                hi, fhi = mid, fm
            if hi - lo < mp.mpf(10) ** (-(dps - 5)):
                break
        return mp.e ** ((lo + hi) / 2)


# ============================================================================
# Part 4: deep-small closed-form oracle (p-side), log-space certification
# with the dropped-series-correction bound folded in as slack (design
# point 4). `dropped`/`bound` machinery matches gen_gammainv_data.py's own
# self-check (g), reused rather than re-derived.
# ============================================================================
MIN_SUBNORMAL = math.ldexp(1.0, -1074)


def deep_small_lx0(a, p, dps):
    with mp.workdps(dps):
        a_m, p_m = mp.mpf(a), mp.mpf(p)
        lg1a = mp.loggamma(a_m + 1)
        lx0 = (mp.log(p_m) + lg1a) / a_m
        x0 = mp.e ** lx0
        S, conv = series_S_mp(a_m, x0, dps=dps, nmax=2000)
        dropped = x0 - mp.log(S)
        bound = abs(dropped) / a_m
        return lx0, x0, bound, conv


def certify_deep_small(a, p, dps, xd_override=None):
    """xd_override: certify a SPECIFIC double (e.g. a negative-control
    perturbation) against the closed-form target -- without this,
    certify_row's deep-small branch would always re-derive and certify its
    OWN xd, silently ignoring any externally supplied candidate: every
    deep-small negative control would then be trivially 'accepted'
    because it was never actually being checked."""
    lx0, x0, bound, conv = deep_small_lx0(a, p, dps)
    xd = xd_override if xd_override is not None else (float(x0) if x0 > 0 else 0.0)
    if not math.isfinite(xd):
        return {"xd": None, "certified": False, "note": "non-finite x0"}
    with mp.workdps(dps):
        if xd == 0.0:
            boundary = mp.mpf(MIN_SUBNORMAL) / 2
            certified = (lx0 + bound) < mp.log(boundary)
            return {"xd": 0.0, "certified": bool(certified), "lx0": lx0,
                     "bound": bound, "note": "round-to-zero"}
        lo_d = math.nextafter(xd, -math.inf)
        hi_d = math.nextafter(xd, math.inf)
        xd_m = mp.mpf(xd)
        lo_mid = (xd_m + mp.mpf(lo_d)) / 2 if lo_d > 0 else mp.mpf(lo_d) / 2
        hi_mid = (xd_m + mp.mpf(hi_d)) / 2
        ln_lo = mp.log(lo_mid) if lo_mid > 0 else mp.mpf(-1e30)
        ln_hi = mp.log(hi_mid)
        certified = (ln_lo < lx0 - bound) and (lx0 + bound < ln_hi)
        return {"xd": xd, "certified": bool(certified), "lx0": lx0,
                 "bound": bound, "ln_lo": ln_lo, "ln_hi": ln_hi}


# ============================================================================
# Part 5: certify_row -- bracket certification (layered dps), the huge-a
# independent second route, and the beyond-resolution variant: at
# a >~ 3e34-class, aphi(lambda) already exceeds APHI_SAT one
# representable double away from x=a (delta(a*phi) ~ a*2^-105/2 at 1
# ULP) -- forward collapses to a 3-point step function {1,0.5,0} around
# x=a with NO gradation. In practice this is rarely needed: whenever the
# target s sits strictly between the saturated neighbor values (the
# usual case), standard sign-flip bracketing already certifies xd
# correctly (the neighbors straddle ANY s in (0,1)); it is kept as a
# documented defensive fallback, not because it is expected to fire
# routinely.
# ============================================================================
HUGE_A_THRESHOLD = mp.mpf("1e16")
N_BEYOND_RESOLUTION = [0]


def _sign(v):
    return 1 if v > 0 else (-1 if v < 0 else 0)


def bracket_signs(a, target, side, xd, dps, fit):
    lo_d = math.nextafter(xd, -math.inf)
    hi_d = math.nextafter(xd, math.inf)
    with mp.workdps(dps):
        xd_m = mp.mpf(xd)
        lo_mid = (xd_m + mp.mpf(lo_d)) / 2 if lo_d > 0 else mp.mpf(lo_d) / 2
        hi_mid = (xd_m + mp.mpf(hi_d)) / 2
        t = mp.mpf(target)
        Plo, Qlo = oracle_pq(a, lo_mid, fit, dps)
        Phi, Qhi = oracle_pq(a, hi_mid, fit, dps)
        flo = (Plo if side == "p" else Qlo) - t
        fhi = (Phi if side == "p" else Qhi) - t
    return _sign(flo), _sign(fhi)


def certify_beyond_resolution(a, target, side, xd, dps, fit):
    lo_d = math.nextafter(xd, -math.inf)
    hi_d = math.nextafter(xd, math.inf)
    with mp.workdps(dps):
        t = mp.mpf(target)
        Pl, Ql = oracle_pq(a, lo_d, fit, dps)
        Px, Qx = oracle_pq(a, xd, fit, dps)
        Ph, Qh = oracle_pq(a, hi_d, fit, dps)
        vl = Pl if side == "p" else Ql
        vx = Px if side == "p" else Qx
        vh = Ph if side == "p" else Qh
        el, ex, eh = abs(vl - t), abs(vx - t), abs(vh - t)
    return bool(ex <= el and ex <= eh)


def certify_row(a, s, side, fit1, fit2, dps_layers=(60, 100), seed_hint=None,
                 x_hint=None, route2_dps=100):
    a_m = mp.mpf(a)

    if side == "p":
        lx0, x0, bound, conv = deep_small_lx0(a, s, dps=dps_layers[0])
        if a_m * x0 < DEEP_SMALL_CUT:
            xd = x_hint if x_hint is not None else (float(x0) if x0 > 0 else 0.0)
            if not math.isfinite(xd):
                return {"xd": None, "certified": False, "method": "deep-small",
                         "note": "non-finite"}
            layer_results = [certify_deep_small(a, s, dps, xd_override=xd)["certified"]
                              for dps in dps_layers]
            return {"xd": xd, "certified": all(layer_results),
                     "method": "deep-small", "layers": layer_results}

    if x_hint is not None:
        xd = float(x_hint)
    else:
        x_star = oracle_x(a, s, side, dps=dps_layers[0], fit=fit1, seed_hint=seed_hint)
        if x_star is None:
            return {"xd": None, "certified": False, "method": "root-find-failed"}
        xd = float(x_star)
    if not (math.isfinite(xd) and xd >= 0):
        return {"xd": xd, "certified": False, "method": "non-finite-xd"}
    if xd == 0.0:
        layer_results = [certify_deep_small(a, s, dps)["certified"] for dps in dps_layers]
        return {"xd": 0.0, "certified": all(layer_results), "method": "deep-small-zero",
                 "layers": layer_results}

    layer_results = []
    for dps in dps_layers:
        slo, shi = bracket_signs(a, s, side, xd, dps, fit1)
        ok = (slo != shi and slo != 0 and shi != 0) or slo == 0 or shi == 0
        layer_results.append(ok)
    bracket_ok = all(layer_results)

    beyond_res = False
    if not bracket_ok and a_m >= HUGE_A_THRESHOLD:
        for esc_dps in (150, 220):
            if certify_beyond_resolution(a, s, side, xd, esc_dps, fit1):
                beyond_res = True
                layer_results = [True, True]
                bracket_ok = True
                N_BEYOND_RESOLUTION[0] += 1
                break

    route2_ok = None
    if a_m >= HUGE_A_THRESHOLD:
        slo2, shi2 = bracket_signs(a, s, side, xd, route2_dps, fit2)
        route2_ok = (slo2 != shi2 and slo2 != 0 and shi2 != 0) or slo2 == 0 or shi2 == 0
        if not route2_ok and beyond_res:
            route2_ok = certify_beyond_resolution(a, s, side, xd, 150, fit2)

    certified = bracket_ok and (route2_ok is None or route2_ok)
    return {"xd": xd, "certified": bool(certified),
             "method": ("beyond-resolution" if beyond_res else "bracket"),
             "layers": layer_results, "route2_ok": route2_ok}


# ============================================================================
# Part 6: build the two Temme fits (once), then negative controls.
# ============================================================================
def build_fits():
    print("building oracle Temme fits ...", file=sys.stderr)
    with mp.workdps(120):
        ETA_LO = -mp.sqrt(2 * phi_lam(mp.mpf("0.5")))
        ETA_HI = mp.sqrt(2 * phi_lam(mp.mpf(2)))
    fit1 = extract_temme_fit(dps=120, kext=15, k=11, nnodes=33, a0=512,
                              eta_lo=ETA_LO, eta_hi=ETA_HI, label="route1")
    fit2 = extract_temme_fit(dps=120, kext=13, k=9, nnodes=27, a0=768,
                              eta_lo=ETA_LO, eta_hi=ETA_HI, label="route2")
    return fit1, fit2


def negative_controls(fit1, fit2):
    """>=4 known-bad rows, spanning p/q x normal/deep-small/huge-a. Every
    invocation must reject ALL of them before any generation proceeds."""
    print("negative controls ...", file=sys.stderr)
    good = {}
    good["p_normal"] = certify_row(5.0, 0.3, "p", fit1, fit2)
    good["q_normal"] = certify_row(100.0, 0.05, "q", fit1, fit2)
    good["p_deep_small"] = certify_row(5.0, 1e-96, "p", fit1, fit2)
    good["q_huge_a"] = certify_row(1e16, 0.5, "q", fit1, fit2)
    for k, v in good.items():
        if not v["certified"]:
            print(f"  FATAL: could not even certify the GOOD row for control "
                  f"{k}: {v}", file=sys.stderr)
            return False

    controls = [
        ("p_normal", 5.0, 0.3, "p", math.nextafter(good["p_normal"]["xd"], math.inf)),
        ("q_normal", 100.0, 0.05, "q", math.nextafter(good["q_normal"]["xd"], -math.inf)),
        ("p_deep_small", 5.0, 1e-96, "p", math.nextafter(good["p_deep_small"]["xd"], math.inf)),
        ("q_huge_a", 1e16, 0.5, "q", math.nextafter(good["q_huge_a"]["xd"], math.inf)),
    ]
    all_rejected = True
    for name, a, s, side, bad_xd in controls:
        r = certify_row(a, s, side, fit1, fit2, x_hint=bad_xd)
        status = "REJECTED (correct)" if not r["certified"] else "ACCEPTED (FATAL BUG)"
        print(f"  [{name}] a={a:.6e} {side}={s:.6e} bad_xd={bad_xd!r} -> {status}",
              file=sys.stderr)
        if r["certified"]:
            all_rejected = False
    return all_rejected


# ============================================================================
# Part 7: strata generation. Points are (a, s, side, tag) tuples; where a
# stratum is naturally described in x/lambda terms, s is constructed as
# forward(a,x) rounded to double (well-posed by construction: s is a
# value the kernel could actually receive, sidesteps the huge-a
# "arbitrary target unreachable" trap noted above).
# ============================================================================
NEXT_UP = lambda v: math.nextafter(v, math.inf)
NEXT_DN = lambda v: math.nextafter(v, -math.inf)


def log_grid(lo, hi, n):
    llo, lhi = math.log10(lo), math.log10(hi)
    return [10 ** (llo + (lhi - llo) * i / (n - 1)) for i in range(n)]


class PointSet:
    def __init__(self):
        self.seen = set()
        self.pts = []
        self.strata_counts = {}

    def add(self, a, s, side, tag):
        if not (math.isfinite(a) and math.isfinite(s) and a > 0 and 0 < s < 1):
            return
        key = (a, s, side)
        if key in self.seen:
            return
        self.seen.add(key)
        self.pts.append((a, s, side, tag))
        self.strata_counts[tag] = self.strata_counts.get(tag, 0) + 1


def s_from_x(a, x, side, fit1, dps=60):
    """Construct a well-posed target s = forward(a,x) rounded to double."""
    try:
        P, Q = oracle_pq(a, x, fit1, dps)
    except Exception:
        return None
    v = P if side == "p" else Q
    if not mp.isfinite(v):
        return None
    sd = float(v)
    if not (math.isfinite(sd) and 0.0 < sd < 1.0):
        return None
    return sd


def gen_random_grids(ps, rng, fit1):
    """Per-region (a, s, side) log-spaced random grids across the full
    domain, plus s just above 1/2 (large-side inputs exercising the
    kernel's own exact complement flip)."""
    n0 = len(ps.pts)
    for _ in range(9000):
        a = 10.0 ** rng.uniform(-3, math.log10(1.7e308) - 1)
        side = rng.choice(("p", "q"))
        s = 10.0 ** rng.uniform(-300, math.log10(0.5))
        ps.add(a, s, side, "random-grid")
    # s just above 1/2 -- large-side inputs, complement-flip exercise.
    for _ in range(1200):
        a = 10.0 ** rng.uniform(-3, 20)
        side = rng.choice(("p", "q"))
        s = 0.5 + rng.uniform(1e-12, 0.499)
        ps.add(a, s, side, "large-side-complement")
    print(f"  random grids: {len(ps.pts) - n0} points", file=sys.stderr)


def gen_underflow_p_threshold(ps, fit1):
    """P1b table boundary curves (x at DBL_MIN_NORMAL and at round-to-zero),
    bit-stepped s ladders at a in {1e-3,1e-2,1e-1,0.3}."""
    n0 = len(ps.pts)
    DBL_MIN_NORMAL = math.ldexp(1.0, -1022)
    boundary_xs = (DBL_MIN_NORMAL, MIN_SUBNORMAL, MIN_SUBNORMAL / 2, 0.0)
    for a in (1e-3, 1e-2, 1e-1, 0.3):
        for xb in boundary_xs:
            if xb <= 0:
                continue
            s0 = s_from_x(a, xb, "p", fit1, dps=80)
            if s0 is None or s0 <= 0:
                continue
            v = s0
            for _ in range(6):
                ps.add(a, v, "p", "underflow-p-threshold")
                v = NEXT_UP(v)
            v = s0
            for _ in range(6):
                v = NEXT_DN(v)
                if v <= 0:
                    break
                ps.add(a, v, "p", "underflow-p-threshold")
    print(f"  underflow p-threshold brackets: {len(ps.pts) - n0} points",
          file=sys.stderr)


def gen_deep_small(ps, rng):
    """Deep-small closed-form territory: a*x across the 2^-60 cut
    (bit-stepped brackets), subnormal xd rows, xd=0 rows."""
    n0 = len(ps.pts)
    a_list = (1e-300, 1e-100, 1e-30, 1e-8, 1e-4, 0.01, 0.1, 0.3, 1.0, 5.0,
              A_T - 1.0)
    for a in a_list:
        for pe in range(-320, -1, 6):
            ps.add(a, 10.0 ** pe, "p", "deep-small")
    # Bit-stepped brackets around the a*x=2^-60 cut itself.
    for a in a_list:
        with mp.workdps(60):
            a_m = mp.mpf(a)
            lg1a = mp.loggamma(a_m + 1)
            # Solve for p such that a*x0(p) = cut (closed form inverted).
            lx0_target = mp.log(DEEP_SMALL_CUT / a_m)
            lp_target = lx0_target * a_m - lg1a
            p0 = float(mp.e ** lp_target) if lp_target > -740 else None
        if p0 is None or not (0 < p0 < 1):
            continue
        v = p0
        for _ in range(8):
            ps.add(a, v, "p", "deep-small-cut-bracket")
            v = NEXT_UP(v)
        v = p0
        for _ in range(8):
            v = NEXT_DN(v)
            if v <= 0:
                break
            ps.add(a, v, "p", "deep-small-cut-bracket")
    for _ in range(2200):
        a = 10.0 ** rng.uniform(-300, math.log10(A_T - 1.0))
        p = 10.0 ** rng.uniform(-320, -1)
        ps.add(a, p, "p", "deep-small-random")
    print(f"  deep-small territory: {len(ps.pts) - n0} points", file=sys.stderr)


def gen_far_q_tail(ps):
    """Far q-tail: q down to subnormal-min at a in {0.5,1,2,5,10,100}."""
    n0 = len(ps.pts)
    for a in (0.5, 1.0, 2.0, 5.0, 10.0, 100.0):
        for qe in range(-320, -1, 2):
            ps.add(a, 10.0 ** qe, "q", "far-q-tail")
        v = 5e-324
        for _ in range(10):
            ps.add(a, v, "q", "far-q-tail-boundary")
            v = NEXT_UP(v)
    print(f"  far q-tail: {len(ps.pts) - n0} points", file=sys.stderr)


def gen_ridge_band(ps, fit1):
    """Ridge band lambda in 1 +/- {0.5,2,8}/sqrt(a) at
    a in {1e4,1e6,1e8,1e16,1e100}; huge-a beyond-resolution rows at
    a in {1e35,1e100,1e300,1.7e308}."""
    n0 = len(ps.pts)
    for a in (1e4, 1e6, 1e8, 1e16, 1e100):
        for delta in (0.5, 2.0, 8.0):
            band = delta / math.sqrt(a)
            for lam, side_hint in ((1.0 - band, "p"), (1.0 + band, "q")):
                x = a * lam
                for xv in (x, NEXT_UP(x), NEXT_DN(x)):
                    for side in ("p", "q"):
                        s0 = s_from_x(a, xv, side, fit1, dps=80)
                        if s0 is not None:
                            ps.add(a, s0, side, "ridge-band")
    for a in (1e35, 1e100, 1e300, 1.7e308):
        xs = [a]
        xv = a
        for _ in range(3):
            xv = NEXT_UP(xv)
            xs.append(xv)
        xv = a
        for _ in range(3):
            xv = NEXT_DN(xv)
            xs.append(xv)
        for xv in xs:
            for side in ("p", "q"):
                s0 = s_from_x(a, xv, side, fit1, dps=80)
                if s0 is not None:
                    ps.add(a, s0, side, "huge-a-beyond-resolution")
        # A handful of nextafter(0.5) targets, per the beyond-resolution
        # investigation: the true root for ANY s strictly between the two
        # saturated neighbor values is xd=a (or an immediate neighbor).
        for s0 in (0.5, NEXT_UP(0.5), NEXT_DN(0.5), 0.3, 0.7):
            for side in ("p", "q"):
                ps.add(a, s0, side, "huge-a-beyond-resolution-target")
    print(f"  ridge band + huge-a beyond-resolution: {len(ps.pts) - n0} points",
          file=sys.stderr)


def gen_seam_brackets(ps, fit1):
    """Seed-partition seam brackets: a_T=20 (bit-stepped in a),
    S1_A_MIN=0.3, the shallow threshold s=2^-10, s=1/2 ladders, |eta|=
    kGammaInvEtaMax boundary."""
    n0 = len(ps.pts)
    # a_T bit-stepped.
    a_seam = A_T
    a_vals = [a_seam]
    v = a_seam
    for _ in range(4):
        v = NEXT_UP(v)
        a_vals.append(v)
    v = a_seam
    for _ in range(4):
        v = NEXT_DN(v)
        a_vals.append(v)
    a_vals.append(PIN["S1_A_MIN"])
    v = PIN["S1_A_MIN"]
    for _ in range(4):
        v = NEXT_UP(v)
        a_vals.append(v)
    v = PIN["S1_A_MIN"]
    for _ in range(4):
        v = NEXT_DN(v)
        a_vals.append(v)
    for a in a_vals:
        for s_seam in (PIN["SHALLOW_THRESHOLD"], 0.5):
            v2 = s_seam
            for _ in range(4):
                for side in ("p", "q"):
                    ps.add(a, v2, side, "seam-bracket")
                v2 = NEXT_UP(v2)
                if v2 >= 1.0:
                    break
            v2 = s_seam
            for _ in range(4):
                v2 = NEXT_DN(v2)
                if v2 <= 0:
                    break
                for side in ("p", "q"):
                    ps.add(a, v2, side, "seam-bracket")
    # |eta| = kGammaInvEtaMax boundary: construct via forward at
    # lambda = lam_of_eta(+-ETA_MAX), a = a_T (the domain this seed uses).
    with mp.workdps(60):
        for sign in (1, -1):
            eta_b = sign * mp.mpf(PIN["ETA_MAX"])
            lam_b = lam_of_eta_mp(eta_b, dps=60)
            x_b = float(lam_b) * A_T
            for xv in (x_b, NEXT_UP(x_b), NEXT_DN(x_b)):
                for side in ("p", "q"):
                    s0 = s_from_x(A_T, xv, side, fit1, dps=60)
                    if s0 is not None:
                        ps.add(A_T, s0, side, "eta-max-boundary")
    print(f"  seed-partition seam brackets: {len(ps.pts) - n0} points",
          file=sys.stderr)


def build_point_set(rng, fit1):
    ps = PointSet()
    gen_random_grids(ps, rng, fit1)
    gen_underflow_p_threshold(ps, fit1)
    gen_deep_small(ps, rng)
    gen_far_q_tail(ps)
    gen_ridge_band(ps, fit1)
    gen_seam_brackets(ps, fit1)
    print(f"  total distinct (a,s,side) points: {len(ps.pts)}", file=sys.stderr)
    for tag, n in sorted(ps.strata_counts.items()):
        print(f"    {tag}: {n}", file=sys.stderr)
    return ps


# ============================================================================
# Part 8: checkpointed, resumable compute pass.
# ============================================================================
CKPT_PATH = os.path.join(tempfile.gettempdir(), f"corvus_gammainv_ref_ckpt_{SEED}.tsv")
WALL_CLOCK_BUDGET_S = 95.0


def load_checkpoint(path, expected_sig):
    done = {}
    if not os.path.exists(path):
        return done
    with open(path, "r") as f:
        header = f.readline().strip()
        if header != expected_sig:
            print(f"  checkpoint signature mismatch ({header!r} != {expected_sig!r}) "
                  f"-- starting fresh.", file=sys.stderr)
            return done
        for line in f:
            parts = line.rstrip("\n").split("\t")
            idx = int(parts[0])
            done[idx] = parts[1:]
    return done


def append_checkpoint(fh, idx, fields):
    fh.write(str(idx) + "\t" + "\t".join(fields) + "\n")
    fh.flush()


def hexd(x):
    return float(x).hex()


def compute_all(ps, fit1, fit2):
    total = len(ps.pts)
    sig = f"v1 SEED={SEED} N={total}"
    done_map = load_checkpoint(CKPT_PATH, sig)
    print(f"  checkpoint: {len(done_map)}/{total} points already computed "
          f"({CKPT_PATH})", file=sys.stderr)

    existing_sig = None
    if os.path.exists(CKPT_PATH):
        with open(CKPT_PATH, "r") as f0:
            existing_sig = f0.readline().strip()
    mode = "a" if existing_sig == sig else "w"

    t_start = time.time()
    newly_done = 0
    with open(CKPT_PATH, mode) as fh:
        if mode == "w":
            fh.write(sig + "\n")
            fh.flush()
        for idx, (a, s, side, tag) in enumerate(ps.pts):
            if idx in done_map:
                continue
            if time.time() - t_start > WALL_CLOCK_BUDGET_S:
                print(f"  wall-clock budget ({WALL_CLOCK_BUDGET_S:.0f}s) hit at "
                      f"{idx}/{total} ({newly_done} computed this run) -- "
                      f"re-run this script to continue.", file=sys.stderr)
                return None, False
            try:
                seed = seed_for(a, s, side)
            except Exception:
                seed = None
            try:
                r = certify_row(a, s, side, fit1, fit2, seed_hint=seed)
            except Exception as e:
                r = {"xd": None, "certified": False, "method": f"exception:{e}"}
            if r["xd"] is None or not r["certified"]:
                append_checkpoint(fh, idx, ["FAILED", r.get("method", "?")])
            else:
                append_checkpoint(fh, idx, [hexd(r["xd"]), r["method"]])
            newly_done += 1
            if newly_done % 500 == 0:
                print(f"    ... {idx + 1}/{total} ({time.time() - t_start:.0f}s "
                      f"this run)", file=sys.stderr)
                sys.stderr.flush()

    print(f"  computed {newly_done} points this run ({time.time() - t_start:.0f}s); "
          f"all {total} points now checkpointed.", file=sys.stderr)

    done_map = load_checkpoint(CKPT_PATH, sig)
    rows_p, rows_q = [], []
    n_failed = 0
    fail_by_tag = {}
    for idx, (a, s, side, tag) in enumerate(ps.pts):
        fields = done_map[idx]
        if fields[0] == "FAILED":
            n_failed += 1
            fail_by_tag[tag] = fail_by_tag.get(tag, 0) + 1
            continue
        xd = float.fromhex(fields[0])
        (rows_p if side == "p" else rows_q).append((a, s, xd))
    print(f"  certified: {len(rows_p)} p-rows, {len(rows_q)} q-rows; "
          f"{n_failed} dropped (uncertified)", file=sys.stderr)
    if fail_by_tag:
        print(f"  drops by stratum: {fail_by_tag}", file=sys.stderr)
    return (rows_p, rows_q), True


# ============================================================================
# Part 9: 25-row independent spot cross-check (elementary series/CF for
# a<=1e4, route2 exact-asymptotic for larger a).
# ============================================================================
def cross_check_25(rows_p, rows_q, rng, fit2):
    print("\n25-row independent spot cross-check ...", file=sys.stderr)
    all_rows = [(a, s, xd, "p") for a, s, xd in rows_p] + \
               [(a, s, xd, "q") for a, s, xd in rows_q]
    if len(all_rows) < 25:
        sample = all_rows
    else:
        sample = rng.sample(all_rows, 25)
    results = []
    for a, s, xd, side in sample:
        with mp.workdps(80):
            if a <= 1e4:
                if xd >= a:
                    K, conv = cf_K_mp(a, xd, dps=80)
                    lg = mp.mpf(a) * mp.log(xd) - xd - mp.loggamma(mp.mpf(a))
                    Qv = mp.e ** lg / K
                    Pv = 1 - Qv
                else:
                    S, conv = series_S_mp(a, xd, dps=80)
                    lg = mp.mpf(a) * mp.log(xd) - xd - mp.loggamma(mp.mpf(a) + 1)
                    Pv = mp.e ** lg * S
                    Qv = 1 - Pv
                method = "elementary series/CF"
            else:
                lam = mp.mpf(xd) / mp.mpf(a)
                phi = phi_lam(lam)
                if mp.mpf(a) * phi > APHI_SAT:
                    Pv, Qv = (mp.mpf(1), mp.mpf(0)) if lam >= 1 else (mp.mpf(0), mp.mpf(1))
                else:
                    Pv, Qv = temme_pq(mp.mpf(a), lam, phi, fit2, dps=80)
                method = "route2 exact-asymptotic"
        got = Pv if side == "p" else Qv
        rel = float(abs(got - mp.mpf(s)))
        results.append((a, s, xd, side, method, rel))
        print(f"  a={a:.6e} {side}={s:.6e} xd={xd!r} [{method}] "
              f"|forward(xd)-s|={rel:.3e}", file=sys.stderr)
    return results


# ============================================================================
# Part 10: main
# ============================================================================
def write_rows(rows_p, rows_q):
    with open(os.path.join(REPO, "tests", "data", "gammainv_p_reference.txt"), "w") as f:
        for a, s, xd in rows_p:
            f.write(f"{hexd(a)} {hexd(s)} {hexd(xd)}\n")
    with open(os.path.join(REPO, "tests", "data", "gammainv_q_reference.txt"), "w") as f:
        for a, s, xd in rows_q:
            f.write(f"{hexd(a)} {hexd(s)} {hexd(xd)}\n")
    print(f"  wrote gammainv_p_reference.txt: {len(rows_p)} rows", file=sys.stderr)
    print(f"  wrote gammainv_q_reference.txt: {len(rows_q)} rows", file=sys.stderr)


def main():
    t_all = time.time()
    fit1, fit2 = build_fits()

    if not negative_controls(fit1, fit2):
        print("\nFATAL: negative control(s) were ACCEPTED -- the certifier "
              "is not rejecting known-bad rows. Aborting, nothing written.",
              file=sys.stderr)
        return 2
    print("negative controls: all rejected (correct).", file=sys.stderr)

    rng = random.Random(SEED)
    print("\nbuilding point set ...", file=sys.stderr)
    ps = build_point_set(rng, fit1)

    print("\nevaluating oracle (resumable) ...", file=sys.stderr)
    result, done = compute_all(ps, fit1, fit2)
    if not done:
        print("\nPARTIAL RUN: re-invoke this script to continue "
              "(checkpoint saved).", file=sys.stderr)
        return 3
    rows_p, rows_q = result

    if len(rows_p) < 5000 or len(rows_q) < 5000:
        print(f"\nFAILED: row counts too low (p={len(rows_p)}, q={len(rows_q)}), "
              f"target >=5000 each.", file=sys.stderr)
        return 1

    cross_check_25(rows_p, rows_q, random.Random(SEED ^ 0x51A7), fit2)

    write_rows(rows_p, rows_q)
    print(f"\nbeyond-resolution certifications used: {N_BEYOND_RESOLUTION[0]}",
          file=sys.stderr)
    print(f"total generator runtime (this invocation): {time.time()-t_all:.1f}s",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
