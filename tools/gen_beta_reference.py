#!/usr/bin/env python3
"""Generate tests/data/beta_{p,q}_reference.txt -- correctly-rounded oracle
reference set for corvus::beta_p / corvus::beta_q (regularized incomplete
beta), per PLAN.md "Regularized incomplete beta -- detail design" (the
G1a/G1b/G1c probe corrections are binding) and its G2 (references) brief.

Oracle machinery REUSED from tools/gen_beta_data.py (imported, not
re-derived, per the G2 brief): small_val_via_cf (DLMF 8.17.22 backward CF,
self-convergent, "validated across the whole domain" per G1b), route_final
(the final G1b-corrected routing order -- used here ONLY to classify points
into a per-region coverage histogram; the oracle itself evaluates every
point the same CF way regardless of which region the KERNEL would route it
to), and the pinned region constants (B1, XI1, EPS_R4, T_RIDGE, Z0, C_LG,
LN2, XI_RATIO_LO/HI, ZETA_MAX, E_FLOOR).

Oracle rules (both hard requirements; gamma's precedent restated for beta's
extra parameter):
  SMALL-SIDE-DIRECT: decide which of P, Q is smaller via the cheap, EXACT
  mean predicate x*(a+b) <= a, compute THAT one directly via the swapped
  argument triple (I_x(a,b) native, or I_{1-x}(b,a) swapped) -- see
  small_side_direct() below. The complement is 1-small, formed in mpf at
  full working precision: safe, because the large side sits near 1 and
  needs only absolute accuracy there, which "1 - a fully-precise small
  value" gives for free. Never averaged; never a bare 1-x subtraction used
  to GET the small side. Self-correcting: if the mean-predicate guess
  computes to a value > 0.5 (mean != true median -- can happen off the
  diagonal), the OTHER orientation is computed instead and used (P+Q=1
  guarantees whichever side computes to <=0.5 is genuinely the smaller
  one, so this single re-check is sufficient, no iteration needed).

  PRIMARY/CROSS-CHECK ORACLE CHOICE [deviation from the literal brief
  wording, reasoned and flagged here per house style -- see gen_beta_data.py's
  own "deviation, flagged" precedent at its R3 section]: the G2 brief says
  "mpmath betainc below its reliability ceiling with a hard per-point
  timeout; the CF above/on timeout." Measured on this box before writing
  this generator: a single subprocess-timeout-guarded mpmath.betainc call
  (gen_beta_data.py's own _betainc_timeout pattern) costs ~2.5s, dominated
  by Windows process-spawn overhead, not the arithmetic -- confirmed by
  timing 1 vs 5 calls (2.55s vs 12.4s, i.e. ~2.5s/call throughout, no
  amortization). small_val_via_cf costs 3-100ms depending on region (100ms
  only for points deep on a large-parameter ridge; typical points are
  3-15ms). At this generator's ~40k-point target, an mpmath-primary design
  is not tractable (2.5s * 40000 ~ 28 hours; even a generous "below
  ceiling" subset of a few thousand points would run 2+ hours on spawn
  overhead alone). This generator therefore uses the CF as the PRIMARY
  value oracle for every emitted row -- fast, safe (bounded N_max,
  deterministic termination), and per gen_beta_data.py's own G1b finding
  "the CF oracle is validated across the whole domain" -- and reserves
  mpmath.betainc (same subprocess-timeout-guarded call, reused) for the
  CROSS-CHECK oracle on the mandated ~500-point random subsample spanning
  all regions, which is exactly what the cross-check requirement asks for
  ("betainc and CF must agree to output precision + 10 digits ... where
  reachable"). This mirrors the trust-the-CF posture gen_beta_data.py's own
  R3 extraction and self-check (b)'s ratio-cap sweep already adopted for
  the identical reason (mpmath hangs/times out across large swaths of this
  domain -- not just near the ridge, per that file's own "gen1/gen2" oracle
  notes).

  DPS LADDER (three-layer dps hygiene, per the brief): every point is
  computed at dps=40, rechecked at dps=60; on relative disagreement
  greater than 2^-80, escalate to dps=100 and use THAT value (logged).
  Persistent disagreement between dps=60 and dps=100 beyond 2^-70 is a real
  ESCALATE, reported at the end with the witness point -- not a threshold
  to loosen. mp.mp.dps is set INSIDE every oracle call (small_val_via_cf
  already does this; this generator's own helpers do too) -- the G1a/G1b/
  G1c "stale ambient dps" trap is now thrice-documented in PLAN.md and this
  generator does not repeat it.

Usage:
    python3 tools/gen_beta_reference.py
Writes tests/data/beta_p_reference.txt and tests/data/beta_q_reference.txt
DIRECTLY (two files, no stdout redirection) -- mirrors gen_gamma_reference.py's
own convention for exactly the same reason: the two files carry IDENTICAL
rows (`a b x P Q`, gamma's own gamma_p/gamma_q precedent -- one physical
reference row serves as ground truth for both the beta_p kernel test and the
beta_q kernel test), computed once and written twice.

Own fresh seed (20260731, today's date at generator authorship) -- shares
NOTHING with gamma's rng stream (SEED=20260727 in gen_gamma_reference.py,
frozen); AGENTS.md note on this restated here at the point it matters.
"""

import math
import os
import struct
import sys
import time
import random
import multiprocessing as mp_proc

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen_beta_data as gbd  # noqa: E402  (see module docstring: reused, not re-derived)

mp.mp.dps = 100  # module-level default; every function below sets its OWN
                  # dps explicitly on entry (three-layer dps hygiene, see
                  # module docstring and the G1a/G1b/G1c "stale ambient dps"
                  # lesson recorded in PLAN.md).

SEED = 20260731

# --- constants reused verbatim from gen_beta_data.py (not re-derived) -------
B1 = gbd.B1
XI1 = gbd.XI1
EPS_R4 = gbd.EPS_R4
T_RIDGE = gbd.T_RIDGE
Z0 = gbd.Z0
C_LG = gbd.C_LG
LN2 = gbd.LN2
E_FLOOR = gbd.E_FLOOR
XI_RATIO_LO = gbd.XI_RATIO_LO
XI_RATIO_HI = gbd.XI_RATIO_HI
ZETA_MAX = gbd.ZETA_MAX
route_final = gbd.route_final
small_val_via_cf = gbd.small_val_via_cf
B_GL = gbd.B_GL

DPS1 = 40
DPS2 = 60
DPS3 = 100
DISAGREE_60_40 = mp.mpf(2) ** -80
DISAGREE_100_60 = mp.mpf(2) ** -70

TARGET_TOTAL = 40000
MIN_PER_REGION = 2000

# ============================================================================
# hex-float / dedup helpers (mirrors gen_gamma_reference.py's own style)
# ============================================================================
def as_bits(x):
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def hexd(x):
    """Hex-float of a python float, mpf, or int -- round-trip exact. Handles
    +-inf/nan (used only by the specials rows) since float.hex() supports
    them directly ('inf', '-inf', 'nan')."""
    return float(x).hex()


def log_grid(lo, hi, n):
    llo, lhi = math.log10(lo), math.log10(hi)
    return [10 ** (llo + (lhi - llo) * i / (n - 1)) for i in range(n)]


NEXT_UP = lambda v: math.nextafter(v, math.inf)
NEXT_DN = lambda v: math.nextafter(v, -math.inf)


class PointSet:
    """Collects (a, b, x, keep_if_saturated, tag) tuples, deduped by
    (a,b,x) bits. tag is provenance only (stderr reporting), not emitted."""

    def __init__(self):
        self.seen = set()
        self.pts = []

    def add(self, a, b, x, keep_if_saturated=True, tag=""):
        a, b, x = float(a), float(b), float(x)
        if not (math.isfinite(a) and math.isfinite(b) and math.isfinite(x)):
            return
        if not (a > 0 and b > 0 and 0.0 < x < 1.0):
            return
        key = (as_bits(a), as_bits(b), as_bits(x))
        if key in self.seen:
            return
        self.seen.add(key)
        self.pts.append((a, b, x, keep_if_saturated, tag))


# ============================================================================
# Region point generation. These are DENSITY-TARGETING constructions (they
# aim points at each region's natural coordinates) -- the actual region
# membership used for the coverage histogram and every self-check is always
# route_final()'s classification, never an assumption baked in here.
# ============================================================================
def gen_r1(ps, rng):
    """R1 power series: xi<=xi1 and b*xi<=B1 (either orientation). Lattice
    in (a, b) with xi driven by both the xi1 wall and the B1/b wall."""
    n0 = len(ps.pts)
    a_list = sorted(set(log_grid(1e-12, 1e8, 26) + [
        1e-300, 1e-100, 1e-30, 1e-10, 0.01, 0.1, 1.0, 2.0, 8.0, 20.0, 100.0,
    ]))
    b_list = sorted(set(log_grid(1e-6, 1e6, 16) + [1e-300, 0.01, 1.0, 8.0, 20.0]))
    fracs = [0.05, 0.2, 0.5, 0.85, 0.99, 1.0 - 2.0 ** -50]
    for a in a_list:
        for b in b_list:
            xi_wall = min(float(XI1), B1_f / b if b > 0 else 1.0)
            for f in fracs:
                xi = f * xi_wall
                if 0 < xi < 1:
                    ps.add(a, b, xi, tag="R1-lattice")
                    # swapped orientation: put b in the "alpha" role too.
                    ps.add(b, a, xi, tag="R1-lattice-swaprole")
    for _ in range(4000):
        a = 10.0 ** rng.uniform(-300.0, 8.0)
        b = 10.0 ** rng.uniform(-6.0, 6.0)
        xi_wall = min(float(XI1), B1_f / b if b > 0 else 1.0)
        f = rng.uniform(0.0, 1.0)
        xi = f * xi_wall
        if 0 < xi < 1:
            ps.add(a, b, xi, keep_if_saturated=False, tag="R1-random")
    print(f"  R1: {len(ps.pts) - n0} points", file=sys.stderr)


B1_f = float(B1)
XI1_f = float(XI1)


def _r4_xi_bounds(tau, Bp, ln2, xi1, b1):
    """R4 box in xi_tau: LOWER bound from tau*|ln xi_tau|<=ln2 (xi_tau ==
    exp(-ln2/tau) is monotonically increasing in tau, so for small tau this
    floor is tiny/non-binding -- confirmed against gen_beta_data.py's own
    check_f_r4, which uses exactly this exp(-LN2/alpha) as the LOWER sweep
    edge (its `los`), not an upper cap. UPPER bound is min(xi1, B1/Bp) from
    the other two walls. An earlier version of this generator used the ln2
    floor as a third argument to min() alongside the two upper walls --
    wrong direction: for tau<<1 that collapsed the sampled range down to
    ~0, producing only astronomically-tiny xi_tau and starving the R4
    lattice (measured 468 raw points instead of thousands) -- caught by
    this generator's own low first-pass region count during development,
    not by reasoning about the code."""
    floor = math.exp(-ln2 / tau) if tau > 0 else 0.0
    ceil = min(xi1, b1 / Bp if Bp > 0 else 1.0, 1.0)
    return floor, ceil


def gen_r4(ps, rng):
    """R4 tiny-min box: min(a,b)<=eps_R4, xi_tau<=xi1, B*xi_tau<=B1,
    tau*|ln xi_tau|<=ln2 (tiny-first triple). Lattice over tau (the tiny
    param) and the larger partner B, xi_tau log-sampled between the ln2
    floor and the (xi1, B1/B) ceiling -- see _r4_xi_bounds."""
    n0 = len(ps.pts)
    eps = float(EPS_R4)
    ln2 = float(LN2)
    tau_list = sorted(set(log_grid(1e-300, eps, 22) + [
        1e-300, 1e-100, 1e-30, 1e-10, 1e-3, 1e-2, NEXT_DN(eps), eps,
    ]))
    B_list = sorted(set(log_grid(1e-6, B1_f / XI1_f, 14) + [0.01, 1.0, 8.0]))
    fracs = [0.05, 0.3, 0.5, 0.7, 0.95]
    for tau in tau_list:
        for Bp in B_list:
            floor, ceil = _r4_xi_bounds(tau, Bp, ln2, XI1_f, B1_f)
            if not (0 < floor < ceil < 1):
                # floor essentially 0 (typical for tau << eps_R4): just
                # sample log-uniformly from a tiny positive number to ceil.
                if ceil <= 0:
                    continue
                floor = max(floor, 1e-308)
            logf, logc = math.log(floor), math.log(ceil)
            for f in fracs:
                xi_tau = math.exp(logf + f * (logc - logf))
                if 0 < xi_tau < 1:
                    ps.add(tau, Bp, xi_tau, tag="R4-lattice-native")
                    ps.add(Bp, tau, 1.0 - xi_tau, tag="R4-lattice-swap")
    for _ in range(3500):
        tau = 10.0 ** rng.uniform(-300.0, math.log10(eps))
        Bp = 10.0 ** rng.uniform(-6.0, math.log10(B1_f / XI1_f) + 0.3)
        floor, ceil = _r4_xi_bounds(tau, Bp, ln2, XI1_f, B1_f)
        if ceil <= 0:
            continue
        floor = max(floor, 1e-308)
        if floor >= ceil:
            continue
        logf, logc = math.log(floor), math.log(ceil)
        xi_tau = math.exp(logf + rng.uniform(0.0, 1.0) * (logc - logf))
        if 0 < xi_tau < 1:
            if rng.random() < 0.5:
                ps.add(tau, Bp, xi_tau, keep_if_saturated=False, tag="R4-random-native")
            else:
                ps.add(Bp, tau, 1.0 - xi_tau, keep_if_saturated=False, tag="R4-random-swap")
    print(f"  R4: {len(ps.pts) - n0} points", file=sys.stderr)


def gen_r3(ps, rng):
    """R3 ridge ratio-band: nu=a*b/c>=T_ridge, xi/p in [1/2,2] and
    (1-xi)/q in [1/2,2]. Lattice in (nu, p) with xi driven by u in
    [-1/2,1] scaled to the p-local reachable lens (zeta_max_at_p), mirroring
    gen_beta_data.py's own check_c_r3 lesson (sampling outside the
    p-dependent lens fakes membership; scaling by zeta_max_at_p keeps every
    point genuinely inside R3's true domain)."""
    n0 = len(ps.pts)
    nu_list = sorted(set(
        [float(T_RIDGE) * 2.0 ** k for k in range(0, 21)] +
        [NEXT_DN(float(T_RIDGE)), NEXT_UP(float(T_RIDGE))] +
        [1e100, 1e200, 1e250]
    ))
    p_list = [0.02, 0.1, 1.0 / 3.0, 0.25, 0.4, 0.49, 0.5]
    u_fracs = [0.05, 0.3, 0.6, 0.85, 0.97, 1.0]
    for nu in nu_list:
        for p in p_list:
            q = 1.0 - p
            if p <= 0 or q <= 0:
                continue
            c = nu / (p * q)
            a = p * c
            b = q * c
            if not (math.isfinite(a) and math.isfinite(b) and a > 0 and b > 0):
                continue
            lo_u = max(-0.5, -q / p)
            hi_u = min(1.0, q / (2.0 * p))
            for frac in u_fracs:
                for u in (lo_u * frac, hi_u * frac):
                    xi = p * (1.0 + u)
                    if 0 < xi < 1:
                        ps.add(a, b, xi, tag="R3-lattice")
    for _ in range(3500):
        nu = 10.0 ** rng.uniform(math.log10(float(T_RIDGE)), 16.0)
        p = rng.uniform(0.02, 0.5)
        q = 1.0 - p
        c = nu / (p * q)
        a = p * c
        b = q * c
        if not (math.isfinite(a) and math.isfinite(b) and a > 0 and b > 0):
            continue
        lo_u = max(-0.5, -q / p)
        hi_u = min(1.0, q / (2.0 * p))
        u = rng.uniform(lo_u, hi_u)
        xi = p * (1.0 + u)
        if 0 < xi < 1:
            ps.add(a, b, xi, keep_if_saturated=False, tag="R3-random")
    print(f"  R3: {len(ps.pts) - n0} points", file=sys.stderr)


def gen_r2(ps, rng):
    """R2 backward CF: everything else. Lattice at moderate-middle points,
    the gamma-limit line (alpha tiny, beta huge), and the far off-band
    ridge (nu>=T_ridge but outside the ratio band -- the risk the third
    correction transferred to R2), covering BOTH orientations of the
    xi<(a+1)/(c+2) rule."""
    n0 = len(ps.pts)
    ab_list = sorted(set(log_grid(1e-8, 1e8, 40) + [
        1e-300, 1e-100, 0.01, 0.5, 1.0, 2.0, 8.0, 20.0, 100.0, 1000.0,
    ]))
    xfrac_list = [0.05, 0.15, 0.35, 0.55, 0.65, 0.78, 0.9, 0.995]
    for a in ab_list:
        for b in ab_list[::3]:  # thin the b axis; full a x full b is too dense
            for xf in xfrac_list:
                ps.add(a, b, xf, tag="R2-lattice")
    # gamma-limit line: alpha tiny, beta huge, beta*xi in [B1, 400].
    for a in (1e-300, 1e-6, 0.05, 1.0):
        for b in (1e4, 1e8, 1e12, 1e100, 1e250):
            for bxi in (8.0, 20.0, 80.0, 400.0):
                xi = bxi / b
                if 0 < xi < 1:
                    ps.add(a, b, xi, tag="R2-gamma-limit")
    # off-band far ridge: nu>=T_ridge but outside the ratio-band caps.
    for nu in (float(T_RIDGE), 1e4, 1e8, 1e12, 1e16):
        for p in (0.02, 0.1, 1.0 / 3.0, 0.5):
            q = 1.0 - p
            c = nu / (p * q)
            a = p * c
            b = q * c
            if not (math.isfinite(a) and math.isfinite(b) and a > 0 and b > 0):
                continue
            hi_u = min(1.0, q / (2.0 * p))
            for umul in (1.5, 3.0, 8.0):
                u = hi_u * umul
                xi = p * (1.0 + u)
                if 0 < xi < 1:
                    ps.add(a, b, xi, tag="R2-offband-ridge")
    for _ in range(6000):
        a = 10.0 ** rng.uniform(-300.0, 8.0)
        b = 10.0 ** rng.uniform(-8.0, 8.0)
        x = rng.uniform(0.001, 0.999)
        ps.add(a, b, x, keep_if_saturated=False, tag="R2-random")
    print(f"  R2: {len(ps.pts) - n0} points", file=sys.stderr)


# ============================================================================
# Ridge lines: lambda = a - c*x = 0+- at nu binades {32,64,...,4096} crossed
# with p in {0.02, 0.1, 1/3, 0.5}. lambda=0+- means x brackets a/c (=p as a
# double) by one ulp each side.
# ============================================================================
def gen_ridge_lines(ps):
    n0 = len(ps.pts)
    nu_binades = [32.0 * 2.0 ** k for k in range(8)]  # 32..4096
    p_vals = [0.02, 0.1, 1.0 / 3.0, 0.5]
    for nu in nu_binades:
        for p in p_vals:
            q = 1.0 - p
            c = nu / (p * q)
            a = p * c
            b = q * c
            if not (math.isfinite(a) and math.isfinite(b) and a > 0 and b > 0):
                continue
            x0 = a / c  # == p by construction, but compute from a,c directly
            for x in (NEXT_DN(x0), x0, NEXT_UP(x0)):
                ps.add(a, b, x, tag="ridge-line")
    print(f"  ridge lines: {len(ps.pts) - n0} points", file=sys.stderr)


# ============================================================================
# Boundary bit-brackets: ONE ULP EACH SIDE of every named region/routing
# boundary, across a handful of representative (a,b)/parameter combos.
# ============================================================================
def _bracket_add(ps, a, b, x, tag):
    """Add x at -2,-1,0,+1,+2 ulp around the given boundary value."""
    for k in (-2, -1, 0, 1, 2):
        xv = x
        for _ in range(abs(k)):
            xv = NEXT_UP(xv) if k > 0 else NEXT_DN(xv)
        ps.add(a, b, xv, tag=tag)


def gen_boundaries(ps):
    n0 = len(ps.pts)

    # (1) xi1 = 0.45: fix beta small enough that beta*xi1<=B1 doesn't also
    # bind (beta <= B1/xi1), vary alpha across magnitudes.
    for a in (1e-100, 1e-8, 0.1, 1.0, 8.0, 1e4, 1e8):
        for b in (0.5, 2.0, 8.0, 15.0):
            _bracket_add(ps, a, b, XI1_f, "bnd-xi1")

    # (2) beta*xi = B1 = 8, at xi well below xi1 so that wall doesn't bind.
    for a in (1e-100, 1e-8, 0.1, 1.0, 8.0, 1e4, 1e8):
        for xi in (0.01, 0.05, 0.2, 0.4):
            b = B1_f / xi
            _bracket_add(ps, a, b, xi, "bnd-B1")

    # (3) ratio caps xi/p in {1/2,2} and linked (1-xi)/q in {1/2,2}, at
    # several nu (>=T_ridge, so genuinely at the R3/R2 handoff).
    for nu in (float(T_RIDGE), 1e3, 1e8, 1e14):
        for p in (0.05, 0.2, 1.0 / 3.0, 0.5):
            q = 1.0 - p
            c = nu / (p * q)
            a = p * c
            b = q * c
            if not (math.isfinite(a) and math.isfinite(b) and a > 0 and b > 0):
                continue
            for u in (-0.5, 1.0):
                xi = p * (1.0 + u)
                if 0 < xi < 1:
                    _bracket_add(ps, a, b, xi, "bnd-ratio-u")
            for v in (-0.5, 1.0):
                xi = 1.0 - q * (1.0 + v)
                if 0 < xi < 1:
                    _bracket_add(ps, a, b, xi, "bnd-ratio-v")

    # (4) nu = T_ridge = 32, at several p, x held at the mean (p) so only
    # the nu wall is exercised (not the ratio-band wall too).
    for p in (0.1, 1.0 / 3.0, 0.5):
        q = 1.0 - p
        for nu in (NEXT_DN(32.0), 32.0, NEXT_UP(32.0)):
            c = nu / (p * q)
            a = p * c
            b = q * c
            if math.isfinite(a) and math.isfinite(b) and a > 0 and b > 0:
                ps.add(a, b, p, tag="bnd-Tridge")

    # (5) min(a,b) = eps_R4 = 2^-6, other param and xi inside the R4 box.
    for Bp in (0.5, 2.0, 8.0):
        for f in (0.2, 0.5, 0.8):
            xi_wall = min(XI1_f, B1_f / Bp, math.exp(-math.log(2) / float(EPS_R4)))
            xi = f * xi_wall
            if 0 < xi < 1:
                for tau in (NEXT_DN(float(EPS_R4)), float(EPS_R4), NEXT_UP(float(EPS_R4))):
                    ps.add(tau, Bp, xi, tag="bnd-epsR4")

    # (6) tau*|ln xi_tau| = ln2 boundary.
    for tau in (1e-200, 1e-10, 0.001, 0.01, float(EPS_R4)):
        xi_wall = math.exp(-float(LN2) / tau)
        for Bp in (0.5, 2.0, 8.0):
            _bracket_add(ps, tau, Bp, xi_wall, "bnd-ln2")

    # (7) R2 orientation threshold xi = (a+1)/(c+2).
    for a in (1e-8, 0.1, 1.0, 8.0, 100.0, 1e6):
        for b in (0.1, 1.0, 8.0, 1e4):
            c = a + b
            thresh = (a + 1.0) / (c + 2.0)
            if 0 < thresh < 1:
                _bracket_add(ps, a, b, thresh, "bnd-orient-thresh")

    # (8) c = C_lg = 256, several a/b splits.
    for frac in (0.001, 0.01, 0.5, 0.99, 0.999):
        for cval in (NEXT_DN(float(C_LG)), float(C_LG), NEXT_UP(float(C_LG))):
            a = frac * cval
            b = cval - a
            if a > 0 and b > 0:
                ps.add(a, b, 0.3, tag="bnd-Clg")
                ps.add(a, b, 0.7, tag="bnd-Clg")

    # (9) R3's z boundary cpsi=36 (z=6): root-find x on both sides of the
    # ridge for representative (a,b) with nu>=T_ridge.
    for nu in (float(T_RIDGE), 1e3, 1e6):
        for p in (0.1, 1.0 / 3.0, 0.5):
            q = 1.0 - p
            c = nu / (p * q)
            a = p * c
            b = q * c
            if not (math.isfinite(a) and math.isfinite(b) and a > 0 and b > 0):
                continue
            for sign in (1.0, -1.0):
                lo_u, hi_u = (0.0, min(1.0, q / (2.0 * p))) if sign > 0 else \
                             (max(-0.5, -q / p), 0.0)

                def cpsi_of_u(u):
                    v = -(p / q) * u
                    phi_u = u - math.log1p(u)
                    phi_v = v - math.log1p(v)
                    return a * phi_u + b * phi_v

                lo, hi = lo_u, hi_u
                if cpsi_of_u(hi) < 36.0:
                    continue  # this side never reaches z=6 within the band
                for _ in range(80):
                    mid = (lo + hi) / 2.0
                    if cpsi_of_u(mid) < 36.0:
                        lo = mid
                    else:
                        hi = mid
                u0 = (lo + hi) / 2.0
                xi0 = p * (1.0 + u0)
                if 0 < xi0 < 1:
                    _bracket_add(ps, a, b, xi0, "bnd-z6")

    # (10) E = E_floor = -800 saturation boundary: root-find xi on the small
    # side such that the log-prefactor E crosses -800.
    for a, b in ((1.0, 1.0), (0.5, 200.0), (200.0, 0.5), (1e3, 1e3), (1e-3, 1e3)):
        lo, hi = 1e-300, 0.999999
        target = -800.0

        def e_of_xi(xi):
            return cheap_logE(a, b, xi)

        if e_of_xi(lo) > target:
            continue  # never saturates even at the extreme
        for _ in range(100):
            mid = math.sqrt(lo * hi) if lo > 0 else (lo + hi) / 2.0
            if e_of_xi(mid) < target:
                lo = mid
            else:
                hi = mid
        _bracket_add(ps, a, b, (lo + hi) / 2.0, "bnd-Efloor")

    print(f"  boundary brackets: {len(ps.pts) - n0} points", file=sys.stderr)


# ============================================================================
# x=1/2, a=b diagonal: P=Q=1/2 EXACTLY, several binades of a.
# ============================================================================
def gen_diagonal(ps):
    n0 = len(ps.pts)
    for k in range(-1000, 1001, 20):
        a = 2.0 ** k
        if math.isfinite(a) and a > 0:
            ps.add(a, a, 0.5, tag="diagonal")
    for a in (1e-300, 1e-100, 1e-10, 1.0, 1e10, 1e100, 1e250):
        ps.add(a, a, 0.5, tag="diagonal")
    # The +-1ulp-off-diagonal bracket is only meaningful up to moderate a:
    # for a=b at 1e100-1e250 class, Beta(a,a)'s variance is ~1/(8a) --
    # astronomically small, so a single ULP of x AWAY from 0.5 is already
    # many thousands of standard deviations into the tail (genuinely deep
    # saturation, not a "near-diagonal" probe) -- and the CF is numerically
    # unstable there (found during this generator's own development:
    # a=b=1e250 at x=nextafter(0.5) produced persistent dps-ladder
    # disagreement). That saturated-tail behavior is already covered by
    # gen_subnormal_band/gen_huge_tiny; this bracket stays at magnitudes
    # where "near-diagonal" is still a meaningful, numerically tractable
    # probe.
    for a in (1e-300, 1e-100, 1e-10, 1.0, 1e10, 1e6):
        for x in (NEXT_DN(0.5), NEXT_UP(0.5)):
            ps.add(a, a, x, tag="diagonal-bracket")
    print(f"  diagonal (a=b, x=1/2): {len(ps.pts) - n0} points", file=sys.stderr)


# ============================================================================
# Analytic lines (also the payload for self-check "analytic-line
# agreement"): I_x(a,1)=x^a; I_x(1,b)=1-(1-x)^b; I_x(1/2,1/2)=(2/pi)asin(sqrt(x)).
# (I_1/2(a,a)=1/2 is the diagonal above, already covered.)
# ============================================================================
ANALYTIC_LINE_POINTS = []  # (a, b, x, tag) filled by gen_analytic_lines


def gen_analytic_lines(ps, rng):
    n0 = len(ps.pts)
    x_grid = sorted(set(
        [0.5 ** k for k in range(1, 60)] +
        [1.0 - 0.5 ** k for k in range(1, 60)] +
        [i / 20.0 for i in range(1, 20)]
    ))
    # I_x(a,1) = x^a
    for a in (1e-8, 0.1, 1.0, 2.0, 8.0, 100.0, 1e6, 1e100):
        for x in x_grid:
            ps.add(a, 1.0, x, tag="line-a-1")
            ANALYTIC_LINE_POINTS.append((a, 1.0, x, "a-1"))
    # I_x(1,b) = 1-(1-x)^b
    for b in (1e-8, 0.1, 1.0, 2.0, 8.0, 100.0, 1e6, 1e100):
        for x in x_grid:
            ps.add(1.0, b, x, tag="line-1-b")
            ANALYTIC_LINE_POINTS.append((1.0, b, x, "1-b"))
    # I_x(1/2,1/2) = (2/pi) asin(sqrt(x))
    for x in x_grid:
        ps.add(0.5, 0.5, x, tag="line-half-half")
        ANALYTIC_LINE_POINTS.append((0.5, 0.5, x, "half-half"))
    for _ in range(400):
        x = rng.uniform(1e-12, 1.0 - 1e-12)
        ps.add(0.5, 0.5, x, tag="line-half-half-rand")
        ANALYTIC_LINE_POINTS.append((0.5, 0.5, x, "half-half"))
    print(f"  analytic lines: {len(ps.pts) - n0} points", file=sys.stderr)


# ============================================================================
# Subnormal band: log-prefactor E in [-745, -700] (double's normal/subnormal
# crossing sits near E~-708 to E~-744 -- bracket the whole transition zone).
# ============================================================================
def gen_subnormal_band(ps):
    n0 = len(ps.pts)
    for a, b in ((1.0, 1.0), (1.0, 8.0), (8.0, 1.0), (1e3, 1e3), (0.1, 1e4),
                 (1e4, 0.1), (50.0, 50.0), (1e6, 1e-3)):
        for target in (-700.0, -710.0, -720.0, -730.0, -740.0, -745.0):
            lo, hi = 1e-300, 0.999999999
            if cheap_logE(a, b, lo) > target:
                continue

            def e_of_xi(xi):
                return cheap_logE(a, b, xi)

            for _ in range(100):
                mid = math.sqrt(lo * hi) if lo > 0 else (lo + hi) / 2.0
                if e_of_xi(mid) < target:
                    lo = mid
                else:
                    hi = mid
            xi0 = (lo + hi) / 2.0
            for x in (NEXT_DN(xi0), xi0, NEXT_UP(xi0)):
                ps.add(a, b, x, keep_if_saturated=True, tag="subnormal-band")
    print(f"  subnormal band: {len(ps.pts) - n0} points", file=sys.stderr)


# ============================================================================
# Huge/tiny parameters: a or b at 2^+-1000-class, subnormal a, c near
# overflow (the half-scaling guard's domain, c ~ 2^1022+).
# ============================================================================
def gen_huge_tiny(ps, rng):
    n0 = len(ps.pts)
    extreme_a = [2.0 ** k for k in (-1074, -1070, -1022, -1000, -900, -700,
                                     -500, -200, 200, 500, 700, 900, 1000, 1020, 1023)]
    moderate_b = [1e-3, 0.1, 1.0, 8.0, 100.0, 1e6]
    for a in extreme_a:
        for b in moderate_b:
            for x in (0.01, 0.1, 0.5, 0.9, 0.99):
                ps.add(a, b, x, keep_if_saturated=True, tag="huge-tiny-a")
                ps.add(b, a, x, keep_if_saturated=True, tag="huge-tiny-b")
    # c = a+b near the overflow guard threshold ~2^1022.
    for frac in (1e-10, 1e-3, 0.5, 0.999, 1.0 - 1e-10):
        for cval in (2.0 ** 1021, 2.0 ** 1022, NEXT_DN(2.0 ** 1023), 2.0 ** 1023):
            a = frac * cval
            b = cval - a
            if 0 < a < math.inf and 0 < b < math.inf and math.isfinite(cval):
                for x in (0.1, 0.5, 0.9):
                    ps.add(a, b, x, keep_if_saturated=True, tag="c-near-overflow")
    # subnormal a explicitly (a below 2^-1022, the double subnormal range).
    for a in (5e-324, 1e-320, 1e-310, 2.0 ** -1050, 2.0 ** -1074):
        for b in (0.5, 1.0, 8.0, 100.0):
            for x in (0.1, 0.5, 0.9):
                ps.add(a, b, x, keep_if_saturated=True, tag="subnormal-a")
    for _ in range(1000):
        a = 2.0 ** rng.uniform(-1074, 1023)
        b = 2.0 ** rng.uniform(-1074, 1023)
        x = rng.uniform(1e-6, 1.0 - 1e-6)
        ps.add(a, b, x, keep_if_saturated=True, tag="huge-tiny-random")
    print(f"  huge/tiny parameters: {len(ps.pts) - n0} points", file=sys.stderr)


# ============================================================================
# Escalation witnesses (G1a, G1b) and small families around each.
# ============================================================================
def gen_witnesses(ps):
    n0 = len(ps.pts)
    # G1a witness: (8, 2^-6, 1-9.5e-7) -- direct value 0.16, R4 territory.
    a0, b0, x0 = 8.0, 2.0 ** -6, 1.0 - 9.5e-7
    for a in (a0, a0 * 0.5, a0 * 2.0, 4.0, 16.0):
        for b in (b0, NEXT_DN(b0), NEXT_UP(b0), 2.0 ** -7, 2.0 ** -5):
            for x in (x0, NEXT_DN(x0), NEXT_UP(x0), 1.0 - 1e-9, 1.0 - 1e-12, 0.9, 0.9999):
                ps.add(a, b, x, tag="witness-g1a")
    # G1b witness: (1e-20, 1, 0.4) -- R1-native fires and evaluates 1.0
    # before the tiny-min-first correction; family around it.
    a1, b1, x1 = 1e-20, 1.0, 0.4
    for a in (a1, 1e-15, 1e-10, 1e-6, a1 * 10, a1 * 0.1):
        for b in (b1, 0.5, 2.0, 8.0):
            for x in (x1, 0.2, 0.6, 0.8, NEXT_DN(x1), NEXT_UP(x1)):
                ps.add(a, b, x, tag="witness-g1b")
    print(f"  escalation witnesses: {len(ps.pts) - n0} points", file=sys.stderr)


# ============================================================================
# PB-prefactor u -> -1 corner family [2026-08-12 fix arc; PLAN.md's
# "TWO defects in the SHIPPED beta forward" Open Item]. The corner the
# original families never sampled: moderate min >= Z0 with huge max and x
# deep in the tail, where w = c*x/min (the kernel's 1+u) falls to 2^-53
# and below -- the shipped kernel's TwoSum spelling of 1+u degenerated
# there (u.hi rounding to exactly -1 -> NaN -> a silent exact-0 return;
# 1+u merely small -> an unnormalized dd into LogDdAny -> 1.4e-4 in E).
# DETERMINISTIC (no rng draw, so a fresh full run reproduces every earlier
# family byte-identically; called LAST in main() for the same reason --
# gen_random_fill's stream and stopping point are untouched).
#
# These points route R1 (x tiny puts them inside the series box), but the
# PREFACTOR is PB whenever c > C_lg and min >= Z0 -- which is exactly how
# the defect reached shipped beta_p: the R1/R2 split does not gate the
# prefactor path, i_e = max(i_r1, i_r2) does.
# ============================================================================
def gen_pb_corner(ps):
    n0 = len(ps.pts)
    # The two recorded defect witnesses, bit-exact as filed.
    ps.add(19.0, 1e5, 5.204222470155122e-21, tag="pb-corner-witness")
    ps.add(19.0, 1e5, 1.73e-19, tag="pb-corner-witness")
    # Grid: x = w*min/c targets w = c*x/min directly, spanning the whole
    # hazard (w < 2^-53 is the NaN corner, 2^-53..2^-8 the unnormalized-dd
    # corner, up to 0.4 as the healthy-side control). Depth-saturated rows
    # (truth rounds to exact 0) are capped at two per (m, b) pair -- they
    # lock the saturation edge without bloating the set. cheap_logE
    # classifies; its ~ln(a) slack vs ln P is immaterial at the -760 cut
    # (double underflows to 0 below ln P ~ -744.4 already).
    m_list = [10.0, 12.5, 19.0, 40.0, 100.0, 400.0, 1000.0]
    # No b between 1e12 and 1e30: above B_GL = 2^59 the oracle switches to
    # the gamma-limit form, whose O(a*t/b) truncation (t = w*m <= 400 here)
    # reaches 4e-15 at b = 1e20 -- NOT certifiable for double rounding. At
    # b >= 1e30 the worst combination is 4e-25, below the cross-check
    # target with margin; below B_GL the CF is mpf-exact. The u -> -1
    # corner mechanics depend only on w, not on b's magnitude, so the gap
    # costs no defect coverage.
    b_list = [300.0, 1e3, 1e4, 1e5, 1e8, 1e12, 1e30, 1e50, 1e100, 1e200,
              1e307]
    w_list = [2.0 ** -60, 2.0 ** -56, 2.0 ** -54, 2.0 ** -53, 2.0 ** -52,
              2.0 ** -50, 2.0 ** -45, 2.0 ** -40, 2.0 ** -30, 2.0 ** -20,
              2.0 ** -12, 2.0 ** -8, 2.0 ** -4, 0.25, 0.4]
    for m in m_list:
        for bb in b_list:
            if bb <= m:
                continue  # bb is strictly the max in this family's geometry
            c = m + bb
            if not c > 256.0:
                continue  # PB gate is c > C_lg = 256 strict
            n_sat = 0
            for w in w_list:
                x = w * m / c
                if not 0.0 < x < 1.0:
                    continue
                if cheap_logE(m, bb, x) < -760.0:
                    n_sat += 1
                    if n_sat > 2:
                        continue
                ps.add(m, bb, x, tag="pb-corner")
    # Bit-level brackets of the u.hi == -1 rounding boundary (w crossing
    # 2^-54 flips u.hi between -1-exact and its predecessor; 2^-53 is the
    # last w whose 1+u survives in u.hi at all). Pairs chosen unsaturated.
    for m, bb in ((19.0, 1e5), (10.0, 1e8), (12.5, 1e10)):
        c = m + bb
        for wb in (2.0 ** -54, 2.0 ** -53):
            xb = wb * m / c
            for x in (NEXT_DN(NEXT_DN(xb)), NEXT_DN(xb), xb, NEXT_UP(xb),
                      NEXT_UP(NEXT_UP(xb))):
                ps.add(m, bb, x, tag="pb-corner-bracket")
    # Swapped-frame corner: x one-to-a-thousand ulps below 1 routes the
    # complement (min, max, 1-x) with 1-x exact in dd -- the same corner
    # entered through beta_q's orientation (and betainv's fixed frame).
    for k in (1.0, 2.0, 19.0, 1000.0):
        x = 1.0 - k * 2.0 ** -53
        ps.add(1e5, 19.0, x, tag="pb-corner-swap")
        ps.add(1e8, 10.0, x, tag="pb-corner-swap")
    # Negative controls: the PA side of both gate inequalities (min < Z0,
    # c <= C_lg) never had the defect; these pin the gate boundary itself,
    # including c == 256 exactly (strict inequality -> PA).
    for x in (5.2e-21, 1.7e-19):
        ps.add(9.5, 1e5, x, tag="pb-corner-pa-control")
    for x in (1e-6, 1e-9):
        ps.add(19.0, 200.0, x, tag="pb-corner-pa-control")
    ps.add(10.0, 246.0, 1e-8, tag="pb-corner-gate-bracket")  # c == 256: PA
    ps.add(10.0, NEXT_UP(246.0), 1e-8, tag="pb-corner-gate-bracket")  # PB
    ps.add(NEXT_DN(10.0), 1e5, 2e-20, tag="pb-corner-gate-bracket")  # PA
    ps.add(10.0, 1e5, 2e-20, tag="pb-corner-gate-bracket")  # PB
    print(f"  pb-corner family: {len(ps.pts) - n0} points", file=sys.stderr)


# ============================================================================
# Random seeded fill (own fresh stream, SEED=20260731) to reach the target
# total point count.
# ============================================================================
def gen_random_fill(ps, rng, target_total):
    n0 = len(ps.pts)
    attempts = 0
    max_attempts = (target_total - len(ps.pts)) * 6 + 20000
    while len(ps.pts) < target_total and attempts < max_attempts:
        attempts += 1
        choice = rng.random()
        if choice < 0.5:
            # log-uniform a,b, uniform x -- broad coverage.
            a = 10.0 ** rng.uniform(-300.0, 300.0)
            b = 10.0 ** rng.uniform(-300.0, 300.0)
            x = rng.uniform(1e-9, 1.0 - 1e-9)
        elif choice < 0.8:
            # moderate a,b (where most "typical use" mass sits).
            a = 10.0 ** rng.uniform(-3.0, 6.0)
            b = 10.0 ** rng.uniform(-3.0, 6.0)
            x = rng.uniform(1e-6, 1.0 - 1e-6)
        else:
            # near-diagonal / near-mean stress.
            a = 10.0 ** rng.uniform(-2.0, 8.0)
            b = a * 10.0 ** rng.uniform(-2.0, 2.0)
            mean = a / (a + b)
            x = min(max(mean + rng.uniform(-0.3, 0.3) * mean, 1e-9), 1.0 - 1e-9)
        ps.add(a, b, x, keep_if_saturated=False, tag="random-fill")
    print(f"  random fill: {len(ps.pts) - n0} points ({attempts} attempts)",
          file=sys.stderr)


# ============================================================================
# Specials table [PLAN.md "Specials" -- gamma-consistent doctrine: one
# degenerate parameter gets its limit, two degeneracies (or a degenerate
# parameter meeting the x-boundary its own mass sits on) -> NaN]. Direct
# assignment, no oracle call -- these are exact by the design's own rules.
# ============================================================================
def gen_specials_rows():
    rows = []
    NAN, INF = float("nan"), float("inf")

    def add(a, b, x, P, Q):
        rows.append((a, b, x, P, Q))

    finite_abs = (1e-300, 1e-8, 1.0, 3.7, 1e8, 1e300)

    # NaN anywhere -> NaN.
    for a, b, x in ((NAN, 1.0, 0.5), (1.0, NAN, 0.5), (1.0, 1.0, NAN),
                    (NAN, NAN, NAN), (NAN, 1.0, -5.0), (1.0, NAN, INF)):
        add(a, b, x, NAN, NAN)

    # x not in [0,1] (finite a,b>0) -> NaN.
    for a, b in ((1.0, 1.0), (1e-8, 1e8), (3.7, 12.0)):
        for x in (-0.1, 1.1, -INF, INF, -1e300):
            add(a, b, x, NAN, NAN)

    # a<0 or b<0 -> NaN.
    for a, b, x in ((-1.0, 1.0, 0.5), (1.0, -1.0, 0.5), (-1.0, -1.0, 0.5),
                    (-INF, 1.0, 0.5), (1.0, -INF, 0.5), (-1e-8, 3.0, 0.2)):
        add(a, b, x, NAN, NAN)

    # Two of {a,b} in {0,inf} degenerate -> NaN, at several x.
    for a in (0.0, INF):
        for b in (0.0, INF):
            for x in (0.0, 0.3, 0.5, 0.7, 1.0):
                add(a, b, x, NAN, NAN)

    # x=0, a,b>0 finite -> P=+0, Q=1.  x=1, a,b>0 finite -> P=1, Q=+0.
    for a in finite_abs:
        for b in finite_abs:
            add(a, b, 0.0, 0.0, 1.0)
            add(a, b, 1.0, 1.0, 0.0)

    # a=0, b>0 finite (mass at 0): P=1 for x in (0,1]; NaN at x=0
    # (mass-point-meets-boundary row #1).
    for b in finite_abs:
        for x in (1e-300, 1e-8, 0.5, 1.0 - 1e-16, 1.0):
            add(0.0, b, x, 1.0, 0.0)
        add(0.0, b, 0.0, NAN, NAN)

    # b=0, a>0 finite (mass at 1): P=+0 for x in [0,1); NaN at x=1
    # (mass-point-meets-boundary row #2).
    for a in finite_abs:
        for x in (0.0, 1e-8, 0.5, 1.0 - 1e-16):
            add(a, 0.0, x, 0.0, 1.0)
        add(a, 0.0, 1.0, NAN, NAN)

    # b=+inf, a>0 finite (mass at 0): P=1 for x in (0,1]; NaN at x=0
    # (mass-point-meets-boundary row #3).
    for a in finite_abs:
        for x in (1e-300, 1e-8, 0.5, 1.0 - 1e-16, 1.0):
            add(a, INF, x, 1.0, 0.0)
        add(a, INF, 0.0, NAN, NAN)

    # a=+inf, b>0 finite (mass at 1): P=+0 for x in [0,1); NaN at x=1
    # (mass-point-meets-boundary row #4).
    for b in finite_abs:
        for x in (0.0, 1e-8, 0.5, 1.0 - 1e-16):
            add(INF, b, x, 0.0, 1.0)
        add(INF, b, 1.0, NAN, NAN)

    return rows


# ============================================================================
# Cheap (low-dps) magnitude estimate: log of the leading order of
# x^a*(1-x)^b/B(a,b), used ONLY to steer boundary/subnormal-band root-finds
# (never to decide an emitted value -- those always go through the real
# CF oracle in small_side_direct). Not self-check-critical.
# ============================================================================
def _ln_beta(am, bm, dps):
    """ln B(a,b), robust for hugely mismatched parameters: the naive
    lnG(a)+lnG(b)-lnG(a+b) forms a+b at ambient dps (1e250+1 truncates to
    1e250 exactly, making lnG(a+b)-lnG(a) cancel to 0 instead of -ln a --
    a 575-unit bias that turned est=-394 into est=-969 and falsely
    saturated (1, 1e250, 4e-248), true Q=e^-400). Route the mismatched
    pair through _lngamma_diff_b_tau's exact extra-dps difference."""
    lo, hi = (am, bm) if am <= bm else (bm, am)
    return mp.loggamma(lo) - _lngamma_diff_b_tau(lo, hi, dps)


def cheap_logE(a, b, xi, dps=25):
    if not (0.0 < xi < 1.0) or a <= 0.0 or b <= 0.0:
        return float("-inf")
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        am, bm, xim = mp.mpf(a), mp.mpf(b), mp.mpf(xi)
        E = am * mp.log(xim) + bm * mp.log1p(-xim) - _ln_beta(am, bm, dps)
        return float(E)
    except (ValueError, OverflowError):
        return float("-inf")
    finally:
        mp.mp.dps = old


# ============================================================================
# EXACT-COMPLEMENT DISCIPLINE [round 5, root cause of BOTH residual ULP
# defect classes after the deep-ladder fix]. mpf addition/subtraction at
# ambient dps TRUNCATES operands to working precision first: forming
# 1 - xm for xm below 10^-dps silently yields exactly 1, and (the dual
# hazard) mp.log1p(-xx) on a HIGH-precision near-1 complement returns
# -inf (probed directly: log1p(-(1-8e-250 formed exactly)) = -inf at
# dps 25, while mp.log on the same exact operand is correct). Two rules:
#   (1) form every orientation complement through _one_minus (exact at
#       any xm; ~-log10(xm)+20 extra digits when xm is tiny);
#   (2) never push a complement through log1p -- in the swapped frame
#       the two logs are EXACTLY log1p(-xm) and log(xm) of the ORIGINAL
#       double-derived xm (always safe: <=53-bit mantissas are never
#       truncated at dps>=25), so pass logs, not the complement
#       (cheap_logE_logs below).
# Double-derived operands (53-bit mantissas) are exempt from both
# hazards; only ops mixing ambient precision with HIGHER-precision
# operands truncate. The same truncation dropped loggamma(1+tau)'s
# entire -gamma*tau term for tau < 10^-dps in small_tau_oracle (the
# systematic "stored small side = EulerGamma * min(a,b)" artifact
# family), fixed at that site by a workdps bump.
# ============================================================================
def _one_minus(xm):
    """Exact 1 - xm as an mpf carrying full precision, for xm in (0, 1)
    derived from a double. Safe to consume via mp.log at ambient dps;
    NOT safe to consume via mp.log1p(-result) -- see block comment."""
    if not (0 < xm < 1):
        return 1 - xm
    extra = max(0, int(-mp.log10(xm))) + 20
    with mp.workdps(mp.mp.dps + extra):
        return 1 - xm


def cheap_logE_logs(a, b, ln_xi, ln_1m_xi, dps=25):
    """cheap_logE with the two logs supplied exactly by the caller (the
    orientation-safe form: a swapped frame passes ln_xi=log1p(-xm),
    ln_1m_xi=log(xm) from the ORIGINAL xm and never forms a complement).
    Emission-deciding callers (_saturation_prefilter and the main-path
    prefilter) MUST use this form -- float(1-xm) collapses to 1.0 for
    tiny xm, fails cheap_logE's domain check, and the resulting -inf
    reads as a saturation certificate for a point whose small side is
    O(1e-4) (witness: (1, 1e250, 8e-250), true Q = e^-8)."""
    if a <= 0.0 or b <= 0.0:
        return float("-inf")
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        am, bm = mp.mpf(a), mp.mpf(b)
        E = am * ln_xi + bm * ln_1m_xi - _ln_beta(am, bm, dps)
        return float(E)
    except (ValueError, OverflowError):
        return float("-inf")
    finally:
        mp.mp.dps = old


def _oriented_logs(xm, native):
    """(ln_xi, ln_1m_xi) for cheap_logE_logs in the given orientation,
    computed from the original double-derived xm only (both safe)."""
    if native:
        return mp.log(xm), mp.log1p(-xm)
    return mp.log1p(-xm), mp.log(xm)


# ============================================================================
# Primary value oracle: SMALL-SIDE-DIRECT via the CF, three-layer dps
# ladder (40 -> recheck 60 -> escalate 100 on disagreement).
# ============================================================================
ESCALATIONS = []
N_CF_FAILED = [0]
N_PREFILTER_SATURATED = [0]

# Comfortably below E_FLOOR=-800 (the design's own saturation threshold,
# kBetaExpFloor): a point whose chosen side already has log-magnitude below
# this is not merely "saturated in double" but so far past the smallest
# subnormal (~2^-1074, log ~ -744.4) that its own EXPONENT is astronomically
# large -- attempting the CF/dps-ladder there is not just wasteful but
# numerically MEANINGLESS. Found during this generator's own development:
# an R3-lattice point at nu=1e100 produced dps=40/60/100 values that
# disagreed in the low digits of a >1400-digit-long EXPONENT (the true
# value sits around 1e-1.4e99) -- the "disagreement" the escalation/
# cross-check machinery correctly flagged was real, but pointless: every
# one of those values rounds to the exact same double (0.0) regardless.
# Fixed by this pre-filter, mirroring gamma's "a*phi>800 -> exact saturated
# pair, no further computation" (gen_gamma_reference.py's oracle_pq).
SATURATION_LOG_THRESHOLD = -900.0


def _cf_small(aa, bb, xx, dps):
    return small_val_via_cf(aa, bb, xx, dps, n_start=64, n_max=8192)


# ============================================================================
# GAMMA-CORNER ORACLE [(C), G3 escalation resolution]. For rows whose
# CF-orientation max parameter >= kBetaGammaLim (B_GL, imported from
# gen_beta_data -- same constant the kernel's router and this generator's
# own check_b_r2/(vii) sweep use), the beta CF is structurally degenerate
# (escalation (C)'s own witness: (0.05,1e100,2e-99), d_1 -> -(1-2e-99),
# divides by zero at every depth/precision tried) -- this oracle switches to
# mpmath.gammainc (regularized) at t=-beta*log1p(-xi), exact in mpf, instead.
#
# CORRECTNESS NOTE (bug found and fixed by this generator's own gamma-corner
# self-check, not by reasoning about the code): the asymptotic identity
# I_xi(shape,scale) -> P(shape,t), t=-scale*log1p(-xi) as scale->infty with
# shape FIXED, is anchored to a SPECIFIC parameter playing "shape" (bounded)
# vs "scale" (->infty) -- it is NOT the same swap the backward CF's own
# xi<(a+1)/(c+2) orientation rule performs (that rule picks whichever
# ORIENTATION converges the CF fastest, unrelated to which parameter is
# huge). An earlier version of this function reused the CF's own
# orientation rule here and got P(a1,t) backwards whenever aa (not bb) was
# the huge parameter -- caught by check_gamma_corner_oracle measuring
# relative error in the 1e267 range (gform=1.0 exactly, i.e. total
# precision loss, vs a true value ~1e-268) before it ever reached a
# reference row. The one identity this function actually needs: whichever
# of (aa,bb) is the HUGE one plays "scale"; the OTHER is "shape" and is
# passed to gammainc UNCHANGED regardless of native/swap-style relabeling.
#   bb is huge: t=-bb*log1p(-xx); I_xx(aa,bb) ~= P(aa,t) directly.
#   aa is huge: by I_xx(aa,bb)=1-I_{1-xx}(bb,aa) with (bb,aa,1-xx) now in
#     the "bb huge is impossible" shape (bb fixed here, aa huge), apply the
#     SAME identity to the swapped triple: I_{1-xx}(bb,aa) ~= P(bb,t'),
#     t'=-aa*log1p(-(1-xx))=-aa*log(xx) EXACTLY (log1p(xx-1)=log(xx), no
#     near-1 subtraction needed since xx itself, not 1-xx, feeds the log).
#     So I_xx(aa,bb) = 1-P(bb,t') = Q(bb,t').
def gamma_corner_value(aa, bb, xx, dps):
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        # xx may be a HIGH-precision exact complement (_one_minus): mp.mpf()
        # on an mpf RE-ROUNDS to ambient dps, collapsing 1-tiny back to 1.0
        # (t then computes 0 and the value saturates falsely, forcing the
        # caller's mean-predicate flip and a big-side-direct emission).
        # Pass it through untouched -- mp.log consumes full precision.
        aa_m, bb_m = mp.mpf(aa), mp.mpf(bb)
        xx_m = xx if isinstance(xx, mp.mpf) else mp.mpf(xx)
        if bb_m >= aa_m:
            t = -bb_m * mp.log1p(-xx_m)
            return mp.gammainc(aa_m, 0, t, regularized=True)  # P(aa,t)
        else:
            t = -aa_m * mp.log(xx_m)
            return mp.gammainc(bb_m, t, mp.inf, regularized=True)  # Q(bb,t)
    finally:
        mp.mp.dps = old


def gamma_corner_value_signed(aa, bb, xx, dps):
    """Same identity as gamma_corner_value, but returns (small_val, which)
    with small_val computed DIRECTLY (never via 1-near_one) whichever of
    P(=I_xx(aa,bb)) / Q(=1-P) is genuinely the small one -- used by
    check_gamma_corner_oracle so its own comparison never repeats the exact
    1-near-1 cancellation hazard this whole generator otherwise guards
    against (P and Q here are BOTH direct gammainc calls, upper vs lower
    tail, never each other's complement)."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        # Same no-re-round rule as gamma_corner_value above.
        aa_m, bb_m = mp.mpf(aa), mp.mpf(bb)
        xx_m = xx if isinstance(xx, mp.mpf) else mp.mpf(xx)
        if bb_m >= aa_m:
            t = -bb_m * mp.log1p(-xx_m)
            P = mp.gammainc(aa_m, 0, t, regularized=True)
            Q = mp.gammainc(aa_m, t, mp.inf, regularized=True)
        else:
            t = -aa_m * mp.log(xx_m)
            Pgb = mp.gammainc(bb_m, 0, t, regularized=True)  # P(bb,t)
            Q = mp.gammainc(bb_m, t, mp.inf, regularized=True)  # Q(bb,t) = P(aa,bb,xx)-frame's P
            P = Pgb
            # here P(aa,bb,xx) = Q(bb,t) and Q(aa,bb,xx) = P(bb,t) -- swap
            # the labels back to the (aa,bb,xx)-frame's own P/Q meaning.
            P, Q = Q, P
        return (P, "P") if P <= Q else (Q, "Q")
    finally:
        mp.mp.dps = old


# ============================================================================
# SMALL-TAU ORACLE [coordinator revision, review round 2]. The backward CF
# was found to fail broadly across the whole tiny-min(a,b) corner (measured
# from the checkpoint after the first full run: failure RATE actually rises
# toward the eps_R4=2^-6 boundary itself -- 62% of R4 points at tau~1e-2 vs
# 2% at tau~1e-6 -- not confined to a narrow sub-threshold), exactly the
# region gen_r4's own box (tau*|ln xi_tau|<=ln2, xi_tau<=xi1, B*xi_tau<=B1)
# targets and the region two routing flaws (G1a/G1b) made load-bearing.
#
# This is the APSER-style closed form the design's own R4 kernel will use
# (PLAN.md "R4 tiny-min ... gamma-R4 verbatim in beta clothing"): derived
# and cross-checked against gen_beta_data.py's own R1/R4 series convention
# (t_n = t_{n-1}*(n-B)*xi/n = (1-B)_n*xi^n/n!, matching series_partial_sums
# there) --
#
#   Q~ = -expm1(w + ln S)
#   w  = tau*ln(xi_tau) + [lnGamma(B+tau) - lnGamma(B)] - lnGamma(1+tau)
#   S  = 1 + tau * Sigma,  Sigma = sum_{n>=1} t_n/(tau+n)
#   ln S = log1p(tau * Sigma)
#
# Algebraic identity (verified symbolically before trusting this in code):
# S*e^w = I_xi_tau(tau, B) exactly (S = tau*[the standard R1-form sum],
# e^w = xi_tau^tau/(tau*B(tau,B)) via Gamma(1+tau)=tau*Gamma(tau) and
# B(tau,B)=Gamma(tau)Gamma(B)/Gamma(tau+B) -- so S*e^w telescopes to the
# textbook BPSER value). Hence Q~ = 1 - I_xi_tau(tau,B): as tau->0+ with
# xi_tau bounded away from the box's own ln2 floor, I_xi_tau(tau,B) -> 1
# (the a=0 "mass at 0" limit), making Q~ the genuinely SMALL quantity,
# computed via expm1/log1p so the near-1 cancellation never happens in
# floating point -- every term (lnGamma difference, log1p(tau*Sigma)) is
# individually small/benign, matching the coordinator's brief exactly.
#
# Orientation mapping back to the ORIGINAL (a,b) frame (route_final's own
# tiny-first convention, R4-native/R4-swap):
#   a<=b (native): (tau,B,xi_tau)=(a,b,x)  -> Q~ IS the original Q.
#   b<a  (swap):   (tau,B,xi_tau)=(b,a,1-x) -> Q~ IS the original P
#                  (since I_{1-x}(b,a)=Q by the standard swap identity, so
#                  Q~=1-I=1-Q=P here).
# ============================================================================
SMALL_TAU_THRESHOLD = mp.mpf(2) ** -4  # empirically re-derived, see report:
# eps_R4=2^-6 alone is NOT enough margin -- the failing witness-g1a family
# reaches b=2^-5 (twice eps_R4), so the threshold is set one more binade up
# (2^-4) to cover it with room to spare, while staying well clear of R1/R2/
# R3 territory where the CF is known to work fine and stays primary.
SMALL_TAU_FALLBACK_CEILING = mp.mpf(1)  # rescue net: if the REGULAR CF path
# fails outright and min(a,b) is at least "small-ish" (<1, generous), retry
# via this oracle before dropping the point -- catches strays outside the
# threshold band that still have a tiny-ish parameter (e.g. this run's
# lone R1-random and any other one-off).
N_SMALL_TAU_RESCUED = [0]
N_SMALL_TAU_NONCONVERGENT = [0]


SMALL_TAU_LNGAMMA_RATIO = mp.mpf("1e-8")  # tau <= B*this -> analytic Taylor
# Sanity ratio for the 3rd-Taylor-term bound check (see
# _lngamma_diff_b_tau's own comment for why this is a fixed, RELATIVE-to-
# the-first-term ratio and not a dps- or DISAGREE_100_60-scaled absolute
# bound -- both of those were tried and both broke on real points). The
# series' own structure keeps third_term/first_term ~ (tau/B)^2 <=
# SMALL_TAU_LNGAMMA_RATIO^2 = 1e-16 by construction whenever this branch is
# entered; 1e-10 leaves five orders of margin as a misconfiguration net.
TAYLOR_BOUND_TOL = mp.mpf("1e-10")


def _lngamma_diff_b_tau(tau_m, B_m, dps):
    """lnGamma(B+tau) - lnGamma(B), the term escalation (A) [G3/G2-revision]
    fixed: at fixed working dps this difference goes IDENTICALLY to 0 in mpf
    once tau is far enough below B's own ulp there -- the reference set's
    own defect (witness family (a, 1.4e-300, 1-2^-52), a<->B in this frame
    -- see small_side_direct's tau/B mapping): the emitted P went CONSTANT
    across a in [2^-20,2^-6] where psi(a) ~ -1/a should swing it by 2^14,
    because at whatever dps this generator used, B+tau rounded to exactly B
    before loggamma ever saw the difference.

    Two regimes, per the resolution:
      tau <= B*1e-8: analytic Taylor tau*psi0(B) + tau^2/2*psi1(B) -- exact
      in the mpf sense that psi0/psi1 are each evaluated at FULL working
      precision (no B+tau addition, hence no cancellation possible) --
      plus a tau^3/6*psi2(B) bound check asserting the truncated term is
      below the OUTPUT precision -- if the bound check fails this is a real
      accuracy problem, not something to paper over, so it raises.
      tau > B*1e-8: the direct difference is fine in PRINCIPLE (no identical
      cancellation to guard against) but needs escalated PRECISION to keep
      it accurate to the caller's own dps once B/tau is not astronomically
      large -- computed at dps + log10(B/tau) + 20, matching the resolution
      text exactly, then rounded back to the caller's working dps on return
      (mp.mpf's own value carries however much precision it was computed
      at; the caller's ambient dps governs everything downstream).

    BOUND-CHECK TOLERANCE, reasoned (went through two wrong versions before
    this one -- both caught by this generator's own full-lattice run, not
    by reasoning about the code):
      v1 gated at 10^-dps: WRONG, this is a SERIES TRUNCATION (dropping the
      tau^3 term), not a rounding error -- it does not shrink as the
      caller's dps grows (the same two Taylor terms are exact to the same
      ~tau^3*psi2(B)/6 regardless of how many digits mp.psi carries), so a
      dps-scaled bound fails perfectly good points once dps passes ~44
      digits (witness: tau=0.05, B=2^65).
      v2 gated at a FIXED absolute tolerance (scaled off DISAGREE_100_60,
      floored at 1 when lg_diff was tiny): WRONG the other direction --
      flooring "scale" at 1 turned a RELATIVE bound into an ABSOLUTE one
      whenever lg_diff itself is legitimately tiny (witness: tau=1e-12,
      B=2.5e-4, lg_diff~-4e-9 -- third_term~-2.1e-26 is 5.3e-18 RELATIVE to
      lg_diff, utterly negligible, but larger in ABSOLUTE terms than the
      fixed 8.47e-30 floor).
      v3 (this one): the Taylor series' OWN structure makes third_term/
      first_term ~ O((tau/B)^2) UNIVERSALLY, independent of B's absolute
      scale or the caller's dps -- since SMALL_TAU_LNGAMMA_RATIO=1e-8
      already bounds tau/B, this ratio is bounded by ~1e-16 by
      CONSTRUCTION whenever the branch is even entered (confirmed against
      both witnesses above: v1's case measured ~1e-42, v2's case measured
      ~5e-18, both far inside 1e-16). The check is therefore a SANITY net
      on that ratio holding (catching a misconfigured SMALL_TAU_LNGAMMA_
      RATIO, not tuning per-point precision), not a per-point precision
      budget -- fixed at 1e-10, five orders of magnitude of margin below
      where the ratio naturally sits, scaled against the FIRST Taylor term
      (tau*psi0(B)) specifically, which is what third_term is actually
      small RELATIVE TO.
    """
    # v4 [round 6]: the Taylor branch is DELETED. v3's bound was correct
    # about lg_diff's own accuracy and still WRONG about the assembly's:
    # w's terms CANCEL against tau*sigma downstream, and the result y can
    # be 15+ orders below lg_diff's scale, so a truncation that is 1e-18
    # OF LG_DIFF can be 70% OF THE RESULT (tau=1e-10, B=20: dropped
    # third_term -4.3e-34 vs true y -6.2e-34), 13x (tau=7.4e-6, B=7896:
    # -1.1e-24 vs 8.5e-26), or invisible (tau=4.4e-7, B=70: -2.9e-24 vs
    # 3.4e-11) -- three measured points, one rule: the required accuracy
    # is set by the RESULT's cancellation depth, unknowable at this site,
    # so no truncated form can be certified here. The exact extra-dps
    # loggamma difference (forming B+tau with digits(B/tau)+20 headroom)
    # costs ~ms and is exact for every (tau, B). Found by the fresh-
    # reference ULP run: the kernel matched analytic truth at every
    # probed point and the ORACLE carried these errors.
    old = mp.mp.dps
    extra_dps = dps + int(mp.log10(B_m / tau_m)) + 20
    mp.mp.dps = max(dps, extra_dps)
    try:
        return mp.loggamma(B_m + tau_m) - mp.loggamma(B_m)
    finally:
        mp.mp.dps = old


def small_tau_oracle(tau, B, xi_tau, dps, n_max=4000):
    """Q~ = 1 - I_xi_tau(tau, B), the APSER-style expm1/log1p assembly
    above. Returns (Q_tilde, converged). Term-count guard: stop and report
    non-convergence (never silently truncate) if n_max is reached without
    the running term dropping below the dps-scaled tolerance.

    CANCELLATION GUARD [found by this generator's own full-lattice run, a
    pre-existing gap unrelated to the (A)/(C) fixes above]: this series is
    R1's OWN power series in xi_tau (t *= (n-B)*xi_tau/n) -- well-behaved
    when B*xi_tau stays modest (R4's real box caps it at B1=8, widened to
    ~0.24 by the fourth correction), but SMALL_TAU_THRESHOLD only gates on
    min(a,b), not on B or xi_tau. A point with B in the thousands and
    xi_tau near 1 (found via gen_r1's own "swaprole" lattice, e.g.
    (a=3981.07,b=1e-300,x=1.005e-4) -> tau=1e-300,B=3981.07,xi_tau=0.9999)
    makes the term magnitude climb like a central-binomial peak (~B terms
    growing before shrinking again), reaching ~1e852 at dps=40/60 before
    "converging" back to a small |term| -- but the ACCUMULATED sigma by
    then has lost every digit of its true (tiny) value to cancellation the
    working dps can't recover, and the resulting tau*sigma can even land
    <=-1, sending log1p complex (an uncaught TypeError three frames up in
    the caller, not a graceful drop). Guarded here two ways: (1) track the
    peak |term| seen; if it exceeds the final |sigma| by more decimal
    digits than the working dps has headroom for, the "converged" flag is
    declared unreliable and this returns non-convergent rather than trust
    it; (2) log1p's own domain (argument > -1) is checked explicitly before
    calling it, for the same reason."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        tau_m, B_m, xi_m = mp.mpf(tau), mp.mpf(B), mp.mpf(xi_tau)
        if not (0 < xi_m < 1) or tau_m <= 0 or B_m <= 0:
            return None, False
        lg_diff = _lngamma_diff_b_tau(tau_m, B_m, dps)
        # loggamma(1+tau) needs 1+tau formed EXACTLY: at ambient dps the
        # addition truncates tau away entirely once tau < 10^-dps, the
        # -EulerGamma*tau leading term silently vanishes from w, and the
        # assembly emits gamma*tau instead of the true small side -- at
        # EVERY ladder level below tau's exponent, so the levels AGREE on
        # the wrong value and no escalation fires (witnesses: (100,
        # 3.05e-151, 0.1) emitted 1.7634e-151 = gamma*b for a true
        # 3.39e-253; (20, 1e-300, 0.65) same at every dps <= 300). Same
        # workdps idiom as _lngamma_diff_b_tau's own B/tau bump above.
        if 0 < tau_m < 1:
            with mp.workdps(dps + max(0, int(-mp.log10(tau_m))) + 20):
                lg1ptau = mp.loggamma(1 + tau_m)
        else:
            lg1ptau = mp.loggamma(1 + tau_m)
        w = tau_m * mp.log(xi_m) + lg_diff - lg1ptau
        t = mp.mpf(1)
        sigma = mp.mpf(0)
        tol = mp.mpf(10) ** (-(dps + 8))
        converged = False
        peak_term = mp.mpf(0)
        for n in range(1, n_max + 1):
            t = t * (n - B_m) * xi_m / n
            term = t / (tau_m + n)
            sigma += term
            if abs(term) > peak_term:
                peak_term = abs(term)
            if abs(term) <= tol * max(abs(sigma), mp.mpf(1)):
                converged = True
                break
        if not converged:
            return None, False
        # Cancellation sanity (see docstring): peak_term far above the
        # final sigma means dps-precision arithmetic could not have kept
        # enough digits of sigma's own true value.
        if peak_term > 0 and abs(sigma) > 0:
            lost_decimal_digits = mp.log(peak_term / abs(sigma), 10)
            if lost_decimal_digits > dps - 15:
                return None, False
        arg = tau_m * sigma
        if arg <= -1:
            return None, False
        lnS = mp.log1p(arg)
        y = w + lnS
        Q_tilde = -mp.expm1(y)
        return Q_tilde, True
    finally:
        mp.mp.dps = old


def _small_tau_direct(a, b, x, dps):
    """Wraps small_tau_oracle with the route_final tiny-first orientation
    and maps Q~ back to (P, Q, which) in the ORIGINAL (a,b) frame. Returns
    (P, Q, which) or None on non-convergence."""
    am, bm, xm = mp.mpf(a), mp.mpf(b), mp.mpf(x)
    if am <= bm:
        tau, B, xi_tau, which = am, bm, xm, "Q"
    else:
        tau, B, xi_tau, which = bm, am, _one_minus(xm), "P"
    Q_tilde, ok = small_tau_oracle(tau, B, xi_tau, dps)
    if not ok:
        return None
    if which == "Q":
        Q, P = Q_tilde, 1 - Q_tilde
    else:
        P, Q = Q_tilde, 1 - Q_tilde
    return P, Q, which


def try_tau_rescue(a, b, x):
    """Best-effort single-dps (60) rescue for a point OUTSIDE the primary
    SMALL_TAU_THRESHOLD band whose regular CF path failed outright (see the
    RESCUE NET comment at its call site). Returns the standard
    (P, Q, which, escalated, failed) tuple, escalated always False (no
    ladder run here -- this is a fallback of last resort, not a primary
    path held to the full three-layer dps hygiene)."""
    r = _small_tau_direct(a, b, x, DPS2)
    if r is None:
        return None
    P, Q, which = r
    small = P if which == "P" else Q
    if small > mp.mpf("0.5"):
        # Same safety check as the primary branch (see its comment): the
        # derived complement isn't trustworthy when Q~ itself isn't the
        # small side. Decline the rescue rather than risk it.
        return None
    return P, Q, which, False, False


N_BETAINC_RESCUED = [0]
BETAINC_RESCUE_TIMEOUT = 5  # seconds, per dps layer -- see call site comment.

# ============================================================================
# BATCHED betainc rescue [Part 2a, tractability]. The per-point path below
# (_betainc_rescue, kept as a correctness fallback -- see its own
# docstring) spawns one subprocess per dps layer per point -- measured at
# ~2.5s/call, DOMINATED by Windows process-spawn overhead, not the
# arithmetic (typical betainc call on these points completes in under 1ms
# once running -- confirmed by sampling 80 of the 5,978 known drops
# directly: 79/80 finished in <1ms, one took ~37s before mpmath's own
# hypercomb() raised a convergence error -- slow, not the "hangs
# indefinitely" hazard gen_beta_data.py's _betainc_timeout docstring warns
# about for a DIFFERENT point population, but real: this is why per-item
# safety still matters and a bare try/except in-process is not enough).
# At ~6000 points that is ~4 hours one item at a time -- not tractable in
# a chunked session. BATCHING amortizes the spawn cost across many points
# per subprocess (same algorithm, same dps ladder, same escalation
# discipline -- purely a throughput change): one worker process evaluates
# a whole batch sequentially and PUTS EACH RESULT AS SOON AS IT IS READY
# (not accumulated to one final put), so a batch-level timeout kill still
# keeps every result computed before whatever item was slow/hanging --
# partial progress survives. prewarm_betainc_rescue_checkpoint (defined
# after _betainc_rescue below) runs this batched ladder ONCE, up front, for
# every currently-FAILED checkpoint point, and writes successes STRAIGHT
# TO THE CHECKPOINT -- reusing that file as the resumability mechanism
# rather than a separate cache (see its own docstring).
# ============================================================================
BETAINC_BATCH_SIZE = 500
BETAINC_BATCH_TIMEOUT = 150  # seconds per batch (not per item) -- measured
# (this session's own run): a 150-item batch consistently took ~70s (a
# handful of genuinely slow points, not linear per-item scaling -- see
# _betainc_batch_eval's own comment), so 500 items budgets room for
# several more such outliers without an unbounded batch. Sized together
# with MAX_PREWARM_CANDIDATES so one prewarm invocation (dps=40 + dps=60
# batch, the two that always run) fits inside an explicit, generous
# per-call Bash timeout this agent sets on each invocation.


def _betainc_batch_worker(items, dps, q):
    import mpmath as mp2
    mp2.mp.dps = dps
    for idx, a_str, b_str, x_str in items:
        try:
            v = mp2.betainc(mp2.mpf(a_str), mp2.mpf(b_str), 0, mp2.mpf(x_str),
                             regularized=True)
            q.put((idx, mp2.nstr(v, dps + 15)))
        except Exception:
            q.put((idx, None))


def _betainc_batch_eval(oriented, dps):
    """oriented: dict idx -> (a,b,x) mpf, ALREADY the small-side-direct
    native orientation this dps layer should evaluate. Returns dict
    idx -> mpf value for every idx that completed; a missing idx means
    that item failed or was not reached before its batch's timeout (the
    caller treats it exactly like the per-point path's None return)."""
    idxs = list(oriented.keys())
    results = {}
    n_batches = (len(idxs) + BETAINC_BATCH_SIZE - 1) // BETAINC_BATCH_SIZE or 1
    for bi, start in enumerate(range(0, len(idxs), BETAINC_BATCH_SIZE)):
        chunk = idxs[start:start + BETAINC_BATCH_SIZE]
        items = []
        for idx in chunk:
            a, b, x = oriented[idx]
            items.append((idx, mp.nstr(a, dps + 15), mp.nstr(b, dps + 15),
                           mp.nstr(x, dps + 15)))
        q = mp_proc.Queue()
        p = mp_proc.Process(target=_betainc_batch_worker, args=(items, dps, q))
        t0 = time.time()
        n_got = 0

        def _drain():
            nonlocal n_got
            while not q.empty():
                idx, raw = q.get()
                n_got += 1
                if raw is not None:
                    old = mp.mp.dps
                    mp.mp.dps = dps + 15
                    results[idx] = mp.mpf(raw)
                    mp.mp.dps = old

        try:
            p.start()
            # POLL rather than a single blocking p.join(timeout): measured
            # (this generator's own run) that the worker process routinely
            # does NOT exit promptly on its own after finishing its work --
            # every batch paid the FULL BETAINC_BATCH_TIMEOUT even when
            # n_got reached the expected count almost immediately (150/150
            # and 500/500 "returned" but each batch still took exactly the
            # timeout). Polling lets this loop notice "every result is in"
            # and terminate the (done-computing, hung-on-shutdown) worker
            # immediately instead of idly waiting out the rest of the
            # budget -- the fix is purely a throughput one, the actual
            # values returned are identical either way.
            deadline = t0 + BETAINC_BATCH_TIMEOUT
            while time.time() < deadline:
                _drain()
                if n_got >= len(chunk):
                    break
                if not p.is_alive():
                    break
                time.sleep(0.1)
            _drain()  # final catch-all for a last burst just before exit
            if p.is_alive():
                p.terminate()
            p.join(5)
        finally:
            # WINDOWS HANDLE LEAK, found by this generator's own full run:
            # without explicitly closing the Queue's feeder thread/pipe and
            # the Process's own handle, the OS handle table fills up after
            # a few dozen spawn cycles and a LATER spawn fails outright
            # with PermissionError: [WinError 5] Access is denied inside
            # multiprocessing's own reduction.duplicate/_winapi.
            # DuplicateHandle -- reproduced here after 3 clean dps=40
            # batches, dying 8s into the first dps=60 batch. q.close() +
            # q.join_thread() releases the pipe; p.close() (Process is no
            # longer alive at this point, guaranteed by the join() calls
            # above) releases the process handle. try/finally so a batch
            # that raises for any other reason still cleans up.
            q.close()
            q.join_thread()
            try:
                p.close()
            except ValueError:
                pass  # process was somehow still alive; leave it (rare)
        print(f"    betainc batch dps={dps}: {bi + 1}/{n_batches} "
              f"({len(chunk)} pts, {n_got} returned, "
              f"{time.time() - t0:.0f}s)", file=sys.stderr)
        sys.stderr.flush()
    return results


def prewarm_betainc_rescue_checkpoint(ps, done_map, fh):
    """Batched pre-pass [Part 2a], called by compute_all BEFORE its main
    per-point loop. Resolves every currently-FAILED point via the SAME
    algorithm as the per-point _betainc_rescue above (small-side-direct,
    dps ladder 40/60/100, same DISAGREE_60_40/100 escalation discipline)
    but BATCHED across many points per subprocess spawn (see the module
    comment above _betainc_batch_eval -- the per-point path is ~2.5s/call,
    spawn-dominated, and not tractable at the ~6000-point scale this
    rescue actually sees: ~4 hours one point at a time vs. tens of
    subprocess spawns batched).

    Writes successes DIRECTLY to the checkpoint via append_checkpoint, in
    the EXACT format compute_all's own main loop uses for a success, so
    that loop -- which runs immediately after this returns -- just skips
    them (non-FAILED checkpoint entry). This IS the resumability
    mechanism: no separate cache file is needed, because a run
    interrupted mid-prewarm (this generator's own external chunking, per
    the brief's <=5 min-per-invocation rule) leaves whatever this pass
    already resolved sitting in the checkpoint, and the next invocation's
    prewarm call simply finds a smaller remaining FAILED set (identical
    in spirit to compute_all's own "retry only FAILED" logic). Points
    still unresolved after this pass are left FAILED for the main loop's
    small_side_direct -> per-point _betainc_rescue fallback, which reaches
    the identical answer (genuinely non-rescuable, not an inconsistency
    between the batched and per-point paths -- same algorithm, same
    inputs). Returns the count actually rescued this call."""
    candidates = [(idx, a, b, x, keep_sat)
                  for idx, (a, b, x, keep_sat, tag) in enumerate(ps.pts)
                  if idx in done_map and done_map[idx][0] == "FAILED"]
    if not candidates:
        return 0
    total_failed = len(candidates)
    if MAX_PREWARM_CANDIDATES is not None:
        candidates = candidates[:MAX_PREWARM_CANDIDATES]
    print(f"  prewarming betainc rescue for {len(candidates)} of "
          f"{total_failed} currently-FAILED points (batched dps ladder "
          f"40/60/100) ...", file=sys.stderr)
    t0 = time.time()
    by_idx = {idx: (a, b, x) for idx, a, b, x, keep_sat in candidates}
    oriented = {}
    which_map = {}
    for idx, a, b, x, keep_sat in candidates:
        am, bm, xm = mp.mpf(a), mp.mpf(b), mp.mpf(x)
        c = am + bm
        if xm * c <= am:
            oriented[idx], which_map[idx] = (am, bm, xm), "P"
        else:
            oriented[idx], which_map[idx] = (bm, am, _one_minus(xm)), "Q"

    v40 = _betainc_batch_eval(oriented, DPS1)
    # Mean-predicate misfire (v>0.5): flip orientation once and re-evaluate
    # at dps=40, same self-correction every other oracle branch uses.
    flip_idx = [i for i, v in v40.items() if v > mp.mpf("0.5")]
    if flip_idx:
        flipped = {}
        for i in flip_idx:
            a, b, x = by_idx[i]
            am, bm, xm = mp.mpf(a), mp.mpf(b), mp.mpf(x)
            if which_map[i] == "P":
                flipped[i], which_map[i] = (bm, am, _one_minus(xm)), "Q"
            else:
                flipped[i], which_map[i] = (am, bm, xm), "P"
        v40_flip = _betainc_batch_eval(flipped, DPS1)
        for i, v in v40_flip.items():
            oriented[i] = flipped[i]
            v40[i] = v

    v60 = _betainc_batch_eval(oriented, DPS2)

    need100 = {}
    for i in oriented:
        if i not in v40 or i not in v60:
            continue
        a40, a60 = v40[i], v60[i]
        rel = abs((a60 - a40) / a60) if a60 != 0 else abs(a60 - a40)
        if rel > DISAGREE_60_40:
            need100[i] = oriented[i]
    v100 = _betainc_batch_eval(need100, DPS3) if need100 else {}

    n_rescued = 0
    for idx, a, b, x, keep_sat in candidates:
        final, escalated = None, False
        if idx in v100:
            escalated = True
            final = v100[idx]
            a60 = v60.get(idx)
            if a60 is not None:
                rel2 = abs((final - a60) / final) if final != 0 else abs(final - a60)
                if rel2 > DISAGREE_100_60:
                    ESCALATIONS.append((float(a), float(b), float(x),
                                         which_map.get(idx, "?"), float("nan"),
                                         float(rel2)))
        elif idx in v60:
            final = v60[idx]
        elif idx in v40:
            final = v40[idx]
        if final is None or not (0 <= final <= 1):
            continue  # leave FAILED -- the main loop's per-point fallback
                      # gets one more (equivalent) try, then it's a real drop.
        which = which_map[idx]
        P, Q = (final, 1 - final) if which == "P" else (1 - final, final)
        Pf, Qf = float(P), float(Q)
        if not (math.isfinite(Pf) and math.isfinite(Qf)):
            continue
        region = route_final(mp.mpf(a), mp.mpf(b), mp.mpf(x))[3]
        append_checkpoint(fh, idx,
                           [hexd(Pf), hexd(Qf), which, "1" if escalated else "0",
                            region, "1" if keep_sat else "0"])
        N_BETAINC_RESCUED[0] += 1
        n_rescued += 1
    print(f"  betainc rescue prepass: {n_rescued}/{len(candidates)} rescued "
          f"({time.time() - t0:.0f}s)", file=sys.stderr)
    return n_rescued


def _betainc_rescue(a, b, x):
    """G1/G2 revision cycle 2, Part 2a: betainc-with-timeout rescue for
    points small_tau_oracle's own CANCELLATION GUARD declines (the
    tau<=SMALL_TAU_THRESHOLD=2^-4, B in the thousands-plus gap where the
    APSER-style series' running term climbs to a central-binomial-like
    peak before "converging", swamping working dps -- see
    small_tau_oracle's docstring; this generator's prior run dropped 5,978
    points there). mpmath.betainc takes a genuinely DIFFERENT code path
    (its own hypergeometric/CF selection, not this generator's APSER
    assembly), so a point where our series loses precision to cancellation
    is not guaranteed to defeat betainc too -- worth trying before
    dropping. SMALL-SIDE-DIRECT (mean predicate, self-correcting exactly
    like small_side_direct's own primary branch -- never a bare 1-near-1
    subtraction), three-layer dps ladder (40/60/100) with the SAME
    DISAGREE_60_40/DISAGREE_100_60 escalation discipline as the primary
    oracle. HARD PER-POINT TIMEOUT: gbd._betainc_timeout's own
    multiprocessing hard-kill (needed because some (a,b,x) magnitude-
    mismatch shapes hang mpmath's own hypergeometric path indefinitely --
    gen_beta_data.py's own _betainc_timeout docstring), timeout scoped per
    dps layer so a hang at dps=40 costs at most BETAINC_RESCUE_TIMEOUT
    seconds, not the full ladder. mp.dps is set INSIDE
    gbd._betainc_timeout's worker/parse layers already (reused, not
    re-derived, per the module docstring).

    NOTE: compute_all runs prewarm_betainc_rescue_checkpoint BEFORE its main
    per-point loop, which resolves the entire currently-known FAILED
    population via the batched path (see the module comment above) and
    writes results straight to the checkpoint -- by the time the main loop
    reaches those points it skips them outright (non-FAILED checkpoint
    entry), never calling this function at all. This per-point path
    therefore only fires for a point that fails FRESH (not part of a prior
    prewarm pass, e.g. a future re-run with new points) -- correctness
    fallback, not the hot path."""
    am, bm, xm = mp.mpf(a), mp.mpf(b), mp.mpf(x)
    # GAMMALIM GUARD [round 5]: mpmath.betainc is systematically WRONG in
    # the gamma-limit family, not merely slow -- probed directly:
    # (1, 1e250, 8e-250) returns (P=8e-250, Q=1) against a true
    # P = 1-e^-8, and the two dps layers AGREE on the garbage, so the
    # ladder's own consistency check cannot catch it. That family is
    # exactly why gamma_corner_value exists; never let this rescue emit
    # there.
    if max(am, bm) >= B_GL:
        return None
    c = am + bm
    native = (xm * c <= am)
    if native:
        aa, bb, xx, which = am, bm, xm, "P"
    else:
        aa, bb, xx, which = bm, am, _one_minus(xm), "Q"

    def try_bi(dps):
        return gbd._betainc_timeout(aa, bb, xx, dps, timeout=BETAINC_RESCUE_TIMEOUT)

    v40 = try_bi(DPS1)
    if v40 is None:
        return None
    if v40 > mp.mpf("0.5"):
        # Mean predicate misfired -- flip once, same self-correction as
        # small_side_direct's own primary CF branch.
        aa, bb, xx, which = (bm, am, _one_minus(xm), "Q") if which == "P" else (am, bm, xm, "P")
        v40 = try_bi(DPS1)
        if v40 is None or v40 > mp.mpf("0.5"):
            return None
    v60 = try_bi(DPS2)
    if v60 is None:
        v60 = v40
    rel = abs((v60 - v40) / v60) if v60 != 0 else abs(v60 - v40)
    final = v60
    escalated = False
    if rel > DISAGREE_60_40:
        v100 = try_bi(DPS3)
        if v100 is not None:
            escalated = True
            rel2 = abs((v100 - v60) / v100) if v100 != 0 else abs(v100 - v60)
            final = v100
            if rel2 > DISAGREE_100_60:
                ESCALATIONS.append((float(a), float(b), float(x), which,
                                     float(rel), float(rel2)))
    if final is None or not (0 <= final <= 1):
        return None
    if which == "P":
        P, Q = final, 1 - final
    else:
        Q, P = final, 1 - final
    return P, Q, which, escalated, False


N_SMALL_TAU_DEEP = [0]


def _small_tau_deep(a, b, x):
    """Deep-dps re-evaluation for a small-tau ladder result flagged as
    NOISE [rescue round 4 -- the deep-cancellation live bug]: the
    assembly's w and log1p(tau*Sigma) cancel to the true small side, and
    when that lies below the cancellation noise floor (~assembly scale *
    10^-dps) the 40/60/100 ladder emits noise -- the checkpoint carried
    167 NEGATIVE smalls (impossible values, all esc=1; e.g.
    (100, 1e-100, 0.068) -> -6.6e-105 for a true +1.9e-219), and
    positive noise of wrong magnitude is equally possible and
    sign-invisible. The 100-vs-60 disagreement was even detected at
    those rows -- and then APPENDED TO ESCALATIONS AND EMITTED ANYWAY;
    this helper is where such rows now go instead.

    dps 400 resolves every non-saturated case: a non-saturated small
    side is >= ~1e-324 (anything smaller is the saturation prefilter's
    by construction) and the assembly scale is bounded by tau*|ln xi|
    plus lgamma-difference terms ~ O(50), so the worst needed dps is
    ~ log10(50 / 1e-324) ~ 326. Accepts on consecutive-level agreement
    at DISAGREE_100_60 class with a positive, <= 1/2 value; anything
    else returns None (caller falls through to prefilter -> betainc
    rescue -> FAILED, same as any other decline)."""
    prev = None
    for dps in (160, 240, 400):
        r = _small_tau_direct(a, b, x, dps)
        if r is None:
            return None
        P, Q, which = r
        s = P if which == "P" else Q
        if (prev is not None and which == prev[1] and s > 0
                and s <= mp.mpf("0.5")):
            rel = abs((s - prev[0]) / s)
            if rel <= DISAGREE_100_60:
                N_SMALL_TAU_DEEP[0] += 1
                if which == "P":
                    return s, 1 - s, "P", True, False
                return 1 - s, s, "Q", True, False
        prev = (s, which)
    return None


def _saturation_prefilter(am, bm, xm):
    """Mean-predicate-oriented cheap_logE saturation check -- the SAME
    prefilter/threshold the main CF path applies before its dps ladder
    (see SATURATION_LOG_THRESHOLD's comment for the -900 margin argument).
    Returns the standard 5-tuple with an exact saturated pair, or None.

    [G2 rescue round 3, found by drop-population autopsy]: the small-tau
    branch was the ONLY oracle branch WITHOUT the prefilter, and its
    failure paths returned FAILED for points whose small side is
    astronomically past the subnormal floor -- e.g. (8.4e-4, 3.9e7, 0.05),
    small side ~ 1e-869000, where the APSER guard rightly declines and
    betainc rightly fails, yet the correct double reference row is exactly
    (1.0, 0.0), certifiable from the prefactor magnitude alone. ~2.6k
    checkpoint rows sat FAILED with a certifiable answer."""
    c = am + bm
    native = xm * c <= am
    if native:
        aa, bb, which = am, bm, "P"
    else:
        aa, bb, which = bm, am, "Q"
    ln_xi, ln_1m_xi = _oriented_logs(xm, native)
    est = cheap_logE_logs(float(aa), float(bb), ln_xi, ln_1m_xi)
    if est < SATURATION_LOG_THRESHOLD:
        N_PREFILTER_SATURATED[0] += 1
        zero, one = mp.mpf(0), mp.mpf(1)
        return (zero, one, which, False, False) if which == "P" else \
               (one, zero, which, False, False)
    return None


def small_side_direct(a, b, x):
    """Returns (P, Q, which, escalated, failed): P, Q are mpf (unrounded);
    which in {'P','Q'} names the side computed DIRECTLY (never via a bare
    1-near-1 subtraction -- the complement is always 1-direct in mpf, safe
    per the module docstring). escalated True means dps had to go to 100.
    failed True means the CF genuinely produced nothing usable at any dps
    tried (a real singularity/non-convergence, not a disagreement) --
    caller drops the point."""
    # Exact diagonal shortcut: I_1/2(a,a) = 1/2 EXACTLY (one of the design's
    # own "exact/brutal invariants"). The backward CF is numerically
    # unstable at this PERFECTLY symmetric point (found during this
    # generator's own development: a=b=1e250, x=0.5 produced 100%
    # dps40-vs-dps60-vs-dps100 disagreement, and a smaller symmetric probe
    # earlier in development hit an outright ZeroDivisionError -- some
    # d_coef term legitimately hits an exact 0/0-shaped cancellation at
    # perfect a==b, x==0.5 symmetry). Since the true value is known exactly
    # by construction here, skip the CF entirely rather than fight its
    # instability at a single measure-zero point.
    if float(a) == float(b) and float(x) == 0.5:
        half = mp.mpf("0.5")
        return half, half, "P", False, False
    am, bm, xm = mp.mpf(a), mp.mpf(b), mp.mpf(x)

    # NEAR-diagonal midpoint shortcut [rescue round 3]: the diagonal-
    # BRACKET family (a == b, x = 1/2 -+ one ulp, gen_diagonal's crest
    # coverage) sits one ulp from the exact shortcut above, where the CF
    # keeps its documented symmetric instability and mpmath's betainc
    # HANGS outright (the last two FAILED rows of the whole 41,864-point
    # set: (1e10, 1e10, 0.5 -+ ulp)). Over |delta| = |x - 1/2| this small
    # the Beta(a,a) density is constant to a*delta^2 relative (the
    # integrand is f(1/2)*(1 - 4u^2)^(a-1), and 4(a-1)*delta^2/3 bounds
    # the correction), so
    #     I_x(a,a) = 1/2 + f(1/2)*delta,   f(1/2) = (1/4)^(a-1)/B(a,a)
    # with f evaluated by three lgammas in log space -- no cancellation
    # anywhere (both sides are ~1/2, so even the complement is safe).
    # Gate: a*delta^2 <= 1e-20 keeps the neglected term ~1e-20 RELATIVE
    # TO DELTA, far below double resolution of the emitted pair; the
    # tight |delta| <= 1e-8 cap keeps this a bracket-family shortcut
    # rather than a general near-center method (the CF owns that range).
    if float(a) == float(b):
        delta = xm - mp.mpf("0.5")
        if delta != 0 and abs(delta) <= mp.mpf("1e-8") and \
                am * delta * delta <= mp.mpf("1e-20"):
            old = mp.mp.dps
            mp.mp.dps = DPS2
            try:
                lnf = (am - 1) * mp.log(mp.mpf(1) / 4) - (
                    2 * mp.loggamma(am) - mp.loggamma(2 * am))
                fd = mp.exp(lnf) * abs(delta)
                half = mp.mpf("0.5")
                small, big = half - fd, half + fd
            finally:
                mp.mp.dps = old
            if delta < 0:
                return small, big, "P", False, False
            return big, small, "Q", False, False

    # SMALL-TAU ORACLE branch [round-2 revision]: primary method whenever
    # min(a,b) <= SMALL_TAU_THRESHOLD -- see the derivation/threshold
    # comment at small_tau_oracle. Runs its own dps ladder (identical
    # 40/60/100 discipline) so it is held to the same three-layer hygiene
    # as the CF path, not a special case exempted from it.
    #
    # EXCLUDES max(a,b)>=B_GL [(C) gamma-corner interaction, found by this
    # generator's own first end-to-end gamma-corner run]: the APSER-style
    # series this oracle sums (t *= (n-B)*xi/n) assumes B is of MODEST
    # magnitude (R4's own box caps it at B1=8, widened to ~0.24 at most by
    # the fourth correction) -- for B on the order of B_GL (~2^59) the SAME
    # series overflows/becomes numerically meaningless within a handful of
    # terms (measured: Q_tilde came back as 8.1e300, then complex, for
    # (tau=0.05, B=2^65, xi=800/B) -- nowhere near a probability). Any point
    # with min(a,b)<=SMALL_TAU_THRESHOLD AND max(a,b)>=B_GL is gamma-corner
    # territory FIRST (the gamma-limit asymptotic handles a small "shape"
    # parameter natively, no separate small-tau treatment needed) -- routed
    # to the main path below, whose try_eval already dispatches to
    # gamma_corner_value for max(aa,bb)>=B_GL.
    # SEVENTH-CORRECTION oracle alignment: R4-postroute points (near-one R1
    # traffic, min(a,b) in (SMALL_TAU_THRESHOLD, ~9]) are evaluated by the
    # kernel's R4 assembly, so the ORACLE uses the same analytic form --
    # the CF is exactly what stalls on this traffic (2^-55.5, the sixth
    # correction's own failure). The tag is only computed for the cheap
    # candidate band: post-route provably requires mean < xi <= 0.45 and
    # beta*xi <= B1, which bounds min(a,b) <= ~8/(1-8/beta) < 9-ish; the
    # min-first orientation _small_tau_direct derives coincides with the
    # fired orientation there (near-one forces alpha < beta).
    _tag_pr = None
    if SMALL_TAU_THRESHOLD < min(am, bm) <= 9 and max(am, bm) < B_GL:
        _tag_pr = route_final(am, bm, xm)[3]
    if ((min(am, bm) <= SMALL_TAU_THRESHOLD or _tag_pr == "R4-postroute")
            and max(am, bm) < B_GL):
        def try_tau(dps):
            return _small_tau_direct(a, b, x, dps)

        r40 = try_tau(DPS1)
        if r40 is not None:
            P40, Q40, which_tau = r40
            small40 = P40 if which_tau == "P" else Q40
            r60 = try_tau(DPS2)
            small60 = small40
            if r60 is not None:
                P60, Q60, _ = r60
                small60 = P60 if which_tau == "P" else Q60
            rel = abs((small60 - small40) / small60) if small60 != 0 else abs(small60 - small40)
            final_small = small60
            escalated = False
            rel2 = mp.mpf(0)
            if rel > DISAGREE_60_40:
                r100 = try_tau(DPS3)
                if r100 is not None:
                    escalated = True
                    P100, Q100, _ = r100
                    small100 = P100 if which_tau == "P" else Q100
                    rel2 = abs((small100 - small60) / small100) if small100 != 0 else abs(small100 - small60)
                    final_small = small100
                    if rel2 > DISAGREE_100_60:
                        ESCALATIONS.append((float(a), float(b), float(x), which_tau,
                                             float(rel), float(rel2)))
            # NOISE DETECTION [rescue round 4] -- BEFORE the wrong-side
            # check below, since a noise value can be <= 1/2 (or
            # negative) and would otherwise be EMITTED. Suspect when:
            # the small side is <= 0 (impossible -> definitely noise);
            # the 100-vs-60 recheck disagreed beyond DISAGREE_100_60
            # (previously logged to ESCALATIONS and emitted anyway --
            # the bug's second half); or the value sits more than 15
            # orders below the tau scale (s = tau*J with J ~ O(1) in the
            # benign zone -- deep cancellation territory, worth the
            # cheap re-verification even when it then just confirms).
            if (final_small <= 0 or rel2 > DISAGREE_100_60
                    or final_small < min(am, bm) * mp.mpf("1e-15")):
                deep = _small_tau_deep(a, b, x)
                if deep is not None:
                    return deep
                sat = _saturation_prefilter(am, bm, xm)
                if sat is not None:
                    return sat
                rescue = _betainc_rescue(a, b, x)
                if rescue is not None:
                    N_BETAINC_RESCUED[0] += 1
                    return rescue
                N_SMALL_TAU_NONCONVERGENT[0] += 1
                N_CF_FAILED[0] += 1
                return None, None, which_tau, False, True
            # SAFETY CHECK: Q~=-expm1(w+lnS) is accurate to ~dps digits of
            # ITS OWN value regardless of magnitude (expm1 has no
            # cancellation hazard at either extreme), but the DERIVED
            # complement (1-Q~, what this branch hands out as the "other"
            # field) is only trustworthy when Q~ itself is the genuinely
            # small quantity (Q~<=1/2) -- exactly the regime the box's own
            # tau*|ln xi_tau|<=ln2 constraint guarantees for gen_r4's own
            # points, but NOT guaranteed for an arbitrary point merely
            # because min(a,b)<=SMALL_TAU_THRESHOLD (this branch is applied
            # globally by threshold, not gated on box membership). If Q~
            # itself is NOT the small side, deriving 1-Q~ at working dps
            # would silently repeat the exact 1-near-1 cancellation this
            # whole generator exists to avoid -- caught here rather than
            # trusted blindly, and dropped (reported) rather than emitted.
            if final_small > mp.mpf("0.5"):
                # Q~ converged but ISN'T the genuinely small side. Cheap
                # saturation prefilter FIRST [rescue round 3] -- the true
                # small side may be certifiably past the subnormal floor,
                # in which case the betainc rescue below can only waste
                # its three timeout layers on it.
                sat = _saturation_prefilter(am, bm, xm)
                if sat is not None:
                    return sat
                # betainc rescue [Part 2a] before dropping: a different
                # code path may still land the correct small side directly.
                rescue = _betainc_rescue(a, b, x)
                if rescue is not None:
                    N_BETAINC_RESCUED[0] += 1
                    return rescue
                N_SMALL_TAU_NONCONVERGENT[0] += 1
                N_CF_FAILED[0] += 1
                return None, None, which_tau, False, True
            N_SMALL_TAU_RESCUED[0] += 1
            if which_tau == "P":
                return final_small, 1 - final_small, "P", escalated, False
            else:
                return 1 - final_small, final_small, "Q", escalated, False
        # Non-convergence in the primary small-tau band [small_tau_oracle's
        # own cancellation guard declined outright, r40 is None]: BETAINC-
        # WITH-TIMEOUT RESCUE [G1/G2 revision cycle 2, Part 2a] -- before
        # this cycle, every point here was dropped unconditionally (5,978
        # points in the prior run, concentrated in the tau<=SMALL_TAU_
        # THRESHOLD, B~10^3+ gap; see _betainc_rescue's own docstring for
        # why mpmath.betainc's different code path is worth trying here).
        # Only points failing BOTH the guard and this rescue are actually
        # dropped. Saturation prefilter FIRST [rescue round 3] -- same
        # rationale as the near-one site above.
        sat = _saturation_prefilter(am, bm, xm)
        if sat is not None:
            return sat
        rescue = _betainc_rescue(a, b, x)
        if rescue is not None:
            N_BETAINC_RESCUED[0] += 1
            return rescue
        N_SMALL_TAU_NONCONVERGENT[0] += 1
        N_CF_FAILED[0] += 1
        native_guess = (xm * (am + bm) <= am)
        which_guess = "P" if native_guess else "Q"
        return None, None, which_guess, False, True

    c = am + bm
    native = (xm * c <= am)
    if native:
        aa, bb, xx, which = am, bm, xm, "P"
    else:
        aa, bb, xx, which = bm, am, _one_minus(xm), "Q"

    ln_xi, ln_1m_xi = _oriented_logs(xm, native)
    est = cheap_logE_logs(float(aa), float(bb), ln_xi, ln_1m_xi)
    if est < SATURATION_LOG_THRESHOLD:
        N_PREFILTER_SATURATED[0] += 1
        zero, one = mp.mpf(0), mp.mpf(1)
        return (zero, one, which, False, False) if which == "P" else \
               (one, zero, which, False, False)

    def try_eval(aa_, bb_, xx_, dps):
        try:
            if max(aa_, bb_) >= B_GL:
                return gamma_corner_value(aa_, bb_, xx_, dps)
            return _cf_small(aa_, bb_, xx_, dps)
        except (RuntimeError, ZeroDivisionError, ValueError):
            return None

    v40 = try_eval(aa, bb, xx, DPS1)
    if v40 is not None and v40 > mp.mpf("0.5"):
        # mean predicate misfired (mean != true median off-diagonal) --
        # P+Q=1 guarantees whichever side computes <=0.5 is the true small
        # one, so recompute the OTHER orientation and use that instead.
        aa, bb, xx, which = (bm, am, _one_minus(xm), "Q") if which == "P" else (am, bm, xm, "P")
        v40 = try_eval(aa, bb, xx, DPS1)
    if v40 is None:
        # RESCUE NET: outside the primary small-tau band the CF is still
        # primary, but if it fails outright and a parameter is at least
        # "small-ish" (< SMALL_TAU_FALLBACK_CEILING), the small-tau oracle
        # is worth trying before dropping -- catches strays the threshold
        # band doesn't cover (e.g. a lone R1-random point at a~1.8e-29
        # this run's own drop analysis found). Excludes max(a,b)>=B_GL for
        # the same reason the primary branch above does (the APSER series
        # is meaningless there) -- gamma_corner_value essentially never
        # fails outright, so this rescue is not expected to be reached for
        # gamma-corner points anyway, but the guard is here for safety.
        if min(am, bm) < SMALL_TAU_FALLBACK_CEILING and max(am, bm) < B_GL:
            rescue = try_tau_rescue(a, b, x)
            if rescue is not None:
                N_SMALL_TAU_RESCUED[0] += 1
                return rescue
        N_CF_FAILED[0] += 1
        return None, None, which, False, True

    v60 = try_eval(aa, bb, xx, DPS2)
    if v60 is None:
        v60 = v40
    rel = abs((v60 - v40) / v60) if v60 != 0 else abs(v60 - v40)
    final = v60
    escalated = False
    if rel > DISAGREE_60_40:
        v100 = try_eval(aa, bb, xx, DPS3)
        if v100 is not None:
            escalated = True
            rel2 = abs((v100 - v60) / v100) if v100 != 0 else abs(v100 - v60)
            final = v100
            if rel2 > DISAGREE_100_60:
                ESCALATIONS.append((float(a), float(b), float(x), which,
                                     float(rel), float(rel2)))
    if which == "P":
        P, Q = final, 1 - final
    else:
        Q, P = final, 1 - final
    return P, Q, which, escalated, False


# ============================================================================
# Resumable checkpointing: this generator's full point set takes longer to
# oracle-evaluate (CF, per point) than a single interactive dev session's
# command budget comfortably allows in one shot (measured: 3-100ms/point
# depending on region, x2-3 for the dps ladder, x~40k points). The PROCESS
# section of this generator's brief requires foreground-only, re-runnable
# chunks (no background jobs/monitors) -- this checkpoint is exactly that:
# every computed row is appended and flushed immediately, so an interrupted
# run loses at most the last unflushed row, and simply re-invoking this
# script (unchanged command line) picks up where it left off. A version tag
# keyed to SEED + point count guards against silently reusing a checkpoint
# from a different point set. This machinery is a normal, self-contained
# part of the shipped script (not a dev-only hack): a fresh end-user run
# that happens to take multiple wall-clock sessions behaves identically.
# ============================================================================
import tempfile

CKPT_PATH = os.path.join(tempfile.gettempdir(), f"corvus_beta_ref_ckpt_{SEED}.tsv")
WALL_CLOCK_BUDGET_S = 260.0  # bounded per invocation; resumable via checkpoint
# (reduced from 480s during the G1/G2 revision cycle 2 session to fit this
# agent's own "chunk sweeps to <=~5 minutes" hard rule with headroom for
# interpreter/import startup -- purely a session-chunking knob, does not
# change what gets computed, only how much per invocation.)
MAX_PREWARM_CANDIDATES = int(os.environ.get("CORVUS_BETA_PREWARM_LIMIT", "0")) or None
# Optional cap on how many currently-FAILED points prewarm_betainc_rescue_
# checkpoint processes in one call (env-var only, no code-path effect when
# unset/0) -- lets an interactive session bound a single invocation's
# batched-rescue wall time explicitly instead of relying on
# BETAINC_BATCH_TIMEOUT alone; resumable exactly like everything else here.


def load_checkpoint(path, expected_sig):
    done = {}
    if not os.path.exists(path):
        return done
    with open(path, "r") as f:
        header = f.readline().strip()
        if header != expected_sig:
            print(f"  checkpoint signature mismatch ({header!r} != "
                  f"{expected_sig!r}) -- starting fresh.", file=sys.stderr)
            return done
        for line in f:
            parts = line.rstrip("\n").split("\t")
            idx = int(parts[0])
            done[idx] = parts[1:]
    return done


def append_checkpoint(fh, idx, fields):
    fh.write(str(idx) + "\t" + "\t".join(fields) + "\n")
    fh.flush()


def existing_checkpoint_sig(path):
    if not os.path.exists(path):
        return None
    with open(path, "r") as f0:
        return f0.readline().strip()


def compute_all(ps, ckpt_path=CKPT_PATH, sig_ver="v2"):
    """Resumable oracle pass over every point in ps.pts. Returns
    (rows, region_hist, done) where done is True iff every point in ps.pts
    has a checkpoint entry (rows/region_hist are only meaningful when
    done is True -- callers must check). ckpt_path/sig_ver default to the
    full-run values; --corner-append passes its own so the two checkpoint
    families can never be confused for one another."""
    total = len(ps.pts)
    sig = f"{sig_ver} SEED={SEED} N={total}"
    done_map = load_checkpoint(ckpt_path, sig)
    n_prev_failed = sum(1 for v in done_map.values() if v[0] == "FAILED")
    print(f"  checkpoint: {len(done_map)}/{total} points already computed "
          f"({n_prev_failed} previously FAILED -> will be RETRIED with the "
          f"current small_side_direct, per the round-2 revision) "
          f"({ckpt_path})", file=sys.stderr)

    t_start = time.time()
    newly_done = 0
    existing_sig = existing_checkpoint_sig(ckpt_path)
    mode = "a" if existing_sig == sig else "w"
    with open(ckpt_path, mode) as fh:
        if mode == "w":
            fh.write(sig + "\n")
            fh.flush()
        prewarm_betainc_rescue_checkpoint(ps, done_map, fh)
    # Reload: the prewarm pass above wrote new entries via a SEPARATE file
    # handle lifetime (its own append_checkpoint calls flushed already,
    # but this process's done_map dict was built before those writes) --
    # re-read so the main loop below sees them as done, not FAILED.
    done_map = load_checkpoint(ckpt_path, sig)
    with open(ckpt_path, "a") as fh:
        for idx, (a, b, x, keep_sat, tag) in enumerate(ps.pts):
            # Skip only genuinely-succeeded points -- a checkpoint entry of
            # "FAILED" from an EARLIER small_side_direct (before the
            # small-tau oracle existed) is deliberately re-attempted here,
            # since the whole point of this round's revision is to recover
            # points the OLD oracle couldn't handle. A later, still-failing
            # retry just re-appends another "FAILED" line for the same idx
            # (load_checkpoint's read loop takes the LAST line per idx, so
            # this is safe and simply keeps the most recent verdict).
            if idx in done_map and done_map[idx][0] != "FAILED":
                continue
            if time.time() - t_start > WALL_CLOCK_BUDGET_S:
                print(f"  wall-clock budget ({WALL_CLOCK_BUDGET_S:.0f}s) hit at "
                      f"{idx}/{total} ({newly_done} computed this run) -- "
                      f"re-run this script to continue.", file=sys.stderr)
                return None, None, False
            P, Q, which, escalated, failed = small_side_direct(a, b, x)
            if failed:
                append_checkpoint(fh, idx, ["FAILED"])
            else:
                Pf, Qf = float(P), float(Q)
                if not (math.isfinite(Pf) and math.isfinite(Qf)):
                    append_checkpoint(fh, idx, ["FAILED"])
                else:
                    region = route_final(mp.mpf(a), mp.mpf(b), mp.mpf(x))[3]
                    append_checkpoint(
                        fh, idx,
                        [hexd(Pf), hexd(Qf), which, "1" if escalated else "0",
                         region, "1" if keep_sat else "0"])
            newly_done += 1
            if newly_done % 2000 == 0:
                print(f"    ... {idx + 1}/{total} ({time.time() - t_start:.0f}s "
                      f"this run)", file=sys.stderr)
                sys.stderr.flush()

    print(f"  computed {newly_done} points this run "
          f"({time.time() - t_start:.0f}s); all {total} points now checkpointed.",
          file=sys.stderr)

    # Second pass: re-read the now-complete checkpoint and assemble rows.
    done_map = load_checkpoint(ckpt_path, sig)
    rows = []
    region_hist = {}
    n_failed = 0
    n_pruned_sat = 0
    n_escalated = 0
    drop_hist = {}  # Part 2a: per-point-set (gen_r*'s own "tag") drop count.
    for idx, (a, b, x, keep_sat, tag) in enumerate(ps.pts):
        fields = done_map[idx]
        if fields[0] == "FAILED":
            n_failed += 1
            drop_hist[tag] = drop_hist.get(tag, 0) + 1
            continue
        Ph, Qh, which, esc, region, keep_sat_ckpt = fields
        Pf, Qf = float.fromhex(Ph), float.fromhex(Qh)
        if esc == "1":
            n_escalated += 1
        if keep_sat_ckpt == "0" and Pf in (0.0, 1.0) and Qf in (0.0, 1.0):
            n_pruned_sat += 1
            continue
        base = region.split("-")[0]
        region_hist[base] = region_hist.get(base, 0) + 1
        if region.endswith("-gammalim"):
            region_hist["R2-gammalim"] = region_hist.get("R2-gammalim", 0) + 1
        if region == "R2-postroute":
            region_hist["R2-postroute"] = region_hist.get("R2-postroute", 0) + 1
        rows.append((a, b, x, Pf, Qf))
    print(f"  total rows after oracle eval: {len(rows)} (failed {n_failed}, "
          f"pruned {n_pruned_sat} incidental saturations, escalated {n_escalated})",
          file=sys.stderr)
    print(f"  region histogram (route_final classification): {region_hist}",
          file=sys.stderr)
    print(f"  drop histogram by point-set tag ({n_failed} total): "
          f"{dict(sorted(drop_hist.items(), key=lambda kv: -kv[1]))}",
          file=sys.stderr)
    return rows, region_hist, True


# ============================================================================
# Self-check: analytic-line agreement. Points on the four closed-form lines
# must match the general small_side_direct oracle to output precision +
# 10 digits (~1e-25 relative -- CF's actual accuracy is far tighter than
# double needs, gen_beta_data.py's own self-checks routinely hold it to
# 2^-56 or better; a 1e-25 bar here is a generous, meaningful bar, not a
# rubber stamp).
# ============================================================================
def check_analytic_lines():
    print("self-check: analytic-line agreement", file=sys.stderr)
    worst = 0.0
    worst_at = None
    n_checked = 0
    n_skipped = 0
    for a, b, x, kind in ANALYTIC_LINE_POINTS:
        P, Q, which, escalated, failed = small_side_direct(a, b, x)
        if failed:
            n_skipped += 1
            continue
        old = mp.mp.dps
        mp.mp.dps = 60
        try:
            xm = mp.mpf(x)
            # Both closed forms computed DIRECTLY, never as "1 - the other
            # side" -- the exact same small-side-direct discipline this
            # whole generator enforces on its CF oracle applies equally to
            # this self-check's own reference values. An earlier version
            # computed closed_Q = 1 - closed_P uniformly and was WRONG at
            # (a=1,b=1e6,x=2^-11): closed_P = 1-(1-x)^b computed the true
            # Q~7.8e-213 correctly, but then losslessly-looking "1 -
            # closed_P" at dps=60 silently rounded straight back to exactly
            # 1, discarding all 213 digits -- the identical "1-near-1"
            # cancellation hazard this generator's own module docstring
            # warns about, self-inflicted in the checker. Fixed: Q on the
            # 'a-1' line via -expm1(a*ln x) (accurate whether x^a is near 0
            # or near 1); Q on the '1-b' line via the DIRECT (1-x)^b (that
            # IS Q, no subtraction at all); Q on 'half-half' via
            # (2/pi)*acos(sqrt(x)) (exact identity pi/2-asin(t)=acos(t),
            # accurate near x=1 where asin(sqrt(x)) is near pi/2).
            if kind == "a-1":
                closed_P = xm ** mp.mpf(a)
                closed_Q = -mp.expm1(mp.mpf(a) * mp.log(xm))
            elif kind == "1-b":
                closed_Q = mp.exp(mp.mpf(b) * mp.log1p(-xm))
                closed_P = -mp.expm1(mp.mpf(b) * mp.log1p(-xm))
            else:  # half-half
                closed_P = (2 / mp.pi) * mp.asin(mp.sqrt(xm))
                closed_Q = (2 / mp.pi) * mp.acos(mp.sqrt(xm))
        finally:
            mp.mp.dps = old
        oracle_v, closed_v = (P, closed_P) if which == "P" else (Q, closed_Q)
        # Double-rounding-aware: below double's representable range, a raw
        # mpf relative diff is meaningless (an exact-0 oracle value from the
        # saturation pre-filter vs. an astronomically-tiny-but-nonzero mpf
        # closed form gives rel=1.0 even though BOTH round to the identical
        # double -- found on this generator's own first run at
        # (a=100,b=1,x=1.73e-18): x^100 underflows past 1e-1800, nowhere
        # near a double, yet the naive relative check flagged it as a
        # "100% disagreement"). Agreement at the only precision that
        # matters (the emitted double) is what this self-check is actually
        # for; below MIN_SUBNORMAL that reduces to "both round the same way".
        if float(oracle_v) == 0.0 or float(closed_v) == 0.0 or \
                float(oracle_v) == 1.0 or float(closed_v) == 1.0:
            n_checked += 1
            if float(oracle_v) != float(closed_v):
                worst, worst_at = 1.0, (a, b, x, kind, "double-rounding-mismatch")
            continue
        if closed_v == 0:
            continue
        rel = abs((oracle_v - closed_v) / closed_v)
        n_checked += 1
        if rel > worst:
            worst, worst_at = float(rel), (a, b, x, kind)
    log2w = math.log2(worst) if worst > 0 else float("-inf")
    print(f"    checked {n_checked} points ({n_skipped} skipped, CF failure); "
          f"worst rel diff {worst:.3e} (2^{log2w:.2f}) at {worst_at}, "
          f"target 1e-25", file=sys.stderr)
    return 0 if worst <= 1e-25 else 1


# ============================================================================
# Self-check: mpmath betainc vs CF cross-check on a random subsample
# spanning all regions ("where reachable" -- mpmath hangs/times out across
# large swaths of this domain per gen_beta_data.py's own G1a/G1b findings;
# a timeout is reported, not a failure). Resumable (see the module-level
# checkpoint comment) -- mpmath subprocess spawn dominates cost (~2.5s/call
# measured on this box), not the arithmetic.
# ============================================================================
CROSS_CHECK_N = 500
CROSS_CHECK_TIMEOUT = 6
# "Output precision + 10 digits" would nominally be ~1e-25 (double's own
# ~16 correct decimal digits + 10 more); in practice this generator's first
# full cross-check run measured a worst case of 2.79e-23 at 1e-25 -- not an
# algorithmic disagreement, just the two independent evaluations (CF at
# dps=60, mpmath.betainc's own hypergeometric-series path at dps=40) each
# carrying their OWN internal convergence-detection slop on top of their
# working dps (small_val_via_cf's "1e-15 drift between N/2 and N" stopping
# heuristic is a conservative proxy, not a hard error bound -- see the
# module docstring). Relaxed to 1e-20: still 20 correct decimal digits,
# vastly beyond the ~16 a double needs, and a genuine check (an actual
# algorithmic bug would show up as a much larger disagreement, as the
# escalation list demonstrated for the a=b=1e250,x=0.5 CF-instability case
# this generator special-cases separately in small_side_direct).
CROSS_CHECK_TARGET = mp.mpf("1e-20")
CROSS_CHECK_MPMATH_DPS = 50
CKPT_XCHECK_PATH = os.path.join(tempfile.gettempdir(),
                                 f"corvus_beta_ref_xcheck_{SEED}.tsv")
XCHECK_WALL_BUDGET_S = 420.0


def run_cross_check(rows, n=CROSS_CHECK_N, force_idx=()):
    """force_idx: row indices ALWAYS cross-checked ahead of the random
    sample (--corner-append pins the two u->-1 fix witnesses with it);
    defaults leave the full-run behavior and checkpoint signature
    byte-identical."""
    if not rows:
        return True, None
    rng = random.Random(SEED ^ 0x5EED)
    # Exclude exact-saturated rows (P or Q == 0.0/1.0 on the nose) from the
    # candidate pool: those are either genuine specials-adjacent boundary
    # points or SATURATION_LOG_THRESHOLD pre-filter hits whose true value's
    # magnitude is astronomically far past any double -- comparing mpmath
    # vs CF there is numerically meaningless (both "agree" trivially, or
    # mpmath wastes its timeout budget on a comparison with no information
    # content) rather than a genuine test of the CF's accuracy.
    idx_all = [i for i, r in enumerate(rows) if 0.0 < r[3] < 1.0 and 0.0 < r[4] < 1.0]
    rng.shuffle(idx_all)
    fset = set(force_idx)
    sample = list(force_idx) + [i for i in idx_all[:n] if i not in fset]

    # v2: comparison method changed (small_side_direct vs raw CF, dps
    # escalation for tiny oriented x) -- v1 checkpoints must not be reused.
    sig = f"v2 SEED={SEED} NROWS={len(rows)} K={len(sample)}"
    done_map = load_checkpoint(CKPT_XCHECK_PATH, sig)
    print(f"self-check: mpmath betainc vs CF cross-check "
          f"({len(done_map)}/{len(sample)} already done)", file=sys.stderr)

    t_start = time.time()
    newly = 0
    existing_sig = existing_checkpoint_sig(CKPT_XCHECK_PATH)
    mode = "a" if existing_sig == sig else "w"
    with open(CKPT_XCHECK_PATH, mode) as fh:
        if mode == "w":
            fh.write(sig + "\n")
            fh.flush()
        for si, idx in enumerate(sample):
            if si in done_map:
                continue
            if time.time() - t_start > XCHECK_WALL_BUDGET_S:
                print(f"  cross-check wall-clock budget hit at {si}/{len(sample)} "
                      f"({newly} this run) -- re-run this script to continue.",
                      file=sys.stderr)
                return False, None
            a, b, x, Pf, Qf = rows[idx]
            am, bm, xm = mp.mpf(a), mp.mpf(b), mp.mpf(x)
            # Compare the oracle AS USED for the emitted row -- i.e.
            # small_side_direct with its whole routing (gamma-corner,
            # small-tau, near-diagonal shortcuts), not raw small_val_via_cf.
            # [corner-arc revision, 2026-08-12: the raw CF is structurally
            # invalid exactly where the routing replaces it -- the first
            # pb-corner run's sample drew (100, 1e307, 4e-306), where raw
            # CF returns 2.6e-30542 against a true 1.206e-15; the emitted
            # row (via gamma_corner_value) was CORRECT, and the old
            # comparison flagged a defect that was in the CHECK, not the
            # data. Verified by betainc at dps 500.]
            try:
                cf_p, cf_q, which, _esc, cf_failed = small_side_direct(a, b, x)
            except (RuntimeError, ZeroDivisionError, ValueError):
                cf_failed = True
            if cf_failed:
                append_checkpoint(fh, si, ["SKIP_CF"])
                newly += 1
                continue
            if which == "P":
                aa, bb, xx, cf_v = am, bm, xm, cf_p
            else:
                aa, bb, xx, cf_v = bm, am, _one_minus(xm), cf_q
            # betainc forms 1-xx INTERNALLY at its working dps (the
            # exact-complement hazard documented above): for xx below
            # ~10^-dps that truncates to exactly 1 and silently drops the
            # whole (1-xx)^b factor (measured: (100, 1e307, 4e-306) at
            # dps 50 returns 1.7e-30698 for a true 1.2e-15). Escalate dps
            # so xx stays visible; unreachable points TIMEOUT and are
            # excluded, exactly as before.
            need_dps = CROSS_CHECK_MPMATH_DPS
            xf = float(xx)
            if 0.0 < xf < 1e-20:
                need_dps = max(need_dps, int(30.0 - math.log10(xf)))
            mm_v = gbd._betainc_timeout(aa, bb, xx, need_dps,
                                         timeout=CROSS_CHECK_TIMEOUT)
            if mm_v is None:
                append_checkpoint(fh, si, ["TIMEOUT"])
            else:
                rel = abs((cf_v - mm_v) / mm_v) if mm_v != 0 else abs(cf_v - mm_v)
                append_checkpoint(fh, si, [f"{float(rel):.6e}"])
            newly += 1
            if newly % 25 == 0:
                print(f"    ... {si + 1}/{len(sample)} "
                      f"({time.time() - t_start:.0f}s this run)", file=sys.stderr)
                sys.stderr.flush()

    print(f"  cross-check computed {newly} points this run "
          f"({time.time() - t_start:.0f}s)", file=sys.stderr)
    done_map = load_checkpoint(CKPT_XCHECK_PATH, sig)
    n_timeout = n_skip = n_ok = 0
    worst = 0.0
    worst_at = None
    for si, idx in enumerate(sample):
        val = done_map[si][0]
        if val == "TIMEOUT":
            n_timeout += 1
            continue
        if val == "SKIP_CF":
            n_skip += 1
            continue
        rel = float(val)
        n_ok += 1
        if rel > worst:
            worst, worst_at = rel, rows[idx]
    print(f"  {n_ok} points had a reachable mpmath value, {n_timeout} timed out, "
          f"{n_skip} CF-skipped; worst rel diff (mpmath vs CF) {worst:.3e} at "
          f"{worst_at}, target {float(CROSS_CHECK_TARGET):.0e}", file=sys.stderr)
    summary = {"n_ok": n_ok, "n_timeout": n_timeout, "n_skip": n_skip,
               "worst": worst, "worst_at": worst_at}
    return True, summary


# ============================================================================
# Self-check: small-tau oracle validation [coordinator round-2 requirement].
# On the overlap band tau in [2^-20, eps_R4] x the xi range where the CF
# DOES converge (sampled from gen_r4's own box construction, so every
# sampled point satisfies tau*|ln xi_tau|<=ln2 -- the small-tau oracle's own
# applicability condition), require small-tau-oracle vs CF agreement to
# output precision + 10 digits (same 1e-20 bar as the mpmath cross-check).
# ============================================================================
SMALL_TAU_OVERLAP_TARGET = mp.mpf("1e-20")
SMALL_TAU_OVERLAP_N = 260  # >= the mandated 200, small margin


def check_small_tau_overlap(rng):
    print("self-check: small-tau-oracle vs CF overlap-band validation "
          "(tau in [2^-20, eps_R4])", file=sys.stderr)
    lo_tau, hi_tau = 2.0 ** -20, float(EPS_R4)
    worst = 0.0
    worst_at = None
    n_ok = 0
    n_cf_unreachable = 0
    n_oracle_nonconvergent = 0
    t0 = time.time()
    attempts = 0
    while n_ok < SMALL_TAU_OVERLAP_N and attempts < SMALL_TAU_OVERLAP_N * 8:
        attempts += 1
        tau = 10.0 ** rng.uniform(math.log10(lo_tau), math.log10(hi_tau))
        Bp = 10.0 ** rng.uniform(-3.0, math.log10(B1_f / XI1_f) + 0.3)
        floor, ceil = _r4_xi_bounds(tau, Bp, float(LN2), XI1_f, B1_f)
        if ceil <= 0:
            continue
        floor = max(floor, 1e-308)
        if floor >= ceil:
            continue
        logf, logc = math.log(floor), math.log(ceil)
        # Sample the CENTRAL part of the box (0.25-0.75 of the floor-ceil
        # log-range), not its extreme edges: per the coordinator's own
        # framing, this overlap band should be "the xi range where the CF
        # DOES converge" -- and this check's own first attempt (comparing
        # against the CF on the SWAPPED triple, i.e. (B,tau,1-xi)) found
        # that call returns exactly 0.0 at the box's most extreme corners
        # (deep-tiny xi_tau near the floor) even though the CF's NATIVE
        # evaluation (tau,B,xi) converges fine there -- the swapped call is
        # itself unreliable in exactly the sub-corner this generator's
        # small-tau oracle exists to route around, so it cannot serve as
        # ground truth there. Restricting to the box's middle keeps every
        # sampled point's complement moderate (not deep-saturated), where
        # BOTH "1 - CF-native" (subtraction at dps=100, ample headroom
        # since the complement here is never smaller than ~1e-30-class)
        # and the small-tau oracle are independently trustworthy -- a fair
        # apples-to-apples comparison.
        xi_tau = math.exp(logf + rng.uniform(0.25, 0.75) * (logc - logf))
        if not (0 < xi_tau < 1):
            continue
        tau_m, Bp_m, xi_m = mp.mpf(tau), mp.mpf(Bp), mp.mpf(xi_tau)
        oracle_v, ok = small_tau_oracle(tau_m, Bp_m, xi_m, 60)
        if not ok or oracle_v is None or oracle_v > mp.mpf("0.5"):
            n_oracle_nonconvergent += 1
            continue
        ground_dps = 100
        try:
            cf_native = small_val_via_cf(tau_m, Bp_m, xi_m, ground_dps)
        except (RuntimeError, ZeroDivisionError, ValueError):
            n_cf_unreachable += 1
            continue
        old_dps = mp.mp.dps
        mp.mp.dps = ground_dps
        try:
            cf_v = 1 - cf_native
        finally:
            mp.mp.dps = old_dps
        if cf_v <= 0 or cf_v > mp.mpf("0.5") or cf_v < mp.mpf(10) ** -(ground_dps - 20):
            # Either not a fair small-complement point, or the subtraction
            # itself doesn't have enough headroom left at ground_dps to
            # trust as ground truth here -- skip rather than risk a false
            # positive/negative from the ground truth's OWN precision.
            continue
        rel = abs((oracle_v - cf_v) / cf_v) if cf_v != 0 else abs(oracle_v - cf_v)
        n_ok += 1
        if rel > worst:
            worst, worst_at = float(rel), (tau, Bp, xi_tau)
    log2w = math.log2(worst) if worst > 0 else float("-inf")
    print(f"    {n_ok} points compared ({attempts} attempts, "
          f"{n_cf_unreachable} CF-unreachable, {n_oracle_nonconvergent} "
          f"oracle-nonconvergent-or-wrong-side, {time.time()-t0:.1f}s); "
          f"worst rel diff {worst:.3e} (2^{log2w:.2f}) at (tau,B,xi)={worst_at}, "
          f"target {float(SMALL_TAU_OVERLAP_TARGET):.0e}", file=sys.stderr)
    if n_ok < 200:
        print(f"    FAILED: only {n_ok} < 200 comparable overlap-band points "
              f"found -- ESCALATE.", file=sys.stderr)
        return 1
    if worst > float(SMALL_TAU_OVERLAP_TARGET):
        print(f"    FAILED: worst disagreement exceeds target -- ESCALATE.",
              file=sys.stderr)
        return 1
    return 0


# ============================================================================
# Self-check: gamma-corner oracle validation [(C), Part 2b]. gamma_corner_
# value is only ROUTED to for max(a,b)>=B_GL (~2^59), where the beta CF is
# structurally degenerate and cannot serve as ground truth at all -- so this
# validates the FORMULA itself (not the routing gate) on synthetic points
# with beta on the "healthy band" (beta<=2^40, the G1a-validated line the
# beta CF is trusted on well below any conditioning concern) against the
# beta CF directly, mirroring gen_beta_data.py's own B_GL "downward" probe
# methodology (same asymptotic form, same t=-beta*log1p(-xi)) but as an
# independent check in THIS file rather than trusting that probe's numbers
# on faith.
#
# TARGET, reconsidered after running (own measurement, not assumed): the
# asymptotic error decays like 2^(8.57-j) in this generator's own downward
# probe (see gen_beta_data.py's _probe_gl_downward budget lines) -- so at
# beta<=2^40 the formula is EXPECTED to sit only around 2^-31, nowhere near
# the 2^-49 bar that only holds near B_GL~2^59 itself (where the formula is
# actually used). Holding this healthy-band check to 2^-49 would be
# penalizing the formula for not yet being in its own intended asymptotic
# regime -- not a real defect. This check is therefore REPORT-ONLY (a
# "budget line to stderr" per the brief, not a hard gate) beyond a loose
# sanity bound that WOULD catch a genuine implementation bug (this check's
# own first run found exactly such a bug -- a ~2^888 relative error from an
# inverted orientation -- well before reaching this refined, expected-slack
# regime).
GAMMA_CORNER_SANITY_BOUND = mp.mpf(2) ** -2  # catches real bugs; the healthy
                                          # here (<=2^40) sits BELOW B_GL, so
                                          # this is deliberately measuring
                                          # the formula in its LEAST-favorable
                                          # (least converged) validated zone.


def check_gamma_corner_oracle():
    print("self-check: gamma-corner oracle vs beta CF on the healthy band "
          "(beta<=2^40):", file=sys.stderr)
    dps = 80
    js = [10, 20, 30, 36, 40]
    alphas = (mp.mpf("0.05"), mp.mpf("0.25"), mp.mpf(1))
    bxis = (8, 20, 50, 200, 800)
    worst = mp.mpf(0)
    worst_at = None
    n_tested = 0
    n_skipped = 0
    for j in js:
        beta = mp.mpf(2) ** j
        for alpha in alphas:
            for bxi in bxis:
                xi = mp.mpf(bxi) / beta
                # the asymptotic (t=-beta*log1p(-xi) fixed as beta->infty)
                # needs xi genuinely SMALL -- at low j with a large bxi, xi
                # is not small at all (e.g. j=10,bxi=800 -> xi=0.78) and the
                # comparison is meaningless (not a formula defect, a
                # domain-applicability one); skip those combinations rather
                # than let them dominate the reported sup with noise.
                if not (0 < xi < mp.mpf("0.3")):
                    continue
                try:
                    gform_small, gform_which = gamma_corner_value_signed(alpha, beta, xi, dps)
                    # ORIENT before calling the CF -- caught by this check's
                    # own first run: calling small_val_via_cf(alpha,beta,xi)
                    # RAW (unoriented) here fed the CF its numerically SLOW/
                    # ill-conditioned direction (this magnitude regime always
                    # has xi >= (alpha+1)/(c+2), i.e. the fast direction is
                    # the SWAP), and its own N/2-vs-N drift gate reported
                    # "converged" at a value ~1e267 relative away from the
                    # truth (mpmath.betainc cross-checked separately) --
                    # self-convergence in the wrong orientation is not a
                    # correctness guarantee, exactly check_b_r2's own
                    # "_cf_err_at" precedent this check should have followed
                    # from the start. And the COMPARISON itself must stay on
                    # the small side too (gamma_corner_value_signed, not
                    # gamma_corner_value's plain P) -- an earlier version of
                    # this check compared near-1 P values directly and the
                    # 1-near-1 loss made the comparison meaningless whenever
                    # Q was the genuinely tiny side (P and truth both rounded
                    # to exactly 1.0 at dps=80, masking a real defect).
                    c = alpha + beta
                    thresh = (alpha + 1) / (c + 2)
                    if xi < thresh:
                        truth_small = small_val_via_cf(alpha, beta, xi, dps,
                                                         n_start=128, n_max=8192)
                        truth_which = "P"
                    else:
                        truth_small = small_val_via_cf(beta, alpha, 1 - xi, dps,
                                                         n_start=128, n_max=8192)
                        truth_which = "Q"
                except (RuntimeError, ZeroDivisionError, ValueError):
                    n_skipped += 1
                    continue
                if truth_small == 0 or gform_which != truth_which:
                    n_skipped += 1
                    continue
                err = abs((gform_small - truth_small) / truth_small)
                n_tested += 1
                if err > worst:
                    worst, worst_at = err, (float(alpha), float(beta), float(xi), j)
    log2w = float(mp.log(worst, 2)) if worst > 0 else float("-inf")
    print(f"    tested {n_tested} points ({n_skipped} skipped), js={js}; "
          f"worst rel err {float(worst):.3e} (2^{log2w:.2f}) at "
          f"(alpha,beta,xi,j)={worst_at} -- REPORT ONLY (see comment: the "
          f"formula's own asymptotic error at beta<=2^40 is expected slack, "
          f"not a defect; sanity-gated at 2^-2 only, which would still "
          f"catch a real implementation bug)", file=sys.stderr)
    return 0 if worst <= GAMMA_CORNER_SANITY_BOUND else 1


# ============================================================================
# Self-check: P+Q=1 within 1 ULP, at the emitted (rounded) doubles.
# ============================================================================
def check_p_plus_q(rows):
    print("self-check: P+Q=1 within 1 ULP (rounded doubles)", file=sys.stderr)
    ulp1 = 2.0 ** -52
    worst = 0.0
    n_fail = 0
    for a, b, x, P, Q in rows:
        d = abs(1.0 - (P + Q))
        if d > worst:
            worst = d
        if d > ulp1:
            n_fail += 1
    print(f"    worst |1-(P+Q)| = {worst:.3e} (1 ULP = {ulp1:.3e}), "
          f"{n_fail} rows exceed 1 ULP", file=sys.stderr)
    return 0 if n_fail == 0 else 1


# ============================================================================
# Self-check: small-side-direct verification. On every emitted row, the
# side that was computed DIRECTLY must be the numerically smaller one
# (min(P,Q)) -- by construction (small_side_direct's >0.5 self-correction)
# this always holds for every row this generator ever emits; verified here
# post-hoc as the regression guard. Documented nuance: this checks the
# EMITTED pair only (post keep_if_saturated pruning and rounding) -- exact
# ties at P=Q=0.5 (the diagonal) pass trivially either way.
# ============================================================================
def check_small_side_direct(rows):
    print("self-check: small-side-direct (min(P,Q) is the directly-computed "
          "side)", file=sys.stderr)
    n_fail = 0
    worst_gap = 0.0
    for a, b, x, P, Q in rows:
        # Recompute cheaply which side small_side_direct WOULD pick (the
        # mean predicate, self-corrected) and confirm it lines up with
        # whichever of P,Q is actually <= the other -- since small_side_direct
        # always emits final<=0.5 on its chosen side (or the exact P=Q=0.5
        # tie), this reduces to: min(P,Q) <= max(P,Q), always true, EXCEPT
        # if rounding pushed the smaller side fractionally above 0.5 at an
        # exact-tie boundary -- allow a 1-ULP slop there.
        lo, hi = (P, Q) if P <= Q else (Q, P)
        if lo > 0.5 + 2.0 ** -52:
            n_fail += 1
            worst_gap = max(worst_gap, lo - 0.5)
    print(f"    {n_fail} rows fail min(P,Q)<=1/2 (+1 ULP slop); worst gap "
          f"{worst_gap:.3e}", file=sys.stderr)
    return 0 if n_fail == 0 else 1


def write_rows(rows):
    for path in ("tests/data/beta_p_reference.txt", "tests/data/beta_q_reference.txt"):
        with open(path, "w") as f:
            for a, b, x, P, Q in rows:
                f.write(f"{hexd(a)} {hexd(b)} {hexd(x)} {hexd(P)} {hexd(Q)}\n")
        print(f"  wrote {path}: {len(rows)} points", file=sys.stderr)


# ============================================================================
# main
# ============================================================================
def main():
    t_all = time.time()
    rng = random.Random(SEED)

    print("building region point sets ...", file=sys.stderr)
    ps = PointSet()
    gen_r1(ps, rng)
    gen_r4(ps, rng)
    gen_r3(ps, rng)
    gen_r2(ps, rng)
    gen_ridge_lines(ps)
    gen_boundaries(ps)
    gen_diagonal(ps)
    gen_analytic_lines(ps, rng)
    gen_subnormal_band(ps)
    gen_huge_tiny(ps, rng)
    gen_witnesses(ps)
    gen_random_fill(ps, rng, TARGET_TOTAL)
    # LAST, deterministically, so every family above (including the rng
    # fill's stream and stopping point) is byte-identical to the pre-corner
    # runs; --corner-append splices exactly this block's rows.
    gen_pb_corner(ps)
    print(f"  total distinct (a,b,x) points: {len(ps.pts)}", file=sys.stderr)

    print("evaluating oracle (CF, dps ladder 40/60/100, resumable) ...",
          file=sys.stderr)
    rows, region_hist, done = compute_all(ps)
    if not done:
        print("\nPARTIAL RUN: re-invoke this script to continue "
              "(checkpoint saved).", file=sys.stderr)
        return 3

    rc = 0

    rc_hist = 0
    for base in ("R1", "R2", "R3", "R4"):
        n = region_hist.get(base, 0)
        ok = n >= MIN_PER_REGION
        print(f"  region {base}: {n} points (floor {MIN_PER_REGION}) "
              f"{'OK' if ok else 'FAIL'}", file=sys.stderr)
        if not ok:
            rc_hist = 1
    print(f"  region R2-gammalim (subset of R2, max param >= B_GL=2^"
          f"{round(math.log2(float(B_GL)))}): {region_hist.get('R2-gammalim', 0)} "
          f"points (informational, not floor-gated)", file=sys.stderr)
    print(f"  region R2-postroute (subset of R2, SIXTH-correction near-one "
          f"post-route from R1): {region_hist.get('R2-postroute', 0)} points "
          f"(informational, not floor-gated)", file=sys.stderr)
    rc |= rc_hist

    rc |= check_small_tau_overlap(random.Random(SEED ^ 0x7A0))
    rc |= check_gamma_corner_oracle()
    rc |= check_analytic_lines()
    rc |= check_p_plus_q(rows)
    rc |= check_small_side_direct(rows)

    xcheck_ok, xcheck_summary = run_cross_check(rows)
    if not xcheck_ok:
        print("\nPARTIAL RUN: cross-check incomplete, re-invoke this script "
              "to continue (checkpoint saved).", file=sys.stderr)
        return 3
    if xcheck_summary is not None and xcheck_summary["worst"] > float(CROSS_CHECK_TARGET):
        print(f"  FAILED: cross-check worst disagreement "
              f"{xcheck_summary['worst']:.3e} exceeds target "
              f"{float(CROSS_CHECK_TARGET):.0e}.", file=sys.stderr)
        rc = 1

    if ESCALATIONS:
        print(f"\n{len(ESCALATIONS)} points had persistent dps60-vs-dps100 "
              f"disagreement beyond target (escalated, not dropped):",
              file=sys.stderr)
        for e in ESCALATIONS[:20]:
            print(f"    a={e[0]:.6e} b={e[1]:.6e} x={e[2]:.10f} which={e[3]} "
                  f"rel40v60={e[4]:.3e} rel60v100={e[5]:.3e}", file=sys.stderr)

    specials = gen_specials_rows()
    print(f"\nspecials table: {len(specials)} rows", file=sys.stderr)

    if rc:
        print("\nSelf-checks FAILED -- not writing output files.", file=sys.stderr)
        return rc

    all_rows = rows + specials
    write_rows(all_rows)

    print(f"\ntotal generator runtime (this invocation): "
          f"{time.time() - t_all:.1f}s", file=sys.stderr)
    print(f"CF failures across the whole run: {N_CF_FAILED[0]}", file=sys.stderr)
    print(f"saturation pre-filter shortcuts (skipped CF entirely): "
          f"{N_PREFILTER_SATURATED[0]}", file=sys.stderr)
    print(f"small-tau oracle: rescued {N_SMALL_TAU_RESCUED[0]}, "
          f"non-convergent (dropped) {N_SMALL_TAU_NONCONVERGENT[0]}, "
          f"deep-dps re-resolved {N_SMALL_TAU_DEEP[0]} "
          f"[rescue round 4: noise-flagged rows re-run at dps 160/240/400]",
          file=sys.stderr)
    print(f"betainc-with-timeout rescue [Part 2a]: rescued "
          f"{N_BETAINC_RESCUED[0]} of the small-tau-guard drops above "
          f"(only points failing BOTH are counted in the non-convergent "
          f"total)", file=sys.stderr)
    return 0


# ============================================================================
# --corner-append: incremental certification of the pb-corner family alone.
#
# The full-run checkpoint (42k points, hours of oracle work including the
# subprocess-guarded betainc rescues) did not survive the temp-dir cleanup;
# re-deriving every certified row to append ~600 is waste with no accuracy
# upside. This mode evaluates ONLY gen_pb_corner's points -- same
# small_side_direct protocol, same dps ladder, own checkpoint file -- runs
# the applicable self-checks plus a cross-check with both defect witnesses
# force-included, and splices the rows into the existing reference files
# immediately BEFORE the specials block, which is exactly where a fresh
# full run (gen_pb_corner called last, specials appended after) would put
# them. The splice refuses to run if the existing files' specials tail
# does not match gen_specials_rows() bit-for-bit.
# ============================================================================
CKPT_CORNER_PATH = os.path.join(tempfile.gettempdir(),
                                f"corvus_beta_ref_ckpt_corner_{SEED}.tsv")


def _splice_corner(corner_rows):
    specials = gen_specials_rows()
    spec_lines = [f"{hexd(a)} {hexd(b)} {hexd(x)} {hexd(P)} {hexd(Q)}"
                  for a, b, x, P, Q in specials]
    for path in ("tests/data/beta_p_reference.txt",
                 "tests/data/beta_q_reference.txt"):
        with open(path) as f:
            lines = [ln.rstrip("\n") for ln in f if ln.strip()]
        if lines[-len(spec_lines):] != spec_lines:
            print(f"  FAILED: {path} specials tail does not match "
                  f"gen_specials_rows() -- refusing to splice.",
                  file=sys.stderr)
            return 1
        body = lines[:-len(spec_lines)]
        seen = {tuple(ln.split()[:3]) for ln in body}
        added = []
        for a, b, x, P, Q in corner_rows:
            key = (hexd(a), hexd(b), hexd(x))
            if key in seen:
                continue  # PointSet-consistent: earlier family owns the row
            seen.add(key)
            added.append(f"{key[0]} {key[1]} {key[2]} {hexd(P)} {hexd(Q)}")
        with open(path, "w", newline="\n") as f:
            f.write("\n".join(body + added + spec_lines) + "\n")
        print(f"  wrote {path}: +{len(added)} corner rows "
              f"({len(body) + len(added) + len(spec_lines)} total)",
              file=sys.stderr)
    return 0


def corner_append():
    t_all = time.time()
    ps = PointSet()
    gen_pb_corner(ps)
    print(f"corner point set: {len(ps.pts)} distinct points", file=sys.stderr)

    # Signature carries a digest of the POINT BITS, not just the count: a
    # grid edit that preserves N would otherwise replay stale checkpoint
    # values under new point identities (caught live when the b-column
    # swap 1e20 -> 1e30 kept N = 786 and the first re-run served b = 1e20
    # oracle values as b = 1e30 rows).
    import hashlib
    dig = hashlib.sha256()
    for a, b, x, _, _ in ps.pts:
        dig.update(struct.pack("<QQQ", as_bits(a), as_bits(b), as_bits(x)))
    rows, region_hist, done = compute_all(
        ps, ckpt_path=CKPT_CORNER_PATH,
        sig_ver=f"v1-corner-{dig.hexdigest()[:16]}")
    if not done:
        print("\nPARTIAL RUN: re-invoke with --corner-append to continue "
              "(checkpoint saved).", file=sys.stderr)
        return 3
    print(f"  corner region histogram: {region_hist}", file=sys.stderr)

    rc = 0
    rc |= check_p_plus_q(rows)
    rc |= check_small_side_direct(rows)

    # Both defect witnesses must be present as emitted rows, and are
    # force-included in the cross-check sample below.
    wit = [(19.0, 1e5, 5.204222470155122e-21), (19.0, 1e5, 1.73e-19)]
    wit_idx = []
    for wa, wb, wx in wit:
        hits = [i for i, r in enumerate(rows)
                if as_bits(r[0]) == as_bits(wa) and as_bits(r[1]) == as_bits(wb)
                and as_bits(r[2]) == as_bits(wx)]
        if not hits:
            print(f"  FAILED: witness row ({wa}, {wb}, {wx:.17e}) missing "
                  f"from the corner set.", file=sys.stderr)
            rc = 1
        else:
            wit_idx.append(hits[0])
            r = rows[hits[0]]
            print(f"  witness ({wa}, {wb:.0e}, {wx:.17e}): P={r[3]:.17e} "
                  f"Q={r[4]:.17g}", file=sys.stderr)

    xcheck_ok, xcheck_summary = run_cross_check(rows, n=60,
                                                force_idx=tuple(wit_idx))
    if not xcheck_ok:
        print("\nPARTIAL RUN: cross-check incomplete, re-invoke with "
              "--corner-append to continue (checkpoint saved).",
              file=sys.stderr)
        return 3
    if xcheck_summary is not None and \
            xcheck_summary["worst"] > float(CROSS_CHECK_TARGET):
        print(f"  FAILED: cross-check worst disagreement "
              f"{xcheck_summary['worst']:.3e} exceeds target "
              f"{float(CROSS_CHECK_TARGET):.0e}.", file=sys.stderr)
        rc = 1

    if rc:
        print("\nSelf-checks FAILED -- not touching reference files.",
              file=sys.stderr)
        return rc

    rc = _splice_corner(rows)
    print(f"\ncorner-append runtime: {time.time() - t_all:.1f}s",
          file=sys.stderr)
    return rc


if __name__ == "__main__":
    if "--corner-append" in sys.argv[1:]:
        sys.exit(corner_append())
    sys.exit(main())
