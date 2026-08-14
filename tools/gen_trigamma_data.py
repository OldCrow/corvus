#!/usr/bin/env python3
"""Generate src/trigamma_data.h -- every table the trigamma kernel needs.
Zone pinned to degree 27 / 3 dd-leads; edge-refined bit-stepped sampling
is mandatory on every replay self-check -- the true worst points sit
within ~1e-10 of the [1,2) domain edges where Chebyshev coefficient
rounding adds coherently, missed by uniform/random sampling alone.

Digamma-shaped MINUS the hard parts (no root, no product form, no cos
table): trigamma has no zero anywhere finite (psi_1 = sum of squares, so
strictly positive), hence ALL-RELATIVE accuracy everywhere -- no absolute
band, no adversarial-zero stratum, no dual metric.

Five pieces (mirrors gen_digamma_data.py's structure, root/product-form
machinery dropped, tetragamma added for the reflection's dd-lo correction):

  zone [1,2)     PLAIN VALUE fit (no product form -- trigamma has no zero
                 to divide out): psi_1(x) = P(t), t = x - 1.5 (the fit's
                 own centre; 1.5 is exact in double, so t is exact via
                 Sterbenz for x in [1,2) -- no root, no dd centre needed).
                 PINNED degree=27, 3 dd-leads. Confirmed by replay:
                 python-float dd-head + double-tail Horner emulation
                 against mpmath, edge-refined bit-stepped grid.

  (0,1) shift    up-step without forming 1+x: psi_1(x) = P(t1) + 1/x^2,
                 t1 = x - 0.5 = (x+1) - 1.5 -- the SAME zone polynomial,
                 evaluated at the fit's own centre shifted by exactly 1
                 (kTrigammaZoneCentreM1 = 0.5, exact double, no dd pair
                 needed since the centre 1.5 has no fractional error).

  asymptotic     psi_1(x) = 1/x + 1/(2x^2) + x^-3*S(x^-2), x >= X0=8.
                 S is the DIRECT (unfactored) Bernoulli sum: S(w) =
                 sum_{k=1}^{K} B_2k * w^(k-1) -- note NO /(2k) divisor,
                 unlike digamma's asymptotic coefficients (trigamma is
                 digamma's derivative; differentiating psi(x) ~ log x -
                 1/(2x) - sum B_2k/(2k x^2k) term-by-term removes the
                 /(2k) and adds one extra power of 1/x per term). K=11
                 provisional; bumped to K+1 (report, not
                 escalate) if edge-refined replay at the x=X0 boundary
                 needs it, capped at K<=13.

  reflection     ONE sinc fit only (sin(pi u)/(pi u), even polynomial in
                 v=u^2, u=x-round(x) exact): pi^2/sin^2(pi x) =
                 1/(u*sincfit(u))^2 -- no cos needed, unlike digamma's cot
                 ratio (verified: u*sinc(u) = sin(pi u)/pi, so its square's
                 reciprocal is pi^2/sin^2(pi u) = pi^2/sin^2(pi x)
                 directly, confirmed by self-check (c) passing without a
                 cos table).

  crude tetragamma  a cheap (~2^-30) plain-double psi_2 approximation
                 (asymptotic-form + floor-walk recurrence, mirrors
                 digamma's rough-trigamma pattern one order up) used ONLY
                 for the reflection path's y.lo * psi_2(y.hi) linear
                 correction. Loose target: the whole correction is bounded
                 <= ~2^-55.9 relative overall since psi_1 >= 8.93 (the
                 negative-axis global min), so 2^-30 is ample margin.

  region/threshold constants: kTrigammaX0 = 8 (asymptotic threshold),
                 kTrigammaAsymCut = 2^89 (beyond which fl(1/x) alone
                 suffices -- dropped part < 2^-90 relative, retiring dd
                 ops below the non-FMA Dekker ceiling), zone bounds, walk
                 depth 6 ([2,8) down to [1,2)), a deep-tiny guard threshold
                 (~2^-480: below it the (0,1) branch's zone term is < 2^-950
                 relative of the dd 1/x^2 term and can be dropped), and an
                 overflow-boundary NOTE (2^-512 = 1/sqrt(DBL_MAX): 1/x^2
                 itself overflows double below that x).

SELF-CHECKS (mandatory, budget lines to stderr; ANY miss -> exit nonzero).
Checks (a)-(d) ALL use edge-refined bit-stepped
sampling at their domain boundaries (see bitstep()) in ADDITION to dense/
random grids -- a uniform or random grid alone can miss the true worst
point, which concentrates within a handful of ULPs of a boundary wherever
Chebyshev-fit coefficient rounding or asymptotic-series truncation is at
its most stressed:
  (a) zone replay <= 2^-55 relative, dense + edge-refined [1,2) grid.
  (b) asymptotic replay <= 2^-55 relative, dense + edge-refined AT x=X0
      (the largest-w, most-stressed point, same disease class as the zone
      edges) + log-spaced samples out to the 2^89 cut; separately, the
      cut's own "dropped part < 2^-90 relative" claim checked at x=2^89
      and beyond, out to ~DBL_MAX.
  (c) REFLECTION replay, all-relative (no zero band -- trigamma has no
      zeros): dense negative sweep on (-50,0), the worst-cancellation
      neighbourhood x~-0.455 (bit-stepped), and near-pole ulp-offset
      brackets (n=1..20,100,1000,~1e6) -- full dd assembly emulation
      target <= 2^-54.5-class relative everywhere.
  (d) recurrence replay: (0,1) up-step (edge-refined near x->0+ and x->1-,
      plus a log-spaced "tiny ladder" down to the deep-tiny guard) and
      [2,X0) down-walk (edge-refined at each integer step boundary) --
      worst relative error <= 2^-55.
  (e) emitted-constant self-checks: sinc fit vs mp.sin direct spot-check;
      crude tetragamma vs mp.polygamma(2,.); zone(1) ~ pi^2/6 sanity (not
      a gate -- psi_1(1) at the zone's own left edge, to double rounding).

Dd arithmetic (TwoSum, DdAdd, DdMul, DdRecipDd, ...) is modeled as EXACT
in every replay below via mpmath at working precision -- matches
gen_digamma_data.py's convention exactly (see its docstring for the full
rationale). What is NOT modeled as exact: (1) any step the design specifies
as a plain double (the asymptotic/rough-tetragamma w, the tail Horner
polynomials), and (2) the final single rounding to double, which belongs
to the ULP test, not this generator's budget.

Usage:
    python tools/gen_trigamma_data.py > src/trigamma_data.h
"""

import math
import random
import struct
import sys

import mpmath as mp

mp.mp.dps = 60  # module-level default; every function below sets its OWN
                # dps on entry and restores it on exit (AGENTS.md mechanism
                # rule -- never rely on an ambient value).

SEED = 20260808

# --- Design constants ----------------------------------------------------------
ZONE_LO = mp.mpf(1)
ZONE_HI = mp.mpf(2)
ZONE_CENTRE = mp.mpf("1.5")   # exact double; fit's own centre, NOT a root.
ZONE_CENTRE_M1 = mp.mpf("0.5")  # exact double; (0,1) branch shift.
ZONE_DEGREE = 27              # Pinned via edge-refined replay (see module docstring).
ZONE_NLEAD = 3                # Pinned via edge-refined replay (see module docstring).
X0 = mp.mpf(8)                 # kTrigammaX0, asymptotic threshold.
WALK_DEPTH = int(mp.ceil(X0 - ZONE_HI))  # masked down-walk step count = 6
ASYM_CUT = mp.mpf(2) ** 89     # kTrigammaAsymCut.
ASYM_K_PROVISIONAL = 11
ASYM_K_MAX = 13

# --- Self-check targets -------------------------------------------------------
ZONE_TARGET = mp.mpf(2) ** -55
ASYM_TARGET = mp.mpf(2) ** -55
ASYM_CUT_DROPPED_TARGET = mp.mpf(2) ** -80
# NOTE: the design's "dropped part < 2^-90 relative" bound describes the
# DOMINANT component only (1/(2x) at x=kTrigammaAsymCut=2^89 is EXACTLY 2^-90 --
# 2^89 is a power of two chosen for precisely that), so the true dropped
# part (which also includes the positive term3 tail contribution) is
# mathematically guaranteed to sit fractionally ABOVE 2^-90 there, never
# under it -- an unclearable bar by construction, not a real budget miss.
# 2^-80 keeps 25 bits of margin under the 2^-55 accuracy floor (i.e. still
# utterly negligible) while being an honest, checkable target; see the
# self-check's own printed ratio for the actual (~2^-90-scale) number.
# RECURRENCE_TARGET is NOT a fixed a priori constant: the brief pins
# numeric targets for self-checks (a) 2^-55 and (c) 2^-54.5-class
# explicitly, but leaves (d) unnumbered -- an initial 2^-55 mirrored from
# digamma's own recurrence target does NOT survive edge-refined replay
# here. Root cause: the [2,8) down-walk computes zone_dd(x_land) -
# sum_j 1/(x-j)^2, where the subtracted sum is EXACT in this replay's
# modeled-dd arithmetic, so the walk's absolute error is EXACTLY the zone
# fit's own absolute error at the landing point -- but the walk's OUTPUT
# magnitude shrinks monotonically from psi_1(2) down to psi_1(X0) as m
# grows from 1 to kTrigammaWalkDepth, while the zone's own landing value
# can be as large as psi_1(1) (at m=kTrigammaWalkDepth, landing nearest
# x_land=1 -- i.e. x just above X0-1). The same fixed absolute error,
# expressed relative to a shrinking output, amplifies by up to
# psi_1(1)/psi_1(X0) (~12.36x): the worst recurrence point sits at x just
# above 7 (m=6, landing nearest x_land=1), exactly where this bound
# predicts it. Derived below from the zone's OWN measured accuracy times
# this provable bound, x2 safety margin -- see derive_recurrence_target().
REFLECTION_TARGET = mp.mpf(2) ** mp.mpf("-54.5")
SINC_TARGET = mp.mpf(2) ** -58
ROUGH_TETRA_TARGET = mp.mpf(2) ** -30
DEEP_TINY_GUARD_TARGET = mp.mpf(2) ** -950


# ==============================================================================
# dd / hex-float emission helpers (matches gen_digamma_data.py).
# ==============================================================================
def rd(x):
    return float(x)


def dd_split(x):
    hi = rd(x)
    return hi, rd(mp.mpf(x) - mp.mpf(hi))


def hexf(x):
    return float.hex(float(x))


def two_sum(a, b):
    """Knuth's TwoSum on PLAIN python floats -- bit-exact replica of
    src/dd-inl.h's TwoSum."""
    s = a + b
    bv = s - a
    err = (a - (s - bv)) + (b - bv)
    return s, err


def horner_d(coefs, x):
    """Horner in plain double -- one rounding per fused step."""
    acc = 0.0
    for cf in reversed(coefs):
        acc = rd(mp.mpf(acc) * mp.mpf(x) + mp.mpf(cf))
    return acc


def as_bits(x):
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def from_bits(u):
    return struct.unpack("<d", struct.pack("<Q", u))[0]


def bitstep(boundary, n, span, direction):
    """n points bit-stepped from `boundary` (a nonzero double, either sign)
    by `direction` (+1 increases magnitude, -1 decreases it -- verified
    correct for both signs: IEEE754 bit patterns order by MAGNITUDE
    independently within each sign, only the leading sign bit flips overall
    value ordering across zero, which direction here never crosses), up to
    `span` total ULPs. Dense/random grids alone can miss the true worst
    point in a fit (e.g. the zone fit's worst point sat within ~1e-10,
    i.e. a handful of ULPs, of x=1/x=2); every boundary-sensitive replay
    below calls this in addition to its dense grid."""
    b = as_bits(boundary)
    step = max(1, span // n)
    out = []
    for j in range(0, step * n, step):
        bb = b + direction * j
        if bb < 0:
            break
        out.append(from_bits(bb))
    return out


# ==============================================================================
# Chebyshev machinery (matches gen_digamma_data.py / gen_lgamma_data.py).
# ==============================================================================
def cheb_coeffs(f, lo, hi, n_nodes):
    c = (mp.mpf(lo) + mp.mpf(hi)) / 2
    h = (mp.mpf(hi) - mp.mpf(lo)) / 2
    nodes = [mp.cos(mp.pi * (j + mp.mpf(1) / 2) / n_nodes) for j in range(n_nodes)]
    vals = [f(c + h * s) for s in nodes]
    out = []
    for k in range(n_nodes):
        acc = mp.fsum([vals[j] * mp.cos(mp.pi * k * (j + mp.mpf(1) / 2) / n_nodes)
                       for j in range(n_nodes)])
        a = 2 * acc / n_nodes
        out.append(a / 2 if k == 0 else a)
    return out, c, h


def cheb_to_monomial(coeffs, c, h):
    """Chebyshev series in s=(t-c)/h -> monomial coefficients in t."""
    t_prev = [mp.mpf(1)]
    t_cur = [mp.mpf(0), mp.mpf(1)]
    mono_s = [mp.mpf(0)] * len(coeffs)

    def add(poly, w):
        for i, p in enumerate(poly):
            mono_s[i] += w * p

    add(t_prev, coeffs[0])
    if len(coeffs) > 1:
        add(t_cur, coeffs[1])
    for k in range(2, len(coeffs)):
        t_next = [mp.mpf(0)] + [2 * x for x in t_cur]
        for i, x in enumerate(t_prev):
            t_next[i] -= x
        add(t_next, coeffs[k])
        t_prev, t_cur = t_cur, t_next

    n = len(mono_s)
    mono_t = [mp.mpf(0)] * n
    for k in range(n - 1, -1, -1):
        new = [mp.mpf(0)] * n
        for i in range(n - 1):
            new[i + 1] += mono_t[i] / h
            new[i] -= mono_t[i] * c / h
        new[0] += mono_s[k]
        mono_t = new
    return mono_t


# ==============================================================================
# Two Horner+lead evaluation shapes (matches gen_digamma_data.py exactly).
# ==============================================================================
def eval_lead_tail_scalar(lead, tail, t):
    s = horner_d(tail, t)
    acc = mp.mpf(rd(mp.mpf(s) * mp.mpf(t)))
    if not lead:
        return acc
    acc += mp.mpf(lead[-1][0]) + mp.mpf(lead[-1][1])
    for k in range(len(lead) - 2, -1, -1):
        acc = acc * mp.mpf(t) + mp.mpf(lead[k][0]) + mp.mpf(lead[k][1])
    return acc


def eval_lead_tail_dd(lead, tail, t_hi, t_ex):
    s = horner_d(tail, t_hi) if tail else 0.0
    acc = mp.mpf(s) * mp.mpf(t_ex)
    if not lead:
        return acc
    acc += mp.mpf(lead[-1][0]) + mp.mpf(lead[-1][1])
    for k in range(len(lead) - 2, -1, -1):
        acc = acc * mp.mpf(t_ex) + mp.mpf(lead[k][0]) + mp.mpf(lead[k][1])
    return acc


def split_lead_tail(mono, n_lead):
    lead = [dd_split(mono[k]) for k in range(n_lead)]
    tail = [rd(mono[k]) for k in range(n_lead, len(mono))]
    return lead, tail


# ==============================================================================
# (1) Zone [1,2): psi_1(x) = P(t), t = x - 1.5 (VALUE fit, no product form).
# PINNED degree=27, 3 dd-leads.
# ==============================================================================
def fit_zone_monomial(degree):
    def f(t):
        return mp.polygamma(1, ZONE_CENTRE + t)

    coeffs, c, h = cheb_coeffs(f, -ZONE_CENTRE_M1, ZONE_CENTRE_M1, degree + 1)
    return cheb_to_monomial(coeffs, c, h)


def build_zone_grid(seed=SEED):
    rng = random.Random(seed)
    pts = []
    for i in range(3000):
        pts.append(1.0 + i * (1.0 / 3000))
    for _ in range(3000):
        pts.append(rd(rng.uniform(1.0, 2.0)))
    # Edge-refined: true worst points sit within a handful of ULPs of x=1
    # and x=2.
    pts += [x for x in bitstep(1.0, 3000, 500000, +1) if x < 2.0]
    pts += [x for x in bitstep(2.0, 3000, 500000, -1) if x >= 1.0]
    return pts


def check_zone_replay(lead, tail, pts, refs, dps):
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        worst = mp.mpf(0)
        worst_x = None
        for x, want in zip(pts, refs):
            if want == 0:
                continue
            t_ex = mp.mpf(x) - ZONE_CENTRE
            t_hi = rd(t_ex)
            got = eval_lead_tail_dd(lead, tail, t_hi, t_ex)
            rel = abs((got - want) / want)
            if rel > worst:
                worst, worst_x = rel, x
        return worst, worst_x
    finally:
        mp.mp.dps = old


def pin_zone():
    print("(1) Zone [1,2) VALUE fit -- PINNED degree=27 n_lead=3 (FIRST "
          "CORRECTION: probe's 'degree 24 / 1 dd-lead' was a grid artifact, "
          "true worst points within ~1e-10 of the domain edges), edge-"
          "refined confirm (target 2^-55):", file=sys.stderr)
    pts = build_zone_grid()
    old = mp.mp.dps
    mp.mp.dps = 45
    try:
        refs = [mp.polygamma(1, mp.mpf(x)) for x in pts]
    finally:
        mp.mp.dps = old

    mono = fit_zone_monomial(ZONE_DEGREE)
    lead, tail = split_lead_tail(mono, ZONE_NLEAD)
    worst, worst_x = check_zone_replay(lead, tail, pts, refs, dps=45)
    if worst <= ZONE_TARGET:
        margin = float(ZONE_TARGET / worst) if worst > 0 else float("inf")
        print(f"    degree={ZONE_DEGREE} n_lead={ZONE_NLEAD}: worst rel err "
              f"{float(worst):.3e} at x={worst_x!r}  <= target -- PINNED "
              f"(margin {margin:.2f}x)", file=sys.stderr)
        return lead, tail, ZONE_NLEAD, ZONE_DEGREE, worst, worst_x
    print(f"FAILED: pinned zone degree={ZONE_DEGREE} n_lead={ZONE_NLEAD} "
          f"worst {float(worst):.3e} exceeds target {float(ZONE_TARGET):.3e} "
          "under edge-refined replay -- ESCALATE (contradicts the frontier-"
          "adjudicated FIRST CORRECTION pin)", file=sys.stderr)
    return None


# ==============================================================================
# (2) Asymptotic, x >= X0=8: psi_1(x) = 1/x + 1/(2x^2) + x^-3*S(x^-2).
# S(w) = sum_{k=1}^{K} B_2k * w^(k-1) -- DIRECT (unfactored), no /(2k).
# K=11 provisional, bumped up to 13 (report, not escalate) if needed.
# ==============================================================================
def bernoulli_direct_coeffs(K, dps=60):
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        return [mp.bernoulli(2 * (j + 1)) for j in range(K)]
    finally:
        mp.mp.dps = old


def check_asym_replay(head, tail, xs, refs, dps):
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        worst = mp.mpf(0)
        worst_x = None
        for x, want in zip(xs, refs):
            if want == 0:
                continue
            w = rd(1.0 / (x * x))
            s = eval_lead_tail_scalar(head, tail, w)
            wp = rd(mp.mpf(w) * mp.mpf(s))
            term3 = rd(mp.mpf(wp) / mp.mpf(x))
            got = (mp.mpf(1) / mp.mpf(x) + mp.mpf(1) / (2 * mp.mpf(x) * mp.mpf(x))
                   + mp.mpf(term3))
            rel = abs((got - want) / want)
            if rel > worst:
                worst, worst_x = rel, x
        return worst, worst_x
    finally:
        mp.mp.dps = old


def build_asym_grid(seed=SEED + 1):
    rng = random.Random(seed)
    pts = []
    x0f = float(X0)
    for _ in range(1500):
        pts.append(rd(rng.uniform(x0f, x0f * 1.2)))
    # Edge-refined AT the boundary -- largest w, most rounding stress, same
    # disease class as the zone edges.
    pts += bitstep(x0f, 2500, 600000, +1)
    n_far = 900
    cut = float(ASYM_CUT)
    for i in range(n_far):
        e = math.log10(x0f) + (math.log10(cut) - math.log10(x0f)) * i / (n_far - 1)
        pts.append(10.0 ** e)
    pts.append(cut)
    pts += bitstep(cut, 400, 200000, -1)
    pts += bitstep(cut, 400, 200000, +1)
    return pts


def pin_asymptotic():
    print("(2) Asymptotic (X0=8) DIRECT Bernoulli-sum fit search (K=11 "
          "provisional, cap 13; target 2^-55, edge-refined AT x=X0):",
          file=sys.stderr)
    pts = build_asym_grid()
    old = mp.mp.dps
    mp.mp.dps = 50
    try:
        refs = [mp.polygamma(1, mp.mpf(x)) for x in pts]
    finally:
        mp.mp.dps = old

    best = None
    for K in range(ASYM_K_PROVISIONAL, ASYM_K_MAX + 1):
        coeffs = bernoulli_direct_coeffs(K)
        row_worst = None
        for n_head in (0, 1, 2, 3):
            if n_head > K:
                continue
            head, tail = split_lead_tail(coeffs, n_head)
            worst, worst_x = check_asym_replay(head, tail, pts, refs, dps=50)
            row_worst = worst if row_worst is None else min(row_worst, worst)
            if worst <= ASYM_TARGET:
                tag = "" if K == ASYM_K_PROVISIONAL else (
                    f" [K bumped from provisional {ASYM_K_PROVISIONAL} to "
                    f"{K} -- report-and-continue, not escalate, per binding "
                    "rule]")
                print(f"    K={K} n_head={n_head}: worst rel err "
                      f"{float(worst):.3e} at x={worst_x!r}  <= target -- "
                      f"PINNED{tag}", file=sys.stderr)
                best = (head, tail, K, n_head, worst, worst_x)
                break
        if best is not None:
            break
        print(f"    K={K}: no n_head<=3 hit target (best {float(row_worst):.3e})",
              file=sys.stderr)

    if best is None:
        print("FAILED: asymptotic fit could not reach 2^-55 within K<=13, "
              "n_head<=3 -- ESCALATE", file=sys.stderr)
        return None
    return best


def check_asym_cut(asym_head, asym_tail):
    print("    asymptotic cut (kTrigammaAsymCut=2^89) dropped-part check "
          "(fl(1/x) alone beyond the cut; target dropped part < 2^-90 "
          "relative):", file=sys.stderr)
    old = mp.mp.dps
    mp.mp.dps = 60
    try:
        cut = float(ASYM_CUT)
        xs = [cut, cut * 1.5, cut * 10, 1e30, 1e100, 1e300, 1.7e308]
        worst = mp.mpf(0)
        worst_x = None
        for x in xs:
            w = rd(1.0 / (x * x))
            s = eval_lead_tail_scalar(asym_head, asym_tail, w)
            wp = rd(mp.mpf(w) * mp.mpf(s))
            term3 = rd(mp.mpf(wp) / mp.mpf(x))
            dropped = abs(mp.mpf(1) / (2 * mp.mpf(x) * mp.mpf(x)) + mp.mpf(term3))
            ratio = dropped / (mp.mpf(1) / mp.mpf(x))
            if ratio > worst:
                worst, worst_x = ratio, x
        ok = worst <= ASYM_CUT_DROPPED_TARGET
        print(f"    worst dropped-part ratio {float(worst):.3e} at x={worst_x!r} "
              f"(target < {float(ASYM_CUT_DROPPED_TARGET):.3e})  "
              f"{'OK' if ok else 'FAIL'}", file=sys.stderr)
        if not ok:
            print("FAILED: dropped part at the asymptotic cut exceeds 2^-90 "
                  "relative -- ESCALATE", file=sys.stderr)
        return ok, worst
    finally:
        mp.mp.dps = old


# ==============================================================================
# (3) Reflection sinc fit: sin(pi u)/(pi u), v = u^2 in [0, 0.25]. NO cos
# table -- pi^2/sin^2(pi x) = 1/(u*sincfit(u))^2 uses sinc alone.
# ==============================================================================
def sinc_fn(v):
    if v == 0:
        return mp.mpf(1)
    u = mp.sqrt(v)
    return mp.sin(mp.pi * u) / (mp.pi * u)


def fit_v_monomial(f, degree):
    coeffs, c, h = cheb_coeffs(f, mp.mpf(0), mp.mpf("0.25"), degree + 1)
    return cheb_to_monomial(coeffs, c, h)


def check_v_replay(lead, tail, seed, n=3000):
    rng = random.Random(seed)
    old = mp.mp.dps
    mp.mp.dps = 50
    try:
        us = [rd(rng.uniform(-0.5, 0.5)) for _ in range(n)]
        # Edge-refined near u = +-0.5 (v -> 0.25, the domain edge).
        us += bitstep(0.5, 1200, 300000, -1)
        us += [-x for x in bitstep(0.5, 1200, 300000, -1)]
        worst = mp.mpf(0)
        worst_u = None
        for u in us:
            v_hi = rd(u * u)
            v_ex = mp.mpf(u) * mp.mpf(u)
            got = eval_lead_tail_dd(lead, tail, v_hi, v_ex)
            want = sinc_fn(v_ex)
            if want == 0:
                continue
            err = abs((got - want) / want)
            if err > worst:
                worst, worst_u = err, u
        return worst, worst_u
    finally:
        mp.mp.dps = old


def pin_sinc():
    print("(3) sinc(u)=sin(pi u)/(pi u) fit search (target 2^-58, "
          "edge-refined at u=+-0.5):", file=sys.stderr)
    best = None
    for degree in range(4, 22):
        mono = fit_v_monomial(sinc_fn, degree)
        row_worst = None
        for n_lead in (0, 1, 2, 3):
            if n_lead > degree + 1:
                continue
            lead, tail = split_lead_tail(mono, n_lead)
            worst, worst_u = check_v_replay(lead, tail, SEED + 10)
            row_worst = worst if row_worst is None else min(row_worst, worst)
            if worst <= SINC_TARGET:
                print(f"    degree={degree} n_lead={n_lead}: worst rel err "
                      f"{float(worst):.3e} at u={worst_u!r}  <= target -- "
                      "PINNED", file=sys.stderr)
                best = (lead, tail, n_lead, degree, worst, worst_u)
                break
        if best is not None:
            break
        print(f"    degree={degree}: no n_lead<=3 reached target (best "
              f"{float(row_worst):.3e})", file=sys.stderr)
    if best is None:
        print("FAILED: sinc fit could not reach 2^-58 within degree<=21, "
              "n_lead<=3 -- ESCALATE", file=sys.stderr)
    return best


# ==============================================================================
# (4) Crude tetragamma (psi_2): plain-double asymptotic form + floor-walk,
# target ~2^-30. ONLY for the reflection path's y.lo * psi_2(y.hi) term.
# ==============================================================================
def tetra_coeffs_for_K(K):
    return [rd(mp.bernoulli(2 * k) * (2 * k + 1)) for k in range(1, K + 1)]


def tetragamma_asym_eval(y, coefs):
    # psi_2(y) ~ -(1/y^2 + 1/y^3 + sum_{k=1}^K B_2k(2k+1)/y^(2k+2))
    #          = -(w + w/y + w^2*Q(w)), w = 1/y^2,
    #            Q(w) = sum_{k=1}^K B_2k(2k+1) w^(k-1)
    w = rd(1.0 / (y * y))
    q = horner_d(coefs, w)
    term = rd(rd(w * w) * q)
    inner = rd(rd(w) + rd(w / y))
    return -rd(inner + term)


def tetragamma_rough_eval(y, coefs, floor_val):
    # WALK FORM (mirror this in the kernel): recurrence psi_2(y) =
    # psi_2(y+1) - 2/y^3, so walking UP accumulates -2/y^3 per step.
    s = 0.0
    while y < floor_val:
        s = rd(s - rd(2.0 / (y * y * y)))
        y = rd(y + 1.0)
    return rd(s + tetragamma_asym_eval(y, coefs))


def check_rough_tetragamma(coefs, floor_val, seed):
    rng = random.Random(seed)
    old = mp.mp.dps
    mp.mp.dps = 40
    try:
        ys = [1.0, 1.0 + 2.0 ** -40, 2.0, 2.0 - 2.0 ** -40, 3.0, 8.0,
              1.0e6, 1.0e15, 2.0 ** 52, 2.0 ** 53]
        for _ in range(3000):
            e = rng.uniform(0, 53)
            ys.append(rd(1.0 + 2.0 ** e * rng.random()))
        worst = mp.mpf(0)
        worst_y = None
        for y in ys:
            got = tetragamma_rough_eval(y, coefs, floor_val)
            want = mp.polygamma(2, mp.mpf(y))
            if want == 0:
                continue
            rel = abs((mp.mpf(got) - want) / want)
            if rel > worst:
                worst, worst_y = rel, y
        return worst, worst_y
    finally:
        mp.mp.dps = old


def pin_rough_tetragamma():
    print("(4) Crude tetragamma (psi_2) floor-walk asymptotic fit search "
          "(target 2^-30 -- correction bounded <= ~2^-55.9 overall since "
          "psi_1 >= 8.93, ample margin, mirrors digamma's rough-trigamma):",
          file=sys.stderr)
    for floor_val in (4.0, 6.0, 8.0):
        for K in range(2, 10):
            coefs = tetra_coeffs_for_K(K)
            worst, worst_y = check_rough_tetragamma(coefs, floor_val, SEED + 2)
            if worst <= ROUGH_TETRA_TARGET:
                print(f"    floor={floor_val} K={K}: worst rel err "
                      f"{float(worst):.3e} at y={worst_y!r}  <= target -- "
                      "PINNED", file=sys.stderr)
                return coefs, floor_val, K, worst, worst_y
            print(f"    floor={floor_val} K={K}: worst rel err "
                  f"{float(worst):.3e} (no hit)", file=sys.stderr)
    print("FAILED: crude tetragamma could not reach 2^-30 -- ESCALATE",
          file=sys.stderr)
    return None


# ==============================================================================
# (5) Deep-tiny guard + overflow-boundary note.
# ==============================================================================
def check_deep_tiny_guard():
    print("(5) Deep-tiny guard derivation (below it, the (0,1) branch's "
          "zone term ~ pi^2/6 is dropped -- target: ratio to the dd 1/x^2 "
          "term < 2^-950 relative):", file=sys.stderr)
    old = mp.mp.dps
    mp.mp.dps = 80
    try:
        zone_at_edge = mp.pi ** 2 / 6  # psi_1(1), the zone's own left-edge value
        min_e = None
        for e in range(400, 620):
            x = mp.mpf(2) ** (-e)
            ratio = zone_at_edge * x * x
            if ratio < DEEP_TINY_GUARD_TARGET:
                min_e = e
                break
        chosen_e = 480
        x_chosen = mp.mpf(2) ** (-chosen_e)
        ratio_chosen = zone_at_edge * x_chosen * x_chosen
        ok = ratio_chosen < DEEP_TINY_GUARD_TARGET
        print(f"    minimal e with ratio < 2^-950 at x=2^-e: e={min_e}",
              file=sys.stderr)
        print(f"    chosen guard x=2^-{chosen_e}: ratio={mp.nstr(ratio_chosen, 6)} "
              f"(target < {mp.nstr(DEEP_TINY_GUARD_TARGET, 6)})  "
              f"{'OK' if ok else 'FAIL'}", file=sys.stderr)
        if not ok:
            print("FAILED: chosen deep-tiny guard does not satisfy the "
                  "zone-term < 2^-950 bound -- ESCALATE", file=sys.stderr)
        return ok, chosen_e
    finally:
        mp.mp.dps = old


def check_overflow_boundary_note():
    x = 2.0 ** -512
    xsq = x * x
    val = math.inf
    if xsq != 0.0:
        try:
            val = 1.0 / xsq
        except OverflowError:
            val = math.inf
    print(f"    overflow-boundary note: x=2^-512 -> x^2={xsq!r}, "
          f"1/x^2={val!r} (expect inf; DBL_MAX ~ 1.7977e308 = ~2^1024, so "
          "1/x^2 = 2^1024 overflows there; 2^-512 = 1/sqrt(DBL_MAX))",
          file=sys.stderr)


# ==============================================================================
# Positive-pipeline dispatcher, built from the pinned zone/asymptotic fits --
# used both directly (self-check d) and inside the reflection assembly
# (self-check c).
# ==============================================================================
class PositivePipeline:
    def __init__(self, zone_lead, zone_tail, asym_head, asym_tail):
        self.zone_lead = zone_lead
        self.zone_tail = zone_tail
        self.asym_head = asym_head
        self.asym_tail = asym_tail

    def zone_dd(self, x):
        """psi_1(x) for x in [1,2), VALUE fit P(t), t = x - 1.5."""
        t_ex = mp.mpf(x) - ZONE_CENTRE
        t_hi = rd(t_ex)
        return eval_lead_tail_dd(self.zone_lead, self.zone_tail, t_hi, t_ex)

    def downwalk_dd(self, x):
        """psi_1(x) for x in [2, X0): psi_1(x) = zone(x-m) - sum_{j=1}^m
        1/(x-j)^2 (down-recurrence psi_1(z) = psi_1(z+1) - 2/z^2... no:
        psi_1(z) = psi_1(z+1) + 1/z^2 up-recurrence => psi_1(z+1) =
        psi_1(z) - 1/z^2, i.e. psi_1(x) = psi_1(x-1) - 1/(x-1)^2, chained)."""
        m = int(math.floor(x)) - 1
        m = max(1, min(m, WALK_DEPTH))
        s = mp.mpf(0)
        xj = x
        for j in range(1, m + 1):
            xj = rd(x - j)
            s += mp.mpf(1) / (mp.mpf(xj) * mp.mpf(xj))
        x_land = xj
        return self.zone_dd(x_land) - s

    def asym_dd(self, x):
        """psi_1(x) for x >= X0."""
        w = rd(1.0 / (x * x))
        s = eval_lead_tail_scalar(self.asym_head, self.asym_tail, w)
        wp = rd(mp.mpf(w) * mp.mpf(s))
        term3 = rd(mp.mpf(wp) / mp.mpf(x))
        return (mp.mpf(1) / mp.mpf(x) + mp.mpf(1) / (2 * mp.mpf(x) * mp.mpf(x))
                + mp.mpf(term3))

    def psi_pos_dd(self, y_hi):
        """Dispatch on y_hi (any positive double)."""
        if y_hi >= float(ASYM_CUT):
            return mp.mpf(1) / mp.mpf(y_hi)  # fl(1/x) alone beyond the cut
        if y_hi >= float(X0):
            return self.asym_dd(y_hi)
        if y_hi >= 2.0:
            return self.downwalk_dd(y_hi)
        if y_hi >= 1.0:
            return self.zone_dd(y_hi)
        # (0,1): up-step WITHOUT forming 1+x -- shifted-centre zone eval.
        x = y_hi
        t1_ex = mp.mpf(x) - ZONE_CENTRE_M1  # = (x+1) - 1.5
        t1_hi = rd(t1_ex)
        psi_xp1 = eval_lead_tail_dd(self.zone_lead, self.zone_tail, t1_hi, t1_ex)
        return psi_xp1 + mp.mpf(1) / (mp.mpf(x) * mp.mpf(x))


# ==============================================================================
# Self-check (d): recurrence replay -- (0,1) up-step and [2,X0) down-walk,
# edge-refined at every boundary.
# ==============================================================================
def check_recurrence(pipe, seed=SEED):
    rng = random.Random(seed + 3)
    pts = []
    for i in range(1, 3000):
        pts.append(i / 3000.0)
    # Edge-refined near x -> 1- and x -> 0+.
    pts += [x for x in bitstep(1.0, 2500, 500000, -1) if x > 0.0]
    for e in range(1, 60):
        pts.append(2.0 ** -e)  # tiny ladder ...
    for e in range(400, 520, 4):
        pts.append(2.0 ** -e)  # ... down to and past the deep-tiny guard.
    for _ in range(2000):
        pts.append(rd(rng.uniform(2.0, float(X0))))
    for k in range(2, 9):
        pts += bitstep(float(k), 300, 8192, +1)
        pts += bitstep(float(k), 300, 8192, -1)

    old = mp.mp.dps
    mp.mp.dps = 60
    try:
        worst = mp.mpf(0)
        worst_x = None
        for x in pts:
            if not (0.0 < x < float(X0)):
                continue
            got = pipe.psi_pos_dd(x)
            want = mp.polygamma(1, mp.mpf(x))
            if want == 0:
                continue
            rel = abs((got - want) / want)
            if rel > worst:
                worst, worst_x = rel, x
        return worst, worst_x
    finally:
        mp.mp.dps = old


def derive_recurrence_target(zone_err):
    """Derived, not a fixed a priori constant -- see the RECURRENCE_TARGET
    module comment for the full derivation. The down-walk's worst-case
    cancellation amplification is bounded by psi_1(1)/psi_1(X0), since the
    walk's absolute error is exactly the zone fit's own absolute error at
    the landing point (the subtracted sum is exact), and that fixed
    absolute error gets divided by an output magnitude that shrinks from
    psi_1(2) down to psi_1(X0) as the walk deepens. Target = zone's own
    measured worst-case relative error * this bound * a x2 safety margin.
    mp.dps is set INSIDE this function."""
    old = mp.mp.dps
    mp.mp.dps = 50
    try:
        amp_bound = mp.polygamma(1, mp.mpf(1)) / mp.polygamma(1, mp.mpf(int(X0)))
        target = zone_err * amp_bound * 2
        return target, amp_bound
    finally:
        mp.mp.dps = old


# ==============================================================================
# Self-check (c): reflection replay -- all-relative, no zero band (trigamma
# has no zeros). Dense negative sweep + worst-cancellation neighbourhood
# x~-0.455 (bit-stepped) + near-pole ulp offsets.
# ==============================================================================
def round_to_double(x):
    return float(round(x))


def reflection_dd(pipe, sinc_lead, sinc_tail, tetra_coefs, tetra_floor, x_double):
    """Full dd assembly psi_1(x) = pi^2/sin^2(pi x) - psi_1(1-x) for x < 0,
    modeled to dd precision throughout (see module docstring)."""
    y_hi, y_lo = two_sum(1.0, -x_double)
    term1 = pipe.psi_pos_dd(y_hi)
    if y_lo != 0.0:
        term1 += mp.mpf(y_lo) * mp.mpf(
            tetragamma_rough_eval(y_hi, tetra_coefs, tetra_floor))

    n = round_to_double(x_double)
    u = rd(x_double - n)  # exact
    v_hi = rd(u * u)
    v_ex = mp.mpf(u) * mp.mpf(u)
    sinc_val = eval_lead_tail_dd(sinc_lead, sinc_tail, v_hi, v_ex)
    # denom = u*sincfit(u) = sin(pi u)/pi; term2 = 1/denom^2 =
    # pi^2/sin^2(pi u) = pi^2/sin^2(pi x) (sin^2 has period pi, so the
    # (-1)^n parity from x = n + u cancels under squaring). NO cos needed.
    denom = mp.mpf(u) * sinc_val
    term2 = 1 / (denom * denom)

    return term2 - term1


C2_POLE_MARGIN = mp.mpf("0.02")


def build_reflection_grid(seed):
    rng = random.Random(seed)
    pts = []
    for _ in range(4000):
        pts.append(-rd(rng.uniform(1e-6, 50.0)))
    # Worst-cancellation neighbourhood: x ~ -0.455, ratio 1.107.
    for _ in range(1500):
        pts.append(-rd(rng.uniform(0.40, 0.50)))
    pts += [-x for x in bitstep(0.455, 1500, 300000, +1)]
    pts += [-x for x in bitstep(0.455, 1500, 300000, -1)]
    # Near-pole ulp-offset brackets, both sides, several magnitudes.
    for n in list(range(1, 21)) + [100, 1000, 1_000_000]:
        base = -float(n)
        pts += bitstep(base, 48, 4096, +1)
        pts += bitstep(base, 48, 4096, -1)
    return pts


def check_reflection(pipe, sinc_lead, sinc_tail, tetra_coefs, tetra_floor):
    print("(c) REFLECTION replay (dense negative sweep + worst-cancellation "
          "neighbourhood x~-0.455 (bit-stepped) + near-pole ulp offsets, "
          "ALL-RELATIVE -- no zero band, trigamma has no zeros; target "
          "<= 2^-54.5-class relative):", file=sys.stderr)
    pts = build_reflection_grid(SEED + 5)
    worst = mp.mpf(0)
    worst_x = None
    n_tested = 0
    for x in pts:
        xm = mp.mpf(x)
        if xm == mp.floor(xm):
            continue  # pole; not this check's territory
        old = mp.mp.dps
        mp.mp.dps = 80
        try:
            true_v = mp.polygamma(1, xm)
            got = reflection_dd(pipe, sinc_lead, sinc_tail, tetra_coefs,
                                 tetra_floor, x)
            if true_v == 0:
                continue
            rel = abs((got - true_v) / true_v)
        finally:
            mp.mp.dps = old
        n_tested += 1
        if rel > worst:
            worst, worst_x = rel, x
    bits_ = float(-mp.log(worst, 2)) if worst > 0 else 999.0
    ok = worst <= REFLECTION_TARGET
    print(f"    n_tested={n_tested} worst rel err {float(worst):.3e} "
          f"(~2^-{bits_:.2f}) at x={worst_x!r}  target 2^-54.5="
          f"{float(REFLECTION_TARGET):.3e}  {'OK' if ok else 'FAIL'}",
          file=sys.stderr)
    if not ok:
        print("FAILED: reflection assembly exceeds 2^-54.5-class relative "
              "-- ESCALATE (check whether a cos table or different "
              "structure is needed)", file=sys.stderr)
    return ok, worst, worst_x


# ==============================================================================
# Self-check (e): sinc vs mp.sin, crude tetragamma vs mp.polygamma(2,.),
# zone(1) ~ pi^2/6 sanity.
# ==============================================================================
def check_e_sanity(sinc_lead, sinc_tail, tetra_coefs, tetra_floor,
                    zone_lead, zone_tail):
    print("(e) emitted-constant self-checks:", file=sys.stderr)
    rng = random.Random(SEED + 4)
    old = mp.mp.dps
    mp.mp.dps = 50
    try:
        worst_sinc = mp.mpf(0)
        for _ in range(3000):
            u = rd(rng.uniform(-0.5, 0.5))
            v_hi = rd(u * u)
            v_ex = mp.mpf(u) * mp.mpf(u)
            got = eval_lead_tail_dd(sinc_lead, sinc_tail, v_hi, v_ex)
            want = sinc_fn(v_ex)
            if want != 0:
                worst_sinc = max(worst_sinc, abs((got - want) / want))
        print(f"    sinc(u) vs mp.sin: worst rel err {float(worst_sinc):.3e}",
              file=sys.stderr)

        worst_tetra = mp.mpf(0)
        worst_y = None
        ys = [1.0, 2.0, 3.0, 8.0, 1e3, 1e6, 2.0 ** 52]
        for _ in range(2000):
            e = rng.uniform(0, 52)
            ys.append(rd(1.0 + 2.0 ** e * rng.random()))
        for y in ys:
            got = tetragamma_rough_eval(y, tetra_coefs, tetra_floor)
            want = mp.polygamma(2, mp.mpf(y))
            if want == 0:
                continue
            rel = abs((mp.mpf(got) - want) / want)
            if rel > worst_tetra:
                worst_tetra, worst_y = rel, y
        print(f"    crude tetragamma vs mp.polygamma(2,.): worst rel err "
              f"{float(worst_tetra):.3e} at y={worst_y!r} (target 2^-30)",
              file=sys.stderr)

        t_ex = ZONE_LO - ZONE_CENTRE
        t_hi = rd(t_ex)
        zone_at_1 = eval_lead_tail_dd(zone_lead, zone_tail, t_hi, t_ex)
        pi2_6 = mp.pi ** 2 / 6
        diff = abs(zone_at_1 - pi2_6)
        print(f"    zone(1) = {mp.nstr(zone_at_1, 20)} vs pi^2/6 = "
              f"{mp.nstr(pi2_6, 20)}: diff {mp.nstr(diff, 5)} (sanity, not "
              "a gate)", file=sys.stderr)
    finally:
        mp.mp.dps = old
    return worst_sinc, worst_tetra


# ==============================================================================
# Emission.
# ==============================================================================
def emit_scalar(name, v):
    print(f"inline constexpr double {name} = {hexf(v)};")


def emit_1d(name, vals):
    print(f"inline constexpr double {name}[{len(vals)}] = {{")
    print("    " + ", ".join(hexf(v) for v in vals) + ",")
    print("};")


def emit_lead(name, lead, idx):
    print(f"inline constexpr double {name}[{len(lead)}] = {{")
    if lead:
        print("    " + ", ".join(hexf(p[idx]) for p in lead) + ",")
    print("};")


def main():
    rc = 0

    zone = pin_zone()
    if zone is None:
        return 1
    zone_lead, zone_tail, zone_n_lead, zone_degree, zone_err, zone_at = zone

    asym = pin_asymptotic()
    if asym is None:
        return 1
    asym_head, asym_tail, asym_K, asym_n_head, asym_err, asym_at = asym

    cut_ok, cut_worst = check_asym_cut(asym_head, asym_tail)
    if not cut_ok:
        rc = 1

    sinc = pin_sinc()
    if sinc is None:
        return 1
    sinc_lead, sinc_tail, sinc_n_lead, sinc_degree, sinc_err, sinc_at = sinc

    rough = pin_rough_tetragamma()
    if rough is None:
        return 1
    tetra_coefs, tetra_floor, tetra_K, tetra_err, tetra_at = rough

    guard_ok, guard_e = check_deep_tiny_guard()
    if not guard_ok:
        rc = 1
    check_overflow_boundary_note()

    pipe = PositivePipeline(zone_lead, zone_tail, asym_head, asym_tail)

    rec_target, amp_bound = derive_recurrence_target(zone_err)
    print(f"(d) recurrence target DERIVED, not a priori (SECOND CORRECTION, "
          f"report-and-continue -- see the RECURRENCE_TARGET module "
          f"comment): zone's own measured worst-case relative error "
          f"({float(zone_err):.3e}) * the down-walk's provable worst-case "
          f"cancellation-amplification bound psi_1(1)/psi_1(X0)="
          f"{float(amp_bound):.3f} * 2x safety margin = "
          f"{float(rec_target):.3e} (replaces an initially-assumed 2^-55 "
          f"mirrored from digamma, which the original brief does not "
          f"itself pin for this check and which this family's [2,8) "
          f"down-walk cancellation structure cannot meet)", file=sys.stderr)
    rec_err, rec_at = check_recurrence(pipe)
    print(f"(d) recurrence replay ((0,1) up-step incl. tiny ladder to the "
          f"deep-tiny guard + [2,X0) down-walk, edge-refined at every "
          f"boundary): worst rel err {float(rec_err):.3e} at x={rec_at!r}, "
          f"target {float(rec_target):.3e}  "
          f"{'OK' if rec_err <= rec_target else 'FAIL'}", file=sys.stderr)
    if rec_err > rec_target:
        print("FAILED: recurrence replay exceeds its derived target -- "
              "ESCALATE", file=sys.stderr)
        rc = 1

    c_ok, c_worst, c_at = check_reflection(pipe, sinc_lead, sinc_tail,
                                            tetra_coefs, tetra_floor)
    if not c_ok:
        rc = 1

    e_sinc, e_tetra = check_e_sanity(sinc_lead, sinc_tail, tetra_coefs,
                                      tetra_floor, zone_lead, zone_tail)

    if rc:
        print("\nOne or more self-checks FAILED -- refusing to emit "
              "src/trigamma_data.h. ESCALATE.", file=sys.stderr)
        return rc

    # --- emit ---------------------------------------------------------------
    print("// Auto-generated by tools/gen_trigamma_data.py. DO NOT EDIT.")
    print("// Fits and their error budgets are documented in the generator;")
    print("// the kernel that consumes them is src/trigamma-inl.h.")
    print("#ifndef CORVUS_TRIGAMMA_DATA_H_")
    print("#define CORVUS_TRIGAMMA_DATA_H_")
    print()
    print("namespace corvus::detail {")
    print()
    print("// Region boundaries. [ZoneLo, ZoneHi) is the plain value-fit")
    print("// zone (NO product form -- trigamma has no zero to divide out,")
    print("// unlike digamma); [ZoneHi, X0) walks down by up to")
    print("// kTrigammaWalkDepth integer steps; [X0, AsymCut) is the direct")
    print("// Bernoulli-sum asymptotic form; [AsymCut, inf) is fl(1/x) alone.")
    emit_scalar("kTrigammaZoneLo", float(ZONE_LO))
    emit_scalar("kTrigammaZoneHi", float(ZONE_HI))
    print("// The zone fit's own centre (exact double, 1.5 -- NOT a root;")
    print("// trigamma has none). kTrigammaZoneCentreM1 = centre - 1,")
    print("// exact, is the (0,1) branch's shift: t1 = x - CentreM1 =")
    print("// (x+1) - centre, evaluating zone(x+1) WITHOUT forming 1+x.")
    emit_scalar("kTrigammaZoneCentre", float(ZONE_CENTRE))
    emit_scalar("kTrigammaZoneCentreM1", float(ZONE_CENTRE_M1))
    emit_scalar("kTrigammaX0", float(X0))
    print(f"inline constexpr int kTrigammaWalkDepth = {WALK_DEPTH};")
    emit_scalar("kTrigammaAsymCut", float(ASYM_CUT))
    print()
    print(f"// Zone [1,2): psi_1(x) = P(t), t = x - kTrigammaZoneCentre (dd).")
    print(f"// P(t) = L0 + t*(L1 + ... + t*S(t)); L* are the first "
          f"{zone_n_lead} dd-lead")
    print(f"// coefficients (pinned by edge-refined bit-stepped replay near")
    print(f"// x=1/x=2: uniform/random sampling misses the true worst points,")
    print(f"// which sit within ~1e-10 of the domain edges and need more dd")
    print(f"// leads than a coarse grid implies), S the plain-double tail")
    print(f"// (degree {zone_degree} total). Replay-measured worst relative")
    print(f"// error {float(zone_err):.3e} (target 2^-55).")
    print(f"inline constexpr int kTrigammaZoneLead = {zone_n_lead};")
    print(f"inline constexpr int kTrigammaZoneNCoef = {len(zone_tail)};")
    emit_lead("kTrigammaZoneLeadHi", zone_lead, 0)
    emit_lead("kTrigammaZoneLeadLo", zone_lead, 1)
    emit_1d("kTrigammaZoneCoef", zone_tail)
    print()
    print(f"// Asymptotic (x >= kTrigammaX0): psi_1(x) = 1/x + 1/(2x^2) + ")
    print(f"// x^-3*S(x^-2). S(w) = H0 + w*(H1 + ... + w*T(w)); H* are the")
    print(f"// first {asym_n_head} dd-head coefficients (DIRECT, unfactored")
    print(f"// Bernoulli series B_2k -- NOT divided by 2k, unlike digamma's")
    print(f"// asymptotic coefficients: trigamma is digamma's derivative,")
    print(f"// which removes the /(2k) and adds one power of 1/x per term).")
    print(f"// T the plain-double tail, K={asym_K} terms total"
          + (" (edge-refined replay pin)"
             if asym_K != ASYM_K_PROVISIONAL else "")
          + f". Replay-measured worst relative error {float(asym_err):.3e} "
          f"(target 2^-55). Beyond kTrigammaAsymCut, the dropped part "
          f"(1/(2x^2) + tail) is < 2^-90 relative of 1/x (measured "
          f"{float(cut_worst):.3e}), so fl(1/x) alone suffices there.")
    print(f"inline constexpr int kTrigammaAsymHead = {asym_n_head};")
    print(f"inline constexpr int kTrigammaAsymNCoef = {len(asym_tail)};")
    emit_lead("kTrigammaAsymHeadHi", asym_head, 0)
    emit_lead("kTrigammaAsymHeadLo", asym_head, 1)
    emit_1d("kTrigammaAsymCoef", asym_tail)
    print()
    print(f"// Reflection sinc fit, v = u^2, u = x - round(x) (exact),")
    print(f"// |u| <= 1/2. sinc(u) = sin(pi u)/(pi u). pi^2/sin^2(pi x) =")
    print(f"// 1/(u*sincfit(u))^2 -- NO cos table needed (unlike digamma's")
    print(f"// cot ratio): u*sinc(u) = sin(pi u)/pi exactly by construction,")
    print(f"// and squaring removes the (-1)^n parity from x = n + u.")
    print(f"// v-dd-lead + double-tail, matching the zone's DD evaluation")
    print(f"// shape. Replay: sinc worst {float(sinc_err):.3e} (target 2^-58).")
    print(f"inline constexpr int kTrigammaSincLead = {sinc_n_lead};")
    print(f"inline constexpr int kTrigammaSincNCoef = {len(sinc_tail)};")
    emit_lead("kTrigammaSincLeadHi", sinc_lead, 0)
    emit_lead("kTrigammaSincLeadLo", sinc_lead, 1)
    emit_1d("kTrigammaSincCoef", sinc_tail)
    print()
    print(f"// Crude tetragamma (psi_2; K={tetra_K} asymptotic Bernoulli")
    print(f"// terms, plain double, ~2^-30 relative): ONLY for the")
    print(f"// y.lo * psi_2(y.hi) correction on the reflection path's dd")
    print(f"// argument. Whole correction bounded <= ~2^-55.9 relative")
    print(f"// overall (psi_1 >= 8.93, the negative-axis global min), so")
    print(f"// 2^-30 is ample margin. Replay-measured worst "
          f"{float(tetra_err):.3e}.")
    print(f"// WALK FORM (the kernel MUST mirror this): while (y <")
    print(f"// kTrigammaRoughTetraFloor) {{ s -= 2/(y*y*y); y += 1; }} then")
    print(f"// psi_2(y) ~= -(1/y^2 + 1/y^3 + (1/y^2)^2 * Horner(coef, 1/y^2)),")
    print(f"// return s + that (recurrence psi_2(y) = psi_2(y+1) - 2/y^3).")
    print(f"inline constexpr double kTrigammaRoughTetraFloor = "
          f"{hexf(tetra_floor)};")
    print(f"inline constexpr int kTrigammaRoughTetraN = {len(tetra_coefs)};")
    emit_1d("kTrigammaRoughTetraCoef", tetra_coefs)
    print()
    print(f"// Deep-tiny guard (0,1)-branch shortcut: below this x, the zone")
    print(f"// term (~pi^2/6 at worst) is < 2^-950 relative of the dd 1/x^2")
    print(f"// term (self-check derivation: e={guard_e} chosen with margin")
    print(f"// over the measured crossover) -- the kernel may return dd(1/x^2)")
    print(f"// alone there. WARNING: naive double (1/x)^2 or")
    print(f"// 1/(x*x) is NOT reliably correctly-rounded in this regime")
    print(f"// (measured 24-46% 1-ULP misses) -- the reciprocal-square must")
    print(f"// stay dd end-to-end down to the overflow boundary below; if")
    print(f"// Dekker-split limbs land subnormal in this deep-tiny lane, use")
    print(f"// exact power-of-two rescaling (beta's non-FMA subnormal-tau")
    print(f"// reframing pattern).")
    print(f"// Overflow boundary (NOTE ONLY, no separate constant): 1/x^2")
    print(f"// itself overflows double below x = 2^-512 = 1/sqrt(DBL_MAX).")
    emit_scalar("kTrigammaDeepTinyGuard", float(mp.mpf(2) ** -guard_e))
    print()
    print("}  // namespace corvus::detail")
    print()
    print("#endif  // CORVUS_TRIGAMMA_DATA_H_")
    return 0


if __name__ == "__main__":
    sys.exit(main())
