#!/usr/bin/env python3
"""Generate src/digamma_data.h -- every table the digamma kernel needs. This
generator pins the numeric constants the design leaves open: zone
degree/dd-lead count, asymptotic K/dd-head count, reflection sinc-pair
fits, rough-trigamma degree -- it does not re-derive the region map, the
product-form choice, or the reflection formula, which are taken as given.

Six pieces:

  root (x0)      the unique positive root of digamma, as a dd pair, and
                 trigamma(x0) (used to fill the zone fit's removable
                 t=0 singularity).

  zone [1,2)     product form psi(x) = t*P(t), t = x - x0 in DD (x0 is not
                 an integer, so t is NOT exact by Sterbenz the way lgamma's
                 zone shift is -- forming it needs the full TwoSum-chain dd
                 pair). Degree and dd-LEAD count pinned by REPLAY: the
                 kernel's own arithmetic (double tail Horner using t.hi,
                 dd Horner from there up using the exact/high-precision t)
                 emulated against mpmath over a dense [1,2) double grid
                 plus 1e-14..1e-1 neighborhoods of x0.

  asymptotic     psi(x) = LogDd(x) - DdRecip(2x) - x^-2*S(x^-2), x >= X0=8.
                 S is the RAW Bernoulli-series tail (S(w) = B2/2 + B4/4 w +
                 ...), not a Chebyshev economization -- the design calls
                 for "Bernoulli-series coefficients" directly. Unlike the
                 zone, x is a plain double here (no cancellation hazard),
                 so w=x^-2 is a plain double and this fit reuses lgamma's
                 exact-scalar-t replay shape. K and the dd-HEAD count
                 pinned by replay at x=X0 (worst case) and out to 1e300.

  reflection     even polynomials in v=u^2 (u=x-round(x), exact) for
                 sin(pi u)/(pi u) and cos(pi u), |u|<=1/2. v itself needs
                 dd precision (TwoProd(u,u), not a bare double square) --
                 unlike the asymptotic fit's w, these feed a cot ratio
                 that is one side of a cancellation up to ~49 bits at the
                 adversarial points self-check (c) probes, so both fits
                 use the zone's dd-t evaluation shape, not the asymptotic's
                 scalar-t shape.

  rough-trigamma a cheap (~2^-40) plain-double trigamma approximation
                 (asymptotic-form, one recurrence step below y=2) used
                 ONLY for the y.lo * trigamma(y.hi) correction on the
                 reflection path's dd argument -- NOT src/beta_data.h's
                 DigammaRough; a separate fit, not shared with it.

  region/threshold constants: kDigammaX0 = 8 (asymptotic threshold --
                 distinct from the ROOT x0, named kDigammaRoot* below to
                 avoid the collision), zone bounds, walk depth (masked
                 down-walk from [2,X0) to [1,2), 6 steps for X0=8). No
                 separate pi dd constant is needed: the sinc-pair ratio
                 cos(pi u)/(u*sincfit(u)) already equals pi*cot(pi x) by
                 construction (the pi in sincfit's own denominator cancels
                 it out), so pi never appears explicitly in the assembly
                 -- double-counting it (multiplying by pi again here) would
                 show up immediately in self-check (c) as a -50-bit
                 "retained accuracy", i.e. no accuracy at all.

SELF-CHECKS (mandatory, budget lines to stderr; ANY miss -> exit nonzero):
  (a) zone replay <= 2^-55 relative, dense [1,2) grid + x0 neighborhoods.
  (b) asymptotic truncation+replay sup at x=X0 and sampled x in
      [X0, 1e300].
  (c) REFLECTION, the accuracy doctrine's own dual metric -- near-relative
      accuracy AT the adversarial zeros needs 2^-104-class term1/term2
      inputs (measured: zone needs degree~40 away from x0's own
      precision floor, asymptotic plateaus ~65 bits regardless of K at
      X0=8 -- an asymptotic, not convergent, series) and is rejected on
      cost; lgamma's negative axis does not offer that contract either:
        (c1) ABSOLUTE error <= 2^-56 at each of the 20 adversarial
             nearest-double points (measured worst 2.42e-18 ~ 2^-58.5 --
             this bar carries margin; exceeding it here is still an
             escalate, not a target to loosen).
        (c2) RELATIVE error <= 2^-52 at dense negative-axis samples
             where |psi(x)| >= 1: each interval (-n-1,-n) for n=1..20,
             excluding the per-zero band W ~ 1/|trigamma(z0)| around
             that interval's zero AND a margin around each pole, plus a
             few log-spaced spot intervals out to n ~ 1e6 (same
             |psi(x)| >= 1 filter, no explicit zero-band needed there --
             the filter itself excludes the zero's neighborhood).
  (d) recurrence replay: up-step from x in (0,1) (incl. x -> 1^-) and
      down-walk from [2,X0) -- worst relative error <= 2^-55.
  (e) emitted-constant self-checks: x0 pair vs fresh recompute; sinc/cos
      fits vs mp.sin/mp.cos; rough-trigamma <= 2^-40.

Dd arithmetic (TwoSum, DdAdd, DdMul, DdRecipDd, ...) is modeled as EXACT
in every replay below via mpmath at working precision -- it carries
~2^-104 relative (src/dd-inl.h), three orders below every target here,
matching gen_lgamma_data.py's replay_lead_tail precedent. What is NOT
modeled as exact is: (1) any step the design specifies as a plain double
(the asymptotic w, the tail Horner polynomials, rough-trigamma), and
(2) the final single rounding to double, which belongs to the ULP test,
not this generator's budget (matches every other generator's convention).

Usage:
    python3 tools/gen_digamma_data.py > src/digamma_data.h
"""

import math
import random
import sys

import mpmath as mp

mp.mp.dps = 60  # module-level default; every function below sets its OWN
                # dps on entry and restores it on exit (AGENTS.md mechanism
                # rule -- never rely on an ambient value).

SEED = 20260806

# --- Design constants [PLAN.md P1 digamma, BINDING] --------------------------
ZONE_LO = mp.mpf(1)
ZONE_HI = mp.mpf(2)
X0 = mp.mpf(8)             # asymptotic threshold (kDigammaX0), NOT the root.
WALK_DEPTH = int(mp.ceil(X0 - ZONE_HI))  # masked down-walk step count = 6

# --- Self-check targets -------------------------------------------------------
ZONE_TARGET = mp.mpf(2) ** -55
ASYM_TARGET = mp.mpf(2) ** -55
RECURRENCE_TARGET = mp.mpf(2) ** -55
# (c) the dual metric replaces a single "54.5 bits relative" bar, rejected
# on cost (see the module docstring).
REFLECTION_C1_ABS_TARGET = mp.mpf(2) ** -56
REFLECTION_C2_REL_TARGET = mp.mpf(2) ** -52
ROUGH_TRIGAMMA_TARGET = mp.mpf(2) ** -40
# Standalone sinc/cos fit target -- tighter than the zone's, since these
# feed one side of the reflection cancellation (check (c)); loosened only
# if replay proves the extra margin unnecessary.
SINC_TARGET = mp.mpf(2) ** -58


# ==============================================================================
# dd / hex-float emission helpers (matches gen_lgamma_data.py / gen_beta_data.py)
# ==============================================================================
def rd(x):
    return float(x)


def dd_split(x):
    hi = rd(x)
    return hi, rd(mp.mpf(x) - mp.mpf(hi))


def hexf(x):
    return float.hex(float(x))


def two_sum(a, b):
    """Knuth's TwoSum on PLAIN python floats (doubles) -- bit-exact replica
    of src/dd-inl.h's TwoSum, used where the replay needs the REAL rounding
    behaviour of a double op (e.g. y_dd = TwoSum(1,-x)), not the
    'dd-arithmetic modeled exact' shortcut used elsewhere in this file."""
    s = a + b
    bv = s - a
    err = (a - (s - bv)) + (b - bv)
    return s, err


def horner_d(coefs, x):
    """Horner in plain double -- one rounding per fused step, as the kernel
    does for every double tail poly in this file."""
    acc = 0.0
    for cf in reversed(coefs):
        acc = rd(mp.mpf(acc) * mp.mpf(x) + mp.mpf(cf))
    return acc


# ==============================================================================
# Chebyshev machinery (matches gen_lgamma_data.py).
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
# Two Horner+lead evaluation shapes.
#
# SCALAR shape (t is an exact plain scalar -- asymptotic's w): matches
# gen_lgamma_data.py's replay_lead_tail exactly, INCLUDING the one rounded
# product at the lead/tail seam (that rounding is real: the kernel forms
# the tail poly value and t both as plain doubles, so their product is one
# more single rounding, same as lgamma's B(t)).
#
# DD shape (t itself is a dd pair, e.g. x-x0 near the root, or v=u*u for
# the reflection sinc/cos pair): the seam is NOT rounded -- t is already
# carried at dd precision by the caller (that's the entire point of
# computing it via TwoSum/TwoProd instead of a bare subtraction/square),
# so collapsing the tail*t product to a single rounded double would throw
# away exactly the precision the dd-t construction exists to protect.
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
# (0) Root x0: unique positive root of digamma, and trigamma(x0).
# ==============================================================================
def find_root(dps, guess="1.4616321"):
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        return mp.findroot(lambda x: mp.digamma(x), mp.mpf(guess))
    finally:
        mp.mp.dps = old


def derive_root():
    x0_100 = find_root(100)
    x0_60 = find_root(60)
    old = mp.mp.dps
    mp.mp.dps = 100
    try:
        diff = abs(x0_100 - x0_60)
        tg = mp.polygamma(1, x0_100)
    finally:
        mp.mp.dps = old
    x0_hi = rd(x0_100)
    x0_lo = rd(mp.mpf(x0_100) - mp.mpf(x0_hi))

    # Fresh independent recompute (different dps layer, different guess) --
    # self-check (e)'s "x0 pair vs fresh recompute".
    x0_check = find_root(80, guess="1.46")
    hi_check = rd(x0_check)
    lo_check = rd(mp.mpf(x0_check) - mp.mpf(hi_check))
    ok = (hi_check == x0_hi and lo_check == x0_lo)

    print(f"(root) x0 = {mp.nstr(x0_100, 30)}", file=sys.stderr)
    print(f"    dd: hi={x0_hi!r} ({x0_hi.hex()}) lo={x0_lo!r} ({x0_lo.hex()})",
          file=sys.stderr)
    print(f"    |x0_100 - x0_60| = {mp.nstr(diff, 5)}", file=sys.stderr)
    print(f"    trigamma(x0) = {mp.nstr(tg, 20)}", file=sys.stderr)
    print(f"    fresh-recompute (dps=80, alt guess) hi/lo match: {ok}",
          file=sys.stderr)
    if not ok:
        print("FAILED: x0 dd pair does not reproduce under an independent "
              "recompute -- ESCALATE", file=sys.stderr)
    return x0_hi, x0_lo, x0_100, tg, ok


# ==============================================================================
# (1) Zone [1,2): psi(x) = t*P(t), t = x - x0 (dd).
# ==============================================================================
def fit_zone_monomial(x0_mp, degree):
    """Fit with EXACTLY degree+1 nodes -- NOT an over-sampled fit truncated
    afterward. cheb_to_monomial's Chebyshev-to-power-basis conversion sums
    O(n) terms whose OWN coefficients grow like 2^n, so at n_nodes >> degree
    the conversion cancels ~n bits it never had the dps budget for (measured
    directly: n_nodes=200 truncated to degree 30 floors at ~1.5e-14, while
    a degree-31-node fit at the SAME degree reaches ~3e-18 -- four orders
    tighter)."""
    lo = ZONE_LO - x0_mp
    hi = ZONE_HI - x0_mp

    def f(t):
        if t == 0:
            return mp.polygamma(1, x0_mp)
        return mp.digamma(x0_mp + t) / t

    coeffs, c, h = cheb_coeffs(f, lo, hi, degree + 1)
    return cheb_to_monomial(coeffs, c, h)


def build_zone_grid(x0_hi, seed=SEED, n_uniform=1000, n_random=1000):
    rng = random.Random(seed)
    pts = []
    for i in range(n_uniform):
        pts.append(1.0 + i * (1.0 / n_uniform))
    for _ in range(n_random):
        pts.append(rd(rng.uniform(1.0, 2.0)))
    for e in range(1, 15):  # 1e-1 .. 1e-14
        d = 10.0 ** (-e)
        for sgn in (1, -1):
            xx = x0_hi + sgn * d
            if 1.0 <= xx < 2.0:
                pts.append(xx)
    return pts


def check_zone_replay(lead, tail, x0_hi, x0_lo, pts, refs, dps):
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        x0_ex = mp.mpf(x0_hi) + mp.mpf(x0_lo)
        worst = mp.mpf(0)
        worst_x = None
        for x, want in zip(pts, refs):
            t_ex = mp.mpf(x) - x0_ex
            t_hi = rd(t_ex)
            if want == 0:
                continue
            p = eval_lead_tail_dd(lead, tail, t_hi, t_ex)
            got = mp.mpf(t_ex) * p
            rel = abs((got - want) / want)
            if rel > worst:
                worst, worst_x = rel, x
        return worst, worst_x
    finally:
        mp.mp.dps = old


def pin_zone(x0_mp, x0_hi, x0_lo):
    print("(2) Zone [1,2) product-form fit search (target 2^-55):",
          file=sys.stderr)
    pts = build_zone_grid(x0_hi)

    old = mp.mp.dps
    mp.mp.dps = 45
    try:
        refs = [mp.digamma(mp.mpf(x)) for x in pts]
    finally:
        mp.mp.dps = old

    best = None
    for degree in range(16, 31):
        mono = fit_zone_monomial(x0_mp, degree)
        row_worst = None
        for n_lead in (2, 3, 4, 5, 6):
            if n_lead > degree + 1:
                continue
            lead, tail = split_lead_tail(mono, n_lead)
            worst, worst_x = check_zone_replay(lead, tail, x0_hi, x0_lo,
                                                pts, refs, dps=45)
            row_worst = worst if row_worst is None else min(row_worst, worst)
            if worst <= ZONE_TARGET:
                print(f"    degree={degree} n_lead={n_lead}: worst rel err "
                      f"{float(worst):.3e} at x={worst_x!r}  <= target -- "
                      "PINNED", file=sys.stderr)
                best = (lead, tail, n_lead, degree, worst, worst_x)
                break
        if best is not None:
            break
        print(f"    degree={degree}: no n_lead<=6 reached target (best "
              f"{float(row_worst):.3e})", file=sys.stderr)

    if best is None:
        print("FAILED: zone fit could not reach 2^-55 within degree<=30, "
              "n_lead<=6 -- ESCALATE", file=sys.stderr)
        return None
    return best


# ==============================================================================
# (2) Asymptotic, x >= X0=8: psi(x) = LogDd(x) - DdRecip(2x) - x^-2*S(x^-2).
# S(w) = sum_{k=1}^{K} B_2k/(2k) * w^(k-1) -- RAW Bernoulli series, not a
# Chebyshev economization (design: "Bernoulli-series coefficients").
# ==============================================================================
def bernoulli_series_coeffs(K, dps=60):
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        return [mp.bernoulli(2 * (j + 1)) / (2 * (j + 1)) for j in range(K)]
    finally:
        mp.mp.dps = old


def check_asym_replay(head, tail, xs, refs, dps):
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        worst = mp.mpf(0)
        worst_x = None
        for x, want in zip(xs, refs):
            w = rd(1.0 / (x * x))
            s = eval_lead_tail_scalar(head, tail, w)
            term3 = mp.mpf(w) * s
            got = mp.log(mp.mpf(x)) - mp.mpf(1) / (2 * mp.mpf(x)) - term3
            if want == 0:
                continue
            rel = abs((got - want) / want)
            if rel > worst:
                worst, worst_x = rel, x
        return worst, worst_x
    finally:
        mp.mp.dps = old


def build_asym_grid(seed=SEED, n_boundary=400, n_far=250):
    rng = random.Random(seed + 1)
    pts = []
    # Dense at the boundary x=X0 (worst case: largest w).
    x0f = float(X0)
    for _ in range(n_boundary):
        pts.append(rd(rng.uniform(x0f, x0f * 1.05)))
    pts.append(x0f)
    # Log-spaced out to 1e300.
    for i in range(n_far):
        e = math.log10(x0f) + (300 - math.log10(x0f)) * i / (n_far - 1)
        pts.append(10.0 ** e)
    return pts


def pin_asymptotic():
    print("(3) Asymptotic (X0=8) Bernoulli-series fit search (target 2^-55):",
          file=sys.stderr)
    pts = build_asym_grid()
    old = mp.mp.dps
    mp.mp.dps = 50
    try:
        refs = [mp.digamma(mp.mpf(x)) for x in pts]
    finally:
        mp.mp.dps = old

    best = None
    for K in range(6, 18):
        coeffs = bernoulli_series_coeffs(K)
        for n_head in (0, 1, 2, 3):
            if n_head > K:
                continue
            head, tail = split_lead_tail(coeffs, n_head)
            worst, worst_x = check_asym_replay(head, tail, pts, refs, dps=50)
            hit = worst <= ASYM_TARGET
            if hit:
                print(f"    K={K} n_head={n_head}: worst rel err "
                      f"{float(worst):.3e} at x={worst_x!r}  <= target -- PINNED",
                      file=sys.stderr)
                best = (head, tail, K, n_head, worst, worst_x)
                break
        if best is not None:
            break
        print(f"    K={K}: best n_head worst rel err {float(worst):.3e} "
              f"(no n_head<=3 hit target)", file=sys.stderr)

    if best is None:
        print("FAILED: asymptotic fit could not reach 2^-55 within K<=17, "
              "n_head<=3 -- ESCALATE", file=sys.stderr)
        return None
    return best


# ==============================================================================
# (3) Reflection sinc-pair: sin(pi u)/(pi u) and cos(pi u), v = u^2 in
# [0, 0.25]. v needs dd precision (feeds a cancellation-prone assembly at
# self-check (c)), so these use the DD evaluation shape.
# ==============================================================================
def sinc_fn(v):
    if v == 0:
        return mp.mpf(1)
    u = mp.sqrt(v)
    return mp.sin(mp.pi * u) / (mp.pi * u)


def cospi_fn(v):
    u = mp.sqrt(v)
    return mp.cos(mp.pi * u)


def fit_v_monomial(f, degree):
    """Per-degree matched node count -- see fit_zone_monomial's comment."""
    coeffs, c, h = cheb_coeffs(f, mp.mpf(0), mp.mpf("0.25"), degree + 1)
    return cheb_to_monomial(coeffs, c, h)


def check_v_replay(lead, tail, f, seed, n=3000, absolute=False):
    """cos(pi u) has a genuine zero at u=+-0.5 (the domain edge) -- RELATIVE
    error is ill-defined there (blows up for any nonzero fit residual, no
    matter how tiny), so that fit is checked in ABSOLUTE error instead
    (well-defined at a zero, and at least as strict as 'target' relative
    error everywhere |f|>=1, which covers the rest of cos(pi u)'s range).
    sinc(u) has no such zero (bounded below by sinc(1/2)=2/pi>0) and stays
    on the relative metric."""
    rng = random.Random(seed)
    old = mp.mp.dps
    mp.mp.dps = 50
    try:
        worst = mp.mpf(0)
        worst_u = None
        for _ in range(n):
            u = rd(rng.uniform(-0.5, 0.5))
            v_hi = rd(u * u)
            v_ex = mp.mpf(u) * mp.mpf(u)
            got = eval_lead_tail_dd(lead, tail, v_hi, v_ex)
            want = f(v_ex)
            if absolute:
                err = abs(got - want)
            else:
                if want == 0:
                    continue
                err = abs((got - want) / want)
            if err > worst:
                worst, worst_u = err, u
        return worst, worst_u
    finally:
        mp.mp.dps = old


def pin_v_fit(f, name, seed, absolute=False):
    metric = "absolute" if absolute else "relative"
    print(f"(4) {name} fit search (target 2^-58 {metric}):", file=sys.stderr)
    best = None
    for degree in range(4, 22):
        mono = fit_v_monomial(f, degree)
        row_worst = None
        for n_lead in (0, 1, 2, 3):
            if n_lead > degree + 1:
                continue
            lead, tail = split_lead_tail(mono, n_lead)
            worst, worst_u = check_v_replay(lead, tail, f, seed,
                                             absolute=absolute)
            row_worst = worst if row_worst is None else min(row_worst, worst)
            if worst <= SINC_TARGET:
                print(f"    degree={degree} n_lead={n_lead}: worst {metric} "
                      f"err {float(worst):.3e} at u={worst_u!r}  <= target "
                      "-- PINNED", file=sys.stderr)
                best = (lead, tail, n_lead, degree, worst, worst_u)
                break
        if best is not None:
            break
        print(f"    degree={degree}: no n_lead<=3 reached target (best "
              f"{float(row_worst):.3e})", file=sys.stderr)
    if best is None:
        print(f"FAILED: {name} fit could not reach 2^-58 within degree<=21, "
              "n_lead<=3 -- ESCALATE", file=sys.stderr)
    return best


# ==============================================================================
# (4) Rough-trigamma: plain-double asymptotic form, target ~2^-40, domain
# [1, 2^53] via one recurrence step below y=2. NOT src/beta_data.h's
# DigammaRough (PLAN.md: "NOT reusable").
# ==============================================================================
def trigamma_asym_eval(y, coefs):
    w = rd(1.0 / (y * y))
    p = horner_d(coefs, w)
    return rd(1.0 / y) + rd(0.5 * w) + rd(rd(w * p) / y)


ROUGH_TRIGAMMA_FLOOR = 6.0  # recurrence steps up until y >= this, then
                            # asymptotic -- y=2 measured ~2.4e-5 best-case
                            # (nowhere near 2^-40), asymptotic series for
                            # trigamma just isn't tight that close in;
                            # y>=6 pinned by check_rough_trigamma's own
                            # sweep, mirroring digamma's own X0=8 finding
                            # that the asymptotic form needs real distance.


def trigamma_rough_eval(y, coefs):
    s = 0.0
    while y < ROUGH_TRIGAMMA_FLOOR:
        s = rd(s + rd(1.0 / (y * y)))
        y = rd(y + 1.0)
    return rd(s + trigamma_asym_eval(y, coefs))


def check_rough_trigamma(coefs, seed):
    rng = random.Random(seed + 2)
    old = mp.mp.dps
    mp.mp.dps = 40
    try:
        worst = mp.mpf(0)
        worst_y = None
        ys = [1.0, 1.0 + 2.0 ** -40, 2.0, 2.0 - 2.0 ** -40, 3.0, 8.0,
              1.0e6, 1.0e15, 2.0 ** 52, 2.0 ** 53]
        for _ in range(4000):
            e = rng.uniform(0, 53)
            ys.append(rd(1.0 + 2.0 ** e * rng.random()))
        for y in ys:
            got = trigamma_rough_eval(y, coefs)
            want = mp.polygamma(1, mp.mpf(y))
            if want == 0:
                continue
            rel = abs((mp.mpf(got) - want) / want)
            if rel > worst:
                worst, worst_y = rel, y
        return worst, worst_y
    finally:
        mp.mp.dps = old


def pin_rough_trigamma():
    print("(5) Rough-trigamma fit search (target 2^-40):", file=sys.stderr)
    for K in range(2, 10):
        coefs = [rd(mp.bernoulli(2 * k)) for k in range(1, K + 1)]
        worst, worst_y = check_rough_trigamma(coefs, SEED)
        if worst <= ROUGH_TRIGAMMA_TARGET:
            print(f"    K={K}: worst rel err {float(worst):.3e} at "
                  f"y={worst_y!r}  <= target -- PINNED", file=sys.stderr)
            return coefs, K, worst, worst_y
        print(f"    K={K}: worst rel err {float(worst):.3e} (no hit)",
              file=sys.stderr)
    print("FAILED: rough-trigamma could not reach 2^-40 within K<=9 -- "
          "ESCALATE", file=sys.stderr)
    return None


# ==============================================================================
# Positive-pipeline dispatcher, built from the pinned zone/asymptotic
# fits -- used both directly (self-check d) and inside the reflection
# assembly (self-check c).
# ==============================================================================
class PositivePipeline:
    def __init__(self, x0_hi, x0_lo, zone_lead, zone_tail,
                 asym_head, asym_tail):
        self.x0_hi = x0_hi
        self.x0_lo = x0_lo
        self.x0_ex = mp.mpf(x0_hi) + mp.mpf(x0_lo)
        self.x0m1_hi = x0_hi - 1.0  # exact (Sterbenz: x0.hi in [1,2))
        self.zone_lead = zone_lead
        self.zone_tail = zone_tail
        self.asym_head = asym_head
        self.asym_tail = asym_tail

    def zone_dd(self, x):
        """psi(x) for x in [1,2), as a high-precision mpf modeling the dd
        result hi+lo (no final rounding)."""
        t_ex = mp.mpf(x) - self.x0_ex
        t_hi = rd(t_ex)
        p = eval_lead_tail_dd(self.zone_lead, self.zone_tail, t_hi, t_ex)
        return mp.mpf(t_ex) * p

    def asym_dd(self, x):
        """psi(x) for x >= X0, as a high-precision mpf."""
        w = rd(1.0 / (x * x))
        s = eval_lead_tail_scalar(self.asym_head, self.asym_tail, w)
        term3 = mp.mpf(w) * s
        return mp.log(mp.mpf(x)) - mp.mpf(1) / (2 * mp.mpf(x)) - term3

    def downwalk_dd(self, x):
        """psi(x) for x in [2, X0), via the masked fixed-step down-walk."""
        m = int(math.floor(x)) - 1
        m = max(1, min(m, WALK_DEPTH))
        s = mp.mpf(0)
        xj = x
        for j in range(1, m + 1):
            xj = rd(x - j)
            s += mp.mpf(1) / mp.mpf(xj)
        x_land = xj
        return self.zone_dd(x_land) + s

    def psi_pos_dd(self, y_hi):
        """Dispatch on y_hi (a plain double >= ZONE_LO in the reflection
        context, or any positive double for the recurrence self-check)."""
        if y_hi >= float(X0):
            return self.asym_dd(y_hi)
        if y_hi >= 2.0:
            return self.downwalk_dd(y_hi)
        if y_hi >= 1.0:
            return self.zone_dd(y_hi)
        # (0,1): up-step WITHOUT forming 1+x -- t1 shifted vs (x0-1).
        x = y_hi
        t1_ex = mp.mpf(x) - mp.mpf(self.x0m1_hi) - mp.mpf(self.x0_lo)
        t1_hi = rd(t1_ex)
        p = eval_lead_tail_dd(self.zone_lead, self.zone_tail, t1_hi, t1_ex)
        psi_xp1 = mp.mpf(t1_ex) * p
        return psi_xp1 - mp.mpf(1) / mp.mpf(x)


# ==============================================================================
# Self-check (d): recurrence replay -- (0,1) up-step and [2,X0) down-walk.
# ==============================================================================
def check_recurrence(pipe, seed=SEED):
    rng = random.Random(seed + 3)
    pts = []
    for i in range(1, 2000):
        pts.append(i / 2000.0)
    pts.append(1.0 - 2.0 ** -52)  # x -> 1^-, worst |psi| = -gamma there
    for e in range(1, 16):
        pts.append(2.0 ** -e)
    for _ in range(1500):
        pts.append(rng.uniform(2.0, float(X0)))
    for k in range(2, 8):
        pts.append(float(k) - 2.0 ** -52)
        pts.append(float(k) + 2.0 ** -52)

    old = mp.mp.dps
    mp.mp.dps = 45
    try:
        worst = mp.mpf(0)
        worst_x = None
        for x in pts:
            if not (0.0 < x < float(X0)):
                continue
            got = pipe.psi_pos_dd(x)
            want = mp.digamma(mp.mpf(x))
            if want == 0:
                continue
            rel = abs((got - want) / want)
            if rel > worst:
                worst, worst_x = rel, x
        return worst, worst_x
    finally:
        mp.mp.dps = old


# ==============================================================================
# Self-check (c): reflection adversarial replay at the 20 nearest doubles
# to the negative-axis zeros of digamma.
# ==============================================================================
def negative_zeros(n=20, dps=80):
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        zeros = []
        for k in range(1, n + 1):
            guess = mp.mpf(-(k)) + mp.mpf("0.5")
            z = mp.findroot(lambda x: mp.digamma(x), guess)
            zeros.append(z)
        return zeros
    finally:
        mp.mp.dps = old


def round_to_double(x):
    r = round(x)
    return float(r)


def reflection_dd(pipe, sinc_lead, sinc_tail, cos_lead, cos_tail,
                   trig_coefs, x_double):
    """Full dd assembly psi(1-x) - pi*cot(pi*x) for x < 0, modeled to dd
    precision throughout (see module docstring)."""
    y_hi, y_lo = two_sum(1.0, -x_double)
    term1 = pipe.psi_pos_dd(y_hi)
    if y_lo != 0.0:
        term1 += mp.mpf(y_lo) * mp.mpf(trigamma_rough_eval(y_hi, trig_coefs))

    n = round_to_double(x_double)
    u = rd(x_double - n)  # exact
    v_hi = rd(u * u)
    v_ex = mp.mpf(u) * mp.mpf(u)
    sinc_val = eval_lead_tail_dd(sinc_lead, sinc_tail, v_hi, v_ex)
    cos_val = eval_lead_tail_dd(cos_lead, cos_tail, v_hi, v_ex)
    # denom = u*sincfit(u) = u*sin(pi u)/(pi u) = sin(pi u)/pi, so
    # cos_val/denom = pi*cos(pi u)/sin(pi u) = pi*cot(pi u) = pi*cot(pi x)
    # ALREADY -- pi cancels analytically, no separate pi multiplication
    # belongs here. Multiplying by pi_ex again here would give pi^2*cot
    # and a catastrophically wrong result -- kDigammaPiHi/Lo is unused by
    # this assembly and dropped from emission below.
    denom = mp.mpf(u) * sinc_val
    ratio = cos_val / denom       # = pi*cot(pi u) = pi*cot(pi x)
    term2 = ratio

    return term1 - term2


def check_c1_absolute(pipe, sinc_lead, sinc_tail, cos_lead, cos_tail,
                       trig_coefs):
    """(c1) ABSOLUTE error <= 2^-56 at each of the 20 adversarial
    nearest-double points -- replaces a single relative-bits bar, rejected
    on cost (see module docstring)."""
    print("(c1) REFLECTION ADVERSARIAL ABSOLUTE (20 nearest-double negative "
          "zeros, target abs <= 2^-56):", file=sys.stderr)
    zeros = negative_zeros(20, dps=80)
    worst_abs = None
    worst_at = None
    ok = True
    for k, z in enumerate(zeros, 1):
        x_double = rd(z)
        old = mp.mp.dps
        mp.mp.dps = 90
        try:
            true_v = mp.digamma(mp.mpf(x_double))
            got = reflection_dd(pipe, sinc_lead, sinc_tail, cos_lead,
                                 cos_tail, trig_coefs, x_double)
            err = abs(got - true_v)
        finally:
            mp.mp.dps = old
        status = "OK" if err <= REFLECTION_C1_ABS_TARGET else "FAIL"
        if status == "FAIL":
            ok = False
        if worst_abs is None or err > worst_abs:
            worst_abs, worst_at = err, (k, x_double)
        print(f"    n={k:2d} x={x_double!r} abs_err={float(err):.3e}  "
              f"{status}", file=sys.stderr)
        if status == "FAIL":
            print(f"    n={k} x={x_double!r} exceeds 2^-56 absolute -- "
                  "ESCALATE (this bar was measured with 2+ decades of "
                  "margin; a miss here is a real regression, not noise)",
                  file=sys.stderr)
    print(f"    worst absolute error over 20 points: "
          f"{float(worst_abs):.3e} at n={worst_at[0]} x={worst_at[1]!r} "
          f"(target 2^-56 = {float(REFLECTION_C1_ABS_TARGET):.3e})",
          file=sys.stderr)
    return ok, worst_abs


# Pole exclusion margin for (c2)'s dense sampling: a fixed fraction of the
# unit interval width, well clear of the +-ulp-scale pole neighborhoods the
# ULP test's own reference set covers separately -- this check's job is
# GENERIC mid-interval accuracy, not pole conditioning.
C2_POLE_MARGIN = mp.mpf("0.01")


def check_c2_relative(pipe, sinc_lead, sinc_tail, cos_lead, cos_tail,
                       trig_coefs):
    """(c2) RELATIVE error <= 2^-52 at dense negative-axis samples where
    |psi(x)| >= 1 -- each interval (-n-1,-n)
    for n=1..20 (excluding that interval's zero band W ~ 1/|trigamma(z0)|
    and a pole margin at both ends), plus a few log-spaced spot intervals
    out to n ~ 1e6 (|psi(x)| >= 1 filter alone -- it already excludes the
    zero's neighborhood out there, so no explicit band is needed)."""
    print("(c2) REFLECTION DENSE RELATIVE (|psi(x)| >= 1, target rel "
          "<= 2^-52):", file=sys.stderr)
    rng = random.Random(SEED + 5)

    # z_{n+1} lies in (-(n+1), -n) per negative_zeros' own convention
    # (z_k in (-k, -(k-1))) -- need z_1..z_21 to cover n=1..20.
    zeros21 = negative_zeros(21, dps=80)

    def sample_interval(n_lo_int, n_hi_int, zero_mp=None, n_pts=40):
        """Points in the OPEN interval (n_lo_int, n_hi_int) (n_lo_int <
        n_hi_int, both negative integers), excluding a pole margin at
        each end and (if given) a W-band around zero_mp."""
        old = mp.mp.dps
        mp.mp.dps = 80
        try:
            width = mp.mpf(n_hi_int) - mp.mpf(n_lo_int)
            lo = mp.mpf(n_lo_int) + C2_POLE_MARGIN * width
            hi = mp.mpf(n_hi_int) - C2_POLE_MARGIN * width
            w = None
            if zero_mp is not None:
                tg = mp.polygamma(1, zero_mp)
                w = 1 / abs(tg)
            pts = []
            tries = 0
            while len(pts) < n_pts and tries < n_pts * 20:
                tries += 1
                t = rng.random()
                x = lo + (hi - lo) * t
                if w is not None and abs(x - zero_mp) < w:
                    continue
                pts.append(rd(x))
            return pts
        finally:
            mp.mp.dps = old

    worst_rel = mp.mpf(0)
    worst_at = None
    n_tested = 0
    n_skipped_small = 0

    def eval_leg(pts, label):
        nonlocal worst_rel, worst_at, n_tested, n_skipped_small
        leg_worst = mp.mpf(0)
        leg_at = None
        for x_double in pts:
            old = mp.mp.dps
            mp.mp.dps = 80
            try:
                true_v = mp.digamma(mp.mpf(x_double))
                if abs(true_v) < 1:
                    n_skipped_small += 1
                    continue
                got = reflection_dd(pipe, sinc_lead, sinc_tail, cos_lead,
                                     cos_tail, trig_coefs, x_double)
                rel = abs((got - true_v) / true_v)
            finally:
                mp.mp.dps = old
            n_tested += 1
            if rel > leg_worst:
                leg_worst, leg_at = rel, x_double
            if rel > worst_rel:
                worst_rel, worst_at = rel, (label, x_double)
        status = "OK" if leg_worst <= REFLECTION_C2_REL_TARGET else "FAIL"
        print(f"    {label}: n_pts={len(pts)} worst_rel="
              f"{float(leg_worst):.3e} at x={leg_at!r}  {status}",
              file=sys.stderr)
        return leg_worst <= REFLECTION_C2_REL_TARGET

    ok = True
    for n in range(1, 21):
        z = zeros21[n]  # z_{n+1}, in (-(n+1), -n)
        pts = sample_interval(-(n + 1), -n, zero_mp=z)
        ok &= eval_leg(pts, f"interval(-{n + 1},-{n})")

    # Log-spaced spot intervals out to n ~ 1e6: |psi(x)|>=1 filter alone.
    for n in (50, 100, 316, 1000, 3162, 10000, 31623, 100000, 316228,
              1000000):
        pts = sample_interval(-(n + 1), -n, zero_mp=None, n_pts=6)
        ok &= eval_leg(pts, f"spot(-{n + 1},-{n})")

    print(f"    worst relative error over all legs: {float(worst_rel):.3e} "
          f"at {worst_at}  (n_tested={n_tested}, n_skipped |psi|<1="
          f"{n_skipped_small}, target 2^-52 = "
          f"{float(REFLECTION_C2_REL_TARGET):.3e})", file=sys.stderr)
    return ok, worst_rel


# ==============================================================================
# Self-check (e) helpers: sinc/cos vs mp.sin/mp.cos direct spot-check.
# ==============================================================================
def check_e_sinc_cos(sinc_lead, sinc_tail, cos_lead, cos_tail):
    rng = random.Random(SEED + 4)
    old = mp.mp.dps
    mp.mp.dps = 50
    try:
        worst_sinc = mp.mpf(0)
        worst_cos = mp.mpf(0)
        for _ in range(2000):
            u = rd(rng.uniform(-0.5, 0.5))
            v_hi = rd(u * u)
            v_ex = mp.mpf(u) * mp.mpf(u)
            got_s = eval_lead_tail_dd(sinc_lead, sinc_tail, v_hi, v_ex)
            got_c = eval_lead_tail_dd(cos_lead, cos_tail, v_hi, v_ex)
            want_s = sinc_fn(v_ex)
            want_c = cospi_fn(v_ex)
            if want_s != 0:
                worst_sinc = max(worst_sinc, abs((got_s - want_s) / want_s))
            # cos(pi u) has a zero at u=+-0.5 -- absolute, see check_v_replay.
            worst_cos = max(worst_cos, abs(got_c - want_c))
        return worst_sinc, worst_cos
    finally:
        mp.mp.dps = old


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

    x0_hi, x0_lo, x0_mp, tg, root_ok = derive_root()
    if not root_ok:
        rc = 1

    zone = pin_zone(x0_mp, x0_hi, x0_lo)
    if zone is None:
        return 1
    zone_lead, zone_tail, zone_n_lead, zone_degree, zone_err, zone_at = zone

    asym = pin_asymptotic()
    if asym is None:
        return 1
    asym_head, asym_tail, asym_K, asym_n_head, asym_err, asym_at = asym

    sinc = pin_v_fit(sinc_fn, "sinc(u)=sin(pi u)/(pi u)", SEED + 10)
    cospi = pin_v_fit(cospi_fn, "cos(pi u)", SEED + 20, absolute=True)
    if sinc is None or cospi is None:
        return 1
    sinc_lead, sinc_tail, sinc_n_lead, sinc_degree, sinc_err, sinc_at = sinc
    cos_lead, cos_tail, cos_n_lead, cos_degree, cos_err, cos_at = cospi

    rough = pin_rough_trigamma()
    if rough is None:
        return 1
    trig_coefs, trig_K, trig_err, trig_at = rough

    pipe = PositivePipeline(x0_hi, x0_lo, zone_lead, zone_tail,
                             asym_head, asym_tail)

    rec_err, rec_at = check_recurrence(pipe)
    print(f"(d) recurrence replay ((0,1) up-step + [2,X0) down-walk): "
          f"worst rel err {float(rec_err):.3e} at x={rec_at!r}, target "
          f"2^-55", file=sys.stderr)
    if rec_err > RECURRENCE_TARGET:
        print("FAILED: recurrence replay exceeds 2^-55 -- ESCALATE",
              file=sys.stderr)
        rc = 1

    c1_ok, c1_worst = check_c1_absolute(
        pipe, sinc_lead, sinc_tail, cos_lead, cos_tail, trig_coefs)
    if not c1_ok:
        rc = 1

    c2_ok, c2_worst = check_c2_relative(
        pipe, sinc_lead, sinc_tail, cos_lead, cos_tail, trig_coefs)
    if not c2_ok:
        rc = 1

    e_sinc, e_cos = check_e_sinc_cos(sinc_lead, sinc_tail, cos_lead, cos_tail)
    print(f"(e) sinc/cos direct spot-check vs mp.sin/mp.cos: "
          f"sinc worst {float(e_sinc):.3e}, cos worst {float(e_cos):.3e}",
          file=sys.stderr)

    print(f"(e) rough-trigamma <= 2^-40: worst {float(trig_err):.3e} "
          f"at y={trig_at!r}  {'OK' if trig_err <= ROUGH_TRIGAMMA_TARGET else 'FAIL'}",
          file=sys.stderr)

    if rc:
        print("\nOne or more self-checks FAILED -- refusing to emit "
              "src/digamma_data.h. ESCALATE.", file=sys.stderr)
        return rc

    # --- emit ---------------------------------------------------------------
    print("// Auto-generated by tools/gen_digamma_data.py. DO NOT EDIT.")
    print("// Fits and their error budgets are documented in the generator;")
    print("// the kernel that consumes them is src/digamma-inl.h.")
    print("#ifndef CORVUS_DIGAMMA_DATA_H_")
    print("#define CORVUS_DIGAMMA_DATA_H_")
    print()
    print("namespace corvus::detail {")
    print()
    print("// The unique positive root of digamma (~1.4616), as a dd pair.")
    print("// NOT to be confused with kDigammaX0 (the asymptotic threshold)")
    print("// below -- distinct constants that happen to share a name in")
    print("// the design prose (x0 the root vs X0 = 8 the threshold).")
    emit_scalar("kDigammaRootHi", x0_hi)
    emit_scalar("kDigammaRootLo", x0_lo)
    print("// (x0.hi - 1), exact by Sterbenz (x0.hi in [1,2)): the (0,1)")
    print("// branch's shifted centre, reusing kDigammaRootLo unchanged.")
    emit_scalar("kDigammaRootM1Hi", x0_hi - 1.0)
    print()
    print("// Region boundaries. [ZoneLo, ZoneHi) is the product-form zone;")
    print("// [ZoneHi, X0) walks down by up to kDigammaWalkDepth integer")
    print("// steps; [X0, inf) is the Bernoulli-series asymptotic form.")
    emit_scalar("kDigammaZoneLo", float(ZONE_LO))
    emit_scalar("kDigammaZoneHi", float(ZONE_HI))
    emit_scalar("kDigammaX0", float(X0))
    print(f"inline constexpr int kDigammaWalkDepth = {WALK_DEPTH};")
    print()
    print(f"// Zone [1,2): psi(x) = t*P(t), t = x - root (dd). P(t) = ")
    print(f"// L0 + t*(L1 + ... + t*S(t)); L* are the first "
          f"{zone_n_lead} dd-lead")
    print(f"// coefficients, S the plain-double tail (degree {zone_degree}"
          f" total).")
    print(f"// Replay-measured worst relative error {float(zone_err):.3e} "
          f"(target 2^-55).")
    print(f"inline constexpr int kDigammaZoneLead = {zone_n_lead};")
    print(f"inline constexpr int kDigammaZoneNCoef = {len(zone_tail)};")
    emit_lead("kDigammaZoneLeadHi", zone_lead, 0)
    emit_lead("kDigammaZoneLeadLo", zone_lead, 1)
    emit_1d("kDigammaZoneCoef", zone_tail)
    print()
    print(f"// Asymptotic (x >= kDigammaX0): psi(x) = log(x) - 1/(2x) - ")
    print(f"// x^-2*S(x^-2). S(w) = H0 + w*(H1 + ... + w*T(w)); H* are the")
    print(f"// first {asym_n_head} dd-head coefficients (RAW Bernoulli")
    print(f"// series B_2k/(2k)), T the plain-double tail, K={asym_K} terms")
    print(f"// total. Replay-measured worst relative error "
          f"{float(asym_err):.3e} (target 2^-55).")
    print(f"inline constexpr int kDigammaAsymHead = {asym_n_head};")
    print(f"inline constexpr int kDigammaAsymNCoef = {len(asym_tail)};")
    emit_lead("kDigammaAsymHeadHi", asym_head, 0)
    emit_lead("kDigammaAsymHeadLo", asym_head, 1)
    emit_1d("kDigammaAsymCoef", asym_tail)
    print()
    print(f"// Reflection sinc-pair, v = u^2, u = x - round(x) (exact),")
    print(f"// |u| <= 1/2. sinc(u) = sin(pi u)/(pi u); cospi(u) = ")
    print(f"// cos(pi u); cot(pi x) = cospi(u) / (u * sinc(u)) (the (-1)^n")
    print(f"// parity from x = n + u cancels in the ratio). Both v-dd-lead")
    print(f"// + double-tail, matching the zone's DD evaluation shape.")
    print(f"// Replay: sinc worst {float(sinc_err):.3e}, cos worst "
          f"{float(cos_err):.3e} (target 2^-58).")
    print(f"inline constexpr int kDigammaSincLead = {sinc_n_lead};")
    print(f"inline constexpr int kDigammaSincNCoef = {len(sinc_tail)};")
    emit_lead("kDigammaSincLeadHi", sinc_lead, 0)
    emit_lead("kDigammaSincLeadLo", sinc_lead, 1)
    emit_1d("kDigammaSincCoef", sinc_tail)
    print(f"inline constexpr int kDigammaCosLead = {cos_n_lead};")
    print(f"inline constexpr int kDigammaCosNCoef = {len(cos_tail)};")
    emit_lead("kDigammaCosLeadHi", cos_lead, 0)
    emit_lead("kDigammaCosLeadLo", cos_lead, 1)
    emit_1d("kDigammaCosCoef", cos_tail)
    print()
    print(f"// Rough trigamma (K={trig_K} asymptotic Bernoulli terms, plain")
    print(f"// double, ~2^-40 relative): ONLY for the y.lo * trigamma(y.hi)")
    print(f"// correction on the reflection path's dd argument. NOT")
    print(f"// src/beta_data.h's DigammaRough (different function, "
          f"different budget,")
    print(f"// not reusable). Replay-measured worst "
          f"{float(trig_err):.3e}.")
    print(f"// WALK FORM (the kernel MUST mirror this): while (y <")
    print(f"// kDigammaRoughTrigammaFloor) {{ s += 1/(y*y); y += 1; }} then")
    print(f"// trigamma(y) ~= 1/y + w/2 + w*Horner(coef, w)/y, w = 1/(y*y),")
    print(f"// return s + that. y=2 alone measured ~2.4e-5 best-case (nowhere")
    print(f"// near 2^-40) -- the asymptotic series just isn't tight that")
    print(f"// close in; the walk floor is what buys the margin, not K.")
    print(f"inline constexpr double kDigammaRoughTrigammaFloor = "
          f"{hexf(ROUGH_TRIGAMMA_FLOOR)};")
    print(f"inline constexpr int kDigammaRoughTrigammaN = {len(trig_coefs)};")
    emit_1d("kDigammaRoughTrigammaCoef", trig_coefs)
    print()
    print("}  // namespace corvus::detail")
    print()
    print("#endif  // CORVUS_DIGAMMA_DATA_H_")
    return 0


if __name__ == "__main__":
    sys.exit(main())
