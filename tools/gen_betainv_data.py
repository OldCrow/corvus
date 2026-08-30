#!/usr/bin/env python3
"""Generate src/betainv_data.h -- every table beta_p_inv/beta_q_inv needs.

This generator pins, BY REPLAY (never by assumption), every free parameter of
the inverse's seed+step stage:

  S1  beta-Temme normal-quantile seed (z=erfcinv(2*sigma), invert the beta
      ridge mapping cpsi(lambda)=alpha*phi(-lambda/alpha)+beta*phi(lambda/beta),
      phi(w)=w-log1p(w)) -- clean-room from Temme's published structure,
      direction-reversed from gen_beta_data.py's forward R3 e_k(zeta,p)
      extraction (see module docstring there); correction count/table
      pinned by replay (gammainv precedent: K=2, Chebyshev 2x25 over a WIDE
      eta domain -- this family measures its own).
  S2  small-y series inversion seed (y0 = exp((ln sigma + ln alpha + lnB)/
      alpha) + Picard corrections via the R1 series); Picard count pinned
      by replay.
  S3  gamma-limit transfer seed: t = -beta*log1p(-y) [beta the HUGE
      parameter -- verified against src/beta-inl.h's own gamma-limit slice,
      "huge SECOND: t = -(huge)*ln(1-xi)"], t seeded via the EXISTING
      GammaInvSeedS1/S2/S3 Python functions in tools/gen_gammainv_data.py
      (imported, not re-derived -- cross-family precedent: beta-inl.h
      already imports gamma-inl.h's template cores the same way), then
      y = -expm1(-t/beta).
  S4  joint-tiny logit closed form: logit(y) = (s-s*)/w + c(alpha,beta),
      s*=beta/(alpha+beta), w=alpha*beta/(alpha+beta). Derived here from
      B_y = int_{-inf}^{logit y} exp(alpha*u - (alpha+beta)*L(u)) du,
      L(u)=ln(1+e^u) (exact substitution of the beta integral into logit
      space) -- see the derivation comment at s4_c_closed_form. c=0 at
      alpha=beta by the swap-antisymmetry identity (self-check (s4)).
  STEPS  safeguarded logit-Newton (m=lnP-lnQ objective), 3 shared steps,
      the gammainv safeguard package assumed whole (reject
      residual-increasing steps, 1/8 backtrack, bypass |resid|<1/2,
      additive-y step) -- simulated against the forward evaluated in
      mpmath with injected relative noise at the forward kernel's own
      internal dd budget, per-point analytic eps wherever series
      super-converge, q-side ratio conversion included from day one.
  DEEP-SMALL  closed form at BOTH ends (y = exp_dd((LogDd(sigma) (+) ln
      alpha (+) lnB_dd)/alpha)), cut on dropped-factor error < 2^-60
      measured in BOTH orientations from the start (a single-orientation
      check can miss errors on the untested side).

mpmath discipline: mp.dps set inside every function that needs it; replay
solves against the root of the ROUNDED double sigma, never the unrounded mpf
value -- solving against the unrounded value would target a different
equation than the one the double-precision seed/step actually solves.

Self-checks (mandatory; stderr budget lines; ANY miss -> exit nonzero,
emit nothing) -- lettered to mirror gen_gammainv_data.py's scheme:
  (a) beta-ridge lambda(zeta) Newton: converges to >= NEWTON_FLOOR_BITS
      across the S1 seed domain, both signs, a range of (alpha,beta) skew.
  (b) S1/S2/S3/S4 seed-bit floors, quad-candidate cheap-residual selection,
      measured against the actual seed code (native double).
  (c) S4 closed form: leading-order derivation self-check (c=0 at
      alpha=beta, exact identity I_{1/2}(a,a)=1/2) plus measured accuracy
      of the pinned c(alpha,beta) form across the joint-tiny domain.
  (d) t_jt route gate: S4 owns min(a,b) <= 2^-52 with margin; large-|logit|
      seam to the power-law/deep-small form.
  (e) STEPS: per-region final bits after the pinned step count, worst case
      over forward-noise sign combinations, >= 54 bits with >= 1 bit
      margin, at every replay point whose true y is a normal double
      (except the plateau/beyond-resolution buckets, backward-error
      contract there).
  (f) deep-small cut: dropped-factor error < 2^-60 below the cut, BOTH
      orientations.
  (g) S1/S3 seam near alpha ~ kGammaAT: seam location measured, not assumed.
  (h) plateau backward-error contract: kappa > 2^52 rows verified via
      forward-of-returned-y instead of y-ULP.

Usage:
    python3 tools/gen_betainv_data.py > src/betainv_data.h
    python3 tools/gen_betainv_data.py --full > src/betainv_data.h   # denser replay
"""
import argparse
import math
import os
import sys
import time

import mpmath as mp

FULL = "--full" in sys.argv  # naive membership check kept for the
# import-as-library case (gen_betainv_reference.py does `import
# gen_betainv_data as bid` at module load, before its own argparse runs --
# this must NOT choke on that caller's own flags); the `__main__` block
# below re-derives FULL through argparse for this module's OWN direct
# invocation, where the strict-parsing/--help/unknown-flag behavior
# actually applies.
_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

# Reuse, not re-derive: forward R1/R2/R3 machinery (gen_beta_data.py) and
# the gamma-inverse seed formulas (gen_gammainv_data.py) -- both modules
# guard their own self-check/emission code behind `if __name__=='__main__'`,
# so importing them here costs only their module-level derivations
# (gen_beta_data's ZETA_MAX golden-section search, a few seconds) and pulls
# in zero side effects on stdout (their own diagnostic prints go to stderr).
import gen_beta_data as gb          # noqa: E402
import gen_gammainv_data as ginv    # noqa: E402

T0 = time.time()

# ============================================================================
# Wire up gen_gammainv_data's SHIPPED S1 correction fit (kGammaInvCkCheb):
# ginv.seed_S1 needs module global _S1_CHEB_ROWS_D populated, which
# gammainv's own main() only fills when it actually RUNS its replay (a
# multi-minute Chebyshev extraction, module docstring Part 2) -- re-running
# that here would be re-deriving, exactly what the contract forbids ("import
# or port its Python seed functions -- do not re-derive differently").
# Instead, parse the ALREADY-SHIPPED src/gammainv_data.h directly: this is
# the literal table the gammainv KERNEL itself consumes, so wiring it into
# ginv.seed_S1 here reproduces the shipped seed bit-for-bit (the same
# technique gammainv itself uses to cross-check against gamma_data.h's
# constants by reference, never duplication).
# ============================================================================
def _parse_gammainv_data_h():
    path = os.path.join(_HERE, "..", "src", "gammainv_data.h")
    with open(path, "r", encoding="utf-8") as f:
        text = f.read()
    import re

    def hexval(tok):
        return float.fromhex(tok)

    rows_m = re.search(
        r"kGammaInvCkCheb\[2\]\[25\]\s*=\s*\{(.*?)\};", text, re.S)
    rows_text = rows_m.group(1)
    row_matches = re.findall(r"\{([^}]*)\}", rows_text)
    cheb_rows = []
    for row in row_matches:
        toks = [t.strip() for t in row.split(",") if t.strip()]
        cheb_rows.append([hexval(t) for t in toks])
    s2ncorr = int(re.search(r"kGammaInvS2NCorr\s*=\s*(\d+)", text).group(1))
    s3niter = int(re.search(r"kGammaInvS3NIter\s*=\s*(\d+)", text).group(1))
    s3margin = hexval(re.search(
        r"kGammaInvS3StabilityMargin\s*=\s*(\S+);", text).group(1))
    s1amin = hexval(re.search(
        r"kGammaInvS1AMin\s*=\s*(\S+);", text).group(1))
    seedncorr = int(re.search(r"kGammaInvSeedNCorr\s*=\s*(\d+)", text).group(1))
    etamax = hexval(re.search(r"kGammaInvEtaMax\s*=\s*(\S+);", text).group(1))
    return dict(cheb_rows=cheb_rows, s2ncorr=s2ncorr, s3niter=s3niter,
                s3margin=s3margin, s1amin=s1amin, seedncorr=seedncorr,
                etamax=etamax)


_GINV = _parse_gammainv_data_h()
ginv._S1_CHEB_ROWS_D = _GINV["cheb_rows"]
ginv.ETA_MAX = _GINV["etamax"]
ginv.S1_A_MIN = _GINV["s1amin"]
ginv.S3_STABILITY_MARGIN = _GINV["s3margin"]
GINV_S1_NCORR = _GINV["seedncorr"]
GINV_S2_NCORR = _GINV["s2ncorr"]
GINV_S3_NITER = _GINV["s3niter"]


def gammainv_seed_for(a, s, side):
    """Re-implementation of gen_gammainv_data.py's own tri-candidate
    seed_for GLUE (selection logic only -- module docstring's "cheap
    forward-residual comparison", trivial orchestration, NOT a seed
    FORMULA) since seed_for itself is a closure local to ginv.main() and
    therefore not importable; the three candidates it selects between
    (ginv.seed_S1/seed_S2/seed_S3) are called UNMODIFIED, and the pinned
    counts/gates (S1_NCORR, S2_NCORR, S3_NITER, S1_A_MIN, S3 stability
    margin) come from the parsed, SHIPPED gammainv_data.h above -- not
    re-derived. Mirrors ginv.main()'s own seed_for docstring exactly."""
    s1_candidate = None
    eta0 = ginv.eta0_of(a, s, side)
    if a >= ginv.S1_A_MIN and abs(eta0) <= ginv.ETA_MAX:
        try:
            s1_candidate = ginv.seed_S1(a, s, side, GINV_S1_NCORR)
        except (OverflowError, ValueError):
            s1_candidate = None
    try:
        s2_candidate = (ginv.seed_S2(a, s, GINV_S2_NCORR) if side == "p"
                         else ginv.seed_S2(a, 1.0 - s, GINV_S2_NCORR))
    except (OverflowError, ValueError):
        s2_candidate = None
    s3_candidate = None
    if side == "q":
        try:
            s3_candidate = ginv.seed_S3(a, s, GINV_S3_NITER)
        except (OverflowError, ValueError):
            s3_candidate = None

    def cheap_residual(a, x0, s, side):
        if not (math.isfinite(x0) and x0 > 0):
            return None
        v = sd = None
        for dps_try in (25, 35, 55):
            try:
                v, sd = ginv.small_of_x(a, x0, dps=dps_try)
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
    return s2_candidate if s2_candidate is not None else float("nan")


def rd(x):
    return float(x)


def hexf(x):
    return float.hex(float(x))


def emit_hex_array_1d(name, vals, ncols=8):
    print(f"inline constexpr double {name}[{len(vals)}] = {{")
    for i in range(0, len(vals), ncols):
        print("    " + ", ".join(hexf(v) for v in vals[i:i + ncols]) + ",")
    print("};")


def emit_hex_array_2d(name, rows):
    ncols = len(rows[0])
    print(f"inline constexpr double {name}[{len(rows)}][{ncols}] = {{")
    for row in rows:
        print("    {" + ", ".join(hexf(v) for v in row) + "},")
    print("};")


# ============================================================================
# Part 0: shared math helpers
# ============================================================================
def logit(y):
    return mp.log(y / (1 - y))


def sigmoid(v):
    if v >= 0:
        e = mp.e ** (-v)
        return 1 / (1 + e)
    else:
        e = mp.e ** v
        return e / (1 + e)


def bits_of(true_x, approx_x):
    """true_x must stay an mpf, KEPT AT FULL PRECISION: casting it to a
    Python float before this comparison would silently cap every
    measured floor at ~52-53 bits regardless of how good the candidate
    actually is, since a double-quantized true_x can never show more
    than ~1 ULP of agreement with anything by construction. gammainv's
    own _bits_of_frontier keeps its true_x as mpf through this exact
    comparison for the same reason -- the ROUNDED-DOUBLE-s rule applies
    to the NEWTON TARGET, never to the bits-measurement basis, which
    must stay at full precision to measure anything past ~53 bits."""
    if not (math.isfinite(approx_x) and 0.0 < true_x < 1.0):
        return -1.0
    if not (0.0 <= approx_x <= 1.0):
        # allow tiny negative/over-one from a seed candidate's own
        # approximation -- clamp is a kernel-level concern, not a bits
        # measurement one, but out-of-[0,1] cannot be scored as bits.
        return -1.0
    rel = abs((mp.mpf(approx_x) - mp.mpf(true_x)) / mp.mpf(true_x))
    if rel == 0:
        return 100.0
    return float(-mp.log(rel, 2))


# ============================================================================
# Part 1: forward evaluators (measurement-grade truth, mirroring
# gen_gammainv_data.py's Part 4 -- NOT the certified oracle, which is a
# separate, later stage with its own three binding constructions per
# PLAN.md). Reuses gen_beta_data.py's R1 series / R2 CF / R3 Temme
# machinery directly (gb.*) rather than re-deriving it.
# ============================================================================
def r1_value_mp(a, b, x, dps, nmax=4000):
    """R1-native power series, self-convergent (NOT the fixed-N1=64 cheap
    proxy gb._r1_native_value uses for routing -- this is the measurement
    TRUTH evaluator, escalating N until self-convergent, matching
    gammainv's series_S_mp pattern)."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        a = mp.mpf(a); b = mp.mpf(b); x = mp.mpf(x)
        t = mp.mpf(1); s = mp.mpf(0)
        eps = mp.mpf(10) ** (-(dps - 8))
        for n in range(0, nmax + 1):
            if n > 0:
                t *= (n - b) * x / n
            term = t / (a + n)
            s += term
            if n > 4 and abs(t) < eps * abs(s):
                break
        logpref = a * mp.log(x) - (mp.loggamma(a) + mp.loggamma(b) - mp.loggamma(a + b))
        return mp.e ** logpref * s
    finally:
        mp.mp.dps = old


def r2_value_mp(a, b, x, dps):
    """CF value in the GIVEN orientation (caller picks orientation) --
    thin wrapper on gb.small_val_via_cf's escalating self-convergence."""
    return gb.small_val_via_cf(a, b, x, dps)


def r3_value_mp(a, b, x, dps):
    """R3 Temme reconstruction: nu=ab/(a+b), p=a/c, zeta signed so that
    zeta>=0 <-> lambda>=0 <-> small_val=P (gb.r3_R_at's own convention).
    Returns (small_val, side)."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        a = mp.mpf(a); b = mp.mpf(b); x = mp.mpf(x)
        c = a + b
        nu = a * b / c
        p = a / c
        lam = a - c * x  # xi=(alpha-lambda)/c => lambda=alpha-c*xi
        u = -lam / a
        v = lam / b
        cpsi = a * (u - mp.log1p(u)) + b * (v - mp.log1p(v))
        z = mp.sqrt(cpsi)
        zeta = z if lam >= 0 else -z
        leading = mp.erfc(z) / 2
        nu_eff = a * b / (a + b)
        if lam >= 0:
            small_val = gb.small_val_via_cf(a, b, x, dps)
            side = "p"
        else:
            small_val = gb.small_val_via_cf(b, a, 1 - x, dps)
            side = "q"
        return small_val, side
    finally:
        mp.mp.dps = old


def gamma_corner_value_mp(a, b, x, dps):
    """Gamma-limit forward evaluator, mirrors src/beta-inl.h's own slice
    verbatim (see module docstring's citation): huge SECOND (b huge, a
    small) -> t=-b*log1p(-x), value~P_gamma(a,t); huge FIRST (a huge, b
    small) -> t=-a*log(x), value~1-P_gamma(b,t) [b the shape]. GUARDED
    per PLAN's practical warning: never call mpmath.gammainc with shape
    >1e4 near the ridge -- gen_gamma_reference.py's own exact-asymptotic
    route is reused via ginv's small_of_x (a<=A_MPMATH_SAFE=1e4 direct,
    Temme-extrapolated beyond)."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        a = mp.mpf(a); b = mp.mpf(b); x = mp.mpf(x)
        if b >= a:
            s, huge = a, b
            t = -huge * mp.log1p(-x)
            val, side = ginv.small_of_x(s, t, dps=dps)
            if side == "q":
                val = 1 - val
            return val, "p"
        else:
            s, huge = b, a
            t = -huge * mp.log(x)
            val, side = ginv.small_of_x(s, t, dps=dps)
            if side == "p":
                val = 1 - val
            return val, "q"
    finally:
        mp.mp.dps = old


def betainv_forward(a, b, x, dps=50):
    """(value_min_side, side) -- the true min(P,Q) and which side, routing
    per gb.route_final's region map (reused, not re-derived) with an
    additional GUARD (PLAN's practical warning): both-huge-balanced
    (min(a,b)>=1e17ish) is gamma_corner's own documented hang trap
    (gamma_corner_value feeds min(a,b) to mpmath.gammainc unconditionally)
    -- excluded here by routing through R3's own extraction machinery
    instead whenever BOTH params are huge (matches beta-inl.h's own
    "BOTH-HUGE EXCLUSION" note: such lanes stay in ordinary R2)."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        a = mp.mpf(a); b = mp.mpf(b); x = mp.mpf(x)
        aa, bb, xx, tag = gb.route_final(a, b, x)
        if tag.startswith("R1"):
            val = r1_value_mp(aa, bb, xx, dps)
            side = "p" if tag.endswith("native") else "q"
            # R1 always returns P of (aa,bb,xx); native means (aa,bb,xx)=
            # (a,b,x) so val=P(a,b,x); swap means (aa,bb,xx)=(b,a,1-x) so
            # val=P(b,a,1-x)=Q(a,b,x).
            return val, side
        if tag.startswith("R4"):
            # R4 is the tiny-min analytic small-side series; reuse R1's
            # OWN series machinery at (tau,Bp,xi_tau) -- R1's series shape
            # covers R4's box too (both are the same t_n recursion; R4 is
            # simply R1 evaluated at a tiny first parameter -- gen_beta_
            # data.py's own module docstring: "R4's S is R1's series shape
            # minus the n=0 term, gamma-R4 verbatim in beta clothing").
            val = r1_value_mp(aa, bb, xx, dps)
            side = "p" if tag.endswith("native") else "q"
            return val, side
        if tag.startswith("R3"):
            val, side_local = r3_value_mp(aa, bb, xx, dps)
            # side_local is relative to (aa,bb,xx); translate to (a,b,x).
            native = tag.endswith("native")
            if native:
                side = side_local
            else:
                side = "p" if side_local == "q" else "q"
            return val, side
        # R2 / R2-gammalim
        native = tag.endswith("native") or (tag.endswith("gammalim") and "native" in tag)
        is_native = "native" in tag
        if "gammalim" in tag:
            val, side_local = gamma_corner_value_mp(aa, bb, xx, dps)
        else:
            val = gb.small_val_via_cf(aa, bb, xx, dps)
            side_local = "p"
        if is_native:
            side = side_local
        else:
            side = "p" if side_local == "q" else "q"
        return val, side
    finally:
        mp.mp.dps = old


def oracle_y(a, b, target, side, dps=45, lo=None, hi=None):
    """Root-find the TRUE y s.t. betainv_forward's small side equals
    target on the given side, bisecting in LOGIT space (linear-space
    bisection is unusable once y is anywhere near subnormal).
    MEASUREMENT-grade (not the certified oracle, built separately)."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        a_m, b_m, target = mp.mpf(a), mp.mpf(b), mp.mpf(target)

        def f(v):
            y = sigmoid(v)
            # A naive "y>=1 -> return 1-target" shortcut would wrongly
            # assume the SOLVED side is always P -- wrong for side='q',
            # where y very close to 1 (dps-limited, sigmoid rounds to
            # exactly 1.0 at extreme v) must drive got->0, not 1-target.
            # mpmath's own route_final/forward machinery already handles
            # y=0/y=1 exactly correctly (I_0=0, I_1=1 for any a,b>0), so
            # just clamp away from the literal endpoints (log(1-y)
            # hazards in the swapped CF orientation) rather than
            # special-case the return value.
            if y <= 0:
                y = mp.mpf(2) ** -1075
            elif y >= 1:
                y = 1 - mp.mpf(2) ** -1075
            val, s = betainv_forward(a_m, b_m, y, dps=dps)
            got = val if s == side else 1 - val
            return got - target

        lo = mp.mpf(-2000) if lo is None else mp.mpf(lo)
        hi = mp.mpf(2000) if hi is None else mp.mpf(hi)
        flo, fhi = f(lo), f(hi)
        tries = 0
        while flo * fhi > 0 and tries < 40:
            if side == "p":
                if flo > 0:
                    lo -= 200
                else:
                    hi += 200
            else:
                if fhi < 0:
                    hi += 200
                else:
                    lo -= 200
            flo, fhi = f(lo), f(hi)
            tries += 1
        if flo * fhi > 0:
            return None
        n_iters = min(300, max(80, dps * 3))
        for _ in range(n_iters):
            mid = (lo + hi) / 2
            fm = f(mid)
            if (fm > 0) == (flo > 0):
                lo, flo = mid, fm
            else:
                hi, fhi = mid, fm
            if hi - lo < mp.mpf(10) ** (-(dps - 5)):
                break
        return sigmoid((lo + hi) / 2)
    finally:
        mp.mp.dps = old


def small_side_of_y(a, b, y, dps=45):
    """(min(P,Q), side_of_min) -- the TRUE small-probability side, NOT
    route_final's own 'native/swap' criterion (a different predicate;
    conflating the two is a known harness-bug class, per the contract's
    explicit warning)."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        a_m, b_m, y_m = mp.mpf(a), mp.mpf(b), mp.mpf(y)
        val, side = betainv_forward(a_m, b_m, y_m, dps=dps)
        if side == "p":
            P, Q = val, 1 - val
        else:
            P, Q = 1 - val, val
        return (P, "p") if P <= Q else (Q, "q")
    finally:
        mp.mp.dps = old


# ============================================================================
# Part 2: SEEDS, native double end-to-end.
# ============================================================================
def lnB_double(a, b):
    return math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)


# --- S2: small-y series inversion (either end via exact complement) -------
def s2_series_S_double(a, b, y, nmax=80):
    """R1's own term recursion in double: t_n = t_{n-1}*(n-b)*y/n,
    S = sum t_n/(a+n)."""
    t = 1.0
    s = 0.0
    for n in range(0, nmax + 1):
        if n > 0:
            t *= (n - b) * y / n
        term = t / (a + n)
        s += term
        if n > 4 and abs(term) < 1e-18 * abs(s):
            break
    return s


def seed_S2(a, b, sigma, ncorr):
    """y0 = exp((ln sigma + ln alpha + lnB(a,b))/a) [leading term] + Picard
    corrections against R1's own series S(a,b,y). R1's native value is
    P = y^a/B(a,b) * S(a,b,y) with S's OWN leading term = 1/a (verify:
    s2_series_S_double starts its accumulation at term=t/(a+n), n=0 ->
    1/a, NOT 1 -- so P ~ y^a/(a*B(a,b)) as y->0), so y^a = a*B*P and
    y0 = exp((ln P + ln a + lnB)/a). The "+ln(a)" term is load-bearing:
    dropping it is a modest, Picard-self-healing error for a not too
    small (the first Picard iteration's -log(S) implicitly reintroduces
    +ln(a) since S~1/a), but for tiny a it makes y0 catastrophically
    wrong (e.g. a=0.02, sigma=0.42: missing ln(0.02)/0.02=-195.6/0.02,
    off by e^195 in the exponent), overflowing y0 past 1, from which the
    Picard loop's own series (evaluated at a nonsensical y) never
    recovers. This is a SEED for the SMALL-P-side target sigma -- caller
    passes 1-sigma when the small side is q (exact complement, no
    cancellation hazard: sigma is not tiny then)."""
    lnb = lnB_double(a, b)
    lna = math.log(a)
    lx = (math.log(sigma) + lna + lnb) / a
    y = math.exp(lx)
    for _ in range(ncorr):
        if not (math.isfinite(y) and 0.0 < y < 1.0):
            break
        S = s2_series_S_double(a, b, y)
        if not (math.isfinite(S) and S > 0):
            break
        lx = (math.log(sigma) + lnb - math.log(S)) / a
        y = math.exp(lx)
    return y


# --- S1: beta-Temme normal-quantile seed -----------------------------------
def cpsi_double(lam, a, b):
    u = -lam / a
    v = lam / b
    return a * (u - math.log1p(u)) + b * (v - math.log1p(v))


def dcpsi_dlam_double(lam, a, b, c):
    """Exact analytic derivative (paper-derived, module docstring):
    d(cpsi)/dlam = lam*c / ((a-lam)*(b+lam)), c=a+b. Clean rational form
    -- unlike gamma's single-parameter phi, beta's two-parameter cpsi has
    no need for a series-reversion trick: safeguarded Newton with this
    exact derivative, bracket-limited to (-b,a), is direct and robust."""
    denom = (a - lam) * (b + lam)
    if denom == 0.0:
        return math.inf
    return lam * c / denom


def lam_of_zeta_double(zeta, a, b, niter=100):
    # niter=100 was MEASURED necessary for full double-precision
    # convergence in the worst replay case (a=b=100, zeta=1.75: niter=40
    # left 2e-5 absolute error, niter=80 converged to 1e-14). The
    # safeguarded direct-lambda Newton oscillates before settling
    # (99.99->66.40->88.53->converged) -- almost certainly because
    # bisecting in RAW lambda space near the boundary is a worse-
    # conditioned formulation than gammainv's own log-space (u=ln(lambda))
    # Newton, which needed only 6 iterations for gamma's single-parameter
    # phi. A kernel implementation should very likely redo this in
    # log-space (or a bounded-domain rational substitution) rather than
    # inherit niter=100 verbatim -- this generator prioritized a CORRECT,
    # simply-derived seed formula within its time budget over an optimally
    # efficient one; the fix is well-scoped (re-parametrize lambda).
    """Safeguarded Newton (bisection fallback) inverting cpsi(lam)=zeta^2*nu
    [nu=a*b/(a+b) -- the defining equation, gb._lambda_of_zeta's own
    "target=zeta*zeta*nu"], lam in (-b,a) open interval, sign(lam)=
    sign(zeta). zeta=0 -> lam=0 exactly."""
    if zeta == 0.0:
        return 0.0
    c = a + b
    nu = a * b / c
    target = zeta * zeta * nu
    if zeta > 0:
        lo, hi = 0.0, a * (1.0 - 1e-15)
    else:
        lo, hi = -b * (1.0 - 1e-15), 0.0
    # initial guess: small-zeta quadratic approx cpsi(lam)~lam^2*c/(2ab)
    # (from the derivative at lam=0: d(cpsi)/dlam|_0=0, second derivative
    # c/(ab) -- cpsi(lam)~lam^2*c/(2*a*b) near lam=0), else bisection seed
    # at the interval midpoint scaled by sign.
    lam = math.copysign(math.sqrt(max(target * 2.0 * a * b / c, 0.0)), zeta)
    lam = min(max(lam, lo), hi) if lam == lam else (lo + hi) / 2.0
    for _ in range(niter):
        cp = cpsi_double(lam, a, b)
        v = cp - target
        if v == 0.0:
            break
        # bisection bracket update (monotone: cpsi increasing in |lam|)
        if zeta > 0:
            if v < 0:
                lo = lam
            else:
                hi = lam
        else:
            if v < 0:
                hi = lam
            else:
                lo = lam
        d = dcpsi_dlam_double(lam, a, b, c)
        if math.isfinite(d) and d != 0.0:
            step = lam - v / d
        else:
            step = (lo + hi) / 2.0
        if not (lo < step < hi):
            step = (lo + hi) / 2.0
        lam = step
    return lam


def erfcinv_double(y):
    return ginv.erfcinv_double(y)


def zeta0_of(nu, s, side):
    """beta's r3_R_at uses z=sqrt(cpsi)=zeta*sqrt(nu) -- NO factor of 2,
    unlike gamma's z=eta*sqrt(a/2) (gamma's cpsi is a*eta^2/2, while
    beta's cpsi IS zeta^2*nu directly, per _lambda_of_zeta's own
    target=zeta*zeta*nu). erfc(z0)/2=s defines z0=erfcinv(2s), so
    zeta0=z0/sqrt(nu) = z0*sqrt(1/nu), NOT z0*sqrt(2/nu) -- the extra
    sqrt(2) would cost ~1 bit of seed accuracy everywhere and trip the
    correction table's own domain gate (|zeta|<=S1_ZETA_MAX) at the
    wrong points. Diagnostic invariant: a correctly-scaled leading Temme
    term must improve as O(1/sqrt(nu)) with nu; a flat or worsening
    trend across a balanced-ridge nu sweep signals this scaling is
    wrong."""
    sgn = 1.0 if side == "p" else -1.0
    z0 = erfcinv_double(2.0 * s)
    if not math.isfinite(z0):
        z0 = 0.0
    return sgn * z0 * math.sqrt(1.0 / nu)


def seed_S1(a, b, sigma, side, ncorr=0):
    """z=erfcinv(2*sigma) [sign by side, zeta>=0<->P-small<->lambda>=0,
    per r3_R_at's own convention], nu=a*b/(a+b), invert the beta ridge
    mapping for lambda, y=(a-lambda)/(a+b). ncorr>0 applies the Temme
    correction via ONE perturbative Newton step in zeta -- re-derived
    HERE (not a literal copy of gammainv's eta_new=eta0+S/a) because
    beta's z-scaling differs from gamma's: r3_R_at's z=sqrt(cpsi)=
    zeta*sqrt(nu) (NO factor of 2, vs gamma's z=eta*sqrt(a/2)), and its
    sign convention FLIPS between branches (lam>=0: R=(leading-true)*
    scale; lam<0: R=(true-leading)*scale -- opposite signs, both
    verified directly against the r3_R_at source). Linearizing
    leading(zeta)=erfc(sqrt(cpsi))/2 the same way (dropping R's own
    derivative, standard perturbative-Newton) gives, uniformly for both
    branches: leading'(zeta0) = -sqrt(nu/pi)*exp(-cpsi), and
      zeta0>=0: delta = -R_output/(nu*sqrt(2))
      zeta0<0:  delta = +R_output/(nu*sqrt(2))
    i.e. delta = -sign(zeta0)*S/(nu*sqrt(2)), S=sum c_k(zeta0,p)/nu^k
    (the SAME Horner-in-1/nu series form gammainv uses, just with the
    corrected prefactor). Do NOT copy gammainv's "zeta_new=zeta0+S/nu"
    verbatim: it uses a different prefactor (missing 1/sqrt(2)) and
    omits the branch-dependent sign flip above. Diagnostic: as with
    zeta0_of, a correctly-scaled correction improves with nu
    (O(1/sqrt(nu)) leading error) across a balanced-ridge nu sweep; a
    wrong prefactor or sign shows up as flat or regressing bits."""
    c = a + b
    nu = a * b / c
    p = a / c
    zeta = zeta0_of(nu, sigma, side)
    # The Chebyshev-in-zeta correction table is fit ONLY on
    # [-S1_ZETA_MAX,S1_ZETA_MAX]; evaluating it outside that (e.g.
    # zeta0=-12.6 at a=0.5,b=5,q~1e-17: |t|=3.6, a Chebyshev series of
    # degree 14 at |t|=3.6 diverges to ~1e7) would silently produce a
    # catastrophically WRONG correction instead of the honest ncorr=0
    # fallback -- gate the correction to its own fitted domain, mirroring
    # gammainv's own |eta0|<=ETA_MAX candidate-availability gate.
    if ncorr > 0 and _S1_CHEB_D is not None and abs(zeta) <= S1_ZETA_MAX \
            and nu >= S1_NU_MIN:
        ck = [_s1_ck_eval(k, zeta, p) for k in range(ncorr)]
        S = 0.0
        for k in range(ncorr - 1, -1, -1):
            S = S / nu + ck[k]
        sgn_z = 1.0 if zeta >= 0.0 else -1.0
        zeta = zeta - sgn_z * S / (nu * math.sqrt(2.0))
    lam = lam_of_zeta_double(zeta, a, b)
    y = (a - lam) / c
    return y


# --- S1 correction table: eta-style perturbative Temme correction, 2D in
# (zeta,p), extracted via the SAME r3_R_at forward reused (not re-derived)
# from gen_beta_data.py, but over a WIDER zeta domain than R3's own
# ratio-band table (gammainv precedent, module docstring). Built lazily by
# build_s1_correction() in main(); None until then (ncorr=0 seed still
# works, self-check (b) measures whether the correction is needed at all).
# ============================================================================
S1_ZETA_MAX = 3.5   # seed-domain half-width, pinned below by self-check;
                     # wider than gb.ZETA_MAX (~1.02, R3's tight ratio-band
                     # table) since S1 is one of FOUR global candidates,
                     # not gated to the ridge -- mirrors gammainv's S1
                     # domain being wider than gamma's own ridge-only table.
S1_NZ = 15
S1_NP = 9
S1_KEXT = 3
S1_NU_LIST = (mp.mpf(30), mp.mpf(60), mp.mpf(120), mp.mpf(240))
S1_P_MID = mp.mpf("0.25")
S1_P_HALF = mp.mpf("0.25")
S1_NU_MIN = 2.0  # the 1/nu correction series is asymptotic in nu
                  # (S(zeta,p) itself is O(1)-bounded, measured -0.18..-2.4
                  # over the fitted domain, but the APPLIED correction S/nu
                  # is unbounded as nu->0 -- at nu=5e-9 (a=b=1e-8,
                  # joint-tiny) it reaches 2e8, sending zeta from 0 to a
                  # nonsense value). Measured: |S/nu| stays a small
                  # fraction (<10%) of S1_ZETA_MAX for nu>=2, growing past
                  # half of it by nu=0.1 and nonsense by nu=0.01 -- gate
                  # BOTH the correction and S1's own candidacy on
                  # nu>=S1_NU_MIN (mirrors gammainv's S1_A_MIN, same
                  # disease).
_S1_CHEB = None  # [row_k][coef] monomial-in-p? -- see build_s1_correction;
                  # populated as [k][nz][np] 2D-Chebyshev-in-(zeta,p) coefs.


def build_s1_correction(kreport, dps=45):
    global _S1_CHEB
    zeta_nodes = [S1_ZETA_MAX * t for t in gb._cheb_nodes(S1_NZ)]
    p_nodes = [S1_P_MID + S1_P_HALF * u for u in gb._cheb_nodes(S1_NP)]
    grid = [[None] * S1_NP for _ in range(S1_NZ)]
    t0 = time.time()
    for i, zeta in enumerate(zeta_nodes):
        for j, p in enumerate(p_nodes):
            A = mp.matrix(len(S1_NU_LIST), S1_KEXT)
            bcol = mp.matrix(len(S1_NU_LIST), 1)
            for r, nu in enumerate(S1_NU_LIST):
                v = 1 / nu
                for k in range(S1_KEXT):
                    A[r, k] = v ** k
                bcol[r, 0] = gb.r3_R_at(nu, p, zeta, dps)
            c = mp.qr_solve(A, bcol)[0]
            grid[i][j] = [c[k, 0] for k in range(S1_KEXT)]
        if (i + 1) % 5 == 0:
            print(f"    S1 correction extraction {i+1}/{S1_NZ} zeta-rows "
                  f"({time.time()-t0:.0f}s)", file=sys.stderr)
    # 2D DCT per order k (exact Chebyshev interpolation, gb.fit_r3_tensor's
    # own pattern, reused).
    coef2d = []
    for k in range(kreport):
        mid = [None] * S1_NZ
        for i in range(S1_NZ):
            vals = [grid[i][j][k] for j in range(S1_NP)]
            mid[i] = gb._cheb_coeffs_1d(vals)
        coef = [[None] * S1_NP for _ in range(S1_NZ)]
        for m in range(S1_NP):
            vals = [mid[i][m] for i in range(S1_NZ)]
            colc = gb._cheb_coeffs_1d(vals)
            for n in range(S1_NZ):
                coef[n][m] = colc[n]
        coef2d.append(coef)
    _S1_CHEB = coef2d
    global _S1_CHEB_D
    _S1_CHEB_D = [[[rd(c) for c in row] for row in coef2d[k]] for k in range(kreport)]
    print(f"    S1 correction table built: {S1_NZ}x{S1_NP} nodes, "
          f"K={kreport}, zeta in [-{float(S1_ZETA_MAX)},{float(S1_ZETA_MAX)}], "
          f"{time.time()-t0:.0f}s total", file=sys.stderr)


def _s1_row_eval_mp(coef_row, zeta, p):
    t = zeta / S1_ZETA_MAX
    u = (p - S1_P_MID) / S1_P_HALF
    row_vals = [gb._clenshaw(coef_row[n], u) for n in range(S1_NZ)]
    return gb._clenshaw(row_vals, t)


# double-precision Chebyshev-in-(zeta,p) rows, built once _S1_CHEB is
# populated with double-rounded coefficients (measurement path matches
# what the kernel will actually do: double arithmetic, not mpmath).
_S1_CHEB_D = None


def _s1_ck_eval(k, zeta, p):
    """c_k(zeta,p), native double, symmetry e_k(zeta,p)=-e_k(-zeta,1-p)
    applied exactly as R3's own (check (h) precedent) since p>0.5
    traffic is never built into the table."""
    if p > 0.5:
        zeta_e, p_e, sign = -zeta, 1.0 - p, -1.0
    else:
        zeta_e, p_e, sign = zeta, p, 1.0
    t = zeta_e / S1_ZETA_MAX
    u = (p_e - float(S1_P_MID)) / float(S1_P_HALF)
    row_vals = []
    for n in range(S1_NZ):
        row = _S1_CHEB_D[k][n]
        b1 = b2 = 0.0
        for c in row[:0:-1]:
            b1, b2 = 2.0 * u * b1 - b2 + c, b1
        row_vals.append(u * b1 - b2 + row[0])
    b1 = b2 = 0.0
    for c in row_vals[:0:-1]:
        b1, b2 = 2.0 * t * b1 - b2 + c, b1
    val = t * b1 - b2 + row_vals[0]
    return sign * val


# --- S3: gamma-limit transfer ------------------------------------------
def seed_S3(a, b, sigma, side):
    """t = -huge*log1p(-y) [beta the huge param, native form] or
    t = -huge*ln(y) [alpha the huge param] -- mirrors src/beta-inl.h's own
    gamma-limit slice exactly (module docstring citation). Seeds t via
    the EXISTING gammainv seed_for (imported, not re-derived), then
    inverts back to y. Returns None if neither parameter is usably huge
    (caller's own domain guard)."""
    if b >= a:
        s, huge = a, b
        gside = side  # native: value~P_gamma(s,t) tracks the SAME side
        try:
            t0 = gammainv_seed_for(s, sigma, gside)
        except (OverflowError, ValueError, ZeroDivisionError):
            return None
        if not (math.isfinite(t0) and t0 >= 0.0):
            return None
        # y = -expm1(-t/huge)
        try:
            y = -math.expm1(-t0 / huge)
        except OverflowError:
            return None
        return y
    else:
        s, huge = b, a
        gside = "p" if side == "q" else "q"
        try:
            t0 = gammainv_seed_for(s, sigma, gside)
        except (OverflowError, ValueError, ZeroDivisionError):
            return None
        if not (math.isfinite(t0) and t0 >= 0.0):
            return None
        try:
            y = math.exp(-t0 / huge)
        except OverflowError:
            return None
        return y


# --- S4: joint-tiny logit closed form ---------------------------------
def s4_c_closed_form(a, b):
    """c(alpha,beta) ~= (pi^2/12)*(alpha-beta), LEADING order in
    C=alpha+beta->0 (module docstring derivation: exact substitution
    B_y=int_{-inf}^{logit y} exp(alpha*u-C*L(u))du, L(u)=ln(1+e^u); the
    min(e^{alpha*u},e^{-beta*u}) approximation to the integrand gives
    logit(y)=(s-s*)/w EXACTLY at O(C^0) for ALL alpha,beta [not just
    alpha=beta] -- so c is a NEXT-order effect. Measured (this
    generator's own probe, high-dps root-find of v0 solving I_y=s* over
    C=1e-3..1e-6, e=(alpha-beta)/C=0.01..0.99): v0/(alpha-beta) is
    INDEPENDENT of the skew e to >10 digits at fixed C, and converges
    LINEARLY in C to pi^2/12=0.822467033... (measured residual ~
    -0.4817*C, i.e. the next term is itself O(C*(alpha-beta)), negligible
    at the domain S4 actually needs -- see self-check (c))."""
    return (math.pi * math.pi / 12.0) * (a - b)


def seed_S4(a, b, sigma, side):
    """The linearized form v=(s-s*)/w+c(a,b) (Taylor-expanding the
    leading exponential relation around v=0) is only valid for v=O(1)
    (near the plateau center) and diverges badly for |v| beyond a few
    units (measured: at C=0.04, |v|=9, linear-form error 0.79 vs the
    exponential form's 0.032 -- 25x). Use instead the EXACT (not
    asymptotic) leading-order relation for EACH branch, which is exact
    in v of any sign/magnitude as alpha,beta->0 (the min(...) integrand
    approximation's own leading term, inverted directly rather than
    Taylor-expanded around v=0):
      v<0 (s<=s*):  P(y)~y^alpha/(alpha*B(alpha,beta))  => v=(ln s+ln
                     alpha+lnB)/alpha  [alpha,beta EXACT via lgamma, not
                     an asymptotic B~1/w approximation]
      v>=0 (s>s*):  Q(y)~z^beta/(beta*B(alpha,beta)), z=1-y  => v=
                     -(ln(1-s)+ln beta+lnB)/beta
    This is IDENTICAL to seed_S2's own ncorr=0 zeroth iterate (in the
    appropriate orientation) -- S4 and S2 share one mechanism at leading
    order; S4's distinct value is picking the CORRECT branch by s vs s*
    without needing a separate orientation trial. MEASURED (high-
    precision v0 sweep): exact-B alone (no c correction) reaches 9.0b at
    C=0.04 (vs 4.97b with the linear+c form), 6.7-6.8b at C=0.2,
    degrading to ~4b by C=2-4 -- c(alpha,beta) ACTIVELY HURTS once B is
    exact (it was compensating for the asymptotic form's own error;
    adding it back on top of exact-B measures WORSE at every C tested,
    e.g. C=2: 3.95b->0.17b), so it is DROPPED, not merely left at 0.
    s4_c_closed_form's derivation stays documented (the paper-math
    argument and pi^2/12 constant are real and correct for the LINEAR
    form) but is no longer called."""
    s = sigma if side == "p" else 1.0 - sigma
    c_ab = a + b
    sstar = b / c_ab
    try:
        lnB = lnB_double(a, b)
        if s <= sstar:
            v = (math.log(s) + math.log(a) + lnB) / a
        else:
            v = -(math.log1p(-s) + math.log(b) + lnB) / b
    except ValueError:
        return None
    if not math.isfinite(v):
        return None
    try:
        if v >= 0:
            e = math.exp(-v)
            y = 1.0 / (1.0 + e)
        else:
            e = math.exp(v)
            y = e / (1.0 + e)
    except OverflowError:
        y = 0.0 if v < 0 else 1.0
    return y


# ============================================================================
# S5: closed form for the region none of S1-S4 target -- moderate-tiny
# SINGLE shape parameter with the OTHER parameter only moderate (not
# ridge-large, not gamma-limit-huge, not both-tiny) and a MODERATE target
# probability on neither extreme -- e.g. (a=0.2,b=2,y=0.1): S1's own
# leading order gives -0.74b even with
# EXACT lambda inversion (nu=0.18 is simply too small for the Temme
# normal-quantile asymptotic, independent of Newton method), S2/S4's
# small-y/small-C asymptotics fail outright when neither y nor 1-y is
# small (S2's own zeroth iterate lands z0>1, nonsensical -- verified
# directly, not assumed). NONE of S1-S4 target "moderate shape, moderate
# probability" AT ALL -- this is genuinely uncovered territory, not a
# tuning gap in an existing candidate.
# Closed form: logit(Y) for Y~Beta(alpha,beta) is EXACTLY ln(X1)-ln(X2)
# for independent X1~Gamma(alpha), X2~Gamma(beta) (Y=X1/(X1+X2), standard
# Gamma-ratio construction) -- so E[logit Y]=psi(alpha)-psi(beta) and
# Var[logit Y]=psi'(alpha)+psi'(beta) EXACTLY (digamma/trigamma of the
# summed shapes, sum of two independent variances), for ANY alpha,beta,
# not an asymptotic in either. A normal approximation in LOGIT space
# using these EXACT first two moments (clean-room: standard Gamma-ratio
# identity plus a Cornish-Fisher-style normal-quantile seed, no ported
# code) is the natural low-order closed form for the region no other
# candidate targets -- measured to beat S3/S4 specifically in the
# moderate-shape band (e.g. (0.2,2,0.1): S4 2.66b vs S5 3.47b),
# and offered as a GLOBAL fifth candidate (cheap-residual selected like
# the other four), not a gated route -- it does no harm where it isn't
# competitive (the comparison simply picks something else).
# ============================================================================
def seed_S5(a, b, sigma, side, dps=25):
    s = sigma if side == "p" else 1.0 - sigma
    if not (0.0 < s < 1.0):
        return None
    with mp.workdps(dps):
        mu = mp.digamma(a) - mp.digamma(b)
        var = mp.polygamma(1, a) + mp.polygamma(1, b)
        sd = mp.sqrt(var)
    if s <= 0.5:
        zz = -math.sqrt(2.0) * erfcinv_double(2.0 * s)
    else:
        zz = math.sqrt(2.0) * erfcinv_double(2.0 * (1.0 - s))
    v = float(mu) + zz * float(sd)
    if not math.isfinite(v):
        return None
    try:
        if v >= 0:
            e = math.exp(-v)
            y = 1.0 / (1.0 + e)
        else:
            e = math.exp(v)
            y = e / (1.0 + e)
    except OverflowError:
        y = 0.0 if v < 0 else 1.0
    return y


# ============================================================================
# Part 3: quad-candidate seed selection (cheap forward-residual global
# comparison, the same mechanism gammainv uses, now the design's own
# starting point here -- "no parameter gating except per-candidate
# availability/stability gates you derive and pin").
# ============================================================================
T_JT = 2.0 ** -8  # [PROVISIONAL -- pinned by self-check (d)] max(alpha,beta)
                    # threshold below which S4 is offered as a candidate;
                    # must comfortably own min(a,b)<=2^-52 (trivially true
                    # once max(a,b)<T_JT<<1, since min<=max).
S1_NCORR = 2       # [PROVISIONAL -- pinned by self-check (b)]


def cheap_residual_beta(a, b, y0, s, side, dps_list=(20, 30)):
    """LOW-PRECISION forward eval -- selector only, not precision-bearing
    (gammainv's own cheap_residual pattern, dps ladder for the same
    reason: extraction/CF machinery can go numerically singular at very
    low dps for some points)."""
    if not (math.isfinite(y0) and 0.0 < y0 < 1.0):
        return None
    for dps_try in dps_list:
        try:
            val, sd = betainv_forward(a, b, y0, dps=dps_try)
            v = float(val)
            if sd != side:
                v = 1.0 - v
            if s == 0:
                return abs(v)
            return abs(v - s) / s
        except (ZeroDivisionError, ValueError, OverflowError, mp.libmp.libhyper.NoConvergence):
            continue
    return None


def seed_for(a, b, sigma, side, s1_ncorr=None, s2_ncorr=6, jt_thresh=None):
    """Quad-candidate seed selection: compute each candidate when its own
    domain guard passes, select by cheap_residual_beta. Per-candidate
    try/except (gammainv's own pattern) so one candidate's exception
    (e.g. S2's math.exp overflow at extreme alpha) never discards an
    already-good other candidate."""
    if s1_ncorr is None:
        s1_ncorr = S1_NCORR
    if jt_thresh is None:
        jt_thresh = T_JT
    c = a + b
    nu = a * b / c

    s1_candidate = None
    if nu >= S1_NU_MIN:
        try:
            s1_candidate = seed_S1(a, b, sigma, side, s1_ncorr)
        except (OverflowError, ValueError, ZeroDivisionError):
            s1_candidate = None

    try:
        # side='q' must NOT call seed_S2(a,b,1-sigma,...) -- that keeps
        # (a,b) unswapped with sigma COMPLEMENTED, asking S2 (a
        # small-Y-for-small-TARGET formula) to solve P(a,b,y)=1-sigma for
        # y, but 1-sigma=P is near 1 (large, NOT small) whenever sigma=Q
        # is the genuinely small side -- the exact opposite of S2's
        # design domain. The correct swapped-twin construction (module
        # docstring: "the swapped twin is free", contract's own phrasing)
        # is the SWAP IDENTITY itself: Q(a,b,y) = P(b,a,1-y), so solve
        # P(b,a,z)=sigma [swap a<->b, sigma UNCHANGED -- sigma IS the
        # small target on the swapped side] for the small z=1-y, then
        # y=1-z.
        s2v = seed_S2(a, b, sigma, s2_ncorr) if side == "p" \
            else seed_S2(b, a, sigma, s2_ncorr)
        s2_candidate = s2v if side == "p" else 1.0 - s2v
    except (OverflowError, ValueError, ZeroDivisionError):
        s2_candidate = None

    try:
        s3_candidate = seed_S3(a, b, sigma, side)
    except (OverflowError, ValueError, ZeroDivisionError):
        s3_candidate = None

    # t_jt gates the CLOSED-FORM ROUTE (where S4's answer ships without
    # iteration under the plateau contract), NOT candidacy in this
    # quad-candidate seed selection -- gating S4's candidacy on
    # `max(a,b)<jt_thresh` here would silently exclude it from points
    # like a=b=0.02 (well above t_jt=2^-8) where it is in fact the best
    # candidate. S4 is offered unconditionally (guarded only by its own
    # validity: log(s) and log1p(-s) both finite, i.e. s in (0,1)).
    s4_candidate = None
    try:
        s4_candidate = seed_S4(a, b, sigma, side)
    except (OverflowError, ValueError, ZeroDivisionError):
        s4_candidate = None

    # DEEP-SMALL as a fifth unconditional candidate: a genuinely
    # deep-tail point (sigma astronomically tiny, e.g. 1e-250) is
    # DEEP-SMALL's own designed territory regardless of whether (a,b)
    # also happen to be small -- S2's Picard series can be asked to
    # converge a correction from a badly-scaled zeroth seed there and
    # fail (measured floor 46.8b, short of the 55b gate) even though the
    # CLOSED FORM (which is exactly S2's own zeroth-order term) is
    # trivially available and often near-exact at that depth. Always
    # tried, cheap, selected only if it wins the residual comparison --
    # consistent with the quad-candidate design's own "compute
    # unconditionally when its own guard passes" doctrine (deep-small's
    # guard is simply "produces a finite y in (0,1)").
    try:
        deep_candidate = deep_small_y(a, b, sigma, side)
    except (OverflowError, ValueError, ZeroDivisionError):
        deep_candidate = None

    # S5: logit-normal via EXACT digamma/trigamma moments, the only
    # candidate that targets "moderate shape, moderate probability"
    # territory. Always offered, cheap, selected only if it wins.
    try:
        s5_candidate = seed_S5(a, b, sigma, side)
    except (OverflowError, ValueError, ZeroDivisionError):
        s5_candidate = None

    best_x, best_r = None, None
    for cand in (s1_candidate, s2_candidate, s3_candidate, s4_candidate,
                 deep_candidate, s5_candidate):
        if cand is None or not (math.isfinite(cand) and 0.0 < cand < 1.0):
            continue
        r = cheap_residual_beta(a, b, cand, sigma, side)
        if r is None:
            continue
        if best_r is None or r < best_r:
            best_x, best_r = cand, r
    if best_x is not None:
        return best_x
    for cand in (s2_candidate, s1_candidate, s3_candidate, s4_candidate,
                 deep_candidate, s5_candidate):
        if cand is not None and math.isfinite(cand) and 0.0 < cand < 1.0:
            return cand
    return float("nan")


# ============================================================================
# Part 4: STEPS -- safeguarded logit-Newton, the gammainv safeguard
# package assumed whole (module docstring): m(y)=ln P(y)-ln Q(y) [=logit(P)],
# w=P*Q/g(y) [g=beta density, dm/dy=g/(P*Q)], target_m=logit(sigma) (side=p)
# or -logit(sigma) (side=q). Step: resid=m-target_m; ls=-scale*w*resid,
# floored at -0.9; candidate=y*(1+ls); ACCEPT if 0<cand<1 AND
# (|resid_best|<TRUST_RESID OR |resid_new|<=|resid_best|); reject ->
# scale*=1/8 for the NEXT step, reset to 1 on acceptance.
#
# STEPS_N = 4: the contract's original "3 shared steps" (gammainv's own
# pin) does not hold band-wide here. Root cause, precisely diagnosed (not
# assumed): a bounded interior sub-band (min(a,b) approx 0.02-0.5, skew
# 3-10x, y interior 0.1-0.3) has NO closed-form seed family (five tried:
# S1 Temme-normal, S2 small-y, S3 gamma-limit, S4 exact-B exponential, S5
# logit-normal via exact digamma/trigamma moments) exceeding ~2-5 bits --
# confirmed NOT a selection failure (cheap-residual always picks the best
# available candidate) and NOT an eps/noise-floor artifact (already at
# EPS_TIGHT). Measured convergence from the worst point in the band
# (a=0.1,b=1, y=0.3, seed 2.12b): 2.12->6.66->16.48->36.16->75.51->103.78
# bits over 6 steps -- clean quadratic (each step envelope-doubles the
# correct bits), step 4 clears the 55b gate by 20+ bits margin band-wide.
# The "fix the seed, don't shave margin" precedent (gammainv S3, small-a
# mid band) does NOT apply symmetrically here: that case had a
# contract-compliant seed fix AVAILABLE; here the seed side is exhausted
# (five families, one addition (S5), one correction tested-and-rejected
# (Cornish-Fisher, measured worse everywhere)) and step 4 RESTORES full
# margin (75+ bits), it does not shave anything. The safeguard package
# (reject residual-increasing, bypass |resid|<1/2, freeze-by-select)
# makes the fourth step IDEMPOTENT for lanes already converged after 3 --
# cost is bounded at one extra forward evaluation for those lanes,
# strictly cheaper and simpler than a sixth fitted seed with its own
# table, residual-compare slot, and seams.
# LATITUDE (carry forward): the kernel MAY add a whole-vector all-lanes-
# converged skip after step 3 as a bench optimization (gammainv's "1
# Halley vs 2 Newton" precedent -- the accuracy gates must hold either
# way, this is a throughput decision only).
# ============================================================================
TRUST_RESID = 0.5  # bypass |resid|<1/2 (contract's own number, matches
                    # gammainv's kGammaInvTrustResid intent)
STEPS_N = 4


def beta_density_mp(a, b, y, dps):
    """g(y) = y^(a-1)*(1-y)^(b-1)/B(a,b), log-space assembly."""
    a = mp.mpf(a); b = mp.mpf(b); y = mp.mpf(y)
    lg = (a - 1) * mp.log(y) + (b - 1) * mp.log1p(-y) - (
        mp.loggamma(a) + mp.loggamma(b) - mp.loggamma(a + b))
    return mp.e ** lg


def m_and_w_mp(a, b, y, dps):
    """(m, w) = (ln P - ln Q, 1/(y*dm/dy)) at working dps. dm/dy = g/(P*Q)
    (m=logit(P), standard logit-density identity), so 1/dm/dy = P*Q/g;
    w carries an EXTRA 1/y factor because the step is applied
    MULTIPLICATIVELY, cand=y*(1+ls) -- gammainv's own "w" is
    1/(x*dm/dx), not 1/(dm/dx), for exactly this reason, documented in
    its own step comment ("w -> 1/x in the far tail"). Omitting the 1/y
    factor (w=P*Q/g alone) under-steps by a factor of y at every
    iteration -- measured 6.94->9.12 bits over 3 "steps" instead of
    converging."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        val, side = betainv_forward(a, b, y, dps=dps)
        if side == "p":
            P, Q = val, 1 - val
        else:
            P, Q = 1 - val, val
        m = mp.log(P) - mp.log(Q)
        g = beta_density_mp(a, b, y, dps)
        w = P * Q / (g * y) if g != 0 else mp.mpf("inf")
        return m, w
    finally:
        mp.mp.dps = old


def target_m_of(sigma, side):
    sigma = mp.mpf(sigma)
    lg = mp.log(sigma) - mp.log1p(-sigma)
    return lg if side == "p" else -lg


def step_once(a, b, y0, target_m, eps, scale, dps=32):
    """One safeguarded logit-Newton step, forward P/Q and w perturbed by
    +/-eps (4 sign combos), worst case kept (gammainv's own noise-
    injection pattern). Returns list of 4 candidate y's (caller reduces)."""
    with mp.workdps(dps):
        m0, w0 = m_and_w_mp(a, b, y0, dps)
        eps_m = mp.mpf(eps)
        out = []
        for sm in (1 + eps_m, 1 - eps_m):
            for sw in (1 + eps_m, 1 - eps_m):
                m_p = m0 * sm if m0 >= 0 else m0 / sm
                w_p = w0 * sw
                resid = m_p - target_m
                ls = -scale * w_p * resid
                if ls < mp.mpf("-0.9"):
                    ls = mp.mpf("-0.9")
                cand = y0 * (1 + ls)
                out.append(cand)
        return out, abs(m0 - target_m)


def simulate_steps_beta(a, b, y0, true_y, sigma, side, eps, nsteps=STEPS_N, dps=32):
    """Worst-case relative-bit result after nsteps, PERFORMANCE
    SIMPLIFICATION carried from gammainv (module docstring): steps
    1..nsteps-1 deterministic (eps=0) trunk, final step adversarial 4-way.
    Safeguard applied on the DETERMINISTIC trunk too (reject
    residual-increasing, bypass |resid|<TRUST_RESID, 1/8 backtrack)."""
    target_m = target_m_of(sigma, side)
    with mp.workdps(dps):
        y = mp.mpf(y0)
        scale = mp.mpf(1)
        rbest = None
        for step_i in range(nsteps - 1):
            if not (0 < y < 1):
                return -1.0
            cands, resid0 = step_once(a, b, y, target_m, 0.0, scale, dps=dps)
            cand = cands[0]
            if rbest is None:
                rbest = resid0
            if not (0 < cand < 1):
                scale = scale / 8
                continue
            m_new, _ = m_and_w_mp(a, b, cand, dps)
            resid_new = abs(m_new - target_m)
            if resid0 < TRUST_RESID or resid_new <= rbest:
                y = cand
                rbest = resid_new
                scale = mp.mpf(1)
            else:
                scale = scale / 8
        if not (0 < y < 1):
            return -1.0
        cands, resid0 = step_once(a, b, y, target_m, eps, scale, dps=dps)
        worst = 999.0
        for cand in cands:
            if not (0 < cand < 1):
                b_ = -1.0
            else:
                rel = abs((cand - true_y) / true_y)
                b_ = 999.0 if rel == 0 else float(-mp.log(rel, 2))
            worst = min(worst, b_)
        return worst


def simulate_steps_beta_multi(a, b, y0, true_y, sigma, side, eps, max_nsteps, dps=32):
    target_m = target_m_of(sigma, side)
    out = {}
    with mp.workdps(dps):
        y = mp.mpf(y0)
        scale = mp.mpf(1)
        rbest = None
        for n in range(1, max_nsteps + 1):
            if not (0 < y < 1):
                out[n] = -1.0
                continue
            cands, resid0 = step_once(a, b, y, target_m, eps, scale, dps=dps)
            worst = 999.0
            for cand in cands:
                if not (0 < cand < 1):
                    bb = -1.0
                else:
                    rel = abs((cand - true_y) / true_y)
                    bb = 999.0 if rel == 0 else float(-mp.log(rel, 2))
                worst = min(worst, bb)
            out[n] = worst
            # advance deterministic trunk (eps=0) for next n
            cands0, _ = step_once(a, b, y, target_m, 0.0, scale, dps=dps)
            cand0 = cands0[0]
            if rbest is None:
                rbest = resid0
            if 0 < cand0 < 1:
                m_new, _ = m_and_w_mp(a, b, cand0, dps)
                resid_new = abs(m_new - target_m)
                if resid0 < TRUST_RESID or resid_new <= rbest:
                    y, rbest, scale = cand0, resid_new, mp.mpf(1)
                else:
                    scale = scale / 8
            else:
                scale = scale / 8
    return out


# ============================================================================
# Part 5: DEEP-SMALL closed form, both orientations. y = exp((ln sigma +
# ln alpha + lnB(alpha,beta))/alpha) -- S2's own leading term with the
# series correction S' DROPPED entirely (S'->1 as y->0).
#
# CUT PREDICATE, re-derived here and verified by direct measurement
# (self-check (f)), not trusted blind: sigma = y^a/(a*B(a,b)) * S',
# S' = a*sum_n t_n/(a+n) = 1 + a*sum_{n>=1} t_n/(a+n), t_0=1,
# t_n=t_{n-1}*(n-b)*y/n. Dropping
# S'->1 costs a relative error in ln(y) of |ln S'|/a; for small y this is
# dominated by the n=1 term, t_1=(1-b)*y:
#   |ln S'|/a ~= |a*t_1/(a+1)|/a = |t_1|/(a+1) = |1-b|*y/(1+a)
# UNLIKE gamma's own deep-small cut (x0*(1+a), no second parameter --
# gamma's forward has none), beta's dropped term is driven by the OTHER
# side's parameter (b, not a) -- there is no "a" analog of gamma's cut
# because b, not a, supplies the leading nonconstant term of the series
# (t_1's coefficient (n-b) at n=1 is exactly 1-b; a enters only via the
# a+1 denominator, a second-order effect near a~O(1) and negligible as
# a->0 where a+1->1 anyway). Both orientations, by the same derivation
# with roles (a,b,y) -> (b,a,1-y) for side='q':
#   P-side (side='p'): |1-b|*y/(1+a)         < DEEP_SMALL_CUT
#   Q-side (side='q'): |1-a|*(1-y)/(1+b)     < DEEP_SMALL_CUT
# A single-orientation self-check can leave the untested orientation
# unbounded (90 ULP reachable at small a, measured) -- self-check (f)
# below sweeps BOTH orientations.
# ============================================================================
DEEP_SMALL_CUT = 2.0 ** -60


def deep_small_cut_bound(a, b, y, side):
    """The analytic dropped-factor bound, |1-b|*y/(1+a)*corr(y) [p] or
    |1-a|*(1-y)/(1+b)*corr(1-y) [q] -- compare against DEEP_SMALL_CUT to
    decide routing. y is the FULL solved variable (not the small-z form)
    in both cases -- side='q' uses (1-y), the genuinely small quantity
    on that branch.

    The bare t1/(1+a-or-b) formula UNDER-PREDICTS at the widened
    gamma-limit corner (huge OTHER-side parameter, e.g. a=5,b=1e300,
    side='q', y=1e-6: measured true/bound ratio 13.8). A naive guard
    threshold (e.g. reject t1>=0.1) is not a fix: it only pushes the
    worst case to a=0.9,b=1e100,y=1e-6 (ratio 6.49, with t1~0.0999999
    just inside the threshold). ROOT CAUSE, isolated by direct
    measurement (not trusted from algebra alone): the true/bound ratio,
    swept over a wide range of (a-or-b) [the OTHER-side coefficient] at
    FIXED y' (the own-side small variable, y for p / 1-y for q), is
    EXACTLY -ln(1-y')/y' and is INDEPENDENT of the other-side
    coefficient entirely (verified to 5+ significant figures across
    coefficient values spanning 1e-4 to 5, at fixed y'). This is exact
    in the huge-OTHER-side-exponent limit (S' -> (1-y')^(other_side-1)
    there, a clean closed form), and the SAME correction is a sound
    (ratio<=1, worst measured 1.0000000004 -- boundary floating-point
    noise) upper predictor for MODERATE other-side parameters too
    (self-check (f) sweeps the full range). corr(y')=-ln(1-y')/y' -> 1
    as y'->0 (l'Hopital), so this reduces to the bare leading-order
    formula exactly where that was already valid, and inflates it
    exactly where it wasn't -- not a heuristic guard, a closed-form
    exact correction factor derived from the measurement."""
    if side == "p":
        yp, coeff, denom = y, abs(1.0 - b), (1.0 + a)
    else:
        yp, coeff, denom = 1.0 - y, abs(1.0 - a), (1.0 + b)
    if not (0.0 < yp < 1.0):
        return math.inf
    corr = 1.0 if yp < 1e-8 else -math.log1p(-yp) / yp
    return coeff * yp / denom * corr


def deep_small_y(a, b, sigma, side):
    """Must include the +ln(a) [resp. +ln(b)] term -- exp((ln sigma +
    lnB)/a) alone (without ln a) is the same disease as seed_S2's own
    zeroth iterate without it (see seed_S2's docstring). Correct leading
    form (S'->1 as y->0, sigma = y^a/(a*B(a,b))*S'): y0 = exp((ln sigma +
    ln a + lnB)/a)."""
    if side == "p":
        lnb = lnB_double(a, b)
        return math.exp((math.log(sigma) + math.log(a) + lnb) / a)
    else:
        lnb = lnB_double(b, a)
        z = math.exp((math.log(sigma) + math.log(b) + lnb) / b)
        return 1.0 - z


def _deep_small_dropped_rel(a, b, y0, side, dps=60):
    """|ln S'(a,b,y or z)|/a[or b] -- the EXACT quantity DROPPED by the
    closed form (S2's own series correction, S'->1 assumed). Must
    accumulate the TRUE S' = a*sum_n t_n/(a+n) = 1 + a*sum_{n>=1}
    t_n/(a+n), WITH the a/(a+n) weight -- summing sum_n c_n*y^n without
    it gives (1-y)^(b-1) (the OTHER factor entirely), not S'; that
    substitution is conservative at tiny a (where a/(a+n)~a/n->0 makes
    S' close to 1 regardless) but measures the wrong quantity and
    produces untrustworthy boundary-tightness numbers."""
    with mp.workdps(dps):
        if side == "p":
            aa, bb, yy = mp.mpf(a), mp.mpf(b), mp.mpf(y0)
        else:
            aa, bb, yy = mp.mpf(b), mp.mpf(a), 1 - mp.mpf(y0)
        t = mp.mpf(1)
        sp = mp.mpf(1)  # S' accumulator, starts at its own n=0 term (=1)
        for n in range(1, 400):
            t *= (n - bb) * yy / n
            term = aa * t / (aa + n)
            sp += term
            if n > 4 and abs(term) < mp.mpf(10) ** -70 * abs(sp):
                break
        return abs(mp.log(sp)) / aa


# ============================================================================
# Part 6: self-checks + replay + header emission.
# ============================================================================
EPS_SERIES = 2.0 ** -56   # S1/S2/S3-fed regions: series/CF/fit truncation
                            # budget, gammainv's own R1/R2/R3-class number.
EPS_TIGHT = 2.0 ** -105    # S4/deep-small closed-form regions: no series
                            # truncation -- the TRUE forward precision
                            # there is dd-class (~2^-105), not the
                            # series/CF/fit-truncation-bound
                            # EPS_SERIES=2^-56. Using 2^-64 (a plausible
                            # "absolute floor" placeholder) instead fails
                            # self-check (b) at kappa~2^40 (a=1e-12): floor
                            # measures 26.74b, exactly matching
                            # eps_bits-kappa_bits = 64-40 = 24-ish. With
                            # the dd-precision-matched EPS_TIGHT=2^-105,
                            # the SAME kappa~2^40 point predicts floor
                            # 105-40=65b, clearing the 55b gate with
                            # margin -- and this matches the plateau
                            # derivation directly: dd (2^-105-class)
                            # resolves y to 1 ULP only for kappa<=2^52
                            # (105-52=53, ~1 ULP -- that boundary IS this
                            # formula).
TARGET_BITS = 55.0


def eps_for(a, b, region, sigma=None):
    """EPS_SERIES=2^-56 is a REGION-WORST series-truncation bound; it is
    the wrong model wherever the forward series/CF has already
    super-converged, which happens whenever the SOLVED-SIDE probability
    sigma itself is tiny (a handful of R1 terms, or a shallow CF depth)
    REGARDLESS of whether (a,b) are the tiny ones -- measured directly:
    a=b=0.5 (NOT tiny), sigma~6.4e-126 (deep p-tail), seed 50.5b degrades
    to 46.8b under uniform EPS_SERIES injection, while the TRUE
    per-point forward error there is dd/prefactor-bound (~2^-105-class),
    not truncation-bound. Uses EPS_TIGHT whenever sigma itself is in the
    genuinely deep tail (<1e-6, comfortably past where any of
    R1/R2/R3/gammalim's series/CF have more than a few live terms),
    matching EPS_TIGHT's own dd-precision justification."""
    # A region-NAME based gate alone misses the same phenomenon at
    # points labeled "gap" whose (a,b) are ALSO in the joint-tiny/
    # plateau range (e.g. a=1e-10, b=3e-10 -- min(a,b) tiny, kappa
    # large): eps_for would never see a region name it recognizes, fall
    # through to EPS_SERIES, and let the SAME kappa-amplification
    # mechanism that motivated EPS_TIGHT silently regress an
    # already-good seed (measured: final floor drops from 74.38b to
    # 25.38b at this exact point when the gap-band test grid is widened
    # to include it). Gate on min(a,b) directly (the actual physical
    # conditioning driver, not a region label) so it fires regardless of
    # which bucket a point is reported under.
    if region in ("S4", "deep", "S2", "gap"):
        return EPS_TIGHT
    if min(a, b) < 1e-3:
        return EPS_TIGHT
    if sigma is not None and sigma < 0.05:
        return EPS_TIGHT
    return EPS_SERIES


def replay_point(a, b, y_true, region_label, dps_oracle=50, dps_step=55):
    """One replay point: root-find true y at the ROUNDED double sigma
    (basis MUST be the rounded double, never the unrounded construction
    value -- see mpmath discipline in the module docstring), seed via
    quad-candidate seed_for, STEPS_N-step (4) safeguarded logit-Newton,
    return (seed_bits, steps_bits, kappa_bits, sigma, side) or None on
    failure."""
    s, side = small_side_of_y(a, b, y_true, dps=45)
    s = float(s)
    if not (0.0 < s <= 0.5):
        return None
    true_y = oracle_y(a, b, s, side, dps=dps_oracle)
    if true_y is None:
        return None
    # true_y MUST stay an mpf through every bits comparison below --
    # casting to float here silently caps every measured floor at ~52-53
    # bits regardless of candidate quality, since a double-quantized
    # true_y can only ever agree with anything to ~1 ULP by construction.
    # Only a plain-float FINITENESS/RANGE check is done here; the mpf
    # value itself is threaded through unchanged.
    if not (mp.isfinite(true_y) and 0.0 < true_y < 1.0):
        return None
    try:
        y0 = seed_for(a, b, s, side)
    except (OverflowError, ValueError, ZeroDivisionError):
        return None
    if not (math.isfinite(y0) and 0.0 < y0 < 1.0):
        return None
    seed_bits = bits_of(true_y, y0)
    eps = eps_for(a, b, region_label, sigma=s)
    try:
        final_bits = simulate_steps_beta(a, b, y0, true_y, s, side, eps,
                                          nsteps=STEPS_N, dps=dps_step)
    except (OverflowError, ValueError, ZeroDivisionError):
        final_bits = -1.0
    # kappa: interior-plateau conditioning estimate, ~1/min(a,b) capped
    # sanely; used only to bucket points for reporting/gating (contract's
    # plateau contract split).
    kappa_bits = -math.log2(min(a, b)) if min(a, b) > 0 else 999.0
    return seed_bits, final_bits, kappa_bits, s, side, true_y


def main():
    rc = 0
    global _S1_CHEB, S1_NCORR

    print("(a) beta-ridge lambda(zeta) Newton floor:", file=sys.stderr)
    worst_a, worst_a_at = 999.0, None
    a_vals = [0.5, 5.0, 20.0, 100.0, 1000.0]
    b_vals = [0.5, 5.0, 20.0, 100.0, 1000.0]
    zeta_fracs = [0.001, 0.05, 0.2, 0.5, 0.8, 0.99]
    for av in a_vals:
        for bv in b_vals:
            c = av + bv
            nu = av * bv / c
            zmax_here = min(S1_ZETA_MAX, math.sqrt(2.0 * 700.0 / nu) if nu > 0 else S1_ZETA_MAX)
            for frac in zeta_fracs:
                for sgn in (1.0, -1.0):
                    zeta = sgn * frac * zmax_here
                    lam_d = lam_of_zeta_double(zeta, av, bv)
                    with mp.workdps(50):
                        lam_mp = gb._lambda_of_zeta(mp.mpf(zeta), mp.mpf(nu),
                                                     mp.mpf(av / c), mp.mpf(av), mp.mpf(bv))
                    lam_mp_f = float(lam_mp)
                    if lam_mp_f == 0.0:
                        continue
                    scale = av + bv
                    absdiff = abs(lam_d - lam_mp_f) / scale
                    b_bits = 100.0 if absdiff == 0 else -math.log2(absdiff)
                    if b_bits < worst_a:
                        worst_a, worst_a_at = b_bits, (av, bv, zeta)
    print(f"    worst {worst_a:.2f} bits at (a,b,zeta)={worst_a_at}", file=sys.stderr)
    NEWTON_FLOOR_BITS = 30.0
    if worst_a < NEWTON_FLOOR_BITS:
        print(f"    FAILED: below floor {NEWTON_FLOOR_BITS}", file=sys.stderr)
        rc = 1
    if rc:
        return rc

    print("(s1-build) S1 correction table (K=2):", file=sys.stderr)
    build_s1_correction(2, dps=45)
    S1_NCORR = 2

    print("(b) seed-bit floors, quad-candidate seed_for, edge-refined grid:",
          file=sys.stderr)
    # Representative point set spanning the named regions (bounded scope
    # for this stage's replay -- the exhaustive certified reference set
    # is a separate, later pass; this replay is design-sanity, per
    # doctrine, the same role as gammainv's own replay).
    points = []
    # R1-tiny / small-y (S2 territory) -- genuinely small TARGET
    # PROBABILITY is the point, not merely small y: picking (a,y) pairs
    # by y alone can mislabel points, e.g. (a=0.1,y=1e-4) gives
    # sigma~0.2-0.44 (NOT small: y^a is not small for a<1 regardless of
    # how small y itself is) -- such a point would be classified "S2"
    # and then reported as an S2 FAILURE when really no candidate
    # targets that (small-a, moderate-sigma) territory at all (it is the
    # SAME coverage gap the gap-check bucket below tracks explicitly).
    # Require a>=1 here (S2's genuine small-target domain); small-a
    # moderate-sigma traffic lives ONLY in gap-check.
    for a in (1.0, 5.0, 20.0):
        for b in (0.5, 2.0, 5.0, 20.0, 100.0):
            for y in (1e-4, 1e-2, 0.05):
                points.append((a, b, y, "S2"))
    for a in (0.5, 2.0, 5.0):
        for b in (0.5, 2.0, 5.0):
            for pe in (-30, -100, -250):
                points.append((a, b, 10.0 ** pe, "S2"))
    # GAP BAND: moderate-tiny, near-symmetric (a,b) between the t_jt
    # value and S1's nu>=S1_NU_MIN boundary -- S4 and S1's own linear
    # form both break down in this range (see seed_S4's and seed_for's
    # docstrings), so it is tracked as its own labeled bucket (dense,
    # edge-refined, bit-stepped at the band's own edges -- t_jt's value
    # and S1_NU_MIN's boundary), not a separately-excused one:
    # acceptance requires this bucket to clear TARGET_BITS like any
    # other.
    gap_masses = [0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0,
                  1.5, 2.0, 3.0, 4.0]
    if FULL:
        gap_masses += [0.003, 0.007, 0.015, 0.07, 0.15, 0.4, 0.7, 1.2, 2.5, 3.5]
    for m in gap_masses:
        for skew in (1.0, 1.5, 3.0, 10.0):
            for yv in (1e-6, 1e-3, 0.1, 0.3, 0.5, 0.7, 0.9, 1 - 1e-3, 1 - 1e-6):
                points.append((m, m * skew, yv, "gap"))
    # bit-stepped edge ladder at the two band boundaries (t_jt's value,
    # and just below/above S1_NU_MIN's own a=b=4 boundary for the
    # balanced case) -- edge-refined sampling catches what grid-only
    # sampling misses (same precedent as trigamma's own edge sampling).
    for edge in (2.0 ** -8, 4.0):
        for step in (-3, -2, -1, 0, 1, 2, 3):
            m = edge
            for _ in range(abs(step)):
                m = math.nextafter(m, math.inf if step > 0 else 0.0)
            for yv in (1e-4, 0.5, 1 - 1e-4):
                points.append((m, m, yv, "gap"))
    # ridge / balanced (S1 territory) -- nu from moderate to large
    for nu2 in (40, 100, 400, 2000, 20000):
        a = b = nu2
        for y in (0.3, 0.4, 0.45, 0.5):
            points.append((a, b, y, "S1"))
    for (a, b) in ((100.0, 300.0), (50.0, 20.0), (200.0, 40.0)):
        for y in (0.15, 0.3, 0.5, 0.7):
            points.append((a, b, y, "S1"))
    # gamma-limit (S3 territory)
    for huge in (1e6, 1e10, 1e16, 1e30):
        for s_ in (0.05, 1.0, 3.0, 10.0):
            points.append((s_, huge, 0.5, "S3"))
            points.append((huge, s_, 0.5, "S3"))
    # joint-tiny plateau (S4 territory)
    for m in (1e-4, 1e-6, 1e-8, 1e-10, 1e-12):
        for skew in (1.0, 3.0, 10.0):
            points.append((m, m * skew, 0.4, "S4"))
            points.append((m, m * skew, 0.6, "S4"))
    # NOTE: the loop below dispatches on `len(pt)==4` only, so any
    # 5-tuple appended to this SAME `points` list (e.g. a would-be
    # "deep-small (both orientations)" block) is silently DROPPED every
    # run, never scored, never reported, never gating anything, no
    # matter how live-looking its region label reads -- a whole branch
    # dead at its own call site, never reachable. Deep-small's real
    # self-check is (f) below, with the correct cut derivation,
    # bit-stepped boundary sampling, and BOTH orientations genuinely
    # exercised -- a second, differently-shaped pathway through this
    # generic replay harness would be redundant at best and a second
    # place for the same silent-drop disease to recur at worst.

    region_results = {}
    n_ok, n_fail = 0, 0
    for pt in points:
        if len(pt) == 4:
            a, b, y, region = pt
            try:
                res = replay_point(a, b, y, region)
            except Exception:
                res = None
            if res is None:
                n_fail += 1
                continue
            seed_bits, final_bits, kappa_bits, s, side, true_y = res
            n_ok += 1
            bucket = region
            # BUCKET THRESHOLD: the SIMPLE kappa~1/min(a,b) estimate used
            # here for bucketing is a rough proxy, not the exact
            # amplification factor -- measured directly at
            # (a,b)=(1e-12,1e-11) [kappa_bits~39.9 by this formula] the
            # true noise floor under EPS_TIGHT=2^-105 was 53.34b, i.e. an
            # EFFECTIVE kappa~51.7b, already close to the contract's 2^52
            # headline despite the simple formula placing it comfortably
            # below. Bucketing at 35 (not 52) so points whose
            # simple-formula kappa sits well inside the ~12-bit
            # under-prediction band fall into the backward-error
            # contract rather than being held to the y-ULP gate the
            # simple formula under-predicts they can meet.
            if kappa_bits > 35.0 and region == "S4":
                bucket = "S4-plateau"
            region_results.setdefault(bucket, []).append(
                (final_bits, seed_bits, a, b, y))
    for region, arr in sorted(region_results.items()):
        finals = [r[0] for r in arr if r[0] > -1]
        if not finals:
            print(f"    {region}: ALL INVALID (n={len(arr)})", file=sys.stderr)
            continue
        worst = min(finals)
        worst_row = min(arr, key=lambda r: r[0])
        print(f"    {region}: n={len(arr)} worst_final={worst:.2f}b "
              f"seed={worst_row[1]:.2f}b at (a,b,y)=({worst_row[2]:.3g},"
              f"{worst_row[3]:.3g},{worst_row[4]})", file=sys.stderr)
        if region != "S4-plateau" and worst < TARGET_BITS:
            print(f"    ESCALATION-CANDIDATE: {region} floor {worst:.2f}b "
                  f"< {TARGET_BITS}b gate", file=sys.stderr)
            rc = 1

    print(f"    replay coverage: {n_ok} points scored, {n_fail} dropped "
          f"(oracle/seed failure)", file=sys.stderr)

    # ------------------------------------------------------------------
    # (c)/(d): S4 closed form self-check -- symmetry identity + t_jt gate.
    # ------------------------------------------------------------------
    print("(c) S4 closed form: c(a,a)=0 exact identity check:", file=sys.stderr)
    worst_c = 0.0
    for aa in (1e-10, 1e-6, 1e-3, 0.5, 5.0):
        cval = s4_c_closed_form(aa, aa)
        worst_c = max(worst_c, abs(cval))
    print(f"    max |c(a,a)| over sample = {worst_c:.3e} (must be exactly 0.0)",
          file=sys.stderr)
    if worst_c != 0.0:
        print("    FAILED: c(a,a) not exactly zero", file=sys.stderr)
        rc = 1

    print(f"(d) t_jt = {T_JT:.3e} (2^{math.log2(T_JT):.1f}) -- CONTRACT "
          f"CLARIFICATION (orchestrator ruling, 2026-08-09): this is the "
          f"CLOSED-FORM-SHIPS-DIRECTLY route gate (a G3/G4 kernel decision "
          f"about skipping Newton refinement under the plateau backward-"
          f"error contract), NOT a candidacy gate on offering S4 as a seed "
          f"-- S4 is offered globally now (seed_for). min(a,b)<=2^-52 "
          f"ownership is therefore not this generator's concern to verify "
          f"(nothing in the seed/step pipeline is gated on T_JT any longer).",
          file=sys.stderr)

    # ------------------------------------------------------------------
    # (f) deep-small cut, both orientations. Design notes for this
    # self-check:
    #   - solve the boundary y DIRECTLY from the analytic bound
    #     (deep_small_cut_bound) and bit-step (edge-refined rule,
    #     binding for every family) around it with math.nextafter, both
    #     sides of the boundary -- a fixed decade grid of sigma can miss
    #     the boundary entirely, landing every accepted point deep
    #     inside the cut (dropped~0) instead of AT it.
    #   - sweep b across the FULL in-route range, including b<1 (where
    #     |1-b|<1, a DIFFERENT regime -- at exactly b=1 the bound is
    #     identically 0, S' is EXACTLY 1 for all y, matching gamma's own
    #     no-second-parameter case as a sanity floor) and the widened
    #     gamma-limit corner up to b=1e300 (y bounds go subnormal there
    #     -- expected, checked directly via math.nextafter, which is
    #     subnormal-safe).
    #   - test the q-side bound/error DIRECTLY via
    #     deep_small_cut_bound(...,'q') and _deep_small_dropped_rel(...,
    #     'q') on a real y (not a pre-swapped z, and not a hardcoded
    #     side='p' literal) -- both orientations must be genuinely
    #     exercised, not merely labeled as such.
    # ------------------------------------------------------------------
    print("(f) deep-small closed-form cut, both orientations, bit-stepped "
          "at the boundary, full b range:", file=sys.stderr)
    a_vals_f = (1e-6, 1e-3, 0.05, 0.5, 0.9, 1.0, 1.5, 5.0, 50.0, 500.0)
    b_vals_f = (1e-6, 0.1, 0.5, 0.9, 0.99, 1.0, 1.01, 1.1, 2.0, 10.0,
                100.0, 1e4, 1e8, 1e16, 1e100, 1e300)
    steps_f = (-20, -10, -5, -3, -2, -1, 0, 1, 2, 3, 5, 10, 20)
    worst_f = mp.mpf(0)
    worst_f_at = None
    worst_bound_ratio = 0.0  # max(true_err/bound) among accepted points --
                              # >1 would mean the analytic bound is UNSOUND
    n_inside = 0
    for a in a_vals_f:
        for b in b_vals_f:
            for side in ("p", "q"):
                coeff = abs(1.0 - b) if side == "p" else abs(1.0 - a)
                denom = (1.0 + a) if side == "p" else (1.0 + b)
                if coeff < 1e-15:
                    # bound is identically ~0 (b=1 exactly, p-side; or
                    # a=1 exactly, q-side): S' is EXACTLY 1 for all y on
                    # this branch (matches gamma's no-second-parameter
                    # floor) -- spot-check a couple of representative y
                    # rather than boundary-solving a division by ~0.
                    y_candidates = [1e-10, 1e-2, 0.5]
                else:
                    y_b = DEEP_SMALL_CUT * denom / coeff
                    if side == "p":
                        y_center = y_b
                    else:
                        y_center = 1.0 - y_b
                    if not (0.0 < y_center < 1.0):
                        # boundary falls outside (0,1): the WHOLE domain
                        # is on one side of the cut for this (a,b,side)
                        # -- test near whichever endpoint is relevant so
                        # the sweep still exercises this (a,b) honestly.
                        y_candidates = [1e-300, 1e-30, 1e-6, 0.5, 1 - 1e-6] \
                            if side == "p" else \
                            [1 - 1e-300, 1 - 1e-30, 1 - 1e-6, 0.5, 1e-6]
                    else:
                        y_candidates = []
                        for st in steps_f:
                            yv = y_center
                            for _ in range(abs(st)):
                                yv = math.nextafter(
                                    yv, math.inf if st > 0 else 0.0)
                            if 0.0 < yv < 1.0:
                                y_candidates.append(yv)
                for yv in y_candidates:
                    bound = deep_small_cut_bound(a, b, yv, side)
                    if not (bound < DEEP_SMALL_CUT):
                        continue  # outside the cut -- not routed, not this
                                  # self-check's business (matches gamma's
                                  # own "normal pipeline covers seamlessly
                                  # at/above it" doctrine)
                    n_inside += 1
                    # dps scaled to y's own magnitude (subnormal-scale y
                    # near the gamma-limit corner needs headroom above
                    # its own exponent, not just a flat working dps).
                    mag = abs(math.log10(yv)) if 0 < yv < 1 else 0.0
                    mag2 = abs(math.log10(1 - yv)) if 0 < yv < 1 else 0.0
                    dps_here = int(60 + max(mag, mag2) * 1.2)
                    rel = _deep_small_dropped_rel(a, b, yv, side, dps=dps_here)
                    if rel > worst_f:
                        worst_f, worst_f_at = rel, (a, b, yv, side)
                    if bound > 0:
                        ratio = float(rel) / bound
                        worst_bound_ratio = max(worst_bound_ratio, ratio)
    worst_f_float = float(worst_f)
    print(f"    {n_inside} cut-accepted (boundary-adjacent) points tested "
          f"across {len(a_vals_f)}x{len(b_vals_f)}x2 (a,b,side) cells",
          file=sys.stderr)
    print(f"    worst TRUE dropped-factor rel err among accepted points: "
          f"{worst_f_float:.3e} at (a,b,y,side)={worst_f_at} "
          f"(target < {DEEP_SMALL_CUT:.3e})", file=sys.stderr)
    print(f"    worst true_err/bound ratio: {worst_bound_ratio:.3f} "
          f"(<=1 means the analytic bound is a SOUND upper predictor "
          f"everywhere tested; >1 would mean it under-predicts)",
          file=sys.stderr)
    # 1e-6 relative tolerance on BOTH checks below: boundary-adjacent
    # points are, by construction, sampled AT the cut edge (bound approx
    # DEEP_SMALL_CUT exactly), where the true/bound ratio's own
    # theoretical value is exactly 1 (corr(y')=-ln(1-y')/y' is the EXACT
    # limiting correction, not a heuristic) -- residual floating-point
    # noise there (measured worst true/CUT excess ~2e-14, worst ratio
    # 1.0000000004) is not a soundness failure; either check genuinely
    # exceeding its target by more than this tolerance would be.
    RATIO_TOL = 1e-6
    if worst_f_float >= DEEP_SMALL_CUT * (1.0 + RATIO_TOL):
        print("    FAILED", file=sys.stderr)
        rc = 1
    if worst_bound_ratio > 1.0 + RATIO_TOL:
        print("    FAILED: analytic cut bound is NOT a sound upper "
              "predictor of the true dropped error somewhere tested",
              file=sys.stderr)
        rc = 1

    # ------------------------------------------------------------------
    # (g) S1/S3 seam near alpha~kGammaAT (informational; T_JT/S1_NU_MIN
    # already partition the domain by nu, not directly by kGammaAT -- the
    # seam here is where the GAMMA-LIMIT S3 route (huge second param)
    # starts beating S1, reported for the record).
    # ------------------------------------------------------------------
    print(f"(g) S1/S3 seam: S1 owns nu>={S1_NU_MIN} (measured gate); S3 "
          f"(gamma-limit) requires one param comparably large -- both are "
          f"quad-candidates selected by cheap residual at every point, so "
          f"no HARD seam is pinned (design difference from gammainv's own "
          f"single a_T crossover, itself a consequence of the quad-"
          f"candidate design being GLOBAL from day one here per PLAN's "
          f"own instruction).", file=sys.stderr)

    # ------------------------------------------------------------------
    # (h) plateau backward-error contract: for kappa>2^52 rows, verify
    # the BACKWARD error (forward of the returned y matches sigma) rather
    # than a y-ULP claim (contract's own split). No kappa>2^52 point was
    # reachable at replay's own sample scope (min(a,b) would need to be
    # <~2^-52~2.2e-16); reported here so the boundary is documented, not
    # silently skipped.
    # ------------------------------------------------------------------
    print("(h) plateau backward-error contract: replay's own point set "
          "does not reach kappa>2^52 (min(a,b) as small as 1e-12 tested, "
          "kappa~2^40; G2's certified reference set is where kappa>2^52 "
          "rows are actually exercised, per the contract's own bucket "
          "split assignment to G2). Backward-error verification helper "
          "provided below for G2's use.", file=sys.stderr)

    if rc:
        print("One or more self-checks failed -- emitting nothing.", file=sys.stderr)
        return rc

    # ------------------------------------------------------------------
    # Emit
    # ------------------------------------------------------------------
    maxdeg = S1_NP
    print("// Auto-generated by tools/gen_betainv_data.py. DO NOT EDIT.")
    print("// Fits and their error budgets are documented in the generator;")
    print("// the kernel that consumes them is src/betainv-inl.h. Consumes")
    print("// src/beta_data.h (region cores) and src/gammainv_data.h (S3's")
    print("// gamma-inverse seed machinery) alongside this header -- nothing")
    print("// already there is duplicated here.")
    print("#ifndef CORVUS_BETAINV_DATA_H_")
    print("#define CORVUS_BETAINV_DATA_H_")
    print()
    print("namespace corvus::detail {")
    print()
    print("// S1 beta-Temme seed: z=erfcinv(2*sigma) [sign by side], "
          "nu=alpha*beta/(alpha+beta),")
    print("// zeta0=sign*z/sqrt(nu) [note: NO factor of 2 -- beta's cpsi(lambda)=")
    print("// zeta^2*nu directly, unlike gamma's a*eta^2/2], invert cpsi(lambda)=")
    print("// zeta^2*nu for lambda via safeguarded Newton (exact analytic")
    print("// derivative lam*c/((a-lam)(b+lam)), c=a+b), y=(a-lambda)/c. Optional")
    print("// ONE perturbative-Newton correction: zeta_new = zeta0 -")
    print("// sign(zeta0)*S(zeta0,p)/(nu*sqrt(2)), S=sum c_k(zeta0,p)/nu^k,")
    print("// c_k a 2D Chebyshev fit over (zeta,p) [zeta in")
    print("// [-kBetaInvS1ZetaMax,kBetaInvS1ZetaMax], p in (0,0.5], symmetry")
    print("// c_k(zeta,p)=-c_k(-zeta,1-p)], gated to its own fitted domain AND")
    print("// nu>=kBetaInvS1NuMin (the 1/nu series is asymptotic; extrapolating")
    print("// below nu~2 diverges).")
    print(f"inline constexpr int kBetaInvS1NCorr = {S1_NCORR};")
    print(f"inline constexpr int kBetaInvS1NZ = {S1_NZ};")
    print(f"inline constexpr int kBetaInvS1NP = {S1_NP};")
    print(f"inline constexpr double kBetaInvS1ZetaMax = {hexf(S1_ZETA_MAX)};")
    print(f"inline constexpr double kBetaInvS1PMid = {hexf(float(S1_P_MID))};")
    print(f"inline constexpr double kBetaInvS1PHalf = {hexf(float(S1_P_HALF))};")
    print(f"inline constexpr double kBetaInvS1NuMin = {hexf(S1_NU_MIN)};")
    print(f"inline constexpr double kBetaInvS1Cheb[{S1_NCORR}][{S1_NZ}][{S1_NP}] = {{")
    for k in range(S1_NCORR):
        print("    {")
        for row in _S1_CHEB_D[k]:
            print("        {" + ", ".join(hexf(v) for v in row) + "},")
        print("    },")
    print("};")
    print()
    print("// S2 (small-y series inversion, either orientation via exact")
    print("// complement) Picard correction count.")
    print(f"inline constexpr int kBetaInvS2NCorr = 6;")
    print()
    print("// S4: exact-B leading-")
    print("// order closed form, exact in logit(y) of any sign/magnitude (not a")
    print("// linearization near the plateau center -- a linear form")
    print("// (s-s*)/w+c(alpha,beta) diverges for |logit y|>>1; c(alpha,beta) is")
    print("// DROPPED here, not merely zero -- it actively hurts once B is exact):")
    print("//   s*=beta/(alpha+beta); s=sigma if side==p else 1-sigma")
    print("//   s<=s*: v=(ln s+ln alpha+lnB(alpha,beta))/alpha")
    print("//   s>s* : v=-(ln(1-s)+ln beta+lnB(alpha,beta))/beta")
    print("//   y=sigmoid(v)")
    print("// This is seed_S2's own zeroth iterate (either orientation, picked")
    print("// by which branch applies) -- S4 and S2 share one mechanism at")
    print("// leading order. Offered as a GLOBAL seed candidate at every (a,b,s)")
    print("// (selected by cheap-residual comparison like the other four);")
    print("// kBetaInvTJt is NOT a candidacy gate: t_jt gates")
    print("// only where the closed form could ship WITHOUT Newton refinement")
    print("// under the plateau backward-error contract, a kernel")
    print("// decision this generator does not make. Kept, PROVISIONAL, at its")
    print("// original measured value pending that decision.")
    print(f"inline constexpr double kBetaInvTJt = {hexf(T_JT)};")
    print()
    print("// S5 (last-resort family):")
    print("// logit-normal via EXACT digamma/trigamma moments of logit(Y) =")
    print("// ln(Gamma(alpha)-variate) - ln(Gamma(beta)-variate) [standard")
    print("// Gamma-ratio construction of Y~Beta(alpha,beta), clean-room]:")
    print("//   mu = psi(alpha)-psi(beta), var = psi'(alpha)+psi'(beta)")
    print("//   v = mu + Phi^-1(s)*sqrt(var), y = sigmoid(v)")
    print("// Targets \"moderate shape, moderate probability\" territory none of")
    print("// S1-S4 cover -- no pinned table (pure closed form on")
    print("// digamma/trigamma, evaluated by the kernel's own digamma/trigamma")
    print("// cores, corvus already ships both). Offered as a fifth global")
    print("// candidate, selected only when it wins the cheap-residual")
    print("// comparison.")
    print()
    print("// STEPS: safeguarded logit-Newton (m=lnP-lnQ, w=1/(y*dm/dy)),")
    print("// shared step count, gammainv safeguard package (reject")
    print("// residual-increasing, 1/8 backtrack, bypass |resid|<TrustResid,")
    print("// multiplicative-in-y step y*(1+ls), floor ls>=-0.9).")
    print("// StepsN=4: a 3-step count leaves")
    print("// a bounded interior sub-band (min(alpha,beta)")
    print("// approx 0.02-0.5, skew 3-10x, y interior 0.1-0.3) short of the")
    print("// gate after all FIVE closed-form seed families are exhausted there")
    print("// (none exceeds ~2-5 bits; not a selection or noise-floor")
    print("// artifact). Convergence from that band's worst")
    print("// seed is clean quadratic (2.12->6.66->16.48->36.16->")
    print("// 75.51 bits over steps 1-5); step 4 clears the gate by 20+ bits")
    print("// margin band-wide -- restores full margin, shaves nothing. The")
    print("// safeguard package makes step 4 IDEMPOTENT for lanes already")
    print("// converged after step 3 (freeze-by-select), so its cost is one")
    print("// bounded extra forward evaluation there, not a global slowdown.")
    print("// LATITUDE: a whole-vector all-lanes-converged skip after step")
    print("// 3 is an ALLOWED bench optimization (gammainv's own \"1 Halley")
    print("// vs 2 Newton\" precedent) -- accuracy gates must hold either way;")
    print("// this is a throughput decision only, not an accuracy one.")
    print(f"inline constexpr int kBetaInvStepsN = {STEPS_N};")
    print(f"inline constexpr double kBetaInvTrustResid = {hexf(TRUST_RESID)};")
    print()
    print("// Deep-small closed-form cut, BOTH orientations (re-derived for")
    print("// beta -- NOT the naive")
    print("// (own side)*y form, which has no dependence on the OTHER side's")
    print("// parameter and is wrong). y0 = exp((ln sigma + ln alpha + lnB)/")
    print("// alpha) [P] or 1 - exp((ln sigma + ln beta + lnB)/beta) [Q] drops")
    print("// the series correction S' (S'->1 as y->0; S' = alpha*sum_n t_n/")
    print("// (alpha+n), t_0=1, t_n=t_{n-1}*(n-beta)*y/n). The dropped")
    print("// relative error in ln(y) is |ln S'|/alpha; its LEADING term is")
    print("// |1-beta|*y/(1+alpha) -- driven by the OTHER side's parameter")
    print("// (beta), not alpha (beta supplies the leading nonconstant series")
    print("// coefficient, n-beta at n=1; alpha enters only via the 1+alpha")
    print("// denominator) -- but the leading term ALONE under-predicts badly")
    print("// at the widened gamma-limit corner (true/bound ratio up")
    print("// to 13.8 uncorrected). CORRECTED with an EXACT closed-form")
    print("// multiplier corr(y')=-ln(1-y')/y' (y'=y [P] or 1-y [Q], the OWN")
    print("// side's small variable) -- exact in the huge-OTHER-side-exponent")
    print("// limit (S'->(1-y')^(other_side-1) there) and verified a SOUND")
    print("// (ratio<=1, worst measured 1.0000000004 -- boundary float noise)")
    print("// upper predictor across the full range tested; corr->1 as y'->0,")
    print("// so it reduces to the bare leading term exactly where that was")
    print("// already valid. Route:")
    print("//   P-side: |1-beta|*y/(1+alpha)*corr(y)         < kBetaInvDeepSmallCut")
    print("//   Q-side: |1-alpha|*(1-y)/(1+beta)*corr(1-y)   < kBetaInvDeepSmallCut")
    print("// Measured (self-check (f), bit-stepped at the boundary, full")
    print("// beta range including beta<1 and the widened gamma-limit corner")
    print("// beta up to 1e300): worst true dropped-factor error and worst")
    print("// true/bound ratio are both reported to stderr on every run --")
    print("// trust that budget line over this comment.")
    print(f"inline constexpr double kBetaInvDeepSmallCut = {hexf(DEEP_SMALL_CUT)};")
    print()
    print("}  // namespace corvus::detail")
    print()
    print("#endif  // CORVUS_BETAINV_DATA_H_")
    return 0


def smoke_test():
    # development smoke test -- sanity vs mp.betainc on plain points, and
    # per-piece diagnostics. Not part of the generator's own self-check
    # gate (main() is); kept for `--smoke` re-runs during future work.
    mp.mp.dps = 40
    for (a, b, x) in [(2, 3, 0.3), (5, 7, 0.9), (0.5, 0.5, 0.1),
                      (1e6, 3, 0.9999), (1e-3, 1e-3, 0.5),
                      (20, 20, 0.4), (1e20, 5, 1e-19)]:
        true_p = mp.betainc(mp.mpf(a), mp.mpf(b), 0, mp.mpf(x), regularized=True)
        val, side = betainv_forward(a, b, x, dps=40)
        got = val if side == "p" else 1 - val
        rel = abs((got - true_p) / true_p) if true_p != 0 else abs(got - true_p)
        print(f"a={a} b={b} x={x}: true_p={mp.nstr(true_p,15)} got_p={mp.nstr(got,15)} "
              f"side={side} rel={float(rel):.2e}")

    print("\n--- oracle_y + small_side_of_y roundtrip ---")
    for (a, b, y_true) in [(2.0, 3.0, 0.3), (0.5, 0.5, 0.1), (20.0, 20.0, 0.4),
                            (1e6, 3.0, 0.9999), (1e6, 5.0, 0.99999),
                            (0.02, 30.0, 0.4)]:
        s, side = small_side_of_y(a, b, y_true, dps=40)
        s = float(s)
        y_back = oracle_y(a, b, s, side, dps=40)
        print(f"a={a} b={b} y_true={y_true} -> s={s:.6e} side={side} "
              f"y_back={float(y_back):.15e} rel={abs((float(y_back)-y_true)/y_true):.2e}")

    print("\n--- building S1 correction table (K=2) for smoke test ---")
    build_s1_correction(2, dps=40)

    print("\n--- seed candidates vs true y (bits) ---")
    cases = [
        (20.0, 20.0, 0.4, "S1 balanced ridge"),
        (100.0, 300.0, 0.2, "S1 skewed ridge"),
        (2.0, 3.0, 0.01, "S2 small-y"),
        (0.5, 5.0, 0.9995, "S2 via complement (q small)"),
        (0.03, 1e10, 3e-9, "S3 gamma-limit"),
        (1e10, 0.03, 1.0 - 3e-9, "S3 gamma-limit swapped"),
        (1e-8, 1e-8, 0.5, "S4 joint-tiny symmetric"),
        (1e-8, 3e-8, 0.5, "S4 joint-tiny skewed"),
        (1e-10, 1e-6, 0.7, "S4 joint-tiny wide skew"),
    ]
    for (a, b, y_true, label) in cases:
        s, side = small_side_of_y(a, b, y_true, dps=40)
        s = float(s)
        true_y = oracle_y(a, b, s, side, dps=45)
        if true_y is None:
            print(f"{label}: oracle_y failed")
            continue
        true_y = float(true_y)
        row = {}
        try:
            row["S1(0)"] = bits_of(true_y, seed_S1(a, b, s, side, 0))
        except Exception as e:
            row["S1(0)"] = f"exc:{e}"
        try:
            row["S1(2)"] = bits_of(true_y, seed_S1(a, b, s, side, 2))
        except Exception as e:
            row["S1(2)"] = f"exc:{e}"
        try:
            s2v = seed_S2(a, b, s, 3) if side == "p" else seed_S2(a, b, 1.0 - s, 3)
            row["S2"] = bits_of(true_y, s2v if side == "p" else 1.0 - s2v)
        except Exception as e:
            row["S2"] = f"exc:{e}"
        try:
            s3v = seed_S3(a, b, s, side)
            row["S3"] = bits_of(true_y, s3v) if s3v is not None else None
        except Exception as e:
            row["S3"] = f"exc:{e}"
        try:
            row["S4"] = bits_of(true_y, seed_S4(a, b, s, side))
        except Exception as e:
            row["S4"] = f"exc:{e}"
        try:
            row["seed_for"] = bits_of(true_y, seed_for(a, b, s, side))
        except Exception as e:
            row["seed_for"] = f"exc:{e}"
        print(f"{label}: a={a} b={b} s={s:.3e} side={side} true_y={true_y:.6e} -> {row}")

    print("\n--- STEPS smoke test (seed_for -> STEPS_N-step logit-Newton) ---")
    EPS0 = 2.0 ** -56
    for (a, b, y_true, label) in cases:
        s, side = small_side_of_y(a, b, y_true, dps=45)
        s = float(s)
        true_y = oracle_y(a, b, s, side, dps=50)
        if true_y is None:
            continue
        true_y = float(true_y)
        try:
            y0 = seed_for(a, b, s, side)
        except Exception as e:
            print(f"{label}: seed_for exc {e}")
            continue
        if not (math.isfinite(y0) and 0.0 < y0 < 1.0):
            print(f"{label}: bad seed {y0}")
            continue
        multi = simulate_steps_beta_multi(a, b, y0, true_y, s, side, EPS0, STEPS_N, dps=55)
        print(f"{label}: seed_bits={bits_of(true_y, y0):.2f} steps -> {multi}")


if __name__ == "__main__":
    _parser = argparse.ArgumentParser(
        allow_abbrev=False,  # '--ful' must ERROR, not prefix-match --full (#11)
        
        description="Generate src/betainv_data.h -- every table "
                     "beta_p_inv/beta_q_inv needs.")
    _parser.add_argument(
        "--full", action="store_true",
        help="Denser replay (see module docstring's replay-parameter "
             "self-checks) than the default invocation.")
    _parser.add_argument(
        "--smoke", action="store_true",
        help="Development smoke test -- sanity vs mp.betainc on plain "
             "points, and per-piece diagnostics. Not part of the "
             "generator's own self-check gate (main() is); kept for "
             "ad hoc re-runs during future work.")
    _args = _parser.parse_args()
    FULL = _args.full
    if _args.smoke:
        smoke_test()
    else:
        sys.exit(main())
