#!/usr/bin/env python3
"""Generate tests/data/betainv_{p,q}_reference.txt -- certified reference
set for corvus::beta_p_inv / corvus::beta_q_inv. Reuses the pinned
seed/step constants and module-level machinery from the checked-in
src/betainv_data.h and tools/gen_betainv_data.py (this generator does
NOT re-run that module's own replay/self-check pipeline, only consumes
its exported machinery and pinned constants).

ORACLE, three binding constructions:

  1. FAST-PATH forward evaluator for R1-tiny/joint-tiny certification
     traffic: bid.r1_value_mp -- a plain, self-convergent mpf power
     series at target dps (gen_betainv_data.py's own "measurement-grade
     truth" evaluator, NOT the fixed-N cheap routing proxy), bypassing
     gen_beta_reference.py's small_side_direct escalation ladder: that
     ladder ranges 0-5ms on ordinary points but 400-524ms on CF-heavy/
     ridge traffic -- infeasible at R1-tiny/joint-tiny volume. Validated
     fast-vs-full on a stratum sample (reported), then used for BOTH
     root-finding and bracket certification in the R1-tiny/joint-tiny
     strata specifically.

  2. GUARD on the reused gamma-corner route, at the enforcement site
     (this file, not gen_beta_reference.py/gen_beta_data.py):
     small_side_direct's own try_eval() calls gamma_corner_value(aa,bb,
     xx,dps) whenever max(aa,bb)>=B_GL, and gamma_corner_value ALWAYS
     feeds min(aa,bb)=min(a,b) to mpmath.gammainc as a shape argument
     (whichever of aa,bb is NOT picked as the huge "scale" side) --
     when BOTH a,b >= B_GL, that shape argument is itself huge (hang
     risk); gen_betainv_data.py's own betainv_forward is ALSO unsafe
     there (its R3 branch calls the raw backward CF unconditionally,
     which raises RuntimeError "CF not converged" at a=b=1e18,x=0.5 --
     not a hang, but not a value either). GUARD: whenever min(a,b) >=
     B_GL (kBetaGammaLim) AND the skew ratio max(a,b)/min(a,b) <=
     SKEW_SAFE_CAP, route through beta_temme_value() below -- a
     dual-anchored R3-Temme extraction built from gen_beta_data.py's
     own extract_e_monomial/r3_R_at machinery (gamma_ck itself is that
     module's GAMMA-side anchor cross-check; the R3 extraction
     apparatus it validates against is what this generator actually
     calls). This same route is also route-2 (the huge-nu stratum's
     independent second certification).
     SCOPE: the anchor-ladder extraction holds p=a/(a+b) FIXED across
     anchors at a MODEST nu (so alpha=nu/(1-p), beta=nu/p at each
     anchor stay CF-safe) -- this is safe exactly when the skew ratio
     is bounded (extreme skew forces an anchor parameter to blow up
     even at modest nu). This construction is scoped to
     "both-huge-BALANCED" traffic; SKEW_SAFE_CAP bounds how far from
     balanced this generator extends it, verified empirically (see
     both_huge_balanced) rather than assumed.

  3. Plateau rows: kappa computed per row (kappa = sigma/(y*f(y)), exact
     mpf, f = beta density). kappa <= 2^52 -> normal half-ulp bracket
     certification (the y-ULP gate). kappa > 2^52 -> BACKWARD-ERROR
     certification: forward of the STORED y at dps 100, required within
     a ~2-ulp-in-sigma contract -- no y-bracket exists to certify there
     (dd precision cannot resolve it). Deep-small rows (both
     orientations): log-space certification against the subnormal/zero
     boundary midpoints, with an exact dropped-term bound
     (bid._deep_small_dropped_rel) folded in as certification slack.

Everything else follows the same CERTIFICATION CORE as the gammainv
reference generator: root-find y* seeded by bid.seed_for, round to
double yd, certify sign(value-target) flips across the two half-ulp
midpoints of yd as exact mpf, layered dps 60 -> 100 (never lower). Rows
the certifier cannot prove: DECLINED and counted, not guessed. NEGATIVE
CONTROLS (>=4, 1-ULP-perturbed known-good rows) must be REJECTED on
every invocation, checked FIRST, exit 2 otherwise.

File format: five hex tokens per row: a b sigma yd marker. marker in
{N, P, B} (Normal bracket-gate row / Plateau backward-error row /
Beyond-resolution row) -- a marker token is used instead of separate
per-bucket files so the swap-identity orientation bookkeeping stays in
one place per side.

BINDING for any future ULP test consuming these files: the marker
column carries CERTIFICATION semantics ONLY (which construction/
contract proved the row -- bracket / backward-error / beyond-
resolution). It is NOT a routing or dilution-avoidance label. A ULP
test MUST bucket huge-nu statistics BY FORMULA computed from (a, b)
alone (the same predicate the kernel itself uses to select its huge-nu
path), independent of which marker a row happens to carry. In
particular: rows in the huge-nu collapse-ONSET band (nu ~ 1e32-1e35,
where the achievable y-transition has begun to narrow but has not yet
collapsed to <=1 ULP) are marked N (they pass ordinary bracket
certification) but are TRIVIALLY SATISFIABLE by nearly any kernel
answer in that neighbourhood -- diluting a same-bucket-as-N-elsewhere
ULP statistic the same way unbucketed beyond-resolution rows would
dilute an unbucketed gamma_inv ULP test. Existing N-marked rows in that
band are not relabeled -- a consuming test must derive its own bucket
boundary from (a,b) directly, using the huge_nu_beyond_resolution()
criterion below (or the kernel's own equivalent routing predicate)
rather than trusting the marker column to do that job.

Usage:
    python3 tools/gen_betainv_reference.py     # resumable; re-run until
                                                  # it reports DONE
"""
import hashlib
import math
import os
import random
import struct
import sys
import tempfile
import time

import mpmath as mp

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "tools"))
import gen_beta_data as gb          # noqa: E402  region cores, B_GL, T_RIDGE, ZETA_MAX
import gen_beta_reference as gbr    # noqa: E402  small_side_direct (audited oracle)
import gen_betainv_data as bid      # noqa: E402  seed_for, betainv_forward, r1_value_mp, ...
from refgen_common import round_to_double  # noqa: E402  single-rounding mpf->double (#13 N14)

SEED = 20260810  # fixed seed for reproducible point-set generation --
                  # changing it changes the RNG draw sequence, which
                  # invalidates any existing checkpoint (see CKPT_PATH,
                  # keyed off SEED) and forces a full regeneration

STATUS_PATH = os.path.join(
    r"C:\Users\gdwol\AppData\Local\Temp\claude\C--Users-gdwol-Development-corvus"
    r"\e81b05d8-c230-46b2-8caa-e48c35f168d2\scratchpad\betainv_g2b", "G2B-STATUS.md")


def status(line):
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    msg = f"- [{ts}] {line}\n"
    try:
        with open(STATUS_PATH, "a", encoding="utf-8") as f:
            f.write(msg)
    except OSError:
        pass
    print(f"STATUS: {line}", file=sys.stderr)


# ============================================================================
# Part 0: pinned thresholds. B_GL/T_RIDGE/ZETA_MAX come directly from the
# gb module's own (already-derived-at-import) values -- not re-derived.
# ============================================================================
B_GL = float(gb.B_GL)                    # kBetaGammaLim, ~2^59
SKEW_SAFE_CAP = 2.0e6                    # construction #2 scope bound -- covers
                                          # the beyond-resolution stratum's own
                                          # skewed corner (p~1e-6, skew~1e6);
                                          # validated by direct route1-vs-route2
                                          # agreement at skew=1e6/nu=1e20 (~30
                                          # decimal digits, ample for double-
                                          # precision certification).
HUGE_NU_THRESHOLD = 1.0e16               # route-2 dual-certification trigger (gammainv's own)
BEYOND_RESOLUTION_THRESHOLD = 3.0e34     # measured empirically (see huge_nu_beyond_resolution below)
PLATEAU_KAPPA_CUT = 2.0 ** 52            # contract split between bracket and backward-error certification
DEEP_SMALL_CUT = bid.DEEP_SMALL_CUT      # 2^-60, both orientations


GUARD_MIN_THRESHOLD = 1.0e10   # min(a,b)>=B_GL alone is NOT a sufficient
                                # guard boundary: small_side_direct's
                                # gamma_corner_value -> mp.gammainc can hang
                                # (no NoConvergence raised, no return within
                                # 90s under a hard-timeout probe) for a
                                # shape argument as "low" as 1.636e17 (well
                                # under B_GL~5.76e17) whenever the OTHER
                                # side is huge (>=B_GL) and the point sits
                                # near gammainc's OWN ridge/turning-point
                                # (a direct shape-only calibration up to
                                # 1e16 was fast when x was NOT near the
                                # ridge -- the hazard is ridge proximity,
                                # not magnitude alone, but this generator
                                # has no cheap way to predict ridge
                                # proximity ahead of the call, so the safe
                                # mitigation is routing more broadly through
                                # this generator's OWN extraction, which is
                                # exactly the ridge-region method).


def both_huge_balanced(a, b):
    """GUARD predicate (construction #2's enforcement condition): fires
    whenever the LARGER of (a,b) is >= B_GL (small_side_direct's own
    gamma_corner_value trigger) AND the smaller is >= GUARD_MIN_THRESHOLD
    (empirically far below where mp.gammainc's own ridge-proximity hazard
    was observed to hang) AND skew stays within the validated extraction
    range."""
    lo, hi = (a, b) if a <= b else (b, a)
    if hi < B_GL or lo < GUARD_MIN_THRESHOLD:
        return False
    return (hi / lo) <= SKEW_SAFE_CAP


# ============================================================================
# Part 1: construction #2 -- dual-anchored R3-Temme extraction. Reuses
# gb.extract_e_monomial + gb.r3_R_at VERBATIM (imported, not re-derived);
# only the anchor LADDERS (route 1 vs route 2) are this generator's own
# parameterization, independently chosen so route-1/route-2 agreement is
# a genuine cross-check, not a repeated computation.
# ============================================================================
def _ladder(route, n_anchor=None):
    if route == 1:
        # gb's OWN R3_NLADDER convention (T_RIDGE*2^(j/3)) -- the same
        # ladder the shipped kernel's R3 table extraction already uses
        # and gen_beta_data.py's own self-checks validated.
        n = n_anchor or 30
        return [float(gb.T_RIDGE) * 2.0 ** (j / 3.0) for j in range(n)]
    else:
        # INDEPENDENT second derivation: different base, different
        # spacing/growth, different count -- deliberately distinct from
        # route 1 so route1-vs-route2 agreement is a genuine cross-check,
        # not a repeated computation (the gammainv reference generator's
        # own route-2 precedent).
        n = n_anchor or 22
        return [45.0 * (1.6 ** j) for j in range(n)]


def beta_temme_value(a, b, x, dps, route=1, korder=None):
    """Dual-anchored R3-Temme reconstruction of (small_val, side) at
    ARBITRARY (a,b,x) with a,b possibly astronomically huge -- valid
    when both_huge_balanced(a,b) holds (anchor ladder shares the query's
    OWN p=a/(a+b), so anchor alpha/beta stay CF-safe only when the skew
    ratio is bounded -- see module docstring construction #2). korder
    defaults distinguish the two routes (15 vs 11 -- oversampled against
    the ladder length, matching gen_beta_data's own K_EXT-vs-ladder-
    length margin)."""
    if korder is None:
        korder = 15 if route == 1 else 11
    with mp.workdps(dps):
        a_m, b_m, x_m = mp.mpf(a), mp.mpf(b), mp.mpf(x)
        c = a_m + b_m
        nu = a_m * b_m / c
        p = a_m / c
        lam = a_m - c * x_m
        u = -lam / a_m
        v = lam / b_m
        cpsi = a_m * (u - mp.log1p(u)) + b_m * (v - mp.log1p(v))
        z = mp.sqrt(cpsi)
        zeta = z if lam >= 0 else -z
        leading = mp.erfc(z) / 2
        if cpsi > mp.mpf(800):
            # Saturated: the erfc leading term alone is exact to working
            # precision (R's own scale is a CORRECTION relative to the
            # exp(cpsi) amplification -- at cpsi>800 that correction is
            # already below dps=100 resolution, matching gamma's own
            # APHI_SAT=800 doctrine, reused here since ZETA_MAX's own
            # derivation used the identical cpsi<=800 bound).
            small_val = leading
        else:
            ladder = _ladder(route)
            e_poly, _ = gb.extract_e_monomial(ladder, p, zeta, korder, dps)
            S = mp.mpf(0)
            for k in range(korder - 1, -1, -1):
                S = S / nu + e_poly[k]
            W = mp.sqrt(2 * mp.pi * nu) * mp.exp(cpsi)
            if lam >= 0:
                small_val = leading - S / W
            else:
                small_val = leading + S / W
        side = "p" if lam >= 0 else "q"
        return small_val, side


def guarded_temme_pq(a, b, x, dps, route=1):
    """(P, Q) via beta_temme_value, either orientation."""
    small_val, side = beta_temme_value(a, b, x, dps, route=route)
    if side == "p":
        return small_val, 1 - small_val
    return 1 - small_val, small_val


# ============================================================================
# Part 2: forward evaluators dispatched by region -- FAST (root-finding,
# every stratum) vs FULL (bracket certification, guard-scoped per
# construction #1/#2). Every one of these is guard-checked FIRST so the
# both-huge-balanced hang/non-convergence never reaches small_side_direct
# or betainv_forward's raw CF.
# ============================================================================
def leading_only_pq(a, b, x, dps):
    """CHEAP guard-route approximation for BISECTION only: the erfc
    LEADING term alone (no R correction), a single mp.erfc call -- no
    ladder extraction. The FULL beta_temme_value (extract_e_monomial +
    QR solve, ~400-800ms per call near the ridge) is too expensive
    inside bisection's ~150-200-iteration inner loop. The R correction
    is a RELATIVE O(1/sqrt(nu)) effect (~1e-9 at nu~1e18); for
    root-finding purposes (not the final certified value) the leading
    term alone lands within a handful of ULPs of the true root, and
    certify_row's own local-nudge refinement (using the FULL evaluator,
    a bounded few calls) closes the gap before certifying."""
    with mp.workdps(dps):
        a_m, b_m, x_m = mp.mpf(a), mp.mpf(b), mp.mpf(x)
        c = a_m + b_m
        lam = a_m - c * x_m
        u = -lam / a_m
        v = lam / b_m
        cpsi = a_m * (u - mp.log1p(u)) + b_m * (v - mp.log1p(v))
        z = mp.sqrt(cpsi)
        leading = mp.erfc(z) / 2
        if lam >= 0:
            return leading, 1 - leading
        return 1 - leading, leading


def guarded_fast_forward(a, b, x, dps=60):
    """(P, Q) -- fast, used for BISECTION root-finding everywhere (many
    calls per row; must be cheap). Guard first, else gen_betainv_data's
    own router-based evaluator (R1 series / R2 CF / R3 CF-reconstruction
    / gamma-corner -- all fast, no escalation ladder). RESCUE (same
    defect as full_forward's own docstring -- bid.betainv_forward's
    gamma_corner_value_mp is equally unguarded against NoConvergence):
    on any exception, fall back to the erfc-leading-only approximation
    (cheap, always defined) rather than propagate -- root-finding only
    needs an APPROXIMATE value (certification is the accuracy gate), so
    this degrades gracefully instead of crashing the whole run."""
    if both_huge_balanced(a, b):
        return leading_only_pq(a, b, x, dps)
    try:
        val, side = bid.betainv_forward(a, b, x, dps=dps)
        if side == "p":
            return val, 1 - val
        return 1 - val, val
    except Exception:
        return leading_only_pq(a, b, x, dps)


def full_forward(a, b, x):
    """(P, Q, which, escalated, failed) -- mirrors small_side_direct's own
    signature exactly, so certify_row can treat both uniformly. Guard
    first (dps=100, matching the certifier's own top dps layer -- the
    extraction's cost is region-driven, not dps-driven, so paying dps=100
    unconditionally here is cheap for the common saturated case and only
    expensive on the already-rare near-transition rows).

    RESCUE WRAPPER (huge-nu scaling): small_side_direct's own try_eval
    calls gamma_corner_value -> mp.gammainc UNGUARDED against
    mp.libmp.libhyper.NoConvergence (it only catches RuntimeError/
    ZeroDivisionError/ValueError) -- a shape argument as 'low' as ~1e17
    (well UNDER kBetaGammaLim=2^59~5.76e17, outside this generator's own
    both_huge_balanced guard scope) can still make mpmath's hyp1f1
    series give up with NoConvergence, propagating as an UNCAUGHT
    exception out of small_side_direct itself. This is a robustness gap
    in the shipped, audited oracle that this generator cannot fix by
    editing gen_beta_reference.py (out of scope here) -- guarded
    EXTERNALLY instead: on NoConvergence (or any other exception
    small_side_direct doesn't itself catch), fall back to this
    generator's OWN dual-anchored extraction (construction #2's route)
    regardless of whether both_huge_balanced's own narrower threshold
    applies -- the extraction is mathematically valid at any nu, just
    more expensive, so this is a safe general rescue, not a correctness
    compromise."""
    if both_huge_balanced(a, b):
        P, Q = guarded_temme_pq(a, b, x, dps=100, route=1)
        which = "P" if P <= Q else "Q"
        return P, Q, which, False, False
    try:
        return gbr.small_side_direct(a, b, x)
    except Exception as e:
        status(f"RESCUE: small_side_direct raised {type(e).__name__}: {e} "
               f"at a={a!r} b={b!r} x={x!r} -- falling back to the "
               f"dual-anchored extraction (route1).")
        try:
            P, Q = guarded_temme_pq(a, b, x, dps=100, route=1)
            which = "P" if P <= Q else "Q"
            return P, Q, which, False, False
        except Exception:
            return None, None, "P", False, True


def fast_series_forward(a, b, x, dps, side):
    """Construction #1: plain self-convergent mpf series (bid.r1_value_mp)
    at target dps. DEPRECATED call shape (kept for the sample validator
    at line ~312, which always wants NATIVE P at a caller-guaranteed-
    small x) -- ROOT-FINDING call sites use fast_series_value/
    fast_series_pq below instead, which pick the orientation
    DYNAMICALLY by which argument is actually small."""
    return bid.r1_value_mp(a, b, x, dps)


def fast_series_value(a, b, y, dps):
    """(value, which) -- R1 series evaluated at whichever orientation has
    the SMALL 3rd argument (fast, reliable convergence): native (a,b,y)
    if y<=0.5 (value=P), swapped (b,a,1-y) if y>0.5 (value=Q). This
    orientation must be chosen DYNAMICALLY by which argument is actually
    small, not by which of P/Q the caller wants: for joint-tiny rows
    whose true y sits on the far side (e.g. solving P(a,b,y)=sigma with
    the true y near 1), a fixed-by-side orientation forces
    bid.r1_value_mp to converge a series in x close to 1, which is slow
    (multi-second) and for some points never converges within the
    4000-term cap at all -- silently missing the sign flip, which
    manifests as a spurious root-find-failed decline rather than merely
    a speed problem."""
    if y <= 0.5:
        return bid.r1_value_mp(a, b, y, dps), "p"
    return bid.r1_value_mp(b, a, 1 - y, dps), "q"


def fast_series_pq(a, b, y, dps):
    """(P, Q) via fast_series_value, complementing whichever side wasn't
    computed directly (the complement is a plain 1-v in mpf -- both
    values share the same tiny-side-direct guarantee small_side_direct
    itself relies on, since exactly one of P,Q is ever computed via a
    slow-converging orientation and it is never the one returned as
    'value')."""
    v, which = fast_series_value(a, b, y, dps)
    if which == "p":
        return v, 1 - v
    return 1 - v, v


def fast_vs_full_validate(rng, n=40):
    """Validates fast-vs-full agreement on a stratum sample (reports the
    sample size and worst disagreement) per construction #1's mandate.
    Samples R1-tiny/joint-tiny-shaped points, compares
    bid.r1_value_mp (fast) against gbr.small_side_direct (full, audited)
    at dps=80, reports worst relative disagreement."""
    worst = 0.0
    worst_at = None
    n_done = 0
    t0 = time.time()
    for _ in range(n):
        a = 10.0 ** rng.uniform(-6, 2)
        b = 10.0 ** rng.uniform(-6, 2)
        x = 10.0 ** rng.uniform(-300, -1)
        try:
            P, Q, which, esc, failed = gbr.small_side_direct(a, b, x)
        except Exception:
            continue
        if failed:
            continue
        true_small = P if which == "P" else Q
        try:
            fast_p = fast_series_forward(a, b, x, dps=80, side="p")
        except Exception:
            continue
        fast_small = fast_p if which == "P" else (1 - fast_p if which == "Q" else None)
        # fast_series_forward always returns NATIVE P(a,b,x); compare
        # against whichever of P/Q the audited oracle called "small".
        got = fast_p if which == "P" else (1 - fast_p)
        ref = true_small
        rel = float(abs((got - ref) / ref)) if ref != 0 else float(abs(got))
        n_done += 1
        if rel > worst:
            worst, worst_at = rel, (a, b, x)
    dt = time.time() - t0
    status(f"construction #1 fast-vs-full validation: n={n_done}/{n} sampled "
           f"({dt:.1f}s), worst relative disagreement {worst:.3e} at {worst_at}")
    return n_done, worst, worst_at


# ============================================================================
# Part 3: root-finder -- bisection in logit(y) space (linear-space
# bisection is unusable near y=0/1 -- gen_betainv_data hit this same
# hazard; avoided here from the start), seeded via bid.seed_for for
# speed only (never for correctness -- falls back to the full default
# bracket whenever the seeded one fails to bracket).
# ============================================================================
def oracle_y(a, b, target, side, dps, use_fast_series=False, seed_hint=None,
             root_dps=None, bracket_halfwidth=80.0, seed_only=False):
    """bracket_halfwidth/seed_only: generic bisection-hardening params --
    the logit-space half-width around seed_hint (default 80) and, when
    seed_only=True, a HARD REFUSAL to fall back to the wide/expanding
    ((-2000,2000), then +-200-per-try) global search when the seeded
    bracket doesn't contain a sign change (an honest decline instead of
    letting bisection wander to a numerically-plausible-but-WRONG point
    far from the seed). For the gamma-limit-seam stratum, a tighter
    bracket around this function's own (cheap, misrouting-prone at
    extreme skew) evaluator is not enough -- see oracle_y_audited below,
    which is what certify_row's gamma-limit-seam branch calls instead.
    These params stay as general-purpose bisection hardening for any
    future caller that only needs a tighter/non-wandering search against
    THIS function's cheap evaluator.

    root_dps: root-FINDING precision, decoupled from dps (the CALLER's
    certification precision). Running root-finding's ~150-190-iteration
    bisection loop AT THE FULL CERTIFICATION dps is expensive: for
    fast_series (bid.r1_value_mp), series truncation eps scales with dps
    (eps=10^-(dps-8)), so every one of those ~180 series evaluations
    would pay dps=60's full term count, including the many iterations
    that land far from the true root (interior y, not small x -- where
    the series is genuinely slow, needing up to its 4000-term cap;
    measured outlier 9.1s/row at dps=60). Root-find instead at a much
    lower dps (default 28, ~1e-20-class eps -- ample to land within a
    handful of ULPs of the true double, since bisection's own stopping
    tolerance is dps-scaled too), THEN certify at the caller's real
    dps_layers (60->100, a handful of calls, not ~180) -- the same
    guarded_fast_forward/full_forward cheap-root-find-vs-expensive-
    certify split used elsewhere in this generator, applied to
    fast_series too."""
    if root_dps is None:
        root_dps = dps
    a_m, b_m = mp.mpf(a), mp.mpf(b)
    with mp.workdps(root_dps):
        target_m = mp.mpf(target)

        def f(v):
            y = bid.sigmoid(v)
            if y <= 0:
                y = mp.mpf(2) ** -1075
            elif y >= 1:
                y = 1 - mp.mpf(2) ** -1075
            if use_fast_series:
                # dynamic orientation by which argument is small (see
                # fast_series_pq's own docstring for why -- fixed-by-
                # side orientation silently missed the bracket entirely
                # for far-side roots).
                P, Q = fast_series_pq(a_m, b_m, y, root_dps)
                got = P if side == "p" else Q
            else:
                P, Q = guarded_fast_forward(a_m, b_m, y, dps=min(root_dps, 60))
                got = P if side == "p" else Q
            return got - target_m

        lo = mp.mpf(-2000)
        hi = mp.mpf(2000)
        seeded_ok = False
        if seed_hint is not None and math.isfinite(seed_hint) and 0.0 < seed_hint < 1.0:
            v_seed = bid.logit(mp.mpf(seed_hint))
            slo, shi = v_seed - bracket_halfwidth, v_seed + bracket_halfwidth
            try:
                if f(slo) * f(shi) <= 0:
                    lo, hi = slo, shi
                    seeded_ok = True
            except (ValueError, OverflowError):
                pass
        if seed_only and not seeded_ok:
            return None
        flo, fhi = f(lo), f(hi)
        tries = 0
        while (not seed_only) and flo * fhi > 0 and tries < 40:
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
        n_iters = min(250, max(80, root_dps * 3))
        for _ in range(n_iters):
            mid = (lo + hi) / 2
            fm = f(mid)
            if (fm > 0) == (flo > 0):
                lo, flo = mid, fm
            else:
                hi, fhi = mid, fm
            if hi - lo < mp.mpf(10) ** (-(root_dps - 5)):
                break
        return bid.sigmoid((lo + hi) / 2)


# ============================================================================
# Part 4: deep-small closed-form certification (construction #3, both
# orientations), reusing: bid.deep_small_cut_bound (routing decision,
# native double), bid.deep_small_y (candidate y0, native double),
# bid._deep_small_dropped_rel (EXACT mpf dropped-term bound, folded in
# as certification slack, the gammainv pattern).
# ============================================================================
def deep_small_ly0(a, b, sigma, side, dps):
    """log-space target + exact bound, either orientation. Returns
    (ln_target, yd_mpf_rounded, bound_mpf) where ln_target is ln(y) [p]
    or ln(1-y) [q] at the closed-form leading order, and
    yd_mpf_rounded is the double CORRECTLY ROUNDED from ln_target at
    full working dps. bid.deep_small_y's own NATIVE-DOUBLE seed formula
    must NOT be used here directly: that function is a fast SEED
    candidate for the kernel's own root search, built from several
    chained double-precision log/exp calls, and its compounded rounding
    error is enough to land OUTSIDE the half-ulp bracket around the true
    high-precision value, causing every deep-small row to spuriously
    fail certification. The oracle's own yd must come from rounding the
    HIGH-PRECISION ln_t, exactly the gammainv deep_small_lx0 pattern --
    native-double formulas are seeds, never the certified value."""
    with mp.workdps(dps):
        a_m, b_m, sigma_m = mp.mpf(a), mp.mpf(b), mp.mpf(sigma)
        if side == "p":
            lnB = mp.loggamma(a_m) + mp.loggamma(b_m) - mp.loggamma(a_m + b_m)
            ln_t = (mp.log(sigma_m) + mp.log(a_m) + lnB) / a_m
        else:
            lnB = mp.loggamma(b_m) + mp.loggamma(a_m) - mp.loggamma(a_m + b_m)
            ln_t = (mp.log(sigma_m) + mp.log(b_m) + lnB) / b_m
        y0_mpf = mp.e ** ln_t
        # round_to_double (#13 N14): y0_mpf can land in SUBNORMAL territory
        # (this IS the deep-small closed form's own target range) --
        # float(mpf) there double-rounds (53-bit mantissa round, then a
        # second ldexp round onto the 2^-1074 grid), up to 1 ulp off
        # correctly rounded. Single-rounding matters here specifically
        # because this yd_rounded IS the certified value (yd_override
        # aside), not a seed.
        if side == "p":
            yd_rounded = round_to_double(y0_mpf) if y0_mpf > 0 else 0.0
        else:
            yd_rounded = 1.0 - round_to_double(y0_mpf) if y0_mpf > 0 else 1.0
        bound = bid._deep_small_dropped_rel(a_m, b_m, y0_mpf, side, dps=dps)
        return ln_t, yd_rounded, bound


MIN_SUBNORMAL = math.ldexp(1.0, -1074)


def certify_deep_small(a, b, sigma, side, dps, yd_override=None):
    """Certify a SPECIFIC double yd against the closed-form log-space
    target, at half-ulp midpoints, bound folded in as slack (gammainv
    pattern). yd_override must be actually threaded through here, or
    negative controls on deep-small rows are silently never checked."""
    ln_t, yd_from_ln, bound = deep_small_ly0(a, b, sigma, side, dps)
    yd = yd_override if yd_override is not None else yd_from_ln
    if not math.isfinite(yd):
        return {"yd": None, "certified": False, "note": "non-finite"}
    yp = yd if side == "p" else (1.0 - yd)
    with mp.workdps(dps):
        if yp <= 0.0:
            # p-side: yd rounds to 0.0 (round-to-zero boundary, compare
            # against MIN_SUBNORMAL/2). q-side: yd rounds to 1.0 (the
            # SMALL quantity z=1-y rounds below ulp(1.0)/2=2^-53) -- the
            # q-side boundary differs from the p-side one by ~700
            # decades, so the two must NOT be conflated.
            boundary = (mp.mpf(MIN_SUBNORMAL) / 2 if side == "p"
                        else mp.mpf(2) ** -53)
            certified = (ln_t + bound) < mp.log(boundary)
            return {"yd": yd, "certified": bool(certified),
                    "note": "round-to-zero" if side == "p" else "round-to-one"}
        lo_d = math.nextafter(yp, -math.inf)
        hi_d = math.nextafter(yp, math.inf)
        yp_m = mp.mpf(yp)
        lo_mid = (yp_m + mp.mpf(lo_d)) / 2 if lo_d > 0 else mp.mpf(lo_d) / 2
        hi_mid = (yp_m + mp.mpf(hi_d)) / 2
        ln_lo = mp.log(lo_mid) if lo_mid > 0 else mp.mpf(-1e30)
        ln_hi = mp.log(hi_mid)
        certified = (ln_lo < ln_t - bound) and (ln_t + bound < ln_hi)
        return {"yd": yd, "certified": bool(certified), "ln_t": ln_t, "bound": bound,
                "ln_lo": ln_lo, "ln_hi": ln_hi}


# ============================================================================
# Part 5: standard bracket certification (half-ulp midpoints of yd, sign
# flip of forward-target). fwd_kind selects the evaluator: "fast_series"
# (construction #1, R1-tiny/joint-tiny), "guard" (construction #2,
# both-huge-balanced), "full" (default, gbr.small_side_direct).
# ============================================================================
def _sign(v):
    return 1 if v > 0 else (-1 if v < 0 else 0)


def bracket_signs(a, b, target, side, yd, dps, fwd_kind):
    lo_d = math.nextafter(yd, -math.inf)
    hi_d = math.nextafter(yd, math.inf)
    with mp.workdps(dps):
        yd_m = mp.mpf(yd)
        lo_mid = (yd_m + mp.mpf(lo_d)) / 2 if lo_d > 0 else mp.mpf(lo_d) / 2
        hi_mid = (yd_m + mp.mpf(hi_d)) / 2
        t = mp.mpf(target)

        def val_at(y):
            if fwd_kind == "fast_series":
                P, Q = fast_series_pq(a, b, y, dps)
                return P if side == "p" else Q
            if fwd_kind == "guard":
                P, Q = guarded_temme_pq(a, b, y, dps, route=1)
                return P if side == "p" else Q
            if fwd_kind == "route2":
                P, Q = guarded_temme_pq(a, b, y, dps, route=2)
                return P if side == "p" else Q
            P, Q, which, esc, failed = full_forward(a, b, y)
            if failed:
                return None
            return P if side == "p" else Q

        vlo = val_at(lo_mid)
        vhi = val_at(hi_mid)
        if vlo is None or vhi is None:
            return None, None
        flo = vlo - t
        fhi = vhi - t
    return _sign(flo), _sign(fhi)


def certify_bracket(a, b, target, side, yd, dps_layers, fwd_kind):
    layer_results = []
    for dps in dps_layers:
        slo, shi = bracket_signs(a, b, target, side, yd, dps, fwd_kind)
        if slo is None:
            layer_results.append(False)
            continue
        ok = (slo != shi and slo != 0 and shi != 0) or slo == 0 or shi == 0
        layer_results.append(ok)
    return all(layer_results), layer_results


# ============================================================================
# Part 6: plateau backward-error certification (construction #3, kappa >
# 2^52 branch) -- forward of the STORED yd at dps=100 must land within
# a ~2-ulp-in-sigma contract: |forward(yd) - sigma|
# <= 2 * ulp(sigma) in the solved side's own probability space.
# ============================================================================
def ulp_of(v):
    if v == 0.0:
        return MIN_SUBNORMAL
    return math.ulp(v)


def certify_plateau(a, b, sigma, side, yd, dps=100):
    with mp.workdps(dps):
        P, Q, which, esc, failed = full_forward(a, b, yd)
        if failed:
            return False
        got = P if side == "p" else Q
        diff = abs(float(got) - sigma)
    return diff <= 2.0 * ulp_of(sigma)


def compute_kappa(a, b, yd, sigma, dps=60):
    """kappa = sigma/(y*f(y)), exact mpf (the plateau formula).
    beta_density_mp's log-space assembly
    lg=(a-1)*ln(y)+(b-1)*ln(1-y)-lnB(a,b) is a CANCELLATION of terms
    scaling with a,b themselves (~1e35 at the huge-nu stratum's own
    scale) down to an O(1)-ish true log-density -- at dps=60 (60
    significant decimal digits) that cancellation loses ALL its
    resolving digits once a,b exceed ~1e40ish, and the naive symptom is
    f rounding to EXACTLY 0, kappa->inf, misrouting a huge-nu row
    (which should have kappa NEAR ZERO -- density is enormous, not
    tiny, at that scale) into the plateau-backward branch (wrong
    contract; the row is really either an ordinary bracket
    certification or a genuine beyond-resolution one). FIX: scale dps
    with max(a,b)'s own decimal magnitude so the cancellation always
    has working digits left over, independent of the caller's
    certification-layer dps."""
    scale_dps = 60 + int(math.log10(max(float(a), float(b), 10.0))) + 20
    use_dps = max(dps, scale_dps)
    with mp.workdps(use_dps):
        f = bid.beta_density_mp(a, b, yd, use_dps)
        if f == 0:
            return mp.mpf("inf")
        return mp.mpf(sigma) / (mp.mpf(yd) * f)


# ============================================================================
# Part 7: beyond-resolution certification (nearest-neighbor, escalated
# dps) -- the gammainv reference generator's own construction, reused
# for the huge-nu stratum where the transition collapses below 1 ULP of y.
# ============================================================================
def certify_beyond_resolution(a, b, target, side, yd, dps, fwd_kind):
    lo_d = math.nextafter(yd, -math.inf)
    hi_d = math.nextafter(yd, math.inf)
    with mp.workdps(dps):
        t = mp.mpf(target)

        def val_at(y):
            if fwd_kind == "guard":
                P, Q = guarded_temme_pq(a, b, y, dps, route=1)
                return P if side == "p" else Q
            if fwd_kind == "route2":
                P, Q = guarded_temme_pq(a, b, y, dps, route=2)
                return P if side == "p" else Q
            P, Q, which, esc, failed = full_forward(a, b, y)
            if failed:
                return None
            return P if side == "p" else Q

        vl, vx, vh = val_at(lo_d), val_at(yd), val_at(hi_d)
        if vl is None or vx is None or vh is None:
            return False
        el, ex, eh = abs(vl - t), abs(vx - t), abs(vh - t)
    return bool(ex <= el and ex <= eh)


# ============================================================================
# Part 8: certify_row -- the top-level dispatcher tying all three
# constructions together, per the CERTIFICATION CORE doctrine (gammainv
# pattern): seed -> root-find -> round -> half-ulp bracket certify at
# layered dps 60->100, with construction-specific routing.
# ============================================================================
N_BEYOND_RESOLUTION = [0]
N_ROUTE2_CHECKED = [0]
N_ROUTE2_DISAGREE = [0]


def refine_guard_root(a, b, target, side, yd, dps=60, max_nudge=12):
    """Bounded local ULP-nudge refinement for the both-huge-balanced
    guard route: bisection's own seed there is the CHEAP leading-erfc-
    only approximation (construction #2's own cost fix, see
    leading_only_pq) -- the dropped R correction is O(1/sqrt(nu))
    relative in P, and at nu>=B_GL~5.76e17 that maps to a Y-space error
    orders of magnitude BELOW 1 ULP (density there scales ~sqrt(nu), so
    Delta-y ~ Delta-P/sqrt(nu) ~ nu^-1 relative to ulp(y)~2^-53 --
    verified empirically), so this search is a defensive bound, not the
    expected common path."""
    with mp.workdps(dps):
        t = mp.mpf(target)

        def f(y):
            P, Q = guarded_temme_pq(a, b, y, dps, route=1)
            return float((P if side == "p" else Q) - t)

        best_y, best_abs = yd, abs(f(yd))
        if best_abs == 0.0:
            return yd
        for direction in (1, -1):
            cur = yd
            fcur = best_abs if direction == 1 else best_abs  # recompute below
            fcur = f(cur)
            for _ in range(max_nudge):
                nxt = math.nextafter(cur, math.inf if direction > 0 else -math.inf)
                fnxt = f(nxt)
                if abs(fnxt) < best_abs:
                    best_y, best_abs = nxt, abs(fnxt)
                if (fnxt > 0) != (fcur > 0):
                    break
                if abs(fnxt) > abs(fcur):
                    break
                cur, fcur = nxt, fnxt
        return best_y


def oracle_y_audited(a, b, target, side, seed_hint, dps=60,
                      bracket_halfwidth=60.0, max_expand=6, n_iters=90):
    """Certifies the gamma-limit-seam stratum by bisecting directly
    against the AUDITED evaluator, replacing an insufficient earlier
    approach (a tight bracket around the CHEAP evaluator, kept as
    oracle_y's seed_only/bracket_halfwidth params -- still useful
    generic machinery, just not sufficient for this stratum alone).

    ROOT CAUSE: the CHEAP evaluator oracle_y bisects against
    (guarded_fast_forward -> bid.betainv_forward -> gb.route_final)
    silently MISROUTES a neighbourhood just past the true root into
    gb.route_final's 'R3-native' tag, whose evaluator
    (gb.small_val_via_cf) is NOT valid at this skew (one shape param
    ~20-300, the other 1e100-1e250) and returns GARBAGE with no
    exception -- witness a=41.4216 b=1.69580e+111: at
    y=6.140121697366854e-111 (tag R2-native-gammalim, the CORRECT
    route) small_val_via_cf-chain gives P=3.470392478062253e-13,
    matching sigma=3.4751741894519584e-13 to ~1e-3 relative; at
    y=2.4425909667268353e-110 (only ~1.38 logit-units away -- INSIDE
    any bracket tolerant of the seed's own honest uncertainty) tag
    flips to R3-native and the same call chain returns
    1.0045e-4541 against a TRUE audited value of 0.5207
    (gbr.small_side_direct) -- catastrophically, silently wrong. A
    tight bracket around a good seed cannot fix this: the false
    crossing sits INSIDE any bracket wide enough to tolerate the
    seed's own uncertainty (the seed is measured 15-52 bits accurate,
    NOT ULP-accurate -- a real bisection is still required, just not
    against a function that lies).

    FIX: bisect directly against the AUDITED evaluator (full_forward,
    the SAME one certify_bracket checks against below) rather than the
    cheap misrouting one -- no second untrustworthy function is ever
    consulted. Bracket widens in BOUNDED doublings around the seed
    (never unbounded -- a genuinely unprovable point should DECLINE,
    not wander); gamma-limit-seam's single-huge-parameter shape keeps
    full_forward itself cheap here (this is the intended fast native
    path for small_side_direct, not its CF-heavy/ridge cost class)."""
    if seed_hint is None or not (math.isfinite(seed_hint) and 0.0 < seed_hint < 1.0):
        return None
    with mp.workdps(dps):
        t = mp.mpf(target)
        v_seed = bid.logit(mp.mpf(seed_hint))

        def f(v):
            y = bid.sigmoid(v)
            if y <= 0:
                y = mp.mpf(2) ** -1075
            elif y >= 1:
                y = 1 - mp.mpf(2) ** -1075
            P, Q, which, esc, failed = full_forward(a, b, y)
            if failed:
                return None
            got = P if side == "p" else Q
            return got - t

        hw = mp.mpf(bracket_halfwidth)
        lo = hi = None
        flo = fhi = None
        for _ in range(max_expand):
            slo, shi = v_seed - hw, v_seed + hw
            flo, fhi = f(slo), f(shi)
            if flo is None or fhi is None:
                return None
            if flo * fhi <= 0:
                lo, hi = slo, shi
                break
            hw *= 2
        if lo is None:
            return None
        for _ in range(n_iters):
            mid = (lo + hi) / 2
            fm = f(mid)
            if fm is None:
                return None
            if (fm > 0) == (flo > 0):
                lo, flo = mid, fm
            else:
                hi, fhi = mid, fm
            if hi - lo < mp.mpf(10) ** (-(dps - 5)):
                break
        # round_to_double (#13 N14): this IS the certified double the
        # gamma-limit-seam stratum ships -- single rounding, not float()'s
        # potential double-rounding in the subnormal band.
        return round_to_double(bid.sigmoid((lo + hi) / 2))


def certify_row(a, b, sigma, side, tag, dps_layers=(60, 100), yd_override=None):
    is_r1tiny = tag in ("r1-tiny", "joint-tiny", "r1-tiny-seam")
    huge_bal = both_huge_balanced(a, b)
    fwd_kind = "fast_series" if is_r1tiny else ("guard" if huge_bal else "full")

    # --- deep-small routing decision (BOTH orientations) -- cheap
    # native-double check first. yd_override (negative
    # controls) bypasses the routing decision -- the override IS the
    # hand-perturbed candidate to certify/reject regardless of which
    # bucket it would naturally fall in.
    try:
        y0_native = bid.deep_small_y(a, b, sigma, side)
        db = bid.deep_small_cut_bound(a, b, y0_native, side)
    except (OverflowError, ValueError):
        y0_native, db = None, math.inf
    if db < DEEP_SMALL_CUT and yd_override is None:
        layer_results = [certify_deep_small(a, b, sigma, side, dps)["certified"]
                          for dps in dps_layers]
        yd = certify_deep_small(a, b, sigma, side, dps_layers[-1])["yd"]
        return {"yd": yd, "certified": all(layer_results), "method": "deep-small",
                "marker": "N", "kappa": None}
    if db < DEEP_SMALL_CUT and yd_override is not None:
        layer_results = [certify_deep_small(a, b, sigma, side, dps,
                                             yd_override=yd_override)["certified"]
                          for dps in dps_layers]
        return {"yd": yd_override, "certified": all(layer_results), "method": "deep-small",
                "marker": "N", "kappa": None}

    if yd_override is not None:
        y_star = mp.mpf(yd_override)
    elif tag == "gamma-limit-seam":
        # Seed with seed_S3 (gamma-transfer -- it owns this territory)
        # and root-find via oracle_y_audited -- see that function's own
        # docstring for the root cause (the cheap evaluator's routing
        # silently returns garbage near this stratum's own seeds) and
        # the fix (bisect directly against the SAME audited evaluator
        # certification checks against, bounded bracket doubling around
        # the seed, never an unbounded search).
        try:
            seed = bid.seed_S3(a, b, sigma, side)
        except Exception:
            seed = None
        if seed is None or not (math.isfinite(seed) and 0.0 < seed < 1.0):
            try:
                seed = bid.seed_for(a, b, sigma, side)
            except Exception:
                seed = None
        y_star = oracle_y_audited(a, b, sigma, side, seed)
        if y_star is not None:
            y_star = mp.mpf(y_star)
    else:
        # --- seed + root-find ---
        try:
            seed = bid.seed_for(a, b, sigma, side)
        except Exception:
            seed = None
        y_star = oracle_y(a, b, sigma, side, dps=dps_layers[0],
                           use_fast_series=is_r1tiny, seed_hint=seed,
                           root_dps=(40 if is_r1tiny else None))
    if y_star is None:
        return {"yd": None, "certified": False, "method": "root-find-failed",
                "marker": "N", "kappa": None}
    # round_to_double (#13 N14): the oracle root becomes the double output
    # here -- y_star can land in subnormal territory (underflow-threshold/
    # subnormal-y strata target exactly that range), where float(mpf)
    # double-rounds. Single rounding, matching deep_small_ly0/
    # oracle_y_audited's own sites.
    yd = round_to_double(y_star)
    if not (math.isfinite(yd) and 0.0 <= yd <= 1.0):
        return {"yd": yd, "certified": False, "method": "non-finite-yd",
                "marker": "N", "kappa": None}
    if yd in (0.0, 1.0):
        # boundary collapse -- treat as deep-small-zero on the appropriate
        # side (mirrors gammainv's own xd==0 branch).
        layer_results = [certify_deep_small(a, b, sigma, side, dps)["certified"]
                          for dps in dps_layers]
        return {"yd": yd, "certified": all(layer_results), "method": "deep-small-boundary",
                "marker": "N", "kappa": None}

    # --- guard-route refinement (bounded, see refine_guard_root) --
    # SKIPPED when yd_override is given (negative controls): refinement
    # would walk a deliberately-wrong candidate back toward the correct
    # root, defeating the control.
    if fwd_kind == "guard" and yd_override is None:
        yd = refine_guard_root(a, b, sigma, side, yd, dps=dps_layers[0])

    # --- kappa / plateau routing ---
    kappa = compute_kappa(a, b, yd, sigma, dps=dps_layers[0])
    if kappa > PLATEAU_KAPPA_CUT:
        ok = certify_plateau(a, b, sigma, side, yd, dps=100)
        if ok:
            return {"yd": yd, "certified": True, "method": "plateau-backward",
                     "marker": "P", "kappa": float(kappa)}
        # FALL THROUGH (robustness, huge-nu escalation): a huge-kappa
        # reading that fails ITS OWN backward-error contract is not
        # necessarily a genuinely unprovable row -- at extreme a,b the
        # density cancellation compute_kappa now dps-scales against can
        # still leave a borderline/wrong kappa classification in rare
        # cases; the ordinary bracket path (with its own beyond-
        # resolution fallback below) is a strictly more general
        # certifier, so give it a chance rather than declining outright
        # on a single failed contract.

    # --- standard bracket certification ---
    ok, layer_results = certify_bracket(a, b, sigma, side, yd, dps_layers, fwd_kind)
    marker = "N"
    method = "bracket"
    if not ok and min(a, b) >= 1.0e15:
        for esc_dps in (150, 220):
            if certify_beyond_resolution(a, b, sigma, side, yd, esc_dps, fwd_kind):
                ok = True
                marker = "B"
                method = "beyond-resolution"
                N_BEYOND_RESOLUTION[0] += 1
                break

    route2_ok = None
    if min(a, b) >= HUGE_NU_THRESHOLD:
        N_ROUTE2_CHECKED[0] += 1
        ok2, _ = certify_bracket(a, b, sigma, side, yd, (100,), "route2")
        if not ok2:
            ok2 = certify_beyond_resolution(a, b, sigma, side, yd, 150, "route2")
        route2_ok = ok2
        if not ok2:
            N_ROUTE2_DISAGREE[0] += 1

    certified = ok and (route2_ok is None or route2_ok)
    return {"yd": yd, "certified": bool(certified), "method": method,
            "marker": marker, "kappa": float(kappa), "route2_ok": route2_ok}


# ============================================================================
# Part 9: negative controls (>=4, spanning normal/deep-small/plateau/
# huge-nu; MUST be REJECTED on every invocation; checked FIRST, exit 2
# otherwise -- the beta/gammainv doctrine).
# ============================================================================
def s_from_y(a, b, y, side, dps=80, tag="moderate"):
    """Construct a well-posed sigma = forward(a,b,y) rounded to double
    (gammainv's own s_from_x pattern) -- guarantees the (a,b,sigma,side)
    triple is REACHABLE, sidestepping the ill-posed-guess pathologies a
    blind log-uniform sigma draw hits (an unreachable (a,b,sigma) triple
    silently fails root-finding)."""
    try:
        if tag in ("r1-tiny", "joint-tiny", "r1-tiny-seam"):
            P, Q = fast_series_pq(a, b, y, dps)
            v = P if side == "p" else Q
        else:
            P, Q = guarded_fast_forward(a, b, y, dps=min(dps, 60))
            v = P if side == "p" else Q
    except Exception:
        return None
    if v is None or not mp.isfinite(v):
        return None
    sd = float(v)
    if not (math.isfinite(sd) and 0.0 < sd < 1.0):
        return None
    return sd


def negative_controls():
    print("negative controls ...", file=sys.stderr)
    status("negative controls: starting")
    cases = []
    # normal bracket-gate row
    s = s_from_y(5.0, 3.0, 0.4, "p")
    if s is not None:
        cases.append(("normal", 5.0, 3.0, s, "p", "moderate", 1))
    # deep-small row (both orientations) -- sigma picked DIRECTLY (a
    # representable double), since the deep-small route computes yd from
    # sigma via the closed form, not the other way around (a forward-
    # from-y construction underflows sigma to exactly 0.0 here --
    # sigma~1e-100 maps to y~1e-500-scale internal to the closed form,
    # far below double range, which is exactly deep-small's own point).
    cases.append(("deep-small-p", 5.0, 3.0, 1e-100, "p", "moderate", 1))
    cases.append(("deep-small-q", 3.0, 5.0, 1e-100, "q", "moderate", 1))
    # plateau row (large perturbation -- the backward-error contract
    # tolerates ~2 ulp by design, so a 1-ulp control would not be a
    # valid negative control there; use a coarse, clearly-outside
    # perturbation instead).
    s = s_from_y(1e-8, 1.5e-8, 0.4, "p")
    if s is not None:
        cases.append(("plateau", 1e-8, 1.5e-8, s, "p", "moderate", 4096))
    # huge-nu / both-huge-balanced row
    s = s_from_y(1e18, 1e18, 0.5 + 1e-9, "p")
    if s is not None:
        cases.append(("huge-nu", 1e18, 1e18, s, "p", "moderate", 1))

    if len(cases) < 4:
        print(f"  FATAL: could only construct {len(cases)}/4+ good negative-control "
              f"rows (forward construction failed on the rest)", file=sys.stderr)
        return False

    all_rejected = True
    for name, a, b, s, side, tag, nulp in cases:
        good = certify_row(a, b, s, side, tag)
        if not good["certified"] or good["yd"] is None:
            print(f"  FATAL: could not even certify the GOOD row for control "
                  f"{name}: {good}", file=sys.stderr)
            return False
        bad_yd = good["yd"]
        for _ in range(nulp):
            bad_yd = math.nextafter(bad_yd, math.inf)
        r = certify_row(a, b, s, side, tag, yd_override=bad_yd)
        status_str = "REJECTED (correct)" if not r["certified"] else "ACCEPTED (FATAL BUG)"
        print(f"  [{name}] a={a:.6e} b={b:.6e} s={s:.6e} good_yd={good['yd']!r} "
              f"bad_yd={bad_yd!r} ({nulp} ulp) -> {status_str}", file=sys.stderr)
        status(f"negative control [{name}]: {status_str}")
        if r["certified"]:
            all_rejected = False
    return all_rejected


# ============================================================================
# Part 10: strata generation. Every point is (a, b, sigma, side, tag).
# sigma is constructed FROM a chosen y via forward evaluation wherever the
# region is hard (R1-tiny/joint-tiny/ridge/gammalim/huge-nu/seams) --
# gammainv's own s_from_x pattern, adopted here after log-uniform sigma
# sampling in these regions was found to hit ill-posed-guess pathologies.
# The swap identity halves orientation coverage (one of (a,b)/(b,a) per
# logical point) EXCEPT near-diagonal/plateau rows, where both are kept.
# ============================================================================
NEXT_UP = lambda v: math.nextafter(v, math.inf)
NEXT_DN = lambda v: math.nextafter(v, -math.inf)


class PointSet:
    def __init__(self):
        self.seen = set()
        self.pts = []
        self.strata_counts = {}

    def add(self, a, b, sigma, side, tag):
        if not (math.isfinite(a) and math.isfinite(b) and math.isfinite(sigma)
                and a > 0 and b > 0 and 0.0 < sigma < 1.0):
            return
        key = (a, b, sigma, side)
        if key in self.seen:
            return
        self.seen.add(key)
        self.pts.append((a, b, sigma, side, tag))
        self.strata_counts[tag] = self.strata_counts.get(tag, 0) + 1

    def add_from_y(self, a, b, y, side, tag, dps=60):
        s = s_from_y(a, b, y, side, dps=dps, tag=tag)
        if s is not None:
            self.add(a, b, s, side, tag)

    def add_from_y_smallside(self, a, b, y, tag, dps=60, both=False):
        """Constructs sigma from y, choosing the SIDE by which of P/Q is
        actually <=0.5 (the 'solve against s=min(p,1-p)' input contract)
        rather than a random coin flip -- a random side with y not
        already known to be on that side's small end can construct an
        out-of-contract row (sigma>0.5 for the requested side), which
        the certifier's deep-small/boundary machinery is not designed to
        handle (it assumes the REQUESTED side's sigma is the one
        approaching 0, not 1) and which manifests as spurious declines.
        both=True (near-diagonal/plateau strata, the swap-maps-s<->1-s
        rule) adds BOTH orientations regardless."""
        is_r1 = tag in ("r1-tiny", "joint-tiny", "r1-tiny-seam")
        try:
            if is_r1:
                P, Q = fast_series_pq(a, b, y, dps)
            else:
                P, Q = guarded_fast_forward(a, b, y, dps=min(dps, 60))
        except Exception:
            return
        try:
            if not (mp.isfinite(P) and mp.isfinite(Q)):
                return
        except TypeError:
            return
        sp, sq = float(P), float(Q)
        if both:
            if math.isfinite(sp) and 0.0 < sp < 1.0:
                self.add(a, b, sp, "p", tag)
            if math.isfinite(sq) and 0.0 < sq < 1.0:
                self.add(a, b, sq, "q", tag)
            return
        if sp <= sq:
            if math.isfinite(sp) and 0.0 < sp < 1.0:
                self.add(a, b, sp, "p", tag)
        else:
            if math.isfinite(sq) and 0.0 < sq < 1.0:
                self.add(a, b, sq, "q", tag)


def gen_r1_tiny(ps, rng, n=800):
    """R1-tiny: small y both orientations, moderate (a,b), fast-series
    territory (constructions #1's own certification traffic)."""
    n0 = len(ps.pts)
    for _ in range(n):
        a = 10.0 ** rng.uniform(-2, 2)
        b = 10.0 ** rng.uniform(-2, 2)
        y = 10.0 ** rng.uniform(-300, -3)
        side = rng.choice(("p", "q"))
        if side == "p":
            ps.add_from_y(a, b, y, "p", "r1-tiny")
        else:
            ps.add_from_y(a, b, 1.0 - y, "q", "r1-tiny")
    status(f"  r1-tiny: {len(ps.pts) - n0} points")


def gen_r1_tiny_qside(ps, rng, n=3000):
    """gen_r1_tiny's own p/q imbalance (measured 17:1, 1872:110 in the
    shipped rows) is STRUCTURALLY EXPECTED per the design's own y_med
    argument (tiny a puts the P<=1/2 crossover at y_med ~ exp(-ln2/a),
    so log-uniform y sampling over a WIDE decade range lands almost
    entirely on the p-side of that crossover) -- but is ALSO compounded
    by a genuine, independently-verifiable defect: gen_r1_tiny's
    q-branch forms the near-1 y argument as NATIVE-FLOAT '1.0 - y' with
    y log-uniform over 10**[-300,-3] -- for y below ~2^-53 (measured:
    1930/2000 = 96.5% of that draw range) this collapses EXACTLY to 1.0
    in double arithmetic (1.0 - 1e-150 == 1.0, hex witness
    0x1.0000000000000p+0), so s_from_y then computes Q(a,b,1.0)=0.0
    EXACTLY and add()'s own 0.0<sigma<1.0 check silently drops the row
    -- nearly every q-side draw was rejected before certification ever
    saw it, independent of the y_med argument. NOT fixed by rewriting
    gen_r1_tiny itself (out of scope here) -- fixed via a direct
    construction instead: build q-side rows DIRECTLY, sigma picked
    log-uniform ON THE Q SIDE (a target probability, never a
    native-float near-1 y), seeded via bid.seed_for and certified
    through the SAME inversion-first machinery (certify_row's
    fast_series/oracle_y path, tag='r1-tiny') every other direct-sigma
    stratum in this file already uses (the huge-nu beyond-resolution
    branch is the precedent). (a,b) drawn from the SAME distribution as
    gen_r1_tiny's own p-side for a genuinely comparable population."""
    n0 = len(ps.pts)
    for _ in range(n):
        a = 10.0 ** rng.uniform(-2, 2)
        b = 10.0 ** rng.uniform(-2, 2)
        sigma = 10.0 ** rng.uniform(-300, -3)
        ps.add(a, b, sigma, "q", "r1-tiny")
    status(f"  r1-tiny q-side (direct construction, RULING 3): "
           f"{len(ps.pts) - n0} points")


def gen_joint_tiny(ps, rng, n=400):
    """Joint-tiny plateau band: both a,b tiny, y interior near s* --
    SEPARATE BUCKET (kappa split happens at certification time)."""
    n0 = len(ps.pts)
    for _ in range(n):
        a = 10.0 ** rng.uniform(-16, -1)
        b = 10.0 ** rng.uniform(-16, -1)
        sstar = b / (a + b)
        # interior y NEAR s* (the actual high-kappa territory) plus a
        # spread further out (still interior, lower kappa -- both
        # orientations needed here per the swap-maps-s<->1-s rule).
        width = rng.choice((1e-6, 1e-3, 1e-1))
        y = min(max(sstar + rng.uniform(-width, width), 1e-6), 1 - 1e-6)
        # plateau/near-diagonal rule (PLAN): BOTH orientations needed.
        ps.add_from_y_smallside(a, b, y, "joint-tiny", both=True)
    status(f"  joint-tiny plateau: {len(ps.pts) - n0} points")


def gen_ridge(ps, rng, n=600):
    """Ridge balanced + skewed sub-bands + the S1/S3 skew seam -- MODERATE
    magnitude only (both-huge-balanced is its own, separate, expensive
    stratum below)."""
    n0 = len(ps.pts)
    for _ in range(n):
        nu = 10.0 ** rng.uniform(1, 8)
        skew = 10.0 ** rng.uniform(-2, 2)
        a = nu * (1 + skew)
        b = a / skew
        c = a + b
        delta = rng.uniform(-6, 6) / math.sqrt(nu)
        y = max(min(a / c + delta, 1 - 1e-15), 1e-15)
        # near-diagonal (balanced, skew close to 1) needs BOTH
        # orientations per the swap-maps-s<->1-s rule; skewed sub-bands
        # get the single small-side orientation (swap identity halves
        # coverage there).
        near_diag = 0.5 <= skew <= 2.0
        ps.add_from_y_smallside(a, b, y, "ridge", both=near_diag)
    status(f"  ridge: {len(ps.pts) - n0} points")


def gen_gammalim_seam(ps, rng, n=400):
    """Gamma-limit dense at the alpha~kGammaAT=20 seam: one param huge
    (gammalim), the other dense around 20 -- moderate cost (single-huge,
    small_side_direct's own guarded gamma_corner_value path). y
    constructed from the ACTUAL gamma-corner transition mapping (t ~
    shape-param scale, gen_beta_data.py's own gamma_corner_value
    identity inverted): a handful of O(1) constants nudged by a bare
    small/huge ratio would land in the saturated regime for all but
    1/400 draws once huge_param reaches 1e100+ (t=shape*O(1) requires
    y within ~shape/huge of the 0/1 boundary, astronomically narrower
    than an O(1)-scale guess for huge_param this large)."""
    n0 = len(ps.pts)
    for _ in range(n):
        small_param = 20.0 * (2.0 ** rng.uniform(-4, 4))
        huge_param = 10.0 ** rng.uniform(2, 250)
        b_is_huge = rng.random() < 0.5
        a, b = (small_param, huge_param) if b_is_huge else (huge_param, small_param)
        t = small_param * (10.0 ** rng.uniform(-1, 1))  # near the gamma transition
        if b_is_huge:
            # t = -b*log1p(-y)  =>  y = -expm1(-t/b)
            y = -math.expm1(-t / huge_param)
        else:
            # t = -a*log(y)     =>  y = exp(-t/a)
            y = math.exp(-t / huge_param)
        y = max(min(y, 1 - 1e-300), 1e-300)
        ps.add_from_y_smallside(a, b, y, "gamma-limit-seam", dps=80)
    status(f"  gamma-limit seam: {len(ps.pts) - n0} points")


def gen_underflow(ps, rng, n=300):
    """Underflow thresholds both ends across the widened
    a*(b)~1/(1074-log2(b)) boundary -- bit-stepped brackets at the
    subnormal/zero y boundary, several (a,b)."""
    n0 = len(ps.pts)
    boundary_ys = (math.ldexp(1.0, -1022), MIN_SUBNORMAL)
    ab_list = [(1e-3, 1.0), (0.1, 0.5), (1.0, 1e-3)]
    for a, b in ab_list:
        for yb in boundary_ys:
            if yb <= 0:
                continue
            s0 = s_from_y(a, b, yb, "p", dps=80)
            if s0 is None:
                continue
            v = s0
            for _ in range(4):
                ps.add(a, b, v, "p", "underflow-threshold")
                v = NEXT_UP(v)
            v = s0
            for _ in range(4):
                v = NEXT_DN(v)
                if v <= 0:
                    break
                ps.add(a, b, v, "p", "underflow-threshold")
    for _ in range(n):
        a = 10.0 ** rng.uniform(-4, 1)
        b = 10.0 ** rng.uniform(-2, 3)
        y = 10.0 ** rng.uniform(-320, -300)
        ps.add_from_y(a, b, y, "p", "underflow-threshold")
    status(f"  underflow thresholds: {len(ps.pts) - n0} points")


def gen_underflow_qside(ps, rng, n=300):
    """q-side direct construction, underflow-threshold
    (single-sided by original construction -- gen_underflow only ever
    built p rows, no q-branch existed at all). Uses the swap identity
    Q(a,b,y) = P(b,a,1-y) AT THE SIGMA LEVEL -- sigma_q = P(b,a, yb)
    computed directly at (b,a) roles and shipped as row
    (a, b, sigma_q, 'q', ...); certify_row solves Q(a,b,y)=sigma_q for
    y independently, and by the swap identity the true y lands near 1
    at the same boundary yb approaches near 0 for P(b,a,.). This NEVER
    forms '1 - y_tiny' in native float (the same collapse-to-1.0
    hazard r1-tiny's q-branch hit, see gen_r1_tiny_qside's docstring),
    since no near-1 y value is ever constructed -- only sigma is."""
    n0 = len(ps.pts)
    boundary_ys = (math.ldexp(1.0, -1022), MIN_SUBNORMAL)
    ab_list = [(1e-3, 1.0), (0.1, 0.5), (1.0, 1e-3)]
    for a, b in ab_list:
        for yb in boundary_ys:
            if yb <= 0:
                continue
            s0 = s_from_y(b, a, yb, "p", dps=80)
            if s0 is None:
                continue
            v = s0
            for _ in range(4):
                ps.add(a, b, v, "q", "underflow-threshold")
                v = NEXT_UP(v)
            v = s0
            for _ in range(4):
                v = NEXT_DN(v)
                if v <= 0:
                    break
                ps.add(a, b, v, "q", "underflow-threshold")
    for _ in range(n):
        a = 10.0 ** rng.uniform(-4, 1)
        b = 10.0 ** rng.uniform(-2, 3)
        y = 10.0 ** rng.uniform(-320, -300)
        s0 = s_from_y(b, a, y, "p", dps=80)
        if s0 is not None:
            ps.add(a, b, s0, "q", "underflow-threshold")
    status(f"  underflow thresholds q-side (direct, RULING 3): "
           f"{len(ps.pts) - n0} points")


def gen_subnormal_y(ps, rng, n=200):
    """Subnormal-y both ends."""
    n0 = len(ps.pts)
    for _ in range(n):
        a = 10.0 ** rng.uniform(-4, 2)
        b = 10.0 ** rng.uniform(-4, 2)
        y = math.ldexp(rng.uniform(1, 2), rng.randint(-1074, -1023))
        side = rng.choice(("p", "q"))
        if side == "p":
            ps.add_from_y(a, b, y, "p", "subnormal-y")
        else:
            ps.add_from_y(a, b, 1.0 - y, "q", "subnormal-y")
    status(f"  subnormal-y: {len(ps.pts) - n0} points")


def gen_subnormal_y_qside(ps, rng, n=400):
    """q-side direct construction, subnormal-y. gen_subnormal_y's
    own q-branch forms '1.0 - y_subnormal' in native float, which is
    EXACTLY 1.0 for EVERY subnormal y (100% collapse -- worse than
    r1-tiny's partial collapse, since y here is always << 2^-1022, hex
    witness 1.0 - 5e-324 == 1.0 -> 0x1.0000000000000p+0), so that
    branch constructed ZERO valid q rows despite a 50/50 coin flip in
    the source. Fixed the same way as underflow-threshold's q-side:
    sigma_q = P(b,a, y_sub) computed directly (swap identity at the
    SIGMA level, no near-1 y ever formed), shipped as
    (a, b, sigma_q, 'q', 'subnormal-y')."""
    n0 = len(ps.pts)
    for _ in range(n):
        a = 10.0 ** rng.uniform(-4, 2)
        b = 10.0 ** rng.uniform(-4, 2)
        y = math.ldexp(rng.uniform(1, 2), rng.randint(-1074, -1023))
        s0 = s_from_y(b, a, y, "p", dps=80)
        if s0 is not None:
            ps.add(a, b, s0, "q", "subnormal-y")
    status(f"  subnormal-y q-side (direct, RULING 3): {len(ps.pts) - n0} points")


# ============================================================================
# huge_nu_beyond_resolution's derivation.
#
# TRAP AVOIDED: testing "resolvable" by evaluating the FORWARD
# probability at a SINGLE z-probe point (z=1 or z=3) and checking
# whether it is exactly 0.0/1.0 in double tests saturation of P/Q in
# SIGMA-space, not collapse of the SOLVED y in Y-space -- the two are
# unrelated at the balanced point: for skew=1 (a=b), any z=O(1) probe
# sits near the MEAN, where P=Q=0.5 by symmetry -- erfc(z)/2 for z=O(1)
# is never within orders of magnitude of 0 or 1 at ANY nu, so a
# single-probe test reports "resolvable" regardless of nu (the balanced
# case degenerates entirely: there is no nu at which a bounded-z probe
# ever saturates, because the probe point never leaves the O(1)
# neighbourhood of the mean).
#
# CORRECT CRITERION: what the kernel's contract actually needs is
# whether the ENTIRE achievable P/Q transition -- from the leftmost
# resolvable z to the rightmost -- maps into <=1 ULP of y. Z_MAX is the
# z (in the SAME convention _huge_nu_y_mpf already uses: z=sqrt(cpsi),
# the direct argument to erfc, i.e. the standard-normal-quantile
# convention divided by sqrt(2)) at which the leading erfc(z)/2 term
# itself underflows the smallest positive double (2^-1074): beyond it no
# sigma in (0,1) can ever select a MORE extreme point, so z=+-Z_MAX are
# the outermost distinguishable probe points. Derived by mpf bisection
# (dps=60, 200 iterations): erfc(Z_MAX)/2 = 2^-1074 at
# Z_MAX = 27.2005633665362563777... (the same point in the OTHER,
# sqrt(2)-scaled Phi-quantile convention: 27.2006*sqrt(2) = 38.465).
#
# BEYOND-RESOLUTION(a,b) iff y(z=+Z_MAX) and y(z=-Z_MAX), each
# independently constructed at mpf precision then double-rounded ONCE
# (matching _huge_nu_y_mpf's own one-rounding protocol), differ by <=1
# ULP of y -- the entire achievable transition collapses inside one ulp,
# so no sigma in (0,1) can ever select more than one double for this
# (a,b) pair.
#
# SANITY ANCHORS (order-of-magnitude cross-checks): balanced
# CENTRAL-band collapse (sigma=0.3 rounding onto the exact double mean)
# measured between nu=1e31 (False) and nu=1e32 (True). Balanced FULL
# collapse (this criterion) measured at nu* ~ 6.0e34 (skew=1); skewed
# sub-bands (skew 1e2/1e4/1e6) at nu* ~ 1.1e35/1.3e35/4.7e34 -- all in
# the same mid-10^34-to-10^35 decade.
#
# USE: per-(a,b), NOT a nu* lookup table (nu* depends on skew, as
# measured above) -- huge_nu_beyond_resolution() below is called
# directly wherever this generator needs to classify a candidate (a,b)
# pair for the huge-nu strata. Existing shipped B rows (nu>=1e35) are
# unchanged by this criterion -- it governs NEW row construction only.
# ============================================================================
_ZMAX_STR = ("27.2005633665362563777429614681942258114152388192067686497169"
             "338076019593011258022472543490099964173")


def _ulp_distance_pos(x, y):
    """ULP distance between two FINITE, POSITIVE (or +0.0) doubles -- exact
    bit-pattern difference. y in (0,1) always here, so no sign handling is
    needed (unlike a general ulp-distance utility)."""
    xi = struct.unpack('<Q', struct.pack('<d', x))[0]
    yi = struct.unpack('<Q', struct.pack('<d', y))[0]
    return abs(xi - yi)


def huge_nu_beyond_resolution(a, b, dps=80):
    """The correct per-(a,b) criterion: True iff the entire achievable
    y-transition (probed at the two outermost distinguishable z,
    +-Z_MAX) collapses inside <=1 ULP of y."""
    zmax = mp.mpf(_ZMAX_STR)
    y_pos = _huge_nu_y_mpf(a, b, zmax, dps=dps)
    y_neg = _huge_nu_y_mpf(a, b, -zmax, dps=dps)
    return _ulp_distance_pos(y_pos, y_neg) <= 1


def _huge_nu_y_mpf(a, b, target_z, dps=80):
    """Construct y at mpf precision so that cpsi = target_z^2 EXACTLY
    (to mpf precision) at the point BEFORE double-rounding, THEN round
    to double. TWO hazards to avoid:
    (1) building y = a/c + delta in NATIVE PYTHON FLOAT arithmetic: its
        ~1e-16 relative rounding error in the 'mean' translates (via
        lam = a - c*y) to an ABSOLUTE error in lam of order c*1e-16 --
        at a,b ~ 1e40+ that alone drives cpsi up to ~1e19-scale even at
        the INTENDED exact-mean point.
    (2) gb.r3_setup's own 'zeta' parameter is NOT target_z directly --
        r3_setup's docstring/_lambda_of_zeta target is
        cpsi = zeta^2 * nu (the RIDGE-normalized deviation, zeta =
        z/sqrt(nu)), not cpsi = zeta^2. Passing target_z straight through
        as zeta silently requests cpsi = target_z^2 * nu ~ nu itself
        (e.g. target_z=1.0 at nu=1e15 produces cpsi=999999999999999.9,
        not 1) -- deep in saturation regardless of target_z. FIX:
        zeta = target_z / sqrt(nu), so cpsi = zeta^2*nu = target_z^2
        exactly, as intended."""
    with mp.workdps(dps):
        a_m, b_m = mp.mpf(a), mp.mpf(b)
        c = a_m + b_m
        nu = a_m * b_m / c
        p = a_m / c
        zeta = mp.mpf(target_z) / mp.sqrt(nu)
        alpha, beta, cc, xi, lam = gb.r3_setup(nu, p, zeta, dps)
        # r3_setup's (alpha,beta,c) are the CANONICAL (nu,p)-derived
        # values, which equal (a,b,c) exactly by construction (nu,p come
        # from a,b) -- xi is the y this generator wants, in the ORIGINAL
        # (a,b) frame already (r3_setup solves in the same frame its
        # nu,p,c came from).
        y = float(xi)
        return max(min(y, 1.0 - 1e-300), 1e-300)


def _huge_nu_mean_ulp(a, b, k, dps=80):
    """The double-rounded mean a/(a+b), stepped k ULPs (k=0 is the mean
    itself). mpf-precision mean, THEN one double rounding, THEN native
    nextafter stepping -- the ONLY robust way to probe this close to the
    resolution boundary (native-float mean arithmetic reintroduces the
    same c*eps error _huge_nu_y_mpf's own docstring documents)."""
    with mp.workdps(dps):
        a_m, b_m = mp.mpf(a), mp.mpf(b)
        y = float(a_m / (a_m + b_m))
    for _ in range(abs(k)):
        y = math.nextafter(y, math.inf if k > 0 else -math.inf)
    return max(min(y, 1.0 - 1e-300), 1e-300)


def gen_huge_nu(ps, rng, n=70, n_beyond=90):
    """Huge-nu: TWO regimes, CLASSIFIED BY the per-(a,b) criterion
    (huge_nu_beyond_resolution) rather than a nu* lookup table (nu*
    depends on skew, see that function's own derivation).
    (a) resolvable GUARD territory (construction #2's own route -- z=O(1)
    points land on distinguishable doubles) via _huge_nu_y_mpf, a
    mpf-precision z-targeted construction (validated: r1-vs-r2 agreement
    to ~30+ decimal digits even at skew=1e6/nu=1e20). Candidates are
    DRAWN from a range comfortably below the measured threshold band
    (~6e34-1.3e35 across skew 1..1e6) and VERIFIED per-pair via
    huge_nu_beyond_resolution -- any draw that turns out to already be
    beyond-resolution is redrawn (rejection sampling), so this bucket is
    a genuine ordinary-bracket-certifiable population, not an
    approximation.
    (b) genuine BEYOND-RESOLUTION territory: nu high enough the WHOLE
    transition collapses within <=1 ULP of the mean -- constructed via
    _huge_nu_mean_ulp (mpf-precision mean, THEN double-rounded, THEN
    ULP-stepped -- the z-targeted construction is NOT usable here since
    at this nu ANY z=O(1..10) point rounds back to one of the SAME
    3-5 doubles nearest the mean), with sigma picked to deliberately
    NOT be exactly the trivial 0.5 (mid-ulp) so the certifier's
    nearest-neighbor contract is genuinely exercised rather than the
    trivial straddle case gammainv's own 'dilution lesson' warns about.
    Candidates are drawn from a range comfortably ABOVE the measured
    threshold band and VERIFIED per-pair via huge_nu_beyond_resolution
    (redrawn if not yet collapsed) -- classifying candidates per-(a,b)
    this way, rather than trusting a fixed nu range, is what makes both
    buckets genuinely reliable populations."""
    n0 = len(ps.pts)
    zmax = mp.mpf(_ZMAX_STR)

    # (a) resolvable guard territory -- rejection-sampled against the
    # huge_nu_beyond_resolution criterion.
    n_ok = 0
    tries = 0
    while n_ok < n and tries < n * 6:
        tries += 1
        skewed = rng.random() < 0.5
        if skewed:
            nu = 10.0 ** rng.uniform(math.log10(B_GL) - 1, 34.0)
            skew = 10.0 ** rng.uniform(2, math.log10(SKEW_SAFE_CAP))
        else:
            nu = 10.0 ** rng.uniform(math.log10(B_GL) - 1, 34.0)
            skew = 10.0 ** rng.uniform(0, 1)
        a, b = (nu, nu * skew) if skewed else (nu * (1 + skew), nu * (1 + skew) / skew)
        if huge_nu_beyond_resolution(a, b):
            continue  # already collapsed -- belongs in bucket (b), redraw
        delta_z = rng.uniform(-4, 4)
        y = _huge_nu_y_mpf(a, b, delta_z)
        before = len(ps.pts)
        ps.add_from_y_smallside(a, b, y, "huge-nu", dps=80)
        if len(ps.pts) > before:
            n_ok += 1

    # (b) genuine beyond-resolution: sigma picked DIRECTLY (the
    # gammainv reference generator's own huge-a-beyond-resolution-target
    # pattern) -- NOT forward-constructed from y here: forward(mean +/-
    # k ulps) is ITSELF already saturated to exactly {0,0.5,1} at this
    # collapse depth, so s_from_y's well-posedness filter would reject
    # essentially every attempt -- there is no 'intermediate' reachable
    # sigma to construct from y at true collapse depth; the well-posed
    # INPUT contract here is simply sigma in (0,1), the same as any
    # kernel call, and certify_row's own nearest-of-{neighbors} contract
    # is what proves the answer, not a forward round-trip.
    beyond_sigmas = (0.5, NEXT_UP(0.5), NEXT_DN(0.5), 0.3, 0.7, 0.1, 0.9,
                     1e-3, 1.0 - 1e-3)
    n_beyond_ok = 0
    tries = 0
    while n_beyond_ok < n_beyond and tries < n_beyond * 6:
        tries += 1
        skewed = rng.random() < 0.5
        if skewed:
            nu = 10.0 ** rng.uniform(35.0, 42.0)
            skew = 10.0 ** rng.uniform(2, math.log10(SKEW_SAFE_CAP))
            a, b = nu, nu * skew
        else:
            nu = 10.0 ** rng.uniform(35.0, 42.0)
            a = b = nu * 2
        if not huge_nu_beyond_resolution(a, b):
            continue  # not actually collapsed yet -- belongs in (a), redraw
        sigma = rng.choice(beyond_sigmas)
        side = rng.choice(("p", "q"))
        before = len(ps.pts)
        ps.add(a, b, sigma, side, "huge-nu")
        if len(ps.pts) > before:
            n_beyond_ok += 1
    status(f"  huge-nu (both-huge-balanced guard + beyond-resolution): "
           f"{len(ps.pts) - n0} points ({n_ok} guard / {n_beyond_ok} beyond, "
           f"criterion-verified)")


def gen_seam_bracket(ps, rng, n=150):
    """a_T-seam bit-stepped bracket (S1/S3 seed seam near kGammaAT=20)."""
    n0 = len(ps.pts)
    a_seam = 20.0
    a_vals = [a_seam]
    v = a_seam
    for _ in range(3):
        v = NEXT_UP(v)
        a_vals.append(v)
    v = a_seam
    for _ in range(3):
        v = NEXT_DN(v)
        a_vals.append(v)
    for a in a_vals:
        for b in (1.0, 20.0, 100.0):
            for y in (0.1, 0.5, 0.9):
                for side in ("p", "q"):
                    ps.add_from_y(a, b, y, side, "at-seam-bracket")
    for _ in range(n):
        a = a_seam * (1 + rng.uniform(-1e-6, 1e-6))
        b = 10.0 ** rng.uniform(-1, 2)
        y = rng.uniform(0.05, 0.95)
        ps.add_from_y_smallside(a, b, y, "at-seam-bracket")
    status(f"  a_T-seam bracket: {len(ps.pts) - n0} points")


def find_y_cut(a, b, side):
    """Bisect (in log-y) for the y where bid.deep_small_cut_bound(a,b,y,
    side) crosses DEEP_SMALL_CUT. A blind log-uniform sigma grid (10^pe,
    pe in -320..-1) has NO regard for whether that sigma's corresponding
    y is anywhere near the actual cut for the (a,b) pair drawn -- for a
    near a_T=19 the cut sits at y~2^-60/19~4.6e-20, so almost the ENTIRE
    pe range constructs a sigma whose true y is nowhere near deep-small
    territory, correctly falling through to ordinary root-finding, which
    then legitimately fails for many of those ill-matched (a,b,sigma)
    combinations (measured 415/540, 77%, drop rate) -- a
    GENERATOR-QUALITY gap, not a certifier defect. FIX (inversion-first,
    matching gen_gammalim_seam's own fix): locate y_cut precisely per
    (a,b,side) via bisection on the ACTUAL routing predicate (not a
    re-derived approximation of it), then sample y at chosen MULTIPLES
    of y_cut (both inside and outside the deep-small regime -- boundary
    coverage is itself valuable, per the design's own
    'deep-small-cut-bracket' precedent in gammainv), and construct
    sigma = forward(y) (in-band, reachable BY CONSTRUCTION). Returns
    z_cut, the SMALL variable's cut location (z=y for side='p', z=1-y
    for side='q' -- bid.deep_small_cut_bound's own yp convention) --
    caller converts to y."""
    def bound_at_z(z):
        y = z if side == "p" else (1.0 - z)
        return bid.deep_small_cut_bound(a, b, y, side)

    lo, hi = 1e-300, 0.5
    if bound_at_z(hi) < DEEP_SMALL_CUT:
        return hi
    if bound_at_z(lo) > DEEP_SMALL_CUT:
        return lo
    for _ in range(80):
        mid = math.sqrt(lo * hi)
        v = bound_at_z(mid)
        if v < DEEP_SMALL_CUT:
            lo = mid
        else:
            hi = mid
    return lo


def gen_deep_small_both(ps, rng, n=300):
    """Deep-small closed-form territory, BOTH orientations from the
    start -- y constructed from the (a,b,side)-SPECIFIC cut location
    (find_y_cut), sigma=forward(y) (in-band by construction)."""
    n0 = len(ps.pts)
    a_list = (1e-300, 1e-30, 1e-4, 0.1, 1.0, 19.0)
    b_list = (0.5, 3.0, 1e5, 1e300)
    factors = (1e-10, 1e-4, 1e-1, 1.0, 3.0, 30.0, 1e3)
    for a in a_list:
        for b in b_list:
            for side in ("p", "q"):
                try:
                    zc = find_y_cut(a, b, side)
                except (OverflowError, ValueError, ZeroDivisionError):
                    continue
                if not (0 < zc < 1):
                    continue
                for factor in factors:
                    zv = min(max(zc * factor, 1e-320), 0.5)
                    yv = zv if side == "p" else (1.0 - zv)
                    ps.add_from_y(a, b, yv, side, "deep-small", dps=80)
    for _ in range(n):
        a = 10.0 ** rng.uniform(-300, 1.3)
        b = 10.0 ** rng.uniform(-2, 300)
        side = rng.choice(("p", "q"))
        try:
            zc = find_y_cut(a, b, side)
        except (OverflowError, ValueError, ZeroDivisionError):
            continue
        if not (0 < zc < 1):
            continue
        zv = min(max(zc * (10.0 ** rng.uniform(-8, 3)), 1e-320), 0.5)
        yv = zv if side == "p" else (1.0 - zv)
        ps.add_from_y(a, b, yv, side, "deep-small", dps=80)
    status(f"  deep-small (both orientations): {len(ps.pts) - n0} points")


def build_point_set(rng):
    ps = PointSet()
    # Per-stratum n values below are sized against measured per-row
    # cost, not just a target row count: a MINORITY of joint-tiny rows
    # (sigma within a few ULPs of 0 or 1) drive bisection into many
    # extra iterations of bid.r1_value_mp's own up-to-4000-term series
    # (measured outlier cost 9.1s vs a ~150-300ms typical row), so
    # bumping joint-tiny's n has an outsized effect on total runtime.
    # The direct q-side strata (gen_r1_tiny_qside, gen_underflow_qside,
    # gen_subnormal_y_qside) sidestep the native-float 1-y collapse bug
    # entirely, so their yield is close to 100% (see each gen_*_qside
    # docstring) and need less oversampling than the y-constructed
    # strata to hit their target count.
    # Design targets (both sides combined): r1-tiny 4-6k, ridge 3-4k,
    # gammalim 2-3k, joint-tiny 1.5-2.5k SEPARATE, underflow 1-1.5k,
    # subnormal-y 0.8-1.2k, huge-nu-B 1-1.5k SEPARATE, a_T-seam
    # 0.5-0.8k.
    gen_r1_tiny(ps, rng, n=6000)
    gen_r1_tiny_qside(ps, rng, n=3000)          # direct q-side construction (near-100% yield)
    gen_joint_tiny(ps, rng, n=1000)             # both=True -> ~2k, within the 1.5-2.5k design band
    gen_ridge(ps, rng, n=5500)                  # yield ~58%; sized toward the 3-4k target
    gen_gammalim_seam(ps, rng, n=9000)          # construction yield ~28%; sized toward the 2-3k target
    gen_underflow(ps, rng, n=600)
    gen_underflow_qside(ps, rng, n=600)         # direct q-side construction (near-100% yield)
    gen_subnormal_y(ps, rng, n=1300)            # bumped for p/q balance against the qside stratum's much higher yield
    gen_subnormal_y_qside(ps, rng, n=800)       # direct q-side construction (near-100% yield)
    gen_huge_nu(ps, rng, n=300, n_beyond=1400)  # n_beyond sized up: B-marked-row
                                                 # conversion is well under 50%;
                                                 # criterion-verified, sized
                                                 # toward the 1-1.5k SEPARATE target
    gen_seam_bracket(ps, rng, n=600)            # within the 0.5-0.8k target
    gen_deep_small_both(ps, rng, n=400)
    status(f"total distinct (a,b,sigma,side) points: {len(ps.pts)}")
    for tag, n in sorted(ps.strata_counts.items()):
        status(f"  stratum {tag}: {n}")
    return ps


# ============================================================================
# Part 11: checkpointed, resumable compute pass (beta-reference /
# gammainv-reference precedent: re-run until it writes the files).
# ============================================================================
CKPT_PATH = os.path.join(tempfile.gettempdir(), f"corvus_betainv_ref_ckpt_{SEED}.tsv")
WALL_CLOCK_BUDGET_S = 470.0


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


def as_bits(x):
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def compute_all(ps):
    total = len(ps.pts)
    # Point-bits digest signature (NUMERICAL-DOCTRINE.md's binding rule: a
    # checkpoint signature must bind to the POINT IDENTITIES, not just the
    # count -- an edit that preserves N would otherwise replay stale
    # oracle values under new point identities; gen_beta_reference.py's
    # r4huge_append is the pattern). Digests the FULL dedup identity
    # PointSet.add keys on -- (a, b, sigma, side), side as one byte; tag
    # is display-only and deliberately excluded.
    dig = hashlib.sha256()
    for a, b, sigma, side, _tag in ps.pts:
        dig.update(struct.pack("<QQQB", as_bits(a), as_bits(b),
                               as_bits(sigma), 0 if side == "p" else 1))
    sig = f"v1-{dig.hexdigest()[:16]} SEED={SEED} N={total}"
    done_map = load_checkpoint(CKPT_PATH, sig)
    status(f"checkpoint: {len(done_map)}/{total} points already computed ({CKPT_PATH})")

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
        for idx, (a, b, sigma, side, tag) in enumerate(ps.pts):
            if idx in done_map:
                continue
            if time.time() - t_start > WALL_CLOCK_BUDGET_S:
                status(f"wall-clock budget ({WALL_CLOCK_BUDGET_S:.0f}s) hit at "
                       f"{idx}/{total} ({newly_done} computed this run) -- "
                       f"re-run to continue.")
                return None, False
            try:
                r = certify_row(a, b, sigma, side, tag)
            except Exception as e:
                r = {"yd": None, "certified": False, "method": f"exception:{e}",
                     "marker": "N", "kappa": None}
            if r["yd"] is None or not r["certified"]:
                append_checkpoint(fh, idx, ["FAILED", r.get("method", "?")])
            else:
                append_checkpoint(fh, idx, [hexd(r["yd"]), r["marker"], r["method"]])
            newly_done += 1
            if newly_done % 25 == 0:
                status(f"  ... {idx + 1}/{total} ({time.time() - t_start:.0f}s this run)")

    status(f"computed {newly_done} points this run ({time.time() - t_start:.0f}s); "
           f"all {total} points now checkpointed.")

    done_map = load_checkpoint(CKPT_PATH, sig)
    rows_p, rows_q = [], []
    n_failed = 0
    fail_by_tag = {}
    marker_counts = {}
    for idx, (a, b, sigma, side, tag) in enumerate(ps.pts):
        fields = done_map[idx]
        if fields[0] == "FAILED":
            n_failed += 1
            fail_by_tag[tag] = fail_by_tag.get(tag, 0) + 1
            continue
        yd = float.fromhex(fields[0])
        marker = fields[1]
        marker_counts[marker] = marker_counts.get(marker, 0) + 1
        (rows_p if side == "p" else rows_q).append((a, b, sigma, yd, marker))
    status(f"certified: {len(rows_p)} p-rows, {len(rows_q)} q-rows; "
           f"{n_failed} dropped (uncertified)")
    if fail_by_tag:
        status(f"drops by stratum: {fail_by_tag}")
    status(f"marker distribution: {marker_counts}")
    return (rows_p, rows_q), True


# ============================================================================
# Part 11b: huge-parameter corner append. Gates the betainv-side half of
# the Dekker-ceiling audit (src/betainv-inl.h's exact power-of-two
# prescale on the E prefactor's DdMulD(lrxi, ra) / DdMulD(lryv, rb) sites,
# the twin of beta forward's BetaR4Tiny fix): rows like
# beta_p_inv(2, 1e307, 0.5) exercise ra/rb above ops::ProdLow's 2^996
# non-FMA Dekker ceiling with the b*ln(1-y) product O(1) and load-bearing.
#
# min(a,b) here is MODERATE (0.5..100), not huge -- both_huge_balanced's
# own guard (min(a,b) >= GUARD_MIN_THRESHOLD = 1e10) never fires, so the
# gamma-corner hang risk Part 1/2's guard exists for does not apply; the
# standard "full" forward route (gbr.small_side_direct) is safe and fast
# here, exactly as the task brief anticipated.
#
# What DOES bite at this corner: certify_row's DEFAULT root-find path
# (oracle_y -> guarded_fast_forward -> bid.betainv_forward -> gb.route_final)
# SILENTLY MISROUTES at this skew -- measured witness: at (a=2, b=1e307),
# bid.betainv_forward(a, b, y) returns a pure {0.0, 1.0} STEP over
# y in [2e-307, 1e-306] where the true P rises smoothly 0.594 -> 0.997 ->
# 0.9995, driving the root-find to converge on a numerically-plausible but
# WRONG root -- the module docstring's own "SHARED-MACHINERY CAVEAT",
# also independently documented by oracle_y_audited's own docstring for
# the gamma-limit-seam stratum (same disease, different witness). FIX:
# exactly what oracle_y_audited exists for -- bisect directly against
# full_forward (the audited small_side_direct route), then certify the
# result through certify_row's own yd_override path (which also calls
# full_forward, never the buggy cheap router, for the actual bracket
# check).
#
# A SECOND, distinct hazard shows up only with the huge parameter FIRST
# (a huge, the exponent-on-y side): the true quantile sits within O(1/a)
# of y = 1, which for a >= 2^900 is hundreds of decimal digits closer to 1
# than any double distinguishes from 1.0 -- oracle_y_audited's bracket
# search correctly finds NO sign change (both endpoints are already
# saturated at working dps) and returns None. This is not a decline: it
# is the SAME deep-small-boundary collapse certify_row's own yd-in-{0,1}
# branch exists for, just never reached because the audited bisection
# never lands on the boundary value it would recognize. FALLBACK:
# certify_deep_small directly at both dps layers (60, 100) -- the same
# call certify_row itself would make had oracle_y/oracle_y_audited handed
# it yd in {0.0, 1.0}. Rows that clear NEITHER route are DECLINED and
# counted, never guessed (checked empirically: ~34% of the huge-FIRST
# orientation genuinely needs this fallback and a further fraction still
# declines -- the quantile there is neither deep-small-provable nor
# bracket-representable at standard dps, a real DECLINE, not a bug).
# ============================================================================
def _nextafter_pos_inf(v):
    return math.nextafter(v, math.inf)


def gen_betainv_huge_corner():
    """(a, b, sigma, side, tag) attempted points. One parameter log-spaced
    across [2^900, 1.7e308) bracketing the non-FMA Dekker ceiling 2^996;
    the other in {0.5, 1, 2, 5, 20, 100} (task spec, verbatim). sigma in
    {1e-6, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99}; side by the natural
    sigma<=0.5 -> p / sigma>0.5 -> q convention (s>1/2 needs no separate
    complement construction -- unlike beta forward's xi, sigma is used
    directly as either function's own probability argument, so the s>1/2
    half of the list already exercises the q-side calls). BOTH parameter
    orders (huge first and second)."""
    huge_list = [_nextafter_pos_inf(2.0 ** 900), 2.0 ** 996, 1.4e308]
    other_list = [0.5, 1.0, 2.0, 5.0, 20.0, 100.0]
    sigma_list = [1e-6, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99]
    pts = []
    seen = set()

    def add(a, b, sigma, tag):
        side = "p" if sigma <= 0.5 else "q"
        key = (a, b, sigma, side)
        if key in seen:
            return
        seen.add(key)
        pts.append((a, b, sigma, side, tag))

    for hh in huge_list:
        for oo in other_list:
            for sigma in sigma_list:
                add(hh, oo, sigma, "betainv-huge-corner-first")
                add(oo, hh, sigma, "betainv-huge-corner-second")
    add(2.0, 1e307, 0.5, "betainv-huge-corner-witness")  # task-example witness
    return pts


def certify_huge_corner_row(a, b, sigma, side, tag):
    """Three-tier certification (see module comment above): audited
    bisection + standard bracket (covers the huge-SECOND orientation and
    most huge-FIRST rows whose quantile is not deep-small); deep-small-
    direct fallback (the huge-FIRST boundary-collapse corner); DECLINE
    (counted, never guessed) if neither certifies."""
    try:
        seed = bid.seed_S3(a, b, sigma, side)
    except Exception:
        seed = None
    if seed is None or not (math.isfinite(seed) and 0.0 < seed < 1.0):
        try:
            seed = bid.seed_for(a, b, sigma, side)
        except Exception:
            seed = None

    if seed is not None:
        yd_aud = oracle_y_audited(a, b, sigma, side, seed, dps=60)
        if yd_aud is not None:
            r = certify_row(a, b, sigma, side, tag, yd_override=yd_aud)
            if r["certified"]:
                return r

    layer_results = [certify_deep_small(a, b, sigma, side, dps) for dps in (60, 100)]
    if (layer_results[0]["certified"] and layer_results[1]["certified"]
            and layer_results[0]["yd"] == layer_results[1]["yd"]):
        return {"yd": layer_results[1]["yd"], "certified": True,
                "method": "deep-small-fallback", "marker": "N"}

    return {"yd": None, "certified": False, "method": "declined", "marker": "N"}


CKPT_HUGE_CORNER_PATH = os.path.join(
    tempfile.gettempdir(), f"corvus_betainv_ref_ckpt_hugecorner_{SEED}.tsv")
HUGE_CORNER_WALL_BUDGET_S = 260.0


def compute_huge_corner(pts):
    total = len(pts)
    # Point-bits digest signature -- same rule/shape as compute_all above
    # (full (a, b, sigma, side) identity, side as one byte).
    dig = hashlib.sha256()
    for a, b, sigma, side, _tag in pts:
        dig.update(struct.pack("<QQQB", as_bits(a), as_bits(b),
                               as_bits(sigma), 0 if side == "p" else 1))
    sig = f"v1-{dig.hexdigest()[:16]} SEED={SEED} N={total}"
    done_map = load_checkpoint(CKPT_HUGE_CORNER_PATH, sig)
    status(f"huge-corner checkpoint: {len(done_map)}/{total} already computed "
           f"({CKPT_HUGE_CORNER_PATH})")
    existing_sig = None
    if os.path.exists(CKPT_HUGE_CORNER_PATH):
        with open(CKPT_HUGE_CORNER_PATH, "r") as f0:
            existing_sig = f0.readline().strip()
    mode = "a" if existing_sig == sig else "w"
    t_start = time.time()
    newly = 0
    with open(CKPT_HUGE_CORNER_PATH, mode) as fh:
        if mode == "w":
            fh.write(sig + "\n")
            fh.flush()
        for idx, (a, b, sigma, side, tag) in enumerate(pts):
            if idx in done_map:
                continue
            if time.time() - t_start > HUGE_CORNER_WALL_BUDGET_S:
                status(f"huge-corner wall-clock budget hit at {idx}/{total} "
                       f"({newly} this run) -- re-run with --huge-corner-append "
                       f"to continue.")
                return None, False
            try:
                r = certify_huge_corner_row(a, b, sigma, side, tag)
            except Exception as e:
                r = {"yd": None, "certified": False, "method": f"exception:{e}"}
            if r["yd"] is None or not r["certified"]:
                append_checkpoint(fh, idx, ["FAILED", r.get("method", "?")])
            else:
                append_checkpoint(fh, idx, [hexd(r["yd"]), r["marker"], r["method"]])
            newly += 1
    status(f"huge-corner: computed {newly} points this run "
           f"({time.time() - t_start:.0f}s); all {total} points now checkpointed.")

    done_map = load_checkpoint(CKPT_HUGE_CORNER_PATH, sig)
    rows_p, rows_q = [], []
    fail_by_key = {}
    for idx, (a, b, sigma, side, tag) in enumerate(pts):
        fields = done_map[idx]
        if fields[0] == "FAILED":
            key = (tag, fields[1] if len(fields) > 1 else "?")
            fail_by_key[key] = fail_by_key.get(key, 0) + 1
            continue
        yd = float.fromhex(fields[0])
        marker = fields[1]
        (rows_p if side == "p" else rows_q).append((a, b, sigma, yd, marker))
    status(f"huge-corner certified: {len(rows_p)} p-rows, {len(rows_q)} q-rows "
           f"(of {total} attempted)")
    if fail_by_key:
        status(f"huge-corner declines by (tag, method): {fail_by_key}")
    return (rows_p, rows_q), True


def _splice_huge_corner(rows_p, rows_q):
    for path, rows in (
        (os.path.join(REPO, "tests", "data", "betainv_p_reference.txt"), rows_p),
        (os.path.join(REPO, "tests", "data", "betainv_q_reference.txt"), rows_q)):
        with open(path) as f:
            lines = [ln.rstrip("\n") for ln in f if ln.strip()]
        seen = {tuple(ln.split()[:3]) for ln in lines}
        added = []
        for a, b, sigma, yd, marker in rows:
            key = (hexd(a), hexd(b), hexd(sigma))
            if key in seen:
                continue  # PointSet-consistent: an earlier family owns the row
            seen.add(key)
            added.append(f"{key[0]} {key[1]} {key[2]} {hexd(yd)} {marker}")
        with open(path, "w", newline="\n") as f:
            f.write("\n".join(lines + added) + "\n")
        status(f"wrote {path}: +{len(added)} huge-corner rows "
               f"({len(lines) + len(added)} total)")


def huge_corner_append():
    t_all = time.time()
    if not negative_controls():
        print("\nFATAL: negative control(s) were ACCEPTED -- the certifier "
              "is not rejecting known-bad rows. Aborting, nothing written.",
              file=sys.stderr)
        return 2
    status("negative controls: all rejected (correct).")

    pts = gen_betainv_huge_corner()
    status(f"huge-corner point set: {len(pts)} attempted points")

    result, done = compute_huge_corner(pts)
    if not done:
        status("PARTIAL RUN: re-invoke with --huge-corner-append to continue "
               "(checkpoint saved).")
        return 3
    rows_p, rows_q = result
    if not rows_p and not rows_q:
        print("\nFAILED: zero rows certified.", file=sys.stderr)
        return 1
    _splice_huge_corner(rows_p, rows_q)
    status(f"huge-corner-append runtime: {time.time() - t_all:.1f}s")
    return 0


# ============================================================================
# Part 12: write rows. Format (ratified deviation, see module docstring):
# five hex tokens per row: a b sigma yd marker.
# ============================================================================
def write_rows(rows_p, rows_q):
    with open(os.path.join(REPO, "tests", "data", "betainv_p_reference.txt"), "w") as f:
        for a, b, sigma, yd, marker in rows_p:
            f.write(f"{hexd(a)} {hexd(b)} {hexd(sigma)} {hexd(yd)} {marker}\n")
    with open(os.path.join(REPO, "tests", "data", "betainv_q_reference.txt"), "w") as f:
        for a, b, sigma, yd, marker in rows_q:
            f.write(f"{hexd(a)} {hexd(b)} {hexd(sigma)} {hexd(yd)} {marker}\n")
    status(f"wrote betainv_p_reference.txt: {len(rows_p)} rows")
    status(f"wrote betainv_q_reference.txt: {len(rows_q)} rows")


# ============================================================================
# Part 12b: x=0 stratum (#14 gap). NEITHER betainv reference file has a row
# whose correctly-rounded output is exactly 0.0, so the ULP gate's x=0
# cross bucket has been silently empty since it was written. Constructed
# OUTSIDE compute_all/PointSet: no root-finding is needed (the output is
# the fixed constant +0.0), so these rows carry NO checkpoint entry --
# deterministic, code-defined points carry no staleness hazard (nothing
# for the point-bits digest in N14.1 to protect against).
#
# CERTIFICATION (mandatory per row): the correctly rounded root is +0.0
# iff the true (infinite-precision) root is strictly < 2^-1075 -- exactly
# half the smallest positive subnormal double; the exact tie 2^-1075
# itself also rounds to 0.0 under ties-to-even (0 is the even neighbor),
# so a strict "<" test is the whole predicate.
#   p-side: P(a,b,y) is increasing in y, so P(a,b,2^-1075) > s implies the
#     true root y* (solving P(a,b,y*)=s) satisfies y* < 2^-1075.
#   q-side: Q(a,b,y) = 1-P(a,b,y) is DEcreasing in y, so Q(a,b,2^-1075) < s,
#     i.e. P(a,b,2^-1075) > 1-s, implies the true root y* satisfies
#     y* < 2^-1075 by the same increasing-P argument.
# Evaluated via bid.r1_value_mp -- the file's own validated tiny-x
# machinery (the fast_series/R1 route fast_vs_full_validate exercises;
# this call sits exactly in its R1-tiny regime), at dps=80 (>=60 floor).
# ============================================================================
_X0_Y_TINY = mp.mpf(2) ** -1075
_X0_DPS = 80
_X0_P_SIDE_AB = [(a, b) for a in (0.05, 0.25, 0.5, 0.9) for b in (0.5, 2.0, 37.5)]
_X0_Q_SIDE_AB = [(a, b) for a in (0.01, 0.02, 0.04) for b in (0.5, 2.0, 8.0)]


def _x0_certified(a, b, s, side, dps=_X0_DPS):
    """True iff the true root for (a,b,s,side) is provably < 2^-1075 (so
    the correctly-rounded double is +0.0) -- see the module comment above
    for the two per-side inequalities."""
    p_tiny = bid.r1_value_mp(a, b, _X0_Y_TINY, dps)
    s_m = mp.mpf(s)
    if side == "p":
        return bool(p_tiny > s_m)
    return bool((1 - p_tiny) < s_m)


def gen_x_zero_rows():
    """Returns (rows_p, rows_q, err). err is None on success; otherwise a
    string naming the failure ('predicate-broken' from the negative
    control, or 'candidate-failed' from a per-row certification miss) and
    rows_p/rows_q are None -- caller (main()) turns this into the exit
    code (mirrors negative_controls()'s own contract: nothing is written
    on failure)."""
    # NEGATIVE CONTROL FIRST (house exit-2 pattern, see negative_controls()
    # above): (a=0.5, b=2.0, s=0.5) has an ORDINARY root (P(a,b,2^-1075) is
    # astronomically smaller than 0.5) -- the predicate MUST reject it.
    if _x0_certified(0.5, 2.0, 0.5, "p"):
        print("  [x0-control] a=5.000000e-01 b=2.000000e+00 s=5.000000e-01 "
              "(p) -> ACCEPTED (FATAL BUG)", file=sys.stderr)
        status("x0 negative control: ACCEPTED (FATAL BUG)")
        return None, None, "predicate-broken"
    print("  [x0-control] a=5.000000e-01 b=2.000000e+00 s=5.000000e-01 "
          "(p) -> REJECTED (correct)", file=sys.stderr)
    status("x0 negative control: REJECTED (correct)")

    rows_p, rows_q = [], []
    for a, b in _X0_P_SIDE_AB:
        e = min(1074, math.ceil(a * 1085))
        s = 2.0 ** -e
        if not _x0_certified(a, b, s, "p"):
            print(f"FATAL: x=0 p-side candidate failed certification: "
                  f"a={a!r} b={b!r} s={s!r} (E={e}) -- true root does not "
                  f"provably round to 0.0 (E formula may be wrong).",
                  file=sys.stderr)
            return None, None, "candidate-failed"
        rows_p.append((a, b, s, 0.0, "N"))
    for a, b in _X0_Q_SIDE_AB:
        s = math.nextafter(1.0, 0.0)
        if not _x0_certified(a, b, s, "q"):
            print(f"FATAL: x=0 q-side candidate failed certification: "
                  f"a={a!r} b={b!r} s={s!r} -- true root does not provably "
                  f"round to 0.0 (E formula may be wrong).", file=sys.stderr)
            return None, None, "candidate-failed"
        rows_q.append((a, b, s, 0.0, "N"))
    status(f"x0 stratum: {len(rows_p)} p-rows, {len(rows_q)} q-rows certified")
    return rows_p, rows_q, None


# ============================================================================
# Part 13: main.
# ============================================================================
def main():
    t_all = time.time()

    if not negative_controls():
        print("\nFATAL: negative control(s) were ACCEPTED -- the certifier "
              "is not rejecting known-bad rows. Aborting, nothing written.",
              file=sys.stderr)
        return 2
    status("negative controls: all rejected (correct).")

    rng = random.Random(SEED)
    status("building point set ...")
    ps = build_point_set(rng)

    status("evaluating oracle (resumable) ...")
    result, done = compute_all(ps)
    if not done:
        status("PARTIAL RUN: re-invoke this script to continue (checkpoint saved).")
        return 3
    rows_p, rows_q = result

    if len(rows_p) < 300 or len(rows_q) < 300:
        print(f"\nFAILED: row counts too low (p={len(rows_p)}, q={len(rows_q)}).",
              file=sys.stderr)
        return 1

    n_done, worst, worst_at = fast_vs_full_validate(random.Random(SEED ^ 0xB17A), n=40)
    # Hard gate (#13 N14.2): fast_vs_full_validate was print-only, and its
    # own try/except-continue silently shrinks the sample below n=40 on
    # any exception -- neither failure mode was enforced before this. The
    # fast route (bid.r1_value_mp) only SEEDS root-finding in the
    # r1-tiny/joint-tiny strata; bracket certification against the
    # audited oracle is independent of it, so this gate catches gross
    # route divergence, not a certification defect. 1e-20 is ~4 orders of
    # magnitude below double resolution (~1e-16) and far above the dps-80
    # mpf route agreement measured empirically -- do not loosen this
    # threshold without re-deriving that margin.
    if n_done < 30 or worst > 1e-20:
        print(f"\nFAILED: fast-vs-full validation gate tripped "
              f"(n_done={n_done}/40, worst={worst:.3e} at {worst_at}). "
              f"Aborting, nothing written.", file=sys.stderr)
        return 1

    zero_p, zero_q, zero_err = gen_x_zero_rows()
    if zero_err == "predicate-broken":
        print("\nFATAL: x=0 stratum negative control was ACCEPTED -- the "
              "certification predicate is not rejecting a known-ordinary "
              "root. Aborting, nothing written.", file=sys.stderr)
        return 2
    if zero_err is not None:
        print("\nFAILED: an x=0 stratum candidate failed certification. "
              "Aborting, nothing written.", file=sys.stderr)
        return 1
    rows_p = rows_p + zero_p
    rows_q = rows_q + zero_q

    write_rows(rows_p, rows_q)
    status(f"beyond-resolution certifications used: {N_BEYOND_RESOLUTION[0]}")
    status(f"route2 dual-checks: {N_ROUTE2_CHECKED[0]}, disagreements: {N_ROUTE2_DISAGREE[0]}")
    status(f"total generator runtime (this invocation): {time.time()-t_all:.1f}s")
    print(f"\nDONE. p-rows={len(rows_p)} q-rows={len(rows_q)} "
          f"beyond-resolution={N_BEYOND_RESOLUTION[0]} "
          f"route2-disagree={N_ROUTE2_DISAGREE[0]}/{N_ROUTE2_CHECKED[0]} "
          f"fast-vs-full worst={worst:.3e} (n={n_done})", file=sys.stderr)
    return 0


if __name__ == "__main__":
    if "--huge-corner-append" in sys.argv[1:]:
        sys.exit(huge_corner_append())
    sys.exit(main())
