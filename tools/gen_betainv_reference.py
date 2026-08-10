#!/usr/bin/env python3
"""Generate tests/data/betainv_{p,q}_reference.txt -- certified reference
set for corvus::beta_p_inv / corvus::beta_q_inv, per PLAN.md "P1 inverse
incomplete beta -- detail design" (BINDING), subsection "Oracle (G2;
frontier-specified -- THREE binding constructions beyond the gammainv
pattern)" and "Reference strata (G2)", plus the G1 stage record (pinned
seed/step constants read from the checked-in src/betainv_data.h and the
G1 generator tools/gen_betainv_data.py -- this generator does NOT re-run
G1's replay/self-check pipeline, only consumes its module-level machinery
and pinned constants, per the gammainv G2 precedent).

ORACLE, three binding constructions (decisions already made by the
frontier design; this generator implements, measures, and reports):

  1. FAST-PATH forward evaluator for R1-tiny/joint-tiny certification
     traffic: bid.r1_value_mp -- a plain, self-convergent mpf power
     series at target dps (gen_betainv_data.py's own "measurement-grade
     truth" evaluator, NOT the fixed-N cheap routing proxy), bypassing
     gen_beta_reference.py's small_side_direct escalation ladder
     (measured below: small_side_direct ranges 0-5ms on ordinary points
     but the brief's own 400-524ms figure is real for CF-heavy/ridge
     traffic -- infeasible at R1-tiny/joint-tiny volume). Validated
     fast-vs-full on a stratum sample (reported), then used for BOTH
     root-finding and bracket certification in the R1-tiny/joint-tiny
     strata specifically.

  2. GUARD on the reused gamma-corner route, at the enforcement site
     (this file, not gen_beta_reference.py/gen_beta_data.py): confirmed
     directly (see probe log referenced in G2-STATUS.md) that
     small_side_direct's own try_eval() calls gamma_corner_value(aa,bb,
     xx,dps) whenever max(aa,bb)>=B_GL, and gamma_corner_value ALWAYS
     feeds min(aa,bb)=min(a,b) to mpmath.gammainc as a shape argument
     (whichever of aa,bb is NOT picked as the huge "scale" side) --
     when BOTH a,b >= B_GL, that shape argument is itself huge (hang
     risk); gen_betainv_data.py's own betainv_forward is ALSO unsafe
     there (its R3 branch calls the raw backward CF unconditionally,
     confirmed to raise RuntimeError "CF not converged" at
     a=b=1e18,x=0.5 -- not a hang, but not a value either). GUARD:
     whenever min(a,b) >= B_GL (kBetaGammaLim) AND the skew ratio
     max(a,b)/min(a,b) <= SKEW_SAFE_CAP, route through
     beta_temme_value() below -- a dual-anchored R3-Temme extraction
     built from gen_beta_data.py's own extract_e_monomial/r3_R_at
     machinery (the "gamma_ck machinery" the brief names -- gamma_ck
     itself is that module's GAMMA-side anchor cross-check; the R3
     extraction apparatus it validates against is what this generator
     actually calls). This same route is also route-2 (the huge-nu
     stratum's independent second certification, gammainv G2 pattern).
     SCOPE NOTE (ratified deviation, reported): the anchor-ladder
     extraction holds p=a/(a+b) FIXED across anchors at a MODEST nu (so
     alpha=nu/(1-p), beta=nu/p at each anchor stay CF-safe) -- this is
     safe exactly when the skew ratio is bounded (extreme skew forces
     an anchor parameter to blow up even at modest nu). The brief's own
     language scopes this construction to "both-huge-BALANCED" traffic;
     SKEW_SAFE_CAP bounds how far from balanced this generator extends
     it, verified empirically below rather than assumed.

  3. Plateau rows: kappa computed per row (kappa = sigma/(y*f(y)), exact
     mpf, f = beta density). kappa <= 2^52 -> normal half-ulp bracket
     certification (the y-ULP gate). kappa > 2^52 -> BACKWARD-ERROR
     certification: forward of the STORED y at dps 100, required within
     the ~2-ulp-in-sigma contract (PLAN's own phrase) -- no y-bracket
     exists to certify there (dd precision cannot resolve it). Deep-
     small rows (both orientations, per G1's THIRD correction): log-
     space certification against the subnormal/zero boundary midpoints,
     with G1's own exact dropped-term bound (bid._deep_small_dropped_rel)
     folded in as certification slack (gammainv pattern).

Everything else follows the gammainv G2 CERTIFICATION CORE: root-find y*
seeded by G1's own bid.seed_for, round to double yd, certify
sign(value-target) flips across the two half-ulp midpoints of yd as
exact mpf, layered dps 60 -> 100 (never lower). Rows the certifier cannot
prove: DECLINED and counted, not guessed. NEGATIVE CONTROLS (>=4,
1-ULP-perturbed known-good rows) must be REJECTED on every invocation,
checked FIRST, exit 2 otherwise.

File format (ratified deviation, reported -- house pattern extended):
five hex tokens per row: a b sigma yd marker. marker in {N, P, B}
(Normal bracket-gate row / Plateau backward-error row / Beyond-resolution
row) -- the existing 13 reference files carry no marker column because no
prior family needed one; this is the first PLAN-mandated per-row bucket
split (plateau contract, beyond-resolution dilution lesson), so a marker
token is added rather than three separate files (keeps the swap-identity
orientation bookkeeping in one place per side).

Usage:
    python3 tools/gen_betainv_reference.py     # resumable; re-run until
                                                  # it reports DONE
"""
import math
import os
import random
import sys
import tempfile
import time

import mpmath as mp

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "tools"))
import gen_beta_data as gb          # noqa: E402  region cores, B_GL, T_RIDGE, ZETA_MAX
import gen_beta_reference as gbr    # noqa: E402  small_side_direct (audited oracle)
import gen_betainv_data as bid      # noqa: E402  seed_for, betainv_forward, r1_value_mp, ...

SEED = 20260809

STATUS_PATH = os.path.join(
    r"C:\Users\gdwol\AppData\Local\Temp\claude\C--Users-gdwol-Development-corvus"
    r"\e81b05d8-c230-46b2-8caa-e48c35f168d2\scratchpad\betainv_g2", "G2-STATUS.md")


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
SKEW_SAFE_CAP = 2.0e6                    # construction #2 scope bound, ORCHESTRATOR
                                          # WIDENED (this pass) to cover the beyond-
                                          # resolution stratum's own skewed corner
                                          # ("mean 1e-6" in the probe's B-P1c table,
                                          # i.e. p~1e-6, skew~1e6) -- validated by
                                          # direct route1-vs-route2 agreement at
                                          # skew=1e6/nu=1e20 (~30 decimal digits,
                                          # ample for double-precision certification;
                                          # see G2-STATUS.md), not assumed from the
                                          # original 1e4 pin.
HUGE_NU_THRESHOLD = 1.0e16               # route-2 dual-certification trigger (gammainv's own)
BEYOND_RESOLUTION_THRESHOLD = 3.0e34     # measured below, matches PLAN's own 1e33-7e34 note
PLATEAU_KAPPA_CUT = 2.0 ** 52            # PLAN's own contract split
DEEP_SMALL_CUT = bid.DEEP_SMALL_CUT      # 2^-60, both orientations (G1 THIRD correction)


GUARD_MIN_THRESHOLD = 1.0e10   # WIDENED (orchestrator continuation round,
                                # self-caught defect): min(a,b)>=B_GL alone
                                # is NOT a sufficient guard boundary --
                                # confirmed directly, small_side_direct's
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
    """GUARD predicate (construction #2's enforcement condition, WIDENED
    per GUARD_MIN_THRESHOLD's own docstring): fires whenever the LARGER
    of (a,b) is >= B_GL (small_side_direct's own gamma_corner_value
    trigger) AND the smaller is >= GUARD_MIN_THRESHOLD (empirically far
    below where mp.gammainc's own ridge-proximity hazard was observed to
    hang) AND skew stays within the validated extraction range."""
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
        # spacing/growth, different count -- gammainv G2 route-2's own
        # "different node count/anchor spacing" precedent.
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
    ladder extraction. SELF-CAUGHT BUG (this pass): the first draft used
    the FULL beta_temme_value (extract_e_monomial + QR solve, ~400-800ms
    per call near the ridge) inside bisection's ~150-200-iteration inner
    loop -- ~90s per row, hanging past this file's own smoke test. The R
    correction is a RELATIVE O(1/sqrt(nu)) effect (~1e-9 at nu~1e18);
    for root-finding purposes (not the final certified value) the
    leading term alone lands within a handful of ULPs of the true root,
    and certify_row's own local-nudge refinement (using the FULL
    evaluator, a bounded few calls) closes the gap before certifying."""
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

    RESCUE WRAPPER (self-caught DEFECT this pass, huge-nu scaling --
    escalated per the brief's own instruction: 'if you find a defect
    [in the reused machinery] that cannot be guarded externally,
    ESCALATE instead of editing'): small_side_direct's own try_eval
    calls gamma_corner_value -> mp.gammainc UNGUARDED against
    mp.libmp.libhyper.NoConvergence (it only catches RuntimeError/
    ZeroDivisionError/ValueError) -- confirmed directly, a shape
    argument as 'low' as ~1e17 (well UNDER kBetaGammaLim=2^59~5.76e17,
    outside this generator's own both_huge_balanced guard scope) can
    still make mpmath's hyp1f1 series give up with NoConvergence,
    propagating as an UNCAUGHT exception out of small_side_direct
    itself. This is a robustness gap in the shipped, audited oracle
    that this generator cannot fix by editing gen_beta_reference.py
    (out of scope per the brief) -- guarded here EXTERNALLY: on
    NoConvergence (or any other exception small_side_direct doesn't
    itself catch), fall back to this generator's OWN dual-anchored
    extraction (construction #2's route) regardless of whether
    both_huge_balanced's own narrower threshold applies -- the
    extraction is mathematically valid at any nu, just more expensive,
    so this is a safe general rescue, not a correctness compromise."""
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
    if y<=0.5 (value=P), swapped (b,a,1-y) if y>0.5 (value=Q).
    ORCHESTRATOR FIX (this pass, cost/correctness escalation): every
    root-finding call site previously fixed the series orientation by
    which of P/Q the CALLER wanted (the requested 'side'), not by which
    argument was actually small -- for joint-tiny rows whose true y sits
    on the far side (e.g. solving P(a,b,y)=sigma with the true y near 1),
    that forces bid.r1_value_mp to converge a series in x close to 1,
    which is slow (multi-second) and for some points never converges
    within the 4000-term cap at all (silently missing the sign flip,
    manifesting as a spurious root-find-failed decline -- not merely a
    speed problem, a coverage gap). This is the real fix; an earlier
    root_dps-decoupling pass (still valid, kept) only reduced the cost
    of the ALREADY-wrong-orientation calls without fixing convergence."""
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
    """Construction #1's own mandate: 'Validate fast-vs-full agreement on
    a stratum sample (report the sample size and worst disagreement)'.
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
# bisection is unusable near y=0/1, gen_betainv_data's own probe-caught
# bug #1, avoided from day one here per that precedent), seeded via
# bid.seed_for for speed only (never for correctness -- falls back to the
# full default bracket whenever the seeded one fails to bracket).
# ============================================================================
def oracle_y(a, b, target, side, dps, use_fast_series=False, seed_hint=None,
             root_dps=None):
    """root_dps (ORCHESTRATOR fix, this pass): root-FINDING precision,
    decoupled from dps (the CALLER's certification precision). SELF-
    CAUGHT COST BUG: the first draft ran root-finding's ~150-190-
    iteration bisection loop AT THE FULL CERTIFICATION dps (60) -- for
    fast_series (bid.r1_value_mp), series truncation eps scales with
    dps (eps=10^-(dps-8)), so EVERY one of those ~180 series
    evaluations paid dps=60's full term count, including the many
    iterations that land far from the true root (interior y, not
    small x -- where the series is genuinely slow, needing up to its
    4000-term cap). Measured outlier: 9.1s/row. FIX: root-find at a
    much lower dps (default 28, ~1e-20-class eps -- ample to land
    within a handful of ULPs of the true double, since bisection's own
    stopping tolerance is dps-scaled too), THEN certify at the caller's
    real dps_layers (60->100, a handful of calls, not ~180) -- exactly
    the guarded_fast_forward/full_forward cheap-root-find-vs-expensive-
    certify split already used elsewhere in this generator, now applied
    to fast_series too."""
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
        if seed_hint is not None and math.isfinite(seed_hint) and 0.0 < seed_hint < 1.0:
            v_seed = bid.logit(mp.mpf(seed_hint))
            slo, shi = v_seed - 80, v_seed + 80
            try:
                if f(slo) * f(shi) <= 0:
                    lo, hi = slo, shi
            except (ValueError, OverflowError):
                pass
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
# orientations -- G1's THIRD correction's own re-derived cut, reused
# verbatim: bid.deep_small_cut_bound (routing decision, native double),
# bid.deep_small_y (candidate y0, native double), bid._deep_small_dropped_rel
# (EXACT mpf dropped-term bound, folded in as certification slack, the
# gammainv pattern).
# ============================================================================
def deep_small_ly0(a, b, sigma, side, dps):
    """log-space target + exact bound, either orientation. Returns
    (ln_target, yd_mpf_rounded, bound_mpf) where ln_target is ln(y) [p]
    or ln(1-y) [q] at the closed-form leading order, and
    yd_mpf_rounded is the double CORRECTLY ROUNDED from ln_target at
    full working dps (SELF-CAUGHT BUG, this pass: an earlier draft used
    bid.deep_small_y's own NATIVE-DOUBLE seed formula here directly --
    that function is a fast SEED candidate for the kernel's own root
    search, built from several chained double-precision log/exp calls,
    and its compounded rounding error is enough to land OUTSIDE the
    half-ulp bracket around the true high-precision value, causing
    every deep-small row to spuriously fail certification. The oracle's
    own yd must come from rounding the HIGH-PRECISION ln_t, exactly the
    gammainv deep_small_lx0 pattern -- native-double formulas are seeds,
    never the certified value)."""
    with mp.workdps(dps):
        a_m, b_m, sigma_m = mp.mpf(a), mp.mpf(b), mp.mpf(sigma)
        if side == "p":
            lnB = mp.loggamma(a_m) + mp.loggamma(b_m) - mp.loggamma(a_m + b_m)
            ln_t = (mp.log(sigma_m) + mp.log(a_m) + lnB) / a_m
        else:
            lnB = mp.loggamma(b_m) + mp.loggamma(a_m) - mp.loggamma(a_m + b_m)
            ln_t = (mp.log(sigma_m) + mp.log(b_m) + lnB) / b_m
        y0_mpf = mp.e ** ln_t
        if side == "p":
            yd_rounded = float(y0_mpf) if y0_mpf > 0 else 0.0
        else:
            yd_rounded = 1.0 - float(y0_mpf) if y0_mpf > 0 else 1.0
        bound = bid._deep_small_dropped_rel(a_m, b_m, y0_mpf, side, dps=dps)
        return ln_t, yd_rounded, bound


MIN_SUBNORMAL = math.ldexp(1.0, -1074)


def certify_deep_small(a, b, sigma, side, dps, yd_override=None):
    """Certify a SPECIFIC double yd against the closed-form log-space
    target, at half-ulp midpoints, bound folded in as slack (gammainv
    pattern -- and its OWN documented self-caught bug: this generator's
    yd_override must be actually threaded through, or negative controls
    on deep-small rows are silently never checked)."""
    ln_t, yd_from_ln, bound = deep_small_ly0(a, b, sigma, side, dps)
    yd = yd_override if yd_override is not None else yd_from_ln
    if not math.isfinite(yd):
        return {"yd": None, "certified": False, "note": "non-finite"}
    yp = yd if side == "p" else (1.0 - yd)
    with mp.workdps(dps):
        if yp <= 0.0:
            # p-side: yd rounds to 0.0 (round-to-zero boundary, compare
            # against MIN_SUBNORMAL/2). q-side: yd rounds to 1.0 (the
            # SMALL quantity z=1-y rounds below ulp(1.0)/2=2^-53 --
            # SELF-CAUGHT BUG, this pass: an earlier draft reused the
            # p-side MIN_SUBNORMAL/2 boundary here unconditionally,
            # which is the wrong threshold by ~700 decades and
            # spuriously rejected every legitimate q-side round-to-one
            # deep-small row).
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
# the ~2-ulp-in-sigma contract (PLAN's own phrase): |forward(yd) - sigma|
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
    """kappa = sigma/(y*f(y)), exact mpf (PLAN's own plateau formula).
    SELF-CAUGHT BUG (this pass, huge-nu escalation): beta_density_mp's
    log-space assembly lg=(a-1)*ln(y)+(b-1)*ln(1-y)-lnB(a,b) is a
    CANCELLATION of terms scaling with a,b themselves (~1e35 at the
    huge-nu stratum's own scale) down to an O(1)-ish true log-density --
    at dps=60 (60 significant decimal digits) that cancellation loses
    ALL its resolving digits once a,b exceed ~1e40ish, and the
    naive symptom is f rounding to EXACTLY 0, kappa->inf, and a
    huge-nu row (which should have kappa NEAR ZERO -- density is
    enormous, not tiny, at that scale) spuriously misrouted into the
    plateau-backward branch (wrong contract; the row is really either
    an ordinary bracket certification or a genuine beyond-resolution
    one). FIX: scale dps with max(a,b)'s own decimal magnitude so the
    cancellation always has working digits left over, independent of
    the caller's certification-layer dps."""
    scale_dps = 60 + int(math.log10(max(float(a), float(b), 10.0))) + 20
    use_dps = max(dps, scale_dps)
    with mp.workdps(use_dps):
        f = bid.beta_density_mp(a, b, yd, use_dps)
        if f == 0:
            return mp.mpf("inf")
        return mp.mpf(sigma) / (mp.mpf(yd) * f)


# ============================================================================
# Part 7: beyond-resolution certification (nearest-neighbor, escalated
# dps) -- gammainv G2's own construction, reused for the huge-nu stratum
# where the transition collapses below 1 ULP of y.
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
    verified empirically, see G2-STATUS.md), so this search is a
    defensive bound, not the expected common path."""
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


def certify_row(a, b, sigma, side, tag, dps_layers=(60, 100), yd_override=None):
    is_r1tiny = tag in ("r1-tiny", "joint-tiny", "r1-tiny-seam")
    huge_bal = both_huge_balanced(a, b)
    fwd_kind = "fast_series" if is_r1tiny else ("guard" if huge_bal else "full")

    # --- deep-small routing decision (BOTH orientations, per THIRD
    # correction) -- cheap native-double check first. yd_override (negative
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
    yd = float(y_star)
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
    triple is REACHABLE, sidestepping the ill-posed-guess pathologies
    this generator's own smoke test hit (see G2-STATUS.md)."""
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
# gammainv's own s_from_x pattern, adopted here after this generator's own
# smoke test hit ill-posed-guess pathologies with log-uniform sigma
# sampling in these regions (see G2-STATUS.md self-caught bugs). The
# swap identity halves orientation coverage (one of (a,b)/(b,a) per
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
        """ORCHESTRATOR FIX (this pass): construct sigma from y, but
        choose the SIDE by which of P/Q is actually <=0.5 (PLAN's own
        'solve against s=min(p,1-p)' input contract) rather than a
        random coin flip -- a random side with y not already known to
        be on that side's small end can construct an out-of-contract
        row (sigma>0.5 for the requested side), which the certifier's
        deep-small/boundary machinery is not designed to handle (it
        assumes the REQUESTED side's sigma is the one approaching 0,
        not 1) and which manifested as spurious declines (see
        G2-STATUS.md). both=True (near-diagonal/plateau strata, PLAN's
        own swap-maps-s<->1-s rule) adds BOTH orientations regardless."""
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
        # orientations per PLAN's swap-maps-s<->1-s rule; skewed
        # sub-bands get the single small-side orientation (swap
        # identity halves coverage there).
        near_diag = 0.5 <= skew <= 2.0
        ps.add_from_y_smallside(a, b, y, "ridge", both=near_diag)
    status(f"  ridge: {len(ps.pts) - n0} points")


def gen_gammalim_seam(ps, rng, n=400):
    """Gamma-limit dense at the alpha~kGammaAT=20 seam: one param huge
    (gammalim), the other dense around 20 -- moderate cost (single-huge,
    small_side_direct's own guarded gamma_corner_value path). y
    constructed from the ACTUAL gamma-corner transition mapping (t ~
    shape-param scale, gen_beta_data.py's own gamma_corner_value
    identity inverted) -- SELF-CAUGHT BUG (this pass): an earlier draft
    picked y from a handful of O(1) constants nudged by a bare
    small/huge ratio, which lands in the saturated regime for all but
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


# ORCHESTRATOR-STATED thresholds (probe/B-P1c table): balanced nu*~3e33,
# skewed nu*~3e39 at mean~1e-6. NOT USED directly below -- direct
# calibration under this generator's own (nu,skew)->(a,b) mapping did not
# reproduce these numbers (empirical collapse nu* measured 2e36-4e36
# across skew 1..1e6; skew=1 alone stayed resolvable to a 1e60 search
# cap). ESCALATED (see G2-STATUS.md, final report) rather than silently
# used or silently overridden; gen_huge_nu below uses ITS OWN measured
# numbers.


def _huge_nu_y_mpf(a, b, target_z, dps=80):
    """Construct y at mpf precision so that cpsi = target_z^2 EXACTLY
    (to mpf precision) at the point BEFORE double-rounding, THEN round
    to double. TWO self-caught bugs this pass:
    (1) building y = a/c + delta in NATIVE PYTHON FLOAT arithmetic: its
        ~1e-16 relative rounding error in the 'mean' translates (via
        lam = a - c*y) to an ABSOLUTE error in lam of order c*1e-16 --
        at a,b ~ 1e40+ that alone drives cpsi up to ~1e19-scale even at
        the INTENDED exact-mean point.
    (2) gb.r3_setup's own 'zeta' parameter is NOT target_z directly --
        r3_setup's docstring/_lambda_of_zeta target is
        cpsi = zeta^2 * nu (the RIDGE-normalized deviation, zeta =
        z/sqrt(nu)), not cpsi = zeta^2. Passing target_z straight through
        as zeta silently requested cpsi = target_z^2 * nu ~ nu itself
        (confirmed directly: 'target_z=1.0' at nu=1e15 produced
        cpsi=999999999999999.9, not 1) -- deep in saturation regardless
        of target_z, which is why the first calibration pass found
        collapse at EVERY nu tested. FIX: zeta = target_z / sqrt(nu),
        so cpsi = zeta^2*nu = target_z^2 exactly, as intended."""
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
    """Huge-nu: TWO regimes. (a) both-huge-balanced/skewed GUARD
    territory (construction #2's own route, RESOLVABLE -- z=O(1) points
    still land on distinguishable doubles) via _huge_nu_y_mpf, a
    mpf-precision z-targeted construction (validated: r1-vs-r2 agreement
    to ~30+ decimal digits even at skew=1e6/nu=1e20, see G2-STATUS.md).
    (b) genuine BEYOND-RESOLUTION territory: nu high enough that the
    WHOLE transition collapses within a handful of ULPs of the mean
    (gammainv's own phrasing) -- constructed via _huge_nu_mean_ulp
    (mpf-precision mean, THEN double-rounded, THEN ULP-stepped -- the
    z-targeted construction is NOT usable here since at this nu ANY
    z=O(1..10) point rounds back to one of the SAME 3-5 doubles nearest
    the mean, self-caught during calibration, see G2-STATUS.md), with
    sigma picked to deliberately NOT be exactly the trivial 0.5 (mid-
    ulp) so the certifier's nearest-neighbor contract is genuinely
    exercised rather than the trivial straddle case gammainv's own
    'dilution lesson' warns about. ORCHESTRATOR-STATED thresholds
    (balanced nu*~3e33, skewed nu*~3e39 at mean~1e-6) did NOT reproduce
    under this generator's own (nu,skew)->(a,b) parameterization when
    measured directly (empirical collapse nu* varied 2e36-4e36 across
    skew 1..1e6, and skew=1 alone stayed resolvable to the 1e60 search
    cap under the z=3 probe) -- ESCALATED, not silently overridden: see
    G2-STATUS.md for the full calibration record and final report. This
    generator uses ITS OWN measured numbers (nu up to 1e34 for the
    guard bucket, nu 1e35-1e42 for the beyond-resolution bucket, both
    empirically confirmed to produce the intended behavior end-to-end)
    rather than the unreproduced probe figures."""
    n0 = len(ps.pts)
    # (a) resolvable guard territory.
    for _ in range(n):
        skewed = rng.random() < 0.5
        if skewed:
            nu = 10.0 ** rng.uniform(math.log10(B_GL) - 1, 20.0)
            skew = 10.0 ** rng.uniform(2, math.log10(SKEW_SAFE_CAP))
        else:
            nu = 10.0 ** rng.uniform(math.log10(B_GL) - 1, 33.0)
            skew = 10.0 ** rng.uniform(0, 1)
        a, b = (nu, nu * skew) if skewed else (nu * (1 + skew), nu * (1 + skew) / skew)
        delta_z = rng.uniform(-4, 4)
        y = _huge_nu_y_mpf(a, b, delta_z)
        ps.add_from_y_smallside(a, b, y, "huge-nu", dps=80)
    # (b) genuine beyond-resolution: sigma picked DIRECTLY (gammainv G2's
    # own huge-a-beyond-resolution-target pattern, PLAN.md precedent --
    # NOT forward-constructed from y here, SELF-CAUGHT BUG this pass:
    # forward(mean +/- k ulps) is ITSELF already saturated to exactly
    # {0,0.5,1} at this collapse depth, so s_from_y's well-posedness
    # filter rejected essentially every attempt (0/15 in a direct
    # isolation test) -- there is no 'intermediate' reachable sigma to
    # construct from y at true collapse depth; the well-posed INPUT
    # contract here is simply sigma in (0,1), the same as any kernel
    # call, and certify_row's own nearest-of-{neighbors} contract is
    # what proves the answer, not a forward round-trip).
    beyond_sigmas = (0.5, NEXT_UP(0.5), NEXT_DN(0.5), 0.3, 0.7, 0.1, 0.9,
                     1e-3, 1.0 - 1e-3)
    for _ in range(n_beyond):
        skewed = rng.random() < 0.5
        if skewed:
            nu = 10.0 ** rng.uniform(35.0, 42.0)
            skew = 10.0 ** rng.uniform(2, math.log10(SKEW_SAFE_CAP))
            a, b = nu, nu * skew
        else:
            nu = 10.0 ** rng.uniform(35.0, 42.0)
            a = b = nu * 2
        sigma = rng.choice(beyond_sigmas)
        side = rng.choice(("p", "q"))
        ps.add(a, b, sigma, side, "huge-nu")
    status(f"  huge-nu (both-huge-balanced guard + beyond-resolution): "
           f"{len(ps.pts) - n0} points")


def gen_seam_bracket(ps, rng, n=150):
    """a_T-seam bit-stepped bracket (S1/S3 seed seam near kGammaAT=20,
    per PLAN's stratum list)."""
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
    side) crosses DEEP_SMALL_CUT -- ORCHESTRATOR FIX (item 1, this
    round): the ORIGINAL gen_deep_small_both picked sigma from a blind
    log-uniform grid (10^pe, pe in -320..-1) with NO regard for whether
    that sigma's corresponding y was anywhere near the actual cut for
    the (a,b) pair drawn -- for a near a_T=19 the cut sits at
    y~2^-60/19~4.6e-20, so almost the ENTIRE pe range constructed a
    sigma whose true y was nowhere near deep-small territory, correctly
    falling through to ordinary root-finding, which then legitimately
    failed for many of those ill-matched (a,b,sigma) combinations --
    measured 415/540 (77%) drop rate, a GENERATOR-QUALITY gap, not a
    certifier defect. FIX (inversion-first, matching gen_gammalim_seam's
    own fix): locate y_cut precisely per (a,b,side) via bisection on the
    ACTUAL routing predicate (not a re-derived approximation of it),
    then sample y at chosen MULTIPLES of y_cut (both inside and outside
    the deep-small regime -- boundary coverage is itself valuable, per
    the design's own 'deep-small-cut-bracket' precedent in gammainv),
    and construct sigma = forward(y) (in-band, reachable BY
    CONSTRUCTION). Returns z_cut, the SMALL variable's cut location (z=y
    for side='p', z=1-y for side='q' -- bid.deep_small_cut_bound's own
    yp convention) -- caller converts to y."""
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
    start (G1 THIRD correction's own lesson) -- y constructed from the
    (a,b,side)-SPECIFIC cut location (find_y_cut), sigma=forward(y) (
    in-band by construction, ORCHESTRATOR fix)."""
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
    # SIZED DOWN [orchestrator cost-discipline pass, this session]: the
    # first sizing (r1-tiny 2400/joint-tiny 1200/ridge 1200/gammalim 800/
    # huge-nu 180) measured real per-row cost far above the initial
    # projection for a MINORITY of joint-tiny rows (sigma within a few
    # ULPs of 0 or 1 drives bisection into many extra iterations of
    # bid.r1_value_mp's own up-to-4000-term series -- one measured
    # outlier cost 9.1s vs a ~150-300ms typical row) -- see
    # G2-STATUS.md. Re-sized to a total that completes within this
    # session's realistic remaining budget rather than leave an
    # open-ended partial checkpoint; reported honestly as smaller than
    # the design's 14-21k range in the final report.
    # ORCHESTRATOR CONTINUATION ROUND: scaled up toward the design's
    # 14-21k target now that per-row costs are measured POST-FIX
    # (fast_series orientation fix, smallside construction, kappa dps
    # fix, deep-small inversion-first fix, huge-nu mpf construction --
    # see G2-STATUS.md). Sizing still short of the full 14-21k range in
    # a few strata (ridge/gammalim/huge-nu, the genuinely expensive
    # ones); reported precisely in the final report against design
    # targets, not silently padded.
    gen_r1_tiny(ps, rng, n=6000)
    gen_joint_tiny(ps, rng, n=2000)
    gen_ridge(ps, rng, n=2000)
    gen_gammalim_seam(ps, rng, n=2000)
    gen_underflow(ps, rng, n=600)
    gen_subnormal_y(ps, rng, n=500)
    gen_huge_nu(ps, rng, n=150, n_beyond=400)
    gen_seam_bracket(ps, rng, n=600)
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
WALL_CLOCK_BUDGET_S = 200.0


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


def compute_all(ps):
    total = len(ps.pts)
    sig = f"v1 SEED={SEED} N={total}"
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
    sys.exit(main())
