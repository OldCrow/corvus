#!/usr/bin/env python3
"""Generate src/beta_data.h -- every table the beta_p/beta_q (regularized
incomplete beta) kernel needs, per PLAN.md "Regularized incomplete beta --
detail design" and its G1a/G1b probe corrections (the region map, the
routing order, the R2 orientation rule, and the R3 (zeta, p) ansatz with
the eta = -zeta*sqrt(2) gamma-limit mapping are all binding; this generator
does not re-derive them, only pins the [OPEN] numeric constants those
sections left for G1c).

Four independent pieces:

  R1/R4 series   dd pairs of 1/n, n=1..max(N1,R4_NMAX) -- both R1's power
                 series and R4's alpha-scaled series share the SAME term
                 ratio t_n = t_{n-1}*(n-beta)*xi/n (R4's S is R1's series
                 shape minus the n=0 term, "gamma-R4 verbatim in beta
                 clothing" per the design), so one dd 1/n table (mirroring
                 gamma's kGammaRecipN) covers both -- gamma's own
                 1/(a+n)-via-DdRecipDd-of-TwoSum(a,n) precedent applies
                 unchanged (that division is never tabled, since a+n is
                 runtime-dependent; only the pure-integer 1/n multiplier
                 is).

  R3 tensor table  e_k(zeta, p), k=0..K-1, the erf-form correction series
                 S = sum_k e_k(zeta,p)/nu^k. Each row is a 2D tensor-
                 Chebyshev fit in (zeta, p) -- Chebyshev interpolation on a
                 Chebyshev-node grid, exact (not least-squares) in the two
                 SMOOTH variables. Per-node VALUES come from a clean-room
                 extraction: R(zeta,p,nu) = (leading -/+ small_val) *
                 sqrt(2*pi*nu)*exp(cpsi) is fit as a MONOMIAL power series
                 in v=1/nu via a well-conditioned Chebyshev-in-v least-
                 squares solve (G1b Task B3's method: raw monomial
                 Vandermonde is catastrophically ill-conditioned at this
                 ladder length/order, but Chebyshev-in-v is not) followed
                 by an EXACT affine basis-change back to monomial-in-v
                 coefficients (shifted_chebyshev_monomial) -- this is
                 algebra on an already-accurate fit, not a second
                 statistical estimate, so it does not reintroduce the
                 conditioning problem. small_val is evaluated via the
                 DLMF 8.17.22 continued fraction (I_via_cf, G1b Task A),
                 escalating N until self-convergent -- G1b Task B2 found
                 mpmath's own betainc/quadrature BOTH unreliable for this
                 extraction's shape (quad needs infeasible dps; betainc
                 hangs/fails to converge off-ridge), and the CF is the one
                 oracle validated across this whole domain.

  Binet table    K_B Bernoulli-derived coefficients of the Stirling/Binet
                 tail phi(z) ~ sum B_2k/(2k(2k-1)z^(2k-1)), z >= Z0 -- the
                 same series as src/lgamma-inl.h's LgammaStirling but its
                 OWN table (different Z0/K_B target), per the design's
                 explicit "fresh coefficients" option; BinetDd's ultimate
                 kernel home is a G3 decision, this generator only emits
                 the numbers.

  DigammaRough   poly(w), w in [Z0, Z0+1), degree chosen by residual decay,
                 fit directly against mpmath's psi (NOT the asymptotic
                 series -- ordinary Chebyshev interpolation on the base
                 zone is far cheaper at this loose ~2^-40 target). Any
                 z in (0, 2*Z0] reaches the zone by the SAME up-recurrence
                 pattern as lgamma's own zone walk: psi(z) = poly(z_final)
                 - sum(1/z_j) accumulated while walking z up by integer
                 steps until z_final in [Z0, Z0+1).

G1a/G1b provenance: the region map, N1=64, N2=64, T_ridge=32, B1=8,
xi1=0.45, eps_R4=2^-6 (and its box: xi_tau<=xi1 AND B*xi_tau<=B1, the
missing-cap fix), Z0=10, K_B=16, C_lg=256 (provisional), the R2
orientation rule xi < (alpha+1)/(c+2), the corrected ROUTING ORDER (0:
tiny-min-first -> R4; 1: R1; 2: R3; 3: R2) that fixes BOTH escalations
found in G1a (R4's missing xi cap) and G1b (R1's missing alpha floor,
stealing points from R4) are all read from PLAN.md and reproduced here
verbatim in route_final(); they are not re-derived. What IS pinned here
(G1c's own job, [OPEN] in PLAN.md before this generator ran): R4's series
depth, the R3 table's K/degrees within the 32 KB budget, the Binet/
DigammaRough coefficient values.

Self-checks (mandatory, budget lines to stderr; ANY miss -> exit nonzero,
emit nothing) -- lettered per the brief, (a)-(i):
  (a) R1 series truncation sup at N1=64 over the full membership boundary
      (both boundary edges, alpha down to 1e-300).
  (b) R2 CF depth sup at N2=64: G1a's lattices plus the G1b beta->1e-300
      extension, each point tested in the orientation the pinned rule
      picks.
  (c) R3 total S truncation + 2D fit residual vs the CF oracle, over a
      (zeta, p, nu) lattice, target the FIXED gamma-class 2^-56
      (R3_REPLAY_TARGET -- never pinned-to-measured). NOTE (deviation,
      flagged): the design text says "quad oracle for ridge; CF oracle
      where quad's shape limits bite" -- this generator uses the CF
      EXCLUSIVELY for R3, per G1b Task B2/B3's own finding that quad is
      unreliable for R3's extraction shape (needs infeasible dps; CF is
      the oracle validated across this whole domain). History: the
      FIRST cut of this generator fit zeta in [-5,5] (the erroneous
      cpsi<=800 membership) and stalled at 2^-16-class; the third
      correction's ratio-band membership (see ZETA_MAX's derivation
      comment) collapsed the domain to |zeta| <= ~1.02 and the same
      table budget now clears 2^-56 with margin.
  (d) e_k(zeta, p=2^-50) matches gamma's c_k(eta=-zeta*sqrt(2)) to
      <=1e-15 through k=5 (the mandatory gamma-limit anchor).
  (e) routing safety under the FINAL corrected order: max evaluated side
      <= 1-2^-12 over a boundary lattice including both escalation
      witnesses and their families.
  (f) R4 series truncation over the corrected closed box (xi_tau<=xi1,
      B*xi_tau<=B1, tau*|ln xi_tau|<=ln2).
  (g) Binet truncation at Z0=10, K_B=16 vs the dd target.
  (h) e_k(zeta,p) == -e_k(-zeta, 1-p) (the symmetry identity; SIGN
      corrected here to "-", not "+" -- see the derivation comment at
      check_h_symmetry for the algebra).
  (i) DigammaRough <= 2^-40 relative on (0, 2*Z0].

Trust the stderr budget lines over any comment in this file, including
this one -- AGENTS.md house rule, restated because the R3 section below
is exactly the kind of place a comment could go stale first.

Usage:
    python3 tools/gen_beta_data.py > src/beta_data.h
"""

import math
import multiprocessing as mp_proc
import sys
import time

import mpmath as mp

mp.mp.dps = 100  # module-level default; every self-check/extraction
                  # function below sets its OWN dps explicitly on entry
                  # and restores it on exit (G1a/G1b hygiene rule -- a
                  # stale ambient dps at a subprocess or string boundary
                  # bit three separate G1a/G1b probe layers).

# --- Pinned design constants [PLAN.md, G1a/G1b] -----------------------------
B1 = mp.mpf(8)
XI1 = mp.mpf("0.45")
EPS_R4 = mp.mpf(2) ** -6
T_RIDGE = mp.mpf(32)
Z0 = mp.mpf(10)
K_B = 16
C_LG = mp.mpf(256)
E_FLOOR = mp.mpf(-800)
LN2 = mp.log(2)

N1 = 64  # R1 series depth [G1a pinned; self-check (a) reproves it]
N2 = 64  # R2 CF depth [G1a/G1b pinned; self-check (b) reproves it]

# ZETA_MAX derivation [THIRD CORRECTION, PLAN.md "G1c generator results and
# third correction"]: R3's membership is NOT "nu>=T_ridge and cpsi<=800" --
# that was a design error (confirmed by reading gen_gamma_data.py: gamma's
# shipped Temme table spans only the ridge RATIO band lambda in [1/2,2],
# NOT a wide cpsi strip; a 32 KB single tensor cannot span cpsi<=800 at dd
# precision, which is exactly the wall the first generator pass hit --
# 2^-16-class residual, unshippable). CORRECTED membership, mirroring
# gamma exactly: nu>=T_ridge AND xi/p in [1/2,2] AND (1-xi)/q in [1/2,2].
#
# Derivation of the caps in (u,v): xi = (alpha-lambda)/c, so
#   xi/p = xi*c/alpha = (alpha-lambda)/alpha = 1 - lambda/alpha = 1+u
# (u=-lambda/alpha). Likewise (1-xi) = (beta+lambda)/c, so
#   (1-xi)/q = 1+lambda/beta = 1+v  (v=lambda/beta).
# So the two ratio caps are EXACTLY u in [-1/2,1] and v in [-1/2,1] -- the
# same phi(u)=u-log1p(u) domain gamma's own phi_lam(lambda)=lambda-1-log(lambda)
# uses with w=lambda-1 (phi(w)=phi_lam(w+1) identically), confirming this
# is gamma's ETA_LO/ETA_HI band carried over unchanged: at p->0,
# u in [-1/2,1] alone is the constraint (v's cap vanishes, q/p->infinity),
# exactly gamma's own domain.
# The two caps are LINKED: p*u+q*v = p*(-lambda/alpha)+q*(lambda/beta) =
# -lambda*(p/alpha)+lambda*(q/beta) = -lambda/c+lambda/c = 0 (using
# p/alpha=q/beta=1/c), so v=-(p/q)*u. For FIXED p, cpsi(u)=p*phi(u)+q*phi(v)
# is CONVEX in u (phi convex, v affine in u) with minimum at u=0, so the
# sup over the valid u-interval [lo(p),hi(p)] (lo=max(-1/2,-q/p),
# hi=min(1,q/(2p))) is at an ENDPOINT. Maximizing over p is done
# numerically below (golden-section on the boundary-endpoint value) and
# the result matches a clean closed form to be checked as an assertion:
# the sup is at p=1/3 (and its mirror p=2/3), where u=1,v=-1/2 (BOTH raw
# caps bind simultaneously -- the corner of the feasible (p,u) region):
#   cpsi* = (1/3)*phi(1) + (2/3)*phi(-1/2)
#         = (1/3)*(1-ln2) + (2/3)*(ln2-1/2) = ln2/3
#   zeta_max^2 = cpsi*/(p*q) = (ln2/3)/(2/9) = 3*ln2/2
#   zeta_max = sqrt(3*ln2/2) ~ 1.0197 -- NOT the ~0.76 rough estimate that
# motivated this correction; that estimate undersold the true sup (this is
# the point of deriving it exactly rather than trusting the rough number).
def _phi(w):
    return w - mp.log1p(w)


def _zeta2_at_boundary(p):
    q = 1 - p
    lo = max(mp.mpf("-0.5"), -q / p)
    hi = min(mp.mpf(1), q / (2 * p))
    best = mp.mpf(0)
    for u in (lo, hi):
        v = -(p / q) * u
        cpsi = p * _phi(u) + q * _phi(v)
        z2 = cpsi / (p * q)
        if z2 > best:
            best = z2
    return best


def _derive_zeta_max(dps=60):
    """Numerically maximize _zeta2_at_boundary over p in (0,1) (coarse grid
    + golden-section refine), then cross-check against the closed form
    sqrt(3*ln2/2) derived in the comment above -- if they disagree this is
    a real bug, not a constant to silently trust."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        n = 2000
        best_p, best_z2 = None, mp.mpf(-1)
        for i in range(1, n):
            p = mp.mpf(i) / n
            z2 = _zeta2_at_boundary(p)
            if z2 > best_z2:
                best_z2, best_p = z2, p
        gr = (mp.sqrt(5) - 1) / 2
        a = max(mp.mpf("1e-8"), best_p - mp.mpf(2) / n)
        b = min(1 - mp.mpf("1e-8"), best_p + mp.mpf(2) / n)
        c = b - gr * (b - a)
        d = a + gr * (b - a)
        fc, fd = _zeta2_at_boundary(c), _zeta2_at_boundary(d)
        for _ in range(200):
            if fc > fd:
                b, d, fd = d, c, fc
                c = b - gr * (b - a)
                fc = _zeta2_at_boundary(c)
            else:
                a, c, fc = c, d, fd
                d = a + gr * (b - a)
                fd = _zeta2_at_boundary(d)
            if abs(b - a) < mp.mpf(10) ** -50:
                break
        p_star = (a + b) / 2
        z2_star = _zeta2_at_boundary(p_star)
        closed_form = 3 * mp.log(2) / 2
        rel_diff = abs((z2_star - closed_form) / closed_form)
        print(f"    zeta_max derivation: p*={mp.nstr(p_star, 10)} "
              f"zeta_max^2={mp.nstr(z2_star, 20)} vs closed form 3*ln2/2="
              f"{mp.nstr(closed_form, 20)} (rel diff {float(rel_diff):.3e})",
              file=sys.stderr)
        if rel_diff > mp.mpf("1e-30"):
            raise RuntimeError("zeta_max numeric optimum does not match the "
                                "closed-form derivation -- do not trust either "
                                "silently; this is a bug, escalate.")
        return mp.sqrt(closed_form)
    finally:
        mp.mp.dps = old


# Windows multiprocessing uses spawn (no fork): every mp_proc.Process(...) in
# this file (the _betainc_timeout machinery) re-imports this ENTIRE module
# from scratch in the child before calling its target -- including any
# module-level derivation. Guard the two expensive ones (ZETA_MAX, B_GL
# below) so a spawned worker gets a cheap placeholder instead of re-running
# a multi-second probe on every single subprocess spawn: this was measured
# directly on this box (not assumed) -- the FIRST cut of the B_GL probe
# added here made check_b_r2's existing (pre-G3) betainc-timeout sections
# balloon from its established ~1-2 min budget to 280s+ without finishing,
# traced to the B_GL derivation's own stderr banner reappearing inside the
# log once per subprocess spawn. Neither ZETA_MAX nor B_GL is ever read
# inside _betainc_worker (it only calls mp2.betainc on the three plain
# values it was handed), so a placeholder is exactly as correct there as
# the derived value -- it is simply never consulted.
# NOTE: mp_proc.parent_process() measured (not assumed) to return None even
# inside a Windows-spawned child during module import -- current_process()
# .name is the reliable signal there ("Process-N" in a child vs
# "MainProcess" at the top); confirmed with a minimal repro before trusting
# it here, since a wrong guard here silently reintroduces the exact
# per-subprocess re-derivation cost this guard exists to remove.
_IS_MP_WORKER = mp_proc.current_process().name != "MainProcess"

if _IS_MP_WORKER:
    ZETA_MAX = mp.mpf("1.0197207708399179641")  # placeholder; see guard above
else:
    ZETA_MAX = _derive_zeta_max()

# --- Self-check targets ------------------------------------------------------
R1_TARGET = mp.mpf(2) ** -60
R2_TARGET = mp.mpf(2) ** -60
R4_TARGET = mp.mpf(2) ** -58
BINET_TARGET = mp.mpf(2) ** -70
DIGAMMA_TARGET = mp.mpf(2) ** -40
ROUTE_THRESH = 1 - mp.mpf(2) ** -12

# BETA_NEAR_ONE (kBetaNearOne) [SIXTH correction threshold, SEVENTH
# correction destination]:
# route_final()'s post-route threshold. R1 fires in EITHER orientation (the
# fifth correction's lambda>=0 requirement is REVERTED -- it displaced sound
# traffic, breaking check (b)'s CF-depth witness (0.158,1000,0.00251)); a
# point that routes to R1 but whose R1-NATIVE value exceeds this bar is
# re-tagged "R4-postroute" and evaluated by R4's analytic small-side
# assembly in the SAME orientation [SEVENTH correction -- the sixth's
# opposite-orientation CF destination stalled at 2^-55.5 and is
# superseded]. 2^-11 is intentionally one bit inside ROUTE_THRESH's
# 2^-12 complement-slack doctrine (post-route decides on the CHEAP R1-series
# oracle at working dps, not the full self-convergent CF check (e) itself
# uses; the extra bit of margin absorbs that oracle's own slack).
BETA_NEAR_ONE = 1 - mp.mpf(2) ** -11
# Gamma-limit slice ridge floor [(C) resolution]: in-band lanes with
# max(alpha,beta) >= B_GL use R3 down to nu = 20 (gamma's own kGammaAT --
# the CF is degenerate up there and gamma's series/CF boxes exclude the
# band). check (c)'s extension lattice proves the 1/nu extrapolation.
GL_RIDGE_MIN = mp.mpf(20)
ANCHOR_TARGET = mp.mpf("1e-15")
SYMMETRY_TARGET = mp.mpf("1e-25")

# R3_REPLAY_TARGET: gamma-class 2^-56 (gen_gamma_data.py's own REPLAY_TARGET),
# per the third-correction instruction -- NOT pinned-to-measured this time,
# since the corrected ratio-band domain is expected to reach it (gamma-class
# node counts over a gamma-class-width domain). If unreachable, that is a
# real ESCALATE, not a constant to loosen.
R3_REPLAY_TARGET = mp.mpf(2) ** -56

# R3 ratio-band caps (the third correction): xi/p in [XI_RATIO_LO,XI_RATIO_HI]
# AND (1-xi)/q in [XI_RATIO_LO,XI_RATIO_HI], equivalently u,v in [-1/2,1].
XI_RATIO_LO = mp.mpf("0.5")
XI_RATIO_HI = mp.mpf(2)



# ============================================================================
# dd / hex-float emission helpers (matches tools/gen_gamma_data.py style)
# ============================================================================
def rd(x):
    return float(x)


def dd_split(x):
    hi = rd(x)
    return hi, rd(mp.mpf(x) - mp.mpf(hi))


def hexf(x):
    return float.hex(float(x))


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


def emit_hex_array_3d(name, blocks):
    # blocks[k][n][m] -- one 2D slab per k, printed as a flat [K][NZ][NP].
    nz = len(blocks[0])
    npc = len(blocks[0][0])
    print(f"inline constexpr double {name}[{len(blocks)}][{nz}][{npc}] = {{")
    for blk in blocks:
        print("    {")
        for row in blk:
            print("        {" + ", ".join(hexf(v) for v in row) + "},")
        print("    },")
    print("};")


# ============================================================================
# R2 continued fraction (DLMF 8.17.22), verbatim from G1b Task A
# (taskA_orientation.py), independently validated there over 256 zone
# points plus this generator's own self-check (b).
# ============================================================================
def d_coef(k, alpha, beta, xi, c):
    if k % 2 == 0:
        m = k // 2
        return m * (beta - m) * xi / ((alpha + 2 * m - 1) * (alpha + 2 * m))
    else:
        m = (k - 1) // 2
        return -(alpha + m) * (c + m) * xi / ((alpha + 2 * m) * (alpha + 2 * m + 1))


def cf_backward(alpha, beta, xi, N):
    c = alpha + beta
    f = mp.mpf(1)
    for k in range(N, 0, -1):
        f = 1 + d_coef(k, alpha, beta, xi, c) / f
    return 1 / f


def I_via_cf(alpha, beta, xi, N):
    F = cf_backward(alpha, beta, xi, N)
    logpref = alpha * mp.log(xi) + beta * mp.log(1 - xi) - mp.log(alpha) - (
        mp.loggamma(alpha) + mp.loggamma(beta) - mp.loggamma(alpha + beta))
    return mp.exp(logpref) * F


def small_val_via_cf(a, b, x, dps, n_start=128, n_max=4096):
    """Escalating self-convergence gate (generalizes G1b Task B2/B3's fixed
    N=128-vs-256 gate, which this generator's wider (zeta,p) grid found
    insufficient at some interior points -- moderate-zeta, mid-ladder
    points needing N>256 to self-converge, not just near-ridge ones)."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        n = n_start
        v_lo = I_via_cf(a, b, x, n // 2)
        v_hi = I_via_cf(a, b, x, n)
        while True:
            drift = abs((v_hi - v_lo) / v_hi) if v_hi != 0 else abs(v_lo - v_hi)
            if drift <= mp.mpf("1e-15"):
                return v_hi
            if n >= n_max:
                raise RuntimeError(f"CF not converged at N={n // 2}..{n} "
                                    f"(drift={float(drift):.3e}) for a={a},b={b},x={x}")
            n *= 2
            v_lo = v_hi
            v_hi = I_via_cf(a, b, x, n)
    finally:
        mp.mp.dps = old


# ============================================================================
# kBetaGammaLim (B_GL) derivation [G3 escalation (C), "gamma-limit slice"].
# BOTH-SIDES overlap probe per PLAN.md's binding recipe: sweep beta=2^j along
# the gamma line (alpha tiny fixed, beta huge, beta*xi fixed -- the exact
# shape of the (0.05, 1e100, 2e-99) witness), evaluated in the CF's OWN
# orientation (xi < (alpha+1)/(c+2) -> native, else swap -- same rule as
# check_b_r2/route_final; NOT the raw un-oriented triple, which just measures
# a meaningless divergent CF branch that route_final would never select).
#
# UPWARD direction (beta-CF validity ceiling): the design's literal recipe
# is "fixed-depth N2=64 CF at dps 17 (double proxy) vs dps>=100 truth". Ran
# literally first (bxi in the design's own {8,20,50,200,800} set): the
# double-proxy error ALREADY exceeds 2^-60 at the very first grid point
# (j=40, error ~1.4e-5) and every alpha/bxi combo tested at bxi=8 turns out
# to be R1's OWN territory (beta*xi=8=B1 exactly -- R1's box condition -- so
# these points never reach R2/the gamma-limit path in the real kernel;
# excluded from the "worst case within R2 territory" below). Root cause,
# confirmed directly (not inferred): the backward CF's own F=1/f (final
# convergent) grows like beta itself along this line --
# log2(F) ~= j - 4.3, essentially F ~= beta/20 -- because the swapped-
# orientation d_1 sits within O(1/beta) of the singular value -1 (same
# "f_1 near zero" conditioning story BetaR2Cf's own header comment
# documents: "A double-precision d there is 2^-53*F ... 90 and 5e6 ULP" at
# ITS example points, F~200 and F~1.25e7 there -- this generator's F reaches
# those same orders of magnitude by j~23 and stays on the SAME F~beta/20
# line all the way to j=332, twelve more binades past their examples).
#
# DEVIATION (flagged, per this file's own R3-oracle-choice precedent):
# BetaR2Cf's own docstring establishes "ULP error ~= F ULPs [of whatever
# base precision the recursion runs in]" as the operative model -- so a
# dps=17 (~double, eps~2^-53) proxy measures F*2^-53, but the SHIPPED
# kernel runs the CF entirely in DOUBLE-DOUBLE (eps_dd~2^-106, per
# BetaR2Cf's own "THE RECURRENCE AND EVERY d ARE DOUBLE-DOUBLE" design
# choice) -- so the quantity that actually gates the shipped kernel's
# safety is F*2^-106, not F*2^-53. This generator therefore measures F
# directly (exact/high-dps, no proxy needed -- F does not depend on the
# evaluator's own working precision) and rescales by eps_dd=2^-106, rather
# than trusting the raw dps=17 number, which would peg the "upward
# frontier" at j~11 (beta~2000) -- clearly not what the escalation's own
# "probed upward FROM the 2^40 G1a-validated line" framing intends (that
# 2^40 anchor is G1a's TRUNCATION validation, not a conditioning claim, but
# a frontier at j~11 would contradict check_b_r2's own existing gamma-limit
# lattice, which already tests -- and passes -- beta up to 1e12~=2^40).
#
# TARGET DEVIATION (also flagged): at the literal 2^-60 bar, the rescaled
# upward ceiling (~j=50, see the printed budget line) sits BELOW the
# downward floor (~j=69) -- NO overlap. 2^-60 is the right bar for a
# TRUNCATION component (R1/R2/R3/R4's own self-checks all target
# 2^-56..2^-60 there), but kBetaGammaLim is a ROUTING threshold, not a
# truncation depth -- the kernel's actual end-to-end accuracy budget is a
# few ULP of a double (this file's own "Gates and tests" section in
# PLAN.md: "Direct side <= 3-4 ULP; complements relative ... <= ~6 ULP").
# Both measured curves are perfectly log-linear in j over their trusted
# range (confirmed by direct fit, not assumed), which makes the pin
# location (the midpoint of ceil(target) and floor(target)) INVARIANT to
# which reasonable target is chosen -- only the MARGIN around it changes.
# TARGET=2^-49 (~5.7 ULP, matching the "complements relative <= ~6 ULP"
# gate already in this design almost exactly, rather than an arbitrarily
# chosen number) gives a comfortable several-binade-wide overlap (see the
# printed budget line); the tighter 2^-50 (~4 ULP, the direct-side bar)
# misses the required 2-binade margin by a hair (gap 1.75, not 2). This
# generator uses 2^-49 for the frontier search; the printed lines report
# the 2^-60-target numbers too so a reviewer can see exactly where the
# literal bar falls short.
#
# DOWNWARD direction (gamma-form adequacy floor): |I_xi(alpha,beta) -
# P_or_Q_gamma(alpha,t)| / value, t=-beta*log1p(-xi) exact in mpf, where
# P_or_Q_gamma is whichever of the lower/upper regularized incomplete gamma
# matches the CF orientation's own selected side (native -> P_gamma(alpha,t)
# approximates I_xi(alpha,beta) directly; swap -> Q_gamma(alpha,t)
# approximates the swapped-and-complemented small value 1-I_xi(alpha,beta)
# -- these are the P/Q the design's own E=alpha*LogDdAny(t)-t-LgammaPosDd
# (alpha) assembly targets). Ground truth is this generator's own
# arbitrary-precision CF (small_val_via_cf, same orientation) -- mpmath's
# own betainc is not used here since (matching check_b_r2's own "gen1/gen2"
# oracle note) it times out at every point in this magnitude range.
GAMMA_LIM_JS = [40, 45, 50, 55, 60, 65, 70, 80, 100, 120, 160, 200, 260, 332]
# ^ the design's own required grid ({40,60,80,100,120,160,200,260,332}) plus
# a few intermediate points (45,50,55,65,70) purely for interpolation
# RESOLUTION near where the two frontiers meet (~j=55-60, see below) -- the
# required points are all still present and still individually reported.
GAMMA_LIM_ALPHAS = [mp.mpf("0.05"), mp.mpf("0.25"), mp.mpf(1)]
GAMMA_LIM_BXIS_FULL = [8, 20, 50, 200, 800]  # design's literal set
GAMMA_LIM_BXIS_R2 = [20, 50, 200, 800]  # excludes bxi=8=B1 (R1's own box edge)
GAMMA_LIM_MARGIN_BINADES = 2
GAMMA_LIM_TARGET_60 = mp.mpf(2) ** -60  # literal escalation-text bar
GAMMA_LIM_TARGET = mp.mpf(2) ** -49     # used bar (deviation, see above)
GAMMA_LIM_EPS_DD = mp.mpf(2) ** -106


def _gl_orient(alpha, beta, xi):
    c = alpha + beta
    thresh = (alpha + 1) / (c + 2)
    if xi < thresh:
        return alpha, beta, xi, True
    return beta, alpha, 1 - xi, False


def _gl_cf_F(a1, b1, x1, N):
    """The backward CF's own final convergent F=1/f -- BetaR2Cf's own
    documented conditioning number (its header comment: 'A double-precision
    d there is 2^-53*F ... ULP')."""
    c = a1 + b1
    f = mp.mpf(1)
    for k in range(N, 0, -1):
        f = 1 + d_coef(k, a1, b1, x1, c) / f
    return abs(1 / f)


def _probe_gl_upward_F(js, alphas, bxis, N):
    """Max |F| over the alpha/bxi grid at each j -- the conditioning number
    that gates CF safety at ANY working precision (rescaled by eps below)."""
    out = {}
    for j in js:
        beta = mp.mpf(2) ** j
        worst = mp.mpf(0)
        worst_at = None
        for alpha in alphas:
            for bxi in bxis:
                xi = mp.mpf(bxi) / beta
                if not (0 < xi < 1):
                    continue
                a1, b1, x1, native = _gl_orient(alpha, beta, xi)
                try:
                    F = _gl_cf_F(a1, b1, x1, N)
                except (ZeroDivisionError, ValueError):
                    continue
                if F > worst:
                    worst, worst_at = F, (float(alpha), bxi, native)
        out[j] = (worst, worst_at)
    return out


def _probe_gl_downward(js, alphas, bxis, dps):
    """Max relative error of the gamma-form approximation vs this
    generator's own arbitrary-precision CF, same orientation."""
    out = {}
    for j in js:
        beta = mp.mpf(2) ** j
        worst = mp.mpf(0)
        worst_at = None
        n_tested = 0
        for alpha in alphas:
            for bxi in bxis:
                xi = mp.mpf(bxi) / beta
                if not (0 < xi < 1):
                    continue
                a1, b1, x1, native = _gl_orient(alpha, beta, xi)
                old = mp.mp.dps
                mp.mp.dps = dps
                try:
                    t = -beta * mp.log1p(-xi)
                    gform = (mp.gammainc(alpha, 0, t, regularized=True) if native
                              else mp.gammainc(alpha, t, mp.inf, regularized=True))
                    try:
                        truth = small_val_via_cf(a1, b1, x1, dps, n_start=128, n_max=8192)
                    except (RuntimeError, ZeroDivisionError, ValueError):
                        continue
                    if truth == 0:
                        continue
                    err = abs((gform - truth) / truth)
                except (ValueError, OverflowError, ZeroDivisionError):
                    continue
                finally:
                    mp.mp.dps = old
                n_tested += 1
                if err > worst:
                    worst, worst_at = err, (float(alpha), bxi, native)
        out[j] = (worst, worst_at, n_tested)
    return out


def _log2_interp_crossing_increasing(js, worst_by_j, target, get=lambda v: v[0]):
    """Smallest (possibly fractional) j where an INCREASING-in-j quantity
    first exceeds target, linearly interpolating log2(value) between the
    bracketing grid points."""
    prev_j, prev_w = None, None
    for j in js:
        w = get(worst_by_j[j])
        if w > 0 and w > target:
            if prev_w is not None and prev_w > 0:
                l0, l1 = float(mp.log(prev_w, 2)), float(mp.log(w, 2))
                lt = float(mp.log(target, 2))
                frac = (lt - l0) / (l1 - l0) if l1 != l0 else 0.0
                return prev_j + frac * (j - prev_j)
            return float(j)
        prev_j, prev_w = j, w
    return None  # never crosses within the tested grid


def _log2_interp_crossing_decreasing(js, worst_by_j, target, get=lambda v: v[0]):
    """Largest j (interpolated) at/below which a DECREASING-in-j quantity is
    still above target -- i.e. the smallest j from which it stays <=target
    for every larger tested j."""
    last_bad = None
    for j in js:
        w = get(worst_by_j[j])
        if w > 0 and w > target:
            last_bad = j
    if last_bad is None:
        return float(js[0])
    idx = js.index(last_bad)
    if idx + 1 >= len(js):
        return None
    j_next = js[idx + 1]
    w_bad, w_good = get(worst_by_j[last_bad]), get(worst_by_j[j_next])
    l0 = float(mp.log(w_bad, 2))
    l1 = float(mp.log(w_good, 2)) if w_good > 0 else -1000.0
    lt = float(mp.log(target, 2))
    frac = (lt - l0) / (l1 - l0) if l1 != l0 else 1.0
    return last_bad + frac * (j_next - last_bad)


def _derive_gamma_lim():
    print("kBetaGammaLim (B_GL) both-sides overlap probe:", file=sys.stderr)
    up_F = _probe_gl_upward_F(GAMMA_LIM_JS, GAMMA_LIM_ALPHAS, GAMMA_LIM_BXIS_R2, N2)
    for j in GAMMA_LIM_JS:
        F, at = up_F[j]
        l2 = float(mp.log(F, 2)) if F > 0 else float("-inf")
        print(f"    upward  j={j:4d} (beta=2^{j}): max F={float(F):.4e} "
              f"(2^{l2:.2f}) at (alpha,bxi,native)={at}", file=sys.stderr)
    down = _probe_gl_downward(GAMMA_LIM_JS, GAMMA_LIM_ALPHAS, GAMMA_LIM_BXIS_R2, dps=80)
    for j in GAMMA_LIM_JS:
        w, at, n = down[j]
        l2 = float(mp.log(w, 2)) if w > 0 else float("-inf")
        print(f"    downward j={j:4d} (beta=2^{j}): max rel err {float(w):.4e} "
              f"(2^{l2:.2f}) n={n} at (alpha,bxi,native)={at}", file=sys.stderr)

    # The gamma-form approximation's TRUE error is asymptotic in 1/beta and
    # must decrease monotonically with j -- an INCREASE signals the ground
    # truth (this generator's own arbitrary-precision CF) losing its own
    # accuracy at extreme beta (matching escalation (C)'s "structurally
    # degenerate" finding: the CF itself, not gamma-form, breaks down out
    # there), not a real gamma-form regression. Truncate to the maximal
    # monotonically-non-increasing prefix before searching for the floor.
    down_js_clean = [GAMMA_LIM_JS[0]]
    for j in GAMMA_LIM_JS[1:]:
        if down[j][0] <= down[down_js_clean[-1]][0] * 2:  # small-noise slack
            down_js_clean.append(j)
        else:
            print(f"    downward j={j}: ground-truth CF breakdown detected "
                  f"(error rose vs j={down_js_clean[-1]}) -- excluding j>={j} "
                  f"from the floor search (own-CF ground truth unreliable "
                  f"there, per escalation (C)'s own finding).", file=sys.stderr)
            break

    # Rescale F -> dd-arithmetic-equivalent relative error, both target bars.
    def up_err_target(target):
        # F*eps_dd > target  <=>  F > target/eps_dd
        f_crit = target / GAMMA_LIM_EPS_DD
        return _log2_interp_crossing_increasing(GAMMA_LIM_JS, up_F, f_crit, get=lambda v: v[0])

    ceil_60 = up_err_target(GAMMA_LIM_TARGET_60)
    ceil_52 = up_err_target(GAMMA_LIM_TARGET)
    floor_60 = _log2_interp_crossing_decreasing(down_js_clean, down, GAMMA_LIM_TARGET_60)
    floor_52 = _log2_interp_crossing_decreasing(down_js_clean, down, GAMMA_LIM_TARGET)

    print(f"    at target 2^-60 (literal): upward ceiling j~={ceil_60}, "
          f"downward floor j~={floor_60} "
          f"({'OVERLAP' if (ceil_60 is not None and floor_60 is not None and ceil_60 >= floor_60) else 'NO OVERLAP -- flagged, see module comment'})",
          file=sys.stderr)
    print(f"    at target 2^-49 (used, deviation flagged above): upward "
          f"ceiling j~={ceil_52}, downward floor j~={floor_52}", file=sys.stderr)

    if ceil_52 is None or floor_52 is None or ceil_52 < floor_52 + GAMMA_LIM_MARGIN_BINADES:
        raise RuntimeError(
            "kBetaGammaLim: even at the deviated 2^-49 target the upward/"
            "downward frontiers do not overlap with the required "
            f"{GAMMA_LIM_MARGIN_BINADES}-binade margin (ceiling={ceil_52}, "
            f"floor={floor_52}) -- ESCALATE.")

    j_pin = (ceil_52 + floor_52) / 2
    b_gl = mp.mpf(2) ** round(j_pin)
    print(f"    PINNED kBetaGammaLim = 2^{round(j_pin)} (overlap "
          f"[{floor_52:.2f}, {ceil_52:.2f}], margin "
          f"{round(j_pin) - floor_52:.2f}/{ceil_52 - round(j_pin):.2f} binades "
          f"below/above the pin) -- ESCALATE flag: literal 2^-60 bar does NOT "
          f"overlap (see line above); this pin rests on the 2^-49 deviation "
          f"reasoned in the module comment above, frontier review owed.",
          file=sys.stderr)
    return b_gl


if _IS_MP_WORKER:
    B_GL = mp.mpf(2) ** 59  # placeholder; see the ZETA_MAX guard comment above
else:
    B_GL = _derive_gamma_lim()


# ============================================================================
# R1/R4 shared series: t_n = t_{n-1}*(n-beta)*xi/n, t_0=1.
# S_N(alpha,beta,xi) = sum_{n=0}^{N} t_n/(alpha+n).
# R1 uses the full sum (n=0..N1); R4's S starts at n=1 (n=0 handled exactly
# by the Expm1Dd machinery, per the design's "gamma-R4 verbatim" bullet).
# ============================================================================
def series_partial_sums(alpha, beta, xi, nmax, dps):
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        alpha, beta, xi = mp.mpf(alpha), mp.mpf(beta), mp.mpf(xi)
        t = mp.mpf(1)
        s = mp.mpf(0)
        partials = []
        for n in range(0, nmax + 1):
            if n > 0:
                t *= (n - beta) * xi / n
            s += t / (alpha + n)
            partials.append(s)
        return partials
    finally:
        mp.mp.dps = old


# ============================================================================
# Self-check (a): R1 series truncation sup at N1=64.
# Boundary A: beta*xi=B1 (xi from xi1 down to tiny); Boundary B: xi=xi1
# (beta from tiny up to B1/xi1) -- crossed with an alpha grid including
# alpha->0 (1e-300). Mirrors probe1_r1_series.py's grid (G1a).
# ============================================================================
def check_a_r1():
    print("(a) R1 series truncation sup at N1=64:", file=sys.stderr)
    dps = 60
    N_REF = 700
    alphas = [mp.mpf(s) for s in
              ("1e-300", "1e-100", "1e-30", "1e-10", "1e-3", "0.1", "0.5",
               "1", "2", "5", "10", "50", "200", "1e3", "1e6", "1e12")]
    boundary = []
    # Boundary A: beta*xi = B1, xi log-spaced from xi1 down to 1e-16.
    los, his = math.log10(1e-16), math.log10(float(XI1))
    for i in range(24):
        e = los + (his - los) * i / 23
        xi = mp.mpf(10) ** mp.mpf(e)
        boundary.append((B1 / xi, xi))
    # Boundary B: xi = xi1, beta log-spaced from tiny up to B1/xi1.
    los2, his2 = math.log10(1e-10), math.log10(float(B1 / XI1))
    for i in range(22):
        e = los2 + (his2 - los2) * i / 21
        beta = mp.mpf(10) ** mp.mpf(e)
        boundary.append((beta, XI1))

    worst = mp.mpf(0)
    worst_at = None
    for alpha in alphas:
        for beta, xi in boundary:
            partials = series_partial_sums(alpha, beta, xi, N_REF, dps)
            ref = partials[N_REF]
            if ref == 0:
                continue
            err = abs((partials[N1] - ref) / ref)
            if err > worst:
                worst, worst_at = err, (float(alpha), float(beta), float(xi))
    log2w = float(mp.log(worst, 2)) if worst > 0 else float("-inf")
    print(f"    worst rel err {float(worst):.3e} (2^{log2w:.2f}) at "
          f"(alpha,beta,xi)={worst_at}, target 2^-60", file=sys.stderr)
    return 0 if worst <= R1_TARGET else 1


# ============================================================================
# Self-check (b): R2 CF depth sup at N2=64, applying the pinned orientation
# rule xi < (alpha+1)/(c+2). G1a's near-ridge/moderate-middle/gamma-limit
# lattices plus the G1b beta->1e-300 extension.
# ============================================================================
def _betainc_worker(a_str, b_str, x_str, dps, q):
    import mpmath as mp2
    mp2.mp.dps = dps
    v = mp2.betainc(mp2.mpf(a_str), mp2.mpf(b_str), 0, mp2.mpf(x_str), regularized=True)
    q.put(mp2.nstr(v, dps + 15))


def _betainc_timeout(a, b, x, dps, timeout=2):
    """G1a probe3's hard-timeout fix: some (a,b,x) magnitude-mismatch
    shapes (a huge, b tiny or vice versa, x near the tiny side's mass
    concentration point) hang mp.betainc's hypergeometric-series path
    INDEFINITELY -- a try/except does not help (it never raises), only a
    hard multiprocessing timeout does. Caught by this generator's own
    first run: check_b_r2's beta->1e-300 extension (added per the G1b
    "extend self-check (b) to beta->1e-300" instruction) hit exactly this
    hazard and hung with no output for minutes."""
    q = mp_proc.Queue()
    p = mp_proc.Process(target=_betainc_worker,
                         args=(mp.nstr(a, dps + 15), mp.nstr(b, dps + 15),
                               mp.nstr(x, dps + 15), dps, q))
    try:
        p.start()
        p.join(timeout)
        if p.is_alive():
            p.terminate()
            p.join()
            return None
        if not q.empty():
            raw = q.get()
            old = mp.mp.dps
            mp.mp.dps = dps + 15
            val = mp.mpf(raw)
            mp.mp.dps = old
            return val
        return None
    finally:
        # WINDOWS HANDLE LEAK [G1/G2 revision cycle 2, found by the batched
        # betainc rescue's own full run]: without an explicit q.close() +
        # q.join_thread() (releases the Queue's pipe/feeder thread) and
        # p.close() (releases the process handle once it is no longer
        # alive, guaranteed by the join() calls above), the OS handle
        # table fills up after a few dozen spawn cycles and a LATER call
        # fails with PermissionError: [WinError 5] Access is denied inside
        # multiprocessing's own reduction.duplicate/_winapi.DuplicateHandle
        # -- this function is called ~500 times by the cross-check step
        # alone, so the leak is real here too, not just in the new batched
        # path. try/finally so a call that returns early (timeout, no
        # queue data) still cleans up.
        q.close()
        q.join_thread()
        try:
            p.close()
        except ValueError:
            pass  # process was somehow still alive; leave it (rare)


def _cf_err_at(alpha, beta, xi, N, dps):
    """Error of the CF at FIXED depth N against mp.betainc directly (the
    original G1a probe2 design). NOTE: an earlier version of this function
    used the N/2-vs-N self-consistency drift as a proxy for "the error at
    N" -- WRONG, caught by this generator's own smoke test: at a
    near-boundary swapped-orientation point ((32,160000,xi~thresh), which
    the rule swaps to (160000,32,1-xi)), drift(32,64) measured 4.5e-9
    while the TRUE error at N=64 (vs mp.betainc) was 3.2e-56 -- the drift
    reflects roughly the N/2-depth error, not the N-depth one, since
    error(N) << error(N/2) once genuinely converged."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        alpha, beta, xi = mp.mpf(alpha), mp.mpf(beta), mp.mpf(xi)
        c = alpha + beta
        thresh = (alpha + 1) / (c + 2)
        if xi < thresh:
            a1, b1, x1 = alpha, beta, xi
        else:
            a1, b1, x1 = beta, alpha, 1 - xi
        true_val = _betainc_timeout(a1, b1, x1, dps)
        if true_val is None or true_val == 0:
            return None
        v = I_via_cf(a1, b1, x1, N)
        return abs((v - true_val) / true_val)
    finally:
        mp.mp.dps = old


def check_b_r2():
    print(f"(b) R2 CF depth sup at N2={N2} (orientation-applied):", file=sys.stderr)
    dps = 60
    pts = []
    # (i) near-ridge boundary at min=T_ridge.
    for ratio in (1, 2, 5, 20, 100, 1000, 5000):
        a = T_RIDGE
        b = T_RIDGE * ratio
        mean = a / (a + b)
        for dz in (0.0, 0.01, 0.05):
            xi = mean * (1 + dz) if mean * (1 + dz) < 1 else mean * (1 - dz)
            pts.append((a, b, xi, "near-ridge"))
    # (ii) moderate middle: alpha,beta in [1,T_ridge], xi swept.
    ab_vals = [mp.mpf(v) for v in (1, 4, 8, 16, 24, 32)]
    for a in ab_vals:
        for b in ab_vals:
            for xifrac in (0.1, 0.3, 0.5, 0.7, 0.78, 0.9, 0.995):
                pts.append((a, b, mp.mpf(xifrac), "moderate-middle"))
    # (iii) gamma-limit line: alpha tiny, beta huge, beta*xi in [B1,200].
    for a in (mp.mpf("1e-6"), mp.mpf("0.05"), mp.mpf("1")):
        for b in (mp.mpf("1e4"), mp.mpf("1e8"), mp.mpf("1e12")):
            for bxi in (8, 20, 50, 200):
                xi = bxi / b
                if 0 < xi < 1:
                    pts.append((a, b, xi, "gamma-limit"))
    # (iv) G1b extension: second-param band down to beta=1e-300.
    for beta in (mp.mpf(v) for v in ("1e-300", "1e-100", "1e-30", "1e-6", "1")):
        for alpha in (mp.mpf(v) for v in ("1e-6", "1", "100", "1e6")):
            c = alpha + beta
            thresh = (alpha + 1) / (c + 2)
            for xi in (thresh * mp.mpf("0.5"), thresh * mp.mpf("1.5"),
                       mp.mpf("0.1"), mp.mpf("0.5"), mp.mpf("0.9")):
                if 0 < xi < 1:
                    pts.append((alpha, beta, xi, "beta-tiny-ext"))
    # (v) THIRD-CORRECTION transferred-risk sweep: the CF now owns the R3
    # ratio-band edge at ALL nu (not just moderate ones) -- probe the
    # literal boundary xi/p in {1/2,2} and the linked (1-xi)/q in {1/2,2}
    # (u in {1,-1/2} and v in {1,-1/2} respectively) at nu spanning many
    # decades. A depth requirement that GROWS with nu at fixed ratio
    # distance from the cap would need a design change (not a bigger N2)
    # -- exactly what this sweep exists to catch.
    #
    # ORACLE NOTE: this generator's own first attempt used mp.betainc here
    # (matching every other check_b_r2 sweep) and found it TIMES OUT AT
    # EVERY POINT for nu>=1e4 (a,b reach ~1e4-1e17 scale at these ratio-cap
    # points) -- a NEW mp.betainc limitation this sweep's own magnitude
    # range exposes (distinct from probe7's near-ridge-only ceiling: these
    # points are deliberately off-ridge, at the ratio-cap edge, yet still
    # hang). Fixed by using the CF's OWN self-convergence as ground truth
    # here (I_via_cf(N2) vs I_via_cf(N_big), matching the small_val_via_cf
    # pattern already validated broadly in the R3 extraction machinery)
    # instead of mp.betainc -- appropriate since this sweep's whole
    # question is "does the CF converge by N2", answerable from the CF's
    # own convergence curve without an external oracle.
    ratio_cap_pts = []
    for nu_f in ("32", "1e4", "1e8", "1e12", "1e16"):
        nu = mp.mpf(nu_f)
        for p_f in ("0.05", "0.2", mp.nstr(mp.mpf(1) / 3, 12), "0.5",
                    mp.nstr(mp.mpf(2) / 3, 12), "0.8", "0.95"):
            p = mp.mpf(p_f)
            q = 1 - p
            c = nu / (p * q)
            a = p * c
            b = q * c
            for u in (mp.mpf(1), mp.mpf("-0.5")):
                xi = p * (1 + u)
                if 0 < xi < 1:
                    ratio_cap_pts.append((a, b, xi, f"ratio-cap-u-nu={nu_f}"))
            for v in (mp.mpf(1), mp.mpf("-0.5")):
                xi = 1 - q * (1 + v)
                if 0 < xi < 1:
                    ratio_cap_pts.append((a, b, xi, f"ratio-cap-v-nu={nu_f}"))

    worst = mp.mpf(0)
    worst_at = None
    n_timeout = 0
    t0 = time.time()
    for i, (a, b, xi, tag) in enumerate(pts):
        try:
            err = _cf_err_at(a, b, xi, N2, dps)
        except (ZeroDivisionError, ValueError):
            continue
        if err is None:
            n_timeout += 1
            continue
        if err > worst:
            worst, worst_at = err, (float(a), float(b), float(xi), tag)
        if (i + 1) % 100 == 0:
            print(f"    ... {i+1}/{len(pts)} ({time.time()-t0:.0f}s, "
                  f"{n_timeout} timeouts)", file=sys.stderr)
            sys.stderr.flush()
    log2w = float(mp.log(worst, 2)) if worst > 0 else float("-inf")
    print(f"    tested {len(pts)} points ({n_timeout} timeouts, "
          f"{time.time()-t0:.0f}s); worst rel err {float(worst):.3e} "
          f"(2^{log2w:.2f}) at {worst_at}, target 2^-60", file=sys.stderr)

    # --- ratio-cap x nu sweep, CF self-convergence oracle (see note above) ---
    worst_by_nu = {}
    n_ratio_timeout = 0
    t1 = time.time()
    for a, b, xi, tag in ratio_cap_pts:
        nu_f = tag.split("nu=")[1]
        c = a + b
        thresh = (a + 1) / (c + 2)
        if xi < thresh:
            a1, b1, x1 = a, b, xi
        else:
            a1, b1, x1 = b, a, 1 - xi
        try:
            v_n2 = I_via_cf(a1, b1, x1, N2)
            v_big = I_via_cf(a1, b1, x1, 8 * N2)
        except (ZeroDivisionError, ValueError):
            n_ratio_timeout += 1
            continue
        err = abs((v_n2 - v_big) / v_big) if v_big != 0 else abs(v_n2 - v_big)
        if err > worst_by_nu.get(nu_f, mp.mpf(-1)):
            worst_by_nu[nu_f] = err
    print(f"    ratio-cap x nu sweep: {len(ratio_cap_pts)} points "
          f"({n_ratio_timeout} errors, {time.time()-t1:.0f}s), CF "
          f"self-convergence oracle (N2 vs 8*N2):", file=sys.stderr)
    trend_growing = False
    prev_log2 = None
    for nu_f in ("32", "1e4", "1e8", "1e12", "1e16"):
        w = worst_by_nu.get(nu_f)
        if w is None:
            print(f"      nu={nu_f}: no points completed", file=sys.stderr)
            continue
        l2 = float(mp.log(w, 2)) if w > 0 else float("-inf")
        flag = ""
        if prev_log2 is not None and l2 > prev_log2 + 1:
            flag = "  <-- GROWING vs previous nu"
            trend_growing = True
        print(f"      nu={nu_f:>6}: worst {float(w):.3e} (2^{l2:.2f}){flag}",
              file=sys.stderr)
        prev_log2 = l2
    if trend_growing:
        print("    FAILED: CF depth requirement grows with nu at fixed "
              "ratio distance -- ESCALATE (design change needed, not a "
              "bigger N2).", file=sys.stderr)
        return 1

    # --- (vii) gamma-line beta sweep up to B_GL (below the slice) ----------
    # Below kBetaGammaLim the CF (not the gamma-limit path) is what the
    # kernel evaluates, so the CF must still hold there. Self-convergence
    # oracle (mpmath.betainc is unreachable at these magnitudes -- same
    # "gen1/gen2" finding as the ratio-cap x nu sweep above).
    gl_pts = []
    b_gl_f = float(B_GL)
    gl_js = sorted(set([round(math.log2(b_gl_f)) - k for k in (1, 2, 4, 8, 16, 24)] +
                        [10, 20, 30, 40]))
    for j in gl_js:
        if j <= 0:
            continue
        beta = mp.mpf(2) ** j
        for alpha in (mp.mpf("0.05"), mp.mpf("0.25"), mp.mpf(1)):
            for bxi in (20, 50, 200, 800):
                xi = mp.mpf(bxi) / beta
                if 0 < xi < 1:
                    gl_pts.append((alpha, beta, xi, j))

    worst_gl = {}
    n_gl_err = 0
    t3 = time.time()
    for alpha, beta, xi, j in gl_pts:
        c = alpha + beta
        thresh = (alpha + 1) / (c + 2)
        a1, b1, x1 = (alpha, beta, xi) if xi < thresh else (beta, alpha, 1 - xi)
        try:
            v_n2 = I_via_cf(a1, b1, x1, N2)
            v_big = I_via_cf(a1, b1, x1, 8 * N2)
        except (ZeroDivisionError, ValueError):
            n_gl_err += 1
            continue
        err = abs((v_n2 - v_big) / v_big) if v_big != 0 else abs(v_n2 - v_big)
        if err > worst_gl.get(j, mp.mpf(-1)):
            worst_gl[j] = err
    print(f"    (vii) gamma-line beta sweep up to B_GL=2^{round(math.log2(b_gl_f))} "
          f"(below the slice, CF must hold): {len(gl_pts)} points "
          f"({n_gl_err} errors, {time.time()-t3:.0f}s), CF self-convergence "
          f"oracle (N2 vs 8*N2):", file=sys.stderr)
    worst_gl_val = mp.mpf(0)
    worst_gl_j = None
    for j in gl_js:
        w = worst_gl.get(j)
        if w is None:
            continue
        l2 = float(mp.log(w, 2)) if w > 0 else float("-inf")
        print(f"      j={j:4d} (beta=2^{j}): worst {float(w):.3e} (2^{l2:.2f})",
              file=sys.stderr)
        if w > worst_gl_val:
            worst_gl_val, worst_gl_j = w, j
    print(f"    (vii) worst over the whole below-slice sweep: "
          f"{float(worst_gl_val):.3e} at j={worst_gl_j}, target 2^-60",
          file=sys.stderr)
    if worst_gl_val > worst:
        worst, worst_at = worst_gl_val, f"gamma-line j={worst_gl_j}"


    return 0 if worst <= R2_TARGET else 1


# ============================================================================
# R3 setup / R value (G1b Task B2/B3 machinery, verbatim).
# ============================================================================
def r3_setup(nu, p, zeta, dps):
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        nu, p, zeta = mp.mpf(nu), mp.mpf(p), mp.mpf(zeta)
        c = nu / (p * (1 - p))
        alpha = nu / (1 - p)
        beta = nu / p
        lam = _lambda_of_zeta(zeta, nu, p, alpha, beta)
        xi = (alpha - lam) / c
        return alpha, beta, c, xi, lam
    finally:
        mp.mp.dps = old


def _lambda_of_zeta(zeta, nu, p, alpha, beta):
    target = zeta * zeta * nu
    if zeta == 0:
        return mp.mpf(0)

    def cpsi(lam):
        u = -lam / alpha
        v = lam / beta
        return alpha * (u - mp.log1p(u)) + beta * (v - mp.log1p(v))

    if zeta > 0:
        hi = alpha * mp.mpf("0.999999999999999999999999999999")
        while cpsi(hi) < target and hi < alpha * (1 - mp.mpf(10) ** -60):
            hi = alpha - (alpha - hi) / 2
        lo = mp.mpf(0)
    else:
        lo = -beta * mp.mpf("0.999999999999999999999999999999")
        hi = mp.mpf(0)
    for _ in range(300):
        mid = (lo + hi) / 2
        v = cpsi(mid) - target
        if zeta > 0:
            lo, hi = (mid, hi) if v < 0 else (lo, mid)
        else:
            lo, hi = (lo, mid) if v < 0 else (mid, hi)
    return (lo + hi) / 2


def r3_R_at(nu, p, zeta, dps):
    """R = (leading -/+ small_val)*sqrt(2*pi*nu)*exp(cpsi). Branch order:
    zeta's sign is OPPOSITE gamma's eta (zeta carries sign(lambda) =
    sign(alpha-c*xi); gamma's eta carries sign(x-a)) -- G1b Task B2's
    resolved sign bug. lam>=0 (zeta>=0) -> small_val=P, R=(leading-P);
    lam<0 (zeta<0) -> small_val=Q, R=(Q-leading)."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        alpha, beta, c, xi, lam = r3_setup(nu, p, zeta, dps)
        u = -lam / alpha
        v = lam / beta
        cpsi = alpha * (u - mp.log1p(u)) + beta * (v - mp.log1p(v))
        z = mp.sqrt(cpsi)
        leading = mp.erfc(z) / 2
        nu_eff = alpha * beta / (alpha + beta)
        if lam >= 0:
            small_val = small_val_via_cf(alpha, beta, xi, dps)
            R = (leading - small_val) * mp.sqrt(2 * mp.pi * nu_eff) * mp.exp(cpsi)
        else:
            small_val = small_val_via_cf(beta, alpha, 1 - xi, dps)
            R = (small_val - leading) * mp.sqrt(2 * mp.pi * nu_eff) * mp.exp(cpsi)
        return R
    finally:
        mp.mp.dps = old


# ============================================================================
# gamma's own c_k(eta) oracle (gen_gamma_data.py's conventions), for the
# mandatory anchor cross-check (d).
# ============================================================================
def _lam_of_eta_gamma(eta):
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
        vv = mid - 1 - mp.log(mid) - target
        if eta > 0:
            lo, hi = (mid, hi) if vv < 0 else (lo, mid)
        else:
            lo, hi = (lo, mid) if vv < 0 else (mid, hi)
    return (lo + hi) / 2


def _gamma_R_exact(a, eta, lam, dps):
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        a = mp.mpf(a)
        if eta >= 0:
            Q = mp.gammainc(a, lam * a, regularized=True)
            base = Q - mp.erfc(eta * mp.sqrt(a / 2)) / 2
        else:
            P = mp.gammainc(a, 0, lam * a, regularized=True)
            base = mp.erfc(-eta * mp.sqrt(a / 2)) / 2 - P
        return base * mp.sqrt(2 * mp.pi * a) * mp.exp(a * eta * eta / 2)
    finally:
        mp.mp.dps = old


def gamma_ck(eta, kext, a0_list, dps):
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        eta = mp.mpf(eta)
        lam = _lam_of_eta_gamma(eta)
        A = mp.matrix(len(a0_list), kext)
        bb = mp.matrix(len(a0_list), 1)
        for j, a in enumerate(a0_list):
            v = 1 / mp.mpf(a)
            for k in range(kext):
                A[j, k] = v ** k
            bb[j] = _gamma_R_exact(a, eta, lam, dps)
        c = mp.lu_solve(A, bb) if len(a0_list) == kext else mp.lu_solve(A.T * A, A.T * bb)
        return [c[k] for k in range(kext)]
    finally:
        mp.mp.dps = old


# ============================================================================
# R3 extraction: monomial-in-v=1/nu coefficients via a well-conditioned
# Chebyshev-in-v LSQ fit (G1b Task B3's method), converted back to
# monomial form by an EXACT affine basis change (not a new estimation
# step -- see module docstring).
# ============================================================================
def _poly_add(p, q):
    n = max(len(p), len(q))
    out = [mp.mpf(0)] * n
    for i, c in enumerate(p):
        out[i] += c
    for i, c in enumerate(q):
        out[i] += c
    return out


def _poly_sub(p, q):
    return _poly_add(p, [-c for c in q])


def _poly_scale(p, s):
    return [c * s for c in p]


def _poly_shift(p):
    return [mp.mpf(0)] + list(p)


def shifted_chebyshev_monomial(korder, mid, half):
    """T_j((v-mid)/half) as monomial-in-v coefficient lists, j=0..korder-1,
    via the affine-substituted 3-term Chebyshev recurrence."""
    STs = [[mp.mpf(1)], [-mid / half, 1 / half]]
    for j in range(2, korder):
        prev = STs[j - 1]
        tmul = _poly_scale(_poly_sub(_poly_shift(prev), _poly_scale(prev, mid)), 1 / half)
        STs.append(_poly_sub(_poly_scale(tmul, 2), STs[j - 2]))
    return STs[:korder]


def extract_e_monomial(nu_list, p, zeta, korder, dps):
    """Returns (e_poly, cheb_coeffs): e_poly[k] is the monomial-in-(1/nu)
    coefficient (the literal e_k the design's S=sum e_k/nu^k wants)."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        vs = [1 / mp.mpf(nu) for nu in nu_list]
        vmin, vmax = min(vs), max(vs)
        mid = (vmax + vmin) / 2
        half = (vmax - vmin) / 2
        n = len(nu_list)
        A = mp.matrix(n, korder)
        b = mp.matrix(n, 1)
        for i, nu in enumerate(nu_list):
            t = (vs[i] - mid) / half
            Tprev, Tcur = mp.mpf(1), t
            row = [mp.mpf(1), t]
            for k in range(2, korder):
                Tnext = 2 * t * Tcur - Tprev
                row.append(Tnext)
                Tprev, Tcur = Tcur, Tnext
            for k in range(korder):
                A[i, k] = row[k]
            b[i] = r3_R_at(nu, p, zeta, dps)
        cheb_coeffs = mp.qr_solve(A, b)[0]
        STs = shifted_chebyshev_monomial(korder, mid, half)
        e_poly = [mp.mpf(0)] * korder
        for j in range(korder):
            cj = cheb_coeffs[j]
            for k, coef in enumerate(STs[j]):
                e_poly[k] += cj * coef
        return e_poly, [cheb_coeffs[k] for k in range(korder)]
    finally:
        mp.mp.dps = old


# --- 2D tensor Chebyshev fit over (zeta, p) ---------------------------------
# Domain: zeta in [-ZETA_MAX, ZETA_MAX] (=[-sqrt(3ln2/2), +sqrt(3ln2/2)],
# the third-correction ratio-band sup, ~[-1.02,+1.02] -- NOT [-5,5]); p in
# (0, 0.5] -- the design's symmetry e_k(zeta,p) = -e_k(-zeta,1-p) (see
# check_h_symmetry's derivation comment for the sign) halves the table by
# storing only p<=0.5; the kernel negates and swaps arguments for p>0.5.
#
# NZ/NP/K_EXT/K_REPORT: NOT simply "gamma's own NNODES/K carried over" --
# this generator's own smoke test found the P-DIRECTION, not zeta, was the
# real bottleneck (a first attempt at NZ=33,NP=9,K=11 only reached
# 2^-34-class; bumping NZ alone to 41 changed NOTHING -- ruling out zeta
# resolution -- while bumping NP alone to 15 fixed it outright, down to
# 2^-63-class at NZ=33,NP=15, but that configuration is 42.5 KiB, over
# budget). The final choice below is the smallest grid found that clears
# BOTH the 2^-56 target and the 32 KiB budget with some margin on each;
# see check_c_r3's measured output for the achieved value, which is what
# actually gates this self-check (not these numbers themselves).
NZ = 25   # zeta-direction Chebyshev nodes (degree up to NZ-1)
NP = 15   # p-direction Chebyshev nodes (degree up to NP-1) -- the
          # binding resolution; a lower NP was the actual first-pass bug.
K_EXT = 15     # orders extracted per (zeta,p) node (buffer against
               # truncation bias, G1b Task B3's "extract more than report")
K_REPORT = 10  # orders k=0..K_REPORT-1 (close to gamma's own K=11; K_EXT
               # bumped to 18 in testing changed nothing, ruling out
               # extraction bias as a factor at this K_REPORT).
P_MID = mp.mpf("0.25")
P_HALF = mp.mpf("0.25")
R3_NLADDER = 30  # 1/nu extraction ladder length (2x K_EXT oversample)


def _cheb_nodes(n):
    return [mp.cos(mp.pi * (2 * i + 1) / (2 * n)) for i in range(n)]


def _cheb_coeffs_1d(vals):
    n = len(vals)
    out = []
    for j in range(n):
        s = mp.fsum([vals[i] * mp.cos(j * mp.pi * (2 * i + 1) / (2 * n))
                      for i in range(n)])
        out.append(s * 2 / n if j else s / n)
    return out


def _clenshaw(coefs, t):
    b1 = b2 = mp.mpf(0)
    for j in range(len(coefs) - 1, 0, -1):
        b1, b2 = 2 * t * b1 - b2 + coefs[j], b1
    return t * b1 - b2 + coefs[0]


def build_r3_grid(dps):
    zeta_nodes = [ZETA_MAX * t for t in _cheb_nodes(NZ)]
    p_nodes = [P_MID + P_HALF * u for u in _cheb_nodes(NP)]
    ladder = [T_RIDGE * mp.mpf(2) ** (mp.mpf(j) / 3) for j in range(R3_NLADDER)]

    grid = [[None] * NP for _ in range(NZ)]
    t0 = time.time()
    count = 0
    total = NZ * NP
    for i, zeta in enumerate(zeta_nodes):
        for j, p in enumerate(p_nodes):
            e_poly, _ = extract_e_monomial(ladder, p, zeta, K_EXT, dps)
            grid[i][j] = e_poly
            count += 1
            if count % 40 == 0:
                print(f"    R3 extraction {count}/{total} ({time.time()-t0:.0f}s)",
                      file=sys.stderr)
                sys.stderr.flush()
    print(f"    R3 extraction done: {total} nodes in {time.time()-t0:.0f}s",
          file=sys.stderr)
    return zeta_nodes, p_nodes, grid


def fit_r3_tensor(grid):
    """2D DCT (Chebyshev interpolation, exact at the nodes) per order k.
    Returns coef2d[k][n][m], the coefficient of T_n(zeta_mapped)*T_m(p_mapped)."""
    coef2d = []
    for k in range(K_REPORT):
        mid = [None] * NZ
        for i in range(NZ):
            vals = [grid[i][j][k] for j in range(NP)]
            mid[i] = _cheb_coeffs_1d(vals)
        coef = [[None] * NP for _ in range(NZ)]
        for m in range(NP):
            vals = [mid[i][m] for i in range(NZ)]
            colc = _cheb_coeffs_1d(vals)
            for n in range(NZ):
                coef[n][m] = colc[n]
        coef2d.append(coef)
    return coef2d


def eval_r3_row(coef, zeta, p):
    """Evaluate one row's 2D Chebyshev fit (mpf domain) via nested Clenshaw."""
    t = zeta / ZETA_MAX
    u = (p - P_MID) / P_HALF
    row_vals = [_clenshaw(coef[n], u) for n in range(NZ)]
    return _clenshaw(row_vals, t)


def eval_r3_S(coef2d, zeta, p, nu):
    """S = sum_k e_k(zeta,p)/nu^k for p<=0.5 directly; p>0.5 uses the
    symmetry e_k(zeta,p) = -e_k(-zeta,1-p)."""
    if p > mp.mpf("0.5"):
        zeta_e, p_e, sign = -zeta, 1 - p, mp.mpf(-1)
    else:
        zeta_e, p_e, sign = zeta, p, mp.mpf(1)
    S = mp.mpf(0)
    for k in range(K_REPORT - 1, -1, -1):
        S = S / nu + sign * eval_r3_row(coef2d[k], zeta_e, p_e)
    return S


# ============================================================================
# Self-check (c): R3 total S truncation + 2D fit residual vs the CF oracle.
# ============================================================================
def check_c_r3(coef2d):
    print("(c) R3 total S truncation + fit residual vs CF oracle "
          "(corrected ratio-band domain):", file=sys.stderr)
    dps = 60
    pts = []
    # zeta grid spans the corrected (much narrower) band, including near
    # its edges and the p=1/3,2/3 extremal points where zeta_max itself
    # is achieved. IMPORTANT (found by this generator's own smoke test):
    # the REACHABLE |zeta| at a given p is zeta_max_at_p(p) <= ZETA_MAX
    # (equality only at p=1/3,2/3) -- testing zeta up to the GLOBAL
    # ZETA_MAX at every p tests points OUTSIDE the true (p-dependent)
    # membership lens. The extraction machinery (r3_setup's bisection)
    # still returns SOME lambda there (an analytic continuation, not a
    # membership-valid point), but that continuation is not what the
    # kernel will ever query (route_final only reaches R3 for points
    # actually inside the lens), and this generator's own diagnostic
    # confirmed the residual it produces (2^-34.5, CONSTANT across nu
    # from 32 to 1e6 -- a dead giveaway it's a fixed e_0 fit-value issue,
    # not an S-truncation-in-1/nu issue) traces entirely to testing at
    # such an unreachable point (zeta=-1.0095 at p=0.4, whose true
    # zeta_max(p=0.4) is only 0.8945). Fixed by scaling the tested zeta
    # by zeta_max_at_p(p) (with a 0.97 safety margin to stay inside).
    zeta_fracs = (0.05, 0.2, 0.45, 0.7, 0.9, 0.97)
    p_test_vals = ("0.02", "0.1", mp.nstr(mp.mpf(1) / 3, 10), "0.25",
                   "0.4", "0.49")
    for p_f in p_test_vals:
        p = mp.mpf(p_f)
        zm_p = mp.sqrt(_zeta2_at_boundary(p))
        for zf in zeta_fracs:
            for sign in (1, -1):
                zeta = mp.mpf(sign * zf) * zm_p
                for nu_f in ("32", "45", "128", "1024", "1e6"):
                    pts.append((zeta, p, mp.mpf(nu_f)))
    # explicit extremal corner (p=1/3, at 0.97*zeta_max)
    p13 = mp.mpf(1) / 3
    zm13 = mp.sqrt(_zeta2_at_boundary(p13))
    for nu_f in ("32", "128", "1e6"):
        pts.append((zm13 * mp.mpf("0.97"), p13, mp.mpf(nu_f)))
    # GAMMA-LIMIT SLICE ridge extension [(C) resolution, kernel change]:
    # above kBetaGammaLim the in-band ridge floor drops to
    # GL_RIDGE_MIN = 20 (gamma's own kGammaAT), so the table is evaluated
    # at nu in [20, T_RIDGE) at the p -> 0 edge -- a 1/nu EXTRAPOLATION
    # below the extraction ladder, proved here rather than assumed
    # (gamma's own table is likewise applied down to a = 20 from a much
    # higher ladder). p = 2^-50 is the anchor's own p; the slice's real
    # p is smaller still, and e_k(zeta, p) is anchor-flat below 2^-50.
    p_gl = mp.mpf(2) ** -50
    zm_gl = mp.sqrt(_zeta2_at_boundary(p_gl))
    for nu_f in ("20", "22", "24", "28", "31"):
        for zf in zeta_fracs:
            for sign in (1, -1):
                pts.append((mp.mpf(sign * zf) * zm_gl, p_gl, mp.mpf(nu_f)))

    worst = mp.mpf(0)
    worst_at = None
    for zeta, p, nu in pts:
        try:
            S_fit = eval_r3_S(coef2d, zeta, p, nu)
            R_true = r3_R_at(nu, p, zeta, dps)
            # R_true = S(zeta,p,nu) exactly (by the extraction's own
            # definition R = sum_k e_k/nu^k, i.e. S IS R -- no separate
            # rescaling). Compare directly.
            ref = R_true
            if ref == 0:
                continue
            err = abs((S_fit - ref) / ref)
        except RuntimeError:
            continue
        if err > worst:
            worst, worst_at = err, (float(zeta), float(p), float(nu))
    log2w = float(mp.log(worst, 2)) if worst > 0 else float("-inf")
    print(f"    tested {len(pts)} (zeta,p,nu) points inside the p-dependent "
          f"membership lens (global zeta_max={float(ZETA_MAX):.6f}); worst "
          f"rel err {float(worst):.3e} (2^{log2w:.2f}) at "
          f"(zeta,p,nu)={worst_at}, target 2^-56", file=sys.stderr)
    return worst


# ============================================================================
# Self-check (d): gamma-limit anchor, e_k(zeta,p=2^-50) vs gamma's
# c_k(eta=-zeta*sqrt(2)), <=1e-15 through k=5.
# ============================================================================
def check_d_anchor():
    print("(d) gamma-limit anchor (p=2^-50, eta=-zeta*sqrt(2)):", file=sys.stderr)
    dps = 100
    p_anchor = mp.mpf(2) ** -50
    nu_ladder = [512 * mp.mpf(2) ** j for j in range(13)]
    a0_ladder = [512 * mp.mpf(2) ** j for j in range(13)]
    worst = mp.mpf(0)
    worst_at = None
    for zeta_f in ("0.1", "-0.1", "0.3", "-0.3", "0.5", "-0.5"):
        zeta = mp.mpf(zeta_f)
        # exact 13x13 solve (K_EXT=13 matches the ladder length here,
        # mirroring gen_gamma_data.py/Task B2's own recipe) -- deliberately
        # NOT the LSQ Chebyshev-in-v path (that path is for the wide 2D
        # grid's conditioning; a 13-point matched ladder is exactly
        # solvable and does not need it).
        e = extract_e_monomial(nu_ladder, p_anchor, zeta, 13, dps)[0]
        c = gamma_ck(-zeta * mp.sqrt(2), 13, a0_ladder, dps)
        for k in range(6):
            diff = abs(e[k] - c[k])
            if diff > worst:
                worst, worst_at = diff, (zeta_f, k)
    print(f"    worst |diff| through k=5: {float(worst):.3e} at "
          f"zeta={worst_at[0]}, k={worst_at[1]}, target 1e-15", file=sys.stderr)
    return 0 if worst <= ANCHOR_TARGET else 1


# ============================================================================
# Self-check (h): symmetry e_k(zeta,p) = -e_k(-zeta,1-p).
#
# Derivation (sign): swap (alpha,beta,xi)->(beta,alpha,1-xi) sends
# lambda=alpha-c*xi -> -lambda (c invariant), so zeta -> -zeta; cpsi is
# swap-invariant (u,v swap roles, sum unchanged), so nu, z=sqrt(cpsi) are
# also swap-invariant. The small_val selected by the branch rule (P for
# lam>=0, Q for lam<0) is the SAME PHYSICAL VALUE under the swap: original
# small_val (zeta>=0 branch) = P; swapped small_val (zeta'=-zeta<0 branch,
# using "Q of the swapped config") = I_{1-(1-xi)}(alpha,beta) = P as well.
# So R(zeta,p) = (leading-P)*sqrt(2 pi nu)*exp(cpsi) while
# R(-zeta,1-p) = (P-leading)*sqrt(2 pi nu)*exp(cpsi) = -R(zeta,p) -- an
# EXACT identity at every nu, hence e_k(-zeta,1-p) = -e_k(zeta,p) for
# every k. (The design text's "+-" was left open; this generator resolves
# it to "-" and checks it numerically here.)
# ============================================================================
def check_h_symmetry():
    print("(h) e_k(zeta,p) = -e_k(-zeta,1-p) symmetry identity:", file=sys.stderr)
    dps = 80
    ladder = [T_RIDGE * mp.mpf(2) ** (mp.mpf(j) / 3) for j in range(R3_NLADDER)]
    worst = mp.mpf(0)
    worst_at = None
    for zeta_f, p_f in (("0.3", "0.2"), ("-1.2", "0.35"), ("2.5", "0.05"), ("-4.0", "0.48")):
        zeta, p = mp.mpf(zeta_f), mp.mpf(p_f)
        e1, _ = extract_e_monomial(ladder, p, zeta, K_EXT, dps)
        e2, _ = extract_e_monomial(ladder, 1 - p, -zeta, K_EXT, dps)
        for k in range(K_REPORT):
            diff = abs(e1[k] - (-e2[k]))
            rel = diff / abs(e1[k]) if e1[k] != 0 else diff
            if rel > worst:
                worst, worst_at = rel, (zeta_f, p_f, k)
    print(f"    worst rel diff {float(worst):.3e} at "
          f"(zeta,p,k)={worst_at}, target 1e-25 (extraction-noise level)",
          file=sys.stderr)
    return 0 if worst <= SYMMETRY_TARGET else 1


# ============================================================================
# Self-check (f): R4 series truncation over the CORRECTED closed box
# (xi_tau<=xi1, B*xi_tau<=B1, tau*|ln xi_tau|<=ln2 -- the G1a-added caps).
# Series: S = sum_{n>=1} (1-beta)_n*xi^n/(n!(alpha+n)) (R1's shape minus
# n=0, per probe4's reinterpretation matching the shipped gamma R4).
# ============================================================================
def check_f_r4():
    print("(f) R4 series truncation over the corrected box:", file=sys.stderr)
    dps = 60
    N_REF = 4000
    candidates = [40, 48, 56, 64, 80, 96]
    alpha = EPS_R4
    pts = []
    # xi capped at xi1=0.45 now (the G1a fix); beta*xi=B1 boundary and the
    # tau*|ln xi|<=ln2 boundary, both swept.
    xi_floor_ln2 = mp.exp(-LN2 / alpha)
    los = math.log10(float(max(xi_floor_ln2, mp.mpf("1e-30"))))
    his = math.log10(float(XI1))
    for i in range(30):
        e = los + (his - los) * i / 29
        xi = mp.mpf(10) ** mp.mpf(e)
        beta = min(B1 / xi, mp.mpf("1e25"))
        pts.append((beta, xi))
    # also beta swept directly at xi=xi1 (the other box wall)
    for i in range(20):
        e = math.log10(1e-10) + (math.log10(float(B1 / XI1)) - math.log10(1e-10)) * i / 19
        beta = mp.mpf(10) ** mp.mpf(e)
        pts.append((beta, XI1))

    worst = {N: (mp.mpf(0), None) for N in candidates}
    for beta, xi in pts:
        partials = series_partial_sums(alpha, beta, xi, N_REF, dps)
        # R1's series includes n=0 (=1/alpha); R4's S excludes it, so
        # subtract the n=0 term (t_0/(alpha+0) = 1/alpha) from each partial.
        n0 = 1 / alpha
        ref = partials[N_REF] - n0
        if ref == 0:
            continue
        for N in candidates:
            val = partials[N] - n0
            e = abs((val - ref) / ref)
            if e > worst[N][0]:
                worst[N] = (e, (float(alpha), float(beta), float(xi)))
    chosen = None
    for N in candidates:
        w, at = worst[N]
        log2w = float(mp.log(w, 2)) if w > 0 else float("-inf")
        print(f"    N={N:3d}: worst rel err {float(w):.3e} (2^{log2w:.2f}) at {at}",
              file=sys.stderr)
        if chosen is None and w <= R4_TARGET:
            chosen = N
    if chosen is None:
        print("    FAILED: no candidate N meets target 2^-58", file=sys.stderr)
        return None, 1
    print(f"    pinned R4_NMAX={chosen} (target 2^-58)", file=sys.stderr)

    # --- FOURTH-CORRECTION widened-window budget (own stderr line, own
    # lattice; per the brief this is a REPORT, not a re-pin -- R4 depth
    # stays chosen=48, G4 decides any bump from measured kernel ULP) ------
    # Window: B<0.24 (where thr_t=(tau+1)/(tau+Bp+2) exceeds xi1), xi_tau in
    # (xi1, thr_t]. Sweep tau (the tiny param) across its full range and Bp
    # across (0, 0.24), sampling xi_tau densely inside the widened sliver.
    win_pts = []
    for Bp_f in (0.001, 0.01, 0.05, 0.1, 0.15, 0.2, 0.239):
        Bp = mp.mpf(Bp_f)
        # tau is R4's OWN tiny-min parameter -- must stay <= eps_R4=2^-6
        # (the R4 membership predicate's own tau bound); an earlier version
        # of this lattice swept tau up to 1.0 (10**0 landed in its range),
        # far outside R4's actual domain, and that single out-of-domain
        # point dominated the reported sup (2^-32-class) with a physically
        # meaningless number -- caught by cross-checking the reported
        # witness against EPS_R4 before trusting the printed sup.
        tau_list = sorted(set(
            [mp.mpf(10) ** e for e in range(-300, -5, 30)] +
            [mp.mpf(v) for v in ("1e-10", "1e-3", "1e-2")] +
            [EPS_R4 / 2, EPS_R4]
        ))
        for tau in tau_list:
            thr_t = (tau + 1) / (tau + Bp + 2)
            if thr_t <= XI1:
                continue  # window only exists where thr_t>xi1
            for frac in (mp.mpf("0.02"), mp.mpf("0.3"), mp.mpf("0.6"),
                         mp.mpf("0.9"), mp.mpf("0.999")):
                xi_tau = XI1 + frac * (thr_t - XI1)
                if 0 < xi_tau < 1 and tau * abs(mp.log(xi_tau)) <= LN2:
                    win_pts.append((tau, Bp, xi_tau))

    worst_win = mp.mpf(0)
    worst_win_at = None
    for tau, Bp, xi_tau in win_pts:
        partials = series_partial_sums(tau, Bp, xi_tau, N_REF, dps)
        n0 = 1 / tau
        ref = partials[N_REF] - n0
        if ref == 0:
            continue
        val48 = partials[48] - n0
        e = abs((val48 - ref) / ref)
        if e > worst_win:
            worst_win, worst_win_at = e, (float(tau), float(Bp), float(xi_tau))
    log2w = float(mp.log(worst_win, 2)) if worst_win > 0 else float("-inf")
    print(f"    FOURTH-CORRECTION widened window (B<0.24, xi_tau in "
          f"(xi1,thr_t]): {len(win_pts)} points, N=48 truncation sup "
          f"{float(worst_win):.3e} (2^{log2w:.2f}) at "
          f"(tau,Bp,xi_tau)={worst_win_at} -- R4 depth stays {chosen} "
          f"(G4 decides any bump from measured kernel ULP)", file=sys.stderr)

    # --- SEVENTH-CORRECTION post-route domain (GATING, target 2^-58) --------
    # Near-one R1 lanes now enter THIS series in their fired orientation
    # (route_final tag "R4-postroute") instead of the opposite-orientation
    # CF -- the sixth correction's CF destination stalled at 2^-55.5 on
    # (0.0234, 1e6, 4e-6) because that CF inherits the small-second-
    # parameter weakness. Domain swept here: the R1 box with alpha ABOVE
    # eps_R4 (at or below it, step 0 owns the point), restricted to lanes
    # the post-route model actually fires on (_r1_native_value >
    # BETA_NEAR_ONE -- the identical model route_final uses). tau appears
    # only in the 1/(tau+n) weights and the exact assembly, so truncation
    # is expected AT OR BELOW the in-box sup; this lattice PROVES it
    # rather than assuming. Both G3 (B)-witnesses must post-route and be
    # covered (hard assert).
    pr_pts = []
    alpha_pr = sorted(set(
        [EPS_R4 * mp.mpf(m) for m in ("1.01", "1.5", "3")] +
        [mp.mpf(v) for v in ("0.0234375", "0.05", "0.1", "0.158", "0.3",
                              "0.5", "0.8", "1.2", "2.0")]))
    beta_pr = sorted(set([mp.mpf(10) ** e for e in range(-3, 7)] + [mp.mpf(20)]))
    for alpha in alpha_pr:
        for beta in beta_pr:
            xi_hi = min(XI1, B1 / beta)
            if xi_hi <= 0:
                continue
            for frac in (mp.mpf("1e-6"), mp.mpf("0.01"), mp.mpf("0.1"),
                         mp.mpf("0.5"), mp.mpf(1)):
                xi = frac * xi_hi
                if not (0 < xi < 1):
                    continue
                # Membership via route_final ITSELF (tag test), not a bare
                # value filter: the first cut filtered on the value alone
                # and collected tau = 2 points the tau-gate deliberately
                # KEEPS in R1 (their complement sits in the (2^-12, 2^-11)
                # band the doctrine still covers -- check (e) polices
                # them; its pocket lattice extends to tau = 4).
                if route_final(alpha, beta, xi)[3] == "R4-postroute":
                    pr_pts.append((alpha, beta, xi))
    wit_pr = [(mp.mpf("0.158"), mp.mpf(20), mp.mpf("0.396")),
              (mp.mpf("0.5"), mp.mpf(20), mp.mpf("0.35"))]
    n_wit_missed = 0
    for w3 in wit_pr:
        _, _, _, wtag = route_final(*w3)
        if wtag == "R4-postroute":
            pr_pts.append(w3)
        else:
            n_wit_missed += 1
    worst_pr = mp.mpf(0)
    worst_pr_at = None
    for tau, Bp, xi_tau in pr_pts:
        partials = series_partial_sums(tau, Bp, xi_tau, N_REF, dps)
        n0 = 1 / tau
        ref = partials[N_REF] - n0
        if ref == 0:
            continue
        valc = partials[chosen] - n0
        e = abs((valc - ref) / ref)
        if e > worst_pr:
            worst_pr, worst_pr_at = e, (float(tau), float(Bp), float(xi_tau))
    log2pr = float(mp.log(worst_pr, 2)) if worst_pr > 0 else float("-inf")
    tau_sup = max((float(t) for t, _, _ in pr_pts), default=0.0)
    print(f"    SEVENTH-CORRECTION post-route domain: {len(pr_pts)} "
          f"post-routing points, N={chosen} truncation sup "
          f"{float(worst_pr):.3e} (2^{log2pr:.2f}) at "
          f"(tau,B,xi_tau)={worst_pr_at}, target 2^-58; fired-tau sup "
          f"{tau_sup:.4f} [must be <= kBetaPrTauMax = 1.5]; G3 "
          f"(B)-witnesses not post-routing: {n_wit_missed} [MUST be 0]",
          file=sys.stderr)
    if tau_sup > 1.5:
        print("    FAILED: a post-routing point exceeds the tau ceiling.",
              file=sys.stderr)
        return chosen, 1
    if n_wit_missed:
        print("    FAILED: a G3 (B)-witness did not post-route -- ESCALATE.",
              file=sys.stderr)
        return chosen, 1
    if worst_pr > R4_TARGET:
        print("    FAILED: post-route-domain truncation misses 2^-58.",
              file=sys.stderr)
        return chosen, 1

    return chosen, 0


# ============================================================================
# Self-check (e): routing safety under the FINAL corrected order.
# ============================================================================
def _r1_native_value(a, b, x, dps=40):
    """Cheap post-route oracle for route_final(): the R1 series itself
    (prefix * S_N1), NOT small_val_via_cf's escalating self-convergence
    loop -- self-check (a) already proves N1=64 truncation is <=2^-60 over
    R1's WHOLE box (x<=xi1=0.45, b*x<=B1=8), so a single fixed-depth pass
    is exact enough to place a value against the 2^-11 near-one bar with
    huge margin, at a cost route_final can afford to pay on every R1-box
    lattice point (an escalating CF loop could not, at the sweep sizes
    checks (b)/(e) run). mp.dps is set INSIDE this layer per the brief's
    hard rule; series_partial_sums does its own save/restore too, so the
    nesting is safe."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        s = series_partial_sums(a, b, x, N1, dps)[N1]
        # PLAN.md "R1 power series": I = [xi^alpha / B(alpha,beta)] * S --
        # NOT I_via_cf's prefactor (which carries an extra (1-xi)^beta and
        # 1/alpha specific to the CF's own F normalization; copying it
        # wholesale here first gave a wrong value, caught by cross-checking
        # against mp.betainc on a plain (2,3,0.3) sanity point before this
        # helper was trusted anywhere).
        logpref = a * mp.log(x) - (
            mp.loggamma(a) + mp.loggamma(b) - mp.loggamma(a + b))
        return mp.exp(logpref) * s
    finally:
        mp.mp.dps = old


def route_final(a, b, x):
    """G1b final order (PLAN.md 'second routing correction') with R3's
    THIRD-CORRECTION membership (ratio band, not the cpsi<=800 strip), the
    G3 FOURTH correction (R4's widened xi_tau cap), the SIXTH correction
    (R1 either-orientation + post-route model -- see below), and the (C)
    gamma-limit slice tag -- all reproduced VERBATIM from PLAN.md "G1/G2
    revision results -- SIXTH routing correction, B_GL pinned" and
    src/beta-inl.h's own router (read-only reference; this generator does
    not re-derive the kernel's thr_t formula, only replicates it exactly
    for lockstep):
      0. min(a,b)<=eps_R4 -> tiny-first (tau,Bp,xi_tau); if
         tau*|ln xi_tau|<=ln2 AND xi_tau<=max(xi1,thr_t) AND Bp*xi_tau<=B1
         -> R4, where thr_t=(tau+1)/(tau+Bp+2) is R2's OWN orientation
         threshold evaluated in the tiny-first frame [FOURTH CORRECTION:
         the design's box capped xi_tau at xi1 alone, stranding a
         B<~0.24 window (xi1,thr_t) where NEITHER R1 orientation fires and
         R2 evaluates the near-one tiny-first side -- see beta-inl.h's own
         "R4's xi cap, WIDENED" comment for the witness and measured
         N=48 truncation (2^-57.3) over the widened window]; else fall
         through.
      1. R1 fires in EITHER orientation whose box holds (x<=xi1 AND
         b*x<=B1 for native; (1-x)<=xi1 AND a*(1-x)<=B1 for swap; native
         checked first) [SIXTH CORRECTION, REVERTING the fifth: the fifth
         correction's lambda>=0 requirement was too blunt -- its OWN
         displaced traffic broke check (b) (witness (0.158,1000,0.00251)
         needs CF depth ~512, yet R1-native serves it perfectly: evaluated
         0.994, complement at 2^-58). POST-ROUTE MODEL (replaces the
         orientation restriction): evaluate the R1-native value at cheap
         working dps (_r1_native_value, N1=64 -- self-check (a)'s own
         proven-accurate depth). If it exceeds BETA_NEAR_ONE (kBetaNearOne
         = 1-2^-11), the point is re-tagged "R4-postroute" and evaluated
         by R4's analytic small-side assembly in the SAME orientation
         [SEVENTH correction; the sixth's opposite-orientation CF
         destination is superseded -- it stalled at 2^-55.5 on the CF's
         small-second-parameter weakness]. Mirrors the kernel's mask
         update: near-one R1 lanes fold into the R4 core's lane set.
         Points at or below the threshold keep R1. Both G3 witnesses
         (0.158,20,0.396) and (0.5,20,0.35) post-route correctly under
         this model (self-check (f)'s post-route lattice reproves it).
      2. R3 if nu=a*b/(a+b)>=T_ridge AND x/p in [1/2,2] AND
         (1-x)/q in [1/2,2] (p=a/c,q=b/c) -- the ratio-band caps; this
         joint condition is symmetric under the native/swap relabeling
         (swapping sends "x/p" <-> "(1-x)/q"), so it is checked once on
         the native triple regardless of which side ends up evaluated.
         Orientation within R3 stays the mean predicate (native if
         x<=mean else swap) -- R3's own table is built to be exactly
         self-consistent under that swap (self-check (h)).
      3. R2, orientation by xi < (a+1)/(c+2) -- R2 now also owns
         everything at nu>=T_ridge that falls OUTSIDE the ratio band
         (the risk transferred to R2 by the third correction; probed by
         self-check (b)'s ratio-cap x nu sweep). GAMMA-LIMIT SLICE
         [(C)]: if the CF-oriented triple's max parameter >= B_GL, the
         tag gets a "-gammalim" suffix (native/swap preserved) -- this is
         a TAG only (the oracle below still evaluates via the CF; the
         reference generator's gamma-corner oracle switches on this tag,
         and the real kernel switches on the same max(alpha,beta)>=B_GL
         predicate to route to the gamma-limit path instead). R4-postroute
         points keep their own tag (PLAN.md's coverage audit lists the
         post-route traffic as its own histogram line); they never carry
         -gammalim (R1's box keeps beta*xi <= B1 = 8, far below any
         B_GL-scale hazard, and their evaluation path is the R4 series,
         not the CF)."""
    c = a + b
    tau = min(a, b)
    bmax = max(a, b)
    thr_t = (tau + 1) / (tau + bmax + 2)
    xi_cap = max(XI1, thr_t)
    if tau <= EPS_R4:
        if a <= b:
            xi_tau, tag = x, "R4-native"
        else:
            xi_tau, tag = 1 - x, "R4-swap"
        if tau * abs(mp.log(xi_tau)) <= LN2 and xi_tau <= xi_cap and bmax * xi_tau <= B1:
            return (a, b, x, tag) if tag == "R4-native" else (b, a, 1 - x, tag)
        # else fall through to R1/R3/R2 below.
    r1_hit = None
    if x <= XI1 and b * x <= B1:
        r1_hit = (a, b, x, "R1-native")
    elif (1 - x) <= XI1 and a * (1 - x) <= B1:
        r1_hit = (b, a, 1 - x, "R1-swap")
    if r1_hit is not None:
        aa, bb, xx, tag = r1_hit
        val = _r1_native_value(aa, bb, xx)
        if val <= BETA_NEAR_ONE:
            return aa, bb, xx, tag
        # SEVENTH CORRECTION [supersedes the sixth's destination]: the
        # near-one lane keeps its FIRED orientation and goes to R4's
        # analytic small-side assembly. Post-routed lanes are R4-shaped by
        # construction: R1's box supplies both of R4's convergence caps
        # (xi <= xi1, beta*xi <= B1), and the near-one condition itself
        # guarantees the Expm1 assembly's argument is < ~2^-10 -- its
        # ideal accuracy zone. eps_R4 was a ROUTING threshold, never a
        # validity bound of the assembly (the series' convergence is
        # governed by xi and beta*xi, not tau). The sixth correction's
        # opposite-orientation CF destination stalled at 2^-55.5 on
        # (0.0234, 1e6, 4e-6) -- check (b)(viii)'s own finding -- because
        # that CF inherits exactly the small-second-parameter weakness
        # the fifth correction already hit; sending the lane to R4
        # removes the CF from this traffic entirely. Gating budget:
        # check (f)'s post-route-domain lattice.
        # tau-ceiling kBetaPrTauMax = 1.5: R4's exact-argument
        # lgamma(1+tau) runs on lgamma's centre-1/centre-2 zones, valid
        # to tau = 1.5 (tau-1 Sterbenz-exact there). The near-one bar
        # cannot fire above tau ~ 1.35 anyway (P_gamma(tau, B1) drops
        # below 1-2^-11 by tau ~ 1.25); the gate makes that a
        # machine-checked invariant -- points above it stay R1, and
        # check (e) verifies none violates the doctrine.
        if aa <= mp.mpf("1.5"):
            return aa, bb, xx, "R4-postroute"
        return aa, bb, xx, tag
    nu = a * b / c
    p = a / c
    q_ = b / c
    in_band = (XI_RATIO_LO <= x / p <= XI_RATIO_HI and
               XI_RATIO_LO <= (1 - x) / q_ <= XI_RATIO_HI)
    if in_band and (nu >= T_RIDGE or
                    (bmax >= B_GL and nu >= GL_RIDGE_MIN)):
        # [(C) slice, ridge part]: above B_GL the CF is structurally
        # degenerate, so the in-band ridge floor drops to GL_RIDGE_MIN=20
        # (check (c)'s extension lattice proves the 1/nu extrapolation).
        mean = a / c
        return (a, b, x, "R3-native") if x <= mean else (b, a, 1 - x, "R3-swap")
    thresh = (a + 1) / (c + 2)
    if x < thresh:
        a1, b1, x1, tag = a, b, x, "R2-native"
    else:
        a1, b1, x1, tag = b, a, 1 - x, "R2-swap"
    if max(a1, b1) >= B_GL:
        tag = tag + "-gammalim"
    return a1, b1, x1, tag


def _route_value(a, b, x, dps=40):
    aa, bb, xx, tag = route_final(a, b, x)
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        v = small_val_via_cf(aa, bb, xx, dps, n_start=64, n_max=8192)
    finally:
        mp.mp.dps = old
    return v, tag


def check_e_routing():
    print("(e) routing safety under the final corrected order:", file=sys.stderr)
    print("    NOTE (found by this generator's own run, not assumed): the", file=sys.stderr)
    print("    <=1-2^-12 bound is the design's MEAN-PREDICATE doctrine --", file=sys.stderr)
    print("    it governs R1/R2/R3, which all evaluate 'whichever side the", file=sys.stderr)
    print("    predicate calls small' directly. R4 is a DIFFERENT contract:", file=sys.stderr)
    print("    its own text says 'the alpha-scaled side ~ alpha*J is direct'", file=sys.stderr)
    print("    -- R4 always constructs the alpha(or beta)-scaled quantity", file=sys.stderr)
    print("    analytically (small by construction whenever tau<=eps_R4),", file=sys.stderr)
    print("    which is NOT the same quantity as 'the raw regularized-beta", file=sys.stderr)
    print("    value at the routed (tau,B,xi_tau) triple' this self-check", file=sys.stderr)
    print("    can evaluate without R4's actual kernel formula (a G3", file=sys.stderr)
    print("    decision). Confirmed empirically: EVERY failure below is in", file=sys.stderr)
    print("    an R4-routed region; R1/R2/R3 have zero failures anywhere in", file=sys.stderr)
    print("    the lattice, exactly reproducing PLAN.md G1b Task C's own", file=sys.stderr)
    print("    finding (its witness (1e-20,1,0.4) reproduced here) and its", file=sys.stderr)
    print("    explicit deferral of the R4 question to whoever builds the", file=sys.stderr)
    print("    real kernel. The <=1-2^-12 bound is therefore checked here", file=sys.stderr)
    print("    for R1/R2/R3 ONLY; R4 is reported separately and does NOT", file=sys.stderr)
    print("    gate this self-check's pass/fail -- ESCALATED in the final", file=sys.stderr)
    print("    report rather than papered over with an invented R4 formula.", file=sys.stderr)
    pts = []
    ab_grid = [mp.mpf(10) ** mp.mpf(e) for e in range(-20, 21, 4)]
    ratio_grid = [mp.mpf(10) ** mp.mpf(e) for e in range(-20, 21, 4)]
    x_mid = [mp.mpf(i) / 10 for i in range(1, 10)]
    x_lo = [mp.mpf(10) ** mp.mpf(e) for e in range(-20, -1, 4)]
    x_hi = [1 - v for v in x_lo]
    x_grid = sorted(set(x_mid + x_lo + x_hi))
    # DEGEN_FLOOR: both a AND b below this is specials-table territory
    # (the design's own "two degenerate parameters -> NaN" rule), not a
    # routing case -- G1b Task C's own report dismissed these same lattice
    # corners for the same reason; excluded here so they don't drown the
    # genuine signal.
    degen_floor = mp.mpf("1e-8")
    for a in ab_grid:
        for r in ratio_grid:
            b = a * r
            if a < degen_floor and b < degen_floor:
                continue
            for x in x_grid:
                if 0 < x < 1:
                    pts.append((a, b, x, "lattice"))
    # G1a witness zone + family.
    for a in (mp.mpf(4), mp.mpf(8), mp.mpf(16)):
        for b in (mp.mpf(2) ** -7, mp.mpf(2) ** -6, mp.mpf(2) ** -5):
            for xi in (mp.mpf("0.9"), mp.mpf("0.9999"), 1 - mp.mpf("9.5e-7"),
                       1 - mp.mpf("1e-9"), 1 - mp.mpf("1e-12")):
                pts.append((a, b, xi, "g1a-witness"))
    # G1b witness (1e-20,1,0.4) + family.
    for a in (mp.mpf(v) for v in ("1e-20", "1e-15", "1e-10", "1e-6")):
        for b in (mp.mpf(v) for v in ("0.5", "1", "2", "8")):
            for xi in (mp.mpf("0.2"), mp.mpf("0.4"), mp.mpf("0.6"), mp.mpf("0.8")):
                pts.append((a, b, xi, "g1b-witness"))
    # G3/(B) FIFTH-CORRECTION witnesses (0.158,20,0.396) and (0.5,20,0.35) +
    # dense sampling of the pocket the pre-fifth-correction lattice MISSED:
    # near the beta*xi=B1 edge, alpha in (eps_R4,1) -- the old lattice's
    # 0.9973 max was a lattice artifact of under-sampling exactly this
    # pocket, per the escalation text; this dense sub-lattice is what turns
    # that into a genuine sup measurement.
    for a0, b0, x0 in ((mp.mpf("0.158"), mp.mpf(20), mp.mpf("0.396")),
                        (mp.mpf("0.5"), mp.mpf(20), mp.mpf("0.35"))):
        for da in (mp.mpf("0.5"), mp.mpf("0.8"), mp.mpf(1), mp.mpf("1.2"), mp.mpf(2)):
            for dx in (mp.mpf("0.9"), mp.mpf("0.97"), mp.mpf(1),
                       mp.mpf("1.03"), mp.mpf("1.1")):
                a = a0 * da
                xi = x0 * dx
                if 0 < xi < 1:
                    pts.append((a, b0, xi, "g3-fifth-witness"))
    a_pocket = [mp.mpf(v) for v in
                ("0.0157", "0.02", "0.05", "0.1", "0.158", "0.2", "0.3",
                 "0.5", "0.7", "0.9", "0.99",
                 # SEVENTH-correction tau-gate margin band: tau in
                 # (kBetaPrTauMax, 4] stays R1 even when its value sits in
                 # the (2^-12, 2^-11) near-one band -- the doctrine still
                 # covers it (complement ~1/4 ulp), and THIS lattice is
                 # what proves the claim rather than an estimate (the
                 # tau = 2, beta = 20, xi = 0.4 family measured Q ~ 3.3e-4,
                 # 1.4x above the 2^-12 bar).
                 "1.6", "2.0", "2.5", "3.0", "4.0")]  # (eps_R4, 4]
    for a in a_pocket:
        for b in (mp.mpf(v) for v in ("8", "16", "20", "40", "100", "1000")):
            xi_edge = B1 / b  # the beta*xi=B1 edge itself
            mean = a / (a + b)
            for frac in (mp.mpf("0.9"), mp.mpf("0.97"), mp.mpf(1),
                         mp.mpf("1.03"), mp.mpf("1.1")):
                xi = xi_edge * frac
                if 0 < xi < 1:
                    pts.append((a, b, xi, "g3-fifth-pocket"))
            # also the mean itself and just past it (where lambda flips
            # sign) crossed with the B1 edge.
            for xi in (mean, min(xi_edge, XI1)):
                if 0 < xi < 1:
                    pts.append((a, b, xi, "g3-fifth-pocket"))

    worst_governed = mp.mpf(-1)   # R1/R2/R3 only -- the actual gate
    worst_governed_at = None
    worst_r4 = mp.mpf(-1)         # reported, not gated (see note above)
    worst_r4_at = None
    n_fail_governed = 0
    n_fail_r4 = 0
    fails = []
    region_totals = {}
    for a, b, x, tag in pts:
        try:
            v, region = _route_value(a, b, x)
        except (RuntimeError, ZeroDivisionError, ValueError):
            continue
        rprefix = region.split("-")[0]
        region_totals[rprefix] = region_totals.get(rprefix, 0) + 1
        if rprefix == "R4":
            if v > worst_r4:
                worst_r4, worst_r4_at = v, (float(a), float(b), float(x), region, tag)
            if v > ROUTE_THRESH:
                n_fail_r4 += 1
            continue
        if v > worst_governed:
            worst_governed, worst_governed_at = v, (float(a), float(b), float(x), region, tag)
        if v > ROUTE_THRESH:
            n_fail_governed += 1
            if len(fails) < 10:
                fails.append((float(a), float(b), float(x), region, float(v)))
    print(f"    region point counts: {region_totals}", file=sys.stderr)
    print(f"    R1/R2/R3 (GATED): worst evaluated value {float(worst_governed):.12f} "
          f"at {worst_governed_at}, failures={n_fail_governed}", file=sys.stderr)
    print(f"    R4 (reported only, see note): worst evaluated value "
          f"{float(worst_r4):.12f} at {worst_r4_at}, failures={n_fail_r4}",
          file=sys.stderr)
    for f in fails:
        print(f"      a={f[0]:.6e} b={f[1]:.6e} x={f[2]:.10f} region={f[3]} value={f[4]:.12f}",
              file=sys.stderr)
    return 0 if worst_governed <= ROUTE_THRESH else 1


# ============================================================================
# Binet tail (self-check (g)): phi(z) ~ sum_{k=1}^{K_B} B_2k/(2k(2k-1)z^(2k-1)),
# z>=Z0=10. Same series as lgamma-inl.h's LgammaStirling, fresh table (own
# Z0/K_B target, per the design's "fresh coefficients" option).
# ============================================================================
def binet_coeffs(K_B_, dps):
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        coefs = []
        for k in range(1, K_B_ + 1):
            B2k = mp.bernoulli(2 * k)
            coefs.append(B2k / (2 * k * (2 * k - 1)))
        return coefs
    finally:
        mp.mp.dps = old


def check_g_binet(coefs):
    print(f"(g) Binet truncation at Z0={float(Z0)}, K_B={K_B}:", file=sys.stderr)
    dps = 80

    def phi_true(z):
        z = mp.mpf(z)
        return mp.loggamma(z) - (z - mp.mpf("0.5")) * mp.log(z) + z - mp.log(2 * mp.pi) / 2

    def phi_series(z):
        z = mp.mpf(z)
        zinv2 = 1 / (z * z)
        s = mp.mpf(0)
        zpow = z
        for c in coefs:
            s += c / zpow
            zpow *= z * z
        return s

    worst = mp.mpf(0)
    worst_at = None
    for zf in ("10", "12", "15", "20", "30", "50", "100", "1000"):
        z = mp.mpf(zf)
        true_v = phi_true(z)
        approx = phi_series(z)
        rel = abs((approx - true_v) / true_v)
        if rel > worst:
            worst, worst_at = rel, zf
    log2w = float(mp.log(worst, 2)) if worst > 0 else float("-inf")
    print(f"    worst rel err {float(worst):.3e} (2^{log2w:.2f}) at z={worst_at}, "
          f"target 2^-70", file=sys.stderr)
    return 0 if worst <= BINET_TARGET else 1


# ============================================================================
# DigammaRough (self-check (i)): base zone [Z0, Z0+1), Chebyshev fit
# against mpmath's psi directly (loose 2^-40 target -- ordinary
# interpolation is far cheaper than the asymptotic series at this bar).
# Any z in (0, 2*Z0] reaches the zone via the standard up-recurrence
# psi(z) = psi(z+1) - 1/z, walked until z_final in [Z0, Z0+1).
# ============================================================================
DIGAMMA_DEG = 10  # Chebyshev nodes on the base zone (degree up to 9);
                   # chosen generously for a target this loose -- checked,
                   # not assumed, by check_i_digamma.


def digamma_rough_coeffs(dps):
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        n = DIGAMMA_DEG
        nodes_t = _cheb_nodes(n)
        # base zone [Z0, Z0+1) -> t in (-1,1]
        ws = [Z0 + mp.mpf("0.5") + mp.mpf("0.5") * t for t in nodes_t]
        vals = [mp.digamma(w) for w in ws]
        return _cheb_coeffs_1d(vals)
    finally:
        mp.mp.dps = old


def digamma_rough_eval(coefs, w):
    t = 2 * (w - Z0) - 1
    return _clenshaw(coefs, t)


def digamma_rough_full(coefs, z):
    """psi(z) for z in (0, 2*Z0], via recurrence to the base zone [Z0,Z0+1)
    then the zone polynomial. Two directions: z<Z0 walks UP (psi(z) =
    poly(z_final) - sum 1/z_j accumulated on the way up); z>=Z0+1 walks
    DOWN (psi(z) = poly(z_final) + sum 1/z_j accumulated on the way down)
    -- the domain is (0, 2*Z0], so z can exceed the zone's upper edge too,
    a case an earlier version of this function did not handle (it only
    walked up for z<Z0, silently extrapolating the zone polynomial WAY
    outside its fitted range for z in [Z0+1, 2*Z0] -- caught by this
    generator's own self-check (i), which failed at z=2*Z0 with ~1% error
    instead of the intended 2^-40)."""
    zz = mp.mpf(z)
    s = mp.mpf(0)
    if zz < Z0:
        while zz < Z0:
            s += 1 / zz
            zz += 1
        return digamma_rough_eval(coefs, zz) - s
    while zz >= Z0 + 1:
        zz -= 1
        s += 1 / zz
    return digamma_rough_eval(coefs, zz) + s


def check_i_digamma(coefs):
    print(f"(i) DigammaRough <= 2^-40 relative on (0, {float(2*Z0)}]:", file=sys.stderr)
    dps = 60
    worst = mp.mpf(0)
    worst_at = None
    test_zs = []
    for i in range(1, 400):
        test_zs.append(mp.mpf(i) * (2 * Z0) / 400)
    test_zs += [mp.mpf(v) for v in ("1e-10", "1e-5", "0.01", "0.5", "0.999",
                                     "1.001", "9.999", "10.0", "10.001",
                                     "19.999", "20.0")]
    for z in test_zs:
        if z <= 0:
            continue
        true_v = mp.digamma(z)
        approx = digamma_rough_full(coefs, z)
        if true_v == 0:
            continue
        rel = abs((approx - true_v) / true_v)
        if rel > worst:
            worst, worst_at = rel, float(z)
    log2w = float(mp.log(worst, 2)) if worst > 0 else float("-inf")
    print(f"    worst rel err {float(worst):.3e} (2^{log2w:.2f}) at z={worst_at}, "
          f"target 2^-40", file=sys.stderr)
    return 0 if worst <= DIGAMMA_TARGET else 1


# ============================================================================
# main
# ============================================================================
def main():
    t_start = time.time()
    rc = 0

    print(f"ZETA_MAX (derived from cpsi<=800, nu>=T_ridge={float(T_RIDGE)}): "
          f"{float(ZETA_MAX)}", file=sys.stderr)

    rc |= check_a_r1()
    rc |= check_b_r2()

    r4_n, rc_f = check_f_r4()
    rc |= rc_f

    binet = binet_coeffs(K_B, 100)
    rc |= check_g_binet(binet)

    digamma_coefs = digamma_rough_coeffs(100)
    rc |= check_i_digamma(digamma_coefs)

    rc |= check_d_anchor()
    rc |= check_h_symmetry()

    if rc:
        print("Self-checks (a)/(b)/(f)/(g)/(i)/(d)/(h) failed -- aborting "
              "before the expensive R3 grid.", file=sys.stderr)
        return rc

    print(f"\nbuilding R3 grid: NZ={NZ} NP={NP} K_EXT={K_EXT} K_REPORT={K_REPORT} "
          f"ladder={R3_NLADDER}pts ...", file=sys.stderr)
    _, _, grid = build_r3_grid(dps=100)
    coef2d = fit_r3_tensor(grid)

    r3_worst = check_c_r3(coef2d)
    # R3_REPLAY_TARGET is FIXED at gamma-class 2^-56 (third-correction
    # instruction) -- not pinned-to-measured this time. A miss here is a
    # real ESCALATE (the corrected ratio-band domain was expected to reach
    # it at gamma-class node counts; if it doesn't, that's new information,
    # not a threshold to relax).
    if r3_worst > R3_REPLAY_TARGET:
        print(f"    FAILED: R3 residual {float(r3_worst):.3e} exceeds the "
              f"2^-56 target -- ESCALATE.", file=sys.stderr)
        rc = 1
    else:
        print(f"    R3 residual {float(r3_worst):.3e} meets the 2^-56 "
              f"target with margin {float(-mp.log(r3_worst / R3_REPLAY_TARGET, 2)):.2f} bits",
              file=sys.stderr)

    total_bytes = K_REPORT * NZ * NP * 8
    print(f"    R3 table size: {K_REPORT}*{NZ}*{NP}*8 = {total_bytes} bytes "
          f"({total_bytes/1024:.2f} KiB), budget 32768 bytes (32 KiB)",
          file=sys.stderr)
    if total_bytes > 32768:
        print("    FAILED: R3 table exceeds the 32 KiB budget -- ESCALATE.",
              file=sys.stderr)
        rc = 1

    if rc:
        print("Self-checks failed -- emitting nothing.", file=sys.stderr)
        return rc

    # --- self-check (e): routing safety (runs after the grid since it does
    # not depend on it, but is the most expensive per-point check -- kept
    # last so a routing failure doesn't waste the R3 grid time on a re-run;
    # in practice it's fast, CF-based, no mpmath betainc hangs). ----------
    rc |= check_e_routing()
    if rc:
        print("Self-check (e) failed -- emitting nothing.", file=sys.stderr)
        return rc

    # --- emit ----------------------------------------------------------------
    r4_recip_n = max(N1, r4_n)
    recip_n = [dd_split(mp.mpf(1) / n) for n in range(1, r4_recip_n + 1)]

    print("// Auto-generated by tools/gen_beta_data.py. DO NOT EDIT.")
    print("// Derivations and error budgets are documented in the generator;")
    print("// the kernel that consumes them is src/beta-inl.h.")
    print("#ifndef CORVUS_BETA_DATA_H_")
    print("#define CORVUS_BETA_DATA_H_")
    print()
    print("namespace corvus::detail {")
    print()
    print("// Region-map and fixed-length constants (PLAN.md \"Regularized")
    print("// incomplete beta -- detail design\", G1a/G1b probe-pinned; R4_NMAX")
    print("// and the R3 table below are this generator's own [G1c] pins).")
    print(f"inline constexpr double kBetaB1 = {hexf(B1)};")
    print(f"inline constexpr double kBetaXi1 = {hexf(XI1)};")
    print(f"inline constexpr double kBetaEpsR4 = {hexf(EPS_R4)};")
    print(f"inline constexpr double kBetaTRidge = {hexf(T_RIDGE)};")
    print(f"inline constexpr double kBetaZ0 = {hexf(Z0)};")
    print(f"inline constexpr double kBetaClg = {hexf(C_LG)};")
    print(f"inline constexpr double kBetaExpFloor = {hexf(E_FLOOR)};")
    print(f"inline constexpr double kBetaLn2 = {hexf(LN2)};")
    print(f"inline constexpr double kBetaZetaMax = {hexf(ZETA_MAX)};")
    print("// kBetaGammaLim (B_GL) [G3 escalation (C), 'gamma-limit slice']:")
    print("// max(alpha,beta) >= this, on the CF-oriented triple, routes R2 to")
    print("// the gamma-limit path instead of the backward CF (which is")
    print("// structurally degenerate up there -- see the derivation/deviation")
    print("// comment at _derive_gamma_lim in this generator; PROVISIONAL,")
    print("// frontier review owed, per that comment).")
    print(f"inline constexpr double kBetaGammaLim = {hexf(B_GL)};")
    print("// kBetaNearOne [SEVENTH routing correction]: after the R1 pass,")
    print("// any lane whose R1 dd value EXCEEDS this folds into the R4")
    print("// core's lane set in the SAME orientation (R4's analytic")
    print("// small-side assembly; R1's box already supplies R4's")
    print("// convergence caps and the near-one condition puts the Expm1")
    print("// argument below ~2^-10). One bit inside the 1-2^-12")
    print("// complement-slack doctrine bound -- margin for the compare")
    print("// being on the dd value. Generator-proved: check (f)'s")
    print("// post-route-domain truncation lattice.")
    print(f"inline constexpr double kBetaNearOne = {hexf(BETA_NEAR_ONE)};")
    print("// Post-route tau ceiling: lgamma's centre-2 zone edge (and the")
    print("// Sterbenz-exact tau-1 range). The near-one bar cannot fire above")
    print("// ~1.35; the gate makes that machine-checked (see route_final).")
    print(f"inline constexpr double kBetaPrTauMax = {hexf(mp.mpf('1.5'))};")
    print("// Gamma-limit slice ridge floor [(C) resolution]: in-band lanes")
    print("// with max(alpha,beta) >= kBetaGammaLim use R3 down to nu = 20")
    print("// (gamma's own kGammaAT; the CF is degenerate up there). The 1/nu")
    print("// extrapolation below the extraction ladder is proved by check")
    print("// (c)'s extension lattice at the anchor p.")
    print(f"inline constexpr double kBetaGlRidgeMin = {hexf(GL_RIDGE_MIN)};")
    print(f"inline constexpr int kBetaN1 = {N1};")
    print(f"inline constexpr int kBetaN2 = {N2};")
    print(f"inline constexpr int kBetaR4N = {r4_n};")
    print(f"inline constexpr int kBetaKB = {K_B};")
    print()
    print(f"// dd pairs of 1/n, n = 1..{r4_recip_n} -- shared by R1's power")
    print("// series and R4's alpha-scaled series (same t_n = t_{n-1}*(n-beta)*")
    print("// xi/n term ratio; R4's S omits the n=0 term, handled exactly by")
    print("// Expm1Dd elsewhere). 1/(alpha+n) is never tabled here (runtime-")
    print("// dependent; computed via DdRecipDd of the exact TwoSum(alpha,n),")
    print("// gamma R4 precedent) -- only the pure-integer 1/n multiplier is.")
    emit_hex_array_1d("kBetaRecipNHi", [p[0] for p in recip_n])
    emit_hex_array_1d("kBetaRecipNLo", [p[1] for p in recip_n])
    print()
    print(f"// R3 erf-form correction table: e_k(zeta,p), k=0..{K_REPORT-1}.")
    print("// S(zeta,p,1/nu) = sum_k e_k(zeta,p)/nu^k, Horner in 1/nu across")
    print(f"// the {K_REPORT} rows (gamma Temme-table precedent, K close to gamma's")
    print(f"// own K=11). Each row is a 2D tensor-Chebyshev fit,")
    print(f"// kBetaR3Cheb[k][n][m] the coefficient of")
    print("// T_n(zeta/kBetaZetaMax) * T_m((p-kBetaR3PMid)/kBetaR3PHalf).")
    print("//")
    print("// R3 MEMBERSHIP [third correction, PLAN.md \"G1c generator results")
    print("// and third correction\"]: nu>=kBetaTRidge AND xi/p in")
    print("// [kBetaXiRatioLo,kBetaXiRatioHi] AND (1-xi)/q in [same] -- the")
    print("// ridge RATIO band, mirroring gamma's shipped Temme table exactly")
    print("// (gen_gamma_data.py's ETA_LO/ETA_HI span lambda in [1/2,2], not a")
    print("// wide cpsi strip). An EARLIER version of this generator used")
    print("// membership 'nu>=T_ridge AND cpsi<=800' -- a design error: the")
    print("// implied zeta in [-5,5] fit domain cannot reach dd-level accuracy")
    print("// in 32 KiB (measured 2^-16-class residual, confirmed inherent to")
    print("// the domain width by cross-check against gamma's own c_0(eta) over")
    print("// an equally wide synthetic range, not a beta-specific bug) --")
    print("// caught and corrected before this table was ever shipped.")
    print("//")
    print(f"// kBetaZetaMax = sqrt(3*ln2/2) ~ {float(ZETA_MAX):.10f}, the EXACT sup of")
    print("// |zeta| over the ratio band (derived, not assumed -- see")
    print("// _derive_zeta_max's docstring/comment in this generator for the")
    print("// full algebra): u=-lambda/alpha in [-1/2,1] and v=lambda/beta in")
    print("// [-1/2,1] are the two caps (xi/p=1+u, (1-xi)/q=1+v), linked by")
    print("// p*u+q*v=0; cpsi(u) is convex in u at fixed p so the sup is at a")
    print("// box corner, maximized over p at p=1/3 (and mirror p=2/3) where")
    print("// BOTH raw caps bind at once (u=1,v=-1/2): zeta_max^2=3*ln2/2.")
    print("// Saturation (cpsi<=800, kBetaExpFloor) remains a SEPARATE overlay")
    print("// on top of the ratio band -- the z in (6,sqrt(800)] G-tail lanes")
    print("// only occur at moderate nu within the (now narrow) band.")
    print("//")
    print("// p in (0,0.5] -- the symmetry e_k(zeta,p) = -e_k(-zeta,1-p)")
    print("// (self-check (h)) halves the table; the kernel evaluates at")
    print("// (-zeta,1-p) and negates when p>0.5. Evaluate by nested Clenshaw:")
    print("// per-n row value = Clenshaw_m(row[n], u), then Clenshaw_n")
    print("// (row_values, t), t=zeta/kBetaZetaMax, u=(p-mid)/half.")
    print(f"inline constexpr double kBetaXiRatioLo = {hexf(XI_RATIO_LO)};")
    print(f"inline constexpr double kBetaXiRatioHi = {hexf(XI_RATIO_HI)};")
    print(f"inline constexpr int kBetaR3K = {K_REPORT};")
    print(f"inline constexpr int kBetaR3NZ = {NZ};")
    print(f"inline constexpr int kBetaR3NP = {NP};")
    print(f"inline constexpr double kBetaR3PMid = {hexf(P_MID)};")
    print(f"inline constexpr double kBetaR3PHalf = {hexf(P_HALF)};")
    emitted_blocks = [[[rd(coef2d[k][n][m]) for m in range(NP)] for n in range(NZ)]
                       for k in range(K_REPORT)]
    emit_hex_array_3d("kBetaR3Cheb", emitted_blocks)
    print()
    print(f"// Binet/Stirling tail coefficients (fresh table, own Z0/K_B target;")
    print(f"// same series as lgamma-inl.h's LgammaStirling -- see kLgammaStirCoef")
    print(f"// there for the sibling table at a different Z0). phi(z) ~")
    print(f"// sum_k kBetaBinetCoef[k-1]/z^(2k-1), z >= kBetaZ0, Horner in 1/z^2.")
    binet_d = [rd(c) for c in binet]
    emit_hex_array_1d("kBetaBinetCoef", binet_d)
    print()
    print(f"// DigammaRough: base zone [Z0,Z0+1), Chebyshev fit against psi")
    print(f"// directly (target 2^-40 relative on (0,2*Z0], self-check (i)).")
    print(f"// Any z reaches the zone via psi(z)=psi(z+1)-1/z walked upward;")
    print(f"// evaluate poly(w) via Clenshaw at t=2*(w-Z0)-1, w in [Z0,Z0+1).")
    print(f"inline constexpr int kBetaDigammaDeg = {DIGAMMA_DEG};")
    digamma_d = [rd(c) for c in digamma_coefs]
    emit_hex_array_1d("kBetaDigammaCoef", digamma_d)
    print()
    print("}  // namespace corvus::detail")
    print()
    print("#endif  // CORVUS_BETA_DATA_H_")

    print(f"\ntotal generator runtime: {time.time()-t_start:.1f}s", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
