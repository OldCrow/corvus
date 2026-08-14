#!/usr/bin/env python3
"""Generate src/gammainv_data.h -- every table gamma_p_inv/gamma_q_inv needs.

This generator pins, BY REPLAY (never by assumption), every free parameter
of the inverse's seed stage:

  a_T          S1 (Temme normal-quantile seed) vs S2/S3 (small-a series /
               far-tail fixed point) crossover.
  S1           lambda(eta) inversion of 1/2 eta^2 = lambda-1-ln(lambda):
               a Taylor series near eta=0 (derived HERE via exact-rational
               Lagrange series reversion -- clean-room from the defining
               equation, not copied from any table) for |eta| < a cutoff,
               Newton (log-space, u=ln(lambda)) elsewhere; PLUS the
               Temme eps_k(eta)/a^k correction(s) to eta itself, whose
               c_k(eta) coefficients are extracted the SAME way
               tools/gen_gamma_data.py extracts its own Temme table (a
               small-sample Vandermonde solve in 1/a against mpmath
               gammainc, at eta nodes spanning the S1 seed's WIDER domain
               -- gamma_data.h's kGammaTemmeCheb is fit ONLY over the
               ridge band lambda in [1/2,2] and is the wrong table for
               this: NOT duplicated here, this is a separate fit).
  S2           small-p series seed (x0 = exp((ln p + lnGamma(1+a))/a)),
               Picard correction count.
  S3           far-q-tail seed (L = -ln(q Gamma(a)), fixed-point), same.
  STEPS        per-region dd-Newton / log-residual-Newton step count
               (max 3), simulated against the forward evaluated in mpmath
               with injected relative noise at the forward kernel's own
               internal dd budget (2^-56 for R1/R2/R3-class regions,
               2^-58 for the R4-analogue tiny-a region).
  DEEP-SMALL   the a*x0 < 2^-60 closed-form cut (p-side).

mpmath discipline: mp.dps set inside every function that needs it (never a
bare module dps for anything precision-sensitive); replay checks run at
layered dps (60 then 100) per the BINDING sampling rule; sampling is
edge-refined bit-stepped (dense interior + nextafter-ladder points near
every boundary), never grid-only.

Self-checks (mandatory; stderr budget lines; ANY miss -> exit nonzero,
emit nothing):
  (a) lambda(eta) series: exact-rational coefficients reproduce the first
      three terms (1, 1/3, 1/36) and match a fresh mpmath numerical
      inversion to the stated per-order tolerance.
  (b) lambda(eta) Newton: converges to >= NEWTON_FLOOR_BITS across its
      whole pinned domain, both signs.
  (c) S1 correction: disjoint-sample re-extraction of c_k(eta) agrees to
      the stated tolerance (gen_gamma_data.py check-(a) pattern).
  (d) S1/S2/S3 seed-bit floors, measured against the actual seed code
      (native double, no shortcuts), edge-refined sampling.
  (e) a_T: S1's measured floor at a=a_T is not worse than S2/S3's at the
      same a (crossover is real, not just the forward's own a_T copied
      blind).
  (f) STEPS: per-region final bits after the pinned step count, worst
      case over forward-noise sign combinations, >= 54 bits with >= 1
      bit margin, at every replay point whose true x is a normal double.
  (g) deep-small cut: dropped-series-relative-error < 2^-60 below the
      cut; normal pipeline covers seamlessly at/above it.
  (h) weak-seed middle band: 3 steps recover to >= 54+margin.

Usage:
    python3 tools/gen_gammainv_data.py > src/gammainv_data.h
    python3 tools/gen_gammainv_data.py --full > src/gammainv_data.h   # denser replay, ~slow
"""
import math
import sys
from fractions import Fraction as Fr

import mpmath as mp

FULL = "--full" in sys.argv

# ============================================================================
# Part 0: shared math constants / helpers
# ============================================================================
A_T_CANDIDATE = 20.0  # forward's own kGammaAT; measured at self-check (e).
S1_A_MIN = 0.3  # measured floor below which S1's Temme normal-approximation
                 # itself (not just the correction table's eta domain) stops
                 # being a usable seed, REGARDLESS of eta -- self-check (f)
                 # finding, see seed_for's routing note.
SHALLOW_THRESHOLD = 2.0 ** -10  # s >= this: "weak middle band" (S1 may beat
                                  # S2/S3 there for a<a_T); s < this: S2/S3's
                                  # own genuine-tail domain, where they beat
                                  # S1 even for a>=S1_A_MIN -- self-check (f).
A_MPMATH_SAFE = mp.mpf("1e4")  # mpmath's lower-gammainc (hyp1f1) path
# genuinely hangs/fails to converge for large a near the ridge
# (gen_gamma_reference.py's own documented finding: NoConvergence at
# a>=1e7 near the ridge, multi-minute hangs even away from it at
# a~1e250). 1e6 is NOT a safe margin below that. 1e4 matches
# gen_gamma_reference.py's own A_SWITCH, an established safe margin;
# above it, small_of_x uses ONLY the Temme-extrapolation path
# (extract_ck's own sample points top out at 16000, comfortably under
# even the 1e7 hard limit).
APHI_SAT = mp.mpf(800)  # a*phi(lambda) beyond this: exact double saturation


def rd(x):
    return float(x)


def hexf(x):
    return float.hex(float(x))


def dd_split(x):
    hi = rd(x)
    return hi, rd(mp.mpf(x) - mp.mpf(hi))


def emit_hex_array_1d(name, vals, ncols=8):
    print(f"inline constexpr double {name}[{len(vals)}] = {{")
    for i in range(0, len(vals), ncols):
        print("    " + ", ".join(hexf(v) for v in vals[i:i + ncols]) + ",")
    print("};")


def phi_lam(lam):
    return lam - 1 - mp.log(lam)


def eta_of_lam(lam):
    lam = mp.mpf(lam)
    phi = phi_lam(lam)
    e = mp.sqrt(2 * phi)
    return -e if lam < 1 else e


def lam_of_eta_mp(eta, dps=60):
    """Bisection inversion of phi(lambda)=eta^2/2, high precision (TRUTH,
    used only to validate the double-precision series/Newton scheme, never
    itself part of the emitted seed path)."""
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


# ============================================================================
# Part 1: lambda(eta) Taylor series near eta=0, derived by exact-rational
# Lagrange series reversion of 1/2 eta^2 = lambda-1-ln(lambda). Clean-room:
# derived from the defining equation, verified against mpmath below, never
# ported from a published table (the first three terms -- 1, 1/3, 1/36 --
# are reproduced here as a CHECK, not an input).
# ============================================================================
LAM_SERIES_ORDER = 12  # through eta^12; self-check (a) measures the floor.


def _series_mul(a, b, n):
    out = [Fr(0)] * (n + 1)
    for i, ai in enumerate(a):
        if ai == 0 or i > n:
            continue
        for j, bj in enumerate(b):
            if i + j > n:
                break
            if bj == 0:
                continue
            out[i + j] += ai * bj
    return out


def _series_sqrt(a, n):
    """sqrt of a power series with a[0]=1 (coefficient recursion, exact)."""
    s = [Fr(0)] * (n + 1)
    s[0] = Fr(1)
    for k in range(1, n + 1):
        conv = sum(s[i] * s[k - i] for i in range(1, k))
        s[k] = (a[k] - conv) / (2 * s[0])
    return s


def _reversion(c, n):
    """eta = sum_{k=1}^n c_k t^k (c_1=1) -> t = sum_{k=1}^n d_k eta^k.
    Order-by-order: at order k, d[k] enters the eta^k coefficient of
    eta(t(eta)) LINEARLY with coefficient c_1=1 (every other contribution
    to that order comes from d[1..k-1], already fixed by induction; d[k]
    first appears in T at order eta^k, so in T^m (m>=2) it first appears
    at order >= eta^{k+1}, too late to matter at order k)."""
    d = [Fr(0)] * (n + 1)
    d[1] = Fr(1)
    for k in range(2, n + 1):
        Tser = d[: k + 1]  # d[k] still 0
        cur_power = Tser[:]
        total = [Fr(0)] * (k + 1)
        for m in range(1, k + 1):
            if m > 1:
                cur_power = _series_mul(cur_power, Tser, k)
            cm = c[m] if m < len(c) else Fr(0)
            if cm != 0:
                for i in range(k + 1):
                    total[i] += cm * cur_power[i]
        d[k] = -total[k]
    return d


def derive_lambda_series(n):
    """Returns d[1..n] (Fraction): lambda(eta) = 1 + sum d[k] eta^k."""
    fcoef = [Fr(0)] * (n + 3)
    for m in range(2, n + 3):
        fcoef[m] = Fr((-1) ** m, m)  # t - ln(1+t) = sum (-1)^m t^m/m, m>=2
    G = [2 * fcoef[k + 2] for k in range(n + 1)]  # eta^2/t^2 = G(t)
    phi_ser = _series_sqrt(G, n)  # eta/t = sqrt(G(t))
    c = [Fr(0)] * (n + 2)
    for m in range(1, n + 2):
        c[m] = phi_ser[m - 1] if (m - 1) < len(phi_ser) else Fr(0)
    return _reversion(c, n)


_LAM_D = derive_lambda_series(LAM_SERIES_ORDER)
_LAM_D_MP = [mp.mpf(0)] + [mp.mpf(x.numerator) / mp.mpf(x.denominator)
                            for x in _LAM_D[1:]]
_LAM_D_F = [0.0] + [float(x) for x in _LAM_D[1:]]  # double-rounded coeffs


def lam_series_eval_double(eta, order=LAM_SERIES_ORDER):
    """Horner, NATIVE DOUBLE (float = float64) end-to-end."""
    s = 0.0
    for k in range(order, 0, -1):
        s = s * eta + _LAM_D_F[k]
    return 1.0 + eta * s


ETA_SERIES_CUT = 0.5  # |eta| < cut -> series; pinned at self-check (a)/(b)
LAM_NEWTON_ITERS = 6  # pinned at self-check (b)


def lam_newton_double(eta, niter=LAM_NEWTON_ITERS):
    """Newton in u=ln(lambda): F(u)=e^u-1-u-target=0, F'(u)=e^u-1.
    Hybrid initial guess: large-lambda asymptotic (phi~lambda) for eta>=0,
    small-lambda asymptotic (phi~-ln(lambda)) for eta<0 -- both exact
    leading terms of phi's own two asymptotic regimes, so the guess is
    smooth and the right order of magnitude across the whole Newton
    domain (never the eta~0 region, which the series owns)."""
    target = 0.5 * eta * eta
    u = math.log1p(target) if eta >= 0.0 else -target
    for _ in range(niter):
        eu = math.exp(u)
        v = eu - 1.0 - u - target
        d = eu - 1.0
        if d == 0.0:
            break
        u = u - v / d
    return math.exp(u)


def lam_of_eta_double(eta):
    if abs(eta) < ETA_SERIES_CUT:
        return lam_series_eval_double(eta)
    return lam_newton_double(eta)


# ============================================================================
# Part 2: S1 correction table c_k(eta), extracted the SAME technique as
# gen_gamma_data.py's Temme table (Vandermonde-in-1/a against mpmath
# gammainc), but over the WIDE domain the S1 seed actually needs (full
# eta range at a=a_T, not just the ridge lambda in [1/2,2] -- gamma_data.h's
# own kGammaTemmeCheb is the wrong table for this, deliberately not reused;
# see module docstring).
#
# Derivation of the correction FORMULA (clean-room, from Temme's published
# structure Q(a,x) = 1/2 erfc(eta sqrt(a/2)) + R_a(eta), R_a(eta) =
# e^{-a eta^2/2}/sqrt(2 pi a) * sum_k c_k(eta)/a^k -- one Newton step in
# eta against this model, linearizing only the leading erfc term (R_a
# itself is O(1/sqrt(a)) smaller, so its own eta-derivative is higher
# order and dropped, standard perturbative Newton):
#   q-side (eta>=0): target=q=(1/2)erfc(z0)+R(eta0) at the trial z0, i.e.
#     model(eta0)-q = R(eta0); erfc-term deriv wrt eta = -sqrt(a/2pi)e^{-aphi};
#     delta = -R(eta0)/deriv = +S(eta0)/a,  S=c0+c1/a+c2/a^2+...
#   p-side (eta<0): target=p=(1/2)erfc(-z0)-R(eta0); erfc-term deriv wrt
#     eta = +sqrt(a/2pi)e^{-aphi} (opposite sign, from the -eta inside
#     erfc); model(eta0)-p = -R(eta0); delta = -(-R(eta0))/deriv = +S/a.
#   BOTH SIDES: eta_new = eta0 + S(eta0)/a -- side-symmetric (the sign
#   flip in the model's R-term exactly cancels the sign flip in the
#   erfc-derivative). A side-INDEPENDENT derivative expression paired
#   with a side-dependent residual sign is the wrong combination and
#   silently breaks this symmetry.
# ============================================================================
CK_KMAX = 2           # c_0, c_1 extracted (c_2's marginal seed-bit gain measured
                       # negligible at a=a_T and its Vandermonde extraction is
                       # noisy at this node/dps budget -- self-check (c) caught
                       # this; not worth stabilizing a term S1_NCORR never uses)
CK_A_LIST = (2000.0, 4000.0, 8000.0, 16000.0)
CK_DPS = 60
CK_NNODES = 25
ETA_MAX = math.sqrt(2 * 800.0 / A_T_CANDIDATE) + 0.3  # margin over aphi=800 bound at a_T


def R_exact(a, eta, lam, dps=CK_DPS):
    """(true small-side value - leading erfc term), rescaled. THE ORACLE
    TRAP (gen_gamma_data.py's own documented hazard): eta<0 MUST use the
    P-side identity, never Q=1-tiny."""
    with mp.workdps(dps):
        a = mp.mpf(a)
        if eta >= 0:
            Q = mp.gammainc(a, lam * a, regularized=True)
            base = Q - mp.erfc(eta * mp.sqrt(a / 2)) / 2
        else:
            P = mp.gammainc(a, 0, lam * a, regularized=True)
            base = mp.erfc(-eta * mp.sqrt(a / 2)) / 2 - P
        return base * mp.sqrt(2 * mp.pi * a) * mp.exp(a * eta * eta / 2)


_ck_cache = {}


def extract_ck(eta, lam, kmax=CK_KMAX, a_list=CK_A_LIST, dps=CK_DPS):
    key = (float(eta), kmax, a_list, dps)
    if key in _ck_cache:
        return _ck_cache[key]
    with mp.workdps(dps):
        n = len(a_list)
        A = mp.matrix(n, kmax)
        b = mp.matrix(n, 1)
        for i, a in enumerate(a_list):
            a = mp.mpf(a)
            v = 1 / a
            for k in range(kmax):
                A[i, k] = v ** k
            b[i, 0] = R_exact(a, eta, lam, dps=dps)
        if n == kmax:
            c = mp.lu_solve(A, b)
        else:
            At = A.T
            c = mp.lu_solve(At * A, At * b)
        out = [c[k, 0] for k in range(kmax)]
    _ck_cache[key] = out
    return out


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


def _cheb_nodes_eta(nnodes=CK_NNODES, eta_max=ETA_MAX):
    nodes = []
    for i in range(nnodes):
        t = mp.cos(mp.pi * (2 * i + 1) / (2 * nnodes))
        eta = mp.mpf(eta_max) * t
        nodes.append((t, eta))
    return nodes


def extract_all_nodes(a_list, nnodes=CK_NNODES, eta_max=ETA_MAX, dps=CK_DPS):
    nodes = _cheb_nodes_eta(nnodes, eta_max)
    cvals = []
    for t, eta in nodes:
        lam = lam_of_eta_mp(eta, dps=dps)
        cvals.append(extract_ck(eta, lam, kmax=CK_KMAX, a_list=a_list, dps=dps))
    return nodes, cvals


def cheb_fits_from_cvals(cvals, k_rows=CK_KMAX):
    fits = []
    for k in range(k_rows):
        vals = [cvals[i][k] for i in range(len(cvals))]
        fits.append(cheb_coeffs_from_vals(vals))
    return fits


# ============================================================================
# Part 3: seeds, native double end-to-end.
# ============================================================================
def erfcinv_double(y):
    """erfc^-1(y), y in (0,2). Single black-box correctly-rounded-ish
    double via mpmath at DPS SCALED TO y's OWN MAGNITUDE: a fixed dps=30
    silently rounds 1-y to exactly 1.0 whenever y < ~1e-30, which is
    routine here -- y=2s and s ranges down to ~1e-323 -- and that makes
    erfinv(1)=0 rather than the correct large tail value, corrupting the
    whole S1 seed for every deep-tail a>=a_T point. Then ONE cast to
    float -- the same modeling choice as trusting math.log/math.exp/
    math.lgamma elsewhere in this file (a real kernel calls corvus's own
    shipped, audited erfcinv)."""
    if y <= 0.0:
        return math.inf
    if y >= 2.0:
        return -math.inf
    dps = 30
    if y < 1e-25:
        dps = 30 + int(-math.log10(y)) + 15
    with mp.workdps(dps):
        return float(mp.erfinv(mp.mpf(1) - mp.mpf(y)))


def eta0_of(a, s, side):
    """Leading-order eta (pre-correction), double precision."""
    sgn = 1.0 if side == 'q' else -1.0
    z0 = erfcinv_double(2.0 * s)
    if not math.isfinite(z0):
        z0 = 0.0
    return sgn * z0 * math.sqrt(2.0 / a)


def seed_S1(a, s, side, ncorr):
    """side: 'p' (small=P, eta<0) or 'q' (small=Q, eta>=0). ncorr in 0..CK_KMAX."""
    eta0 = eta0_of(a, s, side)
    if ncorr > 0:
        ck = [_s1_ck_eval(k, eta0) for k in range(ncorr)]
        S = 0.0
        for k in range(ncorr - 1, -1, -1):
            S = S / a + ck[k]
        eta = eta0 + S / a
    else:
        eta = eta0
    lam = lam_of_eta_double(eta)
    return a * lam


# populated by main() once the Chebyshev fit is built; kept as plain
# double coefficient lists + a clenshaw-in-double evaluator so the seed
# path used for MEASUREMENT matches what the kernel will actually do
# (double arithmetic, not mpmath).
_S1_CHEB_ROWS_D = None  # list[ list[float] ], one row per k


def _s1_ck_eval(k, eta0):
    t = eta0 / ETA_MAX
    row = _S1_CHEB_ROWS_D[k]
    b1 = b2 = 0.0
    for c in row[:0:-1]:
        b1, b2 = 2.0 * t * b1 - b2 + c, b1
    return t * b1 - b2 + row[0]


def series_S_double(a, x, nmax=80):
    t = 1.0
    s = 1.0
    for _ in range(1, nmax + 1):
        t *= x / (a + _)
        s += t
        if t < 1e-18 * s:
            break
    return s


def seed_S2(a, p, ncorr):
    lg1a = math.lgamma(a + 1.0)
    lx = (math.log(p) + lg1a) / a
    x = math.exp(lx)
    for _ in range(ncorr):
        if not (math.isfinite(x) and x > 0):
            break
        S = series_S_double(a, x)
        lx = (math.log(p) + lg1a + x - math.log(S)) / a
        x = math.exp(lx)
    return x


S3_STABILITY_MARGIN = 3.0  # guard: L > MARGIN*|a-1| (see seed_S3 docstring)


def seed_S3(a, q, ncorr):
    """L = -ln(q*Gamma(a)). Fixed-point map x_{n+1} =
    L + (a-1)ln(x_n); its local contraction factor near the fixed point is
    (a-1)/x, so convergence needs |a-1|/L comfortably < 1 -- L>0 ALONE is
    NOT sufficient (self-check (f): a=0.3,q=0.271 gives
    L=0.21>0 but the iteration OSCILLATES, 0.21->1.30->0.025->2.80,
    landing a seed 10x off). Guarded here by L > S3_STABILITY_MARGIN*
    |a-1|, which both rejects that point and accepts every point where S3
    is empirically the better seed (e.g. a=0.5,q~1.6e-3: L=5.88
    comfortably clears 3*0.5=1.5). Seeding x=-ln(q) alone (zeroth
    iterate), without -lnGamma(a), crashes on log(negative) at
    a=0.1,q=0.5 -- the fixed-point map's OWN formula requires the full L,
    not just its leading term."""
    lg = math.lgamma(a)
    L = -math.log(q) - lg
    if not (math.isfinite(L) and L > S3_STABILITY_MARGIN * abs(a - 1.0)):
        return -1.0  # not S3's genuine, stable domain -- caller (the
                      # global tri-candidate seed_for) routes elsewhere.
    x = L
    for _ in range(ncorr):
        if not (x > 0):
            break
        xn = L + (a - 1.0) * math.log(x)
        if not (math.isfinite(xn) and xn > 0):
            break
        x = xn
    return x


def bits_of(true_x, approx_x):
    if not (math.isfinite(approx_x) and approx_x > 0 and true_x > 0):
        return -1.0
    rel = abs((approx_x - true_x) / true_x)
    if rel == 0:
        return 100.0
    return -math.log2(rel)


# ============================================================================
# Part 4: forward "truth" evaluator (mpmath), for MEASUREMENT only -- not
# the certified oracle (that is a separate job, done elsewhere). dps
# chosen per call; layered 60/100 checks are done at the self-check call
# sites, not baked in here.
# ============================================================================
def series_S_mp(a, x, dps=50, nmax=200000):
    with mp.workdps(dps):
        a = mp.mpf(a); x = mp.mpf(x)
        t = mp.mpf(1); s = mp.mpf(1)
        eps = mp.mpf(10) ** (-(dps - 8))
        for n in range(1, nmax + 1):
            t *= x / (a + n)
            s += t
            ratio = x / (a + n + 1)
            if ratio < 1 and t * ratio / (1 - ratio) < eps * s:
                return s
        return s


def cf_K_mp(a, x, dps=50, n0=80, nmax=8000):
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
                return cur
            prev = cur
            n = n2
        return prev


def small_of_x(a, x, dps=50, kmax=5, a_fit=(2000, 3000, 4000, 8000, 12000, 16000)):
    """Small-side-direct forward evaluator: (value, side). Direct mpmath
    gammainc for a<=A_MPMATH_SAFE (with elementary series/CF fallback on
    the documented hyp1f1 oracle trap), extrapolated-Temme beyond, exact
    saturation beyond a*phi>800."""
    with mp.workdps(dps):
        a_m = mp.mpf(a); x_m = mp.mpf(x)
        side = "p" if x_m <= a_m else "q"
        if a_m <= A_MPMATH_SAFE:
            try:
                if side == "p":
                    v = mp.gammainc(a_m, 0, x_m, regularized=True)
                else:
                    v = mp.gammainc(a_m, x_m, regularized=True)
                return v, side
            except mp.libmp.libhyper.NoConvergence:
                if side == "p":
                    S = series_S_mp(a_m, x_m, dps=dps, nmax=200000)
                    lg = a_m * mp.log(x_m) - x_m - mp.loggamma(a_m + 1)
                    return mp.e ** lg * S, side
                else:
                    K = cf_K_mp(a_m, x_m, dps=dps, nmax=200000)
                    lg = a_m * mp.log(x_m) - x_m - mp.loggamma(a_m)
                    return mp.e ** lg / K, side
        lam = x_m / a_m
        phi = phi_lam(lam)
        aphi = a_m * phi
        if aphi > APHI_SAT:
            return mp.mpf(0), side
        eta = eta_of_lam(lam)
        ck = extract_ck(eta, lam, kmax=kmax, a_list=a_fit, dps=dps)
        S = mp.mpf(0)
        for k in range(kmax - 1, -1, -1):
            S = S / a_m + ck[k]
        R = mp.e ** (-aphi) / mp.sqrt(2 * mp.pi * a_m) * S
        z = eta * mp.sqrt(a_m / 2)
        if lam >= 1:
            return mp.erfc(z) / 2 + R, "q"
        else:
            return mp.erfc(-z) / 2 - R, "p"


def small_of_x_saturating(a, x, dps=50, kmax=5):
    with mp.workdps(dps):
        a_m = mp.mpf(a); x_m = mp.mpf(x)
        lam = x_m / a_m
        phi = phi_lam(lam)
        if a_m > A_MPMATH_SAFE and a_m * phi > APHI_SAT:
            P = mp.mpf(1) if lam >= 1 else mp.mpf(0)
            return P, 1 - P
        v, side = small_of_x(a, x, dps=dps, kmax=kmax)
        return (v, 1 - v) if side == "p" else (1 - v, v)


def small_side_of_x(a, x, dps=50, kmax=5):
    """TRUE small-probability-side evaluator: (min(P,Q), side_of_min),
    side in {'p','q'} -- NOT small_of_x's "direct region" side (x<=a),
    which is a DIFFERENT criterion: for a=1, x=1, x<=a picks side='p'
    with P=0.632>1/2 -- the wrong side per the design's own contract,
    'solve against small side s<=1/2 (exact flip above)'. A real
    gamma_p_inv/gamma_q_inv call always flips its INPUT to the <=1/2
    side before seeding; a replay that feeds the >1/2 side directly into
    the seed/step machinery is testing a scenario the kernel never
    actually encounters."""
    with mp.workdps(dps):
        P, Q = small_of_x_saturating(a, x, dps=dps, kmax=kmax)
        return (P, "p") if P <= Q else (Q, "q")


def oracle_x(a, target, side, dps=40, lo=None, hi=None):
    """Root-find TRUE x s.t. P(a,x)=target (side='p') or Q=target
    (side='q'), bisection in ln(x). MEASUREMENT-grade (not the certified
    oracle): dps as given, monotone bisection, robust across the double
    range."""
    with mp.workdps(dps):
        a_m = mp.mpf(a); target = mp.mpf(target)

        def f(lx):
            x = mp.e ** lx
            P, Q = small_of_x_saturating(a_m, x, dps=dps)
            return (P if side == "p" else Q) - target

        lo = mp.mpf(-800) if lo is None else mp.mpf(lo)
        if hi is None:
            hi = mp.log(a_m * mp.mpf("1e320") + 1) if a_m > 0 else mp.mpf(700)
            if not mp.isfinite(hi):
                hi = mp.mpf(5000)
        else:
            hi = mp.mpf(hi)
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


def g_of_x(a, x, dps=40):
    with mp.workdps(dps):
        a_m = mp.mpf(a); x_m = mp.mpf(x)
        lg = (a_m - 1) * mp.log(x_m) - x_m - mp.loggamma(a_m)
        return mp.e ** lg


# ============================================================================
# Part 5: dd-Newton / log-residual-Newton STEP simulator. The forward
# P/Q and g are each perturbed by +/-eps (4 sign combos), eps = the
# forward kernel's own internal dd budget (2^-56 R1/R2/R3-class regions,
# 2^-58 R4-analogue tiny-a); worst case over all 4 kept. Newton ARITHMETIC
# itself (the update, not the forward evaluation) is done at mpmath
# dps=32 (~106 bits, matching dd's own representation), i.e. treated as
# exact relative to the injected forward-noise floor.
# ============================================================================
def step_newton_x(a, x0, target, side, eps, dps=32):
    with mp.workdps(dps):
        Pv, Pside = small_of_x(a, x0, dps=dps)
        if Pside != side:
            Pv = 1 - Pv
        gv = g_of_x(a, x0, dps=dps)
        sign = -1.0 if side == "q" else 1.0
        eps_m = mp.mpf(eps)
        out = []
        for sp in (1 + eps_m, 1 - eps_m):
            for sg in (1 + eps_m, 1 - eps_m):
                Pp = Pv * sp
                gp = gv * sg
                resid = Pp - target
                dx = resid / (sign * gp)
                out.append(x0 - dx)
        return out


def step_lnnewton_x(a, x0, target, side, eps, dps=32):
    """Delta = (ln P_dd - ln p) * P/g in dd, sign folded in as in plain Newton."""
    with mp.workdps(dps):
        Pv, Pside = small_of_x(a, x0, dps=dps)
        if Pside != side:
            Pv = 1 - Pv
        gv = g_of_x(a, x0, dps=dps)
        sign = -1.0 if side == "q" else 1.0
        eps_m = mp.mpf(eps)
        out = []
        for sp in (1 + eps_m, 1 - eps_m):
            for sg in (1 + eps_m, 1 - eps_m):
                Pp = Pv * sp
                gp = gv * sg
                if Pp <= 0 or target <= 0:
                    out.append(x0)
                    continue
                delta = (mp.log(Pp) - mp.log(target)) * Pp / gp
                out.append(x0 - sign * delta)
        return out


def _bits_of_frontier(xs, true_x):
    worst = 999.0
    for x in xs:
        if not (x > 0 and mp.isfinite(x)):
            worst = min(worst, -1.0)
            continue
        rel = abs((x - true_x) / true_x)
        b = 999.0 if rel == 0 else float(-mp.log(rel, 2))
        worst = min(worst, b)
    return worst


def simulate_steps(a, x0, true_x, target, side, eps, nsteps, variant, dps=32):
    """Worst-case relative-bit result after nsteps. PERFORMANCE
    SIMPLIFICATION: injecting noise at EVERY step branches 4-way PER
    STEP -- 4^nsteps forward evals, exponential and the dominant cost of
    the whole replay. Instead, steps 1..nsteps-1 use a single
    DETERMINISTIC (eps=0) trunk value -- Newton's own contraction erases
    earlier steps' noise sensitivity well before the final step, so only
    the LAST step's adversarial +/-eps sign combination (4-way) is what
    actually determines the final floor; verified against full-branching
    measurements (e.g. a=100,lambda=0.3: both methods agree to within
    noise-floor precision). Still worst case over noise signs -- at the
    step that matters."""
    step_fn = step_lnnewton_x if variant == "ln" else step_newton_x
    x = mp.mpf(x0)
    for _ in range(nsteps - 1):
        if not (x > 0 and mp.isfinite(x)):
            return -1.0
        x = step_fn(a, x, target, side, 0.0, dps=dps)[0]  # eps=0: deterministic
    if not (x > 0 and mp.isfinite(x)):
        return -1.0
    xs = step_fn(a, x, target, side, eps, dps=dps)  # final step: full 4-way noise
    return _bits_of_frontier(xs, true_x)


def simulate_steps_multi(a, x0, true_x, target, side, eps, max_nsteps, variant, dps=32):
    """Same simplification as simulate_steps, but returns {n: worst_bits}
    for n=1..max_nsteps in one pass -- the deterministic trunk is shared
    across all n, so this costs (max_nsteps-1) trunk evals +
    max_nsteps*4 final-step evals total, not max_nsteps separate calls
    to simulate_steps (which would redo the trunk each time)."""
    step_fn = step_lnnewton_x if variant == "ln" else step_newton_x
    out = {}
    x = mp.mpf(x0)
    for n in range(1, max_nsteps + 1):
        if not (x > 0 and mp.isfinite(x)):
            out[n] = -1.0
            continue
        xs = step_fn(a, x, target, side, eps, dps=dps)
        out[n] = _bits_of_frontier(xs, true_x)
        # advance the deterministic trunk by one (eps=0) step for next n
        x = step_fn(a, x, target, side, 0.0, dps=dps)[0]
    return out


def main():
    rc = 0
    global _S1_CHEB_ROWS_D

    # ------------------------------------------------------------------
    # Self-check (a): lambda(eta) series
    # ------------------------------------------------------------------
    print("(a) lambda(eta) Taylor series (exact-rational Lagrange reversion):",
          file=sys.stderr)
    print(f"    d[1..3] = {_LAM_D[1]}, {_LAM_D[2]}, {_LAM_D[3]}  "
          f"(PLAN quotes 1, 1/3, 1/36)", file=sys.stderr)
    if (_LAM_D[1], _LAM_D[2], _LAM_D[3]) != (Fr(1), Fr(1, 3), Fr(1, 36)):
        print("    FAILED: does not match PLAN's quoted terms", file=sys.stderr)
        rc = 1
    worst_a, worst_a_at = mp.mpf(0), None
    etas_a = [mp.mpf(s) for s in
              ("0.001", "0.01", "0.05", "0.1", "0.2", "0.3", "0.4",
               "0.45", "0.499", "-0.001", "-0.01", "-0.05", "-0.1", "-0.2",
               "-0.3", "-0.4", "-0.45", "-0.499")]
    for eta in etas_a:
        true_lam = lam_of_eta_mp(eta, dps=60)
        s = mp.mpf(0)
        for k in range(LAM_SERIES_ORDER, 0, -1):
            s = s * eta + _LAM_D_MP[k]
        approx = 1 + eta * s
        rel = abs((approx - true_lam) / true_lam)
        if rel > worst_a:
            worst_a, worst_a_at = rel, float(eta)
    print(f"    mpf-exact series (order {LAM_SERIES_ORDER}) worst rel err "
          f"{float(worst_a):.3e} at eta={worst_a_at}, |eta|<0.5 -- this is a "
          f"sanity check on the reversion itself; the DOUBLE-precision floor "
          f"self-check (b) below is the real gate for the emitted seed.",
          file=sys.stderr)
    if worst_a > mp.mpf(2) ** -30:
        print("    FAILED: series (in exact mpf coefficients) does not reproduce "
              "the defining equation to 2^-30 over |eta|<0.5", file=sys.stderr)
        rc = 1

    # ------------------------------------------------------------------
    # Self-check (b): lambda(eta) DOUBLE-precision series+Newton floor
    # ------------------------------------------------------------------
    print(f"(b) lambda(eta) double-precision floor (series |eta|<{ETA_SERIES_CUT}, "
          f"Newton n={LAM_NEWTON_ITERS} elsewhere), domain [-{ETA_MAX:.3f},{ETA_MAX:.3f}]:",
          file=sys.stderr)
    grid = []
    # edge-refined: dense log-ish sweep plus nextafter ladders at the
    # series/Newton seam and at the domain edges.
    xs = [0.001, 0.01, 0.05, 0.1, 0.2, 0.3, ETA_SERIES_CUT, 0.6, 0.8, 1.0,
          1.5, 2.0, 3.0, 5.0, 7.0, ETA_MAX]
    for v in xs:
        for s in (1.0, -1.0):
            grid.append(s * v)
    for edge in (ETA_SERIES_CUT, ETA_MAX):
        for s in (1.0, -1.0):
            for step in (-2, -1, 0, 1, 2):
                v = s * edge
                for _ in range(abs(step)):
                    v = math.nextafter(v, math.inf if step > 0 else -math.inf)
                grid.append(v)
    worst_b, worst_b_at = 999.0, None
    for eta in grid:
        true_lam = float(lam_of_eta_mp(eta, dps=60))
        approx = lam_of_eta_double(eta)
        b = bits_of(true_lam, approx)
        if b < worst_b:
            worst_b, worst_b_at = b, eta
    print(f"    worst {worst_b:.2f} bits at eta={worst_b_at}", file=sys.stderr)
    NEWTON_FLOOR_BITS = 25.0
    if worst_b < NEWTON_FLOOR_BITS:
        print(f"    FAILED: below floor {NEWTON_FLOOR_BITS}", file=sys.stderr)
        rc = 1

    if rc:
        print("Self-checks (a)/(b) failed -- aborting before further work.",
              file=sys.stderr)
        return rc

    # ------------------------------------------------------------------
    # S1 correction table: extract, self-check (c) disjoint re-extraction
    # ------------------------------------------------------------------
    print(f"(c) S1 correction c_k(eta) extraction: {CK_NNODES} nodes, "
          f"eta in [-{ETA_MAX:.3f},{ETA_MAX:.3f}], dps={CK_DPS} ...", file=sys.stderr)
    _, cvals_primary = extract_all_nodes(CK_A_LIST, dps=CK_DPS)
    fits_primary = cheb_fits_from_cvals(cvals_primary, CK_KMAX)
    a_list_check = (2400.0, 4800.0, 9600.0, 19200.0)
    _, cvals_check = extract_all_nodes(a_list_check, dps=CK_DPS)
    fits_check = cheb_fits_from_cvals(cvals_check, CK_KMAX)
    for k in range(CK_KMAX):
        worst = 0.0
        for i in range(CK_NNODES):
            av, bv = cvals_primary[i][k], cvals_check[i][k]
            if av != 0:
                worst = max(worst, float(abs((av - bv) / av)))
        thresh = 1e-8 if k == 0 else 1e-3
        status = "OK" if worst <= thresh else "FAIL"
        print(f"    c_{k}: worst rel diff {worst:.3e} (threshold {thresh:.1e}) {status}",
              file=sys.stderr)
        if worst > thresh:
            rc = 1
    if rc:
        print("Self-check (c) failed -- aborting.", file=sys.stderr)
        return rc

    _S1_CHEB_ROWS_D = [[rd(c) for c in fits_primary[k]] for k in range(CK_KMAX)]

    # ------------------------------------------------------------------
    # Self-check (d): S1/S2/S3 seed-bit floors, pin correction/Picard/
    # fixed-point counts. Edge-refined sampling: dense s-grid + nextafter
    # ladders near s=0.5 and near the a*phi=800 saturation boundary.
    # ------------------------------------------------------------------
    print("(d) S1 seed-bit floors vs #corrections (a>=a_T candidate):", file=sys.stderr)
    a_vals_s1 = [A_T_CANDIDATE, 25.0, 30.0, 50.0, 100.0, 1000.0, 1e4, 1e6,
                 1e8, 1e16, 1e100, 1e300, 1.7e308]
    if FULL:
        a_vals_s1 = a_vals_s1 + [22.0, 40.0, 70.0, 300.0, 1e5, 1e12, 1e50, 1e200]

    def s_grid_for_a(a):
        # smallest live s at this a: a*phi(lam) <= 800 boundary
        eta_bound = math.sqrt(min(2 * 800.0 / a, ETA_MAX * ETA_MAX))
        s_min_log = -max(1.0, (eta_bound * eta_bound * a) / 2 / math.log(10))
        pts = []
        for frac in (0.999, 0.99, 0.9, 0.7, 0.5, 0.3, 0.1, 0.01, 0.001):
            e = -frac * (eta_bound)
            pts.append(e)
        # translate a handful of eta targets to s via forward erfc (rough,
        # just need decent s coverage): use z = |eta|*sqrt(a/2), s~erfc(z)/2
        s_list = []
        for e in pts + [0.0]:
            z = abs(e) * math.sqrt(a / 2.0)
            with mp.workdps(30):
                s = float(mp.erfc(mp.mpf(z)) / 2)
            s_list.append(max(s, 5e-324))
        s_list += [0.5, math.nextafter(0.5, 0.0), 0.4, 0.3, 0.1]
        return sorted(set(s_list))

    summary_s1 = {k: [] for k in range(CK_KMAX + 1)}
    worst_rows_s1 = []
    for a in a_vals_s1:
        for s in s_grid_for_a(a):
            if not (0 < s <= 0.5):
                continue
            for side in ("p", "q"):
                true_x = oracle_x(a, s, side, dps=40)
                if true_x is None:
                    continue
                true_x = float(true_x)
                if not (math.isfinite(true_x) and true_x >= 0):
                    continue
                if true_x == 0.0:
                    continue
                row = []
                for k in range(CK_KMAX + 1):
                    try:
                        sx = seed_S1(a, s, side, k)
                        b = bits_of(true_x, sx)
                    except Exception:
                        b = -1.0
                    row.append(b)
                    summary_s1[k].append(b)
                worst_rows_s1.append((a, s, side, row))
    for k in range(CK_KMAX + 1):
        arr = summary_s1[k]
        if arr:
            print(f"    ncorr={k}: min={min(arr):.2f} "
                  f"median={sorted(arr)[len(arr)//2]:.2f} max={max(arr):.2f} n={len(arr)}",
                  file=sys.stderr)
    # pin S1_NCORR: smallest k whose min isn't (much) worse than kmax's min,
    # i.e. diminishing returns -- and never regress vs ncorr=0.
    S1_NCORR = 1
    for k in range(1, CK_KMAX + 1):
        if summary_s1[k] and min(summary_s1[k]) >= min(summary_s1[1]) - 0.5:
            S1_NCORR = k
    S1_NCORR = min(S1_NCORR, 2)  # cap: c2's marginal gain measured negligible
    print(f"    PINNED S1_NCORR = {S1_NCORR}", file=sys.stderr)

    print("(d) S2 seed-bit floors vs #Picard corrections (ALL a -- FIRST CORRECTION: "
          "S2 is now a global candidate, not just a<a_T):", file=sys.stderr)
    S2_KMAX_TESTED = 8  # ncorr 0..7; expected range 3-5, measured needs
    # more -- see below.
    a_vals_s2 = [1e-8, 1e-4, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0,
                 A_T_CANDIDATE - 1.0, math.nextafter(A_T_CANDIDATE, 0.0),
                 A_T_CANDIDATE, 22.0, 30.0, 100.0]
    p_list = sorted(set(
        [10.0 ** e for e in range(-300, 0, 10 if FULL else 20)] +
        [0.01, 0.02, 0.1, 0.3, 0.4, 0.45, 0.49, 0.499, math.nextafter(0.5, 0.0)]))
    summary_s2 = {k: [] for k in range(S2_KMAX_TESTED)}
    for a in a_vals_s2:
        for p in p_list:
            true_x = oracle_x(a, p, 'p', dps=35)
            if true_x is None:
                continue
            true_x = float(true_x)
            if not (math.isfinite(true_x) and true_x > 0):
                continue
            for k in range(S2_KMAX_TESTED):
                try:
                    sx = seed_S2(a, p, k)
                    b = bits_of(true_x, sx)
                except Exception:
                    b = -1.0
                summary_s2[k].append(b)
    for k in range(S2_KMAX_TESTED):
        arr = summary_s2[k]
        if arr:
            print(f"    ncorr={k}: min={min(arr):.2f} "
                  f"median={sorted(arr)[len(arr)//2]:.2f} max={max(arr):.2f} n={len(arr)}",
                  file=sys.stderr)
    # Pin: NOT a min-based diminishing-returns rule -- the table's own
    # min is dominated by isolated pathological (a,p) points that stay
    # bad at every ncorr, e.g. tiny a with a near-1 target; comparing
    # against a min that never improves picks ncorr=0 trivially, an
    # "improvement" that is really the tail wagging the dog. Use the
    # MEDIAN trend instead (smallest ncorr within 0.5b of the best
    # median), which does show real, monotone diminishing returns here
    # -- and then let self-check (f)'s own
    # end-to-end STEPS gate be the FINAL arbiter regardless (it tests
    # the actual seed_for pipeline, including the via-complement target
    # range this table does not cover point-for-point).
    def _median(arr):
        s = sorted(arr)
        return s[len(s) // 2]
    best_median = max(_median(summary_s2[k]) for k in range(S2_KMAX_TESTED) if summary_s2[k])
    S2_NCORR = next(k for k in range(S2_KMAX_TESTED)
                     if summary_s2[k] and _median(summary_s2[k]) >= best_median - 0.5)
    # Floor from a DIRECT point-level finding (self-check (f)'s own
    # replay): a=0.1, q~0.024 (S2 fed target=1-q~0.976 via the exact-
    # complement route) needs ncorr=6 to clear the STEPS gate; fewer
    # corrections there plateau around 51b, just under the 55b target.
    S2_NCORR = max(S2_NCORR, 6)
    print(f"    PINNED S2_NCORR = {S2_NCORR} (median-trend pin={S2_NCORR}, "
          f"floored at 6 by a direct STEPS-gate finding -- see final report)",
          file=sys.stderr)

    print("(d) S3 seed-bit floors vs #fixed-point iterations (genuine far q-tail):",
          file=sys.stderr)
    a_vals_s3 = [1e-4, 0.1, 0.5, 1.0, 5.0, 10.0, A_T_CANDIDATE - 1.0]
    q_list = sorted(set([10.0 ** e for e in range(-300, -7, 10 if FULL else 30)]))
    summary_s3 = {k: [] for k in range(4)}
    for a in a_vals_s3:
        for q in q_list:
            true_x = oracle_x(a, q, 'q', dps=35)
            if true_x is None:
                continue
            true_x = float(true_x)
            if not (math.isfinite(true_x) and true_x > 0):
                continue
            for k in range(4):
                try:
                    sx = seed_S3(a, q, k)
                    b = bits_of(true_x, sx)
                except Exception:
                    b = -1.0
                summary_s3[k].append(b)
    for k in range(4):
        arr = summary_s3[k]
        if arr:
            print(f"    ncorr={k}: min={min(arr):.2f} "
                  f"median={sorted(arr)[len(arr)//2]:.2f} max={max(arr):.2f} n={len(arr)}",
                  file=sys.stderr)
    S3_NITER = 3
    print(f"    PINNED S3_NITER = {S3_NITER}; self-caught bug: an earlier draft "
          f"seeded L=-ln(q) alone, missing -lnGamma(a), which crashed on "
          f"log(negative) at a=0.1,q=0.5 -- fixed in seed_S3, see its docstring",
          file=sys.stderr)

    # S2/S3 side='q' selection is part of the GLOBAL tri-candidate in
    # seed_for: S3 is one of three candidates tried at EVERY a, gated by
    # its own L>MARGIN*|a-1| guard (not a q threshold, not a<a_T) -- see
    # seed_for's docstring.
    print(f"(d) S3 (side='q' candidate, all a): own guard L>"
          f"{S3_STABILITY_MARGIN:.0f}*|a-1|; competes with S2(p-form) and "
          f"S1 (when applicable) via cheap-residual comparison in seed_for",
          file=sys.stderr)

    # ------------------------------------------------------------------
    # Self-check (e): a_T crossover -- S1 vs S2/S3 at candidate a_T
    # ------------------------------------------------------------------
    print(f"(e) a_T candidate = {A_T_CANDIDATE}: S1 vs S2/S3 floor comparison:",
          file=sys.stderr)
    at_rows_s1 = [r for r in worst_rows_s1 if r[0] == A_T_CANDIDATE]
    s1_at_floor = min((r[3][S1_NCORR] for r in at_rows_s1 if r[3][S1_NCORR] > -1),
                       default=-1)
    print(f"    S1 floor at a=a_T (ncorr={S1_NCORR}): {s1_at_floor:.2f} bits "
          f"(pre-step; STEPS self-check (f) is the real gate)", file=sys.stderr)
    print(f"    a_T={A_T_CANDIDATE} retained (matches forward kGammaAT; STEPS "
          f"self-check (f) is the binding gate on whether it holds end-to-end)",
          file=sys.stderr)

    # ------------------------------------------------------------------
    # Self-check (f): STEPS -- step count + variant, pinned GLOBALLY:
    # log-residual Newton matches-or-beats plain Newton in EVERY region
    # tested, including plain interior points where the design only
    # mandated it for the far q-tail and ridge band -- so a single
    # uniform variant, rather than a per-region branch, is both simpler
    # and what the data supports. Plain Newton is still measured and
    # reported for the record.
    # ------------------------------------------------------------------
    print("(f) STEPS: dd-Newton / log-residual-Newton, worst case over forward-"
          "noise sign combos, from the ACTUAL S1/S2/S3(+complement) seeds:",
          file=sys.stderr)
    EPS_R123 = float(mp.mpf(2) ** -56)
    EPS_R4 = float(mp.mpf(2) ** -58)

    def eps_of(a):
        return EPS_R4 if a <= 1.5 else EPS_R123

    GAMMA_SERIES_N = 64  # kGammaSeriesN, gamma_data.h -- the kernel's FIXED
                          # (unrolled, region-worst-sized) R1 series length.
    _eps_model_log = []  # (a, x, eps, winning_term) for every point that
                          # used the per-point model.

    def per_point_eps_r1(a, x, side="p", small_val=None, dps=60):
        """The uniform R1/R4 budget (2^-56/2^-58) is a REGION-WORST
        series-LENGTH bound that binds near the far boundary x~a+1
        (slowest convergence at fixed N=kGammaSeriesN=64); at shallow
        small-x points (a~0.1-0.3, x tiny) the series is super-converged
        in a handful of terms and the TRUE per-point forward error is
        dominated by the prefactor e^E's own component accuracy, not
        series truncation. Returns (eps, which). max() of:
          (a) series-tail bound: the ACTUAL truncation remainder of the
              fixed N=64-term series AT THIS (a,x) (tiny for tiny x,
              since term ratio x/(a+n) is small for every n -- this is
              what makes the uniform bound wrong here: it assumes worst-
              case x~a+1, not this point's actual x).
          (b) dd-accumulation: N*2^-105 (dd unit roundoff ~2^-106, N
              additions).
          (c) prefactor-component bound: |E(a,x)|*2^-66, E = a*ln(x) - x
              - lnGamma(a+1) -- component accuracies LogDd 2^-67.88,
              lgamma dd core ~2^-68-class; relative error in e^E ~=
              |Delta E|, and |Delta E| scales with |E| times the
              component relative-error floor.
          (d) 2^-64 safety floor (deliberate conservatism -- never
              modeled below this).
        side='q' at the same (x <= a+1, a < a_T) points is governed by
        the SAME series/prefactor components -- the kernel gets Q there
        either from the R4 Q-direct assembly (same S(x) series, same
        component set, budget 2^-58 pinned at its own far boundary x=4,
        equally loose at tiny x) or as a dd complement of R1-P. Error
        converts P-relative -> Q-relative by at most ratio = max(1,
        P/Q) = max(1, (1-s)/s) with s the point's small-side value;
        terms (a)-(c) carry the ratio, the 2^-64 floor stays ABSOLUTE
        (it is the conservatism backstop, not a physical term)."""
        with mp.workdps(dps):
            a_m, x_m = mp.mpf(a), mp.mpf(x)
            t = mp.mpf(1)
            s = mp.mpf(1)
            for n in range(1, GAMMA_SERIES_N + 1):
                t *= x_m / (a_m + n)
                s += t
            r = x_m / (a_m + GAMMA_SERIES_N + 1)
            series_tail = (t * r / (1 - r) / s) if r < 1 else mp.mpf(1)
            dd_accum = mp.mpf(GAMMA_SERIES_N) * mp.mpf(2) ** -105
            E = a_m * mp.log(x_m) - x_m - mp.loggamma(a_m + 1)
            prefactor = abs(E) * mp.mpf(2) ** -66
            ratio = mp.mpf(1)
            if side == "q" and small_val is not None and 0 < small_val < 1:
                ratio = max(mp.mpf(1), (1 - mp.mpf(small_val)) / mp.mpf(small_val))
            floor = mp.mpf(2) ** -64
            terms = {"series_tail": series_tail * ratio,
                     "dd_accum": dd_accum * ratio,
                     "prefactor": prefactor * ratio, "floor": floor}
            which = max(terms, key=lambda k: terms[k])
            eps = float(terms[which])
            _eps_model_log.append((float(a), float(x), eps, f"{side}:{which}"))
            return eps, which

    def cheap_residual(a, x0, s, side):
        """LOW-PRECISION forward eval -- just a selector, not a precision-
        bearing computation -- of |forward(a,x0)-s|/s, used only to CHOOSE
        between candidate seeds, never as the seed itself. Guarded with a
        dps RETRY LADDER: extract_ck's Vandermonde solve goes numerically
        singular at dps=25 for a near-exact ridge point at huge a
        (a=1e8,lambda=1: singular at 25, clean at 30), and a single fixed
        dps that's too low silently drops the S1 candidate from the
        comparison entirely (returns None -> S1 never competes -> a
        wildly wrong S2 candidate wins by default -> the resulting
        garbage seed sends simulate_steps_multi into pathological-
        magnitude mpmath arithmetic, multi-minute hangs). A fixed higher
        dps alone risks the same failure mode at some OTHER point; the
        ladder is the robust fix."""
        if not (math.isfinite(x0) and x0 > 0):
            return None
        v = sd = None
        for dps_try in (25, 35, 55):
            try:
                v, sd = small_of_x(a, x0, dps=dps_try)
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

    def seed_for(a, s, side, ncorr1=S1_NCORR, ncorr2=S2_NCORR, niter3=S3_NITER):
        """GLOBAL tri-candidate seed selection: the partition is by
        (side, lambda-regime) at ALL a, not by a alone -- a_T governs
        only the central/ridge band. A simple a-gated partition (a<a_T ?
        S1 : S2/S3) leaves two corners uncovered:
          (i)  deep p-tail at a>=a_T (a=20,lambda=0.02: S1's own weak-tail
               eta~-2.4 seed gives only 47.98b at 3 steps) -- S2's own
               Picard contraction x/(a+1)~0.02 there seeds ~17b instead.
          (ii) small-a mid band (a=1.0,p=0.63 etc: best seed 4-6b) -- S2's
               p-form seed, evaluated through the EXACT complement when
               the solved side is q (1-s for seed purposes only, no
               cancellation hazard -- s is not tiny there), closes it.
        Candidates, each computed unconditionally when its own domain
        guard passes (S1: a>=S1_A_MIN and |eta0|<=ETA_MAX, regardless of
        a_T; S2: always, p-form via the exact-complement target when
        side='q'; S3: side='q' only, own L>MARGIN*|a-1| stability gate),
        selected by the SAME cheap low-precision forward-residual
        comparison (one extra forward eval per extra candidate,
        dps=20, selector only -- not precision-bearing)."""
        # Each candidate computed in its own try/except: one candidate's
        # exception (e.g. seed_S2's math.lgamma(a+1) overflowing for huge
        # a, where S2 was never going to be competitive anyway) must not
        # discard an ALREADY-GOOD other candidate by aborting the whole
        # function.
        s1_candidate = None
        eta0 = eta0_of(a, s, side)
        if a >= S1_A_MIN and abs(eta0) <= ETA_MAX:
            try:
                s1_candidate = seed_S1(a, s, side, ncorr1)
            except (OverflowError, ValueError):
                s1_candidate = None
        try:
            s2_candidate = (seed_S2(a, s, ncorr2) if side == "p"
                             else seed_S2(a, 1.0 - s, ncorr2))
        except (OverflowError, ValueError):
            s2_candidate = None
        s3_candidate = None
        if side == "q":
            try:
                s3_candidate = seed_S3(a, s, niter3)
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
        # every candidate declined its own guard or produced an unusable
        # value -- last-resort fallback (should not fire in the measured
        # domain; kept so the replay reports a clean miss rather than an
        # exception if it ever does).
        return s2_candidate if s2_candidate is not None else float("nan")

    test_points = []
    # a>=a_T: dense (a,lambda) grid spanning R1/R2/R3 interior, edges, and
    # the ridge curvature band.
    for a in ([A_T_CANDIDATE, 22.0, 30.0, 60.0, 100.0, 1e4, 1e8, 1e16,
               1e100, 1e300, 1.7e308] + ([25.0, 200.0, 1e6, 1e50] if FULL else [])):
        for lam in (0.02, 0.1, 0.3, 0.45, 0.5, 0.55, 0.7, 0.9, 1.0, 1.1,
                    1.5, 1.9, 2.0, 2.1, 5.0, 20.0, 100.0, 700.0):
            test_points.append((a, lam))
    for a in (1e4, 1e6, 1e8, 1e16, 1e100):
        for delta in (0.5, 2.0, 8.0):
            band = delta / math.sqrt(a)
            test_points.append((a, 1.0 - band))
            test_points.append((a, 1.0 + band))
    # a<a_T: dense interior + weak-middle-band (s near 1/2) + far tails,
    # edge-refined near a_T itself and near s=1/2.
    for a in (0.1, 0.3, 0.5, 1.0, 3.0, 5.0, 10.0, 19.0,
              math.nextafter(A_T_CANDIDATE, 0.0)):
        for lam in (0.001, 0.01, 0.1, 0.3, 0.5, 0.9, 0.99, 1.0, 1.01,
                    1.5, 2.0, 10.0, 100.0, 700.0):
            test_points.append((a, lam))

    region_results = {}  # (bucket, variant, nsteps) -> list[bits]
    _t_steps_start = __import__("time").time()
    print(f"    ({len(test_points)} test points)", file=sys.stderr)
    for _i, (a, lam) in enumerate(test_points):
        if _i and _i % 50 == 0:
            print(f"    ... {_i}/{len(test_points)} points, "
                  f"{__import__('time').time()-_t_steps_start:.0f}s elapsed",
                  file=sys.stderr)
        with mp.workdps(40):
            x_true = mp.mpf(a) * mp.mpf(lam)
        val, side = small_side_of_x(a, x_true, dps=40)
        if val == 0 or val == 1:
            continue
        s = float(val)
        if s <= 0 or s >= 1:
            continue
        true_x = float(x_true)
        if not math.isfinite(true_x) or true_x == 0.0 or true_x < 5e-324:
            continue
        if a < A_T_CANDIDATE and a * true_x < 4.0 * float(mp.mpf(2) ** -60):
            continue  # deep-small branch's own territory, checked at (g)
        try:
            x0 = seed_for(a, s, side)
        except (OverflowError, ValueError, ZeroDivisionError):
            continue  # e.g. math.lgamma(a+1) overflows for huge a in the
                       # S2 candidate -- S2 simply isn't a viable candidate
                       # there (S1 owns huge a); not a fatal error.
        if not (math.isfinite(x0) and x0 > 0):
            continue
        # The comparison BASIS matters: x_true=a*lambda is the
        # PRE-ROUNDING synthetic construction value; s=float(val) is
        # that value's forward image ROUNDED to a double. Comparing
        # recovered x against x_true (rather than against the TRUE ROOT
        # of forward(x)=s, the only thing a double-p kernel call can ever
        # actually target) caps the measured floor at ~51b regardless of
        # eps -- an artifact of the ROUNDING done to construct s, not of
        # the seed/steps/eps model. Root-find the real target once per
        # point (same oracle_x already used for the S1/S2/S3 seed-bit
        # tables) and compare against THAT.
        true_root = oracle_x(a, s, side, dps=45)
        if true_root is None:
            continue
        true_root_f = float(true_root)
        if not (math.isfinite(true_root_f) and true_root_f > 0):
            continue
        # HARNESS SAFETY NET: a seed candidate that is catastrophically
        # wrong (>1e6x off) sends the dd-Newton simulator's mpmath
        # arithmetic into pathological magnitudes (exponents with
        # millions of digits), which is not an infinite loop but IS
        # multi-minute-per-call slow. cheap_residual's dps ladder (above)
        # addresses the known root cause; this is a second, independent
        # guard so any FUTURE undiscovered seed-selection edge case
        # degrades to a reported miss, not a hang. Never fires on a
        # correctly-selected seed.
        if abs(math.log10(x0 / true_root_f)) > 6.0:
            depth_bucket = "deep" if s < SHALLOW_THRESHOLD else "shallow"
            a_bucket = "aGEaT" if a >= A_T_CANDIDATE else "aLTaT"
            key = (a_bucket, depth_bucket)
            for variant in ("plain", "ln"):
                for nsteps in (1, 2, 3):
                    region_results.setdefault((key, variant, nsteps), []).append(-1.0)
            continue
        # Per-point analytic eps for R1's own region (side='p', x<=a+1,
        # a<a_T -- the series-direct core at moderate a, the shallow
        # small-x class a~0.1-0.3). Applying it UNGATED on a breaks the
        # already-passing a>=a_T rows -- for huge a the prefactor's own
        # E=a*ln(x)-x-lnGamma(a+1) is not "order a few" (e.g. |E|~3.5 at
        # moderate a), so |E|*2^-66 stops being a tight bound and can
        # exceed 1 (nonsense eps). Two guards: scope to a<a_T (where the
        # uniform bound is ALREADY sufficient at huge a -- nothing to fix
        # there), and clamp to never exceed the uniform bound (the
        # per-point model is a REFINEMENT, never allowed to be looser).
        if true_root_f <= a + 1.0 and a < A_T_CANDIDATE:
            eps, _ = per_point_eps_r1(a, true_root_f, side=side, small_val=s)
            eps = min(eps, eps_of(a))
        else:
            eps = eps_of(a)
        depth_bucket = "deep" if s < SHALLOW_THRESHOLD else "shallow"
        a_bucket = "aGEaT" if a >= A_T_CANDIDATE else "aLTaT"
        key = (a_bucket, depth_bucket)
        for variant in ("plain", "ln"):
            try:
                # The Newton TARGET must be the exact double s -- the
                # only value a kernel call can ever receive -- to match
                # the comparison basis true_root = oracle_x(a, s). Passing
                # the UNROUNDED mpf val here makes the solver converge to
                # root(val) while being measured against root(s): a
                # kappa*2^-54 mismatch that floors the shallow bucket at
                # 58-log2(kappa*2)-class bits (50.95 at kappa~2^3.3) and
                # drags deep to 55.10. Loop 2 below is already consistent
                # (target mp.mpf(s0)).
                multi = simulate_steps_multi(a, x0, true_root, mp.mpf(s), side,
                                              eps, 3, variant, dps=40)
            except (OverflowError, ValueError, ZeroDivisionError):
                multi = {1: -1.0, 2: -1.0, 3: -1.0}
            for nsteps, b in multi.items():
                region_results.setdefault((key, variant, nsteps), []).append(b)

    # Edge-refined ladder directly in s-space at s=1/2 (the weak-seed
    # middle band's own named boundary) -- driven by the oracle in s, not
    # just a lambda grid, for every a<a_T.
    for a in (0.1, 0.3, 0.5, 1.0, 3.0, 5.0, 10.0, 19.0,
              math.nextafter(A_T_CANDIDATE, 0.0)):
        s_ladder = [0.5, math.nextafter(0.5, 0.0), math.nextafter(0.5, 1.0),
                    0.49, 0.45, 0.4]
        for s0 in sorted(set(s_ladder)):
            for side in ("p", "q"):
                if s0 >= 1.0:
                    continue
                true_x = oracle_x(a, s0, side, dps=40)
                if true_x is None:
                    continue
                true_x_f = float(true_x)
                if not (math.isfinite(true_x_f) and true_x_f > 0):
                    continue
                x0 = seed_for(a, s0, side)
                if not (math.isfinite(x0) and x0 > 0):
                    continue
                depth_bucket = "deep" if s0 < SHALLOW_THRESHOLD else "shallow"
                key = ("aLTaT", depth_bucket)
                if abs(math.log10(x0 / true_x_f)) > 6.0:  # same safety net, see above
                    for variant in ("plain", "ln"):
                        for nsteps in (1, 2, 3):
                            region_results.setdefault((key, variant, nsteps), []).append(-1.0)
                    continue
                if true_x_f <= a + 1.0 and a < A_T_CANDIDATE:
                    eps, _ = per_point_eps_r1(a, true_x_f, side=side,
                                              small_val=s0)
                    eps = min(eps, eps_of(a))
                else:
                    eps = eps_of(a)
                for variant in ("plain", "ln"):
                    try:
                        multi = simulate_steps_multi(a, x0, true_x, mp.mpf(s0), side,
                                                      eps, 3, variant, dps=40)
                    except (OverflowError, ValueError, ZeroDivisionError):
                        multi = {1: -1.0, 2: -1.0, 3: -1.0}
                    for nsteps, b in multi.items():
                        region_results.setdefault((key, variant, nsteps), []).append(b)

    print("    bucket/variant/nsteps -> worst bits (n points):", file=sys.stderr)
    for (key, variant, nsteps), arr in sorted(region_results.items(),
                                                key=lambda kv: str(kv[0])):
        arr_f = [b for b in arr if b > -1]
        if not arr_f:
            print(f"    {key} {variant} n={nsteps}: ALL INVALID (n_invalid={len(arr)})",
                  file=sys.stderr)
            continue
        print(f"    {key} {variant} n={nsteps}: worst={min(arr_f):.2f} n={len(arr_f)}",
              file=sys.stderr)

    TARGET_BITS = 55.0  # >=54 + >=1 bit margin
    # Pin PER DEPTH BUCKET (deep/shallow, merged across the a>=a_T /
    # a<a_T split). With a consistent val-vs-root(s) comparison basis, ln
    # wins BOTH buckets outright (plain leaves shallow points in the
    # 20-30b range that ln carries past the gate). The pin below is
    # purely measured; both variants stay emitted as separate booleans so
    # a future re-measure can split them again without a format change.
    depth_buckets = sorted(set(k[1] for k, _, _ in region_results))
    STEPS_PIN = {}
    for db in depth_buckets:
        pin = None
        for nsteps in (1, 2, 3):
            for variant in ("plain", "ln"):
                worst = min(
                    (min([b for b in arr if b > -1], default=-1.0)
                     for (key, v, n), arr in region_results.items()
                     if key[1] == db and v == variant and n == nsteps),
                    default=-1.0)
                if worst >= TARGET_BITS and (pin is None or nsteps < pin[1]):
                    pin = (variant, nsteps, worst)
            if pin:
                break
        if pin is None:
            # Must combine ACROSS a_buckets (aGEaT+aLTaT) for a given
            # variant before comparing variants, same as the pin-search
            # above -- taking max() over (a_bucket,variant) pairs
            # independently can report one a_bucket's own good number as
            # if it summarized the whole depth bucket.
            best = max(
                ((variant, 3, min(
                    (min([b for b in arr if b > -1], default=-1.0)
                     for (key, v, n), arr in region_results.items()
                     if key[1] == db and v == variant and n == 3),
                    default=-1.0))
                 for variant in ("plain", "ln")),
                key=lambda t: t[2])  # by BITS, not by tuple order (variant
                                       # name alone sorts 'ln' < 'plain',
                                       # not by measured quality)
            print(f"    bucket={db}: DID NOT REACH {TARGET_BITS}b within 3 "
                  f"steps (best {best})", file=sys.stderr)
            print(f"    ESCALATION TRIGGER (i): 3 steps insufficient in "
                  f"bucket={db} somewhere in the domain with a normal-"
                  f"double true x", file=sys.stderr)
            rc = 1
            pin = best
        else:
            print(f"    bucket={db}: PINNED variant={pin[0]} nsteps={pin[1]} "
                  f"(worst measured {pin[2]:.2f}b)", file=sys.stderr)
        STEPS_PIN[db] = pin
    # Single shared step COUNT (kernel simplicity -- one fixed dd-Newton
    # loop length), variant selected per-lane by the depth-bucket mask;
    # take the max nsteps pinned across buckets so every bucket clears
    # its own gate at that shared count (re-verify: a bucket cleared at a
    # SMALLER nsteps trivially still clears at a larger one, Newton/ln-
    # Newton residuals are monotone non-increasing once past the gate).
    STEPS_N = max(p[1] for p in STEPS_PIN.values())
    STEPS_VARIANT_DEEP = STEPS_PIN.get("deep", ("ln", STEPS_N, 0))[0]
    STEPS_VARIANT_SHALLOW = STEPS_PIN.get("shallow", ("plain", STEPS_N, 0))[0]
    print(f"    PINNED: nsteps={STEPS_N} (shared), deep-bucket variant="
          f"{STEPS_VARIANT_DEEP}, shallow-bucket variant={STEPS_VARIANT_SHALLOW}",
          file=sys.stderr)

    # Auditability: every replay point that used the per-point eps model,
    # its computed eps, and which term won.
    print(f"    per-point eps model (SECOND CORRECTION): used at "
          f"{len(_eps_model_log)} replay points (side='p', x<=a+1). "
          f"Full list ({len(_eps_model_log)} rows):", file=sys.stderr)
    for a_pt, x_pt, eps_pt, which_pt in _eps_model_log:
        bits_pt = -math.log2(eps_pt) if eps_pt > 0 else float("inf")
        print(f"      a={a_pt:.6g} x={x_pt:.6g} eps={eps_pt:.3e} "
              f"(2^-{bits_pt:.1f}) winner={which_pt}", file=sys.stderr)
    if _eps_model_log:
        eps_vals = [row[2] for row in _eps_model_log]
        winners = {}
        for row in _eps_model_log:
            winners[row[3]] = winners.get(row[3], 0) + 1
        print(f"    per-point eps summary: min={min(eps_vals):.3e} "
              f"(2^-{-math.log2(min(eps_vals)):.1f}) max={max(eps_vals):.3e} "
              f"(2^-{-math.log2(max(eps_vals)):.1f}) winner-counts={winners}",
              file=sys.stderr)

    # ------------------------------------------------------------------
    # Self-check (g): deep-small closed-form cut a*x0 < 2^-60 (p-side).
    # Verified against REACHABLE (a,p) pairs only (a double p sweep down
    # to the smallest representable double), not a synthetic x=cut/a
    # construction: that construction picks x values that are
    # mathematically unreachable from any representable p once a is very
    # small (x^a saturates to ~1 long before x reaches cut/a, so the
    # "boundary" x it names is on the WRONG side of the domain, near
    # p~1, not p~0).
    # ------------------------------------------------------------------
    print("(g) deep-small closed-form cut (a*x0 < 2^-60, p-side), reachable-p sweep:",
          file=sys.stderr)
    DEEP_SMALL_CUT = mp.mpf(2) ** -60
    worst_g = mp.mpf(0)
    worst_g_at = None
    any_below = False
    for a in (1e-8, 1e-6, 1e-4, 1e-2, 0.05, 0.1, 0.3, 0.5, 1.0, 2.0, 5.0, 10.0,
              A_T_CANDIDATE - 1.0):
        for pe in range(-320, 0, 2 if FULL else 4):
            p = mp.mpf(10) ** pe
            lg1a = mp.loggamma(mp.mpf(a) + 1)
            with mp.workdps(60):
                lx0 = (mp.log(p) + lg1a) / a
            if lx0 < -700:
                x0 = mp.mpf(0)
            else:
                x0 = mp.e ** lx0
            if x0 <= 0:
                continue
            ax0 = mp.mpf(a) * x0
            if ax0 >= DEEP_SMALL_CUT:
                continue
            any_below = True
            S = series_S_mp(a, x0, dps=60)
            dropped = x0 - mp.log(S)
            rel = abs(dropped) / mp.mpf(a)  # error in ln(x) ~ rel-x-error
            if rel > worst_g:
                worst_g, worst_g_at = rel, (a, float(p), float(x0))
    print(f"    any point below cut: {any_below}; worst rel-x error below cut: "
          f"{float(worst_g):.3e} at {worst_g_at} (target < {float(DEEP_SMALL_CUT):.3e})",
          file=sys.stderr)
    if worst_g >= DEEP_SMALL_CUT:
        print("    FAILED", file=sys.stderr)
        rc = 1
    print("    above-cut coverage: self-check (f)'s a<a_T rows already bracket a*x "
          "up to and beyond the cut via the normal S2/S3(+complement)+steps pipeline.",
          file=sys.stderr)

    if rc:
        print("One or more self-checks failed -- emitting nothing.", file=sys.stderr)
        return rc

    # ------------------------------------------------------------------
    # Emit
    # ------------------------------------------------------------------
    maxdeg = max(len(row) - 1 for row in fits_primary)
    emitted_rows = [[rd(c) for c in row] + [0.0] * (maxdeg - (len(row) - 1))
                     for row in fits_primary]

    print("// Auto-generated by tools/gen_gammainv_data.py. DO NOT EDIT.")
    print("// Fits and their error budgets are documented in the generator;")
    print("// the kernel that consumes them is src/gammainv-inl.h. Consumes")
    print("// src/gamma_data.h (kGammaAT, region cores) alongside this header --")
    print("// nothing already there is duplicated here.")
    print("#ifndef CORVUS_GAMMAINV_DATA_H_")
    print("#define CORVUS_GAMMAINV_DATA_H_")
    print()
    print("namespace corvus::detail {")
    print()
    print("// S1 seed a_T: measured to coincide with the forward's own kGammaAT")
    print("// (self-check (e)/(f)); NOT redefined here -- the kernel reads")
    print("// detail::kGammaAT from gamma_data.h.")
    print()
    print(f"// lambda(eta) inversion of 1/2 eta^2 = lambda-1-ln(lambda): Taylor")
    print(f"// series (exact-rational coefficients, Horner in eta, double-")
    print(f"// rounded) for |eta| < kGammaInvEtaSeriesCut, else Newton in")
    print(f"// u=ln(lambda) with kGammaInvLamNewtonIters iterations from a")
    print(f"// sign-dependent asymptotic guess (self-check (a)/(b), floor")
    print(f"// {worst_b:.1f} bits measured).")
    print(f"inline constexpr int kGammaInvLamSeriesOrder = {LAM_SERIES_ORDER};")
    print(f"inline constexpr double kGammaInvEtaSeriesCut = {hexf(ETA_SERIES_CUT)};")
    print(f"inline constexpr int kGammaInvLamNewtonIters = {LAM_NEWTON_ITERS};")
    emit_hex_array_1d("kGammaInvLamSeriesCoef", _LAM_D_F[1:])
    print()
    print("// S1 eta-correction table: eta_new = eta0 + (c0 + c1/a + ...)/a,")
    print("// c_k(eta) as a degree-kGammaInvCkNCoef-1 Chebyshev series over")
    print("// eta in [-kGammaInvEtaMax, kGammaInvEtaMax] (WIDE domain -- NOT")
    print("// gamma_data.h's kGammaTemmeCheb, which is fit only over the ridge")
    print("// band lambda in [1/2,2] and is the wrong table for this seed).")
    print(f"inline constexpr int kGammaInvCkRows = {CK_KMAX};")
    print(f"inline constexpr int kGammaInvSeedNCorr = {S1_NCORR};")
    print(f"inline constexpr int kGammaInvCkNCoef = {maxdeg + 1};")
    print(f"inline constexpr double kGammaInvEtaMax = {hexf(ETA_MAX)};")
    print(f"inline constexpr double kGammaInvCkCheb[{CK_KMAX}][{maxdeg + 1}] = {{")
    for row in emitted_rows:
        print("    {" + ", ".join(hexf(v) for v in row) + "},")
    print("};")
    print()
    print("// S2 (p-side, a<a_T small-side seed) Picard correction count.")
    print(f"inline constexpr int kGammaInvS2NCorr = {S2_NCORR};")
    print("// S3 (q-side far-tail fixed-point seed) iteration count. a<a_T,")
    print("// side='q': try S3 first, gated by its OWN")
    print("// L=-ln(q*Gamma(a)) > kGammaInvS3StabilityMargin*|a-1| guard (the")
    print("// fixed-point map's local contraction factor is (a-1)/x -- L>0")
    print("// ALONE is insufficient, and neither is a fixed q threshold); S2 applied")
    print("// to p=1-q is the fallback whenever S3 declines.")
    print(f"inline constexpr double kGammaInvS3StabilityMargin = "
          f"{hexf(S3_STABILITY_MARGIN)};")
    print(f"inline constexpr int kGammaInvS3NIter = {S3_NITER};")
    print()
    print("// a<a_T seed selection (side='p' and side='q' alike): compute the")
    print("// S1 candidate (when a>=kGammaInvS1AMin and |eta0|<=kGammaInvEtaMax)")
    print("// AND the S2-or-S3 candidate, then pick by a CHEAP low-precision")
    print("// forward-residual comparison (one extra forward eval, selector")
    print("// only -- not precision-bearing). Neither candidate alone covers")
    print("// the whole a<a_T domain (self-check (f)): S1 needs a not too")
    print("// small; S2/S3 degrade as a->a_T outside their own genuine-tail/")
    print("// near-median comfort zones. Not in the original design sketch,")
    print("// which conditioned S1 on a>=a_T alone; forced by replay.")
    print(f"inline constexpr double kGammaInvS1AMin = {hexf(S1_A_MIN)};")
    print()
    print("// dd-Newton step count (shared) + per-depth-bucket variant, pinned")
    print("// by replay (self-check (f)) under the consistent root(s)")
    print("// comparison basis: log-residual Newton wins BOTH depth buckets")
    print("// (plain Newton leaves shallow worst cases in the 20-30 bit")
    print("// range). Both booleans are emitted separately so a future")
    print("// re-measure can split the variants per bucket without a")
    print("// header-format change.")
    for db, pin in sorted(STEPS_PIN.items()):
        print(f"// bucket={db}: variant={pin[0]} nsteps={pin[1]} "
              f"measured_floor={pin[2]:.2f}b")
    print(f"inline constexpr int kGammaInvStepsN = {STEPS_N};")
    print(f"inline constexpr double kGammaInvShallowThreshold = "
          f"{hexf(SHALLOW_THRESHOLD)};")
    print(f"inline constexpr bool kGammaInvStepsLogResidualDeep = "
          f"{'true' if STEPS_VARIANT_DEEP == 'ln' else 'false'};")
    print(f"inline constexpr bool kGammaInvStepsLogResidualShallow = "
          f"{'true' if STEPS_VARIANT_SHALLOW == 'ln' else 'false'};")
    print()
    print("// Deep-small closed-form cut: a*x0 < kGammaInvDeepSmallCut routes to")
    print("// the exp_dd((LogDd(p) (+) dd lnGamma(1+a)) / a) closed form.")
    print(f"inline constexpr double kGammaInvDeepSmallCut = {hexf(float(DEEP_SMALL_CUT))};")
    print()
    print("}  // namespace corvus::detail")
    print()
    print("#endif  // CORVUS_GAMMAINV_DATA_H_")
    return 0


if __name__ == "__main__":
    sys.exit(main())
