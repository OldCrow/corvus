"""INDEPENDENT verification harness for tests/data/beta_{p,q}_reference.txt.

Shares NO code with gen_beta_reference's assemblies (no import of the
oracle module). Methods:
  series     -- DLMF 8.17.8 in the small-side left-tail frame:
                  small = x'^a'(1-x')^b'/(a' B) * sum_k (a'+b')_k/(a'+1)_k x'^k
                ALL terms positive -- no cancellation is possible, so
                accuracy tracks working precision directly. Prefactor in
                log space. Converges at rate x' (used for x' <= 0.75).
  half-split -- for x' > 0.75 (incl. the x'->1 exact-complement family
                where the series cannot converge):
                  small = series(a',b',1/2)
                        + int_{ln d}^{-ln 2} e^{b'v}(1-e^v)^{a'-1} dv / B
                with v = ln(1-t) and d = 1-x' recovered EXACTLY (small-
                result direction of the complement; exact at any prec
                >= the operand's own mantissa). Log coordinates flatten
                the u^{b'-1} endpoint singularity that defeats tanh-sinh
                in u-space, and for huge a' the integrand localizes at
                v ~ ln d + O(1) where split points are placed. Both
                pieces positive -- no cancellation.
  LAYERING   -- every independent value is computed at dps 160 AND 220
                and must agree to 1e-30, else the row is a harness
                failure (never silently trusted).
  identity   -- exact closed forms on the analytic lines (a=1, b=1,
                a=b=1/2, a=b & x=1/2).
  betainc300 -- mpmath betainc at dps 300, only inside its trusted zone
                (moderate parameters, off-diagonal).
  log-bound  -- rigorous inequality certification for rows stored as
                exact {0,1}: ln(small) <= a'ln x' - ln a'
                + min(0,(b'-1)ln(1-x')) - lnB; if the bound is below
                half the min subnormal, the stored 0 is correctly
                rounded. No assembly, just a monotone bound.

Independence claim, stated honestly: the series shares its mathematical
family with the oracle's R1 region but none of its code, dps ladder, or
assembly; what this harness exists to catch is assembly/complement/
truncation defects, guarded here by brute precision + layering, not by
algorithmic novelty. Cross-method anchors: betainc300 self-cert on the
normal sample, and a NEGATIVE CONTROL run on every invocation -- four
known-bad rows must FAIL against their old values and PASS against
their corrected ones, else the harness refuses to certify anything.

ACCEPTANCE (#61): a stored double s passes only if it is the CORRECT
ROUNDING of the independent value v: v must lie strictly inside s's
round-to-nearest interval (midpoints to its two neighbours, exact in
mpf), shrunk on each side by the independent value's own uncertainty
INDEP_REL*|v| (the layering bound). That is |v - s| < ulp/2 - INDEP_REL*|v|,
i.e. ~2^-53 relative, but exact where a symmetric relative bound is
wrong: at a power of two the lower neighbour is half an ulp closer, and
subnormals have a fixed absolute ulp. v within INDEP_REL*|v| of a midpoint
is a NEAR-TIE: the harness cannot decide the rounding, so the row is
DECLINED -- never counted as certified, and every one is printed by
bits. It does not fail the run: such a stored value is still within
ulp/2*(1+INDEP_REL) of the truth, the error any correct rounding has
there, and exact midpoints do occur on the identity lines (e.g.
1-x^2 at x = 1-2^-53 is 2^-52 - 2^-106, the midpoint below 2^-52, so
either neighbour is a half-ulp rounding). The gates this backs
are pinned at 0-1 ULP, so a 2-ULP acceptance (the old 2^-51 + 1e-22)
could not see a 1-ULP reference error.
1-ULP NEGATIVE CONTROL: the four corrected NEG_CONTROL values nudged one
ulp either way (math.nextafter) must be rejected before any row is
judged; after the run, every certified stored row is re-judged with its
stored value nudged one ulp either way, and every nudge must be
rejected. Either failure is exit 2.
FRESH GAMMA-LIMIT SAMPLING (#54/#61): the gammalim stratum (max param >=
kBetaGammaLim, read from src/beta_data.h; min param below it) is sampled
per (s-decade, t/s side) cell so no populated sub-regime is starved, and
cells with no stored rows are NAMED. Independently of stored rows,
fresh seeded points across the slice's domain are evaluated by the
independent route and compared to mpmath betainc (dps 300) to 1e-25,
so the evaluator is certified where the reference set is thin.

Stage B invariants over the full file: P+Q=1 (<=1 ulp), monotonicity of
P in x within every (a,b) group.

Exit non-zero on any failure. Strata without coverage are NAMED."""
import math
import os
import random
import re
import sys
import mpmath
from mpmath import mp

SEED = 987654
N_PER_STRATUM = 60
MIN_SUBN_HALF = mp.mpf(2) ** -1075
# Relative uncertainty of an independent value: indep_small's dps-160
# and dps-220 layers must agree to 1e-30 (~2^-99.7), so the dps-220 value
# is trusted to that; identity values run at dps 80. The acceptance test
# shrinks the rounding interval by this much on each side.
INDEP_REL = mp.mpf("1e-30")
SELF_TOL = mp.mpf("1e-25")   # indep-vs-betainc self-certification
N_GL_FRESH = 12              # fresh gamma-limit points (plus anchors)
GL_FRESH_HUGE_MAX_J = 64     # fresh huge param in [2^47, 2^64): see below


def read_gamma_lim():
    """kBetaGammaLim from src/beta_data.h (a data read, not shared code):
    the stratum boundary must track the kernel's router, which a
    hardcoded copy did not (it said 2^59 after the pin moved to 2^47)."""
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                        "src", "beta_data.h")
    with open(path) as fh:
        m = re.search(r"kBetaGammaLim\s*=\s*([0-9a-fA-FxXpP.+-]+)\s*;",
                      fh.read())
    if not m:
        raise RuntimeError("kBetaGammaLim not found in src/beta_data.h")
    return float.fromhex(m.group(1))

# Negative control: four known-bad oracle values paired with their
# corrected values. Every invocation must REJECT the bad value and
# ACCEPT the good one, or the harness refuses to certify anything.
NEG_CONTROL = [
    (20.0, float.fromhex("0x1.0efe7e62615a0p-24"),
     float.fromhex("0x1.9999999999994p-2"),
     "5.6079293314231234e-17", "5.6079293424256988e-17"),
    (20.0, float.fromhex("0x1.b7cdfd9d7bdbbp-34"),
     float.fromhex("0x1.9999999999994p-2"),
     "8.8879673214171994e-20", "8.8879673214172437e-20"),
    (float.fromhex("0x1.b48888fc5f595p+4"),
     float.fromhex("0x1.b9d9f72e264a7p-26"), 0.65,
     "1.9940066017325325e-14", "1.9940066017329276e-14"),
    (float.fromhex("0x1.fb2a73489786bp+3"),
     float.fromhex("0x1.5798ee2308c3ap-27"),
     float.fromhex("0x1.87ae147ae147bp-2"),
     "2.3949587672500623e-16", "2.3949587672571292e-16"),
]


def one_minus_exact(x):
    """Exact 1-x for a double-derived x in (0,1)."""
    xm = mp.mpf(x)
    if not (0 < xm < 1):
        return 1 - xm
    with mpmath.workdps(mp.dps + max(0, int(-mp.log10(xm))) + 20):
        return 1 - xm


def ln_beta(a, b):
    """ln B(a,b), forming a+b at enough digits that nothing truncates."""
    am, bm = mp.mpf(a), mp.mpf(b)
    lo, hi = (am, bm) if am <= bm else (bm, am)
    extra = max(0, int(mp.log10(hi / lo))) + 20
    with mpmath.workdps(mp.dps + extra):
        return mp.loggamma(lo) + mp.loggamma(hi) - mp.loggamma(hi + lo)


def small_frame(a, b, x, P, Q):
    """(a', b', x', d', stored_small): frame where the small side is the
    left tail int_0^{x'} t^{a'-1}(1-t)^{b'-1} dt. d' = 1-x' is CARRIED,
    never recomputed downstream: recovering it from x' via subtraction
    re-rounds x' to working precision first (mp.mpf() on an mpf rounds
    to ambient dps), which collapses d' to 0 for x' = 1-2e-216 at dps
    200 -- a complement-collapse hazard."""
    if P <= Q:
        return mp.mpf(a), mp.mpf(b), mp.mpf(x), one_minus_exact(x), P
    return mp.mpf(b), mp.mpf(a), one_minus_exact(x), mp.mpf(x), Q


SERIES_CAP = 400000


def series_regularized(ap, bp, xp, dps):
    """I_x'(a',b') via the all-positive-terms series (DLMF 8.17.8):
      x'^a'(1-x')^b' / (a' B(a',b')) * sum_k (a'+b')_k/(a'+1)_k x'^k
    Raises if the cap is hit before convergence -- never a silent
    partial sum, which would silently understate the tail probability."""
    with mpmath.workdps(dps):
        am, bm, xm = mp.mpf(ap), mp.mpf(bp), mp.mpf(xp)
        # (1-x)^b in log form via log1p, NEVER via a formed complement:
        # forming the complement (omx = 1 - xm at dps+40) truncates to
        # exactly 1 for x below 10^-(dps+40), silently DROPPING the whole
        # b*ln(1-x) term. When b*x reaches ~400 with x ~ 4e-305, this
        # makes the "independent" value come back e^400 too large (an
        # impossible P = 4.75e34), and the same truncation at the
        # 160-layer's 10^-200 horizon produces pure layer disagreements
        # at b = 1e200. log1p consumes xm directly: an original double
        # is exact at any dps >= 25, and the exact carried complement
        # only reaches this branch when <= 3/4, where log1p is
        # cancellation-free and its argument truncation is a relative
        # 10^-dps. (This follows the harness's own small_frame
        # discipline.)
        lnpre = (am * mp.log(xm) + bm * mp.log1p(-xm)
                 - mp.log(am) - ln_beta(am, bm))
        s = mp.mpf(1)
        term = mp.mpf(1)
        eps = mp.mpf(10) ** (-dps - 5)
        k = 0
        while True:
            k += 1
            term *= (am + bm + k - 1) / (am + k) * xm
            s += term
            # safe stop: past the term peak (ratio < 1) and negligible
            if term < s * eps and (am + bm + k) * xm < (am + k):
                break
            if k > SERIES_CAP:
                raise ArithmeticError(f"series cap {SERIES_CAP} hit "
                                      f"(x'={mp.nstr(xm, 8)})")
        return mp.exp(lnpre + mp.log(s))


def tail_piece_logquad(ap, bp, d, dps):
    """int_{1/2}^{1-d} t^(a'-1)(1-t)^(b'-1) dt / B  via v = ln(1-t):
    = int_{ln d}^{-ln 2} e^(b'v) (1-e^v)^(a'-1) dv / B.
    v <= -ln 2 => e^v <= 1/2 => log1p(-e^v) is cancellation-free.
    For huge a' the factor (1-e^v)^(a'-1) ~ e^(-a' e^v) kills the
    integrand for v > -ln a'; split points cover both localizations."""
    with mpmath.workdps(dps):
        am, bm = mp.mpf(ap), mp.mpf(bp)
        lnB = ln_beta(am, bm)
        lnd = mp.log(d)
        vhi = -mp.log(2)

        def lgv(v):
            return bm * v + (am - 1) * mp.log1p(-mp.exp(v))

        pts = [lnd]
        for step in (1, 2, 4, 8, 16, 32, 64, 128, 256, 512):
            c = lnd + step
            if c < vhi - mp.mpf("0.5"):
                pts.append(c)
        if am > 2:
            lna = mp.log(am)
            for off in (-2, 0, 2, 5):
                c = -lna + off
                if lnd < c < vhi - mp.mpf("0.5"):
                    pts.append(c)
        pts.append(vhi)
        pts = sorted(set(pts))
        # mp.quad's convergence logic degrades when the integral's
        # magnitude sits far below working epsilon (e.g. 1e-241 at
        # dps 160 comes back as sqrt(eps)-scale noise). Normalize to
        # O(1) via the peak log-magnitude over the split points, then
        # scale back -- the split points bracket every localization the
        # integrand can have (endpoint at ln d, ridge at -ln a').
        lmax = max(lgv(p) for p in pts)
        val = mp.quad(lambda v: mp.exp(lgv(v) - lmax), pts)
        return val * mp.exp(lmax - lnB)


def indep_small(ap, bp, xp, d):
    """Layered independent value of the small side: dps 160 and 220
    must agree to 1e-30 or the row is a harness failure. d = 1-x' is
    the exact carried complement from small_frame."""
    if not (d > 0) and xp < 1:
        raise ArithmeticError("carried complement d invalid")
    vals = []
    for dps in (160, 220):
        if xp <= mp.mpf("0.75"):
            vals.append(series_regularized(ap, bp, xp, dps))
        else:
            vals.append(series_regularized(ap, bp, mp.mpf("0.5"), dps)
                        + tail_piece_logquad(ap, bp, d, dps))
    if vals[1] == 0:
        if vals[0] == 0:
            return vals[1]
        raise ArithmeticError("layer disagreement (zero vs nonzero)")
    if abs((vals[0] - vals[1]) / vals[1]) > mp.mpf("1e-30"):
        raise ArithmeticError(
            f"layer disagreement {mp.nstr(vals[0], 20)} vs "
            f"{mp.nstr(vals[1], 20)}")
    return vals[1]


def identity_value(a, b, x):
    """Exact closed form if the point sits on an analytic line."""
    am, bm, xm = mp.mpf(a), mp.mpf(b), mp.mpf(x)
    # BOTH sides computed directly -- a "1 - P" here at ambient dps is
    # the exact complement-collapse this harness exists to catch.
    if am == 1:
        t = bm * mp.log1p(-xm)
        return -mp.expm1(t), mp.exp(t)
    if bm == 1:
        t = am * mp.log(xm)
        return mp.exp(t), -mp.expm1(t)
    if am == bm == mp.mpf("0.5"):
        r = mp.sqrt(xm)
        return (2 / mp.pi) * mp.asin(r), (2 / mp.pi) * mp.acos(r)
    if am == bm and xm == mp.mpf("0.5"):
        return mp.mpf("0.5"), mp.mpf("0.5")
    return None


def betainc300_ok(a, b, x):
    return (max(a, b) <= 1e6 and min(a, b) >= 1e-6
            and not (a == b and abs(x - 0.5) < 1e-6))


def log_bound_small(ap, bp, xp, d):
    """Rigorous upper bound on ln(small side) in the left-tail frame.
    d = 1-x' carried exactly from small_frame."""
    with mpmath.workdps(80):
        lb = ap * mp.log(xp) - mp.log(ap) - ln_beta(ap, bp)
        if bp < 1:
            lb += (bp - 1) * mp.log(d)
        return lb


def rounding_verdict(s, v, rel=INDEP_REL):
    """Judge stored double s against independent value v (mpf).
    'pass': v lies strictly inside s's round-to-nearest interval, shrunk
            by rel*|v| on each side -- s is certainly RN(v).
    'fail': v lies outside that interval widened by rel*|v| -- s is
            certainly not RN(v).
    'tie' : v is within rel*|v| of an interval edge (a rounding midpoint)
            -- undecidable at this precision; the caller DECLINES the row
            (not certified, printed, not a failure: |v - s| is still
            <= ulp/2*(1+rel)).
    The interval [(prev+s)/2, (s+next)/2] is exact in binary mpf at 200
    bits and is asymmetric at powers of two; for s = 0 or subnormal s it
    is the fixed half-subnormal-ulp interval."""
    s = float(s)
    with mpmath.workprec(200):
        sm = mp.mpf(s)
        lo = (mp.mpf(math.nextafter(s, -math.inf)) + sm) / 2
        hi = (sm + mp.mpf(math.nextafter(s, math.inf))) / 2
        v = mp.mpf(v)
        dv = rel * abs(v)
        if lo + dv < v < hi - dv:
            return "pass"
        if v < lo - dv or v > hi + dv:
            return "fail"
        return "tie"


def negative_control():
    """Detector validation: the evaluator must separate the four known-
    bad values from their corrected values (gross control), and must
    reject each corrected value nudged one ulp either way (1-ULP
    control) -- the gates are pinned at 0-1 ULP, so a detector that
    passes a 1-ULP error certifies nothing they rely on."""
    ok = True
    n_ulp = 0
    for a, b, x, bad, good in NEG_CONTROL:
        iv = indep_small(mp.mpf(a), mp.mpf(b), mp.mpf(x), one_minus_exact(x))
        where = f"a={a:.6g} b={b:.6g} x={x:.6g}"
        if rounding_verdict(float(bad), iv) != "fail":
            print(f"NEGATIVE CONTROL FAILED: bad value accepted at {where}",
                  file=sys.stderr)
            ok = False
        g = float(good)
        if rounding_verdict(g, iv) != "pass":
            print(f"NEGATIVE CONTROL FAILED: good value not certified at "
                  f"{where} (indep={mp.nstr(iv, 20)})", file=sys.stderr)
            ok = False
        for nudged in (math.nextafter(g, math.inf),
                       math.nextafter(g, -math.inf)):
            n_ulp += 1
            if rounding_verdict(nudged, iv) != "fail":
                print(f"1-ULP CONTROL FAILED: good value nudged to "
                      f"{nudged.hex()} not rejected at {where}",
                      file=sys.stderr)
                ok = False
    return ok, n_ulp


GL_S_BUCKETS = ("s<1", "1<=s<10", "10<=s<100", "100<=s<1e3", "s>=1e3")


def gl_cell(a, b, x, b_gl):
    """(s-bucket, t/s side) of a gamma-limit point: s is the parameter
    below kBetaGammaLim, t the gamma-limit argument, t = -huge*ln(1-x)
    for huge second and -huge*ln(x) for huge first (the unshifted form
    is enough to bin). Side 'lo' (t < s) and 'hi' (t > s) put the small
    side in opposite tails."""
    hf = a >= b_gl
    s, huge = (b, a) if hf else (a, b)
    t = -huge * (math.log(x) if hf else math.log1p(-x))
    k = (min(max(int(math.floor(math.log10(s))) + 1, 0), len(GL_S_BUCKETS) - 1)
         if s > 0 else 0)
    return GL_S_BUCKETS[k], ("lo" if t < s else "hi")


def gammalim_fresh(rng, b_gl, fails):
    """Fresh gamma-limit points, no stored rows needed: independent route
    vs mpmath betainc at dps 300, agreement to SELF_TOL. Domain is the
    slice's own (gen_beta_data.py GAMMA_LIM_SLICE_*): one parameter
    huge >= kBetaGammaLim, the other s in [0.05, 2500], t/s off the
    ratio band (r in [0.15, 0.45] or [2.05, 3.5]), both orientations,
    small side above 1e-300 (saturated points are skipped and redrawn).
    The huge parameter is capped at 2^64: mpmath betainc is known WRONG
    in this family at huge ~1e250 (gen_beta_reference.py's GAMMALIM
    GUARD: (1, 1e250, 8e-250), its dps layers agreeing on the garbage),
    so betainc is never trusted alone -- the certificate is two unrelated
    methods agreeing. Anchors: the #54 witness (2000, 2^59, ...) and a
    huge-first point. Returns (points checked, cells hit)."""
    pts = [(2000.0, 2.0 ** 59, 1.3877787807814447e-15),
           (2.0 ** 50 * 1.37, 0.25, 1.0 - 2.0 ** -48)]
    tries = 0
    while len(pts) < N_GL_FRESH + 2 and tries < 50 * N_GL_FRESH:
        tries += 1
        huge = float(2.0 ** rng.uniform(math.log2(b_gl), GL_FRESH_HUGE_MAX_J))
        s = float(10.0 ** rng.uniform(math.log10(0.05), math.log10(2500.0)))
        r = rng.uniform(0.15, 0.45) if rng.random() < 0.5 else rng.uniform(2.05, 3.5)
        t = r * s
        if rng.random() < 0.5:
            a, b, x = s, huge, -math.expm1(-t / huge)           # huge second
        else:
            # huge first: x = 1 - k*2^-53 is exact (x in [1/2, 1)); the
            # realized t moves by the rounding of k, so r is re-derived.
            k = max(1, round(t / huge * 2.0 ** 53))
            a, b, x = huge, s, 1.0 - k * 2.0 ** -53
            t = -huge * math.log(x)
            r = t / s
        # off the ratio band, and small side ~ s ln t - t - lgamma(s+1)
        # (crude, both tails) not below ~1e-295
        if 0.45 < r < 2.05 or s * math.log(t) - t - math.lgamma(s + 1) < -680:
            continue
        pts.append((a, b, x))
    n_done = 0
    cells = set()
    for a, b, x in pts:
        if not (0.0 < x < 1.0):
            continue
        hf = a >= b_gl
        s, huge = (b, a) if hf else (a, b)
        cell = gl_cell(a, b, x, b_gl)
        # Small side by tail: huge second -> P_gamma-like, small P when
        # t < s; huge first -> 1 - P_gamma-like, small P when t > s. The
        # mean is only a first guess (for s < 1 the gamma median sits far
        # below it), so a value > 1/2 flips to the other frame once.
        p_small = (cell[1] == "lo") != hf
        fr_p = (mp.mpf(a), mp.mpf(b), mp.mpf(x), one_minus_exact(x))
        fr_q = (mp.mpf(b), mp.mpf(a), one_minus_exact(x), mp.mpf(x))
        r = (a, b, x, float("nan"), float("nan"))
        try:
            fr = fr_p if p_small else fr_q
            iv = indep_small(*fr)
            if iv > mp.mpf("0.5"):
                fr = fr_q if p_small else fr_p
                iv = indep_small(*fr)
        except Exception as e:
            fails.append(("gammalim-fresh:indep-error", r, float("nan"), repr(e)))
            continue
        if iv > mp.mpf("0.5"):
            fails.append(("gammalim-fresh:ORIENTATION", r, float("nan"),
                          mp.nstr(iv, 17)))
            continue
        if iv < mp.mpf("1e-300"):
            continue
        with mpmath.workdps(300):
            bv = mp.betainc(fr[0], fr[1], 0, fr[2], regularized=True)
        if bv == 0 or abs((iv - bv) / bv) > SELF_TOL:
            fails.append(("gammalim-fresh:DISAGREE", r, mp.nstr(iv, 20),
                          mp.nstr(bv, 20)))
        n_done += 1
        cells.add(cell)
    return n_done, cells


def main():
    mp.dps = 40
    b_gl = read_gamma_lim()
    nc_ok, nc_ulp = negative_control()
    if not nc_ok:
        print("harness cannot certify itself -- aborting before "
              "judging any reference row.", file=sys.stderr)
        return 2
    print(f"negative control passed (4 known-bad rows rejected, "
          f"4 corrected values certified as correctly rounded, "
          f"{nc_ulp} one-ulp nudges rejected)", file=sys.stderr)
    rows = []
    beta_p_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "tests", "data", "beta_p_reference.txt")
    with open(beta_p_path) as fh:
        for ln in fh:
            t = ln.split()
            if len(t) == 5:
                rows.append(tuple(float.fromhex(v) for v in t))
    print(f"rows loaded: {len(rows)}", file=sys.stderr)

    # ---- Stage B invariants (full file) ----
    bad_pq = 0
    groups = {}
    for a, b, x, P, Q in rows:
        s = P + Q
        if abs(1 - s) > 2.3e-16:
            bad_pq += 1
        groups.setdefault((a, b), []).append((x, P))
    bad_mono = []
    for (a, b), pts in groups.items():
        if len(pts) < 3:
            continue
        pts.sort()
        for i in range(1, len(pts)):
            if pts[i][1] < pts[i - 1][1]:
                bad_mono.append((a, b, pts[i - 1], pts[i]))
    print(f"invariants: P+Q>1ulp rows={bad_pq}; "
          f"monotonicity violations={len(bad_mono)} "
          f"over {sum(1 for g in groups.values() if len(g) >= 3)} groups",
          file=sys.stderr)
    for v in bad_mono[:10]:
        print(f"    MONO: a={v[0]:.6g} b={v[1]:.6g} "
              f"x={v[2][0]:.9g}->P={v[2][1]:.9g} then "
              f"x={v[3][0]:.9g}->P={v[3][1]:.9g}", file=sys.stderr)

    # ---- Strata ----
    def stratum(a, b, x, P, Q):
        mn = min(a, b)
        small = min(P, Q)
        import math
        if (not (0.0 < x < 1.0) or not math.isfinite(a)
                or not math.isfinite(b) or P != P or Q != Q):
            return "boundary-special"
        if (P, Q) in ((1.0, 0.0), (0.0, 1.0)):
            return "saturated"
        if identity_value(a, b, x) is not None:
            return "identity"
        if max(a, b) >= b_gl and mn < b_gl:
            return "gammalim"
        if a == b and abs(x - 0.5) < 1e-9:
            return "near-diagonal"
        if small != 0 and small < 2.3e-308:
            return "subnormal"
        if mn <= 2.0 ** -4:
            if small != 0 and small < mn * 1e-6:
                return "smalltau-deep"
            return "smalltau-typ"
        if max(a, b) > 1e4:
            return "large-param"
        return "normal"

    strata = {}
    for r in rows:
        strata.setdefault(stratum(*r), []).append(r)
    rng = random.Random(SEED)
    print("strata sizes: " +
          ", ".join(f"{k}={len(v)}" for k, v in sorted(strata.items())),
          file=sys.stderr)

    fails = []
    declined = []     # near-tie rows: rounding undecidable, not certified
    certified = []    # (stored, indep, tag) for the post-run 1-ULP control
    checked = {}

    def check_value(r, indep, tag):
        a, b, x, P, Q = r
        stored = min(P, Q)
        if indep == 0 and stored == 0:
            return True
        verdict = rounding_verdict(stored, indep)
        if verdict == "pass":
            certified.append((stored, indep, tag))
            return True
        if verdict == "tie":
            declined.append((tag, r, stored, mp.nstr(indep, 25)))
            return False
        fails.append((tag, r, stored, mp.nstr(indep, 17)))
        return False

    # identity stratum: exact, check every row (cheap); dps 80 keeps the
    # closed forms far inside INDEP_REL even where t = b*log1p(-x) is
    # large.
    n_id = 0
    for r in strata.get("identity", []):
        a, b, x, P, Q = r
        with mpmath.workdps(80):
            iv = identity_value(a, b, x)
        Pi, Qi = iv
        indep = Pi if P <= Q else Qi
        check_value(r, indep, "identity")
        n_id += 1
    checked["identity"] = n_id

    # gammalim: sample PER (s-bucket, t/s side) cell, capped per cell, so a
    # populated sub-regime is never starved by a large one; name the
    # cells the stored rows leave empty (#54 sat in such a gap).
    gl_cells = {}
    for r in strata.get("gammalim", []):
        gl_cells.setdefault(gl_cell(r[0], r[1], r[2], b_gl), []).append(r)
    all_cells = [(sb, side) for sb in GL_S_BUCKETS for side in ("lo", "hi")]
    per_cell = max(4, N_PER_STRATUM // max(1, len(gl_cells)))
    gl_rng = random.Random(SEED + 1)   # own stream: other strata unchanged
    gl_sample = []
    for c in all_cells:
        pool = gl_cells.get(c, [])
        gl_sample += pool if len(pool) <= per_cell else gl_rng.sample(pool, per_cell)
    print(f"  gammalim (max param >= kBetaGammaLim = 2^"
          f"{round(math.log2(b_gl))}, min below): stored/sampled per cell: " +
          ", ".join(f"{c[0]}/{c[1]}={len(gl_cells.get(c, []))}/"
                    f"{min(len(gl_cells.get(c, [])), per_cell)}"
                    for c in all_cells), file=sys.stderr)
    gl_empty = [c for c in all_cells if c not in gl_cells]
    if gl_empty:
        print("  gammalim cells with NO stored rows (fresh sampling only): " +
              ", ".join(f"{c[0]}/{c[1]}" for c in gl_empty), file=sys.stderr)

    # series/half-split-verifiable strata
    for name in ("normal", "smalltau-typ", "smalltau-deep", "large-param",
                 "gammalim", "subnormal"):
        pool = strata.get(name, [])
        if name == "gammalim":
            sample = gl_sample
        else:
            sample = pool if len(pool) <= N_PER_STRATUM else rng.sample(pool, N_PER_STRATUM)
        n_done = 0
        self_pairs = 0
        for r in sample:
            a, b, x, P, Q = r
            ap, bp, xp, dp, stored = small_frame(a, b, x, mp.mpf(P), mp.mpf(Q))
            try:
                iv = indep_small(ap, bp, xp, dp)
            except Exception as e:
                fails.append((f"{name}:indep-error", r, float(stored), repr(e)))
                continue
            if name == "normal" and betainc300_ok(a, b, x) and self_pairs < 20:
                with mpmath.workdps(300):
                    bv = mp.betainc(ap, bp, 0, xp, regularized=True)
                if bv != 0 and abs((iv - bv) / bv) > SELF_TOL:
                    fails.append(("SELF-CERT", r, mp.nstr(iv, 17), mp.nstr(bv, 17)))
                self_pairs += 1
            check_value(r, iv, name)
            n_done += 1
        checked[name] = n_done

    # fresh gamma-limit points: certify the evaluator across the slice
    # domain regardless of where stored rows happen to sit.
    n_fresh, fresh_cells = gammalim_fresh(random.Random(SEED + 2), b_gl, fails)
    checked["gammalim-fresh"] = n_fresh
    print(f"  gammalim-fresh: {n_fresh} fresh points (indep vs betainc300 "
          f"<= {mp.nstr(SELF_TOL, 2)}) over cells " +
          ", ".join(f"{c[0]}/{c[1]}" for c in sorted(fresh_cells)),
          file=sys.stderr)

    # 1-ULP negative control on REAL stored rows: every certified row,
    # stored value nudged one ulp either way, must be rejected.
    ulp_leaks = []
    for stored, iv, tag in certified:
        for nudged in (math.nextafter(stored, math.inf),
                       math.nextafter(stored, -math.inf)):
            if rounding_verdict(nudged, iv) != "fail":
                ulp_leaks.append((tag, stored, nudged))
    print(f"  1-ULP control: {2 * len(certified)} nudges of "
          f"{len(certified)} certified stored rows, {len(ulp_leaks)} "
          f"not rejected", file=sys.stderr)

    # saturated stratum: bound certification
    pool = strata.get("saturated", [])
    sample = pool if len(pool) <= N_PER_STRATUM else rng.sample(pool, N_PER_STRATUM)
    n_cert = n_uncert = 0
    for r in sample:
        a, b, x, P, Q = r
        ap, bp, xp, dp, stored = small_frame(a, b, x, mp.mpf(P), mp.mpf(Q))
        lb = log_bound_small(ap, bp, xp, dp)
        if lb < mp.log(MIN_SUBN_HALF):
            n_cert += 1
        else:
            # bound alone can't certify; evaluate directly
            try:
                iv = indep_small(ap, bp, xp, dp)
                if iv < MIN_SUBN_HALF:
                    n_cert += 1
                else:
                    fails.append(("saturated:NOT-JUSTIFIED", r, 0.0, mp.nstr(iv, 17)))
            except Exception as e:
                n_uncert += 1
    checked["saturated"] = n_cert
    if n_uncert:
        print(f"  saturated: {n_uncert} sample rows UNCOVERED (bound "
              f"inconclusive, evaluator failed)", file=sys.stderr)

    # near-diagonal: NAMED as covered by its own quad validation at
    # creation; not re-derivable by this harness's quad.
    nd = len(strata.get("near-diagonal", []))
    print(f"  near-diagonal: {nd} rows NOT covered here (spike width "
          f"~1/sqrt(a) below quad node resolution; certified separately "
          f"by the round-3 bracket derivation)", file=sys.stderr)

    # boundary-special (#13 N14.2): three row shapes share this label.
    # (1) DEFINITIONAL rows -- finite a,b > 0 at x = 0 or x = 1: I_0 = 0
    #     and I_1 = 1 by the definition of the regularized incomplete
    #     beta, no limit involved. Checked EXACTLY, every row.
    # (2) Infinite/degenerate-parameter rows (a or b in {0, inf}, NaN
    #     params, out-of-domain x): the stored value is corvus's API
    #     contract for a convention-dependent limit (e.g. x = 1 with
    #     a = inf stores NaN while finite a stores 1). No oracle exists
    #     but the documented contract itself; NOT judged here.
    # (3) All rows in the stratum are exact-or-NaN by construction, so a
    #     structural invariant holds regardless of shape: (P, Q) must be
    #     one of (0,1), (1,0), (NaN,NaN). Any other pair is corruption.
    n_def = 0
    n_contract = 0
    for r in strata.get("boundary-special", []):
        a, b, x, P, Q = r
        pq_ok = ((P, Q) in ((0.0, 1.0), (1.0, 0.0))
                 or (P != P and Q != Q))
        if not pq_ok:
            fails.append(("boundary-special:STRUCTURE", r, float("nan"),
                          "expected (0,1), (1,0), or (NaN,NaN)"))
            continue
        a_fin_pos = math.isfinite(a) and a > 0
        b_fin_pos = math.isfinite(b) and b > 0
        if a_fin_pos and b_fin_pos and x in (0.0, 1.0):
            want = (0.0, 1.0) if x == 0.0 else (1.0, 0.0)
            if (P, Q) != want:
                fails.append(("boundary-special:DEFINITIONAL", r,
                              float(P), f"want P={want[0]} Q={want[1]}"))
            n_def += 1
        else:
            n_contract += 1
    checked["boundary-special"] = n_def
    print(f"  boundary-special: {n_def} definitional x=0/x=1 rows checked "
          f"exactly; {n_contract} contract rows structure-checked only "
          f"(convention-dependent limits -- oracle is corvus.h itself)",
          file=sys.stderr)

    print("\ncoverage: " + ", ".join(f"{k}={v}" for k, v in sorted(checked.items())),
          file=sys.stderr)
    print(f"certified correctly rounded: {len(certified)}; "
          f"DECLINED near-tie (not certified, within ulp/2*(1+"
          f"{mp.nstr(INDEP_REL, 1)}), not failures): {len(declined)}",
          file=sys.stderr)
    for tag, r, stored, indep in declined:
        print(f"    DECLINED {tag}: a={r[0].hex()} b={r[1].hex()} "
              f"x={r[2].hex()} stored={stored.hex()} indep={indep}",
              file=sys.stderr)
    print(f"FAILURES: {len(fails)} (+{len(bad_mono)} monotonicity, "
          f"{bad_pq} P+Q)", file=sys.stderr)
    hist = {}
    for tag, *_ in fails:
        hist[tag] = hist.get(tag, 0) + 1
    if hist:
        print("failure histogram: " +
              ", ".join(f"{k}={v}" for k, v in sorted(hist.items())),
              file=sys.stderr)
    dump = os.environ.get("VERIFY_BETA_DUMP")
    if dump and (fails or bad_mono):
        with open(dump, "w") as fh:
            for tag, r, stored, indep in fails:
                fh.write(f"{tag}\t{r[0].hex()}\t{r[1].hex()}\t{r[2].hex()}\t"
                         f"{r[3].hex()}\t{r[4].hex()}\t{stored!r}\t{indep}\n")
            for a, b, p0, p1 in bad_mono:
                fh.write(f"MONO\t{a.hex()}\t{b.hex()}\t{p0[0].hex()}\t"
                         f"{p0[1].hex()}\t{p1[0].hex()}\t{p1[1].hex()}\n")
        print(f"full failure dump -> {dump}", file=sys.stderr)
    for tag, r, stored, indep in fails[:25]:
        a, b, x = r[0], r[1], r[2]
        print(f"    {tag}: ({a:.6g},{b:.6g},{x:.6g}) stored={stored!r} "
              f"indep={indep}", file=sys.stderr)
    if ulp_leaks:
        for tag, stored, nudged in ulp_leaks[:10]:
            print(f"    1-ULP LEAK {tag}: stored={stored.hex()} "
                  f"nudged={nudged.hex()} accepted", file=sys.stderr)
        print("harness 1-ULP control failed -- certification void.",
              file=sys.stderr)
        return 2
    return 1 if (fails or bad_mono or bad_pq) else 0


if __name__ == "__main__":
    sys.exit(main())
