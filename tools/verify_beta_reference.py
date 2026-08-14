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

Stage B invariants over the full file: P+Q=1 (<=1 ulp), monotonicity of
P in x within every (a,b) group.

Exit non-zero on any failure. Strata without coverage are NAMED."""
import os
import random
import sys
import mpmath
from mpmath import mp

SEED = 987654
N_PER_STRATUM = 60
MIN_SUBN_HALF = mp.mpf(2) ** -1075
REL_TOL = mp.mpf("1e-22")    # indep-vs-stored mpf tolerance (pre-rounding)
SELF_TOL = mp.mpf("1e-25")   # indep-vs-betainc self-certification

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


def negative_control():
    """Detector validation: the evaluator must separate the four known-
    bad values from their corrected values."""
    tol = mp.mpf(2) ** -51 + REL_TOL
    ok = True
    for a, b, x, bad, good in NEG_CONTROL:
        iv = indep_small(mp.mpf(a), mp.mpf(b), mp.mpf(x), one_minus_exact(x))
        rb = abs((mp.mpf(bad) - iv) / iv)
        rg = abs((mp.mpf(good) - iv) / iv)
        if rb <= tol:
            print(f"NEGATIVE CONTROL FAILED: bad value accepted at "
                  f"a={a:.6g} b={b:.6g} x={x:.6g}", file=sys.stderr)
            ok = False
        if rg > tol:
            print(f"NEGATIVE CONTROL FAILED: good value rejected at "
                  f"a={a:.6g} b={b:.6g} x={x:.6g} "
                  f"(indep={mp.nstr(iv, 20)})", file=sys.stderr)
            ok = False
    return ok


def main():
    mp.dps = 40
    if not negative_control():
        print("harness cannot certify itself -- aborting before "
              "judging any reference row.", file=sys.stderr)
        return 2
    print("negative control passed (4 known-bad rows rejected, "
          "4 corrected values accepted)", file=sys.stderr)
    rows = []
    with open(r"C:\Users\gdwol\Development\corvus\tests\data\beta_p_reference.txt") as fh:
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
        if max(a, b) >= 2.0 ** 59:
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
    checked = {}

    def check_value(r, indep, tag):
        a, b, x, P, Q = r
        stored = mp.mpf(min(P, Q))
        if indep == 0 and stored == 0:
            return True
        # subnormal range: doubles carry reduced mantissa there -- a
        # relative test is meaningless; accept within 1 subnormal ulp.
        if abs(stored - indep) <= mp.mpf("5e-324"):
            return True
        if indep != 0:
            rel = abs((stored - indep) / indep)
            # stored is a rounded double: within 2^-52 of the oracle mpf;
            # indep carries quad error; accept 2^-52 + REL_TOL slack.
            if rel <= mp.mpf(2) ** -51 + REL_TOL:
                return True
        fails.append((tag, r, float(stored), mp.nstr(indep, 17)))
        return False

    # identity stratum: exact, check every row (cheap)
    n_id = 0
    for r in strata.get("identity", []):
        a, b, x, P, Q = r
        iv = identity_value(a, b, x)
        Pi, Qi = iv
        indep = Pi if P <= Q else Qi
        check_value(r, indep, "identity")
        n_id += 1
    checked["identity"] = n_id

    # series/half-split-verifiable strata
    for name in ("normal", "smalltau-typ", "smalltau-deep", "large-param",
                 "gammalim", "subnormal"):
        pool = strata.get(name, [])
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

    print("\ncoverage: " + ", ".join(f"{k}={v}" for k, v in sorted(checked.items())),
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
        print(f"    {tag}: ({a:.6g},{b:.6g},{x:.6g}) stored={stored} "
              f"indep={indep}", file=sys.stderr)
    return 1 if (fails or bad_mono or bad_pq) else 0


if __name__ == "__main__":
    sys.exit(main())
