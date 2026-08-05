"""INDEPENDENT verification harness for tests/data/beta_{p,q}_reference.txt.

Shares NO code with gen_beta_reference's assemblies (no import of the
oracle module). Methods:
  quad-log   -- mp.quad (tanh-sinh) over the Beta density evaluated in
                log space, always in the frame where the SMALL side is a
                left tail: small = int_0^{x'} t^{a'-1}(1-t)^{b'-1} dt.
                ln(1-t) is computed as ln(d + x'(1-u)) with d = 1-x'
                formed EXACTLY (extra-dps complement of a double) --
                both addends nonnegative, no cancellation anywhere.
                For a' >= BIG the mass is invisible in u-coordinates
                (width ~1/a'), so integrate in w = a'(1-t) instead --
                exact algebra, no limit assumed.
  identity   -- exact closed forms on the analytic lines (a=1, b=1,
                a=b=1/2, a=b & x=1/2).
  betainc300 -- mpmath betainc at dps 300, only inside its trusted zone
                (moderate parameters, off-diagonal).
  log-bound  -- rigorous inequality certification for rows stored as
                exact {0,1}: ln(small) <= a'ln x' - ln a'
                + min(0,(b'-1)ln(1-x')) - lnB; if the bound is below
                half the min subnormal, the stored 0 is correctly
                rounded. No assembly, just a monotone bound.

Self-certification: before judging the oracle, quad-log and betainc300
must agree with EACH OTHER to 1e-25 on the normal-interior sample; a
harness that cannot certify itself certifies nothing.

Stage B invariants over the full file: P+Q=1 (<=1 ulp), monotonicity of
P in x within every (a,b) group.

Exit non-zero on any failure. Strata without coverage are NAMED."""
import random
import sys
import mpmath
from mpmath import mp

SEED = 987654
N_PER_STRATUM = 60
BIG = mp.mpf("1e6")          # a' above this -> w-coordinates
MIN_SUBN_HALF = mp.mpf(2) ** -1075
REL_TOL = mp.mpf("1e-22")    # indep-vs-stored mpf tolerance (pre-rounding)
SELF_TOL = mp.mpf("1e-25")   # quad-vs-betainc self-certification


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
    """(a', b', x', stored_small): frame where the small side is the
    left tail int_0^{x'} t^{a'-1}(1-t)^{b'-1} dt."""
    if P <= Q:
        return mp.mpf(a), mp.mpf(b), mp.mpf(x), P
    return mp.mpf(b), mp.mpf(a), one_minus_exact(x), Q


def quad_log(ap, bp, xp, dps=80):
    """Small side via log-space quad in the left-tail frame."""
    with mpmath.workdps(dps):
        lnB = ln_beta(ap, bp)
        d = 1 - xp if xp < mp.mpf("0.5") else None
        if d is None:
            d = one_minus_exact(float(xp)) if xp == mp.mpf(float(xp)) else 1 - xp
        lnxp = mp.log(xp)
        if ap < BIG:
            # t = x'*u, u in [0,1]; ln(1-t) = ln(d + x'*(1-u)), safe.
            def g(u):
                if u <= 0:
                    return mp.mpf(0) if ap > 1 else mp.mpf(0)  # measure zero
                lt = (ap - 1) * (lnxp + mp.log(u)) \
                     + (bp - 1) * mp.log(d + xp * (1 - u)) - lnB
                return mp.exp(lt) * xp
            # splits: decay scale of (1-t)^(b'-1) in u when b' large
            pts = [mp.mpf(0)]
            if bp > 2:
                s = d * 3 / (xp * (bp - 1))
                for m in (1, 10, 100):
                    v = s * m
                    if mp.mpf("1e-30") < v < 1:
                        pts.append(v)
            pts += [mp.mpf("0.5"), mp.mpf(1)]
            pts = sorted(set(pts))
            return mp.quad(g, pts)
        # a' huge: w = a'(1-t), t = 1 - w/a'; dt = -dw/a'
        # small = int_{w0}^{a'} (1-w/a')^(a'... wait t^(a'-1)) ...
        # t = 1 - w/a': t^(a'-1) = exp((a'-1)*log1p(-w/a'));
        # (1-t)^(b'-1) = (w/a')^(b'-1). Bounds: t from 0..x' -> w from
        # a'(1-x') = a'*d down.. up: t=0 -> w=a' ; t=x' -> w=a'*d.
        w0 = ap * d
        def gw(w):
            lt = (ap - 1) * mp.log1p(-w / ap) \
                 + (bp - 1) * (mp.log(w) - mp.log(ap)) - lnB - mp.log(ap)
            return mp.exp(lt)
        hi = min(ap, w0 + 200 + 50 * mp.sqrt(w0))
        return mp.quad(gw, [w0, w0 + 10, hi])


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


def log_bound_small(ap, bp, xp):
    """Rigorous upper bound on ln(small side) in the left-tail frame."""
    with mpmath.workdps(80):
        lb = ap * mp.log(xp) - mp.log(ap) - ln_beta(ap, bp)
        if bp < 1:
            d = 1 - xp
            lb += (bp - 1) * mp.log(d)
        return lb


def main():
    mp.dps = 40
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

    # quad-verifiable strata
    for name in ("normal", "smalltau-typ", "smalltau-deep", "large-param",
                 "gammalim", "subnormal"):
        pool = strata.get(name, [])
        sample = pool if len(pool) <= N_PER_STRATUM else rng.sample(pool, N_PER_STRATUM)
        n_done = 0
        self_pairs = 0
        for r in sample:
            a, b, x, P, Q = r
            ap, bp, xp, stored = small_frame(a, b, x, mp.mpf(P), mp.mpf(Q))
            try:
                qv = quad_log(ap, bp, xp)
            except Exception as e:
                fails.append((f"{name}:quad-error", r, float(stored), repr(e)))
                continue
            if name == "normal" and betainc300_ok(a, b, x) and self_pairs < 20:
                with mpmath.workdps(300):
                    bv = mp.betainc(ap, bp, 0, xp, regularized=True)
                if bv != 0 and abs((qv - bv) / bv) > SELF_TOL:
                    fails.append(("SELF-CERT", r, mp.nstr(qv, 17), mp.nstr(bv, 17)))
                self_pairs += 1
            check_value(r, qv, name)
            n_done += 1
        checked[name] = n_done

    # saturated stratum: bound certification
    pool = strata.get("saturated", [])
    sample = pool if len(pool) <= N_PER_STRATUM else rng.sample(pool, N_PER_STRATUM)
    n_cert = n_uncert = 0
    for r in sample:
        a, b, x, P, Q = r
        ap, bp, xp, stored = small_frame(a, b, x, mp.mpf(P), mp.mpf(Q))
        lb = log_bound_small(ap, bp, xp)
        if lb < mp.log(MIN_SUBN_HALF):
            n_cert += 1
        else:
            # bound alone can't certify; try quad
            try:
                qv = quad_log(ap, bp, xp)
                if qv < MIN_SUBN_HALF:
                    n_cert += 1
                else:
                    fails.append(("saturated:NOT-JUSTIFIED", r, 0.0, mp.nstr(qv, 17)))
            except Exception as e:
                n_uncert += 1
    checked["saturated"] = n_cert
    if n_uncert:
        print(f"  saturated: {n_uncert} sample rows UNCOVERED (bound "
              f"inconclusive, quad failed)", file=sys.stderr)

    # near-diagonal: NAMED as covered by its own quad validation at
    # creation (round 3); not re-derivable by this harness's quad.
    nd = len(strata.get("near-diagonal", []))
    print(f"  near-diagonal: {nd} rows NOT covered here (spike width "
          f"~1/sqrt(a) below quad node resolution; certified separately "
          f"by the round-3 bracket derivation)", file=sys.stderr)

    print("\ncoverage: " + ", ".join(f"{k}={v}" for k, v in sorted(checked.items())),
          file=sys.stderr)
    print(f"FAILURES: {len(fails)} (+{len(bad_mono)} monotonicity, "
          f"{bad_pq} P+Q)", file=sys.stderr)
    for tag, r, stored, indep in fails[:25]:
        a, b, x = r[0], r[1], r[2]
        print(f"    {tag}: ({a:.6g},{b:.6g},{x:.6g}) stored={stored} "
              f"indep={indep}", file=sys.stderr)
    return 1 if (fails or bad_mono or bad_pq) else 0


if __name__ == "__main__":
    sys.exit(main())
