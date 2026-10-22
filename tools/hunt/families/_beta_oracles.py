"""Bounded oracles for the beta family hunt strata (#82).

Never calls mpmath.betainc (it hangs on magnitude-mismatched and
near-diagonal shapes and returns garbage in the gamma-limit band). Forms:

  series   DLMF 8.17.8, all terms positive, in the frame whose argument is
           <= 0.75:  I = x^a (1-x)^b / (a B) * sum_k (a+b)_k/(a+1)_k x^k.
           Capped; raises on the cap.
  split    x' > 0.75: series at 1/2 plus a log-coordinate quadrature of
           the density from 1/2 to x' (v = ln(1-t)); both pieces positive.
  quadC    endpoint-relative tanh-sinh quadrature of the small side below
           the mode (shapes > 1): I = x^(a-1) y^(b-1)/B *
           int_0^x (1-s/x)^(a-1) (1+s/y)^(b-1) ds; cost independent of
           the shapes, so it carries the large-shape (R3, PB, gamma-limit
           ridge) traffic.
  cf       modified-Lentz continued fraction (DLMF 8.17.22), capped,
           raises on the cap. Used only as the independent second form.
  satur    rigorous upper bound on the small side (series ratio bound);
           below 2^-1100 the side is returned as exact 0 (its complement 1).

The small side is always computed directly; the complement is formed at
working precision only when the computed side is <= ~1/2.
"""
from __future__ import annotations

import math

import mpmath as mp
from mpmath import mpf

SERIES_CAP = 200000
CF_CAP = 40000
SAT_LN = -1100 * math.log(2)


class OracleError(ArithmeticError):
    pass


# ---- ln B -------------------------------------------------------------------
_LNB_CACHE: dict = {}


def ln_beta(a, b, dps):
    """ln B(a, b) to ~dps digits relative to max(|lnB|, 1), cancellation-aware.
    Cached on (a, b, dps): the inverse calls it per Newton step."""
    key = (a, b, dps)
    v = _LNB_CACHE.get(key)
    if v is None:
        if len(_LNB_CACHE) > 4096:
            _LNB_CACHE.clear()
        v = _LNB_CACHE[key] = _ln_beta(a, b, dps)
    return v


def _ln_beta(a, b, dps):
    a, b = _mp(a), _mp(b)
    m, M = (a, b) if a <= b else (b, a)
    if M <= 1:
        # lnB = ln((m+M)/(mM)) + lg1p(m) + lg1p(M) - lg1p(m+M), lg1p = lgamma(1+t)
        with mp.workdps(dps + 20):
            c = mp.fadd(m, M, exact=True)
            return (mp.log(c) - mp.log(m) - mp.log(M) + mp.loggamma(1 + m)
                    + mp.loggamma(1 + M) - mp.loggamma(1 + c))
    with mp.workdps(dps + 10):
        if m < mpf(10) ** (-(dps // 2 + 8)):
            # lgamma(M+m) - lgamma(M) = m psi(M) + m^2/2 psi1(M) + O(m^3), M > 1
            lgm = mp.loggamma(m)
            diff = m * mp.digamma(M) + m * m / 2 * mp.psi(1, M)
            return lgm - diff
    # absolute accuracy 10^-dps: cover the a+b cancellation (M/m) AND the
    # magnitude of lgamma(M) ~ M ln M (huge near-diagonal pairs)
    extra = max(0, int(mp.log10(M / m))) + 15
    if M > 10:
        extra += int(mp.log10(M * mp.log(M))) + 1
    with mp.workdps(dps + extra):
        return mp.loggamma(m) + mp.loggamma(M) - mp.loggamma(mp.fadd(M, m, exact=True))


# ---- forms ------------------------------------------------------------------
def _series(al, be, xi, lnB, dps):
    with mp.workdps(dps + 10):
        lnpre = al * mp.log(xi) + be * mp.log1p(-xi) - mp.log(al) - lnB
        s = mpf(1)
        term = mpf(1)
        eps = mpf(10) ** (-dps - 8)
        c = al + be
        k = 0
        while True:
            k += 1
            term *= (c + k - 1) / (al + k) * xi
            s += term
            if term < s * eps and (c + k) * xi < (al + k):
                break
            if k > SERIES_CAP:
                raise OracleError(f"series cap (x'={mp.nstr(xi, 6)})")
        return mp.exp(lnpre + mp.log(s))


def _split_pts(lo, hi, rate):
    """Split points from lo upward: lo + h*2^j with h = 1/(8 rate)."""
    pts = [lo]
    h = 1 / (8 * max(abs(rate), mpf(1) / 64))
    p = lo + h
    while p < hi:
        pts.append(p)
        h *= 2
        p = lo + h
    pts.append(hi)
    return pts


def _tail_logquad(ap, bp, d, lnB, dps):
    """int_{1/2}^{1-d} t^(a'-1)(1-t)^(b'-1) dt / B via v = ln(1-t)."""
    with mp.workdps(dps + 10):
        lnd = mp.log(d)
        vhi = -mp.log(2)

        def lgv(v):
            return bp * v + (ap - 1) * mp.log1p(-mp.exp(v))

        ev = mp.exp(lnd)
        rate = bp - (ap - 1) * ev / (1 - ev)
        # absolute-error quad: integrate in w = (v - lnd)/h0, h0 the decay
        # scale at the lower end, so the integral is O(1).
        h0 = min(1 / max(abs(rate), mpf(1) / 64), vhi - lnd)
        wpts = [(p - lnd) / h0 for p in _split_pts(lnd, vhi, rate)]
        if ap > 2:
            lna = mp.log(ap)
            for off in (-3, -1, 0, 1, 3, 6):
                c = -lna + off
                if lnd < c < vhi:
                    wpts.append((c - lnd) / h0)
        wpts = sorted(set(wpts))
        lmax = max(lgv(lnd + h0 * w) for w in wpts)
        val = h0 * mp.quad(lambda w: mp.exp(lgv(lnd + h0 * w) - lmax), wpts)
        return val * mp.exp(lmax - lnB)


def _at(dps, fn):
    with mp.workdps(dps):
        return fn()


def _quadC(al, be, xi, y, lnB, dps):
    """Small side below the mode (al, be > 1, k > 0)."""
    # The integrand's two terms are each ~(al-1) s/xi and cancel to k s:
    # carry log10((al-1)/(xi k)) extra digits through it (large only near
    # the mean of a huge-nu pair). k itself needs log10(al+be) digits.
    exc = int(mp.log10(al + be)) + 5
    with mp.workdps(dps + 10 + exc):
        a1, b1 = al - 1, be - 1
        k = a1 / xi - b1 / y
        if not k > 0:
            raise OracleError("quadC above the mode")
        ex = max(0, int(mp.log10(a1 / xi / k))) + 5
    with mp.workdps(dps + 10 + ex):
        def lg(s):
            r = s / xi
            if r >= 1:
                return mpf("-inf")
            return a1 * mp.log1p(-r) + b1 * mp.log1p(s / y)

        cut = (dps + 10) * 2.31 + 30
        # The integrand decays like exp(-k s - (c/(2 x y)) s^2). mp.quad's
        # error test is ABSOLUTE, so integrate in w = s/h0 (h0 the decay
        # scale): the integral is then O(1) whatever the shapes.
        sig = mp.sqrt(xi * y / (al + be))
        h0 = min(1 / k, sig)
        wmax = xi / h0
        h = mpf(1) / 8
        pts = [mpf(0)]
        for _ in range(3000):
            w = pts[-1] + h
            if w >= wmax:
                pts.append(wmax)
                break
            pts.append(w)
            if lg(h0 * w) < -cut:
                break
            h *= 2
        else:
            raise OracleError("quadC split")
        with mp.workdps(dps + 10):
            # quadrature at the base precision; lg raises its own
            integral = h0 * mp.quad(lambda w: _at(dps + 10 + ex, lambda: mp.exp(lg(h0 * w))), pts)
    # prefactor carries the huge-shape cancellation: elevated precision
    with mp.workdps(dps + 10 + exc):
        pre = a1 * mp.log(xi) + b1 * mp.log1p(-xi) - lnB
        return mp.exp(pre) * integral


def cf_value(al, be, xi, y, lnB, dps, cap=CF_CAP):
    with mp.workdps(dps + 10):
        eps = mpf(2) ** (-mp.mp.prec + 4)
        tiny = mpf(10) ** (-(dps + 300))
        c = mpf(1)
        d = 1 - (al + be) * xi / (al + 1)
        d = tiny if d == 0 else d
        d = 1 / d
        h = d
        for m in range(1, cap):
            aa = m * (be - m) * xi / ((al + 2 * m - 1) * (al + 2 * m))
            d = 1 + aa * d
            d = tiny if d == 0 else d
            c = 1 + aa / c
            c = tiny if c == 0 else c
            d = 1 / d
            h *= d * c
            aa = -(al + m) * (al + be + m) * xi / ((al + 2 * m) * (al + 2 * m + 1))
            d = 1 + aa * d
            d = tiny if d == 0 else d
            c = 1 + aa / c
            c = tiny if c == 0 else c
            d = 1 / d
            de = d * c
            h *= de
            if abs(de - 1) < eps:
                break
        else:
            raise OracleError("cf cap")
    ex = int(mp.log10(al + be)) + 5
    with mp.workdps(dps + 10 + ex):
        return mp.exp(al * mp.log(xi) + be * mp.log1p(-xi) - lnB) / al * h


def _side(al, be, xi, y, lnB, dps):
    """I_xi(al, be) by series or split (small/moderate shapes)."""
    if xi <= mpf(0.75):
        return _series(al, be, xi, lnB, dps)
    with mp.workdps(dps + 10):
        half = _series(al, be, mpf(0.5), lnB, dps)
        return half + _tail_logquad(al, be, y, lnB, dps)


def _sat_bound(al, be, xi, lnB, dps):
    """ln of a rigorous upper bound on I_xi(al, be) (frame below the mean),
    or None when the series ratio bound is not < 1."""
    with mp.workdps(dps + 10):
        r = max((al + be) * xi / (al + 1), xi)
        if not r < 1 - mpf(10) ** -6:
            return None
        return (al * mp.log(xi) + be * mp.log1p(-xi) - mp.log(al) - lnB
                - mp.log1p(-r))


def _one_minus(x, dps):
    """1 - x, exact for any double-derived x; an astronomically small mpf x
    (an inverse iterate) gets 1 - x at dps + 60 instead."""
    if x < mpf(2) ** -1200:
        with mp.workdps(dps + 60):
            return 1 - x
    return mp.fsub(1, x, exact=True)


def _mp(v):
    """mpf without rounding an mpf argument to the ambient precision."""
    return v if isinstance(v, mp.mpf) else mpf(v)


def betainc_pq(a, b, x, dps):
    """(P, Q) = (I_x(a,b), 1 - I_x(a,b)), each to ~dps digits relative."""
    a, b, x = _mp(a), _mp(b), _mp(x)
    if x == 0:
        return mpf(0), mpf(1)
    if x == 1:
        return mpf(1), mpf(0)
    y = _one_minus(x, dps)
    lnB = ln_beta(a, b, dps + 10)
    with mp.workdps(dps + 10):
        lam = a * y - b * x
    nat = lam >= 0   # x at or below the mean: P is (nearly always) small
    # saturation shortcut on the mean-frame small side
    lb = _sat_bound(a, b, x, lnB, dps) if nat else _sat_bound(b, a, y, lnB, dps)
    if lb is not None and lb < SAT_LN:
        return (mpf(0), mpf(1)) if nat else (mpf(1), mpf(0))
    if min(a, b) >= 30:
        with mp.workdps(dps + 10):
            kn = (a - 1) / x - (b - 1) / y
        if kn == 0 and a == b:
            return mpf(0.5), mpf(0.5)
        if kn > 0:
            v = _quadC(a, b, x, y, lnB, dps)
            with mp.workdps(dps + 20):
                return +v, 1 - v
        v = _quadC(b, a, y, x, lnB, dps)
        with mp.workdps(dps + 20):
            return 1 - v, +v
    if nat:
        v = _side(a, b, x, y, lnB, dps)
        if v <= 0.5:
            with mp.workdps(dps + 20):
                return +v, 1 - v
        w = _side(b, a, y, x, lnB, dps)
        _agree(v, w, dps)
        return v, w
    w = _side(b, a, y, x, lnB, dps)
    if w <= 0.5:
        with mp.workdps(dps + 20):
            return 1 - w, +w
    v = _side(a, b, x, y, lnB, dps)
    _agree(v, w, dps)
    return v, w


def _agree(v, w, dps):
    with mp.workdps(dps + 10):
        if abs(v + w - 1) > mpf(10) ** (-(dps - 8)):
            raise OracleError(f"P+Q-1 = {mp.nstr(v + w - 1, 5)}")


def betainc_pq_cf(a, b, x, dps):
    """Second form: Lentz CF on the small side (frame picked from a cheap
    primary evaluation), working precision raised by log10(a+b) for the
    CF's conditioning. Where the CF cannot converge in its cap (tiny-shape
    frames above the CF's convergence edge, huge near-diagonal shapes) it
    falls back to the primary forms at dps + 30 -- NOT independent there."""
    a, b, x = _mp(a), _mp(b), _mp(x)
    if x == 0:
        return mpf(0), mpf(1)
    if x == 1:
        return mpf(1), mpf(0)
    p0, q0 = betainc_pq(a, b, x, 25)
    if min(p0, q0) == 0:
        return p0, q0
    y = _one_minus(x, dps)
    nat = p0 <= q0
    al, be, xi, yy = (a, b, x, y) if nat else (b, a, y, x)
    ex = int(mp.log10(a + b + 1)) + 10
    try:
        with mp.workdps(dps + ex):
            if not xi < (al + 1) / (al + be + 2) or min(al, be) < mpf(10) ** -6:
                raise OracleError("cf edge")
        lnB = ln_beta(a, b, dps + ex)
        v = cf_value(al, be, xi, yy, lnB, dps + ex, cap=CF_CAP)
    except OracleError:
        return betainc_pq(a, b, x, dps + 30)
    with mp.workdps(dps + 20):
        w = 1 - v
    return (v, w) if nat else (w, v)


# ---- float approximation (generators only; never an oracle) -----------------
def approx_I(a, b, x):
    """Double-precision I_x(a,b) (NR betai), ~1e-10 for moderate shapes;
    normal approximation for huge shapes. For placing points only."""
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    try:
        if a + b > 1e7:
            m = a / (a + b)
            s = math.sqrt(a * b / ((a + b) ** 2 * (a + b + 1)))
            return 0.5 * math.erfc(-(x - m) / s / math.sqrt(2))
        lbt = (math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
               + a * math.log(x) + b * math.log1p(-x))
        if x < (a + 1) / (a + b + 2):
            return math.exp(lbt) * _cf_float(a, b, x) / a
        return 1.0 - math.exp(lbt) * _cf_float(b, a, 1 - x) / b
    except (OverflowError, ValueError):
        return 0.5


def _cf_float(a, b, x):
    tiny = 1e-300
    qab, qap, qam = a + b, a + 1, a - 1
    c, d = 1.0, 1 - qab * x / qap
    d = tiny if abs(d) < tiny else d
    d = 1 / d
    h = d
    for m in range(1, 3000):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1 + aa * d
        d = tiny if abs(d) < tiny else d
        c = 1 + aa / c
        c = tiny if abs(c) < tiny else c
        d = 1 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1 + aa * d
        d = tiny if abs(d) < tiny else d
        c = 1 + aa / c
        c = tiny if abs(c) < tiny else c
        d = 1 / d
        de = d * c
        h *= de
        if abs(de - 1) < 1e-15:
            break
    return h


def approx_root(a, b, target):
    """x with approx_I(a,b,x) ~ target (bisection in doubles)."""
    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mid <= lo or mid >= hi:
            break
        if approx_I(a, b, mid) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ---- inverse ----------------------------------------------------------------
def _norm_quantile(t):
    """z >= 0 with Phi(-z) = t for 0 < t <= 1/2 (float; seeds only)."""
    if t >= 0.5:
        return 0.0
    if t < 1e-300:
        L = -2 * math.log(t)
        return math.sqrt(L - math.log(L) - math.log(2 * math.pi))
    lo, hi = 0.0, 40.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if 0.5 * math.erfc(mid / math.sqrt(2)) > t:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def beta_inv_root(a, b, p, side, dps, trace=False):
    """Exact x with P(x) = p (side 'p') or Q(x) = p (side 'q'), as mpf.

    Solves in the frame whose root is <= ~1/2 (I_x(a,b) = 1 - I_{1-x}(b,a))
    and on the side whose target is <= 1/2, Newton in u = ln(xi) on
    g = ln S(xi) - ln T with a maintained bracket (bisection fallback).
    Cheap low-precision iterations first, then full precision."""
    a, b, p = mpf(a), mpf(b), mpf(p)
    if p == 0:
        return mpf(0) if side == "p" else mpf(1)
    if p == 1:
        return mpf(1) if side == "p" else mpf(0)
    # frame: the root is <= 1/2 iff P(1/2) >= p (side p) / Q(1/2) <= q (side q)
    # (float approximation when it is unambiguous; exact otherwise)
    low = None
    if a + b <= 1e7:
        ph = approx_I(float(a), float(b), 0.5)
        pf = float(p) if side == "p" else 1.0 - float(p)
        if 0 < ph < 1 and abs(ph - pf) > 1e-6 * min(pf, 1 - pf) + 1e-300:
            low = ph >= pf
    if low is None:
        P_half, Q_half = betainc_pq(a, b, mpf(0.5), 30)
        low = (P_half >= p) if side == "p" else (Q_half <= p)
    if low:
        al, be, flip = a, b, False          # xi = x:   P'(xi) = P(x)
        tgt_on_P = side == "p"
    else:
        al, be, flip = b, a, True           # xi = 1-x: P'(xi) = Q(x)
        tgt_on_P = side == "q"
    # solve on the side whose target is <= 1/2 (1 - p is exact for p >= 1/2)
    if p <= 0.5:
        use_small, Tt = tgt_on_P, p         # use_small: the P' side
    else:
        use_small, Tt = not tgt_on_P, mpf(1) - p
    with mp.workdps(dps + 20):
        lnT = mp.log(Tt)
    lnB = ln_beta(al, be, dps + 10)

    def g(u, wd):
        with mp.workdps(wd + 20):
            xi = mp.exp(u)
            if xi >= 1:
                return mpf(1), mpf(0)
        P_, Q_ = betainc_pq(al, be, xi, wd)
        S = P_ if use_small else Q_
        with mp.workdps(wd + 10):
            if S == 0:
                return mpf(-10) ** 9, None
            gv = mp.log(S) - lnT
            lnf = ((al - 1) * u + (be - 1) * mp.log1p(-xi) - lnB)
            dS = mp.exp(lnf + u) / S          # d ln S / du magnitude
            return gv, (dS if use_small else -dS)

    # bracket in u: g increasing in u if use_small, else decreasing
    sgn = 1 if use_small else -1
    lo, hi = None, mp.log(mpf(0.5)) + mpf(1e-9)
    # initial guess: leading power term for small targets at moderate shapes,
    # a normal approximation (in the frame) for large shapes
    nu = float(al * be / (al + be))
    wd1 = 28 + max(0, int(math.log10(nu) / 2)) if nu > 1 else 28
    with mp.workdps(30):
        umin = -1100 * mp.log(2)
        if min(al, be) > 1e3:
            c = al + be
            mean = al / c
            sd = mp.sqrt(al * be / (c * c * (c + 1)))
            z = _norm_quantile(float(Tt))
            z = -z if use_small else z       # P' = Tt: below the mean
            xi0 = mean + z * sd
            if not 0 < xi0 <= 0.5:
                xi0 = min(mpf(0.5), mean) if xi0 > 0 else mean / 4
            u = mp.log(xi0)
        elif use_small and Tt < mpf(0.1):
            u = (lnT + mp.log(al) + lnB) / al
        else:
            u = mp.log(min(al / (al + be), mpf(0.25)))
        u = min(u, hi)
        u = max(u, umin)
    def underflow():
        """True when the root lies below 2^-1100 in this frame (answer 0)."""
        gv0, _ = g(umin, wd1)
        return sgn * gv0 >= 0

    def underflow_result():
        xi = mp.exp(umin - 10)
        if not flip:
            return xi
        with mp.workdps(dps + 400):
            return 1 - xi

    checked = False
    if u < umin + 40:
        checked = True
        if underflow():
            return underflow_result()
    lo = umin
    uhi = mp.log(mpf(0.5)) + mpf(1e-9)
    # stage 2 accepts a Newton step <= 10^-(dps/2 + 2): convergence is
    # quadratic there, so the remaining error is ~step^2 << 10^-dps.
    for wd, tol, nit in ((wd1, mpf(10) ** -(wd1 - 7), 300),
                         (dps + 5, mpf(10) ** -(dps // 2 + 2), 200)):
        verified = wd == wd1
        if wd > wd1:
            # A local bracket around the stage-1 root, verified on the first
            # bisection fallback (the low-precision bracket may be off).
            delta = mpf(10) ** -(wd1 - 12)
            lo, hi = u - delta, min(u + delta, uhi)
        for _ in range(nit):
            gv, dg = g(u, wd)
            s = sgn * gv
            if s < 0:
                lo = u if lo is None or u > lo else lo
            else:
                hi = u if u < hi else hi
            if dg is None or dg == 0:
                un = None
            else:
                with mp.workdps(wd + 10):
                    un = u - gv / dg
            if un is None or not (lo <= un <= hi):
                if not verified:
                    verified = True
                    glo, _ = g(lo, wd)
                    ghi, _ = g(hi, wd)
                    if not (sgn * glo < 0 <= sgn * ghi):
                        lo, hi = umin, uhi
                if not checked and lo == umin and u - umin < 1:
                    checked = True
                    if underflow():
                        return underflow_result()
                un = (lo + hi) / 2
            step = abs(un - u)
            if trace:
                print(wd, mp.nstr(u, 20), mp.nstr(gv, 5), mp.nstr(dg, 5) if dg else dg, mp.nstr(step, 5))
            u = un
            if step <= tol * max(1, abs(u) * mpf(10) ** -6):
                break
            if hi - lo <= mpf(10) ** -(wd - 7):
                break
        else:
            raise OracleError(f"inverse: no convergence at dps {wd}")
    with mp.workdps(dps + 30):
        xi = mp.exp(u)
        if not flip:
            return xi
        return 1 - xi


# ---- harness hooks -----------------------------------------------------------
# hunt.py imports every module in families/ when --families is empty; this
# helper module contributes no strata.
SEAMS: list = []


def strata():
    return []
