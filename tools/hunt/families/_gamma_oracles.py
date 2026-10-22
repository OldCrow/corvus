"""Capped exact oracles for the gamma and gammainv hunt families (#82).

Forward: pq(a, x, dps) -> (P, Q), each relatively accurate. The small side is
computed DIRECTLY and the large side as 1 - small at the working precision:
  a <= 3/2, x <= 4   Q by the R4-style small-a identity, P by the series
                     (both direct; P's series has no cancellation there)
  x <  a             P by the power series
  x >= a             Q by Legendre's continued fraction (modified Lentz)
  a*phi > 790        saturated: the small side is below e^-790 (Chernoff
                     bound P or Q <= e^{-a*phi}) and is returned as the
                     leading asymptotic estimate; it rounds to +0 and its
                     real error against any value below 2^-1100 is < 2^-26 ulp.
  a > 1e7            only the saturated band and the exact diagonal x == a
                     for a >= 2^100 (Q = 1/2 - 1/(3 sqrt(2 pi a)) + O(a^-3/2): the median
                     sits below the mean)
                     are covered; anything else raises (row uncertified).
Independent forms (pq2): Gauss's continued fraction for gamma(a,x) (DLMF
8.9.2) for P, a backward fixed-depth Legendre CF with doubling depth for Q,
and mpmath's UPPER gammainc for the small-a box.

Every loop carries an iteration cap and raises OracleError when it is hit.

Inverse: inv_root(kind, a, s, dps, forward) -> exact root as mpf. Solved on
the small side (s > 1/2 is flipped to 1 - s, exact), bracketed in double bit
space from a scipy float guess (guess only; never the answer), bisected to
adjacent doubles, then one secant step between them for the real root.
"""
from __future__ import annotations

import math

import mpmath as mp
from mpmath import mpf

CAP = 1_500_000           # series / CF iteration cap
A_MAX_EXACT = 1e7         # above this only saturated / diagonal rows
SAT_APHI = 790            # Chernoff: small side <= e^-790 < 2^-1139


class OracleError(RuntimeError):
    pass


# ---- helpers ----------------------------------------------------------------
def _guard(a, x) -> int:
    """Guard digits for the prefactor a*log x - x - lgamma(a) cancellation."""
    m = max(1.0, abs(a * math.log(x)) if x > 0 else 1.0, x, a)
    return 12 + int(math.log10(m))


def _aphi(a, x):
    """a*phi(x/a), phi(l) = l - 1 - ln l, at the current precision."""
    a, x = mpf(a), mpf(x)
    u = (x - a) / a
    if abs(u) < 0.5:
        return a * (u - mp.log1p(u))
    return a * (x / a - 1 - mp.log(x / a))


def _lgamma1p(a):
    """lgamma(1+a), relative-accurate for tiny a (the caller adds guard)."""
    a = mpf(a)
    if a > mpf(2) ** -20:
        with mp.extradps(10):
            return mp.loggamma(1 + a)
    # ln Gamma(1+a) = -gamma*a + sum_{k>=2} (-1)^k zeta(k) a^k / k
    s = -mp.euler * a
    term_a = a
    eps = mpf(2) ** (-mp.mp.prec - 8) * abs(s)
    for k in range(2, 200):
        term_a *= a
        t = (-1) ** k * mp.zeta(k) * term_a / k
        s += t
        if abs(t) < eps:
            return s
    raise OracleError("lgamma1p series cap")


def _p_series(a, x):
    a, x = mpf(a), mpf(x)
    eps = mpf(2) ** (-mp.mp.prec - 6)
    t = s = mpf(1)
    n = 0
    while True:
        n += 1
        t *= x / (a + n)
        s += t
        r = a + n + 1
        if x < r and t * x / (r - x) < eps * s:
            break
        if n > CAP:
            raise OracleError(f"P series cap a={a} x={x}")
    return mp.exp(a * mp.log(x) - x - mp.loggamma(a + 1)) * s


def _p_hyp(a, x):
    """The same series via mpmath's fixed-point hypsum (20-30x faster at
    large a): P = x^a e^-x / Gamma(a+1) * 1F1(1; a+1; x). Raises
    NoConvergence at the term cap."""
    a, x = mpf(a), mpf(x)
    s = mp.hyp1f1(1, a + 1, x, maxterms=CAP)
    return mp.exp(a * mp.log(x) - x - mp.loggamma(a + 1)) * s


def _q_cf_lentz(a, x):
    # Only ever called with x >= a: below the ridge forward Lentz
    # false-converges (measured: 1e6 - 2 sqrt(1e6) agrees to 1e-23 only).
    a, x = mpf(a), mpf(x)
    eps = mpf(2) ** (-mp.mp.prec - 4)
    tiny = mpf(2) ** (-4 * mp.mp.prec)
    b = x + 1 - a
    c = 1 / tiny
    d = 1 / b
    h = d
    n = 0
    while True:
        n += 1
        an = -n * (n - a)
        b += 2
        d = an * d + b
        if d == 0:
            d = tiny
        c = b + an / c
        if c == 0:
            c = tiny
        d = 1 / d
        delta = d * c
        h *= delta
        if abs(delta - 1) < eps and n > 2:
            break
        if n > CAP:
            raise OracleError(f"Q CF cap a={a} x={x}")
    return mp.exp(a * mp.log(x) - x - mp.loggamma(a)) * h


def _q_r4(a, x):
    """a <= 3/2, x <= 4:  Gamma(1+a) Q = (Gamma(1+a)-1) - (x^a-1) - a x^a S,
    S = sum_{n>=1} (-x)^n / (n! (a+n)); every bracket relative-accurate."""
    a, x = mpf(a), mpf(x)
    with mp.extradps(12):
        g1 = mp.expm1(_lgamma1p(a))
        xa1 = mp.expm1(a * mp.log(x))
        eps = mpf(2) ** (-mp.mp.prec - 4)
        u = mpf(1)
        s = mpf(0)
        n = 0
        while True:
            n += 1
            u *= -x / n
            t = u / (a + n)
            s += t
            if n > x and abs(t) < eps * abs(s):
                break
            if n > 10000:
                raise OracleError("R4 series cap")
        num = g1 - xa1 - a * (xa1 + 1) * s
        return num / (g1 + 1)


def _saturated_small(a, x, aphi):
    a, x = mpf(a), mpf(x)
    lam = x / a
    return mp.exp(-aphi) / (mp.sqrt(2 * mp.pi * a) * abs(lam - 1))


def pq(a, x, dps):
    """(P, Q) at dps digits, each relatively accurate."""
    a, x = float(a), float(x)
    if not (a > 0 and x > 0) or math.isinf(a) or math.isinf(x):
        raise OracleError("special input")
    with mp.workdps(dps + 10):
        aphi = _aphi(a, x)
        if aphi > SAT_APHI:
            sm = _saturated_small(a, x, aphi)
            return (sm, 1 - sm) if x < a else (1 - sm, sm)
    if a > A_MAX_EXACT:
        if x == a and a >= 2.0 ** 100:
            with mp.workdps(dps + 10):
                q = mpf(0.5) - 1 / (3 * mp.sqrt(2 * mp.pi * mpf(a)))
                return 1 - q, q
        raise OracleError(f"a={a:g} > {A_MAX_EXACT:g} off saturation: not covered")
    with mp.workdps(dps + _guard(a, x)):
        if a <= 1.5 and x <= 4:
            return +_p_series(a, x), +_q_r4(a, x)
        if x < a:
            try:
                p = _p_hyp(a, x)
            except mp.libmp.NoConvergence as exc:
                raise OracleError(f"hyp1f1 cap a={a} x={x}") from exc
            return p, 1 - p
        q = _q_cf_lentz(a, x)
        return 1 - q, q


# ---- independent forms ------------------------------------------------------
def _p_gauss_cf(a, x):
    """DLMF 8.9.2: gamma(a,x) = x^a e^-x / (a - a x/(a+1 + x/(a+2 - (a+1)x/
    (a+3 + 2x/(a+4 - ...))))), evaluated backward at doubling depth."""
    a, x = mpf(a), mpf(x)
    eps = mpf(2) ** (-mp.mp.prec + 8)

    def cf(n):
        v = a + n
        for j in range(n, 0, -1):
            m = (j + 1) // 2
            k = -(a + m - 1) * x if j % 2 else m * x
            v = (a + j - 1) + k / v
        return 1 / v

    n, prev = 64, None
    while n <= CAP:
        cur = cf(n)
        if prev is not None and abs(cur - prev) < eps * abs(cur):
            return mp.exp(a * mp.log(x) - x - mp.loggamma(a)) * cur
        prev, n = cur, 2 * n
    raise OracleError("Gauss CF cap")


def _q_cf_backward(a, x):
    a, x = mpf(a), mpf(x)
    eps = mpf(2) ** (-mp.mp.prec + 8)

    def cf(n):
        k = x + 2 * n + 1 - a
        for j in range(n, 0, -1):
            k = (x + 2 * j - 1 - a) - j * (j - a) / k
        return k

    n, prev = 64, None
    while n <= CAP:
        cur = cf(n)
        if prev is not None and abs(cur - prev) < eps * abs(cur):
            return mp.exp(a * mp.log(x) - x - mp.loggamma(a)) / cur
        prev, n = cur, 2 * n
    raise OracleError("backward CF cap")


def _q_box_alt(a, x):
    """Q = -expm1(a ln x - lgamma(1+a) + log1p(a*S4)) with e^-x*sum = 1+a*S4,
    S4 = sum_{n>=1} (-x)^n/(n!(a+n)): a second assembly of the box (shares
    S4 and lgamma1p with _q_r4; used only where mpmath is too slow)."""
    a, x = mpf(a), mpf(x)
    with mp.extradps(12):
        eps = mpf(2) ** (-mp.mp.prec - 4)
        u, s, n = mpf(1), mpf(0), 0
        while True:
            n += 1
            u *= -x / n
            t = u / (a + n)
            s += t
            if n > x and abs(t) <= eps * abs(s):
                break
            if n > 10000:
                raise OracleError("box alt series cap")
        return -mp.expm1(a * mp.log(x) - _lgamma1p(a) + mp.log1p(a * s))


def pq2(a, x, dps):
    a, x = float(a), float(x)
    with mp.workdps(dps + 10):
        aphi = _aphi(a, x)
        if aphi > SAT_APHI or a > A_MAX_EXACT:
            return pq(a, x, dps)  # no second form there (documented)
    with mp.workdps(dps + _guard(a, x)):
        if a <= 1.5 and x <= 4:
            if a < 1e-30 and x < 1e-20:
                q = _q_box_alt(a, x)  # mpmath is seconds per call here
            else:
                q = mp.gammainc(mpf(a), mpf(x), mp.inf, regularized=True)
            # Gauss's CF cancels like 1 - P for P near 1 (tiny a): there P
            # comes from the independent Q instead.
            return (+(1 - q) if q < 0.25 else +_p_gauss_cf(a, x)), +q
        if x < a:
            p = _p_gauss_cf(a, x)
            return p, 1 - p
        q = _q_cf_backward(a, x)
        return 1 - q, q


# ---- inverse ----------------------------------------------------------------
def _bits(d: float) -> int:
    import struct
    return struct.unpack("<Q", struct.pack("<d", d))[0]


def _from_bits(u: int) -> float:
    import struct
    return struct.unpack("<d", struct.pack("<Q", u))[0]


_INF_BITS = 0x7FF0000000000000


def _guess(kind, a, s):
    try:
        import scipy.special as sc
        g = sc.gammaincinv(a, s) if kind == "p" else sc.gammainccinv(a, s)
        g = float(g)
        if g > 0 and math.isfinite(g):
            return g
    except Exception:
        pass
    return None


def inv_root(kind, a, s, dps, forward=pq, max_evals=160):
    """Exact root x of P(a,x) = s (kind 'p') or Q(a,x) = s (kind 'q')."""
    a, s = float(a), float(s)
    if not (0 < s < 1 and a > 0):
        raise OracleError("special input")
    # Beyond resolution: the whole transition sits inside half an ulp of a.
    if a >= 2.0 ** 130:
        with mp.workdps(dps + 10):
            # leading order x = a + sqrt(a) * z, z the normal quantile of P
            # z only places the root inside its half-ulp; the offset is
            # < 2^-20 ulp, so a double normal quantile is ample.
            import scipy.special as sc
            z = float(sc.ndtri(s)) if kind == "p" else -float(sc.ndtri(s))
            return mpf(a) + mp.sqrt(mpf(a)) * z
    # Solve on the small side: g(x) increasing in x, root at g = 0.
    flip = s > 0.5
    side = kind if not flip else ("q" if kind == "p" else "p")
    with mp.workdps(dps + 10):
        t = mpf(s) if not flip else 1 - mpf(s)
    extra = 0
    if side == "p" and a < 1:
        extra = int(-math.log10(a)) + 2   # kappa = d ln P / d ln x ~ a
    wdps = dps + 12 + extra
    evals = [0]

    def g(xd):
        if xd == 0.0:
            return mpf(-1)
        if math.isinf(xd):
            return mpf(1)
        evals[0] += 1
        if evals[0] > max_evals:
            raise OracleError("inverse eval cap")
        P, Q = forward(a, xd, wdps)
        with mp.workdps(wdps):
            return (P - t) if side == "p" else (t - Q)

    x0 = _guess(kind, a, s)
    # Newton in log x from the float guess (log of the small side, whose
    # slope is x f / F with x f = x^a e^-x / Gamma(a)): a poor tail guess
    # otherwise costs ~40 bracket evaluations.
    for _ in range(4):
        if x0 is None:
            break
        P, Q = forward(a, x0, wdps)
        evals[0] += 1
        F = P if side == "p" else Q
        if not F > 0:
            break
        with mp.workdps(wdps):
            xf = mp.exp(a * mp.log(x0) - x0 - mp.loggamma(a))
            slope = xf / F if side == "p" else -xf / F
            if slope == 0:
                break
            step = (mp.log(t) - mp.log(F)) / slope
            if abs(step) > 50:
                break
            x1 = float(mpf(x0) * mp.exp(step))
        if not (x1 > 0 and math.isfinite(x1)) or x1 == x0:
            break
        x0 = x1
        if abs(step) < 2.0 ** -55:
            break
    if x0 is None:
        # No float guess (scipy returns 0 below the subnormals): test the
        # bottom of the range first, which settles the zero-result band.
        lo, hi = (0, 1) if g(_from_bits(1)) >= 0 else (1, _INF_BITS)
    else:
        b0 = _bits(x0)
        g0 = g(x0)
        if g0 == 0:
            return mpf(x0)
        step = 1
        if g0 < 0:
            lo, hi = b0, None
            while hi is None:
                nb = min(b0 + step, _INF_BITS)
                if g(_from_bits(nb)) >= 0:
                    hi = nb
                else:
                    lo = nb
                    step *= 4
        else:
            lo, hi = None, b0
            while lo is None:
                nb = max(b0 - step, 0)
                if g(_from_bits(nb)) < 0:
                    lo = nb
                else:
                    hi = nb
                    step *= 4
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if g(_from_bits(mid)) < 0:
            lo = mid
        else:
            hi = mid
    xl, xh = _from_bits(lo), _from_bits(hi)
    if math.isinf(xh):
        raise OracleError("root beyond DBL_MAX")
    if xh < 2.0 ** -1000:
        # Root below 2^-1000: P = x^a/Gamma(1+a) * (1 - O(x)) there, so the
        # deep-small closed form is exact to ~2^-1000 relative (and, unlike a
        # secant across a subnormal ulp, it respects P's x^a curvature).
        with mp.workdps(wdps):
            lnp = mp.log(t) if side == "p" else mp.log1p(-t)
            return mp.exp((lnp + _lgamma1p(mpf(a))) / a)
    gl, gh = g(xl), g(xh)
    with mp.workdps(wdps):
        if gh == 0:
            return mpf(xh)
        return mpf(xl) + (-gl) * (mpf(xh) - mpf(xl)) / (gh - gl)


def gamma_value(fn, a, x, dps, two=False):
    P, Q = (pq2 if two else pq)(a, x, dps)
    return P if fn == "gamma_p" else Q


def gammainv_value(fn, a, s, dps, two=False):
    kind = "p" if fn == "gamma_p_inv" else "q"
    return inv_root(kind, a, s, dps, forward=pq2 if two else pq)


# hunt.py imports every module in families/ when --families is empty.
SEAMS: list = []


def strata():
    return []
