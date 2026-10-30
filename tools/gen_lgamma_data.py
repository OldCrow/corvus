#!/usr/bin/env python3
"""Generate src/lgamma_data.h -- every table the lgamma kernel needs.

Three fits, all Chebyshev interpolation converted to monomial form, and one
Taylor table:

  zone polys   B_c(t) = lgamma(c + t) / t  for c = 1 and c = 2, t in [-1/2, 1/2].
               lgamma's only positive-axis zeros are at x = 1 and x = 2, and
               dividing them out is what makes RELATIVE accuracy there
               possible: the kernel forms t*B(t) with t exact (Sterbenz), so
               the zero is reproduced by construction rather than by
               cancellation. Padded to one degree so the kernel selects
               coefficients per lane and runs a single Horner pass.

  Stirling     psi(w) = x*phi(x) with w = 1/x^2, phi the Stirling remainder
               lgamma(x) - [(x-1/2)log x - x + log(2pi)/2], fitted on
               w in [0, 1/X0^2].

  sinpi        K(v) = log(sin(pi*u)/(pi*u)) / v with v = u^2, u in [-1/2, 1/2].
               The reflection needs log(pi/|sin pi x|); writing it as
               -log|u| - v*K(v) cancels pi analytically and keeps the u -> 0
               limit finite, where forming sin(pi*x) first would underflow.

  band table   lgamma(x0 + r) on (5/2, X0) (#31): x0 = 5/2 + i*h, h = 1/32,
               r = x - x0 exact, |r| <= h/2, and
                 c0 + c1*r + r^2*(c2 + c3*r + ... + c8*r^6),
               c0 = lgamma(x0) and c1 = psi(x0) as double-double, cn =
               psi^(n-1)(x0)/n! in plain double. Taylor, not a fit: the
               coefficients are exact values, so there is nothing to
               converge, and the truncation (the r^9 term) is a closed
               bound. Stored interleaved, one row per node, so a lane's
               gathers stay in one or two cache lines. SELF-CHECK replays
               the kernel's arithmetic, fused and unfused, at bit-stepped
               cell edges (the doctrine's replay rules).

The leading coefficients of the zone and sinpi fits are emitted as
double-double pairs (kLgammaZoneLead*, kLgammaSinLead*): the kernel evaluates
only the tail in plain double, so the tail's rounding is attenuated by the
t^3 (resp. v^2) it is multiplied by. How many leading terms that takes is not
a guess -- SELF-CHECK below replays the kernel's exact arithmetic (double
where the kernel uses double, exact where it uses dd) against mpmath and
refuses to emit a table that misses its budget.

Usage:
    python3 tools/gen_lgamma_data.py > src/lgamma_data.h
"""

import math
import random
import struct
import sys

import mpmath as mp

mp.mp.dps = 60

# --- Design constants -------------------------------------------------------
# Stirling takes over at X0; (5/2, X0) is the band table. X0 = 8 predates the
# table: it was the cheapest point of the old recurrence that met budget, and
# moving it would change Stirling's fit, so it stays.
X0 = mp.mpf(8)
ZONE_LO = mp.mpf("0.5")    # below this, one log shift instead
ZONE_MID = mp.mpf("1.5")   # centre-1 / centre-2 split
ZONE_HI = mp.mpf("2.5")    # above this, the recurrence

N_NODES = 160              # EVEN on purpose: an odd count puts a node at t = 0,
                           # where mpmath's 1 + t rounds to 1 and B(t) collapses
                           # to 0/0 -- a single bad node, and the Chebyshev
                           # coefficients then plateau instead of decaying.
TAIL_CUT = mp.mpf("1e-21")
ZONE_LEAD = 3              # dd leading coefficients of B(t)
SIN_LEAD = 2               # dd leading coefficients of K(v)
TAB_INV_H = 32             # band table: nodes per unit x
TAB_DEG = 8                # r^8 is the last term
# Relative to lgamma on the band. Replay measures ~2^-64.8; the tail's
# coefficient rounding, not the truncation, sets that floor.
TAB_TARGET = mp.mpf(2) ** -63.5

# Gates, in relative error of the replayed kernel evaluation. The floor is the
# tail Horner's own rounding, not the fit residual (the mpmath-side fit is
# orders of magnitude finer). ~1 ULP of a double is 1.1e-16.
ZONE_TARGET = 3.0e-17
STIR_TARGET = 1.0e-18
SIN_TARGET = 4.0e-17       # absolute, on log(sin(pi u)/(pi u))

SEED = 20260725


# --- Chebyshev machinery ----------------------------------------------------
def cheb_coeffs(f, lo, hi, n_nodes=N_NODES):
    c = (mp.mpf(lo) + mp.mpf(hi)) / 2
    h = (mp.mpf(hi) - mp.mpf(lo)) / 2
    nodes = [mp.cos(mp.pi * (j + mp.mpf(1) / 2) / n_nodes) for j in range(n_nodes)]
    vals = [f(c + h * s) for s in nodes]
    out = []
    for k in range(n_nodes):
        # A list, not a generator: mpmath's fsum silently mis-sums a generator.
        acc = mp.fsum([vals[j] * mp.cos(mp.pi * k * (j + mp.mpf(1) / 2) / n_nodes)
                       for j in range(n_nodes)])
        a = 2 * acc / n_nodes
        out.append(a / 2 if k == 0 else a)
    return out, c, h


def truncate(coeffs):
    d = len(coeffs) - 1
    while d > 0 and abs(coeffs[d]) < TAIL_CUT:
        d -= 1
    return coeffs[: d + 1]


def cheb_to_monomial(coeffs, c, h):
    """Chebyshev series in s = (t - c)/h -> monomial coefficients in t."""
    # Work in s first, then compose with the affine map.
    t_prev = [mp.mpf(1)]
    t_cur = [mp.mpf(0), mp.mpf(1)]
    mono_s = [mp.mpf(0)] * len(coeffs)

    def add(poly, w):
        for i, p in enumerate(poly):
            mono_s[i] += w * p

    add(t_prev, coeffs[0])
    if len(coeffs) > 1:
        add(t_cur, coeffs[1])
    for k in range(2, len(coeffs)):
        t_next = [mp.mpf(0)] + [2 * x for x in t_cur]
        for i, x in enumerate(t_prev):
            t_next[i] -= x
        add(t_next, coeffs[k])
        t_prev, t_cur = t_cur, t_next

    # s = (t - c)/h: expand sum_k a_k ((t-c)/h)^k into powers of t.
    n = len(mono_s)
    mono_t = [mp.mpf(0)] * n
    for k in range(n - 1, -1, -1):  # Horner in (t - c)/h
        new = [mp.mpf(0)] * n
        for i in range(n - 1):
            new[i + 1] += mono_t[i] / h
            new[i] -= mono_t[i] * c / h
        new[0] += mono_s[k]
        mono_t = new
    return mono_t


# --- doubles ----------------------------------------------------------------
def rd(x):
    return float(x)


def dd_split(x):
    hi = rd(x)
    return hi, rd(mp.mpf(x) - mp.mpf(hi))


def horner_d(coefs, x):
    """Horner in plain double -- one rounding per fused step, as the kernel does."""
    acc = 0.0
    for cf in reversed(coefs):
        acc = rd(mp.mpf(acc) * mp.mpf(x) + mp.mpf(cf))
    return acc


# --- the functions being fitted ---------------------------------------------
def zone_B(centre, t):
    """lgamma(centre + t)/t, with the removable singularity at t = 0 filled."""
    if t == 0:
        return mp.digamma(centre)
    return mp.loggamma(centre + t) / t


def stir_psi(w):
    """x*phi(x) with w = 1/x^2; phi is the Stirling remainder."""
    if w == 0:
        return mp.mpf(1) / 12
    x = 1 / mp.sqrt(w)
    phi = mp.loggamma(x) - ((x - mp.mpf(1) / 2) * mp.log(x) - x
                            + mp.log(2 * mp.pi) / 2)
    return x * phi


def sin_K(v):
    """log(sin(pi u)/(pi u))/v with v = u^2; K(0) = -pi^2/6."""
    if v == 0:
        return -(mp.pi ** 2) / 6
    u = mp.sqrt(v)
    return mp.log(mp.sin(mp.pi * u) / (mp.pi * u)) / v


# --- fit + replay -----------------------------------------------------------
def fit_lead_tail(f, lo, hi, n_lead):
    """Fit f, then split it into n_lead dd coefficients plus a double tail."""
    coeffs, c, h = cheb_coeffs(f, lo, hi)
    mono = cheb_to_monomial(truncate(coeffs), c, h)
    lead = [dd_split(mono[k]) for k in range(n_lead)]
    tail = [rd(mono[k]) for k in range(n_lead, len(mono))]
    return lead, tail


def replay_lead_tail(lead, tail, t):
    """Exactly what the kernel computes: double tail Horner, dd from there up.

    dd steps are modelled as exact -- they carry ~2^-105, three orders below
    anything else here, so what this measures is the double tail's rounding
    plus the fit residual, which is the pair the gate is about.
    """
    s = horner_d(tail, t)
    acc = mp.mpf(rd(mp.mpf(s) * mp.mpf(t)))       # the one rounded product
    acc += mp.mpf(lead[-1][0]) + mp.mpf(lead[-1][1])
    for k in range(len(lead) - 2, -1, -1):
        acc = acc * mp.mpf(t) + mp.mpf(lead[k][0]) + mp.mpf(lead[k][1])
    return acc


def check_zone(centre, lead, tail, lo, hi, rng, n=4000):
    worst = mp.mpf(0)
    worst_t = 0.0
    for _ in range(n):
        t = rd(rng.uniform(float(lo), float(hi)))
        if t == 0.0:
            continue
        # No final rounding here: that one is unavoidable and belongs to the
        # ULP test, not to the table's budget. What is measured is how much
        # error the approximation hands to it.
        got = mp.mpf(t) * replay_lead_tail(lead, tail, t)
        want = mp.loggamma(mp.mpf(centre) + mp.mpf(t))
        if want == 0:
            continue
        rel = abs((got - want) / want)
        if rel > worst:
            worst, worst_t = rel, t
    return float(worst), worst_t


def check_stir(coefs, rng, n=4000):
    """Replay the whole Stirling assembly: phi in double, the rest exact."""
    worst = mp.mpf(0)
    worst_x = 0.0
    for _ in range(n):
        # Sweep 1/x uniformly so the X0 end -- the tightest one -- is dense.
        x = rd(1.0 / rng.uniform(1e-4, float(1 / X0)))
        w = rd(1.0 / rd(mp.mpf(x) * mp.mpf(x)))
        phi = rd(mp.mpf(horner_d(coefs, w)) / mp.mpf(x))
        xm = mp.mpf(x)
        got = ((xm - mp.mpf("0.5")) * mp.log(xm) - xm
               + mp.log(2 * mp.pi) / 2 + mp.mpf(phi))
        want = mp.loggamma(xm)
        rel = abs((got - want) / want)  # pre-rounding, as in check_zone
        if rel > worst:
            worst, worst_x = rel, x
    return float(worst), worst_x


def check_sin(lead, tail, rng, n=4000):
    worst = mp.mpf(0)
    worst_u = 0.0
    for _ in range(n):
        u = rd(rng.uniform(-0.5, 0.5))
        if u == 0.0:
            continue
        v = rd(mp.mpf(u) * mp.mpf(u))
        got = mp.mpf(v) * replay_lead_tail(lead, tail, v)
        want = mp.log(mp.sin(mp.pi * mp.mpf(u)) / (mp.pi * mp.mpf(u)))
        err = abs(got - want)
        if err > worst:
            worst, worst_u = err, u
    return float(worst), worst_u


def max_finite_arg():
    """Largest double x with lgamma(x) finite. Above it the kernel returns inf."""
    dbl_max = mp.mpf(float.fromhex("0x1.fffffffffffffp+1023"))
    # Bisect rather than Newton: lgamma is ~x log x here, so a secant step from
    # any plausible guess overshoots into the overflow region and the solver
    # gives up.
    lo, hi = mp.mpf("1e300"), mp.mpf("1e307")
    for _ in range(400):
        mid = (lo + hi) / 2
        if mp.loggamma(mid) > dbl_max:
            hi = mid
        else:
            lo = mid
    xf = rd(lo)
    while mp.loggamma(mp.mpf(xf)) > dbl_max:  # step back onto a finite double
        bits = struct.unpack("<Q", struct.pack("<d", xf))[0] - 1
        xf = struct.unpack("<d", struct.pack("<Q", bits))[0]
    return xf


# --- emit -------------------------------------------------------------------
def tab_row(i):
    """[c0.hi, c0.lo, c1.hi, c1.lo, c8, c7, ..., c2] at node i: Horner order
    for the tail, so a row reads front to back as the kernel uses it."""
    x0 = mp.mpf(5) / 2 + mp.mpf(i) / TAB_INV_H
    c0 = dd_split(mp.loggamma(x0))
    c1 = dd_split(mp.digamma(x0))
    tail = [rd(mp.polygamma(n - 1, x0) / mp.factorial(n))
            for n in range(TAB_DEG, 1, -1)]
    return [c0[0], c0[1], c1[0], c1[1]] + tail


def replay_tab(row, x, fused):
    """The kernel's table evaluation at double x, as an mpf: Horner, r^2 and
    the tail in double (fused or not), the c0 + c1*r assembly exact (its dd
    rounding is ~2^-104, far below this budget)."""
    i = int(mp.nint((mp.mpf(x) - mp.mpf(5) / 2) * TAB_INV_H))
    x0 = float(mp.mpf(5) / 2 + mp.mpf(i) / TAB_INV_H)
    r = x - x0  # exact
    assert mp.mpf(r) == mp.mpf(x) - mp.mpf(x0)
    p = row[i][4]
    for c in row[i][5:]:
        p = rd(mp.mpf(p) * r + c) if fused else p * r + c
    tail = p * (r * r)
    return (mp.mpf(row[i][0]) + row[i][1]
            + (mp.mpf(row[i][2]) + row[i][3]) * r + tail)


def check_tab(rows, rng):
    """Bit-stepped points on both sides of every cell edge and around every
    node, plus random points per cell, both FMA models. One worst-case figure,
    plus the per-arm counts. Coefficients are recomputed at twice the digits
    first: a generator that disagrees with itself emits nothing."""
    with mp.workdps(2 * mp.mp.dps):
        for i in (0, 1, len(rows) // 2, len(rows) - 2, len(rows) - 1):
            again = tab_row(i)
            if again != rows[i]:
                return mp.inf, 0.0, 0, 0, 0
    n_edge = n_node = n_rand = 0
    worst, wx = mp.mpf(0), 0.0
    h = 1.0 / TAB_INV_H
    lo, hi = float(ZONE_HI), float(X0)

    def probe(x):
        nonlocal worst, wx
        if not lo < x < hi:
            return 0
        v = mp.loggamma(mp.mpf(x))
        for fused in (True, False):
            e = abs(replay_tab(rows, x, fused) - v) / abs(v)
            if e > worst:
                worst, wx = e, x
        return 1

    for i in range(len(rows)):
        x0 = lo + i * h
        for edge in (x0 - h / 2, x0 + h / 2):
            a = edge
            b = edge
            for _ in range(4):
                n_edge += probe(a) + probe(b)
                a = math.nextafter(a, 0.0)
                b = math.nextafter(b, 9.0)
        n_node += probe(math.nextafter(x0, 0.0)) + probe(math.nextafter(x0, 9.0))
        for _ in range(8):
            n_rand += probe(x0 + rng.uniform(-h / 2, h / 2))
    return worst, wx, n_edge, n_node, n_rand


def hexes(name, vals, per_line=1):
    print(f"inline constexpr double {name} = {{")
    for v in vals:
        print(f"    {float.hex(v)},")
    print("};")


def main():
    rng = random.Random(SEED)
    rc = 0

    # Zone polynomials. Centre 1 also serves x in (0, 1/2) via
    # lgamma(x) = lgamma(1 + x) - log x, where t = x is exact, so its argument
    # range is the union [-1/2, 1/2]; centre 2 likewise covers both the direct
    # band and every landing point of the recurrence, which is always in
    # (3/2, 5/2].
    zones = []
    for centre in (1, 2):
        lead, tail = fit_lead_tail(lambda t, c=centre: zone_B(c, t),
                                   -0.5, 0.5, ZONE_LEAD)
        err, wt = check_zone(centre, lead, tail, -0.5, 0.5, rng)
        print(f"zone centre {centre}: tail degree {len(tail) - 1 + ZONE_LEAD}, "
              f"max rel err {err:.3e} at t={wt:.17g}", file=sys.stderr)
        if err > ZONE_TARGET:
            print(f"FAILED: zone {centre} exceeds {ZONE_TARGET:g}", file=sys.stderr)
            rc = 1
        zones.append((lead, tail))
    n_zone = max(len(t) for _, t in zones)

    # Stirling remainder.
    st_coeffs, c, h = cheb_coeffs(stir_psi, 0, 1 / (X0 * X0), 81)
    st_mono = [rd(v) for v in cheb_to_monomial(truncate(st_coeffs), c, h)]
    err, wx = check_stir(st_mono, rng)
    print(f"stirling X0={float(X0)}: psi degree {len(st_mono) - 1}, "
          f"max rel err {err:.3e} at x={wx:.17g}", file=sys.stderr)
    if err > STIR_TARGET:
        print(f"FAILED: stirling exceeds {STIR_TARGET:g}", file=sys.stderr)
        rc = 1

    # sinpi log-correction.
    sin_lead, sin_tail = fit_lead_tail(sin_K, 0, 0.25, SIN_LEAD)
    err, wu = check_sin(sin_lead, sin_tail, rng)
    print(f"sinpi: K degree {len(sin_tail) - 1 + SIN_LEAD}, "
          f"max abs err {err:.3e} at u={wu:.17g}", file=sys.stderr)
    if err > SIN_TARGET:
        print(f"FAILED: sinpi exceeds {SIN_TARGET:g}", file=sys.stderr)
        rc = 1

    # Band table.
    n_tab = int((X0 - ZONE_HI) * TAB_INV_H) + 1
    tab = [tab_row(i) for i in range(n_tab)]
    err, wx, n_edge, n_node, n_rand = check_tab(tab, rng)
    print(f"band table h=1/{TAB_INV_H} D={TAB_DEG}: {n_tab} nodes, {n_edge} edge, "
          f"{n_node} node, {n_rand} random points, fused and unfused: "
          f"max rel err {err:.3e} (2^{float(mp.log(err, 2)):.2f}) at x={wx:.17g}",
          file=sys.stderr)
    if err > TAB_TARGET:
        print(f"FAILED: band table exceeds 2^{float(mp.log(TAB_TARGET, 2)):.1f}",
              file=sys.stderr)
        rc = 1

    if rc:
        return rc

    hl2pi = dd_split(mp.log(2 * mp.pi) / 2)
    max_arg = max_finite_arg()

    print("// Auto-generated by tools/gen_lgamma_data.py. DO NOT EDIT.")
    print("// Fits and their error budgets are documented in the generator;")
    print("// the kernel that consumes them is src/lgamma-inl.h.")
    print("#ifndef CORVUS_LGAMMA_DATA_H_")
    print("#define CORVUS_LGAMMA_DATA_H_")
    print()
    print("namespace corvus::detail {")
    print()
    print("// Region boundaries. (0, LO) shifts by one log; [LO, MID) and")
    print("// [MID, HI] are the zero-centred polynomials; (HI, X0) is the band")
    print("// table; [X0, inf) is Stirling.")
    print(f"inline constexpr double kLgammaZoneLo = {float.hex(float(ZONE_LO))};")
    print(f"inline constexpr double kLgammaZoneMid = {float.hex(float(ZONE_MID))};")
    print(f"inline constexpr double kLgammaZoneHi = {float.hex(float(ZONE_HI))};")
    print(f"inline constexpr double kLgammaX0 = {float.hex(float(X0))};")
    print()
    print("// Largest argument with a finite result; above it lgamma overflows.")
    print(f"inline constexpr double kLgammaMaxArg = {float.hex(max_arg)};")
    print()
    print("// B(t) = L0 + t*(L1 + t*(L2 + t*S(t))), t = x - centre, index 0 is")
    print("// centre 1 and index 1 is centre 2. L* are double-double; S is the")
    print("// plain-double tail, zero-padded to one degree for a single Horner.")
    print(f"inline constexpr int kLgammaZoneNCoef = {n_zone};")
    for part, idx in (("Hi", 0), ("Lo", 1)):
        print(f"inline constexpr double kLgammaZoneLead{part}[2][{ZONE_LEAD}] = {{")
        for lead, _ in zones:
            print("    {" + ", ".join(float.hex(p[idx]) for p in lead) + "},")
        print("};")
    print(f"inline constexpr double kLgammaZoneCoef[2][{n_zone}] = {{")
    for centre, (_, tail) in zip((1, 2), zones):
        print(f"    {{  // centre {centre}")
        for cf in tail + [0.0] * (n_zone - len(tail)):
            print(f"        {float.hex(cf)},")
        print("    },")
    print("};")
    print()
    print("// phi(x) = psi(w)/x with w = 1/x^2: the Stirling remainder. Plain")
    print("// double is enough -- phi is ~1/(12x), so even a full ULP of it is")
    print("// ~2^-60 of lgamma(X0), the tightest point in the region.")
    print(f"inline constexpr int kLgammaStirNCoef = {len(st_mono)};")
    hexes(f"kLgammaStirCoef[{len(st_mono)}]", st_mono)
    print(f"inline constexpr double kLgammaHalfLog2PiHi = {float.hex(hl2pi[0])};")
    print(f"inline constexpr double kLgammaHalfLog2PiLo = {float.hex(hl2pi[1])};")
    print()
    print("// Band table: lgamma(x0 + r), x0 = HI + i/kLgammaTabInvH, one row per")
    print("// node: c0.hi, c0.lo, c1.hi, c1.lo, then the plain-double tail")
    print(f"// c{TAB_DEG} .. c2 in Horner order.")
    print(f"inline constexpr int kLgammaTabN = {n_tab};")
    print(f"inline constexpr int kLgammaTabStride = {len(tab[0])};")
    print(f"inline constexpr double kLgammaTabInvH = {float.hex(float(TAB_INV_H))};")
    print(f"alignas(64) inline constexpr double kLgammaTab[{n_tab * len(tab[0])}] = {{")
    for row in tab:
        print("    " + ", ".join(float.hex(v) for v in row) + ",")
    print("};")
    print()
    print("// log(sin(pi u)/(pi u)) = v*(M0 + v*(M1 + v*T(v))), v = u^2.")
    print(f"inline constexpr int kLgammaSinNCoef = {len(sin_tail)};")
    print("inline constexpr double kLgammaSinLeadHi[" + str(SIN_LEAD) + "] = {"
          + ", ".join(float.hex(p[0]) for p in sin_lead) + "};")
    print("inline constexpr double kLgammaSinLeadLo[" + str(SIN_LEAD) + "] = {"
          + ", ".join(float.hex(p[1]) for p in sin_lead) + "};")
    hexes(f"kLgammaSinCoef[{len(sin_tail)}]", sin_tail)
    print()
    print("}  // namespace corvus::detail")
    print()
    print("#endif  // CORVUS_LGAMMA_DATA_H_")
    return 0


if __name__ == "__main__":
    sys.exit(main())
