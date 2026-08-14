#!/usr/bin/env python3
"""Generate src/erfinv_data.h -- every table the erfinv/erfcinv kernel needs.

Both public functions route onto two shared cores:

  Core C (central, |x| <= ~0.4769): x = y*Pc(y^2), a direct polynomial fit,
  NO Newton/Halley step -- the central condition number is ~1, so refining
  against corvus's own (already ~1-ULP) erfc would pass that error straight
  through with no benefit. Pc is fitted the way lgamma's zone polynomials
  are: Chebyshev on v = y^2 in [0, 1/4], converted to monomial, then split
  into a few dd LEADING coefficients (the low-degree, dominant terms) plus a
  plain-double tail (the high-degree terms, whose rounding is attenuated by
  the v^n_lead they are multiplied by).

  Core T (tail, s in (0, 1/2), i.e. x in (~0.4769, ~27.217)): a cheap seed
  polynomial (2^-19 relative, three intervals selected by t = sqrt(-log s))
  followed by ONE dd Halley step. Two Halley formulations, selected by the
  SAME threshold the erfc kernel already uses (x = 6, i.e. t >= kErfinvTFar):
    mid  (t <  kErfinvTFar): residual space, f(x) = erfc(x) - s, using the
         erfc kernel's own pre-rounding compensated dd (ErfcCoreDd) so the
         step is not floored by double-rounding the public erfc/erf.
    far  (t >= kErfinvTFar): log space, F(x) = log erfc(x) - log s, written
         via corvus's own tail model erfc(x) = e^{-x^2} G(1/x)/x (REUSING
         erfc_tail_data.h -- no new tail fit needed) so it needs no
         exponential at all and stays accurate for arbitrarily small
         (subnormal) s, where residual space would underflow to zero.
  Both reduce to a Halley step x1 = x0 + delta with a closed form (derived
  from f'' = -2x*f' for the plain-erfc form, and F'' = -2x*F' - F'^2 for the
  log form -- see the derivation comment in src/erfinv-inl.h):
    mid:  delta = -f / (fp + x0*f)
    far:  delta = -2*F / (2*Fp + F*(2*x0 + Fp))

SELF-CHECK (house rule, same as gen_lgamma_data.py/gen_exp_table.py): this
generator does not just fit curves, it REPLAYS the kernel's actual visible
double-precision arithmetic (t = sqrt(w.hi), the seed Horner, x0 = t*seed,
the Halley formulas above) against an mpmath oracle, and refuses to emit a
table unless:
  * each seed polynomial alone is accurate to 2^-19 relative (x0/t vs the
    true x/t), and
  * the full seed+Halley replay lands within 2^-56 relative of the true
    root (pre-rounding -- the final round to a double is the ULP test's
    job, not this budget's).
Per the project's convention (see gen_lgamma_data.py), already-audited
sub-results are modeled as their IDEAL value in this replay: ErfcCoreDd
stands in for mpmath's erfc, and the fitted G stands in for the true
G(u) = a*e^{a^2}*erfc(a). Both are already accuracy-audited elsewhere (erfc
ULP gates; gen_erfc_tail_poly.py's own self-check); what THIS check exists
to catch is a transcription or formula bug in the new seed/Halley code, not
to re-litigate bounds that already have their own gate. dd steps (LogDdAny,
the DdAdd/DdMulD assembly) are likewise modeled as exact -- they carry
~2^-104, three-plus orders below the 2^-56 budget here.

Two documented mpmath traps, both from gen_lgamma_data.py, apply again:
fsum needs a LIST (a generator silently mis-sums), and the Chebyshev node
count must be EVEN (an odd count puts a node at the domain's exact center,
which is fine here since none of these fits are centered on a singularity,
but the node count is kept even anyway for consistency).

Usage:
    python3 tools/gen_erfinv_data.py > src/erfinv_data.h
"""

import random
import sys

import mpmath as mp

mp.mp.dps = 60

# --- Design constants --------------------------------------------------
Y_SPLIT = mp.mpf("0.5")     # |y| <= Y_SPLIT routes to the central core C
C_LEAD = 3                  # dd leading coefficients of Pc(v)
N_NODES_C = 64               # EVEN (see module docstring)
N_NODES_SEED = 64
TAIL_CUT = mp.mpf("1e-20")

# T region boundaries, all derived, not guessed.
T_LO = mp.sqrt(mp.log(2))                       # t at s = 1/2
T_SPLIT = mp.mpf(2)                             # seed interval 0/1 split
T_FAR = mp.sqrt(-mp.log(mp.erfc(mp.mpf(6))))    # t at x = 6: mid/far split,
                                                 # and reuses erfc's own
                                                 # core/tail boundary.
T_HI = mp.sqrt(mp.mpf(1074) * mp.log(2))        # t at s = 2^-1074 (the
                                                 # smallest subnormal double)

# Gates (relative error), no margin beyond what the design calls for.
C_TARGET = mp.mpf(2) ** -55
SEED_TARGET = mp.mpf(2) ** -19
HALLEY_TARGET = mp.mpf(2) ** -56

SEED = 20260725


# --- Chebyshev machinery (reused verbatim from gen_lgamma_data.py) -----
def cheb_coeffs(f, lo, hi, n_nodes):
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


def truncate(coeffs, cut=TAIL_CUT):
    d = len(coeffs) - 1
    while d > 0 and abs(coeffs[d]) < cut:
        d -= 1
    return coeffs[: d + 1]


def cheb_to_monomial(coeffs, c, h):
    """Chebyshev series in s = (t - c)/h -> monomial coefficients in t."""
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


# --- doubles -------------------------------------------------------------
def rd(x):
    return float(x)


def dd_split(x):
    hi = rd(x)
    return hi, rd(mp.mpf(x) - mp.mpf(hi))


def horner_d(coefs, x):
    """Horner in plain double -- one rounding per fused step."""
    acc = 0.0
    for cf in reversed(coefs):
        acc = rd(mp.mpf(acc) * mp.mpf(x) + mp.mpf(cf))
    return acc


def fit_lead_tail(f, lo, hi, n_lead, n_nodes):
    coeffs, c, h = cheb_coeffs(f, lo, hi, n_nodes)
    mono = cheb_to_monomial(truncate(coeffs), c, h)
    lead = [dd_split(mono[k]) for k in range(n_lead)]
    tail = [rd(mono[k]) for k in range(n_lead, len(mono))]
    return lead, tail


def replay_lead_tail(lead, tail, v):
    """Same shape as gen_lgamma_data.py's replay_lead_tail: double tail
    Horner, dd from there up modeled as exact (see module docstring)."""
    s = horner_d(tail, v)
    acc = mp.mpf(rd(mp.mpf(s) * mp.mpf(v)))
    acc += mp.mpf(lead[-1][0]) + mp.mpf(lead[-1][1])
    for k in range(len(lead) - 2, -1, -1):
        acc = acc * mp.mpf(v) + mp.mpf(lead[k][0]) + mp.mpf(lead[k][1])
    return acc


# --- erfcinv oracle (mpmath has erfinv but not erfcinv) -------------------
def erfcinv_mp(s):
    """x with erfc(x) = s, s in (0, 1/2), via root-find in log space."""
    s = mp.mpf(s)
    w = -mp.log(s)
    x0 = mp.sqrt(max(w - mp.log(mp.pi * w) / 2, mp.mpf(1) / 4))
    return mp.findroot(lambda x: mp.log(mp.erfc(x)) - mp.log(s), x0)


def tail_G(a):
    a = mp.mpf(a)
    return a * mp.e ** (a * a) * mp.erfc(a)


# --- Core C: Pc(v) = erfinv(y)/y, v = y^2, y = sqrt(v) --------------------
def central_f(v):
    if v == 0:
        return 2 / mp.sqrt(mp.pi)  # erfinv'(0)
    y = mp.sqrt(v)
    return mp.erfinv(y) / y


def check_central(lead, tail, rng, n=6000):
    worst, worst_y = mp.mpf(0), 0.0
    for _ in range(n):
        y = rd(rng.uniform(-0.5, 0.5))
        if y == 0.0:
            continue
        v = rd(mp.mpf(y) * mp.mpf(y))
        got = mp.mpf(y) * replay_lead_tail(lead, tail, v)
        want = mp.erfinv(mp.mpf(y))
        if want == 0:
            continue
        rel = abs((got - want) / want)
        if rel > worst:
            worst, worst_y = rel, y
    return worst, worst_y


# --- Core T seeds: x/t as a function of t (mid0/mid1) or u = 1/t (far) ---
def x_over_t(t):
    t = mp.mpf(t)
    return erfcinv_mp(mp.e ** (-t * t)) / t


def fit_seed(lo, hi, in_u, n_nodes=N_NODES_SEED):
    f = (lambda u: x_over_t(1 / mp.mpf(u))) if in_u else x_over_t
    coeffs, c, h = cheb_coeffs(f, lo, hi, n_nodes)
    # Seeds only need 2^-19 (one dd Halley step corrects the rest), so cut
    # much more aggressively than the fits that carry the accuracy claim --
    # a cheaper seed is a cheaper per-lane Horner in the kernel.
    mono = cheb_to_monomial(truncate(coeffs, cut=mp.mpf("3e-7")), c, h)
    return [rd(m) for m in mono]


def check_seed(coefs, lo, hi, in_u, rng, n=800):
    worst, worst_v = mp.mpf(0), 0.0
    for _ in range(n):
        v = rd(rng.uniform(float(lo), float(hi)))
        t = (1 / mp.mpf(v)) if in_u else mp.mpf(v)
        want = x_over_t(t)
        got = horner_d(coefs, v)
        rel = abs((mp.mpf(got) - want) / want)
        if rel > worst:
            worst, worst_v = rel, v
    return worst, worst_v


# --- Full seed+Halley replay, exactly the kernel's arithmetic -------------
def replay_halley(s):
    """Returns (x0, delta) -- NOT x0+delta rounded to a double. The final
    x1 = fl(x0 + delta) is the kernel's one unavoidable rounding (same
    convention as gen_lgamma_data.py's check_zone): it belongs to the ULP
    test's budget, not to this pre-rounding self-check.
    """
    s_mp = mp.mpf(s)
    w_true = -mp.log(s_mp)
    w_hi = rd(w_true)
    t = rd(mp.sqrt(mp.mpf(w_hi)))

    if t < float(T_SPLIT):
        coefs, var = SEED0, t
    elif t < float(T_FAR):
        coefs, var = SEED1, t
    else:
        coefs, var = SEED2, rd(1.0 / t)
    seed_val = horner_d(coefs, var)
    x0 = rd(t * seed_val)

    if t < float(T_FAR):
        # mid: residual space against the IDEAL erfc (ErfcCoreDd modeled
        # exact -- its own ~1 ULP approximation quality is a separate,
        # already-audited bound; see module docstring).
        f = rd(mp.erfc(mp.mpf(x0)) - s_mp)
        fp = rd(-2 / mp.sqrt(mp.pi) * mp.e ** (-mp.mpf(x0) * mp.mpf(x0)))
        denom = rd(fp + rd(x0 * f))
        delta = rd(-f / denom)
    else:
        # far: log space, reusing the IDEAL tail model G (also
        # already-audited via gen_erfc_tail_poly.py's own self-check).
        x0sq_true = mp.mpf(x0) * mp.mpf(x0)  # exact at this precision
        g_val = tail_G(x0)
        log_term = rd(mp.log(mp.mpf(x0)) - mp.log(g_val))
        F = rd(w_true - x0sq_true - mp.mpf(log_term))
        Fp = rd(-2 * mp.mpf(x0) / (mp.sqrt(mp.pi) * g_val))
        denom = rd(2 * Fp + F * rd(2 * x0 + Fp))
        delta = rd(-2 * F / denom)
    return x0, delta


def check_halley(rng, n=20000):
    worst, worst_s = mp.mpf(0), 0.0

    def probe(s):
        nonlocal worst, worst_s
        s = rd(s)
        if s <= 0.0 or s >= 0.5:
            return
        x0, delta = replay_halley(s)
        x1 = mp.mpf(x0) + mp.mpf(delta)  # pre-rounding: no final fl()
        want = erfcinv_mp(s)
        rel = abs((x1 - want) / want)
        if rel > worst:
            worst, worst_s = rel, s

    # Log-spaced s down to the smallest subnormal double, so the far branch
    # (the only one that ever sees a subnormal s) is densely covered.
    for _ in range(n):
        probe(mp.mpf(2) ** rng.uniform(-1074.0, -1.0))

    # Bit/value neighbourhoods of every region boundary: t = T_SPLIT (seed
    # interval 0/1) and t = T_FAR (seed interval 1/2 AND the mid/far Halley
    # split), where a one-ULP shift in t can route two adjacent lanes to
    # different seeds or different Halley formulas entirely.
    for t_b in (T_SPLIT, T_FAR):
        for k in range(-2048, 2049):
            probe(mp.e ** (-(t_b + mp.mpf(k) * mp.mpf("1e-6")) ** 2))
    # s itself right at the smallest/largest ends of T's domain.
    for k in range(-64, 65):
        probe(mp.mpf(2) ** (-1074 + k * 0.001))
        probe(mp.mpf("0.5") - mp.mpf(2) ** -52 * (k + 65))

    return worst, worst_s


# --- emit ------------------------------------------------------------------
def hexes(name, vals):
    print(f"inline constexpr double {name} = {{")
    for v in vals:
        print(f"    {float.hex(v)},")
    print("};")


SEED0 = SEED1 = SEED2 = None  # populated in main(), read by replay_halley


def main():
    global SEED0, SEED1, SEED2
    rng = random.Random(SEED)
    rc = 0

    # --- Core C ---
    lead, tail = fit_lead_tail(central_f, mp.mpf(0), mp.mpf("0.25"), C_LEAD,
                               N_NODES_C)
    err, worst_y = check_central(lead, tail, rng)
    print(f"central: tail degree {len(tail) - 1 + C_LEAD}, "
          f"max rel err {float(err):.3e} at y={worst_y:.17g}", file=sys.stderr)
    if err > C_TARGET:
        print(f"FAILED: central exceeds {float(C_TARGET):g}", file=sys.stderr)
        rc = 1

    # --- Core T seeds ---
    SEED0 = fit_seed(T_LO, T_SPLIT, in_u=False)
    SEED1 = fit_seed(T_SPLIT, T_FAR, in_u=False)
    SEED2 = fit_seed(1 / T_HI, 1 / T_FAR, in_u=True)
    n_seed = max(len(SEED0), len(SEED1), len(SEED2))
    SEED0 += [0.0] * (n_seed - len(SEED0))
    SEED1 += [0.0] * (n_seed - len(SEED1))
    SEED2 += [0.0] * (n_seed - len(SEED2))

    for name, coefs, lo, hi, in_u in (
        ("seed0 [t_lo,2]", SEED0, T_LO, T_SPLIT, False),
        ("seed1 [2,t_far]", SEED1, T_SPLIT, T_FAR, False),
        ("seed2 [1/t_hi,1/t_far] (u)", SEED2, 1 / T_HI, 1 / T_FAR, True),
    ):
        d = next(i for i in range(len(coefs) - 1, -1, -1) if coefs[i] != 0.0)
        err, worst_v = check_seed(coefs, lo, hi, in_u, rng)
        print(f"{name}: degree {d}, max rel err {float(err):.3e} "
              f"at v={worst_v:.17g}", file=sys.stderr)
        if err > SEED_TARGET:
            print(f"FAILED: {name} exceeds {float(SEED_TARGET):g}",
                  file=sys.stderr)
            rc = 1

    # --- Full seed+Halley replay ---
    err, worst_s = check_halley(rng)
    print(f"halley end-to-end: max rel err {float(err):.3e} at s={worst_s:.17g}",
          file=sys.stderr)
    if err > HALLEY_TARGET:
        print(f"FAILED: halley replay exceeds {float(HALLEY_TARGET):g}",
              file=sys.stderr)
        rc = 1

    if rc:
        return rc

    print("// Auto-generated by tools/gen_erfinv_data.py. DO NOT EDIT.")
    print("// Fits and their error budgets are documented in the generator;")
    print("// the kernel that consumes them is src/erfinv-inl.h.")
    print("#ifndef CORVUS_ERFINV_DATA_H_")
    print("#define CORVUS_ERFINV_DATA_H_")
    print()
    print("namespace corvus::detail {")
    print()
    print("// |y| <= kErfinvYSplit routes erfinv to the central core C.")
    print(f"inline constexpr double kErfinvYSplit = {float.hex(float(Y_SPLIT))};")
    print()
    print("// Central: Pc(v) = L0 + v*(L1 + v*(L2 + v*S(v))), v = y^2. L* dd,")
    print("// S the plain-double tail, single Horner pass.")
    print(f"inline constexpr int kErfinvCLead = {C_LEAD};")
    print(f"inline constexpr int kErfinvCNCoef = {len(tail)};")
    print(f"inline constexpr double kErfinvCLeadHi[{C_LEAD}] = {{"
          + ", ".join(float.hex(p[0]) for p in lead) + "};")
    print(f"inline constexpr double kErfinvCLeadLo[{C_LEAD}] = {{"
          + ", ".join(float.hex(p[1]) for p in lead) + "};")
    hexes(f"kErfinvCCoef[{len(tail)}]", tail)
    print()
    print("// T region boundaries (t = sqrt(-log s)). kErfinvTSplit divides")
    print("// the two mid seed intervals; kErfinvTFar is both the far/mid")
    print("// Halley split AND reuses erfc's own core/tail split (x = 6).")
    print(f"inline constexpr double kErfinvTSplit = {float.hex(float(T_SPLIT))};")
    print(f"inline constexpr double kErfinvTFar = {float.hex(float(T_FAR))};")
    print()
    print("// Seed polynomials for x/t, index 0/1 in t on [t_lo,2]/[2,t_far],")
    print("// index 2 in u=1/t on [1/t_hi,1/t_far]. Padded to one degree so")
    print("// the kernel selects per-lane and runs a single Horner pass.")
    print(f"inline constexpr int kErfinvSeedNCoef = {n_seed};")
    print(f"inline constexpr double kErfinvSeedCoef[3][{n_seed}] = {{")
    for name, coefs in (("t in [t_lo, 2]", SEED0), ("t in [2, t_far]", SEED1),
                       ("u = 1/t in [1/t_hi, 1/t_far]", SEED2)):
        print(f"    {{  // {name}")
        for cf in coefs:
            print(f"        {float.hex(cf)},")
        print("    },")
    print("};")
    print()
    print("}  // namespace corvus::detail")
    print()
    print("#endif  // CORVUS_ERFINV_DATA_H_")
    return 0


if __name__ == "__main__":
    sys.exit(main())
