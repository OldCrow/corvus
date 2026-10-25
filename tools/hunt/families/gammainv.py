"""gamma_p_inv / gamma_q_inv hunt strata (#82).

Bound: 1 ULP over the whole domain (docs/ACCURACY.md, gamma_p_inv and
gamma_q_inv). The deep-small closed form (x*(1+a) < 2^-60) and the
subnormal/zero-result rows are documented as correctly rounded (0 ULP), but
docs-correction issue #76 is open, so they are judged at 1 ULP here with
tie_exempt = 0 (a 1-ULP miss there is reported only as at-risk; read the ulp
column against the documented 0).

Oracle: the exact root, _gamma_oracles.inv_root -- solved on the small side
(s > 1/2 flipped exactly), bracketed in double bit space from a scipy float
guess, bisected to adjacent doubles, then one secant step (closed form below
2^-1000). oracle2 repeats the root search with the independent forward forms.
Not covered: 1e7 < a < 2^130 (series/CF too slow, the beyond-resolution
approximation not yet valid) -- this leaves kGammaInvBracketAMax = 1e16 and
the ~3e34 beyond-resolution onset without strata. For a >= 2^130 the root is
a + sqrt(a)*z (transition inside half an ulp of a).
"""
import math

from core import Stratum, full_mantissa, log_uniform, neighbours, subnormal

from . import _gamma_oracles as G

FNS = ("gamma_p_inv", "gamma_q_inv")
BAND_SCALE = 1.0  # sampled strata; bit-stepped strata are fixed-size

DEEP_CUT = 2.0 ** -60
ETA_MAX = float.fromhex("0x1.27d1135318e42p+3")
CR_NOTE = "documented CR (0 ULP); judged at 1 ULP pending #76"

SEAMS = [
    {"fn": "gamma_p/q_inv", "name": "aT", "coord": "a", "value": 20.0, "source": "kGammaAT: forward E split, seed family a < a_T"},
    {"fn": "gamma_p/q_inv", "name": "s1_amin", "coord": "a", "value": float.fromhex("0x1.3333333333333p-2"), "source": "kGammaInvS1AMin (S1 seed domain)"},
    {"fn": "gamma_p/q_inv", "name": "r4_amax", "coord": "a", "value": 1.5, "source": "kGammaR4AMax: lgamma(1+a) zone form vs lga+lna"},
    {"fn": "gamma_p/q_inv", "name": "zone_half", "coord": "a", "value": 0.5, "source": "kLgammaZoneLo: lgamma(1+a) centre switch"},
    {"fn": "gamma_p/q_inv", "name": "tiny_a", "coord": "a", "value": 2.0 ** -900, "source": "kGammaInvTinyA (deep-small scaled numerator, #83)"},
    {"fn": "gamma_p/q_inv", "name": "lnct_log1p", "coord": "t = min(s, 1-s)", "value": 2.0 ** -7, "source": "ln(1 - t) as Log1pDd below 2^-7 (#83)"},
    {"fn": "gamma_p/q_inv", "name": "cf_amax", "coord": "a", "value": 2.0 ** 996, "source": "kGammaInvCfAMax (CF scrub; beyond-resolution band)"},
    {"fn": "gamma_p/q_inv", "name": "bracket_amax", "coord": "a", "value": 1e16, "source": "kGammaInvBracketAMax: NOT COVERED (oracle gap 1e7 < a < 2^130)"},
    {"fn": "gamma_p/q_inv", "name": "beyond_resolution", "coord": "a", "value": 3e34, "source": "derived onset: NOT COVERED (oracle gap)"},
    {"fn": "gamma_p/q_inv", "name": "flip", "coord": "s", "value": 0.5, "source": "GammaInvVec m_flip (s > 1/2 -> 1 - s)"},
    {"fn": "gamma_p/q_inv", "name": "deep_small", "coord": "x*(1+a)", "value": DEEP_CUT, "source": "kGammaInvDeepSmallCut (closed form)"},
    {"fn": "gamma_p/q_inv", "name": "s3_gate", "coord": "L/|a-1|", "value": 3.0, "source": "kGammaInvS3StabilityMargin, L = -ln(q Gamma(a))"},
    {"fn": "gamma_p/q_inv", "name": "eta_max", "coord": "|eta0|", "value": ETA_MAX, "source": "kGammaInvEtaMax (S1 guard), eta0 = erfcinv(2t) sqrt(2/a)"},
    {"fn": "gamma_p/q_inv", "name": "eta_series_cut", "coord": "|eta|", "value": 0.5, "source": "kGammaInvEtaSeriesCut (lambda-of-eta series vs Newton)"},
    {"fn": "gamma_p/q_inv", "name": "out_min_normal", "coord": "x", "value": 2.0 ** -1022, "source": "derived: subnormal results"},
    {"fn": "gamma_p/q_inv", "name": "out_zero", "coord": "x", "value": 2.0 ** -1075, "source": "derived: zero results (kGammaInvQuotMax clamp lies inside this band)"},
    # not mappable / unreachable
    {"fn": "gamma_p/q_inv", "name": "bracket_entry", "coord": "|Newton step|", "value": 2.0 ** -7, "source": "kGammaInvBracketEntry/Tol/RStall: iteration-state thresholds (unmapped)"},
    {"fn": "gamma_p/q_inv", "name": "bracket_ylim", "coord": "log x", "value": 709.0, "source": "kGammaInvBracketYMin/YMax: unreachable outside deep-small with a < 1e16 (unmapped)"},
    {"fn": "gamma_p/q_inv", "name": "shallow", "coord": "-", "value": 2.0 ** -10, "source": "kGammaInvShallowThreshold: emitted by the generator, unused by the kernel"},
]


def _oracle(fn):
    return lambda args, dps: G.gammainv_value(fn, args[0], args[1], dps)


def _oracle2(fn):
    return lambda args, dps: G.gammainv_value(fn, args[0], args[1], dps, two=True)


def _bound(args, v):
    return 1


def _ok(r):
    return (r is not None and r[1] is not None and r[0] > 0 and math.isfinite(r[0])
            and 0 < r[1] < 1)


def _take(rng, rows, n):
    rows = [r for r in rows if _ok(r)]
    return rows if n >= len(rows) else rng.sample(rows, n)


def _many(rng, n, f, tries=30):
    out = []
    for _ in range(n * tries):
        if len(out) >= n:
            break
        r = f(rng)
        if _ok(r):
            out.append(r)
    return out


def _sub(r):
    """Subnormal, half uniform in mantissa bits, half uniform in exponent."""
    return subnormal(r) if r.random() < 0.5 else subnormal(r, 52)


def s_mixed(r):
    u = r.random()
    if u < 0.35:
        return log_uniform(r, 1e-300, 0.5)
    if u < 0.45:
        return _sub(r)
    if u < 0.75:
        return full_mantissa(r, 1e-3, 0.999)
    return 1.0 - log_uniform(r, 1e-16, 0.5)


def _target(fn, a, x):
    """The double s whose exact root lies at x (P side for gamma_p_inv)."""
    try:
        P, Q = G.pq(a, x, 25)
    except G.OracleError:
        return None
    return float(P if fn == "gamma_p_inv" else Q)


S_SET = (1e-300, 1e-8, 0.1, 0.5, 0.9, 1.0 - 1e-12)


def strata():
    out = []

    def add(name, gen_for, n, seam=None, note="", per_fn=False):
        if not (name.endswith(":bits") or name.startswith("witness")):
            n = max(1, int(n * BAND_SCALE))
        for fn in FNS:
            gen = gen_for(fn) if per_fn else gen_for
            out.append(Stratum(fn, name, gen, n, _oracle(fn), _bound,
                               oracle2=_oracle2(fn), seam=seam, note=note))

    # 1. a seams.
    tiny = 2.0 ** -900
    a_seams = [
        ("aT", 20.0, S_SET),
        ("s1_amin", float.fromhex("0x1.3333333333333p-2"), S_SET),
        ("r4_amax", 1.5, S_SET),
        ("zone_half", 0.5, S_SET),
        ("tiny_a", tiny, (6e-323, 1e-300, 1e-272, 1e-270, 1e-269, 0.3)),
        ("cf_amax", 2.0 ** 996, (1e-300, 0.3, 0.5, 0.9)),
    ]
    for name, t, sset in a_seams:
        add(f"seam:{name}:bits",
            lambda r, n, t=t, sset=sset: _take(r, [(a, s) for a in neighbours(t, 64) for s in sset], n),
            129 * len(sset), seam=name, note=CR_NOTE if name == "tiny_a" else "")
        if name == "cf_amax":
            continue
        for side, tag in ((-1, "lo"), (1, "hi")):
            lo, hi = (t * 0.95, t) if side < 0 else (t, t * 1.05)
            if name == "tiny_a":
                samp = lambda q, lo=lo, hi=hi: (lambda a: (a, log_uniform(q, 5e-324, a * 800)))(full_mantissa(q, lo, hi))
            else:
                samp = lambda q, lo=lo, hi=hi: (full_mantissa(q, lo, hi), s_mixed(q))
            add(f"seam:{name}:band_{tag}", lambda r, n, f=samp: _many(r, n, f), 800, seam=name)

    # 2. The s = 1/2 flip.
    a_flip = (0.05, 1.0, 7.5, 20.0, 1e3, 1e6)
    add("seam:flip:bits",
        lambda r, n: _take(r, [(a, s) for a in a_flip for s in neighbours(0.5, 64)], n),
        129 * len(a_flip), seam="flip")
    for lo, hi, tag in ((0.45, 0.5, "lo"), (0.5, 0.55, "hi")):
        add(f"seam:flip:band_{tag}",
            lambda r, n, lo=lo, hi=hi: _many(r, n, lambda q: (log_uniform(q, 1e-3, 1e6), full_mantissa(q, lo, hi))),
            800, seam="flip")

    # 3. Deep-small cut x*(1+a) = 2^-60 (targets placed so the root sits there).
    a_ds = (1e-300, 1e-10, 1e-3, 0.3, 1.0, 4.0)

    def ds_bits(fn):
        def gen(r, n):
            rows = []
            for a in a_ds:
                s0 = _target(fn, a, DEEP_CUT / (1 + a))
                if s0 and 0 < s0 < 1:
                    rows += [(a, s) for s in neighbours(s0, 64)]
            return _take(r, rows, n)
        return gen
    add("seam:deep_small:bits", ds_bits, 129 * len(a_ds), seam="deep_small", note=CR_NOTE, per_fn=True)

    def ds_band(fn):
        return lambda r, n: _many(r, n, lambda q: (lambda a: (a, _target(fn, a, 2.0 ** (-60 + q.uniform(-4, 4)) / (1 + a))))(
            log_uniform(q, 1e-300, 8.0)))
    add("seam:deep_small:band", ds_band, 1200, seam="deep_small", note=CR_NOTE, per_fn=True)

    # 4. Subnormal / zero results: root x in [2^-1080, 2^-1015].
    def sub_band(fn):
        return lambda r, n: _many(r, n, lambda q: (lambda a: (a, _target(fn, a, 2.0 ** -q.uniform(1015, 1080))))(
            log_uniform(q, 1e-300, 1.2) if q.random() < 0.7 else subnormal(q)))
    add("subnormal_x:band", sub_band, 1500, seam="out_min_normal", note=CR_NOTE, per_fn=True)
    a_z = (0.3, 0.5, 1.0)

    def zero_bits(fn):
        def gen(r, n):
            rows = []
            for a in (a_z if fn == "gamma_p_inv" else (1e-300, 1e-100, 1e-20)):
                for xs in (2.0 ** -1075, 2.0 ** -1022):
                    s0 = _target(fn, a, xs)
                    if s0 and 0 < s0 < 1:
                        rows += [(a, s) for s in neighbours(s0, 64)]
            return _take(r, rows, n)
        return gen
    add("seam:out_zero_min_normal:bits", zero_bits, 129 * 2 * len(a_z), seam="out_zero", note=CR_NOTE, per_fn=True)

    # 5. Subnormal targets (#72: a >~ 4500, p <= 2^-1054 is known bad).
    add("subnormal_s:a_20_1e6",
        lambda r, n: [(log_uniform(r, 20.0, 1e6), _sub(r)) for _ in range(n)], 3000,
        note="#72 band (a >~ 4500, p <= 2^-1054)")
    add("subnormal_s:a_lt_20",
        lambda r, n: [(log_uniform(r, 1e-3, 20.0), _sub(r)) for _ in range(n)], 1000)

    # 6. Seed-selection seams (S3 gate, S1 eta0 guard, lambda-of-eta series cut).
    def s3_q(a, f=1.0):
        return math.exp(-3.0 * abs(a - 1.0) * f - math.lgamma(a))

    def s3_gen(fn):
        def conv(a, q):
            return (a, q) if fn == "gamma_q_inv" else (a, 1.0 - q)
        return lambda r, n: _many(r, n, lambda g: (lambda a: conv(a, s3_q(a, g.uniform(0.9, 1.1))))(
            full_mantissa(g, 0.01, 0.7) if g.random() < 0.3 else full_mantissa(g, 1.3, 20.0)))
    add("seam:s3_gate:band", s3_gen, 1500, seam="s3_gate", per_fn=True,
        note="gamma_p_inv rows carry p = 1 - q (exact flip only where q > 2^-53)")
    a_s3 = (0.05, 0.4, 2.5, 7.0, 15.0, 19.5)
    add("seam:s3_gate:bits",
        lambda r, n: _take(r, [(a, s) for a in a_s3 for s in neighbours(s3_q(a), 64)], n),
        129 * len(a_s3), seam="s3_gate")
    for name, c, arange in (("eta_max", ETA_MAX, (0.3, 17.0)), ("eta_series_cut", 0.5, (0.3, 3e3))):
        add(f"seam:{name}:band",
            lambda r, n, c=c, arange=arange: _many(r, n, lambda q: (lambda a: (a, 0.5 * math.erfc(c * math.sqrt(a / 2) * q.uniform(0.9, 1.1))))(
                log_uniform(q, *arange))),
            1000, seam=name)
        a_b = (0.3, 1.0, 4.0, 10.0)
        add(f"seam:{name}:bits",
            lambda r, n, c=c, a_b=a_b: _take(r, [(a, s) for a in a_b
                                                 for s in neighbours(0.5 * math.erfc(c * math.sqrt(a / 2)), 64)], n),
            129 * len(a_b), seam=name)

    # 7. Crossed extremes.
    add("extreme:beyond_resolution",
        lambda r, n: [(math.ldexp(1.0 + r.random(), r.randint(130, 1022)), s_mixed(r)) for _ in range(n)], 600)
    add("extreme:big_a_ridge",
        lambda r, n: [(log_uniform(r, 1e6, 1e7), log_uniform(r, 1e-10, 0.5) if r.random() < 0.5
                       else 1 - log_uniform(r, 1e-10, 0.5)) for _ in range(n)], 300,
        note="a in [1e6, 1e7] (costliest oracle rows)")
    add("extreme:subnormal_a",
        lambda r, n: _many(r, n, lambda q: (lambda a: (a, log_uniform(q, 5e-324, a * 800) if q.random() < 0.7 else s_mixed(q)))(
            subnormal(q))), 800, note=CR_NOTE)
    add("extreme:tiny_a_q_tail",
        lambda r, n: _many(r, n, lambda q: (lambda a: (a, log_uniform(q, 5e-324, a * 800)))(log_uniform(q, 1e-300, 1e-5))),
        1200)

    # 8. Generic coverage.
    add("generic:small_a",
        lambda r, n: [(log_uniform(r, 1e-5, 20.0), s_mixed(r)) for _ in range(n)], 3000)
    add("generic:temme",
        lambda r, n: [(log_uniform(r, 20.0, 1e6), s_mixed(r)) for _ in range(n)], 3000)

    # 9. Known witnesses.
    wit = [(float.fromhex("0x1.1b1c756ce6137p+12"), float.fromhex("0x0.0000000000007p-1022")),
           (7.151707517755251, 0.6176099837213445), (0.1, 0.99), (20.0, 1e-300)]
    add("witness", lambda r, n: list(wit), len(wit), note="#72 (p=7 denorm, a~4508), #56, #48, #49")

    # 10. Subnormal shape AND subnormal target (#82 T7, the route census):
    # a = i * 2^-1074, s = j * 2^-1074. On the q side the forward's Q rounds
    # to 0 above the root, the residual is NaN, and the bracket re-proposes
    # one midpoint until its cap: roots 4-35% off (#83).
    D = 2.0 ** -1074
    aq = [(i * D, j * D) for i in (1, 2, 3, 5, 8, 13, 20, 34, 39, 46, 65, 77, 84, 120, 164, 240)
          for j in (1, 2, 3, 5, 8, 20)]
    add("witness:subnormal_aq", lambda r, n: list(aq), len(aq),
        note="gammainv/bracket/cap and bracket/reject_nan; q side wrong (#83)")
    # 11. Route witnesses: inputs that reach a route no other stratum reaches.
    rw = [(float.fromhex("0x1.35e480e2883cep-3"), float.fromhex("0x1.f8dbcaf539fbdp-1")),
          (float.fromhex("0x1.d0b296b5b1b25p-4"), float.fromhex("0x1.d8620a769c840p-7"))]
    add("witness:route_reject_nan", lambda r, n: list(rw), len(rw),
        note="gammainv/bracket/reject_nan at ordinary shapes (correctly rounded)")
    # 12. T2 completion (#82): seam buckets the audit found empty.
    cfa = 6.696928794914171e+299
    add("seam:cf_amax:band",
        lambda r, n: [(full_mantissa(r, cfa * 0.95, cfa * 1.05), s_mixed(r)) for _ in range(n)], 400,
        seam="cf_amax")
    add("seam:beyond_resolution:bits",
        lambda r, n: [(a, s) for a in neighbours(3e34, 16) for s in (0.3, 1e-5, 0.7)], 3 * 33,
        seam="beyond_resolution")
    add("seam:beyond_resolution:band",
        lambda r, n: [(full_mantissa(r, 3e34 * 0.95, 3e34 * 1.05), s_mixed(r)) for _ in range(n)], 400,
        seam="beyond_resolution")
    # 13. ln(1 - t) switches to Log1pDd at t = 2^-7 (#83). The deep-small rows
    # at tiny a amplify its error by |log x|, so both sides carry tiny-a rows.
    t7 = 2.0 ** -7
    a_t7 = (1e-20, 1e-3, 0.5, 7.5, 1e3)
    add("seam:lnct_log1p:bits",
        lambda r, n: _take(r, [(a, s) for a in a_t7 for s in neighbours(t7, 64)], n),
        129 * len(a_t7), seam="lnct_log1p")
    for lo, hi, tag in ((t7 * 0.95, t7, "lo"), (t7, t7 * 1.05, "hi")):
        add(f"seam:lnct_log1p:band_{tag}",
            lambda r, n, lo=lo, hi=hi: _many(r, n, lambda q: (log_uniform(q, 1e-300, 1e6),
                                                              full_mantissa(q, lo, hi))),
            400, seam="lnct_log1p")
    return out


# ---- seam-audit fields (#82 T2; tools/hunt/seam_audit.py docstring) ----------
def _q_of(r):
    """The Q-side target of a reference row (gammainv_q files store q)."""
    return r["s"] if r["file"].endswith("_q") else 1.0 - r["s"]


def _s3_ratio(r):
    a, q = r["a"], _q_of(r)
    if not (0 < q and a > 0 and a != 1.0):
        return None
    L = -(math.log(q) + math.lgamma(a))
    return L / abs(a - 1.0)


def _eta0(r):
    a, q = r["a"], _q_of(r)
    if not (0 < q < 1 and a > 0):
        return None
    import mpmath as mp
    with mp.workdps(30):
        return abs(float(mp.erfinv(1 - 2 * mp.mpf(q)) * mp.sqrt(2 / mp.mpf(a))))


_AUDIT = {
    "s3_gate": {"audit_f": _s3_ratio},
    "eta_max": {"audit_f": _eta0},
    "eta_series_cut": {"audit_f": _eta0},  # the S1 seed's eta, as the kernel cuts
    "out_zero": {"sides": "above"},
    # t = 1 - p on a p row near 1 moves in steps of ulp(p) = 2^-53, far
    # coarser than 64 ulp of 2^-7: the bits window is 64 of those steps.
    "lnct_log1p": {"audit_f": lambda r: min(r["s"], 1.0 - r["s"]), "window": 64 * 2.0 ** -53},
    "bracket_entry": {"route": "gammainv/bracket"},
    "shallow": {"internal": "emitted by the generator, not read by the kernel"},
    "bracket_ylim": {"internal": "the bracket's log-x limits; the bracket is barred for "
                                 "a >= 1e16, and a root near e^709 needs a ~ 1e308"},
}
for _s in SEAMS:
    _s.update(_AUDIT.get(_s["name"], {}))
