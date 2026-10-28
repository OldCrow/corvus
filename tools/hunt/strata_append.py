"""Strata appendix for the reference files (#82, T3).

Appends to a public-function reference file a block of rows drawn from the
hunt strata (tools/hunt/families), so the populations the hunt designed --
seam bit-neighbourhoods, full-mantissa bands on each side of a seam,
subnormal inputs and outputs, dense fills -- become permanent gate rows.
Two strata are added here:
  near_tie   rows whose exact value lies within 1/32 ulp of a rounding
             midpoint, filtered from the family's non-seam strata (a 0-ULP
             cell cannot fail without them);
  cross      the cartesian product of each argument's extremes (subnormal,
             tiny, kernel thresholds, huge) for the multi-argument families.
The existing rows are kept byte-identical. Each appended stratum starts with
a `# stratum <name>` line; the rows before the first one are stratum `base`.
The per-stratum counts go to tests/data/ref_counts.txt, which LoadRef
(tests/ulp_utils.h) checks exactly.

Every row is certified: the first oracle at 40 digits gives the correctly
rounded value; the second form (or the first at 80 digits) must round to the
same double, with residuals agreeing to 1e-4 ulp and inside the rounding cell
by that margin. A row that fails is dropped
and counted. The kernel is then run through tests/hunt_probe on every tier
in core.TIERS; a row whose integer ULP error exceeds its stratum's published
bound on any tier moves to `# stratum deferred:?`. LoadRef fails on `?`:
name the issue that owns the defect (`deferred:#84`). A --replace run keeps
every row of a `deferred:#NN` section there, rows moved in by hand included
(some fail only a gate's own rule, such as an exact-zero check). Deferred rows are
counted but not gated unless CORVUS_REF_DEFERRED=1, so the fix row has its
fail-first gate ready (docs/NUMERICAL-DOCTRINE.md, test doctrine).

  python tools/hunt/strata_append.py KEY ... [--workers 6] [--dry-run] [--replace]

--replace removes an existing appendix (every `# stratum` section) first.
--count FILE ... records the counts of files this tool does not append to
(the dd files, betainv_backward) or after a generator rewrote one.
Rows are computed through a checkpoint under corvus-scratch/strata/, so an
interrupted run resumes. Run the kernel builds in core.TIERS first.
"""
from __future__ import annotations

import argparse
import importlib
import itertools
import json
import math
import multiprocessing as mp_
import random
import sys
import time
from pathlib import Path

import mpmath as mp
from mpmath import mpf

sys.path.insert(0, str(Path(__file__).resolve().parent))
import annotate  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import gen_bessel_reference  # noqa: E402  the #89 half-x rule
from core import (LowMemory, MAX_WORKERS, MIN_FREE_MB, REPO, TASKS_PER_WORKER,  # noqa: E402
                  TIERS, evaluate, free_mb, plan_workers, quantum, round_double, ulp_dist)

DATA = REPO / "tests" / "data"
COUNTS = DATA / "ref_counts.txt"
CKPT_DIR = REPO.parent / "corvus-scratch" / "strata"
SEED = 20261023
AGREE = mpf("1e-4")   # two-form agreement that certifies a row, in ulps
TIE_RESOLVE = mpf("1e-12")  # |r| this close to 1/2: re-evaluate at TIE_DPS
TIE_DPS = 400
TIE = 1 / 32          # near_tie: distance to the rounding midpoint, in ulps
TIE_ROWS, TIE_POOL = 128, 1024
CAPS = {"bits": 48, "band": 80, "other": 120}  # rows per hunt stratum (Stratum.ref_cap overrides)

# file key -> the probe function whose strata feed it (FILES lists outputs).
PRIMARY = {k: k for k in annotate.FILES}
PRIMARY.update({"erfcinv": "erfcinv", "gammainv_p": "gamma_p_inv",
                "gammainv_q": "gamma_q_inv", "betainv_p": "beta_p_inv",
                "betainv_q": "beta_q_inv"})

# ---- cross: argument extremes -----------------------------------------------
_TINY = [5e-324, 2.0 ** -1022, 1e-300, 1e-20]
SHAPE = _TINY + [2.0 ** -6, 0.5, 1.0, 20.0, 1e7, 2.0 ** 47, 1e15, 1e100, 1e300]
X01 = _TINY + [2.0 ** -11, 0.5, 1.0 - 2.0 ** -11, 1.0 - 2.0 ** -53]
POS = _TINY + [0.5, 1.0, 20.0, 1e7, 2.0 ** 47, 1e100, 1.7976931348623157e308]
PROB = _TINY + [1e-3, 0.5, 1.0 - 2.0 ** -11, 1.0 - 2.0 ** -53]
CROSS = {
    "gamma_p": [SHAPE, POS], "gamma_q": [SHAPE, POS],
    "gamma_p_inv": [SHAPE, PROB], "gamma_q_inv": [SHAPE, PROB],
    "beta_p": [SHAPE, SHAPE, X01], "beta_q": [SHAPE, SHAPE, X01],
    # beta inverses stop at 1e15: past it, off the median, the root solver's
    # precision growth takes gigabytes per row (a = 1e300, p = 1e-20 stopped
    # the pool), and the normal-median closed form covers the median only.
    "beta_p_inv": [SHAPE[:-2], SHAPE[:-2], PROB], "beta_q_inv": [SHAPE[:-2], SHAPE[:-2], PROB],
    "lbeta": [SHAPE + [1.7976931348623157e308], SHAPE + [1.7976931348623157e308]],
}


def _kind(name: str) -> str:
    if name.endswith(":bits"):
        return "bits"
    if name.startswith("seam:"):
        return "band"
    return "other"


def _rng(key, name):
    return random.Random(f"{SEED}|{key}|{name}")


# ---- worker: certified exact values -----------------------------------------
def _exact(fn, args):
    """(cr, r, note) for one output, or (None, None, why)."""
    o1, o2 = annotate._ORACLES[fn]
    with mp.workdps(40):
        v1 = +o1(args, 40)
    with mp.workdps(90):
        cr = round_double(v1) if mp.isfinite(v1) else float(v1)
        r1 = annotate._residual(v1, cr)
    v2 = None
    if o2 is not None:
        try:
            with mp.workdps(60):
                v2 = +o2(args, 60)
            form = "two forms"
        except Exception:  # the second form has a narrower domain (bessel)
            v2 = None
    if v2 is None:
        with mp.workdps(80):
            v2 = +o1(args, 80)
        form = "one form, 40/80 digits"
    with mp.workdps(90):
        r2 = annotate._residual(v2, cr)
    if r1 is None or r2 is None:
        return None, None, "special-value disagreement"
    # Certified when the forms agree to AGREE ulp and the residual stays
    # inside the rounding cell by that margin: then the side of the midpoint
    # is known. (The two Temme fits for large a agree to ~3e-6 ulp, 1e-22
    # relative, which the old 1e-6 cut dropped.)
    tol = abs(r1 - r2)
    if tol > AGREE or abs(r2) + tol > mpf(0.5) + annotate.HALF_SLACK:
        return None, None, f"forms disagree: r {mp.nstr(r1, 8)} vs {mp.nstr(r2, 8)}"
    # A residual this close to 1/2 does not say which side of the midpoint
    # the value is on: an algebraic near-tie such as I_x(1, 1/2) at
    # x = 2^-1074 (x/2 + x^2/8, the excess 2^-1076 relative) reads exactly
    # 0.5 at 40-90 digits, and was stored rounded to even on the wrong side
    # (row 12h). Re-evaluate at TIE_DPS and round from there; a residual
    # still unresolved there is dropped, never guessed.
    if abs(abs(r2) - mpf(0.5)) <= TIE_RESOLVE:
        with mp.workdps(TIE_DPS):
            v3 = +o1(args, TIE_DPS)
            if not mp.isfinite(v3):
                return None, None, "near-tie: high-precision form not finite"
            cr = round_double(v3)
            r3 = annotate._residual(v3, cr)
            # Inside the 400-digit block: the excess can be ~1e-324, which
            # the ambient 15 digits would round to an exact tie.
            unresolved = r3 is None or abs(abs(r3) - mpf(0.5)) <= mpf(10) ** -(TIE_DPS - 40)
        if unresolved:
            return None, None, "near-tie unresolved at %d digits" % TIE_DPS
        return cr, float(r3), form + f", tie resolved at {TIE_DPS} digits"
    return cr, float(r2), form


def _half_x(fn, args):
    """(cr, r) by gen_bessel_reference's #89 rule for i1 / i1e at |x| <
    2^-1021, where x/2 is a midpoint for odd x / 2^-1074 and the excess
    (~2^-2148 relative for i1) is beyond any oracle precision; else None."""
    if fn not in ("i1", "i1e"):
        return None
    cr = gen_bessel_reference.half_x_rule(1, fn == "i1e", args[0])
    if cr is None:
        return None
    with mp.workdps(30):
        r = float((mpf(args[0]) / 2 - mpf(cr)) / quantum(cr)) if cr else 0.0
    return cr, r


def _apply_half_x(fns, args, out):
    """A row's [(cr hex, r)] with the #89 rule applied (cached rows too)."""
    if out is None:
        return out
    fixed = []
    for fn, o in zip(fns, out):
        hx = _half_x(fn, args)
        fixed.append((hx[0].hex(), hx[1]) if hx else o)
    return fixed


def _work(task):
    """task = (outputs, args) -> (args, [(cr, r)] or None, note, v1 of first)."""
    fns, args = task
    out = []
    try:
        for fn in fns:
            hx = _half_x(fn, args)
            if hx is not None:
                out.append((hx[0].hex(), hx[1]))
                note = "#89 half-x rule"
                continue
            cr, r, note = _exact(fn, args)
            if cr is None:
                return args, None, note
            out.append((cr.hex() if math.isfinite(cr) else repr(cr), r))
        return args, out, note
    except Exception as exc:  # oracle refused: the row is dropped, counted
        return args, None, f"oracle: {type(exc).__name__}: {exc}"[:200]


def _pool(families, fns, args_list, ckpt: Path, workers, batch=256, batch_seconds=1200):
    done = {}
    if ckpt.exists():
        for line in ckpt.read_text().splitlines():
            d = json.loads(line)
            done[d["k"]] = (d["out"], d["note"])
    k_of = lambda a: ",".join(x.hex() for x in a)  # noqa: E731
    todo = [a for a in dict.fromkeys(args_list) if k_of(a) not in done]
    print(f"  {len(args_list) - len(todo)} rows from checkpoint, {len(todo)} to compute",
          file=sys.stderr, flush=True)
    if not todo:
        return done, k_of
    ckpt.parent.mkdir(parents=True, exist_ok=True)
    fh = ckpt.open("a")
    workers = plan_workers(workers)
    pool = mp_.Pool(workers, initializer=annotate._init, initargs=(families,),
                    maxtasksperchild=TASKS_PER_WORKER)

    def keep(args, out, note):
        done[k_of(args)] = (out, note)
        fh.write(json.dumps({"k": k_of(args), "out": out, "note": note}) + "\n")

    t0 = time.time()
    try:
        for i in range(0, len(todo), batch):
            pend = [(a, pool.apply_async(_work, ((fns, a),))) for a in todo[i:i + batch]]
            deadline = time.time() + batch_seconds
            while pend:
                still = []
                for a, pr in pend:
                    if pr.ready():
                        keep(*pr.get())
                    else:
                        still.append((a, pr))
                pend = still
                if not pend:
                    break
                fm = free_mb()
                if fm is not None and fm < MIN_FREE_MB:
                    pool.terminate()
                    fh.flush()
                    raise LowMemory(f"free memory {fm} MB < {MIN_FREE_MB} MB; pool stopped")
                if time.time() > deadline:
                    for a, _ in pend:
                        keep(a, None, "wall-clock timeout")
                    pool.terminate()
                    pool = mp_.Pool(workers, initializer=annotate._init, initargs=(families,),
                                    maxtasksperchild=TASKS_PER_WORKER)
                    break
                time.sleep(0.1)
            fh.flush()
            n = min(i + batch, len(todo))
            el = time.time() - t0
            print(f"  {n}/{len(todo)}  {el:.0f}s  eta {el / n * (len(todo) - n):.0f}s  "
                  f"free {free_mb()} MB", file=sys.stderr, flush=True)
    finally:
        pool.terminate()
        fh.close()
    return done, k_of


# ---- file handling -----------------------------------------------------------
def _split(raw: str):
    """(base lines, eol) with any existing appendix removed."""
    lines = raw.splitlines(keepends=True)
    eol = "\r\n" if lines and lines[0].endswith("\r\n") else "\n"
    for i, line in enumerate(lines):
        if line.startswith("# stratum "):
            return lines[:i], eol, True
    return lines, eol, False


def _owners(raw: str, nargs):
    """args -> issue for the rows of existing `deferred:#NN` sections, so a
    --replace run keeps the triage (an unassigned row stays `?`)."""
    out, issue = {}, None
    for line in raw.splitlines():
        if line.startswith("# stratum "):
            name = line[10:].strip()
            issue = name[9:] if name.startswith("deferred:#") else None
            continue
        tok = annotate._data_tokens(line)
        if tok and issue:
            out[tuple(float.fromhex(t) for t in tok[:nargs])] = issue
    return out


def _existing_args(lines, nargs):
    seen = set()
    for line in lines:
        tok = annotate._data_tokens(line)
        if tok:
            seen.add(tuple(float.fromhex(t) for t in tok[:nargs]))
    return seen


def _valid(args) -> bool:
    # Special inputs (NaN, infinities) belong to the smoke tests' contracts.
    return all(math.isfinite(a) for a in args)


def _draw(key, fn, strata, nargs, seen):
    """[(stratum name, [args])] from the hunt strata, the cross product, and
    the near_tie candidate pool (returned separately, unfiltered)."""
    plan, pool_pts = [], []
    for s in strata:
        cap = min(s.n, s.ref_cap or CAPS[_kind(s.name)])
        pts = [tuple(float(a) for a in p) for p in s.gen(_rng(key, s.name), cap)]
        rng = _rng(key, s.name + "|take")
        if len(pts) > cap:
            pts = [pts[i] for i in sorted(rng.sample(range(len(pts)), cap))]
        plan.append((s.name, pts))
        if _kind(s.name) == "other":
            pool_pts += [tuple(float(a) for a in p)
                         for p in s.gen(_rng(key, s.name + "|tie"), s.n)]
    if fn in CROSS:
        pts = [tuple(p) for p in itertools.product(*CROSS[fn])]
        if fn in ("beta_p_inv", "beta_q_inv"):
            # Both shapes tiny: I_x is flat to O(a + b) across the interior,
            # and the root solver's two forms agreed on x = 1 for
            # beta_q_inv(1e-300, 1e-300, 1/2), whose root is 1/2 by symmetry.
            pts = [p for p in pts if max(p[0], p[1]) > 1e-20]
        plan.append(("cross", pts))
    rng = _rng(key, "near_tie")
    if len(pool_pts) > TIE_POOL:
        pool_pts = rng.sample(pool_pts, TIE_POOL)
    # Deduplicate against the base rows and across strata, first claim wins.
    out = []
    for name, pts in plan:
        keep = []
        for p in pts:
            if len(p) == nargs and _valid(p) and p not in seen:
                seen.add(p)
                keep.append(p)
        out.append((name, keep))
    tie = [p for p in dict.fromkeys(pool_pts) if len(p) == nargs and _valid(p) and p not in seen]
    return out, tie


def _tie_dist(cr_hex, r):
    """Distance of the exact value to the rounding midpoint, in ulps."""
    return abs(abs(r) - 0.5)


def _violations(fn_bound, fns, args_list, results, k_of):
    """Indices of rows whose kernel output breaks its published bound on any
    tier, with the worst (tier, output, ulps)."""
    lines = [f"{f} " + " ".join(a.hex() for a in args) for args in args_list for f in fns]
    outs = evaluate(lines, TIERS)
    bad = {}
    for i, args in enumerate(args_list):
        res, _ = results[k_of(args)]
        for j, f in enumerate(fns):
            cr = float.fromhex(res[j][0]) if res[j][0].startswith(("0x", "-0x")) else float(res[j][0])
            bound, tol, tie = fn_bound[f](args, cr, res[j][1])
            # A CR cell's documented near-tie exception (#76): a 1-ULP result
            # whose stored residual is within the stratum's tie_exempt of 1/2
            # is a tie, not a violation -- the hunt scorer's and the gates' rule.
            tie_ok = bound == 0 and tie and 0.5 - abs(res[j][1]) < tie
            for tier in TIERS:
                o = outs[tier][i * len(fns) + j]
                if tol is not None:
                    with mp.workdps(60):
                        exact = mpf(cr) + mpf(res[j][1]) * quantum(cr)
                        if abs(mpf(o) - exact) > tol:
                            bad[i] = (tier, f, "abs")
                    continue
                u = ulp_dist(o, cr)
                if tie_ok and u == 1:
                    continue
                if u > bound and (i not in bad or bad[i][2] == "abs" or u > bad[i][2]):
                    bad[i] = (tier, f, u)
    return bad


def _bounds(strata_by_fn):
    """fn -> (args, cr, r) -> (integer bound, absolute tolerance or None,
    tie_exempt), from the family's first stratum for that fn (its
    region-aware bound)."""
    out = {}
    for f, s in strata_by_fn.items():
        def b(args, cr, r, s=s):
            with mp.workdps(60):
                v = mpf(cr) + mpf(r) * quantum(cr) if math.isfinite(cr) else mpf(cr)
                tol = s.abs_tol(args, v) if s.abs_tol else None
                return s.bound(args, v), tol, getattr(s, "tie_exempt", 0.0) or 0.0
        out[f] = b
    return out


def append(key, workers, dry_run, replace):
    families, outs = annotate.FILES[key]
    fns = [f for f, _, _ in outs]
    nargs = len(outs[0][1])
    fn = PRIMARY[key]
    path = annotate.path_of(key)
    raw = path.read_bytes().decode("ascii")
    base, eol, had = _split(raw)
    if had and not replace:
        print(f"{key}: already has an appendix (use --replace)", file=sys.stderr)
        return 1
    annotate._init(families)  # parent needs the strata too
    all_strata = [s for fam in families
                  for s in importlib.import_module(f"families.{fam}").strata()]
    strata = [s for s in all_strata if s.fn == fn]
    first = {}
    for s in all_strata:
        first.setdefault(s.fn, s)
    seen = _existing_args(base, nargs)
    plan, tie_pool = _draw(key, fn, strata, nargs, seen)
    n_plan = sum(len(p) for _, p in plan)
    print(f"{key}: {len(strata)} hunt strata, {n_plan} rows planned, "
          f"{len(tie_pool)} near-tie candidates", file=sys.stderr, flush=True)
    if dry_run:
        for name, pts in plan:
            print(f"  {name:40s} {len(pts)}")
        return 0
    args_all = [p for _, pts in plan for p in pts] + tie_pool
    results, k_of = _pool(families, fns, args_all, CKPT_DIR / f"{key}.ckpt.jsonl", workers)
    # Checkpoints written before #89 hold ties-to-even values for these rows.
    for p in args_all:
        out, note = results[k_of(p)]
        results[k_of(p)] = (_apply_half_x(fns, p, out), note)

    sections, dropped = [], {}
    for name, pts in plan:
        ok = [p for p in pts if results[k_of(p)][0] is not None]
        for p in pts:
            if results[k_of(p)][0] is None:
                dropped[results[k_of(p)][1][:60]] = dropped.get(results[k_of(p)][1][:60], 0) + 1
        sections.append((name, ok))
    fb = _bounds(first)

    def is_tie(p):
        # A near-tie row on a relative cell: rows judged by an absolute band
        # (zeros of lgamma, digamma) have no rounding midpoint to be near.
        res = results[k_of(p)][0]
        if res is None:
            return False
        for f, (c, r) in zip(fns, res):
            cr = float.fromhex(c) if c.startswith(("0x", "-0x")) else float(c)
            if not math.isfinite(cr) or fb[f](p, cr, r)[1] is not None:
                return False
        return min(_tie_dist(*o) for o in res) <= TIE

    ties = [p for p in tie_pool if is_tie(p)]
    sections.append(("near_tie", ties[:TIE_ROWS]))

    flat = [(name, p) for name, pts in sections for p in pts]
    bad = _violations(fb, fns, [p for _, p in flat], results, k_of)
    # A row already deferred under an issue stays deferred: some fail only a
    # gate's own rule (an exact-zero check, a gate tighter than the published
    # bound), which the probe check above does not model.
    owners = _owners(raw, nargs)
    for i, (_, p) in enumerate(flat):
        if p in owners and i not in bad:
            bad[i] = ("gate", fns[0], "kept")
    deferred = [(flat[i][0], flat[i][1], bad[i]) for i in sorted(bad)]
    keep = {}
    for i, (name, p) in enumerate(flat):
        if i not in bad:
            keep.setdefault(name, []).append(p)

    marker = " N" if key.startswith("betainv") else ""

    def row(p):
        res = results[k_of(p)][0]
        cols = [a.hex() for a in p] + [c for c, _ in res]
        return " ".join(cols) + marker + " r=" + ",".join(f"{r:.6f}" for _, r in res) + eol

    body, counts = list(base), [("base", sum(1 for l in base if annotate._data_tokens(l)))]
    for name, _ in sections:
        rows = keep.get(name, [])
        if rows:
            body.append(f"# stratum {name}{eol}")
            body += [row(p) for p in rows]
            counts.append((name, len(rows)))
    for issue in sorted({owners.get(p, "?") for _, p, _ in deferred},
                        key=lambda s: (s == "?", s)):
        part = [d for d in deferred if owners.get(d[1], "?") == issue]
        body.append(f"# stratum deferred:{issue}{eol}")
        for name, p, (tier, f, u) in part:
            body.append(f"# {name}: {f} {u} ulp on {tier}{eol}")
            body.append(row(p))
        counts.append((f"deferred:{issue}", len(part)))
    path.write_bytes("".join(body).encode("ascii"))
    _write_counts(path.name, counts)
    n_new = sum(c for n, c in counts[1:])
    print(f"{key}: appended {n_new} rows in {len(counts) - 1} strata "
          f"({len(ties[:TIE_ROWS])} near-tie, {len(deferred)} deferred); "
          f"dropped {sum(dropped.values())}", file=sys.stderr)
    for why, n in sorted(dropped.items(), key=lambda kv: -kv[1])[:8]:
        print(f"  dropped {n:5d}  {why}", file=sys.stderr)
    for name, p, (tier, f, u) in deferred[:20]:
        print(f"  DEFERRED {name}: {f}({', '.join(a.hex() for a in p)}) {u} ulp on {tier}",
              file=sys.stderr)
    return 0


def _write_counts(fname, counts):
    """Replace fname's entries in ref_counts.txt, keeping the file sorted."""
    head = ["# Exact row counts per reference file and stratum (#82 T3).",
            "# file|stratum|rows -- LoadRef (tests/ulp_utils.h) fails on any mismatch.",
            "# Written by tools/hunt/strata_append.py; `base` is every row before",
            "# the first `# stratum` line."]
    entries = []
    if COUNTS.exists():
        entries = [l for l in COUNTS.read_text().splitlines()
                   if l and not l.startswith("#") and l.split("|")[0] != fname]
    entries += [f"{fname}|{n}|{c}" for n, c in counts]
    entries.sort(key=lambda l: l.split("|")[0])
    COUNTS.write_text("\n".join(head + entries) + "\n")


def count_file(path: Path):
    """[(stratum, rows)] as LoadRef counts them."""
    counts = [["base", 0]]
    for line in path.read_text().splitlines():
        if line.startswith("# stratum "):
            counts.append([line[10:].rstrip(), 0])
        elif line.strip() and not line.lstrip().startswith("#"):
            counts[-1][1] += 1
    if counts[0][1] == 0:
        counts.pop(0)
    return [tuple(c) for c in counts]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("keys", nargs="*")
    ap.add_argument("--workers", type=int, default=MAX_WORKERS)
    ap.add_argument("--dry-run", action="store_true", help="print the plan only")
    ap.add_argument("--replace", action="store_true", help="rebuild an existing appendix")
    ap.add_argument("--count", nargs="+", metavar="FILE", help="record row counts only")
    a = ap.parse_args()
    if a.count:
        for f in a.count:
            _write_counts(Path(f).name, count_file(Path(f)))
        return 0
    if not a.keys:
        print(" ".join(annotate.FILES))
        return 0
    rc = 0
    for key in a.keys:
        rc |= append(key, min(a.workers, MAX_WORKERS), a.dry_run, a.replace)
    return rc


if __name__ == "__main__":
    sys.exit(main())
