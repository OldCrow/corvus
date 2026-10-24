"""Residual annotation for the reference files (#82, T1).

Appends `r=<residual>` to every row of a public-function reference file:
residual = (exact - stored CR value) / ulp(stored value's binade), so a
correct row has |r| <= 1/2. A file with two output columns (gamma, beta: P
and Q) carries `r=<rP>,<rQ>`. The gates read it through LoadRef and
RealErr (tests/ulp_utils.h) and gate REAL error -- |result - exact| in ulps
-- instead of only the integer ULP distance, which cannot see a value drift
from 1.1 to 1.45 under a bound of 1.

The exact values come from the hunt families' oracles (tools/hunt/families).
A row whose |r| exceeds 1/2, or whose oracle fails, is recomputed at higher
precision and with the independent oracle2 where one exists; if the stored
value is still not the correctly rounded one, the row is a FINDING (a
generator defect), listed on stderr, and the file is not rewritten.

Where the hunt oracles do not reach, other routes take over:
- beta inverses where the hunt root finder fails: forward bisection inside
  the stored value's cell (the gamma Temme route and the beta normal-median
  form, once workarounds here, live in the oracles since #82 T3);
- special rows (NaN, infinite or zero-shape inputs): r = 0, since the gates
  match their exact special values bit for bit.
betainv B rows (beyond resolution) certify a neighbour of the root, so
|r| < 3/2 passes there. KNOWN_BAD lists stored values known not to be
correctly rounded, with the issue that tracks them.

Rows stay byte-identical apart from the token; re-running replaces an
existing token. Comment and blank lines are kept as they are. A row whose
stored value changes is recomputed (the checkpoint records the value).

The internal double-double reference files (exp_dd, log_dd, expm1_dd,
dd_special, log1pmx_ulo) are out of scope: their gates already measure a
continuous dd relative error. betainv_backward holds certified intervals,
not CR values.

  python tools/hunt/annotate.py [--workers 6] [--check] FILE_KEY ...
  python tools/hunt/annotate.py --list

Long families (gamma, beta, the inverses) checkpoint every result under
corvus-scratch/annotate/, so an interrupted run resumes.
"""
from __future__ import annotations

import argparse
import importlib
import json
import math
import multiprocessing as mp_
import sys
import time
from pathlib import Path

import mpmath as mp
from mpmath import mpf

sys.path.insert(0, str(Path(__file__).resolve().parent))
from core import (LowMemory, MAX_WORKERS, MIN_FREE_MB, REPO, TASKS_PER_WORKER,  # noqa: E402
                  free_mb, quantum)

DATA = REPO / "tests" / "data"
CKPT_DIR = REPO.parent / "corvus-scratch" / "annotate"

# file key -> (family modules, [(fn, argument columns, output column), ...]).
# A file listed with several outputs gets one residual per output, in order.
_U = lambda fam, fn: ([fam], [(fn, (0,), 1)])  # noqa: E731
FILES = {
    "erf": _U("erf", "erf"),
    "erfc": _U("erfc", "erfc"),
    "erfinv": _U("erfinv", "erfinv"),
    "erfcinv": _U("erfinv", "erfcinv"),
    "lgamma": _U("lgamma", "lgamma"),
    "digamma": _U("digamma", "digamma"),
    "trigamma": _U("trigamma", "trigamma"),
    "i0": _U("bessel", "i0"),
    "i1": _U("bessel", "i1"),
    "i0e": _U("bessel", "i0e"),
    "i1e": _U("bessel", "i1e"),
    "exp": _U("elementary", "exp"),
    "log": _U("elementary", "log"),
    "log1p": _U("elementary", "log1p"),
    "cos": _U("elementary", "cos"),
    "sin": _U("elementary", "sin"),
    "lbeta": (["lbeta"], [("lbeta", (0, 1), 2)]),
    "gamma_p": (["gamma"], [("gamma_p", (0, 1), 2), ("gamma_q", (0, 1), 3)]),
    "gamma_q": (["gamma"], [("gamma_p", (0, 1), 2), ("gamma_q", (0, 1), 3)]),
    "beta_p": (["beta"], [("beta_p", (0, 1, 2), 3), ("beta_q", (0, 1, 2), 4)]),
    "beta_q": (["beta"], [("beta_p", (0, 1, 2), 3), ("beta_q", (0, 1, 2), 4)]),
    "gammainv_p": (["gammainv"], [("gamma_p_inv", (0, 1), 2)]),
    "gammainv_q": (["gammainv"], [("gamma_q_inv", (0, 1), 2)]),
    "betainv_p": (["betainv"], [("beta_p_inv", (0, 1, 2), 3)]),
    "betainv_q": (["betainv"], [("beta_q_inv", (0, 1, 2), 3)]),
}
# Rows whose stored value is known not to be correctly rounded, with the
# issue that tracks the generator defect. They keep their true residual,
# so the gates measure real error against the exact value regardless.
KNOWN_BAD = {
    ("betainv_q", "0x1.0000000000000p+0 0x1.0000000000000p-1 "
                  "0x1.6a09e667f3bcdp-1 0x1.ffffffffffffep-2"):
        "a=1, b=1/2: root is 1/2 - 1.23 ulp, stored 1/2 - 2 ulp; generator, #87",
}
HALF_SLACK = mpf("1e-9")  # |r| <= 1/2 + this is a correct row
with mp.workprec(64):  # 2^1024 - 2^970 needs 54 bits
    OVERFLOW = mpf(2 ** 54 - 1) * mpf(2) ** 970  # rounds to infinity


def path_of(key: str) -> Path:
    return DATA / f"{key}_reference.txt"


# ---- worker -----------------------------------------------------------------
_ORACLES: dict = {}


def _init(family_names):
    global _ORACLES
    _ORACLES = {}
    for name in family_names:
        for s in importlib.import_module(f"families.{name}").strata():
            o1, o2 = _ORACLES.get(s.fn, (None, None))
            _ORACLES[s.fn] = (o1 or s.oracle, o2 or s.oracle2)


def _residual(v, cr: float):
    """(v - cr) / ulp(cr), or a special-value verdict."""
    if math.isnan(cr):
        return mpf(0) if mp.isnan(v) else None
    if math.isinf(cr):
        ok = mp.isinf(v) or abs(v) >= OVERFLOW
        return mpf(0) if ok and mp.sign(v) == math.copysign(1, cr) else None
    if not mp.isfinite(v):
        return None
    return (v - mpf(cr)) / quantum(cr)


def _beta_inv_fallback(fn, args, cr):
    """Residual of a beta_p_inv / beta_q_inv row where the hunt root finder
    fails or its two forms disagree (2026-10-22: a, b at 1e-94 / 1e-189
    with the root within e^-1e131 of 1): bisection on the forward
    betainc_pq inside cr +- 3/4 ulp, which needs the root in that cell (a
    bracket, not a guess)."""
    from families._beta_oracles import betainc_pq
    a, b, s = (mpf(v) for v in args)
    side_p = fn == "beta_p_inv"
    q = quantum(cr)
    with mp.workdps(60):
        lo = max(mpf(0), mpf(cr) - q * mpf(0.75))
        hi = min(mpf(1), mpf(cr) + q * mpf(0.75))

        def g(x):  # increasing in x, root at 0
            P, Q = betainc_pq(a, b, x, 50)
            return P - s if side_p else s - Q

        glo, ghi = g(lo), g(hi)
        if not (glo <= 0 <= ghi):
            raise RuntimeError("root not inside cr +- 3/4 ulp")
        for _ in range(48):
            mid = (lo + hi) / 2
            if g(mid) <= 0:
                lo = mid
            else:
                hi = mid
        return float(((lo + hi) / 2 - mpf(cr)) / q), "forward bisection"


def _work(task):
    """task = (fn, args, cr) -> (fn, args, r or None, certified, note)."""
    fn, args, cr = task
    zero_shape = fn in ("beta_p", "beta_q") and 0.0 in args[:2]
    if math.isnan(cr) or zero_shape or not all(math.isfinite(v) for v in args):
        # Special rows (NaN, infinite, zero-shape or out-of-domain inputs): the stored
        # value is the contract's exact special value; gates match it bit
        # for bit and take no real error from it.
        return fn, args, 0.0, True, "special: exact by contract"
    res = _work_oracles(task)
    if fn in ("beta_p_inv", "beta_q_inv") and (res[2] is None or not res[3]):
        try:
            r, how = _beta_inv_fallback(fn, args, cr)
            return fn, args, r, True, f"fallback: {how}"
        except Exception as exc:
            return res[:4] + (f"{res[4]}; fallback: {exc}"[:240],)
    return res


def _work_oracles(task):
    fn, args, cr = task
    o1, o2 = _ORACLES[fn]
    note = ""
    try:
        with mp.workdps(40):
            v = +o1(args, 40)
        with mp.workdps(90):
            r = _residual(v, cr)
        if r is not None and abs(r) <= mpf(0.5) + HALF_SLACK:
            return fn, args, float(r), True, ""
    except Exception as exc:  # oracle refused: certify the hard way below
        note = f"oracle: {type(exc).__name__}: {exc}"[:160]
    # Out of range or failed: higher precision, plus the independent form.
    try:
        with mp.workdps(80):
            v2 = +o1(args, 80)
        v3 = None
        if o2 is not None:
            with mp.workdps(60):
                v3 = +o2(args, 60)
        with mp.workdps(90):
            r2 = _residual(v2, cr)
            agree = v3 is None or (mp.isfinite(v2) and abs(v3 - v2) <= abs(v2) * mpf(2) ** -70
                                   + mpf(2) ** -1100) or (v3 == v2)
        if r2 is None:
            return fn, args, None, bool(agree), f"special-value mismatch {mp.nstr(v2, 20)}"
        return fn, args, float(r2), bool(agree), note
    except Exception as exc:
        return fn, args, None, False, f"{note}; cert: {type(exc).__name__}: {exc}"[:240]


# ---- driver -----------------------------------------------------------------
def _ckey(fn, args):
    return fn + "|" + ",".join(a.hex() for a in args)


def _pool(family_names, tasks, ckpt: Path, workers, batch=256, batch_seconds=1200):
    done = {}
    if ckpt.exists():
        for line in ckpt.read_text().splitlines():
            d = json.loads(line)
            done[d["key"]] = (d["r"], d["cert"], d["note"], d.get("cr"))
    # A residual is relative to the stored value: a regenerated row whose
    # value changed is computed again.
    todo = [t for t in tasks if _ckey(t[0], t[1]) not in done
            or done[_ckey(t[0], t[1])][3] not in (None, t[2].hex())]
    print(f"  {len(tasks) - len(todo)} rows from checkpoint, {len(todo)} to compute",
          file=sys.stderr, flush=True)
    if not todo:
        return done
    ckpt.parent.mkdir(parents=True, exist_ok=True)
    fh = ckpt.open("a")
    pool = mp_.Pool(workers, initializer=_init, initargs=(family_names,),
                    maxtasksperchild=TASKS_PER_WORKER)

    cr_of = {_ckey(t[0], t[1]): t[2].hex() for t in todo}

    def keep(fn, args, r, cert, note):
        k = _ckey(fn, args)
        done[k] = (r, cert, note, cr_of[k])
        fh.write(json.dumps({"key": k, "r": r, "cert": cert, "note": note,
                             "cr": cr_of[k]}) + "\n")

    t0 = time.time()
    try:
        for i in range(0, len(todo), batch):
            pend = [(t, pool.apply_async(_work, (t,))) for t in todo[i:i + batch]]
            deadline = time.time() + batch_seconds
            while pend:
                still = []
                for t, pr in pend:
                    if pr.ready():
                        keep(*pr.get())
                    else:
                        still.append((t, pr))
                pend = still
                if not pend:
                    break
                fm = free_mb()
                if fm is not None and fm < MIN_FREE_MB:
                    pool.terminate()
                    fh.flush()
                    raise LowMemory(f"free memory {fm} MB < {MIN_FREE_MB} MB; pool stopped")
                if time.time() > deadline:
                    for t, _ in pend:
                        keep(t[0], t[1], None, False, "wall-clock timeout")
                    pool.terminate()
                    pool = mp_.Pool(workers, initializer=_init, initargs=(family_names,),
                                    maxtasksperchild=TASKS_PER_WORKER)
                    break
                time.sleep(0.1)
            fh.flush()
            n = min(i + batch, len(todo))
            if (i // batch) % 20 and n < len(todo):
                continue
            el = time.time() - t0
            print(f"  {n}/{len(todo)}  {el:.0f}s  eta {el / n * (len(todo) - n):.0f}s  "
                  f"free {free_mb()} MB", file=sys.stderr, flush=True)
    finally:
        pool.terminate()
        fh.close()
    return done


def _split_lines(raw: bytes):
    """Lines with their terminators, so CRLF/LF is preserved exactly."""
    return raw.decode("ascii").splitlines(keepends=True)


def _data_tokens(line: str):
    s = line.strip()
    if not s or s.startswith("#"):
        return None
    return [t for t in s.split() if not t.startswith("r=")]


def annotate(key: str, workers: int, check_only: bool) -> int:
    families, outs = FILES[key]
    path = path_of(key)
    lines = _split_lines(path.read_bytes())
    tasks, seen = [], set()
    for line in lines:
        tok = _data_tokens(line)
        if tok is None:
            continue
        for fn, acols, ocol in outs:
            args = tuple(float.fromhex(tok[c]) for c in acols)
            k = _ckey(fn, args)
            if k not in seen:
                seen.add(k)
                tasks.append((fn, args, float.fromhex(tok[ocol])))
    print(f"{key}: {len(tasks)} oracle evaluations", file=sys.stderr, flush=True)
    done = _pool(families, tasks, CKPT_DIR / f"{families[0]}.ckpt.jsonl", workers)

    findings, out_lines, worst = [], [], 0.0
    for ln, line in enumerate(lines, 1):
        tok = _data_tokens(line)
        if tok is None:
            out_lines.append(line)
            continue
        rs = []
        for fn, acols, ocol in outs:
            args = tuple(float.fromhex(tok[c]) for c in acols)
            r, cert, note, _ = done[_ckey(fn, args)]
            # betainv B rows (beyond resolution) certify a neighbour of the
            # root, not its rounding: |r| < 3/2 is their contract.
            limit = 1.5 if key.startswith("betainv") and tok[4] == "B" else 0.5
            known = KNOWN_BAD.get((key, " ".join(tok[:4])))
            if known:
                print(f"  known: {path.name}:{ln}: r={r} ({known})", file=sys.stderr)
            elif r is None or not cert or abs(r) > limit + 1e-9:
                findings.append(f"{path.name}:{ln}: {fn} r={r} certified={cert} {note}")
            if r is not None:
                worst = max(worst, abs(r))
            rs.append("nan" if r is None else f"{r:.6f}")
        body = line.rstrip("\r\n")
        eol = line[len(body):]
        body = " ".join(t for t in body.split() if not t.startswith("r="))
        out_lines.append(f"{body} r={','.join(rs)}{eol}")
    print(f"{key}: max |r| = {worst:.6f}, {len(findings)} findings", file=sys.stderr)
    for f in findings[:40]:
        print("  FINDING " + f, file=sys.stderr)
    if findings:
        return 1
    if not check_only:
        path.write_bytes("".join(out_lines).encode("ascii"))
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("keys", nargs="*")
    ap.add_argument("--workers", type=int, default=MAX_WORKERS)
    ap.add_argument("--check", action="store_true", help="compute and report; do not rewrite")
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()
    if a.list or not a.keys:
        print(" ".join(FILES))
        return 0
    workers = min(a.workers, MAX_WORKERS)
    rc = 0
    for key in a.keys:
        rc |= annotate(key, workers, a.check)
    return rc


if __name__ == "__main__":
    sys.exit(main())
