"""Hunt harness core (#82, T5).

A stratum is a named input population for one public function, an exact
oracle, and the integer ULP bound that applies to each row. The runner
evaluates every stratum on every built tier through tests/hunt_probe, then
scores each row by REAL error (exact - result, in ulps of the exact value)
as well as by the integer ULP distance to the correctly rounded value.

Why real error: an integer bound B only fails once real error passes B + 0.5,
which a finite population hits rarely (i0's M-3 rows: ~2 per 1e6). The real
error's maximum over a stratum reaches the danger zone with a few thousand
rows. A row is
  VIOLATION  integer ULP > bound            (certified by a second oracle
                                             evaluation before it counts)
  AT-RISK    real error > bound + 0.5 - margin
Oracles must not hang: every series/CF carries an iteration cap and raises
on non-convergence. The pool still bounds each batch by wall clock and
replaces the pool if a worker wedges.
"""
from __future__ import annotations

import json
import math
import multiprocessing as mp_
import random
import struct
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional, Sequence

import mpmath as mp
from mpmath import mpf

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
# Point-selection helpers are shared with the reference generators (#82 T3).
from refgen_common import (bits, from_bits, full_mantissa, log_uniform,  # noqa: E402,F401
                           neighbours, ordered, subnormal, cap_process_memory, worker_mb)

REPO = Path(__file__).resolve().parents[2]
TIERS = {  # tier label -> build directory (Windows clang-cl builds on Z)
    "AVX3_ZEN4": "build-clangcl",
    "SSE4": "build-guards-sse4",
    "SSE2": "build-guards-sse2",
}
MARGIN = 0.25  # AT-RISK margin below the integer bound's real-error limit
MAX_WORKERS = 6  # cap-oracle-sweep-workers rule


# ---- double helpers ---------------------------------------------------------
def ulp_dist(a: float, b: float) -> int:
    if math.isnan(a) or math.isnan(b):
        return 0 if (math.isnan(a) and math.isnan(b)) else 1 << 62
    return abs(ordered(a) - ordered(b))


def quantum(v) -> mpf:
    """ulp of the binade holding |v| (subnormal grid below 2^-1022)."""
    av = abs(mpf(v))
    if av == 0:
        return mpf(2) ** -1074
    if not mp.isfinite(av):
        return mpf(2) ** 971  # ulp of the top binade
    _, e = mp.frexp(av)  # av = m * 2^e, m in [1/2, 1): exact, any precision
    return mp.ldexp(1, max(int(e) - 1, -1022) - 52)


def round_double(v) -> float:
    """Single rounding of an exact mpf to binary64, subnormals included."""
    v = mpf(v)
    if v == 0 or not mp.isfinite(v):
        return float(v)
    if abs(v) >= mpf(2) ** 1024 * (1 - mpf(2) ** -54):
        return math.copysign(math.inf, float(mp.sign(v)))
    q = quantum(v)
    r = float(mp.nint(v / q) * q)
    # An mpf zero is unsigned: a negative value that rounds to zero must come
    # back as -0 (log1pmx's underflow, #45), as refgen_common.round_to_double
    # does.
    return math.copysign(r, -1.0) if r == 0 and v < 0 else r


# ---- strata -----------------------------------------------------------------
Args = tuple  # tuple of floats


@dataclass
class Stratum:
    fn: str                                   # probe function name
    name: str                                 # unique within fn
    gen: Callable[[random.Random, int], list] # -> list of Args
    n: int
    oracle: Callable[[Args, int], mpf]        # exact value at dps
    bound: Callable[[Args, mpf], int]         # integer ULP bound for the row
    oracle2: Optional[Callable[[Args, int], mpf]] = None  # independent form
    abs_tol: Optional[Callable[[Args, mpf], Optional[mpf]]] = None
    # abs_tol: when it returns a value t, the row is judged by absolute error
    # |out - exact| <= t instead of ULPs (documented absolute bands).
    tie_exempt: float = 0.0
    # tie_exempt: a documented near-tie exception for a CR (bound 0) cell, in
    # ulps: a 1-ULP miss whose exact value lies within this distance of a
    # rounding midpoint is reported as "tie", not as a violation.
    seam: Optional[str] = None                # registry key it targets (T2)
    ref_cap: Optional[int] = None
    # ref_cap: rows strata_append.py takes into the reference file, where the
    # default per-kind cap (120 for a population) is too few to measure the
    # stratum's not-CR rate (#79's adversarial strata).
    note: str = ""


# ---- probe ------------------------------------------------------------------
def probe_path(build: str) -> Path:
    exe = "hunt_probe.exe" if sys.platform == "win32" else "hunt_probe"
    return REPO / build / "tests" / exe


def evaluate(lines: Sequence[str], tiers=TIERS) -> dict[str, list[float]]:
    out = {}
    payload = "\n".join(lines) + "\n"
    for tier, build in tiers.items():
        p = subprocess.run([str(probe_path(build))], input=payload,
                           capture_output=True, text=True, check=True)
        vals = p.stdout.split()
        if len(vals) != len(lines) or "?" in vals:
            bad = vals.index("?") if "?" in vals else -1
            raise RuntimeError(f"{tier}: probe rejected line {bad}: {lines[bad] if bad >= 0 else ''}")
        out[tier] = [from_bits(int(v, 16)) for v in vals]
    return out


# ---- scoring ----------------------------------------------------------------
_STRATA: dict = {}


def _init(family_names):
    """Workers rebuild the strata by importing the family modules: strata
    hold lambdas, which do not pickle under Windows' spawn start method.
    Each worker first caps its own memory (refgen_common.cap_process_memory)."""
    global _STRATA
    import importlib
    cap_process_memory()
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    _STRATA = {}
    for name in family_names:
        for s in importlib.import_module(f"families.{name}").strata():
            _STRATA[(s.fn, s.name)] = s


def _score(task):
    key, args, outs = task
    s = _STRATA[key]
    try:
        with mp.workdps(40):
            v = +s.oracle(args, 40)
    except Exception as exc:  # oracle refused (cap hit): row is uncertified
        return key, args, None, f"oracle: {type(exc).__name__}: {exc}"[:200]
    res = _judge(s, args, outs, v)
    if res["violation"]:
        # Certify: higher precision, and the independent form if one exists.
        try:
            with mp.workdps(80):
                v2 = +s.oracle(args, 80)
            v3 = None
            if s.oracle2 is not None:
                with mp.workdps(60):
                    v3 = +s.oracle2(args, 60)
            with mp.workdps(90):
                agree = abs(v2 - v) <= abs(v2) * mpf(2) ** -70 + mpf(2) ** -1100
                if v3 is not None:
                    agree = agree and abs(v3 - v2) <= abs(v2) * mpf(2) ** -70 + mpf(2) ** -1100
            res = _judge(s, args, outs, v2)
            res["certified"] = bool(agree)
            res["forms"] = 2 if v3 is not None else 1
        except Exception as exc:
            res["certified"] = False
            res["cert_error"] = str(exc)[:200]
    return key, args, res, None


def _judge(s: Stratum, args, outs, v):
    # Judge at 90 digits: at mpmath's default 53 bits, mpf(v), v/q and the
    # midpoint distance would re-round the oracle value and mis-pick the CR
    # double near ties.
    with mp.workdps(90):
        return _judge_hp(s, args, outs, v)


def _judge_hp(s: Stratum, args, outs, v):
    cr = round_double(v)
    q = quantum(v)
    bound = s.bound(args, v)
    tol = s.abs_tol(args, v) if s.abs_tol else None
    real, ulps = [], []
    for o in outs:
        if math.isnan(o) or math.isinf(o) or not mp.isfinite(v):
            real.append(0.0 if (o == float(v) or (math.isnan(o) and mp.isnan(v))) else math.inf)
            ulps.append(0 if real[-1] == 0.0 else 1 << 62)
            continue
        if tol is not None:
            real.append(float(abs(mpf(o) - v) / tol))  # in units of the tolerance
            ulps.append(0)
        else:
            real.append(float(abs(mpf(o) - v) / q))
            ulps.append(ulp_dist(o, cr))
    if tol is not None:
        viol = any(r > 1.0 for r in real)
        risk = any(r > 1.0 - MARGIN for r in real)
    else:
        viol = any(u > bound for u in ulps)
        # A CR cell (bound 0) has no real-error margin to erode: any miss is a
        # violation, so at-risk there only marks exempted ties (below).
        risk = bound > 0 and any(r > bound + 0.5 - MARGIN for r in real)
    tie = None
    if tol is None and mp.isfinite(v) and v != 0:
        f = abs(v / q)
        tie = float(abs(f - mp.floor(f) - mpf(0.5)))  # distance to midpoint, ulps
        if viol and bound == 0 and max(ulps) <= 1 and tie < s.tie_exempt:
            viol = False
            risk = True
    return {"real": real, "ulp": ulps, "bound": bound, "abs": tol is not None,
            "tie_dist": tie,
            "cr": cr.hex() if math.isfinite(cr) else str(cr),
            "violation": viol, "at_risk": risk}


# mpmath's caches grow with every new precision and argument, so a worker's
# memory creeps upward for its whole life: 6 workers took Z from 3.6 GB free to
# 1.5 GB in 14k beta rows. Recycling workers caps it.
TASKS_PER_WORKER = 400
MIN_FREE_MB = 1536  # stop the pool before the host starves (Z, 16 GB, froze 2026-10-22)


def free_mb() -> Optional[int]:
    """Available physical memory in MB, or None where unsupported."""
    if sys.platform != "win32":
        try:
            with open("/proc/meminfo") as f:
                for line in f:
                    if line.startswith("MemAvailable:"):
                        return int(line.split()[1]) // 1024
        except OSError:
            return None
        return None
    import ctypes

    class MemStatus(ctypes.Structure):
        _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                    ("ullTotalPhys", ctypes.c_ulonglong), ("ullAvailPhys", ctypes.c_ulonglong),
                    ("ullTotalPageFile", ctypes.c_ulonglong), ("ullAvailPageFile", ctypes.c_ulonglong),
                    ("ullTotalVirtual", ctypes.c_ulonglong), ("ullAvailVirtual", ctypes.c_ulonglong),
                    ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]
    st = MemStatus()
    st.dwLength = ctypes.sizeof(MemStatus)
    ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(st))
    return st.ullAvailPhys // 2**20


class LowMemory(RuntimeError):
    pass


def plan_workers(requested: int) -> int:
    """Workers that fit: each may reach its cap (worker_mb), and MIN_FREE_MB
    stays free on top. Raises LowMemory when not even one fits."""
    n = min(requested, MAX_WORKERS)
    fm = free_mb()
    if fm is None:
        return n
    fit = (fm - MIN_FREE_MB) // worker_mb()
    if fit < 1:
        raise LowMemory(f"free memory {fm} MB leaves no room for one {worker_mb()} MB "
                        f"worker above the {MIN_FREE_MB} MB headroom")
    if fit < n:
        print(f"  workers {n} -> {fit} (free {fm} MB, cap {worker_mb()} MB each, "
              f"headroom {MIN_FREE_MB} MB)", file=sys.stderr, flush=True)
    return min(n, fit)


def _key(k, args):
    return f"{k[0]}|{k[1]}|" + ",".join(a.hex() for a in args)


def _pool_map(family_names, tasks, ckpt: Path, workers=MAX_WORKERS, batch=256,
              batch_seconds=900):
    """Score tasks on a worker pool. Every result is appended to the checkpoint
    as it arrives, so an interrupted run resumes where it stopped. A batch has
    a wall-clock limit (a wedged worker costs one batch), and the pool is torn
    down if free memory drops below MIN_FREE_MB."""
    done = {}
    if ckpt.exists():
        for line in ckpt.read_text().splitlines():
            r = json.loads(line)
            done[r["key"]] = (tuple(r["k"]), tuple(float.fromhex(a) for a in r["args"]),
                              r["res"], r["err"])
    todo = [t for t in tasks if _key(t[0], t[1]) not in done]
    print(f"  resume: {len(tasks) - len(todo)} rows from checkpoint, {len(todo)} to score",
          file=sys.stderr, flush=True)
    results = list(done.values())
    workers = plan_workers(workers)
    pool = mp_.Pool(workers, initializer=_init, initargs=(family_names,),
                    maxtasksperchild=TASKS_PER_WORKER)
    fh = ckpt.open("a")

    def keep(r):
        results.append(r)
        fh.write(json.dumps({"key": _key(r[0], r[1]), "k": list(r[0]),
                             "args": [a.hex() for a in r[1]], "res": r[2], "err": r[3]}) + "\n")

    try:
        for i in range(0, len(todo), batch):
            chunk = todo[i:i + batch]
            pend = [(t, pool.apply_async(_score, (t,))) for t in chunk]
            deadline = time.time() + batch_seconds
            while pend:
                still = []
                for t, pr in pend:
                    if pr.ready():
                        keep(pr.get())
                    else:
                        still.append((t, pr))
                pend = still
                if not pend:
                    break
                fm = free_mb()
                if fm is not None and fm < MIN_FREE_MB:
                    pool.terminate()
                    fh.flush()
                    raise LowMemory(f"free memory {fm} MB < {MIN_FREE_MB} MB; pool stopped. "
                                    f"Next pending row: {pend[0][0][0]} {[a.hex() for a in pend[0][0][1]]}")
                if time.time() > deadline:
                    for t, _ in pend:
                        keep((t[0], t[1], None, "oracle: wall-clock timeout"))
                    pool.terminate()
                    pool = mp_.Pool(workers, initializer=_init, initargs=(family_names,),
                    maxtasksperchild=TASKS_PER_WORKER)
                    break
                time.sleep(0.2)
            fh.flush()
            print(f"  scored {len(results)}/{len(tasks)}  free {free_mb()} MB",
                  file=sys.stderr, flush=True)
    finally:
        pool.terminate()
        fh.close()
    return results


def run(strata: list[Stratum], family_names: list[str], seed: int, out_json: Path,
        tiers=TIERS, workers=MAX_WORKERS):
    rng = random.Random(seed)
    rows, lines = [], []
    for s in strata:
        pts = s.gen(random.Random(rng.getrandbits(64)), s.n)
        for args in pts:
            rows.append(((s.fn, s.name), tuple(float(a) for a in args)))
            lines.append(s.fn + " " + " ".join(float(a).hex() for a in args))
    t0 = time.time()
    outs = evaluate(lines, tiers)
    tasks = [(k, a, [outs[t][i] for t in tiers]) for i, (k, a) in enumerate(rows)]
    print(f"evaluated {len(rows)} rows on {len(tiers)} tiers in {time.time()-t0:.0f}s",
          file=sys.stderr)
    scored = _pool_map(family_names, tasks, out_json.with_suffix(".ckpt.jsonl"), workers)
    return summarise(strata, scored, list(tiers), out_json)


def summarise(strata, scored, tiers, out_json: Path):
    by = {}
    for key, args, res, err in scored:
        b = by.setdefault(key, {"n": 0, "uncert": 0, "max_real": [0.0] * len(tiers),
                                "max_ulp": [0] * len(tiers), "ncr": [0] * len(tiers),
                                "violations": [], "at_risk": [], "errors": []})
        b["n"] += 1
        if res is None:
            b["uncert"] += 1
            if len(b["errors"]) < 3:
                b["errors"].append({"args": [a.hex() for a in args], "err": err})
            continue
        for k in range(len(tiers)):
            b["max_real"][k] = max(b["max_real"][k], res["real"][k])
            b["max_ulp"][k] = max(b["max_ulp"][k], res["ulp"][k])
            b["ncr"][k] += res["ulp"][k] > 0
        rec = {"args": [a.hex() for a in args], **res}
        if res["violation"]:
            b["violations"].append(rec)
        elif res["at_risk"] and len(b["at_risk"]) < 20:
            b["at_risk"].append(rec)
    report = {"tiers": tiers, "strata": []}
    for s in strata:
        b = by.get((s.fn, s.name))
        if b is None:
            continue
        b["violations"].sort(key=lambda r: -max(r["real"]))
        report["strata"].append({"fn": s.fn, "name": s.name, "seam": s.seam,
                                 "note": s.note, **b})
    out_json.write_text(json.dumps(report, indent=1))
    return report


def print_table(report):
    tiers = report["tiers"]
    print(f"{'fn':12s}{'stratum':34s}{'n':>7s}{'unc':>5s}  "
          + "  ".join(f"{t[:9]:>9s}" for t in tiers) + "   max ULP   flags")
    for s in report["strata"]:
        flags = []
        nv = sum(1 for v in s["violations"] if v.get("certified", False))
        nu = len(s["violations"]) - nv
        if nv:
            flags.append(f"VIOLATION x{nv}")
        if nu:
            flags.append(f"uncertified-viol x{nu}")
        if s["at_risk"]:
            flags.append("at-risk")
        print(f"{s['fn']:12s}{s['name'][:33]:34s}{s['n']:7d}{s['uncert']:5d}  "
              + "  ".join(f"{r:9.3f}" for r in s["max_real"])
              + f"   {str(s['max_ulp']):10s}{' '.join(flags)}")
