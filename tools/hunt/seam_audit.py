"""Seam audit (#82, T2): does every registered kernel seam have reference rows?

For each SEAMS entry in tools/hunt/families/, evaluate its coordinate on every
row of the function's checked-in reference files and count
  bits   rows within 64 ulp of the seam,
  below  rows in [t(1-d), t) outside the bit neighbourhood,
  above  rows in (t, t(1+d)] outside the bit neighbourhood   (d = 5%).
A seam with an empty bucket is a population gap: the review's misses (#70's
a*phi band in [34, 36), for one) were exactly such gaps.

A SEAMS entry may carry audit-only fields (#82 T2 completion):
  audit_f      row -> coordinate (None where undefined), for coordinates the
               generic COORDS table cannot express (router orientation,
               region membership, a ratio to an expression-valued seam)
  audit_value  the numeric threshold (or a list of them) when `value` is an
               expression or approximate
  window       the bits window when the coordinate is coarser than 64 ulp of
               the threshold (a difference of a double and an exact value)
  sides        "below" / "above" when the other side is outside the domain
               (a bucket narrower than one double is also skipped)
  route        a census route (tools/hunt/route_census.py) that covers a
               threshold on a loop-internal quantity no row can evaluate
  internal     the reason no row coordinate exists, when not even a route
               marks the threshold
Entries with `route` or `internal` are listed as CENSUS / INTERNAL, not as
gaps; any other entry without a usable coordinate is UNMAPPED.

  python tools/hunt/seam_audit.py [--delta 0.05] [--only-gaps]
"""
from __future__ import annotations

import argparse
import importlib
import math
import pkgutil
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import families  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "tests" / "data"


def phi(lam):
    return lam - 1 - math.log(lam) if lam > 0 else math.inf


# fn (as SEAMS names it) -> (reference files, column names; outputs included)
FILES = {
    "erf": (["erf"], ["x", "y"]), "erfc": (["erfc"], ["x", "y"]),
    "erfinv": (["erfinv"], ["y", "x"]), "erfcinv": (["erfcinv"], ["z", "x"]),
    "lgamma": (["lgamma"], ["x", "y"]), "digamma": (["digamma"], ["x", "y"]),
    "trigamma": (["trigamma"], ["x", "y"]),
    "i0": (["i0"], ["x", "y"]), "i1": (["i1"], ["x", "y"]),
    "i0e": (["i0e"], ["x", "y"]), "i1e": (["i1e"], ["x", "y"]),
    "exp": (["exp"], ["x", "y"]), "log": (["log"], ["x", "y"]), "log1p": (["log1p"], ["x", "y"]),
    "cos": (["cos"], ["x", "y"]), "sin": (["sin"], ["x", "y"]),
    "gamma_p/q": (["gamma_p", "gamma_q"], ["a", "x", "P", "Q"]),
    "gamma_p": (["gamma_p"], ["a", "x", "P", "Q"]), "gamma_q": (["gamma_q"], ["a", "x", "P", "Q"]),
    "gamma_p/q_inv": (["gammainv_p", "gammainv_q"], ["a", "s", "root x"]),
    "beta": (["beta_p", "beta_q"], ["a", "b", "x", "P", "Q"]),
    "beta_p": (["beta_p"], ["a", "b", "x", "P", "Q"]), "beta_q": (["beta_q"], ["a", "b", "x", "P", "Q"]),
    "betainv": (["betainv_p", "betainv_q"], ["a", "b", "s", "root x"]),
    "lbeta": (["lbeta"], ["a", "b", "lnB"]),
}

# coordinate expression -> function of the named columns (None when undefined)
COORDS = {
    "x": lambda r: r.get("x", r.get("root x")), "y": lambda r: r.get("y"),
    "z": lambda r: r.get("z"),
    "a": lambda r: r.get("a"), "b": lambda r: r.get("b"), "s": lambda r: r.get("s"),
    "p": lambda r: r.get("s"),
    "root x": lambda r: r.get("root x"),
    "1-x": lambda r: 1 - r["x"],
    "x/a": lambda r: r["x"] / r["a"] if r["a"] else None,
    "a*phi": lambda r: r["a"] * phi(r["x"] / r["a"]) if r["a"] > 0 and r["x"] > 0 else None,
    "x*(1+a)": lambda r: r["root x"] * (1 + r["a"]),
    "log x": lambda r: math.log(r["root x"]) if r["root x"] > 0 else None,
    "a+b": lambda r: r["a"] + r["b"],
    "b*x": lambda r: r["b"] * r.get("x", r.get("root x", math.nan)),
    "min(a,b)": lambda r: min(r["a"], r["b"]),
    "max(a,b)": lambda r: max(r["a"], r["b"]),
    "p=a/(a+b)": lambda r: r["a"] / (r["a"] + r["b"]),
}


def load_rows(name, cols):
    rows = []
    path = DATA / f"{name}_reference.txt"
    if not path.exists():
        return rows
    for line in path.read_text().splitlines():
        if not line or line[0] in "#;":
            continue
        f = line.split()
        try:
            vals = [float.fromhex(t) if t.startswith(("0x", "-0x")) else float(t)
                    for t in f[:len(cols)]]
        except ValueError:
            continue
        if len(vals) == len(cols):
            row = dict(zip(cols, vals))
            row["file"] = name
            rows.append(row)
    return rows


def ulp(v):
    v = abs(v)
    return math.ulp(v) if v else 5e-324


def census_routes():
    pat = re.compile(r'CORVUS_ROUTE\s*\(\s*d\w*\s*,\s*"([^"]+)"')
    names = set()
    for p in (REPO / "src").glob("*.h"):
        for line in p.read_text(encoding="utf-8", errors="replace").splitlines():
            if not line.lstrip().startswith("//"):
                names.update(pat.findall(line))
    return names


def buckets(rows, f, t, delta, win=None):
    nb = lo = hi = 0
    win = win if win is not None else 64 * ulp(t)
    for r in rows:
        try:
            c = f(r)
        except (KeyError, ZeroDivisionError, ValueError, OverflowError, TypeError):
            c = None
        if c is None or not math.isfinite(c):
            continue
        if abs(c - t) <= win:
            nb += 1
        elif t - abs(t) * delta <= c < t:
            lo += 1
        elif t < c <= t + abs(t) * delta:
            hi += 1
    return nb, lo, hi


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--delta", type=float, default=0.05)
    ap.add_argument("--only-gaps", action="store_true")
    a = ap.parse_args()
    cache = {}
    routes = census_routes()
    gaps = unmapped = covered = total = 0
    print(f"{'fn':14s}{'seam':28s}{'coord':16s}{'value':>14s}{'bits':>7s}{'below':>7s}{'above':>7s}")
    for m in pkgutil.iter_modules(families.__path__):
        if m.name.startswith("_"):
            continue
        for s in importlib.import_module(f"families.{m.name}").SEAMS:
            fn, coord = s["fn"], s["coord"]
            label = f"{fn:14s}{s['name'][:27]:28s}{coord[:15]:16s}"
            if "route" in s or "internal" in s:
                total += 1
                covered += 1
                if "route" in s:
                    ok = s["route"] in routes
                    gaps += not ok
                    if not a.only_gaps or not ok:
                        print(f"{label}{'':>14s}   CENSUS {s['route']}"
                              + ("" if ok else "   ROUTE MISSING"))
                elif not a.only_gaps:
                    print(f"{label}{'':>14s}   INTERNAL: {s['internal'][:60]}")
                continue
            f = s.get("audit_f") or COORDS.get(coord)
            vals = s.get("audit_value", s.get("value"))
            vals = vals if isinstance(vals, list) else [vals]
            spec = FILES.get(fn)
            for t in vals:
                total += 1
                if f is None or spec is None or not isinstance(t, (int, float)) or not math.isfinite(t):
                    unmapped += 1
                    if not a.only_gaps:
                        print(f"{label}{'':>14s}   UNMAPPED")
                    continue
                files, cols = spec
                rows = cache.setdefault(fn, [r for name in files for r in load_rows(name, cols)])
                nb, lo, hi = buckets(rows, f, t, a.delta, s.get("window"))
                sides = s.get("sides", "both")
                # A band narrower than one double holds no row outside the bits.
                need_lo = sides in ("both", "below") and abs(t) * a.delta > 64 * ulp(t)
                need_hi = sides in ("both", "above") and abs(t) * a.delta > 64 * ulp(t)
                gap = (nb == 0) or (need_lo and lo == 0) or (need_hi and hi == 0)
                gaps += gap
                if gap or not a.only_gaps:
                    note = "" if sides == "both" else f"  ({sides} only)"
                    print(f"{label}{t:14.6g}{nb:7d}{lo:7d}{hi:7d}" + ("   GAP" if gap else "") + note)
    print(f"\n{total} seam thresholds: {gaps} with an empty bucket, {unmapped} unmapped, "
          f"{covered} covered by the census or internal")
    return 1 if gaps or unmapped else 0


if __name__ == "__main__":
    sys.exit(main())
