"""Seam audit (#82, T2): does every registered kernel seam have reference rows?

For each SEAMS entry in tools/hunt/families/, evaluate its coordinate on every
row of the function's checked-in reference files and count
  bits   rows within 64 ulp of the seam,
  below  rows in [t(1-d), t) outside the bit neighbourhood,
  above  rows in (t, t(1+d)] outside the bit neighbourhood   (d = 5%).
A seam with an empty bucket is a population gap: the review's misses (#70's
a*phi band in [34, 36), for one) were exactly such gaps. Coordinates the audit
cannot evaluate from a reference row are listed as UNMAPPED -- those entries
still need a machine-readable coordinate.

  python tools/hunt/seam_audit.py [--delta 0.05] [--only-gaps]
"""
from __future__ import annotations

import argparse
import importlib
import math
import pkgutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import families  # noqa: E402

DATA = Path(__file__).resolve().parents[2] / "tests" / "data"


def phi(lam):
    return lam - 1 - math.log(lam) if lam > 0 else math.inf


# fn (as SEAMS names it) -> (reference files, column names)
FILES = {
    "erf": (["erf"], ["x"]), "erfc": (["erfc"], ["x"]),
    "erfinv": (["erfinv"], ["y"]), "erfcinv": (["erfcinv"], ["z"]),
    "lgamma": (["lgamma"], ["x"]), "digamma": (["digamma"], ["x"]),
    "trigamma": (["trigamma"], ["x"]),
    "i0": (["i0"], ["x"]), "i1": (["i1"], ["x"]), "i0e": (["i0e"], ["x"]), "i1e": (["i1e"], ["x"]),
    "exp": (["exp"], ["x"]), "log": (["log"], ["x"]), "log1p": (["log1p"], ["x"]),
    "cos": (["cos"], ["x"]), "sin": (["sin"], ["x"]),
    "gamma_p/q": (["gamma_p", "gamma_q"], ["a", "x"]),
    "gamma_p": (["gamma_p"], ["a", "x"]), "gamma_q": (["gamma_q"], ["a", "x"]),
    "gamma_p/q_inv": (["gammainv_p", "gammainv_q"], ["a", "s", "root x"]),
    "beta": (["beta_p", "beta_q"], ["a", "b", "x"]),
    "beta_p": (["beta_p"], ["a", "b", "x"]), "beta_q": (["beta_q"], ["a", "b", "x"]),
    "betainv": (["betainv_p", "betainv_q"], ["a", "b", "s", "root x"]),
    "lbeta": (["lbeta"], ["a", "b"]),
}

# coordinate expression -> function of the named columns (None when undefined)
COORDS = {
    "x": lambda r: r.get("x"), "y": lambda r: r.get("y"), "z": lambda r: r.get("z"),
    "a": lambda r: r.get("a"), "b": lambda r: r.get("b"), "s": lambda r: r.get("s"),
    "p": lambda r: r.get("s"),
    "root x": lambda r: r.get("root x"),
    "1-x": lambda r: 1 - r["x"],
    "x/a": lambda r: r["x"] / r["a"] if r["a"] else None,
    "a*phi": lambda r: r["a"] * phi(r["x"] / r["a"]) if r["a"] > 0 and r["x"] > 0 else None,
    "x-(a+1)": lambda r: r["x"] - (r["a"] + 1),
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
            vals = [float.fromhex(t) if t.startswith(("0x", "-0x")) else float(t) for t in f[:len(cols)]]
        except ValueError:
            continue
        if len(vals) == len(cols):
            rows.append(dict(zip(cols, vals)))
    return rows


def ulp(v):
    v = abs(v)
    return math.ulp(v) if v else 5e-324


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--delta", type=float, default=0.05)
    ap.add_argument("--only-gaps", action="store_true")
    a = ap.parse_args()
    cache = {}
    gaps = unmapped = total = 0
    print(f"{'fn':14s}{'seam':28s}{'coord':16s}{'value':>14s}{'bits':>7s}{'below':>7s}{'above':>7s}")
    for m in pkgutil.iter_modules(families.__path__):
        if m.name.startswith("_"):
            continue
        for s in importlib.import_module(f"families.{m.name}").SEAMS:
            total += 1
            fn, coord, t = s["fn"], s["coord"], s.get("value")
            f = COORDS.get(coord)
            spec = FILES.get(fn)
            if f is None or spec is None or not isinstance(t, (int, float)) or not math.isfinite(t):
                unmapped += 1
                if not a.only_gaps:
                    print(f"{fn:14s}{s['name'][:27]:28s}{coord[:15]:16s}{'':>14s}   UNMAPPED")
                continue
            files, cols = spec
            rows = cache.setdefault(fn, [r for name in files for r in load_rows(name, cols)])
            nb = lo = hi = 0
            win = 64 * ulp(t)
            for r in rows:
                try:
                    c = f(r)
                except (KeyError, ZeroDivisionError, ValueError, OverflowError):
                    c = None
                if c is None or not math.isfinite(c):
                    continue
                if abs(c - t) <= win:
                    nb += 1
                elif t - abs(t) * a.delta <= c < t:
                    lo += 1
                elif t < c <= t + abs(t) * a.delta:
                    hi += 1
            gap = (nb == 0) or (lo == 0) or (hi == 0)
            gaps += gap
            if gap or not a.only_gaps:
                print(f"{fn:14s}{s['name'][:27]:28s}{coord[:15]:16s}{t:14.6g}{nb:7d}{lo:7d}{hi:7d}"
                      + ("   GAP" if gap else ""))
    print(f"\n{total} seams: {gaps} with an empty bucket, {unmapped} unmapped")


if __name__ == "__main__":
    main()
