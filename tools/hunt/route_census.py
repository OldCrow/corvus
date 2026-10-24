"""Route census (#82, T7): is every kernel route exercised, and does every
row take the route the hunt's router replica predicts?

Kernels mark routes with CORVUS_ROUTE(d, "family/route", mask) and region
blocks with CORVUS_ROUTE_SCOPE (src/route_census.h). The route list comes
from the source, so a route no row reaches is still listed. A census build
(-DCORVUS_ROUTE_CENSUS=ON) of tests/hunt_probe prints the routes lane 0 took
for each row; every row of every public reference file -- deferred rows
included, they are real inputs -- is run through it, once per output
function of its file.

Reported:
  GRID-ONLY    a non-special route only the special-input grid reaches: it
               runs, but no gate checks its result; it needs a stratum.
  UNEXERCISED  a route no reference row takes. Either add a stratum that
               reaches it (a hunt family, then strata_append.py --replace),
               or argue it unreachable in a `// ROUTE-UNREACHABLE: <why>`
               comment on the line above the marker; the census lists those
               separately and does not count them as gaps.
  MISMATCH     a row whose top-level region differs from the replica's
               (tools/hunt/families: gamma.region, beta.route). A mismatch
               is either a replica bug or a router defect (#84's class).
Exit status 1 when either list is non-empty.

  python tools/hunt/route_census.py [--build build-census] [--keys KEY ...]
                                    [--json out.json]
"""
from __future__ import annotations

import argparse
import collections
import json
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import annotate  # noqa: E402
from core import REPO  # noqa: E402
from families import beta as B, gamma as G  # noqa: E402

SRC = REPO / "src"
WITNESSES = Path(__file__).resolve().parent / "route_witnesses.txt"
MARK = re.compile(r'CORVUS_ROUTE\s*\(\s*d\w*\s*,\s*"([^"]+)"')
UNREACH = re.compile(r"//\s*ROUTE-UNREACHABLE:\s*(.*)")


def source_routes():
    """route -> (file:line, unreachable argument or None)."""
    out = {}
    for p in sorted(SRC.glob("*.h")):
        lines = p.read_text(encoding="utf-8", errors="replace").splitlines()
        for i, line in enumerate(lines):
            if line.lstrip().startswith("//"):
                continue  # documentation examples (route_census.h, ops-inl.h)
            for name in MARK.findall(line):
                # The comment block directly above the marker; an argument
                # starts at `ROUTE-UNREACHABLE:` and runs to the block's end.
                j = i
                while j > 0 and lines[j - 1].lstrip().startswith("//"):
                    j -= 1
                block = [ln.lstrip()[2:].strip() for ln in lines[j:i]]
                why = None
                for k, text in enumerate(block):
                    m = UNREACH.search("//" + text)
                    if m:
                        why = " ".join([m.group(1).strip()] + block[k + 1:]).strip()
                out.setdefault(name, (f"{p.name}:{i + 1}", why))
    return out


# ---- replicas: fn -> (route prefix, args -> acceptable top-level labels) ----
def _gamma(args):
    r = G.region(*args)
    # The replica does not model R4 (a subset of R1 u R2: the small-a box).
    return {r, "R4"} if r in ("R1", "R2") else {r}


def _beta(args):
    return {reg for reg, _ in B.route(*args)}


REPLICAS = {
    "gamma_p": ("gamma/", _gamma), "gamma_q": ("gamma/", _gamma),
    "beta_p": ("beta/", _beta), "beta_q": ("beta/", _beta),
}
TOP = re.compile(r"^(R\d\w*|GL)$")


# Special inputs: the reference files hold none (the smoke tests own those
# contracts), so the census adds this grid per argument for the
# special-value routes. Census rows need no expected value.
SPECIALS = [float("nan"), float("inf"), float("-inf"), 0.0, -0.0, 1.0, -1.0,
            0.5, 2.0, -2.5, 5e-324, 1e300]


def specials_of(nargs):
    import itertools
    return [tuple(p) for p in itertools.product(SPECIALS, repeat=nargs)]


def witnesses():
    """[(route, fn, args, why)] from route_witnesses.txt."""
    out = []
    for line in WITNESSES.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        route, fn, args, why = line.split("|", 3)
        out.append((route, fn, tuple(float.fromhex(v) if "x" in v else float(v)
                                     for v in args.split()), why))
    return out


def rows_of(key):
    """[(args tuple)] for every data row of the file, deferred rows included."""
    _, outs = annotate.FILES[key]
    nargs = len(outs[0][1])
    rows = []
    for line in annotate.path_of(key).read_text().splitlines():
        tok = annotate._data_tokens(line)
        if tok:
            rows.append(tuple(float.fromhex(t) if "x" in t or t in ("inf", "-inf", "nan")
                              else float(t) for t in tok[:nargs]))
    return rows


def run_probe(exe: Path, lines):
    p = subprocess.run([str(exe)], input="\n".join(lines) + "\n",
                       capture_output=True, text=True, check=True)
    out = p.stdout.splitlines()
    if len(out) != len(lines):
        raise RuntimeError(f"probe returned {len(out)} lines for {len(lines)}")
    routes = []
    for line in out:
        parts = line.split(" ", 1)
        if len(parts) != 2:
            raise RuntimeError(f"not a census build? line '{line}'")
        routes.append([] if parts[1].strip() == "-" else parts[1].strip().split(","))
    return routes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", default="build-census")
    ap.add_argument("--keys", nargs="*", default=list(annotate.FILES))
    ap.add_argument("--json")
    a = ap.parse_args()
    exe = REPO / a.build / "tests" / ("hunt_probe.exe" if sys.platform == "win32" else "hunt_probe")
    src = source_routes()
    hits = collections.Counter()
    ref_hits = collections.Counter()  # from reference rows, which the gates check
    examples = {}
    mismatches = []
    done = set()  # (fn, row): the P/Q files share rows and functions
    for key in a.keys:
        _, outs = annotate.FILES[key]
        for fn, _, _ in outs:
            ref = set(rows_of(key))
            rows = [r for r in dict.fromkeys(list(ref) + specials_of(len(outs[0][1])))
                    if (fn, r) not in done]
            done.update((fn, r) for r in rows)
            if not rows:
                continue
            lines = [fn + " " + " ".join(v.hex() for v in r) for r in rows]
            routes = run_probe(exe, lines)
            for r, rs in zip(rows, routes):
                for name in set(rs):
                    hits[name] += 1
                    if r in ref:
                        ref_hits[name] += 1
                    examples.setdefault(name, (fn, r))
                if fn in REPLICAS and all(map(lambda v: v == v and abs(v) != float("inf"), r)):
                    prefix, rep = REPLICAS[fn]
                    got = {n[len(prefix):] for n in rs
                           if n.startswith(prefix) and TOP.match(n[len(prefix):])}
                    if not got:
                        continue  # a special row: no region taken
                    try:
                        want = rep(r)
                    except (ZeroDivisionError, ValueError, OverflowError):
                        continue  # outside the replica's domain (zero shape)
                    if not (got & want):
                        mismatches.append((fn, r, sorted(got), sorted(want)))
            print(f"{key}: {fn} {len(rows)} rows", file=sys.stderr, flush=True)

    # Witnesses: each must still take its route.
    wit = witnesses()
    lost = []
    for route, fn, args, _ in wit:
        rs = run_probe(exe, [fn + " " + " ".join(v.hex() for v in args)])[0]
        if route in rs:
            hits[route] += 1
            examples.setdefault(route, (fn, args))
        else:
            lost.append((route, fn, args))
    witnessed = {w[0] for w in wit}

    unknown = sorted(set(hits) - set(src))
    gaps = sorted(n for n, (_, why) in src.items() if hits[n] == 0 and not why)
    argued = sorted(n for n, (_, why) in src.items() if hits[n] == 0 and why)
    reached_argued = sorted(n for n, (_, why) in src.items() if hits[n] and why)
    # A non-special route reached only by the special-input grid runs, but no
    # gate checks its result: it needs a reference row (a stratum).
    grid_only = sorted(n for n in src if hits[n] and not ref_hits[n] and "/special/" not in n
                       and n not in witnessed)

    print(f"\n{len(src)} routes in the source; {sum(1 for n in src if hits[n])} exercised")
    for n in sorted(src):
        print(f"  {hits[n]:8d}  {n:44s} {src[n][0]}")
    for n in gaps:
        print(f"UNEXERCISED {n}  ({src[n][0]})")
    for route, fn, args in lost:
        print(f"WITNESS LOST {route}: {fn}({', '.join(v.hex() for v in args)}) no longer takes it")
    for n in grid_only:
        print(f"GRID-ONLY {n}  ({src[n][0]}): no reference row; e.g. {examples[n]}")
    for n in argued:
        print(f"argued unreachable: {n}  ({src[n][0]}): {src[n][1]}")
    for n in reached_argued:
        print(f"REACHED BUT ARGUED UNREACHABLE {n}: {hits[n]} rows, e.g. {examples[n]}")
    for n in unknown:
        print(f"recorded but not in the source list (parser gap): {n}")
    for fn, r, got, want in mismatches[:40]:
        print(f"MISMATCH {fn}({', '.join(v.hex() for v in r)}): kernel {got}, replica {want}")
    print(f"{len(gaps)} unexercised, {len(grid_only)} grid-only, {len(mismatches)} mismatches, "
          f"{len(argued)} argued unreachable")
    if a.json:
        Path(a.json).write_text(json.dumps({
            "hits": dict(hits), "ref_hits": dict(ref_hits), "gaps": gaps,
            "grid_only": grid_only, "argued": argued,
            "mismatches": [[f, [v.hex() for v in r], g, w] for f, r, g, w in mismatches],
            "examples": {k: [f, [v.hex() for v in r]] for k, (f, r) in examples.items()},
        }, indent=1))
    return 1 if gaps or grid_only or lost or mismatches or reached_argued or unknown else 0


if __name__ == "__main__":
    sys.exit(main())
