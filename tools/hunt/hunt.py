"""Run hunt strata on every built tier and report real-error profiles (#82, T5).

  python tools/hunt/hunt.py [--families erf,gamma] [--scale 1.0] [--seed N]
                            [--out hunt-report.json] [--only-flagged]

Each module in tools/hunt/families/ exports
  SEAMS: list[dict]   the T2 seam registry for its functions
                      (keys: fn, name, coord, value, source)
  strata() -> list[core.Stratum]
Build hunt_probe in each tier's build directory first (core.TIERS).
"""
from __future__ import annotations

import argparse
import importlib
import pkgutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import core  # noqa: E402
import families  # noqa: E402


def load(names):
    mods = []
    for m in pkgutil.iter_modules(families.__path__):
        if m.name.startswith("_") or (names and m.name not in names):  # _ = private helper
            continue
        mods.append(importlib.import_module(f"families.{m.name}"))
    return mods


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--families", default="")
    ap.add_argument("--scale", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=20261022)
    ap.add_argument("--out", default="hunt-report.json")
    ap.add_argument("--only-flagged", action="store_true")
    ap.add_argument("--workers", type=int, default=core.MAX_WORKERS)
    a = ap.parse_args()
    names = [n for n in a.families.split(",") if n]
    strata = []
    mods = load(names)
    for m in mods:
        for s in m.strata():
            s.n = max(1, int(s.n * a.scale))
            strata.append(s)
    seen = set()
    for s in strata:
        if (s.fn, s.name) in seen:
            sys.exit(f"duplicate stratum {s.fn}/{s.name}")
        seen.add((s.fn, s.name))
    print(f"{len(strata)} strata, {sum(s.n for s in strata)} rows", file=sys.stderr)
    rep = core.run(strata, [m.__name__.split(".")[-1] for m in mods], a.seed, Path(a.out),
                   workers=min(a.workers, core.MAX_WORKERS))
    if a.only_flagged:
        rep["strata"] = [s for s in rep["strata"] if s["violations"] or s["at_risk"]]
    core.print_table(rep)


if __name__ == "__main__":
    main()
