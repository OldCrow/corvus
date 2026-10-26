"""Regenerate tests/data/real_pins.txt from measured real error (#82, T1).

Runs every ULP gate (ctest -R _ulp) in each build with
CORVUS_REAL_PIN_REPORT=1, collects the `PIN gate|region|class|measured|
contract` lines, and writes one pin per gate, region and tier class:

  pin = measured maximum over the class's tiers + DELTA, rounded up to 0.01

DELTA absorbs compiler-to-compiler differences inside a tier class (the same
tier built by clang-cl, GCC and Apple Clang). Extra logs (--log, e.g. CI's
ULP report for NEON and AVX2) join the maximum. Entries the run did not
measure are kept, and so is every entry's issue field.

A pin above its contract ceiling must name the issue that tracks the defect
(see CheckReal in tests/ulp_utils.h): such entries are listed as NEEDS ISSUE
and the file is written with an `?` in the issue field, which the gate
rejects until a person names the issue.

  python tools/hunt/pins.py [--builds DIR ...] [--log FILE ...] [--write]
"""
from __future__ import annotations

import argparse
import math
import os
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
PINS = REPO / "tests" / "data" / "real_pins.txt"
DELTA = 0.02
LINE = re.compile(r"PIN (.+)\|([-0-9.e+inf]+)\|([0-9.]+)\s*$")
HEADER = ("# Real-error regression pins (#82 T1), read by CheckReal in tests/ulp_utils.h.\n"
          "# gate|region|tier class|pin[|issue]. Regenerate with tools/hunt/pins.py;\n"
          "# a pin above the contract ceiling names the issue that tracks it.\n")


def collect(text: str, meas: dict, contract: dict):
    for line in text.splitlines():
        m = LINE.search(line)
        if m:
            key, v, c = m.group(1), float(m.group(2)), float(m.group(3))
            meas[key] = max(meas.get(key, 0.0), v)
            contract[key] = c


def load():
    out = {}
    if PINS.exists():
        for line in PINS.read_text().splitlines():
            if line and not line.startswith("#"):
                f = line.split("|")
                out["|".join(f[:3])] = (float(f[3]), f[4] if len(f) > 4 else "")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--builds", nargs="*",
                    default=["build-clangcl", "build-guards-sse4", "build-guards-sse2"])
    ap.add_argument("--log", nargs="*", default=[])
    ap.add_argument("--write", action="store_true")
    a = ap.parse_args()
    meas, contract = {}, {}
    env = dict(os.environ, CORVUS_REAL_PIN_REPORT="1")
    for b in a.builds:
        p = subprocess.run(["ctest", "--test-dir", str(REPO / b), "-R", "_ulp", "-V"],
                           capture_output=True, text=True, env=env)
        collect(p.stdout, meas, contract)
        print(f"{b}: {sum(1 for l in p.stdout.splitlines() if 'PIN ' in l)} pin lines",
              file=sys.stderr)
    for f in a.log:
        collect(Path(f).read_text(errors="replace"), meas, contract)
    old = load()
    rows, need = dict(old), []
    for key, v in meas.items():
        pin = math.ceil((v + DELTA) * 100 - 1e-9) / 100
        if v <= contract[key]:
            pin = min(pin, contract[key])  # inside the contract: it is the cap
        issue = old.get(key, (0, ""))[1]
        if pin > contract[key] and not issue:
            issue = "?"
            need.append(f"{key}: measured {v:.4f} > contract {contract[key]:.2f}")
        rows[key] = (pin, issue)
    for key in sorted(rows):
        o = old.get(key)
        if key in meas and (o is None or abs(o[0] - rows[key][0]) > 1e-9):
            print(f"  {key}: {o[0] if o else 'none'} -> {rows[key][0]:.2f}", file=sys.stderr)
    for n in need:
        print("NEEDS ISSUE " + n, file=sys.stderr)
    if a.write:
        body = "".join(f"{k}|{p:.2f}" + (f"|{i}" if i else "") + "\n"
                       for k, (p, i) in sorted(rows.items()))
        PINS.write_text(HEADER + body)
    return 0


if __name__ == "__main__":
    sys.exit(main())
