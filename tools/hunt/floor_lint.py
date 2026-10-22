"""Floor lint (#82, T4): every Dekker-product site states why it is exact.

ops::ProdLow / SquareLow -- and so TwoProd, DdMul, DdMulD and the other dd
products built on them -- are exact only between a ceiling (|a|, |b| < 2^996)
and a floor (|a*b| >~ 2^-968); see src/ops-inl.h. Below the floor the no-FMA
tiers lose the low product (M-5/#73: erfinv's subnormal inputs, 2 ULP on
SSE4/SSE2). The rule in ops-inl.h -- "every caller must sit between the two
or prescale" -- had no check, so this lint makes the argument visible:

  // FLOOR: <why the operands stay inside [floor, ceiling], or the prescale>

either on the call's line or in the comment lines directly above it, or once
at the head of the enclosing function (between the previous function's end and
this function's signature) to cover every site in it. The lint lists every
uncovered site; --strict exits non-zero when any remain.

The rule, what counts as an argument, and the phasing live in
docs/NUMERICAL-DOCTRINE.md (hazard rules, "Every Dekker-product site states
its floor argument"). In short: row 12-T argues the sites that can reach the
floor, each fix row argues the functions it changes, and #86 (v1.2.0) argues
the rest. The lint stays report-only until it reaches zero; then CI runs it
--strict.

  python tools/hunt/floor_lint.py [--strict] [--summary]
"""
from __future__ import annotations

import argparse
import re
import sys
from collections import Counter
from pathlib import Path

SRC = Path(__file__).resolve().parents[2] / "src"
CALL = re.compile(r"\b(TwoProd|TwoSquare|ProdLow|SquareLow|DdMul|DdMulD|DdMulDd|DdSquare|DdSqr)\(")
DEFINES = {"dd-inl.h", "ops-inl.h"}  # the primitives themselves
FN_HEAD = re.compile(r"^(template\b|HWY_(NOINLINE|INLINE|ATTR)|CORVUS_OUTLINE|inline\b|static\b)")
MARK = "FLOOR:"


def scan(path: Path):
    lines = path.read_text(encoding="utf-8").splitlines()
    head_mark = False  # a FLOOR: in the current function's leading comment
    pending = False    # FLOOR: seen in comments since the last code line
    for i, line in enumerate(lines, 1):
        s = line.strip()
        if FN_HEAD.match(line):
            head_mark = pending
        if s.startswith("//"):
            pending = pending or MARK in s
            continue
        if s == "}" and not line.startswith(" "):
            head_mark = False
        m = CALL.search(line)
        if m and not s.startswith("#"):
            covered = MARK in line or head_mark or pending
            yield i, m.group(1), covered, s
        if s:
            pending = False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--strict", action="store_true")
    ap.add_argument("--summary", action="store_true")
    a = ap.parse_args()
    total = uncovered = 0
    per_file = Counter()
    for path in sorted(SRC.glob("*.h")):
        if path.name in DEFINES:
            continue
        for ln, prim, covered, text in scan(path):
            total += 1
            if not covered:
                uncovered += 1
                per_file[path.name] += 1
                if not a.summary:
                    print(f"src/{path.name}:{ln}: {prim}: no FLOOR: argument   {text[:80]}")
    if a.summary:
        for f, n in per_file.most_common():
            print(f"{f:28s}{n:4d} uncovered")
    print(f"{total} Dekker-product sites, {uncovered} without a FLOOR: argument")
    return 1 if (a.strict and uncovered) else 0


if __name__ == "__main__":
    sys.exit(main())
