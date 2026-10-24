#!/usr/bin/env python3
"""Claims check (#82 T6): do the docs state the bounds the gates enforce?

tests/data/bounds.txt holds every gated region as
  gate|region|unit|gate|published|claim[|issue]
The gates must pass exactly its gate field (tests/ulp_utils.h CheckBound).
This script holds the docs to its published field:

- Table: each claim id has one published value and unit across its rows,
  and published >= gate on every row unless the row names an issue.
- Markers: an HTML comment `<!-- claim: <id> [<id> ...] [#NN] -->` (also
  inside `///` blocks in corvus.h, where Doxygen drops it) sits directly
  before the fragment that states a bound. The FIRST bound phrase after the
  marker (within the same row for a table cell) must state every listed
  claim's published value:
    N ULP (N > 0)       "N ULP" / "N-ULP"
    0 ULP               "correctly rounded" or "0 ULP" / "0-ULP"
    abs53 1             "2^-53 absolute" / "absolute, of order 2^-53"
    abs53 0.5           "half-ulp absolute" / "half of 2^-53" / "0.500·2^-53 absolute"
  A marker that names an issue records a known wrong statement (corrected
  by that issue's row): reported, not failed. `claim: none` covers
  explanatory text that uses a bound phrase without stating a bound.
- Coverage: every claim id is marked in docs/ACCURACY.md, and every
  "max N ULP" phrase in README.md, docs/USER-GUIDE.md and corvus.h is
  covered by a marker.

  python tools/claims_check.py        (exit 1 on any failure; ctest claims_check)
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
TABLE = REPO / "tests" / "data" / "bounds.txt"
DOCS = ["README.md", "docs/USER-GUIDE.md", "docs/ACCURACY.md", "include/corvus/corvus.h"]
CONTRACT_DOCS = {"README.md", "docs/USER-GUIDE.md", "include/corvus/corvus.h"}
MARK = re.compile(r"<!--\s*claim:\s*([^>]*?)\s*-->")
WINDOW = 600

# A bound phrase: what it states, as (unit, value).
PHRASES = [
    (re.compile(r"half[- ]ulp absolute|half of 2\^-53", re.I), ("abs53", 0.5)),
    (re.compile(r"(\d*\.\d+)\s*·\s*2\^-53 absolute", re.I), "abs53"),  # value from the match
    (re.compile(r"2\^-53 absolute|absolute,? of order 2\^-53", re.I), ("abs53", 1.0)),
    (re.compile(r"correctly rounded", re.I), ("ulp", 0.0)),
    (re.compile(r"\b(\d+)[- ]ULP\b"), None),  # value from the match
]
MAX_ULP = re.compile(r"\bmax(?:imum)? (\d+) ULP\b", re.I)


def normalise(text: str) -> str:
    text = re.sub(r"^\s*///\s?", "", text, flags=re.M)
    text = text.replace("−", "-").replace("**", "").replace("`", "")
    return re.sub(r"\s+", " ", text)


def first_phrase(text: str):
    """(position, (unit, value)) of the first bound phrase, or None."""
    best = None
    for pat, val in PHRASES:
        m = pat.search(text)
        if m and (best is None or m.start() < best[0]):
            if val is None:
                v = ("ulp", float(m.group(1)))
            elif isinstance(val, str):
                v = (val, float(m.group(1)))
            else:
                v = val
            best = (m.start(), v, m.group(0))
    return best


def load_table():
    claims, errors, known = {}, [], []
    for n, line in enumerate(TABLE.read_text().splitlines(), 1):
        if not line or line.startswith("#"):
            continue
        f = line.split("|")
        if len(f) not in (6, 7):
            errors.append(f"bounds.txt:{n}: expected 6 or 7 fields")
            continue
        gate, region, unit, g, p, claim = f[:6]
        issue = f[6] if len(f) == 7 else ""
        g, p = float(g), float(p)
        if p < g:
            (known if issue else errors).append(
                f"bounds.txt:{n}: {gate}|{region} published {p:g} < gate {g:g}"
                + (f" (known: {issue})" if issue else ""))
        prev = claims.setdefault(claim, (unit, p))
        if prev != (unit, p):
            errors.append(f"bounds.txt:{n}: claim {claim} is {unit} {p:g} here, "
                          f"{prev[0]} {prev[1]:g} elsewhere")
    return claims, errors, known


def main():
    claims, errors, known = load_table()
    marked = {doc: set() for doc in DOCS}
    for doc in DOCS:
        raw = (REPO / doc).read_text(encoding="utf-8")
        marks = list(MARK.finditer(raw))
        covered = []  # (start, end) spans each marker governs, in raw offsets
        for i, m in enumerate(marks):
            end = marks[i + 1].start() if i + 1 < len(marks) else len(raw)
            end = min(end, m.end() + WINDOW)
            bol = raw.rfind("\n", 0, m.start()) + 1
            if raw[bol:m.start()].lstrip().startswith("|"):
                # A table cell: the claim is this row's, never the next one's.
                eol = raw.find("\n", m.end())
                end = min(end, eol if eol >= 0 else len(raw))
            covered.append((m.start(), end))
            line = raw.count("\n", 0, m.start()) + 1
            toks = m.group(1).split()
            issue = next((t for t in toks if t.startswith("#")), "")
            ids = [t for t in toks if not t.startswith("#")]
            window = normalise(raw[m.end():end])
            ph = first_phrase(window)
            if ids == ["none"]:
                continue  # explanatory text that states no bound
            for cid in ids:
                marked[doc].add(cid)
                if cid not in claims:
                    errors.append(f"{doc}:{line}: claim {cid} is not in bounds.txt")
                    continue
                want = claims[cid]
                if ph is None:
                    msg = f"{doc}:{line}: claim {cid}: no bound phrase follows the marker"
                    (known if issue else errors).append(msg + (f" (known: {issue})" if issue else ""))
                    continue
                got = ph[1]
                if (got[0], float(got[1])) != (want[0], float(want[1])):
                    msg = (f"{doc}:{line}: claim {cid}: the text says '{ph[2]}', "
                           f"the table publishes {want[1]:g} {want[0]}")
                    (known if issue else errors).append(msg + (f" (known: {issue})" if issue else ""))
        if doc in CONTRACT_DOCS:
            for m in MAX_ULP.finditer(raw):
                if not any(s <= m.start() < e for s, e in covered):
                    line = raw.count("\n", 0, m.start()) + 1
                    errors.append(f"{doc}:{line}: '{m.group(0)}' carries no claim marker")
    for cid in sorted(claims):
        if cid not in marked["docs/ACCURACY.md"]:
            errors.append(f"docs/ACCURACY.md: claim {cid} is never marked")
    for k in known:
        print(f"KNOWN {k}")
    for e in errors:
        print(f"FAIL {e}")
    n_marks = sum(len(v) for v in marked.values())
    print(f"{len(claims)} claims, {n_marks} marked claim ids across {len(DOCS)} docs: "
          f"{len(errors)} failures, {len(known)} known mismatches")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
