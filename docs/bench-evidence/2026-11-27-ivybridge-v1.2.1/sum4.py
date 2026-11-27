"""Pass 3: median (min-max) over the fresh processes of each of the four
trees, for bench_span1's two tables. Usage: python3 sum4.py (from this
directory)."""
import glob, re, statistics as st
SIDES = ("main", "branch", "masked", "lane0")
def read(name):
    rows = {}
    for f in sorted(glob.glob(f"pass3/{name}_*.txt")):
        sec = None
        for line in open(f):
            if line.startswith("ns per value"): sec = "value"; continue
            if line.startswith("ns per CALL"): sec = "call"; continue
            m = re.match(r"(.{28})\s+([\d.\s]+)$", line.rstrip())
            if sec and m:
                rows.setdefault((sec, m.group(1).strip()), []).append([float(x) for x in m.group(2).split()])
    return rows
R = {s: read(s) for s in SIDES}
n = len(next(iter(R["main"].values())))
for sec, cols in (("value", ("batch", "array", "local")), ("call", ("n=1", "n=2", "n=3", "n=4", "n=7"))):
    print(f"\n{sec}: median of {n} (min-max), " + " / ".join(SIDES))
    for k in R["main"]:
        if k[0] != sec: continue
        for c, col in enumerate(cols):
            cells = []
            for s in SIDES:
                v = [r[c] for r in R[s][k]]
                cells.append(f"{st.median(v):6.1f} ({min(v):.1f}-{max(v):.1f})")
            print(f"{k[1]:26s} {col:5s} " + " | ".join(cells))
