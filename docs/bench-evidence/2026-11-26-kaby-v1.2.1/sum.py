"""Median of the three fresh processes per tree, and branch/main, for
bench_span1's two tables. Usage: python3 sum.py (from this directory)."""
import re, statistics as st
def read(name):
    rows = {}
    for rep in (1, 2, 3):
        sec = None
        for line in open(f"{name}{rep}/quiet_bench_bench_span1.txt"):
            if line.startswith("ns per value"): sec = "value"; continue
            if line.startswith("ns per CALL"): sec = "call"; continue
            m = re.match(r"(.{28})\s+([\d.\s]+)$", line.rstrip())
            if sec and m:
                rows.setdefault((sec, m.group(1).strip()), []).append([float(x) for x in m.group(2).split()])
    return {k: [st.median(c) for c in zip(*v)] for k, v in rows.items()}, \
           {k: [max(c) / min(c) for c in zip(*v)] for k, v in rows.items()}
(m, ms), (b, bs) = read("main"), read("branch")
for sec, cols in (("value", "batch array local"), ("call", "n=1 n=2 n=3 n=4 n=7")):
    print(f"\n{sec}: median of 3, main -> branch (ratio); columns {cols}")
    for k in m:
        if k[0] != sec: continue
        print(f"{k[1]:26s}" + "".join(f"  {x:6.1f}->{y:6.1f} ({y/x:4.2f})" for x, y in zip(m[k], b[k])))
print(f"\nmax/min across the three processes, worst cell: main {max(max(v) for v in ms.values()):.2f}, branch {max(max(v) for v in bs.values()):.2f}")
