"""Median over the fresh processes per tree, and branch/main, for
bench_span1's two tables. Usage: python3 sum.py pass1|pass2 (from this
directory)."""
import glob, re, statistics as st, sys
d = sys.argv[1]
def read(name):
    rows = {}
    files = sorted(glob.glob(f"{d}/{name}_*.txt"))
    for f in files:
        sec = None
        for line in open(f):
            if line.startswith("ns per value"): sec = "value"; continue
            if line.startswith("ns per CALL"): sec = "call"; continue
            m = re.match(r"(.{28})\s+([\d.\s]+)$", line.rstrip())
            if sec and m:
                rows.setdefault((sec, m.group(1).strip()), []).append([float(x) for x in m.group(2).split()])
    return {k: [st.median(c) for c in zip(*v)] for k, v in rows.items()}, \
           {k: [max(c) / min(c) for c in zip(*v)] for k, v in rows.items()}, len(files)
(m, ms, n), (b, bs, _) = read("main"), read("branch")
for sec, cols in (("value", "batch array local"), ("call", "n=1 n=2 n=3 n=4 n=7")):
    print(f"\n{d} {sec}: median of {n}, main -> branch (ratio); columns {cols}")
    for k in m:
        if k[0] != sec: continue
        print(f"{k[1]:26s}" + "".join(f"  {x:6.1f}->{y:6.1f} ({y/x:4.2f})" for x, y in zip(m[k], b[k])))
print("\nmax/min across the processes, per cell (main | branch):")
for k in m:
    print(f"{k[0]:5s} {k[1]:26s}" + " ".join(f"{x:5.2f}" for x in ms[k]) + "  | " + " ".join(f"{x:5.2f}" for x in bs[k]))
