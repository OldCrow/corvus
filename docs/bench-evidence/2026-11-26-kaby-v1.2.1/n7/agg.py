"""Aggregate bench_tail runs. usage: python3 agg.py <dir> <stat> <name>...
stat: med = median over processes of each process's median;
      low = minimum over processes of each process's minimum."""
import sys, glob, statistics as st, collections
d, stat, names = sys.argv[1], sys.argv[2], sys.argv[3:]
v = collections.defaultdict(list)
for nm in names:
    for f in glob.glob(f"{d}/{nm}_*.txt"):
        for l in open(f):
            p = l.split()
            if len(p) == 5 and p[0] != '#':
                v[(nm, p[0], p[1], int(p[2]))].append((float(p[3]), float(p[4])))
keys = sorted({k[1:] for k in v}, key=lambda k: (["erfinv","exp","erf","erfcinv"].index(k[0]), ["wr+rd","wr","rd","none"].index(k[1]), k[2]))
print(f"{'fn':8s}{'shape':6s}{'n':>3s} " + " ".join(f"{n:>7s}" for n in names))
for k in keys:
    row = []
    for nm in names:
        x = v[(nm,)+k]
        row.append(st.median(a for a, _ in x) if stat == "med" else min(b for _, b in x))
    print(f"{k[0]:8s}{k[1]:6s}{k[2]:3d} " + " ".join(f"{r:7.1f}" for r in row))
