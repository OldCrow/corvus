"""Aggregate bench_tail2 runs: lowest level per cell (min over processes of
each process's min). usage: python3 agg2.py <dir> <fn> <name>..."""
import sys, glob, collections
d, fn, names = sys.argv[1], sys.argv[2], sys.argv[3:]
v = collections.defaultdict(list)
for nm in names:
    for f in glob.glob(f"{d}/{nm}_*.txt"):
        for l in open(f):
            p = l.split()
            if len(p) == 6 and p[0] == fn:
                v[(nm, p[1], int(p[3]), int(p[2]))].append(float(p[5]))
for store in ("s8", "mc"):
    for n in (1, 2, 3, 4, 5, 7, 8):
        print(f"\n{fn} store={store} n={n}   off: " + " ".join(f"{o:6d}" for o in range(8)))
        for nm in names:
            print(f"  {nm:6s}             " + " ".join(f"{min(v[(nm, store, n, o)]):6.1f}" for o in range(8)))
