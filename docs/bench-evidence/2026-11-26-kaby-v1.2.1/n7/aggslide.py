"""Aggregate bench_slide runs: per fn and n, the lowest reading at each
slide K, as a compact row. usage: python3 aggslide.py <dir> <name>..."""
import sys, glob, collections
d, names = sys.argv[1], sys.argv[2:]
for nm in names:
    v = collections.defaultdict(list); addr = set()
    for f in glob.glob(f"{d}/{nm}_*.txt"):
        for l in open(f):
            p = l.split()
            if l.startswith("# addr"): addr.add(p[2])
            elif len(p) == 5 and p[0] != '#': v[(p[0], int(p[1]), int(p[2]))].append(float(p[4]))
    print(f"\n== {nm}  alloca block at K=0: {sorted(addr)}")
    ks = sorted({k[2] for k in v})
    for fn in ("erf", "erfcinv", "erfinv", "exp"):
        for n in (1, 4, 7):
            row = [min(v[(fn, n, K)]) for K in ks]
            base = sorted(row)[len(row)//4]
            print(f"{fn:8s} n={n} base {base:6.1f} max {max(row):6.1f}  " + "".join(
                " ." if x < base*1.03 else f"{round((x/base-1)*100):2d}" for x in row))
print("\none column per slide step; '.' = within 3% of the lower-quartile reading, else % above it")
