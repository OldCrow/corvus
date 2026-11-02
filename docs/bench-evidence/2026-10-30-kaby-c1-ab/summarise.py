"""Medians, spreads and A/B per band over quiet A/B passes.

usage: summarise.py [evidence dir]; each pass is a folder pass<i>-<A|B> holding
quiet_bench_<bench>.txt. A missing bench file is skipped. Prints n = 10^6 rows
for bench_lgamma and every row for the others."""
import glob, os, re, statistics as st, sys

D = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
data, order = {}, []


def add(k, v, x):
    if k not in data:
        order.append(k)
    data.setdefault(k, {}).setdefault(v, []).append(x)


for pdir in sorted(glob.glob(os.path.join(D, "pass*"))):
    v = pdir[-1]
    for b in ["bench_lgamma", "bench_gamma", "bench_beta", "bench_lbeta", "bench_gammainv", "bench_betainv"]:
        label = None
        path = os.path.join(pdir, f"quiet_bench_{b}.txt")
        if not os.path.exists(path):
            continue
        for line in open(path, encoding="utf-8", errors="replace"):
            line = line.rstrip()
            if not line or line.startswith("corvus active") or re.match(r"\s*n\s+", line):
                continue
            m = re.match(r"lbeta (.+?)\s*:\s+([\d.]+) ns/el", line)
            if m:
                add((b, m.group(1), "-", "simd"), v, float(m.group(2)))
                continue
            m = re.match(r"\s*(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", line)
            if m and label:
                n = m.group(1)
                if b == "bench_lgamma":
                    add((b, label, n, "corvus"), v, float(m.group(2)))
                    add((b, label, n, "libm"), v, float(m.group(3)))
                elif "inv" in b:
                    add((b, label, n, "p_inv"), v, float(m.group(2)))
                    add((b, label, n, "q_inv"), v, float(m.group(3)))
                else:
                    add((b, label, n, "simd"), v, float(m.group(2)))
            else:
                label = re.split(r":\s*(n=|a in)", line)[0][:40].split(" (")[0]
for k in order:
    a, b = data[k]["A"], data[k]["B"]
    am, bm = st.median(a), st.median(b)
    sa, sb = (max(a) - min(a)) / am, (max(b) - min(b)) / bm
    if k[2] in ("10000",):
        continue
    print(f"{k[0]:14s} {k[1]:40s} {k[2]:>7s} {k[3]:6s} {am:8.2f} {bm:8.2f} {sa:6.1%} {sb:6.1%} {am / bm:5.2f}x")

