"""Summarise the six quiet passes: per bench, band and size, the median of 3
for A (pre-#42) and B (HEAD), the A/B spread, and A/B. Columns: forward
benches report simd ns/el; inverse benches report p_inv and q_inv."""
import glob
import os
import re
import statistics as st

D = os.path.dirname(os.path.abspath(__file__))
benches = ["bench_gamma", "bench_beta", "bench_gammainv", "bench_betainv"]
data = {}
order = []
for pdir in sorted(glob.glob(D + r"\pass*")):
    v = pdir[-1]
    for b in benches:
        label = None
        for line in open(rf"{pdir}\quiet_bench_{b}.txt", encoding="utf-8", errors="replace"):
            line = line.rstrip()
            if not line or line.startswith("corvus active") or re.match(r"\s*n\s+", line):
                continue
            m = re.match(r"\s*(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", line)
            if m and label:
                cols = [("p_inv", m.group(2)), ("q_inv", m.group(3))] if "inv" in b else [("simd", m.group(2))]
                for c, val in cols:
                    k = (b, label, m.group(1), c)
                    if k not in data:
                        order.append(k)
                    data.setdefault(k, {}).setdefault(v, []).append(float(val))
            else:
                label = re.split(r":\s*(n=|a in)", line)[0][:44]
print(f"{'bench':14s} {'band':44s} {'n':>7s} {'col':5s} {'A med':>9s} {'B med':>9s} {'A sprd':>6s} {'B sprd':>6s} {'A/B':>6s}")
for k in order:
    a, b = data[k]["A"], data[k]["B"]
    am, bm = st.median(a), st.median(b)
    sa, sb = (max(a) - min(a)) / am, (max(b) - min(b)) / bm
    print(f"{k[0]:14s} {k[1]:44s} {k[2]:>7s} {k[3]:5s} {am:9.1f} {bm:9.1f} {sa:6.1%} {sb:6.1%} {am / bm:5.2f}x")
