import re, os
D = os.path.dirname(os.path.abspath(__file__))
B = "gamma beta digamma trigamma gammainv betainv bessel lbeta erfinv".split()


def read(p, b):
    label, out = "", {}
    for line in open(f"{D}/pass{p}/quiet_bench_bench_{b}.txt", encoding="utf-8", errors="replace"):
        line = line.rstrip()
        m = re.search(r"([\d.]+)x\s*$", line)
        if not m:
            if line and not line.startswith("corvus") and "upper bnd" not in line and not line.lstrip().startswith("n "):
                label = line[:44]
            continue
        t = line.split()
        if line.startswith("lbeta"):
            out[line[:15]] = (0, float(m.group(1)), line[16:40])
        elif t[0].isdigit():
            n = int(t[0])
            if label not in out or n >= out[label][0]:
                out[label] = (n, float(m.group(1)), " ".join(t[1:-1]))
        else:
            out[label + "/" + t[0]] = (0, float(m.group(1)), " ".join(t[1:-1]))
    return out


for b in B:
    a, c = read(1, b), read(2, b)
    print("==", b, "range", min(v[1] for v in a.values()), max(v[1] for v in a.values()),
          "|", min(v[1] for v in c.values()), max(v[1] for v in c.values()))
    for k in a:
        y = c.get(k, (0, 0, ""))
        f = " <<" if abs(a[k][1] - y[1]) / a[k][1] > 0.1 else ""
        if f or b in ("bessel", "erfinv"):
            print(f"  {k[:46]:46s} {a[k][1]:7.2f} {y[1]:7.2f} {a[k][2][:26]}{f}")
