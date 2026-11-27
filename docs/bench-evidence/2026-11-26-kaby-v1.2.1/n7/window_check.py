"""For an executable and an AVX2 driver symbol, list the constant loads on the
one-value path that come AFTER the first gather and touch page offsets
0x7c0..0x7cf (the window gather_alias2.cpp measures). usage:
python3 window_check.py <symbol-regex> <binary>...   or, with --sweep <file>,
predict literal_sweep.txt's slow rows from one binary linked at E = 0."""
import re, subprocess, sys
OBJ = "/usr/local/opt/llvm/bin/llvm-objdump"
def path(binary, pat):
    nm = subprocess.run(["nm", binary], capture_output=True, text=True).stdout
    sym = [l.split()[2] for l in nm.splitlines() if re.search(pat, l)][0]
    out = subprocess.run([OBJ, "-d", "--no-show-raw-insn", f"--disassemble-symbols={sym}", binary],
                         capture_output=True, text=True).stdout.splitlines()
    res, started = [], False
    for l in out:
        m = re.match(r"\s*([0-9a-f]+):\s+(\S+)\s*(.*)", l)
        if not m: continue
        mn, ops = m.group(2), m.group(3)
        if mn == "jne" and not started: started = True; continue
        if not started: continue
        a = re.search(r"##\s*0x([0-9a-f]+)", ops)
        res.append((mn, int(a.group(1), 16) if a and "(%rip)" in ops and mn != "leaq" else None))
        if mn in ("jmp", "retq"): break
    return res
def hits(p, shift=0):
    g = [i for i, (mn, _) in enumerate(p) if "gather" in mn]
    if not g: return []
    out = []
    for i, (mn, a) in enumerate(p):
        if a is None or i < g[0]: continue
        w = 32 if mn in ("vmovdqa", "vmovapd", "vmovaps", "vpcmpgtq", "vpand", "vandpd") else 8
        lo = (a - shift) & 0xfff
        if lo < 0x7d0 and lo + w > 0x7c0: out.append((i, mn, lo))
    return out
if sys.argv[1] == "--sweep":
    p = path(sys.argv[3], "N_AVX2L7ErfImplE"); ok = n = 0
    for l in open(sys.argv[2]):
        f = l.split()
        if not f or f[0].startswith("#"): continue
        E, ns = int(f[0]), float(f[5]); pred = bool(hits(p, E)); slow = ns > 40
        n += 1; ok += pred == slow
        if pred or slow: print(f"E={E:5d} measured {ns:6.2f} ns  predicted {'SLOW' if pred else 'fast'}  {[(i, hex(a)) for i, _, a in hits(p, E)]}")
    print(f"prediction matches measurement at {ok} of {n} positions")
else:
    for b in sys.argv[2:]:
        h = hits(path(b, sys.argv[1]))
        print(f"{b}: {'IN WINDOW ' + str([(i, mn, hex(a)) for i, mn, a in h]) if h else 'no post-gather constant in 0x7c0..0x7cf'}")
