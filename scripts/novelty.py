import pandas as pd, numpy as np, re

facs = pd.read_csv("her2-data/tz_her2_facs_3classes_524346.csv")
spr  = pd.read_csv("her2-data/trastuzumab_her2_spr_all_768.csv")

provided = set(facs.cdrh3.dropna()) | set(spr.cdrh3.dropna())
provided = {s for s in provided if len(s) == 10}
print(f"provided unique 10-mers: {len(provided)}")

arr = np.array([list(s) for s in sorted(provided)])
WT10 = "WGGDGFYAMD"

designs = []
for blk in open("participant-11/designs/designs.fasta").read().split(">")[1:]:
    ln = blk.strip().split("\n")
    name, seq = ln[0], "".join(ln[1:])
    m = re.search(r"YYCSR(.*?)WG[QKR]G", seq)
    if not m: continue
    loop = m.group(1)
    designs.append((name, loop, loop[:10], loop[10:]))

print(f"\n{'design':28s} {'loop10':12s} {'tail':5s} {'vsWT':>5s} {'minDist':>8s} {'verdict'}")
for name, loop, l10, tail in designs:
    q = np.array(list(l10))
    mind = int((arr != q).sum(axis=1).min())
    wt = sum(a != b for a, b in zip(l10, WT10))
    ok = (mind >= 2) and (wt >= 3)
    flag = "" if tail == "Y" else " TAIL!"
    print(f"{name:28s} {l10:12s} {tail:5s} {wt:5d} {mind:8d}  "
          f"{'PASS' if ok else 'FAIL'}{flag}")