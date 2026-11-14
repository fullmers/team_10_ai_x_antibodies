import json, glob, os

rows = []
for d in sorted(glob.glob("scores/*/")):
    name = os.path.basename(d.rstrip("/"))
    big = [f for f in glob.glob(d + "*.json")
           if os.path.basename(f) != "result.json" and os.path.getsize(f) > 100000]
    if not big:
        continue
    j = json.load(open(big[0]))
    conf = j["confidence_scores"][0]
    ptm = j["ptm_scores"][0]
    rows.append((name, conf, ptm, j["structures"][0]["structure"]))

wt = next((r for r in rows if "wt_" in r[0]), None)
if wt:
    print(f"BASELINE  conf={wt[1]:.4f}  ptm={wt[2]:.4f}\n")

os.makedirs("structures", exist_ok=True)
print(f"{'design':34s} {'conf':>7s} {'dConf':>7s} {'ptm':>7s} {'dPtm':>7s}")
for name, conf, ptm, struct in sorted(rows, key=lambda r: -r[1]):
    dc = conf - wt[1] if wt else 0
    dp = ptm - wt[2] if wt else 0
    mark = " <<<" if wt and conf > wt[1] else ""
    print(f"{name:34s} {conf:7.4f} {dc:+7.4f} {ptm:7.4f} {dp:+7.4f}{mark}")
    open(f"structures/{name}.cif", "w").write(struct)

print(f"\nwrote {len(rows)} structures to structures/")