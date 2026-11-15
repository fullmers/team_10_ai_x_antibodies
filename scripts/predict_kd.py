import pandas as pd, numpy as np
from sklearn.ensemble import RandomForestRegressor

spr = pd.read_csv("her2-data/spr_with_dG.csv")
AA = "ACDEFGHIKLMNPQRSTVWY"

def enc(s):
    v = np.zeros(10 * 20)
    for i, c in enumerate(s[:10]):
        if c in AA:
            v[i * 20 + AA.index(c)] = 1
    return v

X = np.array([enc(s) for s in spr.cdrh3])
y = spr.dG.values
m = RandomForestRegressor(n_estimators=400, random_state=0).fit(X, y)

print(f"trained on {len(X)} measured designs\n")
print(f"{'id':5s} {'cdrh3':12s} {'dG':>7s} {'KD (nM)':>12s}")

for name, loop in [("WT", "WGGDGFYAMD"),
                   ("d12", "WEGSAYYAFE"),
                   ("d04", "IVGFANNHFD"),
                   ("d15", "WYGSAEYSLD"),
                   ("d09", "HYGWDHFYMD")]:
    dg = m.predict([enc(loop)])[0]
    kd = np.exp(dg / 0.593) * 1e9
    print(f"{name:5s} {loop:12s} {dg:7.2f} {kd:12.2f}")

print("\nWT measured: KD 0.152 nM, dG -13.41")
print("Compare the WT row above — the gap is the model's own error.")