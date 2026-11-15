import numpy as np, pandas as pd

RT = 0.593  # kcal/mol at 298 K

df = pd.read_csv("her2-data/trastuzumab_her2_spr_all_768.csv")
df = df[df.kd_nM.notna()].copy()
df["dG"] = RT * np.log(df.kd_nM * 1e-9)

wt = df[df.cdrh3 == "WGGDGFYAMD"]
if len(wt):
    print(f"trastuzumab: KD {wt.kd_nM.iloc[0]} nM -> dG {wt.dG.iloc[0]:.2f} kcal/mol\n")

print(f"n = {len(df)}")
print(f"KD  {df.kd_nM.min():.3f} to {df.kd_nM.max():.0f} nM, median {df.kd_nM.median():.1f}")
print(f"dG  {df.dG.min():.2f} to {df.dG.max():.2f} kcal/mol, spread {df.dG.max()-df.dG.min():.1f}\n")

print("ten tightest binders:")
print(df.nsmallest(10, "dG")[["cdrh3", "kd_nM", "dG"]].to_string(index=False))

df.to_csv("her2-data/spr_with_dG.csv", index=False)
print("\nwrote her2-data/spr_with_dG.csv")