"""
Build team-10.csv in the format abguard expects.
Run from the folder containing participant-11/.
"""
import csv

PRE = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRY"
       "ADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSR")
POST = "WGQGTLVTVSS"
VL = ("DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPS"
      "RFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK")
WT10 = "WGGDGFYAMD"

# (loop, vsWT, minDist, note) — ranked best first
RANKED = [
    ("WEGSAYYAFEY", 6, 3, "AF2-IG iPTM 0.971 vs 0.975 WT; loop RMSD 2.31 A; framework Y kept"),
    ("IVGFANNHFDY", 8, 4, "best loop placement 1.86 A; framework Y kept"),
    ("WYGSAEYSLDS", 7, 3, "AF2-IG iPTM 0.971; highest Boltz-2 confidence; framework Y to S"),
    ("HYGWDHFYMDY", 7, 5, "largest novelty margin; framework Y kept"),
    ("EFGDEYHYFDY", 7, 4, "highest pTM of the set; framework Y kept"),
    ("THAESNYALDW", 7, 3, "near-baseline Boltz-2; framework Y to W"),
    ("ENGHSQFAMDH", 6, 4, "near-baseline Boltz-2; framework Y to H"),
    ("TYGAAYFAMGS", 7, 5, "high novelty margin; framework Y to S"),
    ("WFGDEYHYFDY", 6, 4, "framework Y kept"),
    ("WVGFANNHFDY", 7, 4, "framework Y kept"),
    ("MWGVARHAMDY", 6, 4, "framework Y kept"),
    ("DAGQGEYAMNY", 5, 4, "framework Y kept"),
    ("DVGPSLWALDY", 7, 3, "framework Y kept"),
    ("WVGPSLWALDY", 6, 2, "framework Y kept; minimum novelty margin"),
    ("EHGRAHYDMDN", 6, 4, "framework Y to N"),
    ("WLGGGNLAMDF", 4, 2, "W99 preserved; framework Y to F; minimum novelty margin"),
    ("TLGGGNLAMDF", 5, 3, "framework Y to F"),
    ("LWGNDFYALDF", 5, 3, "framework Y to F"),
    ("AYGSAEYSLDS", 7, 4, "framework Y to S"),
    ("ELGGAWYYLDF", 7, 3, "framework Y to F"),
]

rows = []
for i, (loop, vswt, mind, note) in enumerate(RANKED, 1):
    vh = PRE + loop + POST
    assert len(vh) <= 250 and len(VL) <= 250, f"too long: {loop}"
    assert vh.count("C") % 2 == 0, f"odd cysteines in VH: {loop}"
    assert VL.count("C") % 2 == 0, "odd cysteines in VL"
    assert set(vh) <= set("ACDEFGHIKLMNPQRSTVWY"), f"bad residue: {loop}"
    desc = (f"CDRH3 {loop}; {vswt}/10 positions changed from trastuzumab; "
            f"min {mind} mutations from all 525,042 provided sequences; {note}")
    rows.append({
        "sequence_id": f"t10-{i:02d}",
        "kind": "paired",
        "heavy_sequence": vh,
        "light_sequence": VL,
        "design_description": desc,
    })

with open("team-10.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["sequence_id", "kind", "heavy_sequence",
                                      "light_sequence", "design_description"])
    w.writeheader()
    w.writerows(rows)

print(f"wrote team-10.csv  ({len(rows)} designs)")
print(f"VH length {len(rows[0]['heavy_sequence'])}, VL length {len(VL)}")
print(f"cysteines  VH {rows[0]['heavy_sequence'].count('C')}, VL {VL.count('C')}")
for r in rows[:4]:
    print(f"  {r['sequence_id']}  {r['design_description'][:60]}...")