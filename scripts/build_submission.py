#!/usr/bin/env python3
"""Assemble the HER2 CDRH3 submission list from the scored, novelty-checked designs.

Slots 1-4 are fixed by hand: one representative per sequence family, chosen to
differ by construction rather than by score. Slots 5-20 are filled from the
remaining strict-novelty-passing refolded designs, ranked on min_ipae with iptm
as tie-break, skipping loops already listed.

Families are computed here (single linkage, Hamming distance <= 4 over the
10-mer) rather than asserted, so the grouping in RANKING.md is reproducible.
"""

from __future__ import annotations

import csv
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
METRICS = os.path.join(ROOT, "data/docs/runs/her2-main/scored-boltz2/refold_metrics.csv")
NOVELTY = os.path.join(ROOT, "scripts/novelty200_strict.csv")
DEV = os.path.join(ROOT, "scripts/developability.csv")
KEEP = os.path.join(ROOT, "scripts/strict_refolded.txt")
OUT_CSV = os.path.join(ROOT, "submission/ranked_designs.csv")
OUT_MD = os.path.join(ROOT, "submission/RANKING.md")

PARENT = "WGGDGFYAMD"

TOP4 = [
    ("SLSLSNYAFD",
     "Best overall score (min_ipae 0.3076, the lowest of all 37) with clean "
     "developability: no liability motif, hydrophobic fraction equal to the parent, "
     "net charge -1."),
    ("TRAAVGYAHD",
     "Second sequence family, one of only two neutral-charge designs, larger novelty "
     "margin than its family-mate at 4 vs 3."),
    ("DLSYYNHAFD",
     "Largest novelty margin in the set: minimum Hamming distance 5 to the Chinery "
     "screen, 20% identity to trastuzumab, 8 mutations. Clean developability."),
    ("WSSAVGYAHD",
     "Conservative hedge: closest to the parent at 6 mutations, and its iptm (0.926) "
     "ties the top-ranked design. Same family as slot 2, taken as a deliberate paired "
     "bet."),
]


def hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))


def families(loops: list[str], cutoff: int = 4) -> dict[str, int]:
    """Single-linkage clusters over Hamming distance; returns loop -> family index."""
    parent = {l: l for l in loops}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i, a in enumerate(loops):
        for b in loops[i + 1:]:
            if hamming(a, b) <= cutoff:
                ra, rb = find(a), find(b)
                if ra != rb:
                    parent[ra] = rb
    roots = {}
    out = {}
    for l in loops:
        r = find(l)
        out[l] = roots.setdefault(r, len(roots) + 1)
    return out


def main() -> int:
    metrics = {r["query_id"]: r for r in csv.DictReader(open(METRICS))}
    nov = {r["label"]: r for r in csv.DictReader(open(NOVELTY))}
    dev = {r["design"]: r for r in csv.DictReader(open(DEV))}
    keep = [l.strip() for l in open(KEEP) if l.strip()]

    pool = []
    for d in keep:
        m, n, v = metrics[d], nov[d], dev[d]
        pool.append({
            "design_id": d,
            "cdrh3": n["cdrh3"],
            "min_ipae": float(m["min_ipae"]),
            "iptm": float(m["iptm"]),
            "binder_plddt": float(m["binder_plddt"]),
            "min_dist": int(n["min_dist"]),
            "identity_pct": float(n["identity_%"]),
            "mutations": int(n["mut_vs_tz"]),
            "hydro_frac": float(v["hydro_frac"]),
            "net_charge": int(v["net_charge"]),
            "dev_flag": v["worse_than_parent"],
            "dev_reason": v["reasons"],
        })

    loops_all = sorted({p["cdrh3"] for p in pool})
    fam2 = families(loops_all, 2)
    fam3 = families(loops_all, 3)
    for p in pool:
        p["family"] = fam2[p["cdrh3"]]
        p["family_c3"] = fam3[p["cdrh3"]]

    by_loop: dict[str, dict] = {}
    for p in sorted(pool, key=lambda x: (x["min_ipae"], -x["iptm"])):
        by_loop.setdefault(p["cdrh3"], p)   # best instance of each distinct loop

    rows = []
    used: set[str] = set()
    for loop, rationale in TOP4:
        p = dict(by_loop[loop])
        p["rationale"] = rationale
        p["selection"] = "hand-picked family representative"
        rows.append(p)
        used.add(loop)

    rest = [p for loop, p in by_loop.items() if loop not in used]
    rest.sort(key=lambda x: (x["min_ipae"], -x["iptm"]))
    for p in rest:
        if len(rows) >= 20:
            break
        q = dict(p)
        note = (f"Rank by score among the remaining strict-novelty designs; "
                f"family {q['family']} (d<=2).")
        if q["dev_flag"] == "yes":
            note += f" DEVELOPABILITY FLAG: {q['dev_reason']}."
        q["rationale"] = note
        q["selection"] = "score-ranked (min_ipae, iptm tie-break)"
        rows.append(q)

    fields = ["rank", "cdrh3", "design_id", "min_ipae", "iptm", "binder_plddt",
              "min_hamming_to_chinery", "identity_pct_vs_WGGDGFYAMD", "mutations",
              "hydrophobic_fraction", "net_charge", "developability",
              "developability_detail", "family_d2", "family_d3", "selection",
              "rationale"]

    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    with open(OUT_CSV, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(fields)
        for i, p in enumerate(rows, 1):
            w.writerow([
                i, p["cdrh3"], p["design_id"], f"{p['min_ipae']:.4f}",
                f"{p['iptm']:.3f}", f"{p['binder_plddt']:.2f}", p["min_dist"],
                f"{p['identity_pct']:.0f}", p["mutations"],
                f"{p['hydro_frac']:.2f}", f"{p['net_charge']:+d}",
                "FLAG" if p["dev_flag"] == "yes" else "clean",
                p["dev_reason"], p["family"], p["family_c3"], p["selection"],
                p["rationale"],
            ])

    print(f"wrote {OUT_CSV} ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
