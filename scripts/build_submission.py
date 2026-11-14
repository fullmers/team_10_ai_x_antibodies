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


DENOVO_RUN = os.path.join(ROOT, "data/docs/runs/denovo")
DENOVO_PICKS = ["her2_vhh_denovo_rep117", "her2_vhh_denovo_rep106"]  # by min_ipae
NOVELTY_SCOPE_REDESIGN = ("Chinery FACS+SPR screen (525,042 length-10 CDRH3), strict rule; "
                          "trastuzumab CDRH3 variants only")
NOVELTY_SCOPE_DENOVO = ("SAbDab2 ML subset ONLY (15,641 structures, 5,301 unique antibodies; "
                        "NOT full SAbDab); CDR LCS identity <70%; chance-level baseline")


def load_denovo() -> list[dict]:
    """The two de novo VHHs: refold scores, epitope footprint, SAbDab2-subset identity."""
    import json
    met = {r["query_id"]: r for r in csv.DictReader(
        open(os.path.join(DENOVO_RUN, "scored-boltz2/refold_metrics.csv")))}
    epi = {r["design"]: r for r in csv.DictReader(
        open(os.path.join(ROOT, "scripts/epitope_check_top5.csv")))}
    nov_lines = [l for l in open(os.path.join(ROOT, "scripts/novelty_sabdab2.csv"))
                 if not l.startswith('"#')]
    nov = {r["design"]: r for r in csv.DictReader(nov_lines)}
    q = json.load(open(os.path.join(DENOVO_RUN, "prepared/holo/queries.json")))["queries"]
    seqs = {x["query_id"]: next(c["sequence"] for c in x["chains"] if c["chain_id"] == "H")
            for x in q}
    out = []
    for d in DENOVO_PICKS:
        m, e, n = met[d], epi[d], nov[d]
        out.append({
            "design_id": d, "cdr3": n["cdr3"], "seq": seqs[d],
            "min_ipae": float(m["min_ipae"]), "iptm": float(m["iptm"]),
            "binder_plddt": float(m["binder_plddt"]),
            "hot": f"{e['hotspots_within_4.5A']}/4",
            "fp": f"{e['her2_residues_contacted']} (of which {e['contacted_in_tz_epitope']} in trastuzumab epitope)",
            "fp_pct": 100 * float(e["frac_footprint_in_tz_epitope"]),
            "id3": float(n["id_cdr3"]), "id123": float(n["id_cdr123"]),
            "nn3": f"{n['nn_cdr3']} ({n['nn_cdr3_type']}, {n['nn_cdr3_pdb'][-4:].upper()})",
            "rationale": (
                "De novo VHH on the caplacizumab framework against HER2 domain IV, hotspots "
                "from the trastuzumab epitope. In the refolded Boltz-2 pose it sits on that "
                f"epitope ({e['hotspots_within_4.5A']}/4 hotspot atoms within 4.5 A, "
                f"{100 * float(e['frac_footprint_in_tz_epitope']):.0f}% of its footprint inside "
                "trastuzumab's). Passes <70% CDR identity against the SAbDab2 ML subset only; "
                "the nearest-neighbour identity is at chance level (shuffled-loop control median "
                "0.61). `pass` is empty (missing pose-RMSD evidence); scores are single-sample "
                "and not comparable with the CDRH3 track."),
        })
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

    fields = ["rank", "track", "cdrh3", "design_id", "min_ipae", "iptm", "binder_plddt",
              "min_hamming_to_chinery", "identity_pct_vs_WGGDGFYAMD", "mutations",
              "hydrophobic_fraction", "net_charge", "developability",
              "developability_detail", "family_d2", "family_d3",
              "epitope_hotspots_within_4.5A", "epitope_footprint_residues",
              "epitope_footprint_in_trastuzumab_pct", "sabdab2_cdr3_identity",
              "sabdab2_cdr123_identity", "sabdab2_nearest_cdr3", "novelty_scope",
              "selection", "rationale", "vhh_sequence"]

    denovo = load_denovo()

    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    with open(OUT_CSV, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(fields)
        for i, p in enumerate(rows, 1):
            w.writerow([
                i, "cdrh3_redesign", p["cdrh3"], p["design_id"], f"{p['min_ipae']:.4f}",
                f"{p['iptm']:.3f}", f"{p['binder_plddt']:.2f}", p["min_dist"],
                f"{p['identity_pct']:.0f}", p["mutations"],
                f"{p['hydro_frac']:.2f}", f"{p['net_charge']:+d}",
                "FLAG" if p["dev_flag"] == "yes" else "clean",
                p["dev_reason"], p["family"], p["family_c3"],
                "", "", "", "", "", "", NOVELTY_SCOPE_REDESIGN,
                p["selection"], p["rationale"], "",
            ])
        for j, d in enumerate(denovo, len(rows) + 1):
            w.writerow([
                j, "de_novo", d["cdr3"], d["design_id"], f"{d['min_ipae']:.4f}",
                f"{d['iptm']:.3f}", f"{d['binder_plddt']:.2f}", "n/a", "n/a", "n/a",
                "not assessed", "not assessed", "not assessed", "", "", "",
                d["hot"], d["fp"], f"{d['fp_pct']:.0f}", f"{d['id3']:.3f}",
                f"{d['id123']:.3f}", d["nn3"], NOVELTY_SCOPE_DENOVO,
                "de novo VHH, appended after the CDRH3 list; ordered by min_ipae",
                d["rationale"], d["seq"],
            ])

    print(f"wrote {OUT_CSV} ({len(rows) + len(denovo)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
