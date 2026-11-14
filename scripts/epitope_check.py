#!/usr/bin/env python3
"""Does the refolded de novo VHH sit on trastuzumab's HER2 epitope?

For each design, reads the Boltz-2 refolded complex (chains H = VHH, C and D =
HER2 domain IV fragments) and reports:
  * the minimum distance from ANY VHH atom to each of the four hotspot atoms
    D/593/NZ, C/573/CD1, D/591/CG2, C/572/CD,
  * the VHH footprint on HER2 (HER2 residues with a heavy atom within CONTACT A
    of the VHH),
  * overlap of that footprint with trastuzumab's own epitope, defined from the
    1N8Z coordinates (HER2 residues within CONTACT A of the trastuzumab Fv).

Distances are within one complex, so they do not depend on the coordinate frame
of the prediction. The Boltz-2 model renumbers HER2 residues 1..N within each
chain; hotspots are mapped by residue ORDER using the residue_ids that were
given to Boltz-2, and the script ASSERTS that the residue names at the mapped
positions match the original 1N8Z residues, so a wrong mapping fails loudly.

Caveat: the input contains only a domain-IV fragment of HER2 (489-607). "Drifted
elsewhere" can therefore only mean elsewhere ON that fragment.
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import math
import os
from collections import defaultdict

CONTACT = 4.5
HOT = [("D", 593, "NZ"), ("C", 573, "CD1"), ("D", 591, "CG2"), ("C", 572, "CD")]


def read_atoms(path: str):
    """Yield dicts with chain, seq (auth_seq_id), comp, atom, xyz, element."""
    cols, start = [], None
    lines = open(path).readlines()
    for i, l in enumerate(lines):
        s = l.strip()
        if s.startswith("_atom_site."):
            cols.append(s.split(".", 1)[1].split()[0])
            start = i + 1
        elif cols and s.startswith(("ATOM", "HETATM")):
            start = i
            break
    ix = {c: i for i, c in enumerate(cols)}
    out = []
    for l in lines[start:]:
        f = l.split()
        if not f or f[0] not in ("ATOM", "HETATM"):
            break
        x = (float(f[ix["Cartn_x"]]), float(f[ix["Cartn_y"]]), float(f[ix["Cartn_z"]]))
        if any(math.isnan(v) for v in x):
            continue
        out.append({
            "chain": f[ix["auth_asym_id"]],
            "seq": int(f[ix["auth_seq_id"]]),
            "comp": f[ix["auth_comp_id"]],
            "atom": f[ix["auth_atom_id"]],
            "el": f[ix["type_symbol"]],
            "xyz": x,
        })
    return out


def heavy(atoms):
    return [a for a in atoms if a["el"] != "H" and not a["atom"].startswith("H")]


def footprint(binder, target, cutoff):
    """Target residues (chain, orig_index_key) within cutoff of any binder atom."""
    hit = set()
    for t in target:
        for b in binder:
            dx = t["xyz"][0] - b["xyz"][0]
            if abs(dx) > cutoff:
                continue
            if math.dist(t["xyz"], b["xyz"]) <= cutoff:
                hit.add((t["chain"], t["seq"]))
                break
    return hit


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", default="data/docs/runs/denovo")
    ap.add_argument("--ref", default="data/docs/Data/12-her2/her2_fv_1n8z.cif")
    ap.add_argument("--designs", help="comma-separated design ids (default: all scored)")
    ap.add_argument("--top", type=int, default=5)
    ap.add_argument("--csv")
    args = ap.parse_args()

    # reference: original numbering, residue names, trastuzumab epitope
    ref = read_atoms(args.ref)
    ref_her2 = [a for a in ref if a["chain"] in ("C", "D")]
    ref_fv = heavy([a for a in ref if a["chain"] in ("A", "B")])
    ref_name = {(a["chain"], a["seq"]): a["comp"] for a in ref_her2}
    tz_epitope = footprint(ref_fv, heavy(ref_her2), CONTACT)

    # mapping model index -> original residue number, from the Boltz-2 query
    queries = json.load(open(os.path.join(args.run, "prepared/holo/queries.json")))["queries"]
    resids = {}
    for c in queries[0]["chains"]:
        if c["chain_id"] in ("C", "D"):
            resids[c["chain_id"]] = [int(r) for r in c["residue_ids"]]

    metrics = list(csv.DictReader(open(os.path.join(args.run, "scored-boltz2/refold_metrics.csv"))))
    key = "query_id" if "query_id" in metrics[0] else "design"
    ranked = sorted(metrics, key=lambda r: float(r["min_ipae"]))
    if args.designs:
        want = args.designs.split(",")
        ranked = [r for r in ranked if r[key] in want]
    else:
        ranked = ranked[:args.top]

    rows = []
    for r in ranked:
        d = r[key]
        path = glob.glob(os.path.join(args.run, f"boltz2/{d}/seed_*/{d}_seed_*_model.cif"))[0]
        atoms = read_atoms(path)
        vhh = heavy([a for a in atoms if a["chain"] == "H"])
        her2 = heavy([a for a in atoms if a["chain"] in ("C", "D")])

        # model index -> original residue number
        model_to_orig = {}
        for ch in ("C", "D"):
            idxs = sorted({a["seq"] for a in her2 if a["chain"] == ch})
            if len(idxs) != len(resids[ch]):
                raise SystemExit(f"{d}: chain {ch} has {len(idxs)} residues, "
                                 f"query listed {len(resids[ch])}")
            for i, orig in zip(idxs, resids[ch]):
                model_to_orig[(ch, i)] = orig
        # assert residue identities match the original structure
        for a in her2:
            orig = model_to_orig[(a["chain"], a["seq"])]
            if ref_name.get((a["chain"], orig)) != a["comp"]:
                raise SystemExit(f"{d}: residue mismatch at {a['chain']}{orig}: "
                                 f"model {a['comp']} vs 1N8Z {ref_name.get((a['chain'], orig))}")

        hot_d = {}
        for ch, orig, atom in HOT:
            tgt = [a for a in her2 if a["chain"] == ch
                   and model_to_orig[(ch, a["seq"])] == orig and a["atom"] == atom]
            if len(tgt) != 1:
                raise SystemExit(f"{d}: hotspot {ch}{orig}/{atom} found {len(tgt)} times")
            hot_d[f"{ch}/{orig}/{atom}"] = min(math.dist(tgt[0]["xyz"], v["xyz"]) for v in vhh)

        fp_model = footprint(vhh, her2, CONTACT)
        fp = {(ch, model_to_orig[(ch, i)]) for ch, i in fp_model}
        overlap = fp & tz_epitope
        rows.append({
            "design": d,
            "min_ipae": float(r["min_ipae"]), "iptm": float(r["iptm"]),
            **{f"d_{k}": v for k, v in hot_d.items()},
            "hotspots_within_4.5A": sum(1 for v in hot_d.values() if v <= CONTACT),
            "her2_residues_contacted": len(fp),
            "tz_epitope_residues": len(tz_epitope),
            "contacted_in_tz_epitope": len(overlap),
            "frac_footprint_in_tz_epitope": (len(overlap) / len(fp)) if fp else 0.0,
            "frac_tz_epitope_covered": len(overlap) / len(tz_epitope),
        })

    print(f"# contact cutoff {CONTACT} A; trastuzumab epitope in 1N8Z = "
          f"{len(tz_epitope)} HER2 residues")
    hk = [f"{c}/{o}/{a}" for c, o, a in HOT]
    hdr = (f"{'design':26}{'ipae':>6}{'iptm':>6}" + "".join(f"{k:>12}" for k in hk)
           + f"{'hot<=4.5':>9}{'contact':>8}{'inTz':>6}{'%fp':>6}{'%tz':>6}")
    print(hdr)
    print("-" * len(hdr))
    for x in rows:
        print(f"{x['design']:26}{x['min_ipae']:6.2f}{x['iptm']:6.2f}"
              + "".join(f"{x['d_' + k]:12.1f}" for k in hk)
              + f"{x['hotspots_within_4.5A']:9d}{x['her2_residues_contacted']:8d}"
              + f"{x['contacted_in_tz_epitope']:6d}"
              + f"{100 * x['frac_footprint_in_tz_epitope']:6.0f}"
              + f"{100 * x['frac_tz_epitope_covered']:6.0f}")
    print("# hot<=4.5: hotspot atoms with a VHH atom within 4.5 A (of 4); contact: HER2 "
          "residues touched; inTz: of those, in trastuzumab's epitope; %fp: share of the "
          "VHH footprint inside the trastuzumab epitope; %tz: share of the trastuzumab "
          "epitope covered")

    if args.csv:
        with open(args.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            for x in rows:
                w.writerow({k: (f"{v:.3f}" if isinstance(v, float) else v) for k, v in x.items()})
        print(f"# wrote {args.csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
