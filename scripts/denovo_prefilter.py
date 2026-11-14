#!/usr/bin/env python3
"""Submission pre-filters for the de novo VHHs (chain H sequences).

1. Extract chain H from the *_generated.cif in design/generation/inference_samples/
   and cross-check it against the sequence Boltz-2 was given (prepared/holo/queries.json).
2. Hard filters: even cysteine count, total length < 250, canonical residues only.
3. Developability motif screen (same definitions as scripts/developability.py)
   run on the WHOLE chain, with every hit labelled as lying in a designed region
   or in the fixed caplacizumab framework.  The tools were written for a
   10-residue CDRH3 and flank it with chain-B context, so they are applied here
   by position instead of by loop.  The "parent" comparison used for the CDRH3
   track does not exist for a VHH (parent is a 10-mer in a different framework);
   the fixed-framework hits, identical in both designs, are the baseline.
4. Stability proxy counts (scripts/stability_proxy.py definitions): Gly, Pro,
   longest G/S run.  PROXY ONLY, not a Tm prediction.

"Designed region" = the three placeholder segments of the framework (13, 9-10
and 12-21 residues).  They are wider than the IMGT CDRs: each carries a few
flanking framework positions that the model also redesigns.
"""

from __future__ import annotations

import argparse
import glob
import importlib.util
import json
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CANON = set("ACDEFGHIKLMNPQRSTVWY")
THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
         "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
         "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
FRAME = re.compile(
    r"^(?P<f1>EVQLVESGGGLVQPGGSLRLSC)(?P<s1>.{13})(?P<f2>WFRQAPGKGRELVA)(?P<s2>.{9,10})"
    r"(?P<f3>YPDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYC)(?P<s3>.{12,21})(?P<f4>WGQGTQVTVS)$")


def load(name: str):
    spec = importlib.util.spec_from_file_location(name, os.path.join(ROOT, "scripts", name + ".py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def chain_h_from_cif(path: str) -> str:
    cols, start, lines = [], None, open(path).readlines()
    for i, l in enumerate(lines):
        s = l.strip()
        if s.startswith("_atom_site."):
            cols.append(s.split(".", 1)[1].split()[0])
            start = i + 1
        elif cols and s.startswith(("ATOM", "HETATM")):
            start = i
            break
    ix = {c: i for i, c in enumerate(cols)}
    res = {}
    for l in lines[start:]:
        f = l.split()
        if not f or f[0] not in ("ATOM", "HETATM"):
            break
        if f[ix["auth_asym_id"]] == "H":
            res.setdefault(int(f[ix["auth_seq_id"]]), f[ix["auth_comp_id"]])
    nums = sorted(res)
    if nums != list(range(nums[0], nums[-1] + 1)):
        raise SystemExit(f"{path}: chain H numbering has gaps")
    return "".join(THREE.get(res[n], "X") for n in nums)


def hits(seq: str, dev) -> list[tuple[str, int, int]]:
    """(motif, start, end) 0-based half-open spans, using developability.py's definitions."""
    out = []
    for m in ("NG", "DG", "DP"):
        i = seq.find(m)
        while i >= 0:
            out.append((m, i, i + 2))
            i = seq.find(m, i + 1)
    out += [("N-x-S/T", x.start(), x.start() + 3) for x in re.finditer(r"(?=N[^P][ST])", seq)]
    out += [("C", i, i + 1) for i, ch in enumerate(seq) if ch == "C"]
    hyd = "[" + "".join(sorted(dev.HYDROPHOBIC)) + "]+"
    out += [(f"hydrophobic run {len(x.group())}", x.start(), x.end())
            for x in re.finditer(hyd, seq) if len(x.group()) >= dev.HYDRO_RUN_MIN]
    out += [(f"K/R run {len(x.group())}", x.start(), x.end())
            for x in re.finditer(r"[KR]+", seq) if len(x.group()) >= dev.KR_RUN_MIN]
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", default=os.path.join(ROOT, "data/docs/runs/denovo"))
    ap.add_argument("--designs", default="her2_vhh_denovo_rep117,her2_vhh_denovo_rep106")
    ap.add_argument("--fasta", default=os.path.join(ROOT, "submission/denovo_vhh_sequences.fasta"))
    ap.add_argument("--csv", default=os.path.join(ROOT, "scripts/denovo_prefilter.csv"))
    args = ap.parse_args()

    dev = load("developability")
    stab = load("stability_proxy")
    q = json.load(open(os.path.join(args.run, "prepared/holo/queries.json")))["queries"]
    boltz = {x["query_id"]: next(c["sequence"] for c in x["chains"] if c["chain_id"] == "H") for x in q}

    records, rows = [], []
    for d in args.designs.split(","):
        cif = glob.glob(os.path.join(args.run, f"design/generation/inference_samples/{d}/{d}_generated.cif"))
        if len(cif) != 1:
            raise SystemExit(f"{d}: generated cif not found")
        seq = chain_h_from_cif(cif[0])
        same = seq == boltz[d]
        m = FRAME.match(seq)
        if not m:
            raise SystemExit(f"{d}: does not match the caplacizumab frame pattern")
        # designed positions (0-based)
        regions, pos = {}, 0
        for g in ("f1", "s1", "f2", "s2", "f3", "s3", "f4"):
            n = len(m.group(g))
            if g.startswith("s"):
                regions[g] = (pos, pos + n)
            pos += n
        designed = {i for a, b in regions.values() for i in range(a, b)}
        dseq = "".join(seq[i] for i in sorted(designed))

        n_c = seq.count("C")
        fixed_c = sum(m.group(g).count("C") for g in ("f1", "f2", "f3", "f4"))
        c_pos = [i + 1 for i, ch in enumerate(seq) if ch == "C"]
        hard = {
            "cys_total": n_c, "cys_even": n_c % 2 == 0, "cys_in_designed_regions": n_c - fixed_c,
            "length": len(seq), "length_ok": len(seq) < 250,
            "canonical_only": set(seq) <= CANON,
            "chain_H_matches_boltz2_query": same,
        }

        found = hits(seq, dev)
        lab = []
        for name, a, b in found:
            if name == "C":
                continue
            where = "designed region" if any(i in designed for i in range(a, b)) else "fixed framework"
            lab.append((name, a + 1, b, where, seq[a:b]))

        def comp(s):
            return {"hydro_frac": sum(1 for c in s if c in dev.HYDROPHOBIC) / len(s),
                    "net_charge": sum(1 for c in s if c in dev.POSITIVE) - sum(1 for c in s if c in dev.NEGATIVE),
                    "n_his": s.count("H")}

        per_region = {}
        for g, nm in (("s1", "CDR1 region"), ("s2", "CDR2 region"), ("s3", "CDR3 region")):
            a, b = regions[g]
            s = seq[a:b]
            per_region[nm] = {"seq": s, **comp(s), **stab.metrics(s)}
        per_region["all designed residues"] = {"seq": dseq, **comp(dseq), **stab.metrics(dseq)}
        per_region["whole chain H"] = {"seq": seq, **comp(seq), **stab.metrics(seq)}

        records.append((d, seq, hard, lab, per_region, c_pos))

    # ---- fasta ----
    os.makedirs(os.path.dirname(args.fasta), exist_ok=True)
    with open(args.fasta, "w") as fh:
        for d, seq, *_ in records:
            fh.write(f">{d} chain_H VHH caplacizumab-framework de_novo HER2_domain_IV len={len(seq)}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")

    # ---- report ----
    for d, seq, hard, lab, per, c_pos in records:
        print(f"\n=== {d} ===\n{seq}")
        print(f"  [1] cysteines: {hard['cys_total']} at positions {c_pos}  -> "
              f"{'EVEN' if hard['cys_even'] else 'ODD: HARD-FILTER REJECTION, CANNOT BE SUBMITTED'}"
              f"  ({hard['cys_in_designed_regions']} inside designed regions)")
        print(f"  [2] length: {hard['length']} aa -> {'OK (<250)' if hard['length_ok'] else 'TOO LONG'}")
        print(f"      canonical residues only: {hard['canonical_only']};  "
              f"matches the sequence given to Boltz-2: {hard['chain_H_matches_boltz2_query']}")
        print("  [3] motif hits (numbering = position in chain H, 1-based):")
        if not lab:
            print("      none")
        for name, a, b, where, s in lab:
            print(f"      {name:<20} {a}-{b}  {s:<8} {where}")
        print("  [3b] composition / [4] stability proxy by region "
              "(gly, pro, longest G/S run):")
        for nm, r in per.items():
            print(f"      {nm:<22} len {len(r['seq']):3d}  hydrophobic {r['hydro_frac']:.2f}  "
                  f"charge {r['net_charge']:+d} (His {r['n_his']})  Gly {r['gly_count']}  "
                  f"Pro {r['pro_count']}  G/S run {r['flex_run']}")

    import csv as _csv
    with open(args.csv, "w", newline="") as fh:
        w = _csv.writer(fh)
        w.writerow(["# De novo VHH pre-filters. Motif definitions from developability.py; "
                    "Gly/Pro/flex-run from stability_proxy.py (PROXY ONLY, not a Tm prediction). "
                    "No parent comparison: the CDRH3 parent is a 10-mer in another framework."])
        w.writerow(["design", "length", "cys_total", "cys_even", "cys_in_designed_regions",
                    "motif", "start", "end", "residues", "location"])
        for d, seq, hard, lab, per, c_pos in records:
            if not lab:
                w.writerow([d, hard["length"], hard["cys_total"], hard["cys_even"],
                            hard["cys_in_designed_regions"], "none", "", "", "", ""])
            for name, a, b, where, s in lab:
                w.writerow([d, hard["length"], hard["cys_total"], hard["cys_even"],
                            hard["cys_in_designed_regions"], name, a, b, s, where])
    print(f"\n# wrote {args.fasta}\n# wrote {args.csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
