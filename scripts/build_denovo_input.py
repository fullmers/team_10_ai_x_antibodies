#!/usr/bin/env python3
"""Build a de novo VHH-vs-HER2 design input for RFD4-Proteina's antibody model.

The antibody model takes a single mmCIF that encodes the whole design spec in
the `_atom_site` loop. This script grafts two existing inputs into one:

  framework : chain H of Data/03-vhh/pmvk_caplacizumab.cif.zst -- the
              caplacizumab VHH framework (86 resolved residues) carrying three
              HETATM/UNK placeholders with NaN coordinates whose comp_id gives
              the CDR length range (H23 "13-13", H38 "9-10", H76 "12-21") and
              `condition_sequence` "<M>". PMVK (chain A) is dropped.
  antigen   : chains C (489-580) and D (591-607) of
              Data/12-her2/her2_fv_1n8z.cif -- HER2 domain IV, split at the
              unresolved 581-590 loop because validation rejects a numbering
              gap inside one chain. The trastuzumab Fv (chains A, B) is dropped.

Column conventions copied from the caplacizumab file:
  antigen backbone  cath `x.x.x.x`, CoordinateGroup 0 (N, CA, C, O only)
  antigen sidechain cath `x.x.x.x`, CoordinateGroup -1; hotspot atoms are sidechain atoms
  framework atom    cath `2.60.40.x`, CoordinateGroup 1, sequence = residue
  designed CDR      comp_id = length range, sequence `<M>`, group -1, NaN xyz

Hotspots are the four HER2 atoms nearest trastuzumab's CDRH3 in 1N8Z, taken
verbatim from Data/12-her2/her2_cdrh3.yaml, so the VHH is aimed at the same
epitope as trastuzumab.

Usage: build_denovo_input.py --out runs/denovo/input/her2_vhh_denovo.cif
"""

from __future__ import annotations

import argparse
import os
import subprocess

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

# chain, auth_seq_id, atom name  -- from Data/12-her2/her2_cdrh3.yaml C_HOT
HOTSPOTS = {("D", 593, "NZ"), ("C", 573, "CD1"), ("D", 591, "CG2"), ("C", 572, "CD")}

ANTIGEN_CHAINS = ("C", "D")


def read_atom_site(path: str):
    """Return (header_text, columns, rows, tail_text) for an mmCIF."""
    with open(path) as fh:
        lines = fh.readlines()
    cols, start = [], None
    for i, line in enumerate(lines):
        s = line.strip()
        if s.startswith("_atom_site."):
            cols.append(s.split(".", 1)[1].split()[0])
            start = i + 1
        elif cols and s.startswith(("ATOM", "HETATM")):
            start = i
            break
    rows, end = [], len(lines)
    for i in range(start, len(lines)):
        s = lines[i].strip()
        if not s.startswith(("ATOM", "HETATM")):
            end = i
            break
        rows.append(s.split())
    return "".join(lines[:start]), cols, rows, "".join(lines[end:])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument("--framework-zst",
                    default=os.path.join(root, "data/docs/Data/03-vhh/pmvk_caplacizumab.cif.zst"))
    ap.add_argument("--antigen",
                    default=os.path.join(root, "data/docs/Data/12-her2/her2_fv_1n8z.cif"))
    ap.add_argument("--out", required=True)
    ap.add_argument("--antigen-mode", choices=["absolute", "relative"], default="absolute",
                    help="absolute: HER2 backbone fixed on coordinate group 0 (as the PMVK "
                         "example). relative: HER2 backbone is internal-geometry group 2 and "
                         "the framework group 1, with no absolute group (as the working "
                         "her2_cdrh3.yaml run).")
    args = ap.parse_args()

    # --- framework ------------------------------------------------------
    fw_cif = subprocess.run(["zstd", "-d", "-c", args.framework_zst],
                            capture_output=True, check=True).stdout.decode()
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    tmp = args.out + ".framework.tmp"
    with open(tmp, "w") as fh:
        fh.write(fw_cif)
    head, cols, fw_rows, tail = read_atom_site(tmp)
    os.remove(tmp)

    idx = {c: i for i, c in enumerate(cols)}
    fw_rows = [r for r in fw_rows if r[idx["auth_asym_id"]] == "H"]
    if not fw_rows:
        raise SystemExit("no chain H atoms in framework")

    # --- antigen --------------------------------------------------------
    _h2, a_cols, a_rows, _t2 = read_atom_site(args.antigen)
    aidx = {c: i for i, c in enumerate(a_cols)}
    ag_rows = [r for r in a_rows if r[aidx["auth_asym_id"]] in ANTIGEN_CHAINS]
    if not ag_rows:
        raise SystemExit("no antigen atoms found")

    # entity ids: 1 = HER2 C, 3 = HER2 D, 2 = VHH framework (keeps the
    # framework's existing entity id so its rows stay untouched)
    entity_of = {"C": "1", "D": "3"}

    out_rows: list[list[str]] = []
    for r in ag_rows:
        ch = r[aidx["auth_asym_id"]]
        seq = r[aidx["auth_seq_id"]]
        comp = r[aidx["auth_comp_id"]]
        atom = r[aidx["auth_atom_id"]] if "auth_atom_id" in aidx else r[aidx["label_atom_id"]]
        hot = (ch, int(seq), atom) in HOTSPOTS
        row = ["."] * len(cols)
        row[idx["group_PDB"]] = "ATOM"
        row[idx["type_symbol"]] = r[aidx["type_symbol"]]
        row[idx["label_atom_id"]] = r[aidx["label_atom_id"]]
        row[idx["label_alt_id"]] = "."
        row[idx["label_comp_id"]] = comp
        row[idx["label_asym_id"]] = ch
        row[idx["label_entity_id"]] = entity_of[ch]
        row[idx["label_seq_id"]] = seq
        row[idx["pdbx_PDB_ins_code"]] = "."
        row[idx["auth_seq_id"]] = seq
        row[idx["auth_comp_id"]] = comp
        row[idx["auth_asym_id"]] = ch
        row[idx["auth_atom_id"]] = atom
        row[idx["B_iso_or_equiv"]] = r[aidx["B_iso_or_equiv"]] if "B_iso_or_equiv" in aidx else "0.0"
        row[idx["occupancy"]] = "1.0"
        row[idx["pdbx_formal_charge"]] = "?"
        row[idx["annotation_expsegmin_1_residue"]] = "-1"
        row[idx["annotation_expsegmax_1_residue"]] = "-1"
        row[idx["condition_sequence_1_residue"]] = comp
        row[idx["condition_hotspot_1_atom"]] = "True" if hot else "False"
        row[idx["condition_cath-label_1_residue"]] = "x.x.x.x"
        # Only backbone atoms are coordinate-conditioned (group 0 = antigen);
        # side chains, including hotspot atoms, are -1, as in the PMVK example.
        ag_group = "0" if args.antigen_mode == "absolute" else "2"
        row[idx["condition_CoordinateGroup_1_atom"]] = (
            ag_group if atom in ("N", "CA", "C", "O") else "-1")
        row[idx["Cartn_x"]] = r[aidx["Cartn_x"]]
        row[idx["Cartn_y"]] = r[aidx["Cartn_y"]]
        row[idx["Cartn_z"]] = r[aidx["Cartn_z"]]
        row[idx["pdbx_PDB_model_num"]] = "1"
        out_rows.append(row)

    n_hot = sum(1 for r in out_rows if r[idx["condition_hotspot_1_atom"]] == "True")
    if n_hot != len(HOTSPOTS):
        raise SystemExit(f"expected {len(HOTSPOTS)} hotspot atoms, marked {n_hot}")

    # The PMVK example lists the antibody chain first, then the antigen.
    out_rows = list(fw_rows) + out_rows

    # renumber ids and array indices across the merged file
    for i, r in enumerate(out_rows):
        r[idx["id"]] = str(i + 1)
        r[idx["atom_array_index"]] = str(i)

    # --- entity blocks ---------------------------------------------------
    def chain_residues(rows, chain, i):
        seen = {}
        for r in rows:
            if r[i["auth_asym_id"]] == chain:
                seen.setdefault(int(r[i["auth_seq_id"]]), r[i["auth_comp_id"]])
        return [(n, seen[n]) for n in sorted(seen)]

    her2 = {c: chain_residues(ag_rows, c, aidx) for c in ANTIGEN_CHAINS}
    fw_res = chain_residues(fw_rows, "H", idx)

    def one_letter(res):
        return "".join(THREE_TO_ONE.get(c, "X") for _, c in res)

    def wrap(s, n=80):
        return "\n".join(s[i:i + n] for i in range(0, len(s), n))

    ent_lines = ["data_her2_vhh_denovo", "#",
                 "_entry.id       her2_vhh_denovo.cif",
                 "_entry.author   RFD4-Proteina",
                 "_entry.date     2026-11-14",
                 "#", "loop_", "_entity.id ", "_entity.type ",
                 "_entity.pdbx_description ",
                 "1 polymer .", "2 polymer .", "3 polymer .", "#",
                 "loop_", "_entity_poly.entity_id ", "_entity_poly.type ",
                 "_entity_poly.nstd_linkage ", "_entity_poly.nstd_monomer ",
                 "_entity_poly.pdbx_seq_one_letter_code ",
                 "_entity_poly.pdbx_seq_one_letter_code_can ",
                 "_entity_poly.pdbx_strand_id ",
                 "_entity_poly.pdbx_target_identifier "]

    fw_seq = ""
    for n, comp in fw_res:
        fw_seq += "X" if comp not in THREE_TO_ONE else THREE_TO_ONE[comp]

    for ent, chain, seq, nstd in (("1", "C", one_letter(her2["C"]), "no"),
                                  ("2", "H", fw_seq, "yes"),
                                  ("3", "D", one_letter(her2["D"]), "no")):
        ent_lines += [f"{ent} polypeptide(L) no {nstd} ", ";" + wrap(seq), ";", " ",
                      ";" + wrap(seq), ";", f" {chain} ?"]
    ent_lines += ["#", "loop_", "_entity_poly_seq.entity_id ",
                  "_entity_poly_seq.hetero ", "_entity_poly_seq.mon_id ",
                  "_entity_poly_seq.num "]
    for ent, res in (("1", her2["C"]), ("2", fw_res), ("3", her2["D"])):
        for num, comp in res:
            ent_lines.append(f"{ent} n {comp:5s} {num}")
    ent_lines.append("#")

    # chem_comp / struct_conn / chem_comp_bond blocks come from the framework
    # file unchanged: the placeholder comp ids and CDR linkage records live
    # there, and every antigen residue is a standard amino acid already listed.
    keep = head.split("#\n")
    tail_blocks = []
    for block in keep:
        if block.lstrip().startswith(("loop_\n_chem_comp.", "loop_\n_struct_conn.",
                                      "loop_\n_chem_comp_bond.")):
            tail_blocks.append(block)
    body = "#\n".join(tail_blocks)

    atom_hdr = "loop_\n" + "".join(f"_atom_site.{c} \n" for c in cols)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write("\n".join(ent_lines) + "\n")
        if body:
            fh.write(body.rstrip() + "\n#\n")
        fh.write(atom_hdr)
        for r in out_rows:
            fh.write(" ".join(r) + "\n")
        fh.write("#\n")

    print(f"wrote {args.out}")
    print(f"  HER2 C: {len(her2['C'])} residues ({her2['C'][0][0]}-{her2['C'][-1][0]})")
    print(f"  HER2 D: {len(her2['D'])} residues ({her2['D'][0][0]}-{her2['D'][-1][0]})")
    print(f"  VHH  H: {len(fw_res)} residues, {sum(1 for _, c in fw_res if c not in THREE_TO_ONE)} CDR placeholders")
    print(f"  atoms: {len(out_rows)}, hotspots marked: {n_hot}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
