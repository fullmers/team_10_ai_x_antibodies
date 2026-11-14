#!/usr/bin/env python3
"""Extract the designed CDRH3 loop from RFD4-Proteina `*_generated.cif` files.

The CDRH3 loop of the trastuzumab Fv in this challenge is chain B (auth), author
residue numbers 99-108 inclusive -- 10 positions. For each design we read the
`_atom_site` loop, collect one residue name per (chain, auth_seq_id), and emit
the one-letter sequence over that window.

Usage:
  extract_cdrh3.py RUNDIR [RUNDIR ...]          # a design/ or inference_samples/ dir
  extract_cdrh3.py --chain B --first 99 --last 108 ...

Output is one line per design: "<loop>  <design_id>", which is exactly the input
format novelty_check.py takes via --loops.
"""

from __future__ import annotations

import argparse
import os
import re
import sys

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

NUM_RE = re.compile(r"^-?\d+$")


def parse_atom_site_columns(path: str) -> tuple[list[str], int]:
    """Return the ordered _atom_site column names and the line index where rows start."""
    cols: list[str] = []
    with open(path) as fh:
        lines = fh.readlines()
    start = None
    for i, line in enumerate(lines):
        s = line.strip()
        if s.startswith("_atom_site."):
            cols.append(s.split(".", 1)[1].split()[0])
            start = i + 1
        elif cols and (s.startswith("ATOM") or s.startswith("HETATM")):
            start = i
            break
    if not cols or start is None:
        raise SystemExit(f"{path}: no _atom_site loop found")
    return cols, start


def extract_loop(path: str, chain: str, first: int, last: int) -> str:
    cols, start = parse_atom_site_columns(path)
    try:
        i_comp = cols.index("auth_comp_id")
    except ValueError:
        i_comp = cols.index("label_comp_id")
    try:
        i_asym = cols.index("auth_asym_id")
    except ValueError:
        i_asym = cols.index("label_asym_id")
    try:
        i_seq = cols.index("auth_seq_id")
    except ValueError:
        i_seq = cols.index("label_seq_id")

    residues: dict[int, str] = {}
    with open(path) as fh:
        for _ in range(start):
            fh.readline()
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                if s.startswith("#"):
                    break
                continue
            if not (s.startswith("ATOM") or s.startswith("HETATM")):
                break
            f = s.split()
            if max(i_comp, i_asym, i_seq) >= len(f):
                continue
            if f[i_asym] != chain:
                continue
            if not NUM_RE.match(f[i_seq]):
                continue
            n = int(f[i_seq])
            if first <= n <= last:
                residues.setdefault(n, f[i_comp])

    missing = [n for n in range(first, last + 1) if n not in residues]
    if missing:
        raise SystemExit(
            f"{path}: chain {chain} missing residues {missing} in window {first}-{last}"
        )
    return "".join(THREE_TO_ONE.get(residues[n], "X") for n in range(first, last + 1))


def find_designs(root: str) -> list[tuple[str, str]]:
    """Return (design_id, cif_path), sorted by the trailing integer when present."""
    hits: list[tuple[str, str]] = []
    for dirpath, _dirnames, filenames in os.walk(root):
        for fn in filenames:
            if fn.endswith("_generated.cif"):
                hits.append((fn[: -len("_generated.cif")], os.path.join(dirpath, fn)))

    def key(item: tuple[str, str]):
        m = re.search(r"(\d+)$", item[0])
        return (0, int(m.group(1))) if m else (1, 0)

    return sorted(hits, key=key)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("roots", nargs="+", help="directory holding *_generated.cif files")
    ap.add_argument("--chain", default="B", help="auth chain id (default: B)")
    ap.add_argument("--first", type=int, default=99, help="first auth_seq_id (default: 99)")
    ap.add_argument("--last", type=int, default=108, help="last auth_seq_id (default: 108)")
    args = ap.parse_args()

    n = 0
    for root in args.roots:
        for design_id, path in find_designs(root):
            print(f"{extract_loop(path, args.chain, args.first, args.last)}  {design_id}")
            n += 1
    if n == 0:
        print("no *_generated.cif found", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
