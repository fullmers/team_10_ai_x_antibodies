#!/usr/bin/env python3
"""CDR-level novelty of the de novo VHH designs against the SAbDab2 ML subset.

SCOPE -- read first.  The reference set is the SAbDab2 *ML dataset* (Zenodo
20083995, CC-BY-4.0): a curated subset of 15,641 antibody structures from 8,641
PDB entries, 5,301 unique antibodies.  It is NOT full SAbDab2 (21,237 distinct
instances) and nowhere near all known antibodies (no OAS, no patents, no
therapeutics outside the PDB).  "Novel" below means "novel against this subset".

Why CDRs, not the full sequence.  Every design keeps the caplacizumab VHH
framework (83 fixed residues of 118-127), so identity of the full sequence to
caplacizumab-like VHHs is ~65-70% by construction before any CDR is considered.
A full-length "<70% to any known antibody" test cannot be passed by this design
strategy, so identity is measured on the designed CDRs.

Region definitions (design side):
  CDR3  residues between Cys104 and the J-region Trp: exactly the IMGT CDR3 and
        exactly the designed segment (placeholder 12-21).  Validated below
        against the database's own CDRH3 annotation.
  CDR1  placeholder is 13 residues = IMGT 24-40 minus the 4 CDR1 gap positions,
        so IMGT CDR1 (8 residues, 27-30 + 35-38) is segment[3:11].  Validated
        below against the database.
  CDR2  placeholder is 9-10 residues = IMGT 55 + CDR2 + 66, so CDR2 is
        segment[1:-1].  INFERRED from IMGT gap conventions; validated only
        indirectly (FR3 length), so treat the CDR1+2+3 result as softer than
        the CDR3 result.

Identity = LCS length / max(len a, len b).  LCS allows gaps and ignores
position, so it can only over-report similarity, which is the conservative
direction for a novelty claim.  Computed with a bit-parallel LCS that is
self-tested against a plain DP before use.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import random
import re
import sys
from collections import Counter

csv.field_size_limit(sys.maxsize)

FRAME = re.compile(
    r"^EVQLVESGGGLVQPGGSLRLSC(?P<s1>.{13})WFRQAPGKGRELVA(?P<s2>.{9,10})"
    r"YPDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYC(?P<c3>.{12,21})WGQGTQVTVS$")


def masks(a: str) -> dict:
    m: dict = {}
    for i, ch in enumerate(a):
        m[ch] = m.get(ch, 0) | (1 << i)
    return m


def lcs_len(a_masks: dict, m: int, b: str) -> int:
    """Bit-parallel LCS length (Hyyro / Crochemore et al.)."""
    full = (1 << m) - 1
    v = full
    for ch in b:
        mk = a_masks.get(ch, 0)
        u = v & mk
        v = ((v + u) | (v & ~mk)) & full
    return m - bin(v).count("1")


def lcs_dp(a: str, b: str) -> int:
    prev = [0] * (len(b) + 1)
    for x in a:
        cur = [0]
        for j, y in enumerate(b, 1):
            cur.append(prev[j - 1] + 1 if x == y else max(prev[j], cur[j - 1]))
        prev = cur
    return prev[-1]


def selftest() -> None:
    rnd = random.Random(7)
    aa = "ACDEFGHIKLMNPQRSTVWY"
    for _ in range(3000):
        a = "".join(rnd.choice(aa) for _ in range(rnd.randint(1, 25)))
        b = "".join(rnd.choice(aa) for _ in range(rnd.randint(1, 25)))
        assert lcs_len(masks(a), len(a), b) == lcs_dp(a, b), (a, b)


def identity(a: str, a_masks: dict, b: str) -> float:
    return lcs_len(a_masks, len(a), b) / max(len(a), len(b))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", required=True, help="splits_final/ab_split.csv from the SAbDab2 ML dataset")
    ap.add_argument("--run", default="data/docs/runs/denovo")
    ap.add_argument("--csv", default="scripts/novelty_sabdab2.csv")
    ap.add_argument("--threshold", type=float, default=0.70)
    args = ap.parse_args()

    selftest()
    print("# LCS self-test vs plain DP: 3000 random pairs OK")

    # ---- database ----
    rows = []
    with open(args.db, newline="") as fh:
        for r in csv.DictReader(fh):
            if r["CDRH3"]:
                rows.append(r)
    n_all = sum(1 for _ in csv.DictReader(open(args.db, newline="")))
    print(f"# database rows: {n_all}; with a heavy chain / CDRH3: {len(rows)}; "
          f"unique SAbDab2 heavy ids: {len({r['HEAVY_ID'] for r in rows})}; "
          f"types: {dict(Counter(r['type'] for r in rows))}")

    # ---- validation of region definitions against the database ----
    ok3 = tot3 = 0
    ok1 = tot1 = 0
    fr3_len = Counter()
    for r in rows:
        h = r["Hseq"]
        c3 = r["CDRH3"]
        i = h.find("C" + c3)
        if i >= 0:
            tot3 += 1
            if h[i + 1 + len(c3):i + 2 + len(c3)] == "W":
                ok3 += 1
        c1, c2 = r["CDRH1"], r["CDRH2"]
        if len(c1) == 8:
            tot1 += 1
            if re.search(r"C.{3}" + re.escape(c1) + r".{2}W[FVILY][RQ]", h):
                ok1 += 1
        j2 = h.find(c2)
        if c2 and i >= 0 and j2 >= 0 and j2 < i:
            fr3_len[i - (j2 + len(c2))] += 1   # residues from CDR2 end up to Cys104 (exclusive)
    print(f"# validation CDR3 = 'C'+CDRH3+'W': {ok3}/{tot3} database heavy chains")
    print(f"# validation CDR1 (8-mer) at C+3 ... +2+W[FVILY][RQ]: {ok1}/{tot1}")
    print(f"# validation CDR2: residues between CDR2 end and Cys104, most common = "
          f"{fr3_len.most_common(3)} (design predicts 37)")

    # ---- designs ----
    queries = json.load(open(os.path.join(args.run, "prepared/holo/queries.json")))["queries"]
    seqs = {q["query_id"]: next(c["sequence"] for c in q["chains"] if c["chain_id"] == "H")
            for q in queries}
    metrics = {r["query_id"]: r for r in csv.DictReader(
        open(os.path.join(args.run, "scored-boltz2/refold_metrics.csv")))}
    designs = []
    for d, s in seqs.items():
        m = FRAME.match(s)
        if not m:
            raise SystemExit(f"{d}: sequence does not match the caplacizumab frame pattern")
        s1, s2, c3 = m.group("s1"), m.group("s2"), m.group("c3")
        designs.append({"design": d, "full": s, "cdr1": s1[3:11], "cdr2": s2[1:-1],
                        "cdr3": c3, "min_ipae": float(metrics[d]["min_ipae"]),
                        "iptm": float(metrics[d]["iptm"])})
    n_fixed = sum(len(x) for x in ("EVQLVESGGGLVQPGGSLRLSC", "WFRQAPGKGRELVA",
                                   "YPDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYC", "WGQGTQVTVS"))
    lens = [len(d["full"]) for d in designs]
    print(f"# {len(designs)} designs; fixed framework residues {n_fixed}; "
          f"design lengths {min(lens)}-{max(lens)} => framework is "
          f"{100 * n_fixed / max(lens):.0f}-{100 * n_fixed / min(lens):.0f}% of each sequence")

    # ---- compare ----
    uniq3 = {}
    uniq123 = {}
    for r in rows:
        uniq3.setdefault(r["CDRH3"], r)
        uniq123.setdefault(r["CDRH1"] + r["CDRH2"] + r["CDRH3"], r)
    print(f"# unique CDRH3 in database: {len(uniq3)}; unique CDRH1+2+3: {len(uniq123)}")

    out = []
    for d in designs:
        c123 = d["cdr1"] + d["cdr2"] + d["cdr3"]
        best3 = (0.0, None)
        a3 = masks(d["cdr3"])
        for seq, r in uniq3.items():
            idn = identity(d["cdr3"], a3, seq)
            if idn > best3[0]:
                best3 = (idn, r)
        best123 = (0.0, None)
        a123 = masks(c123)
        for seq, r in uniq123.items():
            idn = identity(c123, a123, seq)
            if idn > best123[0]:
                best123 = (idn, r)
        d.update({
            "cdr123": c123,
            "id_cdr3": best3[0], "nn_cdr3": best3[1]["CDRH3"], "nn_cdr3_pdb": best3[1]["PDB_ID"],
            "nn_cdr3_type": best3[1]["type"],
            "id_cdr123": best123[0],
            "nn_cdr123": best123[1]["CDRH1"] + best123[1]["CDRH2"] + best123[1]["CDRH3"],
            "nn_cdr123_pdb": best123[1]["PDB_ID"], "nn_cdr123_type": best123[1]["type"],
        })
        d["pass_cdr3"] = d["id_cdr3"] < args.threshold
        d["pass_cdr123"] = d["id_cdr123"] < args.threshold
        out.append(d)
    out.sort(key=lambda x: x["min_ipae"])

    fields = ["design", "min_ipae", "iptm", "cdr1", "cdr2", "cdr3", "id_cdr3", "nn_cdr3",
              "nn_cdr3_pdb", "nn_cdr3_type", "pass_cdr3", "id_cdr123", "nn_cdr123_pdb",
              "nn_cdr123_type", "pass_cdr123"]
    with open(args.csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([f"# Reference = SAbDab2 ML subset (15,641 structures, 5,301 unique antibodies), "
                    f"NOT full SAbDab. pass = max CDR LCS identity < {args.threshold:.0%}. "
                    f"CDR2 boundaries inferred; CDR3 validated."])
        w.writerow(fields)
        for d in out:
            w.writerow([f"{d[f]:.3f}" if isinstance(d[f], float) else d[f] for f in fields])

    n = len(out)
    print(f"\n# threshold: max identity < {args.threshold:.0%}")
    for t in (0.5, 0.6, 0.7, 0.8, 0.9):
        print(f"#   CDR3   max identity < {t:.0%}: {sum(1 for d in out if d['id_cdr3'] < t):2d}/{n}"
              f"    CDR1+2+3 < {t:.0%}: {sum(1 for d in out if d['id_cdr123'] < t):2d}/{n}")
    print(f"\n{'design':24}{'ipae':>6}{'iptm':>6}  {'CDR3':22}{'id3':>6}  nearest CDR3 (type, pdb){'':14}{'id123':>6}")
    for d in out[:10]:
        print(f"{d['design']:24}{d['min_ipae']:6.2f}{d['iptm']:6.2f}  {d['cdr3']:22}{d['id_cdr3']:6.2f}  "
              f"{d['nn_cdr3']:22} {d['nn_cdr3_type']:5} {d['nn_cdr3_pdb'][-4:]}{d['id_cdr123']:9.2f}")
    print(f"\n# wrote {args.csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
