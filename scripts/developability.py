#!/usr/bin/env python3
"""Developability triage for designed CDRH3 loops, against the parent loop.

Scope and caveats, stated up front because they bound what this means:

  * This is a sequence-motif screen over a 10-residue loop. It is a cheap
    prefilter, not a developability assessment -- no aggregation prediction,
    no structural context, no experimental surrogate.
  * Motifs are scanned in FRAMEWORK CONTEXT, not on the loop in isolation.
    In this construct chain B reads ...Y Y C S R [loop 99-108] Y W G Q G...
    so the loop sits between `SR` and `YW`. Both flanking Y and W are
    hydrophobic, so a loop ending in hydrophobic residues extends a run across
    the junction -- invisible if you scan the loop alone. Motif hits are
    reported for the flanked context; `*_loop_only` columns give the isolated
    count so junction-created hits are identifiable as (context - loop_only).
  * Only loop positions can differ between designs; flank hits are shared by
    every design including the parent, which is why comparison is to the
    parent run through the same pipeline rather than to an absolute threshold.

Liabilities scanned (Lu et al. / Jain et al. style chemical liabilities):
  NG            deamidation hotspot
  DG, DP        isomerisation / acid-labile cleavage
  N-x-S/T       N-glycosylation sequon (x != P)
  C             free cysteine -- the conserved C96 is framework; any cysteine
                in the loop is unpaired and is flagged unconditionally
  5+ hydrophobic run   aggregation / HIC retention risk (AILMFVWY)
  3+ K/R run           patch of positive charge, polyreactivity risk

Also reported: hydrophobic fraction over the loop, and net charge at pH 7.4
as (K+R) - (D+E), histidine excluded and counted separately.

"Worse than parent" is per metric:
  * any liability count above the parent's count for that motif
  * hydrophobic fraction above the parent's
  * net charge further from neutral than the parent's, or positive when the
    parent is negative (positive CDR charge tracks with polyreactivity and
    fast clearance; the parent is -2)

Usage:
  developability.py --loops FILE [--subset FILE] [--csv OUT]
"""

from __future__ import annotations

import argparse
import csv
import re
import sys

PARENT = "WGGDGFYAMD"
FLANK_N = "SR"      # chain B 97-98, immediately before loop position 99
FLANK_C = "YW"      # chain B 109-110, immediately after loop position 108

HYDROPHOBIC = set("AILMFVWY")
POSITIVE = set("KR")
NEGATIVE = set("DE")

HYDRO_RUN_MIN = 5
KR_RUN_MIN = 3


def longest_run(seq: str, alphabet: set[str]) -> int:
    best = cur = 0
    for ch in seq:
        cur = cur + 1 if ch in alphabet else 0
        best = max(best, cur)
    return best


def count_runs(seq: str, alphabet: set[str], minimum: int) -> int:
    """Number of maximal runs of at least `minimum` consecutive members."""
    return sum(1 for m in re.finditer(f"[{''.join(sorted(alphabet))}]+", seq)
               if len(m.group()) >= minimum)


def sequons(seq: str) -> int:
    """N-x-S/T with x != P, counted with overlap allowed."""
    return len(re.findall(r"(?=N[^P][ST])", seq))


def scan(seq: str) -> dict:
    return {
        "NG": seq.count("NG"),
        "DG": seq.count("DG"),
        "DP": seq.count("DP"),
        "sequon": sequons(seq),
        "free_C": seq.count("C"),
        "hydro_run": count_runs(seq, HYDROPHOBIC, HYDRO_RUN_MIN),
        "KR_run": count_runs(seq, POSITIVE, KR_RUN_MIN),
    }


MOTIFS = ["NG", "DG", "DP", "sequon", "free_C", "hydro_run", "KR_run"]

PRETTY = {
    "NG": "NG", "DG": "DG", "DP": "DP", "sequon": "N-x-S/T",
    "free_C": "free C", "hydro_run": f"{HYDRO_RUN_MIN}+ hydrophobic run",
    "KR_run": f"{KR_RUN_MIN}+ K/R run",
}


def profile(loop: str) -> dict:
    ctx = FLANK_N + loop + FLANK_C
    ctx_hits = scan(ctx)
    loop_hits = scan(loop)
    # A loop cysteine is unpaired regardless of context; framework C96 is outside
    # the scanned window, so every C seen here is a loop cysteine.
    ctx_hits["free_C"] = loop_hits["free_C"] = loop.count("C")

    n_hydro = sum(1 for ch in loop if ch in HYDROPHOBIC)
    pos = sum(1 for ch in loop if ch in POSITIVE)
    neg = sum(1 for ch in loop if ch in NEGATIVE)

    return {
        "loop": loop,
        "ctx": ctx_hits,
        "loop_only": loop_hits,
        "hydro_frac": n_hydro / len(loop),
        "net_charge": pos - neg,
        "n_his": sum(1 for ch in loop if ch == "H"),
        "longest_hydro_run_ctx": longest_run(ctx, HYDROPHOBIC),
    }


def compare(p: dict, par: dict) -> list[str]:
    """Reasons this loop is worse than the parent; empty list means it is not."""
    why = []
    for m in MOTIFS:
        d = p["ctx"][m] - par["ctx"][m]
        if d > 0:
            junction = (p["ctx"][m] - p["loop_only"][m]) > (par["ctx"][m] - par["loop_only"][m])
            why.append(f"{PRETTY[m]} +{d}" + (" (junction)" if junction else ""))
    if p["hydro_frac"] > par["hydro_frac"]:
        why.append(f"hydrophobic {p['hydro_frac']:.0%} vs {par['hydro_frac']:.0%}")
    if abs(p["net_charge"]) > abs(par["net_charge"]):
        why.append(f"charge {p['net_charge']:+d} vs {par['net_charge']:+d} (further from neutral)")
    elif p["net_charge"] > 0 >= par["net_charge"]:
        why.append(f"charge {p['net_charge']:+d} (positive)")
    return why


def parse_loops(path: str) -> list[tuple[str, str]]:
    out = []
    with open(path) as fh:
        for line in fh:
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            f = line.split()
            out.append((f[0].upper(), f[1] if len(f) > 1 else f[0].upper()))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--loops", required=True, help="loop file, '<loop>  <label>' per line")
    ap.add_argument("--subset", help="file of labels to keep, one per line")
    ap.add_argument("--csv", help="write results here")
    ap.add_argument("--parent", default=PARENT)
    args = ap.parse_args()

    designs = parse_loops(args.loops)
    if args.subset:
        keep = {l.strip() for l in open(args.subset) if l.strip()}
        designs = [(s, lb) for s, lb in designs if lb in keep]
        missing = keep - {lb for _, lb in designs}
        if missing:
            print(f"# warning: {len(missing)} subset labels not in loop file: "
                  f"{sorted(missing)[:5]}", file=sys.stderr)
    if not designs:
        ap.error("no designs selected")

    par = profile(args.parent.upper())

    header = (["design", "cdrh3"]
              + [f"{m}" for m in MOTIFS]
              + [f"{m}_loop_only" for m in MOTIFS]
              + ["longest_hydro_run_ctx", "hydro_frac", "net_charge", "n_his",
                 "hydro_frac_vs_parent", "net_charge_vs_parent",
                 "worse_than_parent", "reasons"])

    rows = []
    for loop, label in sorted(set(designs), key=lambda d: d[1]):
        p = profile(loop)
        why = compare(p, par)
        rows.append([label, loop]
                    + [p["ctx"][m] for m in MOTIFS]
                    + [p["loop_only"][m] for m in MOTIFS]
                    + [p["longest_hydro_run_ctx"],
                       f"{p['hydro_frac']:.2f}", p["net_charge"], p["n_his"],
                       f"{p['hydro_frac'] - par['hydro_frac']:+.2f}",
                       f"{p['net_charge'] - par['net_charge']:+d}",
                       "yes" if why else "no", "; ".join(why)])

    print(f"# parent {args.parent}: " + ", ".join(
        f"{PRETTY[m]}={par['ctx'][m]}" for m in MOTIFS))
    print(f"# parent hydrophobic fraction {par['hydro_frac']:.0%}, "
          f"net charge {par['net_charge']:+d}, "
          f"longest hydrophobic run in context {par['longest_hydro_run_ctx']}")
    print(f"# context scanned: {FLANK_N}[loop]{FLANK_C}\n")

    show = ["design", "cdrh3", "hydro_frac", "net_charge", "worse_than_parent", "reasons"]
    idx = [header.index(c) for c in show]
    disp = [[str(r[i]) for i in idx] for r in rows]
    widths = [max(len(x) for x in [show[j]] + [d[j] for d in disp]) for j in range(len(show))]
    fmt = "  ".join(f"{{:<{w}}}" for w in widths)
    print(fmt.format(*show))
    print(fmt.format(*("-" * w for w in widths)))
    for d in disp:
        print(fmt.format(*d))

    n_worse = sum(1 for r in rows if r[header.index("worse_than_parent")] == "yes")
    print(f"\n# {n_worse}/{len(rows)} flagged worse than parent on at least one metric")
    print(f"# {len(rows) - n_worse}/{len(rows)} clean")

    if args.csv:
        with open(args.csv, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(header)
            w.writerows(rows)
        print(f"# wrote {args.csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
