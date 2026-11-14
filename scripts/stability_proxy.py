#!/usr/bin/env python3
"""Backbone-flexibility proxy for designed CDRH3 loops.  NOT a Tm prediction.

What this is: three crude counts on the 10-residue loop, plus the refolded
binder pLDDT, each shown next to the parent loop.  Glycine adds backbone
flexibility and so raises the entropic cost of ordering the loop; proline
restricts it.  That is a heuristic from loop biophysics, nothing here has been
calibrated against measured thermal stability, and a loop-only count says
nothing about framework packing, VH-VL pairing or the antigen interface.

Metrics (loop residues only; the flanking `SR`/`YW` do not touch the loop's
glycines or serines, so scanning the loop alone is exact for these counts):
  gly_count        number of G
  pro_count        number of P
  flex_run         longest run of consecutive G or S (the classic flexible
                   linker alphabet)
  gly_run          longest run of consecutive G
  binder_plddt     mean pLDDT of the whole refolded binder from
                   refold_metrics.csv -- a global confidence score, not a
                   loop-local one, and a weak stability proxy at best

Direction used for the "worse than parent" flag:
  gly_count  higher than parent  -> worse
  flex_run   longer than parent  -> worse
  pro_count  reported only: proline can rigidify a loop but also impose kinks
             or cis/trans heterogeneity, so a difference is not scored as
             better or worse.
The parent was not refolded in this campaign, so it has no binder_plddt;
pLDDT is compared to the median of the designs instead.
"""

from __future__ import annotations

import argparse
import csv
import re
import statistics

PARENT = "WGGDGFYAMD"
FLEX = set("GS")


def longest_run(seq: str, alphabet: set[str]) -> int:
    best = cur = 0
    for ch in seq:
        cur = cur + 1 if ch in alphabet else 0
        best = max(best, cur)
    return best


def metrics(loop: str) -> dict:
    return {
        "gly_count": loop.count("G"),
        "pro_count": loop.count("P"),
        "flex_run": longest_run(loop, FLEX),
        "gly_run": longest_run(loop, {"G"}),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metrics", default="data/docs/runs/her2-main/scored-boltz2/refold_metrics.csv")
    ap.add_argument("--loops", default="scripts/designs200.txt")
    ap.add_argument("--subset", default="scripts/strict_refolded.txt")
    ap.add_argument("--csv", default="scripts/stability_proxy.csv")
    ap.add_argument("--picks", default="SLSLSNYAFD,TRAAVGYAHD,DLSYYNHAFD,WSSAVGYAHD")
    args = ap.parse_args()

    loops = {}
    for line in open(args.loops):
        s, lab = line.split()[:2]
        loops[lab] = s
    plddt = {r["query_id"]: float(r["binder_plddt"])
             for r in csv.DictReader(open(args.metrics))}
    keep = [l.strip() for l in open(args.subset) if l.strip()]
    picks = args.picks.split(",")

    par = metrics(PARENT)
    med_plddt = statistics.median(plddt[d] for d in keep)

    rows = []
    for d in keep:
        loop = loops[d]
        m = metrics(loop)
        why = []
        if m["gly_count"] > par["gly_count"]:
            why.append(f"glycine {m['gly_count']} > {par['gly_count']}")
        if m["flex_run"] > par["flex_run"]:
            why.append(f"flexible G/S run {m['flex_run']} > {par['flex_run']}")
        rows.append({
            "design": d, "cdrh3": loop, "pick": "yes" if loop in picks else "",
            **m,
            "gly_vs_parent": m["gly_count"] - par["gly_count"],
            "pro_vs_parent": m["pro_count"] - par["pro_count"],
            "flex_run_vs_parent": m["flex_run"] - par["flex_run"],
            "binder_plddt": f"{plddt[d]:.2f}",
            "plddt_vs_design_median": f"{plddt[d] - med_plddt:+.2f}",
            "flexibility_worse_than_parent": "yes" if why else "no",
            "reasons": "; ".join(why),
        })
    rows.sort(key=lambda r: (r["pick"] != "yes", -float(r["binder_plddt"])))

    fields = list(rows[0].keys())
    with open(args.csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["# PROXY ONLY, NOT A Tm PREDICTION. Loop-level glycine/proline/flexible-run "
                    "counts plus whole-binder refold pLDDT. Uncalibrated against any stability "
                    "measurement."])
        w.writerow([f"# parent {PARENT}: gly={par['gly_count']} pro={par['pro_count']} "
                    f"flex_run(G/S)={par['flex_run']} gly_run={par['gly_run']}; "
                    f"parent has no refold pLDDT (not refolded); "
                    f"pLDDT compared with the median of the {len(keep)} designs = {med_plddt:.2f}"])
        w.writerow(fields)
        for r in rows:
            w.writerow([r[f] for f in fields])

    print(f"# PROXY, not a Tm prediction.  parent {PARENT}: gly={par['gly_count']} "
          f"pro={par['pro_count']} G/S-run={par['flex_run']} gly-run={par['gly_run']}")
    print(f"# design median binder_plddt = {med_plddt:.2f} (parent not refolded)\n")
    hdr = f"{'design':16}{'CDRH3':12}{'gly':>4}{'pro':>4}{'G/S run':>8}{'gly run':>8}{'pLDDT':>8}{'vs med':>8}  verdict"
    print(hdr); print("-" * len(hdr))
    for r in rows:
        if r["pick"] == "yes":
            print(f"{r['design']:16}{r['cdrh3']:12}{r['gly_count']:>4}{r['pro_count']:>4}"
                  f"{r['flex_run']:>8}{r['gly_run']:>8}{r['binder_plddt']:>8}"
                  f"{r['plddt_vs_design_median']:>8}  "
                  f"{'WORSE: ' + r['reasons'] if r['reasons'] else 'not worse'}")
    n_worse = sum(1 for r in rows if r["flexibility_worse_than_parent"] == "yes")
    n_gly = sum(1 for r in rows if r["gly_vs_parent"] > 0)
    print(f"\n# all {len(rows)} designs: {n_worse} worse than parent on flexibility "
          f"({n_gly} with more glycines)")
    print(f"# glycine counts: " + ", ".join(
        f"{g}G x{sum(1 for r in rows if r['gly_count'] == g)}"
        for g in sorted({r['gly_count'] for r in rows})))
    print(f"# wrote {args.csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
