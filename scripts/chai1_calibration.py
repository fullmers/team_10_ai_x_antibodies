#!/usr/bin/env python3
"""ipTM calibration from the Champloo VHH-antigen benchmark (Chai-1 replicas).

CROSS-MODEL WARNING -- read before using any threshold this produces.
The confidences here are **Chai-1**. The HER2 CDRH3 designs in this project were
refolded and scored with **Boltz-2**. ipTM is not calibrated across predictors:
`06-evaluate.md` and `16-refolding.md` both say to compare only within one
predictor/protocol, and that published cut-offs shift between models. A
threshold derived here is a sanity check on how a known non-binder scores in a
different model, NOT a decision boundary for our Boltz-2 numbers. The Boltz-2
version of this archive is 1.7 TB and does not fit on this VM.

Data: s3://apheris-oslo-hackathon/data/sigma2-ab-ag-champloo/
      chai1_all_per_replica_confidence_metrics.zip (md5 0b41d2974deffc6bfc34928607e038d2)
Layout: a flat zip of .npz files named
      system_{VHH_idx}_{antigen_idx}_{VHH_PDB}_{antigen_PDB}_trunk_{t}_scores.model_idx_{m}.npz
A pair is REAL when the VHH and antigen indices match (and, verified, the two
PDB ids match too); every other index combination is a SHUFFLED pair, which is
the study's non-binder control ("dreamy" in the AF3 archive).

Each .npz is itself a zip of .npy arrays; we read `iptm.npy` only. numpy is not
installed in this environment, so the scalar float32 is unpacked from the .npy
header + payload directly.

Writes a tidy CSV: one row per summary statistic, plus the per-pair aggregates
needed to reproduce them.
"""

from __future__ import annotations

import argparse
import ast
import csv
import io
import os
import re
import struct
import sys
import zipfile
from collections import defaultdict
from multiprocessing import Pool

NAME_RE = re.compile(
    r"^system_(\d+)_(\d+)_([0-9A-Za-z]+)_([0-9A-Za-z]+)"
    r"_trunk_(\d+)_scores\.model_idx_(\d+)\.npz$"
)

ARCHIVE = None  # per-worker handle


def read_npy_scalar(raw: bytes) -> float:
    """Return the single float in a .npy holding a shape-(1,) array."""
    if raw[:6] != b"\x93NUMPY":
        raise ValueError("not a .npy payload")
    major = raw[6]
    if major == 1:
        hlen = struct.unpack("<H", raw[8:10])[0]
        start = 10 + hlen
        header = raw[10:start]
    else:
        hlen = struct.unpack("<I", raw[8:12])[0]
        start = 12 + hlen
        header = raw[12:start]
    meta = ast.literal_eval(header.decode("latin1").strip())
    descr = meta["descr"]
    if descr not in ("<f4", "<f8"):
        raise ValueError(f"unexpected dtype {descr}")
    fmt = "<f" if descr == "<f4" else "<d"
    return struct.unpack(fmt, raw[start:start + struct.calcsize(fmt)])[0]


def _init(path: str) -> None:
    global ARCHIVE
    ARCHIVE = zipfile.ZipFile(path)


def _iptm_for(names: list[str]) -> list[tuple[str, float]]:
    out = []
    for n in names:
        try:
            blob = ARCHIVE.open(n).read()
            inner = zipfile.ZipFile(io.BytesIO(blob))
            out.append((n, read_npy_scalar(inner.open("iptm.npy").read())))
        except Exception as exc:  # keep going; report at the end
            out.append((n, float("nan")))
            print(f"# warn: {n}: {exc}", file=sys.stderr)
    return out


def quantile(sorted_vals: list[float], q: float) -> float:
    """Linear-interpolated quantile of an already-sorted list."""
    if not sorted_vals:
        return float("nan")
    if len(sorted_vals) == 1:
        return sorted_vals[0]
    pos = q * (len(sorted_vals) - 1)
    lo = int(pos)
    hi = min(lo + 1, len(sorted_vals) - 1)
    frac = pos - lo
    return sorted_vals[lo] * (1 - frac) + sorted_vals[hi] * frac


def mean(v: list[float]) -> float:
    return sum(v) / len(v) if v else float("nan")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("archive", help="chai1_all_per_replica_confidence_metrics.zip")
    ap.add_argument("--csv", default="scripts/calibration.csv")
    ap.add_argument("--pairs-csv", help="optional per-pair aggregate dump")
    ap.add_argument("--jobs", type=int, default=min(16, os.cpu_count() or 4))
    args = ap.parse_args()

    z = zipfile.ZipFile(args.archive)
    names = [n for n in z.namelist() if n.endswith(".npz")]
    z.close()
    print(f"# {len(names):,} replica files in archive", file=sys.stderr)

    chunks = [names[i::args.jobs * 4] for i in range(args.jobs * 4)]
    vals: list[tuple[str, float]] = []
    with Pool(args.jobs, initializer=_init, initargs=(args.archive,)) as pool:
        for i, part in enumerate(pool.imap_unordered(_iptm_for, chunks), 1):
            vals.extend(part)
            print(f"# chunk {i}/{len(chunks)}  ({len(vals):,} files)", file=sys.stderr)

    # group replicas by pair
    per_pair_real: dict[tuple, list[float]] = defaultdict(list)
    per_pair_shuf: dict[tuple, list[float]] = defaultdict(list)
    samples_real: list[float] = []
    samples_shuf: list[float] = []
    skipped = 0
    for n, v in vals:
        m = NAME_RE.match(n)
        if not m or v != v:  # unparsed name or NaN
            skipped += 1
            continue
        vhh, ag, p1, p2 = m.group(1), m.group(2), m.group(3), m.group(4)
        key = (vhh, ag, p1, p2)
        if vhh == ag:
            per_pair_real[key].append(v)
            samples_real.append(v)
        else:
            per_pair_shuf[key].append(v)
            samples_shuf.append(v)

    samples_real.sort()
    samples_shuf.sort()
    pair_mean_real = sorted(mean(v) for v in per_pair_real.values())
    pair_mean_shuf = sorted(mean(v) for v in per_pair_shuf.values())
    pair_max_shuf = sorted(max(v) for v in per_pair_shuf.values())
    pair_max_real = sorted(max(v) for v in per_pair_real.values())

    def frac_above(sorted_vals: list[float], t: float) -> float:
        lo, hi = 0, len(sorted_vals)
        while lo < hi:
            mid = (lo + hi) // 2
            if sorted_vals[mid] > t:
                hi = mid
            else:
                lo = mid + 1
        return (len(sorted_vals) - lo) / len(sorted_vals) if sorted_vals else float("nan")

    rows = []

    def add(population: str, n_items: str, vals_sorted: list[float]):
        rows.append({
            "population": population,
            "unit": n_items,
            "n": len(vals_sorted),
            "mean": f"{mean(vals_sorted):.4f}",
            "median": f"{quantile(vals_sorted, 0.50):.4f}",
            "p75": f"{quantile(vals_sorted, 0.75):.4f}",
            "p90": f"{quantile(vals_sorted, 0.90):.4f}",
            "p95": f"{quantile(vals_sorted, 0.95):.4f}",
            "p99": f"{quantile(vals_sorted, 0.99):.4f}",
            "max": f"{vals_sorted[-1]:.4f}" if vals_sorted else "",
            "frac_above_0.7": f"{frac_above(vals_sorted, 0.7):.4f}",
            "frac_above_0.8": f"{frac_above(vals_sorted, 0.8):.4f}",
            "frac_above_0.9": f"{frac_above(vals_sorted, 0.9):.4f}",
        })

    add("shuffled (non-binder)", "replica sample", samples_shuf)
    add("shuffled (non-binder)", "pair mean over replicas", pair_mean_shuf)
    add("shuffled (non-binder)", "pair max over replicas", pair_max_shuf)
    add("real (binder)", "replica sample", samples_real)
    add("real (binder)", "pair mean over replicas", pair_mean_real)
    add("real (binder)", "pair max over replicas", pair_max_real)

    fields = list(rows[0].keys())
    with open(args.csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["# MODEL: Chai-1. Our HER2 designs were scored with Boltz-2."
                    " ipTM does not transfer across predictors -- treat these"
                    " thresholds as indicative only."])
        w.writerow([f"# source: {os.path.basename(args.archive)};"
                    f" real pairs={len(per_pair_real)}, shuffled pairs={len(per_pair_shuf)},"
                    f" replica files={len(vals)}, skipped={skipped}"])
        w.writerow(fields)
        for r in rows:
            w.writerow([r[f] for f in fields])
        w.writerow([])
        w.writerow(["threshold excluding 95% of shuffled pairs (unit = pair mean)",
                    f"{quantile(pair_mean_shuf, 0.95):.4f}"])
        w.writerow(["threshold excluding 95% of shuffled replicas (unit = replica sample)",
                    f"{quantile(samples_shuf, 0.95):.4f}"])
        w.writerow(["threshold excluding 95% of shuffled pairs (unit = pair max)",
                    f"{quantile(pair_max_shuf, 0.95):.4f}"])
        w.writerow(["real pairs retained at the pair-mean 95% threshold",
                    f"{frac_above(pair_mean_real, quantile(pair_mean_shuf, 0.95)):.4f}"])

    if args.pairs_csv:
        with open(args.pairs_csv, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["pair_type", "vhh_idx", "antigen_idx", "vhh_pdb", "antigen_pdb",
                        "n_replicas", "iptm_mean", "iptm_median", "iptm_max"])
            for label, d in (("real", per_pair_real), ("shuffled", per_pair_shuf)):
                for k, v in sorted(d.items()):
                    sv = sorted(v)
                    w.writerow([label, *k, len(v), f"{mean(sv):.4f}",
                                f"{quantile(sv, 0.5):.4f}", f"{sv[-1]:.4f}"])

    print(f"# real pairs {len(per_pair_real)}, shuffled pairs {len(per_pair_shuf)},"
          f" replicas {len(vals):,}, skipped {skipped}", file=sys.stderr)
    print(f"# wrote {args.csv}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
