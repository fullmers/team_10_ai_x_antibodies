#!/usr/bin/env python3
"""Check CDRH3 designs against the HER2 challenge novelty filter.

For each input CDRH3 loop, report:
  * identity % to trastuzumab's CDRH3 (WGGDGFYAMD)
  * the minimum Hamming distance to any CDRH3 in the provided challenge data
    (the FACS table + the SPR table, 525k unique sequences)
  * pass / fail against the novelty rule

The rule as stated in data/docs/12-her2.md ("Challenge filter"):

    A design passes if its CDRH3 differs at 2 or more of its 10 positions from
    every sequence in the provided data AND is at most 70% identical to
    trastuzumab's CDRH3 (at least 3 of 10 positions changed).

That is the default here (--rule docs). A stricter reading -- identity strictly
under 70% and strictly more than 2 mutations from every provided sequence -- is
available as --rule strict. Both verdicts are printed in every run so the
difference is visible; --rule only selects which one the PASS/FAIL column and
the exit code follow.

Usage:
  novelty_check.py LOOP [LOOP ...]
  novelty_check.py --loops designs.txt
  cat designs.txt | novelty_check.py

Input files may be one loop per line, optionally followed by whitespace and a
label (e.g. "DLSLSDFSFD  her2_cdrh3_0"). Blank lines and #-comments are ignored.
"""

from __future__ import annotations

import argparse
import csv
import io
import os
import sys
import zipfile

TRASTUZUMAB_CDRH3 = "WGGDGFYAMD"

DEFAULT_DATA_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "data",
    "hackathon-data",
    "her2-trastuzumab",
)
FACS_BASENAME = "tz_her2_facs_3classes_524346.csv"
SPR_BASENAME = "trastuzumab_her2_spr_all_768.csv"

# Documented rule: identity <= 70%, min distance >= 2.
# Strict rule:     identity <  70%, min distance >  2.
RULES = {
    "docs": {"max_identity": 70.0, "identity_inclusive": True,
             "min_distance": 2, "distance_inclusive": True},
    "strict": {"max_identity": 70.0, "identity_inclusive": False,
               "min_distance": 2, "distance_inclusive": False},
}


def hamming(a: str, b: str) -> int:
    """Hamming distance between two equal-length strings."""
    return sum(x != y for x, y in zip(a, b))


def min_hamming(query: str, pool: list[str], floor: int = 0) -> int:
    """Smallest Hamming distance from query to any sequence in pool.

    Uses the running best as an early-exit bound, and stops entirely once the
    distance cannot get any lower.
    """
    best = len(query)
    for seq in pool:
        d = 0
        for x, y in zip(query, seq):
            if x != y:
                d += 1
                if d >= best:
                    break
        else:
            best = d
            if best <= floor:
                return best
    return best


def _read_cdrh3_column(handle) -> list[str]:
    reader = csv.DictReader(handle)
    if reader.fieldnames is None or "cdrh3" not in reader.fieldnames:
        raise SystemExit(f"no 'cdrh3' column in {getattr(handle, 'name', handle)}")
    return [row["cdrh3"].strip().upper() for row in reader if row.get("cdrh3")]


def load_provided_data(data_dir: str) -> tuple[list[str], dict]:
    """Load every CDRH3 from the two provided challenge tables.

    Returns (unique sequences, stats). The FACS table is read from the plain
    CSV if it has been unzipped, otherwise straight out of the .zip.
    """
    seqs: list[str] = []
    stats: dict = {}

    facs_csv = os.path.join(data_dir, FACS_BASENAME)
    facs_zip = facs_csv + ".zip"
    if os.path.exists(facs_csv):
        with open(facs_csv, newline="") as fh:
            facs = _read_cdrh3_column(fh)
        stats["facs_source"] = facs_csv
    elif os.path.exists(facs_zip):
        with zipfile.ZipFile(facs_zip) as zf:
            name = next(n for n in zf.namelist()
                        if n.endswith(".csv") and not n.startswith("__MACOSX"))
            with zf.open(name) as raw:
                facs = _read_cdrh3_column(io.TextIOWrapper(raw, "utf-8", newline=""))
        stats["facs_source"] = facs_zip
    else:
        raise SystemExit(
            f"FACS table not found in {data_dir}.\n"
            "Fetch it with:\n"
            "  aws s3 cp s3://apheris-oslo-hackathon/data/her2-trastuzumab/ "
            f"{data_dir}/ --recursive"
        )
    stats["facs_rows"] = len(facs)
    seqs.extend(facs)

    spr_csv = os.path.join(data_dir, SPR_BASENAME)
    if os.path.exists(spr_csv):
        with open(spr_csv, newline="") as fh:
            spr = _read_cdrh3_column(fh)
        stats["spr_source"] = spr_csv
        stats["spr_rows"] = len(spr)
        seqs.extend(spr)
    else:
        stats["spr_source"] = None
        stats["spr_rows"] = 0

    unique = sorted(set(seqs))
    stats["unique"] = len(unique)
    return unique, stats


def parse_loops(args_loops: list[str], loops_file: str | None) -> list[tuple[str, str]]:
    """Return a list of (loop, label). Reads stdin when nothing else is given."""
    out: list[tuple[str, str]] = []

    def add_line(line: str) -> None:
        line = line.split("#", 1)[0].strip()
        if not line:
            return
        parts = line.split()
        loop = parts[0].strip().upper()
        label = parts[1] if len(parts) > 1 else loop
        out.append((loop, label))

    for loop in args_loops:
        add_line(loop)
    if loops_file:
        with open(loops_file) as fh:
            for line in fh:
                add_line(line)
    if not out and not sys.stdin.isatty():
        for line in sys.stdin:
            add_line(line)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Novelty check for HER2 CDRH3 designs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument("loops", nargs="*", help="CDRH3 loops, e.g. DLSLSDFSFD")
    ap.add_argument("--loops", dest="loops_file", metavar="FILE",
                    help="file of loops, one per line, optional label after whitespace")
    ap.add_argument("--data-dir", default=DEFAULT_DATA_DIR,
                    help=f"folder holding the two provided tables (default: {DEFAULT_DATA_DIR})")
    ap.add_argument("--rule", choices=sorted(RULES), default="docs",
                    help="which reading of the filter drives PASS/FAIL (default: docs)")
    ap.add_argument("--reference", default=TRASTUZUMAB_CDRH3,
                    help=f"reference CDRH3 (default: {TRASTUZUMAB_CDRH3})")
    ap.add_argument("--csv", metavar="FILE", help="also write results as CSV")
    ap.add_argument("--fail-on-reject", action="store_true",
                    help="exit 1 if any design fails the rule (default: always exit 0)")
    args = ap.parse_args()

    designs = parse_loops(args.loops, args.loops_file)
    if not designs:
        ap.error("no CDRH3 loops given (pass them as arguments, with --loops, or on stdin)")

    ref = args.reference.strip().upper()
    rule = RULES[args.rule]

    pool_all, stats = load_provided_data(args.data_dir)
    ref_len = len(ref)
    pool = [s for s in pool_all if len(s) == ref_len]
    skipped = len(pool_all) - len(pool)

    print(f"# reference CDRH3       : {ref} ({ref_len} residues)")
    print(f"# FACS table            : {stats['facs_source']} ({stats['facs_rows']:,} rows)")
    print(f"# SPR table             : {stats['spr_source']} ({stats['spr_rows']:,} rows)")
    print(f"# unique provided CDRH3 : {stats['unique']:,}"
          f"  ({len(pool):,} of length {ref_len}, {skipped:,} other lengths skipped)")
    print(f"# rule in force         : {args.rule} -- identity "
          f"{'<=' if rule['identity_inclusive'] else '<'} {rule['max_identity']:.0f}% "
          f"AND min distance {'>=' if rule['distance_inclusive'] else '>'} "
          f"{rule['min_distance'] if rule['distance_inclusive'] else rule['min_distance']}")
    print()

    header = ("label", "cdrh3", "mut_vs_tz", "identity_%", "min_dist",
              "pass_docs", "pass_strict", "verdict")
    rows = []
    for loop, label in designs:
        if len(loop) != ref_len:
            rows.append((label, loop, "-", "-", "-", "-", "-", "SKIP (length)"))
            continue
        mut = hamming(loop, ref)
        identity = (ref_len - mut) / ref_len * 100.0
        dist = min_hamming(loop, pool)

        ok_id_docs = identity <= RULES["docs"]["max_identity"]
        ok_id_strict = identity < RULES["strict"]["max_identity"]
        ok_d_docs = dist >= RULES["docs"]["min_distance"]
        ok_d_strict = dist > RULES["strict"]["min_distance"]
        pass_docs = ok_id_docs and ok_d_docs
        pass_strict = ok_id_strict and ok_d_strict

        chosen = pass_docs if args.rule == "docs" else pass_strict
        if chosen:
            verdict = "PASS"
        else:
            why = []
            ok_id = ok_id_docs if args.rule == "docs" else ok_id_strict
            ok_d = ok_d_docs if args.rule == "docs" else ok_d_strict
            if not ok_id:
                why.append(f"identity {identity:.0f}%")
            if not ok_d:
                why.append("in data" if dist == 0 else f"dist {dist}")
            verdict = "FAIL (" + ", ".join(why) + ")"

        rows.append((label, loop, str(mut), f"{identity:.0f}", str(dist),
                     "yes" if pass_docs else "no",
                     "yes" if pass_strict else "no", verdict))

    widths = [max(len(str(r[i])) for r in (header,) + tuple(rows)) for i in range(len(header))]
    fmt = "  ".join(f"{{:<{w}}}" for w in widths)
    print(fmt.format(*header))
    print(fmt.format(*("-" * w for w in widths)))
    for r in rows:
        print(fmt.format(*r))

    scored = [r for r in rows if r[7] != "SKIP (length)"]
    n_docs = sum(1 for r in scored if r[5] == "yes")
    n_strict = sum(1 for r in scored if r[6] == "yes")
    print()
    print(f"# {n_docs}/{len(scored)} pass the documented rule; "
          f"{n_strict}/{len(scored)} pass the strict rule")

    if args.csv:
        with open(args.csv, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(header)
            w.writerows(rows)
        print(f"# wrote {args.csv}")

    if args.fail_on_reject:
        n_failed = sum(1 for r in scored
                       if r[5 if args.rule == "docs" else 6] == "no")
        return 1 if n_failed else 0
    return 0


if __name__ == "__main__":
    sys.exit(main())
