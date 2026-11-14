"""
Oslo AI x Antibodies hackathon - CDRH3 design filter.

Reads AntiFold FASTA output, extracts CDRH3 loops, applies the optimisation-track
hard filters, and writes a ranked, diverse set of VH sequences for folding.

Still TODO (hard filter, blocking): novelty condition 2 - every design must be
>2 point mutations from all 524,346 Chinery CDRH3 variants.
"""

import glob, re, os
from collections import Counter

WT = "WGGDGFYAMDY"                 # trastuzumab CDRH3, 11 residues
CONTACTS = [1, 5, 6, 7]            # PyMOL 1N8Z: W99, G103, F104, Y105 contact HER2 (chain C, 4.5 A)
INVARIANT = [3, 10, 11]            # 100% conserved across AntiFold samples
MIN_MUT = 4                        # <70% identity on an 11-mer
N_DIVERSE = 15
N_W99 = 5

# trastuzumab VH framework flanking CDRH3
PRE = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRY"
       "ADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSR")
POST = "WGQGTLVTVSS"


def cdrh3(seq):
    """Pull the CDRH3 loop out of a full-length heavy chain."""
    m = re.search(r"YYC(.*?)WG[QKR]G", seq)
    if not m:
        return None
    loop = m.group(1)
    return loop[2:] if len(loop) > 3 else None   # drop the SR/AR prefix


def ham(a, b):
    """Hamming distance, or None if lengths differ."""
    return sum(1 for x, y in zip(a, b) if x != y) if len(a) == len(b) else None


# ---------- 1. parse every AntiFold batch ----------
seen, keep = set(), []
for f in sorted(glob.glob("results/*/*.txt")):
    if os.path.getsize(f) < 10000:
        continue
    txt = open(f, errors="ignore").read()
    if not txt.startswith(">"):
        continue
    batch = f.split("/")[1]
    for block in txt.split(">")[1:]:
        lines = block.strip().split("\n")
        seq = "".join(lines[1:]).strip()
        h3 = cdrh3(seq)
        if not h3 or h3 in seen:
            continue
        seen.add(h3)
        keep.append((batch, h3, ham(h3, WT), len(h3)))

print(f"unique CDRH3 parsed: {len(keep)}")
print("mutation counts:", dict(sorted(Counter(k[2] for k in keep).items(),
                                      key=lambda x: (x[0] is None, x[0]))))

# ---------- 2. hard filters ----------
pool = [k for k in keep if k[2] and k[2] >= MIN_MUT and "C" not in k[1]]
print(f"pool after >={MIN_MUT} mutations and no cysteine: {len(pool)}")

# ---------- 3. per-position tolerance ----------
print("\nper-position residue distribution (pool):")
for i in range(11):
    tag = []
    if i + 1 in CONTACTS:
        tag.append("HER2-contact")
    if i + 1 in INVARIANT:
        tag.append("invariant")
    print(f"  {i+1:2d} {WT[i]}  {Counter(h[1][i] for h in pool).most_common(5)}"
          f"  {' '.join(tag)}")

n_keep_w = sum(1 for k in pool if k[1][0] == "W")
print(f"\ndesigns preserving W99 (a HER2 contact): {n_keep_w} / {len(pool)}")

# ---------- 4. hypothesis A: model-chosen, maximally diverse ----------
pool.sort(key=lambda x: (-x[2], x[1]))
sel = [pool[0]]
while len(sel) < N_DIVERSE:
    rest = [c for c in pool if c not in sel]
    if not rest:
        break
    sel.append(max(rest, key=lambda c: min(ham(c[1], s[1]) for s in sel)))

# ---------- 5. hypothesis B: restore the W99 contact ----------
hyp_b, seen_b = [], set(k[1] for k in sel)
for k in sel:
    cand = "W" + k[1][1:]
    if cand in seen_b or ham(cand, WT) < MIN_MUT or "C" in cand:
        continue
    seen_b.add(cand)
    hyp_b.append((k[0] + "+W99", cand, ham(cand, WT), 11))
    if len(hyp_b) == N_W99:
        break

final = sel + hyp_b

# ---------- 6. report ----------
print(f"\n--- final set: {len(sel)} diverse + {len(hyp_b)} W99-restored = {len(final)} ---")
for b, h3, d, _ in final:
    broken = [i for i in INVARIANT if h3[i - 1] != WT[i - 1]]
    flag = f"breaks invariant {broken}" if broken else ""
    print(f"  {h3}  mut={d}  {b:24s} {flag}")

# ---------- 7. write VH sequences for folding ----------
os.makedirs("designs", exist_ok=True)
with open("designs/designs.fasta", "w") as fh:
    fh.write(f">wt_trastuzumab_VH\n{PRE + WT + POST}\n")
    for i, (b, h3, d, _) in enumerate(final):
        fh.write(f">d{i:02d}_{h3}_mut{d}\n{PRE + h3 + POST}\n")

with open("designs/designs.csv", "w") as fh:
    fh.write("id,cdrh3,mutations,batch,breaks_invariant\n")
    for i, (b, h3, d, _) in enumerate(final):
        broken = ";".join(str(x) for x in INVARIANT if h3[x - 1] != WT[x - 1])
        fh.write(f"d{i:02d},{h3},{d},{b},{broken}\n")

print(f"\nwrote designs/designs.fasta ({len(final)} designs + WT, "
      f"VH length {len(PRE + WT + POST)} aa)")
print("wrote designs/designs.csv")
print("\nNOT YET CHECKED: novelty condition 2 (>2 mutations from all 524k Chinery variants)")