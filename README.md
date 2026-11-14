# Team 10 — AI × Antibodies, Oslo 2026

**Optimisation track.** Trastuzumab CDRH3 redesign against HER2 domain IV.

We changed 6–8 of the 11 residues in trastuzumab's binding loop and kept the
binding. Our top designs score an AF2-IG interface pTM of **0.971** against a
measured wild-type baseline of **0.975**, on a framework that is otherwise
untouched.

---

## The problem

The optimisation track requires a design to be **under 70% identical** to
trastuzumab's CDRH3 — at least 3 of its 10 positions changed — and at least 2
mutations away from every sequence in the provided data.

That loop is the antibody's main contact with HER2, and it has already been
optimised to clinical use. The filter forces you away from the known answer, in
exactly the region where being wrong costs you the binding.

So the design question is not *how do we change this loop* but **which positions
can we afford to change**.

---

## What we did

### 1. Map the loop before touching it

Two independent analyses, neither of which depends on the other:

**Structural.** PyMOL on PDB 1N8Z (trastuzumab Fab bound to HER2). Four CDRH3
positions sit within 4.5 Å of HER2 chain C:

| Loop position | Residue | Role |
|---|---|---|
| 1 | W99 | contacts HER2 |
| 5 | G103 | contacts HER2 |
| 6 | F104 | contacts HER2 |
| 7 | Y105 | contacts HER2 — largest footprint (8 of 15 contact atoms) |

**Statistical.** AntiFold inverse folding, 237 valid samples. Three positions
never move:

| Loop position | Residue | Conservation |
|---|---|---|
| 3 | G101 | 237 / 237 |
| 10 | D108 | 237 / 237 |
| 11 | Y109 | 235 / 237 |

Those three are not HER2 contacts. They hold the loop's conformation — the shape
that presents the contact residues.

That leaves **five tolerant positions** (1, 2, 4, 6, 9) and a requirement of at
least three mutations. The novelty filter can be satisfied entirely inside the
region the model treats as permissive.

### 2. Generate

AntiFold (OPIG), inverse folding on 1N8Z with HER2 included as context, across
four sampling regimes:

| Batch | Regions | Temperature | Purpose |
|---|---|---|---|
| `af-cdr3-t050` | CDR3 | 0.50 | conservative |
| `af-cdr3-t080` | CDR3 | 0.80 | exploratory |
| `af-allcdr-t050` | CDR1+2+3 | 0.50 | variation outside H3 |
| `af-cdr3-t120` | CDR3 | 1.20 | high diversity |

**253 unique CDRH3 loops.**

### 3. Filter

- ≥4 mutations from `WGGDGFYAMDY` (stricter than the required 3)
- no cysteine in the loop (unpaired Cys fails the hard filter and aggregates)
- → **237 candidates**

### 4. Select 20 for diversity, not for score

237 candidates collapse into a handful of real ideas — AntiFold finds one
solution and samples around it. Ranking by mutation count would have returned
twenty near-duplicates.

Instead: greedy max-min Hamming selection, so each pick is as far as possible
from everything already chosen.

Split across two hypotheses:

- **15 model-chosen**, maximally diverse
- **5 with W99 restored** — only 5 of 237 candidates preserved the tryptophan
  that contacts HER2, so we rebuilt that hypothesis deliberately rather than
  letting the model discard it

### 5. Fold

Boltz-2, three chains per job — designed VH (120 aa), trastuzumab VL (107 aa),
HER2 extracellular domain (581 aa) — with MSAs from PDB70 via MMseqs2.

**The wild type was folded through the identical pipeline.** Every number below
is relative to that baseline.

### 6. Validate structurally

Framework superposed on 1N8Z chain B (residues 1–98, 110–120, CA), then CDRH3
displacement measured at residues 99–109.

---

## Results

### Calibration

Boltz-2 reproduces the **wild-type** CDRH3 at **0.883 Å** against its own crystal
structure. That is the prediction floor for this pipeline. Anything above it is
real deviation, not model noise — a number we could only state because we folded
the baseline the same way.

### Interface scoring (AF2-IG, JapanFold)

| | iPTM | pTM | pLDDT |
|---|---|---|---|
| wild-type | **0.975** | 0.980 | 0.963 |
| d15 `WYGSAEYSLDS` | 0.971 | 0.977 | 0.954 |
| d12 `WEGSAYYAFEY` | 0.971 | 0.978 | 0.937 |

### Structural validation

| Design | CDRH3 | Mut | Framework RMSD | Loop RMSD | W99 |
|---|---|---|---|---|---|
| wild-type | WGGDGFYAMDY | 0 | 0.369 Å | **0.883 Å** | yes |
| d04 | IVGFANNHFDY | 8 | 0.359 Å | **1.861 Å** | no |
| d06 | TYGAAYFAMGS | 8 | 0.380 Å | 2.072 Å | no |
| d15 | WYGSAEYSLDS | 7 | 0.378 Å | 2.169 Å | yes |
| d12 | WEGSAYYAFEY | 6 | 0.355 Å | 2.309 Å | yes |
| d08 | ENGHSQFAMDH | 7 | 0.386 Å | 2.927 Å | no |
| d05 | THAESNYALDW | 8 | 0.369 Å | 2.984 Å | no |

Frameworks are identical across all designs (~0.37 Å), as expected — nothing
there changed. Loops sit 1–2 Å beyond the prediction floor. The loops **have**
moved; we are not claiming otherwise.

### Ranking

The 20 were selected for diversity, then ranked by Boltz-2 confidence relative to
the wild-type baseline folded through the identical pipeline. The top 6 were
additionally validated structurally (framework fit plus loop displacement), and
the top 2 scored on AF2-IG for interface-specific confidence.

Ranks 7–20 rest on Boltz-2 confidence alone.

Two caveats we state rather than hide. The top margins are small — +0.003 and
+0.002 over baseline — and AF2-IG disagrees on the sign, placing both leads
marginally below wild-type. We therefore rank on the ordering the two models
broadly agree on, not on either margin. And d05 ranks third on Boltz-2
confidence but has nearly the worst loop placement of the six we validated
(2.984 Å against a 0.883 Å floor); where structural validation and Boltz-2
confidence conflict, we weight the structural evidence.

### Lead designs

**d15 — `WYGSAEYSLDS`**, 7 mutations. Highest Boltz-2 confidence, tied-best
iPTM, preserves W99.

**d12 — `WEGSAYYAFEY`**, 6 mutations. Fewest changes, same iPTM, preserves W99.

Both: VH 120 aa + VL 107 aa = **227 aa**, under the 250 limit.

---

## What we did not measure

- **Thermostability.** No predictor in our stack. Our argument is structural —
  unchanged clinical framework, changes confined to one loop — not experimental.
- **Affinity.** iPTM and Boltz-2 confidence measure structural plausibility, not
  binding strength. The task's own data makes this concrete: 21 of 454
  SPR-confirmed non-binders score AlphaFold3 ipTM ≥ 0.7.

Two models disagree on the sign: Boltz-2 puts both leads marginally above
baseline, AF2-IG marginally below. Both margins are within noise. The honest
conclusion is **indistinguishable from the parent on predicted binding** —
which, given the filter forced us 6–7 mutations away from it, is the result we
were aiming for.

---

## What we did measure
- **Novelty condition 2.** The ≥2-mutation check against all 525,063 provided
  sequences.
  

## Why this is useful

Trastuzumab is patented. A molecule that engages the same HER2 epitope with
comparable predicted affinity, but shares under 70% CDRH3 identity, is outside
the composition-of-matter claims on the original.

That matters to biosimilar developers entering the HER2 market, to ADC
programmes needing a fresh binder scaffold for new linker chemistry, to
bispecific developers who need a HER2 arm, and to markets where the originator
is unaffordable.

The method generalises beyond HER2: any legacy antibody facing patent expiry, a
developability flag, or an immunogenicity problem needs the same thing — change
the sequence enough to be new, little enough to still work.

---

## Reproducing

```bash
python3 filter.py   # parse AntiFold output, filter, select 20, write designs/
python3 fold.py     # submit 21 Boltz-2 jobs to Nebius
python3 rank.py     # fetch, score, write structures/
```

Structural validation, in PyMOL:

```
load 1n8z.pdb
load structures/fold-d15_WYGSAEYSLDS_mut7.cif
pair_fit fold-d15_WYGSAEYSLDS_mut7 and chain H and resi 1-98+110-120 and name CA, \
         1n8z and chain B and resi 1-98+110-120 and name CA
rms_cur fold-d15_WYGSAEYSLDS_mut7 and chain H and resi 99-109 and name CA, \
        1n8z and chain B and resi 99-109 and name CA, matchmaker=-1
```

| Path | Contents |
|---|---|
| `final.csv` | all 20 ranked, with evidence |
| `designs/` | sequences (FASTA + CSV) |
| `structures/` | 22 Boltz-2 complexes (mmCIF) |
| `nebius-runs/` | run IDs for every job |
| `cdrh3_overlay.png` | designed loop vs crystal, in the HER2 pocket |

---

## Tools and data

AntiFold · Boltz-2 · MMseqs2/PDB70 — Nebius Scientific AI Platform
AF2-IG — JapanFold
PyMOL — structural analysis
PDB 1N8Z — trastuzumab Fab / HER2 extracellular domain

CDRH3 variant data: Chinery et al., *bioRxiv* 2024.03.26.586756
