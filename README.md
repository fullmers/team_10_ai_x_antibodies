# Team 10 — AI × Antibodies, Oslo 2026

HER2 binder design. **30 submitted designs: 28 optimisation, 2 de novo.**

Every optimisation design is trastuzumab's framework with a redesigned CDRH3.
All 30 were checked against the full provided dataset — 525,042 unique CDRH3
sequences — and pass both novelty conditions. The top 16 were additionally
screened for cross-reactivity against EGFR, HER3 and HER4.

Submission: [`team-10.csv`](team-10.csv)

---

## The problem

The optimisation track requires a design to be **under 70% identical** to
trastuzumab's CDRH3 — at least 3 of its 10 positions changed — and at least 2
mutations away from every sequence in the provided data.

That loop is the antibody's main contact with HER2, and it has already been
optimised to clinical use. The filter forces you away from the known answer, in
exactly the region where being wrong costs you the binding.

The design question is therefore not *how do we change this loop* but **which
positions can we afford to change, and what does the constraint cost us**.

---

## What we did

### 1. Map the loop before touching it

Two independent analyses, neither depending on the other.

**Structural.** PyMOL on PDB 1N8Z (trastuzumab Fab bound to HER2). Four CDRH3
positions sit within 4.5 Å of HER2 chain C:

| Loop position | Residue | Note |
|---|---|---|
| 1 | W99 | contacts HER2 |
| 5 | G103 | contacts HER2 |
| 6 | F104 | contacts HER2 |
| 7 | Y105 | contacts HER2 — 8 of 15 contact atoms |

**Statistical.** 237 valid AntiFold samples. Three positions never move:

| Loop position | Residue | Conservation |
|---|---|---|
| 3 | G101 | 237 / 237 |
| 10 | D108 | 237 / 237 |
| 11 | Y109 | 235 / 237 |

Those three are *not* HER2 contacts. They hold the loop's conformation — the
shape that presents the contact residues.

Five positions (1, 2, 4, 6, 9) are tolerant, and the filter needs three
mutations. **The novelty requirement can be satisfied entirely inside the region
the model treats as permissive.**

### 2. Generate with two independent models

| Model | Training | Output |
|---|---|---|
| AntiFold (OPIG) | antibody structures; never trained on binding | 253 unique CDRH3 loops, four sampling temperatures |
| RFD4-Proteina `her2-novel-2026-11-09` | Apheris finetune on 524,346 measured HER2 binders | 100 designs, top 20 refolded and scored |

Using one generator that has seen binding data and one that has not is
deliberate — it is what lets us separate *model limitation* from *constraint
limitation* in section 5.

### 3. Filter

Implemented from the published rules in `scripts/novelty.py`:

- ≥3 of 10 CDRH3 positions changed from trastuzumab
- ≥2 mutations from **every** one of the 525,042 provided CDRH3 sequences
  (FACS screen + SPR panel combined)
- no cysteine in the loop; ≤250 residues; even cysteine count

All 30 submitted designs pass. Trastuzumab itself fails, as it must — it is in
the provided data. That is the control proving the check works.

We also flag designs that alter the framework tyrosine immediately following
CDRH3, conserved across all 525,042 provided sequences and held fixed in the
reference design setup. Those were deprioritised: the change sits outside the
loop the filter measures and carries uncharacterised risk.

### 4. Score, three independent ways

**Structural confidence.** Boltz-2 (Nebius), three chains — designed VH, VL,
HER2 ECD — with PDB70 MSAs. AF2-IG (JapanFold) for interface pTM.

**Interface quality.** RFD4-Proteina `refold-score` on template-guided Boltz-2
folds: interface PAE, binder pLDDT, ipTM. Best design 0.328 Å interface PAE
against the Apheris antibody pass criterion of ≤2.0 Å.

**Cross-reactivity.** Top candidates screened against EGFR, HER3 and HER4
(`candidates_merged.csv`). HER2 sits in a receptor family with high structural
homology; a binder that also engages EGFR or HER3 is a liability regardless of
how tightly it binds the target. Selectivity ranges 0.491 to 0.277 across the
panel. **The 16 panel-screened designs are ranked by selectivity, ahead of
designs scored on ipTM alone.**

### 5. Calibrate everything against a baseline

Wild-type trastuzumab was run through every pipeline alongside the designs.

| Measurement | Wild type | Why it matters |
|---|---|---|
| Boltz-2 CDRH3 refold | **0.883 Å** vs 1N8Z | our prediction floor |
| AF2-IG iPTM | **0.975** | designs land at 0.971 — comparable |
| Measured K_D (SPR) | **0.152 nM**, ΔG −13.41 | the affinity bar |

Apheris advise against ranking on CDRH3 refold agreement, since Protenix does
not reproduce trastuzumab's own loop. We tested whether that caveat applied to
our pipeline rather than assuming either way: Boltz-2 reproduces the wild-type
loop at 0.883 Å, so our loop displacements are measured against a known floor.

---

## Results

### Lead designs

| Rank | CDRH3 | Changed | Selectivity | Worst off-target | Basis |
|---|---|---|---|---|---|
| t10-01 | WHSLDEYYFDY | 7/10 | **0.491** | 0.358 | panel |
| t10-02 | SLSLGNYAFDY | 6/10 | 0.439 | **0.288** | panel |
| t10-03 | WHSLYEYYFDY | 7/10 | 0.422 | 0.372 | panel |
| t10-17 | SLSLSNYAFDY | 7/10 | not screened | — | ipTM 0.926 |

`SLSLSNYAFDY` has the highest ipTM in the entire set and still ranks 17th,
because it has no cross-reactivity data. An untested design is not ranked above
one shown to be selective.

### De novo track

Two VHHs, 120 aa, interface PAE 1.64 and 1.70 Å, ipTM 0.887 and 0.888.

### The constraint finding

This is the result we consider most useful, and it is a negative one.

Of 190 SPR-measured designs, **every sub-nanomolar binder sits within 1–4
mutations of trastuzumab**, and nearly all retain the `WGG` loop start. The
novelty filter requires ≥3 mutations. That leaves a one-residue window.

We ran two generators across that boundary. A random-forest ΔG model fitted to
the 190 measured K_D values — which reproduces wild type at −13.69 kcal/mol
against a measured −13.41 — predicts:

| Source | Best predicted K_D |
|---|---|
| AntiFold | 18.2 nM |
| HER2-finetuned RFD4-Proteina | 46.6 nM |
| Trastuzumab (measured) | 0.152 nM |

**Both generators hit the same ceiling**, including the one trained on 524,000
measured binders. The bound is the novelty constraint, not the model. Any team
reporting a predicted sub-nanomolar design under this filter should be asked
how.

### Epitope: a check that failed usefully

Our first de novo VHH was explicitly targeted at HER2 domain IV. We superposed
the predicted complex onto 1N8Z (HER2 RMSD 1.103 Å over 3,477 atoms) and
compared contact sets:

- Trastuzumab epitope: residues 557, 558, 560, 561, 570, 571, 572, 573, 579,
  591, 593
- VHH footprint: residues 90–196 — **domain I/II**
- **Shared: 0 of 12**

The design intent did not survive generation. Epitope preservation is judged by
competitive binding against trastuzumab, so that molecule was excluded. The two
submitted VHHs carry the same caveat in their row descriptions, flagged rather
than assumed away.

---

## What we did not measure

- **Thermostability.** No predictor in our stack. Our argument is structural —
  unchanged clinical framework, changes confined to one loop — not experimental.
- **Affinity, directly.** ipTM and interface PAE measure structural
  plausibility. The task's own data makes the gap concrete: 21 of 454
  SPR-confirmed non-binders score AlphaFold3 ipTM ≥ 0.7. Our ΔG predictions
  extrapolate beyond their training distribution and are reported as such.
- **De novo epitope.** Unconfirmed for both submitted VHHs.

Boltz-2 and AF2-IG disagree on the sign of our designs' difference from wild
type; both margins are within noise. The honest conclusion is **indistinguishable
from the parent on predicted binding** — which, given the filter forced us 6–8
mutations away, is the outcome we were aiming for.

---

## Why this is useful

Trastuzumab is patented. A molecule that engages the same HER2 epitope with
comparable predicted affinity, but under 70% CDRH3 identity, falls outside the
composition-of-matter claims on the original — and our panel data says it does
so without picking up EGFR or HER3 cross-reactivity.

That matters to biosimilar developers entering the HER2 market, to ADC
programmes needing a fresh binder scaffold for new linker chemistry, and to
bispecific developers who need a selective HER2 arm.

The method generalises: any legacy antibody facing patent expiry, a
developability flag, or an immunogenicity problem needs the same thing — change
the sequence enough to be new, little enough to still work, and know in advance
what that trade costs.

---

## Repository

| Path | Contents |
|---|---|
| `team-10.csv` | **the submission** — 30 designs, both tracks |
| `candidates_merged.csv` | off-target panel: HER2 vs EGFR / HER3 / HER4 |
| `ranked_designs.csv` | HER2-finetuned RFD4-Proteina designs, scored |
| `antifold-pipeline/` | generation, filtering, folding, ranking |
| `antifold-pipeline/structures/` | 22 Boltz-2 complexes (mmCIF) |
| `antifold-pipeline/nebius-runs/` | run IDs for every job |
| `scripts/novelty.py` | novelty check against 525,042 provided sequences |
| `scripts/predict_kd.py` | ΔG model fitted to the 190 SPR measurements |
| `scripts/make_submission.py` | builds the final CSV |
| `rfd4/` | Apheris RFD4-Proteina documentation (organiser-provided) |

Every number above traces to a file here: a run ID, a structure, or a row in the
provided data.

### Reproducing

```bash
python3 antifold-pipeline/filter.py   # parse, filter, select
python3 antifold-pipeline/fold.py     # submit Boltz-2 jobs
python3 antifold-pipeline/rank.py     # fetch, score, write structures
python3 scripts/novelty.py            # novelty vs all provided sequences
python3 scripts/predict_kd.py         # predicted dG / K_D
python3 scripts/make_submission.py    # build team-10.csv
```

Structural validation, in PyMOL:

```
load antifold-pipeline/1n8z.pdb
load antifold-pipeline/structures/fold-d15_WYGSAEYSLDS_mut7.cif
pair_fit fold-d15_WYGSAEYSLDS_mut7 and chain H and resi 1-98+110-120 and name CA, \
         1n8z and chain B and resi 1-98+110-120 and name CA
rms_cur fold-d15_WYGSAEYSLDS_mut7 and chain H and resi 99-109 and name CA, \
        1n8z and chain B and resi 99-109 and name CA, matchmaker=-1
```

---

## Tools and data

AntiFold · Boltz-2 · MMseqs2 / PDB70 — Nebius Scientific AI Platform
RFD4-Proteina (`her2-novel-2026-11-09`) · refold-score — Apheris Foundry
AF2-IG — JapanFold
PyMOL — structural analysis
PDB 1N8Z — trastuzumab Fab / HER2 extracellular domain

CDRH3 variant data: Chinery et al., *bioRxiv* 2024.03.26.586756