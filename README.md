# HER2 CDRH3 redesign — team 10 pipeline

Optimisation track: redesign trastuzumab's CDRH3 loop against HER2, plus a de novo VHH
side track. Everything here was built during the hackathon.

## The finding

The design this pipeline ranked 5th of 60 on interface confidence (`WGPPGFYAFD`) sits at
Hamming distance **zero** from the Chinery screen. It is in the training data, labelled a
high binder, and the model sampled that exact loop **five separate times** out of 200.
Two more of the top ten are 80% identical to trastuzumab at distance 1.

The model is most confident about sequences it has already seen. Ranking on `min_ipae`
alone would have submitted a known, ineligible sequence. Hence: filter before ranking.

## Calibration

We parsed 735,172 replica predictions from the Champloo VHH-antigen benchmark (121
confirmed binders, 14,601 shuffled non-binders) to find a confidence threshold.

| population | median ipTM | above 0.7 |
| --- | --- | --- |
| non-binders | 0.172 | 0.13% |
| real binders | 0.187 | 5.8% |

The distributions overlap. A threshold excluding 95% of non-binders keeps only **14% of
real binders**. Our designs all score above 0.73 on Boltz-2, where only 5.8% of confirmed
binders reach 0.7 on Chai-1 — evidence the scales differ, not that our designs are better.
We used no borrowed threshold and ranked relatively within one model.

Note: the task brief describes 454 negatives. The archive holds 121 real and 14,601
shuffled. We report what we found.

## Pipeline

    cd data/docs
    # design: RFD4-Proteina her2-novel-2026-11-09, N=200, seed 2
    # (run as container directly; run-her2.sh hardcodes --seed 1)
    N=200 TOPN=60 OUT=runs/her2-main bash Scripts/run-her2.sh prepare
    N=200 TOPN=60 OUT=runs/her2-main bash Scripts/run-her2.sh predict
    N=200 TOPN=60 OUT=runs/her2-main bash Scripts/run-her2.sh score

    python3 -I scripts/extract_cdrh3.py data/docs/runs/her2-main/design > scripts/designs200.txt
    python3 -I scripts/novelty_check.py --loops scripts/designs200.txt --rule strict --csv scripts/novelty200_strict.csv
    python3 -I scripts/developability.py --loops scripts/designs200.txt --subset scripts/strict_refolded.txt --csv scripts/developability.csv
    python3 -I scripts/stability_proxy.py
    python3 -I scripts/build_submission.py

Exact commands, image ids, weight checkpoints, seeds and timings are in `log.md`,
including the attempts that failed.

## Filters

**Novelty** — both readings of the rule against 525,042 provided sequences. 172/200 pass
the documented rule, 142/200 strict. Validated: our Hamming distance agrees with the FACS
table's own `edit_distance` column on 524,346/524,346 rows.

**Developability** — six liability motifs plus hydrophobic fraction and net charge,
scanned in framework context (`SR[loop]YW`), not on the bare loop. Three designs flagged
for hydrophobic runs turned out to be junction artifacts from the flanking tryptophan.
19 of 37 survivors clean.

**Stability proxy** — glycine/proline counts. Not a Tm prediction. All 37 designs have
fewer glycines than the parent; none has more, so this does not discriminate between
designs, only between designs and the parent.

## Selection

Slots 1-4 are hand-picked to differ by construction, not by score. The top four by
`min_ipae` would be near-identical variations on one binding mode.

1. `SLSLSNYAFD` — best score (0.3076), clean
2. `TRAAVGYAHD` — second family, neutral charge, distance 4
3. `DLSYYNHAFD` — largest novelty margin (distance 5, 20% identity)
4. `WSSAVGYAHD` — conservative hedge, 6 mutations, iptm ties the lead

Slots 2 and 4 are the same family: a deliberate paired bet, and their order is arbitrary
(0.3234 vs 0.3238). Slots 1 and 3 also cluster at d<=2, so the four span three families,
not four. Both cutoffs are in the CSV.

## Results

### Submitted ranking (20 designs)

| # | CDRH3 | design id | min_ipae | iptm | b_pLDDT | dist | id% | mut | hyd | chg | dev |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | SLSLSNYAFD | her2_cdrh3_140 | 0.3076 | 0.926 | 98.59 | 3 | 30 | 7 | 0.50 | -1 | clean |
| 2 | TRAAVGYAHD | her2_cdrh3_98 | 0.3234 | 0.903 | 98.32 | 4 | 30 | 7 | 0.50 | +0 | clean |
| 3 | DLSYYNHAFD | her2_cdrh3_57 | 0.3559 | 0.885 | 97.94 | 5 | 20 | 8 | 0.50 | -2 | clean |
| 4 | WSSAVGYAHD | her2_cdrh3_31 | 0.3238 | 0.926 | 98.16 | 3 | 40 | 6 | 0.50 | -1 | clean |
| 5 | SLYVNNYSYD | her2_cdrh3_39 | 0.3288 | 0.902 | 98.13 | 4 | 20 | 8 | 0.50 | -1 | FLAG N-x-S/T |
| 6 | DLFFPYYAFD | her2_cdrh3_27 | 0.3454 | 0.892 | 95.99 | 3 | 30 | 7 | 0.70 | -2 | FLAG hydrophobic |
| 7 | TTAAVGYAHD | her2_cdrh3_37 | 0.3481 | 0.879 | 98.29 | 4 | 30 | 7 | 0.50 | -1 | clean |
| 8 | YMSLYYYEFD | her2_cdrh3_184 | 0.3504 | 0.909 | 95.28 | 3 | 20 | 8 | 0.70 | -2 | FLAG hydrophobic |
| 9 | YMSLYYYAFD | her2_cdrh3_134 | 0.3578 | 0.882 | 96.62 | 3 | 30 | 7 | 0.80 | -1 | FLAG hydrophobic |
| 10 | DLSSSDHAFD | her2_cdrh3_192 | 0.3618 | 0.856 | 98.22 | 4 | 20 | 8 | 0.30 | -3 | FLAG charge |
| 11 | SLSLSDYAFD | her2_cdrh3_43 | 0.3646 | 0.879 | 97.52 | 3 | 30 | 7 | 0.50 | -2 | clean |
| 12 | DLSSSDFAFD | her2_cdrh3_147 | 0.3650 | 0.882 | 98.06 | 4 | 20 | 8 | 0.40 | -3 | FLAG charge |
| 13 | SLSLSNFAFD | her2_cdrh3_9 | 0.3663 | 0.836 | 98.60 | 4 | 20 | 8 | 0.50 | -1 | clean |
| 14 | WRSLSEYYFD | her2_cdrh3_148 | 0.3692 | 0.877 | 97.18 | 3 | 30 | 7 | 0.50 | -1 | clean |
| 15 | DLSLHNHAFD | her2_cdrh3_181 | 0.3726 | 0.878 | 97.41 | 4 | 20 | 8 | 0.40 | -2 | clean |
| 16 | WASLHDYYHD | her2_cdrh3_193 | 0.3735 | 0.893 | 95.63 | 3 | 30 | 7 | 0.50 | -2 | clean |
| 17 | YASSSDYYFD | her2_cdrh3_114 | 0.3785 | 0.877 | 97.20 | 3 | 20 | 8 | 0.50 | -2 | clean |
| 18 | WKSTTQYYFD | her2_cdrh3_22 | 0.3839 | 0.875 | 96.85 | 4 | 30 | 7 | 0.40 | +0 | clean |
| 19 | AARPSSRSGSDDWDY | rep117 (de novo VHH) | 1.6377 | 0.887 | 94.64 | n/a | n/a | n/a | — | — | clean |
| 20 | AARPWWKAGSDAYDY | rep106 (de novo VHH) | 1.7012 | 0.888 | 97.62 | n/a | n/a | n/a | — | — | FLAG N-x-S/T |

Parent: `WGGDGFYAMD`, hydrophobic fraction 0.50, net charge -2.
`min_ipae` lower is better; `iptm` and `b_pLDDT` higher; `dist` is minimum Hamming
distance to 525,042 measured sequences, higher is more novel.

**The spread is narrow.** `min_ipae` runs 0.308 to 0.384 across all 18 CDRH3 designs.
The scores do not separate candidates as cleanly as a ranking implies, which is why the
filters do the real work. De novo rows are a different molecule in a different system
and are not comparable with the CDRH3 rows.

### Excluded

| design | why |
|---|---|
| WGPPGFYAFD (her2_cdrh3_70) | Distance 0. Present verbatim in the FACS screen, labelled `high`. A known binder, sampled 5x |
| WGSPGFYAMD (her2_cdrh3_11, _28) | 80% identity to trastuzumab, distance 1 |
| YMNLYYYAFD, DLSLSEHAFD | Cut to meet the 20-design limit; both developability-flagged |

### Filter funnel

| stage | count |
|---|---|
| designs generated | 200 |
| distinct loops | 164 |
| refolded with Boltz-2 | 60 (52 unique) |
| pass strict novelty | 37 |
| clean developability | 19 |
| submitted | 18 CDRH3 + 2 VHH |

### Stability proxy

| | glycines | prolines | longest G/S run |
|---|---|---|---|
| parent WGGDGFYAMD | 3 | 0 | 2 |
| SLSLSNYAFD | 0 | 0 | 1 |
| TRAAVGYAHD | 1 | 0 | 1 |
| DLSYYNHAFD | 0 | 0 | 1 |
| WSSAVGYAHD | 1 | 0 | 2 |

Across all 37: 33 have zero glycines, 4 have one, none has more than the parent.

### De novo epitope check

| design | min_ipae | hotspots within 4.5 A | footprint on trastuzumab epitope | verdict |
|---|---|---|---|---|
| rep106 | 1.70 | 3/4 | 73% (11 of 15) | on epitope |
| rep117 | 1.87 | 3/4 | 79% (11 of 14) | on epitope |
| rep48 | ~3 | 2/4 | 100% (9 of 9) | on epitope |
| rank 4 | 5.54 | — | 53% | drifted |
| rank 5 | 5.97 | — | partial | drifted |

### De novo sequences

    rep117  EVQLVESGGGLVQPGGSLRLSCAASGRSFSDYHMGWFRQAPGKGRELVASISSSGSTWYPDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYCAARPSSRSGSDDWDYWGQGTQVTVS
    rep106  EVQLVESGGGLVQPGGSLRLSCAASGRTFTDYAMGWFRQAPGKGRELVANITTSGSTWYPDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYCAARPWWKAGSDAYDYWGQGTQVTVS

## De novo

A VHH built on the caplacizumab framework against HER2 domain IV. The input had to be
constructed by grafting (`scripts/build_denovo_input.py`); validation rejected three
attempts before we found that hotspot atoms must sit on coordinate group 0.

Consequence: hotspots steer rather than constrain. Three of the top five land on
trastuzumab's epitope (`scripts/epitope_check.py`); two drifted. 27/30 pass novelty
against the SAbDab2 ML subset (5,301 unique antibodies — NOT full SAbDab), with a
shuffled control at 0.61 median confirming the threshold is meaningful.

rep117 and rep106 submitted. rep106 carries an N-x-S/T sequon in CDR2, flagged.

## Not done

- **Thermostability.** AbMelt needs MD at three temperatures per design. ThermoMPNN's own
  authors showed stability models fail on double mutations; ours carry 6-8. We report a
  sequence proxy and label it one.
- **Affinity.** No KD prediction. The 524k FACS labels and 768 SPR points were used for
  novelty, not as training data. That is the obvious next step.
- **Off-target.** Only a charge-patch proxy, no cross-reactivity screen.
- `pass` is empty for all 60 rows (missing pose-RMSD); ranking uses `min_ipae` directly.
- One diffusion sample per design, so no error bars.

## Files

| path | contents |
| --- | --- |
| `log.md` | every command, version, seed, timing, and the failures |
| `submission/team10.csv` | the submitted 20, abguard format |
| `submission/ranked_designs.csv` | full columns per design |
| `submission/RANKING.md` | readable rationale and caveats |
| `scripts/` | nine scripts, each with scope and caveats in its docstring |

Run outputs under `data/docs/runs/` are gitignored for size; the scored CSVs are here.
