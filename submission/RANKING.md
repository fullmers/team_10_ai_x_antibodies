# HER2 CDRH3 submission — ranked designs

Generated 2026-11-14 from `runs/her2-main` (RFD4-Proteina `her2-novel-2026-11-09`, 200 designs, seed 2; top 60 refolded with Boltz-2).

Parent CDRH3: `WGGDGFYAMD` (trastuzumab).

## How this list was built

**Slots 1-4 were chosen to differ by construction, not by score.** One representative per sequence family was taken, rather than the top four by `min_ipae`, so that a failure of one binding mode does not take the whole submission with it.

**Slots 5-20** are the remaining strict-novelty-passing refolded designs ranked on `min_ipae` with `iptm` as tie-break, skipping exact duplicates of loops already listed.

All 20 pass the **strict** novelty reading (identity < 70% to trastuzumab AND minimum Hamming distance > 2 to all 525,042 length-10 CDRH3s in the Chinery FACS + SPR screen).

## Sequence families

Families are single-linkage clusters over Hamming distance on the 10-mer, computed in `scripts/build_submission.py`. **The grouping depends on the cutoff, and the headline claim is only partly supported** — both columns are in the CSV so the reader can judge:

| claim | at d<=2 | at d<=3 |
| --- | --- | --- |
| `SLSLSNYAFD` and `DLSYYNHAFD` same family (slots 1, 3) | yes | yes |
| `TRAAVGYAHD`/`TTAAVGYAHD` same family (slots 2, 7) | yes | yes |
| `WSSAVGYAHD` same family as slot 2 (slots 2, 4) | **no** — singleton, 3 substitutions from `TTAAVGYAHD` | yes |

So **slots 1 and 3 are same-family at every cutoff tested**, which weakens the one-per-family framing: the four picks span three families at d<=2 and two at d<=3, not four. They still differ on the axes that motivated them (score, charge, novelty margin, conservatism).

**Slots 2 and 4 are a deliberate paired bet on the same family**, and their order is arbitrary: `TRAAVGYAHD` (0.3234) and `WSSAVGYAHD` (0.3238) are separated by 0.0004 in `min_ipae`, far below any meaningful resolution of the metric.

## Excluded by deliberate choice

`SLYVNNYSYD` (her2_cdrh3_39) scores 5th on `min_ipae` (0.3288) and has a novelty distance of 4, but carries an **N-x-S/T glycosylation sequon** (`NYS`). It is kept at slot 5 rather than promoted into the top four. Whether the sequon matters depends on the expression system: if the assay is cell-free without glycosylation machinery it is not a real liability, but **this could not be verified** — the Adaptyv document in the task materials contains no expression protocol, and the word 'glycosylation' appears nowhere in the documentation set.

`WGPPGFYAFD` (her2_cdrh3_70) scored 5th overall by `min_ipae` (0.3141) across all 60 refolded designs but is **excluded**: it is present verbatim in the FACS screen, labelled `high`. It is a known binder, not a novel design. Its rediscovery is a sanity check on the model, not a candidate.

## Caveats

- `pass` is empty for every row in `refold_metrics.csv` (`pass_note`: missing `refolded_target_aligned_binder_rmsd.rmsd`), so the module's own automated verdict is unavailable for all 60. Ranking uses `min_ipae` directly.

- Developability is a sequence-motif prefilter only, scanned in framework context (`SR[loop]YW`), not an assessment. See `scripts/developability.py`.

- One diffusion sample per design (`--num-diffusion-samples 1`, seed 42). No replicate spread.

- Boltz-2 `iptm` here is not comparable to the Chai-1 calibration in `scripts/calibration.csv`; ipTM does not transfer across predictors.

## The list

| # | CDRH3 | design id | min_ipae | iptm | b_pLDDT | dist | id % | mut | hyd | chg | dev |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | SLSLSNYAFD | her2_cdrh3_140 | 0.3076 | 0.926 | 98.59 | 3 | 30 | 7 | 0.50 | -1 | clean |
| 2 | TRAAVGYAHD | her2_cdrh3_98 | 0.3234 | 0.903 | 98.32 | 4 | 30 | 7 | 0.50 | +0 | clean |
| 3 | DLSYYNHAFD | her2_cdrh3_57 | 0.3559 | 0.885 | 97.94 | 5 | 20 | 8 | 0.50 | -2 | clean |
| 4 | WSSAVGYAHD | her2_cdrh3_31 | 0.3238 | 0.926 | 98.16 | 3 | 40 | 6 | 0.50 | -1 | clean |
| 5 | SLYVNNYSYD | her2_cdrh3_39 | 0.3288 | 0.902 | 98.13 | 4 | 20 | 8 | 0.50 | -1 | FLAG |
| 6 | DLFFPYYAFD | her2_cdrh3_27 | 0.3454 | 0.892 | 95.99 | 3 | 30 | 7 | 0.70 | -2 | FLAG |
| 7 | TTAAVGYAHD | her2_cdrh3_37 | 0.3481 | 0.879 | 98.29 | 4 | 30 | 7 | 0.50 | -1 | clean |
| 8 | YMSLYYYEFD | her2_cdrh3_184 | 0.3504 | 0.909 | 95.28 | 3 | 20 | 8 | 0.70 | -2 | FLAG |
| 9 | YMSLYYYAFD | her2_cdrh3_134 | 0.3578 | 0.882 | 96.62 | 3 | 30 | 7 | 0.80 | -1 | FLAG |
| 10 | DLSSSDHAFD | her2_cdrh3_192 | 0.3618 | 0.856 | 98.22 | 4 | 20 | 8 | 0.30 | -3 | FLAG |
| 11 | SLSLSDYAFD | her2_cdrh3_43 | 0.3646 | 0.879 | 97.52 | 3 | 30 | 7 | 0.50 | -2 | clean |
| 12 | DLSSSDFAFD | her2_cdrh3_147 | 0.3650 | 0.882 | 98.06 | 4 | 20 | 8 | 0.40 | -3 | FLAG |
| 13 | SLSLSNFAFD | her2_cdrh3_9 | 0.3663 | 0.836 | 98.60 | 4 | 20 | 8 | 0.50 | -1 | clean |
| 14 | WRSLSEYYFD | her2_cdrh3_148 | 0.3692 | 0.877 | 97.18 | 3 | 30 | 7 | 0.50 | -1 | clean |
| 15 | DLSLSEHAFD | her2_cdrh3_35 | 0.3700 | 0.838 | 98.18 | 4 | 20 | 8 | 0.40 | -3 | FLAG |
| 16 | DLSLHNHAFD | her2_cdrh3_181 | 0.3726 | 0.878 | 97.41 | 4 | 20 | 8 | 0.40 | -2 | clean |
| 17 | WASLHDYYHD | her2_cdrh3_193 | 0.3735 | 0.893 | 95.63 | 3 | 30 | 7 | 0.50 | -2 | clean |
| 18 | YMNLYYYAFD | her2_cdrh3_149 | 0.3783 | 0.877 | 96.05 | 3 | 30 | 7 | 0.80 | -1 | FLAG |
| 19 | YASSSDYYFD | her2_cdrh3_114 | 0.3785 | 0.877 | 97.20 | 3 | 20 | 8 | 0.50 | -2 | clean |
| 20 | WKSTTQYYFD | her2_cdrh3_22 | 0.3839 | 0.875 | 96.85 | 4 | 30 | 7 | 0.40 | +0 | clean |

## Rationales

**1. `SLSLSNYAFD` (her2_cdrh3_140)** — Best overall score (min_ipae 0.3076, the lowest of all 37) with clean developability: no liability motif, hydrophobic fraction equal to the parent, net charge -1.
**2. `TRAAVGYAHD` (her2_cdrh3_98)** — Second sequence family, one of only two neutral-charge designs, larger novelty margin than its family-mate at 4 vs 3.
**3. `DLSYYNHAFD` (her2_cdrh3_57)** — Largest novelty margin in the set: minimum Hamming distance 5 to the Chinery screen, 20% identity to trastuzumab, 8 mutations. Clean developability.
**4. `WSSAVGYAHD` (her2_cdrh3_31)** — Conservative hedge: closest to the parent at 6 mutations, and its iptm (0.926) ties the top-ranked design. Same family as slot 2, taken as a deliberate paired bet.
5. `SLYVNNYSYD` (her2_cdrh3_39) — Rank by score among the remaining strict-novelty designs; family 5 (d<=2). DEVELOPABILITY FLAG: N-x-S/T +1.
   - Developability flag: N-x-S/T +1
6. `DLFFPYYAFD` (her2_cdrh3_27) — Rank by score among the remaining strict-novelty designs; family 2 (d<=2). DEVELOPABILITY FLAG: hydrophobic 70% vs 50%.
   - Developability flag: hydrophobic 70% vs 50%
7. `TTAAVGYAHD` (her2_cdrh3_37) — Rank by score among the remaining strict-novelty designs; family 7 (d<=2).
8. `YMSLYYYEFD` (her2_cdrh3_184) — Rank by score among the remaining strict-novelty designs; family 13 (d<=2). DEVELOPABILITY FLAG: hydrophobic 70% vs 50%.
   - Developability flag: hydrophobic 70% vs 50%
9. `YMSLYYYAFD` (her2_cdrh3_134) — Rank by score among the remaining strict-novelty designs; family 13 (d<=2). DEVELOPABILITY FLAG: 5+ hydrophobic run +1; hydrophobic 80% vs 50%.
   - Developability flag: 5+ hydrophobic run +1; hydrophobic 80% vs 50%
10. `DLSSSDHAFD` (her2_cdrh3_192) — Rank by score among the remaining strict-novelty designs; family 3 (d<=2). DEVELOPABILITY FLAG: charge -3 vs -2 (further from neutral).
   - Developability flag: charge -3 vs -2 (further from neutral)
11. `SLSLSDYAFD` (her2_cdrh3_43) — Rank by score among the remaining strict-novelty designs; family 3 (d<=2).
12. `DLSSSDFAFD` (her2_cdrh3_147) — Rank by score among the remaining strict-novelty designs; family 3 (d<=2). DEVELOPABILITY FLAG: charge -3 vs -2 (further from neutral).
   - Developability flag: charge -3 vs -2 (further from neutral)
13. `SLSLSNFAFD` (her2_cdrh3_9) — Rank by score among the remaining strict-novelty designs; family 3 (d<=2).
14. `WRSLSEYYFD` (her2_cdrh3_148) — Rank by score among the remaining strict-novelty designs; family 9 (d<=2).
15. `DLSLSEHAFD` (her2_cdrh3_35) — Rank by score among the remaining strict-novelty designs; family 3 (d<=2). DEVELOPABILITY FLAG: charge -3 vs -2 (further from neutral).
   - Developability flag: charge -3 vs -2 (further from neutral)
16. `DLSLHNHAFD` (her2_cdrh3_181) — Rank by score among the remaining strict-novelty designs; family 3 (d<=2).
17. `WASLHDYYHD` (her2_cdrh3_193) — Rank by score among the remaining strict-novelty designs; family 8 (d<=2).
18. `YMNLYYYAFD` (her2_cdrh3_149) — Rank by score among the remaining strict-novelty designs; family 13 (d<=2). DEVELOPABILITY FLAG: 5+ hydrophobic run +1; hydrophobic 80% vs 50%.
   - Developability flag: 5+ hydrophobic run +1; hydrophobic 80% vs 50%
19. `YASSSDYYFD` (her2_cdrh3_114) — Rank by score among the remaining strict-novelty designs; family 12 (d<=2).
20. `WKSTTQYYFD` (her2_cdrh3_22) — Rank by score among the remaining strict-novelty designs; family 10 (d<=2).

---
Full columns in [`ranked_designs.csv`](ranked_designs.csv).
