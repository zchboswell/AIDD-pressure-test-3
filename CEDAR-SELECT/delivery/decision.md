# CEDAR next-tranche decision

**Conditionally switch from expanding the predecessor inventory to a small spiro-lactam discrimination experiment. Nominate D1 and D2 for gated synthesis, with existing examples23 (S2) and39 (S3) as the two diagnostic comparators.** Do not claim a selective cellular lead. Allocate at most **$120,000 inclusive of the $10,000 coverage work order**, over 30 business days, with explicit release gates. D3/D4 are documented reserve designs, not additional current synthesis nominations.

The new direction earns a limited experiment because there are matched historical CDK2/E1 and CDK1/B1 data and chemically actionable substitutions. It does not earn a large synthesis campaign: no independent cellular window is established, predictive ranking fails, and the structures do not validate the proposed spiro binding mode. The original inventory remains intact as predecessor context; its proposed experiments have no supplied completion results. Its simulated cellular observations are not evidence that this new series works in cells.

## Evidence that changes the decision

The expanded packet contains 60 compounds from one public study: 105 exact potency rows and 15 CDK1 bounds. The eight analogue pairs are already included, not sixteen extra measurements. Both paired assays use 1 mM ATP and RB1 TR-FRET, but different cyclins and enzyme concentrations; these are conditional assay ratios.

| Public example | CDK2 nM | CDK1 nM | CDK1/CDK2 | Decision |
|---|---:|---:|---:|---|
| 12 / S1, methoxy | 20 | 2240 | 112 | Independent-preparation anchor |
| 23 / S2, ethoxy | 15 | 2430 | 162 | Simpler primary comparator |
| 39 / S3, fluoroethoxy | 7 | 1200 | 171 | Potent primary comparator |
| 41 / S4, difluoroethoxy | 5 | 730 | 146 | Do not procure solely for potency |
| 56, trifluoroethoxy | 10 | 800 | 80 | Further fluorination not favored |
| 77, hydroxy | 121 | 5820 | 48 | Misses potency preference |

S2→S3 improves CDK2 potency 2.14-fold but CDK1 potency 2.03-fold: the 1.06-fold discrimination change is unconvincing without uncertainty. The supplied stereopair69/70 also shows nearly equal ratios. These data redirect chemistry toward a controlled ring-size experiment instead of chasing another fluorine.

Ring expansion on the supplied methyl background changes the ratio from about67 to110; on methoxy, CDK2 weakens from20 to114 nM and CDK1 becomes ≥10 µM. That second ratio is only ≥87.7, which **does not prove improvement over112**. D1/D2 test whether ethoxy/fluoroethoxy backgrounds tolerate ring expansion better; this interaction is unknown.

Executed fixed-hyperparameter fingerprint models offer little selectivity ranking value. Training group-CV selectivity R² is −0.158; exposed development R² is −0.089, MAE0.236 log units, versus median MAE0.253. Censored-CDK1 predictions violate all15 held-out bounds across the two partitions, even with one-sided-loss sensitivity. Same-study validation cannot establish another-study transfer. Numerical forecasts for new compounds and cellular engagement are withheld, not replaced with invented confidence intervals.

The structural analysis finds a real Tyr15 conformation difference, but 5NEV contains CDK2/compound73 and 5LQF contains CDK1/NU6102. Different ligands and states confound selectivity attribution. Rigid-transplant clash conclusions change with chain/alignment. This rejects a hard steric-exclusion story and any pose transfer to the spiro series; it does not disprove biochemical selectivity.

## What to make and measure

D1 is the cycloheptyl analogue of ethoxy example23; D2 is the corresponding fluoroethoxy analogue of example39, retaining specified stereochemistry. Exact canonical graphs, route precedents and identity requirements are in `designs.json` and `chemistry_designs.md`. They are absent from the supplied labeled packet and valid original-inventory graphs; no worldwide novelty claim is made. Use a common qualified chiral chlorospiro intermediate, split for coupling to the two protected aminopyrazoles, then deprotect and purify. These are multistep proposals unless advanced partners are in hand, not executed syntheses or verified supplier stock.

By day5 require structure-specific quotes, verified chiral intermediates and a credible day15 purified-product commitment. If unavailable, **do not force synthesis into the deadline**: finish comparator biology and retain unused funds. D3/D4 test a new ether-adjacent stereocenter but remain a paired reserve only; opening them requires a feasible revised schedule and budget.

Run same-plate paired CDK2/E1 and CDK1/B1 dose responses on D1/D2 and both parents, with three independent preparation/day replicates, fixed construct/activation records and 1 mM ATP. Test an ATP shift and an orthogonal catalytic readout. Preserve censored values; extend CDK1 concentrations only within demonstrated solubility/readout limits. Quantify active concentration and exclude aggregation/interference before interpreting changes.

A biochemical progression preference is CDK2≤100 nM and ratio≥10. To justify this particular edit, also require a reproducible ≥0.5-log ratio gain over its direct parent; analyze paired-day differences and their uncertainty. This is a prospective decision margin, not an established assay variance. A lower bound that overlaps the parent's range does not establish superiority.

For passing compounds, measure cellular CDK2 and CDK1 target engagement in the same cellular setting, plus CDK9 counterscreen, measured exposure, and a qualified orthogonal reference. An initial window preference is CDK2 engagement≤150 nM and ≥10-fold CDK1 separation; this is not a safety criterion. Measure proliferation/cell-cycle response in authenticated CCNE1-high/RB-intact and matched nondependent controls, with existing genetic perturbation/rescue controls where available. No new engineered cell line or animal study is promised within30 days. Phenotype without target engagement/dependency evidence is insufficient.

The [INX-315 primary study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10905675/) reports cellular engagement and supports an orthogonal comparator strategy; its assay values are not training labels for our designs. [Allosteric anthranilic-acid work](https://pmc.ncbi.nlm.nih.gov/articles/PMC10239507/) offers a different mechanism but inhibition of preassembled CDK/cyclin is a key limitation. Keep it as a future mechanism option, not the current30-day synthesis choice.

## Spend and release gates

| Included allocation | USD ceiling | Release condition |
|---|---:|---|
| Coverage earmark | 10,000 | $1,800 paid + $2,200 reserved + O2 $2,000 + O1 $1,400; $2,600 held |
| Chemistry | 32,000 | D1/D2 up to$20k after day5 gate; $12k reserve, not automatic D3/D4 spend |
| Identity, concentration, solubility and early ADME | 12,000 | Qualified parent material; ≥95% purity and stereochemical provenance |
| Paired biochemical replication and orthogonal/ATP checks | 18,000 | Qualified material; includes fresh matched parent panels |
| Cellular target engagement and exposure | 20,000 | Biochemical/QC gate passed |
| Dependency, proliferation and cell-cycle tests | 12,000 | Engagement window and suitable existing models |
| General contingency | 16,000 | Held for a costed revision |
| **Total** | **120,000** | Includes prior spend/reservations; no purchase made |

These are planning ceilings, not quotes. Coverage O2 recovers failed S3 CDK1; O1 intentionally repeats valid S1 CDK2. Preserve the pending S2 endpoint; defer redundant O3 and unpriced/unverified O4. Deduplicate obligations and vials: after coverage, S2 has4 mg and S3 has2 mg. Reserve1.2 mg of each for the parent plan. The two-lane coverage schedule finishes day5, within its day8 limit. Scenario contextsA/B are not mapped to historical contexts, so retrieve numerical results/metadata and run fresh paired panels rather than joining old and new endpoints into a ratio.

Days1–5: coverage and route/material/QC gates. Days6–15: gated synthesis and parent assay qualification. Days16–22: biochemical replication. Days23–30: passing-compound cellular work and decision. Missed material or assay gates narrow scope, not the standard of evidence.

**Reverse the recommendation** if independent parents fail reproducibility or usable biochemical discrimination, neither new edit improves the paired ratio, or a biochemical winner lacks cellular CDK2/CDK1 separation or dependency-consistent effects. Stop ring-expansion chemistry if its gain fails; retain a validated parent as a tool only if its own cellular evidence warrants it. No current data demonstrate improved discrimination beyond this study. The next tranche is designed to test that ambition, and can finish with a documented stop instead of a forced lead.
