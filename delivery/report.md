# CEDAR next-cycle decision

**Advance CED-02, CED-03, CED-08, CED-13, CED-17 and CED-22 into a decision cycle; add CED-23 as the sole diagnostic.** These are provisional experimental priorities, not validated leads. The proposal costs **$17,500**, uses **seven unique vials**, leaves **$500**, and completes by **business day 7** under the supplied parallel scheduling terms. No new experiments have been run. All current engagement, CDK1, ADME and inventory values are **simulated program observations attached to real chemical structures**; historical ChEMBL measurements and 1H1Q coordinates are distinct public evidence.

## Why these compounds

| Compound | Role and decision it resolves | Main risk |
|---|---|---|
| **CED-22** | Leading observed anchor: cellular engagement 130 nM, HLM half-life 22 min, solubility 70 µM. Establish its matched biochemical window and confirm stability. | CDK1 400 nM is not a demonstrated selective window; stability barely exceeds preference. |
| **CED-02** | Second anchor with a lactam-bearing chemical variation: 130 nM engagement, 22 min HLM, 18 µM solubility. | Same borderline profile; a chemical alternative rather than independent replication of CED-22. |
| **CED-17** | Focused cyclopropyl-versus-tert-butyl comparison with CED-22, with 130 nM engagement, 22 min HLM and 120 µM solubility. Historical IC50 is 700 versus 54.9 nM for CED-22. | Same reported current profile does not prove reproducibility. Historical/current rank changes reflect different assays, not an established contradiction. |
| **CED-08** | High-value unknown: engagement unmeasured, CDK1 1,800 nM, HLM 38 min, solubility 18 µM. Two new engagement runs resolve the central gap. | Favorable ADME/CDK1 context cannot predict cellular CDK2 potency. |
| **CED-03** | Second unknown: both activities missing, HLM 38 min, solubility 35 µM. Check VIAL-03 first, then two engagement runs. | STOCK-28 has the same label with a different graph. Quarantine that export; do not infer that all CED-03 stock is verified. |
| **CED-13** | Conditional tradeoff candidate: engagement 320 nM, CDK1 1,800 nM, HLM 38 min. Tests whether useful matched separation compensates for weaker engagement. | Misses the approximately 150 nM preference; 8 µM solubility can limit interpretation of upper assay concentrations. |
| **CED-23** | Diagnostic: well-soluble structural-reference family, HLM 38 min, unmeasured biology. Tests a different chemotype rather than assuming that a crystal structure validates the main series. | Not an established biological positive control. Inventory and deposited ligand encode different purine tautomers; no cellular or CDK1 window follows from the structure. |

CED-22 is the first observed-profile choice, followed by CED-02; CED-17 earns its slot through the named substitution comparison. CED-03/08 are information priorities, not model-ranked winners. We omit CED-01 as a diagnostic because its existing 65 nM engagement, 90 nM CDK1 and 14 min HLM already describe an unfavorable profile; measuring the unknowns twice changes more decisions. CED-18 illustrates why the old shortcut fails: historical IC50 22.5 nM does not overcome current engagement 320 nM and solubility 3 µM. CED-13 has the same supplied engagement/CDK1/HLM values and higher solubility. CED-19 remains a stability/window fallback, but 850 nM engagement is weak. CED-24 has unmeasured endpoints and 3 µM solubility versus CED-23's 120 µM.

CED-04/12 arrive on day 10 and cannot enter this cycle. CED-05's 72% purity would require IDENTITY before interpreting new biology; its 35/45 nM engagement/CDK1 profile and 8 min HLM do not justify that spend. Duplicate shipment rows add no material. Every inventory request is accounted for in [selection.csv](selection.csv), including invalid STOCK-27 and identity-conflicted STOCK-28. CED-03 is included only through STOCK-03/VIAL-03, with identity-first work. LC-MS distinguishes the gross label/graph discrepancy; it does not establish every stereochemical or protonation detail.

## What to order and how to decide

Run **CDK2-E1, CDK1-B and ENGAGE on all seven compounds**. Order **two separate ENGAGE experiments for CED-03 and CED-08**, **IDENTITY for CED-03**, and **MIC for CED-02 and CED-22**. These are 24 unique menu combinations and 26 charged experimental runs. The extra MIC runs test a consequential 22 min value near the 20 min preference. Existing SOL data inform concentration feasibility; repeating SOL everywhere would displace higher-value endpoint measurements. Inspect nominal versus achieved concentrations and precipitation/curve behavior as assay-quality requirements; this plan does not pretend the menu includes an extra orthogonal assay.

| Compound | Total USD, including handling | Mass used / stock mg | Completion day |
|---|---:|---:|---:|
| CED-02 | 2,500 | 0.45 / 1.30 | 7 |
| CED-03 | 3,300 | 0.55 / 1.60 | 6 |
| CED-08 | 2,950 | 0.50 / 1.60 | 7 |
| CED-13 | 2,050 | 0.35 / 1.60 | 6 |
| CED-17 | 2,050 | 0.35 / 1.30 | 7 |
| CED-22 | 2,500 | 0.45 / 1.30 | 6 |
| CED-23 | 2,150 | 0.35 / 1.60 | 7 |
| **Total** | **17,500** | No duplicate-vial stock counted | **7** |

Handling is charged once per physical compound: $1,650 total; assays cost $15,850. Arrival plus identity first, where present, plus longest parallel assay determines completion. Both ENGAGE repeats are booked prospectively in parallel; a sequential five-day confirmation after seeing a result would not fit. If VIAL-03 fails identity/purity, stop its biological interpretation and retain the failure. There is no silently added replacement or assumed assay refund. The exact order file is [experiments.csv](experiments.csv); independent accounting is in [budget_material_schedule.csv](support/budget_material_schedule.csv).

Use these predeclared review rules, as **team triage choices rather than validated clinical thresholds**:

1. Reject uninterpretable curves, unresolved identity, or concentrations beyond demonstrated assay feasibility. Preserve qualifiers and failures. Identity failure cannot be repaired by averaging biological results.
2. Prefer cellular engagement below approximately 150 nM and HLM above 20 min, as the program requests. For CED-03/08 require the two new engagement runs to agree on that decision; a greater than threefold discrepancy or opposite sides of the preference is inconclusive, not a fabricated confidence interval. Borderline values require judgment, not rounding into a pass.
3. Compute **CDK1-B IC50 / CDK2-E1 IC50** from the new matched 1 mM ATP panel. Adopt **at least tenfold** as a provisional useful biochemical-window goal for this cycle. This is a new operational criterion; report sensitivity at fivefold and twentyfold. Carry one-sided limits through the ratio. Do not calculate selectivity by dividing historical assays or cellular CDK2 engagement into biochemical CDK1 activity.
4. Fund another focused cycle if at least one original-series candidate has credible useful engagement, a useful matched biochemical window and acceptable stability/assay behavior. CED-13 is a fallback only if its prospective result improves materially or offers a compelling measured tradeoff; the current 320 nM value does not meet preference. If only CED-23 succeeds, consider a separate chemotype direction with identity/state verification. If no compound meets the combined objective, pause expansion rather than retaining the old biochemical threshold.

The menu lacks cellular CDK1 engagement and cell viability. Thus even a successful cycle cannot establish a cellular safety window, CDK2-dependent antiproliferative mechanism, or in-vivo efficacy. Those would be the next evidentiary gap before a therapeutic-selectivity claim, not uncosted work hidden in this order. A lower-spend alternative omitting CED-17 costs $15,450 but forfeits the substitution question; removing CED-13 too costs $13,400 and forfeits the potency/window fallback. The chosen spend buys both questions; filling all eight slots is unnecessary.

## Evidence and model audit

The repaired [assay_math.py](../analysis/assay_math.py) preserves IC50/Ki/Kd, supports nM/uM/µM/μM/mM/M, reverses inequalities under the negative logarithm, and rejects unsupported, nonfinite or nonpositive inputs. Five regression-test methods cover the legacy failure modes. For example, ACT-9004 ≤0.01 µM becomes ≤10 nM and pIC50 ≥8.000; ACT-9008 Ki ≥0.3 µM remains Ki, ≥300 nM and pKi ≤6.522879. These fictional QC records are not training labels. Zero, negative, unsupported ng/mL and invalid-graph records are quarantined.

All **266 activity rows** are accounted for: **122 eligible, 138 context-only, 2 duplicate and 4 quarantined**. All 258 source-linked records agree with the 256 supplied originals, including the two reexports. Duplicate accounting uses nonempty source_activity_id, not the potential_duplicate flag; chemistry grouping uses canonical isomeric graphs with charge/fragments/stereochemistry preserved. The 122 exact historical labels comprise 117 chemical groups. Different source observations of the same graph remain measurements, grouped together in validation. No salt, tautomer or assay equivalence is silently imposed.

The declared modeling endpoint is **exact historical IC50 in CHEMBL5736732**, not present cellular engagement. Its multi-kinase description lists CDK2/cyclin A2 at 15 µM ATP and CDK2/cyclin E1 at 20 µM, among other conditions; it does not uniquely assign each row to the new E1/1 mM ATP assay. CHEMBL661125 has incompletely specified context; CHEMBL3705366 describes CDK2/cyclin A at 20 µM ATP through a related protocol. They and censored thresholds remain context. No assay correction or cross-assay pooling is justified by this snapshot.

**Do not use the vendor model to choose compounds.** Its original R² **0.960** is reproducible but relies on a retrospective activity-QC index unavailable in inventory, shares nine chemical identities across its random split, mixes four assays/two endpoint types, and treats 49 nonexact responses as labels. The index correlates **0.999609** with corrected historical labels. On a shared curated holdout, R² drops from **0.989** with the index to **0.370** with fingerprints only. Honest historical chemical-group CV gives MAE **0.330 log units**, R² **0.344**; scaffold CV gives **0.357**, **0.249**, with four of five within-fold R² values nonpositive. The fingerprint model beats the median baseline on pooled MAE but has uneven, same-document evidence. This earns chemical-similarity context, not a current-endpoint forecast.

CED-17–22 match historical training chemicals; those labels are observations, not prediction successes. CED-23/24 have maximum fingerprint similarity approximately 0.192/0.200 to the historical set, below the domain tested in either validation. We abstain from forecasts and invented intervals. [predictions.csv](predictions.csv) instead reports the **18 supplied cellular observations**, **6 explicit abstentions**, **2 duplicates** and **2 quarantines**, with assay attribution and no unsupported uncertainty bounds. [model_audit.md](model_audit.md) and its saved held-out predictions retain methods, failures and quantitative evidence.

## What the molecular views establish

The deposited source is **1H1Q, NU6094 / component 2A6, author chain A residue 1298**, bound to CDK2 chain A with cyclin A2 chain B, at 2.5 Å resolution. Complete heavy-atom graph mapping, without fitting, gives original-frame RMSD **0.000 Å** for the reference export and **8.774964 Å** for the priority export. The latter is translated **(+8, −3, +2) Å**. Row-order comparison falsely gives 6.152 Å for the correct reference. Fitting the priority returns approximately zero but conceals its broken receptor relationship.

The installed RDKit `CalcRMS` with an explicit atom map and `RenumberAtoms` resolve the old workflow obstacle; an independent coordinate calculation confirms the result and coordinate immutability. **The priority export is not usable for supplied-frame contacts.** Separately fitted/reordered artifacts are labeled demonstrations of reference recovery, not newly supported candidate poses. Source polar proximities to Glu81 and Leu83 are consistent with a hinge-binding discussion, but distances alone do not establish hydrogen bonds, protonation, binding affinity or CDK1 selectivity. CED-23 shares the source tautomer family; CED-24 is a related graph without a supplied pose. Neither authorizes pose transfer to the main series.

Use [the deposited pocket](figures/deposited_pocket.png), [frame audit](figures/structure_audit.png), [selected chemical structures](figures/selected_chemical_structures.png), [decision comparison](figures/decision_comparison.png) and [model diagnostic](figures/model_diagnostic.png). Numerical mappings and interface evidence accompany [structure_checks.json](structure_checks.json). No docking score or molecular movie is needed to answer this coordinate-validity question.

## Execution and handoff

Scientific analysis, converter checks, molecular rendering and a clean-output replay are completed; experimental outcomes remain prospective. [REPRODUCE.md](REPRODUCE.md) gives exact commands and limits; [NEXT_BATCH.md](NEXT_BATCH.md) carries reusable findings with reassessment conditions. The supplied input package remains byte-identical locally. Automatic approval review restricted remote capture of supplied inputs and administrative files; the supervisor directed generated-analysis-only publication. The private repository therefore needs the separately supplied input package for replay. This operational limitation is not a scientific success/failure metric or an AIDD capture claim. Exact exclusions and rejections are recorded in `admin/CAPTURE_LIMITATION.json`.

Source attribution: supplied ChEMBL_37 official API snapshot (CC BY-SA 3.0), PDB 1H1Q (CC0), and OpenFF benchmark lig_1h1q/lig_1oiu revision fd88824f9114244f95a14b485e6d6c96c1de716d (CC BY 4.0); see [ATTRIBUTION.md](ATTRIBUTION.md) and local `inputs/project/PROVENANCE.md`. All curation, analysis, derived depictions and operational recommendations here are adaptations; public redistribution is not authorized.
