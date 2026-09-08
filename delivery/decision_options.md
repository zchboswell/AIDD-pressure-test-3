# Independent decision review

This is a proposed next-cycle design, not executed experimental work. All inventory, budget and recent program measurements below are supplied **simulated** records attached to real chemical graphs. The public biochemical archive and deposited structure provide separate, endpoint-specific evidence; neither measures the current cellular objective.

## Recommended provisional slate

| Role | Compound | Decision value and principal risk |
|---|---|---|
| Advance | CED-22 | Cellular engagement 130 nM, solubility 70 µM, microsomal half-life 22 min. Balanced starting anchor, but CDK1 400 nM and borderline stability need scrutiny. |
| Advance | CED-17 | Same supplied engagement/CDK1/stability values as CED-22, solubility 120 µM; chemically distinct substitution tests whether the preferred profile survives prospective matched assays. Repeated values are not independent evidence of reproducibility. |
| Advance | CED-02 | 130 nM engagement, 22 min half-life, 18 µM solubility; distinct lactam-bearing analog gives another chemically concrete test of the leading observed profile. |
| Advance | CED-08 | Unmeasured engagement, CDK1 1,800 nM, half-life 38 min, solubility 18 µM: high decision value from closing the cellular evidence gap. No cellular potency prediction is warranted from these properties alone. |
| Advance | CED-03 | Both primary activities missing, but half-life 38 min and solubility 35 µM support an informative prospective test. Use only valid STOCK-03/VIAL-03. Quarantine STOCK-28/VIAL-CONFLICT, whose graph conflicts with its CED-03/source label; do not add its mass to VIAL-03. |
| Advance | CED-13 | 320 nM engagement is weaker than preference, but half-life 38 min and CDK1 1,800 nM justify a conditional exploration of the potency/window tradeoff. Solubility 8 µM limits confidence at high assay concentrations. |
| Diagnostic | CED-01 | 65 nM engagement but CDK1 90 nM and half-life 14 min. Prospective comparator tests the failure mode of prioritizing potency alone; not a validated selective positive control. |
| Diagnostic | CED-23 | Structurally anchored, different chemical graph with both program activities missing, half-life 38 min and solubility 120 µM. Tests whether the deposited-ligand reference contributes useful current biology. Its deposited structure does not validate other compounds' poses or current selectivity. |

These are conditional advancement candidates into experiments, not six validated leads. CED-13 is deliberately a weaker-engagement fallback. CED-19 is an alternative fallback with 55 min half-life and CDK1 6,000 nM, but engagement 850 nM is substantially worse; its apparent cross-assay separation cannot override that shortfall. CED-07 is a reasonable reserve but repeats the 130/400/22 profile with only 8 µM solubility. Do not use stock that arrives on day 10 (CED-04/12) for an eight-day cycle. CED-05 is 72% pure and would require IDENTITY before biological interpretation; its existing 35/45 nM activities and 8 min half-life do not justify this limited cycle's slot. CED-INVALID and the conflicting STOCK-28 are quarantined; STOCK-25/26 duplicate VIAL-01/02 and add no mass.

## Concrete experiment option

Run CDK2-E1, CDK1-B and ENGAGE once on every selected compound. CDK2-E1/CDK1-B are a prospective matched-ATP biochemical panel; ENGAGE measures the requested cellular endpoint. Run IDENTITY on VIAL-03 before its biological assays as a precaution against the conflicting inventory export, while retaining the quarantine on VIAL-CONFLICT. Add a second ENGAGE experiment on CED-08 because its missing cellular result is the strongest plausible decision-changing gap. Schedule both ENGAGE experiments prospectively in parallel, rather than promising a sequential five-day confirmation after seeing a five-day result.

| Compound | Handling USD | Assays USD | Total USD | Mass used / available mg | Arrival + identity + longest assay, days |
|---|---:|---:|---:|---:|---:|
| CED-01 | 200 | 1,850 | 2,050 | 0.35 / 1.00 | 1 + 0 + 5 = 6 |
| CED-02 | 250 | 1,850 | 2,100 | 0.35 / 1.30 | 2 + 0 + 5 = 7 |
| CED-03 | 300 | 2,050 | 2,350 | 0.40 / 1.60 | 0 + 1 + 5 = 6 |
| CED-08 | 150 | 2,800 | 2,950 | 0.50 / 1.60 | 2 + 0 + 5 = 7 |
| CED-13 | 200 | 1,850 | 2,050 | 0.35 / 1.60 | 1 + 0 + 5 = 6 |
| CED-17 | 200 | 1,850 | 2,050 | 0.35 / 1.30 | 2 + 0 + 5 = 7 |
| CED-22 | 250 | 1,850 | 2,100 | 0.35 / 1.30 | 1 + 0 + 5 = 6 |
| CED-23 | 300 | 1,850 | 2,150 | 0.35 / 1.60 | 2 + 0 + 5 = 7 |
| Total | 1,850 | 15,950 | **17,800** | 8 unique physical compounds | **latest day 7** |

The proposal leaves $200 unspent. The experimental combination ENGAGE/CED-08 has replicates=2; every other combination has replicates=1. No purity correction creates extra usable material. If VIAL-03 identity fails, stop interpreting its downstream assays and preserve the failure; an unplanned replacement cannot silently exceed the eight-compound cap. No optional extra experimental cost or delayed conditional assay is assumed to fit.

## Decision gates and interpretation

1. Identity and assay validity first. The CED-03 check must support the expected graph and adequate purity; no result resolves VIAL-CONFLICT by association. Review concentration-response quality and solubility limits before interpreting potency or ratios. The menu has no orthogonal cellular CDK1 or viability assay, so this cycle cannot establish a cellular safety window or cytotoxicity mechanism.
2. Prefer engagement below approximately 150 nM and microsomal half-life above 20 min, exactly as project preferences rather than hard clinical thresholds. The supplied 22 min values are close to preference and have no supplied uncertainty. CED-13 needs a material improvement or compelling prospective window to earn the subsequent cycle; its current result does not meet the engagement preference.
3. Calculate prospective biochemical CDK1-B IC50 / CDK2-E1 IC50 under the matching 1 mM ATP assays. A working goal of at least tenfold separation would be a new team decision threshold, not a supplied or clinically validated criterion. Retain raw values and censoring bounds; a ratio bounded above/below must stay bounded. Seek agreement with cellular engagement without declaring the biochemical ratio a cellular CDK1 window.
4. Existing CDK1-B / cellular CDK2-engagement quotients are descriptive cross-assay comparisons only. For example, 400/130 does not establish threefold intrinsic selectivity. ATP, compartment, permeability, protein binding and endpoint differences prevent that inference. No compound is currently proven to have the desired matched window.
5. Continue a focused chemistry cycle if at least one credible candidate combines reproducible useful engagement with a convincing matched biochemical window and acceptable stability/assay behavior. If engagement is weak in CED-08/03 and current anchors remain narrowly separated, pause expansion and resolve the cellular CDK1/viability gap before claiming therapeutic selectivity. If only CED-23 succeeds, consider a separate chemotype direction rather than transferring its pose to the original series.

The vendor model can earn at most a historical endpoint-specific supporting role after leakage-resistant validation; it cannot replace the proposed cellular measurement or infer unmeasured CDK1 results. A favorable random-row historical score is not evidence for the changed program endpoint. The proposed six-candidate selection follows current simulated evidence and missing-information value, not a biochemical IC50 <100 nM shortcut.

## Preferred refinement: seven compounds, greater replication value

After independent challenge, prefer the parent's revised seven-compound design over the eight-compound option above. Remove CED-01: its existing unfavorable stability and cross-assay activity profile make repeat testing lower value, and it is not an independently validated assay control. Keep the same six provisional advances and diagnostic CED-23. Perform the matched three-assay panel on all seven, IDENTITY first on CED-03, two ENGAGE experiments each for CED-03 and CED-08, and one new MIC experiment each on CED-02 and CED-22. These MIC experiments are prospective repeats of supplied measurements, charged once each; the two new ENGAGE experiments each are charged twice. No hidden technical replicate is included.

| Compound | Total USD including handling | Total material mg | Available mg | Completion business day |
|---|---:|---:|---:|---:|
| CED-02 | 2,500 | 0.45 | 1.30 | 7 |
| CED-03 | 3,300 | 0.55 | 1.60 | 6 |
| CED-08 | 2,950 | 0.50 | 1.60 | 7 |
| CED-13 | 2,050 | 0.35 | 1.60 | 6 |
| CED-17 | 2,050 | 0.35 | 1.30 | 7 |
| CED-22 | 2,500 | 0.45 | 1.30 | 6 |
| CED-23 | 2,150 | 0.35 | 1.60 | 7 |
| **Total** | **17,500** | | | **7** |

This leaves $500 and one compound slot unused. Totals were independently checked from the supplied CSVs with Python Decimal arithmetic, one thread on leased CPU 9. The runtime launcher did not support the first stdin `-` invocation and failed before running analysis; the subsequent `-c` invocation succeeded. No result of the failed invocation was used. The check used no external data and changed no input files.

CED-22 ranks first among the three compounds with supplied 130 nM engagement and 22 min half-life because it offers a practical 70 µM solubility anchor. CED-02 ranks next because the lactam-containing chemistry adds a different structural comparison; its 18 µM solubility is adequate as a relative starting point but does not establish assay validity. CED-17 ranks third for standalone advancement value, but retaining it produces a focused N-cyclopropyl versus N-tert-butyl comparison with CED-22. This is a chemical comparison, not a biological replicate. The parent's historical activity comparison (700 versus 54.9 nM) should be supported by the final curation table and kept endpoint-specific: identical supplied current 130 nM engagement does not contradict different historical biochemical measurements. It illustrates why a historical 100 nM cutoff is not an appropriate current cellular eligibility rule.

CED-03 and CED-08 deserve the extra ENGAGE experiments because both have unmeasured cellular engagement and plausible stability. Results there can change whether a useful lead exists, whereas additional confirmation of CED-01's unfavorable profile is less likely to change the decision. Run their two experiments in parallel after any identity prerequisite; two experiments do not by themselves establish a calibrated confidence interval, and a disagreement should trigger qualification of the result rather than an invented precise mean. New MIC data for CED-02 and CED-22 test the fragility of the 22 min values near the 20 min preference. CED-17's stability remains supported only by its supplied measurement; the pair comparison cannot substitute for a new MIC result on that compound.

If a smaller cycle is preferred, dropping CED-17 yields **five advances plus CED-23**, costs **$15,450**, and preserves repeated engagement on both unknowns and both new MIC measurements. The explicit lost question is the cyclopropyl/tert-butyl substitution comparison. Dropping CED-13 as well yields **four advances plus CED-23**, costs **$13,400**, and focuses entirely on observed near-preference or unmeasured engagement cases, but sacrifices the conditional 320 nM engagement/38 min stability fallback. Neither saved-budget option should be called scientifically inferior by default. The seven-compound option is defensible because its two extra compounds answer named comparison and tradeoff questions; it is not necessary to fill eight slots or exhaust $18,000.
