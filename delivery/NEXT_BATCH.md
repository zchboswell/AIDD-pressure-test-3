# Portable CEDAR handoff

**Decision now:** six provisional advances CED-02/03/08/13/17/22 plus diagnostic CED-23; $17,500, seven unique vials, latest day7. Exact orders and reasons are in experiments.csv and selection.csv. CED03 identity precedes biology; its conflicted export stays quarantined. All repeat ENGAGE runs are booked in parallel. No experiments have yet been executed.

| Reusable finding | Supporting evidence | Applies when | Reassess when / practical action |
|---|---|---|---|
| Count physical vials and source observations, not export rows. | STOCK25/26 duplicate VIAL01/02; ACT9001/9002 reproduce existing source_activity_id. Identity-conflicted STOCK28 is a separate quarantine. | Same vial/source-ID conventions with consistent chemistry and provenance. | New IDs, salts, stereochemistry, inconsistent mass/purity or labels require reconciliation. Never sum duplicated stock. |
| Historical potency is assay-specific context for the changed objective. | CED18 historical22.5nM/current TE320nM; CED17 historical700nM/current TE130nM. Primary prose is multi-kinase; new panel is1mM ATP. | Current question is cellular engagement plus CDK1 window. | New matched assays may support a new relationship. Do not reinstate the <100nM shortcut without endpoint-specific validation. |
| Labels and inequality direction must survive normalization. | ACT9004 gives pIC50≥8; ACT9008 remains pKi≤6.522879; invalid values quarantined. Five regression methods pass. | Supported concentration endpoints and units in the repaired callable. | New unit/endpoint/qualifier needs explicit support and tests. Never infer nM for unknown units or clip nonpositive values. |
| Retrospective QC features are not prospective predictors. | vendor index r=.999609 to corrected label; vendor random R².960 vs modest grouped historical signal. | Feature timing follows the response, or field absent at decision time. | Independently acquired pre-assay features and external endpoint-matched validation could justify a new model; never deploy this vendor score by deleting its feature after training. |
| Atom reordering and receptor-frame validity are different questions. | CalcRMS+explicit graph map: reference0.000Å, priority8.774964Å; priority translation(+8,-3,+2)Å; fitting hides this. | Comparing chemically equivalent exports with a common receptor frame. | Changed chemistry/stereo/atom completeness or legitimate receptor transform needs explicit registration provenance. Prefer deposited source; do not use original priority contacts. |
| Tautomer grouping is not bound-state evidence. | CED23 and 1H1Q ligand canonical tautomers match but encoded graphs differ. | Identity linkage across representations. | Hydrogen-bond assignments, new poses or protonation-sensitive chemistry need state-aware evidence. No transfer of the 2A6 pose to the main series. |
| Buy decision-changing measurements, not a full menu for every stock. | Seven-compound plan replaces known-poor comparator01 with repeated unknown TE03/08 and borderline MIC02/22; exact accounting passes. | Same cycle budget, stock and schedule. | Changed prices, arrival, material or assay durations require recomputation. No sequential five-day rescue is promised within eight days. |

Next review: evaluate identity, curve validity, repeat agreement, approximately150nM engagement preference, >20min HLM preference, and prospective matched CDK1-B/CDK2-E1 biochemical ratio. A tenfold ratio is this team's provisional goal; show fivefold/twentyfold sensitivity. It does not establish cellular selectivity. Stop expansion if no credible compound combines the required properties; if only CED23 works, discuss a separate chemotype. Cellular CDK1 and viability are unmeasured and unavailable in this cycle's menu.

Keep useful history (vial accounting) and retire only the context-specific shortcuts (biochemical cutoff, random-row model, fitted/row-order export integrity). These are project findings with explicit applicability, not company-wide policy or external guidance.

## Portable assembly

Create a small code/report bundle from the final generated-analysis commit; the supplied input package is intentionally not republished due the recorded approval limitation:

```sh
git archive --format=tar.gz --output=/tmp/cedar-generated-replay.tar.gz HEAD analysis delivery admin/PREFLIGHT.json admin/RUNTIME_CHECKS.json admin/INPUT_VERIFICATION.json admin/CAPTURE_LIMITATION.json
```

Unpack in a new directory, provide the separately supplied `inputs/` folder, and run the command in REPRODUCE.md with a currently approved runtime/CPU. The pinned manifest SHA256 is `77934c089aa941e191e08c2459d12d6c522596ba4614f81c7fb7850aba6a7094`. The bundle does not contain an installed environment, credentials, private traces, AIDD material or the source input package. Without that original package this is a portable analysis/handoff, not a standalone reproducible dataset. A new data revision requires a deliberate manifest update and scientific reassessment, not bypassing the hash guard.
