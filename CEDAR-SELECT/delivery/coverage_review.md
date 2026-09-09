# Coverage and executable work-order decision

Propose **O2 ($2,000) and O1 ($1,400)**, retain the already authorized S2 CDK1 endpoint, and defer O3/O4. This is a recommendation within a simulated scenario, not a purchase, assay execution, or newly measured result. Its scientific purpose is to recover a failed selectivity endpoint and assess independent-preparation reproducibility before spending on a small potency gain that may weaken the CDK1 margin.

## Evidence that changes the priority

The supplied `inputs/coverage/analogue_evidence.json` attributes eight analogues to US20240360137A1, with exact source-reported values in a single paired biochemical study. The figures below are calculated directly from its original micromolar values, not from the retrieval responses or work-order states.

| Slot / public example | Alkoxy change | CDK2 IC50 nM | CDK1 IC50 nM | CDK1/CDK2 |
|---|---|---:|---:|---:|
| S1 / 12 | Methoxy | 20 | 2240 | 112 |
| S2 / 23 | Ethoxy | 15 | 2430 | 162 |
| S3 / 39 | 2-Fluoroethoxy | 7 | 1200 | 171.4 |
| S4 / 41 | 2,2-Difluoroethoxy | 5 | 730 | 146 |
| — / 56 | 2,2,2-Trifluoroethoxy | 10 | 800 | 80 |
| — / 69 | Difluorocyclopropoxy stereoisomer | 6 | 800 | 133.3 |
| — / 70 | Opposite side-chain stereoisomer | 6 | 820 | 136.7 |
| — / 77 | Hydroxy | 121 | 5820 | 48.1 |

Adding a first fluorine to S2 improves reported CDK2 potency 2.14-fold but also improves unwanted CDK1 potency 2.03-fold: discrimination changes only 1.058-fold. A second fluorine improves CDK2 potency 1.4-fold while CDK1 potency improves 1.64-fold, reducing the ratio to 0.852 of S3. Those small ratio differences are not established improvements without replicates. Thus prioritize a prospective S2/S3 comparison, retain S1 as a reproducibility anchor, and do not buy S4 merely because it is the most potent. Hydroxy loses the 100 nM preference; the side-chain stereoisomers are a potential stereochemical diagnostic but equal rounded CDK2 values do not prove equivalent binding.

These comparisons use CDK2/E1 versus CDK1/B1, both 1 mM ATP and the same RB1 TR-FRET readout, with different enzyme concentrations (100 versus 10 pM). They establish conditional source-assay ratios, not cellular selectivity. Work-order context-A/B names have no supplied mapping to those public contexts; historical values cannot fill missing work-order result cells or be paired with a new endpoint to fabricate a new ratio. Request raw curves, QC, numerical estimates and construct/ATP metadata for completed and pending work before making a paired comparison. An O2-only CDK1 result does not establish a fresh S3 ratio without a documented comparable CDK2 measurement.

## Retrieval reconciliation

| Response | What the archive supports | Remaining action / forbidden inference |
|---|---|---|
| R1 | Complete, scope-matched fictional CDK2 IC50 index hit | Obtain primary measurement and assay context; index has no numerical potency. |
| R2 | Executed unsupported `target` key; returned CDK1 P06493 | Repeat using `target_accession=P24941`, verify returned target. HTTP 200 does not validate intent. |
| R3 | An assay entity, numeric ID 771 | Repeat for compound namespace; compound:771 and assay:771 are different identities. |
| R4 | Page 1 of 2; archive lacks next page | Retrieve/reconcile page 2 before claiming coverage; one returned record is not exhaustive. |
| R5 | Complete empty compound:999 lookup | Scoped index absence only; neither inactivity nor universal evidence absence. |
| R6 | HTTP 503; null records | Retrieval failed and coverage is unknown; retry or independently source-check. |
| R7 | Empty form-keyed mechanism lookup | Mechanisms are parent-indexed; inspect the explicit parent link, not an assay-label transfer. |
| R8 | Parent:771 simulated CDK2 annotation | Resolves parent annotation, without numerical assay or automatic salt-form equivalence. |
| R9 | Wrong endpoint: percent inhibition at 10 µM | Recover actual IC50 concentration-response record; never convert this metadata row to an IC50. |

The contract permits R7/R8's separate entity-keyed mechanism operation; their keys are not erroneous activity filters. No envelope is eligible numerical IC50 evidence. These are inert records: no claim of a real provider defect, live search capability, or successful service retry is made.

## Exact budget and material reconciliation

Financial identity is `obligation_id`, not export row. F1/F2 are one paid OB-17 ($1,800); F3/F4 are one reserved OB-23 ($2,200). The **$10,000 earmark therefore has $6,000 free**, not $2,000 or $10,000. The earmark is included in the proposed $120,000 tranche; $110,000 lies outside it.

| Charge | USD | Treatment |
|---|---:|---|
| OB-17 | 1,800 | Paid once; includes prior S3 procurement |
| OB-23 | 2,200 | Preserve open reservation; do not recharge or silently cancel |
| O2 | 2,000 | New proposed QC-controlled failed CDK1 repeat on S3 |
| O1 | 1,400 | New proposed deliberate independent preparation/day S1 CDK2 repeat |
| Unallocated earmark | 2,600 | Hold pending actual results; no invented assay-price commitment |
| **Total earmark** | **10,000** | **$7,400 paid/reserved/proposed + $2,600 unallocated** |

O1 deliberately repeats valid W1 to answer reproducibility, so completed status alone is no reason to reject it. Conversely W3 was attempted but failed QC, so it has no eligible CDK1 result and O2 fills a real endpoint gap. W2 has only its CDK2 child complete; the CDK1 child remains pending to day 3. O3 duplicates that pending work and is not justified now. O4 has unavailable cost and duration and no verified material; it cannot enter an executable allocation merely because a catalogue entry or retrosynthesis exists.

| Vial | Unique total usable mg | Reserved within total mg | Free at cutoff mg | Planned usage mg | Remaining after plan mg |
|---|---:|---:|---:|---:|---:|
| V1 / S1 | 5 | 0 | 5 | O1: 1 | 4 |
| V2 / S2 | 6 | 2 | 4 | Pending W2: 2 | 4 |
| V3 / S3 | 4 | 0 | 4 | O2: 2 | 2 |

I2/I3 are duplicate exports of V2, not 12 mg. Reserved material is part of the stated 6 mg, not additional material. M3's received 4 mg is V3 already counted, not another vial; its cost is already OB-17. M4 records one simulated historical S1 preparation (41% yield, simulated NMR/LCMS identity checks); it establishes neither current stock nor S4 synthesis. M1 is a proposed S4 route with no experimental procedure/yield; M2 is only a fictional listing with quantity, price, lead time and availability unverified. None adds real chemical feasibility evidence.

## Schedule, decision gates and replay

Assume conservatively that existing pending W2 occupies lane 1 on days 0–3. Run O2 in lane 2 on days 0–3; run O1 in lane 1 on days 3–5. This respects two concurrent jobs and completes by day 5 versus the day-8 limit, with three days slack. New work starts only on prospective authorization; the schedule is a feasible simulated proposal, not a provider booking. If the as-of day has advanced or existing lane occupancy differs, reschedule before authorization rather than declaring these dates guaranteed.

At day 3, inspect pending S2 and recovered S3 numerical curves/QC. If either fails, stop chemical escalation and use the held funds only after a costed revised plan. At day 5, assess S1 day/preparation reproducibility against W1 after retrieving its actual numeric result and metadata. A repeat alone without the original numerical result cannot assess reproducibility. A greater-than-threefold shift is a pragmatic investigation trigger, not a claimed statistical threshold. The decisive prospective comparison remains matched CDK2/CDK1 dose responses for S2 and S3 with independent preparations and replicate uncertainty; fund any missing measurements explicitly elsewhere in the $120,000 parent plan. Reverse preference for S3 if its paired selectivity falls below 10-fold, if biochemical CDK2 exceeds 100 nM, or if S2 supplies a reproducible superior selectivity/cellular-window result. Similar S2/S3 ratios favor the simpler ethoxy reference rather than implying a fluorination advantage.

Run `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/zbos/.local/share/cedar-runtime/20260908/cedar-python CEDAR-SELECT/coverage_work/reconcile.py` from the project. It checks duplicate consistency, exact costs, material constraints, option durations, lane occupancy and deadline, and recreates `delivery/coverage_reconciliation.json` with all nine interpretations and the original-value calculations. Execution passed; computation used standard-library Python, one bounded foreground process, no installs, no service calls and no delegates. That process exited; no scientific jobs remain from this subtask.
