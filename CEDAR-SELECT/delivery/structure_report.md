# Structural qualification and falsification test

**Decision contribution:** retain the old purine structures as mechanism controls; do not rank the spiro series by docking into this receptor pair or claim its CDK1 selectivity mechanism is established. Coordinate analysis confirms different Tyr15 conformations, but fails to establish robust steric exclusion of compound 73 from CDK1.

## Identity, suitability and provenance

The supplied [5NEV](https://www.rcsb.org/structure/5NEV) is CDK2/cyclin A with **compound 73, ligand 72L**, not CDK1 and not compound 72. It supersedes obsolete 5LQE. [5LQF](https://www.rcsb.org/structure/5LQF) is CDK1/cyclin B/CKS2 with **NU6102, compound 3, ligand 4SP**. These are different ligands and different cyclin complexes. Neither structure contains the coverage packet's spiro-lactam aminopyrimidine/pyrazole series. The CCD describes a 7H tautomer for 72L whereas Coxon names compound 73 as 9H; the same connectivity does not validate a particular protonation assignment at 2.97 Å.

The deposited 5NEV model has resolution 2.97 Å, Rfree 0.2612. Its two 72L instances (author chains A/C) have validation RSCC 0.86/0.88, RSR 0.17/0.14 and all 32 dictionary heavy atoms present. Calculated ligand mean B factors are 115.8/129.9 Å² versus 66.6/79.5 Å² for heavy atoms of contacting residues. Thus the pose is informative at scaffold level; exact torsions, short distances and proposed hydrogen geometry are less secure. Occupancies recorded as full do not eliminate disorder uncertainty. In 5LQF, resolution is 2.06 Å, Rfree 0.2539, 4SP RSCC is 0.91 in both copies, and mean ligand B is 45.3/48.5 Å². Do not treat the two receptors as equally well resolved. Validation statistics come from supplied wwPDB full reports, section 6.4; paper Table 4 describes obsolete 5LQE and differs from revised 5NEV statistics.

## Executed coordinate test

[Replay script](../structure_work/analyze_structures.py), [numerical results](../structure_work/results.json), [contact table](../structure_work/contacts.csv), [transplant overlap table](../structure_work/transplant_clashes.csv), and [generated view](../structure_work/structure_comparison.png).

The script analyzes both kinase copies per crystal using author chain numbering. Heavy-atom contacts are defined at ≤4.0 Å. Observed protein sequences are globally aligned, avoiding the CDK1/CDK2 numbering offset after residue 100. CDK1 is fitted to CDK2 with conserved Cα atoms, either across residues 26–284 (168–170 anchors) or around the pocket at residues 30–90 and 120–145 (65 anchors). It transplants the observed CDK2 ligand coordinates without optimization and counts heavy-atom overlaps below 2.2 Å. This threshold identifies conspicuous geometric overlap; it is not a force-field energy or donor/acceptor hydrogen-bond test. Symmetry mates, solvent thermodynamics, relaxation and alternate conformer populations are not modeled. Crystal copies are sensitivity checks, not independent biological replicates.

* The expected hinge contacts recur: 72L N3–Glu81 O is 2.82/3.09 Å; N5–Leu83 O is 2.72/2.63 Å; sulfonamide OH–Asp86 N is 2.89/3.25 Å. Distances support the general purine pose, not selective binding. Protein atoms contacted around His84/Gln85 are often backbone atoms; these sequence differences from CDK1 Ser84/Met85 do not establish a selective side-chain interaction.
* Tyr15 OH moves **10.60–11.46 Å** between the fitted structures across all eight chain/fit combinations; loop Cα RMSD is 1.91–2.37 Å. The loop conformation difference is real in these deposited snapshots and qualitatively agrees with Coxon Figure 2 (visually inspected).
* The stronger hypothesis, that the ligand is geometrically excluded from CDK1, **fails robustness**. Only 3/8 transplant combinations have a <2.2 Å overlap (one atom pair each), versus zero for every native 72L control. The overlap shifts between distal phenyl C22/Ile10 O and sulfonamide OH/Asp86 OD2 as the fitting method changes. Domain fitting has RMSD 2.08–2.24 Å; pocket fitting 1.54–1.59 Å. These substantial differences make rigid-transplant penalties highly alignment dependent. The eight comparisons cannot be counted as eight independent tests.

This is a useful negative computational result: do not use a single fixed CDK1 structure to assert that the bulky purine substituent is impossible there. Conformational free-energy preference remains possible and is not measured by this calculation. The comparison is further confounded by distinct ligands, cyclin partners, crystal packing and resolution. No docking affinity, MM/GBSA, free-energy calculation, redocking qualification or new structural observation is claimed.

## Falsifiable next use

Coxon et al. reported 73 at 44 nM CDK2 and 86 µM CDK1 and proposed an indirectly encoded loop conformational preference; they could not obtain a CDK1–73 complex. These are published results and a mechanism proposal, not current CEDAR measurements. The paper also reports minimal/no cellular growth effect for 73 at the maximum tested 30 µM in its panel, demonstrating why isolated-enzyme discrimination does not establish a cellular window. [Primary article](https://doi.org/10.1021/acs.jmedchem.6b01254).

For an old-series diagnostic, compare authentic 73 and NU6102 in matched CDK2/A, CDK2/E and CDK1/B assays with documented ATP concentrations and inhibition curves, alongside solubility and cellular target engagement. A same-ligand CDK2/CDK1 structural or biophysical comparison would distinguish a transferable binding hypothesis from ligand/context confounding. Absence of a CDK1 crystal is not evidence of zero binding.

For the proposed spiro-series switch, biochemical matched analogues and the observed stereo pair should drive selection before structure-based claims. If one earns progression, obtain its CDK2 co-complex and pursue the CDK1 comparison on that same compound. In the meantime, keep hydroxycyclohexyl configuration fixed when testing the pyrazole O-substituent change; a separately tested stereochemical inversion is a diagnostic perturbation. No evidence here assigns the spiro substituent to the purine C6 pocket, predicts a particular Tyr15 response, or supports a tenfold selectivity improvement. New binding data showing no loop difference for a selective spiro compound would refute transfer of the old mechanism while leaving that compound's biochemical value intact.

## Reproduction, exposure and resources

From repository root:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLCONFIGDIR=CEDAR-SELECT/structure_work/.mplconfig /home/zbos/.local/share/cedar-runtime/20260908/cedar-python CEDAR-SELECT/structure_work/analyze_structures.py
```

Environment: pinned cedar-python; Biopython 1.88, numpy and matplotlib from provided runtime. No installation or environment alteration. Runs are single-process and use declared one-thread numerical limits, not claimed kernel-enforced memory caps. Final process usage is in results.json (~seconds wall/CPU and <200 MiB RSS); two preliminary coordinate runs and a figure-generation rerun also completed, each under two seconds. Read-only pdftotext extraction and two article-page renders were performed. No persistent scientific job or delegate was launched; all processes are stopped at handoff.

Exposure: read only this project's AGENTS/admin rules, common request, supplied structural sources and analogue evidence identities, and ordinary PDF skill. Public retrieval on 2026-09-09 opened only the RCSB 5NEV and 5LQF pages to verify current accession identity; no private data was transmitted. Source file hashes/provenance remain in inputs/structural_sources/SOURCE_METADATA.json. No peer, evaluator or AIDD sources accessed. This is exposed historical-source analysis, not blind validation. Extracted full article/validation text and article-page images are scratch excluded from delivery; the generated molecular view is original coordinate visualization.
