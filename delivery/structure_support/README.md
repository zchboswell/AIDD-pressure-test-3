# Structure audit and reusable handoff

The source is PDB **1H1Q**, an X-ray model at 2.5 Å of Thr160-phosphorylated CDK2 / cyclin A, bound to **NU6094**, chemical component **2A6**. The extracted ligand is author chain A, residue 1298; the associated kinase is A and cyclin is B (a second pair is C/D). PDB header and supplied entry JSON establish these identities. This is not a CDK1 complex or a main-series candidate co-crystal.

`structure_checks.json` concerns the original supplied exports. The source SDF's 24 heavy-atom positions are checked against the corresponding PDB 2A6 A1298 atoms. All four complete graph isomorphisms are enumerated, and the smallest original-frame RMSD selected. RDKit `CalcRMS` is passed the explicit graph atom map and checked against a NumPy calculation and unchanged input coordinates. `RenumberAtoms` provides reference-ordered outputs without a fit. Installed interface documentation is retained in `installed_interface.txt`.

The reference export preserves the deposited positions exactly: mapped RMSD **0.000 Å**. Comparing atom rows would incorrectly report **6.152 Å**. The priority export has mapped RMSD **8.775 Å** and is a pure translation of **(+8, −3, +2) Å**. Its original position is not usable for binding-site contacts. A separate alignment of a copy gives essentially zero RMSD; this illustrates why fitting cannot test preservation of the receptor frame. Fitted demonstration SDFs are clearly separate outputs and do not validate a novel pose. Prefer the original deposited coordinates for any future discussion. Require a correctly registered export with recorded transformation provenance before using a vendor pose for contacts.

The informative view is an orthographic projection of actual 3D ligand coordinates, using a shared PCA projection for both panels. Its left panel shows source polar-atom proximity: Glu81 O–N9 **2.572 Å**, Leu83 N–N3 **3.046 Å**, and Leu83 O–N2 **2.587 Å** (exact distances in CSV). These are geometrical measurements consistent with a hinge-binding hypothesis, not assigned hydrogen bonds: no hydrogens, angles or energetic analysis establish them. No selectivity or cellular engagement follows from this single structure. The right panel directly shows the original export displacement. Atom correspondence and all ≤4 Å heavy-atom contacts to chain A are supplied as CSV.

CED-23 is supplied as the OpenFF lig_1h1q reference. The inventory and SDF encode different purine tautomers; their RDKit canonical tautomers match, but exact canonical SMILES do not. Do not erase this chemical-state distinction when assigning hydrogen donors/acceptors. CED-24 is a related purine analog with a sulfonamide; no coordinates for CED-24 are supplied here. For main-series compounds, test activity/selectivity empirically; a structure-guided binding hypothesis would require new candidate-specific evidence.

Replay from the repository root into a fresh output directory:

```sh
taskset -c 7 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLCONFIGDIR=/tmp/cedar3-mpl /home/zbos/.local/share/cedar-runtime/20260908/cedar-python analysis/structure_audit.py --out /tmp/cedar3-structure-replay
```

Use only a currently leased CPU when replaying. RDKit 2026.03.6, NumPy 2.5.2, matplotlib 3.11.1 were used. First render completed in 1.82 s wall time, peak RSS 146972 KiB, one CPU and no GPU. Input SDF produced a warning that it was tagged as 2D despite nonzero Z coordinates; RDKit treated it as 3D, and coordinates were verified against PDB. No input was changed. No commercial tooling or scientific internet was used.

Applicability: atom mapping is required whenever chemically identical exports reorder atoms. Fixed-frame RMSD is required whenever preservation of a receptor relationship is the claim. Canonical tautomer grouping aids identity reconciliation but does not establish the bound tautomer. Reassess mappings if chemistry/stereochemistry changes, if hydrogens are added, or if deposited identifiers/coordinates change.
