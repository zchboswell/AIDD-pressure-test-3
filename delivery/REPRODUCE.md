# Reproduce CEDAR generated analysis

Run from the repository root with the separately supplied original `inputs/` package. The private remote contains generated analysis, not the supplied package or administrative protocol files: automatic approval review excluded those files and the supervisor instructed us not to repeat that operation. `admin/CAPTURE_LIMITATION.json` lists every excluded path. Do not reconstruct missing source data from candidate identifiers or the internet.

Input identity: `inputs/MANIFEST.json` SHA-256 **77934c089aa941e191e08c2459d12d6c522596ba4614f81c7fb7850aba6a7094**; all **30** listed file hashes must match. `analysis/reproduce.py` verifies the manifest and every listed input before and after execution. Source files are read only by these scripts.

Runtime used: Python **3.12.14**, RDKit **2026.03.6**, NumPy **2.5.2**, scipy **1.18.0**, scikit-learn **1.9.0**, pandas **3.0.5**, matplotlib **3.11.1**, PyMOL **3.2.0a**. Some imported packages are not needed by the final pipeline. The shared launcher is `/home/zbos/.local/share/cedar-runtime/20260908/cedar-python`; exact import paths are in `admin/RUNTIME_CHECKS.json`. No installation, proprietary software, GPU or scientific internet access is needed. The launcher accepts script filenames and `-c`, but not stdin `-`; that failed attempt was retained in the progress record.

Obtain a valid compute lease first. The original final check used CPU6, one thread and at most 1 GiB RAM; this is historical evidence, not an ongoing resource permission. Substitute the CPU ID and runtime path approved for a new run. Choose an output directory that **does not exist**:

```sh
taskset -c 6 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /home/zbos/.local/share/cedar-runtime/20260908/cedar-python \
  analysis/reproduce.py --inputs inputs --out /tmp/cedar-new-replay
```

Use `--runtime /approved/path/to/cedar-python` if the launcher is relocated. The launcher and installed packages must exist separately; no environment is bundled. The pipeline is sequential and passes one-thread limits to children. It refuses to overwrite an existing output directory. Per-stage stdout/stderr and exit codes are written to its `support/` directory. These are scientific execution logs, not private native traces.

The stages are:

1. Run `analysis/test_assay_math.py` (five regression methods; checks units, bounds, endpoints, invalid values and no clipping).
2. Run `curate.py` (row-level source comparison, chemistry/assay groups and inventory identity).
3. Run `model_audit.py` (untouched vendor execution; retrospective-feature ablation; fixed chemical/scaffold CV and applicability diagnostics).
4. Run `structure_audit.py` (PDB/SDF identity, graph-mapped original-frame checks, independent RMSD, reordering and separately fitted demonstrations).
5. Run `decision.py` (explicit selection policy, observed/abstained prediction accounting, experimental order, exact decimal material/cost/day calculations and chemical/decision figures).
6. Run `render_pocket.py` (headless PyMOL depiction of the deposited source, no candidate fitting/docking).
7. Run `validate_delivery.py` (independent accounting, contract, sentinel numerical and status checks).

Expected outputs include five required CSVs, `structure_checks.json`, five PNG figures, SVG diagnostics, `support/validation.json`, `support/replay_status.json`, `model_artifacts/` and `structure_support/`. Expected numerical results: 266 curated rows with 122/138/2/4 eligible/context/duplicate/quarantine; 28 inventory decisions; 18 observed/6 abstained/2 duplicate/2 quarantined current-endpoint entries; $17,500 for seven compounds; latest day7; RMSD0.000 and8.774964Å. `support/replay_comparison.json` in the delivered directory records comparison of the independently generated replay against final outputs. OOF CSVs and primary numerical files should match exactly in the pinned runtime; SVG timestamps/IDs and execution durations need not be byte-identical.

The decision report, model interpretation, attribution and reusable handoff are authored review artifacts rather than synthesized text from the numerical pipeline. Preserve them alongside replay outputs, and reassess them when evidence changes. The fixed model uses seed7 and 80 trees without tuning; validation metrics are historical, not a coverage guarantee. Export frame metrics concern the original supplied files before fitting. The PyMOL view depicts only 1H1Q. No new biological assays have been executed.

Repair reconstruction: use `analysis/assay_math.py` as the compatible replacement for downstream imports of `normalize(value, unit, relation='=', endpoint='IC50')`. Leave `inputs/legacy/assay_math.py` unchanged. Do not copy the repaired helper over the legacy file when reproducing the original vendor metric: `model_audit.py` intentionally launches the original vendor script with its original helper. No serialized model is required or supplied.

Replay can fail correctly for a changed manifest, unavailable package, invalid output directory or missing input package; keep the error logs. Low-quality predictions, invalid rows and shifted exports are findings, not pipeline errors. Passing the validator does not establish scientific efficacy or successful experimental delivery.
