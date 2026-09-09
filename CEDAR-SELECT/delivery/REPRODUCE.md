# Replay and custody

Use the repository root and pinned launcher `/home/zbos/.local/share/cedar-runtime/20260908/cedar-python`. Python3.12.14, RDKit2026.03.6, numpy2.5.2, scipy1.18.0, scikit-learn1.9.0, pandas3.0.5, matplotlib3.11.1 and Biopython1.88 were available. Runtime is external and deliberately not committed. No package was installed. Other environments may change byte-level outputs; this is a pinned replay, not a cross-version numerical guarantee.

The clean-output assembly copies only the new supplied inputs, explicit Python scripts, and original inventory CSV into a new destination. It then executes prediction, censor sensitivity, design checks, coverage reconciliation, structural checks and tranche arithmetic in sequence, each with a180-second timeout and numerical threads1. It does not overwrite delivered scientific files. Use a destination that does not exist:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/zbos/.local/share/cedar-runtime/20260908/cedar-python CEDAR-SELECT/analysis/replay.py --out /tmp/cedar-select-replay-unique
```

`replay_check.json` is written inside that new destination. All13 delivered CSVs must be byte-identical, three deterministic JSONs equal, and structural scientific results equal; runtime counters are excluded. The completed pre-freeze replay receipt is `evidence/replay_check.json`. The expected RDKit invalid-ring message comes from the preserved invalid original-inventory `C1CC` record. It is excluded, not repaired. No network, hidden outcomes, original predecessor model, paid service, background job or installed environment copy is needed.

For individual runs, scripts assume repository-root working directory and write their documented relative paths; use the isolated assembly when preserving frozen outputs. `predictive.py --out` and `censor_sensitivity.py --base --out` also support separate output paths. Other scripts have fixed assembly-relative paths. Do not run them over the frozen package to append new analyses without owner authorization.

`delivery/SHA256SUMS` binds the active package's inputs, scripts, outputs and handoff documents; it excludes itself, ignored caches/replay scratch and unrelated predecessor files. Run from repository root:

```sh
sha256sum -c CEDAR-SELECT/delivery/SHA256SUMS
```

The original17-file input manifest is separately retained as `admin/R005_INPUT_MANIFEST.json`; input bytes were rechecked at completion. The final commit binds that manifest, current conduct and the entire tracked repository. Full article/validation text extracts and screenshots are ignored scratch; original supplied PDFs remain under inputs. Structural numerical outputs and original generated images are retained, with no duplicated publication text.

The frozen scientific handoff is immutable. Any next experimental results require a separately authorized continuation. This replay demonstrates computational reproduction, not validation of the proposed biological outcome.
