---
name: run-experiment
description: Run one or more pipeline stages for the Kinase-Substrate HGN project. Knows the correct venv, module paths, and execution order for every experiment variant. Invoke as /run-experiment <stage> where stage is: generic-graph, liver-graph, freeze-trials, generic-loo, liver-loo, generic-multi, liver-multi, compare, or "all" for the full pipeline.
disable-model-invocation: true
---

The user wants to run a pipeline stage for the Kinase-Substrate HGN project.
$ARGUMENTS contains the stage(s) to run (e.g. "liver-loo", "compare", "all").

**Always use `.venv312/Scripts/python`** (not bare `python`). Run all commands from the project root.

## Stage map

| Argument | Command | Output |
|---|---|---|
| `generic-graph` | `.venv312/Scripts/python -m src.main` | `outputs/nodes.csv.gz`, `outputs/edges.csv.gz` |
| `liver-graph` | `.venv312/Scripts/python -m src.run_liver` | `outputs/Liver_network_nodes.csv.gz`, `outputs/Liver_network_edges.csv.gz` |
| `freeze-trials` | `.venv312/Scripts/python -m src_prediction.run_freeze_multi_kinase_trials` | `outputs_prediction/frozen_trials_multi_kinase.csv.gz` |
| `generic-loo` | `.venv312/Scripts/python -m src_prediction.run_generic_model` | `outputs_prediction/generic_trained_model_debug/` |
| `liver-loo` | `.venv312/Scripts/python -m src_prediction.run_baseline_similarity` | `outputs_prediction/liver_trained_model_fc_corr_debug/` |
| `generic-multi` | `.venv312/Scripts/python -m src_prediction.run_multi_kinase_generic` | `outputs_prediction/generic_multi_kinase/` |
| `liver-multi` | `.venv312/Scripts/python -m src_prediction.run_multi_kinase_liver` | `outputs_prediction/liver_multi_kinase/` |
| `compare` | `.venv312/Scripts/python -m src_prediction.compare_experiments` | `outputs_prediction/comparison_*` |
| `full-run` | see note below | `outputs_prediction/generic_multi_kinase/`, `outputs_prediction/liver_multi_kinase/` |

**Full pipeline order** (`all`): generic-graph → liver-graph → freeze-trials → generic-loo → liver-loo → generic-multi → liver-multi → compare.

### Full 14,596-trial evaluation (`full-run`)

The multi-kinase stages (`generic-multi`, `liver-multi`) are the bottleneck for the complete LOO evaluation over all 14,596 kinase-site edges. These are intended to run on a university HPC cluster when access is available. The local commands are the same:

```
.venv312/Scripts/python -m src_prediction.run_multi_kinase_generic
.venv312/Scripts/python -m src_prediction.run_multi_kinase_liver
```

Current debug results use 10 LOO trials. The `frozen_trials_multi_kinase.csv.gz` file already contains all 6,909 multi-kinase trials — no new trial set generation is needed. Running `full-run` locally is possible but may take hours; check `progress_notes.txt` for the latest timing estimates before kicking it off on a laptop.

## Instructions

1. Parse $ARGUMENTS to identify which stage(s) to run.
   - If empty or `all`, run every stage in order.
   - If a comma-separated list, run them in the correct pipeline order regardless of user order.
2. Before running, check that prerequisites exist:
   - `liver-loo`, `liver-multi`, `compare` require `outputs/Liver_network_edges.csv.gz` to exist.
   - `generic-loo`, `generic-multi` require `outputs/edges.csv.gz` to exist.
   - `compare` requires at least one of the LOO result directories to exist.
   If a prerequisite is missing, tell the user which upstream stage to run first.
3. Run each stage with Bash. Show the command before running it.
4. After each stage completes, report: exit code, wall-clock time, and where output was written.
5. If a stage fails, show the last 30 lines of output and stop — do not proceed to downstream stages.

## Key constraints (do not violate)

- **node2vec runs once before the LOO loop** — do not add per-trial embedding code.
- **Protein co-expression is disabled** — do not re-enable or reference it.
- The 80th percentile FC correlation threshold is the production setting; do not change it without being asked.
