---
name: run-experiment
description: Run one or more pipeline stages for the Kinase-Substrate HGN project. Knows the correct venv, module paths, and execution order for every experiment variant. Invoke as /run-experiment <stage> where stage is: generic-graph, liver-graph, freeze-trials, generic-multi, control-multi, cancer-multi, liver-multi, spectral-all, ksea, inference-all, result1, result2, compare, or "all" for the full pipeline.
disable-model-invocation: true
---

The user wants to run a pipeline stage for the Kinase-Substrate HGN project.
$ARGUMENTS contains the stage(s) to run (e.g. "liver-multi", "compare", "all").

**Always use `.venv312/Scripts/python`** (not bare `python`). Run all commands from the project root.

Both `src/` and `src_prediction/` are organized into functional subpackages
(`core/`, `embeddings/`, `data_prep/`, `engine/`, `analysis/`, `runners/` —
reorganized 2026-06-25). Every stage below is invoked from `runners/` or
`analysis/`.

## Stage map

| Argument | Command | Output |
|---|---|---|
| `generic-graph` | `.venv312/Scripts/python -m src.runners.main` | `outputs/nodes.csv.gz`, `outputs/edges.csv.gz` |
| `liver-graph` | `.venv312/Scripts/python -m src.runners.run_liver` | `outputs/Liver_network_nodes.csv.gz`, `outputs/Liver_network_edges.csv.gz` |
| `freeze-trials` | `.venv312/Scripts/python -m src_prediction.runners.run_freeze_multi_kinase_trials` | `outputs_prediction/frozen_trials_multi_kinase.csv.gz` |
| `generic-multi` | `.venv312/Scripts/python -m src_prediction.runners.run_multi_kinase_generic` | `outputs_prediction/generic_multi_kinase/` |
| `control-multi` | `.venv312/Scripts/python -m src_prediction.runners.run_multi_kinase_control` | `outputs_prediction/control_multi_kinase/` |
| `cancer-multi` | `.venv312/Scripts/python -m src_prediction.runners.run_multi_kinase_cancer` | `outputs_prediction/cancer_multi_kinase/` |
| `liver-multi` | `.venv312/Scripts/python -m src_prediction.runners.run_multi_kinase_liver` | `outputs_prediction/liver_multi_kinase/` |
| `spectral-all` | `.venv312/Scripts/python -m src_prediction.runners.run_all_spectral` | `outputs_prediction/{network}_multi_kinase_spectral/` (all 4 networks) |
| `compare` | `.venv312/Scripts/python -m src_prediction.analysis.compare_multi_kinase` | `outputs_prediction/comparison_multi_kinase_*` |
| `analyze` | `.venv312/Scripts/python -m src_prediction.analysis.analyze_results` | `outputs_prediction/analysis/` (8-network summary) |
| `ksea` | `.venv312/Scripts/python -m src_prediction.analysis.run_ksea` | `outputs_prediction/ksea/` |
| `inference-all` | `.venv312/Scripts/python -m src_prediction.runners.run_inference_all` | `outputs_prediction/inference_all/{network}_{embedding}/` (12 combos) |
| `result1` | `.venv312/Scripts/python -m src_prediction.analysis.build_rank_heatmap` | `outputs_prediction/analysis/rank_heatmap_*` |
| `result2` | `.venv312/Scripts/python -m src_prediction.analysis.build_rank_ks_test` | `outputs_prediction/analysis/ks_test_*` |
| `full-run` | see note below | all of the above |

**Full pipeline order** (`all`): generic-graph → liver-graph → freeze-trials →
generic-multi → control-multi → cancer-multi → liver-multi → spectral-all →
compare → analyze → ksea → inference-all → result1 → result2.

### Timing notes

- `generic-multi`/`control-multi`/`cancer-multi`/`liver-multi` (Phase 1, node2vec) and
  `spectral-all` (Phase 2, spectral, runs all 4 networks sequentially) are each
  6,909-trial LOO runs — tens of minutes to ~1 hour each. `spectral-all` alone took
  ~186 min for all four networks in the last full run.
- `inference-all` (Phase 3, predict-all) runs all 4 networks x 3 embeddings = 12
  combos, but each combo only needs ~736 masked trials (not 6,909), so the full
  sweep took ~30 min total in the last run — much cheaper than it looks.
- Check `current_state_of_project.txt` for the latest timing numbers before kicking
  off a long stage on a laptop.

## Instructions

1. Parse $ARGUMENTS to identify which stage(s) to run.
   - If empty or `all`, run every stage in order.
   - If a comma-separated list, run them in the correct pipeline order regardless of user order.
2. Before running, check that prerequisites exist:
   - `liver-graph` requires `outputs/edges.csv.gz` (built by `generic-graph`).
   - `control-multi`, `cancer-multi`, `liver-multi`, `spectral-all` require
     `outputs/Liver_network_edges.csv.gz` to exist.
   - `generic-multi`, `control-multi`, `cancer-multi`, `liver-multi`, `spectral-all`
     require `outputs_prediction/frozen_trials_multi_kinase.csv.gz` (built by `freeze-trials`).
   - `compare` and `analyze` require at least one of the LOO result directories to exist.
   - `ksea` requires `data/liver/LiverCancer_ProtExp_Phospho_casecntrl.xlsx` (raw data, not a
     pipeline output) and benefits from `spectral-all` having run first (for the
     three-way KSEA comparison against predicted substrate sets).
   - `inference-all` requires `outputs_prediction/liver_positive_kinase_site_edges.csv.gz`
     (built by `generic-graph` + `liver-graph`, no `freeze-trials` needed — it uses the
     full known-edge set, not the frozen multi-kinase trial set).
   - `result1` and `result2` require `inference-all` to have completed for at least
     one network/embedding combo.
   If a prerequisite is missing, tell the user which upstream stage to run first.
3. Run each stage with Bash. Show the command before running it.
4. After each stage completes, report: exit code, wall-clock time, and where output was written.
5. If a stage fails, show the last 30 lines of output and stop — do not proceed to downstream stages.

## Key constraints (do not violate)

- **node2vec/spectral embeddings run once before the trial loop** (both in LOO and in
  the Phase 3 predict-all engine) — do not add per-trial embedding code.
- **Protein co-expression is disabled** — do not re-enable or reference it.
- The 80th percentile FC correlation threshold is the production setting; do not change it without being asked.
- **Do not `import` `src.runners.run_liver` directly** (e.g. for a smoke test) — it
  calls `run(...)` unconditionally at module level (no `__main__` guard), so importing
  it executes the full liver-graph build pipeline as a side effect. Always invoke it
  via `python -m src.runners.run_liver`.
- Nine files at the root of `src_prediction/` (`compare_experiments.py`,
  `compare_multi_kinase_runs.py`, `run_generic_model.py`, `run_baseline_similarity.py`,
  `run_generic_multi_kinase.py`, `run_liver_multi_kinase.py`, `make_multi_kinase_fold_set.py`,
  `make_multi_kinase_fold_set_10.py`, `migrate_trial_index.py`) are confirmed dead/legacy
  code, left unmoved and with stale imports. Never invoke or import these.
