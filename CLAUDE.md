# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

Research project building a heterogeneous graph network (HGN) for context-aware kinase-substrate association (KSA) prediction in liver cancer. Two independent layers:

- **`src/`** — graph construction (generic HGN + liver extension)
- **`src_prediction/`** — leave-one-out evaluation, node2vec embeddings, logistic regression scoring

## Python environment

**Always use `.venv312/Scripts/python`** (Python 3.12.7), not bare `python`. The `.venv312/` virtual environment holds all dependencies.

```
.venv312/Scripts/python -m <module>
```

## Key commands

All commands run from the **project root**:

| Task | Command |
|---|---|
| Build generic graph | `.venv312/Scripts/python -m src.main` |
| Build liver graph | `.venv312/Scripts/python -m src.run_liver` |
| Freeze LOO trial set (multi-kinase) | `.venv312/Scripts/python -m src_prediction.run_freeze_multi_kinase_trials` |
| Run generic LOO (50-trial debug) | `.venv312/Scripts/python -m src_prediction.run_generic_model` |
| Run liver LOO (50-trial debug) | `.venv312/Scripts/python -m src_prediction.run_baseline_similarity` |
| Run generic multi-kinase (6,909 trials) | `.venv312/Scripts/python -m src_prediction.run_multi_kinase_generic` |
| Run liver multi-kinase (6,909 trials) | `.venv312/Scripts/python -m src_prediction.run_multi_kinase_liver` |
| Compare debug generic vs liver | `.venv312/Scripts/python -m src_prediction.compare_experiments` |
| Compare multi-kinase generic vs liver | `.venv312/Scripts/python -m src_prediction.compare_multi_kinase` |
| Run all supervisor analyses | `.venv312/Scripts/python -m src_prediction.analyze_results` |

**Full pipeline order** (start fresh): generic graph → liver graph → freeze trial set → run generic/liver → compare → analyze.

**Analysis outputs** go to `outputs_prediction/analysis/`: `performance_summary.txt`, `top1_frequency_summary.txt`, `ttk_deepdive.txt`, `ranking_shifts.txt`.

## Data layout

```
data/raw/            <- 7 source datasets (Kinase_Substrate_Dataset.tsv, PPase, PTMsigDB, MsigDB, PPI, gene_protein, PTMcode2)
data/liver/          <- LiverCancer_ProtExp_Phospho_casecntrl.xlsx
data/processed/      <- cached preprocessed files (auto-generated, safe to delete and regenerate)
outputs/             <- graph outputs (nodes.csv.gz, edges.csv.gz, Liver_network_*.csv.gz)
outputs_prediction/  <- frozen trials, candidate kinases, experiment results
```

## Node ID format

All node IDs follow:
- Proteins: `PROTEIN:GENENAME` (e.g., `PROTEIN:AKT1`)
- Phosphosites: `SITE:GENENAME_pSXXX` (e.g., `SITE:AKT1_pS473`)

## Critical architectural decisions

**node2vec runs once before the LOO loop** (not per trial). `leave_one_out.py` pre-computes embeddings from the full graph once, then each trial only masks one training edge and re-trains the logistic regression. Do not move embedding computation inside the trial loop — this was intentionally redesigned in May 2026 for performance (50-trial run: 24 h -> <5 min).

**Protein co-expression edges are disabled.** `liver_network.py` has logic for protein correlation edges but it is turned off — the input data has only 3 usable rows. Do not re-enable without new data.

**80th percentile FC correlation threshold** is the working setting for liver fold-change correlation edges (`site_corr_fc_pos`, `site_corr_fc_neg`). 85th percentile was tested and performed worse.

**Clique size limits** in `src/config.py`: `MAX_SITE_CLIQUE = 50`, `MAX_PROTEIN_CLIQUE = 100`. These guard against exploding edge counts from fully-connected pathway cliques.

## Current experiment scale

- Generic graph: ~13,318 nodes, ~1.1 M edges
- Liver graph: ~16,605 nodes, ~2.1 M edges
- Candidate kinases: 420
- Frozen LOO trial set (multi-kinase): 6,909 trials across 2,459 multi-kinase sites
- 50-trial results: liver improves sparse sites; generic retains advantage for hub sites (AKT1-S473, GSK3B-S9, etc.)

## Progress tracking

Three files record the evolving state of the project — read these before making changes:
- `progress_notes.txt` — technical decisions and benchmarks
- `current_state_of_project.txt` — full experimental results and next steps
- `PROJECT_PIPELINE_EXPLANATION.txt` — detailed walkthrough of every module and how they connect
