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
| Build liver graph (all 3 corr types) | `.venv312/Scripts/python -m src.run_liver` |
| Freeze LOO trial set (multi-kinase) | `.venv312/Scripts/python -m src_prediction.run_freeze_multi_kinase_trials` |
| Run generic LOO (50-trial debug) | `.venv312/Scripts/python -m src_prediction.run_generic_model` |
| Run liver FC LOO (50-trial debug) | `.venv312/Scripts/python -m src_prediction.run_baseline_similarity` |
| Run generic multi-kinase (6,909 trials) | `.venv312/Scripts/python -m src_prediction.run_multi_kinase_generic` |
| Run liver FC multi-kinase (6,909 trials) | `.venv312/Scripts/python -m src_prediction.run_multi_kinase_liver` |
| Run control multi-kinase (6,909 trials) | `.venv312/Scripts/python -m src_prediction.run_multi_kinase_control` |
| Run cancer multi-kinase (6,909 trials) | `.venv312/Scripts/python -m src_prediction.run_multi_kinase_cancer` |
| Compare generic vs liver FC | `.venv312/Scripts/python -m src_prediction.compare_multi_kinase` |
| Run all supervisor analyses (4 networks) | `.venv312/Scripts/python -m src_prediction.analyze_results` |

**Full pipeline order** (start fresh):
1. Build generic graph
2. Build liver graph ← now embeds all 3 correlation types in one file
3. Freeze trial set
4. Run all four LOO experiments (generic, control, cancer, liver FC)
5. Compare + analyze

**Analysis outputs** go to `outputs_prediction/analysis/`:
- `performance_summary.txt` — all four networks side-by-side (held-out, adjusted, best-true)
- `top1_frequency_summary.txt` — most predicted kinases per network
- `per_kinase_topk.txt` / `.csv` — per-kinase top-1/5/10% across all networks
- `ttk_deepdive.txt` — TTK case study
- `ranking_shifts.txt` — biggest rank improvements/worsenings (generic vs liver FC)

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

**node2vec runs once, on a PSP-stripped graph.** `leave_one_out.py` builds the embedding graph with all `phosphorylates` edges removed (these are the PhosphoSitePlus KSA edges = the prediction target). Embeddings are therefore blind to kinase-substrate relationships and only reflect structural/pathway context. Each trial re-trains logistic regression using the full KSA edges as labels. Do not move embedding computation inside the trial loop.

**Embedding strategy is modular.** `src_prediction/embedding_strategy.py` defines `EmbeddingStrategy` (ABC) and `Node2VecStrategy` (current). To swap in graph contrastive learning or spectral methods, implement `EmbeddingStrategy.fit(graph) -> Dict[str, np.ndarray]` and pass the new strategy to `run_leave_one_out(embedding_strategy=...)`.

**Four networks, one liver graph file.** `Liver_network_edges.csv.gz` now contains all three correlation edge types. LOO run scripts select which type via `allowed_relations`. No need to rebuild the graph per experiment.

| Network | Correlation edges used |
|---|---|
| Generic | none — `GENERIC_BASE_RELATIONS` only |
| Control | `site_corr_ctrl_pos/neg` (18 healthy samples) |
| Cancer | `site_corr_cancer_pos/neg` (18 tumor samples) |
| Liver FC | `site_corr_fc_pos/neg` (tumor/mean_control ratio) |

**adjusted_held_out_rank** is the primary per-trial rank metric. Before ranking the held-out kinase, all OTHER known true kinases for that site are removed from the scored list. This prevents other co-kinases (whose edges were NOT held out) from artificially inflating the held-out kinase's rank. The `held_out_kinase_rank` column (unadjusted) is still saved for backward compatibility.

**Protein co-expression edges are disabled.** `liver_network.py` has logic for protein correlation edges but it is turned off — the input data has only 3 usable rows. Do not re-enable without new data.

**80th percentile correlation threshold** applies to all three correlation types (FC, control, cancer). 85th percentile was tested on FC and performed worse.

**Clique size limits** in `src/config.py`: `MAX_SITE_CLIQUE = 50`, `MAX_PROTEIN_CLIQUE = 100`. These guard against exploding edge counts from fully-connected pathway cliques.

## Current experiment scale

- Generic graph: ~13,318 nodes, ~1.1 M edges
- Liver graph: ~16,605 nodes, ~2.1 M edges (contains all 3 correlation edge types)
- Candidate kinases: 420
- Frozen LOO trial set (multi-kinase): 6,909 trials across 2,459 multi-kinase sites
- n_jobs_outer: 8 (ThreadPoolExecutor); Node2Vec Word2Vec workers: 4

## Current run status (2026-05-27)

All four networks complete. Analysis regenerated.

| Experiment | adjusted median rank |
|---|---|
| generic_multi_kinase | 188.0 |
| control_multi_kinase | 166.5 |
| cancer_multi_kinase | 174.0 |
| liver_multi_kinase | 165.0 ← best |

**Critical finding:** Adding correlation edges causes severe hub bias in node2vec.
In Liver FC, ULK1 is predicted top-1 for 100% of 6,909 trials. In Control, VRK2
for 94.8%. Cancer: VRK3 for 95.5%. The improved rank METRICS are real but top-1
predictions are biologically meaningless due to degree concentration.

**Next priority:** Swap `Node2VecStrategy` for a degree-normalised embedding via
the `EmbeddingStrategy` ABC (no other code changes needed).

## Progress tracking

Three files record the evolving state of the project — read these before making changes:
- `progress_notes.txt` — technical decisions and benchmarks
- `current_state_of_project.txt` — full experimental results and next steps
- `PROJECT_PIPELINE_EXPLANATION.txt` — detailed walkthrough of every module and how they connect

---

## Coding guidelines (Karpathy Skills)

Behavioral guidelines to reduce common LLM coding mistakes.

**Tradeoff:** These guidelines bias toward caution over speed. For trivial tasks, use judgment.

### 1. Think Before Coding

**Don't assume. Don't hide confusion. Surface tradeoffs.**

Before implementing:
- State your assumptions explicitly. If uncertain, ask.
- If multiple interpretations exist, present them - don't pick silently.
- If a simpler approach exists, say so. Push back when warranted.
- If something is unclear, stop. Name what's confusing. Ask.

### 2. Simplicity First

**Minimum code that solves the problem. Nothing speculative.**

- No features beyond what was asked.
- No abstractions for single-use code.
- No "flexibility" or "configurability" that wasn't requested.
- No error handling for impossible scenarios.
- If you write 200 lines and it could be 50, rewrite it.

Ask yourself: "Would a senior engineer say this is overcomplicated?" If yes, simplify.

### 3. Surgical Changes

**Touch only what you must. Clean up only your own mess.**

When editing existing code:
- Don't "improve" adjacent code, comments, or formatting.
- Don't refactor things that aren't broken.
- Match existing style, even if you'd do it differently.
- If you notice unrelated dead code, mention it - don't delete it.

When your changes create orphans:
- Remove imports/variables/functions that YOUR changes made unused.
- Don't remove pre-existing dead code unless asked.

The test: Every changed line should trace directly to the user's request.

### 4. Goal-Driven Execution

**Define success criteria. Loop until verified.**

Transform tasks into verifiable goals:
- "Add validation" → "Write tests for invalid inputs, then make them pass"
- "Fix the bug" → "Write a test that reproduces it, then make it pass"
- "Refactor X" → "Ensure tests pass before and after"

For multi-step tasks, state a brief plan:
```
1. [Step] → verify: [check]
2. [Step] → verify: [check]
3. [Step] → verify: [check]
```

Strong success criteria let you loop independently. Weak criteria ("make it work") require constant clarification.

---

**These guidelines are working if:** fewer unnecessary changes in diffs, fewer rewrites due to overcomplication, and clarifying questions come before implementation rather than after mistakes.
