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
| Run generic multi-kinase — node2vec (6,909 trials) | `.venv312/Scripts/python -m src_prediction.run_multi_kinase_generic` |
| Run liver FC multi-kinase — node2vec (6,909 trials) | `.venv312/Scripts/python -m src_prediction.run_multi_kinase_liver` |
| Run control multi-kinase — node2vec (6,909 trials) | `.venv312/Scripts/python -m src_prediction.run_multi_kinase_control` |
| Run cancer multi-kinase — node2vec (6,909 trials) | `.venv312/Scripts/python -m src_prediction.run_multi_kinase_cancer` |
| **Run all 4 — spectral + categories (Phase 2)** | `.venv312/Scripts/python -m src_prediction.run_all_spectral` |
| Smoke test: node2vec + categories (50 trials) | `.venv312/Scripts/python -m src_prediction.run_test_categorized` |
| Smoke test: spectral + categories (50 trials) | `.venv312/Scripts/python -m src_prediction.run_test_spectral` |
| Compare generic vs liver FC | `.venv312/Scripts/python -m src_prediction.compare_multi_kinase` |
| Run all supervisor analyses (8 networks: node2vec + spectral) | `.venv312/Scripts/python -m src_prediction.analyze_results` |
| Run KSEA (kinase activity enrichment) | `.venv312/Scripts/python -m src_prediction.run_ksea` |

**Full pipeline order** (start fresh):
1. Build generic graph
2. Build liver graph ← embeds all 3 correlation types in one file
3. Freeze trial set
4. Run Phase 1 (node2vec): run_multi_kinase_{generic,control,cancer,liver}
5. Run Phase 2 (spectral): run_all_spectral  ← adds category columns to results
6. Analyze: analyze_results.py  ← 8-network comparison (Phase 1 + Phase 2)
7. KSEA: run_ksea.py  ← kinase activity enrichment on liver FC dataset

**Phase 2 output columns** (spectral results, in addition to all Phase 1 columns):
- `held_out_kinase_category` — poor / average / rich (PSP substrate count)
- `adjusted_held_out_rank_in_category` — rank within category (denom ~120-166, not 420)
- `top1_poor` / `top1_average` / `top1_rich` — top-1 kinase within each category

**Analysis outputs** go to `outputs_prediction/analysis/` (8 networks: Phase 1 + Phase 2):
- `performance_summary.txt` — four node2vec networks side-by-side
- `top1_frequency_summary.txt` — most predicted kinases per network (node2vec)
- `per_kinase_topk.txt` / `.csv` — per-kinase top-1/5/10% across all networks (node2vec)
- `ttk_deepdive.txt` — TTK case study (node2vec + spectral)
- `ranking_shifts.txt` — rank improvements/worsenings node2vec generic vs liver FC
- `phase_comparison.txt` — Phase 1 vs Phase 2 median rank / top-20 / improvement%
- `category_ranks.txt` — per-category (poor/average/rich) breakdown (spectral)
- `ranking_shifts_spectral.txt` — rank improvements/worsenings spectral generic vs liver FC

**KSEA outputs** go to `outputs_prediction/ksea/`:
- `ksea_psp.csv` — KSEA scores using known PSP substrates
- `ksea_generic_spectral.csv` — KSEA using generic spectral top-1 predictions
- `ksea_liver_spectral.csv` — KSEA using liver FC spectral top-1 predictions
- `ksea_three_way_comparison.csv` — merged PSP vs generic vs liver scores
- `ksea_comparison.txt` — human-readable report

**KSEA formula (Wiredja et al. 2017):** score = (s_bar - p_bar) / (m * delta).
m in denominator penalizes kinases with many substrates. No significant hits after BH-FDR
(adj_p ~0.5 throughout) — dataset underpowered for this formula. Top-ranked: MAPK12
(n=3, score=1.08), CDKL5 (n=3, score=1.05) — few substrates with very high FC.
CDK1 (n=61) and CSNK2A1 (n=35) penalized by large m despite consistent upregulation.
Predicted-substrate KSEA is dominated by hub kinases (TTK/PIK3R1) — uninformative.
For richer KSEA: need a "predict all" inference pass over all liver phosphosites.

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

## Current run status (2026-06-17)

**Phase 1 complete** — node2vec, 6,909 trials each:

| Experiment | adj median rank |
|---|---|
| generic_multi_kinase | 188.0 |
| control_multi_kinase | 166.5 |
| cancer_multi_kinase | 174.0 |
| liver_multi_kinase | 165.0 |

**Phase 2 complete** — SpectralEmbeddingStrategy + category ranking, 6,909 trials each:

| Experiment | adj median rank | adj top-20 | global top-1 hub |
|---|---|---|---|
| generic_multi_kinase_spectral | 79.0 | 10.6% | TTK 92.9% |
| control_multi_kinase_spectral | 74.0 | 22.1% | CDK1 57% / PIK3R1 43% |
| cancer_multi_kinase_spectral | 74.0 | 22.1% | PIK3R1 99.3% |
| liver_multi_kinase_spectral | 74.0 | 22.0% | PIK3R1 100% |

**Key Phase 2 findings:**
- Spectral halves median rank vs node2vec across all networks (~57% improvement).
- Under spectral, Control = Cancer = Liver FC (all median 74). Correlation TYPE no longer
  differentiates — the signal is in having correlation edges, not in which type.
- Control is unique: top-1 splits CDK1/PIK3R1 (57%/43%) — most biologically plausible.
- Hub bias persists (within-category: PIK3R1→poor, TGFBR2→average, CDK1→rich at 100%),
  but hub kinases shifted from VRK family to cancer-relevant nodes (CDK1, PI3K, TGF-β).

**Next priority:** Draft email to Dr. Ayati with Phase 2 results + KSEA findings. Then consider a "predict all" inference pass for richer KSEA substrate sets.

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
