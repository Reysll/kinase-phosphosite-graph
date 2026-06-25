# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

Research project building a heterogeneous graph network (HGN) for context-aware kinase-substrate association (KSA) prediction in liver cancer. Two independent layers:

- **`src/`** — graph construction (generic HGN + liver extension)
- **`src_prediction/`** — leave-one-out evaluation, node2vec embeddings, logistic regression scoring

Both trees are organized into functional subpackages (reorganized 2026-06-25):

```
src/
├── core/        config, io_utils, ids, normalize
├── graph/       liver_network, preprocess, sanity_checks
├── builders/    per-dataset edge builders (kinase_substrate, ppi, ptmsigdb, msigdb, ptmcode2, ppase_substrate)
└── runners/     main.py (generic graph), run_liver.py (liver graph)

src_prediction/
├── core/        config, io_utils, graph_loader, parallel_utils, relation_filters, experiment_utils, metrics
├── embeddings/  embeddings.py, embedding_strategy.py
├── data_prep/   evaluation_set, truth_utils, negative_sampling, pair_features, kinase_categories,
│                trial_sampling, find_multi_kinase_sites
├── engine/      leave_one_out.py (LOO trial engine), inference_all.py (predict-all engine),
│                model_scoring.py, scoring.py
├── analysis/    analyze_results, compare_multi_kinase, ksea, run_ksea, build_rank_heatmap, build_rank_ks_test
├── runners/     every run_*.py entry point invoked via `python -m`
└── [root]       9 confirmed-dead/legacy files, intentionally left flat and unmoved — see
                 "Dead/legacy code" below. Do not import from these.
```

**Dead/legacy code** (flat at `src_prediction/` root, NOT in any subpackage — confirmed
unused 2026-06-25 by checking whether their inputs/outputs exist): `compare_experiments.py`,
`compare_multi_kinase_runs.py`, `run_generic_model.py`, `run_baseline_similarity.py`,
`run_generic_multi_kinase.py` (pre-leakage-fix duplicate of `runners/run_multi_kinase_generic.py`
— reads a frozen-trials file that no longer exists), `run_liver_multi_kinase.py` (empty stub),
`make_multi_kinase_fold_set.py`, `make_multi_kinase_fold_set_10.py`, `migrate_trial_index.py`.
Their internal imports still use the old flat `src_prediction.X` paths and will fail if run —
left as-is rather than fixed or deleted, pending a decision on whether to remove them.

**`src/runners/run_liver.py` is not import-safe** — it calls `run(...)` unconditionally at
module level (no `if __name__ == "__main__":` guard), so `import src.runners.run_liver` executes
the full liver-graph build pipeline as a side effect. This is pre-existing behavior (the file is
only ever meant to be run via `python -m src.runners.run_liver`), not introduced by the 2026-06-25
reorg — don't `import` it directly in scripts or smoke tests.

## Python environment

**Always use `.venv312/Scripts/python`** (Python 3.12.7), not bare `python`. The `.venv312/` virtual environment holds all dependencies.

```
.venv312/Scripts/python -m <module>
```

## Key commands

All commands run from the **project root**:

| Task | Command |
|---|---|
| Build generic graph | `.venv312/Scripts/python -m src.runners.main` |
| Build liver graph (all 3 corr types) | `.venv312/Scripts/python -m src.runners.run_liver` |
| Freeze LOO trial set (multi-kinase) | `.venv312/Scripts/python -m src_prediction.runners.run_freeze_multi_kinase_trials` |
| Run generic multi-kinase — node2vec (6,909 trials) | `.venv312/Scripts/python -m src_prediction.runners.run_multi_kinase_generic` |
| Run liver FC multi-kinase — node2vec (6,909 trials) | `.venv312/Scripts/python -m src_prediction.runners.run_multi_kinase_liver` |
| Run control multi-kinase — node2vec (6,909 trials) | `.venv312/Scripts/python -m src_prediction.runners.run_multi_kinase_control` |
| Run cancer multi-kinase — node2vec (6,909 trials) | `.venv312/Scripts/python -m src_prediction.runners.run_multi_kinase_cancer` |
| **Run all 4 — spectral + categories (Phase 2)** | `.venv312/Scripts/python -m src_prediction.runners.run_all_spectral` |
| Smoke test: node2vec + categories (50 trials) | `.venv312/Scripts/python -m src_prediction.runners.run_test_categorized` |
| Smoke test: spectral + categories (50 trials) | `.venv312/Scripts/python -m src_prediction.runners.run_test_spectral` |
| Compare generic vs liver FC | `.venv312/Scripts/python -m src_prediction.analysis.compare_multi_kinase` |
| Run all supervisor analyses (8 networks: node2vec + spectral) | `.venv312/Scripts/python -m src_prediction.analysis.analyze_results` |
| Run KSEA (kinase activity enrichment) | `.venv312/Scripts/python -m src_prediction.analysis.run_ksea` |
| **Run predict-all inference pass (Phase 3, 12 combos)** | `.venv312/Scripts/python -m src_prediction.runners.run_inference_all` |
| Result 1: kinase x network rank heatmap/clustergram | `.venv312/Scripts/python -m src_prediction.analysis.build_rank_heatmap` |
| Result 2: per-kinase rank-distribution KS-test | `.venv312/Scripts/python -m src_prediction.analysis.build_rank_ks_test` |

**Full pipeline order** (start fresh):
1. Build generic graph
2. Build liver graph ← embeds all 3 correlation types in one file
3. Freeze trial set
4. Run Phase 1 (node2vec): runners.run_multi_kinase_{generic,control,cancer,liver}
5. Run Phase 2 (spectral): runners.run_all_spectral  ← adds category columns to results
6. Analyze: analysis.analyze_results  ← 8-network comparison (Phase 1 + Phase 2)
7. KSEA: analysis.run_ksea  ← kinase activity enrichment on liver FC dataset
8. Run Phase 3 (predict-all): runners.run_inference_all  ← 4 networks x 3 embeddings, no LOO holdout
9. Result 1 + 2: analysis.build_rank_heatmap, analysis.build_rank_ks_test ← consume Phase 3 output

Note: the two 50-trial debug commands previously listed here (`run_generic_model`,
`run_baseline_similarity`) were confirmed dead/stale during the 2026-06-25 reorg — see the
"Dead/legacy code" list above. They predate the embedding-leakage fix and their results are
invalid; use the Phase 1 multi-kinase runners instead.

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

**KSEA formula (Wiredja et al. 2017):** score = (s_bar - p_bar) * sqrt(m) / delta.
sqrt(m) is in the NUMERATOR (fixed 2026-06-23 — an earlier version had m in the
denominator, which penalized large substrate sets; the correct formula rewards
them via lower SEM). With the fix: 17/64 scoreable kinases reach adj_p < 0.05
(was 0/64 before). Top UP-regulated: MAPK12, CDKL5, MAPK13, CSNK2A1. Top
DOWN-regulated: PRKACA, PRKD1, PRKCD, RPS6KA1, PAK1, PRKG1, RPS6KA3, PRKG2.
None of these 12 kinases ever reach top-1/5/10 in the LOO results — hub bias
crowds them out. Predicted-substrate KSEA (using LOO top-1 predictions) is
dominated by hub kinases (TTK/PIK3R1) — uninformative for that comparison.
See outputs_prediction/ksea/ksea_comparison.txt for the full report.

**Phase 3 — predict-all inference pass** (`inference_all.py` + `run_inference_all.py`,
built 2026-06-24/25): scores every (candidate kinase, site) pair across the
3,228 FC-valid liver sites, for all 4 networks x 3 embeddings (node2vec,
spectral, concat) — no LOO holdout. Needed because LOO only ever evaluates the
kinase that was actually held out, on the 2,459 multi-kinase sites; it cannot
produce every kinase's rank at every site (needed for the heatmap/KS-test
below). Leakage avoided via two-tier scoring: known true (kinase, site) edges
are scored by a model with ONLY that edge masked (~736 such edges, restricted
to the 3,228 sites); everything else is scored by one global model trained on
all known edges (nothing to leak for a pair that was never a positive label).
"generic" uses the small GENERIC_NODES/EDGES file (matches existing LOO
scripts) and only covers 391/3,228 (12%) of sites — its site population is
essentially PSP-annotated sites, a different population than the broader
liver-proteomics-detected set. Do NOT swap "generic" to the LIVER graph
restricted to base relations to try to fix this — that fragments the liver
site population into many disconnected protein-site components (sites whose
only edges are the correlation ones being stripped out), which breaks
SpectralEmbeddingStrategy's eigensolver (near-zero eigenvalue multiplicity).
Output: `outputs_prediction/inference_all/{network}_{embedding}/ranks.csv.gz`
(long format: site_node_id, kinase_node_id, score, rank_in_site, is_true_kinase).

**Result 2 KS-test caveat:** most kinases' `rank_in_site` is nearly CONSTANT
across all 3,228 sites within one network (hub bias — a kinase's own
embedding dominates its score far more than the site's does), so even a
1-rank shift between networks gives `ks_statistic` ~= 1.0 trivially.
`build_rank_ks_test.py` reports each kinase's within-network rank std-dev and
an `informative` flag (std_rank > 1.0 in at least one network); only
informative + significant kinases are listed in the summary's highlight list,
though the per-pair CSVs keep every kinase so degenerate cases stay visible.

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

**node2vec runs once, on a PSP-stripped graph.** `src_prediction/engine/leave_one_out.py` builds the embedding graph with all `phosphorylates` edges removed (these are the PhosphoSitePlus KSA edges = the prediction target). Embeddings are therefore blind to kinase-substrate relationships and only reflect structural/pathway context. Each trial re-trains logistic regression using the full KSA edges as labels. Do not move embedding computation inside the trial loop.

**Embedding strategy is modular.** `src_prediction/embeddings/embedding_strategy.py` defines `EmbeddingStrategy` (ABC), `Node2VecStrategy`, `SpectralEmbeddingStrategy`, and `ConcatEmbeddingStrategy` (concatenates the two — built 2026-06-25 since spectral fixes the GLOBAL hub bias but not the WITHIN-CATEGORY bias). To swap in a new method, implement `EmbeddingStrategy.fit(graph) -> Dict[str, np.ndarray]` and pass the new strategy to `run_leave_one_out(embedding_strategy=...)` or `run_inference_all(embedding_strategy=...)`. `ConcatEmbeddingStrategy.directed=True`: node2vec is fit on the directed graph it's given directly; spectral gets an explicit `.to_undirected()` copy, since each sub-strategy needs its own appropriately-directed graph rather than sharing one.

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

## Current run status (last updated 2026-06-25)

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

**Phase 3 complete (2026-06-25)** — predict-all inference pass, 4 networks x 3 embeddings = 12 combos, 3,228 FC-valid sites, ~736 own-edge-masked trials per combo (much cheaper than the 6,909-trial LOO runs — full sweep took ~30 min total):

| Network | Site coverage | Notes |
|---|---|---|
| generic | 391/3,228 (12%) | small GENERIC graph file; PSP-annotated sites only |
| control | 3,228/3,228 (100%) | |
| cancer | 3,228/3,228 (100%) | |
| liver_fc | 3,228/3,228 (100%) | |

Result 1 (heatmap/clustergram) confirms the documented spectral finding visually: column dendrogram clusters cancer+liver_fc tightest, control next, generic as the clear outlier branch. Result 2 (KS-test) surfaced a new finding (see caveat above) — most kinases' ranks are network-invariant baselines, not site-varying; only the `informative`-flagged subset (~260-410/415 depending on embedding) is a meaningful distributional-shift candidate list.

**Next priority:** Send Result 1 (clustergrams) and Result 2 (KS-test, both pairs: control-vs-liver_fc and cancer-vs-control) to Dr. Ayati. Open questions for her: (1) does the generic-graph site-coverage gap (12%) need addressing, or is it an acceptable/expected asymmetry; (2) which KS-test pair she finds more informative now that both exist; (3) interpretation of the `informative`-flag caveat.

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
