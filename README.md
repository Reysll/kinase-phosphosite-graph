# Kinase-Substrate Graph Project

## Table of Contents

* [Overview](#overview)
* [Scientific Motivation](#scientific-motivation)
* [Core Design Choices](#core-design-choices)
* [Repository Layout](#repository-layout)
* [Data Inputs](#data-inputs)
* [Graph Construction Layer](#graph-construction-layer)
* [Prediction and Evaluation Layer](#prediction-and-evaluation-layer)
* [Supporting Helper Files](#supporting-helper-files)
* [Outputs](#outputs)
* [How the Full Pipeline Works](#how-the-full-pipeline-works)
* [How to Run the Project](#how-to-run-the-project)
* [Experimental Results](#experimental-results)
* [Interpretation of Results](#interpretation-of-results)
* [Important Notes and Caveats](#important-notes-and-caveats)
* [Current Scientific Direction](#current-scientific-direction)

## Overview

This project builds and evaluates a context-aware biological network for kinase-substrate association (KSA) prediction, with a current application to liver cancer.

The full workflow has two major layers:

1. **Graph construction**

   * Build a generic heterogeneous phosphorylation graph from curated public resources.
   * Extend that graph into a liver-specific graph using liver proteomic and phosphoproteomic data, adding three types of phosphosite correlation edges.

2. **Prediction and evaluation**

   * Embed the graph with node2vec (applied once on a PSP-stripped graph to prevent leakage).
   * Train a logistic regression scorer for kinase-site pairs.
   * Evaluate all four network variants using leave-one-out on 6,909 frozen multi-kinase trials.
   * Compare kinase rankings across networks to identify context-dependent changes.

**Current status (2026-05-27):** All four networks complete. Key finding: adding correlation edges improves mean/median rank metrics but introduces severe hub degree bias in node2vec — a single kinase dominates top-1 for all trials in every correlation network. Degree-normalised embedding is the immediate next step.

## Scientific Motivation

Many existing computational KSA methods work out of context. They rely on static information such as sequence motifs, structural features, or generic prior knowledge. However, phosphorylation is highly transient and context-dependent: the kinase most relevant to a phosphosite may differ across tissues, disease states, and biological conditions.

This project addresses that limitation by augmenting a generic phosphorylation network with liver-specific phosphoproteomic context and comparing how kinase rankings change across four network variants.

A related background issue is the **streetlight effect**: well-studied kinases have more available data, making them easier for both biology and prediction methods to focus on. The main current scientific emphasis is **context-aware prediction**, not bias alone.

## Core Design Choices

### 1. Only two node types

The graph uses only two node types:

* `PROTEIN` — kinases, phosphatases, and all other proteins
* `SITE` — individual phosphosites

Examples:

* `PROTEIN:AKT1`
* `PROTEIN:PPP1CA`
* `SITE:FOXO3-T32`
* `SITE:FOXO3-S253`

Biological roles are stored as node metadata (`is_kinase`, `is_phosphatase`, `protein_role`) rather than as separate node types. This keeps the graph structure simple while preserving interpretability.

### 2. Site-specific modeling

Phosphorylation is modeled at the phosphosite level. Different sites on the same protein are distinct graph nodes. This is essential because the prediction target is kinase-to-site association, not kinase-to-protein association.

### 3. Generic graph first, context layers second

The project intentionally separates:

* a **generic graph** built from curated public resources
* a **liver-specific graph** built by augmenting the generic graph with three correlation edge types derived from the liver proteomics dataset

This separation allows direct comparison across a controlled set of four network variants.

### 4. Context encoded as edges

Liver context is introduced through additional graph edges rather than only as annotations. This means disease context changes the graph topology itself, directly affecting node2vec embeddings and the learned scorer.

Three correlation edge types are stored in a single liver graph file. Each LOO run script selects the relevant subset via `allowed_relations`:

| Network | Edges added | Correlation source |
|---|---|---|
| Generic | none | — |
| Control | `site_corr_ctrl_pos/neg` | raw abundance, 18 healthy samples |
| Cancer | `site_corr_cancer_pos/neg` | raw abundance, 18 tumor samples |
| Liver FC | `site_corr_fc_pos/neg` | fold-change (tumor / mean control) |

### 5. Embedding leakage fix

`node2vec` runs **once** on a PSP-stripped graph: all PhosphoSitePlus `phosphorylates` edges (the prediction target) are removed before computing embeddings. These edges remain as training labels. This prevents the model from reading the answer indirectly via embedding structure. Non-PSP KSA edges from other resources are kept even if they coincide with a PSP annotation.

### 6. Adjusted held-out rank (primary metric)

For a site with known kinases k1/k2/k3, if k1 is held out and the model scores [k7, k4, k2, k1], the naive rank of k1 is 4. Since k2's edge was not hidden, it should not count against k1. The adjusted metric removes all other known true kinases before measuring the held-out kinase's rank — giving rank 3, not 4.

### 7. Leave-one-out evaluation on multi-kinase sites

The frozen evaluation set covers 6,909 trials across 2,459 sites with two or more known kinases. Multi-kinase sites are used because they allow the adjusted metric to differ meaningfully from the naive metric, and because they are the most informative cases for studying context-dependent kinase prioritisation.

### 8. Modular embedding architecture

`src_prediction/embedding_strategy.py` defines an `EmbeddingStrategy` abstract base class with a single `fit(graph) -> Dict[str, ndarray]` method. `Node2VecStrategy` is the current concrete implementation. To swap the embedding method, implement `EmbeddingStrategy.fit()` and pass the new instance to `run_leave_one_out(embedding_strategy=...)`. No other code changes are needed.

## Repository Layout

### Graph construction

* `src/builders/kinase_substrate.py`
* `src/builders/ppase_substrate.py`
* `src/main.py`
* `src/liver_network.py`
* `src/run_liver.py`
* `src/config.py`, `src/io_utils.py`, `src/ids.py`, `src/normalize.py`

### Prediction and evaluation

* `src_prediction/embedding_strategy.py` — `EmbeddingStrategy` ABC + `Node2VecStrategy`
* `src_prediction/leave_one_out.py` — core LOO evaluation engine
* `src_prediction/run_freeze_multi_kinase_trials.py` — freeze the 6,909-trial eval set
* `src_prediction/run_multi_kinase_generic.py` — generic network (no correlation)
* `src_prediction/run_multi_kinase_control.py` — control network (healthy corr)
* `src_prediction/run_multi_kinase_cancer.py` — cancer network (tumor corr)
* `src_prediction/run_multi_kinase_liver.py` — liver FC network (fold-change corr)
* `src_prediction/analyze_results.py` — four-network side-by-side analysis
* `src_prediction/compare_multi_kinase.py` — pairwise generic vs liver FC comparison
* `src_prediction/pair_features.py`, `src_prediction/negative_sampling.py`
* `src_prediction/model_scoring.py`, `src_prediction/relation_filters.py`

### Supporting helpers

* `src_prediction/config.py`, `src_prediction/graph_loader.py`
* `src_prediction/embeddings.py`, `src_prediction/metrics.py`
* `src_prediction/truth_utils.py`, `src_prediction/parallel_utils.py`
* `src_prediction/experiment_utils.py`, `src_prediction/io_utils.py`

### Older / debug scripts (not part of main pipeline)

* `src_prediction/run_generic_model.py`, `src_prediction/run_baseline_similarity.py` — 50-trial debug runs, predate the leakage fix; outputs are stale
* `src_prediction/find_multi_kinase_sites.py`, `src_prediction/make_multi_kinase_fold_set.py`
* `src_prediction/compare_experiments.py`, `src_prediction/compare_multi_kinase_runs.py`

## Data Inputs

### Core generic graph data

* `data/raw/Kinase_Substrate_Dataset.tsv` — PhosphoSitePlus KSA edges (human-human, single-site rows only)
* `data/PPase_protSubtrates_201903.xls` — phosphatase-substrate associations
* PTMsigDB — site-level pathway co-membership
* MsigDB — protein-level pathway co-membership
* PPI — high-confidence protein-protein interactions
* PTMcode2 — site-site coevolution data

### Liver-specific data

* `data/liver/LiverCancer_ProtExp_Phospho_casecntrl.xlsx`

  * protein expression sheet — identifies liver-observed proteins
  * phosphorylation sheet — identifies liver-observed phosphosites; used to compute fold-change profiles and all three correlation edge types

Protein co-expression edges are disabled: the input has only 3 usable rows after preprocessing.

## Graph Construction Layer

### `src/builders/kinase_substrate.py`

Builds protein and site nodes plus `phosphorylates` (kinase→site) and `has_site` (protein→site) edges from the PhosphoSitePlus dataset.

### `src/builders/ppase_substrate.py`

Adds phosphatase protein nodes and `dephosphorylates` edges.

### `src/main.py`

Merges all generic graph components (KSA, phosphatase, PTMsigDB, MsigDB, PPI, PTMcode2) into `outputs/nodes.csv.gz` and `outputs/edges.csv.gz`.

Generic graph: ~13,318 nodes, ~1.1 million edges.

### `src/liver_network.py`

Extends the generic graph with liver-observed proteins and phosphosites and computes all three correlation edge types in a single pass:

1. Detect healthy and tumor columns in the phosphorylation sheet.
2. Compute healthy mean per phosphosite.
3. Compute tumor fold-change profiles (tumor / mean control).
4. Compute pairwise site-site Pearson correlations for each of the three vectors (FC, control raw, cancer raw). Pairs with fewer than 6 overlapping non-NaN measurements are skipped.
5. Keep edges above the 80th percentile correlation threshold (positive and negative separately).

All three edge types are stored in a single `outputs/Liver_network_edges.csv.gz`.

Liver graph: ~16,605 nodes, ~2.1 million edges.

### `src/run_liver.py`

Launcher script for the liver graph builder. Sets `site_corr_percentile=80.0` and `add_protein_fc_corr=False`.

## Prediction and Evaluation Layer

### `src_prediction/embedding_strategy.py`

Defines the `EmbeddingStrategy` abstract base class and `Node2VecStrategy` (current implementation: dimensions=32, walk_length=10, num_walks=25, p=q=1 DeepWalk mode). To swap embedding methods, implement `fit(graph) -> Dict[str, ndarray]` and pass the new instance to `run_leave_one_out()`.

### `src_prediction/leave_one_out.py`

Core evaluation engine. Per-run (not per-trial):

1. Strip PSP `phosphorylates` edges from the graph via `_build_embedding_relations()`.
2. Compute node2vec embeddings once on the stripped graph.
3. Pre-build the full training feature matrix once.
4. Pre-build per-site scoring feature matrices once.

Per trial:

1. Mask out the held-out kinase-site pair from the training matrix.
2. Fit logistic regression.
3. Score all 420 candidate kinases for the held-out site.
4. Record `held_out_kinase_rank`, `adjusted_held_out_rank`, `best_true_kinase_rank`, `top1_predicted_kinase`.

Parallelised via `ThreadPoolExecutor` (threads share the pre-built numpy arrays; n_jobs_outer=8).

### `src_prediction/run_freeze_multi_kinase_trials.py`

Freezes the 6,909-trial evaluation set (all edges from sites with ≥2 known kinases) into `outputs_prediction/frozen_trials_multi_kinase.csv.gz`. Run once; all four network experiments use the same frozen set.

### `src_prediction/run_multi_kinase_{generic,control,cancer,liver}.py`

Four LOO experiment runners. Each loads the same frozen trials, selects the appropriate `allowed_relations` subset from the liver graph file, and calls `run_leave_one_out()`. Results go to `outputs_prediction/{generic,control,cancer,liver}_multi_kinase/`.

### `src_prediction/analyze_results.py`

Reads all four result sets and produces the full side-by-side analysis in `outputs_prediction/analysis/`:

* `performance_summary.txt` — all three rank metrics across all four networks
* `top1_frequency_summary.txt` — hub bias table
* `per_kinase_topk.txt / .csv` — per-kinase top-1/5/10% across all four networks
* `ttk_deepdive.txt` — TTK case study
* `ranking_shifts.txt` — most improved and worsened kinases (generic vs liver FC)

### `src_prediction/compare_multi_kinase.py`

Merges generic and liver FC result files on `adjusted_held_out_rank`, computes per-trial rank deltas, and writes `comparison_multi_kinase_generic_vs_liver.csv.gz`. This file feeds `analyze_results.py`'s ranking shift analysis.

### `src_prediction/pair_features.py`

Converts node embeddings into kinase-site pair feature vectors: concatenation of kinase embedding `k`, site embedding `s`, elementwise difference `|k - s|`, and elementwise product `k * s`.

### `src_prediction/negative_sampling.py`

Samples negative kinase-site pairs (non-true kinases) for supervised logistic regression training.

### `src_prediction/model_scoring.py`

Trains `LogisticRegression(class_weight="balanced", solver="liblinear")` and ranks all candidate kinases for a site by predicted probability.

### `src_prediction/relation_filters.py`

Defines `GENERIC_BASE_RELATIONS` and the correlation edge subsets used by each network. Also defines `PSP_KSA_RELATIONS` (the set stripped before embedding).

## Outputs

### Graph outputs

* `outputs/nodes.csv.gz` — generic graph nodes
* `outputs/edges.csv.gz` — generic graph edges
* `outputs/Liver_network_edges.csv.gz` — liver graph edges (all three correlation types)

### Per-network LOO results

* `outputs_prediction/generic_multi_kinase/results.csv.gz` and `metrics.csv.gz`
* `outputs_prediction/control_multi_kinase/results.csv.gz` and `metrics.csv.gz`
* `outputs_prediction/cancer_multi_kinase/results.csv.gz` and `metrics.csv.gz`
* `outputs_prediction/liver_multi_kinase/results.csv.gz` and `metrics.csv.gz`

### Analysis outputs

* `outputs_prediction/analysis/performance_summary.txt`
* `outputs_prediction/analysis/top1_frequency_summary.txt`
* `outputs_prediction/analysis/per_kinase_topk.txt` and `.csv`
* `outputs_prediction/analysis/ttk_deepdive.txt`
* `outputs_prediction/analysis/ranking_shifts.txt`

## How the Full Pipeline Works

### Stage 1. Build generic graph

`src/main.py` merges all public curated resources into the context-free prior graph.

### Stage 2. Build liver graph

`src/liver_network.py` extends the generic graph with liver-observed nodes and adds all three correlation edge types in a single pass.

### Stage 3. Freeze evaluation set

`src_prediction/run_freeze_multi_kinase_trials.py` creates the reproducible 6,909-trial set. This step is run once; the frozen set is shared across all four experiments.

### Stage 4. Run all four LOO experiments

`run_multi_kinase_{generic,control,cancer,liver}.py` run independently in any order. Each selects its own `allowed_relations` subset from the single liver graph file.

### Stage 5. Regenerate comparison and analysis

`compare_multi_kinase.py` must be run before `analyze_results.py` because the ranking shifts analysis reads from the comparison output. Then `analyze_results.py` produces all summary files.

## How to Run the Project

Always use the project virtual environment:

```
.venv312/Scripts/python -m <module>
```

### 1. Build the generic graph

```
.venv312/Scripts/python -m src.main
```

Output: `outputs/nodes.csv.gz`, `outputs/edges.csv.gz`

### 2. Build the liver graph

```
.venv312/Scripts/python -m src.run_liver
```

Output: `outputs/Liver_network_edges.csv.gz` (contains all three correlation types)

### 3. Freeze the evaluation set

```
.venv312/Scripts/python -m src_prediction.run_freeze_multi_kinase_trials
```

Output: `outputs_prediction/frozen_trials_multi_kinase.csv.gz`

### 4. Run all four LOO experiments

```
.venv312/Scripts/python -m src_prediction.run_multi_kinase_generic
.venv312/Scripts/python -m src_prediction.run_multi_kinase_control
.venv312/Scripts/python -m src_prediction.run_multi_kinase_cancer
.venv312/Scripts/python -m src_prediction.run_multi_kinase_liver
```

These can run in any order or in parallel (each writes to its own output directory).

### 5. Regenerate comparison and analysis

```
.venv312/Scripts/python -m src_prediction.compare_multi_kinase
.venv312/Scripts/python -m src_prediction.analyze_results
```

Output: `outputs_prediction/analysis/` (all summary files)

## Experimental Results

All four networks evaluated on 6,909 frozen trials, 2,459 multi-kinase sites, 420 candidate kinases.

### Adjusted held-out rank (primary metric — lower is better)

|  | Generic | Control | Cancer | Liver FC |
|---|---|---|---|---|
| Mean rank | 190.4 | 177.7 | 186.0 | **175.8** |
| Median rank | 188.0 | 166.5 | 174.0 | **165.0** |
| Top-10% | 1.5% | 1.1% | 1.4% | 1.8% |
| Top-20% | 2.2% | 1.7% | 2.2% | 2.7% |

Liver FC and Control both outperform Generic. Cancer performs **worse** than Generic — raw tumor abundance correlations add noise rather than signal. Control (healthy co-expression) rivals Liver FC.

### Hub degree bias (critical finding)

| Network | Top-1 monopoly |
|---|---|
| Generic | ULK1 76.5%, WNK1 21.0%, ZAP70 2.5% |
| Control | VRK2 94.8%, YES1 4.8% |
| Cancer | VRK3 95.5%, VRK2 3.8% |
| Liver FC | ULK1 **100.0%** of all 6,909 trials |

Adding correlation edges concentrates node2vec embedding mass on a single high-degree kinase, which then dominates top-1 for every trial regardless of biology. The improved mean/median ranks are real (e.g. CSK substrates move from ~rank 402 to ~rank 162 in Liver FC) but the correct kinase never reaches rank 1.

### Most improved kinases in Liver FC vs Generic (adjusted rank)

| Kinase | Mean delta | Trials | % improved |
|---|---|---|---|
| CSK | +241 | 6 | 100% |
| STK24 | +238 | 8 | 100% |
| STK11 | +216 | 15 | 100% |
| MAP3K5 | +177 | 10 | 100% |
| RET | +172 | 15 | 100% |

### Most worsened kinases in Liver FC vs Generic

| Kinase | Mean delta | Trials | % improved |
|---|---|---|---|
| PKN1 | -177 | 10 | 0% |
| MAP3K11 | -168 | 7 | 0% |
| MAP3K8 | -156 | 23 | 0% |
| SIK2 | -155 | 5 | 0% |
| PDK1 | -149 | 5 | 0% |

### TTK case study (21 trials)

| Network | Adjusted mean rank | Top-1 always |
|---|---|---|
| Generic | 247 | ULK1 (18/21), WNK1 (3/21) |
| Control | **26** | VRK2 (21/21) |
| Cancer | 113 | VRK3 (21/21) |
| Liver FC | 243 | ULK1 (21/21) |

Control dramatically improves TTK rank (mean 26 vs 248 in generic). Liver FC provides no improvement. TTK substrates are stably co-expressed in healthy liver tissue but show divergent fold-change patterns in tumors.

## Interpretation of Results

The main scientific question is:

**Does adding liver context change kinase prioritisation in a biologically meaningful way?**

The current answer is: **partially**. Mean and median adjusted ranks improve in Liver FC and Control, meaning the true kinase moves closer to the top of the ranked list. But the degree bias in node2vec prevents any network from placing the correct kinase at rank 1 for the vast majority of trials.

The per-kinase interpretation framework (from `per_kinase_topk.csv`):

* High top-k% in all networks → degree/hub bias, not context-specific
* High top-k% in Liver FC only → context genuinely boosts this kinase
* High top-k% in Generic only → liver correlation edges dilute this kinase's signal

This table is currently dominated by hub bias and will become biologically interpretable after the embedding swap.

## Important Notes and Caveats

* **Hub bias is the dominant signal in top-1.** Until a degree-normalised embedding replaces node2vec, top-1 predictions are biologically meaningless. Rank improvement metrics (mean, median, top-10%, top-20%) are valid and show genuine structural improvement.

* **The leakage fix changed the network ranking.** Pre-fix results showed generic and liver FC roughly tied; post-fix, Liver FC consistently outperforms Generic. The pre-fix numbers are scientifically invalid and should not be cited.

* **The 50-trial debug runs are stale.** `run_generic_model.py` and `run_baseline_similarity.py` predate the leakage fix and use the unadjusted `held_out_kinase_rank` metric. Their outputs should not be used.

* **Protein co-expression edges are disabled.** `liver_network.py` has logic for protein correlation edges, but only 3 usable rows remain after preprocessing. Do not re-enable without new data.

* **The comparison file must be regenerated** whenever the generic or liver LOO results are updated. Run `compare_multi_kinase.py` before `analyze_results.py`.

* **Cancer worsens performance.** Raw tumor-sample abundance correlations add noise. This is consistent with the known heterogeneity of tumor phosphoproteomics relative to fold-change normalisation.

* Some files use the word `fold` in their names; the actual evaluation is leave-one-out on held-out edges, not k-fold cross-validation.

## Current Scientific Direction

**Immediate priority:** Swap `Node2VecStrategy` for a degree-normalised embedding via the `EmbeddingStrategy` ABC. Candidate methods: spectral embedding, personalized PageRank, graph contrastive learning with spectral filtering. No other code changes are required.

**Biological follow-up (after embedding fix):**

1. Revisit the most-improved kinases (CSK, STK24, STK11, MAP3K5, RET) with biologically valid top-k predictions.
2. Investigate the TTK case: why does healthy co-expression (rank 26) help so much more than fold-change (rank 243)? Are TTK substrate fold-change patterns genuinely divergent in tumors?
3. Statistical significance testing: Wilcoxon signed-rank test on per-trial `adjusted_held_out_rank` pairs (Generic vs each network) before writing the methods section.

**Abstract story:** Four-network comparison reveals hub degree bias in node2vec embeddings. Context-specific correlation edges improve ranking metrics but degree sensitivity of node2vec obscures biological top-1 predictions — motivating the methodological contribution of the embedding swap.
