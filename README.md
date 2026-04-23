# Kinase-Substrate Graph Project

## Table of Contents

* [Overview](#overview)
* [Scientific Motivation](#scientific-motivation)
* [Core Design Choices](#core-design-choices)
* [Repository Layout](#repository-layout)
* [Data Inputs](#data-inputs)
* [Graph Construction Layer](#graph-construction-layer)

  * [`src/builders/kinase_substrate.py`](#srcbuilderskinase_substratepy)
  * [`src/builders/ppase_substrate.py`](#srcbuildersppase_substratepy)
  * [`src/main.py`](#srcmainpy)
  * [`src/liver_network.py`](#srcliver_networkpy)
  * [`src/run_liver.py`](#srcrun_liverpy)
* [Prediction and Evaluation Layer](#prediction-and-evaluation-layer)

  * [`src_prediction/run_freeze_folds.py`](#src_predictionrun_freeze_foldspy)
  * [`src_prediction/pair_features.py`](#src_predictionpair_featurespy)
  * [`src_prediction/negative_sampling.py`](#src_predictionnegative_samplingpy)
  * [`src_prediction/model_scoring.py`](#src_predictionmodel_scoringpy)
  * [`src_prediction/leave_one_out.py`](#src_predictionleave_one_outpy)
  * [`src_prediction/run_baseline_similarity.py`](#src_predictionrun_baseline_similaritypy)
  * [`src_prediction/compare_experiments.py`](#src_predictioncompare_experimentspy)
  * [`src_prediction/relation_filters.py`](#src_predictionrelation_filterspy)
* [Supporting Helper Files](#supporting-helper-files)
* [Recent Exploratory Analysis Files](#recent-exploratory-analysis-files)
* [Outputs](#outputs)
* [How the Full Pipeline Works](#how-the-full-pipeline-works)
* [How to Run the Project](#how-to-run-the-project)

  * [1. Build the generic graph](#1-build-the-generic-graph)
  * [2. Build the liver-specific graph](#2-build-the-liver-specific-graph)
  * [3. Freeze leave-one-out evaluation cases](#3-freeze-leave-one-out-evaluation-cases)
  * [4. Run the generic prediction experiment](#4-run-the-generic-prediction-experiment)
  * [5. Run the liver prediction experiment](#5-run-the-liver-prediction-experiment)
  * [6. Compare generic vs liver predictions](#6-compare-generic-vs-liver-predictions)
  * [7. Optional exploratory analysis for multi-kinase sites](#7-optional-exploratory-analysis-for-multi-kinase-sites)
* [Interpretation of Results](#interpretation-of-results)
* [Important Notes and Caveats](#important-notes-and-caveats)
* [Current Scientific Direction](#current-scientific-direction)

## Overview

This project builds and evaluates a context-aware biological network for kinase-substrate association, or KSA, prediction, with a current application to liver cancer.

The full workflow has two major layers:

1. **Graph construction**

   * Build a generic heterogeneous phosphorylation graph from curated public resources.
   * Extend that graph into a liver-specific graph using liver proteomic and phosphoproteomic data.

2. **Prediction and evaluation**

   * Use the graph to rank candidate kinases for phosphosites.
   * Compare generic versus liver-specific rankings to identify context-dependent changes.

The full pipeline is:

1. Build the generic graph
2. Build the liver-specific graph
3. Prepare evaluation inputs
4. Freeze held-out leave-one-out cases for reproducibility
5. Run prediction on the generic graph
6. Run prediction on the liver graph
7. Compare ranking differences
8. Inspect biologically interesting changes, especially at sites with multiple reported kinases

## Scientific Motivation

Many existing computational KSA methods work out of context. They often rely on static information such as sequence motifs, structural features, or generic prior knowledge. However, phosphorylation is highly transient and context-dependent. The kinase most relevant to a phosphosite may differ across tissues, disease states, and biological conditions.

This project addresses that limitation by augmenting a generic phosphorylation network with liver-specific phosphoproteomic context and then comparing how kinase rankings change.

A related background issue is the **streetlight effect**. Well-studied kinases tend to have much more available data, which makes them easier for both biological research and prediction methods to focus on. That remains an important background motivation, but the main current scientific emphasis of the repository is **context-aware prediction**, not bias alone.

The broader long-term goal is to build a framework that can support other disease-specific or tissue-specific contexts beyond liver cancer.

## Core Design Choices

### 1. Only two node types

The graph uses only two main node types:

* `protein`
* `site`

Examples:

* `PROTEIN:AKT1`
* `PROTEIN:PPP1CA`
* `PROTEIN:FOXO3`
* `SITE:FOXO3-T32`

Kinases, phosphatases, and other proteins are all represented as `protein` nodes. Their biological roles are stored as metadata instead of being split into separate node types.

Important protein metadata fields include:

* `is_kinase`
* `is_phosphatase`
* `protein_role`

This keeps the graph structure simple while preserving biological interpretability.

### 2. Site-specific modeling

Phosphorylation is modeled at the phosphosite level, not only at the protein level. This means different sites on the same protein are distinct graph nodes.

For example:

* `SITE:FOXO3-T32`
* `SITE:FOXO3-S253`

This is essential because the prediction target is kinase-to-site association, not just kinase-to-protein association.

### 3. Generic graph first, context later

The project intentionally separates:

* a **generic graph** built from curated public resources
* a **liver-specific graph** built by augmenting the generic graph with disease context

This separation allows direct comparison between:

* what the context-free graph prioritizes
* what the liver-specific graph prioritizes

### 4. Context is encoded as edges

Liver context is not only stored as annotations. It is introduced through additional graph edges, especially phosphosite correlation edges derived from tumor fold-change profiles relative to healthy controls.

This means disease context changes the graph topology itself.

### 5. Leave-one-out evaluation, not standard k-fold cross-validation

Some files historically use the word `fold`, but the actual evaluation framework is **leave-one-out** on held-out edges.

A held-out case means:

1. remove one known kinase-site edge
2. train on the remaining data
3. rank candidate kinases for that one held-out site

So the "frozen folds" files are really frozen held-out leave-one-out evaluation cases.

## Repository Layout

### Graph construction files

* `src/builders/kinase_substrate.py`
* `src/builders/ppase_substrate.py`
* `src/main.py`
* `src/liver_network.py`
* `src/run_liver.py`

### Prediction files

* `src_prediction/run_freeze_folds.py`
* `src_prediction/pair_features.py`
* `src_prediction/negative_sampling.py`
* `src_prediction/model_scoring.py`
* `src_prediction/leave_one_out.py`
* `src_prediction/run_baseline_similarity.py`
* `src_prediction/compare_experiments.py`
* `src_prediction/relation_filters.py`

### Supporting helper files

* `src/config.py`
* `src/io_utils.py`
* `src/ids.py`
* `src/normalize.py`
* `src_prediction/config.py`
* `src_prediction/graph_loader.py`
* `src_prediction/embeddings.py`
* `src_prediction/metrics.py`
* `src_prediction/truth_utils.py`
* `src_prediction/parallel_utils.py`
* `src_prediction/experiment_utils.py`

### Recent exploratory analysis files

* `src_prediction/find_multi_kinase_sites.py`
* `src_prediction/make_multi_kinase_fold_set.py`
* `src_prediction/make_multi_kinase_fold_set_10.py`
* `src_prediction/run_generic_multi_kinase.py`
* `src_prediction/run_liver_multi_kinase.py`
* `src_prediction/compare_multi_kinase_runs.py`
* `src_prediction/check_ttk_from_existing_debug_results.py`
* `src_prediction/count_ttk_substrates_in_graph.py`

These newer files support supervisor-driven follow-up analyses and are not all part of the original core pipeline.

## Data Inputs

### Core generic graph data

The generic graph integrates multiple curated biological resources, including:

* kinase-substrate associations
* phosphatase-substrate associations
* PTMsigDB site-level pathway information
* MsigDB protein-level pathway information
* high-confidence protein-protein interaction data
* PTMcode2 site-site coevolution data

### Liver-specific data

The liver-specific network uses a liver proteomic and phosphoproteomic Excel dataset with at least:

* a protein expression sheet
* a phosphorylation sheet

These sheets are used to:

* identify proteins observed in liver data
* identify phosphosites observed in liver data
* compute tumor-relative-to-control fold-change profiles
* derive context-specific phosphosite correlation edges

## Graph Construction Layer

### `src/builders/kinase_substrate.py`

**Purpose**
Builds the core graph components from the kinase-substrate dataset.

**Creates**

* protein nodes for kinases
* protein nodes for substrate proteins
* site nodes for phosphosites
* `has_site` edges from substrate protein to phosphosite
* `phosphorylates` edges from kinase protein to phosphosite

**Why it exists**
This file defines the main biological supervision signal in the project:

`kinase -> phosphosite`

Without it, there is no known kinase-site relationship to learn from.

**Main function**

* `build_kinase_substrate_graph(path: str)`

**Input**

* `data/raw/Kinase_Substrate_Dataset.tsv`

**Output**
Two pandas DataFrames:

* `nodes_df`
* `edges_df`

**Important internal logic**

1. Read the table
2. Extract kinase, substrate, organism, and site information
3. Keep only human-human rows
4. Keep only rows with exactly one phosphosite token
5. Normalize protein names and site labels
6. Construct canonical node IDs
7. Create nodes and edges
8. Drop duplicates

### `src/builders/ppase_substrate.py`

**Purpose**
Adds phosphatase-substrate information to the graph.

**Creates**

* phosphatase protein nodes
* substrate protein nodes
* phosphosite nodes
* `has_site` edges
* `dephosphorylates` edges

**Why it exists**
Phosphorylation regulation includes both addition and removal of phosphate groups. Adding phosphatase information makes the graph more mechanistic.

**Main function**

* `build_ppase_substrate_graph(path: str)`

**Input**

* `data/PPase_protSubtrates_201903.xls`

### `src/main.py`

**Purpose**
Main entry point for building the generic graph.

**Combines**

1. kinase-substrate graph
2. phosphatase-substrate graph
3. PTMsigDB site-site pathway edges
4. MsigDB protein-protein pathway edges
5. PPI edges
6. PTMcode2 site-site coevolution edges

**Writes**

* `outputs/nodes.csv.gz`
* `outputs/edges.csv.gz`

**Why it exists**
Each builder produces one graph component. This file merges everything into one final generic network and adds protein role metadata.

### `src/liver_network.py`

**Purpose**
Builds the liver-specific graph by extending the generic graph with liver cancer data.

**Main idea**
Start with the generic graph, then:

1. add liver-observed proteins
2. add liver-observed phosphosites
3. reconnect those nodes to known resources
4. add liver-specific phosphosite correlation edges

**Most important context step**
Phosphosite fold-change correlation edges are computed from tumor fold-change profiles relative to healthy controls.

**Step-by-step logic**

1. Detect healthy and tumor columns
2. Compute healthy mean per phosphosite
3. Compute tumor fold-change profiles
4. Compute pairwise site-site correlations
5. Keep only strong positive and strong negative correlations above percentile thresholds
6. Add context-specific edges such as:

   * `site_corr_fc_pos`
   * `site_corr_fc_neg`

**Why it matters scientifically**
This is the step that makes the graph context-aware and liver-cancer-specific.

### `src/run_liver.py`

**Purpose**
Small launcher script for liver graph construction.

**Why it exists**
It makes experiment-specific liver-network settings explicit and reproducible.

**Typical settings**

* `site_corr_percentile=80.0`
* `add_protein_fc_corr=False`

## Prediction and Evaluation Layer

### `src_prediction/run_freeze_folds.py`

**Purpose**
Creates a fixed set of held-out leave-one-out evaluation cases.

**Why it exists**
Leave-one-out is expensive. Freezing the held-out cases ensures that experiments are compared on the same evaluation subset.

**Typical output**

* `outputs_prediction/frozen_folds_10.csv.gz`

### `src_prediction/pair_features.py`

**Purpose**
Converts node embeddings into machine-learning features for kinase-site pairs.

**Feature design**
Given:

* kinase embedding `k`
* site embedding `s`

The pair feature vector includes:

1. `k`
2. `s`
3. `|k - s|`
4. `k * s`

### `src_prediction/negative_sampling.py`

**Purpose**
Creates negative kinase-site examples for supervised training.

**Why it exists**
The classifier needs both positive and negative examples.

### `src_prediction/model_scoring.py`

**Purpose**
Trains the supervised scorer and ranks candidate kinases for a site.

**Current model**

* `LogisticRegression`

**Why it exists**
Earlier ranking was based on cosine similarity between embeddings. The current version uses a learned scoring function because raw similarity did not interpret the added liver context reliably.

### `src_prediction/leave_one_out.py`

**Purpose**
Core leave-one-out evaluation engine.

**Main steps per held-out case**

1. Remove one known phosphorylation edge
2. Rebuild the graph view
3. Recompute node2vec embeddings
4. Train the logistic regression model
5. Score all candidate kinases for the held-out site
6. Record results

**Important outputs recorded**

* rank of the specific held-out kinase
* best rank among all known true kinases for the site
* top-1 predicted kinase
* top-1 score

**Why the best-true metric matters**
A site may have more than one known kinase. This metric makes evaluation more biologically faithful.

### `src_prediction/run_baseline_similarity.py`

**Purpose**
Main experiment runner, despite the older name.

**What it controls**

* graph choice: generic or liver
* whether site-context edges are included
* whether protein-context edges are included
* held-out evaluation file
* node2vec parameters
* parallelization settings

**What it does**

1. Load graph
2. Load candidate kinases
3. Load held-out evaluation cases
4. Build the allowed relation set
5. Run leave-one-out evaluation
6. Summarize metrics
7. Write outputs

### `src_prediction/compare_experiments.py`

**Purpose**
Compares generic and liver experiment outputs case by case.

**What it computes**

* whether top-1 changed
* held-out rank delta
* best-true rank delta
* whether liver improved held-out rank
* whether liver improved best-true rank

**Why it matters**
This is the main bridge from prediction results to biological interpretation.

### `src_prediction/relation_filters.py`

**Purpose**
Defines which edge types are allowed in each prediction experiment.

**Current relation groups**

Generic base relations:

* `has_site`
* `phosphorylates`
* `dephosphorylates`
* `site_same_pathway`
* `protein_same_pathway`
* `ppi_high_confidence`
* `site_coevolution`

Site context relations:

* `site_corr_fc_pos`
* `site_corr_fc_neg`

Protein context relations:

* `protein_corr_fc_pos`
* `protein_corr_fc_neg`

## Supporting Helper Files

### `src/ids.py`

Creates standardized node IDs.

### `src/normalize.py`

Normalizes protein names and phosphosite labels.

### `src/io_utils.py`

Provides centralized file reading and writing utilities.

### `src_prediction/graph_loader.py`

Loads node and edge tables into the structures needed by the prediction layer.

### `src_prediction/embeddings.py`

Builds graph objects and computes node2vec embeddings.

### `src_prediction/truth_utils.py`

Builds mappings such as:

* site -> all known true kinases

### `src_prediction/parallel_utils.py`

Handles outer parallel execution across held-out cases.

### `src_prediction/metrics.py`

Summarizes case-level results into experiment-level metrics.

### `src_prediction/experiment_utils.py`

Writes results, metrics, and summary text in a standard format.

## Recent Exploratory Analysis Files

These files support newer biological follow-up and debugging.

### `src_prediction/find_multi_kinase_sites.py`

Finds phosphosites with more than one reported kinase.

### `src_prediction/make_multi_kinase_fold_set.py`

Builds a leave-one-out evaluation set restricted to multi-kinase sites.

### `src_prediction/make_multi_kinase_fold_set_10.py`

Builds a small 10-case debug subset from the multi-kinase evaluation set.

### `src_prediction/run_generic_multi_kinase.py`

Runs the generic model on the multi-kinase subset.

### `src_prediction/run_liver_multi_kinase.py`

Runs the liver model on the multi-kinase subset.

### `src_prediction/compare_multi_kinase_runs.py`

Compares generic and liver results specifically for the multi-kinase-site analysis.

### `src_prediction/check_ttk_from_existing_debug_results.py`

Examines TTK behavior in existing generic-vs-liver comparison outputs.

### `src_prediction/count_ttk_substrates_in_graph.py`

Counts known TTK substrate-site edges in the graph.

## Outputs

### Generic graph outputs

* `outputs/nodes.csv.gz`
* `outputs/edges.csv.gz`

### Liver graph outputs

Typical liver-network runs write liver-specific node and edge tables under the configured output paths.

### Prediction outputs

Each experiment typically writes:

* results file
* metrics file
* summary text file

### Comparison outputs

Comparison scripts typically write files that summarize:

* changed top-1 predictions
* rank deltas
* improved liver cases

### Exploratory outputs

Additional recent outputs include:

* multi-kinase site lists
* restricted evaluation subsets
* TTK-specific summaries
* debug tables for changed liver predictions

## How the Full Pipeline Works

### Stage 1. Build generic graph

Files:

* `src/builders/kinase_substrate.py`
* `src/builders/ppase_substrate.py`
* `src/main.py`

Creates the context-free biological prior graph.

### Stage 2. Build liver-specific graph

Files:

* `src/liver_network.py`
* `src/run_liver.py`

Adds liver-specific proteins, phosphosites, and phosphosite-correlation context.

### Stage 3. Prepare evaluation inputs

File:

* `src_prediction/run_freeze_folds.py`

Creates a reproducible held-out evaluation subset.

### Stage 4. Run trained model on generic graph

Uses the prediction stack on the generic graph.

### Stage 5. Run trained model on liver graph

Uses the same stack on the liver-specific graph.

### Stage 6. Compare generic versus liver

Uses comparison scripts to identify changed kinase rankings.

### Stage 7. Biological follow-up

Current biological follow-up focuses on:

* sites with multiple reported kinases
* changed top-1 predictions between generic and liver
* most frequently predicted liver-specific kinases if multi-kinase cases are not yet strong enough

## How to Run the Project

### 1. Build the generic graph

```bash
python -m src.main
```

Expected main outputs:

```text
outputs/nodes.csv.gz
outputs/edges.csv.gz
```

### 2. Build the liver-specific graph

```bash
python -m src.run_liver
```

Expected outputs:

```text
liver-specific node and edge tables under the configured output paths
```

### 3. Freeze leave-one-out evaluation cases

```bash
python -m src_prediction.run_freeze_folds
```

Expected output:

```text
outputs_prediction/frozen_folds_10.csv.gz
```

### 4. Run the generic prediction experiment

Edit `src_prediction/run_baseline_similarity.py` so that:

```python
graph_choice = "generic"
include_site_corr = False
include_protein_corr = False
```

Then run:

```bash
python -m src_prediction.run_baseline_similarity
```

### 5. Run the liver prediction experiment

Edit `src_prediction/run_baseline_similarity.py` so that:

```python
graph_choice = "liver"
include_site_corr = True
include_protein_corr = False
```

Then run:

```bash
python -m src_prediction.run_baseline_similarity
```

### 6. Compare generic vs liver predictions

```bash
python -m src_prediction.compare_experiments
```

This step identifies:

* top-1 changes
* held-out rank changes
* best-true rank changes
* cases where liver improved ranking

### 7. Optional exploratory analysis for multi-kinase sites

Find sites with multiple known kinases:

```bash
python -m src_prediction.find_multi_kinase_sites
```

Build the multi-kinase evaluation set:

```bash
python -m src_prediction.make_multi_kinase_fold_set
```

Build a small 10-case debug subset:

```bash
python -m src_prediction.make_multi_kinase_fold_set_10
```

Run generic multi-kinase experiment:

```bash
python -m src_prediction.run_generic_multi_kinase
```

Run liver multi-kinase experiment:

```bash
python -m src_prediction.run_liver_multi_kinase
```

Compare multi-kinase runs:

```bash
python -m src_prediction.compare_multi_kinase_runs
```

## Interpretation of Results

The main scientific question is not only whether a global metric improves. The main question is:

**Does adding liver context change kinase prioritization in a biologically meaningful way?**

A typical interpretation pattern is:

* the generic graph ranks one kinase highly for a site
* the liver graph ranks a different kinase highly for the same site
* that changed ranking becomes a biological follow-up question

The strongest future cases are expected to be phosphosites with multiple known kinases, where the generic and liver graphs may prioritize different biologically plausible kinases.

## Important Notes and Caveats

* `run_baseline_similarity.py` is still named after an earlier cosine-similarity stage, but the current pipeline uses a trained logistic scorer.
* Some newer multi-kinase analysis scripts are exploratory and may still require technical refinement.
* A recent broad multi-kinase run exposed a logistic solver limitation for some larger targeted subsets.
* The term `fold` appears in the code, but the evaluation method is leave-one-out on held-out edges, not standard k-fold cross-validation.
* Protein fold-change correlation was explored conceptually, but the current working liver configuration uses phosphosite correlation as the main active context layer.

## Current Scientific Direction

The main contribution of the repository in its current form is:

* **context-aware KSA prediction**
* **generic vs disease-specific kinase ranking comparison**
* **biological interpretation of ranking shifts**

The current biological follow-up strategy is:

1. inspect sites with multiple reported kinases
2. check whether the liver-specific graph prioritizes a different and potentially more liver-relevant kinase
3. use those cases for downstream biological investigation
4. in future work, experimentally validate the strongest candidates

If multi-kinase sites do not yield strong enough examples, another strategy is to analyze the most frequently predicted liver-specific kinases and examine their relationship to liver cancer biology.
