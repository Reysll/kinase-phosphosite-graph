---
name: project-overview
description: High-level architecture of the Kinase-Substrate HGN project — two pipelines, their outputs, and the biological goal
metadata:
  type: project
---

This is a master's-level research project building a heterogeneous graph network (HGN) for kinase-substrate phosphorylation prediction in a liver cancer context.

**Two pipelines:**

1. `src/` — Generic graph construction
   - Builds HGN from: Kinase_Substrate_Dataset, PPase substrates, PTMsigDB, MsigDB, PPI (high_confidence_score), PTMcode2
   - Entry: `python -m src.main`
   - Outputs: `outputs/nodes.csv.gz`, `outputs/edges.csv.gz`
   - Node types: protein (kinase/phosphatase/both/protein), site
   - Edge types: has_site, phosphorylates, dephosphorylates, site_same_pathway, protein_same_pathway, ppi_high_confidence, site_coevolution

2. `src/liver_network.py` — Liver extension
   - Expands generic graph with liver proteomics (LiverCancer_ProtExp_Phospho_casecntrl.xlsx)
   - Adds liver-specific nodes, reconnects via all databases, adds fold-change correlation edges
   - Entry: `python -m src.run_liver`
   - Outputs: `outputs/Liver_network_nodes.csv.gz`, `outputs/Liver_network_edges.csv.gz`
   - Fold-change correlation: healthy sample mean → tumor fold change → pairwise Pearson → 80th percentile threshold
   - Protein correlation disabled (only 3 usable protein rows in data — awaiting professor clarification)

3. `src_prediction/` — LOO-CV prediction pipeline
   - Evaluation set prep: `run_prepare_eval_set.py` → candidate kinases + liver positive edges
   - Freeze trials: `run_freeze_trials.py` → frozen_trials_10.csv.gz (debug: 10 trials)
   - Run LOO: `run_baseline_similarity.py` → node2vec embeddings + logistic regression scorer
   - Compare: `compare_experiments.py` → generic vs liver comparison

**Biological goal:** Identify liver-cancer-relevant kinase activity from phosphoproteomics. Validate that liver-specific graph context improves kinase ranking over generic graph. Ultimately discuss frequently predicted kinases and their liver cancer relevance; future work = experimental validation.

**Key scale numbers:**
- Generic graph: ~13,318 nodes, ~1,095,542 edges
- Liver graph: ~16,605 nodes, ~2,110,777 edges
- Candidate kinases: ~500 (from generic graph)
- Positive LOO edges (liver sites): see prep_summary.txt
- TTK has 68 substrates in the graph

**Current active pipeline (as of 2026-06-25) — supersedes item 3's debug-era description above for current work:**

4. `src_prediction/` — Phase 1: multi-kinase node2vec LOO
   - Embedding leakage fix: PSP `phosphorylates` edges stripped from the graph before fitting node2vec; kept as training labels
   - Four networks (generic / control / cancer / liver FC) selected via `allowed_relations` on a single liver graph file — no rebuild per experiment
   - Entry: `run_multi_kinase_{generic,control,cancer,liver}.py`
   - Trial set: 6,909 multi-kinase trials (sites with 2+ known PSP kinases), frozen via `run_freeze_multi_kinase_trials.py`
   - Primary metric: `adjusted_held_out_rank` (other known true co-kinases dropped from the scored list before ranking)

5. Phase 2: spectral embedding + category ranking
   - `SpectralEmbeddingStrategy` (normalized graph Laplacian, degree-normalizing) added to address node2vec's hub bias (VRK family dominated top-1 in Phase 1)
   - Category-based ranking: kinases split into poor / average / rich by PSP substrate count, ranked only within their own category
   - Entry: `run_all_spectral.py` (runs all 4 networks)
   - Finding: spectral halves median rank vs node2vec, but hub bias persists at the category level (one kinase still wins ~100% within each category)

6. Phase 3 (2026-06-25): predict-all inference pass
   - Generalizes LOO from "evaluate only the held-out edge on multi-kinase sites" to "score every (candidate kinase, site) pair" — needed because LOO can't produce a full kinase x site rank matrix
   - Leakage avoided via two-tier scoring: known true edges get their own edge masked before scoring; everything else uses one global model (nothing to leak for a pair that was never a positive label)
   - `ConcatEmbeddingStrategy` added (concatenates node2vec + spectral vectors) as a third embedding option alongside the two from Phase 1/2
   - Scope: 3,228 FC-valid liver sites x ~420 candidate kinases x 4 networks x 3 embeddings = 12 combos
   - Entry: `run_inference_all.py`
   - Caveat: "generic" network only covers 391/3,228 (12%) of FC-valid sites (its graph is essentially PSP-annotated sites, a much smaller population than the broader liver-proteomics-detected set the other three networks cover at 100%)

7. Result 1 + Result 2 (consume Phase 3 output)
   - Result 1: average rank per kinase per network, clustergram (`build_rank_heatmap.py`) — one per embedding strategy
   - Result 2: per-kinase rank-distribution KS-test between two networks (`build_rank_ks_test.py`), run for both control-vs-liver_fc and cancer-vs-control
   - Caveat found while building Result 2: most kinases' rank is nearly constant across all sites within one network (hub bias), so a naive KS-test over-reports "significant" kinases — fixed by flagging kinases whose rank std-dev shows they actually vary by site (`informative` column)

8. KSEA (kinase activity enrichment, on the liver FC proteomics dataset directly, independent of LOO/inference-all)
   - Formula (Wiredja et al. 2017): `score = (s_bar - p_bar) * sqrt(m) / delta` — sqrt(m) in the numerator rewards (not penalizes) larger substrate sets
   - Entry: `run_ksea.py`
   - 17/64 scoreable kinases reach significance after the 2026-06-23 formula fix (an earlier version had `m` in the denominator instead, which penalized large substrate sets)

9. Result 1 fix — informative filter + delta-vs-Control (2026-06-25)
   - Found the same hub-bias degeneracy documented in item 7's Result 2 caveat also applied to Result 1: under spectral, a kinase's mean rank was nearly identical across Control/Cancer/Liver FC (e.g. TSSK4: 414.0/414.0/414.0), so the clustergram's "Cancer+Liver FC cluster tightest, Generic the outlier" finding mostly just restated that Generic differs structurally, not a real liver-cancer signal
   - Fix in `build_rank_heatmap.py`: (a) clustermaps restricted to `informative` kinases (same std_rank > 1.0 definition as Result 2, computed over control/cancer/liver_fc only); (b) added delta tables (`rank_heatmap_delta_{embedding}.{csv,png}`) — `mean_rank[disease] - mean_rank[Control]` for Cancer/Liver FC, per Dr. Ayati's Option 2 framing (score = score_disease − score_control) — which cancels the hub-dominated baseline and isolates the real network-dependent shift
   - Finding: under spectral (271/415 informative kinases — the least hub-biased embedding), IRAK4/RIPK1/RIPK3/TRPM7/MLKL all gain predicted relevance in both Cancer-vs-Control and Liver-FC-vs-Control — a biologically coherent necroptosis/innate-immune cluster, flagged for Dr. Ayati as a candidate finding (not yet independently validated)

**Directory structure (reorganized 2026-06-25):** both `src/` and `src_prediction/` are split into functional subpackages — `core/` (config, io_utils, and other shared infra), `embeddings/` (`src_prediction` only — embeddings.py, embedding_strategy.py), `data_prep/` (`src_prediction` only — eval set / trial / negative-sampling helpers), `engine/` (`src_prediction` only — leave_one_out.py, inference_all.py), `analysis/` (`src_prediction` only — analyze_results.py, ksea.py, build_rank_heatmap.py, build_rank_ks_test.py, etc.), `runners/` (every `run_*.py` entry point in both trees, plus `main.py`/`run_liver.py` in `src/`). `src/builders/` already existed as its own subpackage and was left untouched. All `python -m` commands now use the longer dotted path, e.g. `python -m src_prediction.runners.run_multi_kinase_generic`.

**Legacy/dead code (intentionally kept, not removed or moved):** nine files sit flat at the root of `src_prediction/`, confirmed unused during the 2026-06-25 reorg but preserved in case they're useful as reference later — `compare_experiments.py`, `compare_multi_kinase_runs.py`, `run_generic_model.py`, `run_baseline_similarity.py` (the original LOO debug scripts referenced in item 3 above — predate the embedding-leakage fix, results invalid), `run_generic_multi_kinase.py` (pre-leakage-fix duplicate of the active `run_multi_kinase_generic.py`), `run_liver_multi_kinase.py` (empty stub), `make_multi_kinase_fold_set.py`, `make_multi_kinase_fold_set_10.py`, `migrate_trial_index.py`. Their internal imports still use the old pre-reorg flat paths and will raise `ImportError` if run — that's expected, since they were already non-functional before the reorg too.

**Updated key scale numbers (2026-06-25):** candidate kinases is actually **420** (the item-3-era "~500" above was an early rough estimate, kept as-is per the no-removal note); 2,459 multi-kinase sites (2+ known PSP kinases) → 6,909 multi-kinase LOO trials; 3,228 phosphosites have valid liver FC data (the scope for Phase 3 / KSEA, a different and stricter filter than the 10,146/13,410 liver-site counts above).
